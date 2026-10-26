#!/usr/bin/env python3
"""
Route processing and metric extraction.

This module handles the analysis of individual routes and target directories,
extracting metrics such as route length, solving status, and aggregated statistics
for benchmark evaluation.
"""

from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Callable

from .utils import (
    canonicalize_smiles,
    extract_reactant_smiles,
    is_finite_number,
    mean_or_nan,
    shorten_notes,
)


class UnreliableReactantCountsError(ValueError):
    """An MCTS state sequence cannot supply a complete reactant distribution."""


def route_length_from_pkl(route_pkl: Path, extract_data: Callable[[str], str]) -> float:
    """
    Extract route length from a pickled route file.

    Analyzes the route structure to determine the number of steps (reactions or transitions).
    Supports two formats:
    - Reaction-based: counts the number of is_rxn=True nodes
    - State-based: counts state nodes minus the root (used by MCTS)

    Args:
        route_pkl: Path to the route pickle file.
        extract_data: Function that deserializes pickle to JSON dict.

    Returns:
        Route length as float, or math.nan if extraction fails.
    """
    raw = extract_data(str(route_pkl))
    data = json.loads(raw)
    node_infos = [node_info for node_info in data.values() if isinstance(node_info, dict)]
    if not node_infos:
        return math.nan

    rxn_count = sum(1 for node_info in node_infos if bool(node_info.get("is_rxn", False)))
    if rxn_count:
        return float(rxn_count)

    # Some route formats store only state nodes (for example MolSetNode chains in MCTS).
    # In that case, the number of transitions is the number of nodes minus the root node.
    state_node_count = sum(
        1
        for node_info in node_infos
        if not bool(node_info.get("is_mol", False)) and not bool(node_info.get("is_rxn", False))
    )
    if state_node_count == len(node_infos) and state_node_count > 1:
        return float(state_node_count - 1)

    return math.nan


def reactants_count_from_pkl(route_pkl: Path, extract_data: Callable[[str], str]) -> list[int]:
    """
    Extract the number of reactants for each reaction step in a route pickle.

    Args:
        route_pkl: Path to the route pickle file.
        extract_data: Function that deserializes pickle routes to JSON.

    Returns:
        One integer per reaction node, in extracted route order.

    Raises:
        UnreliableReactantCountsError: Inferred MCTS transitions are incomplete
            or inconsistent. The caller must omit the method's reactant
            distribution rather than normalize a partial set of counts.
    """
    raw = extract_data(str(route_pkl))
    data = json.loads(raw)
    reactants_counts: list[int] = []
    inferred_reactants_counts: list[int] = []

    for node_info in data.values():
        if not isinstance(node_info, dict):
            continue
        inferred_count = node_info.get("inferred_reactant_count")
        if isinstance(inferred_count, int) and inferred_count > 0:
            inferred_reactants_counts.append(inferred_count)
        if not bool(node_info.get("is_rxn", False)):
            continue

        rxn_smiles = node_info.get("rxn_smiles")
        if not isinstance(rxn_smiles, str):
            rxn_smiles = None
        reactants_counts.append(len(extract_reactant_smiles(rxn_smiles)))

    if reactants_counts:
        return reactants_counts

    state_nodes = [
        node_info for node_info in data.values()
        if isinstance(node_info, dict) and bool(node_info.get("is_state", False))
    ]
    for index, node_info in enumerate(state_nodes):
        if node_info.get("depth") != index:
            raise UnreliableReactantCountsError(
                f"cannot infer reactant counts: state {index} has nonconsecutive depth"
            )
        if index == 0:
            continue
        inferred_count = node_info.get("inferred_reactant_count")
        expanded = node_info.get("inferred_expanded_mol_smiles", [])
        if (
            not isinstance(inferred_count, int)
            or isinstance(inferred_count, bool)
            or inferred_count <= 0
            or len(expanded) != 1
        ):
            raise UnreliableReactantCountsError(
                f"cannot infer reactant counts at state {index}: expected one "
                "expanded molecule and at least one newly introduced reactant"
            )

    return inferred_reactants_counts


def building_blocks_count_from_pkl(route_pkl: Path, extract_data: Callable[[str], str]) -> int:
    """
    Count unique purchasable building blocks in a route pickle.

    Args:
        route_pkl: Path to the route pickle file.
        extract_data: Function that deserializes pickle routes to JSON.

    Returns:
        Number of unique canonical molecule SMILES marked as purchasable.
    """
    raw = extract_data(str(route_pkl))
    data = json.loads(raw)
    building_blocks: set[str] = set()

    for node_info in data.values():
        if not isinstance(node_info, dict):
            continue

        if bool(node_info.get("is_mol", False)) and bool(node_info.get("is_purchasable", False)):
            mol_smiles = node_info.get("mol_smiles")
            if isinstance(mol_smiles, str):
                canonical = canonicalize_smiles(mol_smiles)
                if canonical is not None:
                    building_blocks.add(canonical)

        for mol_info in node_info.get("mols", []):
            if not isinstance(mol_info, dict):
                continue
            if not bool(mol_info.get("is_purchasable", False)):
                continue
            mol_smiles = mol_info.get("smiles")
            if not isinstance(mol_smiles, str):
                continue
            canonical = canonicalize_smiles(mol_smiles)
            if canonical is not None:
                building_blocks.add(canonical)

    return len(building_blocks)


def analyze_method(
    method_dir: Path,
    extract_data: Callable[[str], str],
) -> dict[str, object]:
    """
    Analyze all targets in a method directory and extract aggregated metrics.

    Processes all numeric target directories under the method, computing metrics
    including solving rate, wall-clock time, route length, and route count.

    Args:
        method_dir: Path to a method directory containing numeric target directories.
        extract_data: Function that deserializes pickle routes to JSON.

    Returns:
        Dictionary with aggregated metrics:
            - method: Method name
            - run_dir: Method repository name (the parent subrepository under results root)
            - n_total_targets: Total number of target molecules in run
            - n_targets: Number of targets (excluding inventory)
            - n_inventory_excluded: Count of inventory-already-owned targets
            - inventory_true_pct: Percentage of inventory-excluded targets
            - n_solved: Number of successfully solved targets
            - solving_rate: n_solved / n_targets
            - mean_wall_time_s: Average wall-clock time for solved targets
            - avg_route_length: Average number of reactions per route
            - avg_num_routes: Average number of routes found per target
            - route_lengths: Number of reactions for each route
            - reaction_reactant_counts: Number of reactants for each reaction step
            - route_building_block_counts: Number of unique building blocks for each route
            - note: Concatenated error/warning messages (up to 5)

    Raises:
        RuntimeError: A route cannot be read or a required metric cannot be
            extracted. Includes the file path and preserves the original cause.
            Unreliable MCTS reactant inference instead emits a RuntimeWarning
            and omits that method's entire reactant distribution.
    """

    target_dirs = sorted(
        d for d in method_dir.iterdir() if d.is_dir() and d.name.isdigit()
    )
    if not target_dirs:
        return {
            "method": method_dir.name,
            "run_dir": "",
            "n_total_targets": 0,
            "n_targets": 0,
            "n_inventory_excluded": 0,
            "inventory_true_pct": math.nan,
            "n_solved": 0,
            "solving_rate": math.nan,
            "mean_wall_time_s": math.nan,
            "avg_route_length": math.nan,
            "avg_num_routes": math.nan,
            "route_lengths": [],
            "reaction_reactant_counts": [],
            "route_building_block_counts": [],
            "note": "no target directories",
        }

    repository_name = method_dir.parent.name

    n_total_targets = 0
    n_targets = 0
    n_inventory_excluded = 0
    n_solved = 0
    wall_times: list[float] = []
    route_lengths: list[float] = []  
    route_counts: list[float] = []
    reaction_reactant_counts: list[int] = []
    route_building_block_counts: list[int] = []
    notes: list[str] = []
    reactant_issues: list[str] = []

    for target_dir in target_dirs:
        stats_path = target_dir / "stats.json"
        if not stats_path.exists():
            continue

        n_total_targets += 1
        with stats_path.open("r", encoding="utf-8") as f:
            stats = json.load(f)

        if bool(stats.get("target_in_inventory", False)):
            n_inventory_excluded += 1
            continue

        n_targets += 1
        wall = stats.get("soln_time_wallclock")
        solved = is_finite_number(wall)
        if not solved:
            continue

        n_solved += 1
        wall_times.append(float(wall))
        route_files = sorted(target_dir.glob("route_*.pkl"))
        route_counts.append(float(len(route_files)))
        for route_file in route_files:
            metric = "route length"
            try:
                route_length = route_length_from_pkl(route_file, extract_data)
                if not is_finite_number(route_length):
                    raise ValueError("no valid route length data")
                metric = "reactant counts"
                try:
                    reactant_counts = reactants_count_from_pkl(route_file, extract_data)
                except UnreliableReactantCountsError as exc:
                    reactant_issues.append(f"{route_file}: {exc}")
                    reactant_counts = []
                metric = "building-block count"
                building_block_count = building_blocks_count_from_pkl(route_file, extract_data)
            except Exception as exc:  # noqa: BLE001
                raise RuntimeError(
                    f"Failed to extract {metric} from {route_file}: {exc}"
                ) from exc
            route_lengths.append(float(route_length))
            reaction_reactant_counts.extend(reactant_counts)
            route_building_block_counts.append(building_block_count)

    if reactant_issues:
        # Mixing valid routes with silently dropped transitions biases the
        # normalized distribution. Keep the other metrics, but omit this entire
        # method's reactant distribution and report an example of the problem.
        reaction_reactant_counts = []
        message = (
            f"reactant counts omitted for {method_dir.name}: "
            f"{len(reactant_issues)} route(s) have unreliable MCTS transitions; "
            f"first issue: {reactant_issues[0]}"
        )
        notes.append(message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)

    return {
        "method": method_dir.name,
        "run_dir": repository_name,
        "n_total_targets": n_total_targets,
        "n_targets": n_targets,
        "n_inventory_excluded": n_inventory_excluded,
        "inventory_true_pct": (
            (100.0 * n_inventory_excluded / n_total_targets) if n_total_targets else math.nan
        ),
        "n_solved": n_solved,
        "solving_rate": (n_solved / n_targets) if n_targets else math.nan,
        "mean_wall_time_s": mean_or_nan(wall_times),
        "avg_route_length": mean_or_nan(route_lengths),
        "avg_num_routes": mean_or_nan(route_counts),
        "route_lengths": route_lengths,
        "reaction_reactant_counts": reaction_reactant_counts,
        "route_building_block_counts": route_building_block_counts,
        "note": shorten_notes(notes),
    }
