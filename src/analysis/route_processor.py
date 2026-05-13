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
from pathlib import Path
from typing import Callable

from .utils import (
    canonicalize_smiles,
    extract_reactant_smiles,
    is_finite_number,
    mean_or_nan,
    shorten_notes,
)


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

    return reactants_counts or inferred_reactants_counts


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
    run_dir: Path | None = None,
) -> dict[str, object]:
    """
    Analyze all targets in a method directory and extract aggregated metrics.

    Processes target directories in the latest run (or specified run), computing metrics
    including solving rate, wall-clock time, route length, and route count.

    Args:
        method_dir: Path to the method directory (e.g., "results_6_runs/MEGAN_retro_star").
        extract_data: Function that deserializes pickle routes to JSON.
        run_dir: Specific run directory to analyze. If None, uses latest run.

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
    """

    all_run_dirs = sorted(
        d for d in method_dir.iterdir() if d.is_dir() and d.name.isdigit()
    )
    if not all_run_dirs:
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
            "note": "no run directories",
        }

    selected_run_dir = run_dir or all_run_dirs[-1]
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

    for target_dir in all_run_dirs:
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
            try:
                route_length = route_length_from_pkl(route_file, extract_data)
                reaction_reactant_counts.extend(
                    reactants_count_from_pkl(route_file, extract_data)
                )
                route_building_block_counts.append(
                    building_blocks_count_from_pkl(route_file, extract_data)
                )
            except Exception as exc:  # noqa: BLE001
                notes.append(f"{target_dir.name}/{route_file.name}: parse failed ({exc})")
                continue
            if is_finite_number(route_length):
                route_lengths.append(float(route_length))
            else:
                notes.append(f"{target_dir.name}/{route_file.name}: no valid route length data")

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
