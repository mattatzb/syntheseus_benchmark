#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Callable, Iterator

EXTRACTOR_PATH = Path(
    "/home/tatzber/saturn/oracles/synthesizability/utils/extract_syntheseus_route_data.py"
)


def load_extract_data_func() -> Callable[[str], str]:
    spec = importlib.util.spec_from_file_location("extract_syntheseus_route_data", EXTRACTOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load extractor from {EXTRACTOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    extract_data = getattr(module, "extract_data", None)
    if not callable(extract_data):
        raise RuntimeError(f"'extract_data' not found in {EXTRACTOR_PATH}")
    return extract_data


def iter_target_dirs(run_dir: Path) -> Iterator[Path]:
    def sort_key(path: Path) -> tuple[int, int | str]:
        return (0, int(path.name)) if path.name.isdigit() else (1, path.name)

    for target_dir in sorted(run_dir.iterdir(), key=sort_key):
        if target_dir.is_dir() and target_dir.name.isdigit():
            yield target_dir


def has_route_file(target_dir: Path) -> bool:
    return any(target_dir.glob("route_*.pkl"))


def get_primary_route_file(target_dir: Path) -> Path | None:
    route_files = sorted(target_dir.glob("route_*.pkl"))
    return route_files[0] if route_files else None


def find_solved(run_dir: Path) -> dict[str, str]:
    solved_dict: dict[str, str] = {}
    for target_dir in iter_target_dirs(run_dir):
        solved_dict[target_dir.name] = "solved" if has_route_file(target_dir) else "unsolved"
    return solved_dict


def find_common_solved(run_dir1: Path, run_dir2: Path) -> list[str]:
    solved1 = find_solved(run_dir1)
    solved2 = find_solved(run_dir2)
    return sorted(
        [
            target
            for target, status in solved1.items()
            if target in solved2 and status == "solved" and solved2[target] == "solved"
        ],
        key=int,
    )


def find_first_only_solved(run_dir1: Path, run_dir2: Path) -> list[str]:
    solved1 = find_solved(run_dir1)
    solved2 = find_solved(run_dir2)
    return sorted(
        [
            target
            for target, status in solved1.items()
            if status == "solved" and solved2.get(target) != "solved"
        ],
        key=int,
    )


def find_second_only_solved(run_dir1: Path, run_dir2: Path) -> list[str]:
    solved1 = find_solved(run_dir1)
    solved2 = find_solved(run_dir2)
    return sorted(
        [
            target
            for target, status in solved2.items()
            if status == "solved" and solved1.get(target) != "solved"
        ],
        key=int,
    )


def load_stats(target_dir: Path) -> dict[str, object]:
    stats_path = target_dir / "stats.json"
    if not stats_path.exists():
        return {}
    with stats_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def get_route_info(route_pkl: Path, extractor: Callable[[str], str]) -> dict[str, object]:
    return json.loads(extractor(str(route_pkl)))


def get_reaction_smiles_from_route(route_pkl: Path, extractor: Callable[[str], str]) -> list[str]:
    route_info = get_route_info(route_pkl, extractor)
    reaction_nodes = [
        node_info
        for node_info in route_info.values()
        if isinstance(node_info, dict) and bool(node_info.get("is_rxn", False))
    ]
    reaction_nodes.sort(key=lambda node_info: node_info.get("depth", 0))
    return [
        str(node_info["rxn_smiles"])
        for node_info in reaction_nodes
        if node_info.get("rxn_smiles")
    ]


def compare_common_solved(
    run_dir1: Path,
    run_dir2: Path,
    extractor: Callable[[str], str],
) -> list[dict[str, object]]:
    comparisons: list[dict[str, object]] = []
    for target in find_common_solved(run_dir1, run_dir2):
        target_dir1 = run_dir1 / target
        target_dir2 = run_dir2 / target
        route_file1 = get_primary_route_file(target_dir1)
        route_file2 = get_primary_route_file(target_dir2)
        if route_file1 is None or route_file2 is None:
            continue

        rxns1 = get_reaction_smiles_from_route(route_file1, extractor)
        rxns2 = get_reaction_smiles_from_route(route_file2, extractor)
        if rxns1 == rxns2:
            continue

        stats1 = load_stats(target_dir1)
        stats2 = load_stats(target_dir2)
        comparisons.append(
            {
                "target": target,
                "index": stats1.get("index", stats2.get("index", target)),
                "smiles": stats1.get("smiles", stats2.get("smiles", "")),
                "run1_rxn_smiles": rxns1,
                "run2_rxn_smiles": rxns2,
                "only_in_run1": [rxn for rxn in rxns1 if rxn not in rxns2],
                "only_in_run2": [rxn for rxn in rxns2 if rxn not in rxns1],
            }
        )
    return comparisons


def get_first_only_solved_reactions(
    run_dir1: Path,
    run_dir2: Path,
    extractor: Callable[[str], str],
) -> list[dict[str, object]]:
    first_only: list[dict[str, object]] = []
    for target in find_first_only_solved(run_dir1, run_dir2):
        target_dir = run_dir1 / target
        route_file = get_primary_route_file(target_dir)
        if route_file is None:
            continue

        stats = load_stats(target_dir)
        first_only.append(
            {
                "target": target,
                "index": stats.get("index", target),
                "smiles": stats.get("smiles", ""),
                "run1_rxn_smiles": get_reaction_smiles_from_route(route_file, extractor),
            }
        )
    return first_only


def get_second_only_solved_reactions(
    run_dir1: Path,
    run_dir2: Path,
    extractor: Callable[[str], str],
) -> list[dict[str, object]]:
    second_only: list[dict[str, object]] = []
    for target in find_second_only_solved(run_dir1, run_dir2):
        target_dir = run_dir2 / target
        route_file = get_primary_route_file(target_dir)
        if route_file is None:
            continue

        stats = load_stats(target_dir)
        second_only.append(
            {
                "target": target,
                "index": stats.get("index", target),
                "smiles": stats.get("smiles", ""),
                "run2_rxn_smiles": get_reaction_smiles_from_route(route_file, extractor),
            }
        )
    return second_only


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare reaction SMILES between two Syntheseus run directories. "
            "The script uses route_0.pkl for each solved target."
        )
    )
    parser.add_argument("run_dir1", type=Path, help="Reference run directory.")
    parser.add_argument("run_dir2", type=Path, help="Comparison run directory.")
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path(__file__).resolve().parent / "path_comparison_results.json",
        help="Path to write the comparison results JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    extractor = load_extract_data_func()

    common_differences = compare_common_solved(args.run_dir1, args.run_dir2, extractor)
    first_only_solved = get_first_only_solved_reactions(args.run_dir1, args.run_dir2, extractor)
    second_only_solved = get_second_only_solved_reactions(args.run_dir1, args.run_dir2, extractor)

    result = {
        "run_dir1": str(args.run_dir1),
        "run_dir2": str(args.run_dir2),
        "n_common_solved": len(find_common_solved(args.run_dir1, args.run_dir2)),
        "n_common_solved_different_paths": len(common_differences),
        "n_solved_only_in_run1": len(first_only_solved),
        "n_solved_only_in_run2": len(second_only_solved),
        "common_solved_different_paths": common_differences,
        "solved_only_in_run1": first_only_solved,
        "solved_only_in_run2": second_only_solved,
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(f"Wrote comparison results to {args.output_json}")
    print(f"Common solved targets: {result['n_common_solved']}")
    print(f"Common solved with different paths: {result['n_common_solved_different_paths']}")
    print(f"Solved only in run1: {result['n_solved_only_in_run1']}")
    print(f"Solved only in run2: {result['n_solved_only_in_run2']}")


if __name__ == "__main__":
    main()
