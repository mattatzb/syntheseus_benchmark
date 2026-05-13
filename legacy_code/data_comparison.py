#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from statistics import mean
from typing import Callable, Iterator

EXCLUDED_METHOD_DIRS = {"analysis", "configs"}
INVENTORY_COLUMNS = ["n_inventory_excluded", "inventory_true_pct"]
DEBUG_COLUMNS = ["n_targets", "n_solved", "note"]
EXTRACTOR_PATH = Path(
    "/home/tatzber/saturn/oracles/synthesizability/utils/extract_syntheseus_route_data.py"
)


def is_finite_number(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value)


def mean_or_nan(values: list[float]) -> float:
    return mean(values) if values else math.nan


def format_float(value: object, digits: int = 6) -> str:
    return f"{float(value):.{digits}f}" if is_finite_number(value) else ""


def find_latest_run_dir(method_dir: Path) -> Path | None:
    run_dirs = sorted(d for d in method_dir.iterdir() if d.is_dir())
    return run_dirs[-1] if run_dirs else None


def iter_target_dirs(run_dir: Path) -> Iterator[Path]:
    for target_dir in sorted(run_dir.iterdir(), key=lambda p: p.name):
        if target_dir.is_dir() and target_dir.name.isdigit():
            yield target_dir


def shorten_notes(notes: list[str], limit: int = 5) -> str:
    if not notes:
        return ""
    short = "; ".join(notes[:limit])
    if len(notes) > limit:
        short += f"; ... +{len(notes) - limit} more"
    return short


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


def route_length_from_pkl(route_pkl: Path, extract_data: Callable[[str], str]) -> float:
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


def analyze_method(
    method_dir: Path,
    extract_data: Callable[[str], str],
    run_dir: Path | None = None,
) -> dict[str, object]:
    selected_run_dir = run_dir or find_latest_run_dir(method_dir)
    if selected_run_dir is None:
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
            "note": "no run directory",
        }

    n_total_targets = 0
    n_targets = 0
    n_inventory_excluded = 0
    n_solved = 0
    wall_times: list[float] = []
    route_lengths: list[float] = []  # One entry per extracted route file.
    route_counts: list[float] = []
    notes: list[str] = []

    for target_dir in iter_target_dirs(selected_run_dir):
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
            except Exception as exc:  # noqa: BLE001
                notes.append(f"{target_dir.name}/{route_file.name}: parse failed ({exc})")
                continue
            if is_finite_number(route_length):
                route_lengths.append(float(route_length))
            else:
                notes.append(f"{target_dir.name}/{route_file.name}: no valid route length data")

    return {
        "method": method_dir.name,
        "run_dir": selected_run_dir.name,
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
        "note": shorten_notes(notes),
    }


def build_fieldnames(show_inventory_columns: bool, show_debug_columns: bool) -> list[str]:
    fieldnames = [
        "method",
        "run_dir",
        "n_total_targets",
        "solving_rate",
        "mean_wall_time_s",
        "avg_route_length",
        "avg_num_routes",
    ]
    if show_inventory_columns:
        fieldnames = fieldnames[:4] + INVENTORY_COLUMNS + fieldnames[4:]
    if show_debug_columns:
        fieldnames.extend(DEBUG_COLUMNS)
    return fieldnames


def format_row(
    row: dict[str, object], show_inventory_columns: bool, show_debug_columns: bool
) -> dict[str, object]:
    out: dict[str, object] = {
        "method": row["method"],
        "run_dir": row["run_dir"],
        "n_total_targets": row["n_total_targets"],
        "solving_rate": format_float(row["solving_rate"]),
        "mean_wall_time_s": format_float(row["mean_wall_time_s"]),
        "avg_route_length": format_float(row["avg_route_length"]),
        "avg_num_routes": format_float(row["avg_num_routes"]),
    }
    if show_inventory_columns:
        out["n_inventory_excluded"] = row["n_inventory_excluded"]
        out["inventory_true_pct"] = format_float(row["inventory_true_pct"], 2)
    if show_debug_columns:
        out["n_targets"] = row["n_targets"]
        out["n_solved"] = row["n_solved"]
        out["note"] = row["note"]
    return out


def write_csv(
    rows: list[dict[str, object]],
    output_csv: Path,
    show_inventory_columns: bool,
    show_debug_columns: bool,
) -> None:
    fieldnames = build_fieldnames(show_inventory_columns, show_debug_columns)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(format_row(row, show_inventory_columns, show_debug_columns))


def read_existing_csv(input_csv: Path) -> list[dict[str, object]]:
    if not input_csv.exists():
        return []

    with input_csv.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def merge_rows(
    existing_rows: list[dict[str, object]], new_rows: list[dict[str, object]]
) -> list[dict[str, object]]:
    new_rows_by_method = {str(row["method"]): row for row in new_rows}
    merged_rows: list[dict[str, object]] = []
    seen_methods: set[str] = set()

    for row in existing_rows:
        method = str(row.get("method", ""))
        if method in new_rows_by_method:
            merged_rows.append(new_rows_by_method[method])
            seen_methods.add(method)
        else:
            merged_rows.append(row)

    for method, row in sorted(new_rows_by_method.items()):
        if method not in seen_methods:
            merged_rows.append(row)

    return merged_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate benchmark metrics across methods into a CSV."
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Root directory containing method folders (default: results_6_runs).",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path(__file__).resolve().parent / "results_summary.csv",
        help="Path to output CSV.",
    )
    parser.add_argument(
        "--show-inventory-columns",
        action="store_true",
        help="Include inventory-related columns in CSV output.",
    )
    parser.add_argument(
        "--show-debug-columns",
        action="store_true",
        help="Include n_targets, n_solved, and note columns in CSV output.",
    )
    parser.add_argument(
        "--method-dir",
        dest="method_dirs",
        type=Path,
        action="append",
        default=[],
        help=(
            "Specific method directory to analyze. Can be passed multiple times. "
            "If omitted, all method directories under --results-root are analyzed."
        ),
    )
    parser.add_argument(
        "--merge-with-existing",
        action="store_true",
        help="Merge analyzed rows into the existing output CSV instead of replacing all rows.",
    )
    parser.add_argument(
        "--run-dir",
        dest="run_dirs",
        type=Path,
        action="append",
        default=[],
        help=(
            "Specific run directory to analyze. Can be passed multiple times. "
            "The method name is taken from the parent directory of the run directory."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    extract_data = load_extract_data_func()
    if args.run_dirs:
        run_dirs = sorted(Path(run_dir).resolve() for run_dir in args.run_dirs)
        analyzed_rows = [
            analyze_method(run_dir.parent, extract_data, run_dir=run_dir) for run_dir in run_dirs
        ]
    else:
        if args.method_dirs:
            method_dirs = sorted(Path(method_dir).resolve() for method_dir in args.method_dirs)
        else:
            method_dirs = sorted(
                d
                for d in args.results_root.iterdir()
                if d.is_dir() and d.name not in EXCLUDED_METHOD_DIRS
            )

        analyzed_rows = [analyze_method(method_dir, extract_data) for method_dir in method_dirs]
    rows = analyzed_rows
    if args.merge_with_existing:
        existing_rows = read_existing_csv(args.output_csv)
        rows = merge_rows(existing_rows, analyzed_rows)

    write_csv(
        rows,
        args.output_csv,
        show_inventory_columns=args.show_inventory_columns,
        show_debug_columns=args.show_debug_columns,
    )

    print(f"Wrote {len(rows)} rows to {args.output_csv}")
    for row in analyzed_rows:
        line = (
            f"{row['method']}: solved {row['n_solved']}/{row['n_targets']} "
            f"(rate={format_float(row['solving_rate'], 3)}), "
            f"wall={format_float(row['mean_wall_time_s'], 3)} s, "
            f"route_len={format_float(row['avg_route_length'], 3)}, "
            f"n_routes={format_float(row['avg_num_routes'], 3)}"
        )
        if args.show_inventory_columns:
            line += (
                f", inv_excluded={row['n_inventory_excluded']}, "
                f"inv_true_pct={format_float(row['inventory_true_pct'], 2)}%"
            )
        print(line)


if __name__ == "__main__":
    main()
