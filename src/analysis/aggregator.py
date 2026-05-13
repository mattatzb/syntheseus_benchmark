#!/usr/bin/env python3
"""
Benchmark results aggregation and CSV output.

This module handles CSV I/O operations including reading, writing, and merging
benchmark results across multiple methods and runs.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter
from pathlib import Path
from typing import Any

from .utils import format_float
from .route_processor import analyze_method
from config import Config
from .mcts_extractor import extract_mcts_data
from .retro_star_extractor import extract_data as extract_retro_star_data
from .plots import plot_distribution_csvs



# Column configuration constants
INVENTORY_COLUMNS = ["n_inventory_excluded", "inventory_true_pct"]
DEBUG_COLUMNS = ["n_targets", "n_solved", "note"]
EXCLUDED_METHOD_DIRS = {"analysis", "configs"}
DISTRIBUTION_FIELDNAMES = ["method", "run_dir", "metric", "value", "count", "proportion"]
DISTRIBUTION_SPECS = {
    "route_lengths": ("route_lengths_distribution.csv", "route_length"),
    "route_building_block_counts": (
        "building_blocks_distribution.csv",
        "n_building_blocks",
    ),
    "reaction_reactant_counts": ("reactants_distribution.csv", "n_reactants"),
}


def extractor_for_method(method_dir: Path):
    """Select the route pickle extractor matching a method directory."""
    if method_dir.name.endswith("_mcts"):
        return extract_mcts_data
    return extract_retro_star_data


def build_fieldnames(show_inventory_columns: bool, show_debug_columns: bool) -> list[str]:
    """
    Build the CSV column headers based on output configuration.

    Args:
        show_inventory_columns: Include inventory-related columns.
        show_debug_columns: Include debug columns (n_targets, n_solved, note).

    Returns:
        List of column header names in output order.
    """
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
    """
    Format a single row for CSV output, converting floats to strings.

    Args:
        row: Dictionary containing raw metric values.
        show_inventory_columns: Include inventory columns in output.
        show_debug_columns: Include debug columns in output.

    Returns:
        Formatted row dictionary ready for CSV writing.
    """
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
    """
    Write aggregated metrics to a CSV file.

    Creates parent directories if needed and writes rows with appropriate formatting.

    Args:
        rows: List of metric dictionaries to write.
        output_csv: Path to output CSV file.
        show_inventory_columns: Include inventory columns.
        show_debug_columns: Include debug columns.
    """
    fieldnames = build_fieldnames(show_inventory_columns, show_debug_columns)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(format_row(row, show_inventory_columns, show_debug_columns))


def _integer_distribution_values(values: object) -> list[int]:
    """
    Convert a raw distribution list into integer bin values.

    Non-finite values and non-integer numeric values are ignored.
    """
    if not isinstance(values, list):
        return []

    integer_values: list[int] = []
    for value in values:
        if isinstance(value, bool):
            continue
        if isinstance(value, int):
            integer_values.append(value)
            continue
        if isinstance(value, float) and math.isfinite(value) and value.is_integer():
            integer_values.append(int(value))
    return integer_values


def build_distribution_rows(
    rows: list[dict[str, object]],
    metric_key: str,
    metric_name: str,
) -> list[dict[str, object]]:
    """
    Build value/count/proportion rows for one distribution metric.

    Args:
        rows: Analyzed method rows returned by analyze_method.
        metric_key: Key containing the raw list values in each analyzed row.
        metric_name: Human-readable metric name written to the CSV.

    Returns:
        Rows with method, run_dir, metric, value, count, and proportion columns.
    """
    distribution_rows: list[dict[str, object]] = []

    for row in rows:
        values = _integer_distribution_values(row.get(metric_key))
        total = len(values)
        if total == 0:
            continue

        counts = Counter(values)
        for value, count in sorted(counts.items()):
            distribution_rows.append(
                {
                    "method": row["method"],
                    "run_dir": row["run_dir"],
                    "metric": metric_name,
                    "value": value,
                    "count": count,
                    "proportion": format_float(count / total),
                }
            )

    return distribution_rows


def write_distribution_csv(rows: list[dict[str, object]], output_csv: Path) -> None:
    """
    Write one distribution CSV.

    Args:
        rows: Distribution rows from build_distribution_rows.
        output_csv: Path to the output CSV.
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=DISTRIBUTION_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def write_plot_data_csvs(
    analyzed_rows: list[dict[str, object]], output_dir: Path
) -> dict[str, Path]:
    """
    Write distribution CSVs used by downstream plotting.

    This consumes the existing analyze_method output and does not trigger any
    additional route parsing.

    Args:
        analyzed_rows: Rows returned by analyze_method.
        output_dir: Directory where plot-data CSVs should be written.

    Returns:
        Mapping from metric key to written CSV path.
    """
    written_paths: dict[str, Path] = {}
    for metric_key, (filename, metric_name) in DISTRIBUTION_SPECS.items():
        distribution_rows = build_distribution_rows(
            analyzed_rows,
            metric_key=metric_key,
            metric_name=metric_name,
        )
        output_csv = output_dir / filename
        write_distribution_csv(distribution_rows, output_csv)
        written_paths[metric_key] = output_csv
    return written_paths


def read_existing_csv(input_csv: Path) -> list[dict[str, object]]:
    """
    Read an existing CSV file into a list of dictionaries.

    Args:
        input_csv: Path to CSV file to read.

    Returns:
        List of dictionaries (one per row), or empty list if file doesn't exist.
    """
    if not input_csv.exists():
        return []

    with input_csv.open("r", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def merge_rows(
    existing_rows: list[dict[str, object]], new_rows: list[dict[str, object]]
) -> list[dict[str, object]]:
    """
    Merge new analysis rows into existing CSV data.

    Replaces rows with matching method names or appends new methods.
    Preserves the order of existing rows.

    Args:
        existing_rows: Rows from the existing CSV.
        new_rows: Newly analyzed rows.

    Returns:
        Merged list of rows with updated data for matched methods.
    """
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
        default=Config.RESULTS_ROOT,
        help="Root directory containing method folders.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Config.DEFAULT_OUTPUT_CSV,
        help="Path to output CSV.",
    )
    parser.add_argument(
        "--write-plot-data",
        action="store_true",
        help="Write additional distribution CSVs for plotting without changing the summary CSV.",
    )
    parser.add_argument(
        "--plot-distributions",
        action="store_true",
        help=(
            "Write distribution CSVs and create the route distribution PDF using "
            "the same analyzed data."
        ),
    )
    parser.add_argument(
        "--plot-data-dir",
        type=Path,
        default=Config.OUTPUT_ROOT,
        help="Directory for plot-data CSVs.",
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
    """
    Main entry point for benchmark aggregation.

    Loads extraction function, analyzes method directories (or specific runs),
    optionally merges with existing results, and writes aggregated CSV.
    """
    args = parse_args()

    if args.run_dirs:
        run_dirs = sorted(Path(run_dir).resolve() for run_dir in args.run_dirs)
        analyzed_rows = [
            analyze_method(
                run_dir.parent,
                extractor_for_method(run_dir.parent),
                run_dir=run_dir,
            )
            for run_dir in run_dirs
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

        analyzed_rows = [
            analyze_method(method_dir, extractor_for_method(method_dir))
            for method_dir in method_dirs
        ]

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
    if args.write_plot_data or args.plot_distributions:
        written_paths = write_plot_data_csvs(analyzed_rows, args.plot_data_dir)
        for path in written_paths.values():
            print(f"Wrote plot data to {path}")
        if args.plot_distributions:
            plot_distribution_csvs(args.plot_data_dir, Config.DEFAULT_OUTPUT_PLOT)
            print(f"Wrote plots to {Config.DEFAULT_OUTPUT_PLOT}")

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
