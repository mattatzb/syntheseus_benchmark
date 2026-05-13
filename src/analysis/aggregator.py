#!/usr/bin/env python3
"""
Benchmark results aggregation and CSV output.

This module handles CSV I/O operations including reading, writing, and merging
benchmark results across multiple methods and runs.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

from .utils import format_float
from .route_processor import analyze_method
from config import Config
from .extractor import extract_data



# Column configuration constants
INVENTORY_COLUMNS = ["n_inventory_excluded", "inventory_true_pct"]
DEBUG_COLUMNS = ["n_targets", "n_solved", "note"]
EXCLUDED_METHOD_DIRS = {"analysis", "configs"}


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
