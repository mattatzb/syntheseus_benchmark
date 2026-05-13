#!/usr/bin/env python3
"""Filter out molecules whose stats.json marks target_in_inventory as true."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def collect_smiles_in_inventory(stats_root: Path) -> set[str]:
    """Collect target SMILES where target_in_inventory is true."""
    smiles_to_drop: set[str] = set()

    for stats_path in stats_root.glob("*/stats.json"):
        with stats_path.open("r", encoding="utf-8") as f:
            stats = json.load(f)
        if stats.get("target_in_inventory") is True and isinstance(stats.get("smiles"), str):
            smiles_to_drop.add(stats["smiles"])

    return smiles_to_drop


def filter_csv(input_csv: Path, output_csv: Path, smiles_to_drop: set[str]) -> tuple[int, int]:
    """Write filtered CSV and return (total_rows, kept_rows)."""
    with input_csv.open("r", newline="", encoding="utf-8") as fin:
        reader = csv.DictReader(fin)
        if not reader.fieldnames or "smiles" not in reader.fieldnames:
            raise ValueError(f"Missing required 'smiles' column in {input_csv}")

        rows = list(reader)

    kept_rows = [row for row in rows if row.get("smiles") not in smiles_to_drop]

    with output_csv.open("w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(kept_rows)

    return len(rows), len(kept_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a CSV copy excluding molecules with target_in_inventory=true."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("/home/tatzber/syntheseus_benchmark/pesticides_qed_ge_0_5.csv"),
        help="Input CSV path containing a 'smiles' column.",
    )
    parser.add_argument(
        "--stats-root",
        type=Path,
        default=Path(
            "/home/tatzber/syntheseus_benchmark/results_6_runs/LocalRetro_mcts/"
            "LocalRetro_2026-02-25T17:46:14"
        ),
        help="Directory containing numbered subfolders with stats.json files.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path(
            "/home/tatzber/syntheseus_benchmark/pesticides_qed_ge_0_5_no_inventory.csv"
        ),
        help="Output filtered CSV path.",
    )
    args = parser.parse_args()

    smiles_to_drop = collect_smiles_in_inventory(args.stats_root)
    total_rows, kept_rows = filter_csv(args.input_csv, args.output_csv, smiles_to_drop)
    dropped_rows = total_rows - kept_rows

    print(f"Read {total_rows} molecules from {args.input_csv}")
    print(
        f"Found {len(smiles_to_drop)} unique target SMILES with target_in_inventory=true "
        f"in {args.stats_root}"
    )
    print(f"Wrote {kept_rows} molecules to {args.output_csv} (dropped {dropped_rows})")


if __name__ == "__main__":
    main()
