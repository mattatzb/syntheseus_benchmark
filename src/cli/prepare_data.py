#!/usr/bin/env python3
"""Prepare the benchmark tables and targets from the archived pesticide records."""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from data_prep.filter_pesticides import filter_by_qed
from data_prep.filter_pesticides_by_inventory import clean_inventory_centroids


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--raw-csv", type=Path, default=ROOT / "data/pesticides_raw.csv",
                        help="Archived raw CSV; defaults to the supplied snapshot.")
    source.add_argument("--regenerate-raw", action="store_true",
                        help="Fetch live databases. Not recommended for reproducing the paper.")
    parser.add_argument("--inventory-smiles", type=Path,
                        default=ROOT / "data/enamine_sep_lpdc_blocks.smi")
    parser.add_argument("--output-dir", type=Path,
                        help="Defaults to data/, or data/refreshed/ for live regeneration.")
    args = parser.parse_args(argv)
    if np.__version__ != "1.26.4":
        raise ValueError(
            f"Preparation requires NumPy 1.26.4 for historical cluster ordering; "
            f"found {np.__version__}. Update the environment from environment.yml."
        )
    output_dir = args.output_dir or ROOT / ("data/refreshed" if args.regenerate_raw else "data")
    if args.regenerate_raw and output_dir.resolve() == (ROOT / "data").resolve():
        raise ValueError("Use another --output-dir for live data to preserve the archived raw CSV")
    inventory = {line.strip() for line in args.inventory_smiles.read_text().splitlines() if line.strip()}
    if args.regenerate_raw:
        warnings.warn(
            "Live regeneration uses changing external databases and is not recommended "
            "for reproducing the published benchmark.", UserWarning,
        )
        # Import live-service clients only for the explicitly requested refresh.
        from data_prep.generate_pesticide import generate_raw_dataset
        raw = generate_raw_dataset()
    else:
        raw = pd.read_csv(args.raw_csv)
    filtered = filter_by_qed(raw)
    result = clean_inventory_centroids(filtered, inventory)
    targets = result.summary.head(100)["centroid_smiles"].tolist()

    # Compute everything before writing the three related artifacts.
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.regenerate_raw:
        raw.to_csv(output_dir / "pesticides_raw.csv", index=False)
    result.assignments.to_csv(output_dir / "pesticides_update_filter_clusters.csv", index=False)
    result.summary.to_csv(output_dir / "pesticides_update_filter_cluster_summary.csv", index=False)
    # Preserve the original benchmark's exact SMILES strings and file format.
    (output_dir / "targets_no_bb_top100clust.smi").write_text("\n".join(targets), encoding="utf-8")
    print(f"Raw records: {len(raw)}; valid records with QED >= 0.5: {len(filtered)}")
    for iteration, removal in enumerate(result.removals, start=1):
        print(f"Round {iteration}: removed {removal.centroids_removed} inventory centroids "
              f"({removal.rows_removed} rows); {removal.rows_remaining} records remain")
    print(f"Final records: {len(result.assignments)}; clusters: {len(result.summary)}; targets: {len(targets)}")
    print(f"Outputs written to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
