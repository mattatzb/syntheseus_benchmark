#!/usr/bin/env python3
"""
Visualization of solved and unsolved molecules.

This module generates RDKit grid images of sampled molecules from benchmark results,
allowing visual inspection of which targets were successfully solved vs unsolved.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Iterator

from rdkit import Chem
from rdkit.Chem import Draw

from .utils import is_finite_number, iter_target_dirs
from config import Config


# Constants
EXCLUDED_METHOD_DIRS = {"analysis", "configs"}


def target_record_from_stats(stats: dict[str, object], target_dir: Path) -> dict[str, object]:
    """
    Create a record object from a target's stats.json file.

    Args:
        stats: Dictionary loaded from stats.json.
        target_dir: Path to the target directory.

    Returns:
        Dictionary containing target metadata including index, SMILES, solving status, and timing.
    """
    wallclock = stats.get("soln_time_wallclock")
    solved = is_finite_number(wallclock)
    return {
        "index": stats.get("index", target_dir.name),
        "smiles": stats.get("smiles", ""),
        "target_dir": str(target_dir),
        "solved": solved,
        "soln_time_wallclock": wallclock,
        "rxn_model_calls_used": stats.get("rxn_model_calls_used"),
    }


def load_method_records(method_dir: Path) -> tuple[Path | None, list[dict[str, object]]]:
    """
    Load all target records from a method's latest run directory.

    Reads stats.json files from each target directory and filters out inventory-owned molecules.

    Args:
        method_dir: Path to the method directory.

    Returns:
        Tuple of (run_dir, records) where:
            - run_dir: Path to the latest run directory (or None if not found)
            - records: List of target record dictionaries
    """
    if method_dir is None:
        return None, []

    records: list[dict[str, object]] = []
    for target_dir in iter_target_dirs(method_dir):
        stats_path = target_dir / "stats.json"
        if not stats_path.exists():
            continue

        with stats_path.open("r", encoding="utf-8") as f:
            stats = json.load(f)

        if bool(stats.get("target_in_inventory", False)):
            continue

        records.append(target_record_from_stats(stats, target_dir))

    return method_dir, records


def sample_records(
    records: list[dict[str, object]],
    n_solved: int,
    n_unsolved: int,
    rng: random.Random,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """
    Randomly sample solved and unsolved molecules from records.

    Args:
        records: List of target records.
        n_solved: Number of solved molecules to sample.
        n_unsolved: Number of unsolved molecules to sample.
        rng: Random number generator instance for reproducibility.

    Returns:
        Tuple of (solved_sample, unsolved_sample) lists.
    """
    solved = [r for r in records if bool(r["solved"])]
    unsolved = [r for r in records if not bool(r["solved"])]
    return (
        rng.sample(solved, k=min(n_solved, len(solved))),
        rng.sample(unsolved, k=min(n_unsolved, len(unsolved))),
    )


def save_molecule_grid(
    records: list[dict[str, object]],
    out_path: Path,
    title_prefix: str,
    mols_per_row: int = 5,
    sub_img_size: tuple[int, int] = (300, 240),
) -> int:
    """
    Generate an RDKit grid image from a list of molecule records.

    Renders molecules to a PNG image using RDKit's MolsToGridImage function.
    Skips molecules with invalid SMILES.

    Args:
        records: List of target records containing SMILES strings.
        out_path: Path where PNG image will be saved.
        title_prefix: Prefix for molecule legends (e.g., "Solved" or "Unsolved").
        mols_per_row: Number of molecules per row in grid (default: 5).
        sub_img_size: Size of each molecule image as (width, height) tuple (default: (300, 240)).

    Returns:
        Count of records with invalid SMILES that were skipped.
    """
    mols = []
    legends = []
    invalid_smiles = 0
    for record in records:
        smiles = str(record.get("smiles", ""))
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_smiles += 1
            continue
        mols.append(mol)
        legends.append(f"{title_prefix} idx={record['index']}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not mols:
        out_path.write_text("No valid molecules to render.\n", encoding="utf-8")
        return invalid_smiles

    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=mols_per_row,
        subImgSize=sub_img_size,
        legends=legends,
        useSVG=False,
    )
    img.save(str(out_path))
    return invalid_smiles


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the visualization script.

    Returns:
        Parsed arguments for results root, sample sizes, random seed, and output directory.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Sample molecules from each method's latest run directory and display solved/unsolved sets."
        )
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Config.RESULTS_ROOT,
        help="Root directory containing method folders (default: results_6_runs).",
    )
    parser.add_argument(
        "--n-solved",
        type=int,
        default=10,
        help="Number of solved molecules to sample per method (default: 10).",
    )
    parser.add_argument(
        "--n-unsolved",
        type=int,
        default=15,
        help="Number of unsolved molecules to sample per method (default: 15).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for reproducible sampling (default: 0).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Config.DEFAULT_SAMPLE_DIR,
        help="Directory where RDKit image grids are saved.",
    )
    return parser.parse_args()


def main() -> None:
    """
    Main entry point for molecule visualization.

    For each method in results root, loads target records, samples solved/unsolved molecules,
    and generates RDKit grid images with visual and text output.
    """
    args = parse_args()
    rng = random.Random(args.seed)

    method_dirs = sorted(
        d
        for d in args.results_root.iterdir()
        if d.is_dir() and d.name not in EXCLUDED_METHOD_DIRS
    )

    for method_dir in method_dirs:
        run_dir, records = load_method_records(method_dir)
        print(f"\n=== {method_dir.name} ===")
        if run_dir is None:
            print("No run directory found.")
            continue

        solved_sample, unsolved_sample = sample_records(
            records=records,
            n_solved=args.n_solved,
            n_unsolved=args.n_unsolved,
            rng=rng,
        )

        total_solved = sum(1 for r in records if bool(r["solved"]))
        total_unsolved = len(records) - total_solved

        print(f"Run: {run_dir.name}")
        print(
            f"Available (excluding inventory): solved={total_solved}, unsolved={total_unsolved}, total={len(records)}"
        )
        print(
            f"Sampled: solved={len(solved_sample)} (requested {args.n_solved}), "
            f"unsolved={len(unsolved_sample)} (requested {args.n_unsolved})"
        )

        print("\nSolved molecules:")
        if not solved_sample:
            print("  (none)")
        for i, record in enumerate(solved_sample, start=1):
            print(
                f"  {i:>2}. idx={record['index']}, wall={record['soln_time_wallclock']}, "
                f"smiles={record['smiles']}"
            )

        print("\nUnsolved molecules:")
        if not unsolved_sample:
            print("  (none)")
        for i, record in enumerate(unsolved_sample, start=1):
            print(f"  {i:>2}. idx={record['index']}, smiles={record['smiles']}")

        solved_img_path = args.output_dir / method_dir.name / "solved.png"
        unsolved_img_path = args.output_dir / method_dir.name / "unsolved.png"
        solved_invalid = save_molecule_grid(
            solved_sample,
            solved_img_path,
            title_prefix="Solved",
        )
        unsolved_invalid = save_molecule_grid(
            unsolved_sample,
            unsolved_img_path,
            title_prefix="Unsolved",
        )
        print("\nRDKit visualization files:")
        print(f"  solved:   {solved_img_path}")
        print(f"  unsolved: {unsolved_img_path}")
        if solved_invalid or unsolved_invalid:
            print(
                f"  skipped invalid SMILES: solved={solved_invalid}, unsolved={unsolved_invalid}"
            )


if __name__ == "__main__":
    main()
