#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Iterator

from rdkit import Chem
from rdkit.Chem import Draw

EXCLUDED_METHOD_DIRS = {"analysis", "configs"}


def is_finite_number(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value)


def find_latest_run_dir(method_dir: Path) -> Path | None:
    run_dirs = sorted(d for d in method_dir.iterdir() if d.is_dir())
    return run_dirs[-1] if run_dirs else None


def iter_target_dirs(run_dir: Path) -> Iterator[Path]:
    for target_dir in sorted(run_dir.iterdir(), key=lambda p: p.name):
        if target_dir.is_dir() and target_dir.name.isdigit():
            yield target_dir


def target_record_from_stats(stats: dict[str, object], target_dir: Path) -> dict[str, object]:
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
    run_dir = find_latest_run_dir(method_dir)
    if run_dir is None:
        return None, []

    records: list[dict[str, object]] = []
    for target_dir in iter_target_dirs(run_dir):
        stats_path = target_dir / "stats.json"
        if not stats_path.exists():
            continue

        with stats_path.open("r", encoding="utf-8") as f:
            stats = json.load(f)

        if bool(stats.get("target_in_inventory", False)):
            continue

        records.append(target_record_from_stats(stats, target_dir))

    return run_dir, records


def sample_records(
    records: list[dict[str, object]],
    n_solved: int,
    n_unsolved: int,
    rng: random.Random,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
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
    parser = argparse.ArgumentParser(
        description=(
            "Sample molecules from each method's latest run directory and display solved/unsolved sets."
        )
    )
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
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
        default=Path(__file__).resolve().parent / "molecule_samples",
        help="Directory where RDKit image grids are saved.",
    )
    return parser.parse_args()


def main() -> None:
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
