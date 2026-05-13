#!/usr/bin/env python3
"""Filter pesticide compounds by QED and structural exclusion rules."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED


HEAVY_METAL_ATOMIC_NUMBERS = {
    13,  # Al
    22, 23, 24, 25, 26, 27, 28, 29, 30,  # Transition metals period 4
    31, 32, 33, 34,  # Heavier main-group atoms commonly excluded
    40, 42, 44, 45, 46, 47, 48, 50, 51, 52, 53,
    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
}


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """Return a molecule from SMILES, or None if invalid."""
    if not isinstance(smiles, str) or not smiles or smiles == "Failed":
        return None
    return Chem.MolFromSmiles(smiles)


def compute_qed(mol: Chem.Mol | None) -> float | None:
    """Return QED for a molecule, or None if invalid."""
    if mol is None:
        return None
    return float(QED.qed(mol))


def has_heavy_metal(mol: Chem.Mol | None) -> bool:
    """True when a molecule contains at least one heavy metal atom."""
    if mol is None:
        return False
    return any(atom.GetAtomicNum() in HEAVY_METAL_ATOMIC_NUMBERS for atom in mol.GetAtoms())


def is_short_aliphatic(mol: Chem.Mol | None, min_carbons: int = 5) -> bool:
    """True for non-aromatic molecules with fewer than `min_carbons` carbons."""
    if mol is None:
        return False

    has_aromatic_atom = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    return (not has_aromatic_atom) and carbon_count < min_carbons


def is_single_simple_aromatic_cycle(mol: Chem.Mol | None) -> bool:
    """True when molecule is exactly one aromatic ring with no extra heavy atoms."""
    if mol is None:
        return False

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) != 1:
        return False

    ring_atoms = set(atom_rings[0])
    if not ring_atoms:
        return False

    if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
        return False

    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    return heavy_atom_count == len(ring_atoms)


def main() -> None:
    parser = argparse.ArgumentParser(description="Filter compounds with QED >= threshold.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/home/tatzber/syntheseus_benchmark/pesticides_smiles/pesticides.csv"),
        help="Input CSV path (must contain a 'smiles' column).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/home/tatzber/syntheseus_benchmark/pesticides_qed_ge_0_5.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Minimum QED value to keep (default: 0.5).",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if "smiles" not in df.columns:
        raise ValueError(f"Missing required 'smiles' column in {args.input}")

    df["mol"] = df["smiles"].apply(mol_from_smiles)
    df["qed"] = df["mol"].apply(compute_qed)
    df["has_heavy_metal"] = df["mol"].apply(has_heavy_metal)
    df["is_short_aliphatic"] = df["mol"].apply(is_short_aliphatic)
    df["is_single_simple_aromatic_cycle"] = df["mol"].apply(is_single_simple_aromatic_cycle)

    keep_mask = (
        df["qed"].notna()
        & (df["qed"] >= args.threshold)
        & (~df["has_heavy_metal"])
        & (~df["is_short_aliphatic"])
        & (~df["is_single_simple_aromatic_cycle"])
    )
    filtered = df[keep_mask].drop(columns=["mol"]).copy()
    filtered.to_csv(args.output, index=False)

    print(
        f"Wrote {len(filtered)} / {len(df)} compounds to {args.output} "
        f"(QED >= {args.threshold}, no heavy metals, not short aliphatic, not single aromatic cycle)."
    )


if __name__ == "__main__":
    main()
