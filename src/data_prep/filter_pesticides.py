"""The historical pesticide filter: valid SMILES and QED >= 0.5 only."""

from __future__ import annotations

import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """Return a molecule from SMILES, or None if invalid."""
    if not isinstance(smiles, str) or not smiles or smiles == "Failed":
        return None
    return Chem.MolFromSmiles(smiles)


def compute_qed(mol: Chem.Mol | None) -> float | None:
    """Return QED for a molecule, or None if invalid."""
    return None if mol is None else float(QED.qed(mol))


def filter_by_qed(raw: pd.DataFrame, threshold: float = 0.5) -> pd.DataFrame:
    """Keep qualifying rows in source order, without canonicalizing or deduplicating."""
    if "smiles" not in raw.columns:
        raise ValueError("Missing required 'smiles' column in the raw pesticide data")
    filtered = raw.copy()
    filtered["qed"] = filtered["smiles"].apply(mol_from_smiles).apply(compute_qed)
    return filtered.loc[filtered["qed"] >= threshold].reset_index(drop=True)
