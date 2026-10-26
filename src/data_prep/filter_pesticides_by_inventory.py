"""Reconstruct the historical recursive removal of purchasable centroids."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from rdkit import Chem

from .cluster_pesticides_butina import cluster_molecules


@dataclass(frozen=True)
class RemovalRound:
    centroids_removed: int
    rows_removed: int
    rows_remaining: int


@dataclass
class PreparationResult:
    assignments: pd.DataFrame
    summary: pd.DataFrame
    removals: list[RemovalRound]


def clean_inventory_centroids(
    molecules: pd.DataFrame, inventory: set[str], num_targets: int = 100,
) -> PreparationResult:
    """Recluster until the leading centroids are all outside the inventory.

    Compare canonical centroid SMILES against inventory lines as supplied, matching
    Syntheseus's canonicalize_inventory=False. Remove all rows whose raw SMILES
    equal an offending centroid's raw SMILES. Other inventory records remain in
    the clustering pool, including non-centroids and alternate SMILES spellings.
    """
    if num_targets < 1:
        raise ValueError("num_targets must be positive")
    remaining = molecules.copy()
    removals = []
    while True:
        assignments, summary = cluster_molecules(remaining)
        if len(summary) < num_targets:
            raise ValueError(f"Need {num_targets} clusters, found only {len(summary)}")
        selected = summary.head(num_targets)["centroid_smiles"]
        purchased = {
            smiles for smiles in selected
            if Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) in inventory
        }
        if not purchased:
            return PreparationResult(assignments, summary, removals)
        kept = remaining.loc[~remaining["smiles"].isin(purchased)].reset_index(drop=True)
        removed_rows = len(remaining) - len(kept)
        if removed_rows == 0:
            raise RuntimeError("Inventory cleaning made no progress")
        removals.append(RemovalRound(len(purchased), removed_rows, len(kept)))
        remaining = kept
