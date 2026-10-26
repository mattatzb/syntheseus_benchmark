"""Cluster pesticide records with Morgan fingerprints and Butina clustering."""

from __future__ import annotations

import pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina


def morgan_fp(smiles: str, radius: int, n_bits: int):
    """Build the historical Morgan bit fingerprint, or return None if invalid."""
    if not isinstance(smiles, str) or not smiles or smiles == "Failed":
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # Preserve the historical API without its per-record deprecation warning.
    with rdBase.BlockLogs():
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def butina_cluster(fps, distance_cutoff: float):
    """Return Butina clusters, preserving input order and historical defaults."""
    if not fps:
        return []
    if len(fps) == 1:
        return [(0,)]
    distances = []
    for i in range(1, len(fps)):
        similarities = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        distances.extend(1.0 - similarity for similarity in similarities)
    return Butina.ClusterData(
        distances, len(fps), distance_cutoff, isDistData=True, reordering=False
    )


def cluster_molecules(
    molecules: pd.DataFrame, radius: int = 2, n_bits: int = 2048,
    distance_cutoff: float = 0.4,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return assignments and a summary ranked once by descending population.

    Centroid row indices refer to zero-based positions in the assignments table.
    Duplicate source records contribute separately to cluster populations.
    Quicksort tie order is historical and requires the pinned NumPy version.
    """
    if "smiles" not in molecules.columns or molecules.empty:
        raise ValueError("Clustering requires a nonempty table with a 'smiles' column")
    assignments = molecules.reset_index(drop=True).copy()
    fps = [morgan_fp(smiles, radius, n_bits) for smiles in assignments["smiles"]]
    if any(fp is None for fp in fps):
        raise ValueError("Clustering requires valid SMILES; apply the QED filter first")
    clusters = butina_cluster(fps, distance_cutoff)
    cluster_ids = [-1] * len(assignments)
    cluster_sizes = [0] * len(assignments)
    rows = []
    for cluster_id, cluster in enumerate(clusters):
        for row_index in cluster:
            cluster_ids[row_index] = cluster_id
            cluster_sizes[row_index] = len(cluster)
        centroid_index = cluster[0]
        rows.append({
            "cluster_id": cluster_id,
            "cluster_size": len(cluster),
            "centroid_row_index": centroid_index,
            "centroid_smiles": assignments.at[centroid_index, "smiles"],
        })
    assignments["cluster_id"] = cluster_ids
    assignments["cluster_size"] = cluster_sizes
    assignments["is_valid_smiles"] = True
    # Do not sort again: unstable ties would change the target order.
    summary = pd.DataFrame(rows).sort_values(
        "cluster_size", ascending=False, kind="quicksort", ignore_index=True
    )
    return assignments, summary
