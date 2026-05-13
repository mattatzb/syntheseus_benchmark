#!/usr/bin/env python3
"""Cluster molecules with Butina using Morgan fingerprints (RDKit)."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina


def morgan_fp(smiles: str, radius: int, n_bits: int):
    """Build a Morgan fingerprint from SMILES, or return None if invalid."""
    if not isinstance(smiles, str) or not smiles or smiles == "Failed":
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)


def butina_cluster(fps, distance_cutoff: float):
    """Return Butina clusters from a list of fingerprints."""
    n_fps = len(fps)
    if n_fps == 0:
        return []
    if n_fps == 1:
        return [(0,)]

    dists = []
    for i in range(1, n_fps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1.0 - x for x in sims)

    return Butina.ClusterData(dists, n_fps, distance_cutoff, isDistData=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cluster molecules from CSV with RDKit Morgan fingerprints + Butina."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("/home/tatzber/syntheseus_benchmark/pesticides_smiles/pesticides_update_filtered.csv"),
        help="Input CSV file.",
    )
    parser.add_argument(
        "--smiles-col",
        type=str,
        default="smiles",
        help="Name of SMILES column in input CSV.",
    )
    parser.add_argument(
        "--distance-cutoff",
        type=float,
        default=0.4,
        help="Butina distance cutoff (distance = 1 - Tanimoto).",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=2,
        help="Morgan fingerprint radius.",
    )
    parser.add_argument(
        "--n-bits",
        type=int,
        default=2048,
        help="Morgan fingerprint bit size.",
    )
    parser.add_argument(
        "--out-assignments",
        type=Path,
        default=Path("/home/tatzber/syntheseus_benchmark/pesticides_qed_ge_0_5_clusters.csv"),
        help="Output CSV with cluster assignments.",
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=Path("/home/tatzber/syntheseus_benchmark/pesticides_qed_ge_0_5_cluster_summary.csv"),
        help="Output CSV with per-cluster summary.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    if args.smiles_col not in df.columns:
        raise ValueError(f"Column '{args.smiles_col}' not found in {args.input}")

    fps = []
    valid_row_indices = []
    for row_idx, smiles in enumerate(df[args.smiles_col]):
        fp = morgan_fp(smiles, radius=args.radius, n_bits=args.n_bits)
        if fp is not None:
            valid_row_indices.append(row_idx)
            fps.append(fp)

    clusters = butina_cluster(fps, distance_cutoff=args.distance_cutoff)

    cluster_id_by_valid_fp_idx = {}
    for cluster_id, cluster in enumerate(clusters):
        for fp_idx in cluster:
            cluster_id_by_valid_fp_idx[fp_idx] = cluster_id

    df["cluster_id"] = -1
    for fp_idx, row_idx in enumerate(valid_row_indices):
        df.at[row_idx, "cluster_id"] = cluster_id_by_valid_fp_idx[fp_idx]

    cluster_sizes = (
        df[df["cluster_id"] >= 0].groupby("cluster_id", as_index=True).size().to_dict()
    )
    df["cluster_size"] = df["cluster_id"].map(cluster_sizes).fillna(0).astype(int)
    df["is_valid_smiles"] = df["cluster_id"] >= 0

    args.out_assignments.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_assignments, index=False)

    summary_rows = []
    for cluster_id, cluster in enumerate(clusters):
        centroid_fp_idx = cluster[0]
        centroid_row_idx = valid_row_indices[centroid_fp_idx]
        summary_rows.append(
            {
                "cluster_id": cluster_id,
                "cluster_size": len(cluster),
                "centroid_row_index": centroid_row_idx,
                "centroid_smiles": df.at[centroid_row_idx, args.smiles_col],
            }
        )

    summary_df = pd.DataFrame(summary_rows).sort_values(
        by="cluster_size", ascending=False, ignore_index=True
    )
    summary_df.to_csv(args.out_summary, index=False)

    print(f"Input rows: {len(df)}")
    print(f"Valid SMILES clustered: {len(valid_row_indices)}")
    print(f"Invalid/failed SMILES: {len(df) - len(valid_row_indices)}")
    print(f"Number of clusters: {len(clusters)}")
    print(f"Assignments written to: {args.out_assignments}")
    print(f"Summary written to: {args.out_summary}")


if __name__ == "__main__":
    main()
