"""Protect the recovered preparation procedure and the published target order."""

from __future__ import annotations

import contextlib
import hashlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cli.prepare_data import main
from data_prep.filter_pesticides import filter_by_qed
from data_prep.filter_pesticides_by_inventory import clean_inventory_centroids


class DataPreparationTests(unittest.TestCase):
    def test_filter_uses_only_validity_and_inclusive_qed_preserving_order_and_duplicates(self):
        raw = pd.DataFrame({"smiles": ["CC", "[Zn]", "c1ccccc1", "CC", "CO", "Failed", None]})
        # The first four records pass the boundary, including the structures
        # that the unused metal/short-chain/simple-ring filters excluded.
        def score(mol):
            if mol is None:
                return None
            return 0.49 if Chem.MolToSmiles(mol) == "CO" else 0.5

        with patch("data_prep.filter_pesticides.compute_qed", new=score):
            filtered = filter_by_qed(raw)
        self.assertEqual(filtered.smiles.tolist(), raw.smiles.tolist()[:4])
        self.assertEqual(filtered.qed.tolist(), [0.5] * 4)
        self.assertEqual(raw.columns.tolist(), ["smiles"])

    def test_inventory_cleaning_reclusters_and_removes_only_exact_centroid_rows(self):
        raw = pd.DataFrame({"smiles": ["CCO", "CCO", "OCC", "CCN", "CCCN", "CCCCN", "c1ccccc1"]})
        centroids = iter(["CCO", "CCN", "c1ccccc1"])

        def cluster(frame):
            centroid = next(centroids)
            self.assertIn(centroid, frame.smiles.tolist())
            summary = pd.DataFrame({"centroid_smiles": [centroid]})
            return frame.copy(), summary

        with patch("data_prep.filter_pesticides_by_inventory.cluster_molecules",
                   side_effect=cluster) as clustering:
            result = clean_inventory_centroids(raw, {"CCO", "CCN", "CCCCN"}, num_targets=1)
        self.assertEqual(clustering.call_count, 3)
        # OCC is canonical CCO, but was not the raw centroid string; it stays.
        # CCCCN is in inventory but is never a selected centroid; it stays too.
        self.assertEqual(result.assignments.smiles.tolist(), ["OCC", "CCCN", "CCCCN", "c1ccccc1"])
        self.assertEqual([(r.centroids_removed, r.rows_removed, r.rows_remaining)
                          for r in result.removals], [(1, 2, 5), (1, 1, 4)])

    def test_archived_raw_reproduces_all_targets_and_consistent_csvs_offline(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            log = io.StringIO()
            with patch("requests.sessions.Session.request", side_effect=AssertionError("Network forbidden")), \
                    contextlib.redirect_stdout(log):
                self.assertEqual(main(["--output-dir", directory]), 0)
            for message in ["valid records with QED >= 0.5: 1590",
                            "removed 19 inventory centroids (47 rows); 1543 records remain",
                            "removed 4 inventory centroids (8 rows); 1535 records remain",
                            "removed 2 inventory centroids (5 rows); 1530 records remain",
                            "removed 1 inventory centroids (3 rows); 1527 records remain"]:
                self.assertIn(message, log.getvalue())
            targets = (output / "targets_no_bb_top100clust.smi").read_bytes()
            self.assertEqual(hashlib.sha256(targets).hexdigest(),
                             "f9dcbfd502cc8dacd531bf7906f22543baba9ffc043fb4140ecf41ef644be622")
            self.assertEqual(targets, (ROOT / "data/targets_no_bb_top100clust.smi").read_bytes())
            assignments = pd.read_csv(output / "pesticides_update_filter_clusters.csv")
            summary = pd.read_csv(output / "pesticides_update_filter_cluster_summary.csv")
            self.assertEqual((len(assignments), len(summary)), (1527, 918))
            self.assertEqual(summary.centroid_smiles.iloc[:100].tolist(), targets.decode().splitlines())
            for row in summary.itertuples(index=False):
                centroid = assignments.iloc[row.centroid_row_index]
                self.assertEqual(centroid.smiles, row.centroid_smiles)
                self.assertEqual(centroid.cluster_id, row.cluster_id)
                members = assignments[assignments.cluster_id == row.cluster_id]
                self.assertEqual(len(members), row.cluster_size)
                self.assertTrue((members.cluster_size == row.cluster_size).all())
            self.assertFalse({"has_heavy_metal", "is_short_aliphatic",
                              "is_single_simple_aromatic_cycle"} & set(assignments.columns))
            for name in ["pesticides_update_filter_clusters.csv",
                         "pesticides_update_filter_cluster_summary.csv"]:
                self.assertEqual((output / name).read_bytes(), (ROOT / "data" / name).read_bytes(), name)
            self.assertEqual(set(p.name for p in output.iterdir()), {
                "pesticides_update_filter_clusters.csv",
                "pesticides_update_filter_cluster_summary.csv", "targets_no_bb_top100clust.smi"})

    def test_wrong_numpy_version_fails_before_writing_outputs(self):
        with tempfile.TemporaryDirectory() as directory, \
                patch("cli.prepare_data.np.__version__", "1.24.3"):
            with self.assertRaisesRegex(ValueError, "NumPy 1.26.4"):
                main(["--output-dir", directory])
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_invalid_inputs_fail_before_writing_outputs(self):
        with self.assertRaisesRegex(ValueError, "smiles"):
            filter_by_qed(pd.DataFrame({"name": ["missing structure"]}))
        with self.assertRaisesRegex(ValueError, "100 clusters"):
            clean_inventory_centroids(pd.DataFrame({"smiles": ["CCO"]}), set())
        archived_raw = ROOT / "data/pesticides_raw.csv"
        before = archived_raw.read_bytes()
        with self.assertRaisesRegex(ValueError, "preserve the archived raw CSV"):
            main(["--regenerate-raw", "--output-dir", str(ROOT / "data")])
        self.assertEqual(archived_raw.read_bytes(), before)

    def test_live_refresh_is_explicit_and_keeps_the_archived_raw_unchanged(self):
        raw_path = ROOT / "data/pesticides_raw.csv"
        before = raw_path.read_bytes()
        with tempfile.TemporaryDirectory() as directory, \
                patch("data_prep.generate_pesticide.generate_raw_dataset",
                      return_value=pd.read_csv(raw_path)) as generate, \
                contextlib.redirect_stdout(io.StringIO()):
            with self.assertWarnsRegex(UserWarning, "external databases"):
                self.assertEqual(main(["--regenerate-raw", "--output-dir", directory]), 0)
            generate.assert_called_once_with()
            self.assertEqual(len(pd.read_csv(Path(directory) / "pesticides_raw.csv")), 2564)
        self.assertEqual(raw_path.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
