"""Protect target classification and reproducible molecule sampling."""

from __future__ import annotations

import json
import math
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analysis.visualizer import load_method_records, sample_records


class MoleculeSamplingTests(unittest.TestCase):
    def test_record_loading_excludes_inventory_and_classifies_nonfinite_times(self):
        with tempfile.TemporaryDirectory() as tmp:
            method = Path(tmp)
            cases = [("0", 0.0, False), ("1", math.inf, False),
                     ("2", math.nan, False), ("3", None, False),
                     ("4", 2.0, True), ("analysis", 1.0, False)]
            for name, wall, inventory in cases:
                target = method / name
                target.mkdir()
                (target / "stats.json").write_text(json.dumps({
                    "smiles": "CCO", "soln_time_wallclock": wall,
                    "target_in_inventory": inventory,
                }))
            (method / "5").mkdir()  # Incomplete target without stats.
            returned_path, records = load_method_records(method)
            self.assertEqual(returned_path, method)
            self.assertEqual([r["index"] for r in records], ["0", "1", "2", "3"])
            self.assertEqual([r["solved"] for r in records], [True, False, False, False])

    def test_seeded_sampling_is_repeatable_separate_and_capped_by_availability(self):
        records = [{"index": i, "solved": i < 6} for i in range(9)]
        original = [dict(r) for r in records]
        first = sample_records(records, 3, 99, random.Random(42))
        second = sample_records(records, 3, 99, random.Random(42))
        self.assertEqual(first, second)
        solved, unsolved = first
        self.assertEqual(len(solved), 3)
        self.assertEqual(len({r["index"] for r in solved}), 3)
        self.assertTrue(all(r["solved"] for r in solved))
        self.assertEqual({r["index"] for r in unsolved}, {6, 7, 8})
        self.assertEqual(records, original)


if __name__ == "__main__":
    unittest.main()
