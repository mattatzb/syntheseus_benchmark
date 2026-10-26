"""Regression tests captured before the failure-handling changes.

Run from the repository root with python -B -m unittest discover -s tests -v.
All synthetic input/output files live in temporary directories, never outputs/.
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import Mock, patch

from rdkit import rdBase

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analysis import route_processor as processor


def reaction_route():
    return {
        "root": {"is_mol": True, "mol_smiles": "CCCO", "depth": 0},
        "reaction_1": {"is_rxn": True, "rxn_smiles": "CCO.C>>CCCO", "depth": 1},
        "intermediate": {"is_mol": True, "mol_smiles": "CCO", "depth": 2},
        "block_1": {"is_mol": True, "mol_smiles": "C", "is_purchasable": True, "depth": 2},
        "reaction_2": {"is_rxn": True, "rxn_smiles": "CC.O>>CCO", "depth": 3},
        "block_2": {"is_mol": True, "mol_smiles": "CC", "is_purchasable": True, "depth": 4},
        "block_3": {"is_mol": True, "mol_smiles": "O", "is_purchasable": True, "depth": 4},
    }


def state_route():
    route = {}
    states = [("CCCO",), ("CCO", "C"), ("CC", "O", "C")]
    for depth, smiles in enumerate(states):
        previous = set(states[depth - 1]) if depth else set(smiles)
        route[str(depth)] = {
            "is_state": True,
            "is_mol": False,
            "is_rxn": False,
            "depth": depth,
            "mols": [{"smiles": s, "is_purchasable": s in {"C", "CC", "O"}} for s in smiles],
            "inferred_reactant_count": len(set(smiles) - previous),
            "inferred_expanded_mol_smiles": sorted(previous - set(smiles)),
        }
    return route


def make_target(method, index, wall=2.0, inventory=False, routes=1):
    target = method / str(index)
    target.mkdir(parents=True)
    (target / "stats.json").write_text(json.dumps({
        "soln_time_wallclock": wall,
        "target_in_inventory": inventory,
    }))
    for i in range(routes):
        (target / f"route_{i}.pkl").touch()
    return target


class RouteProcessorRegressionTests(unittest.TestCase):
    def test_reaction_route_metrics(self):
        extractor = Mock(return_value=json.dumps(reaction_route()))
        route = Path("route.pkl")
        self.assertEqual(processor.route_length_from_pkl(route, extractor), 2.0)
        self.assertEqual(processor.reactants_count_from_pkl(route, extractor), [2, 2])
        self.assertEqual(processor.building_blocks_count_from_pkl(route, extractor), 3)

    def test_state_route_metrics(self):
        extractor = Mock(return_value=json.dumps(state_route()))
        route = Path("route.pkl")
        self.assertEqual(processor.route_length_from_pkl(route, extractor), 2.0)
        self.assertEqual(processor.reactants_count_from_pkl(route, extractor), [2, 2])
        self.assertEqual(processor.building_blocks_count_from_pkl(route, extractor), 3)

    def test_method_statistics_keep_existing_weighting_and_inventory_exclusions(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            method = Path(tmp) / "experiment" / "Example_retro_star"
            make_target(method, 0, wall=2.0, routes=2)
            make_target(method, 1, wall=4.0)
            make_target(method, 2, wall=math.inf, routes=0)
            make_target(method, 3, wall=0.0, inventory=True, routes=0)
            row = processor.analyze_method(method, Mock(return_value=json.dumps(reaction_route())))
        expected = {
            "method": "Example_retro_star", "run_dir": "experiment",
            "n_total_targets": 4, "n_targets": 3, "n_inventory_excluded": 1,
            "inventory_true_pct": 25.0, "n_solved": 2,
            "solving_rate": 2 / 3, "mean_wall_time_s": 3.0,
            "avg_route_length": 2.0, "avg_num_routes": 1.5,
            "route_lengths": [2.0, 2.0, 2.0],
            "reaction_reactant_counts": [2] * 6,
            "route_building_block_counts": [3, 3, 3], "note": "",
        }
        self.assertEqual(row, expected)

    def test_empty_extracted_route_has_no_length(self):
        self.assertTrue(math.isnan(processor.route_length_from_pkl(Path("route.pkl"), lambda _: "{}")))

    def test_building_block_identity_uses_canonical_smiles_and_preserves_stereochemistry(self):
        molecules = [("OCC", True), ("CCO", True),
                     ("C[C@H](O)F", True), ("C[C@@H](O)F", True),
                     ("N", False), ("not-a-smiles", True), (None, True)]
        reaction_nodes = {str(i): {"is_mol": True, "mol_smiles": smiles, "is_purchasable": owned}
                          for i, (smiles, owned) in enumerate(molecules)}
        state_nodes = {"state": {"mols": [{"smiles": smiles, "is_purchasable": owned}
                                           for smiles, owned in molecules]}}
        for name, data in [("reaction nodes", reaction_nodes), ("state nodes", state_nodes)]:
            with self.subTest(format=name), rdBase.BlockLogs():
                self.assertEqual(processor.building_blocks_count_from_pkl(
                    Path("unused.pkl"), lambda _: json.dumps(data)), 3)

    def test_route_length_mean_weights_exported_routes_rather_than_targets(self):
        lengths = {("0", "route_0.pkl"): 1, ("0", "route_1.pkl"): 3,
                   ("1", "route_0.pkl"): 5}

        def extractor(path):
            path = Path(path)
            length = lengths[(path.parent.name, path.name)]
            return json.dumps({str(i): {"is_rxn": True, "rxn_smiles": "C.O>>CO"}
                               for i in range(length)})

        with tempfile.TemporaryDirectory() as tmp:
            method = Path(tmp) / "Example_retro_star"
            make_target(method, 0, routes=2)
            make_target(method, 1, routes=1)
            row = processor.analyze_method(method, extractor)
        self.assertEqual(row["route_lengths"], [1, 3, 5])
        self.assertEqual(row["avg_route_length"], 3)
        self.assertEqual(row["avg_num_routes"], 1.5)


class RouteProcessorSafetyTests(unittest.TestCase):
    def test_parse_errors_raise_with_the_route_path_and_original_cause(self):
        for function in ("route_length_from_pkl", "reactants_count_from_pkl", "building_blocks_count_from_pkl"):
            with self.subTest(function=function), tempfile.TemporaryDirectory(dir="/tmp") as tmp:
                method = Path(tmp) / "Example_retro_star"
                target = make_target(method, 0)
                error = ModuleNotFoundError("No module named 'syntheseus'")
                with patch.object(processor, function, side_effect=error):
                    with self.assertRaisesRegex(RuntimeError, "route_0.pkl") as caught:
                        processor.analyze_method(method, Mock(return_value=json.dumps(reaction_route())))
                self.assertIs(caught.exception.__cause__, error)
                self.assertIn(str(target), str(caught.exception))

    def test_invalid_route_length_cannot_silently_drop_a_route(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            method = Path(tmp) / "Example_retro_star"
            make_target(method, 0)
            with self.assertRaisesRegex(RuntimeError, "route_0.pkl"):
                processor.analyze_method(method, Mock(return_value="{}"))

    def test_unreliable_state_reactants_are_rejected_instead_of_returning_a_partial_list(self):
        for defect in ("zero_new_reactants", "two_expanded_molecules", "nonconsecutive_depth"):
            data = state_route()
            if defect == "zero_new_reactants":
                data["2"]["inferred_reactant_count"] = 0
            elif defect == "two_expanded_molecules":
                data["2"]["inferred_expanded_mol_smiles"] = ["CCO", "C"]
            else:
                data["2"]["depth"] = 4
            with self.subTest(defect=defect), self.assertRaisesRegex(ValueError, "reactant"):
                processor.reactants_count_from_pkl(Path("route.pkl"), lambda _: json.dumps(data))

    def test_unreliable_reactants_omit_the_whole_method_distribution_but_preserve_other_metrics(self):
        good = state_route()
        bad = state_route()
        bad["2"]["inferred_reactant_count"] = 0
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            method = Path(tmp) / "Example_mcts"
            make_target(method, 0, routes=2)
            extractor = lambda path: json.dumps(bad if path.endswith("route_1.pkl") else good)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                row = processor.analyze_method(method, extractor)
        self.assertEqual(row["reaction_reactant_counts"], [])
        self.assertEqual(row["route_lengths"], [2.0, 2.0])
        self.assertEqual(row["route_building_block_counts"], [3, 3])
        self.assertEqual(row["solving_rate"], 1.0)
        self.assertEqual(row["mean_wall_time_s"], 2.0)
        self.assertEqual(row["avg_num_routes"], 2.0)
        self.assertIn("omitted", row["note"])
        self.assertIn("Example_mcts", row["note"])
        self.assertIn("route_1.pkl", row["note"])
        self.assertTrue(any("reactant" in str(w.message) and "omitted" in str(w.message) for w in caught))


if __name__ == "__main__":
    unittest.main()
