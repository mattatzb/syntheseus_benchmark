"""Exercise both pickle extractors with small, independently specified routes."""

from __future__ import annotations

import json
import pickle
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analysis.mcts_extractor import extract_mcts_data
from analysis.retro_star_extractor import extract_data
from analysis import route_processor


def molecule(smiles, purchasable=False):
    return SimpleNamespace(smiles=smiles, metadata={"is_purchasable": purchasable})


class PickleExtractorTests(unittest.TestCase):
    def test_branched_retro_star_route_counts_all_reactions(self):
        # Both precursors of the root reaction require another reaction:
        # three reactions in total, two along either longest branch, depth four.
        def mol_node(smiles, depth, purchasable=False):
            return SimpleNamespace(depth=depth, mol=molecule(smiles, purchasable))

        def reaction(smiles, depth):
            return SimpleNamespace(depth=depth, reaction=smiles)

        route = [
            mol_node("O", 4, True),
            reaction("CC.O>>CCO", 3),
            mol_node("CCOCO", 0),
            mol_node("CC", 4, True),
            mol_node("CO", 2),
            reaction("CCO.CO>>CCOCO", 1),
            mol_node("C", 4, True),
            mol_node("CCO", 2),
            reaction("C.O>>CO", 3),
            mol_node("O", 4, True),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "branched.pkl"
            path.write_bytes(pickle.dumps(route))
            nodes = list(json.loads(extract_data(str(path))).values())
            self.assertEqual([n["depth"] for n in nodes], [0, 1, 2, 2, 3, 3, 4, 4, 4, 4])
            self.assertEqual(sum(n["is_mol"] for n in nodes), 7)
            self.assertEqual(route_processor.route_length_from_pkl(path, extract_data), 3)
            self.assertEqual(route_processor.reactants_count_from_pkl(path, extract_data), [2, 2, 2])
            self.assertEqual(route_processor.building_blocks_count_from_pkl(path, extract_data), 3)

    def test_mcts_extracts_ordered_states_and_only_new_reactants(self):
        block = molecule("C", True)
        block.metadata.update({"rdkit_mol": object(), "recorded_at": datetime(2026, 1, 2)})
        root = SimpleNamespace(depth=0, mols=[molecule("CCCO")])
        intermediate = SimpleNamespace(depth=1, mols=[molecule("CCO"), block])
        terminal = SimpleNamespace(depth=2, mols=[molecule("O", True), block, molecule("CC", True)])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "states.pkl"
            path.write_bytes(pickle.dumps([terminal, root, intermediate]))
            nodes = list(json.loads(extract_mcts_data(str(path))).values())
            self.assertEqual([n["depth"] for n in nodes], [0, 1, 2])
            self.assertEqual([m["smiles"] for m in nodes[2]["mols"]], ["C", "CC", "O"])
            self.assertEqual([n["inferred_reactant_smiles"] for n in nodes],
                             [[], ["C", "CCO"], ["CC", "O"]])
            self.assertEqual([n["inferred_expanded_mol_smiles"] for n in nodes],
                             [[], ["CCCO"], ["CCO"]])
            metadata = nodes[2]["mols"][0]["metadata"]
            self.assertNotIn("rdkit_mol", metadata)
            self.assertEqual(metadata["recorded_at"], "2026-01-02T00:00:00")
            self.assertEqual(route_processor.route_length_from_pkl(path, extract_mcts_data), 2)
            self.assertEqual(route_processor.reactants_count_from_pkl(path, extract_mcts_data), [2, 2])
            self.assertEqual(route_processor.building_blocks_count_from_pkl(path, extract_mcts_data), 3)


if __name__ == "__main__":
    unittest.main()
