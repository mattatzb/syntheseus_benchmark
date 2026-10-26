"""Check that the benchmark inputs match all archived searches and configurations."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import yaml
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGNS = [("1_top_results", 1), ("50_top_results", 50)]


class BenchmarkInputTests(unittest.TestCase):
    def test_targets_are_unique_outside_inventory_and_match_every_campaign_in_order(self):
        targets = (ROOT / "data/targets_no_bb_top100clust.smi").read_text().splitlines()
        self.assertEqual(len(targets), 100)
        canonical = []
        for index, smiles in enumerate(targets):
            mol = Chem.MolFromSmiles(smiles)
            self.assertIsNotNone(mol, f"Invalid target {index}: {smiles}")
            canonical.append(Chem.MolToSmiles(mol))
        self.assertEqual(len(set(canonical)), 100)
        # Syntheseus configs use canonicalize_inventory=False: inventory lines
        # are compared with canonical target SMILES as supplied to the search.
        inventory = set((ROOT / "data/enamine_sep_lpdc_blocks.smi").read_text().splitlines())
        self.assertFalse(set(canonical) & inventory)

        for campaign, _ in CAMPAIGNS:
            for method in sorted((ROOT / "syntheseus_outputs" / campaign).iterdir()):
                if not method.is_dir():
                    continue
                with self.subTest(campaign=campaign, method=method.name):
                    indices = {p.name for p in method.iterdir() if p.is_dir() and p.name.isdigit()}
                    self.assertEqual(indices, {str(i) for i in range(100)})
                    for index, smiles in enumerate(targets):
                        stats = json.loads((method / str(index) / "stats.json").read_text())
                        self.assertEqual(stats["smiles"], smiles, f"{method.name}, target {index}")
                        self.assertIs(stats["target_in_inventory"], False)

    def test_search_configs_match_archived_campaigns_and_include_complete_pistachio_model(self):
        config_dir = ROOT / "syntheseus_outputs/configs"
        for campaign, top_n in CAMPAIGNS:
            for method in sorted((ROOT / "syntheseus_outputs" / campaign).iterdir()):
                if not method.is_dir():
                    continue
                with self.subTest(campaign=campaign, method=method.name):
                    suffix = "_top50" if top_n == 50 else ""
                    config = yaml.safe_load((config_dir / f"{method.name}{suffix}.yml").read_text())
                    algorithm = "mcts" if method.name.endswith("_mcts") else "retro_star"
                    model = method.name.removeprefix("Pistachio_").removesuffix(f"_{algorithm}")
                    self.assertEqual(config["model_class"], model)
                    self.assertEqual(config["search_algorithm"], algorithm)
                    self.assertEqual(config["num_top_results"], top_n)
                    self.assertEqual((ROOT / config["results_dir"]).resolve(), method.resolve())
                    self.assertEqual((ROOT / config["search_targets_file"]).resolve(),
                                     ROOT / "data/targets_no_bb_top100clust.smi")
                    self.assertEqual((ROOT / config["inventory_smiles_file"]).resolve(),
                                     ROOT / "data/enamine_sep_lpdc_blocks.smi")
                    if top_n == 50:
                        top_one = yaml.safe_load((config_dir / f"{method.name}.yml").read_text())
                        keys = {"num_top_results", "results_dir"}
                        self.assertEqual({k: v for k, v in config.items() if k not in keys},
                                         {k: v for k, v in top_one.items() if k not in keys})
                    if method.name.startswith("Pistachio_"):
                        model_dir = ROOT / config["model_dir"]
                        self.assertEqual(model_dir.resolve(), ROOT / "MEGAN_backward")
                        for name in ["config.gin", "model_best.pt",
                                     "megan_default_32_bfs_randat/all_actions.json",
                                     "megan_default_32_bfs_randat/prop2oh.json"]:
                            path = model_dir / name
                            self.assertTrue(path.is_file(), str(path))
                            self.assertGreater(path.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
