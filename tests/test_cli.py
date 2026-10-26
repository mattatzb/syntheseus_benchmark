"""Campaign selection must keep top-one and top-fifty outputs separate."""

from __future__ import annotations

import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from analysis import aggregator
from config import Config


class CampaignOutputTests(unittest.TestCase):
    def parse(self, *arguments):
        with patch.object(sys, "argv", ["aggregate_results.py", *arguments]):
            return aggregator.parse_args()

    def test_default_campaign_uses_renamed_top_one_summary(self):
        with patch.object(Config, "RESULTS_ROOT", ROOT / "syntheseus_outputs/1_top_results"):
            args = self.parse()
        self.assertEqual(args.output_csv, Config.OUTPUT_ROOT / "1_top_results_summary.csv")
        self.assertEqual(args.plot_data_dir, Config.OUTPUT_ROOT)

    def test_top_fifty_uses_separate_summary_and_plot_directory(self):
        args = self.parse("--results-root", str(ROOT / "syntheseus_outputs/50_top_results"))
        self.assertEqual(args.output_csv, Config.OUTPUT_ROOT / "50_top_results_summary.csv")
        self.assertEqual(args.plot_data_dir, Config.OUTPUT_ROOT / "50_top_results")

    def test_configured_top_fifty_root_selects_top_fifty_outputs(self):
        with patch.object(Config, "RESULTS_ROOT", ROOT / "syntheseus_outputs/50_top_results"):
            args = self.parse()
        self.assertEqual(args.output_csv.name, "50_top_results_summary.csv")
        self.assertEqual(args.plot_data_dir, Config.OUTPUT_ROOT / "50_top_results")

    def test_method_selection_infers_campaign(self):
        args = self.parse("--method-dir", str(ROOT / "syntheseus_outputs/50_top_results/MEGAN_mcts"))
        self.assertEqual(args.output_csv.name, "50_top_results_summary.csv")
        self.assertEqual(args.plot_data_dir, Config.OUTPUT_ROOT / "50_top_results")

    def test_explicit_output_paths_are_preserved(self):
        args = self.parse("--results-root", "syntheseus_outputs/50_top_results",
                          "--output-csv", "/tmp/custom.csv", "--plot-data-dir", "/tmp/custom-plots")
        self.assertEqual(args.output_csv, Path("/tmp/custom.csv"))
        self.assertEqual(args.plot_data_dir, Path("/tmp/custom-plots"))

    def test_mixed_campaign_method_selection_is_rejected(self):
        with patch("sys.stderr", new=io.StringIO()), self.assertRaises(SystemExit) as caught:
            self.parse("--method-dir", "syntheseus_outputs/1_top_results/MEGAN_mcts",
                       "--method-dir", "syntheseus_outputs/50_top_results/MEGAN_mcts")
        self.assertEqual(caught.exception.code, 2)

    def test_removed_cli_options_are_rejected(self):
        for arguments in [("--run-dir", "/tmp/unused"), ("--merge-with-existing",)]:
            with self.subTest(arguments=arguments):
                with patch("sys.stderr", new=io.StringIO()), self.assertRaises(SystemExit) as caught:
                    self.parse(*arguments)
                self.assertEqual(caught.exception.code, 2)

    def test_custom_results_root_uses_neutral_summary_name(self):
        self.assertEqual(self.parse("--results-root", "/tmp/custom-results").output_csv.name,
                         "benchmark_summary.csv")

    def test_pdf_is_written_beside_selected_distribution_files(self):
        args = self.parse("--results-root", "syntheseus_outputs/50_top_results", "--plot-distributions")
        with patch.object(aggregator, "parse_args", return_value=args), \
             patch.object(Path, "iterdir", return_value=iter([])), \
             patch.object(aggregator, "write_csv"), \
             patch.object(aggregator, "write_plot_data_csvs", return_value={}), \
             patch.object(aggregator, "plot_distribution_csvs") as plotter, \
             patch("sys.stdout", new=io.StringIO()):
            aggregator.main()
        plotter.assert_called_once_with(args.plot_data_dir, args.plot_data_dir / "results_plot.pdf")


if __name__ == "__main__":
    unittest.main()
