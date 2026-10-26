"""Protect the bundled benchmark outputs before changing error handling."""

from __future__ import annotations

import csv
import importlib.util
import io
import shutil
import subprocess
import sys
import tempfile
import unittest
import warnings
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from analysis import aggregator


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


@unittest.skipUnless(importlib.util.find_spec("syntheseus"), "Bundled route pickles require Syntheseus")
class BundledOutputRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = {}
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            for experiment in ("1_top_results", "50_top_results"):
                cls.rows[experiment] = [
                    aggregator.analyze_method(method, aggregator.extractor_for_method(method))
                    for method in sorted((ROOT / "syntheseus_outputs" / experiment).iterdir())
                    if method.is_dir()
                ]
        cls.analysis_warnings = [str(w.message) for w in caught if w.category is RuntimeWarning]

    def test_top_one_summary_matches_saved_csv_byte_for_byte(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            output = Path(tmp) / "summary.csv"
            aggregator.write_csv(self.rows["1_top_results"], output, False, False)
            self.assertEqual(output.read_bytes(), (ROOT / "outputs/1_top_results_summary.csv").read_bytes())

    def test_top_fifty_summary_matches_regenerated_csv(self):
        expected = read_csv(ROOT / "outputs/50_top_results_summary.csv")
        actual = [
            {k: str(v) for k, v in aggregator.format_row(row, False, False).items()}
            for row in self.rows["50_top_results"]
        ]
        self.assertEqual(actual, expected)

    def test_only_affected_top_fifty_methods_omit_reactants_and_warn(self):
        affected = {"MEGAN_mcts": 40, "RetroKNN_mcts": 3}
        for experiment, rows in self.rows.items():
            for row in rows:
                with self.subTest(experiment=experiment, method=row["method"]):
                    if experiment == "50_top_results" and row["method"] in affected:
                        self.assertEqual(row["reaction_reactant_counts"], [])
                        self.assertIn(f"{affected[row['method']]} route(s)", row["note"])
                        self.assertTrue(any(row["note"] == warning for warning in self.analysis_warnings))
                    else:
                        self.assertTrue(row["reaction_reactant_counts"])
                        self.assertEqual(row["note"], "")
        self.assertEqual(len(self.analysis_warnings), 2)

    def test_top_one_distribution_csvs_match_saved_files_byte_for_byte(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            outputs = aggregator.write_plot_data_csvs(self.rows["1_top_results"], Path(tmp))
            for output in outputs.values():
                with self.subTest(filename=output.name):
                    self.assertEqual(output.read_bytes(), (ROOT / "outputs" / output.name).read_bytes())

    @unittest.skipUnless(importlib.util.find_spec("matplotlib") and shutil.which("pdftoppm"),
                         "PDF comparison requires Matplotlib and pdftoppm")
    def test_top_one_pdf_render_matches_saved_figure(self):
        from analysis import plots
        buffer = io.BytesIO()
        pdf_pages = plots.PdfPages
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            aggregator.write_plot_data_csvs(self.rows["1_top_results"], Path(tmp))
            with patch.object(plots, "PdfPages", side_effect=lambda _: pdf_pages(buffer)):
                plots.plot_distribution_csvs(Path(tmp), Path(tmp) / "unused.pdf")
        command = ["pdftoppm", "-r", "72", "-singlefile", "-png", "-"]
        expected = subprocess.run(command, input=(ROOT / "outputs/results_plot.pdf").read_bytes(),
                                  capture_output=True, check=True).stdout
        actual = subprocess.run(command, input=buffer.getvalue(), capture_output=True, check=True).stdout
        self.assertEqual(actual, expected)


class AggregatorSafetyTests(unittest.TestCase):
    def test_distribution_bins_and_proportions_are_computed_separately_per_method(self):
        rows = [
            {"method": "A_mcts", "run_dir": "experiment", "route_lengths": [4, 1, 4]},
            {"method": "B_mcts", "run_dir": "experiment", "route_lengths": [1, 1]},
            {"method": "Empty_mcts", "run_dir": "experiment", "route_lengths": []},
        ]
        distribution = aggregator.build_distribution_rows(rows, "route_lengths", "route_length")
        self.assertEqual([(r["method"], r["value"], r["count"], r["proportion"])
                          for r in distribution],
                         [("A_mcts", 1, 1, "0.333333"), ("A_mcts", 4, 2, "0.666667"),
                          ("B_mcts", 1, 2, "1.000000")])

    def test_analysis_failure_leaves_existing_output_files_untouched(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            root = Path(tmp)
            method = root / "methods" / "Example_mcts"
            method.mkdir(parents=True)
            output = root / "summary.csv"
            output.write_text("existing summary\n")
            args = Namespace(method_dirs=[method],
                             output_csv=output, show_inventory_columns=False, show_debug_columns=False,
                             write_plot_data=True, plot_distributions=True, plot_data_dir=root)
            with patch.object(aggregator, "parse_args", return_value=args), \
                 patch.object(aggregator, "analyze_method", side_effect=RuntimeError("bad route")), \
                 patch.object(aggregator, "write_csv") as summary_writer, \
                 patch.object(aggregator, "write_plot_data_csvs") as data_writer, \
                 patch.object(aggregator, "plot_distribution_csvs") as plotter:
                with self.assertRaisesRegex(RuntimeError, "bad route"):
                    aggregator.main()
            self.assertEqual(output.read_text(), "existing summary\n")
            summary_writer.assert_not_called()
            data_writer.assert_not_called()
            plotter.assert_not_called()

    def test_notes_are_visible_on_stderr_without_debug_columns(self):
        row = {
            "method": "Example_mcts", "run_dir": "experiment", "n_total_targets": 1,
            "n_targets": 1, "n_solved": 1, "solving_rate": 1.0,
            "mean_wall_time_s": 2.0, "avg_route_length": 2.0, "avg_num_routes": 1.0,
            "note": "reactant counts omitted for Example_mcts: unreliable state transition",
        }
        args = Namespace(method_dirs=[Path("Example_mcts")],
                         output_csv=Path("unused.csv"), show_inventory_columns=False, show_debug_columns=False,
                         write_plot_data=False, plot_distributions=False)
        stderr = io.StringIO()
        with patch.object(aggregator, "parse_args", return_value=args), \
             patch.object(aggregator, "analyze_method", return_value=row), \
             patch.object(aggregator, "write_csv"), \
             patch("sys.stdout", new=io.StringIO()), patch("sys.stderr", new=stderr):
            aggregator.main()
        self.assertIn(row["note"], stderr.getvalue())

    def test_omitted_reactants_do_not_create_a_partial_distribution(self):
        rows = [{"method": "Example_mcts", "run_dir": "experiment",
                 "route_lengths": [2.0], "route_building_block_counts": [3],
                 "reaction_reactant_counts": [], "note": "reactant counts omitted"}]
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            paths = aggregator.write_plot_data_csvs(rows, Path(tmp))
            self.assertEqual(read_csv(paths["reaction_reactant_counts"]), [])
            self.assertEqual(len(read_csv(paths["route_lengths"])), 1)
            self.assertEqual(len(read_csv(paths["route_building_block_counts"])), 1)


if __name__ == "__main__":
    unittest.main()
