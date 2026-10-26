"""Missing plotting inputs must fail; intentionally empty distributions are valid."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from analysis import plots


class PlotInputTests(unittest.TestCase):
    def test_each_missing_csv_raises_before_replacing_an_existing_pdf(self):
        for missing in plots.METRIC_FILES.values():
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as directory:
                directory = Path(directory)
                for filename in plots.METRIC_FILES.values():
                    if filename != missing:
                        (directory / filename).write_text(
                            "method,run_dir,metric,value,count,proportion\n", encoding="utf-8"
                        )
                output = directory / "results_plot.pdf"
                output.write_bytes(b"existing figure")
                with patch.object(plots, "plot_distribution_pdf") as render:
                    with self.assertRaises(FileNotFoundError) as error:
                        plots.plot_distribution_csvs(directory, output)
                    self.assertEqual(Path(error.exception.filename), directory / missing)
                    render.assert_not_called()
                self.assertEqual(output.read_bytes(), b"existing figure")

    def test_missing_input_directory_creates_no_output_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            output = directory / "output" / "results_plot.pdf"
            with patch.object(plots, "plot_distribution_pdf") as render:
                with self.assertRaises(FileNotFoundError):
                    plots.plot_distribution_csvs(directory / "missing", output)
                render.assert_not_called()
            self.assertFalse(output.parent.exists())

    def test_header_only_reactants_csv_remains_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for metric, filename in plots.METRIC_FILES.items():
                contents = "method,run_dir,metric,value,count,proportion\n"
                if metric != "n_reactants":
                    contents += f"Example_mcts,50_top_results,{metric},2,1,1.000000\n"
                (directory / filename).write_text(contents, encoding="utf-8")
            output = directory / "results_plot.pdf"
            with patch.object(plots, "plot_distribution_pdf") as render:
                plots.plot_distribution_csvs(directory, output)
            data, actual_output = render.call_args.args
            self.assertEqual(actual_output, output)
            self.assertEqual(data["n_reactants"], {})
            self.assertEqual(data["route_length"], {"Example_mcts": {2: 1.0}})
            self.assertEqual(data["n_building_blocks"], {"Example_mcts": {2: 1.0}})


if __name__ == "__main__":
    unittest.main()
