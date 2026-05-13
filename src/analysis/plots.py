#!/usr/bin/env python3
"""
Plot distribution summaries for Syntheseus benchmark routes.

This module consumes the distribution CSVs produced by aggregator.py and writes a
single-page PDF with MCTS and retro-star experiments as two labeled sections.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Iterable

try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ModuleNotFoundError:
    plt = None
    PdfPages = None


MetricData = dict[str, dict[str, dict[int, float]]]

METRIC_FILES = {
    "route_length": "route_lengths_distribution.csv",
    "n_building_blocks": "building_blocks_distribution.csv",
    "n_reactants": "reactants_distribution.csv",
}
METRIC_TITLES = {
    "route_length": "# Steps",
    "n_building_blocks": "# Building Blocks",
    "n_reactants": "# Reactants",
}
ALGORITHMS = {
    "mcts": "MCTS",
    "retro_star": "Retro Star",
}
ALGORITHM_LABELS = {
    "mcts": "a",
    "retro_star": "b",
}


def _require_matplotlib() -> None:
    if plt is None or PdfPages is None:
        raise RuntimeError(
            "Plotting requires matplotlib. Install it in the active environment "
            "before creating route distribution PDFs."
        )


def _read_distribution_csv(input_csv: Path) -> dict[str, dict[int, float]]:
    """Read one distribution CSV into method -> value -> proportion."""
    distribution: dict[str, dict[int, float]] = defaultdict(dict)
    if not input_csv.exists():
        return distribution

    with input_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            method = row.get("method", "")
            value = row.get("value", "")
            proportion = row.get("proportion", "")
            if not method or not value or not proportion:
                continue
            distribution[method][int(value)] = float(proportion)

    return distribution


def read_distribution_data(plot_data_dir: Path) -> MetricData:
    """
    Read all distribution CSVs from a directory.

    Args:
        plot_data_dir: Directory containing the distribution CSVs.

    Returns:
        Mapping metric name -> method -> value -> proportion.
    """
    return {
        metric: _read_distribution_csv(plot_data_dir / filename)
        for metric, filename in METRIC_FILES.items()
    }


def _method_algorithm(method: str) -> str | None:
    if method.endswith("_retro_star"):
        return "retro_star"
    if method.endswith("_mcts"):
        return "mcts"
    return None


def _display_method(method: str) -> str:
    for suffix in ("_retro_star", "_mcts"):
        if method.endswith(suffix):
            method = method[: -len(suffix)]
            break
    return method.replace("_", " ")


def _methods_for_algorithm(distribution_data: MetricData, algorithm: str) -> list[str]:
    """Return methods for one algorithm, preserving first-seen metric order."""
    methods: list[str] = []
    seen: set[str] = set()
    for metric in METRIC_FILES:
        for method in distribution_data.get(metric, {}):
            if method in seen:
                continue
            if _method_algorithm(method) != algorithm:
                continue
            seen.add(method)
            methods.append(method)
    return methods


def _plasma_colors(methods: Iterable[str]) -> dict[str, tuple[float, float, float, float]]:
    methods = list(methods)
    cmap = plt.get_cmap("plasma")
    if len(methods) <= 1:
        positions = [0.55]
    else:
        positions = [0.15 + (0.7 * i / (len(methods) - 1)) for i in range(len(methods))]
    return {method: cmap(position) for method, position in zip(methods, positions)}


def _axis_values(metric_distribution: dict[str, dict[int, float]], methods: list[str]) -> list[int]:
    values = sorted(
        {
            value
            for method in methods
            for value in metric_distribution.get(method, {})
        }
    )
    return values


def _tick_step(max_value: int) -> int:
    if max_value <= 8:
        return 1
    if max_value <= 20:
        return 5
    return 10


def _plot_metric_group(
    subfig,
    metric: str,
    metric_distribution: dict[str, dict[int, float]],
    methods: list[str],
    colors: dict[str, tuple[float, float, float, float]],
) -> None:
    axes = subfig.subplots(2, 2, sharey=True)
    subfig.suptitle(METRIC_TITLES[metric], fontsize=13)

    values = _axis_values(metric_distribution, methods)
    max_value = max(values) if values else 1
    tick_step = _tick_step(max_value)
    xticks = list(range(0, max_value + 1, tick_step))
    if max_value not in xticks:
        xticks.append(max_value)

    max_proportion = max(
        (
            proportion
            for method in methods
            for proportion in metric_distribution.get(method, {}).values()
        ),
        default=0.0,
    )
    y_max = max(0.1, max_proportion * 1.15)

    flat_axes = list(axes.flat)
    for ax, method in zip(flat_axes, methods):
        method_distribution = metric_distribution.get(method, {})
        xs = sorted(method_distribution)
        ys = [method_distribution[x] for x in xs]
        ax.bar(xs, ys, width=0.75, color=colors[method], edgecolor="black", linewidth=0.6)
        ax.set_title(_display_method(method), fontsize=11)
        ax.set_xlim(-0.5, max_value + 0.5)
        ax.set_ylim(0, y_max)
        ax.set_xticks(xticks)
        ax.tick_params(axis="both", labelsize=9)

    for ax in flat_axes[len(methods):]:
        ax.axis("off")

    for row_index in range(2):
        axes[row_index, 0].set_ylabel("Proportion", fontsize=10)


def plot_distribution_pdf(distribution_data: MetricData, output_pdf: Path) -> None:
    """
    Write the route distribution PDF.

    Args:
        distribution_data: Distribution data from read_distribution_data.
        output_pdf: PDF file to create.
    """
    _require_matplotlib()
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    all_methods = [
        method
        for metric_distribution in distribution_data.values()
        for method in metric_distribution
    ]
    colors = _plasma_colors(dict.fromkeys(all_methods))

    with PdfPages(output_pdf) as pdf:
        fig = plt.figure(figsize=(16, 11.2), constrained_layout=True)
        algorithm_subfigs = fig.subfigures(2, 1, hspace=0.08)

        for algorithm_subfig, (algorithm, algorithm_title) in zip(
            algorithm_subfigs,
            ALGORITHMS.items(),
        ):
            methods = _methods_for_algorithm(distribution_data, algorithm)
            if len(methods) > 4:
                raise ValueError(
                    f"Expected at most 4 methods for {algorithm}, found {len(methods)}."
                )

            algorithm_label = ALGORITHM_LABELS[algorithm]
            algorithm_subfig.suptitle(
                f"{algorithm_label})",
                fontsize=15,
                x=0.01,
                ha="left",
            )
            metric_subfigs = algorithm_subfig.subfigures(1, 3, wspace=0.05)

            for metric_subfig, metric in zip(metric_subfigs, METRIC_FILES):
                _plot_metric_group(
                    metric_subfig,
                    metric,
                    distribution_data.get(metric, {}),
                    methods,
                    colors,
                )

        pdf.savefig(fig)
        plt.close(fig)


def plot_distribution_csvs(plot_data_dir: Path, output_pdf: Path) -> None:
    """
    Read distribution CSVs and write the route distribution PDF.

    Args:
        plot_data_dir: Directory containing distribution CSVs.
        output_pdf: PDF file to create.
    """
    distribution_data = read_distribution_data(plot_data_dir)
    plot_distribution_pdf(distribution_data, output_pdf)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot route-length, building-block, and reactant distributions."
    )
    parser.add_argument(
        "--plot-data-dir",
        type=Path,
        required=True,
        help="Directory containing distribution CSVs from aggregate_results.py.",
    )
    parser.add_argument(
        "--output-pdf",
        type=Path,
        default=None,
        help="Output PDF path. Defaults to route_distributions.pdf in --plot-data-dir.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_pdf = args.output_pdf or (args.plot_data_dir / "route_distributions.pdf")
    plot_distribution_csvs(args.plot_data_dir, output_pdf)
    print(f"Wrote plots to {output_pdf}")


if __name__ == "__main__":
    main()
