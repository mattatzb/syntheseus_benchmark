#!/usr/bin/env python3
"""
Configuration defaults for the Syntheseus benchmark repository.

This module centralizes repository paths and supports environment variable
overrides for portability across different machines and directory layouts.

Modify as needed to fit your local setup, or use environment variables to override defaults.
"""

from __future__ import annotations

import os
from pathlib import Path


def _path_from_env(name: str, default: Path) -> Path:
    value = os.environ.get(name)
    return Path(value) if value else default


class Config:
    PROJECT_ROOT: Path = Path(__file__).resolve().parents[1]
    SRC_ROOT: Path = Path(__file__).resolve().parent
    RESULTS_ROOT: Path = _path_from_env(
        "SYNTHSEUS_BENCHMARK_RESULTS_ROOT",
        PROJECT_ROOT / "syntheseus_outputs" / "1_top_results",
    )
    OUTPUT_ROOT: Path = _path_from_env(
        "SYNTHSEUS_BENCHMARK_OUTPUT_ROOT",
        PROJECT_ROOT / "outputs",
    )
    DATA_DIR: Path = _path_from_env(
        "SYNTHSEUS_BENCHMARK_DATA_DIR",
        PROJECT_ROOT / "data",
    )
    INVENTORY_SMILES_FILE: Path = _path_from_env(
        "SYNTHSEUS_BENCHMARK_INVENTORY_SMILES_FILE",
        DATA_DIR / "enamine_sep_lpdc_blocks.smi",
    )
    DEFAULT_OUTPUT_CSV: Path = OUTPUT_ROOT / "results_summary.csv"
    DEFAULT_SAMPLE_DIR: Path = OUTPUT_ROOT / "molecule_samples"
    DEFAULT_OUTPUT_PLOT: Path = OUTPUT_ROOT / "results_plot.pdf"

    @classmethod
    def validate_path(cls, path: Path, description: str) -> None:
        """
        Validate that a critical path exists.

        Args:
            path: Path to validate.
            description: Human-readable description of the path.

        Raises:
            FileNotFoundError: If the given path does not exist.
        """
        if not path.exists():
            raise FileNotFoundError(f"Missing required path for {description}: {path}")

    @classmethod
    def validate(cls) -> None:
        """
        Validate the most important repository paths.

        This does not validate every possible path, only the paths that are
        usually required for analysis and data processing.
        """
        cls.validate_path(cls.RESULTS_ROOT, "results root")
        cls.validate_path(cls.OUTPUT_ROOT, "output root")
        cls.validate_path(cls.DATA_DIR, "data directory")
