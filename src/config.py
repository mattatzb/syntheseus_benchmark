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
    DEFAULT_SAMPLE_DIR: Path = OUTPUT_ROOT / "molecule_samples"

    @staticmethod
    def campaign_for(path: Path) -> str | None:
        """Identify a bundled campaign, including paths inside a method directory."""
        resolved = path.resolve()
        for candidate in (resolved, *resolved.parents):
            if candidate.name in {"1_top_results", "50_top_results"}:
                return candidate.name
        return None

    @classmethod
    def summary_csv_for(cls, campaign: str | None) -> Path:
        filename = f"{campaign}_summary.csv" if campaign else "benchmark_summary.csv"
        return cls.OUTPUT_ROOT / filename

    @classmethod
    def plot_data_dir_for(cls, campaign: str | None) -> Path:
        # Preserve the published top-one paths; isolate the other campaign.
        if campaign == "50_top_results":
            return cls.OUTPUT_ROOT / campaign
        return cls.OUTPUT_ROOT
