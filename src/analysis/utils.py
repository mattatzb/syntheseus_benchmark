#!/usr/bin/env python3
"""
Utility functions for benchmark analysis.

This module provides common helper functions used across the analysis pipeline,
including number validation, statistics, formatting, and file discovery.
"""

from __future__ import annotations

import math
import rdkit
from pathlib import Path
from statistics import mean
from typing import Iterator
from rdkit import Chem

def is_finite_number(value: object) -> bool:
    """
    Check if a value is a finite number (int or float).

    Args:
        value: The value to check.

    Returns:
        True if value is int or float and is finite (not NaN or inf), False otherwise.
    """
    return isinstance(value, (int, float)) and math.isfinite(value)


def mean_or_nan(values: list[float]) -> float:
    """
    Calculate the mean of a list of values, or return NaN if list is empty.

    Args:
        values: List of numeric values.

    Returns:
        Mean of the values if non-empty, otherwise math.nan.
    """
    return mean(values) if values else math.nan


def format_float(value: object, digits: int = 6) -> str:
    """
    Format a numeric value as a string with specified decimal places.

    Args:
        value: The value to format.
        digits: Number of decimal places (default: 6).

    Returns:
        Formatted string if value is finite number, empty string otherwise.
    """
    return f"{float(value):.{digits}f}" if is_finite_number(value) else ""


def iter_target_dirs(run_dir: Path) -> Iterator[Path]:
    """
    Iterate over target directories in a run, in sorted order.

    Target directories are those with numeric names.

    Args:
        run_dir: Path to the run directory.

    Yields:
        Path to each numeric-named subdirectory (target).
    """
    for target_dir in sorted(run_dir.iterdir(), key=lambda p: p.name):
        if target_dir.is_dir() and target_dir.name.isdigit():
            yield target_dir


def shorten_notes(notes: list[str], limit: int = 5) -> str:
    """
    Shorten a list of notes to a limited number, with ellipsis if truncated.

    Args:
        notes: List of note strings.
        limit: Maximum number of notes to include (default: 5).

    Returns:
        Joined note string with ellipsis if notes exceed limit, or empty string if no notes.
    """
    if not notes:
        return ""
    short = "; ".join(notes[:limit])
    if len(notes) > limit:
        short += f"; ... +{len(notes) - limit} more"
    return short


def canonicalize_smiles(smiles: str | None) -> str | None:
    """Return RDKit canonical SMILES, or None for empty/unparseable input."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)