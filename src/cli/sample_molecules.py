#!/usr/bin/env python3
"""
CLI entry point for molecule visualization.

This script provides a command-line interface for sampling and visualizing solved/unsolved
molecules from benchmark results using RDKit grid images.

Usage:
    python sample_molecules.py --results-root <path> --output-dir <path> [options]
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from analysis.visualizer import parse_args, main

if __name__ == "__main__":
    main()
