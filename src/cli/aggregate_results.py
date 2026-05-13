#!/usr/bin/env python3
"""
CLI entry point for benchmark aggregation.

This script provides a command-line interface for aggregating benchmark results
from multiple method runs into a single CSV file.

Usage:
    python aggregate_results.py --results-root <path> --output-csv <path> [options]
    python aggregate_results.py --run-dir <path> [options]
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from analysis.aggregator import parse_args, main

if __name__ == "__main__":
    main()
