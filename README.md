# syntheseus_benchmark

Post-processing and analysis for Syntheseus benchmark outputs.

This repository does not rerun retrosynthesis experiments. It reads existing benchmark output folders,
extracts summary metrics, and optionally renders solved/unsolved molecules.

## Current Layout

The repository now uses a source-based CLI and analysis package in `src/`:

- `src/cli/aggregate_results.py` - generate aggregated CSV summaries
- `src/cli/sample_molecules.py` - sample and visualize solved/unsolved molecules
- `src/analysis/aggregator.py` - CSV aggregation logic
- `src/analysis/route_processor.py` - target-level metric extraction
- `src/analysis/retro_star_extractor.py` - route pickle extractor for reaction-node routes
- `src/analysis/mcts_extractor.py` - route pickle extractor for MCTS state-node routes
- `src/analysis/plots.py` - route distribution plot generation
- `src/analysis/visualizer.py` - RDKit sampling and grid image generation
- `src/config.py` - centralized path defaults and environment overrides
- `outputs/*.csv` - example CSV and plot outputs generated with `src/cli/aggregate_results.py`
- `outputs/molecule_samples/` - example image outputs generated with `src/cli/sample_molecules.py`

The example files in `outputs/` were generated from Syntheseus result folders under
`syntheseus_outputs/`.

## Expected Results Structure

The current scripts expect the benchmark output root to contain method folders, where each method folder contains numeric target folders:

```text
syntheseus_outputs/
  1_top_results/
    METHOD_NAME/
      0/
        stats.json
        route_*.pkl
      1/
      ...
      99/
```

## Configuration

Path defaults are managed in `src/config.py`. The current defaults are:

- `RESULTS_ROOT`: `PROJECT_ROOT / "syntheseus_outputs" / "1_top_results"`
- `OUTPUT_ROOT`: `PROJECT_ROOT / "outputs" / "test_outputs"`
- `DATA_DIR`: `PROJECT_ROOT / "data"`
- `INVENTORY_SMILES_FILE`: `DATA_DIR / "enamine_sep_lpdc_blocks.smi"`
- `DEFAULT_OUTPUT_CSV`: `OUTPUT_ROOT / "results_summary.csv"`
- `DEFAULT_SAMPLE_DIR`: `OUTPUT_ROOT / "molecule_samples"`
- `DEFAULT_OUTPUT_PLOT`: `OUTPUT_ROOT / "results_plot.pdf"`

### Environment variable overrides

You can override defaults without editing code by exporting environment variables:

```bash
export SYNTHSEUS_BENCHMARK_RESULTS_ROOT=/path/to/results_root
export SYNTHSEUS_BENCHMARK_OUTPUT_ROOT=/path/to/output_root
export SYNTHSEUS_BENCHMARK_DATA_DIR=/path/to/data
export SYNTHSEUS_BENCHMARK_INVENTORY_SMILES_FILE=/path/to/inventory.smi
```

Then run the CLI normally.

### When to edit `src/config.py`

If you want to persistently change the repository defaults for your local setup,
editing `src/config.py` is often easier than passing the same CLI arguments every time.
For example, if you always work with a different results root, update `RESULTS_ROOT`
and `DEFAULT_OUTPUT_CSV` once in `src/config.py`.

> `INVENTORY_SMILES_FILE` is kept in config for inventory-related analysis and future use.

## Usage

From the repository root, use the CLI entrypoints in `src/cli`.

### Aggregate benchmark CSV

```bash
python src/cli/aggregate_results.py
```

This scans all methods under the configured results root and writes a summary CSV to the configured analysis output.
The summary CSV columns are unchanged when plot data is requested.

#### Common options

- `--results-root PATH`
  - Root directory containing method folders.
- `--output-csv PATH`
  - Output CSV file path.
- `--show-inventory-columns`
  - Add inventory-related columns:
    - `n_inventory_excluded`
    - `inventory_true_pct`
- `--show-debug-columns`
  - Add debug columns:
    - `n_targets`
    - `n_solved`
    - `note`
- `--method-dir PATH`
  - Restrict analysis to one or more specified method directories.
- `--merge-with-existing`
  - Merge updated rows into an existing CSV rather than replacing it.
- `--run-dir PATH`
  - Specify individual run directories explicitly.
- `--write-plot-data`
  - Write distribution CSVs for route length, route building blocks, and reaction-step reactant counts.
- `--plot-distributions`
  - Write the distribution CSVs and generate the configured distribution plot PDF.
- `--plot-data-dir PATH`
  - Directory for the distribution CSVs.

### Sample solved/unsolved molecules

```bash
python src/cli/sample_molecules.py
```

This samples solved and unsolved targets from each method's latest benchmark data and generates RDKit image grids.

#### Common options

- `--results-root PATH`
  - Root directory containing method folders.
- `--n-solved INT`
  - Number of solved molecules sampled per method.
- `--n-unsolved INT`
  - Number of unsolved molecules sampled per method.
- `--seed INT`
  - Random seed for reproducible sampling.
- `--output-dir PATH`
  - Output directory for generated PNG images.

## Examples

### Generate a new summary CSV

```bash
python src/cli/aggregate_results.py
```

### Generate a summary with inventory and debug columns

```bash
python src/cli/aggregate_results.py --show-inventory-columns --show-debug-columns
```

### Use a different results root

```bash
python src/cli/aggregate_results.py --results-root syntheseus_outputs/50_top_results
```

### Save results to a custom CSV path

```bash
python src/cli/aggregate_results.py --output-csv outputs/50_top_results_summary.csv
```

### Generate summary CSV and distribution plots

```bash
python src/cli/aggregate_results.py --plot-distributions
```

This writes the normal summary CSV, the three distribution CSVs, and the PDF configured by `DEFAULT_OUTPUT_PLOT`.

### Generate only plot-data CSVs

```bash
python src/cli/aggregate_results.py --write-plot-data
```

### Sample molecules with a custom output directory

```bash
python src/cli/sample_molecules.py --output-dir outputs/samples
```

## Notes

- `src/config.py` centralizes the most important path defaults.
- MCTS routes are stored as state-node paths rather than explicit reaction-node paths. Route length and final building-block counts are extracted from the MCTS state sequence. Reactant counts for MCTS are inferred from newly introduced molecules between consecutive states.
- Retro-star routes use explicit reaction SMILES, so reaction-step reactant counts are extracted directly from `rxn_smiles`.
- For persistent local customization, editing `src/config.py` is often easier than repeating long CLI arguments.
- For temporary changes, use the environment variables listed above.
