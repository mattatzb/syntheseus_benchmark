# Pesticide retrosynthesis benchmark

Reproduce the manuscript's target dataset, summary tables, and figure from the
archived inputs. Run the commands below from the repository root.

## 1. Install the environment

Verified on **Ubuntu 22.04.4 LTS, Linux 5.19.0-41-generic, x86-64**.
Dependencies and the Syntheseus revision are pinned in `environment.yml`.

```bash
conda env create --file environment.yml
conda activate syntheseus_benchmark
```

Preparation and saved-output analysis require no GPU or model checkpoints.
LocalRetro/RetroKNN searches use this environment, with CUDA 11.3 for GPU
execution. MEGAN searches require a separate environment (see below).

## 2. Reproduce data preparation

The frozen raw collection and Enamine September 2025 inventory are included as
`data/pesticides_raw.csv` and `data/enamine_sep_lpdc_blocks.smi`.

```bash
python src/cli/prepare_data.py
```

This regenerates, in `data/`:

- `pesticides_update_filter_clusters.csv`
- `pesticides_update_filter_cluster_summary.csv`
- `targets_no_bb_top100clust.smi`

## 3. Download the archived search outputs

Download the required ZIPs from **[Zenodo: 10.5281/zenodo.23669375](https://doi.org/10.5281/zenodo.23669375)**
and extract them into the repository root, preserving these paths:

| Archive | Extracted directory | Purpose |
| --- | --- | --- |
| `syntheseus_outputs_1_top_results.zip` | `syntheseus_outputs/1_top_results/` | Top-1 tables and the manuscript figure |
| `syntheseus_outputs_50_top_results.zip` | `syntheseus_outputs/50_top_results/` | Top-50 table |
| `MEGAN_backward.zip` | `MEGAN_backward/` | Pistachio MEGAN inference and the complete test suite |

For example: `unzip ~/Downloads/syntheseus_outputs_1_top_results.zip -d .`.
The checkpoint archive is unnecessary for analysing saved searches.

## 4. Reproduce the tables and figure

```bash
python src/cli/aggregate_results.py --results-root syntheseus_outputs/1_top_results --plot-distributions
python src/cli/aggregate_results.py --results-root syntheseus_outputs/50_top_results
```

These write `outputs/1_top_results_summary.csv`,
`outputs/50_top_results_summary.csv`, and `outputs/results_plot.pdf`.
The figure also uses the generated `route_lengths_distribution.csv`,
`building_blocks_distribution.csv`, and `reactants_distribution.csv` in
`outputs/`; summary generation does not depend on these three files.

## 5. Check reproduction

With all three Zenodo archives extracted:

```bash
python -B -m pytest
```

Tests compare regenerated preparation and analysis results with the supplied
references. PDF comparison additionally requires `pdftoppm` (`poppler-utils`
on Ubuntu) and is skipped when that tool is absent.

## Optional: rerun Syntheseus searches

Configurations are in `syntheseus_outputs/configs/`:

- `<Model>_<algorithm>.yml`: top-1, for LocalRetro, RetroKNN, or MEGAN.
- `<Model>_<algorithm>_top50.yml`: corresponding top-50 runs.
- `Pistachio_MEGAN_<algorithm>.yml`: Pistachio-trained MEGAN, top-1 only.

`<algorithm>` is `mcts` or `retro_star`. All use the benchmark target file,
Enamine inventory, and a **180-second** limit. `num_top_results` controls
single-step predictions (1 or 50), not saved routes (default maximum: five).

LocalRetro, RetroKNN, and standard MEGAN use Syntheseus's default **USPTO-50K
checkpoint bundles**, downloaded automatically. Their download links are listed
in the [pinned Syntheseus checkpoint configuration](https://github.com/microsoft/syntheseus/blob/569e9a9dad74217f01b0fb2bb1114f6a8894be95/syntheseus/reaction_prediction/inference/default_checkpoint_links.yml).
Pistachio configurations use `model_dir: MEGAN_backward`: download the **complete
folder** from Zenodo, not just its checkpoint file.

Example, using a dedicated default-model cache and a separate rerun destination:

```bash
export SYNTHESEUS_CACHE_DIR="$PWD/.syntheseus_cache/uspto"
syntheseus search --config syntheseus_outputs/configs/LocalRetro_mcts.yml \
  results_dir=/tmp/syntheseus-reruns/1_top_results/LocalRetro_mcts
```

Fresh searches produce a model/timestamp subdirectory. For analysis, copy its
numeric target directories (`0` through `99`) directly into a method folder
ending in `_mcts` or `_retro_star`, then select it with aggregation's
`--method-dir` option.

### MEGAN environment

**MEGAN's TensorFlow 2.13.0 requires NumPy <= 1.24.3; preparation requires
NumPy 1.26.4.** Keep preparation and analysis in `syntheseus_benchmark` and run
both USPTO-50K and Pistachio MEGAN inference in a separate environment.

<details>
<summary>Create the separate MEGAN environment</summary>

```bash
conda create --name syntheseus_megan --clone syntheseus_benchmark
conda install --name syntheseus_megan numpy=1.24.3
conda activate syntheseus_megan
python -m pip install \
  "numpy==1.24.3" "tensorflow==2.13.0" "gin-config==0.3.0" \
  "GitPython==3.1.45" "torchtext==0.5.0" \
  "syntheseus[megan,viz] @ git+https://github.com/microsoft/syntheseus.git@569e9a9dad74217f01b0fb2bb1114f6a8894be95"
python -m pip check
```

Activate `syntheseus_benchmark` again before preparation or analysis.

</details>

## Reproducibility limits and implementation notes

- **Raw data:** use the supplied snapshot. BCPC, CIR, and PubChem evolve;
  `python src/cli/prepare_data.py --regenerate-raw --output-dir /tmp/pesticide-live`
  collects new data and is not expected to reproduce the paper.
- **Cluster ties:** historical quicksort ordering depends on NumPy/build/CPU
  behavior. Keep the pinned environment, including NumPy 1.26.4; tests check the
  exact 100-target file.
- **Search reruns:** not rerun during these reproduction checks. Hardware,
  timing, and runtime differences can change results under the wall-clock limit.
- **Top-50 reactants:** analysis warns and omits unreliable reactant distributions
  for MEGAN/RetroKNN MCTS. Their summaries and other distributions remain available.
  Optional top-50 plots use `--plot-distributions` and go into `outputs/50_top_results/`.
- **Architecture and metrics:** `src/cli/` provides entry points, `src/data_prep/`
  prepares targets, and `src/analysis/` computes metrics/plots. Retro* route length
  counts all reaction nodes; MCTS counts state transitions. The legacy `run_dir`
  CSV column identifies the campaign. See each CLI's `--help` for further options.
