# Pesticide retrosynthesis benchmark

Reproduce the pesticide target dataset and analyse the archived Syntheseus searches
used in the manuscript. The verified workflow starts from the supplied raw CSV
and saved routes; rerunning live collection or model inference is optional.

Heavy assets are deposited at **[Zenodo: 10.5281/zenodo.23669375](https://doi.org/10.5281/zenodo.23669375)**.
The code, frozen preparation inputs, and lightweight reference outputs remain in
this repository. Run all commands below from the repository root.

## 1. Install the benchmark environment

Reproduction was verified on **Ubuntu 22.04.4 LTS, Linux 5.19.0-41-generic,
x86-64**, with Python 3.9.7, NumPy 1.26.4, pandas 2.1.4, and RDKit 2025.3.5.
The environment pins Syntheseus to commit
`569e9a9dad74217f01b0fb2bb1114f6a8894be95`.

```bash
conda env create --file environment.yml
conda activate syntheseus_benchmark
python -m pip check
```

This environment supports data preparation, analysis of **all** saved searches
(including MEGAN/Pistachio), and LocalRetro/RetroKNN inference. Preparation and
analysis require no GPU or model weights. GPU searches use PyTorch 1.10.2 and
CUDA 11.3 dependencies and require a compatible NVIDIA GPU and driver.
MEGAN inference needs a separate environment; see step 4.

## 2. Reproduce the benchmark targets from the frozen raw data

The required inputs are already in `data/`:

- `pesticides_raw.csv`: the archived BCPC/CIR/PubChem collection.
- `enamine_sep_lpdc_blocks.smi`: the Enamine EU September 2025 inventory used
  by the benchmark.

```bash
python src/cli/prepare_data.py
```

The command regenerates these three files in `data/`:

| File | Contents |
| --- | --- |
| `pesticides_update_filter_clusters.csv` | Final preparation pool with cluster assignments |
| `pesticides_update_filter_cluster_summary.csv` | All clusters ranked by population |
| `targets_no_bb_top100clust.smi` | The 100 benchmark centroids, in search-target order |

To reproduce into another directory while retaining the supplied references:

```bash
python src/cli/prepare_data.py --output-dir /tmp/pesticide-preparation
```

The workflow retains valid SMILES with **QED >= 0.5**, then applies Morgan
fingerprints (radius 2, 2048 bits) and Butina clustering (distance cutoff 0.4,
`reordering=False`). It ranks clusters by descending population and inspects the
first 100 centroids. If a selected centroid is in the inventory, every record
with that exact raw centroid SMILES is removed, and clustering/ranking repeats.
Other inventory members stay in the clustering pool. There are no additional
metal, short-chain, or aromatic-ring exclusion filters.

Inventory checks compare canonical centroid SMILES with the inventory lines as
supplied, matching Syntheseus's `canonicalize_inventory=False`. Source order,
raw SMILES spelling, and duplicate records are preserved: duplicates contribute
to cluster populations but the final 100 targets are distinct.

Expected counts are **2564 raw records -> 1590 QED-qualified records -> 1527 final
records -> 918 clusters -> 100 targets**. Four inventory-cleaning rounds remove
47, 8, 5, and 3 records. The final pool contains 1187 distinct canonical structures;
1527 counts records, not distinct molecules. The summary's first 100 rows match
the `.smi`; `centroid_row_index` refers to zero-based positions in the assignments
CSV, excluding its header.

All 100 targets, including their order, were reproduced exactly. The target-file
SHA-256 is `f9dcbfd502cc8dacd531bf7906f22543baba9ffc043fb4140ecf41ef644be622`.
Ties use the historical NumPy quicksort ordering, not a random seed. Keep the
pinned versions: the command rejects NumPy versions other than 1.26.4. Other
builds/CPU sorting implementations can order equal populations differently;
the regression tests check the exact published files.

### Optional: collect a new raw dataset

**Do not regenerate the raw CSV to reproduce the paper.** BCPC entries and
CIR/PubChem structure resolution can change, so a fresh collection can produce
different records, structures, clusters, and targets.

```bash
python src/cli/prepare_data.py --regenerate-raw --output-dir /tmp/pesticide-live
```

This explicitly invokes live BCPC collection with CAS-based CIR resolution and
PubChem InChIKey fallback, then runs preparation on that new collection. The
new `pesticides_raw.csv` and three derived files go into the chosen directory.
Without `--output-dir`, live regeneration uses `data/refreshed/`; using the
archived `data/` directory as its destination is refused. Normal preparation
never contacts these services. Live collection has been tested with mocked
services, not by collecting a new dataset for the reproduction checks.

## 3. Download and place the heavy assets

Download only the archives needed from
**[the Zenodo deposit](https://doi.org/10.5281/zenodo.23669375)**:

| Archive | Required location after extraction | Needed for |
| --- | --- | --- |
| `syntheseus_outputs_1_top_results.zip` | `syntheseus_outputs/1_top_results/` | Top-1 summaries and the published distribution figure; includes Pistachio MEGAN runs |
| `syntheseus_outputs_50_top_results.zip` | `syntheseus_outputs/50_top_results/` | Top-50 summaries and optional distributions |
| `MEGAN_backward.zip` | `MEGAN_backward/` | Rerunning the Pistachio-trained MEGAN model |

The archives preserve repository-relative paths. For example, after downloading
them to `~/Downloads`, extract the needed files **from the repository root**:

```bash
unzip "$HOME/Downloads/syntheseus_outputs_1_top_results.zip" -d .
unzip "$HOME/Downloads/syntheseus_outputs_50_top_results.zip" -d .
# Needed only for Pistachio MEGAN inference or the complete test suite:
unzip "$HOME/Downloads/MEGAN_backward.zip" -d .
```

The resulting layout must be:

```text
repository-root/
├── data/                              # included frozen inputs and reference CSVs
├── syntheseus_outputs/
│   ├── configs/                       # included search configurations
│   ├── 1_top_results/
│   │   └── METHOD_NAME/
│   │       ├── 0/                     # numeric target directories 0 ... 99
│   │       │   ├── stats.json
│   │       │   └── route_*.pkl         # saved routes for solved targets
│   │       └── ...
│   └── 50_top_results/                 # same method/target layout
├── MEGAN_backward/
│   ├── model_best.pt
│   ├── config.gin
│   └── megan_default_32_bfs_randat/
│       ├── all_actions.json
│       └── prop2oh.json
└── outputs/                           # included lightweight reference results
```

Avoid an extra archive-name directory around these folders. The complete MEGAN
model directory is needed for inference, including its configuration and action
vocabulary; `model_best.pt` alone is insufficient. Saved-route analysis needs
only the relevant campaign archive, not this model directory.

## 4. Inspect the search configurations or optionally rerun searches

**Skip inference and proceed to step 5 to reproduce the archived analyses.**
All configurations are in `syntheseus_outputs/configs/`. In the table below,
`<algorithm>` is either `mcts` or `retro_star`:

| Model/training data | Top-1 configuration | Top-50 configuration |
| --- | --- | --- |
| LocalRetro / USPTO-50K | `LocalRetro_<algorithm>.yml` | `LocalRetro_<algorithm>_top50.yml` |
| RetroKNN / USPTO-50K | `RetroKNN_<algorithm>.yml` | `RetroKNN_<algorithm>_top50.yml` |
| MEGAN / USPTO-50K | `MEGAN_<algorithm>.yml` | `MEGAN_<algorithm>_top50.yml` |
| MEGAN / Pistachio | `Pistachio_MEGAN_<algorithm>.yml` | Not part of this campaign |

Each uses `data/targets_no_bb_top100clust.smi`,
`data/enamine_sep_lpdc_blocks.smi`, `time_limit_s: 180`, and `save_graph: false`.
Top-1 and top-50 set `num_top_results` to 1 and 50 respectively; this controls
single-step predictions, not saved routes. Syntheseus's default
`num_routes_to_plot: 5` extracts up to five routes for each solved target.

### Checkpoints for every model

The non-Pistachio experiments use Syntheseus's default USPTO-50K bundles. The
pinned Syntheseus version automatically downloads them when absent:

| Model | Required bundle |
| --- | --- |
| LocalRetro | [Default USPTO-50K bundle](https://figshare.com/ndownloader/files/42287319): checkpoint, configuration, and template data |
| RetroKNN | [Default USPTO-50K bundle](https://figshare.com/ndownloader/files/45662430): LocalRetro backbone, adapter checkpoint, and nearest-neighbor datastore |
| MEGAN / USPTO-50K | [Default USPTO-50K bundle](https://figshare.com/ndownloader/files/42012732): checkpoint, configuration, and featurizer/action data |
| MEGAN / Pistachio | Complete `MEGAN_backward/` directory from the [Zenodo deposit](https://doi.org/10.5281/zenodo.23669375) |

Keep default USPTO bundles in a dedicated cache:

```bash
export SYNTHESEUS_CACHE_DIR="$PWD/.syntheseus_cache/uspto"
```

Otherwise Syntheseus uses `~/.cache/torch/syntheseus/`. Defaults are cached under
`LocalRetro_backward/`, `RetroKNN_backward/`, and `MEGAN_backward/` within that
cache. The Pistachio configurations instead specify `model_dir: MEGAN_backward`,
referring to the directory at the repository root. See the
[upstream model documentation](https://microsoft.github.io/syntheseus/0.5.0/single_step/)
and the [checkpoint links at the pinned commit](https://github.com/microsoft/syntheseus/blob/569e9a9dad74217f01b0fb2bb1114f6a8894be95/syntheseus/reaction_prediction/inference/default_checkpoint_links.yml).

For LocalRetro or RetroKNN, use `syntheseus_benchmark`. Examples redirect fresh
outputs to `/tmp/` to keep the archived campaigns separate:

```bash
syntheseus search --config syntheseus_outputs/configs/LocalRetro_mcts.yml \
  results_dir=/tmp/syntheseus-reruns/1_top_results/LocalRetro_mcts
syntheseus search --config syntheseus_outputs/configs/RetroKNN_retro_star_top50.yml \
  results_dir=/tmp/syntheseus-reruns/50_top_results/RetroKNN_retro_star
```

### MEGAN inference: create a separate environment

This applies to **both USPTO-50K and Pistachio MEGAN searches**. MEGAN's pinned
package requires TensorFlow 2.13.0, which requires NumPy <= 1.24.3. Benchmark
preparation requires NumPy 1.26.4, so keep these tasks in separate environments.
Saved MEGAN routes can still be analysed in `syntheseus_benchmark` without
TensorFlow or the MEGAN inference package.

Starting from the working benchmark environment, create a separate inference
copy and install the reference MEGAN dependencies:

```bash
conda create --name syntheseus_megan --clone syntheseus_benchmark
conda install --name syntheseus_megan numpy=1.24.3
conda activate syntheseus_megan
python -m pip install \
  "numpy==1.24.3" "tensorflow==2.13.0" "gin-config==0.3.0" \
  "GitPython==3.1.45" "torchtext==0.5.0" \
  "syntheseus[megan,viz] @ git+https://github.com/microsoft/syntheseus.git@569e9a9dad74217f01b0fb2bb1114f6a8894be95"
python -m pip check

export SYNTHESEUS_CACHE_DIR="$PWD/.syntheseus_cache/uspto"
syntheseus search --config syntheseus_outputs/configs/MEGAN_mcts.yml \
  results_dir=/tmp/syntheseus-reruns/1_top_results/MEGAN_mcts
syntheseus search --config syntheseus_outputs/configs/Pistachio_MEGAN_mcts.yml \
  results_dir=/tmp/syntheseus-reruns/1_top_results/Pistachio_MEGAN_mcts
```

The clone retains the benchmark's Python, RDKit, and CUDA/PyTorch versions.
[Upstream installation guidance](https://microsoft.github.io/syntheseus/0.5.0/installation/)
also describes installing a dedicated full Syntheseus environment. Use the pinned
commit above rather than replacing it with the latest release.

### Analyse a fresh search

Syntheseus adds a model/timestamp directory beneath `results_dir`. Aggregation
expects numeric target directories directly inside a method folder ending in
`_mcts` or `_retro_star`. Copy the numeric directories into that layout first.
For the LocalRetro example, replace `LocalRetro_TIMESTAMP` with the actual run
directory printed by Syntheseus:

```bash
conda activate syntheseus_benchmark
mkdir -p /tmp/syntheseus-analysis/1_top_results/LocalRetro_mcts
cp -a /tmp/syntheseus-reruns/1_top_results/LocalRetro_mcts/LocalRetro_TIMESTAMP/[0-9]* \
  /tmp/syntheseus-analysis/1_top_results/LocalRetro_mcts/
python src/cli/aggregate_results.py \
  --results-root /tmp/syntheseus-analysis/1_top_results \
  --output-csv /tmp/syntheseus-analysis/LocalRetro_top1_summary.csv
```

Fresh searches have not been rerun during the repository's reproduction checks.
A 180-second wall-clock budget makes their exact outcomes sensitive to hardware,
timing, and model/runtime differences. The verified results use the frozen
search outputs supplied through Zenodo; checkpoint downloads also depend on
availability of the upstream services. No model retraining is required here.

## 5. Reproduce the archived summaries and figure

After extracting the relevant campaign archives, switch back to the benchmark
environment and run:

```bash
conda activate syntheseus_benchmark
python src/cli/aggregate_results.py \
  --results-root syntheseus_outputs/1_top_results --plot-distributions
python src/cli/aggregate_results.py \
  --results-root syntheseus_outputs/50_top_results
```

The top-1 command writes `outputs/1_top_results_summary.csv`,
`outputs/results_plot.pdf`, and these figure-source tables in `outputs/`:
`route_lengths_distribution.csv`, `building_blocks_distribution.csv`, and
`reactants_distribution.csv`. The top-50 command writes
`outputs/50_top_results_summary.csv`. Summary CSVs are calculated directly from
`stats.json` and route pickles; the distribution CSVs are used for the figure.

To generate top-50 distributions too, add `--plot-distributions` to the second
command. These CSVs and its PDF go into `outputs/50_top_results/`, preserving the
top-1 figure. MEGAN/RetroKNN top-50 MCTS routes have incomplete reactant inference;
analysis warns and omits the affected method's entire reactant distribution.
Route lengths, building-block distributions, and summaries remain available.
Plotting requires all three distribution files to exist; header-only distributions
are accepted when counts are intentionally omitted.

## 6. Verify reproduction

With **all three Zenodo archives extracted**, run:

```bash
python -B -m pytest
```

The 45 tests cover preparation, the exact target file/order, both summary CSVs,
distribution tables, route metrics, configuration/input consistency, and error
handling. The complete suite checks the presence of the Pistachio model files
even though saved-route analysis itself does not need the checkpoint.
The figure comparison checks rendered PDF content, rather than PDF metadata.
It requires the external `pdftoppm` tool (`poppler-utils` on Ubuntu) and is skipped
if that tool is absent. Live collection is mocked and searches are not rerun.

## Additional reference

- **Architecture:** `src/cli/` contains preparation, aggregation, and molecule-grid
  entry points; `src/data_prep/` contains reusable preparation functions;
  `src/analysis/` extracts route metrics and plots distributions; `src/config.py`
  supplies analysis defaults. `data/` and `outputs/` retain lightweight reference
  artifacts; heavy model/search assets are restored from Zenodo.
- **Metric definitions:** Retro* route length counts all reaction nodes, including
  branches; MCTS route length counts state transitions. Average route counts
  describe saved route files among solved targets. Building blocks count unique
  canonical purchasable structures per route. The historical CSV column
  `run_dir` identifies the campaign, such as `1_top_results`.
- **Additional commands:** `python src/cli/sample_molecules.py --seed 0` writes
  solved/unsolved molecule grids to `outputs/molecule_samples/`. CLI `--help`
  lists overrides, including `--method-dir`, `--output-csv`, `--plot-data-dir`,
  `--write-plot-data`, and preparation's `--raw-csv`/`--inventory-smiles`.
  Analysis paths can also be overridden with
  `SYNTHSEUS_BENCHMARK_RESULTS_ROOT` and `SYNTHSEUS_BENCHMARK_OUTPUT_ROOT`.
