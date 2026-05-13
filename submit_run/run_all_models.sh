#!/usr/bin/env bash
set -u -o pipefail

BASE_CONFIG="${1:-/home/tatzber/syntheseus_benchmark/test_run.yml}"
RESULTS_ROOT="${2:-/home/tatzber/syntheseus_benchmark/results_all_models}"
DRY_RUN="${DRY_RUN:-0}"
TARGETS_CSV="${TARGETS_CSV:-/home/tatzber/syntheseus_benchmark/pesticides_qed_ge_0_5_cluster_summary.csv}"
TARGETS_MODE="${TARGETS_MODE:-centroids}"  # centroids | sample
SAMPLE_SIZE="${SAMPLE_SIZE:-50}"
RANDOM_SEED="${RANDOM_SEED:-42}"
RESUME="${RESUME:-1}"
if [[ "${3:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

MODELS=(
  "MEGAN"
  "Chemformer"
  "Graph2Edits"
  "LocalRetro"
  "MHNreact"
  "RetroKNN"
  "RootAligned"
)

SEARCH_ALGOS=(
  "retro_star"
  "mcts"
)

if [[ ! -f "${BASE_CONFIG}" ]]; then
  echo "Base config not found: ${BASE_CONFIG}" >&2
  exit 1
fi
if [[ ! -f "${TARGETS_CSV}" ]]; then
  echo "Targets CSV not found: ${TARGETS_CSV}" >&2
  exit 1
fi

if command -v syntheseus >/dev/null 2>&1; then
  RUNNER=(syntheseus)
else
  CONDA_ENV_NAME="${CONDA_ENV_NAME:-syntheseus-full}"
  if command -v conda >/dev/null 2>&1; then
    RUNNER=(conda run -n "${CONDA_ENV_NAME}" syntheseus)
  else
    echo "Neither 'syntheseus' nor 'conda' was found in PATH." >&2
    exit 1
  fi
fi

mkdir -p "${RESULTS_ROOT}"
RUN_TAG="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${RESULTS_ROOT}/run_${RUN_TAG}.log"
TARGETS_SMI_FILE="${RESULTS_ROOT}/targets_${TARGETS_MODE}_${RUN_TAG}.smi"
CENTROID_MAP_FILE="${RESULTS_ROOT}/centroid_targets_${RUN_TAG}.tsv"

python - <<PY
import pandas as pd
from pathlib import Path

targets_csv = Path("${TARGETS_CSV}")
output_smi = Path("${TARGETS_SMI_FILE}")
centroid_map_file = Path("${CENTROID_MAP_FILE}")
targets_mode = "${TARGETS_MODE}"
sample_size = int("${SAMPLE_SIZE}")
seed = int("${RANDOM_SEED}")

df = pd.read_csv(targets_csv)
if targets_mode == "centroids":
    required_cols = {"cluster_id", "cluster_size", "centroid_smiles"}
    if not required_cols.issubset(df.columns):
        raise SystemExit(
            "TARGETS_MODE=centroids requires columns "
            f"{sorted(required_cols)} in {targets_csv}"
        )
    cdf = df[["cluster_id", "cluster_size", "centroid_smiles"]].copy()
    cdf["centroid_smiles"] = cdf["centroid_smiles"].astype(str)
    cdf = cdf[(cdf["centroid_smiles"] != "") & (cdf["centroid_smiles"] != "Failed")]
    cdf = cdf.drop_duplicates(subset=["cluster_id"], keep="first")
    cdf = cdf.sort_values(by="cluster_id", ascending=True).reset_index(drop=True)
    if len(cdf) == 0:
        raise SystemExit(f"No usable centroid targets found in {targets_csv}")

    output_smi.write_text("\\n".join(cdf["centroid_smiles"].tolist()) + "\\n")
    cdf.to_csv(centroid_map_file, sep="\\t", index=False)
    print(f"Wrote centroid map ({len(cdf)} clusters) -> {centroid_map_file}")
    print(f"Wrote all centroids ({len(cdf)}) -> {output_smi}")
elif targets_mode == "sample":
    if "smiles" not in df.columns:
        raise SystemExit(f"TARGETS_MODE=sample requires 'smiles' column in {targets_csv}")
    smiles_df = pd.DataFrame({"smiles": df["smiles"]}).dropna()
    smiles_df["smiles"] = smiles_df["smiles"].astype(str)
    smiles_df = smiles_df[(smiles_df["smiles"] != "") & (smiles_df["smiles"] != "Failed")]
    smiles_df = smiles_df.drop_duplicates(subset=["smiles"])
    if len(smiles_df) == 0:
        raise SystemExit(f"No usable SMILES found in {targets_csv}")

    n = min(sample_size, len(smiles_df))
    sampled = smiles_df.sample(n=n, random_state=seed)
    output_smi.write_text("\\n".join(sampled["smiles"].tolist()) + "\\n")
    print(f"Sampled {n} SMILES -> {output_smi}")
else:
    raise SystemExit("TARGETS_MODE must be one of: centroids, sample")
PY

echo "Base config : ${BASE_CONFIG}" | tee -a "${LOG_FILE}"
echo "Results root: ${RESULTS_ROOT}" | tee -a "${LOG_FILE}"
echo "Runner      : ${RUNNER[*]}" | tee -a "${LOG_FILE}"
echo "Targets CSV : ${TARGETS_CSV}" | tee -a "${LOG_FILE}"
echo "Targets mode: ${TARGETS_MODE}" | tee -a "${LOG_FILE}"
if [[ "${TARGETS_MODE}" == "centroids" ]]; then
  echo "Centroid map: ${CENTROID_MAP_FILE}" | tee -a "${LOG_FILE}"
fi
if [[ "${TARGETS_MODE}" == "sample" ]]; then
  echo "Sample size : ${SAMPLE_SIZE}" | tee -a "${LOG_FILE}"
  echo "Seed        : ${RANDOM_SEED}" | tee -a "${LOG_FILE}"
fi
echo "Targets SMI : ${TARGETS_SMI_FILE}" | tee -a "${LOG_FILE}"
echo "Dry run     : ${DRY_RUN}" | tee -a "${LOG_FILE}"
echo "Resume mode : ${RESUME}" | tee -a "${LOG_FILE}"
echo "" | tee -a "${LOG_FILE}"

success_count=0
fail_count=0

if [[ "${TARGETS_MODE}" == "centroids" ]]; then
  while IFS=$'\t' read -r cluster_id cluster_size centroid_smiles; do
    if [[ "${cluster_id}" == "cluster_id" ]]; then
      continue
    fi
    cluster_dir="${RESULTS_ROOT}/cluster_${cluster_id}"
    cluster_smi_file="${cluster_dir}/target_cluster_${cluster_id}.smi"
    mkdir -p "${cluster_dir}"
    printf "%s\n" "${centroid_smiles}" > "${cluster_smi_file}"

    echo "== Cluster ${cluster_id} (size=${cluster_size}) ==" | tee -a "${LOG_FILE}"

    for algo in "${SEARCH_ALGOS[@]}"; do
      for model in "${MODELS[@]}"; do
        model_results_dir="${cluster_dir}/${algo}/${model}/"
        mkdir -p "${model_results_dir}"

        if [[ "${RESUME}" == "1" ]]; then
          if [[ -f "${model_results_dir}/.completed" ]] || find "${model_results_dir}" -type f -name "stats.json" | grep -q .; then
            echo "SKIP cluster=${cluster_id} model=${model} algo=${algo} (already completed)" | tee -a "${LOG_FILE}"
            ((success_count++))
            echo "" | tee -a "${LOG_FILE}"
            continue
          fi
        fi

        tmp_override="$(mktemp)"
        cat > "${tmp_override}" <<EOF
model_class: "${model}"
search_algorithm: "${algo}"
results_dir: "${model_results_dir}"
append_timestamp_to_dir: true
search_targets_file: "${cluster_smi_file}"
EOF

        echo ">>> Running cluster=${cluster_id} model=${model} algo=${algo}" | tee -a "${LOG_FILE}"
        echo "${RUNNER[*]} search --config ${BASE_CONFIG} --config ${tmp_override}" >> "${LOG_FILE}"

        if [[ "${DRY_RUN}" == "1" ]]; then
          echo "DRY_RUN=1 -> skipped" | tee -a "${LOG_FILE}"
          ((success_count++))
        else
          if "${RUNNER[@]}" search --config "${BASE_CONFIG}" --config "${tmp_override}" >> "${LOG_FILE}" 2>&1; then
            echo "OK  cluster=${cluster_id} model=${model} algo=${algo}" | tee -a "${LOG_FILE}"
            touch "${model_results_dir}/.completed"
            ((success_count++))
          else
            echo "FAIL cluster=${cluster_id} model=${model} algo=${algo}" | tee -a "${LOG_FILE}"
            ((fail_count++))
          fi
        fi

        rm -f "${tmp_override}"
        echo "" | tee -a "${LOG_FILE}"
      done
    done
  done < "${CENTROID_MAP_FILE}"
else
  for algo in "${SEARCH_ALGOS[@]}"; do
    for model in "${MODELS[@]}"; do
      model_results_dir="${RESULTS_ROOT}/${algo}/${model}/"
      mkdir -p "${model_results_dir}"

      if [[ "${RESUME}" == "1" ]]; then
        if [[ -f "${model_results_dir}/.completed" ]] || find "${model_results_dir}" -type f -name "stats.json" | grep -q .; then
          echo "SKIP model=${model} algo=${algo} (already completed)" | tee -a "${LOG_FILE}"
          ((success_count++))
          echo "" | tee -a "${LOG_FILE}"
          continue
        fi
      fi

      tmp_override="$(mktemp)"
      cat > "${tmp_override}" <<EOF
model_class: "${model}"
search_algorithm: "${algo}"
results_dir: "${model_results_dir}"
append_timestamp_to_dir: true
search_targets_file: "${TARGETS_SMI_FILE}"
EOF

      echo ">>> Running model=${model} algo=${algo}" | tee -a "${LOG_FILE}"
      echo "${RUNNER[*]} search --config ${BASE_CONFIG} --config ${tmp_override}" >> "${LOG_FILE}"

      if [[ "${DRY_RUN}" == "1" ]]; then
        echo "DRY_RUN=1 -> skipped" | tee -a "${LOG_FILE}"
        ((success_count++))
      else
        if "${RUNNER[@]}" search --config "${BASE_CONFIG}" --config "${tmp_override}" >> "${LOG_FILE}" 2>&1; then
          echo "OK  model=${model} algo=${algo}" | tee -a "${LOG_FILE}"
          touch "${model_results_dir}/.completed"
          ((success_count++))
        else
          echo "FAIL model=${model} algo=${algo}" | tee -a "${LOG_FILE}"
          ((fail_count++))
        fi
      fi

      rm -f "${tmp_override}"
      echo "" | tee -a "${LOG_FILE}"
    done
  done
fi

echo "Completed. Success=${success_count} Fail=${fail_count}" | tee -a "${LOG_FILE}"
echo "Full log: ${LOG_FILE}" | tee -a "${LOG_FILE}"
