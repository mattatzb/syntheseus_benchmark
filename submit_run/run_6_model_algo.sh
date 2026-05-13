#!/usr/bin/env bash
set -euo pipefail

BASE_CONFIG="${1:-/home/tatzber/syntheseus_benchmark/test_run.yml}"
RUN_ROOT="${2:-/home/tatzber/syntheseus_benchmark/results_6_runs}"
TARGETS_CSV="${TARGETS_CSV:-/home/tatzber/syntheseus_benchmark/pesticides_qed0_5_tanimotoclust_0_6_top100clust.csv}"
SMILES_COLUMN="${SMILES_COLUMN:-centroid_smiles}"
DRY_RUN="${DRY_RUN:-1}"  # 1 = prepare only, 0 = execute syntheseus
CONDA_ENV_NAME="${CONDA_ENV_NAME:-syntheseus-full}"

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
  if command -v conda >/dev/null 2>&1; then
    RUNNER=(conda run -n "${CONDA_ENV_NAME}" syntheseus)
  else
    echo "Neither 'syntheseus' nor 'conda' found in PATH." >&2
    exit 1
  fi
fi

mkdir -p "${RUN_ROOT}/configs"
TARGETS_SMI="${RUN_ROOT}/targets_top100clust.smi"
LOG_FILE="${RUN_ROOT}/run.log"

python - <<PY
import pandas as pd
from pathlib import Path

csv_path = Path("${TARGETS_CSV}")
out_smi = Path("${TARGETS_SMI}")
col = "${SMILES_COLUMN}"

df = pd.read_csv(csv_path)
if col not in df.columns:
    raise SystemExit(f"Column '{col}' not found in {csv_path}")

smiles = df[col].dropna().astype(str)
smiles = smiles[(smiles != "") & (smiles != "Failed")].drop_duplicates()
if len(smiles) == 0:
    raise SystemExit(f"No usable SMILES found in column '{col}' from {csv_path}")

out_smi.write_text("\\n".join(smiles.tolist()) + "\\n")
print(f"Wrote {len(smiles)} targets to {out_smi}")
PY

echo "Base config : ${BASE_CONFIG}" > "${LOG_FILE}"
echo "Run root    : ${RUN_ROOT}" >> "${LOG_FILE}"
echo "Targets CSV : ${TARGETS_CSV}" >> "${LOG_FILE}"
echo "Targets SMI : ${TARGETS_SMI}" >> "${LOG_FILE}"
echo "SMILES col  : ${SMILES_COLUMN}" >> "${LOG_FILE}"
echo "Runner      : ${RUNNER[*]}" >> "${LOG_FILE}"
echo "Dry run     : ${DRY_RUN}" >> "${LOG_FILE}"
echo "" >> "${LOG_FILE}"

declare -a MODELS=("LocalRetro" "MEGAN" "RetroKNN")
declare -a ALGOS=("mcts" "retro_star")

for model in "${MODELS[@]}"; do
  for algo in "${ALGOS[@]}"; do
    run_name="${model}_${algo}"
    override_file="${RUN_ROOT}/configs/${run_name}.yml"
    run_results_dir="${RUN_ROOT}/${run_name}/"

    cat > "${override_file}" <<EOF
model_class: "${model}"
search_algorithm: "${algo}"
search_targets_file: "${TARGETS_SMI}"
results_dir: "${run_results_dir}"
append_timestamp_to_dir: true
EOF

    echo "Prepared ${override_file}" | tee -a "${LOG_FILE}"
    echo "${RUNNER[*]} search --config ${BASE_CONFIG} --config ${override_file}" >> "${LOG_FILE}"

    if [[ "${DRY_RUN}" == "1" ]]; then
      echo "DRY_RUN=1 -> skipped ${run_name}" | tee -a "${LOG_FILE}"
    else
      if "${RUNNER[@]}" search --config "${BASE_CONFIG}" --config "${override_file}" >> "${LOG_FILE}" 2>&1; then
        echo "OK ${run_name}" | tee -a "${LOG_FILE}"
      else
        echo "FAIL ${run_name}" | tee -a "${LOG_FILE}"
      fi
    fi
    echo "" | tee -a "${LOG_FILE}"
  done
done

echo "Done. Log: ${LOG_FILE}"
