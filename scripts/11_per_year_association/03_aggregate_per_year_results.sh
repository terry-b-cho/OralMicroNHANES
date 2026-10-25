#!/bin/bash

# =============================================================================
# MODULE 11 — STEP 03: aggregate per-year results (scheme-wise FDR + Storey q)
# 03_aggregate_per_year_results.sh
# =============================================================================
# Runs the module-11 copy of aggregate_was_results.R in INDIVIDUAL mode for each
# (cycle x analysis) on the CLR per-year outputs from step 02. TRUE = use_qvalue
# (Storey q-values), exactly as the production aggregation.
#
# FDR family = one (analysis x clr x cycle) directory — the per-year analogue of
# the production scheme-wise correction. Writes ONLY inside the per_year_F/ and
# per_year_G/ trees (aggregated *_complete.rds + *_aggregation_summary.txt land
# next to the per-dep-var .rds files, mirroring production layout).
#
# Run AFTER all step-02 SLURM jobs have finished:
#   squeue -u $USER | grep -E 'm11F_|m11G_'   # must be empty
# =============================================================================

set -u

# PROJECT_ROOT: path of your local clone; export it or run from the repo root.
BASE_DIR="${PROJECT_ROOT:-$(pwd)}"
SCRIPT_DIR="${BASE_DIR}/scripts/11_per_year_association"
OUT_ROOT="${BASE_DIR}/results/analyses_results/11_per_year_association"
LOG_DIR="${OUT_ROOT}/logs/per_year"
# Conda env: create once from envs/nhanes-analysis_for_reviewers.yml
#   (conda env create -f envs/nhanes-analysis_for_reviewers.yml -n nhanes-analysis),
# then export CONDA_ENV or keep the default name below.
CONDA_ENV="${CONDA_ENV:-nhanes-analysis}"

CYCLES=("F" "G")
ANALYSIS_TYPES=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS")
NORM="clr"

mkdir -p "$LOG_DIR"

# Refuse to aggregate while per-year jobs are still running/pending
if squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -qE "^m11F_|^m11G_"; then
    echo "ERROR: m11F_/m11G_ jobs are still in the queue — wait for them to finish."
    exit 1
fi

# Load environment (same as job environment)
module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

echo "MODULE 11 STEP 03: aggregation + scheme-wise BH FDR + Storey q-values"
fail=0
for cycle in "${CYCLES[@]}"; do
    for analysis in "${ANALYSIS_TYPES[@]}"; do
        results_dir="${OUT_ROOT}/per_year_${cycle}"
        in_dir="${results_dir}/${analysis}_out/result_${NORM}"
        log="${LOG_DIR}/aggregate_${cycle}_${analysis}_${NORM}.log"

        if [[ ! -d "$in_dir" ]]; then
            echo "  SKIP (no results dir): $in_dir"
            continue
        fi

        echo "  Aggregating cycle $cycle / $analysis / $NORM ..."
        if Rscript "${SCRIPT_DIR}/aggregate_was_results.R" "$analysis" "$NORM" "$results_dir" TRUE > "$log" 2>&1; then
            echo "    OK  ($(basename "$log"))"
        else
            echo "    FAILED — see $log"
            fail=1
        fi
    done
done

echo ""
if [[ $fail -eq 0 ]]; then
    echo "STEP 03 complete. Aggregated files are in ${OUT_ROOT}/per_year_{F,G}/<analysis>_out/result_${NORM}/"
    echo "NEXT STEPS: Rscript ${SCRIPT_DIR}/04_variable_availability.R  (then 05..08)"
    exit 0
else
    echo "STEP 03 finished with failures — inspect the logs above, then re-run this script."
    exit 1
fi
