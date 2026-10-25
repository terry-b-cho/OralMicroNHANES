#!/bin/bash

# =============================================================================
# MODULE 11 — STEP 02: per-year (single-cycle) WAS submitter, CLR only
# 02_submit_per_year_was_clr.sh
# =============================================================================
# Adapted copy of run_all_was_analyses_flexible_11.sh (itself a path-adapted copy
# of the production submitter). Differences for the per-year task:
#   * Loops over CYCLES=(F G): each combo uses the single-cycle schema from
#     step 01 (schema_per_year/<cycle>/), so universal_was_analysis.R fits
#     single_F / single_G survey designs (unscaled WTMEC2YR, SDMVPSU/SDMVSTRA,
#     nest=TRUE) — identical regressions, just not pooled.
#   * CLR normalization only.
#   * RIGHT-SIZED resources (O2 memory-efficiency notice): 1G per job instead of
#     the production 32G (measured MaxRSS of these jobs is ~0.3-0.45G), with
#     per-analysis time caps derived from measured production wall-times.
#   * Outputs: <OUT_ROOT>/per_year_<cycle>/<analysis>_out/result_clr/
#     Logs:    <OUT_ROOT>/logs/per_year/
#     Test:    <OUT_ROOT>/test_runs_per_year/  (isolated from real outputs)
#   * Job names m11F_* / m11G_*.
# Inherited unchanged: flexible resume (never resubmits an existing .rds),
# triple overwrite protection, per-dependent-variable job granularity, the
# byte-identical engine universal_was_analysis.R, and the conda environment.
#
# Usage: ./02_submit_per_year_was_clr.sh [test_mode]
#   test_mode: 'test' submits 1 dep var per (cycle x analysis), 5 pairs each.
# Resume: just re-run; only missing .rds files are submitted.
# =============================================================================

TEST_MODE=${1:-""}

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $0 [test_mode]"
  echo "  test_mode: 'test' for a smoke test (outputs isolated in test_runs_per_year/)"
  echo "Full run: 2,851 jobs per cycle = 5,702 jobs (1G each, short partition)."
  exit 0
fi

# =============================================================================
# PATH BLOCK
# =============================================================================
# PROJECT_ROOT: path of your local clone; export it or run from the repo root.
BASE_DIR="${PROJECT_ROOT:-$(pwd)}"
SCRIPT_DIR="${BASE_DIR}/scripts/11_per_year_association"
ANALYSIS_R_SCRIPT="${SCRIPT_DIR}/universal_was_analysis.R"

# Read-only input database (same DB as production)
DB_PATH="${BASE_DIR}/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"

# All writes stay under OUT_ROOT
OUT_ROOT="${BASE_DIR}/results/analyses_results/11_per_year_association"
SCHEMA_ROOT="${OUT_ROOT}/schema_per_year"          # produced by step 01
LOG_DIR="${OUT_ROOT}/logs/per_year"
LOG_FILE="${LOG_DIR}/per_year_was_submit_$(date +%Y%m%d_%H%M%S).log"

# Conda env: create once from envs/nhanes-analysis_for_reviewers.yml
#   (conda env create -f envs/nhanes-analysis_for_reviewers.yml -n nhanes-analysis),
# then export CONDA_ENV or keep the default name below.
CONDA_ENV="${CONDA_ENV:-nhanes-analysis}"
# =============================================================================

CYCLES=("F" "G")
ANALYSIS_TYPES=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS")
NORMALIZATION="clr"

# RIGHT-SIZED resources. Memory: measured MaxRSS of these per-year jobs is
# 380-460MB (sacct, 2026-08-30 smoke test). 768M => ~55-60% utilization,
# satisfying the O2 memory-efficiency notice (production wastefully used 32G,
# ~1.3% utilization). If any job is OOM-killed, raise MEM here and re-run this
# script — the flexible resume only resubmits the missing .rds files.
# Time: measured pooled clr wall-times were 1_demoWAS ~2min, 3_exWAS ~1h20,
# 2_oradWAS ~2h12, 5_outWAS ~3h56, 4_pheWAS ~5h31; per-year fits use half the
# data, caps below leave ~2x headroom over the per-year expectation.
get_resources() {
    local analysis_type=$1
    local test_mode=$2

    if [[ "$test_mode" == "test" ]]; then
        echo "768M 0-00:15"
        return
    fi

    # Time caps raised for 3_exWAS/2_oradWAS after the 2026-08-30 full run: 23 jobs
    # on slow/contended nodes hit the original 2h/3h caps (TIME LIMIT, not memory).
    case $analysis_type in
        1_demoWAS) echo "768M 0-00:30" ;;
        3_exWAS)   echo "768M 0-04:00" ;;
        2_oradWAS) echo "768M 0-06:00" ;;
        5_outWAS)  echo "768M 0-05:00" ;;
        4_pheWAS)  echo "768M 0-06:00" ;;
        *)         echo "768M 0-04:00" ;;
    esac
}

check_result_file_exists() {
    local file_path=$1
    if [[ -f "$file_path" ]]; then
        local file_size=$(stat -f%z "$file_path" 2>/dev/null || stat -c%s "$file_path" 2>/dev/null || echo 0)
        if [[ $file_size -gt 100 ]]; then
            return 0
        fi
    fi
    return 1
}

get_missing_dependent_variables() {
    local cycle=$1
    local analysis_type=$2
    local test_mode=$3

    local schema_file="${SCHEMA_ROOT}/${cycle}/${analysis_type}_${NORMALIZATION}_schema_structure.csv"
    local out_base="${OUT_ROOT}/per_year_${cycle}"
    [[ "$test_mode" == "test" ]] && out_base="${OUT_ROOT}/test_runs_per_year/per_year_${cycle}"
    local output_dir="${out_base}/${analysis_type}_out/result_${NORMALIZATION}"

    if [[ ! -f "$schema_file" ]]; then
        echo "ERROR: Schema file not found (run step 01 first): $schema_file" >&2
        return 1
    fi

    local expected_vars_file=$(mktemp)
    cut -d',' -f1 "$schema_file" | tail -n +2 | sort | uniq > "$expected_vars_file"
    local total_expected=$(wc -l < "$expected_vars_file")

    if [[ "$test_mode" == "test" ]]; then
        head -n 1 "$expected_vars_file" > "${expected_vars_file}.test"
        mv "${expected_vars_file}.test" "$expected_vars_file"
        total_expected=1
    fi

    local existing_vars_file=$(mktemp)
    if [[ -d "$output_dir" ]]; then
        find "$output_dir" -maxdepth 1 -name "*.rds" -type f -exec basename {} .rds \; | sort > "$existing_vars_file"
    else
        touch "$existing_vars_file"
    fi
    local total_existing=$(wc -l < "$existing_vars_file")

    echo "   File count: $total_existing/$total_expected exist" >&2

    local missing_vars_file=$(mktemp)
    comm -23 "$expected_vars_file" "$existing_vars_file" > "$missing_vars_file"

    echo "$(date): $cycle $analysis_type - $total_existing/$total_expected exist, $(wc -l < "$missing_vars_file") missing" >> "$LOG_FILE"

    cat "$missing_vars_file"
    rm -f "$expected_vars_file" "$existing_vars_file" "$missing_vars_file"
}

submit_missing_variables_only() {
    local cycle=$1
    local analysis_type=$2
    local missing_vars=$3
    local test_mode=$4
    local memory=$5
    local time=$6

    local schema_file="${SCHEMA_ROOT}/${cycle}/${analysis_type}_${NORMALIZATION}_schema_structure.csv"
    local out_base="${OUT_ROOT}/per_year_${cycle}"
    [[ "$test_mode" == "test" ]] && out_base="${OUT_ROOT}/test_runs_per_year/per_year_${cycle}"
    local output_dir="${out_base}/${analysis_type}_out/result_${NORMALIZATION}"

    if [[ ! -f "$DB_PATH" ]]; then
        echo "ERROR: Database file not found: $DB_PATH" >&2
        return 1
    fi

    mkdir -p "$output_dir" "$LOG_DIR"

    if [[ ! -w "$output_dir" ]]; then
        echo "ERROR: Output directory not writable: $output_dir" >&2
        return 1
    fi

    local counter=0
    local total_missing=0
    if [[ -n "$missing_vars" ]]; then
        total_missing=$(echo "$missing_vars" | wc -l)
    fi

    while IFS= read -r dep_var; do
        [[ -z "$dep_var" ]] && continue
        counter=$((counter + 1))

        local result_file="${output_dir}/${dep_var}.rds"
        if check_result_file_exists "$result_file"; then
            echo "   SAFETY: $result_file appeared since last check - SKIPPING"
            continue
        fi

        local clean_var=$(echo "$dep_var" | sed 's/[^a-zA-Z0-9_]/_/g')
        local job_name="m11${cycle}_${analysis_type}_${clean_var}"
        [[ "$test_mode" == "test" ]] && job_name="${job_name}_test"

        echo "   [$counter/$total_missing] Submitting: $dep_var ($cycle)"
        echo "$(date): Submitting job $job_name for $dep_var" >> "$LOG_FILE"

        sbatch << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$LOG_DIR/${job_name}_%j.out
#SBATCH --error=$LOG_DIR/${job_name}_%j.err
#SBATCH --partition=short
#SBATCH --time=$time
#SBATCH --mem=$memory
#SBATCH --cpus-per-task=2

# Load environment (identical to production job environment)
module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "\$(conda shell.bash hook)"
conda activate $CONDA_ENV

cd $BASE_DIR

# CRITICAL SAFETY: final check before running analysis
RESULT_FILE="$output_dir/$dep_var.rds"
if [[ -f "\$RESULT_FILE" ]]; then
    echo "SAFETY: \$RESULT_FILE already exists - ABORTING to prevent overwrite"
    exit 0
fi

echo "Starting MODULE 11 per-year (single-cycle $cycle) analysis for $dep_var at \$(date)"
echo "Analysis: $analysis_type, Normalization: $NORMALIZATION, Cycle: $cycle"
echo "Resources: $memory memory, $time time"
echo "Output file: \$RESULT_FILE"

if [[ "$test_mode" == "test" ]]; then
    Rscript $ANALYSIS_R_SCRIPT \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$schema_file" \\
        --database_path "$DB_PATH" \\
        --output_path "$output_dir" \\
        --analysis_type "$analysis_type" \\
        --normalization "$NORMALIZATION" \\
        --test
else
    Rscript $ANALYSIS_R_SCRIPT \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$schema_file" \\
        --database_path "$DB_PATH" \\
        --output_path "$output_dir" \\
        --analysis_type "$analysis_type" \\
        --normalization "$NORMALIZATION"
fi

if [[ -f "\$RESULT_FILE" ]]; then
    echo "OK: \$RESULT_FILE created successfully"
else
    echo "ERROR: expected output \$RESULT_FILE was not created"
fi
echo "MODULE 11 per-year analysis for $dep_var ($cycle) completed at \$(date)"
EOF

        sleep 0.2
    done <<< "$missing_vars"

    return 0
}

echo "============================================================================="
echo "MODULE 11 STEP 02: per-year (F and G separately) WAS submission — CLR only"
echo "============================================================================="
echo "  - Cycles: ${CYCLES[*]}   Analyses: ${ANALYSIS_TYPES[*]}"
echo "  - Normalization: $NORMALIZATION"
echo "  - Test mode: ${TEST_MODE:-'No'}"
echo "  - Schema root: $SCHEMA_ROOT"
echo "  - Output root: $OUT_ROOT/per_year_{F,G}/"
echo "  - Log file: $LOG_FILE"
echo "============================================================================="

if [[ ! -f "$ANALYSIS_R_SCRIPT" ]]; then
    echo "ERROR: Analysis R script not found: $ANALYSIS_R_SCRIPT"
    exit 1
fi

mkdir -p "$LOG_DIR"
echo "MODULE 11 per-year submission - $(date)" > "$LOG_FILE"
echo "Test Mode: ${TEST_MODE:-'No'}" >> "$LOG_FILE"

total_submitted=0
total_skipped=0
failed_submissions=0

for cycle in "${CYCLES[@]}"; do
    for analysis_type in "${ANALYSIS_TYPES[@]}"; do
        echo ""
        echo "CHECKING: cycle $cycle, $analysis_type ($NORMALIZATION)"
        echo "-----------------------------------------------------------------------------"

        missing_vars=$(get_missing_dependent_variables "$cycle" "$analysis_type" "$TEST_MODE")
        if [[ $? -ne 0 ]]; then
            echo "ERROR: failed to check missing variables for $cycle $analysis_type"
            failed_submissions=$((failed_submissions + 1))
            continue
        fi

        missing_count=0
        [[ -n "$missing_vars" ]] && missing_count=$(echo "$missing_vars" | wc -l)

        if [[ $missing_count -eq 0 ]]; then
            echo "OK: cycle $cycle $analysis_type: ALL JOBS COMPLETED - SKIPPING"
            total_skipped=$((total_skipped + 1))
            continue
        fi

        resources=$(get_resources "$analysis_type" "$TEST_MODE")
        memory=$(echo $resources | cut -d' ' -f1)
        time=$(echo $resources | cut -d' ' -f2)
        echo "Submitting $missing_count jobs ($memory, $time)"

        if submit_missing_variables_only "$cycle" "$analysis_type" "$missing_vars" "$TEST_MODE" "$memory" "$time"; then
            total_submitted=$((total_submitted + missing_count))
        else
            failed_submissions=$((failed_submissions + 1))
        fi

        sleep 1
    done
done

echo ""
echo "============================================================================="
echo "MODULE 11 PER-YEAR SUBMISSION SUMMARY"
echo "============================================================================="
echo "Total jobs submitted: $total_submitted"
echo "Combos skipped (complete): $total_skipped"
echo "Failed submissions: $failed_submissions"
echo "Log file: $LOG_FILE"
echo ""
echo "MONITORING:  squeue -u \$USER | grep -E 'm11F_|m11G_'"
echo "ERRORS:      find $LOG_DIR -name 'm11*_*.err' -size +0"
echo ""
echo "NEXT STEP after ALL jobs finish:"
echo "  ./scripts/11_per_year_association/03_aggregate_per_year_results.sh"

[[ $failed_submissions -eq 0 ]] && exit 0 || exit 1
