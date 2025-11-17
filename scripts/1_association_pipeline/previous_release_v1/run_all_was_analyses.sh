#!/bin/bash

# =============================================================================
# Complete WAS Analysis Suite - All Individual Regression Analyses
# run_all_was_analyses.sh
# =============================================================================
# This script submits all 18 individual regression analyses:
# - 6 analysis types: 1_demoWAS, 2_oradWAS, 3_exWAS, 4_pheWAS, 5_outWAS, 6_zimWAS
# - 3 normalizations: clr, lognorm, none
# - Total: 18 individual regression analyses (NO aggregation step)
#
# Note: After all analyses complete, run run_aggregation_only.sh for final results
#
# Usage: ./run_all_was_analyses.sh [test_mode] [memory] [time]
# =============================================================================

# Parse command line arguments
TEST_MODE=${1:-""}
MEMORY=${2:-"64G"}
TIME_LIMIT=${3:-"4:00:00"}

# Display help information if requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $0 [test_mode] [memory] [time]"
  echo ""
  echo "This script runs ALL individual regression analyses (18 total):"
  echo "  - 6 analysis types × 3 normalization methods"
  echo "  - Each runs individual regressions only (NO aggregation or FDR correction)"
  echo ""
  echo "IMPORTANT: After all analyses complete, run run_aggregation_only.sh for final results"
  echo ""
  echo "Parameters:"
  echo "  test_mode (optional): If set to 'test', run in test mode for all analyses"
  echo "  memory (optional): Memory allocation (e.g., 64G, 120G). Default: 64G"
  echo "  time (optional): Time limit (e.g., 4:00:00, 12:00:00). Default: 4:00:00"
  echo ""
  echo "Examples:"
  echo "  ./run_all_was_analyses.sh                    # Full analyses"
  echo "  ./run_all_was_analyses.sh test               # Test mode"
  echo "  ./run_all_was_analyses.sh \"\" 120G 12:00:00   # Custom resources"
  exit 0
fi

# Set paths
BASE_DIR="/n/groups/patel/terry/nhanes_oral_mirco_cho"
SCRIPT_DIR="${BASE_DIR}/scripts/1_association_pipeline"
MAIN_SCRIPT="${SCRIPT_DIR}/run_complete_was_analysis.sh"
LOG_FILE="${BASE_DIR}/results/all_was_analyses_$(date +%Y%m%d_%H%M%S).log"

# Analysis types and normalizations
ANALYSIS_TYPES=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS" "6_zimWAS")
NORMALIZATIONS=("clr" "lognorm" "none")

echo "============================================================================="
echo "Complete WAS Individual Regression Suite Submission"
echo "============================================================================="
echo "Submitting ALL 18 individual regression analyses:"
echo "  - Analysis types: ${ANALYSIS_TYPES[*]}"
echo "  - Normalizations: ${NORMALIZATIONS[*]}"
echo "  - Test mode: ${TEST_MODE:-'No'}"
echo "  - Memory: $MEMORY"
echo "  - Time limit: $TIME_LIMIT"
echo "  - Log file: $LOG_FILE"
echo "  - Timestamp: $(date)"
echo "============================================================================="

# Create log file
mkdir -p "$(dirname "$LOG_FILE")"
echo "WAS Analysis Suite Submission Log - $(date)" > "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"

# Function to submit analysis and log result
submit_analysis() {
    local analysis_type=$1
    local normalization=$2
    local count=$3
    local total=$4
    
    echo ""
    echo "[$count/$total] Submitting: $analysis_type with $normalization normalization"
    echo "Command: $MAIN_SCRIPT $analysis_type $normalization $TEST_MODE $MEMORY $TIME_LIMIT"
    
    # Log to file
    echo "[$count/$total] $(date): Submitting $analysis_type $normalization" >> "$LOG_FILE"
    
    # Submit the analysis
    if [[ "$TEST_MODE" == "test" ]]; then
        $MAIN_SCRIPT $analysis_type $normalization test $MEMORY $TIME_LIMIT
    else
        $MAIN_SCRIPT $analysis_type $normalization "" $MEMORY $TIME_LIMIT
    fi
    
    local exit_code=$?
    
    if [[ $exit_code -eq 0 ]]; then
        echo "✓ Successfully submitted: $analysis_type $normalization"
        echo "✓ $(date): SUCCESS - $analysis_type $normalization" >> "$LOG_FILE"
    else
        echo "✗ Failed to submit: $analysis_type $normalization (exit code: $exit_code)"
        echo "✗ $(date): FAILED - $analysis_type $normalization (exit code: $exit_code)" >> "$LOG_FILE"
        return $exit_code
    fi
    
    # Small delay to avoid overwhelming the scheduler
    sleep 2
}

# Submit all analyses
echo ""
echo "Starting submission of all analyses..."

count=0
total=$((${#ANALYSIS_TYPES[@]} * ${#NORMALIZATIONS[@]}))
failed_submissions=0

for analysis_type in "${ANALYSIS_TYPES[@]}"; do
    for normalization in "${NORMALIZATIONS[@]}"; do
        ((count++))
        
        if ! submit_analysis "$analysis_type" "$normalization" "$count" "$total"; then
            ((failed_submissions++))
        fi
    done
done

echo ""
echo "============================================================================="
echo "Submission Summary"
echo "============================================================================="
echo "Total analyses: $total"
echo "Successfully submitted: $((total - failed_submissions))"
echo "Failed submissions: $failed_submissions"
echo "Log file: $LOG_FILE"
echo "Completion time: $(date)"

# Log summary
echo "" >> "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"
echo "Submission Summary - $(date)" >> "$LOG_FILE"
echo "Total analyses: $total" >> "$LOG_FILE"
echo "Successfully submitted: $((total - failed_submissions))" >> "$LOG_FILE"
echo "Failed submissions: $failed_submissions" >> "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"

if [[ $failed_submissions -eq 0 ]]; then
    echo ""
    echo "All $total individual regression analyses submitted successfully!"
    echo ""
    echo "💡 MONITORING:"
    echo "  Check job status: squeue -u \$USER"
    echo "  View individual logs: ls results/*/logs/"
    echo ""
    echo "NEXT REQUIRED STEP - Result Aggregation:"
    echo "  After ALL analyses complete, run aggregation:"
    echo "  ./scripts/1_association_pipeline/run_aggregation_only.sh"
    echo ""
    echo "OPTIONAL - Pipeline Diagnostics:"
    echo "  To check completion status and resume missing analyses:"
    echo "  sbatch scripts/1_association_pipeline/submit_smart_resume.sh"
    echo ""
    echo "📁 INDIVIDUAL RESULTS:"
    echo "  Individual regression files will be saved to:"
    for analysis_type in "${ANALYSIS_TYPES[@]}"; do
        for normalization in "${NORMALIZATIONS[@]}"; do
            echo "    results/${analysis_type}_out/result_${normalization}/"
        done
    done
    echo ""
    echo "FINAL OUTPUT (after aggregation):"
    echo "  - *_tidied_complete.rds files with FDR and Bonferroni correction"
    echo "  - Use p.value.fdr < 0.05 for significance testing"
    
    exit 0
else
    echo ""
    echo "⚠️  Some submissions failed. Check the log file for details:"
    echo "  cat $LOG_FILE"
    echo ""
    echo "You may need to resubmit failed analyses manually."
    
    exit 1
fi 