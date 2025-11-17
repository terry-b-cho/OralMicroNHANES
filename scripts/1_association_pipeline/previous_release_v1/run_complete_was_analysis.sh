#!/bin/bash

# =============================================================================
# WAS Individual Regression Analysis Pipeline
# run_complete_was_analysis.sh
# =============================================================================
# This script runs the individual regression analysis pipeline:
# 1. Submits individual regression jobs for each dependent variable
# 2. Waits for all jobs to complete
# 3. Provides guidance for next steps (aggregation with run_aggregation_only.sh)
#
# Note: Aggregation step is now separate - run run_aggregation_only.sh after this completes
#
# Usage: ./run_complete_was_analysis.sh <analysis_type> <normalization> [test_mode] [memory] [time]
# =============================================================================

# Display help information if requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $0 <analysis_type> <normalization> [test_mode] [memory] [time]"
  echo ""
  echo "This script runs individual regression analyses only:"
  echo "  1. Individual regression analyses for each dependent variable"
  echo "  2. Waits for all regression jobs to complete"
  echo "  3. Provides next step guidance for aggregation"
  echo ""
  echo "IMPORTANT: This script does NOT perform aggregation or multiple comparisons correction."
  echo "After this completes, run: ./scripts/1_association_pipeline/run_aggregation_only.sh"
  echo ""
  echo "Parameters:"
  echo "  analysis_type: 1_demoWAS, 2_oradWAS, 3_exWAS, 4_pheWAS, 5_outWAS, 6_zimWAS"
  echo "  normalization: clr, lognorm, none"
  echo "  test_mode (optional): If set to 'test', run only a subset for testing"
  echo "  memory (optional): Memory allocation (e.g., 64G, 120G). Default: 64G"
  echo "  time (optional): Time limit (e.g., 4:00:00, 12:00:00). Default: 4:00:00"
  echo ""
  echo "Output files (individual results saved in results/<analysis_type>_out/result_<normalization>/):"
  echo "  <dep_var>.rds files (many individual regression results)"
  echo ""
  echo "Next step: Run aggregation script to create final results with FDR correction"
  exit 0
fi

# Check required parameters
if [[ $# -lt 2 ]]; then
  echo "Error: Missing required parameters"
  echo "Usage: $0 <analysis_type> <normalization> [test_mode] [memory] [time]"
  echo "Run with -h or --help for more information"
  exit 1
fi

# Parse command line arguments
ANALYSIS_TYPE=$1
NORMALIZATION=$2
TEST_MODE=${3:-""}
MEMORY=${4:-"64G"}
TIME_LIMIT=${5:-"4:00:00"}

# Set paths
BASE_DIR="/n/groups/patel/terry/nhanes_oral_mirco_cho"
RESULTS_DIR="${BASE_DIR}/results"
SCRIPT_DIR="${BASE_DIR}/scripts/1_association_pipeline"
INDIVIDUAL_SCRIPT="${SCRIPT_DIR}/run_was_analysis_debug.sh"
AGGREGATION_SCRIPT="${SCRIPT_DIR}/aggregate_was_results.R"

echo "============================================================================="
echo "WAS Individual Regression Analysis Pipeline"
echo "============================================================================="
echo "Analysis Type: $ANALYSIS_TYPE"
echo "Normalization: $NORMALIZATION"
echo "Test Mode: ${TEST_MODE:-'No'}"
echo "Memory: $MEMORY"
echo "Time Limit: $TIME_LIMIT"
echo "Timestamp: $(date)"
echo "============================================================================="

# Step 1: Run individual analyses
echo ""
echo "STEP 1: Running individual regression analyses..."
echo "============================================================================="

if [[ "$TEST_MODE" == "test" ]]; then
  echo "Running in test mode..."
  $INDIVIDUAL_SCRIPT $ANALYSIS_TYPE $NORMALIZATION test $MEMORY $TIME_LIMIT
else
  echo "Running full analysis..."
  $INDIVIDUAL_SCRIPT $ANALYSIS_TYPE $NORMALIZATION "" $MEMORY $TIME_LIMIT
fi

if [[ $? -ne 0 ]]; then
  echo "ERROR: Individual analysis submission failed"
  exit 1
fi

echo "Individual analyses submitted successfully"

# Step 2: Wait for all jobs to complete
echo ""
echo "STEP 2: Waiting for all jobs to complete..."
echo "============================================================================="

# Function to check if any jobs are still running
check_running_jobs() {
  local job_pattern="${ANALYSIS_TYPE}_${NORMALIZATION}"
  local running_jobs=$(squeue -u $(whoami) --name="$job_pattern*" --noheader | wc -l)
  echo $running_jobs
}

# Wait for jobs to complete
echo "Monitoring job completion..."
WAIT_COUNT=0
MAX_WAIT=720  # Maximum wait time in minutes (12 hours)

while true; do
  RUNNING_JOBS=$(check_running_jobs)
  
  if [[ $RUNNING_JOBS -eq 0 ]]; then
    echo "All jobs completed!"
    break
  fi
  
  echo "$(date): $RUNNING_JOBS jobs still running..."
  
  # Check if we've exceeded maximum wait time
  if [[ $WAIT_COUNT -ge $MAX_WAIT ]]; then
    echo "ERROR: Maximum wait time exceeded. Some jobs may still be running."
    echo "Please check job status manually with: squeue -u $(whoami)"
    exit 1
  fi
  
  sleep 60  # Wait 1 minute before checking again
  ((WAIT_COUNT++))
done

echo "All individual analyses completed at $(date)"

# Step 3: Individual regressions completed - Next steps information
echo ""
echo "STEP 3: Individual regression analyses completed successfully!"
echo "============================================================================="

OUTPUT_DIR="${RESULTS_DIR}/${ANALYSIS_TYPE}_out/result_${NORMALIZATION}"
echo "Individual result files saved to: $OUTPUT_DIR"

# Count individual result files
if [[ -d "$OUTPUT_DIR" ]]; then
  INDIVIDUAL_FILES=$(find "$OUTPUT_DIR" -name "*.rds" ! -name "*_complete.rds" ! -name "*_models.rds" | wc -l)
  echo "Number of individual result files created: $INDIVIDUAL_FILES"
else
  echo "Result directory created: $OUTPUT_DIR"
fi

echo ""
echo "🔄 IMPORTANT: Aggregation step has been separated for better pipeline control"
echo "============================================================================="
echo ""
echo "NEXT REQUIRED STEP - Result Aggregation:"
echo "  After ALL individual analyses are complete, run:"
echo "  ./scripts/1_association_pipeline/run_aggregation_only.sh"
echo ""
echo "  This step will:"
echo "  ✓ Combine all individual regression results"
echo "  ✓ Apply FDR and Bonferroni multiple comparisons correction"
echo "  ✓ Create final analysis files (*_tidied_complete.rds)"
echo ""
echo "OPTIONAL - Pipeline Diagnostics:"
echo "  To check completion status and resume missing analyses:"
echo "  sbatch scripts/1_association_pipeline/submit_smart_resume.sh"
echo ""
echo "  Smart resume will:"
echo "  ✓ Identify missing individual regressions"
echo "  ✓ Submit targeted jobs for missing variables only"
echo "  ✓ Skip completed analyses (does NOT perform aggregation)"
echo ""
echo "💡 MONITORING:"
echo "  Check job status: squeue -u \$USER"
echo "  View logs: ls $OUTPUT_DIR/logs/"

# Step 4: Final summary
echo ""
echo "STEP 4: Individual Regression Analysis Summary"
echo "============================================================================="
echo "Analysis Type: $ANALYSIS_TYPE"
echo "Normalization: $NORMALIZATION"
echo "Output Directory: $OUTPUT_DIR"
echo "Completion Time: $(date)"
echo ""
echo "Status: SUCCESS - Individual regressions completed"
echo ""
echo "REMEMBER: This pipeline completed ONLY the individual regression step."
echo ""
echo "TO COMPLETE YOUR ANALYSIS:"
echo "  1. Wait for ALL planned analyses to finish running"
echo "  2. Run aggregation: ./scripts/1_association_pipeline/run_aggregation_only.sh"
echo "  3. Use p.value.fdr < 0.05 for reporting significant associations"
echo ""
echo "If some analyses fail, use submit_smart_resume.sh to identify and rerun missing pieces"
echo "=============================================================================" 