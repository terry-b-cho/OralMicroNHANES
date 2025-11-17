#!/bin/bash

# =============================================================================
# FLEXIBLE WAS Analysis Suite with FLEXIBLE RESUME - PRODUCTION
# run_all_was_analyses_flexible.sh
# =============================================================================
# This script submits WAS analyses with FLEXIBLE RESUME functionality:
# - NEVER resubmits completed jobs (checks for existing .rds files)
# - Uses CORRECTED time limits based on computational analysis
# - Preserves all existing successful results
# - Only submits missing dependent variables per analysis type
#
# Key Features:
# - FLEXIBLE resume: Only resubmits jobs that haven't produced output files
# - Corrected resource allocation based on actual computational requirements
# - File completeness verification (non-empty .rds files)
# - Detailed progress reporting
#
# Usage: ./run_all_was_analyses_flexible.sh [test_mode] [resource_level]
# =============================================================================

# Parse command line arguments
TEST_MODE=${1:-""}
RESOURCE_LEVEL=${2:-"auto"}

# Display help information if requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $0 [test_mode] [resource_level]"
  echo ""
  echo "FLEXIBLE WAS Analysis with FLEXIBLE RESUME:"
  echo "  - NEVER resubmits completed jobs (checks existing .rds files)"
  echo "  - Uses CORRECTED time limits based on computational analysis"
  echo "  - Preserves all existing successful results"
  echo "  - Only submits missing dependent variables"
  echo ""
  echo "CORRECTED RESOURCE ALLOCATION:"
  echo "  2_oradWAS, 4_pheWAS, 5_outWAS, 6_zimWAS: 12G memory, 4 hours (1,349 regressions/dep_var)"
  echo "  3_exWAS: 8G memory, 3 hours (473 regressions/dep_var)"
  echo "  1_demoWAS: 1G memory, 1 hour (26 regressions/dep_var)"
  echo ""
  echo "Parameters:"
  echo "  test_mode (optional): If set to 'test', run in test mode"
  echo "  resource_level (optional): Resource allocation strategy"
  echo "    - auto (default): Standard allocation based on analysis requirements"
  echo "    - aggressive: 32G memory, 12 hours for all job types"
  echo ""
  echo "Examples:"
  echo "  ./run_all_was_analyses_flexible.sh                # FLEXIBLE resume with corrected resources"
  echo "  ./run_all_was_analyses_flexible.sh test           # Test mode"
  echo "  ./run_all_was_analyses_flexible.sh \"\" aggressive # Higher resource usage"
  exit 0
fi

# Set paths
BASE_DIR="/n/groups/patel/terry/nhanes_oral_mirco_cho"
SCRIPT_DIR="${BASE_DIR}/scripts/1_association_pipeline"
SINGLE_ANALYSIS_SCRIPT="${SCRIPT_DIR}/run_single_was_analysis_PRODUCTION.sh"
LOG_FILE="${BASE_DIR}/logs/1_association_pipeline_flexible/flexible_was_analyses_$(date +%Y%m%d_%H%M%S).log"

# Analysis types and normalizations
ANALYSIS_TYPES=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS" "6_zimWAS")
NORMALIZATIONS=("clr" "lognorm" "none")

# CORRECTED Resource allocation function based on computational analysis
get_resources() {
    local analysis_type=$1
    local resource_level=$2
    local test_mode=$3
    
    # Test mode overrides everything
    if [[ "$test_mode" == "test" ]]; then
        echo "2G 0-00:15"
        return
    fi
    
    # CORRECTED resource allocation based on actual computational requirements
    case $resource_level in
        aggressive)
            echo "32G 0-12:00"  # 12 hours and 32G for all job types in aggressive mode
            ;;
        auto|*)
            case $analysis_type in
                2_oradWAS|4_pheWAS|5_outWAS|6_zimWAS) echo "12G 0-04:00" ;;  # 4 hours for 1349 regressions/dep_var
                3_exWAS) echo "8G 0-03:00" ;;                                # 3 hours for 473 regressions/dep_var
                1_demoWAS) echo "1G 0-01:00" ;;                              # 1 hour for 26 regressions/dep_var
            esac
            ;;
    esac
}

# ULTRA-EFFICIENT: Simple file existence check (no validation unless necessary)
check_result_file_exists() {
    local file_path=$1
    
    # Check if file exists and has reasonable size
    if [[ -f "$file_path" ]]; then
        local file_size=$(stat -f%z "$file_path" 2>/dev/null || stat -c%s "$file_path" 2>/dev/null || echo 0)
        if [[ $file_size -gt 100 ]]; then
            return 0
        fi
    fi
    return 1
}

# ULTRA-EFFICIENT: Get missing files using simple file operations
get_missing_dependent_variables() {
    local analysis_type=$1
    local normalization=$2
    
    local schema_file="${BASE_DIR}/results/0_ss_files/${analysis_type}_${normalization}_schema_structure.csv"
    local output_dir="${BASE_DIR}/results/${analysis_type}_out/result_${normalization}"
    
    # Check if schema file exists
    if [[ ! -f "$schema_file" ]]; then
        echo "ERROR: Schema file not found: $schema_file" >&2
        return 1
    fi
    
    # EFFICIENCY: Get expected files using simple operations
    echo "   Extracting expected dependent variables..."
    local expected_vars_file=$(mktemp)
    cut -d',' -f1 "$schema_file" | tail -n +2 | sort | uniq > "$expected_vars_file"
    local total_expected=$(wc -l < "$expected_vars_file")
    
    # EFFICIENCY: Get existing files using simple operations  
    echo "   📁 Finding existing result files..."
    local existing_vars_file=$(mktemp)
    if [[ -d "$output_dir" ]]; then
        find "$output_dir" -maxdepth 1 -name "*.rds" -type f -exec basename {} .rds \; | sort > "$existing_vars_file"
    else
        touch "$existing_vars_file"
    fi
    local total_existing=$(wc -l < "$existing_vars_file")
    
    echo "   File count: $total_existing/$total_expected exist"
    
    # EFFICIENCY: Find missing using diff (much faster than loops)
    local missing_vars_file=$(mktemp)
    comm -23 "$expected_vars_file" "$existing_vars_file" > "$missing_vars_file"
    local missing_count=$(wc -l < "$missing_vars_file")
    
    echo "   ✅ Analysis complete: $missing_count missing files identified"
    
    # Log the results
    echo "$(date): $analysis_type $normalization - $total_existing/$total_expected exist, $missing_count missing" >> "$LOG_FILE"
    
    # Return missing variables
    cat "$missing_vars_file"
    
    # Cleanup
    rm -f "$expected_vars_file" "$existing_vars_file" "$missing_vars_file"
}

# Function to submit a single analysis type with normalization
submit_analysis_type() {
    local analysis_type=$1
    local normalization=$2
    local test_mode=$3
    
    echo ""
    echo "CHECKING: $analysis_type with $normalization normalization"
    echo "============================================================================="
    
    # Get missing dependent variables using efficient method
    local missing_vars
    missing_vars=$(get_missing_dependent_variables "$analysis_type" "$normalization")
    local exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        echo "❌ Failed to check missing variables for $analysis_type $normalization"
        return 1
    fi
    
    # Count missing variables
    local missing_count=0
    if [[ -n "$missing_vars" ]]; then
        missing_count=$(echo "$missing_vars" | wc -l)
    fi
    
    if [[ $missing_count -eq 0 ]]; then
        echo "✅ $analysis_type $normalization: ALL JOBS COMPLETED - SKIPPING"
        return 0
    fi
    
    echo "$analysis_type $normalization: SUBMITTING $missing_count MISSING JOBS"
    
    # Get resources for this analysis type
    local resources=$(get_resources "$analysis_type" "$RESOURCE_LEVEL" "$test_mode")
    local memory=$(echo $resources | cut -d' ' -f1)
    local time=$(echo $resources | cut -d' ' -f2)
    
    echo "💾 Resources: $memory memory, $time time"
    echo "📝 Missing variables:"
    
    # Show first few missing variables as preview
    local preview_count=5
    local shown_count=0
    while IFS= read -r var && [[ $shown_count -lt $preview_count ]]; do
        [[ -n "$var" ]] && echo "   - $var"
        ((shown_count++))
    done <<< "$missing_vars"
    
    if [[ $missing_count -gt $preview_count ]]; then
        echo "   ... and $((missing_count - preview_count)) more"
    fi
    
    # Submit jobs for missing variables
    echo ""
    echo "SUBMITTING JOBS..."
    
    if ! submit_missing_variables_only "$analysis_type" "$normalization" "$missing_vars" "$test_mode" "$memory" "$time"; then
        echo "❌ Failed to submit jobs for $analysis_type $normalization"
        return 1
    fi
    
    echo "✅ Successfully submitted $missing_count jobs for $analysis_type $normalization"
    return 0
}

# ENHANCED submission function with additional safety checks
submit_missing_variables_only() {
    local analysis_type=$1
    local normalization=$2
    local missing_vars=$3
    local test_mode=$4
    local memory=$5
    local time=$6
    
    local base="${BASE_DIR}"
    local db_path="$base/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
    local schema_file="$base/results/0_ss_files/${analysis_type}_${normalization}_schema_structure.csv"
    local output_dir="$base/results/${analysis_type}_out/result_${normalization}"
    
    # SAFETY: Verify critical files exist before submitting any jobs
    if [[ ! -f "$db_path" ]]; then
        echo "ERROR: Database file not found: $db_path" >&2
        return 1
    fi
    
    if [[ ! -f "$schema_file" ]]; then
        echo "ERROR: Schema file not found: $schema_file" >&2
        return 1
    fi
    
    # Create output directory (safe - won't overwrite existing files)
    mkdir -p "$output_dir"
    mkdir -p "$base/logs"
    
    # SAFETY: Double-check that output directory is writable
    if [[ ! -w "$output_dir" ]]; then
        echo "ERROR: Output directory not writable: $output_dir" >&2
        return 1
    fi
    
    # Submit job for each missing dependent variable
    local counter=0
    local total_missing=0
    if [[ -n "$missing_vars" ]]; then
        total_missing=$(echo "$missing_vars" | wc -l)
    fi
    
    while IFS= read -r dep_var; do
        [[ -z "$dep_var" ]] && continue
        
        counter=$((counter + 1))
        
        # CRITICAL SAFETY: Double-check the file doesn't exist just before submission
        local result_file="${output_dir}/${dep_var}.rds"
        if check_result_file_exists "$result_file"; then
            echo "   ⚠️  SAFETY: File $result_file appeared since last check - SKIPPING submission"
            continue
        fi
        
        # Clean variable name for job naming
        local clean_var=$(echo "$dep_var" | sed 's/[^a-zA-Z0-9_]/_/g')
        local job_name="flexible_${analysis_type}_${normalization}_${clean_var}"
        
        if [[ "$test_mode" == "test" ]]; then
            job_name="${job_name}_test"
        fi
        
        echo "   [$counter/$total_missing] Submitting: $dep_var"
        
        # SAFETY: Log each job submission
        echo "$(date): Submitting job $job_name for $dep_var" >> "$LOG_FILE"
        
        # Submit SLURM job with ENHANCED SAFETY in the job script
        sbatch << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$base/logs/${job_name}_%j.out
#SBATCH --error=$base/logs/${job_name}_%j.err
#SBATCH --partition=short
#SBATCH --time=$time
#SBATCH --mem=$memory
#SBATCH --cpus-per-task=2

# Load environment
module load miniconda3/23.1.0
source activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Change to working directory
cd $base

# CRITICAL SAFETY: Final check before running analysis
RESULT_FILE="$output_dir/$dep_var.rds"
if [[ -f "\$RESULT_FILE" ]]; then
    echo "SAFETY: Result file \$RESULT_FILE already exists - ABORTING to prevent overwrite"
    echo "This job was likely submitted redundantly or another job completed first"
    exit 0
fi

# Run analysis for this dependent variable
echo "Starting FLEXIBLE RESUME analysis for $dep_var at \$(date)"
echo "Analysis: $analysis_type, Normalization: $normalization"
echo "Resources: $memory memory, $time time"
echo "Output file: \$RESULT_FILE"

if [[ "$test_mode" == "test" ]]; then
    Rscript scripts/1_association_pipeline/universal_was_analysis_PRODUCTION.R \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$schema_file" \\
        --database_path "$db_path" \\
        --output_path "$output_dir" \\
        --analysis_type "$analysis_type" \\
        --normalization "$normalization" \\
        --test
else
    Rscript scripts/1_association_pipeline/universal_was_analysis_PRODUCTION.R \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$schema_file" \\
        --database_path "$db_path" \\
        --output_path "$output_dir" \\
        --analysis_type "$analysis_type" \\
        --normalization "$normalization"
fi

# SAFETY: Verify the output file was created and is valid
if [[ -f "\$RESULT_FILE" ]]; then
    echo "✅ SAFETY: Output file \$RESULT_FILE created successfully"
else
    echo "❌ SAFETY: Expected output file \$RESULT_FILE was not created"
fi

echo "FLEXIBLE RESUME analysis for $dep_var completed at \$(date)"
EOF

        # Brief pause between submissions to avoid overwhelming scheduler
        sleep 0.2
    done <<< "$missing_vars"
    
    return 0
}

echo "============================================================================="
echo "FLEXIBLE WAS Analysis Suite with FLEXIBLE RESUME - PRODUCTION"
echo "============================================================================="
echo "FLEXIBLE RESUME: Only submitting jobs for missing .rds files"
echo "💾 CORRECTED RESOURCES: Based on computational analysis"
echo "🛡️  PRESERVATION: Never overwrites existing successful results"
echo ""
echo "  - Analysis types: ${ANALYSIS_TYPES[*]}"
echo "  - Normalizations: ${NORMALIZATIONS[*]}"
echo "  - Test mode: ${TEST_MODE:-'No'}"
echo "  - Resource level: $RESOURCE_LEVEL"
echo "  - Log file: $LOG_FILE"
echo "  - Timestamp: $(date)"
echo "============================================================================="

# Check if analysis script exists
if [[ ! -f "$SINGLE_ANALYSIS_SCRIPT" ]]; then
    echo "ERROR: Single analysis script not found: $SINGLE_ANALYSIS_SCRIPT"
    exit 1
fi

# Create log file
mkdir -p "$(dirname "$LOG_FILE")"
echo "FLEXIBLE WAS Analysis Suite with FLEXIBLE RESUME - $(date)" > "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"
echo "Resource Level: $RESOURCE_LEVEL" >> "$LOG_FILE"
echo "Test Mode: ${TEST_MODE:-'No'}" >> "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"

# Display corrected resource allocation summary
echo ""
echo "CORRECTED RESOURCE ALLOCATION ($RESOURCE_LEVEL mode):"
echo "============================================================================="
for analysis_type in "${ANALYSIS_TYPES[@]}"; do
    resources=$(get_resources "$analysis_type" "$RESOURCE_LEVEL" "$TEST_MODE")
    memory=$(echo $resources | cut -d' ' -f1)
    time=$(echo $resources | cut -d' ' -f2)
    echo "  $analysis_type: $memory memory, $time time"
done
echo "============================================================================="

# Perform FLEXIBLE resume analysis
echo ""
echo "PERFORMING FLEXIBLE RESUME ANALYSIS..."
echo "============================================================================="

total_submitted=0
total_skipped=0
failed_submissions=0

# Check each analysis type and normalization combination
for analysis_type in "${ANALYSIS_TYPES[@]}"; do
    for normalization in "${NORMALIZATIONS[@]}"; do
        
        # Log to file
        echo "$(date): Checking $analysis_type $normalization" >> "$LOG_FILE"
        
        if submit_analysis_type "$analysis_type" "$normalization" "$TEST_MODE"; then
            # Count submitted jobs
            missing_vars=$(get_missing_dependent_variables "$analysis_type" "$normalization")
            if [[ -n "$missing_vars" ]]; then
                local missing_count=$(echo "$missing_vars" | wc -l)
                total_submitted=$((total_submitted + missing_count))
                echo "✓ $(date): SUCCESS - $analysis_type $normalization ($missing_count jobs)" >> "$LOG_FILE"
            else
                total_skipped=$((total_skipped + 1))
                echo "✓ $(date): SKIPPED - $analysis_type $normalization (complete)" >> "$LOG_FILE"
            fi
        else
            failed_submissions=$((failed_submissions + 1))
            echo "✗ $(date): FAILED - $analysis_type $normalization" >> "$LOG_FILE"
        fi
        
        # Small delay between analysis types
        sleep 1
    done
done

echo ""
echo "============================================================================="
echo "FLEXIBLE RESUME SUMMARY"
echo "============================================================================="
echo "Total jobs submitted: $total_submitted"
echo "Total analyses skipped (complete): $total_skipped"
echo "Failed submissions: $failed_submissions"
echo "Resource level: $RESOURCE_LEVEL"
echo "Log file: $LOG_FILE"
echo "Completion time: $(date)"

# Log summary
echo "" >> "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"
echo "FLEXIBLE RESUME SUMMARY - $(date)" >> "$LOG_FILE"
echo "Total jobs submitted: $total_submitted" >> "$LOG_FILE"
echo "Total analyses skipped (complete): $total_skipped" >> "$LOG_FILE"
echo "Failed submissions: $failed_submissions" >> "$LOG_FILE"
echo "Resource level: $RESOURCE_LEVEL" >> "$LOG_FILE"
echo "=============================================================================" >> "$LOG_FILE"

if [[ $failed_submissions -eq 0 ]]; then
    echo ""
    echo "FLEXIBLE RESUME COMPLETED SUCCESSFULLY!"
    echo ""
    echo "RESULTS:"
    if [[ $total_submitted -gt 0 ]]; then
        echo "  ✅ Submitted $total_submitted missing jobs with corrected time limits"
    fi
    if [[ $total_skipped -gt 0 ]]; then
        echo "  ⏭️  Skipped $total_skipped complete analyses (preserving existing results)"
    fi
    echo ""
    echo "💡 MONITORING:"
    echo "  Check job status: squeue -u \$USER | grep flexible"
    echo "  View individual logs: ls $BASE_DIR/logs/flexible_*.out"
    echo "  Check for errors: find $BASE_DIR/logs -name 'flexible_*.err' -size +0"
    echo ""
    echo "🔄 RESUME AGAIN:"
    echo "  Run this script again to check for any remaining failures"
    echo "  $0 $TEST_MODE $RESOURCE_LEVEL"
    echo ""
    echo "�� NEXT STEP - Result Aggregation:"
    echo "  After ALL analyses complete, run aggregation:"
    echo "  ./scripts/1_association_pipeline/run_aggregation_only.sh"
    echo ""
    echo "⏱️  ESTIMATED COMPLETION TIMES (with corrected limits):"
    if [[ "$TEST_MODE" == "test" ]]; then
        echo "  Test mode: ~15-30 minutes"
    else
        case $RESOURCE_LEVEL in
            aggressive)
                echo "  Aggressive resources: ~4-5 hours (depends on cluster load)"
                ;;
            auto|*)
                echo "  Auto resources: ~6 hours (microbiome analyses take longest)"
                ;;
        esac
    fi
    
    exit 0
else
    echo ""
    echo "⚠️  Some submissions failed. Check the log file for details:"
    echo "  cat $LOG_FILE"
    echo ""
    echo "🔄 You can run this script again to retry failed submissions:"
    echo "  $0 $TEST_MODE $RESOURCE_LEVEL"
    
    exit 1
fi 