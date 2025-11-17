#!/bin/bash

# =============================================================================
# Run Aggregation Only for Failed Aggregations
# run_aggregation_only.sh
# =============================================================================
# This script identifies analyses that have individual results but failed 
# aggregation, and runs ONLY the aggregation step with FDR/Bonferroni correction.
#
# Usage: ./run_aggregation_only.sh [specific_analysis] [normalization]
# Examples:
#   ./run_aggregation_only.sh                    # Check and aggregate all that need it
#   ./run_aggregation_only.sh 3_exWAS clr        # Aggregate specific analysis
# =============================================================================

# Display help information if requested
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $0 [specific_analysis] [normalization]"
  echo ""
  echo "This script runs ONLY the aggregation step for analyses that:"
  echo "  - Have individual result files"
  echo "  - Are missing aggregated files (*_tidied_complete.rds)"
  echo ""
  echo "Parameters:"
  echo "  specific_analysis (optional): Run only for this analysis (e.g., 3_exWAS)"
  echo "  normalization (optional): Run only for this normalization (e.g., clr)"
  echo ""
  echo "Examples:"
  echo "  $0                           # Check all 18 analyses and aggregate those that need it"
  echo "  $0 3_exWAS                   # Aggregate all normalizations of 3_exWAS"
  echo "  $0 3_exWAS clr              # Aggregate only 3_exWAS clr"
  echo "  $0 \"\" clr                   # Aggregate all analyses with clr normalization"
  echo ""
  echo "Output: Creates *_tidied_complete.rds files with FDR/Bonferroni correction"
  exit 0
fi

# Parse command line arguments
SPECIFIC_ANALYSIS=${1:-""}
SPECIFIC_NORMALIZATION=${2:-""}

# Set paths
BASE_DIR="/n/groups/patel/terry/nhanes_oral_mirco_cho"
RESULTS_DIR="${BASE_DIR}/results"
AGGREGATION_SCRIPT="${BASE_DIR}/scripts/1_association_pipeline/aggregate_was_results.R"

# All analysis combinations
ANALYSIS_TYPES=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS" "6_zimWAS")
NORMALIZATIONS=("clr" "lognorm" "none")

echo "============================================================================="
echo "WAS Results Aggregation - Failed Aggregations Only"
echo "============================================================================="
echo "Checking for analyses that need aggregation..."
if [[ -n "$SPECIFIC_ANALYSIS" ]]; then
  echo "Specific analysis: $SPECIFIC_ANALYSIS"
fi
if [[ -n "$SPECIFIC_NORMALIZATION" ]]; then
  echo "Specific normalization: $SPECIFIC_NORMALIZATION"
fi
echo "Timestamp: $(date)"
echo "============================================================================="

# Function to check if aggregation is needed
needs_aggregation() {
  local analysis_type=$1
  local normalization=$2
  
  local result_dir="${RESULTS_DIR}/${analysis_type}_out/result_${normalization}"
  local aggregated_file="${result_dir}/${analysis_type}_${normalization}_tidied_complete.rds"
  
  # Check if result directory exists
  if [[ ! -d "$result_dir" ]]; then
    return 1  # No results directory, doesn't need aggregation
  fi
  
  # Check if aggregated file exists and is not empty
  if [[ -f "$aggregated_file" ]] && [[ -s "$aggregated_file" ]]; then
    local file_size=$(stat -f%z "$aggregated_file" 2>/dev/null || stat -c%s "$aggregated_file" 2>/dev/null || echo 0)
    if [[ $file_size -gt 1000 ]]; then
      return 1  # Aggregation already complete
    fi
  fi
  
  # Check if there are individual result files
  local individual_files=$(find "$result_dir" -name "*.rds" ! -name "*_complete.rds" ! -name "*_models.rds" | wc -l)
  
  if [[ $individual_files -gt 0 ]]; then
    return 0  # Has individual files but no aggregation - needs aggregation
  else
    return 1  # No individual files, doesn't need aggregation
  fi
}

# Function to run aggregation for a specific analysis
run_single_aggregation() {
  local analysis_type=$1
  local normalization=$2
  
  echo ""
  echo "🔧 Running aggregation for ${analysis_type}_${normalization}..."
  
  # Check if aggregation script exists
  if [[ ! -f "$AGGREGATION_SCRIPT" ]]; then
    echo "  ❌ Aggregation script not found: $AGGREGATION_SCRIPT"
    return 1
  fi
  
  # Run aggregation
  local cmd="Rscript $AGGREGATION_SCRIPT $analysis_type $normalization $RESULTS_DIR"
  echo "  Command: $cmd"
  
  if eval "$cmd"; then
    echo "  ✅ Aggregation successful for ${analysis_type}_${normalization}"
    return 0
  else
    echo "  ❌ Aggregation failed for ${analysis_type}_${normalization}"
    return 1
  fi
}

# Main logic
aggregations_needed=0
aggregations_successful=0
aggregations_failed=0

echo ""
echo "Scanning all analyses for aggregation needs..."

for analysis_type in "${ANALYSIS_TYPES[@]}"; do
  # Skip if specific analysis requested and this isn't it
  if [[ -n "$SPECIFIC_ANALYSIS" && "$analysis_type" != "$SPECIFIC_ANALYSIS" ]]; then
    continue
  fi
  
  for normalization in "${NORMALIZATIONS[@]}"; do
    # Skip if specific normalization requested and this isn't it
    if [[ -n "$SPECIFIC_NORMALIZATION" && "$normalization" != "$SPECIFIC_NORMALIZATION" ]]; then
      continue
    fi
    
    echo -n "  Checking ${analysis_type}_${normalization}... "
    
    if needs_aggregation "$analysis_type" "$normalization"; then
      echo "NEEDS AGGREGATION"
      ((aggregations_needed++))
      
      if run_single_aggregation "$analysis_type" "$normalization"; then
        ((aggregations_successful++))
      else
        ((aggregations_failed++))
      fi
    else
      echo "OK (already aggregated or no individual files)"
    fi
  done
done

echo ""
echo "============================================================================="
echo "Aggregation Summary"
echo "============================================================================="
echo "Analyses needing aggregation: $aggregations_needed"
echo "Successful aggregations: $aggregations_successful"
echo "Failed aggregations: $aggregations_failed"
echo "Completion time: $(date)"

if [[ $aggregations_needed -eq 0 ]]; then
  echo ""
  echo "No analyses need aggregation - all are up to date!"
elif [[ $aggregations_failed -eq 0 ]]; then
  echo ""
  echo "All needed aggregations completed successfully!"
  echo ""
  echo "Generated files (with FDR and Bonferroni correction):"
  
  for analysis_type in "${ANALYSIS_TYPES[@]}"; do
    if [[ -n "$SPECIFIC_ANALYSIS" && "$analysis_type" != "$SPECIFIC_ANALYSIS" ]]; then
      continue
    fi
    
    for normalization in "${NORMALIZATIONS[@]}"; do
      if [[ -n "$SPECIFIC_NORMALIZATION" && "$normalization" != "$SPECIFIC_NORMALIZATION" ]]; then
        continue
      fi
      
      local aggregated_file="${RESULTS_DIR}/${analysis_type}_out/result_${normalization}/${analysis_type}_${normalization}_tidied_complete.rds"
      if [[ -f "$aggregated_file" ]]; then
        local file_size=$(du -h "$aggregated_file" | cut -f1)
        echo "  ✓ ${analysis_type}_${normalization}_tidied_complete.rds ($file_size)"
      fi
    done
  done
  
  echo ""
  echo "Use p.value.fdr < 0.05 for significance testing."
else
  echo ""
  echo "⚠️  Some aggregations failed. Check error messages above."
  echo "You may need to:"
  echo "  1. Check if individual result files exist"
  echo "  2. Verify R environment is properly loaded"
  echo "  3. Ensure sufficient disk space"
fi

echo "============================================================================="

exit $aggregations_failed 