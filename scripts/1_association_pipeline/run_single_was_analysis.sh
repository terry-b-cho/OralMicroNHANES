#!/bin/bash

# PER-DEPENDENT-VARIABLE SUBMISSION - O2 COMPATIBLE
# Matches original pipeline: one job per dependent variable
# UPDATED: O2 module compatibility, 4 transformations, enhanced error handling
# Usage: bash run_per_dependent_variable.sh <analysis_type> <normalization> [test]

if [ $# -lt 2 ]; then
    echo "Usage: $0 <analysis_type> <normalization> [test]"
    echo "analysis_type: 1_demoWAS, 2_oradWAS, 3_exWAS, 4_pheWAS, 5_outWAS, 6_zimWAS"
    echo "normalization: clr, lognorm, none, hellinger"  # Updated to include hellinger
    echo "test: optional flag for test mode"
    echo ""
    echo "ENHANCED FEATURES:"
    echo "  • O2 module compatibility (gcc/14.2.0, conda/miniforge3/24.11.3-0)"
    echo "  • 4 transformation methods including Hellinger"
    echo "  • Final technical audit corrections applied"
    echo "  • Enhanced error handling and table/view detection"
    exit 1
fi

ANALYSIS_TYPE=$1
NORMALIZATION=$2
TEST_MODE=${3:-""}

BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho
DB_PATH="$BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
SCHEMA_FILE="$BASE/results/0_ss_files/${ANALYSIS_TYPE}_${NORMALIZATION}_schema_structure.csv"
OUTPUT_DIR="$BASE/results/${ANALYSIS_TYPE}_out/result_${NORMALIZATION}"

echo "PER-DEPENDENT-VARIABLE SUBMISSION - FINAL AUDIT CORRECTED"
echo "=============================================================="
echo "Analysis: $ANALYSIS_TYPE"
echo "Normalization: $NORMALIZATION"
echo "Test mode: $TEST_MODE"
echo "Schema file: $SCHEMA_FILE"
echo "Output dir: $OUTPUT_DIR"
echo ""
echo "ENHANCEMENTS APPLIED:"
echo "  ✅ O2 module compatibility"
echo "  ✅ 4 transformations (none, hellinger, clr, lognorm)"
echo "  ✅ Table/view detection fixes"
echo "  ✅ Effect scale harmonization"
echo "  ✅ Survey-weighted statistics with proper NCHS compliance"

# Validate normalization method
case $NORMALIZATION in
    clr|lognorm|none|hellinger)
        echo "✅ Valid normalization method: $NORMALIZATION"
        ;;
    *)
        echo "❌ Invalid normalization method: $NORMALIZATION"
        echo "Valid options: clr, lognorm, none, hellinger"
        exit 1
        ;;
esac

# Check if schema file exists
if [ ! -f "$SCHEMA_FILE" ]; then
    echo "❌ Schema file not found: $SCHEMA_FILE"
    echo ""
    echo "Available schema files:"
    ls "$BASE/results/0_ss_files/${ANALYSIS_TYPE}_"*"_schema_structure.csv" 2>/dev/null || echo "  None found for $ANALYSIS_TYPE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$BASE/logs"

# Get list of unique dependent variables
echo "Getting list of dependent variables..."
DEPENDENT_VARS=$(cut -d',' -f1 "$SCHEMA_FILE" | tail -n +2 | sort | uniq)
TOTAL_VARS=$(echo "$DEPENDENT_VARS" | wc -l)

echo "Found $TOTAL_VARS unique dependent variables"

if [ "$TEST_MODE" = "test" ]; then
    echo "Test mode: Processing first 3 dependent variables only"
    DEPENDENT_VARS=$(echo "$DEPENDENT_VARS" | head -3)
    TOTAL_VARS=3
fi

# ENHANCED resource allocation based on computational analysis (from technical review)
case $ANALYSIS_TYPE in
    2_oradWAS|4_pheWAS|5_outWAS|6_zimWAS)
        # Microbiome as independent variable: higher computational load
        MEMORY="12G"
        TIME="0-04:00"
        echo "Resource allocation: $MEMORY memory, $TIME time (microbiome as predictor)"
        ;;
    3_exWAS)
        # Microbiome as dependent variable: moderate computational load
        MEMORY="8G"
        TIME="0-03:00"
        echo "Resource allocation: $MEMORY memory, $TIME time (microbiome as outcome)"
        ;;
    1_demoWAS)
        # Demographics to microbiome: lowest computational load
        MEMORY="1G"
        TIME="0-01:00"
        echo "Resource allocation: $MEMORY memory, $TIME time (demographics analysis)"
        ;;
    *)
        # Default allocation
        MEMORY="2G"
        TIME="0-00:30"
        echo "Resource allocation: $MEMORY memory, $TIME time (default)"
        ;;
esac

if [ "$TEST_MODE" = "test" ]; then
    MEMORY="2G"
    TIME="0-00:15"
    echo "Test mode override: $MEMORY memory, $TIME time"
fi

echo ""

# Submit job for each dependent variable
counter=0
for dep_var in $DEPENDENT_VARS; do
    counter=$((counter + 1))
    
    # Clean variable name for job naming
    clean_var=$(echo "$dep_var" | sed 's/[^a-zA-Z0-9_]/_/g')
    job_name="${ANALYSIS_TYPE}_${NORMALIZATION}_${clean_var}"
    
    if [ "$TEST_MODE" = "test" ]; then
        job_name="${job_name}_test"
    fi
    
    echo "[$counter/$TOTAL_VARS] Submitting: $dep_var"
    
    # Submit SLURM job with O2 compatibility
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$BASE/logs/${job_name}_%j.out
#SBATCH --error=$BASE/logs/${job_name}_%j.err
#SBATCH --partition=short
#SBATCH --time=$TIME
#SBATCH --mem=$MEMORY
#SBATCH --cpus-per-task=2

# Load environment with O2 compatibility (CRITICAL: Updated modules)
module purge
module load gcc/14.2.0                 # Updated with O2 update
module load conda/miniforge3/24.11.3-0 # Updated with O2 update

# Initialize conda for bash and activate environment (modern syntax)
eval "\$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Change to working directory
cd $BASE

# Run analysis for this dependent variable with final audit corrections
echo "Starting FINAL AUDIT CORRECTED analysis for $dep_var at \$(date)"
echo "Analysis: $ANALYSIS_TYPE, Normalization: $NORMALIZATION"
echo "Enhancements: Table/view detection, effect harmonization, survey weights"

if [ "$TEST_MODE" = "test" ]; then
    Rscript scripts/1_association_pipeline/universal_was_analysis.R \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$SCHEMA_FILE" \\
        --database_path "$DB_PATH" \\
        --output_path "$OUTPUT_DIR" \\
        --analysis_type "$ANALYSIS_TYPE" \\
        --normalization "$NORMALIZATION" \\
        --test
else
    Rscript scripts/1_association_pipeline/universal_was_analysis.R \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$SCHEMA_FILE" \\
        --database_path "$DB_PATH" \\
        --output_path "$OUTPUT_DIR" \\
        --analysis_type "$ANALYSIS_TYPE" \\
        --normalization "$NORMALIZATION"
fi

echo "Analysis for $dep_var completed at \$(date)"
EOF

    # Brief pause between submissions to avoid overwhelming scheduler
    sleep 0.5
done

echo ""
echo "✅ Submitted $TOTAL_VARS jobs for $ANALYSIS_TYPE $NORMALIZATION"
echo ""
echo "MONITORING:"
echo "  Check running jobs: squeue -u \$USER | grep ${ANALYSIS_TYPE}_${NORMALIZATION}"
echo "  View logs: ls $BASE/logs/${ANALYSIS_TYPE}_${NORMALIZATION}_*.out"
echo "  Check for errors: find $BASE/logs -name '${ANALYSIS_TYPE}_${NORMALIZATION}_*.err' -size +0"
echo ""
echo "EXPECTED OUTPUT:"
echo "  Files: $OUTPUT_DIR/*.rds (one per dependent variable)"
echo "  Enhanced features: Effect scales, survey statistics"
echo ""
echo "NEXT STEPS:"
echo "  1. Wait for all jobs to complete"
echo "  2. Run aggregation: Rscript scripts/1_association_pipeline/aggregate_was_results.R $ANALYSIS_TYPE $NORMALIZATION \$(pwd)/results"
echo "  3. Apply FDR correction separately after aggregation" 
