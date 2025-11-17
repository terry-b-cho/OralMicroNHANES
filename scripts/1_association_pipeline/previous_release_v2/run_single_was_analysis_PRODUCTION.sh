#!/bin/bash

#  PER-DEPENDENT-VARIABLE SUBMISSION
# Matches original pipeline: one job per dependent variable
# Usage: bash run_per_dependent_variable.sh <analysis_type> <normalization> [test]

if [ $# -lt 2 ]; then
    echo "Usage: $0 <analysis_type> <normalization> [test]"
    echo "analysis_type: 1_demoWAS, 2_oradWAS, 3_exWAS, 4_pheWAS, 5_outWAS, 6_zimWAS"
    echo "normalization: clr, lognorm, none"
    echo "test: optional flag for test mode"
    exit 1
fi

ANALYSIS_TYPE=$1
NORMALIZATION=$2
TEST_MODE=${3:-""}

BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho
DB_PATH="$BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
SCHEMA_FILE="$BASE/results/0_ss_files/${ANALYSIS_TYPE}_${NORMALIZATION}_schema_structure.csv"
OUTPUT_DIR="$BASE/results/${ANALYSIS_TYPE}_out/result_${NORMALIZATION}"

echo " PER-DEPENDENT-VARIABLE SUBMISSION"
echo "============================================="
echo "Analysis: $ANALYSIS_TYPE"
echo "Normalization: $NORMALIZATION"
echo "Test mode: $TEST_MODE"
echo "Schema file: $SCHEMA_FILE"
echo "Output dir: $OUTPUT_DIR"

# Check if schema file exists
if [ ! -f "$SCHEMA_FILE" ]; then
    echo "❌ Schema file not found: $SCHEMA_FILE"
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

# Memory allocation based on analysis type
case $ANALYSIS_TYPE in
    3_exWAS)
        MEMORY="8G"
        TIME="0-02:00"
        ;;
    4_pheWAS)
        MEMORY="4G"
        TIME="0-01:00"
        ;;
    *)
        MEMORY="2G"
        TIME="0-00:30"
        ;;
esac

if [ "$TEST_MODE" = "test" ]; then
    MEMORY="2G"
    TIME="0-00:15"
fi

echo "Resource allocation: $MEMORY memory, $TIME time"
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
    
    # Submit SLURM job
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$BASE/logs/${job_name}_%j.out
#SBATCH --error=$BASE/logs/${job_name}_%j.err
#SBATCH --partition=short
#SBATCH --time=$TIME
#SBATCH --mem=$MEMORY
#SBATCH --cpus-per-task=2

# Load environment
module load miniconda3/23.1.0
source activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Change to working directory
cd $BASE

# Run analysis for this dependent variable
echo "Starting analysis for $dep_var at \$(date)"

if [ "$TEST_MODE" = "test" ]; then
    Rscript scripts/1_association_pipeline/universal_was_analysis_PRODUCTION.R \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$SCHEMA_FILE" \\
        --database_path "$DB_PATH" \\
        --output_path "$OUTPUT_DIR" \\
        --analysis_type "$ANALYSIS_TYPE" \\
        --normalization "$NORMALIZATION" \\
        --test
else
    Rscript scripts/1_association_pipeline/universal_was_analysis_PRODUCTION.R \\
        --dependent_var "$dep_var" \\
        --schema_structure_file "$SCHEMA_FILE" \\
        --database_path "$DB_PATH" \\
        --output_path "$OUTPUT_DIR" \\
        --analysis_type "$ANALYSIS_TYPE" \\
        --normalization "$NORMALIZATION"
fi

echo "Analysis for $dep_var completed at \$(date)"
EOF

    # Brief pause between submissions
    sleep 0.5
done

echo ""
echo "✅ Submitted $TOTAL_VARS jobs for $ANALYSIS_TYPE $NORMALIZATION"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER | grep ${ANALYSIS_TYPE}_${NORMALIZATION}"
echo ""
echo "Expected output files:"
echo "  $OUTPUT_DIR/*.rds (one per dependent variable)"
echo ""
echo "After completion, run aggregation script to combine results." 
