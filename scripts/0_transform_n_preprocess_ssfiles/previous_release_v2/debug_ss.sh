#!/bin/bash
# Debug script for  schema structure creation
# Run this manually on O2 to identify the issue

echo " Schema Structure Debug Test"
echo "========================================="

# Set paths
BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho
DB_PATH="$BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
OUT_DIR="$BASE/results/0_ss_files"

echo "BASE: $BASE"
echo "DB_PATH: $DB_PATH"
echo "OUT_DIR: $OUT_DIR"
echo ""

# Check 1: Database exists
echo "1. Checking database..."
if [ -f "$DB_PATH" ]; then
    echo "   ✓ Database exists: $DB_PATH"
    echo "   Size: $(ls -lh "$DB_PATH" | awk '{print $5}')"
else
    echo "   ✗ Database NOT found: $DB_PATH"
    exit 1
fi

# Check 2: Output directory
echo ""
echo "2. Creating output directory..."
mkdir -p "$OUT_DIR"
if [ -d "$OUT_DIR" ]; then
    echo "   ✓ Output directory ready: $OUT_DIR"
else
    echo "   ✗ Cannot create output directory: $OUT_DIR"
    exit 1
fi

# Check 3: Config files
echo ""
echo "3. Checking config files..."
config_file="$BASE/configs/3_exWAS_vars.txt"
if [ -f "$config_file" ]; then
    echo "   ✓ Config exists: $config_file"
    echo "   Variables: $(wc -l < "$config_file")"
    echo "   First 5 variables:"
    head -5 "$config_file" | sed 's/^/     /'
else
    echo "   ✗ Config missing: $config_file"
    exit 1
fi

# Check 4: Load conda environment and run
echo ""
echo "4. Loading conda environment and running test..."
# Load miniconda module and activate environment
module load miniconda3/23.1.0
source activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

echo "   ✓ Conda environment activated"
echo "   ✓ R version: $(R --version | head -1)"

echo "Running  schema structure creation..."
Rscript "$BASE/scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R" \
  --db "$DB_PATH" \
  --otu_F "DADA2RSV_GENUS_RELATIVE_F_clr" \
  --otu_G "DADA2RSV_GENUS_RELATIVE_G_clr" \
  --vars_file "$config_file" \
  --otu_role "dep" \
  --pipeline "3_exWAS" \
  --transform "clr" \
  --out_dir "$OUT_DIR"

exit_code=$?
echo ""
echo "Exit code: $exit_code"

# Check output
expected_output="$OUT_DIR/3_exWAS_clr_schema_structure.csv"
if [ -f "$expected_output" ]; then
    echo "✅ SUCCESS! Output file created:"
    echo "   File: $expected_output"
    echo "   Size: $(ls -lh "$expected_output" | awk '{print $5}')"
    echo "   Rows: $(wc -l < "$expected_output")"
else
    echo "❌ FAILED! Output file not created"
    echo "Expected: $expected_output"
fi
