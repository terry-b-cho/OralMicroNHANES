#!/bin/bash
# LOCAL Schema Structure File Creation Script
# Generates mapping files for all 24 WAS analyses (6 types × 4 transformations)
# LOCAL VERSION - No SLURM dependencies

BASE=$(pwd)
DB_PATH="$BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
OUT_DIR="$BASE/results/0_ss_files"

echo "LOCAL Schema Structure Creation Starting at $(date)"
echo "   Database: $DB_PATH"
echo "   Output: $OUT_DIR"

# Check if database exists
if [ ! -f "$DB_PATH" ]; then
    echo "❌ ERROR: Database not found at $DB_PATH"
    exit 1
fi

# Create output directory
mkdir -p "$OUT_DIR"

# Analysis configurations
declare -a analysis_types=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS" "6_zimWAS")
declare -a transformations=("none" "hellinger" "clr" "lognorm")

# OTU role mapping
declare -A otu_roles
otu_roles["1_demoWAS"]="dep"    # Demographics → Microbiome
otu_roles["2_oradWAS"]="indep"  # Microbiome → Oral Health
otu_roles["3_exWAS"]="dep"      # Exposures → Microbiome
otu_roles["4_pheWAS"]="indep"   # Microbiome → Phenotypes
otu_roles["5_outWAS"]="indep"   # Microbiome → Disease Outcomes
otu_roles["6_zimWAS"]="indep"   # Microbiome → Lab Measurements

echo "Creating 24 schema structure files..."

# Track progress
total_files=24
current_file=0

# Process each combination
for analysis in "${analysis_types[@]}"; do
  for transform in "${transformations[@]}"; do
    
    current_file=$((current_file + 1))
    
    # Check if output file already exists
    output_file="$OUT_DIR/${analysis}_${transform}_schema_structure.csv"
    if [ -f "$output_file" ]; then
        echo "  [$current_file/$total_files] SKIP: $analysis × $transform (already exists)"
        continue
    fi
    
    # Define table names for this transformation
    OTU_F="DADA2RSV_GENUS_RELATIVE_F_${transform}"
    OTU_G="DADA2RSV_GENUS_RELATIVE_G_${transform}"
    
    # Get OTU role for this analysis
    otu_role="${otu_roles[$analysis]}"
    
    # Configuration file
    vars_file="$BASE/configs/${analysis}_vars.txt"
    
    # Check if config file exists
    if [ ! -f "$vars_file" ]; then
        echo "  ❌ [$current_file/$total_files] ERROR: Config file not found: $vars_file"
        continue
    fi
    
    echo "  ⚙️  [$current_file/$total_files] Processing: $analysis × $transform (role: $otu_role)"
    
    # Run schema structure creation
    Rscript "$BASE/scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R" \
      --db "$DB_PATH" \
      --otu_F "$OTU_F" \
      --otu_G "$OTU_G" \
      --vars_file "$vars_file" \
      --otu_role "$otu_role" \
      --pipeline "$analysis" \
      --transform "$transform" \
      --out_dir "$OUT_DIR"
    
    # Check if the file was created successfully
    if [ -f "$output_file" ]; then
        # Get file size and row count
        file_size=$(du -h "$output_file" | cut -f1)
        row_count=$(wc -l < "$output_file")
        echo "  ✅ [$current_file/$total_files] SUCCESS: Created $output_file ($file_size, $row_count rows)"
    else
        echo "  ❌ [$current_file/$total_files] FAILED: $output_file not created"
    fi
    
  done
done

echo ""
echo "Schema Structure Creation Complete!"
echo "   Output directory: $OUT_DIR"
echo "   Expected files: 24"

# Count actual files created
created_files=$(ls -1 "$OUT_DIR"/*_schema_structure.csv 2>/dev/null | wc -l)
echo "   Created files: $created_files"

# List all files
echo ""
echo "Files in $OUT_DIR:"
ls -lh "$OUT_DIR"/*_schema_structure.csv 2>/dev/null || echo "   No schema structure files found"

if [ "$created_files" -eq 24 ]; then
    echo ""
    echo "ALL 24 SCHEMA STRUCTURE FILES CREATED SUCCESSFULLY!"
else
    echo ""
    echo "⚠️  Some files may be missing. Expected: 24, Found: $created_files"
fi 