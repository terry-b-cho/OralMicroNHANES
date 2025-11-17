#!/bin/bash
#  Schema Structure File Creation Script
# Generates mapping files for all 24 WAS analyses (6 types × 4 transformations)
# using  database with no pre-filtering

BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho
DB_PATH="$BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
OUT_DIR="$BASE/results/0_ss_files"

echo " Schema Structure Creation Starting at $(date)"
echo "   Database: $DB_PATH"
echo "   Output: $OUT_DIR"

# Create output directory
mkdir -p "$OUT_DIR"
mkdir -p "$BASE/logs"

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

echo "Submitting 24 schema structure creation jobs..."

# Submit jobs for each combination
for analysis in "${analysis_types[@]}"; do
  for transform in "${transformations[@]}"; do
    
    # Define table names for this transformation (all use RELATIVE naming)
    OTU_F="DADA2RSV_GENUS_RELATIVE_F_${transform}"
    OTU_G="DADA2RSV_GENUS_RELATIVE_G_${transform}"
    
    # Get OTU role for this analysis
    otu_role="${otu_roles[$analysis]}"
    
    # Configuration file
    vars_file="$BASE/configs/${analysis}_vars.txt"
    
    # Job name
    job_name="ss_${analysis}_${transform}"
    
    echo "  Submitting: $job_name"
    
    # Create temporary job script
    job_script="$BASE/logs/temp_${job_name}.sh"
    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$BASE/logs/ss_${analysis}_${transform}_%j.out
#SBATCH --error=$BASE/logs/ss_${analysis}_${transform}_%j.err
#SBATCH --partition=short
#SBATCH --time=0-01:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Load conda environment
module load conda/miniforge3/24.11.3-0
eval "\$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Run  schema structure creation
Rscript $BASE/scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R \
  --db '$DB_PATH' \
  --otu_F '$OTU_F' \
  --otu_G '$OTU_G' \
  --vars_file '$vars_file' \
  --otu_role '$otu_role' \
  --pipeline '$analysis' \
  --transform '$transform' \
  --out_dir '$OUT_DIR'
EOF

    # Submit the job script
    sbatch "$job_script"
    
    # Clean up temporary script after a short delay
    (sleep 5 && rm -f "$job_script") &
    
  done
done

echo ""
echo "All 24 jobs submitted. Monitor with:"
echo "  squeue -u \$USER"
echo ""
echo "Expected output files in: $OUT_DIR"
echo "  • 1_demoWAS_clr_schema_structure.csv"
echo "  • 2_oradWAS_hellinger_schema_structure.csv"
echo "  • ... (24 total files)"
echo ""
echo "After completion, proceed to WAS analyses with  pipeline." 