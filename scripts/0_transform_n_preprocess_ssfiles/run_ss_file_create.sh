#!/bin/bash
#
# Step 3: submit 24 schema-structure creation jobs to SLURM, one per
# (analysis x transformation) combination.
#
# Environment: each SLURM job activates the project conda env at
# envs/.conda/envs/nhanes-analysis. Callers do not need to pre-load anything.

# === USER SETTING =============================================================
# Set PROJECT_ROOT to the absolute path of your local clone of this repository.
PROJECT_ROOT=/n/groups/patel/terry/nhanes_oral_mirco_cho
# ==============================================================================

DB_PATH="$PROJECT_ROOT/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite"
OUT_DIR="$PROJECT_ROOT/results/0_ss_files"

mkdir -p "$OUT_DIR" "$PROJECT_ROOT/logs"

declare -a analysis_types=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS")
declare -a transformations=("none" "hellinger" "clr" "lognorm")

declare -A otu_roles
otu_roles["1_demoWAS"]="dep"
otu_roles["2_oradWAS"]="indep"
otu_roles["3_exWAS"]="dep"
otu_roles["4_pheWAS"]="indep"
otu_roles["5_outWAS"]="indep"

for analysis in "${analysis_types[@]}"; do
  for transform in "${transformations[@]}"; do
    OTU_F="DADA2RSV_GENUS_RELATIVE_F_${transform}"
    OTU_G="DADA2RSV_GENUS_RELATIVE_G_${transform}"
    otu_role="${otu_roles[$analysis]}"
    vars_file="$PROJECT_ROOT/configs/${analysis}_vars.txt"
    job_name="ss_${analysis}_${transform}"
    job_script="$PROJECT_ROOT/logs/temp_${job_name}.sh"

    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=$PROJECT_ROOT/logs/ss_${analysis}_${transform}_%j.out
#SBATCH --error=$PROJECT_ROOT/logs/ss_${analysis}_${transform}_%j.err
#SBATCH --partition=short
#SBATCH --time=0-01:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load conda/miniforge3/24.11.3-0
eval "\$(conda shell.bash hook)"
conda activate $PROJECT_ROOT/envs/.conda/envs/nhanes-analysis

Rscript $PROJECT_ROOT/scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R \\
  --db '$DB_PATH' \\
  --otu_F '$OTU_F' \\
  --otu_G '$OTU_G' \\
  --vars_file '$vars_file' \\
  --otu_role '$otu_role' \\
  --pipeline '$analysis' \\
  --transform '$transform' \\
  --out_dir '$OUT_DIR'
EOF

    sbatch "$job_script"
    (sleep 5 && rm -f "$job_script") &
  done
done
