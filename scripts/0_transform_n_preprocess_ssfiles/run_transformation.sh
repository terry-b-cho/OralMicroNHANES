#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-02:00
#SBATCH --mem=64G
#SBATCH -c 8
#SBATCH -o logs/trans_%j.out
#SBATCH -e logs/trans_%j.err
#
# Step 1: submit the microbiome transformation job to SLURM.
#
# Environment: SLURM module `gcc/14.2.0 R/4.4.2`. The job loads these modules
# itself, so callers do not need to pre-load anything.

# === USER SETTING =============================================================
# Set PROJECT_ROOT to the absolute path of your local clone of this repository.
PROJECT_ROOT=/n/groups/patel/terry/nhanes_oral_mirco_cho
# ==============================================================================

module purge
module load gcc/14.2.0 R/4.4.2

mkdir -p "$PROJECT_ROOT/logs"

Rscript "$PROJECT_ROOT/scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R" \
  --in_db  "$PROJECT_ROOT/data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite" \
  --out_db "$PROJECT_ROOT/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite"
