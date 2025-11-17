#!/bin/bash
#SBATCH -p short      # jobs < 12 h
#SBATCH -t 0-01:00    # hh:mm:ss
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -o logs/trans_%j.out
#SBATCH -e logs/trans_%j.err

module load gcc/9.2.0 R/4.2.1

BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho

Rscript $BASE/scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R \
  --in_db  $BASE/data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite \
  --out_db $BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --prev   0.001
