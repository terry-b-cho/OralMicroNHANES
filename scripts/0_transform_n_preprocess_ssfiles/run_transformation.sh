#!/bin/bash
#SBATCH -p short      # jobs < 12 h
#SBATCH -t 0-02:00    # extended time for full data
#SBATCH --mem=64G     # more memory for larger datasets
#SBATCH -c 8
#SBATCH -o logs/trans_%j.out
#SBATCH -e logs/trans_%j.err

module purge
module load gcc/14.2.0 R/4.4.2

BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho

mkdir -p $BASE/logs

echo " Pipeline Transformation Starting at $(date)"
echo "   • NO prevalence pre-filtering"
echo "   • Four transformations: none, hellinger, clr, lognorm"
echo "   • Universal pseudocount ε = 0.001 for COUNT matrices only"

Rscript $BASE/scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R \
  --in_db  $BASE/data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite \
  --out_db $BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite

echo " transformation completed at $(date)"
echo "Next step: Run missing data filling script" 