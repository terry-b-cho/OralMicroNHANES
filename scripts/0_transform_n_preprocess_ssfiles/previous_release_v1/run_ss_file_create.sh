#!/usr/bin/env bash
# ------------------------------------------------------------------
# Run ss_file_create for all pipelines and transformations
# ------------------------------------------------------------------
module load gcc/9.2.0 R/4.2.1

set -euo pipefail

BASE=/n/groups/patel/terry/nhanes_oral_mirco_cho
DB=$BASE/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
SCRIPT=$BASE/scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R
OUTDIR=$BASE/results/0_ss_files

# List of pipelines (config file prefix)
pipelines=(
  1_demoWAS
  2_oradWAS
  3_exWAS
  4_pheWAS
  5_outWAS
  6_zimWAS
)
# List of transformations
transforms=(none clr lognorm)

for pipeline in "${pipelines[@]}"; do
  vars_file=$BASE/configs/${pipeline}_vars.txt
  # Determine OTU role
  if [ "$pipeline" = "3_exWAS" ]; then
    otu_role=indep
  else
    otu_role=dep
  fi

  for transform in "${transforms[@]}"; do
    otu_F=DADA2RSV_GENUS_RELATIVE_F_${transform}
    otu_G=DADA2RSV_GENUS_RELATIVE_G_${transform}

    echo "[${pipeline}/${transform}] Running schema_structure creation..."
    Rscript "$SCRIPT" \
      --db "$DB" \
      --otu_F "$otu_F" \
      --otu_G "$otu_G" \
      --vars_file "$vars_file" \
      --otu_role "$otu_role" \
      --pipeline "$pipeline" \
      --transform "$transform" \
      --out_dir "$OUTDIR"
    echo "[${pipeline}/${transform}] Done."
    echo
  done
done
