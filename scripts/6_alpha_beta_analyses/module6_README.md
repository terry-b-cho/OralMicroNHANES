# Module 6 — Alpha/beta diversity analyses

## What this module does
Loads NHANES oral-microbiome alpha-diversity metrics and 24+ categorical metadata variables, computes PCoA / centroid distances on three beta-diversity matrices (Bray-Curtis, unweighted UniFrac, weighted UniFrac), runs Kruskal-Wallis + pairwise Wilcoxon tests per categorical variable, and renders an integrated 5×3 per-variable panel (alpha distributions, beta distributions, beta-centroid distributions, PCoA scatter, scree). Also emits a supplementary categorical-variable definitions table.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt`
- `data/00_nhanes_omp_diversity_db/dada2rsv-beta-braycurtis.txt`
- `data/00_nhanes_omp_diversity_db/dada2rsv-beta-unwunifrac.txt`
- `data/00_nhanes_omp_diversity_db/dada2rsv-beta-wunifrac.txt`
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`

## Outputs
- `scripts/6_alpha_beta_analyses/data/` *(intermediate, auto-created by script 1)*:
  - `all_categorical_data.rds`, `alpha_diversity_with_all_categories.rds`, `categorical_variables_info.csv`
- `results/analyses_results/6_alpha_beta_analyses_out/`:
  - `integrated_diversity_plots/Integrated_<var>[_<mode>].pdf` (one per categorical variable, per label mode)
  - `alpha_pairwise_stat_res.csv`, `beta_pairwise_stat_res.csv`, `beta_centroid_PCoA_pairwise_stat_res.csv`
- `scripts/6_alpha_beta_analyses/supplementary_table_categorical_variable_definitions.{csv,md}`

## Scripts
- `1_load_diversity_and_all_categories_HPC.R` — pull alpha + categorical metadata; write intermediate RDS/CSVs
- `2_integrated_diversity_FINAL_HPC_fixed.R` — PCoA + centroid + pairwise tests + integrated plots; sources script 1 internally
- `2_run_all_modes_half_HPC.sh` — HPC wrapper that runs script 2 in three label modes (full P-values, asterisks, no brackets)
- `3_create_categorical_variable_supplementary_table.R` — supplementary categorical-variable definitions table

## Environment
R >= 4.5 with: `dplyr`, `tidyr`, `readr`, `forcats`, `tibble`, `data.table`, `DBI`, `RSQLite`, `dbplyr`, `ggplot2`, `ape`, `egg`, `cowplot`, `ggdist`, `ggsignif`, `knitr`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

On HPC (Slurm / O2):
```bash
module load conda/miniforge3/24.11.3-0
conda activate $PROJECT_ROOT/envs/.conda/envs/nhanes-analysis
```

## How to run
1. Open each `.R` and the `.sh` file and set `PROJECT_ROOT` (at the top) to the absolute path of your local clone of this repository.
2. From the repo root:
   ```bash
   Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories_HPC.R
   Rscript scripts/6_alpha_beta_analyses/2_integrated_diversity_FINAL_HPC_fixed.R
   # or, for all three label modes:
   # bash scripts/6_alpha_beta_analyses/2_run_all_modes_half_HPC.sh
   Rscript scripts/6_alpha_beta_analyses/3_create_categorical_variable_supplementary_table.R
   ```

## Run order
`1 → 2` (or the shell wrapper) `→ 3`. Script 2 sources script 1 internally, so step 1 is optional if you go straight to step 2. Step 3 requires the intermediate `data/` dir, which is created by step 1 (or by step 2 when it sources script 1).
