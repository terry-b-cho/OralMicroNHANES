# Module 8 — Microbial association network analyses

## What this module does
Builds host-microbe association networks (nodes = OTU/host variables, edges = significant WAS associations) and renders them as PDFs. Two parallel pipelines:
- **Baseline** (`*.R`) — writes to `8_network_analyses_out/`
- **Additional** (`*_additional.R`) — writes to `8_network_analyses_out_additional/` (stricter dual-threshold significance gate: BH-FDR ≤ 0.05 AND Storey q ≤ 0.05)

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` *(or `..._complete.sqlite` fallback for the additional pipeline)*
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_*.rds`
- `results/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_out/result_{none,clr}/*_tidied_complete.rds`
- `configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt`
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv`

## Outputs

### Baseline pipeline → `results/analyses_results/8_network_analyses_out/`
- `intermediate/{phyloseq_objects, sample_data_subsets, ubiome_genus_mapping_complete, significant_results, variable_description, otu_pass_prevalence_list, variable_groups, *_case_*}.rds`
- `network_plots_final/<mode>/<transformation>/Network_*.pdf`
- `group_networks/<mode>/<group>/*.pdf`

### Additional pipeline → `results/analyses_results/8_network_analyses_out_additional/`
Same layout as above, plus `supplementary_tables/` with the per-edge association tables.

## Scripts
- `1_load_all_categories.R` / `1_load_all_categories_additional.R` — build intermediates (phyloseq + WAS results + filtered significant tables)
- `2_network.R` / `2_network_additional.R` — render per-variable network PDFs
- `2.1_group_network.R` / `2.1_group_network_additional.R` — render grouped (umbrella) network PDFs

## Environment
R >= 4.5 with: `phyloseq`, `DBI`, `RSQLite`, `dplyr`, `dbplyr`, `tidyr`, `stringr`, `purrr`, `tibble`, `readr`, `igraph`, `ggraph`, `ggplot2`, `gridExtra`, `extrafont`, `ggrepel`, `scales`, `rlang`, `grid`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open each `.R` file and set `PROJECT_ROOT` (top) to your local clone of this repository.
2. From the repo root — choose either pipeline:

   ```bash
   # Baseline pipeline (writes to 8_network_analyses_out/)
   Rscript scripts/8_network_analyses/1_load_all_categories.R
   Rscript scripts/8_network_analyses/2_network.R
   Rscript scripts/8_network_analyses/2.1_group_network.R

   # Additional pipeline (writes to 8_network_analyses_out_additional/)
   Rscript scripts/8_network_analyses/1_load_all_categories_additional.R
   Rscript scripts/8_network_analyses/2_network_additional.R
   Rscript scripts/8_network_analyses/2.1_group_network_additional.R
   ```

## Run order
Within each pipeline: `1 → 2` and `1 → 2.1`. The two plotting scripts (`2` and `2.1`) are independent of each other; either can run after `1`. The baseline and additional pipelines are also independent — they read different inputs in some cases and write to different output dirs.
