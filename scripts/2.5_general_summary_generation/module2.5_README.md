# Module 2.5 — General summary generation

## What this module does
Generates publication-ready supplementary summary tables for the NHANES oral microbiome study.

- `general_summary.R` — main WAS variable summaries, binary distributions with Wilson 95% CIs, alpha-diversity stratification by demographic group, variable-type distribution by analysis.
- `additional_supplementary_generation.R` — full per-variable distribution table, alpha + beta diversity group means (with and without Wilcoxon vs reference + BH-FDR), and a host × OTU regression supplementary table (with a CLR-only manuscript subset).

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite`
- `data/00_nhanes_omp_diversity_db/{dada2rsv-alpha.txt, dada2rsv-beta-{braycurtis,unwunifrac,wunifrac}.txt}`
- `configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt`
- `results/analyses_results/4_otu_host_plot_out/intermediate/{tidied_results_with_otu, taxonomy_annotations, ubiome_variable_mapping, prevalence_filters, binary_case_filters}.rds`

## Outputs (relative to `PROJECT_ROOT`)

`results/analyses_results/2.5_general_summary_generation/out/` *(from `general_summary.R`)*:
- `Table_S1_Complete_WAS_Variable_Summary.csv`, `Table_S2_Binary_Variables_Analysis.csv`, `Table_S3_Alpha_Diversity_Stratification.csv`, `Table_S4_Variable_Type_Distribution.csv`
- `table_s1_complete_was_variable_summary.csv`, `table_s2_was_analysis_summary.csv`, `table_s3_variable_type_distribution.csv`, `table_s4_binary_variable_analysis.csv`

`results/analyses_results/2.5_general_summary_generation/supplementary/` *(from `additional_supplementary_generation.R`)*:
- `full_list_supplementary_variable_distribution.csv`
- `full_list_supplementary_alpha_beta_variable_distribution_per_group_per_level.csv`
- `full_list_supplementary_alpha_beta_variable_distribution_per_group_per_level_compared_to_ref.csv`
- `full_list_supplementary_host_otu_regression_res.csv`
- `full_list_supplementary_host_otu_regression_res_subset_clr.csv`

## Scripts
- `general_summary.R`
- `additional_supplementary_generation.R`

## Environment
R >= 4.5 with: `DBI`, `RSQLite`, `dbplyr`, `dplyr`, `tidyr`, `readr`, `purrr`, `stringr`, `data.table`, `survey`, `knitr`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
Open each `.R` file and set `PROJECT_ROOT` (top) to your local clone of this repository. From the repo root:

```bash
Rscript scripts/2.5_general_summary_generation/general_summary.R
Rscript scripts/2.5_general_summary_generation/additional_supplementary_generation.R
```

## Run order
Both scripts are independent — run in any order.
