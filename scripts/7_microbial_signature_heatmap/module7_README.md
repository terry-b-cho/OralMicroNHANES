# Module 7 — Microbial signature correlation heatmaps

## What this module does
Builds host-host weighted correlation matrices from per-taxon svyglm β-vectors (inverse-variance weighted Pearson over shared microbes; BH-FDR per dataset), then renders two heatmap styles:
- **legacy "sized"** — cell area scales with FDR tier; clustering uses `1 - |r|`
- **optimal-k splits** — silhouette-driven row/column k splits with configurable OTU / FDR / prevalence filters

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite`
- `results/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_out/result_clr/*_clr_tidied_complete.rds`
- `configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt`

## Outputs (under `results/analyses_results/7_microbial_signature_heatmap_out/`)
- `intermediate/microbial_signature/<analysis>_clr_all_results/{signature_matrix,correlation_matrix}.rds`
- `intermediate/microbial_signature/merged_clr_data/{signature_matrix,correlation_matrix}.rds` + `analysis_contributions.csv`, `composition_summary.csv`, `detailed_correlation_summary.csv`
- `intermediate/microbial_signature/heatmap_data/<dataset>_{correlation,fdr,n_effective,n_shared}_matrix.csv` + `<dataset>_variable_mapping.csv` + `<dataset>_summary_stats.csv`
- `intermediate/microbial_signature/{dataset_summary,complete_verification_results}.csv`
- `intermediate/heatmap_significant_results.rds`, `intermediate/tidied_results_clr.rds`
- Legacy sized: `<analysis>_clr_all_results_conventional_heatmap_sized.pdf` + `0_weighted_correlation_legend.pdf` + `dataset_summary.csv`
- Optimal-k: `otu<N>_fdr<X>_prev<Y>_k<K>_<analysis>_clr_siges_heatmap.pdf` + matching `..._weighted_correlation_legend.pdf` + `otu<N>_fdr<X>_prev<Y>_dataset_summary.csv`

## Scripts
- `1_correlation_signature_matrix_intermediate_processing.R` — build signature + correlation intermediates from svyglm outputs
- `2_microbial_signature_heatmap_sized.R` — render legacy sized heatmaps from intermediates
- `2_microbial_signature_heatmap_sized_cutree_optimal_k_sigres.R` — render silhouette-driven optimal-k heatmaps

## Environment
R >= 4.5 with: `data.table`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`, `DBI`, `RSQLite`, `ggplot2`, `ComplexHeatmap`, `circlize`, `RColorBrewer`, `grid`, `cluster`, `extrafont` (Arial).

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open each `.R` file and set `PROJECT_ROOT` (top) to your local clone of this repository.
2. From the repo root:
   ```bash
   Rscript scripts/7_microbial_signature_heatmap/1_correlation_signature_matrix_intermediate_processing.R
   Rscript scripts/7_microbial_signature_heatmap/2_microbial_signature_heatmap_sized.R
   Rscript scripts/7_microbial_signature_heatmap/2_microbial_signature_heatmap_sized_cutree_optimal_k_sigres.R \
     otu_threshold=2 fdr_threshold=0.05 prevalence_threshold=0.10
   ```
   The third script accepts CLI flags: `otu_threshold`, `fdr_threshold`, `prevalence_threshold`, `k_min`, `k_max`.

## Run order
Script 1 must run first — it builds the intermediates that both plotting scripts read. The two plotting scripts can then run in either order, independently of each other.
