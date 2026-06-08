# Module 9 — Abundance/prevalence biplot + categorical heatmaps

## What this module does
1. **Biplot** — top-15 taxa per taxonomic rank (Phylum → Genus) with mean relative abundance vs prevalence, phylum-coloured.
2. **Categorical heatmaps** — for each WAS factor (demoWAS / oradWAS / exWAS / pheWAS / outWAS) and each transformation (`none` z-score, `clr` log2 fold-change), test per-OTU differences across factor levels with Kruskal-Wallis (≥3 levels) or Wilcoxon (2 levels), BH-FDR adjust, and render the top 30 ranked taxa as a pheatmap with phylum row annotations.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite`
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_{counts,relative,relative_none,relative_clr}.rds`
- `configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt`
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv`

## Outputs (under `results/analyses_results/9_biplot_n_categories_heatmap_fixed_kw_wilcoxon_out/`)
- `intermediate/` — `phyloseq_objects.rds`, `heatmap_factorized_datasets.rds`, `sample_data_subsets.rds`, `table_description.rds`, `variable_description.rds`, `ubiome_genus_mapping_complete.rds`
- `fig1a_Top_Taxa_Abundance_With_Prevalence_biplot.pdf`
- `categories_heatmap_fixed/none/{<WAS>_<factor>_Abundance_none.pdf, ALL_none_booklet.pdf}` + per-factor CSV tables under `none/supplementary/`
- `categories_heatmap_fixed/clr/{<WAS>_<factor>_log2FC_clr.pdf, ALL_clr_booklet.pdf}` + per-factor CSV tables under `clr/supplementary/`

## Scripts
- `1_load_and_process_all_data_standalone.R` — pull NHANES metadata + phyloseq objects + factor mappings; write intermediates
- `2_biplot_n_categories_heatmap_fixed_kw_wilcoxon.R` — render biplot + categorical heatmaps + CSV tables

## Environment
R >= 4.5 with: `phyloseq`, `DBI`, `RSQLite`, `dplyr`, `tidyr`, `stringr`, `purrr`, `tibble`, `readr`, `forcats`, `data.table`, `ggplot2`, `egg`, `gridExtra`, `grid`, `pheatmap`, `RColorBrewer`, `extrafont`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open each `.R` file and set `PROJECT_ROOT` (top) to your local clone of this repository.
2. From the repo root:

```bash
Rscript scripts/9_biplot_n_categories_heatmap_fixed_kw_wilcoxon/1_load_and_process_all_data_standalone.R
Rscript scripts/9_biplot_n_categories_heatmap_fixed_kw_wilcoxon/2_biplot_n_categories_heatmap_fixed_kw_wilcoxon.R
```

## Run order
`1 → 2`. Script 1 builds all intermediates that script 2 reads from `intermediate/`.
