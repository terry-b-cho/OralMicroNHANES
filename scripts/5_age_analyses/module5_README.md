# Module 5 — Age-stratified microbiome analyses

## What this module does
Clusters genus-level oral microbiome relative abundances across age groups (Ward.D2 on z-scored mean abundance) and produces clustering diagnostics, taxonomic visualizations, ecology-by-cluster summaries, and a nativity composition check.

## Inputs (relative to `PROJECT_ROOT`)
- `results/intermediate/ubiome_relative_none_updated.rds`
- `data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite` *(used by `5_age_analyses_full.R`)*
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv` *(used by `5_age_analyses_additional.R`)*

## Outputs (relative to `PROJECT_ROOT`)
- `results/analyses_results/5_age_analyses_out/` — outputs from `5_age_analyses_full.R` (PDFs + CSVs)
- `results/analyses_results/5_age_analyses_out_additional/` — outputs from `5_age_analyses_additional.R`:
  - `k_optimization/plots/k_optimization_panels.pdf`, `consensus_stability.pdf`
  - `k_optimization/metrics/{min_size,wk,silhouette,gap,stability,k_decisions}_table.csv`
  - `k_optimization/k_optimization_summary.md`
  - `ecology/` (trait distributions, ternary summaries, manifest)
  - `nativity/nativity_by_age_group_with_counts.pdf` (+ plot-input CSV)
  - `genus_cluster_assignments.csv`
- `results/analyses_results/5_age_analyses_out_additional_k_choices/` — outputs from `5_age_analyses_additional_k_choices.R`: all of the above plus per-k parallel-coordinate (`5.3_*`) and PCA+t-SNE (`5.4_*`) PDFs for k ∈ {3, 7, 9, 12}. The k = 9 PDFs reproduce `5_age_analyses_full.R`'s 5.3 / 5.4. k = 3, 7, 12 use a max-contrast Kelly cluster palette with legends split into companion `5.4_combined_pca_tsne_k{3,7,12}_legend.pdf` files.

## Scripts
- `5_age_analyses_full.R` — clustering pipeline (silhouette + gap), parallel coordinates, PCA/t-SNE, heatmap, integrated centroid plots
- `5_age_analyses_additional.R` — corrected silhouette + gap + bootstrap stability diagnostics, GOLD-trait ecology by cluster, nativity composition by age group
- `5_age_analyses_additional_k_choices.R` — superset of `5_age_analyses_additional.R` that additionally renders 5.3 + 5.4 figures at k = 3, 7, 9, 12. Uses a static t-SNE coordinate fixture (`results/intermediate/historical_tsne_coords_k9.csv`, extracted from the historical 5.4 via `tools/extract_historical_tsne.py`) so the t-SNE pattern is byte-identical across runs.

## Environment
R >= 4.5 with: `phyloseq`, `dplyr`, `tidyr`, `cluster`, `ggplot2`, `egg`, `patchwork`, `scales`, `grid`, `gridExtra`, `data.table`, `stringr`, `GGally`, `Rtsne`, `pheatmap`, `cowplot`, `DBI`, `RSQLite`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open each `.R` file and set `PROJECT_ROOT` (at the top) to the absolute path of your local clone of this repository.
2. From the repo root:

```bash
Rscript scripts/5_age_analyses/5_age_analyses_full.R
Rscript scripts/5_age_analyses/5_age_analyses_additional.R
Rscript scripts/5_age_analyses/5_age_analyses_additional_k_choices.R   # optional: per-k 5.3 + 5.4 at k = 3, 7, 9, 12
```

## Run order
The three scripts are independent and can be run in any order. `5_age_analyses_additional.R` produces the corrected silhouette diagnostics that supersede the silhouette block flagged in `5_age_analyses_full.R`. `5_age_analyses_additional_k_choices.R` is a self-contained superset that reproduces the k-optimization / ecology / nativity outputs plus per-k 5.3 + 5.4 figures.
