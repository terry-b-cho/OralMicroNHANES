# Module 5.5 — Smoking-exposure multicollinearity diagnostics

## What this module does
Diagnostics for a 24-variable smoking/combustion exposure set across four measurement modalities (questionnaire smoking status, urinary PAH metabolites, urinary metals & mercapturates, blood VOCs/combustion analytes). Produces host-level covariance/coverage maps, per-modality PCAs, per-variable significant-taxa counts, taxon-by-variable signed-significance heatmap, and modality-level signature similarity / Jaccard summaries.

## Inputs (relative to `PROJECT_ROOT`)
- `results/intermediate/ubiome_variable_mapping.csv`
- `results/intermediate/ubiome_relative_none_updated.rds`
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_clr.rds`
- `results/3_exWAS_out/result_clr/3_exWAS_clr_tidied_complete.rds`
- `results/3_exWAS_out/result_clr/3_exWAS_clr_glanced_complete.rds`
- `results/3_exWAS_out/result_clr/3_exWAS_clr_rsq_complete.rds`
- `results/3_exWAS_out/result_clr/3_exWAS_clr_aggregation_summary.txt`
- `results/analyses_results/7_microbial_signature_heatmap_out/intermediate/microbial_signature/exWAS_clr_all_results/correlation_matrix.rds`

## Outputs (under `results/analyses_results/5.5_smoking_analyses_out_additional/`)
- `plots/a_signature_similarity_heatmap.pdf`
- `plots/b1_per_variable_sig_taxa_bars.pdf`
- `plots/b2_correlation_heatmap_24vars.pdf`
- `plots/c_modality_modality_jaccard.pdf`
- `plots/d_pca_modality_marginal_PCs.pdf`
- `plots/e_taxon_by_variable_significance_heatmap.pdf`
- `tables/modality_stratified_jaccard.csv`
- `tables/pairwise_complete_case_N.csv`
- `inputs/input_provenance.txt`, `inputs/smoking_variable_modality_mapping.csv`, `inputs/tidied_smoking_subset.rds`, `inputs/clr_tax_metadata.rds`
- `manifest.csv`, `run.log`

## Scripts
- `5.5_smoking_analyses_additional.R` — single self-contained runner

## Environment
R >= 4.5 with: `phyloseq`, `dplyr`, `tidyr`, `ggplot2`, `egg`, `patchwork`, `scales`, `grid`, `gridExtra`, `ComplexHeatmap`, `circlize`, `ggrepel`, `viridis`, `RColorBrewer`, `car`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open `5.5_smoking_analyses_additional.R` and set `PROJECT_ROOT` (at the top) to the absolute path of your local clone of this repository.
2. From the repo root:

```bash
Rscript scripts/5.5_smoking_analyses/5.5_smoking_analyses_additional.R
```

## Run order
Standalone. Depends on upstream `results/3_exWAS_out/` and `results/analyses_results/7_microbial_signature_heatmap_out/` outputs existing — generate those first if not already present.
