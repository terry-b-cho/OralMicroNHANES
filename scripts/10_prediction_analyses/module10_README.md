# Module 10 — Prediction from oral microbiome features

## What this module does
End-to-end demonstration that the NHANES oral microbiome carries predictive signal for participant-level demographics. Three canonical prediction tasks on the same shared sample universe and the same feature backbone (prevalence-filtered CLR-transformed genus abundances):

| Task | Type | Target |
|---|---|---|
| `age` | regression | continuous age in years |
| `age_group` | multi-class classification | 14-19, 20-29, 30-39, 40-49, 50-59, 60-69 |
| `gender` | binary classification | Male vs Female |

Each task is fit twice — **microbiome only** (CLR features) and **microbiome + basic demographics** (excluding the target itself) — with `glmnet` elastic-net regularized regression (α = 0.5, 5-fold CV). Test-set metrics are reported both unweighted and survey-weighted (WTMEC2YR).

## Inputs (relative to `PROJECT_ROOT`)
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_clr.rds` *(CLR feature values)*
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_counts.rds` *(used only to compute prevalence — taxa kept if non-zero in ≥1% of samples)*
- `data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt`
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` *(reads `DEMO_F` + `DEMO_G` only)*

## Outputs (under `results/analyses_results/10_prediction_analyses_out/`)
- `intermediate/` — `clr_features.rds`, `alpha_features.rds`, `covariates_and_targets.rds`, `sample_universe.rds`
- `models/<task>__<feature_set>.rds` — fitted `cv.glmnet` object per task × feature set
- `predictions/<task>__<feature_set>.rds` — per-test-sample y_true + y_pred (+ survey weight)
- `metrics_summary.csv` — tidy one-row-per-(task, feature_set, weighting) results table
- `figures/` — publication-grade visualizations (7 PDFs, all faceted across tasks)
  - `metrics_overview.pdf` — per-task tabular heat-map: rows = (feature_set × weighting), columns = the metrics relevant to that task type, cells coloured by per-metric rank (deeper = better)
  - `roc_curves_classification.pdf` — ROC for every classification task (binary `gender`; multinomial `age_group` one-vs-rest); colour = class (grafify palette), linetype = feature set
  - `precision_recall_curves.pdf` — PR curves, same faceting; dotted horizontal = class-prevalence baseline (preferred over ROC when classes are imbalanced)
  - `calibration_curves.pdf` — decile reliability curves: mean predicted probability vs observed positive frequency per bin; point size = bin n
  - `confusion_matrices.pdf` — row-normalised confusion matrix per (task × feature_set) panel
  - `feature_importance.pdf` — top-20 features per (task × feature_set) by |glmnet coefficient at λ.min| (max over classes for multinomial)
  - `regression_diagnostics.pdf` — 3-panel composite for each regression task: (a) observed vs predicted, (b) residuals vs predicted with loess bias trend, (c) residual Q-Q vs normal

## Scripts
- `01_prepare_features_and_targets.R` — single data-prep script (features, covariates, survey design)
- `02_run_prediction_tasks.R` — single train/eval script (3 tasks × 2 feature sets × 1 model)
- `03_visualize_results.R` — single visualization script (reads metrics + predictions, writes PDFs)

## Environment
R >= 4.5 with: `phyloseq`, `DBI`, `RSQLite`, `dplyr`, `tidyr`, `tibble`, `data.table`, `readr`, `stringr`, `glmnet`, `pROC`, `ggplot2`, `patchwork`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open each `.R` file and set `PROJECT_ROOT` (top) to your local clone of this repository.
2. From the repo root:

```bash
Rscript scripts/10_prediction_analyses/01_prepare_features_and_targets.R
Rscript scripts/10_prediction_analyses/02_run_prediction_tasks.R
Rscript scripts/10_prediction_analyses/03_visualize_results.R
```

## Run order
`1 → 2 → 3`. Script 1 builds intermediates, script 2 fits models and writes metrics + predictions, script 3 renders performance visualizations.
