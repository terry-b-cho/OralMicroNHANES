# Module 11 — Per-cycle re-estimation and cross-cycle validation of the pooled CLR associations

## What this module does
Re-runs every pooled (NHANES 2009–2012) survey-weighted CLR association regression
independently within cycle F (2009–10) and cycle G (2011–12) using the identical engine,
covariates, and coding (single-cycle NCHS designs: unscaled `WTMEC2YR`, `SDMVPSU`/`SDMVSTRA`;
per-cycle p-values on design-df t, df = 16/17), then validates every pooled BH-FDR<0.05
association across the two independent cycles (direction concordance; replication-FDR).
Produces the Figure S16 main figures, accessory figures, core validation tables, the
Supplementary Table S11 workbook, and a self-contained HTML report.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`
- `results/0_ss_files/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_clr_schema_structure.csv`
- `results/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_out/result_clr/*_clr_tidied_complete.rds`
- `results/analyses_results/10_prediction_analyses_out/intermediate/qc/qc_microbiome_prevalence_summary.csv`
- `results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv`

## Outputs (under `results/analyses_results/11_per_year_association/`)
- `schema_per_year/{F,G}/` — single-cycle schema CSVs
- `per_year_{F,G}/<analysis>_out/result_clr/` — per-cycle regression results + aggregates
- `comparisons/` — harmonized master table, F-vs-G and pooled-vs-cycle comparisons, availability
- `tables/` — core tables in workbook order: `table1_cycle_characteristics`,
  `table2_validation_summary`, `table3_replication_tiers`, `table4_validated_core`
- `figures/` — main figures `figS16A_left_*`, `figS16A_right_*`, `figS16B_*`, `figS16C_*`;
  all other figures + accessory tables under `figures/additional_accessory_figures/` (`optional_*`)
- `Table_S11_summary_of_cross_cycle_validation_data_and_results.xlsx` — 7-sheet supplementary workbook
- `reports/module11_per_cycle_validation_report.{md,html}` — self-contained report (figures embedded)
- `logs/per_year/` — SLURM logs

## Scripts
- `00_module11_config.R` — shared config: paths, thresholds, palettes, figure/table routing
- `01_prepare_per_year_schemas.sh` — split the pooled CLR schemas by cycle
- `02_submit_per_year_was_clr.sh` — SLURM submitter, one job per (cycle × dependent variable); resume-safe
- `03_aggregate_per_year_results.sh` — per (cycle × analysis) aggregation, BH-FDR + Storey q
- `04_variable_availability.R` — cross-cycle variable availability; table1 + accessory availability table
- `05_harmonize_and_compare.R` — harmonized master table, validation/replication statistics, core tables 2–4
- `06_figures.R` — all figures (main Figure S16 + accessory)
- `07_report.R` — markdown + self-contained HTML report
- `09_supplementary_excel.R` + `09_make_xlsx.py` — Supplementary Table S11 workbook
- `universal_was_analysis.R`, `aggregate_was_results.R` — regression engine / aggregator (copies of the pooled pipeline's)

## Environment
R >= 4.5 with: `survey`, `qvalue`, `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `purrr`,
`broom`, `DBI`, `RSQLite`, `optparse`, `getopt`, `ggplot2`, `egg`, `extrafont` (Arial),
`jsonlite`, `base64enc`. External tools: `pdftoppm` (poppler), `pandoc`, SLURM (step 02),
Python 3.9 + `openpyxl` (step 09).

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

One-time Python venv for the workbook step:
```bash
python3 -m venv scripts/11_per_year_association/.pyxlsx_venv
scripts/11_per_year_association/.pyxlsx_venv/bin/pip install openpyxl==3.1.5
```

## How to run
1. `export PROJECT_ROOT=/path/to/your/clone` (or simply run everything from the repo root).
2. From the repo root:

```bash
./scripts/11_per_year_association/01_prepare_per_year_schemas.sh
./scripts/11_per_year_association/02_submit_per_year_was_clr.sh    # SLURM; wait for jobs to drain
./scripts/11_per_year_association/03_aggregate_per_year_results.sh
Rscript scripts/11_per_year_association/04_variable_availability.R
Rscript scripts/11_per_year_association/05_harmonize_and_compare.R
Rscript scripts/11_per_year_association/06_figures.R
Rscript scripts/11_per_year_association/09_supplementary_excel.R
Rscript scripts/11_per_year_association/07_report.R
```

## Run order
Strictly 01 → 02 → 03 → 04 → 05 → 06 → 09 → 07. Steps 01–03 are resume-safe and never
overwrite; steps 04–09 are deterministic (seed 42) and regenerate their outputs in place.
Requires the pooled results (`results/<analysis>_out/`) and the module 03/10 intermediates
listed above to exist first.
