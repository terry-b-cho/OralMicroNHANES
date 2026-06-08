# Module 2 — Database preprocessing + phyloseq object creation

## What this module does
1. Derives ~50 comprehensive demographic variables from `DEMO_F`/`DEMO_G` in the NHANES SQLite DB (factors, quartiles, dummy variables, alternate-name mappings) and writes them back into a processed DB.
2. Builds six phyloseq objects (counts, relative, relative_none, relative_clr, relative_lognorm, relative_hellinger) from the processed DB.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`

## Outputs
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` *(processed DB with derived demographics)*
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/`
  - `ubiome_counts.rds`, `ubiome_relative.rds`
  - `ubiome_relative_none.rds`, `ubiome_relative_clr.rds`, `ubiome_relative_lognorm.rds`, `ubiome_relative_hellinger.rds`
  - `nhanes_phyloseq_objects_all.rds` *(combined list of all six)*

## Scripts
- `process_db_n_phyloseq.R` — single end-to-end script (demographic processing → phyloseq creation)

## Environment
R >= 4.5 with: `DBI`, `RSQLite`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`, `forcats`, `fs`, `ape`, `data.tree`, `phyloseq`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open `process_db_n_phyloseq.R` and set `PROJECT_ROOT` (top) to your local clone of this repository.
2. From the repo root:

```bash
Rscript scripts/2_preprocess_db_n_phyloseq/process_db_n_phyloseq.R
```

## Run order
Single script, top-to-bottom. No prior module is required from within this repo.
