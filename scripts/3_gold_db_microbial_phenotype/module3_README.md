# Module 3 — GOLD database microbial phenotype mapping

## What this module does
Aggregates the GOLD microbial phenotype database by genus, parses NHANES oral microbiome genus names, and produces the genus-to-GOLD-feature mapping used downstream.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_GOLDdb/goldData_ubiome_selected.csv`
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite`
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative.rds`

## Outputs (relative to `PROJECT_ROOT`)
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/gold_db_genus.csv`
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv`
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/mapping_summary_stats.csv`

## Scripts
- `gold_db_process_n_map2genus.R` — GOLD aggregation + genus name parsing + merge

## Environment
R >= 4.5 with: `data.table`, `dplyr`, `dbplyr`, `RSQLite`, `readr`, `tidyr`, `tibble`, `purrr`, `phyloseq`. Exact versions in `module3_tool_version_list.txt`.

## How to run
1. Open `gold_db_process_n_map2genus.R` and set `PROJECT_ROOT` (at the top) to the absolute path of your local clone of this repository.
2. From the repo root:

```bash
Rscript scripts/3_gold_db_microbial_phenotype/gold_db_process_n_map2genus.R
```

## Run order
Single script, no order required.
