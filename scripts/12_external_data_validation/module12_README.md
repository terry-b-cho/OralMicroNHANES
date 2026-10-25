# Module 12 — External validation against independent oral microbiome cohorts

## What this module does
Tests whether the NHANES genus-level associations reproduce in independent oral microbiome
cohorts: harmonizes each external cohort's metadata and genus tables, refits the matching
association models, and reports direction agreement, rank agreement, and coverage per
phenotype and cohort (including comparisons that could not be made, and why). Produces the
Supplementary Figure S15 panels and the Supplementary Table S10 workbook. The prepared
external data are committed under
`results/analyses_results/12_external_data_validation/datasets/` — nothing needs to be
downloaded.

## Inputs (relative to `PROJECT_ROOT`, from outside this module)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite`
- `results/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_out/result_clr/*_clr_tidied_complete.rds`
- `results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv`
- `results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/nhanes_phyloseq_objects_all.rds`
- `results/analyses_results/6_alpha_beta_analyses_out/alpha_pairwise_stat_res.csv`

## Outputs (under `results/analyses_results/12_external_data_validation/`)
- `nhanes_genus_reference.tsv`, `genus_overlap_summary.tsv`, `genus_presence_matrix.tsv`,
  `dataset_catalog.tsv`, `phenotype_availability_matrix.tsv`,
  `replication_targets_{all,rows}.tsv`, `pairwise_direction_concordance.tsv`
- `validation/figures/{pdf,png}/` — the Supplementary Figure S15 panels
- validation + audit tables (per-phenotype concordance, opportunity matrix, master validation table)
- `S10_summary_of_external_concordance_validation_data_and_results.xlsx` — Supplementary Table S10
- `datasets/<dataset_id>/` — per-cohort manifests, harmonized metadata, genus tables

## Scripts
- `00_shared/` — `config.sh` (dataset roster), `harmonize_helpers.R` (genus parsing/mapping)
- `config/<dataset_id>.yaml` — one per cohort: source, harmonization rules, disposition
- `04_metadata_harmonization/`, `05_microbiome_processing/`, `08_validation_ready/` — per-cohort
  harmonization to an analysis-ready table
- `09_reports/` — four module-level reference builders (genus reference, overlap matrix, catalog, shortlist)
- `10_validation_analysis/` — numbered steps + `run_all.sh`: extract NHANES effects, refit external
  cohorts, assemble comparisons, draw the S15 panels, write tables
- `11_validation_audit/` — numbered steps + `run_all.sh`: inventories, opportunity matrix, new
  validations, master table, audits, phylum-level and alpha-diversity validation, the S10 workbook
- `notebooks/` — self-contained HTML records of data acquisition and quality-control decisions

## Environment
R >= 4.5 with: `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `purrr`, `broom`, `DBI`,
`RSQLite`, `yaml`, `phyloseq`, `vegan`, `ggplot2`, `ggrepel`, `patchwork`, `grid`,
`markdown`, `base64enc`. Python 3 (standard library only) for the workbook assembler.
Quarto + Pandoc only if re-rendering `notebooks/`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. `export ORALMICRO_ROOT=/path/to/your/clone` (or simply run everything from the repo root).
2. From the repo root:

```bash
S=scripts/12_external_data_validation

# Module-level reference tables, in this order (each depends on the previous)
Rscript $S/09_reports/build_nhanes_genus_reference.R
Rscript $S/09_reports/build_genus_overlap_matrix.R
Rscript $S/09_reports/build_catalog.R
Rscript $S/09_reports/build_shortlist.R

bash $S/10_validation_analysis/run_all.sh   # S15 panels + validation tables
bash $S/11_validation_audit/run_all.sh      # audit tables + Supplementary Table S10
```

## Run order
The four `09_reports/` builders first (in the order shown), then `10_validation_analysis`,
then `11_validation_audit`. Both drivers are safe to re-run into an emptied output
directory; a missing or malformed input trips an assertion rather than passing silently.
Requires the upstream NHANES results listed under Inputs to exist first.
