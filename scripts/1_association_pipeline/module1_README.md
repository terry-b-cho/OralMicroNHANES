# Module 1 — WAS association pipeline (24 schemes)

## What this module does
Per-dependent-variable, survey-weighted regression (NCHS pooled-cycle 2009-2012) for **24 analysis schemes** = 6 WAS types × 4 transformations:

| WAS type | Description | Model |
|---|---|---|
| `1_demoWAS` | Demographics → microbiome | linear |
| `2_oradWAS` | Microbiome → oral health | logistic |
| `3_exWAS` | Exposures → microbiome | linear |
| `4_pheWAS` | Microbiome → phenotypes | linear |
| `5_outWAS` | Microbiome → disease outcomes | logistic |

Transformations: `none` (raw proportions), `hellinger` (√P), `clr` (centred log-ratio), `lognorm` (log₁₀-CPM).

Aggregates the per-dep-var results, applies scheme-wise Benjamini-Hochberg FDR (and optional Storey q-value) on main-effect terms only, and writes combined `*_tidied_complete.rds` files plus 4 supplementary CSV tables.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`
- `results/0_ss_files/<type>_<norm>_schema_structure.csv` — 20 schema files defining the dep-var × indep-var pairs to run

## Outputs

> **The full regression outputs are NOT bundled in this GitHub repo because of their size.** All regression outputs (≈ 11,400 per-dep-var RDS files + aggregated `*_complete.rds` + supplementary CSV tables) are deposited as a single archive on **Zenodo** at **<https://doi.org/10.5281/zenodo.17871009>**. Unpack the archive into `results/` to obtain the layout below, which is also exactly what re-running this pipeline reproduces.

```
results/
├── <type>_out/                       # one per WAS type (1_demoWAS … 5_outWAS)
│   └── result_<norm>/                # one per transformation (none / hellinger / clr / lognorm)
│       ├── <dep_var>.rds             # per-dep-var output: list(pe_tidied, pe_glanced, rsq)
│       ├── <type>_<norm>_tidied_complete.rds      # aggregated coefficient table + FDR columns
│       ├── <type>_<norm>_glanced_complete.rds     # aggregated model-level stats
│       ├── <type>_<norm>_rsq_complete.rds         # aggregated R² values
│       └── <type>_<norm>_aggregation_summary.txt
└── supplementary_tables/
    ├── s_table_aggregation_summary_overall.csv
    ├── s_table_aggregation_summary_significance.csv
    ├── s_table_aggregation_summary_detailed_statistics.csv
    └── s_table_scheme_wise_correction_summary.csv
```

## Scripts
- `universal_was_analysis.R` — core analysis engine (1 dep-var per invocation)
- `run_all_was_analyses_flexible.sh` — SLURM batch submitter; resumes safely (never resubmits a dep-var that already has its RDS)
- `aggregate_was_results.R` — aggregates + FDR-corrects + writes supplementary tables

## Environment
R >= 4.5 with: `optparse`, `DBI`, `RSQLite`, `dplyr`, `readr`, `purrr`, `tibble`, `stringr`, `survey`, `broom`, `getopt`. Optional: `qvalue` (Bioconductor) — Storey q-values; falls back to BH-only if absent.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run

Open `run_all_was_analyses_flexible.sh` and set `PROJECT_ROOT` (top) to your local clone of this repository. From an O2 login or interactive node:

```bash
module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate $PROJECT_ROOT/envs/.conda/envs/nhanes-analysis
cd $PROJECT_ROOT

# 1) Submit all 24 schemes (test mode: only first dep-var per scheme, 15 min)
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh test

# 2) Submit all 24 schemes (production)
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh

# 3) After all per-dep-var jobs finish, aggregate + supplementary tables
Rscript scripts/1_association_pipeline/aggregate_was_results.R --run-all $PROJECT_ROOT/results TRUE
```

## Run order
`run_all_was_analyses_flexible.sh` → wait for all SLURM jobs to finish → `aggregate_was_results.R --run-all`.
