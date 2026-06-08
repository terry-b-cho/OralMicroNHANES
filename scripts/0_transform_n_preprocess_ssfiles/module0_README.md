# Module 0 — Microbiome transformation and schema-structure preprocessing

## What this module does
Three-step preprocessing pipeline:
1. Transform NHANES oral microbiome abundance tables into four normalized versions (`none`, `hellinger`, `clr`, `lognorm`).
2. Complete the transformed SQLite (derive `RIAGENDR_01`, fill missing `variable_names_epcf` rows for derived variables).
3. Build 24 schema-structure CSVs (6 WAS analyses × 4 transformations) used downstream by the association pipeline.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite`
- `configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt`

## Outputs (relative to `PROJECT_ROOT`)

> **Not on GitHub.** The transformed SQLite databases are too large to host on GitHub. The full set of regression outputs from the downstream association pipeline is deposited at **https://doi.org/10.5281/zenodo.17871009** and is structured exactly as it would be reproduced by re-running this pipeline.

Expected tree once all three steps have run locally:

```
data/00_nhanes_omp_transformed_db/
├── nhanes_oral_transformed.sqlite          # Step 1 output
└── nhanes_oral_transformed_complete.sqlite # Step 2 output

results/0_ss_files/                         # Step 3 output (20 CSVs)
├── 1_demoWAS_{clr,hellinger,lognorm,none}_schema_structure.csv
├── 2_oradWAS_{clr,hellinger,lognorm,none}_schema_structure.csv
├── 3_exWAS_{clr,hellinger,lognorm,none}_schema_structure.csv
├── 4_pheWAS_{clr,hellinger,lognorm,none}_schema_structure.csv
└── 5_outWAS_{clr,hellinger,lognorm,none}_schema_structure.csv
```

## Scripts
- `nhanes_omp_transformation.R` — Step 1: build the transformed DB.
- `run_transformation.sh` — SLURM wrapper for Step 1.
- `nhanes_db_filling_missing_data.R` — Step 2: complete the DB.
- `ss_file_create.R` — Step 3 worker: build one schema-structure CSV.
- `run_ss_file_create.sh` — SLURM wrapper that submits the 24 Step 3 jobs.

## Environment
- **Step 1** (`run_transformation.sh`): SLURM module `gcc/14.2.0 R/4.4.2`. The wrapper loads these itself.
- **Step 2** (`nhanes_db_filling_missing_data.R`) and **Step 3** (`run_ss_file_create.sh`): project conda env at `envs/.conda/envs/nhanes-analysis` (R ≥ 4.5). The Step 3 wrapper activates this env in each SLURM job.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`. Exact package versions in `module0_tool_version_list.txt`.

## How to run
1. Open each `.sh` and `.R` file and set `PROJECT_ROOT` (top of file) to the absolute path of your local clone of this repository.
2. From the repo root, run the three steps in order:

```bash
# Step 1
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh

# Step 2 (after Step 1 completes)
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate envs/.conda/envs/nhanes-analysis
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
  --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite

# Step 3 (after Step 2 completes)
bash scripts/0_transform_n_preprocess_ssfiles/run_ss_file_create.sh
```

## Run order
Steps must run in order: 1, then 2, then 3. Step 3 reads the complete DB produced by Step 2; Step 2 reads the transformed DB produced by Step 1.
