# Conda environment

All R scripts in this repository are run under a single conda environment specified by `nhanes-analysis_for_reviewers.yml` (exported from the live env, R 4.5.1, 200+ pinned packages).

## Create

```bash
conda env create -f envs/nhanes-analysis_for_reviewers.yml
```

This creates an env named `nhanes-analysis-for-reviewers`.

## Activate

```bash
conda activate nhanes-analysis-for-reviewers
```

Each script in this repo expects this env to be active. On HPC (Slurm / O2):

```bash
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate nhanes-analysis-for-reviewers
```
