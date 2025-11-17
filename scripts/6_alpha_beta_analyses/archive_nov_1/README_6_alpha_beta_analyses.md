# Alpha and Beta Diversity Analyses

This directory contains scripts for comprehensive alpha and beta diversity analyses of the NHANES oral microbiome dataset.

## Overview

The analysis pipeline consists of three main steps:

1. **Data Loading and Preparation** - Load diversity metrics and categorical variables
2. **Integrated Diversity Analysis** - Generate 4×3 grid plots with statistical tests
3. **Supplementary Table Creation** - Generate categorical variable definitions table

## Files

### Core Analysis Scripts

- **`1_load_diversity_and_all_categories_local.R`** - Data loading script for local execution
- **`1_load_diversity_and_all_categories_HPC.R`** - Data loading script for HPC/SLURM clusters
- **`2_integrated_diversity_FINAL_local.R`** - Integrated analysis script for local execution
- **`2_integrated_diversity_FINAL_HPC.R`** - Integrated analysis script for HPC/SLURM clusters

### Documentation and Utilities

- **`0_data_description_method_result_explanation.md`** - Comprehensive documentation of methods and results
- **`3_create_categorical_variable_supplementary_table.R`** - Generate supplementary table of categorical variables
- **`README_6_alpha_beta_analyses.md`** - This file

### Data Files

- **`data/`** - Directory containing processed data files
- **`supplementary_table_categorical_variable_definitions.csv`** - Categorical variable definitions
- **`supplementary_table_categorical_variable_definitions.md`** - Markdown version of definitions

## Usage

### Local Execution

For local development and testing:

```bash
# Step 1: Load data
Rscript 1_load_diversity_and_all_categories_local.R

# Step 2: Run integrated analysis
Rscript 2_integrated_diversity_FINAL_local.R
```

### HPC/SLURM Execution

For high-performance computing clusters:

```bash
# Step 1: Load data
Rscript 1_load_diversity_and_all_categories_HPC.R

# Step 2: Run integrated analysis
Rscript 2_integrated_diversity_FINAL_HPC.R
```

## Analysis Modes

The integrated diversity analysis supports three modes:

1. **`no_tick_and_bracket_asterisk`** - Remove x-axis labels and show asterisks for significance
2. **`no_tick_and_bracket`** - Remove x-axis labels but show exact p-values
3. **Normal mode** - Standard plots with full labels and p-values

## Output

### Generated Files

- **PDF plots** - 4×3 grid plots for each categorical variable (120+ files total)
- **Statistical results** - CSV files with pairwise comparisons
- **Data files** - Processed datasets in RDS format

### Plot Structure

Each plot contains:
- **Row 1**: Alpha diversity (Observed OTUs, Shannon, Inverse Simpson)
- **Row 2**: Beta diversity centroids (Bray-Curtis, Unweighted UniFrac, Weighted UniFrac)
- **Row 3**: PCoA plots with confidence ellipses
- **Row 4**: Scree plots showing variance explained

## Requirements

### Local Environment
- R 4.3+
- Required packages: dplyr, tidyr, ggplot2, phyloseq, vegan, ape, egg, cowplot

### HPC Environment
- R 4.3+
- Conda environment: `o2-nhanes_oral_env`
- Required packages: dplyr, tidyr, ggplot2, ape, egg, cowplot (vegan and phyloseq not required)

## Data Requirements

- Alpha diversity metrics from `dada2rsv-alpha.txt`
- Beta diversity matrices from `dada2rsv-beta-*.txt`
- Categorical variables from database tables
- Phyloseq object (local only)

## Memory Requirements

- **Local**: 8-16 GB RAM recommended
- **HPC**: 128 GB RAM available for large-scale analysis

## Troubleshooting

### Common Issues

1. **Package dependencies** - Use HPC versions for cluster environments
2. **Memory limitations** - Ensure sufficient RAM for 9,349 samples
3. **File paths** - Verify correct paths for your environment

### Support

For technical issues, refer to:
- `0_data_description_method_result_explanation.md` for detailed methodology
- Check log files for specific error messages
- Verify data file availability and permissions

## Results

The analysis generates comprehensive diversity plots and statistical results for 40+ categorical variables including demographics, oral health, and disease outcomes. All outputs are saved in the `results/analyses_results/6_alpha_beta_analyses_out/` directory.
