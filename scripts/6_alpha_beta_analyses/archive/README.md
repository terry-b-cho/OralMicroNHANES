# Alpha and Beta Diversity Analysis - NHANES Oral Microbiome

## Overview

This directory contains comprehensive end-to-end diversity analyses using pre-computed alpha and beta diversity metrics from the NHANES oral microbiome study.

## Data Sources

### Pre-computed Diversity Metrics
- **Alpha Diversity**: `/data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt` (9,662 samples × 201 metrics)
- **Beta Diversity**:
  - Bray-Curtis: `/data/00_nhanes_omp_diversity_db/dada2rsv-beta-braycurtis.txt` (9,349 × 9,349)
  - Unweighted UniFrac: `/data/00_nhanes_omp_diversity_db/dada2rsv-beta-unwunifrac.txt` (9,349 × 9,349)
  - Weighted UniFrac: `/data/00_nhanes_omp_diversity_db/dada2rsv-beta-wunifrac.txt` (9,349 × 9,349)

### Metadata Sources
- **Database**: `/data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite`
- **Phyloseq Objects**: `/results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/`

## Analysis Pipeline

### 1. Data Loading (`1_load_diversity_data.R`)
**Purpose**: Load and merge pre-computed diversity metrics with categorical variables

**Key Features**:
- Extracts alpha diversity metrics at rarefaction depth 10,000
- Loads demographic data from SQLite database
- Creates age groups (14-19, 20-24, ..., 80-85)
- Ensures proper SEQN matching (character type)
- Uses subset of 1,000 samples for testing

**Outputs**:
- `data/alpha_diversity_with_categories.rds` - Merged alpha diversity + demographics
- `data/demographic_data.rds` - Demographic variables
- `data/beta_diversity_matrices.rds` - Beta diversity distance matrices

**SEQN Matching Logic**:
```r
# Alpha diversity file: SEQN is numeric → convert to character
alpha_diversity$SEQN <- as.character(SEQN)

# Database: SEQN needs character conversion
demo_data$SEQN <- as.character(SEQN)

# Merge uses character SEQN
alpha_with_categories <- left_join(alpha_diversity, demo_data, by = "SEQN")
```

### 2. Alpha Diversity Analysis (`2_alpha_diversity_analysis.R`)
**Purpose**: Comprehensive alpha diversity visualization and statistical testing

**Visualizations**:
1. **Violin plots by Age Group** (3 metrics)
   - Observed OTUs
   - Shannon Diversity
   - Inverse Simpson
   
2. **Violin plots by Gender** (3 metrics)

3. **Faceted plot** - All metrics side by side

**Plot Specifications**:
- Theme: `egg::theme_article(base_size = 5)`
- Font size: 5pt for all text
- Colors: Grafify palette
- Dimensions: 8×6 (individual), 12×4 (faceted)

**Statistical Tests**:
- **Kruskal-Wallis**: For age group comparisons (non-parametric)
- **Wilcoxon**: For gender comparisons (2 groups)
- **FDR Correction**: Benjamini-Hochberg method

**Outputs**:
- PDFs: `results/alpha_diversity_plots/*.pdf` (7 plots)
- CSV: `results/alpha_diversity_statistical_tests.csv`
- CSV: `results/alpha_diversity_summary_by_age.csv`
- CSV: `results/alpha_diversity_summary_by_gender.csv`

**Results (Test Subset, N=995)**:
- Observed OTUs: Significant by Age Group (P < 0.001) and Gender (P < 0.001)
- Shannon Diversity: Significant by Gender (P < 0.001), not by Age (P = 0.21)
- Inverse Simpson: Significant by Gender (P < 0.001), not by Age (P = 0.73)

### 3. Beta Diversity PCoA (PENDING)
**Purpose**: Principal Coordinate Analysis with marginal distributions

**Planned Features**:
- PCoA for all three beta diversity metrics
- Marginal density plots on both axes
- Colored by categorical variables (Age Group, Gender)
- Faceted comparisons
- Confidence ellipses for groups

### 4. Statistical Tests (PENDING)
**Purpose**: PERMANOVA and dispersion tests

**Planned Tests**:
- PERMANOVA (adonis2) for all beta diversity metrics
- Betadisper for within-group dispersion
- Pairwise comparisons
- Effect sizes

## Color Palettes

### Phylum Colors
```r
phylum_colors <- c(
  "Firmicutes"           = "#E69F00",  # orange
  "Bacteroidetes"        = "#56B4E9",  # blue
  "Actinobacteria"       = "#009E73",  # green
  "Proteobacteria"       = "#F0E442",  # yellow
  "Fusobacteria"         = "#CC6677",  # dark blue
  "Spirochaetae"         = "#D55E00",  # reddish-orange
  "Cyanobacteria"        = "#CC79A7",  # purple
  "Acidobacteria"        = "#999999",  # grey
  "Candidate division SR1" = "#AD7700", # brown
  "Planctomycetes"       = "#332288",  # deep indigo
  "Saccharibacteria"     = "#44AA99",  # teal
  "Synergistetes"        = "#88CCEE",  # light blue
  "Tenericutes"          = "#117733",  # dark green
  "unclassified"         = "#DDDDDD",  # very light grey
  "Verrucomicrobia"      = "#0072B2"   # reddish pink
)
```

### Grafify Colors (for categorical variables)
52 colors from Grafify palette for automatic coloring

## Environment Setup

### Conda Environment
```bash
mamba env create -f envs/oral_env.yaml
mamba activate oral_env
```

### Required Packages
- R base ≥4.2
- phyloseq, ggplot2, dplyr, tidyr, readr
- vegan, ape (for diversity analyses)
- cowplot, ggbeeswarm, egg (for visualization)
- DBI, RSQLite (for database access)

## Usage

### Run Complete Pipeline
```bash
cd /Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho

# Step 1: Load data
mamba run -n oral_env Rscript scripts/6_alpha_beta_analyses/1_load_diversity_data.R

# Step 2: Alpha diversity analysis
mamba run -n oral_env Rscript scripts/6_alpha_beta_analyses/2_alpha_diversity_analysis.R

# Step 3: Beta diversity PCoA (when ready)
# mamba run -n oral_env Rscript scripts/6_alpha_beta_analyses/3_beta_diversity_pcoa.R

# Step 4: Statistical tests (when ready)
# mamba run -n oral_env Rscript scripts/6_alpha_beta_analyses/4_statistical_tests.R
```

### Outputs Directory Structure
```
scripts/6_alpha_beta_analyses/
├── data/                                    # Processed data
│   ├── alpha_diversity_with_categories.rds
│   ├── beta_diversity_matrices.rds
│   └── demographic_data.rds
├── results/                                 # Analysis results
│   ├── alpha_diversity_statistical_tests.csv
│   ├── alpha_diversity_summary_by_age.csv
│   ├── alpha_diversity_summary_by_gender.csv
│   └── alpha_diversity_plots/               # Visualization PDFs
│       ├── alpha_observed_otus_by_age_group.pdf
│       ├── alpha_shannon_diversity_by_age_group.pdf
│       ├── alpha_inverse_simpson_by_age_group.pdf
│       ├── alpha_observed_otus_by_gender.pdf
│       ├── alpha_shannon_diversity_by_gender.pdf
│       ├── alpha_inverse_simpson_by_gender.pdf
│       └── alpha_faceted_by_age_group.pdf
└── README.md                                # This file
```

## Key Design Decisions

### 1. Pre-computed Metrics
**Why**: Original diversity calculations took hours and require phylogenetic trees
- Use rarefaction depth 10,000 (most complete coverage)
- Use first resampling (RSV_*_10000_0 columns)
- Bray-Curtis, UniFrac metrics pre-computed in QIIME2/DADA2

### 2. SEQN Handling
**Critical**: SEQN type must be consistent across all datasets
- Alpha/Beta files: SEQN is numeric → convert to character
- Database: SEQN is integer → convert to character  
- Phyloseq rownames: Already character
- **Always use character for merging**

### 3. Age Grouping
**Following 5_age_analysis.Rmd**:
- Breaks: `c(14, seq(20, 85, by = 5))`
- Labels: `"14-19"`, `"20-24"`, ..., `"80-85"`
- Reference group: `"30-39"` (median age)

### 4. Testing Strategy
**Large datasets require careful handling**:
- Full alpha diversity: 9,662 samples × 201 columns
- Full beta diversity: 9,349 × 9,349 matrices (~574MB-991MB each)
- **Solution**: Use 1,000 sample subset for testing
- **Production**: Remove subset limits in scripts

## Troubleshooting

### Error: "Can't join due to incompatible types"
**Solution**: Ensure SEQN is character in both datasets before joining

### Error: "RSQLite package not found"
**Solution**: `mamba install -n oral_env r-rsqlite r-dbi`

### Memory issues with beta diversity
**Solution**: 
- Use subset for testing
- Process one metric at a time
- Consider using `bigstatsr` or `bigmemory` for production

### Plots not displaying properly
**Solution**: Ensure `egg` package is installed and theme is applied with `base_size = 5`

## Next Steps

1. ✅ Load pre-computed diversity data
2. ✅ Alpha diversity visualization & tests
3. ⏳ Beta diversity PCoA with marginal distributions
4. ⏳ PERMANOVA and statistical tests
5. ⏳ Combined visualization dashboard

---

**Last Updated**: October 8, 2025  
**Status**: Alpha diversity analysis complete, Beta diversity in progress  
**Contact**: Terry Cho

