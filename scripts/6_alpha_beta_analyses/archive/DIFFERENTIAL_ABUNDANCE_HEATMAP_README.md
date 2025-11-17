# Differential Abundance Heatmap Analysis - Data Flow & Usage Guide

## Overview

This document explains how categorical variables are selected and processed for differential abundance heatmap generation, and provides a comprehensive guide to using the standalone plotting script.

## Data Input/Output Locations

### 1. Pre-computed Phyloseq Objects

**Location:** `/results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/`

These are pre-computed and saved phyloseq objects with different transformations:

- `ubiome_relative_none.rds` - Raw relative abundance (proportions summing to 1)
- `ubiome_relative_clr.rds` - Centered log-ratio (CLR) transformed abundances
- `ubiome_relative_lognorm.rds` - Log10-normalized counts

Each phyloseq object contains:
- **OTU table**: Taxon abundances per sample
- **Tax table**: Taxonomic annotations (Phylum, Genus, etc.)
- **Sample data**: All host/clinical variables merged from database

### 2. Variable Configuration Files

**Location:** `/config/`

These text files define which variables belong to each analysis type:

- `1_demoWAS_vars.txt` - Demographic variables (age, gender, education, etc.)
- `2_oradWAS_vars.txt` - Oral disease variables (denture, gum disease, etc.)
- `3_exWAS_vars.txt` - Exposure variables (smoking, hepatitis, HPV, etc.)
- `4_pheWAS_vars.txt` - Phenotype variables (BMI, blood pressure, waist circumference, etc.)
- `5_outWAS_vars.txt` - Outcome variables (asthma, diabetes, heart disease, etc.)
- `6_zimWAS_vars.txt` - Zero-inflated model variables

Format: One variable name per line, lines starting with `#` are comments

### 3. Output Location

**Default:** `/results/analyses_results/04_association_phyloseq_out/visualizations/categories_heatmap_standalone/`

Structure:
```
categories_heatmap_standalone/
├── none/                          # Raw abundance analysis
│   ├── demoWAS_Gender_none.pdf
│   ├── demoWAS_Gender_none_table.csv
│   ├── oradWAS_Denture_none.pdf
│   └── ALL_none_booklet.pdf       # Combined booklet
├── clr/                           # CLR transformed analysis
│   ├── demoWAS_Gender_clr.pdf
│   └── ALL_clr_booklet.pdf
└── lognorm/                       # Log-normalized analysis
    ├── demoWAS_Gender_lognorm.pdf
    └── ALL_lognorm_booklet.pdf
```

## Categorical Variable Selection & Factorization

### How Variables are Selected Per Analysis Type

The pipeline follows this logic:

1. **Load raw variables** from phyloseq sample_data
2. **Subset to relevant variables** using config files
3. **Factorize with appropriate reference levels** based on analysis type
4. **Filter for valid factor levels** (≥2 levels, ≥30 samples per level)

### Reference Level Assignment Strategy

#### demoWAS (Demographics)
- **Gender**: "Female" as reference (alphabetically first, biological reference)
- **Age_group**: "30-39" as reference (median age group)
- **Education_level**: "College/AA" as reference (most common category)
- **Ethnicity**: "White" as reference (largest group)
- **US_born**: "US Born" as reference (majority category)

#### oradWAS (Oral Diseases)
- **All binary variables** use "control" as reference
- **Control definition**: Samples with ALL oral disease variables = 0 (no disease)
- Variables: Denture, Gum_disease, Oral_hygiene, Tooth_decay

#### exWAS (Exposures)
- **Smoking_status**: "Never smoker" as reference (unexposed baseline)
- **Hepatitis/HPV**: "Negative" as reference
- **Categorical exposures**: Lowest/unexposed category as reference

#### pheWAS (Phenotypes)
- **BMI_category**: "Healthy weight" as reference (18.5-25 kg/m²)
- **Blood_pressure**: "Normal" as reference (<120/<80 mmHg)
- **Continuous variables** are binned into categories with clinically meaningful cutpoints

#### outWAS (Outcomes)
- **All binary outcomes** use "control" as reference
- **Control definition**: Samples with ALL outcome disease variables = 0
- Variables: Asthma, Bronchitis, Diabetes, Heart_attack, Stroke, etc.

### Factorization Code Pattern

```r
# Example: Binary disease with control reference
Diabetes = case_when(
  DIABETES == 1 ~ "Diabetes",
  SEQN %in% SEQN_outWAS_control ~ "control",
  TRUE ~ NA_character_
) %>% factor(levels = c("control", "Diabetes"))

# Example: Multi-level categorical with specified reference
Education_level = factor(EducationLevel,
  levels = c("< 9th Grade", "9-11th Grade", "High School", 
             "College/AA", "College Graduate")
) %>% fct_relevel("College/AA")
```

## Statistical Testing & Visualization

### Statistical Tests (per taxon, per categorical variable)

1. **Two-level variables** → Wilcoxon rank-sum test (non-parametric)
2. **Multi-level variables** → Kruskal-Wallis test (non-parametric)
3. **Multiple testing correction** → Benjamini-Hochberg FDR ≤ 0.05

### Post-hoc Pairwise Comparisons (Alpha & Beta Diversity)

For integrated diversity analysis, pairwise comparisons are performed when:
- Overall test (Kruskal-Wallis for alpha, PERMANOVA for beta) is significant (P < 0.05)
- Variable has 2-6 levels (for visual clarity)

**Test Method**: Pairwise Wilcoxon rank-sum tests with FDR correction

**Applied to**:
1. **Alpha diversity metrics**: Observed OTUs, Shannon Diversity, Inverse Simpson
2. **Beta diversity centroid distances**: Bray-Curtis, Unweighted UniFrac, Weighted UniFrac

**Visualization**: 
- Significant pairs (FDR-adjusted P<0.05) shown as brackets above violin plots
- **Exact P-values** annotated with 3 significant figures (e.g., P=0.023, P=0.001)
- Dynamic y-axis buffer adjusts to accommodate all brackets
- Brackets ordered by significance (most significant at bottom)

**Saved Results**:
- `alpha_pairwise_stat_res.csv`: All alpha diversity pairwise comparisons
- `beta_centroid_pairwise_stat_res.csv`: All beta centroid pairwise comparisons

**CSV Format**:
| Column | Description |
|--------|-------------|
| Variable | Categorical variable name |
| Metric | Alpha metric (Observed_OTUs, Shannon_Diversity, Inverse_Simpson) or Beta metric (Braycurtis, Unwunifrac, Wunifrac) |
| Group1 | First group in comparison |
| Group2 | Second group in comparison |
| P_value | Raw P-value from Wilcoxon test |
| P_value_FDR | FDR-adjusted P-value (Benjamini-Hochberg) |
| Significant | TRUE if P_value_FDR < 0.05 |

**Example Results** (Full Dataset, n=9,349):
- Alpha diversity: 612 total comparisons, 400 significant (65.4%)
- Beta centroid: 218 total comparisons, 143 significant (65.6%)
- Age group comparisons: Up to 13 significant pairs for alpha diversity
- Ethnicity comparisons: 9 significant pairs for alpha diversity

### Taxa Filtering

Before testing, taxa are filtered:
- **For raw data (none)**: Present in ≥1% of samples
- **For CLR/lognorm**: Non-zero variance across samples
- **All transformations**: Must have non-missing Genus annotation

### Heatmap Display Logic

For each categorical variable:

1. **Identify significant taxa** (FDR ≤ 0.05)
2. **Rank by effect size** (range of group means: max - min)
3. **Select top 30 taxa** (or fewer if <30 significant)
4. **Calculate difference from reference level**:
   - **CLR**: (mean_level - mean_ref) / log(2) → log₂ fold-change
   - **Lognorm**: (mean_level - mean_ref) → log₁₀ difference
   - **None**: z-score across groups

5. **Hierarchical clustering** by Euclidean distance
6. **Annotate by Phylum** with semi-transparent color bars

### Color Scale Interpretation

| Transformation | Breaks | Color Meaning |
|---|---|---|
| **CLR** | -2.5 to +2.5 | **Red**: Higher than reference (up to 8-fold)<br>**Blue**: Lower than reference<br>**White**: No change |
| **Lognorm** | -1 to +1 | **Red**: Higher than reference<br>**Blue**: Lower than reference<br>**White**: No change |
| **None** | -1.5 to +1.5 | **Red**: Above genus average<br>**Blue**: Below genus average<br>**White**: At average (z-score relative to all groups) |

Palette: `#053061` (dark blue) → `#F7F7F7` (white) → `#67001F` (dark red)

## Usage Guide

### Running the Standalone Script

```bash
# Make executable
chmod +x scripts/standalone_differential_abundance_heatmap.R

# Run with default paths
Rscript scripts/standalone_differential_abundance_heatmap.R
```

### Customizing Paths

Edit these variables in the script:

```r
# Line 21-24
base_path <- "/path/to/your/project"
phyloseq_obj_path <- file.path(base_path, "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate")
viz_out_path <- file.path(base_path, "results/analyses_results/04_association_phyloseq_out/visualizations")
config_dir_path <- file.path(base_path, "config")
```

### Customizing Parameters

```r
# Line 519-521 (in run_dataset function)
min_prev = 0.01   # Minimum prevalence (1% of samples)
top_n = 30        # Max taxa to display per heatmap
alpha = 0.05      # FDR threshold
```

### Adding New Categorical Variables

1. **Add variable to config file**:
   ```bash
   echo "NEW_VARIABLE" >> config/1_demoWAS_vars.txt
   ```

2. **Ensure variable exists in phyloseq sample_data**

3. **Add factorization code** (if needed):
   ```r
   # In section 4 of script
   New_category = factor(NEW_VARIABLE, 
     levels = c("ref_level", "other_level")
   ) %>% fct_relevel("ref_level")
   ```

4. **Add to appropriate dataset**:
   ```r
   select(SEQN, ..., New_category)
   ```

### Output Files

For each transformation (none/clr/lognorm) and dataset combination:

1. **Individual PDFs**: One per categorical variable (e.g., `demoWAS_Gender_clr.pdf`)
2. **CSV tables**: Full results including non-significant taxa (e.g., `demoWAS_Gender_clr_table.csv`)
3. **Booklet PDF**: All heatmaps combined (e.g., `ALL_clr_booklet.pdf`)

CSV table columns:
- `OTU`: Taxon ID
- `Genus`: Genus name (or "Uncl <OTU>" if unclassified)
- `Phylum`: Phylum assignment
- `p_value`: Raw p-value from test
- `p_adj_BH`: Benjamini-Hochberg adjusted p-value
- `<level1>`, `<level2>`, etc.: Mean abundance per categorical level
- `in_top_panel`: TRUE if taxon is in top 30 displayed in heatmap

## Key Differences from Full Pipeline

The standalone script differs from the full `4_association_phyloseq_analyses.Rmd` workflow:

### Standalone Script ✓
- Loads **pre-computed** phyloseq objects from disk
- Uses **minimal** factorization (only what's needed for heatmaps)
- **No database queries** (all data already in phyloseq)
- **Fast execution** (~5-10 minutes)
- Self-contained, portable

### Full Pipeline
- Connects to SQLite database
- Retrieves variables from F/G cycle tables
- Updates phyloseq sample_data with new variables
- Computes associations (GLM/linear models)
- Saves updated phyloseq objects
- Generates heatmaps + network plots + other visualizations
- Long execution time (hours)

## Troubleshooting

### Error: "phyloseq object not found"
**Solution**: Check that phyloseq_obj_path points to correct directory with .rds files

### Error: "no factor vars"
**Solution**: Check that variables in config files exist in phyloseq sample_data

### Error: "too few samples"
**Solution**: Some categorical levels have <30 samples. Either:
- Collapse rare levels
- Exclude the variable
- Reduce minimum sample threshold

### Error: "no taxa left"
**Solution**: Minimum prevalence filter is too strict. Reduce min_prev parameter.

### Heatmap shows only a few taxa
**Solution**: Few taxa pass FDR threshold. This is expected for weak associations. Consider:
- Using different transformation (CLR often more powerful)
- Checking if variable is truly associated with microbiome
- Reviewing raw p-values in CSV table

## References

### Key Statistical Methods
- **Wilcoxon rank-sum test**: Non-parametric test for two groups
- **Kruskal-Wallis test**: Non-parametric ANOVA for >2 groups
- **Benjamini-Hochberg procedure**: FDR control for multiple testing
- **CLR transformation**: Aitchison (1982) compositional data analysis
- **Hierarchical clustering**: Euclidean distance, complete linkage

### Color Palette Source
- Diverging blue-white-red: ColorBrewer RdBu scheme
- Phylum colors: Custom palette based on phylogenetic distance

---

**Last Updated**: October 2025  
**Script Version**: 1.0  
**Compatible with**: R ≥ 4.2, phyloseq ≥ 1.42

