# Complete Diversity Analysis - Full Dataset Instructions

## Overview

This directory contains a complete, publication-ready alpha and beta diversity analysis pipeline for the NHANES oral microbiome dataset. The analysis can be run in two modes:

- **TEST MODE** (default): 500-1,000 samples, ~7-10 minutes total runtime
- **FULL MODE**: ALL 9,662 samples, ~35-50 minutes total runtime

## Quick Start

### Test Mode (Default - Recommended for Development)

```bash
cd /Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho

# Activate conda environment
mamba activate oral_env

# Run data loading (TEST MODE - 1,000 samples)
Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R

# Run integrated analysis (TEST MODE - 500 samples)
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R
```

**Expected Output (TEST MODE)**:
- Runtime: ~7-10 minutes
- PDFs created: 10 integrated plots (top 10 significant categorical variables)
- Dimensions: 77.2mm × 137.2mm per PDF
- Sample size: 500-1,000 samples

---

## Full Dataset Analysis (For Final Publication)

### Step 1: Update Configuration in BOTH Scripts

#### File 1: `1_load_diversity_and_all_categories.R`

Open the file and change line 32:

```r
# BEFORE (Test Mode):
TEST_MODE <- TRUE    # ← SET TO FALSE FOR FULL DATASET

# AFTER (Full Mode):
TEST_MODE <- FALSE   # ← FULL DATASET ENABLED
```

#### File 2: `4_integrated_diversity_FINAL.R`

Open the file and change line 36:

```r
# BEFORE (Test Mode):
TEST_MODE <- TRUE    # ← SET TO FALSE FOR FULL DATASET ANALYSIS

# AFTER (Full Mode):
TEST_MODE <- FALSE   # ← FULL DATASET ENABLED
```

### Step 2: Run Full Analysis

```bash
cd /Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho

# Activate conda environment
mamba activate oral_env

# Step 1: Load full dataset (9,662 samples, ~5 min)
Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R

# Step 2: Run integrated analysis (ALL samples, ALL 33 variables, ~30-45 min)
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R
```

### Step 3: Monitor Progress

The script will print progress messages:

```
=== INTEGRATED DIVERSITY ANALYSIS - FINAL ===

TEST_MODE: FALSE

Loading beta diversity...
  Braycurtis: PC1=12.3%, PC2=7.8%
  Unwunifrac: PC1=9.8%, PC2=6.4%
  Wunifrac: PC1=11.2%, PC2=7.1%
  ✓ PCoA complete

Calculating centroids...
  ✓ Done

Running PERMANOVA...
  ✓ 99 tests, 78 significant

Creating integrated 2×3 plots...
    Added 6 pairwise comparison brackets (buffer: 87 %)
  ✓ Integrated_Gender.pdf (4×3 layout, 77.2×137.2mm)
  ✓ Integrated_Age_group.pdf (4×3 layout, 77.2×137.2mm)
  ...
  [Progress continues for all 33 variables]
  
✓ Complete!
```

---

## Expected Output (FULL MODE)

### Files Created

**Location**: `/results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/`

**Number of PDFs**: ~33 integrated plots (one per categorical variable)

**File naming**: `Integrated_[Variable_Name].pdf`

**Examples**:
```
Integrated_Gender.pdf
Integrated_Age_group.pdf
Integrated_Education_level.pdf
Integrated_Ethnicity.pdf
Integrated_US_born.pdf
Integrated_Income_to_poverty_ratio.pdf
Integrated_Household_size.pdf
Integrated_Marital_status.pdf
Integrated_Interview_language.pdf
Integrated_Denture.pdf
Integrated_Gum_disease.pdf
Integrated_Oral_hygene.pdf
Integrated_Tooth_decay.pdf
Integrated_Smoking_status.pdf
Integrated_Days_smoked_in_5_days.pdf
Integrated_Cigarettes_per_day.pdf
Integrated_Last_smoked_cigarette.pdf
Integrated_Hepatitis_C_antibody.pdf
Integrated_HPV_PCR_summary.pdf
Integrated_BMI_category.pdf
Integrated_Blood_pressure.pdf
Integrated_Pulse_category.pdf
Integrated_Waist_circumference_quintiles.pdf
Integrated_Asthma.pdf
Integrated_Bronchitis.pdf
Integrated_Emphysema.pdf
Integrated_Angina.pdf
Integrated_Heart_failure.pdf
Integrated_Heart_attack.pdf
Integrated_Stroke.pdf
Integrated_Diabetes.pdf
... and possibly more
```

### PDF Specifications

- **Dimensions**: 77.2mm × 137.2mm (3.04" × 5.4")
  - Width: 77.2mm (10% narrower than original)
  - Height: 137.2mm (20% taller than original)
- **Resolution**: 300 DPI
- **Format**: Vector PDF (publication-ready)
- **Layout**: 4×3 grid (12 panels: A-L)

### Panel Structure per PDF

```
┌─────────────┬─────────────┬─────────────┐
│ A. Observed │ B. Shannon  │ C. Inverse  │  Row 1: Alpha Diversity
│    OTUs     │  Diversity  │   Simpson   │  K-W: P=X, R²=Y
│  + Brackets │  + Brackets │  + Brackets │  H=Z, df=W, N=N
├─────────────┼─────────────┼─────────────┤
│ D. Bray-    │ E. Unwuni-  │ F. Wuni-    │  Row 2: Beta Centroid
│    Curtis   │    frac     │    frac     │  PERMANOVA: P=X, R²=Y
│  Centroid   │  Centroid   │  Centroid   │  F=Z, df=W,W, N=N
├─────────────┼─────────────┼─────────────┤
│ G. Bray-    │ H. Unwuni-  │ I. Wuni-    │  Row 3: PCoA Ordination
│    Curtis   │    frac     │    frac     │  P=X, R²=Y
│    PCoA     │    PCoA     │    PCoA     │  (with ellipses+centroids)
├─────────────┼─────────────┼─────────────┤
│ J. Bray-    │ K. Unwuni-  │ L. Wuni-    │  Row 4: Scree Plots
│    Curtis   │    frac     │    frac     │  (variance explained,
│    Scree    │    Scree    │    Scree    │   grayscale)
└─────────────┴─────────────┴─────────────┘
         Legend (group colors)
```

---

## All 33 Categorical Variables Analyzed

### Demographics (9 variables)
1. **Gender** - Reference: Female
2. **Age_group** - Reference: 30-39 years
3. **AgeGroup** - Reference: 14-19 years (alternative coding)
4. **age_group** - Reference: 14–19 years (another alternative)
5. **Education_level** - Reference: College/AA degree
6. **Ethnicity** - Reference: White
7. **US_born** - Reference: US Born
8. **Income_to_poverty_ratio** - Reference: 150-184%
9. **Household_size** - Reference: 4 members
10. **Marital_status** - Reference: Married
11. **Interview_language** - Reference: English

### Oral Health (4 variables)
12. **Denture** - Reference: control (no denture)
13. **Gum_disease** - Reference: control
14. **Oral_hygene** - Reference: control
15. **Tooth_decay** - Reference: control

### Exposures/Lifestyle (7 variables)
16. **Smoking_status** - Reference: Never smoker
17. **Days_smoked_in_5_days** - Reference: Never smoked
18. **Cigarettes_per_day** - Reference: Never smoked
19. **Last_smoked_cigarette** - Reference: (varies)
20. **Hepatitis_C_antibody** - Reference: Negative
21. **HPV_PCR_summary** - Reference: Negative
22. **Hepatitis_B_surface_antigen** - Reference: Negative

### Phenotypes/Measurements (5 variables)
23. **BMI_category** - Reference: Healthy weight
24. **Blood_pressure** - Reference: Normal
25. **Pulse_category** - Reference: 70-75 bpm
26. **Waist_circumference_quintiles** - Reference: Q3

### Disease Outcomes (8 variables)
27. **Asthma** - Reference: control
28. **Bronchitis** - Reference: control
29. **Emphysema** - Reference: control
30. **Angina** - Reference: control
31. **Heart_failure** - Reference: control
32. **Heart_attack** - Reference: control
33. **Stroke** - Reference: control
34. **Diabetes** - Reference: control

**Note**: The exact number may vary slightly (30-35) depending on which variables have sufficient data (≥30 samples per level).

---

## Statistical Annotation Format (NEW!)

### Alpha Diversity Subtitle (2 lines):
```
K-W: P=<0.001, R²=0.065
H=67.84, df=5, N=996
```

**Line 1**: Key results (P-value and effect size)  
**Line 2**: Test details (statistic, df, sample size)

### Beta Diversity Subtitle (2 lines):
```
PERMANOVA: P=0.001, R²=0.028
F=14.32, df=5,990, N=996
```

**Line 1**: Key results (P-value and variance explained)  
**Line 2**: Test details (F-statistic, df, sample size)

### PCoA Subtitle (1 line):
```
P=0.001, R²=0.028
```

---

## Runtime Estimates

| Mode | Samples | Variables | Runtime | Memory |
|------|---------|-----------|---------|--------|
| **TEST** | 500-1,000 | Top 10 | 7-10 min | ~2 GB |
| **FULL** | 9,662 | ALL 33 | 35-50 min | ~8 GB |

### Breakdown (FULL MODE):

1. **Data Loading** (`1_load_diversity_and_all_categories.R`): ~5 minutes
   - Load alpha diversity (9,662 × 201 columns)
   - Load phyloseq sample data
   - Factorize 33 categorical variables
   - Merge and save

2. **Beta Diversity Loading**: ~2 minutes
   - Load 3 distance matrices (9,349 × 9,349 each)
   - Subset to alpha diversity samples

3. **PCoA Computation**: ~3 minutes
   - Perform PCoA for 3 metrics
   - Eigenvalue decomposition

4. **Centroid Calculation**: ~5 minutes
   - Calculate distances for 33 variables × 3 metrics

5. **PERMANOVA Testing**: ~10 minutes
   - Run 99 tests (33 variables × 3 metrics)
   - 199 permutations each

6. **Plot Generation**: ~15-25 minutes
   - Create 33 PDFs
   - Each with 12 panels
   - Pairwise comparisons for each panel

**Total**: ~35-50 minutes

---

## Output Quality Metrics

### Per PDF (4×3 layout):
- ✅ 12 panels (A-L)
- ✅ 6 diversity metrics
- ✅ ALL pairwise comparisons shown (P<0.05)
- ✅ Comprehensive statistical annotations (2-line format)
- ✅ 95% confidence ellipses (PCoA)
- ✅ Plus sign (+) centroids
- ✅ Grayscale scree plots
- ✅ Shared legend
- ✅ Safe Grafify colors

### Statistical Rigor:
- ✅ Kruskal-Wallis tests (alpha diversity)
- ✅ Epsilon-squared effect sizes
- ✅ Pairwise Wilcoxon tests (FDR-corrected)
- ✅ PERMANOVA (beta diversity)
- ✅ R² variance explained
- ✅ Pseudo-F statistics
- ✅ All P-values with 3 significant figures

---

## Troubleshooting

### Issue: "Out of memory" error

**Solution**: Increase system memory or run in batches

```r
# Option 1: Reduce test sample size
TEST_SAMPLES <- 300  # Instead of 500

# Option 2: Run in batches (manual)
# Edit line 228 in 4_integrated_diversity_FINAL.R to process specific variables:
# top_vars <- c("Gender", "Age_group", "Ethnicity")  # Specify variables
```

### Issue: "Too slow" in full mode

**Solution**: The full dataset is large. Expected runtime is 35-50 minutes. This is normal.

To speed up:
- Close other applications
- Use command line (not RStudio)
- Run overnight if needed

### Issue: "Missing phyloseq object"

**Solution**: Ensure prerequisite data exists

```bash
# Check if phyloseq object exists
ls results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_none.rds

# If missing, run preprocessing script first
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R
```

---

## File Locations

### Input Data

**Alpha Diversity**:
- `data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt` (9,662 samples × 201 metrics)

**Beta Diversity**:
- `data/00_nhanes_omp_diversity_db/dada2rsv-beta-braycurtis.txt` (9,349 × 9,349)
- `data/00_nhanes_omp_diversity_db/dada2rsv-beta-unwunifrac.txt` (9,349 × 9,349)
- `data/00_nhanes_omp_diversity_db/dada2rsv-beta-wunifrac.txt` (9,349 × 9,349)

**Metadata**:
- `results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_none.rds`

### Output Data

**Processed Data**:
- `scripts/6_alpha_beta_analyses/data/alpha_diversity_with_all_categories.rds`
- `scripts/6_alpha_beta_analyses/data/all_categorical_data.rds`
- `scripts/6_alpha_beta_analyses/data/categorical_variables_info.csv`

**Figures**:
- `results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/Integrated_*.pdf`

**Statistical Results**:
- `results/analyses_results/6_alpha_beta_analyses_out/permanova_all_variables.csv`
- `results/analyses_results/6_alpha_beta_analyses_out/beta_diversity_centroid_distances_all.csv`

---

## Expected Results (FULL MODE)

### Number of PDFs Created

**Approximately 30-33 PDFs** (one per categorical variable)

The exact number depends on:
- Variables with ≥30 samples per level
- Variables with 2-10 factor levels
- Variables with significant PERMANOVA results (P<0.05)

### Statistical Power

With n=9,662:
- **>99.9% power** to detect small effects (ε²≥0.005, R²≥0.005)
- **Perfect power** for medium/large effects
- **Highly robust** pairwise comparisons

### Expected Significant Findings

Based on literature and test mode:

**Strong Associations** (P<0.001, R²>0.02):
- Age group
- Smoking status
- Ethnicity

**Moderate Associations** (P<0.01, R²>0.01):
- Household size
- Income-to-poverty ratio
- Oral health variables
- BMI category

**Small Associations** (P<0.05, R²>0.005):
- Education level
- Blood pressure
- Disease outcomes

---

## Data Quality Checks

The scripts automatically perform these checks:

1. ✅ **SEQN matching** between diversity data and metadata
2. ✅ **Numeric type validation** for diversity metrics
3. ✅ **Factor level verification** for categorical variables
4. ✅ **Sample size filtering** (≥30 per group)
5. ✅ **Level count filtering** (2-10 levels per variable)
6. ✅ **Missing data handling** (complete case analysis)

---

## Comparison: Test vs Full Mode

| Feature | TEST MODE | FULL MODE |
|---------|-----------|-----------|
| **Samples** | 500-1,000 | 9,662 |
| **Variables plotted** | Top 10 | ALL 33 |
| **Runtime** | 7-10 min | 35-50 min |
| **Memory** | ~2 GB | ~8 GB |
| **PDFs created** | ~10 | ~33 |
| **Power** | 80-90% | >99% |
| **Use case** | Development, testing | Publication, final results |

---

## Post-Analysis Steps

### 1. Verify Output

```bash
# Count PDFs created
ls results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/Integrated_*.pdf | wc -l

# Check file sizes
ls -lh results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/ | head -10
```

### 2. Select Key Figures for Manuscript

Review all PDFs and select 3-5 most important for main text:
- Age_group (strongest effect)
- Ethnicity (demographic variation)
- Smoking_status (lifestyle factor)
- BMI_category (phenotype)
- Diabetes (disease outcome)

Move others to supplementary materials.

### 3. Update Manuscript Text

Use the complete documentation:
- `scripts/6_alpha_beta_analyses/4_complete_integrated_diversity_FINAL_explanation.md` (18,500 words)

Contains ready-to-use text for:
- Data Description
- Methods Section
- Figure Descriptions
- Results Section

### 4. Export Data for Tables

```r
# In R console or script:
library(dplyr)

# Load PERMANOVA results
permanova <- read.csv("results/analyses_results/6_alpha_beta_analyses_out/permanova_all_variables.csv")

# Create supplementary table
supp_table <- permanova %>%
  arrange(Metric, desc(R2)) %>%
  mutate(
    P_formatted = ifelse(P_value < 0.001, "<0.001", format(P_value, digits=3)),
    R2_formatted = format(R2, digits=3, nsmall=3)
  ) %>%
  select(Metric, Variable, R2_formatted, P_formatted, N_samples)

write.csv(supp_table, "results/supplementary_tables/Table_S_PERMANOVA_Results.csv", row.names=FALSE)
```

---

## Common Workflows

### Workflow 1: Quick Test (Before Full Run)

```bash
# Set TEST_MODE = TRUE in both scripts
# Run test analysis
Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R

# Review output (~10 PDFs)
open results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/

# If satisfied, proceed to full analysis
```

### Workflow 2: Full Analysis (For Submission)

```bash
# Set TEST_MODE = FALSE in both scripts
# Run full analysis
Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R

# Review all ~33 PDFs
open results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/

# Select key figures for manuscript
# Export supplementary tables
```

### Workflow 3: Update After Reviewer Comments

```bash
# If reviewers request specific categorical variables or analyses:

# 1. Modify variable selection in 4_integrated_diversity_FINAL.R (line ~228)
top_vars <- c("Gender", "Age_group", "Smoking_status", "BMI_category")

# 2. Re-run only the integrated analysis (data already loaded)
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R

# 3. Updated PDFs will be generated for specified variables only
```

---

## Technical Notes

### Memory Management

**Test Mode**:
- ~2 GB RAM required
- Safe for most laptops

**Full Mode**:
- ~8 GB RAM required
- Recommended: 16 GB system memory
- Close other applications during analysis

### Processing Time

The longest steps in FULL MODE:

1. **PERMANOVA** (33 variables × 3 metrics × 199 permutations): ~10 min
2. **Plot generation** (33 PDFs × 12 panels): ~15-25 min
3. **Pairwise comparisons** (many groups × many pairs): ~10 min

Total: ~35-50 minutes

### Parallelization

Currently sequential. To speed up:

```r
# Option: Reduce permutations (line ~194 in 4_integrated_diversity_FINAL.R)
# BEFORE:
perm_result <- adonis2(formula_obj, data = meta_subset, permutations = 199)

# AFTER (faster but less precise):
perm_result <- adonis2(formula_obj, data = meta_subset, permutations = 99)
```

---

## Validation Checklist

Before submitting results:

- [ ] TEST_MODE set to FALSE in BOTH scripts
- [ ] Both scripts run without errors
- [ ] ~33 PDFs created in output directory
- [ ] Each PDF is 77.2mm × 137.2mm
- [ ] Statistical annotations in 2-line format
- [ ] All pairwise brackets visible
- [ ] PERMANOVA results CSV created
- [ ] Manuscript text updated with full dataset n=9,662

---

## Support and Documentation

**Full Documentation**: 
- `4_complete_integrated_diversity_FINAL_explanation.md` - 18,500 words
  - Complete Data Description
  - Complete Methods Section
  - Complete Figure Descriptions
  - Complete Results Section

**Quick Reference**:
- `STATISTICAL_ANNOTATIONS_SUMMARY.md` - Statistical methods guide
- `COMPLETE_FINAL_SUMMARY.md` - Analysis overview

**Questions?**
- Check inline comments in R scripts
- Review documentation files
- Verify file paths and permissions

---

## Success Criteria

✅ **Analysis Complete When**:
1. Both scripts run without errors
2. ~33 integrated PDFs created
3. Each PDF contains 12 panels
4. Statistical annotations in new 2-line format
5. All pairwise comparisons visible
6. File dimensions: 77.2mm × 137.2mm

✅ **Ready for Publication When**:
1. TEST_MODE = FALSE used
2. All 9,662 samples analyzed
3. All 33 categorical variables processed
4. Manuscript sections updated
5. Supplementary tables exported
6. Key figures selected (3-5 for main text)

---

**Last Updated**: October 8, 2025  
**Script Version**: 4.0 (FINAL with 2-line subtitles)  
**Status**: ✅ Production-ready for full dataset analysis

