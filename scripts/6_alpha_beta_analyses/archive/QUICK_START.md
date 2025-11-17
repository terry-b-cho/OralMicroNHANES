# Quick Start Guide - Integrated Diversity Analysis

## Run Test Mode (Fast - 10 minutes)

```bash
cd /Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho
mamba activate oral_env

# Run both scripts (TEST_MODE = TRUE by default)
Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R
```

**Output**: ~10 PDFs in `results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/`

---

## Run Full Mode (Complete - 45 minutes)

### Step 1: Enable Full Mode

Edit both files and change `TEST_MODE <- TRUE` to `TEST_MODE <- FALSE`:

**File 1**: `scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R` (line 32)
**File 2**: `scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R` (line 36)

### Step 2: Run Analysis

```bash
cd /Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho
mamba activate oral_env

# Step 1: Load full dataset (~5 min)
Rscript scripts/6_alpha_beta_analyses/1_load_diversity_and_all_categories.R

# Step 2: Run full analysis (~35-45 min)
Rscript scripts/6_alpha_beta_analyses/4_integrated_diversity_FINAL.R
```

**Output**: ~33 PDFs (one per categorical variable)

---

## What You Get

### Each PDF Contains (4×3 grid, 12 panels):

```
Row 1 (A-C): Alpha Diversity
  - Observed OTUs, Shannon, Inverse Simpson
  - Violin+Box plots with pairwise brackets
  - K-W: P=X, R²=Y
    H=Z, df=W, N=N

Row 2 (D-F): Beta Centroid Distances
  - Bray-Curtis, Unwunifrac, Wunifrac
  - Violin+Box plots
  - PERMANOVA: P=X, R²=Y
    F=Z, df=W,W, N=N

Row 3 (G-I): PCoA Ordination
  - PC1 vs PC2 scatter
  - 95% ellipses + centroids
  - P=X, R²=Y

Row 4 (J-L): Scree Plots
  - Variance explained (grayscale)
  - First 5 PCs shown
```

### PDF Specs:
- **Size**: 77.2mm × 137.2mm (3.04" × 5.4")
- **Format**: Vector PDF, 300 DPI
- **Colors**: Safe Grafify palette (colorblind-friendly)

---

## Variables Analyzed (33 total when FULL MODE)

**Demographics**: Gender, Age, Education, Ethnicity, Income, Household size, etc.  
**Oral Health**: Denture, Gum disease, Hygiene, Tooth decay  
**Lifestyle**: Smoking status, Hepatitis, HPV  
**Phenotypes**: BMI, Blood pressure, Waist circumference  
**Outcomes**: Asthma, Heart disease, Stroke, Diabetes, etc.  

---

## Manuscript Use

**Main Text**: Select 3-5 key PDFs  
**Supplement**: Include all ~33 PDFs  
**Tables**: Export PERMANOVA results CSV  
**Text**: Copy from `4_complete_integrated_diversity_FINAL_explanation.md` (18,500 words)

---

**Need Help?** See `README_FULL_ANALYSIS.md` for complete instructions.

