# NHANES Oral Microbiome Association Analysis Pipeline - Production Documentation

## Overview

This documentation covers the **production-ready NHANES oral microbiome association analysis pipeline** consisting of 6 core scripts designed for Wide Association Studies (WAS) with proper survey design, multiple comparisons correction, and FLEXIBLE resume functionality.

## **Core Scripts**

### **1. Environment Setup**

- **`setup_nhanes_environment.sh`** - Conda environment creation with all required R packages

### **2. Individual Analysis Scripts**

- **`universal_was_analysis_PRODUCTION.R`** - Core statistical analysis engine (survey-weighted regression)
- **`run_single_was_analysis_PRODUCTION.sh`** - Single analysis type submission script

### **3. Batch Processing Scripts**

- **`run_all_was_analyses_flexible.sh`** - FLEXIBLE resume batch submission with corrected resource allocation

### **4. Result Processing Scripts**

- **`aggregate_was_results.R`** - Results aggregation with multiple comparisons correction

### **5. Testing and Validation**

- **`debug_all_pipelines.R`** - Comprehensive pipeline testing and validation

---

## **Analysis Types and Framework**

### **Analysis Types Available:**


| Pipeline      | Description                    | Dependent Variable       | Independent Variable     | Model Type |
| --------------- | -------------------------------- | -------------------------- | -------------------------- | ------------ |
| **1_demoWAS** | Demographics → Microbiome     | Microbiome OTU abundance | Demographic Variables    | Linear     |
| **2_oradWAS** | Microbiome → Oral Health      | Oral health outcomes     | Microbiome OTU abundance | Logistic   |
| **3_exWAS**   | Exposures → Microbiome        | Microbiome OTU abundance | Exposure variables       | Linear     |
| **4_pheWAS**  | Microbiome → Phenotypes       | Phenotype measures       | Microbiome OTU abundance | Linear     |
| **5_outWAS**  | Microbiome → Disease Outcomes | Disease outcomes         | Microbiome OTU abundance | Logistic   |
| **6_zimWAS**  | Microbiome → Lab Measurements | Lab measurements         | Microbiome OTU abundance | Linear     |

### **Normalization Methods:**

- **`clr`** - Uses CLR-transformed microbiome tables (`*_clr`)
- **`lognorm`** - Uses log-normal transformed microbiome tables (`*_lognorm`)
- **`none`** - Uses raw microbiome tables (`*_none`) with pseudocount addition (1e-6); pseudo count added at the point of regression

### **Computational Requirements (Corrected):**

Based on actual regression counts per dependent variable:


| Analysis Type                               | Regressions/dep_var | Memory | Time | Reasoning                          |
| --------------------------------------------- | --------------------: | -------- | ------ | ------------------------------------ |
| **2_oradWAS, 4_pheWAS, 5_outWAS, 6_zimWAS** |               1,349 | 12G    | 4h   | Microbiome as independent variable |
| **3_exWAS**                                 |                 473 | 8G     | 3h   | Microbiome as dependent variable   |
| **1_demoWAS**                               |                  26 | 1G     | 1h   | Microbiome as dependent variable   |

---

## **Quick Start Guide**

### **Step 1: Environment Setup (One-time)**

```bash
# Start interactive session
srun --pty -p interactive -t 12:00:00 --mem=32G bash

# Setup conda environment
bash scripts/1_association_pipeline/setup_nhanes_environment.sh
```

### **Step 2: Validation (Recommended)**

```bash
# Load environment
module purge
module load gcc/9.2.0      # Must be first
module load miniconda3/23.1.0  # Must be last
source activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Run comprehensive validation
Rscript scripts/1_association_pipeline/debug_all_pipelines.R
```

### **Step 3: Run Analysis**

#### **Option A: Flexible (Resume) Batch Processing (Recommended)**

```bash
# Exit interactive session first
exit

# Quick commandline documentation
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh --help

# Run all 18 analyses (6 types × 3 normalizations) with smart resume
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh

# Test mode (15 minutes per job)
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh test

# Aggressive resources (32G, 12h for all jobs)
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh "" aggressive
```

#### **Option B: Individual Analysis**

```bash
# Run single analysis type
bash scripts/1_association_pipeline/run_single_was_analysis_PRODUCTION.sh 1_demoWAS clr

# Test mode
bash scripts/1_association_pipeline/run_single_was_analysis_PRODUCTION.sh 1_demoWAS clr test
```

### **Step 4: Result Aggregation**

```bash
# After ALL individual analyses complete, aggregate results
Rscript scripts/1_association_pipeline/aggregate_was_results.R 1_demoWAS clr /path/to/results
```

---

## **Script Details**

### **1. setup_nhanes_environment.sh**

**Purpose:** Creates conda environment with all required R packages for survey analysis.

**Key Features:**

- Checks for compute node availability
- Creates environment from yml file or individual package installation
- Validates all required packages are installed
- Provides activation instructions

**Required Packages:**

- `survey`, `broom`, `dplyr`, `glue`, `DBI`, `RSQLite`, `getopt`, `logger`, `readr`, `tibble`, `tidyr`, `magrittr`

**Usage:**

```bash
bash scripts/1_association_pipeline/setup_nhanes_environment.sh
```

### **2. universal_was_analysis_PRODUCTION.R**

**Purpose:** Core statistical analysis engine performing survey-weighted regression analysis.

**Key Features:**

- **Survey Design Implementation:** Proper NHANES stratification (`SDMVSTRA`) and clustering (`SDMVPSU`)
- **FLEXIBLE Covariate Selection:** Automatically assesses and uses maximum usable covariates
- **Microbiome Data Handling:** Proper pseudocount addition for 'none' normalization only
- **Robust Error Handling:** Comprehensive error recovery and fallback mechanisms
- **Enhanced Zero Handling:** Adds pseudocount (1e-6) to ALL values for microbiome data with 'none' normalization

**Statistical Methods:**

- **Linear Regression:** For continuous dependent variables (gaussian family)
- **Logistic Regression:** For binary dependent variables (quasibinomial family)
- **Survey Weights:** Uses `WTMEC2YR` with proper scaling
- **R-squared Calculation:** Manual calculation for survey models with robust error handling

**Command Line Arguments:**

```bash
Rscript universal_was_analysis_PRODUCTION.R \
    --dependent_var "RSV_genus1002_relative" \
    --schema_structure_file "schema.csv" \
    --database_path "database.sqlite" \
    --output_path "output_dir" \
    --analysis_type "1_demoWAS" \
    --normalization "clr" \
    [--test]
```

**Output Structure:**

- `pe_tidied` - Coefficient results with metadata
- `pe_glanced` - Model-level statistics
- `rsq` - R-squared values

### **3. run_single_was_analysis_PRODUCTION.sh**

**Purpose:** Submits SLURM jobs for a single analysis type and normalization.

**Key Features:**

- Per-dependent-variable job submission (matches original pipeline structure)
- Automatic resource allocation based on analysis type
- Test mode support with reduced resources
- Comprehensive job monitoring instructions

**Resource Allocation:**

- **3_exWAS:** 8G memory, 2 hours
- **4_pheWAS:** 4G memory, 1 hour
- **Other analyses:** 2G memory, 30 minutes
- **Test mode:** 2G memory, 15 minutes

**Usage:**

```bash
bash run_single_was_analysis_PRODUCTION.sh <analysis_type> <normalization> [test]
```

### **4. run_all_was_analyses_flexible.sh**

**Purpose:** FLEXIBLE resume batch submission system with corrected resource allocation.

**Key Features:**

- **EFFICIENT FLEXIBLE Resume:** Only submits jobs for missing .rds files
- **NEVER Overwrites:** Preserves all existing successful results
- **Corrected Resource Allocation:** Based on actual computational analysis
- **File Detection Speed:** ~5 seconds vs 3+ hours previously
- **Safety Mechanisms:** Race condition protection and comprehensive logging

**FLEXIBLE Resume Logic:**

1. **Efficient Missing File Detection:** Uses Unix utilities (`find`, `cut`, `comm`)
2. **File Existence Check:** Simple file size validation (>100 bytes)
3. **Atomic Job Submission:** Double-checks files don't exist before submission
4. **Comprehensive Logging:** All operations logged to organized directory

**Resource Levels:**

- **auto (default):** Standard allocation based on analysis requirements
- **aggressive:** 32G memory, 12 hours for all job types

**Usage:**

```bash
./run_all_was_analyses_flexible.sh [test_mode] [resource_level]
```

**Log Directory:** `logs/1_association_pipeline_flexible/`

### **5. aggregate_was_results.R**

**Purpose:** Aggregates individual results with multiple comparisons correction.

**Key Features:**

- **FDR Correction:** Applied WITHIN each dependent variable (microbiome OTU)
- **Safe File Reading:** Robust handling of corrupted or incomplete files
- **Multiple Correction Methods:** Both FDR (Benjamini-Hochberg) and Bonferroni
- **Comprehensive Reporting:** Detailed summary statistics and top significant results

**Correction Strategy:**

- FDR correction applied by grouping on dependent variable using `group_by()`
- Only main effect terms (`indep_var`) are corrected
- Maintains FDR < 0.05 as standard significance threshold

**Command Line Usage:**

```bash
Rscript aggregate_was_results.R <analysis_type> <normalization> <results_dir>
```

**Output Files:**

- `*_tidied_complete.rds` - Main results with corrected p-values
- `*_glanced_complete.rds` - Model statistics
- `*_rsq_complete.rds` - R-squared values
- `*_aggregation_summary.txt` - Summary report

### **6. debug_all_pipelines.R**

**Purpose:** Comprehensive testing and validation of all pipeline components.

**Key Features:**

- **Complete Coverage:** Tests all 18 analysis combinations (6 types × 3 normalizations)
- **Schema Validation:** Checks for required schema structure files
- **Pipeline Integration:** Tests full analysis workflow
- **Output Validation:** Verifies result structure and metadata
- **Production Readiness Assessment:** Success rate calculation and recommendations

**Test Components:**

1. Schema file existence and structure
2. Dependent variable extraction
3. Full pipeline execution in test mode
4. Output file creation and validation
5. Result structure verification
6. Metadata column presence

**Success Criteria:**

- **Excellent (≥90%):** Ready for production deployment
- **Good (≥75%):** Mostly ready, minor issues to address
- **Moderate (≥50%):** Significant issues need resolution
- **Poor (<50%):** Major issues prevent deployment

**Usage:**

```bash
Rscript scripts/1_association_pipeline/debug_all_pipelines.R
```

---

## 📁 **File Structure and Outputs**

### **Schema Structure Files:**

```
results/0_ss_files/
├── 1_demoWAS_clr_schema_structure.csv
├── 1_demoWAS_lognorm_schema_structure.csv
├── 1_demoWAS_none_schema_structure.csv
└── ... (18 total files)
```

### **Individual Results:**

```
results/<analysis_type>_out/result_<normalization>/
├── RSV_genus1002_relative.rds    # For 1_demoWAS, 3_exWAS
├── DENTURE_OHAROCDE.rds          # For 2_oradWAS
├── BMXBMI.rds                    # For 4_pheWAS, 6_zimWAS
└── DIQ010.rds                    # For 5_outWAS
```

### **Aggregated Results:**

```
results/<analysis_type>_out/result_<normalization>/
├── <analysis_type>_<normalization>_tidied_complete.rds     ⭐ MAIN RESULTS
├── <analysis_type>_<normalization>_glanced_complete.rds    # Model stats
├── <analysis_type>_<normalization>_rsq_complete.rds        # R-squared
└── <analysis_type>_<normalization>_aggregation_summary.txt # Summary
```

### **Key Result Columns:**

```r
# Main results structure
results <- readRDS("*_tidied_complete.rds")

# Key columns:
# term                - Model term (focus on "indep_var" for main effects)
# estimate            - Effect size (coefficient)
# std.error           - Survey-adjusted standard error  
# p.value             - Original p-value
# p.value.fdr         - FDR-corrected p-value ⭐ USE FOR REPORTING
# p.value.bonferroni  - Bonferroni-corrected p-value
# phenotype           - Dependent variable name
# exposure            - Independent variable name
# dependent_var       - Dependent variable name (same as phenotype)
# independent_var     - Independent variable name (same as exposure)
# n_obs               - Sample size
# formula_used        - Actual formula used in analysis
# n_covariates        - Number of covariates included
```

---

## **Statistical Methods**

### **Survey Design Implementation:**

- **Complex Survey Design:** Proper NHANES stratification (`SDMVSTRA`) and clustering (`SDMVPSU`)
- **Survey Weights:** Uses `WTMEC2YR` with correct scaling across cycles
- **Survey Models:** Uses `svyglm()` from R `survey` package (NOT standard `lm()`/`glm()`)

### **Multiple Comparisons Correction:**

- **FDR:** Applied WITHIN each dependent variable (microbiome OTU)
- **Method:** Benjamini-Hochberg procedure using `p.adjust(method = "fdr")`
- **Implementation:** `group_by(dependent_var) %>% mutate(p.value.fdr = p.adjust(p.value, method = "fdr"))`
- **Significance Threshold:** FDR < 0.05

### **Microbiome Data Handling:**

- **Pseudocount Strategy:** Adds 1e-6 to ALL values for 'none' normalization only
- **Transformation Awareness:** Uses pre-transformed tables for 'clr' and 'lognorm'
- **Zero Handling:** Consistent treatment for both dependent and independent microbiome variables

### **Covariate Selection:**

- **FLEXIBLE Maximum Selection:** Automatically assesses data quality for each covariate
- **Quality Criteria:** ≤30% missing data, sufficient variation, no numerical issues
- **Fallback Strategy:** Progressive reduction with priority order (demographics → education → ethnicity)
- **Full Covariate Set:** `RIDAGEYR`, `AGE_SQUARED`, `RIAGENDR`, `INDFMPIR`, education indicators, ethnicity indicators, `BORN_INUSA`

---

## **Performance and Efficiency**

### **FLEXIBLE Resume Performance:**

- **File Detection:** ~5 seconds (vs 3+ hours previously)
- **Missing File Algorithm:** Unix utilities (`find`, `cut`, `comm`) instead of R validation
- **Memory Usage:** Minimal overhead for status checking
- **Safety:** Never overwrites existing successful results

### **Resource Optimization:**

- **Corrected Time Limits:** Based on actual computational requirements
- **Memory Allocation:** Appropriate for regression complexity
- **Cluster Efficiency:** Reduced job submission overhead

### **Expected Completion Times:**

- **Test Mode:** 15-30 minutes total
- **Production (auto resources):** ~6 hours
- **Production (aggressive resources):** ~4-5 hours

---

## **Troubleshooting**

### **Environment Issues:**

```bash
# If setup fails, check module loading order
module purge
module load gcc/9.2.0      # Must be first
module load miniconda3/23.1.0  # Must be last

# Manual package installation if needed
R -e "install.packages(c('survey', 'broom', 'dplyr'), repos='https://cloud.r-project.org/')"
```

### **Job Monitoring:**

```bash
# Check running jobs
squeue -u $USER | grep flexible

# Check for errors
find logs/ -name 'flexible_*.err' -size +0

# View execution log
tail -f logs/1_association_pipeline_flexible/flexible_was_analyses_*.log
```

### **Common Issues:**

#### **1. Memory/Time Failures:**

```bash
# Use aggressive mode for problematic analyses
./run_all_was_analyses_flexible.sh "" aggressive
```

#### **2. Missing Schema Files:**

```bash
# Check schema structure files exist
ls results/0_ss_files/*_schema_structure.csv
```

#### **3. Database Connection Issues:**

```bash
# Verify database exists and is accessible
ls -la data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
```

---

## ✅ **Success Indicators**

### **Environment Setup Successful:**

- ✅ All R packages load without errors
- ✅ Validation tests pass with >75% success rate
- ✅ Test analysis completes successfully

### **Analysis Successful:**

- ✅ All jobs complete (check with `squeue`)
- ✅ Individual .rds files created for each dependent variable
- ✅ Aggregated files contain expected number of results
- ✅ FDR-corrected p-values present in results
- ✅ Survey-adjusted standard errors reported

### **Results Ready for Publication:**

- ✅ Main results file: `*_tidied_complete.rds`
- ✅ Significant associations: `p.value.fdr < 0.05`
- ✅ Proper survey design implementation
- ✅ multiple comparisons correction applied

---

### **Methods Summary:**

**Survey Design:**

> "Analyses accounted for the complex survey design of NHANES using appropriate stratification, clustering, and survey weights. All models were fitted using survey-weighted generalized linear models (svyglm) from the R survey package."

**Multiple Comparisons:**

> "P-values were adjusted for multiple comparisons using the Benjamini-Hochberg false discovery rate (FDR) procedure applied within each dependent variable. Associations with FDR-adjusted p-values < 0.05 were considered statistically significant."

**Microbiome Data:**

> "Microbiome relative abundance data were analyzed using three normalization approaches: centered log-ratio (CLR) transformation, log-normal transformation, and raw abundances with pseudocount addition (1e-6) for zero values."

**Software:**

> "Analyses were conducted using R with the survey package for complex survey design, broom package for standardized model output extraction, and custom scripts for batch processing and FLEXIBLE resume functionality."

---

#  **Production Deployment Summary**

### **Core Workflow:**

1. **Setup:** `setup_nhanes_environment.sh`
2. **Validate:** `debug_all_pipelines.R`
3. **Analyze:** `run_all_was_analyses_flexible.sh`
4. **Aggregate:** `aggregate_was_results.R`

### **Key Advantages:**

- **FLEXIBLE Resume:** Only processes missing results
- **Corrected Resources:** Based on actual computational requirements
- **Survey Design:** Proper NHANES complex survey implementation
- **Multiple Comparisons:** FDR correction by dependent variable
- **Production Ready:** Comprehensive error handling and logging

### **Main Results:**

- Use `*_tidied_complete.rds` files
- Significance threshold: `p.value.fdr < 0.05`
- Survey-adjusted standard errors included
- Complete metadata for reproducibility
