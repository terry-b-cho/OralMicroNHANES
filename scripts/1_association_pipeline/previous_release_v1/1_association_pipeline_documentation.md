# WAS Analysis Pipeline - Complete File Status and Usage Guide

## **ESSENTIAL FILES - USE THESE**

### **Core Production Files:**
1. **`setup_nhanes_environment.sh`** ✅ **CURRENT**
   - Primary environment setup script
   - **Use this first** for environment creation

2. **`test_survey_implementation.R`** ✅ **CURRENT** 
   - Validation and testing script
   - **Run this after setup** to verify installation

3. **`run_complete_was_analysis.sh`** ✅ **CURRENT** ⭐ **INDIVIDUAL REGRESSION PIPELINE**
   - Individual regression analysis pipeline (NO aggregation)
   - **This runs regressions only** - requires separate aggregation step

4. **`run_aggregation_only.sh`** ✅ **CURRENT** ⭐ **AGGREGATION PIPELINE**
   - Result aggregation with multiple comparisons correction
   - **Run this AFTER individual regressions complete**

5. **`aggregate_was_results.R`** ✅ **CURRENT**
   - Core aggregation script with FDR/Bonferroni correction
   - **Auto-called by run_aggregation_only.sh** (don't run manually)

6. **`universal_was_analysis_debugging_report_information.R`** ✅ **CURRENT**
   - Core analysis script with survey design
   - **Auto-called by regression pipeline** (don't run manually)

7. **`run_all_was_analyses.sh`** ✅ **CURRENT** ⭐ **BATCH REGRESSION SUBMISSION**
   - Submits all 18 individual regression analyses at once (6 types × 3 normalizations)
   - **Use this to submit all regressions** - requires separate aggregation step

8. **`submit_smart_resume.sh`** ✅ **OPTIONAL DIAGNOSTIC**
   - Identifies missing individual regressions and resumes them
   - **Optional diagnostic tool** - does NOT perform aggregation

---

## 🔄 **BACKUP FILES - USE IF NEEDED**

9. **`setup_nhanes_environment_simple.sh`** 🔄 **BACKUP**
   - Alternative environment setup method
   - **Use only if main setup fails**

10. **`run_was_analysis_debug.sh`** 🔄 **BACKUP**
   - Manual individual job submission
   - **Use only for advanced manual control**

---

## 🔄 **EXECUTION FLOW: NEW TWO-STEP PIPELINE**

### **STEP 1: Individual Regressions**

When you run:
```bash
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS clr
```

Here's **exactly** what scripts get executed:

### **STEP-BY-STEP EXECUTION:**

#### **Step 1: Individual Analysis Submission**
```bash
# run_complete_was_analysis.sh calls:
→ run_was_analysis_debug.sh 1_demoWAS clr "" 64G 4:00:00
```

#### **Step 2: SLURM Job Creation (for each dependent variable)**
```bash
# run_was_analysis_debug.sh creates SLURM jobs that call:
→ universal_was_analysis_debugging_report_information.R
```

#### **Step 3: Job Monitoring**
```bash
# run_complete_was_analysis.sh waits for all jobs to complete
→ Monitors with squeue until all jobs finish
```

#### **Step 4: Next Step Guidance**
```bash
# run_complete_was_analysis.sh provides instructions for:
→ Running run_aggregation_only.sh for final results
```

### **STEP 2: Result Aggregation (Separate Step)**

When you run:
```bash
./scripts/1_association_pipeline/run_aggregation_only.sh
```

This script:
```bash
# Identifies analyses needing aggregation and calls:
→ aggregate_was_results.R for each analysis requiring aggregation
```

### **DETAILED BREAKDOWN:**

#### **What `run_complete_was_analysis.sh` Does (NEW):**
1. **Calls** `run_was_analysis_debug.sh` to submit individual jobs
2. **Monitors** job completion automatically  
3. **Provides guidance** for next steps (NO automatic aggregation)
4. **Displays** summary and next step instructions

#### **What `run_was_analysis_debug.sh` Does:**
1. **Reads** schema file to find all dependent variables
2. **Creates** individual SLURM job scripts for each dependent variable
3. **Submits** each job to the cluster
4. **Each job calls** `universal_was_analysis_debugging_report_information.R`

#### **What `universal_was_analysis_debugging_report_information.R` Does:**
1. **Runs** survey-aware regression analysis for one dependent variable
2. **Fits** models for all independent variables associated with that dependent variable
3. **Saves** individual results as `.rds` files

#### **What `run_aggregation_only.sh` Does (NEW):**
1. **Scans** all analyses to identify which need aggregation
2. **Calls** `aggregate_was_results.R` for analyses with individual results but missing aggregated files
3. **Skips** analyses that are already aggregated or have no individual results
4. **Provides** summary of aggregation operations

#### **What `aggregate_was_results.R` Does:**
1. **Reads** all individual `.rds` files
2. **Combines** them into analysis-type-specific files
3. **Applies** FDR and Bonferroni multiple comparisons correction
4. **Creates** final aggregated results

### **CONCRETE EXAMPLE:**

For `1_demoWAS clr` analysis:

```bash
# 1. Main pipeline starts (using CLR-transformed microbiome OTU tables)
run_complete_was_analysis.sh 1_demoWAS clr

# 2. Submits individual jobs (one per microbiome OTU dependent variable)
run_was_analysis_debug.sh 1_demoWAS clr
  ├── Job 1: RSV_genus1002_relative (CLR-transformed) → universal_was_analysis_debugging_report_information.R
  ├── Job 2: RSV_genus1076_relative (CLR-transformed) → universal_was_analysis_debugging_report_information.R  
  ├── Job 3: RSV_genus396_relative (CLR-transformed) → universal_was_analysis_debugging_report_information.R
  └── ... (hundreds more microbiome OTU jobs)

# 3. Waits for all jobs to complete

# 4. Aggregates results
aggregate_was_results.R 1_demoWAS clr results
  ├── Reads: RSV_genus1002_relative.rds
  ├── Reads: RSV_genus1076_relative.rds
  ├── Reads: RSV_genus396_relative.rds
  └── Creates: 1_demoWAS_clr_tidied_complete.rds (with FDR correction)
```

**Key Point**: The `clr` specification tells the pipeline to use the **CLR-transformed RSV tables** (pre-processed microbiome OTU abundance data) rather than raw or log-normal transformed tables.

### **KEY POINTS:**

1. **`universal_was_analysis_debugging_report_information.R`** is the **core analysis script** that does the actual statistical modeling

2. **`aggregate_was_results.R`** is the **post-processing script** that combines results and adds multiple comparisons correction

3. **Both are auto-called** - you don't run them directly

4. **The main pipeline orchestrates everything** - job submission, monitoring, aggregation, and validation

---

## **QUICK START COMMANDS - NEW TWO-STEP WORKFLOW**

### **Step 1: Environment Setup (One-time)**
```bash
# Start interactive session
srun --pty -p interactive -t 12:00:00 --mem=32G bash

# Setup environment
bash scripts/1_association_pipeline/setup_nhanes_environment.sh
```

### **Step 2: Validation (Recommended)**
```bash
# Load modules and test
module load gcc/9.2.0
module load miniconda3/23.1.0
source activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Run validation tests
Rscript scripts/1_association_pipeline/test_survey_implementation.R
```

### **Step 3A: Run Individual Regressions (Choose One Option)**

#### **Option A: Batch Submission (All 18 Regression Analyses) ⭐ RECOMMENDED**
```bash
# Exit interactive session
exit

# Submit ALL individual regression analyses (6 types × 3 normalizations = 18 total)
./scripts/1_association_pipeline/run_all_was_analyses.sh

# For test mode (recommended first time)
./scripts/1_association_pipeline/run_all_was_analyses.sh test

# With custom resources
./scripts/1_association_pipeline/run_all_was_analyses.sh "" 120G 12:00:00
```

#### **Option B: Individual Regression Analysis**
```bash
# Run single regression analysis
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS clr

# For test mode
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS clr test
```

### **Step 3B: Result Aggregation (Required After Regressions Complete)**

#### **After ALL regression analyses complete:**
```bash
# Aggregate all results with FDR/Bonferroni correction
./scripts/1_association_pipeline/run_aggregation_only.sh

# Or aggregate specific analysis only
./scripts/1_association_pipeline/run_aggregation_only.sh 1_demoWAS clr
```

### **Optional: Diagnostic and Resume (If Some Analyses Fail)**
```bash
# Check status and resume missing individual regressions
sbatch scripts/1_association_pipeline/submit_smart_resume.sh
```

**Purpose of `submit_smart_resume.sh`:** Runs the `submit_smart_resume.R` unattended on the cluster.

**What `submit_smart_resume.R` it does:**
- **Skips complete analyses** (saves hours)
- **Runs aggregation** for analyses with all variables but missing final step
- **Submits targeted jobs** for only missing dependent variables (not entire analyses)
- **Flags completely wrong analyses** for full re-runs

---

## **Analysis Types and Commands**

### **Analysis Types Overview:**

| Pipeline | Description | Model Type | Dependent Variable | Independent Variable |
|----------|-------------|------------|-------------------|---------------------|
| **1_demoWAS** | How demographics affect microbiome | Linear | **Microbiome OTU abundance (RSV)** | Demographics |
| **2_oradWAS** | How microbiome affects oral health | Logistic | Oral health outcomes | **Microbiome OTU abundance (RSV)** |
| **3_exWAS** | How exposures affect microbiome | Linear | **Microbiome OTU abundance (RSV)** | Exposure variables |
| **4_pheWAS** | How microbiome affects phenotypes | Linear | Phenotype measures | **Microbiome OTU abundance (RSV)** |
| **5_outWAS** | How microbiome affects disease outcomes | Logistic | Disease outcomes | **Microbiome OTU abundance (RSV)** |
| **6_zimWAS** | How microbiome affects lab measurements | Linear | Lab measurements | **Microbiome OTU abundance (RSV)** |

**Note**: RSV = Ribosomal Sequence Variant tables containing NHANES oral microbiome relative abundance OTU data

### **All Available Analyses:**
```bash
# 1_demoWAS - Demographics → Microbiome
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS clr
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS lognorm
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS none

# 2_oradWAS - Microbiome → Oral Health
./scripts/1_association_pipeline/run_complete_was_analysis.sh 2_oradWAS clr
./scripts/1_association_pipeline/run_complete_was_analysis.sh 2_oradWAS lognorm
./scripts/1_association_pipeline/run_complete_was_analysis.sh 2_oradWAS none

# 3_exWAS - Exposures → Microbiome
./scripts/1_association_pipeline/run_complete_was_analysis.sh 3_exWAS clr
./scripts/1_association_pipeline/run_complete_was_analysis.sh 3_exWAS lognorm
./scripts/1_association_pipeline/run_complete_was_analysis.sh 3_exWAS none

# 4_pheWAS - Microbiome → Phenotypes
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS clr
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS lognorm
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS none

# 5_outWAS - Microbiome → Disease Outcomes
./scripts/1_association_pipeline/run_complete_was_analysis.sh 5_outWAS clr
./scripts/1_association_pipeline/run_complete_was_analysis.sh 5_outWAS lognorm
./scripts/1_association_pipeline/run_complete_was_analysis.sh 5_outWAS none

# 6_zimWAS - Microbiome → Lab Measurements
./scripts/1_association_pipeline/run_complete_was_analysis.sh 6_zimWAS clr
./scripts/1_association_pipeline/run_complete_was_analysis.sh 6_zimWAS lognorm
./scripts/1_association_pipeline/run_complete_was_analysis.sh 6_zimWAS none
```

### **Normalization Methods:**
**IMPORTANT**: These refer to **pre-processed microbiome OTU abundance tables** created in previous pipeline steps. The regression analysis uses different pre-normalized tables based on your selection:

- **`clr`** - Uses **CLR-transformed RSV tables** (Centered log-ratio transformation of microbiome OTU relative abundance)
- **`lognorm`** - Uses **log-normal transformed RSV tables** (Log-normal transformation of microbiome OTU relative abundance)  
- **`none`** - Uses **raw RSV tables** (Untransformed microbiome OTU relative abundance)

**Note**: RSV = Ribosomal Sequence Variant tables containing NHANES oral microbiome relative abundance OTU data. The normalization was applied to these microbiome variables in earlier data processing steps, not during the regression analysis.

### **Database Table Structure:**
The pipeline accesses different pre-processed microbiome tables based on normalization method:

- **`clr`** → Uses tables like `DADA2RSV_GENUS_RELATIVE_F_clr`, `DADA2RSV_GENUS_RELATIVE_G_clr`
- **`lognorm`** → Uses tables like `DADA2RSV_GENUS_RELATIVE_F_lognorm`, `DADA2RSV_GENUS_RELATIVE_G_lognorm`  
- **`none`** → Uses tables like `DADA2RSV_GENUS_RELATIVE_F_none`, `DADA2RSV_GENUS_RELATIVE_G_none`

Each table contains the same microbiome OTU variables (e.g., `RSV_genus1002_relative`) but with different transformations applied.

---

## 📁 **Output Files (What You Get)**

### **Main Result File Naming Patterns by Analysis Type**

| Analysis Type | Microbe Role | File Naming Pattern | Example Files |
|---------------|--------------|---------------------|---------------|
| **1_demoWAS** | Dependent | `RSV_genus*.rds` | `RSV_genus1002_relative.rds` |
| **3_exWAS** | Dependent | `RSV_genus*.rds` | `RSV_genus10_relative.rds` |
| **2_oradWAS** | Independent | `<dep_var>.rds` | `DENTURE_OHAROCDE.rds` |
| **4_pheWAS** | Independent | `<dep_var>.rds` | `BMXBMI.rds` |
| **5_outWAS** | Independent | `<dep_var>.rds` | `DIQ010.rds` |
| **6_zimWAS** | Independent | `<dep_var>.rds` | `LBXGLU.rds` |

### **Main Aggregated Result Files:**
```
results/<analysis_type>_out/result_<normalization>/
├── <analysis_type>_<normalization>_tidied_complete.rds    ⭐ MAIN RESULTS
├── <analysis_type>_<normalization>_glanced_complete.rds   # Model stats
├── <analysis_type>_<normalization>_rsq_complete.rds      # R-squared
└── <analysis_type>_<normalization>_aggregation_summary.txt # Summary
```

### **Key Features of Main Results:**
- ✅ **Simple tibble** (not nested lists)
- ✅ **Multiple comparisons correction** (FDR and Bonferroni)
- ✅ **Survey-adjusted standard errors**
- ✅ **All coefficient results** in one file

### **Main Results Structure:**
```r
# Load main results
results <- readRDS("results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds")

# Key columns:
# term                - Model term (focus on "indep_var" for main effects)
# estimate            - Effect size (coefficient)
# std.error           - Survey-adjusted standard error  
# p.value             - Original p-value
# p.value.fdr         - FDR-corrected p-value ⭐ USE FOR REPORTING
# p.value.bonferroni  - Bonferroni-corrected p-value
# phenotype           - Dependent variable name
# exposure            - Independent variable name
# n_obs               - Sample size
# series              - NHANES cycle (F, G)
```

---

## **Multiple Comparisons Correction**

### **The Problem:**
In WAS analyses, thousands of statistical tests are performed. Without correction:
- **High false positive rate** (Type I error inflation)
- **Many spurious associations** would be reported
- **Unreliable scientific conclusions**

### **Our Solution:**
The pipeline automatically applies **two correction methods**:

1. **FDR (False Discovery Rate)** - Benjamini-Hochberg procedure
   - ✅ **Recommended for reporting**
   - Controls expected proportion of false discoveries
   - Less conservative than Bonferroni
   - **Use `p.value.fdr < 0.05` for significance**

2. **Bonferroni Correction** - Family-wise error rate control
   - Very conservative approach
   - Controls probability of any false discovery
   - Use for very stringent significance thresholds

### **Implementation Details:**
- ✅ Corrections applied only to **main effect terms** (`indep_var`)
- ✅ Intercept and covariate terms retain original p-values
- ✅ Corrections calculated **within each analysis type**
- ✅ Both methods included in results for flexibility

---

## **Statistical Methods**

### **Survey Design Implementation:**
- **Complex Survey Design**: Proper NHANES stratification (`SDMVSTRA`) and clustering (`SDMVPSU`)
- **Weight Scaling**: Correctly scales survey weights across multiple 2-year cycles
- **Survey Models**: Uses `svyglm()` from R `survey` package (NOT standard `lm()`/`glm()`)

### **Model Types:**
- **Linear Regression**: Continuous dependent variables (gaussian family)
- **Logistic Regression**: Binary dependent variables (quasibinomial family)

### **Key Features:**
- ✅ **Complete case analysis** after weight scaling
- ✅ **Proper variance estimation** accounting for survey design  
- ✅ **Binary outcome validation** for logistic regression
- ✅ **Comprehensive error handling** and logging
- ✅ **Automatic resource allocation** based on analysis size

---

## 🔧 **Troubleshooting Quick Reference**

### **Environment Issues:**
```bash
# If main setup fails, try backup:
bash scripts/1_association_pipeline/setup_nhanes_environment_simple.sh

# If still failing, manual R install:
R -e "install.packages(c('survey', 'broom', 'dplyr'), repos='https://cloud.r-project.org/')"
```

### **Job Monitoring:**
```bash
# Check running jobs
squeue -u $USER

# Check specific analysis
squeue -u $USER --name="1_demoWAS_clr*"

# View job logs
cat results/1_demoWAS_out/logs/slurm_*.err
```

### **Resource Adjustment:**
```bash
# For large analyses, increase resources:
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS clr "" 120G 24:00:00
#                                                                    ↑     ↑
#                                                                 memory  time
```

### **Common Issues and Solutions:**

#### **1. Environment Setup Issues**
```bash
# Problem: Package installation fails
# Solution: Try alternative setup methods

# Method 1 failed? Try Method 2:
bash scripts/1_association_pipeline/setup_nhanes_environment_simple.sh

# Still failing? Check available packages:
conda search r-survey
conda search r-broom

# Last resort: Install via R directly
R -e "install.packages(c('survey', 'broom', 'dplyr'), repos='https://cloud.r-project.org/')"
```

#### **2. Job Failures**
```bash
# Problem: Jobs fail due to memory
# Solution: Increase memory allocation
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS clr "" 120G 12:00:00

# Problem: Jobs fail due to time limits  
# Solution: Increase time allocation
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS clr "" 120G 24:00:00

# Problem: Individual job errors
# Solution: Check specific job logs
cat results/4_pheWAS_out/logs/slurm_<job_id>.err
```

#### **3. Module Loading Issues**
```bash
# CRITICAL: Load modules in this exact order
module purge
module load gcc/9.2.0      # Must be first
module load miniconda3/23.1.0  # Must be last
```

#### **4. Aggregation Issues**
```bash
# Problem: Aggregation fails
# Solution: Check individual results exist
ls results/1_demoWAS_out/result_clr/*.rds

# Manual aggregation if needed
Rscript scripts/1_association_pipeline/aggregate_was_results.R 1_demoWAS clr results
```

---

## ✅ **Success Indicators**

### **Setup Successful When:**
- ✅ Environment created without errors
- ✅ All R packages load successfully  
- ✅ Validation tests pass
- ✅ Test analysis completes

### **Analysis Successful When:**
- ✅ All jobs complete (check with `squeue`)
- ✅ Aggregated files created
- ✅ Summary shows expected number of tests
- ✅ FDR/Bonferroni multiple comparisons correction applied
- ✅ Survey-adjusted standard errors reported
- ✅ Results contain p.value.fdr column

---

## 📝 **Citation and Methods for Publications**

### **Methods Section Text:**

**Survey Design:**
> "Analyses accounted for the complex survey design of NHANES using appropriate stratification, clustering, and survey weights scaled across multiple cycles. All models were fitted using survey-weighted generalized linear models (svyglm) from the R survey package."

**Multiple Comparisons:**
> "P-values were adjusted for multiple comparisons using the Benjamini-Hochberg false discovery rate (FDR) procedure. Associations with FDR-adjusted p-values < 0.05 were considered statistically significant."

**Software:**
> "Analyses were conducted using R with the survey package for complex survey design and the broom package for standardized model output extraction."

---

## **Bottom Line - NEW WORKFLOW**

### **For Most Users (Two-Step Process):**
1. **Setup**: `bash scripts/1_association_pipeline/setup_nhanes_environment.sh`
2. **Test**: `Rscript scripts/1_association_pipeline/test_survey_implementation.R`  
3. **Run Regressions**: `./scripts/1_association_pipeline/run_all_was_analyses.sh`
4. **Wait for Completion**: Monitor with `squeue -u $USER`
5. **Aggregate Results**: `./scripts/1_association_pipeline/run_aggregation_only.sh`

### **Main Results Files (After Aggregation):**
- `results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds`
- Use `p.value.fdr < 0.05` for significance
- Survey-adjusted standard errors included

### **If Problems:**
1. Check log files in `results/<analysis_type>_out/logs/`
2. Try backup setup script if environment fails
3. Use test mode first: add `test` to command
4. Increase resources if jobs fail
5. Use `submit_smart_resume.sh` to diagnose and fix missing regressions

### **Key Changes:**
- **Regressions** and **aggregation** are now separate steps
- More control and better error recovery
- Clear workflow with explicit next steps

**This pipeline is production-ready for NHANES oral microbiome X variable Wide Association Studies with proper survey design and multiple comparisons correction!** 