# NHANES Oral Microbiome Data Transformation and Schema Structure Preprocessing Guide

## Overview**

This directory contains the complete preprocessing pipeline for NHANES oral microbiome data analysis. The workflow consists of two main stages:

1. **Data Transformation**: Apply normalization methods to microbiome OTU abundance data
2. **Data Derivation and Completion**: filling in missing data and deriving neccessary variables for analyses
3. **Schema Structure Creation**: Generate mapping files for downstream WAS analyses

## 📁 **Files in This Directory**

### **Core Scripts:**

1. **`nhanes_omp_transformation.R`** - Transforms microbiome abundance data (CLR, log-normal, none)
2. **`nhanes_db_filling_missing_data.R`** - Fills missing metadata and derives additional variables
3. **`ss_file_create.R`** - Creates schema structure files for WAS analyses

### **SLURM Submission Scripts:**

4. **`run_transformation.sh`** - Submits transformation jobs to SLURM
5. **`run_ss_file_create.sh`** - Submits schema structure creation jobs to SLURM

### **Debug Tool Available:**

6. **`debug_ss.sh`** - debug script to test schema structure creation manually

---

## **Complete Workflow - Step by Step**

### **Prerequisites**

- Access to slurm HPC cluster
- R ≥ 4.2.1 with required packages
- Input database: `data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite`
- Configuration files in `configs/` directory: `# These are list of dependent variables to exploit for ssfile creation per analysis specification`

## **Step 1: Data Transformation** **`nhanes_omp_transformation.R`**

#### **What This Does:**

Processes NHANES oral microbiome data with **NO prevalence pre-filtering** (optionally setting this is possible). Creates four normalized versions of each microbiome table while preserving ALL taxa for regression-time filtering decisions.

### **Design Decisions:**

1. **NO Prevalence Pre-filtering**:

   - Keeps ALL taxa regardless of prevalence
   - Filtering decisions deferred to regression-time
   - Preserves analytical flexibility
2. **Universal Pseudocount Strategy** (ε = 0.001):

   ```r
   EPS <- 1e-3  # Universal ε = 0.001
   ```

   - **Fixed pseudocount** across all taxonomic levels for comparability
   - **Applied only to COUNT matrices** before transformation
   - **Prevents log(0) without distorting relative scale**
   - **No pseudocount for "none" and Hellinger** (zeros allowed)
3. **Metadata Synchronization**:

   - Automatically updates `variable_names_epcf` table
   - Creates metadata entries for each transformed table
   - Enables downstream code to resolve tables automatically

#### **Statistical Transformations Applied:**

**1. "none" Transformation:**

```r
# Applied to RELATIVE tables
T_ij = P_ij  (where P_ij = C_ij / n_i)
```

- **Raw relative abundances** (no mathematical transformation)
- **Input**: RELATIVE tables 
- **Preserves original compositional structure**
- **Values remain in [0,1] interval**

**2. Hellinger Transformation:**

```r
# Applied to RELATIVE tables  
T_ij = √P_ij
```

- **Square-root transformation** of relative abundances
- **Input**: RELATIVE tables
- **Reduces impact of dominant taxa**
- **Euclidean distance becomes meaningful**

**3. CLR (Centered Log-Ratio) Transformation:**

```r
# Applied to COUNT tables
# Step 1: Add universal pseudocount
C_pseudo = C_ij + ε  (where ε = 0.001)

# Step 2: Calculate geometric mean per sample
g_i = exp[(1/D) * Σⱼ ln(C_ij + ε)]

# Step 3: CLR transformation
T_ij = ln[(C_ij + ε) / g_i]

where:
- D = number of taxa (columns)
- g_i = geometric mean of (counts + ε) for sample i
```

- **Input**: COUNT tables (not relative abundances)
- **Addresses compositional data constraints** (unit-sum constraint)
- **Centers around zero** (mean of log-ratios = 0)
- **Preserves relative relationships** between taxa
- **Natural log base** for CLR geometry

**4. Log-Normal Transformation:**

```r
# Applied to COUNT tables
# Step 1: Add universal pseudocount
C_pseudo = C_ij + ε  (where ε = 0.001)

# Step 2: Calculate library sizes and mean depth
n_i = rowSums(C_ij)  (raw library sizes for mean calculation)
n̄ = mean(n_i)  (mean library size WITHOUT εD)
lib_size = n_i + ε*D  (library size + pseudocount correction for denominator)

# Step 3: Log-normal transformation  
T_ij = log₁₀[(C_ij + ε)/(n_i + ε*D) × n̄]
```

- **Input**: COUNT tables (not relative abundances)
- **Base-10 logarithm** for interpretability
- **Variance stabilization** for right-skewed abundance distributions
- **Library size normalization** with mean depth scaling
- **Approximates normality** for downstream linear regression

**Universal pseudocount ε = 0.001:**

```r
EPS <- 1e-3  # Fixed across all taxonomic levels
```

**Pseudocount characteristics:**

- **Fixed value**: ε = 0.001 for all transformations
- **Applied only to COUNT matrices** before any log operation
- **log₁₀(0.001) = -3** for log-normal transformation
- **ln(0.001) ≈ -6.9** for CLR transformation

### **Range Description:**


| Transformation | Range                | Distribution       | Input Tables | Use Case                        |
| -------------- | -------------------- | ------------------ | ------------ | ------------------------------- |
| **none**       | [0, 1]               | Right/left-skewed  | RELATIVE     | Raw abundance analysis          |
| **hellinger**  | [0, 1]               | Less skewed        | RELATIVE     | Distance-based analysis         |
| **clr**        | (-∞, +∞)            | Symmetric around 0 | COUNT        | Compositional data analysis     |
| **lognorm**    | (-∞, +∞)            | More normal        | COUNT        | Linear regression               |

**Range Details:**
- **lognorm**: Range depends on count data and library sizes, typically extends to negative values with fixed lower bound influenced by ε = 0.001
- **clr**: Symmetric around zero due to centering; range depends on count variation within samples
- **hellinger**: Bounded [0,1] since it's square-root of relative abundances
- **none**: Bounded [0,1] as original relative abundances

## **Command:**

```bash
# Submit transformation job
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh
```

### **Manual Execution:**

```bash
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R \
  --in_db  data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite
```

### **Input Tables Processed:**

**Taxonomic Levels × Cycles:**

- `DADA2RSV_GENUS_RELATIVE_F/G` (genus level, cycles F & G)
- `DADA2RSV_FAMILY_RELATIVE_F/G` (family level, cycles F & G)
- `DADA2RSV_ORDER_RELATIVE_F/G` (order level, cycles F & G)
- `DADA2RSV_CLASS_RELATIVE_F/G` (class level, cycles F & G)
- `DADA2RSV_PHYLUM_RELATIVE_F/G` (phylum level, cycles F & G)

**Total**: 10 source tables → 40 transformed tables (4 transformations each)

### **Output:**

**Database**: `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite`

**New Tables** (40 total):

```bash
All output tables use RELATIVE naming convention:

From RELATIVE input tables:
DADA2RSV_GENUS_RELATIVE_F → {
  DADA2RSV_GENUS_RELATIVE_F_none,
  DADA2RSV_GENUS_RELATIVE_F_hellinger
}

From COUNT input tables (but output as RELATIVE names):
DADA2RSV_GENUS_COUNT_F → {
  DADA2RSV_GENUS_RELATIVE_F_clr,      ← Uses RELATIVE naming
  DADA2RSV_GENUS_RELATIVE_F_lognorm   ← Uses RELATIVE naming
}
[... repeated for all 10 taxonomic levels]
```

**Metadata Updates**: Corresponding entries added to `variable_names_epcf` with transformation-specific descriptions.

- **Database**: `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite`
- **New Tables**: 40 transformed tables (10 taxonomic levels × 4 transformations)
- **Input Source**: RELATIVE tables for none/hellinger, COUNT tables for clr/lognorm
- **Metadata**: Updated `variable_names_epcf` with transformation-specific entries
- **Taxa Count**: All original taxa preserved (no filtering applied)
- **Pseudocount**: ε = 0.001 applied universally to COUNT matrices only


Expected output:
```bash
(base) Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R

⤷  PIPELINE - NO PRE-FILTERING
⤷ in_db:  data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite
⤷ out_db: data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite
✓ Copied input database to: data/00_nhanes_omp_transformed_db/nhanes_oral_transf
ormed.sqlite
✓ DADA2RSV_GENUS_RELATIVE_F → none / hellinger
✓ DADA2RSV_GENUS_RELATIVE_G → none / hellinger
✓ DADA2RSV_FAMILY_RELATIVE_F → none / hellinger
✓ DADA2RSV_FAMILY_RELATIVE_G → none / hellinger
✓ DADA2RSV_ORDER_RELATIVE_F → none / hellinger
✓ DADA2RSV_ORDER_RELATIVE_G → none / hellinger
✓ DADA2RSV_CLASS_RELATIVE_F → none / hellinger
✓ DADA2RSV_CLASS_RELATIVE_G → none / hellinger
✓ DADA2RSV_PHYLUM_RELATIVE_F → none / hellinger
✓ DADA2RSV_PHYLUM_RELATIVE_G → none / hellinger
✓ DADA2RSV_GENUS_COUNT_F → clr / lognorm
✓ DADA2RSV_GENUS_COUNT_G → clr / lognorm
✓ DADA2RSV_FAMILY_COUNT_F → clr / lognorm
✓ DADA2RSV_FAMILY_COUNT_G → clr / lognorm
✓ DADA2RSV_ORDER_COUNT_F → clr / lognorm
✓ DADA2RSV_ORDER_COUNT_G → clr / lognorm
✓ DADA2RSV_CLASS_COUNT_F → clr / lognorm
✓ DADA2RSV_CLASS_COUNT_G → clr / lognorm
✓ DADA2RSV_PHYLUM_COUNT_F → clr / lognorm
✓ DADA2RSV_PHYLUM_COUNT_G → clr / lognorm

Transformation completed with four scenarios:
   • none (relative)  | hellinger  → from *_RELATIVE_* tables
   • clr  | lognorm   → from *_COUNT_*    tables
   • ε = 0.001 applied only to count data
   • Output DB: data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite
```

## **Step 2: Missing Data Completion and Variable Derivation** **`nhanes_db_filling_missing_data.R`**

### **What This Does:**

Completes the transformed database by deriving additional variables and filling missing metadata entries to ensure all variables used in downstream analyses are properly documented and accessible.

### **Core Operations:**

#### **1. Gender Variable Derivation:**

```r
# In-place SQL transformation for both cycles F and G
ALTER TABLE DEMO_F ADD COLUMN RIAGENDR_01 INTEGER
UPDATE DEMO_F SET RIAGENDR_01 = RIAGENDR - 1

ALTER TABLE DEMO_G ADD COLUMN RIAGENDR_01 INTEGER  
UPDATE DEMO_G SET RIAGENDR_01 = RIAGENDR - 1

# Result: 0 = Male, 1 = Female (from original 1 = Male, 2 = Female)
```

#### **2. Comprehensive Variable Annotation System:**

The script includes manual annotations for **40+ variables** across multiple categories:

**Demographic Derived Variables:**

```r
AGE_SQUARED          → "Age squared (RIDAGEYR²)"
BORN_INUSA           → "Born in USA (DMDBORN4=1)"
RIAGENDR_01          → "Gender (0=Male,1=Female)"
```

**Education Indicators:**

```r
EDUCATION_LESS9      → "Edu: <9th grade (DMDEDUC2=1)"
EDUCATION_9_11       → "Edu: 9–11th grade (DMDEDUC2=2)"  
EDUCATION_AA         → "Edu: AA/some college (DMDEDUC2=4)"
EDUCATION_COLLEGEGRAD → "Edu: college graduate or above (DMDEDUC2=5)"
```

**Ethnicity Indicators:**

```r
ETHNICITY_MEXICAN           → "Mexican American (RIDRETH1=1)"
ETHNICITY_OTHERHISPANIC     → "Other Hispanic (RIDRETH1=2)"
ETHNICITY_OTHER             → "Other race/multi-racial (RIDRETH1=5)"
ETHNICITY_NONHISPANICBLACK  → "Non-Hispanic Black (RIDRETH1=4)"
```

**Laboratory Biomarkers:**

```r
# Prostate markers
LBXP1, LBXP2, LBDP3, LBXPS4 → Various PSA measurements

# Folate measures  
LBDRFO, LBXFOLSI, LBDFOT    → RBC and serum folate levels

# Other biomarkers
LBXSCK, LBXMMASI, URX1DC    → Creatine kinase, arsenic, cysteine
```

**Health Outcomes:**

```r
# Respiratory conditions
ASTHMA, BRONCHITIS, EMPHYSEMA → Self-reported respiratory diseases

# Cardiovascular conditions  
ANGINA, HEART_FAILURE, CHD, CVD → Self-reported cardiovascular diseases

# Cancer history
CANCER_BREAST, CANCER_LUNG, CANCER_PROSTATE → Cancer history variables
```

#### **3. Metadata Completion:**

```r
# Step 1: Variable→Table mapping via field inspection
mapping_algorithm = {
  for each variable v in annotations:
    find all tables T where v ∈ dbListFields(T)
    create mapping: variable → table(s)
}

# Step 2: Anti-join with existing metadata  
new_entries = mapping ∖ existing_metadata_rows

# Step 3: Inherit table-level metadata template
template_fields = {
  Data.File.Name, Data.File.Description,
  Begin.Year, EndYear, Component, Use.Constraints
}

# Step 4: Assemble and append new metadata rows
new_metadata = new_entries + annotations + template_inheritance
```

#### **4. Database Transaction Safety:**

```r
# Atomic metadata updates with rollback protection
dbBegin(con)
tryCatch({
  old_metadata ← existing variable_names_epcf
  new_complete ← bind_rows(old_metadata, new_metadata)  
  dbWriteTable(con, "variable_names_epcf", new_complete, overwrite=TRUE)
  dbCommit(con)
}, error = function(e) {
  dbRollback(con)  # Ensures database integrity
  stop("Error appending metadata: ", e$message)
})
```

### **Variables Added to Metadata:**

The script automatically detects and documents **40+ variables** across:

- **Demographic derivatives** (age, gender, education, ethnicity)
- **Laboratory biomarkers** (PSA, folate, creatine kinase, etc.)
- **Anthropometric measures** (skinfold measurements)
- **Health outcomes** (respiratory, cardiovascular, cancer history)
- **Exposure biomarkers** (arsenic, cysteine derivatives)

### **Command:**

```bash
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
  --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
```

### **Input/Output:**

**Input Database**: `nhanes_oral_transformed.sqlite`

- Transformed microbiome tables
- Incomplete metadata coverage

**Output Database**: `nhanes_oral_transformed_complete.sqlite`

- **All original content preserved** (file copy operation)
- **RIAGENDR_01 derived** in both DEMO_F and DEMO_G tables
- **Complete metadata coverage** for all analysis variables
- **Automatic variable→table mapping** for downstream schema generation


### **Step 3: Schema Structure File Creation**

#### **What This Does:**

- Creates mapping files for each WAS analysis type and transformation method
- Maps microbiome OTU variables to NHANES variables
- Calculates sample sizes for each variable pair
- Generates 24 schema structure files (6 analysis types × 4 transformations)

#### **Analysis Types and Variable Roles:**


| Pipeline      | Description                    | Microbiome Role | Config File                  |
| --------------- | -------------------------------- | ----------------- | ------------------------------ |
| **1_demoWAS** | Demographics → Microbiome     | Dependent       | `configs/1_demoWAS_vars.txt` |
| **2_oradWAS** | Microbiome → Oral Health      | Independent     | `configs/2_oradWAS_vars.txt` |
| **3_exWAS**   | Exposures → Microbiome        | Dependent       | `configs/3_exWAS_vars.txt`   |
| **4_pheWAS**  | Microbiome → Phenotypes       | Independent     | `configs/4_pheWAS_vars.txt`  |
| **5_outWAS**  | Microbiome → Disease Outcomes | Independent     | `configs/5_outWAS_vars.txt`  |
| **6_zimWAS**  | Microbiome → Lab Measurements | Independent     | `configs/6_zimWAS_vars.txt`  |

#### **Command:**

```bash
# Submit all schema structure creation jobs
bash scripts/0_transform_n_preprocess_ssfiles/run_ss_file_create.sh
```

#### **Manual Execution Example:**

```bash
Rscript scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R \
  --db        data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
  --otu_F     DADA2RSV_GENUS_RELATIVE_F_clr \
  --otu_G     DADA2RSV_GENUS_RELATIVE_G_clr \
  --vars_file configs/2_oradWAS_vars.txt \
  --otu_role  indep \
  --pipeline  2_oradWAS \
  --transform clr \
  --out_dir   results/0_ss_files
```

#### **Output:**

- **Schema Files**: `results/0_ss_files/<pipeline>_<transform>_schema_structure.csv`
- **24 Total Files**: e.g., `1_demoWAS_clr_schema_structure.csv`, `2_oradWAS_hellinger_schema_structure.csv`

---

## **Output File Structure**

### **Schema Structure File Format:**

Each CSV contains the following columns:


| Column        | Description                                           |
| --------------- | ------------------------------------------------------- |
| `dep_var`     | Dependent variable name                               |
| `dep_table`   | NHANES table where dependent variable came from       |
| `indep_var`   | Independent variable name                             |
| `indep_table` | NHANES/OTU table where independent variable came from |
| `n`           | Number of non-missing pairs                           |
| `cycle`       | NHANES cycle (F or G)                                 |

### **Example Schema Structure Content:**

```csv
dep_var,dep_table,indep_var,indep_table,n,cycle
RSV_genus1002_relative,DADA2RSV_GENUS_RELATIVE_F_clr,RIDAGEYR,DEMO_F,4679,F
RSV_genus1002_relative,DADA2RSV_GENUS_RELATIVE_F_clr,RIAGENDR,DEMO_F,4679,F
```

---

## 🔧 **Dependencies**

### **R Packages Required:**

```r
# For transformation scripts (nhanes_omp_transformation.R)
optparse, DBI, RSQLite, dplyr, fs, compositions

# For schema structure creation (ss_file_create.R)  
optparse, DBI, RSQLite, dplyr, purrr, readr, fs, glue, stringr
```

### **System Requirements:**

- **SLURM** cluster access
- **Environment**: Uses either R/4.2.1 module OR conda environment
- **Memory**: 64G for transformation, 16G for schema creation
- **Time**: ~2 hours for transformation, ~1 hour per schema job (18 parallel jobs)

---

## **Quick Start Commands**

### **Complete Workflow (Sequential Steps):**

```bash
# Step 1: Transform microbiome data (SLURM job with R module)
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh

# Step 2: Fill missing data (run after Step 1 completes - requires conda)
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
  --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite

# Step 3: Create 24 schema structure files (parallel SLURM jobs with conda)
bash scripts/0_transform_n_preprocess_ssfiles/run_ss_file_create.sh
```

### **Environment Differences:**

- **Step 1**: Uses `module load gcc/14.2.0 R/4.2.1`
- **Steps 2 & 3**: Use conda environment `/n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis`

### **Monitor Progress:**

```bash
# Check SLURM jobs
squeue -u $USER

# Monitor transformation logs
tail -f logs/trans_*.out

# Monitor schema structure jobs (24 jobs named ss_*)
squeue -u $USER | grep ss_

# Check transformation output (should see both databases)
ls -la data/00_nhanes_omp_transformed_db/

# Check schema structure files (should be exactly 24 files)
ls -la results/0_ss_files/*.csv | wc -l
```

---

## **Success Indicators**

### **After Step 1 (Transformation):**

- Database created: `nhanes_oral_transformed.sqlite`
- **40 new transformed tables**: 10 taxonomic levels × 4 transformations each
  - From RELATIVE tables: `DADA2RSV_{GENUS,FAMILY,ORDER,CLASS,PHYLUM}_RELATIVE_{F,G}` → `_none`, `_hellinger`
  - From COUNT tables: `DADA2RSV_{GENUS,FAMILY,ORDER,CLASS,PHYLUM}_COUNT_{F,G}` → `_clr`, `_lognorm`
- **Metadata automatically cloned** for each transformed table
- **ALL taxa retained** (no prevalence filtering)
- **Universal pseudocount** ε = 0.001 applied to COUNT matrices only
- Log shows: "Transformation completed with four scenarios"

### **After Step 2 (Missing Data Filling):**

- Complete database created: `nhanes_oral_transformed_complete.sqlite`
- RIAGENDR_01** columns added to both DEMO_F and DEMO_G tables
- 40+ variable annotations** added to `variable_names_epcf` table
- Message: "Completed: [database_path]"

### **After Step 3 (Schema Structure Creation):**

- **Exactly 24 CSV files** created in `results/0_ss_files/`:
  ```
  1_demoWAS_clr_schema_structure.csv     # Demographics → Microbiome
  1_demoWAS_hellinger_schema_structure.csv
  1_demoWAS_lognorm_schema_structure.csv
  1_demoWAS_none_schema_structure.csv
  2_oradWAS_clr_schema_structure.csv     # Microbiome → Oral Health  
  ... (continuing for all 6 analysis types × 4 transformations = 24 files)
  ```
- Each file contains variable pairs with sample counts and cycle information
- **OTU role correctly assigned**: `dep` for 1_demoWAS & 3_exWAS, `indep` for others

---

## **Troubleshooting**

### **Debug Tool Available:**

```bash
# Use debug script to test schema structure creation manually
bash scripts/0_transform_n_preprocess_ssfiles/debug_ss.sh
```

**What debug_ss.sh actually does:**

1. ✅ **Path verification**: Checks BASE, DB_PATH, OUT_DIR
2. ✅ **Database validation**: Confirms existence and reports file size
3. ✅ **Directory creation**: Creates output directory
4. ✅ **Config file check**: Validates `3_exWAS_vars.txt` exists and shows first 5 variables
5. ✅ **Environment test**: Loads conda environment and checks R version
6. ✅ **Single test run**: Executes schema creation for `3_exWAS clr` specifically
7. ✅ **Output validation**: Checks if expected file created and reports size/row count

### **Common Issues:**

#### **1. Environment Mismatch**

```bash
# Problem: Step 1 uses R module, Steps 2-3 use conda
# Solution: Ensure correct environment for each step
# Step 1: module load gcc/14.2.0 R/4.2.1
# Steps 2-3: eval "$(conda shell.bash hook)" && conda activate /n/groups/.../nhanes-analysis
```

#### **2. Missing Configuration Files**

```bash
# Problem: *_vars.txt files not found in configs/
# Solution: Verify all 6 config files exist
ls -la configs/*_vars.txt
# Should show: 1_demoWAS_vars.txt through 6_zimWAS_vars.txt
```

#### **3. Schema Structure Jobs Failing**

```bash
# Problem: Some of 24 jobs fail
# Solution: Check logs and rerun specific failed jobs
find logs -name "ss_*_*.err" -size +0  # Check for errors
# Rerun individual job if needed:
# sbatch [recreate job script for specific analysis_transform combination]
```

#### **4. Output File Location**

```bash
# Problem: Files not appearing in expected location
# Actual output: results/0_ss_files/
```

---

## **File Structure After Completion**

```
data/00_nhanes_omp_transformed_db/
├── nhanes_oral_transformed.sqlite           # Step 1 output
└── nhanes_oral_transformed_complete.sqlite  # Step 2 output (final)

results/0_ss_files/                          # Step 3 output (NOTE: results/ )
├── 1_demoWAS_clr_schema_structure.csv
├── 1_demoWAS_hellinger_schema_structure.csv
├── 1_demoWAS_lognorm_schema_structure.csv
├── 1_demoWAS_none_schema_structure.csv
├── 2_oradWAS_clr_schema_structure.csv
├── ... (24 total files)
└── 6_zimWAS_none_schema_structure.csv

logs/                                        # All execution logs
├── trans_*.out/.err                         # Transformation logs
├── ss_*_*.out/.err                         # Schema structure logs (24 files)
└── temp_ss_*.sh                            # Temporary job scripts (auto-deleted)
```

---

## 📝 **Next Steps**

After completing this preprocessing workflow:

1. **WAS Analyses**: Use the 24 schema structure files with the association pipeline:

   ```bash
   ./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh
   ```
2. **Quality Control**: Examine generated databases and validate schema files:

   ```bash
   # Check database tables (should see 40 transformed tables with RELATIVE naming)
   sqlite3 data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite ".tables" | grep "_clr\|_lognorm\|_none\|_hellinger"

   # Validate schema structure files (should be exactly 24)
   ls results/0_ss_files/*.csv | wc -l
   
   # Verify naming convention: all output tables use RELATIVE format
   sqlite3 data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite ".tables" | grep "RELATIVE.*_clr\|RELATIVE.*_lognorm\|RELATIVE.*_none\|RELATIVE.*_hellinger" | head -5
   ```
