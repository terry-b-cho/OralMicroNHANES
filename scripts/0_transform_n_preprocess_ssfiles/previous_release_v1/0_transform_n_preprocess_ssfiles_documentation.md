# NHANES Oral Microbiome Data Transformation and Schema Structure Preprocessing Guide

## **Overview**

This directory contains the complete preprocessing pipeline for NHANES oral microbiome data analysis. The workflow consists of two main stages:

1. **Data Transformation**: Apply normalization methods to microbiome OTU abundance data
2. **Schema Structure Creation**: Generate mapping files for downstream WAS analyses

## 📁 **Files in This Directory**

### **Core Scripts:**
1. **`nhanes_omp_transformation.R`** - Transforms microbiome abundance data (CLR, log-normal, none)
2. **`nhanes_db_filling_missing_data.R`** - Fills missing metadata and derives additional variables
3. **`ss_file_create.R`** - Creates schema structure files for WAS analyses

### **SLURM Submission Scripts:**
4. **`run_transformation.sh`** - Submits transformation jobs to SLURM
5. **`run_ss_file_create.sh`** - Submits schema structure creation jobs to SLURM

---

## **Complete Workflow - Step by Step**

### **Prerequisites**
- Access to O2 HPC cluster
- R ≥ 4.2.1 with required packages
- Input database: `data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite`
- Configuration files in `configs/` directory

### **Step 1: Data Transformation**

#### **What This Does:**
- Processes NHANES oral microbiome relative abundance data
- Applies prevalence filtering (default: 0.1% of samples)
- Creates three versions of each microbiome table:
  - **"none"**: Raw abundances
  - **"clr"**: Centered Log-Ratio transformation
  - **"lognorm"**: Log normalization with pseudo-count

#### **Input Tables Processed:**
- `DADA2RSV_GENUS_RELATIVE_F/G`
- `DADA2RSV_FAMILY_RELATIVE_F/G`
- `DADA2RSV_ORDER_RELATIVE_F/G`
- `DADA2RSV_CLASS_RELATIVE_F/G`
- `DADA2RSV_PHYLUM_RELATIVE_F/G`

#### **Command:**
```bash
# Submit transformation job
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh
```

#### **Manual Execution (if needed):**
```bash
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R \
  --in_db  data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --prev   0.001
```

#### **Output:**
- **Database**: `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite`
- **New Tables**: Each original table gets 3 versions (e.g., `DADA2RSV_GENUS_RELATIVE_F_clr`, `DADA2RSV_GENUS_RELATIVE_F_lognorm`, `DADA2RSV_GENUS_RELATIVE_F_none`)

### **Step 2: Missing Data Filling and Variable Derivation**

#### **What This Does:**
- Derives `RIAGENDR_01` (0=Male, 1=Female) from `RIAGENDR`
- Fills missing metadata entries in `variable_names_epcf` table
- Maps new variables to their respective tables

#### **Command:**
```bash
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
  --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
```

#### **Output:**
- **Complete Database**: `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`

### **Step 3: Schema Structure File Creation**

#### **What This Does:**
- Creates mapping files for each WAS analysis type and normalization method
- Maps microbiome OTU variables to NHANES variables
- Calculates sample sizes for each variable pair
- Generates 18 schema structure files (6 analysis types × 3 normalizations)

#### **Analysis Types and Variable Roles:**

| Pipeline | Description | Microbiome Role | Config File |
|----------|-------------|-----------------|-------------|
| **1_demoWAS** | Demographics → Microbiome | Dependent | `configs/1_demoWAS_vars.txt` |
| **2_oradWAS** | Microbiome → Oral Health | Independent | `configs/2_oradWAS_vars.txt` |
| **3_exWAS** | Exposures → Microbiome | Dependent | `configs/3_exWAS_vars.txt` |
| **4_pheWAS** | Microbiome → Phenotypes | Independent | `configs/4_pheWAS_vars.txt` |
| **5_outWAS** | Microbiome → Disease Outcomes | Independent | `configs/5_outWAS_vars.txt` |
| **6_zimWAS** | Microbiome → Lab Measurements | Independent | `configs/6_zimWAS_vars.txt` |

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
- **18 Total Files**: e.g., `1_demoWAS_clr_schema_structure.csv`, `2_oradWAS_lognorm_schema_structure.csv`

---

## **Output File Structure**

### **Schema Structure File Format:**
Each CSV contains the following columns:

| Column | Description |
|--------|-------------|
| `dep_var` | Dependent variable name |
| `dep_table` | NHANES table where dependent variable came from |
| `indep_var` | Independent variable name |
| `indep_table` | NHANES/OTU table where independent variable came from |
| `n` | Number of non-missing pairs |
| `cycle` | NHANES cycle (F or G) |

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
# For transformation scripts
optparse, DBI, RSQLite, dplyr, fs, compositions

# For schema structure creation
optparse, DBI, RSQLite, dplyr, purrr, readr, fs, glue, stringr
```

### **System Requirements:**
- **SLURM** cluster access
- **R** ≥ 4.2.1
- **Memory**: 32G for transformation, 16G for schema creation
- **Time**: ~1 hour for transformation, ~30 minutes for schema creation

---

## **Quick Start Commands**

### **Complete Workflow (Recommended):**
```bash
# Step 1: Transform microbiome data
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh

# Step 2: Fill missing data (run after Step 1 completes)
Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
  --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite

# Step 3: Create schema structure files (run after Step 2 completes)
bash scripts/0_transform_n_preprocess_ssfiles/run_ss_file_create.sh
```

### **Monitor Progress:**
```bash
# Check SLURM jobs
squeue -u $USER

# Check transformation output
ls -la data/00_nhanes_omp_transformed_db/

# Check schema structure files
ls -la results/0_ss_files/
```

---

## ✅ **Success Indicators**

### **After Step 1 (Transformation):**
- ✅ Database created: `nhanes_oral_transformed.sqlite`
- ✅ New tables with suffixes: `_clr`, `_lognorm`, `_none`
- ✅ Log shows successful transformation of all taxonomic levels

### **After Step 2 (Missing Data Filling):**
- ✅ Complete database created: `nhanes_oral_transformed_complete.sqlite`
- ✅ `RIAGENDR_01` column added to DEMO tables
- ✅ Metadata rows appended to `variable_names_epcf`

### **After Step 3 (Schema Structure Creation):**
- ✅ 18 CSV files created in `results/0_ss_files/`
- ✅ Each file contains variable pairs with sample sizes
- ✅ Files ready for downstream WAS analyses

---

## **Troubleshooting**

### **Common Issues:**

#### **1. Missing Input Database**
```bash
# Problem: nhanes_031725.sqlite not found
# Solution: Check path and ensure database exists
ls -la data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite
```

#### **2. R Package Installation**
```bash
# Problem: Missing R packages
# Solution: Install required packages
R -e "install.packages(c('optparse', 'DBI', 'RSQLite', 'dplyr', 'compositions'))"
```

#### **3. Configuration Files Missing**
```bash
# Problem: Config files not found
# Solution: Check configs directory
ls -la configs/*_vars.txt
```

#### **4. Memory Issues**
```bash
# Problem: Jobs fail due to memory
# Solution: Increase memory in SLURM scripts
# Edit run_transformation.sh: #SBATCH --mem=64G
```

---

## 📝 **Next Steps**

After completing this preprocessing workflow, you can proceed to:

1. **WAS Analyses**: Use the schema structure files with the association pipeline
2. **Quality Control**: Examine the generated databases and schema files
3. **Custom Analyses**: Modify configuration files for specific variable sets

**Your NHANES oral microbiome data is now fully preprocessed and ready for Wide Association Studies!** 