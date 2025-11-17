# Chapter 2: Database Preprocessing & phyloseq Object Creation

> **📖 Full Documentation for Stage 2 of the NHANES Oral Microbiome Analysis Pipeline**

## **RECOMMENDED WORKFLOW**

**Follow this sequence for optimal results:**

1. **FIRST:** Run demographic processing to create comprehensive processed database
2. **SECOND:** Run diagnostic tool to verify processing was successful
3. **ANYTIME:** Re-run diagnostic tool to check database status

```bash
# Step 1: Create processed database with all demographic variables
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R

# Step 2: Verify processing was successful
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R --processed

# Step 3: Check original database status anytime
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R
```

---

## Overview

This chapter focuses on database preprocessing, diagnostics, and phyloseq object creation. It provides tools to:

1. **Create comprehensive demographic variables** - Process DEMO_F/DEMO_G tables with 32 derived variables
2. **Diagnose database completeness** - Check what tables/variables exist in original or processed databases
3. **Prevent redundant data creation** - Avoid recreating existing derived tables
4. **Create phyloseq objects** - Convert data for microbiome analysis (future)
5. **Quality control** - Validate data integrity and completeness

---

## **STEP 1: Comprehensive Demographic Processing**

### **Script:** `process_demographics_complete.R`

**Purpose:** One-step comprehensive demographic processing that creates derived variables, calculates quartiles, applies factor levels, and saves to a new processed database with complete metadata registration.

**⚠️ RUN THIS FIRST** - This creates the processed database that contains all demographic variables ready for analysis.

### **Usage:**

```bash
# Test mode first (recommended to verify what will be created)
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R --test

# Create processed database (creates 32 new variables in DEMO_F and DEMO_G)
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R

# Custom database paths (if needed)
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R \
  --input_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
  --output_db data/custom_processed.sqlite
```

### **Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input_db` | Path to input SQLite database | `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite` |
| `--output_db` | Path to output processed SQLite database | `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` |
| `--test` | Test mode - don't create output database | `FALSE` |

### **What It Creates:**

#### **32 Derived Variables Added to DEMO_F and DEMO_G:**

**Core Variables:**
- `cycle` - NHANES Cycle (F=2009-2010, G=2011-2012)
- `RIAGENDR_01` - Gender (0=Male, 1=Female)
- `Gender` - Gender (categorical)
- `AgeGroup` - Age groups (14-19, 20-29, 30-39, 40-49, 50-59, 60-69)
- `Age` - Age in years (capped at 85)

**Education Variables:**
- `EducationLevel` - Education level (5 categories)
- `Education_Level` - Education level (alternative naming)

**Ethnicity & Citizenship:**
- `Ethnicity` - Race/ethnicity (5 categories)
- `BornInUSA` - Born in USA status
- `US_Born` - US born status (alternative naming)
- `US_Citizen` - US citizenship status (alternative naming)

**Income & Poverty:**
- `PIRcat` - Poverty income ratio categories (<1, 1-2, 2-4, >4)
- `Ratio_Family_Income_Poverty` - Family income to poverty ratio (continuous)

**Household Variables:**
- `Household_Size` - Household size (capped at 7)
- `Family_Size` - Family size (capped at 7)
- `Household_Size_Factor` - Household size as factor (1,2,3,4,5,6,7+)

**Marital Status:**
- `Marital_Status` - Marital status (6 categories)

**Interview Variables:**
- `Interview_Language` - Language of interview (English/Spanish)
- `Interview_Proxy` - Proxy used in interview (Yes/No)
- `Interview_Interpreter` - Interpreter used in interview (Yes/No)

**Household Reference Person:**
- `HH_Reference_Gender` - Household reference person gender
- `HH_Reference_Age` - Household reference person age
- `HH_Reference_Education_Level` - Household reference person education
- `HH_Reference_Marital_Status` - Household reference person marital status

**Detailed Income Categories:**
- `Annual_Household_Income` - Annual household income (detailed categories)
- `Annual_Family_Income` - Annual family income (detailed categories)

### **🔄 Alternative Variable Names**

The demographic processing script includes **multiple alternative names** for compatibility with different parts of your pipeline:

#### **Age Variables:**
- `AgeGroup` → "14-19", "20-29", "30-39", "40-49", "50-59", "60-69"
- `age_group` → "14–19", "20–29", "30–39", "40–49", "50–59", "60–69" (en-dash format)
- `Age` → Continuous age (capped at 85)
- `AgeQuartile` → Dynamic quartiles calculated from data

#### **Education Variables:**
- `EducationLevel` → "< 9th Grade", "9-11th Grade", "High School", "College/AA", "College Graduate"
- `Education_Level` → Same as above (alternative naming)
- `EducationLevel_Simple` → Simplified for phyloseq: "High School", "College/AA", "College Graduate"

#### **Gender Variables:**
- `RIAGENDR_01` → Binary (0=Male, 1=Female)
- `Gender` → Categorical ("Male", "Female")

#### **Income Variables:**
- `PIRcat` → Categories ("<1", "1-2", "2-4", ">4")
- `PIRQuartile` → Dynamic quartiles ("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%")
- `Ratio_Family_Income_Poverty` → Continuous PIR value

#### **Sample ID Variables:**
- `SEQN` → Original sample ID
- `sample` → Alternative name for phyloseq compatibility

#### **Household Variables:**
- `Household_Size` → Numeric (1-7)
- `Household_Size_Factor` → Factor ("1", "2", "3", "4", "5", "6", "7+")

#### **Citizenship Variables:**
- `BornInUSA` → "US Born", "non-US Born"
- `US_Born` → Same as above (alternative)
- `US_Citizen` → Same as above (alternative)

#### **Cycle Variables:**
- `cycle` → "F", "G"
- `Release_Cycle` → "2009-2010", "2011-2012"

### **What the Processing Pipeline Does**

The `process_demographics_complete.R` script combines all steps:

1. **Demographic Variable Derivation** - Creates 31 new variables
2. **Age Quartile Calculation** - Data-driven age quartiles
3. **PIR Quartile Calculation** - Data-driven income quartiles  
4. **Factor Level Application** - Proper R factor levels with ordering
5. **Database Creation** - Saves to new processed SQLite file
6. **Metadata Registration** - Complete variable metadata

### **✅ Expected Results**

After running the processing script successfully:

```
COMPREHENSIVE DEMOGRAPHIC PROCESSING SUMMARY
✅ DEMO_F: 10,537 rows, 90 columns (added 31 variables)
✅ DEMO_G: 9,756 rows, 95 columns (added 31 variables)
📝 Total variables in metadata: [updated count]
📂 Processed database saved: data/.../nhanes_oral_transformed_complete_processed.sqlite
```

---

## **STEP 2: Database Diagnostic Verification**

### **Script:** `missing_tables_n_variables_diagnostic.R`

**Purpose:** Comprehensive diagnostic script that checks what tables, variables, and derived data exist in either the original or processed NHANES SQLite database. Use this to verify that Step 1 worked correctly and to check database status anytime.

**⚠️ RUN THIS AFTER STEP 1** - This verifies that the demographic processing was successful.

### **Usage:**

```bash
# Check processed database (after running Step 1)
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R --processed

# Check original database status
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R

# Custom database path with detailed output
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R \
  --db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
  --output_dir results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics/ \
  --format both

# Quick check (CSV output only)
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R \
  --format csv
```

### **Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--db` | Path to SQLite database | `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite` |
| `--processed` | Check processed database instead of original | `FALSE` |
| `--output_dir` | Output directory for reports | `results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics/` |
| `--format` | Output format: `html`, `csv`, or `both` | `html` |

### **What It Checks:**

#### **1. Database Structure**
- Total number of tables
- Metadata table presence (`variable_names_epcf`, `table_names_epcf`)
- Table sizes and column counts
- Table categorization (original, derived, transformed, metadata)

#### **2. Expected Derived Tables**
- **Oral health tables:** `OralDisease_F`, `OralDisease_G`, `d_outcome_oh_f`, `d_outcome_oh_g`
- **Disease outcome tables:** `d_outcome_mcq_f`, `d_outcome_mcq_g`
- **Microbiome transformed tables:** All 30 transformation combinations (5 taxa × 2 cycles × 3 transforms)

#### **3. Expected Derived Variables**
- **Demographics:** `AGE_SQUARED`, `BORN_INUSA`, education categories, ethnicity categories, `RIAGENDR_01`
- **Comprehensive demographics:** All 32 variables created by Step 1
- **Oral health:** `TOOTH_DECAY_OHAROCDT`, `GUM_DISEASE_OHAROCGP`, `ORAL_HYGIENE_OHAROCOH`, `DENTURE_OHAROCDE`
- **Cardiovascular:** `CHD`, `STROKE`, `HEART_ATTACK`, `HEART_FAILURE`, `ANGINA`, `CVD`
- **Respiratory:** `ASTHMA`, `BRONCHITIS`, `EMPHYSEMA`
- **Cancer:** `CANCER_BREAST`, `CANCER_COLON`, `CANCER_LUNG`, `CANCER_ESOPHAGEAL`, `CANCER_PROSTATE`, `CANCER_MOUTH`
- **Metabolic:** `DIABETES`
- **Lab measures:** Various lab measurement variables

#### **4. Microbiome Transformations**
- Checks all 30 expected transformation tables
- Validates row counts and column counts
- Reports missing transformations

#### **5. Demographic Processing Status (when --processed flag used)**
- Verifies all 32 comprehensive demographic variables exist
- Checks proper metadata registration
- Confirms factor levels are applied correctly

### **Expected Output After Step 1:**

#### **Original Database Status:**
```
============================================================
DATABASE DIAGNOSTIC SUMMARY
============================================================
Total tables in database: 1460
Metadata tables present: ✅ YES

EXPECTED TABLES:
  Total expected: 36
  ✅ Exist: 36
  ❌ Missing: 0

EXPECTED VARIABLES:
  Total expected: 45
  ✅ Complete (exist + metadata): 45
  ⚠️ Exist but no metadata: 0
  ❌ Missing: 0

👥 DEMOGRAPHIC PROCESSING (DEMO_F & DEMO_G):
  Total expected variables: 64
  ✅ Complete (exist + metadata): 2
  ⚠️ Exist but no metadata: 0
  ❌ Variables missing: 62
============================================================
```

#### **Processed Database Status (after Step 1):**
```
============================================================
DATABASE DIAGNOSTIC SUMMARY
============================================================
Total tables in database: 1460
Metadata tables present: ✅ YES

EXPECTED TABLES:
  Total expected: 36
  ✅ Exist: 36
  ❌ Missing: 0

EXPECTED VARIABLES:
  Total expected: 77
  ✅ Complete (exist + metadata): 77
  ⚠️ Exist but no metadata: 0
  ❌ Missing: 0

👥 DEMOGRAPHIC PROCESSING (DEMO_F & DEMO_G):
  Total expected variables: 64
  ✅ Complete (exist + metadata): 64
  ⚠️ Exist but no metadata: 0
  ❌ Variables missing: 0
============================================================
```

### **Output Reports:**

#### **Console Summary:**
Real-time summary of database status with clear indicators of what exists and what's missing.

#### **Detailed CSV Reports (when `--format csv` or `both`):**
- `database_table_info.csv` - Complete table inventory
- `expected_tables_status.csv` - Status of all expected derived tables
- `expected_variables_status.csv` - Status of all expected derived variables
- `demographic_processing_status.csv` - Detailed demographic variable status
- `transformation_status.csv` - Status of microbiome transformations
- `action_*.csv` - Action plans for missing items

#### **HTML Report (when `--format html` or `both`):**
- `diagnostic_report.html` - Summary HTML report

### **Status Codes:**

| Status | Meaning |
|--------|---------|
| ✅ COMPLETE | Variable exists in tables AND has metadata |
| ⚠️ EXISTS BUT NO METADATA | Variable exists in tables but missing from `variable_names_epcf` |
| ❌ MISSING | Variable does not exist in any table |
| ✅ EXISTS | Table exists in database |
| ❌ MISSING | Table does not exist in database |

---

## ✅ **CURRENT DATABASE STATUS INFORMATION**

### **ALL EXPECTED TABLES EXIST IN ORIGINAL DATABASE:**

✅ **Oral Health Tables:**
- `OralDisease_F` (8,189 rows) ✅ EXISTS
- `OralDisease_G` (8,956 rows) ✅ EXISTS  
- `d_outcome_oh_f` (8,189 rows) ✅ EXISTS
- `d_outcome_oh_g` (8,956 rows) ✅ EXISTS

✅ **Disease Outcome Tables:**
- `d_outcome_mcq_f` (10,109 rows) ✅ EXISTS
- `d_outcome_mcq_g` (9,364 rows) ✅ EXISTS

✅ **Microbiome Transformation Tables:**
- All 30 transformation combinations ✅ COMPLETE
- 5 taxa levels × 2 cycles × 3 transforms (clr, lognorm, none)

### **ALL EXPECTED VARIABLES EXIST WITH METADATA IN ORIGINAL DATABASE:**

✅ **Oral Health Variables:**
- `TOOTH_DECAY_OHAROCDT` → "Tooth Decay" ✅ COMPLETE
- `GUM_DISEASE_OHAROCGP` → "Gum Disease" ✅ COMPLETE
- `ORAL_HYGIENE_OHAROCOH` → "Oral Hygiene" ✅ COMPLETE
- `DENTURE_OHAROCDE` → "Denture" ✅ COMPLETE

✅ **MCQ/DIQ Variables (Disease Outcomes):**
- `CHD`, `STROKE`, `HEART_ATTACK`, `HEART_FAILURE`, `ANGINA`, `CVD` ✅ COMPLETE
- `ASTHMA`, `BRONCHITIS`, `EMPHYSEMA` ✅ COMPLETE
- `CANCER_BREAST`, `CANCER_COLON`, `CANCER_LUNG`, `CANCER_ESOPHAGEAL`, `CANCER_PROSTATE`, `CANCER_MOUTH` ✅ COMPLETE
- `DIABETES` ✅ COMPLETE

✅ **Basic Demographics Variables:**
- `AGE_SQUARED`, `BORN_INUSA`, `RIAGENDR_01` ✅ COMPLETE
- Education categories: `EDUCATION_LESS9`, `EDUCATION_9_11`, `EDUCATION_AA`, `EDUCATION_COLLEGEGRAD` ✅ COMPLETE
- Ethnicity categories: `ETHNICITY_MEXICAN`, `ETHNICITY_OTHERHISPANIC`, `ETHNICITY_OTHER`, `ETHNICITY_NONHISPANICBLACK` ✅ COMPLETE

### **⚠️ COMPREHENSIVE DEMOGRAPHIC VARIABLES NEED PROCESSING:**

❌ **Missing from Original Database (Created by Step 1):**
- All 32 comprehensive demographic variables (cycle, Gender, AgeGroup, EducationLevel, etc.)
- Data-driven quartiles (AgeQuartile, PIRQuartile)
- Proper factor levels for categorical variables
- Alternative naming conventions for compatibility

**→ This is why Step 1 (demographic processing) is essential!**

---

## 🛡️ **When to Use Each Script**

### **Use Process Demographics Complete (`process_demographics_complete.R`):**
- **FIRST TIME:** To create comprehensive processed database
- **When you need:** All demographic variables with proper factors and quartiles
- **Before:** Any analysis requiring demographic variables
- **If:** You want alternative variable names for compatibility

### **Use Diagnostic Tool (`missing_tables_n_variables_diagnostic.R`):**
- **AFTER Step 1:** To verify demographic processing was successful
- **BEFORE creating new tables/variables:** To check what already exists
- **ANYTIME:** To check database status
- **REPEATEDLY:** As many times as needed for verification
- **WHEN DEBUGGING:** To troubleshoot missing data issues

### **Example Integration in Your Workflow:**
```bash
#!/bin/bash
# Complete Chapter 2 Workflow

echo "Step 1: Creating comprehensive processed database..."
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R

echo "Step 2: Verifying processing was successful..."
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R --processed

echo "Step 3: Checking original database status for comparison..."
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R

echo "✅ Chapter 2 processing complete!"
```

---

## **Using the Processed Database**

After completing Steps 1 and 2, you'll have a processed database with properly factored variables ready for analysis:

```r
# Connect to processed database
con <- dbConnect(SQLite(), "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")

# Load processed demographic data (factors already applied!)
demo_f <- tbl(con, "DEMO_F") %>% collect()
demo_g <- tbl(con, "DEMO_G") %>% collect()

# Combine cycles if needed
demographic_data <- bind_rows(demo_f, demo_g)

# Variables are already properly factored with quartiles calculated!
str(demographic_data$AgeQuartile)    # Ordered factor with quartile labels
str(demographic_data$EducationLevel) # Ordered factor with education levels
str(demographic_data$PIRQuartile)   # Ordered factor with income quartiles

# Ready for analysis - no additional processing needed!
head(demographic_data)
```

---

## **Quick Start Example**

```bash
# Navigate to project root
cd /path/to/nhanes_oral_mirco_cho

# Step 1: Create processed database (REQUIRED FIRST)
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R

# Step 2: Verify processing was successful
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R --processed

# Optional: Check original database for comparison
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R

# Check results
ls data/00_nhanes_omp_transformed_db/
# Should see: nhanes_oral_transformed_complete_processed.sqlite

# View diagnostic results
ls results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics*/
cat results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics_processed/diagnostic_report.html  # or open in browser
```

---

## 🔗 **Integration with Main Pipeline**

This chapter integrates with the overall NHANES analysis pipeline:

```
Chapter 0: Data Transformation
     ↓
Chapter 1: Association Analysis  
     ↓
Chapter 2: Database Processing ← 📍 YOU ARE HERE
     │
     ├── Step 1: Create Processed Database (process_demographics_complete.R)
     └── Step 2: Verify Processing (missing_tables_n_variables_diagnostic.R)
     ↓
Chapter 3: GOLD Database Integration
     ↓  
Chapter 4: Downstream Analysis
```

### **Typical Workflow:**
1. **Complete Chapters 0-1** (transformation and analysis)
2. **Run Step 1** (demographic processing) to create processed database
3. **Run Step 2** (diagnostic verification) to confirm processing success
4. **Use diagnostic tool anytime** to check database status
5. **Proceed to Chapter 3-4** with complete processed database

---

## **Future Extensions**

### **Planned Features:**
- [ ] **phyloseq object creation** from processed database
- [ ] **Data quality checks** (missing values, outliers)
- [ ] **Schema validation** against expected structure
- [ ] **Interactive HTML dashboard** for diagnostics
- [ ] **Integration with SLURM** for HPC environments

### **Contributing:**
To add new expected tables/variables to the diagnostic:

1. Edit the `expected_derived_vars` or `expected_derived_tables` lists in the diagnostic script
2. Add appropriate categorization
3. Update documentation

---

## 🛠️ **Troubleshooting**

### **Common Issues:**

#### **Database Connection Error:**
```
Error: database is locked
```
**Solution:** Ensure no other processes are using the SQLite file.

#### **Missing Packages:**
```
Error: there is no package called 'fs'
```
**Solution:** Install missing packages:
```r
install.packages(c("optparse", "DBI", "RSQLite", "dplyr", "tidyr", "purrr", "knitr", "tibble", "fs", "forcats"))
```

#### **Permission Denied:**
```
Error: cannot create directory 'results/diagnostics'
```
**Solution:** Check write permissions or change output directory:
```bash
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R \
  --output_dir ~/diagnostics/
```

#### **Step 1 Failed:**
```
Error: could not create processed database
```
**Solution:** Check input database exists and you have write permissions:
```bash
# Verify input database exists
ls -la data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite

# Run in test mode first
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R --test
```

#### **Step 2 Shows Missing Variables:**
```
❌ Variables missing: 62
```
**Solution:** This means Step 1 hasn't been run yet or failed. Run Step 1 first:
```bash
# Run Step 1 first
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R

# Then verify with Step 2
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R --processed
```

---

## 📄 **File Structure**

```
scripts/2_preprocess_db_n_phyloseq/
├── process_demographics_complete.R            # STEP 1: Complete demographic processing pipeline
├── missing_tables_n_variables_diagnostic.R    # STEP 2: Database diagnostic verification
├── 2_preprocess_db_n_phyloseq_documentation.md # 📖 This documentation
├── logs/                                       # 📁 Log files directory
│   ├── process_demographics_complete_YYYYMMDD_HHMMSS.log
│   ├── process_demographics_complete_summary_YYYYMMDD_HHMMSS.txt
│   ├── diagnostic_YYYYMMDD_HHMMSS.log
│   └── diagnostic_summary_YYYYMMDD_HHMMSS.txt
└── [future scripts for phyloseq creation]

data/00_nhanes_omp_transformed_db/              # Database files
├── nhanes_oral_transformed_complete.sqlite          # Original database
└── nhanes_oral_transformed_complete_processed.sqlite # Processed database (created by Step 1)

results/analyses_results/02_preprocess_db_n_phyloseq_out/
├── diagnostics/                                # Original database diagnostic outputs
│   ├── database_table_info.csv
│   ├── expected_tables_status.csv
│   ├── expected_variables_status.csv
│   ├── demographic_processing_status.csv
│   ├── transformation_status.csv
│   ├── action_*.csv
│   └── diagnostic_report.html
└── diagnostics_processed/                      # Processed database diagnostic outputs
    ├── database_table_info.csv
    ├── expected_tables_status.csv
    ├── expected_variables_status.csv
    ├── demographic_processing_status.csv
    ├── transformation_status.csv
    ├── action_*.csv
    └── diagnostic_report.html
```

### **📁 Logging System**

All scripts in this chapter automatically save detailed logs to `scripts/2_preprocess_db_n_phyloseq/logs/`:

#### **Log File Types:**

**Step 1 - Demographic Processing Logs:**
- `process_demographics_complete_YYYYMMDD_HHMMSS.log` - Complete processing pipeline log
- `process_demographics_complete_summary_YYYYMMDD_HHMMSS.txt` - Summary of complete processing

**Step 2 - Diagnostic Script Logs:**
- `diagnostic_YYYYMMDD_HHMMSS.log` - Complete timestamped execution log
- `diagnostic_summary_YYYYMMDD_HHMMSS.txt` - Summary report with key statistics

#### **Log Content Examples:**

**Step 1 Processing Summary:**
```
COMPREHENSIVE DEMOGRAPHIC PROCESSING SUMMARY
Timestamp: 2025-06-18 14:26:58
Input database: data/.../nhanes_oral_transformed_complete.sqlite
Output database: data/.../nhanes_oral_transformed_complete_processed.sqlite
Test mode: NO

VARIABLES CREATED:
- cycle (NHANES Cycle)
- RIAGENDR_01 (Binary gender)
- Gender (Categorical gender)
- AgeGroup (Age categories)
[... 28+ more variables ...]

RESULTS:
✅ DEMO_F: 10,537 rows, 90 columns (added 31 variables)
✅ DEMO_G: 9,756 rows, 95 columns (added 31 variables)
📝 Total variables in metadata: [updated count]

STATUS: APPLIED TO PROCESSED DATABASE
```

**Step 2 Diagnostic Summary:**
```
DATABASE DIAGNOSTIC SUMMARY
Timestamp: 2025-06-18 14:27:12
Database: data/.../nhanes_oral_transformed_complete_processed.sqlite
Output directory: results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics_processed/
Format: html

RESULTS:
Total tables: 1460
Expected tables exist: 36/36
Expected variables complete: 77/77
Demographic processing: 64/64 variables complete
Microbiome transformations exist: 30/30

STATUS: ALL DEMOGRAPHIC PROCESSING COMPLETE ✅
```

#### **Benefits of Logging:**

✅ **Complete audit trail** of all operations
✅ **Timestamps** for all steps
✅ **Error tracking** and debugging information  
✅ **Reproducibility** - exact record of what was done
✅ **Summary reports** for quick status checks
✅ **Automatic organization** by timestamp
✅ **Step-by-step verification** of workflow progress

---

## **STEP 3: phyloseq Object Creation**

### **Script:** `ubiome_phyloseq_object_creation_cho.Rmd`

**Purpose:** Creates comprehensive phyloseq objects from the processed NHANES oral microbiome database for downstream microbiome analysis. This step converts the SQLite database format into phyloseq objects optimized for different analysis types.

**⚠️ RUN THIS AFTER STEPS 1 & 2** - This requires the processed database with comprehensive demographics.

### **Usage:**

```bash
# Navigate to the script directory
cd scripts/2_preprocess_db_n_phyloseq/

# Run the RMarkdown file (in R/RStudio)
# Or knit from command line:
Rscript -e "rmarkdown::render('ubiome_phyloseq_object_creation_cho.Rmd')"
```

### **Prerequisites:**

1. ✅ **Processed database exists**: `nhanes_oral_transformed_complete_processed.sqlite` (from Step 1)
2. ✅ **All transformations complete**: 30 microbiome transformation tables (from Chapter 0)
3. ✅ **Comprehensive demographics**: 84 demographic variables (from Step 1)
4. ✅ **Required R packages**: phyloseq, ape, dplyr, tidyr, DBI, RSQLite, etc.

---

## **What phyloseq Objects Are Created**

### **✅ Successfully Created phyloseq Objects**

| Object | Taxa Count | Samples | Description | Use Cases |
|--------|------------|---------|-------------|-----------|
| **`ubiome_counts`** | 1,349 | 9,847 | Raw count data (no filtering) | Alpha diversity, DESeq2 analysis |
| **`ubiome_relative`** | 1,349 | 9,847 | Relative abundance (no filtering) | General compositional analysis |
| **`ubiome_relative_none`** | 659 | 9,847 | 0.1% prevalence filtered, no transformation | Filtered baseline analysis |
| **`ubiome_relative_clr`** | 659 | 9,847 | 0.1% prevalence filtered, CLR transformed | **Beta diversity, compositional analysis (RECOMMENDED)** |
| **`ubiome_relative_lognorm`** | 659 | 9,847 | 0.1% prevalence filtered, log-normal transformed | Standard statistical tests, linear modeling |

### **Key Validation Results**

✅ **Taxonomy Table**: 1,349 taxa with 6 taxonomic levels (Domain → Genus)
✅ **Sample Metadata**: 9,847 samples with **84 comprehensive variables** including survey weights
✅ **Sample Consistency**: All objects have identical sample names
✅ **Data Integrity**: Appropriate negative values only in transformed data (CLR/lognorm)
✅ **Survey Variables**: All essential NHANES design variables included
✅ **Factor Levels**: Proper R factor levels applied to categorical variables

### **Correct Taxa Counts Achieved**

The **union approach** successfully preserved all existing OTU data:

- **Unfiltered data**: 1,349 taxa (counts & relative)
- **Filtered data**: 659 taxa (union of 648 F + 409 G taxa, with 0-filling)
- **No data loss**: All existing taxa preserved through proper union operations

**Technical Details:**
- F tables: 648 taxa → expanded to 659 taxa (11 taxa filled with 0)
- G tables: 409 taxa → expanded to 659 taxa (250 taxa filled with 0)
- Result: Both cycles have identical 659 taxa structure

---

## 👥 **Complete Sample Metadata (84 Variables)**

The sample metadata now includes **all 84 common variables** from DEMO_F and DEMO_G tables instead of just the 23 that were previously selected. This provides access to the **complete demographic dataset** including:

### **Variable Categories:**

#### **🧑‍🤝‍🧑 Demographic Variables (32 variables)**
- **Core Demographics**: `Gender`, `Age`, `AgeGroup`, `AgeQuartile`, `RIAGENDR_01`
- **Ethnicity & Race**: `Ethnicity`, `RIDRETH1`, `RIDRETH3`, ethnic category variables
- **Citizenship**: `BornInUSA`, `US_Born`, `US_Citizen`, `DMDCITZN`, `DMDYRSUS`
- **Education**: `EducationLevel`, `Education_Level`, `EducationLevel_Simple`, `DMDEDUC2`, education categories
- **Income & Poverty**: `PIRcat`, `PIRQuartile`, `Ratio_Family_Income_Poverty`, `INDFMPIR`, income categories
- **Alternative Names**: Multiple naming conventions for compatibility across analysis pipelines

#### **🏠 Household Variables (8 variables)**
- **Household Structure**: `Household_Size`, `Family_Size`, `Household_Size_Factor`, `DMDHHSIZ`, `DMDFMSIZ`
- **Marital Status**: `Marital_Status`, `DMDMARTL`
- **Household Reference Person**: `HH_Reference_Gender`, `HH_Reference_Age`, `HH_Reference_Education_Level`, `HH_Reference_Marital_Status`
- **Income Details**: `Annual_Household_Income`, `Annual_Family_Income`

#### **Survey Design Variables (12 variables)**
- **Essential NHANES Variables**: `WTINT2YR`, `WTMEC2YR`, `SDMVPSU`, `SDMVSTRA`
- **Survey Structure**: `SDDSRVYR`, `RIDSTATR`, `RIDEXMON`
- **Interview Details**: `Interview_Language`, `Interview_Proxy`, `Interview_Interpreter`
- **Cycle Information**: `cycle`, `Release_Cycle`

#### **Other Variables (32 variables)**
- **Original NHANES IDs**: `SEQN`, `sample` (alternative naming)
- **Age Details**: `RIDAGEYR`, `RIDAGEMN`, continuous age variables
- **Additional Demographics**: Various derived and original NHANES demographic variables
- **Quality Control**: Interview and examination status variables

### **🔧 Technical Implementation**

```r
# Sample metadata creation (simplified view)
samples_df_nhanes <- demo_combined %>%
  rename(sample_id = SEQN) %>%
  # Keep ALL 84 variables (no selection/filtering)
  column_to_rownames("sample_id")

# Result: 9,847 samples × 84 variables
```

### **✅ Essential NHANES Survey Design Variables Verified**

All critical survey design variables are present for proper weighted analysis:

- ✅ `WTINT2YR` - Interview weights
- ✅ `WTMEC2YR` - MEC examination weights  
- ✅ `SDMVPSU` - Masked variance pseudo-PSU
- ✅ `SDMVSTRA` - Masked variance pseudo-stratum

---

## **Taxonomy Information**

### **Taxonomic Structure**

- **Levels**: Domain, Phylum, Class, Order, Family, Genus (6 levels)
- **Coverage**: Complete taxonomic annotation for all 1,349 taxa
- **Source**: SILVA 123 database annotations from NHANES metadata
- **Format**: Proper phyloseq `tax_table` format with consistent naming

### **🌳 Phylogenetic Tree**

- **Construction**: Taxonomy-based hierarchical tree
- **Tips**: 1,349 tips matching OTU table
- **Branch Lengths**: Uniform branch lengths for consistent analysis
- **Format**: Standard phylo object compatible with phyloseq

---

## 💾 **Output Files Created**

### **📁 Individual phyloseq Objects**
```
results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/
├── ubiome_counts.rds                    # Raw count data
├── ubiome_relative.rds                  # Relative abundance, no filtering  
├── ubiome_relative_none.rds             # Filtered, no transformation
├── ubiome_relative_clr.rds              # Filtered, CLR transformed
└── ubiome_relative_lognorm.rds          # Filtered, log-normal transformed
```

### **📦 Combined Object File**
```
nhanes_phyloseq_objects_all.rds         # All 5 objects in one list
```

### **🔧 Usage Examples**

#### **Loading Specific Objects:**
```r
# Load individual object
ubiome_clr <- readRDS("results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_clr.rds")

# Load all objects
all_objects <- readRDS("results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/nhanes_phyloseq_objects_all.rds")
ubiome_clr <- all_objects$relative_clr
```

#### **Basic Analysis Workflow:**
```r
# Alpha diversity (use counts data)
alpha_div <- estimate_richness(ubiome_counts)

# Beta diversity (use CLR transformed data - RECOMMENDED)
ordination <- ordinate(ubiome_relative_clr, method='PCoA', distance='euclidean')

# Differential abundance 
# - Use ubiome_counts for DESeq2
# - Use ubiome_relative_clr for compositional analysis
```

#### **Survey-Weighted Analysis:**
```r
# Extract survey design variables (all 84 variables available!)
sample_data <- data.frame(sample_data(ubiome_relative_clr))

# Essential survey variables
weights <- sample_data$WTMEC2YR
strata <- sample_data$SDMVSTRA  
clusters <- sample_data$SDMVPSU

# Create survey design object
library(survey)
svy_design <- svydesign(
  ids = ~clusters,
  strata = ~strata, 
  weights = ~weights,
  data = sample_data,
  nest = TRUE
)
```

---

## **Analysis Recommendations**

### **Use Case Guidelines**

| Analysis Type | Recommended Object | Rationale |
|---------------|-------------------|-----------|
| **Alpha Diversity** | `ubiome_counts` | Raw counts preserve count-based diversity metrics |
| **Beta Diversity** | `ubiome_relative_clr` | CLR transformation handles compositional nature |
| **Differential Abundance** | `ubiome_counts` (DESeq2) or `ubiome_relative_clr` (compositional) | Count-based for DESeq2, CLR for compositional methods |
| **Ordination** | `ubiome_relative_clr` | CLR transformation recommended for PCoA/NMDS |
| **Machine Learning** | `ubiome_relative_clr` or `ubiome_relative_lognorm` | Transformed data for ML algorithms |
| **Standard Statistics** | `ubiome_relative_lognorm` | Log-normal for linear modeling |
| **Baseline Comparison** | `ubiome_relative_none` | Filtered but untransformed data |

### **Transformation Guide**

- **`ubiome_counts`**: Original count data for count-based methods
- **`ubiome_relative`**: Basic relative abundance for general use
- **`ubiome_relative_none`**: Filtered baseline (0.1% prevalence threshold applied)
- **`ubiome_relative_clr`**: **RECOMMENDED** for most microbiome analyses (handles compositionality)
- **`ubiome_relative_lognorm`**: For methods requiring normal-like distributions

### **⚠️ Important Considerations**

1. **Survey Weights**: Always use survey weights for population inference
2. **Prevalence Filtering**: Filtered objects (659 taxa) exclude rare taxa for robust analysis
3. **Negative Values**: CLR and log-normal transformed data contain negative values (expected)
4. **Zero Handling**: Zero values preserved in all objects for proper analysis
5. **Metadata Completeness**: All 84 demographic variables available for comprehensive analysis

---

## 🔗 **Integration with Analysis Pipeline**

### **Workflow Integration**

```
Chapter 0: Data Transformation
     ↓
Chapter 1: Association Analysis  
     ↓
Chapter 2: Database Processing
     │
     ├── Step 1: Demographics Processing ✅
     ├── Step 2: Database Verification ✅  
     └── Step 3: phyloseq Creation ✅ ← 📍 COMPLETED
     ↓
Chapter 3: GOLD Database Integration
     ↓  
Chapter 4: Downstream Analysis
```

### **Next Steps After phyloseq Creation**

1. **Quality Control**: Validate phyloseq objects meet analysis requirements
2. **Exploratory Analysis**: Basic diversity and composition analysis  
3. **Advanced Analysis**: Integration with Chapters 3-4 pipelines
4. **Publication**: Generate publication-ready figures and statistics

---

## **Technical Validation Summary**

### **✅ Data Integrity Checks**

| Check | Status | Details |
|-------|--------|---------|
| **Sample Consistency** | ✅ PASS | All 5 objects have identical 9,847 sample names |
| **Taxa Structure** | ✅ PASS | Unfiltered: 1,349 taxa; Filtered: 659 taxa |
| **Metadata Completeness** | ✅ PASS | All 84 demographic variables included |
| **Survey Variables** | ✅ PASS | All essential NHANES design variables present |
| **Transformation Validation** | ✅ PASS | Appropriate negative values in CLR/lognorm only |
| **Factor Levels** | ✅ PASS | Categorical variables properly factored |
| **Union Approach** | ✅ PASS | No taxa data loss, proper 0-filling |

### **📏 Object Dimensions Summary**

```
phyloseq Object Summary:
  ubiome_counts      : 1349 taxa × 9847 samples × 84 variables
  ubiome_relative    : 1349 taxa × 9847 samples × 84 variables  
  ubiome_relative_none    : 659 taxa × 9847 samples × 84 variables
  ubiome_relative_clr     : 659 taxa × 9847 samples × 84 variables
  ubiome_relative_lognorm : 659 taxa × 9847 samples × 84 variables

Taxonomy: 6 levels (Domain → Genus)
🌳 Phylogeny: 1349 tips, taxonomy-based tree
👥 Metadata: 84 comprehensive demographic variables
Survey Design: Complete NHANES weights and design variables
```

---
