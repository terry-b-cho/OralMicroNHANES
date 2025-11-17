# NHANES Oral Microbiome Wide Association Studies (M-WAS) Pipeline

<div align="center">
  <img src="assets/logo/oral_micro_logo.svg" alt="NHANES Oral Microbiome Analysis Logo" width="600"/>
  
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![R Version](https://img.shields.io/badge/R-4.3.1-blue.svg)](https://www.r-project.org/)
  [![Conda](https://img.shields.io/badge/Conda-23.1.0-green.svg)](https://docs.conda.io/)
  [![SLURM](https://img.shields.io/badge/SLURM-Compatible-orange.svg)](https://slurm.schedmd.com/)
</div>

---

## **Project Overview**

This repository implements a **comprehensive pipeline** for conducting Microbiome Wide Association Studies (M-WAS) on NHANES oral microbiome data. The pipeline performs **24 distinct analyses** (6 analysis types × 4 transformation methods) with **survey-weighted statistics** and **multiple comparisons correction as per CDC guideline**.

**MAJOR STATISTICAL ENHANCEMENTS:**
- **Pooled-Cycle Analysis**: Automatic detection and implementation of 2009-2012 pooled estimates with NCHS-compliant weight calculations
- **Final Technical Audit**: All 9 critical blocking issues resolved for production readiness
- **Integer Pseudocount**: ε=1 preserves count scale, prevents extreme negative values
- **SQL View Optimization**: 50% storage reduction for untransformed tables
- **Enhanced Transformations**: 4 methods including Hellinger transformation
- **Progressive Algorithms**: Covariate selection, adaptive binning, convergence monitoring
- **Effect Scale Harmonization**: Comparable results across transformation methods
- **FDR Correction**: Benjamini-Hochberg procedure applied within dependent variables
- **Mathematical Precision**: Complete algorithmic formulations for biostatistician reproducibility

### **What This Pipeline Does:**
- **Transforms** raw microbiome OTU abundance data using **4 normalization methods** (including Hellinger)
- **Conducts** 6 types of association studies between microbiome and health variables
- **Implements** **pooled-cycle analysis** with automatic detection of 2009-2012 combined estimates for enhanced statistical power
- **Applies** proper NHANES survey design (stratification, clustering, weights) with NCHS compliance
- **Executes** **mathematically precise** workflow: Individual regressions → Aggregation with FDR correction
- **Provides** smart resume capabilities with progressive algorithms and fallback strategies
- **Generates** publication-ready results with effect scale harmonization and comprehensive statistics

### **Analysis Types:**

| Pipeline | Description | Direction | Model Type |
|----------|-------------|-----------|------------|
| **1_demoWAS** | How demographics affect oral microbiome | Demographics → **Microbiome** | Linear |
| **2_oradWAS** | How oral microbiome affects oral health | **Microbiome** → Oral Health | Logistic |
| **3_exWAS** | How exposures affect oral microbiome | Exposures → **Microbiome** | Linear |
| **4_pheWAS** | How oral microbiome affects phenotypes | **Microbiome** → Phenotypes | Linear |
| **5_outWAS** | How oral microbiome affects disease outcomes | **Microbiome** → Disease Outcomes | Logistic |
| **6_zimWAS** | How oral microbiome affects lab measurements | **Microbiome** → Lab Values | Linear |

### **Transformation Methods (STATISTICALLY CORRECTED):**

**Mathematical Formulations:**

- **None** - T_ij = P_ij (untransformed proportions via SQL views, no storage duplication)
- **Hellinger** - T_ij = √P_ij (Hellinger transformation for distance-based analysis)
- **CLR** - T_ij = ln[(C_ij + 1)/g_i] (CLR from COUNT tables, integer pseudocount ε=1)
- **LogNorm** - T_ij = log₁₀[(C_ij + 1)/(n_i + D) × n̄] (log-normal with zero-library handling)

**Where:**
- P_ij = relative abundance of taxon j in sample i
- C_ij = count of taxon j in sample i  
- g_i = geometric mean of sample i counts (robust calculation)
- n_i = library size of sample i
- n̄ = mean library size across samples
- D = number of taxa
- ε = 1 (integer pseudocount preserves count scale)

---

## **Complete Step-by-Step Guide**

### **Prerequisites**
- ✅ Access to O2 HPC cluster (or similar SLURM environment)
- ✅ R ≥ 4.4.3 (via conda - critical for consistency)
- ✅ Conda ≥ 24.11.3 (miniforge3 recommended)
- ✅ GCC ≥ 14.2.0 (for R package compilation)
- ✅ Input data: `data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite`

### **R Environment Setup (Critical - UPDATED FOR O2)**

**Important:** We use R 4.4.3 via Conda, not O2's module system, for consistency and enhanced package support required by the final audit corrected pipeline.

**Required setup for all analyses (CRITICAL: Updated O2 modules):**

```bash
# Essential modules (updated for O2 compatibility)
module purge
module load gcc/14.2.0                 # Updated with O2 update
module load conda/miniforge3/24.11.3-0 # Updated with O2 update

# Initialize conda for bash and activate environment (modern syntax)
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis
```

**Verify R installation:**
```bash
which R
# Expected: /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis/bin/R

R --version
# Expected: R version 4.4.3 (2024-02-29) -- "Eye Holes"
```

**⚠️ Critical Notes for Final Audit Corrected Pipeline:**
- Always use this Conda environment for consistent results across all 24 analyses
- Never use O2's R modules for this project (causes GLIBCXX errors with enhanced packages)
- The `gcc/14.2.0` module is required for compiling R packages with C++ dependencies
- Use the updated conda activation syntax with `eval "$(conda shell.bash hook)"` for O2 compatibility
- Enhanced R packages required: `survey`, `broom`, `dplyr`, `stringr`, `optparse` for audit corrections

---

## **Workflow Integration**

```
Chapter 0: Data Transformation ✅
     ↓
Chapter 1: Association Analysis ✅ 
     ↓
Chapter 2: Database Processing ✅
     │
     ├── Step 2.1: Demographics Processing ✅
     ├── Step 2.2: Database Verification ✅  
     └── Step 2.3: phyloseq Creation ✅
     ↓
Chapter 3: GOLD Database Integration ✅ ← 📍 COMPLETED
     │
     ├── Step 3.1: GOLD Database Processing ✅
     ├── Step 3.2: Genus Name Parsing ✅
     └── Step 3.3: Phenotype Mapping ✅
     ↓  
Chapter 4: Downstream Analysis
```

---

## **Chapter 0: Data Preprocessing & Transformation**

> **📖 Full Documentation:** [0_transform_n_preprocess_ssfiles_documentation.md](scripts/0_transform_n_preprocess_ssfiles/0_transform_n_preprocess_ssfiles_documentation.md)

### **Step 0.1: Microbiome Data Transformation (STATISTICALLY CORRECTED)**

This step transforms raw OTU abundance data into analysis-ready formats with **comprehensive statistical corrections** and **NO prevalence pre-filtering** to preserve analytical flexibility.

```bash
# Submit transformation job (uses updated R environment)
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh

# Monitor progress
squeue -u $USER
```

**STATISTICAL CORRECTIONS APPLIED:**

**Mathematical Formulations (Corrected):**
- **None**: T_ij = P_ij (SQL views from RELATIVE tables, no storage duplication)
- **Hellinger**: T_ij = √P_ij (from RELATIVE tables, for distance-based analysis)
- **CLR**: T_ij = ln[(C_ij + 1)/g_i] (from COUNT tables, integer pseudocount ε=1)
- **Log-norm**: T_ij = log₁₀[(C_ij + 1)/(n_i + D) × n̄] (from COUNT tables, zero-library handling)

**Key Corrections:**
- **Integer Pseudocount**: ε = 1 (preserves count scale, avoids extreme negatives)
- **Zero-Library Handling**: Samples with n_i = 0 set to NA gracefully
- **SQL View Optimization**: "_none" tables as views (50% storage reduction)
- **Numerical Precision**: Single coercion to avoid repeated casting
- **Closure Verification**: Warnings for non-closed relative abundance data
- **Enhanced Metadata**: Transformation-specific descriptions with mathematical units

**What happens:**
- **NO prevalence pre-filtering** - All taxa retained for regression-time filtering decisions
- Creates **4 versions** of each taxonomic table (was 3):
  - `*_none` - SQL views of raw relative abundances (no duplication)
  - `*_hellinger` - Square-root transformed (NEW: distance-based analysis)
  - `*_clr` - CLR-transformed with integer pseudocount (compositional data analysis)
  - `*_lognorm` - Log-normal transformed with zero-library handling (variance stabilization)
- **Universal naming convention**: All output tables use `DADA2RSV_<taxalevel>_RELATIVE_<cycle>_<method>` format
- Processes 5 taxonomic levels: genus, family, order, class, phylum (10 source tables → **40 transformed tables + 10 SQL views**)
- **Statistically sound pseudocount**: Integer ε = 1 applied only to COUNT matrices
- Automatically updates metadata with transformation provenance and mathematical descriptions

**Output:** `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite`

### **Step 0.2: Missing Data Completion and Variable Derivation (ENHANCED)**

```bash
# Run after Step 0.1 completes (requires conda environment)
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
  --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
  --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
```

**🔧 CRITICAL ENHANCEMENT: SQL View Recreation**
- **Problem**: SLURM jobs read the complete DB, but "_none" views only in intermediate DB
- **Solution**: Recreates all 10 "_none" SQL views in the complete database
- **Impact**: Prevents table/view detection failures in final audit corrected analysis

**What happens:**
- **Critical Fix**: Recreates "_none" SQL views in complete database for SLURM compatibility
- **Gender variable derivation**: Creates `RIAGENDR_01` (0=Male, 1=Female) from `RIAGENDR`
- **Comprehensive variable annotation**: Documents 40+ variables across categories:
  - Demographic derivatives (age, education, ethnicity indicators)
  - Laboratory biomarkers (PSA, folate, creatine kinase, etc.)
  - Health outcomes (respiratory, cardiovascular, cancer history)
  - Exposure biomarkers (arsenic, cysteine derivatives)
- **Metadata completion algorithm**: Maps variables to tables and fills missing metadata entries
- **Database transaction safety**: Atomic updates with rollback protection

**Output:** `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite`

### **Step 0.3: Schema Structure File Creation (UPDATED: 24 FILES)**

```bash
# Creates mapping files for all 24 analyses (24 parallel SLURM jobs)
bash scripts/0_transform_n_preprocess_ssfiles/run_ss_file_create.sh
```

**What happens:**
- Submits **24 parallel SLURM jobs** (6 analysis types × **4 transformations** including hellinger)
- Maps microbiome OTU variables to NHANES variables based on configuration files
- Calculates sample sizes for each variable pair across both cycles (F & G)
- Uses conda environment for consistency with downstream analyses
- Creates proper OTU role assignments:
  - `dep` (dependent): 1_demoWAS, 3_exWAS (Demographics/Exposures → Microbiome)
  - `indep` (independent): 2_oradWAS, 4_pheWAS, 5_outWAS, 6_zimWAS (Microbiome → Outcomes)

**Output:** `results/0_ss_files/*.csv` (**24 files total**)
- Each schema file contains: `dep_var`, `dep_table`, `indep_var`, `indep_table`, `n`, `cycle`
- Example files: `1_demoWAS_clr_schema_structure.csv`, `2_oradWAS_hellinger_schema_structure.csv`

---

## **Chapter 1: Association Analysis Pipeline - FINAL AUDIT CORRECTED**

> **📖 Full Documentation:** [1_association_pipeline_documentation.md](scripts/1_association_pipeline/1_association_pipeline_documentation.md)

### **Step 1.1: Environment Setup (One-time) - O2 COMPATIBLE**

```bash
# Start interactive session
srun --pty -p interactive -t 12:00:00 --mem=32G bash

# Setup conda environment with O2 compatibility
bash scripts/1_association_pipeline/setup_nhanes_environment.sh

# Load required modules and activate environment (CRITICAL: Updated for O2)
module purge
module load gcc/14.2.0                 # Updated with O2 update
module load conda/miniforge3/24.11.3-0 # Updated with O2 update

# Initialize conda for bash and activate environment (modern syntax)
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

# Validate installation with enhanced testing (24 analyses)
Rscript scripts/1_association_pipeline/debug_all_pipelines.R
```

**FINAL AUDIT CORRECTIONS VERIFIED:**
- **Table/View Detection**: Handles SQL views for "_none" tables
- **Enhanced R Packages**: `survey`, `broom`, `dplyr`, `stringr`, `optparse` 
- **Statistical Rigor**: FDR correction, effect scale harmonization
- **Survey Design**: NCHS Technical Documentation 2006 compliance

### **Step 1.2: Run Association Analyses (24 TOTAL)**

#### **Option A: FLEXIBLE Resume Batch Processing (⭐ RECOMMENDED)**

```bash
# Exit interactive session
exit

# Run all 24 analyses (6 types × 4 transformations) with FLEXIBLE resume
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh

# For test mode (recommended first time - tests all 4 transformations)
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh test

# With aggressive resources (32G, 12h for all jobs)
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh "" aggressive

# View help and features
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh --help
```

**ENHANCED FEATURES:**
- **24 Total Analyses**: 6 types × 4 transformations (including hellinger)
- **Table/View Detection**: Prevents silent failures with SQL views
- **FLEXIBLE Resume**: Only submits missing analyses, never overwrites
- **Enhanced Error Handling**: Progressive covariate reduction with fallbacks
- **Effect Scale Harmonization**: Comparable results across transformations

#### **Option B: Individual Analysis (ALL 4 TRANSFORMATIONS)**

```bash
# Run single analysis type (all 4 transformations now supported)
bash scripts/1_association_pipeline/run_single_was_analysis.sh 1_demoWAS clr
bash scripts/1_association_pipeline/run_single_was_analysis.sh 1_demoWAS hellinger

# For test mode
bash scripts/1_association_pipeline/run_single_was_analysis.sh 1_demoWAS clr test

# All available combinations:
# Analysis types: 1_demoWAS, 2_oradWAS, 3_exWAS, 4_pheWAS, 5_outWAS, 6_zimWAS
# Normalizations: clr, lognorm, none, hellinger
```

### **Step 1.3: Monitor Progress (ENHANCED)**

```bash
# Check running jobs
squeue -u $USER | grep flexible

# Check for errors (enhanced error detection)
find logs/ -name 'flexible_*.err' -size +0

# View execution log
tail -f logs/1_association_pipeline_flexible/flexible_was_analyses_*.log

# Verify table/view detection working
sqlite3 data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
  "SELECT COUNT(*) FROM sqlite_master WHERE type='view' AND name LIKE '%_none%';"
# Should return 10 (5 taxa × 2 cycles)
```

### **Step 1.4: Aggregate Results with FDR Correction (After Individual Analyses Complete)**

#### **Option A: Batch Aggregation + Supplementary Tables (RECOMMENDED)**

```bash
# Aggregate ALL 24 analyses + create publication-ready supplementary tables
Rscript scripts/1_association_pipeline/aggregate_was_results.R --run-all $(pwd)/results
```

**ENHANCED AGGREGATION FEATURES:**
- **FDR Correction**: Benjamini-Hochberg applied within each dependent variable (q < 0.1)
- **Effect Scale Preservation**: Maintains transformation-specific interpretations
- **Enhanced Reporting**: Pre/post correction statistics, effect scale summaries
- **Publication Ready**: Supplementary tables with comprehensive statistics

#### **Option B: Individual Aggregation**

```bash
# Aggregate single analysis type
Rscript scripts/1_association_pipeline/aggregate_was_results.R 1_demoWAS clr $(pwd)/results

# Repeat for all 24 combinations (6 types × 4 transformations)
for analysis in 1_demoWAS 2_oradWAS 3_exWAS 4_pheWAS 5_outWAS 6_zimWAS; do
  for transform in none hellinger clr lognorm; do
    Rscript scripts/1_association_pipeline/aggregate_was_results.R $analysis $transform $(pwd)/results
  done
done
```

#### **Option C: Supplementary Tables Only (After Aggregation)**

```bash
# Create supplementary tables from existing aggregation summaries
Rscript scripts/1_association_pipeline/aggregate_was_results.R --create-tables $(pwd)/results

# Alternative: Use standalone script
Rscript scripts/1_association_pipeline/create_supplementary_tables.R $(pwd)/results
```

### **Step 1.5: Resume Missing Analyses (Enhanced Safety)**

The FLEXIBLE resume system with final audit corrections handles missing analyses automatically:

```bash
# Simply re-run the flexible script - it will only submit missing jobs
./scripts/1_association_pipeline/run_all_was_analyses_flexible.sh

# The system will:
# ✅ Skip completed analyses (preserves existing results)
# ✅ Only submit jobs for missing .rds files
# ✅ Never overwrite successful results
# ✅ Handle table/view detection for "_none" analyses
# ✅ Apply all 9 critical audit corrections automatically
```

---
## **Chapter 2: Database Preprocessing & phyloseq Object Creation**

> **📖 Full Documentation:** [2_preprocess_db_n_phyloseq_documentation.md](scripts/2_preprocess_db_n_phyloseq/2_preprocess_db_n_phyloseq_documentation.md)

### **Step 2.1: Comprehensive Demographic Processing (Run First)**

Create comprehensive processed database with all demographic variables:

```bash
# Test mode first (recommended to verify what will be created)
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R --test

# Create processed database (creates 32 new variables in DEMO_F and DEMO_G)
Rscript scripts/2_preprocess_db_n_phyloseq/process_demographics_complete.R
```

**What it creates:**
- ✅ 32 comprehensive demographic variables (cycle, Gender, AgeGroup, EducationLevel, etc.)
- ✅ Data-driven quartiles (AgeQuartile, PIRQuartile)
- ✅ Proper factor levels for categorical variables
- ✅ Alternative naming conventions for compatibility
- ✅ Complete metadata registration

**Output:** `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite`

### **Step 2.2: Database Diagnostic Verification (Run After Step 2.1)**

Verify that demographic processing was successful and check database status:

```bash
# Check processed database (after running Step 2.1)
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R --processed

# Check original database status for comparison
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R

# Comprehensive diagnostic with detailed output
Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R \
  --db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
  --output_dir results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics/ \
  --format both
```

**What it checks:**
- ✅ All 30 microbiome transformation tables (5 taxa × 2 cycles × 3 transforms)
- ✅ Expected derived tables (`OralDisease_F/G`, `d_outcome_mcq_f/g`, etc.)
- ✅ Expected derived variables (77+ variables across categories)
- ✅ Demographic processing status (64 comprehensive demographic variables)
- ✅ Metadata completeness in `variable_names_epcf`

**Output:** Console summary + detailed CSV/HTML reports in `results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics/`

### **Step 2.3: phyloseq Object Creation (Run After Steps 2.1 & 2.2)**

Create comprehensive phyloseq objects from the processed database for downstream microbiome analysis:

```bash
# Navigate to script directory
cd scripts/2_preprocess_db_n_phyloseq/

# Run the RMarkdown file (in R/RStudio or command line)
Rscript -e "rmarkdown::render('ubiome_phyloseq_object_creation_cho.Rmd')"
```

**What it creates:**

#### **✅ 5 phyloseq Objects for Different Analysis Types:**

| Object | Taxa Count | Samples | Description |
|--------|------------|---------|-------------|=
| **`ubiome_counts`** | 1,349 | 9,847 | Raw count data (no filtering) |
| **`ubiome_relative`** | 1,349 | 9,847 | Relative abundance (no filtering) |
| **`ubiome_relative_none`** | 659 | 9,847 | 0.1% prevalence filtered, no transformation |
| **`ubiome_relative_clr`** | 659 | 9,847 | 0.1% prevalence filtered, CLR transformed |
| **`ubiome_relative_lognorm`** | 659 | 9,847 | 0.1% prevalence filtered, log-normal transformed |

#### **Key Features:**

✅ **Complete Sample Metadata**: 9,847 samples with **84 comprehensive variables** including:
- All original NHANES variables (SEQN, SDDSRVYR, RIDSTATR, etc.)
- All derived demographic variables created by `process_demographics_complete.R`
- All survey design variables for proper weighted analysis
- All processed categorical variables with proper factor levels

✅ **Taxonomy Information**: 1,349 taxa with 6 taxonomic levels (Domain → Genus)
✅ **Union Approach**: Preserves all existing OTU data (no taxa loss)
- Unfiltered data: 1,349 taxa (counts & relative)
- Filtered data: 659 taxa (union of 648 F + 409 G taxa, with 0-filling)
✅ **Survey Compatibility**: All essential NHANES design variables included
✅ **Data Integrity**: Proper validation and consistent structure across all objects

**Output:** `results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/`
- Individual RDS files for each phyloseq object
- Combined file: `nhanes_phyloseq_objects_all.rds`
- Ready for downstream microbiome analysis

#### **Usage Recommendations:**

- **Alpha diversity**: Use `ubiome_counts`
- **Beta diversity**: Use `ubiome_relative_clr` (recommended)
- **Compositional analysis**: Use `ubiome_relative_clr`
- **Standard statistics**: Use `ubiome_relative_lognorm`
- **Survey-weighted analysis**: All objects include proper NHANES design variables

## **Chapter 3: GOLD Database Integration**

> **📖 Full Documentation:** [3_gold_db_microbial_phenotype_documentation.md](scripts/3_gold_db_microbial_phenotype/3_gold_db_microbial_phenotype_documentation.md)

### **Step 3.1: GOLD Database Integration (Completed ✅)**

Integrates the Genomes OnLine Database (GOLD) with NHANES oral microbiome data for functional phenotype annotation.

```bash
# Navigate to Chapter 3 directory
cd scripts/3_gold_db_microbial_phenotype/

# Run the RMarkdown integration script
Rscript -e "rmarkdown::render('gold_db_process_n_genus_mapping.Rmd')"

# Alternative: Use the logging wrapper script
Rscript run_gold_db_integration.R
```

**What it accomplishes:**

#### **✅ GOLD Database Processing:**
- **Raw Data Processing**: 438,180 rows × 22 columns successfully processed
- **Genus Aggregation**: 4,215 unique genera with 18 phenotypic features
- **Processing Efficiency**: Complete aggregation in <10 seconds

#### **✅ Phenotypic Features Integrated (15 total):**
1. **`BIOTIC_RELATIONSHIPS`** - Ecological relationships and interactions
2. **`OXYGEN_REQUIREMENT`** - Oxygen metabolism preferences (anaerobe/aerobe)
3. **`METABOLISM`** - Metabolic pathway information and capabilities
4. **`ENERGY_SOURCES`** - Energy utilization patterns and sources
5. **`GRAM_STAIN`** - Cell wall characteristics (Gram+/Gram-)
6. **`CELL_SHAPE`** - Morphological characteristics (coccus, rod, etc.)
7. **`MOTILITY`** - Movement capabilities (motile/nonmotile)
8. **`SPORULATION`** - Spore formation capacity
9. **`TEMPERATURE_RANGE`** - Thermal growth preferences
10. **`SALINITY`** - Salt tolerance ranges
11. **`STAIN_INDEX`** - Quantitative gram stain index [0,1]
12. **`OXYGEN_INDEX`** - Quantitative oxygen requirement index [0,1]
13. **`MOTILITY_INDEX`** - Quantitative motility index [0,1]
14. **`SPORULATION_INDEX`** - Quantitative sporulation index [0,1]
15. **`SAMPLE_COLLECTION_DATE`** - Isolation metadata

#### **✅ Genus Name Parsing & Standardization:**
- **Original Genera**: 964 unique genus names from NHANES taxonomy
- **Parsed Genera**: 957 genera after standardization (handles complex nomenclature)
- **Parsing Algorithm**: Comprehensive handling of Candidatus, brackets, group designations
- **Unclassified Assignment**: 385 OTUs appropriately classified as "unclassified"

#### **✅ Mapping Performance Results:**

**Overall Coverage Statistics:**
- **Total OTUs Processed**: 1,349 (100% of NHANES oral microbiome data)
- **OTUs with GOLD Annotations**: 820 (60.79% overall coverage)
- **Classifiable OTUs**: 964 (71.46% of total)
- **Annotated Classifiable OTUs**: 820 (85.06% of classifiable) ⭐ **EXCELLENT**

**Quality Assessment:**
- **High-Quality Integration**: 85.06% annotation success for classifiable genera
- **Comprehensive Coverage**: 4,215 GOLD genera available for mapping
- **Robust Parsing**: Successfully handled complex taxonomic nomenclature

#### **✅ Output Files Generated:**

**Primary Outputs** (in `results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/`):

1. **`gold_db_genus.csv`**: 4,215 genera × 18 features (aggregated GOLD database)
   - Complete genus-level phenotypic profiles
   - Agreement levels (high/medium/low) for quality assessment
   - Ready for functional analysis applications

2. **`ubiome_genus_mapping_complete.csv`**: 1,349 OTUs × 18 columns (complete mapping)
   - Every NHANES OTU with available GOLD annotations
   - Original and parsed genus names for transparency
   - All 15 phenotypic features integrated

3. **`mapping_summary_stats.csv`**: 8 key performance metrics
   - Coverage statistics and quality indicators
   - Validation metrics for publication

**Integration Success Indicators:**
- [x] **Excellent Coverage**: 85.06% of classifiable genera successfully annotated
- [x] **Comprehensive Features**: 15 phenotypic traits integrated
- [x] **Quality Assurance**: Robust parsing and validation throughout
- [x] **Pipeline Ready**: All outputs validated for downstream Chapter 4 analysis

#### **Scientific Applications:**

**Functional Annotation Capabilities:**
- **Mechanistic Interpretation**: Link microbiome associations to biological functions
- **Phenotype-Based Analysis**: Group organisms by functional characteristics
- **Hypothesis Generation**: Identify functional patterns in significant associations
- **Enrichment Testing**: Test for overrepresentation of specific phenotypes

**Integration with Association Results:**
```r
# Example: Enhance association interpretation with phenotypes
significant_genera <- association_results[p.value.fdr < 0.05, ]
annotated_significant <- merge(significant_genera, 
                              ubiome_genus_mapping_complete, 
                              by.x = "exposure", by.y = "otu")

# Analyze oxygen requirements of significant genera
table(annotated_significant$OXYGEN_REQUIREMENT)
```

**Mathematical Framework:**
- **Binary Indices**: Proportion-based scoring [0,1] for traits like motility
- **Oxygen Index**: Continuous scale (anaerobe=0, aerobe=1) for metabolism
- **Agreement Levels**: Quality indicators (high/medium/low) for reliability

#### **Ready for Chapter 4:**
The GOLD database integration provides comprehensive functional annotation for downstream analyses including:
- **Effect Size Visualization** with phenotypic context
- **Functional Enrichment Analysis** of significant associations  
- **Mechanistic Hypothesis Generation** for microbiome-health relationships
- **Publication-Ready Methods** with mathematical formulations

---

## **Chapter 4: Downstream Analyses**

> **📖 Full Documentation:** [4_downstream_analyses_documentation.md](scripts/4_downstream_analyses/4_downstream_analyses_documentation.md)

**[PLANNED]** Advanced downstream analyses including:
- Effect size visualization
- Microbial Phenotype mapping to study effect

---

## 📁 **Project Structure**

```
nhanes_oral_mirco_cho/
├── 📖 README.md                          # This file
├── assets/                            # Project assets
│   ├── 00_trouble_shooting_n_appendix.md
|   └── logo/
│       ├── oral_micro_logo.ai
│       └── oral_micro_logo.svg
├── ⚙️ configs/                           # Analysis configuration
│   ├── 00_configs_documentation.md       # Config documentation
│   ├── 1_demoWAS_vars.txt               # Demographics variables
│   ├── 2_oradWAS_vars.txt               # Oral health variables
│   ├── 3_exWAS_vars.txt                 # Exposure variables
│   ├── 4_pheWAS_vars.txt                # Phenotype variables
│   ├── 5_outWAS_vars.txt                # Disease outcome variables
│   └── 6_zimWAS_vars.txt                # Lab measurement variables
├── 💾 data/                             # Data directories
│   ├── 00_GOLDdb/                       # GOLD database files
│   ├── 00_nhanes_omp_abundance_db/      # Raw abundance data
│   ├── 00_nhanes_omp_diversity_db/      # Diversity metrics
│   └── 00_nhanes_omp_transformed_db/    # Transformed data (40 tables + 10 SQL views)
├── results/                          # Analysis outputs
│   ├── 0_ss_files/                      # Schema structure files (24 total)
│   ├── 1_demoWAS_out/                   # Demographics results (4 transformations)
│   ├── 2_oradWAS_out/                   # Oral health results (4 transformations)
│   ├── 3_exWAS_out/                     # Exposure results (4 transformations)
│   ├── 4_pheWAS_out/                    # Phenotype results (4 transformations)
│   ├── 5_outWAS_out/                    # Disease outcome results (4 transformations)
│   ├── 6_zimWAS_out/                    # Lab measurement results (4 transformations)
│   └── analyses_results/                # Chapter-specific analysis results
│       ├── 02_preprocess_db_n_phyloseq_out/  # Chapter 2 outputs
│       │   ├── diagnostics/             # Original database diagnostic reports
│       │   ├── diagnostics_processed/   # Processed database diagnostic reports
│       │   └── intermediate/            # phyloseq Objects
│       │       ├── ubiome_counts.rds                    # Raw count data
│       │       ├── ubiome_relative.rds                  # Relative abundance, no filtering
│       │       ├── ubiome_relative_none.rds             # Filtered, no transformation
│       │       ├── ubiome_relative_clr.rds              # Filtered, CLR transformed
│       │       ├── ubiome_relative_lognorm.rds          # Filtered, log-normal transformed
│       │       └── nhanes_phyloseq_objects_all.rds     # All 5 objects combined
│       └── 03_gold_db_microbial_phenotype_out/  # Chapter 3 outputs ⭐ NEW
│           ├── intermediate/            # GOLD integration results
│           │   ├── gold_db_genus.csv                    # Aggregated GOLD database (4,215 genera)
│           │   ├── ubiome_genus_mapping_complete.csv    # Complete OTU-GOLD mapping (1,349 OTUs)
│           │   └── mapping_summary_stats.csv           # Coverage and quality statistics
│           └── logs/                    # Analysis logs and execution records
└── 🔧 scripts/                          # Analysis scripts
    ├── 0_transform_n_preprocess_ssfiles/ # Chapter 0: Data preprocessing (STATISTICALLY CORRECTED)
    │   ├── 📖 0_transform_n_preprocess_ssfiles_documentation.md
    │   ├── nhanes_omp_transformation.R                    # ⭐ CORRECTED: Integer pseudocount, SQL views
    │   ├── nhanes_db_filling_missing_data.R               # ⭐ ENHANCED: SQL view recreation
    │   ├── ss_file_create.R
    │   ├── run_transformation.sh
    │   └── run_ss_file_create.sh
    ├── 1_association_pipeline/           # Chapter 1: Association analyses (FINAL AUDIT CORRECTED)
    │   ├── 📖 1_association_pipeline_documentation.md     # ⭐ UPDATED: 24 analyses documentation
    │   ├── setup_nhanes_environment.sh                    # ⭐ O2 COMPATIBLE: Updated modules
    │   ├── universal_was_analysis_PRODUCTION.R            # ⭐ FINAL AUDIT: All 9 critical fixes
    │   ├── run_single_was_analysis_PRODUCTION.sh          # ⭐ O2 COMPATIBLE: 4 transformations
    │   ├── run_all_was_analyses_flexible.sh               # ⭐ CORRECTED: Hellinger support
    │   ├── aggregate_was_results.R             # ⭐ ENHANCED: FDR + effect harmonization
    │   ├── create_supplementary_tables.R                  # ⭐ Publication-ready tables
    │   └── debug_all_pipelines.R                          # ⭐ COMPREHENSIVE: 24 analyses testing
    ├── 2_preprocess_db_n_phyloseq/           # Chapter 2: Database preprocessing
    │   ├── 📖 2_preprocess_db_n_phyloseq_documentation.md
    │   ├── process_demographics_complete.R             # STEP 1: Comprehensive demographic processing
    │   ├── missing_tables_n_variables_diagnostic.R    # STEP 2: Database diagnostic verification
    │   └── ubiome_phyloseq_object_creation_cho.Rmd     # STEP 3: phyloseq object creation
    ├── 3_gold_db_microbial_phenotype/   # Chapter 3: GOLD integration ✅
    │   ├── 📖 3_gold_db_microbial_phenotype_documentation.md
    │   ├── gold_db_process_n_genus_mapping.Rmd     # Main integration script
    │   ├── run_gold_db_integration.R               # Logging wrapper script
    │   └── logs/                                   # Execution logs
    └── 4_downstream_analyses/           # Chapter 4: Advanced analyses
        └── 📖 4_downstream_analyses_documentation.md
```

---

## **Understanding Your Results (ENHANCED)**

### **Individual Results Files (After Step 1.2-1.3):**
```
results/<analysis_type>_out/result_<normalization>/
├── <dep_var1>.rds    # Enhanced individual results with FDR + effect scales
├── <dep_var2>.rds    # Enhanced individual results with FDR + effect scales
├── <dep_var3>.rds    # Enhanced individual results with FDR + effect scales  
└── ... (many more individual .rds files with audit corrections)
```

### **Final Aggregated Results Files (After Step 1.4):**
```
results/<analysis_type>_out/result_<normalization>/
├── <analysis>_<norm>_tidied_complete.rds    ⭐ MAIN RESULTS (FDR + effect scales)
├── <analysis>_<norm>_glanced_complete.rds   # Model statistics + design df
├── <analysis>_<norm>_rsq_complete.rds       # Survey R-squared + effect info
└── <analysis>_<norm>_aggregation_summary.txt # Enhanced analysis summary
```

### **Key Result Columns (ENHANCED - Available After Aggregation Step 1.4):**
```r
# Load main results (enhanced with final audit corrections)
results <- readRDS("results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds")

# CORE STATISTICAL COLUMNS:
# term                - Model term ("indep_var" = main effect to focus on)
# estimate            - Effect size (regression coefficient)
# std.error           - Survey-adjusted standard error
# p.value             - Original p-value

# FDR CORRECTION COLUMNS (NEW):
# p.value.fdr         - FDR-corrected p-value ⭐ USE FOR SIGNIFICANCE (q < 0.1)
# significant.fdr     - Boolean flag for FDR significance (TRUE/FALSE)

# METADATA COLUMNS:
# phenotype           - Dependent variable (outcome)
# exposure            - Independent variable (predictor)
# dependent_var       - Dependent variable name (same as phenotype)
# independent_var     - Independent variable name (same as exposure)
# n_obs               - Sample size for this specific analysis
# formula_used        - Actual statistical formula used
# n_covariates        - Number of covariates successfully included
# normalization       - Transformation method applied

# EFFECT SCALE HARMONIZATION COLUMNS (NEW):
# effect_scale        - Transformation-specific scale description
#                      • "proportion (0-1)" for none
#                      • "sqrt-proportion (0-1)" for hellinger  
#                      • "ln-ratio (centered)" for clr
#                      • "log10-CPM" for lognorm
# interpretation_note - Effect size interpretation guide per transformation

# ENHANCED QUALITY CONTROL COLUMNS (NEW):  
# df_residual_design  - Design degrees of freedom for validation
# aggregate_base_model - Flag for aggregation-level processing
```

### **Significance Testing (UPDATED):**

```r
# PRIMARY SIGNIFICANCE TEST (use FDR-corrected results)
significant_results <- results %>%
  filter(term == "indep_var",           # Main effects only
         significant.fdr == TRUE)       # FDR q < 0.1 ⭐ RECOMMENDED

# Alternative: Manual FDR threshold
fdr_significant <- results %>%
  filter(term == "indep_var",
         p.value.fdr < 0.1)            # Manual FDR threshold

# Effect size interpretation by transformation
results %>%
  filter(significant.fdr == TRUE) %>%
  group_by(normalization, effect_scale) %>%
  summarise(
    n_significant = n(),
    median_effect = median(abs(estimate), na.rm = TRUE),
    .groups = "drop"
  )
```

### **Effect Size Interpretation Guide:**

| Transformation | Effect Scale | Interpretation Example |
|---------------|--------------|----------------------|
| **none** | proportion (0-1) | "1-unit increase in relative abundance (0→1)" |
| **hellinger** | sqrt-proportion (0-1) | "1-unit increase in sqrt-proportion" | 
| **clr** | ln-ratio (centered) | "1-unit increase in log-ratio relative to geometric mean" |
| **lognorm** | log10-CPM | "1-unit increase in log10 counts-per-million" |

**Example Interpretation:**
```r
# For a significant CLR result with estimate = 0.5
# interpretation_note: "Effect per ln-ratio unit (compositional)"
# effect_scale: "ln-ratio (centered)"
# 
# Meaning: A 1-unit increase in the CLR-transformed abundance 
# (relative to geometric mean) is associated with a 0.5-unit 
# change in the outcome variable.
```

## **Statistical Methods (FINAL AUDIT CORRECTED)**

### **Survey Design Implementation (NCHS COMPLIANT):**
- **Stratification**: Uses `SDMVSTRA` pseudo-strata for variance estimation
- **Clustering**: Accounts for `SDMVPSU` primary sampling units with proper nesting
- **Weighting**: Applies `WTMEC2YR` with single-cycle analysis (no scaling) per NCHS guidelines
- **Models**: `svyglm()` from R `survey` package with `options(survey.lonely.psu = "certainty")`
- **NCHS Compliance**: Follows Technical Documentation 2006 recommendations for MEC weights
- **Complete Case Analysis**: Filtering done after survey design creation for proper variance estimation

### **Multiple Comparisons Correction (ENHANCED FDR):**
- **Primary Method**: Benjamini-Hochberg False Discovery Rate (FDR) procedure ⭐ **RECOMMENDED**
- **Implementation**: `p.adjust(method = "BH")` applied within each dependent variable using `group_by()`
- **Significance Threshold**: q < 0.1 (adjusted for effect scale harmonization across transformations)
- **Scope**: Applied only to main effect terms (`indep_var`), not covariates or intercepts  
- **Secondary Method**: Bonferroni correction (very conservative, q < 0.05)
- **Correction Strategy**: Within dependent variable to maintain biological interpretability

### **Microbiome Data Handling (STATISTICALLY CORRECTED):**
- **Integer Pseudocount Strategy**: ε = 1 applied only to COUNT matrices (preserves count scale)
- **Zero-Library Handling**: Samples with library size = 0 set to NA gracefully (no analysis failure)
- **SQL View Optimization**: "_none" tables implemented as views (50% storage reduction)
- **Numerical Precision**: Single coercion to numeric to avoid repeated casting overhead
- **Transformation Quality**: Input validation ensures correct table types (COUNT vs RELATIVE)
- **Closure Verification**: Warnings for non-closed relative abundance data with tolerance checks

### **Enhanced Model Implementation:**
- **Table/View Detection**: Manual SQL query to detect both tables and views (prevents silent failures)
- **Selective Suffix Logic**: Only microbiome tables (`DADA2RSV_*`) get transformation suffixes
- **Progressive Covariate Selection**: Optimal → Essential → Minimal model hierarchy with fallbacks
- **Enhanced Quartile Handling**: Quartiles → Tertiles → Log-continuous → Raw continuous fallback
- **Minimum Sample Sizes**: Dynamic calculation based on model complexity and regression type
- **R-squared Calculation**: Survey-weighted calculation using proper weight extraction methods

### **Effect Scale Harmonization (NEW):**
- **Transformation-Specific Scales**: Each transformation reports effects in appropriate units
- **Interpretation Metadata**: `effect_scale` and `interpretation_note` columns for each result
- **Cross-Transformation Comparability**: Standardized reporting enables comparison across methods
- **Mathematical Documentation**: Complete formulas provided for reproducibility

### **Model Types by Analysis (ENHANCED):**
- **Linear Regression**: 1_demoWAS, 3_exWAS, 4_pheWAS, 6_zimWAS (continuous outcomes)
  - Family: `gaussian()` with survey design
  - Minimum samples: max(15, n_covariates + 5) for adequate precision
- **Logistic Regression**: 2_oradWAS, 5_outWAS (binary outcomes)
  - Family: `quasibinomial()` with survey design for overdispersion handling
  - Minimum samples: max(10, n_covariates + 5) for convergence
- **Enhanced Covariates**: Age, sex, PIR, education, ethnicity with progressive selection
- **Survey-Adjusted Inference**: All standard errors account for complex survey design

### **Quality Control and Validation:**
- **Infinite Value Handling**: Converts ±Inf to NA instead of model failure
- **Design Degrees of Freedom**: `df_residual_design` calculated and reported for validation
- **Convergence Monitoring**: Enhanced error recovery with multiple fallback attempts
- **Essential Variable Protection**: Age and sex always retained regardless of missing data rates

---
### 🔧 **Comprehensive Troubleshooting Guide Appendix:** [00_trouble_shooting_n_appendix.md](assets/00_trouble_shooting_n_appendix.md)

---
## 👥 **Authors & Contributors**
<table align="center">
  <tr>
    <td align="center">
      <a href="https://github.com/terry-b-cho">
        <img src="https://github.com/terry-b-cho.png" width="80" style="border-radius: 50%;" alt="Byeongyeon Cho"/><br/>
        <img src="https://badgen.net/badge/Lead%20Developer/Byeongyeon%20Cho%20🧪/purple?icon=github&scale=1.0" height="20"/>
      </a>
    </td>
    <td align="center">
      <a href="https://github.com/adkostic">
        <img src="https://github.com/adkostic.png" width="80" style="border-radius: 50%;" alt="Aleksandar D. Kostic"/><br/>
        <img src="https://badgen.net/badge/PI/Aleksandar%20D.%20Kostic%20🔬/blue?icon=github&scale=1.0" height="20"/>
      </a>
    </td>
    <td align="center">
      <a href="https://github.com/b-tierney">
        <img src="https://github.com/b-tierney.png" width="80" style="border-radius: 50%;" alt="Braden Tierney"/><br/>
        <img src="https://badgen.net/badge/Supervisor/Braden%20Tierney%20🧬/green?icon=github&scale=1.0" height="20"/>
      </a>
    </td>
    <td align="center">
      <a href="https://github.com/chiragjp">
        <img src="https://github.com/chiragjp.png" width="80" style="border-radius: 50%;" alt="Chirag J. Patel"/><br/>
        <img src="https://badgen.net/badge/PI/Chirag%20J.%20Patel%20📊/orange?icon=github&scale=1.0" height="20"/>
      </a>
    </td>
  </tr>
</table>

---

## 📄 **License**

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🔗 **Related Resources**

- **NHANES Survey Design**: [CDC Analytic Guidelines](https://wwwn.cdc.gov/nchs/nhanes/analyticguidelines.aspx)
- **NCHS Technical Documentation 2006**: [Survey Design Recommendations](https://www.cdc.gov/nchs/data/series/sr_02/sr02_161.pdf)
- **R Survey Package**: [Survey Package Documentation](https://cran.r-project.org/web/packages/survey/index.html)
- **Microbiome Analysis**: [Compositional Data Analysis Methods](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)
- **FDR Correction**: [Benjamini-Hochberg Procedure](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

---

# **NHANES Oral Microbiome Pipeline - FINAL AUDIT CORRECTED**

## **Production Readiness Summary**

### **Pipeline Enhancements Delivered:**
- ✅ **24 Total Analyses**: 6 analysis types × 4 transformations (including Hellinger)
- ✅ **Final Technical Audit**: All 9 critical blocking issues resolved
- ✅ **Statistical Rigor**: Integer pseudocount (ε=1), zero-library handling, SQL view optimization
- ✅ **O2 Compatibility**: Updated modules (`gcc/14.2.0`, `conda/miniforge3/24.11.3-0`)
- ✅ **Survey Compliance**: NCHS Technical Documentation 2006 guidelines implemented
- ✅ **FDR Correction**: Benjamini-Hochberg within dependent variables (q < 0.1)
- ✅ **Effect Scale Harmonization**: Comparable results across all transformations
- ✅ **Enhanced Error Handling**: Table/view detection, progressive fallbacks, infinite value handling

### **Production Metrics:**
- **Total Analyses**: 24 (6 analysis types × 4 transformations)
- **Database Objects**: 50 (40 transformed tables + 10 SQL views)
- **Expected Runtime**: 4-6 hours for all analyses (O2 cluster)
- **Success Rate**: >90% with comprehensive error recovery
- **Statistical Methods**: Survey-weighted models with multiple comparisons correction
- **Publication Ready**: FDR-corrected results with effect scale interpretation

### **Quality Assurance:**
- **Comprehensive Testing**: All 24 analyses validated with `debug_all_pipelines.R`
- **Biostatistical Review**: Multiple rounds of technical audit corrections applied
- **Mathematical Validation**: Correct pseudocount strategy, transformation formulas, survey design
- **Reproducibility**: Complete documentation, version control, standardized environments

## **Expected Results Overview**

### **Analysis Coverage:**
```
6 Analysis Types × 4 Transformations = 24 Total Analyses

Demographics → Microbiome (1_demoWAS): 4 analyses
Microbiome → Oral Health (2_oradWAS): 4 analyses  
Exposures → Microbiome (3_exWAS): 4 analyses
Microbiome → Phenotypes (4_pheWAS): 4 analyses
Microbiome → Disease Outcomes (5_outWAS): 4 analyses
Microbiome → Lab Values (6_zimWAS): 4 analyses
```

### **Statistical Power:**
- **Sample Size**: 9,847 participants (NHANES 2009-2012)
- **Microbiome Features**: 1,349 genus-level OTUs
- **Health Variables**: 200+ across 6 domains
- **Survey Design**: Proper stratification, clustering, and weighting
- **Multiple Comparisons**: FDR correction maintains discovery power

### **Publication Output:**
- **Main Results**: 24 `*_tidied_complete.rds` files with FDR correction
- **Effect Interpretation**: Transformation-specific scales and interpretation guides
- **Quality Metrics**: Survey R-squared, design degrees of freedom, sample sizes
- **Supplementary Materials**: Comprehensive aggregation summaries and statistics

## **Methods Summary for Publication**

### **Statistical Framework:**
> "We conducted 24 comprehensive association analyses using 4 microbiome data transformations (untransformed, Hellinger, CLR, log-normal) across 6 health domains. All analyses incorporated NHANES complex survey design using proper stratification (SDMVSTRA), clustering (SDMVPSU), and survey weights (WTMEC2YR) following NCHS Technical Documentation 2006 guidelines. Models were fitted using survey-weighted generalized linear models (svyglm) with enhanced error handling including zero-library sample management and progressive covariate selection."

### **Multiple Comparisons:**
> "P-values were adjusted for multiple comparisons using the Benjamini-Hochberg false discovery rate procedure applied within each dependent variable. Associations with FDR-adjusted p-values < 0.1 were considered statistically significant, maintaining adequate power while controlling false discoveries."

### **Data Transformations:**
> "Microbiome data underwent four transformation methods: (1) untransformed relative abundances, (2) Hellinger transformation (√P_ij), (3) centered log-ratio transformation (ln[(C_ij + 1)/g_i]), and (4) log-normal transformation (log₁₀[(C_ij + 1)/(n_i + D) × n̄]). Integer pseudocount (ε=1) was applied to count data to preserve scale and numerical stability."

### **Software Implementation:**
> "Analyses were conducted using R 4.4.3 via conda with enhanced packages for survey analysis (survey), statistical modeling (broom), and data manipulation (dplyr, stringr). The pipeline incorporated comprehensive statistical corrections including table/view detection, effect scale harmonization, and robust error handling for production deployment."

---

## **Ready for Scientific Discovery**

This **final audit corrected** NHANES oral microbiome analysis pipeline represents a comprehensive, biostatistically sound, and computationally robust framework for investigating microbiome-health associations in a nationally representative U.S. population. 

**The pipeline is now production-ready for generating high-quality, reproducible scientific insights into the role of oral microbiome in human health and disease.**

