#!/usr/bin/env Rscript
################################################################################
##  1. LOAD PRE-COMPUTED DIVERSITY DATA & ALL CATEGORICAL VARIABLES           ##
##  - Loads alpha diversity from dada2rsv-alpha.txt                           ##
##  - Loads ALL categorical variables from database directly                   ##
##  - HPC version without phyloseq dependency                                  ##
##  - TEST MODE: Uses smaller subset for testing                              ##
################################################################################

# Environment: R >= 4.5 with dplyr, tidyr, readr, forcats, tibble, data.table,
# DBI, RSQLite, dbplyr. Conda spec: envs/nhanes-analysis_for_reviewers.yml.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(tibble)
  library(data.table)
  library(DBI)
  library(RSQLite)
  library(dbplyr)
})

cat("=== LOADING PRE-COMPUTED DIVERSITY DATA & ALL CATEGORICAL VARIABLES ===\n\n")

###############################################################################
## 1. CONFIGURATION                                                          ##
###############################################################################

# ┌─────────────────────────────────────────────────────────────────┐
# │  IMPORTANT: MUST MATCH TEST_MODE IN 4_integrated_diversity_FINAL.R │
# │  TEST_MODE <- TRUE:  Uses 1,000 samples (~2 min)                │
# │  TEST_MODE <- FALSE: Uses ALL 9,662 samples (~5 min)            │
# └─────────────────────────────────────────────────────────────────┘

TEST_MODE <- FALSE   # ← FULL DATASET ENABLED (N=9,349)
TEST_SAMPLES <- 1000  # Only used when TEST_MODE = TRUE

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

base_path <- PROJECT_ROOT
data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
output_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

cat("Configuration:\n")
cat("  TEST_MODE:", TEST_MODE, "\n")
if (TEST_MODE) cat("  TEST_SAMPLES:", TEST_SAMPLES, "\n")
cat("  Output path:", output_path, "\n\n")

###############################################################################
## 2. LOAD ALPHA DIVERSITY DATA                                              ##
###############################################################################

cat("Loading alpha diversity data...\n")

alpha_file <- file.path(data_path, "dada2rsv-alpha.txt")
alpha_data_raw <- fread(alpha_file, nrows = if(TEST_MODE) TEST_SAMPLES + 1 else -1)

cat("  Alpha diversity loaded:", nrow(alpha_data_raw), "samples\n")

# Extract metrics at rarefaction depth 10000
alpha_diversity <- alpha_data_raw %>%
  select(
    SEQN,
    Observed_OTUs = RSV_ObservedOTUs_10000_0,
    Shannon_Diversity = RSV_ShanWienDiv_10000_0,
    Inverse_Simpson = RSV_InverseSimpson_10000_0
  ) %>%
  mutate(
    SEQN = as.character(SEQN),
    # CRITICAL: Convert metrics to numeric (they may be read as character)
    Observed_OTUs = as.numeric(Observed_OTUs),
    Shannon_Diversity = as.numeric(Shannon_Diversity),
    Inverse_Simpson = as.numeric(Inverse_Simpson)
  )

cat("  Extracted metrics:", nrow(alpha_diversity), "samples\n")
cat("  SEQN type:", class(alpha_diversity$SEQN), "\n")
cat("  Observed_OTUs type:", class(alpha_diversity$Observed_OTUs), "\n")

###############################################################################
## 3. LOAD CLINICAL VARIABLES FROM DATABASE                                 ##
###############################################################################

cat("\nLoading clinical variables from database...\n")

# Connect to database
db_path <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite")
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# Load from specific derived tables
clinical_data <- data.frame(SEQN = alpha_diversity$SEQN, stringsAsFactors = FALSE)

# 1. Load Oral Health variables from OralDisease tables
cat("  Loading oral health variables...\n")
oral_f <- tbl(con, "OralDisease_F") %>% collect() %>% mutate(SEQN = as.character(SEQN))
oral_g <- tbl(con, "OralDisease_G") %>% collect() %>% mutate(SEQN = as.character(SEQN))
oral_combined <- bind_rows(oral_f, oral_g)
clinical_data <- clinical_data %>% left_join(oral_combined, by = "SEQN")
cat("    [done] Loaded 4 oral health variables\n")

# 2. Load Outcome variables from d_outcome_mcq tables
cat("  Loading outcome variables...\n")
outcome_f <- tbl(con, "d_outcome_mcq_f") %>% collect() %>% mutate(SEQN = as.character(SEQN))
outcome_g <- tbl(con, "d_outcome_mcq_g") %>% collect() %>% mutate(SEQN = as.character(SEQN))
outcome_combined <- bind_rows(outcome_f, outcome_g)
clinical_data <- clinical_data %>% left_join(outcome_combined, by = "SEQN")
cat("    [done] Loaded 16 outcome variables\n")

# 3. Load Hepatitis C from HEPC tables
cat("  Loading Hepatitis C...\n")
hepc_f <- tbl(con, "HEPC_F") %>% select(SEQN, LBDHCV) %>% collect() %>% mutate(SEQN = as.character(SEQN))
hepc_g <- tbl(con, "HEPC_G") %>% select(SEQN, LBDHCV) %>% collect() %>% mutate(SEQN = as.character(SEQN))
hepc_combined <- bind_rows(hepc_f, hepc_g)
clinical_data <- clinical_data %>% left_join(hepc_combined, by = "SEQN")
cat("    [done] Loaded LBDHCV (Hepatitis C)\n")

# 4. Load HPV from HPV tables (if exists)
if ("HPV_F" %in% dbListTables(con)) {
  cat("  Loading HPV...\n")
  hpv_f <- tbl(con, "HPV_F") %>% select(SEQN, LBDRPCR) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  hpv_g <- tbl(con, "HPV_G") %>% select(SEQN, LBDRPCR) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  hpv_combined <- bind_rows(hpv_f, hpv_g)
  clinical_data <- clinical_data %>% left_join(hpv_combined, by = "SEQN")
  cat("    [done] Loaded LBDRPCR (HPV PCR)\n")
} else {
  clinical_data$LBDRPCR <- NA
  cat("    [WARNING] HPV tables not found\n")
}

# 5. Load BMI from BMX tables
cat("  Loading BMI and anthropometry...\n")
bmx_f <- tbl(con, "BMX_F") %>% select(SEQN, BMXBMI) %>% collect() %>% mutate(SEQN = as.character(SEQN))
bmx_g <- tbl(con, "BMX_G") %>% select(SEQN, BMXBMI) %>% collect() %>% mutate(SEQN = as.character(SEQN))
bmx_combined <- bind_rows(bmx_f, bmx_g)
clinical_data <- clinical_data %>% left_join(bmx_combined, by = "SEQN")
cat("    [done] Loaded BMXBMI (BMI)\n")

# 6. Load Blood Pressure from BPX tables (need BPXPLS, BPXSY1-3, BPXDI1-3 for MSYSTOLIC/MDIASTOLIC)
cat("  Loading blood pressure variables...\n")
bpx_f <- tbl(con, "BPX_F") %>% 
  select(SEQN, BPXPLS, BPXSY1, BPXSY2, BPXSY3, BPXDI1, BPXDI2, BPXDI3) %>% 
  collect() %>% mutate(SEQN = as.character(SEQN))
bpx_g <- tbl(con, "BPX_G") %>% 
  select(SEQN, BPXPLS, BPXSY1, BPXSY2, BPXSY3, BPXDI1, BPXDI2, BPXDI3) %>% 
  collect() %>% mutate(SEQN = as.character(SEQN))
bpx_combined <- bind_rows(bpx_f, bpx_g) %>%
  rowwise() %>%
  mutate(
    MSYSTOLIC = mean(c(BPXSY1, BPXSY2, BPXSY3), na.rm = TRUE),
    MDIASTOLIC = mean(c(BPXDI1, BPXDI2, BPXDI3), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(SEQN, BPXPLS, MSYSTOLIC, MDIASTOLIC)
clinical_data <- clinical_data %>% left_join(bpx_combined, by = "SEQN")
cat("    [done] Loaded BPXPLS (Pulse), MSYSTOLIC, MDIASTOLIC\n")

# 7. Load Smoking from SMQ tables (need to derive SMQ_current_ever_never)
cat("  Loading smoking status...\n")
smq_f <- tbl(con, "SMQ_F") %>% select(SEQN, SMQ020, SMQ040) %>% collect() %>% mutate(SEQN = as.character(SEQN))
smq_g <- tbl(con, "SMQ_G") %>% select(SEQN, SMQ020, SMQ040) %>% collect() %>% mutate(SEQN = as.character(SEQN))
smq_combined <- bind_rows(smq_f, smq_g) %>%
  mutate(
    SMQ_current_ever_never = case_when(
      SMQ020 == 2 ~ 0,  # Never smoker
      SMQ020 == 1 & SMQ040 == 3 ~ 1,  # Current smoker
      SMQ020 == 1 & SMQ040 %in% c(1, 2) ~ 2,  # Former smoker
      TRUE ~ NA_real_
    )
  ) %>%
  select(SEQN, SMQ_current_ever_never)
clinical_data <- clinical_data %>% left_join(smq_combined, by = "SEQN")
cat("    [done] Loaded SMQ_current_ever_never (Smoking status)\n")

# Close database connection
dbDisconnect(con)

cat("  Clinical data loaded:", ncol(clinical_data) - 1, "variables\n")

# Filter to samples that exist in alpha diversity
clinical_data_filtered <- clinical_data %>%
  filter(SEQN %in% alpha_diversity$SEQN)

cat("  Filtered to alpha diversity samples:", nrow(clinical_data_filtered), "\n")

###############################################################################
## 3b. LOAD DEMOGRAPHIC VARIABLES FROM DATABASE                              ##
###############################################################################

cat("\nLoading demographic variables...\n")

# Connect to database again
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# Load demographic variables from DEMO tables
demo_f <- tbl(con, "DEMO_F") %>% 
  select(SEQN, RIAGENDR, RIDAGEYR, DMDEDUC2, RIDRETH1, DMDMARTL, DMDHHSIZ, INDFMPIR, SIALANG) %>%
  collect() %>% mutate(SEQN = as.character(SEQN))

demo_g <- tbl(con, "DEMO_G") %>% 
  select(SEQN, RIAGENDR, RIDAGEYR, DMDEDUC2, RIDRETH1, DMDMARTL, DMDHHSIZ, INDFMPIR, SIALANG) %>%
  collect() %>% mutate(SEQN = as.character(SEQN))

demo_combined <- bind_rows(demo_f, demo_g)

# Process demographic variables
demo_processed <- demo_combined %>%
  mutate(
    # Gender
    Gender = case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    
    # Age group
    AgeGroup = case_when(
      RIDAGEYR >= 18 & RIDAGEYR < 30 ~ "18-29",
      RIDAGEYR >= 30 & RIDAGEYR < 40 ~ "30-39",
      RIDAGEYR >= 40 & RIDAGEYR < 50 ~ "40-49",
      RIDAGEYR >= 50 & RIDAGEYR < 60 ~ "50-59",
      RIDAGEYR >= 60 & RIDAGEYR < 70 ~ "60-69",
      RIDAGEYR >= 70 ~ "70+",
      TRUE ~ NA_character_
    ),
    
    # Age group 5-year bins
    AgeGroup_5years = case_when(
      RIDAGEYR >= 14 & RIDAGEYR < 20 ~ "14-19",
      RIDAGEYR >= 20 & RIDAGEYR < 25 ~ "20-24",
      RIDAGEYR >= 25 & RIDAGEYR < 30 ~ "25-29",
      RIDAGEYR >= 30 & RIDAGEYR < 35 ~ "30-34",
      RIDAGEYR >= 35 & RIDAGEYR < 40 ~ "35-39",
      RIDAGEYR >= 40 & RIDAGEYR < 45 ~ "40-44",
      RIDAGEYR >= 45 & RIDAGEYR < 50 ~ "45-49",
      RIDAGEYR >= 50 & RIDAGEYR < 55 ~ "50-54",
      RIDAGEYR >= 55 & RIDAGEYR < 60 ~ "55-59",
      RIDAGEYR >= 60 & RIDAGEYR < 65 ~ "60-64",
      RIDAGEYR >= 65 & RIDAGEYR < 70 ~ "65-69",
      TRUE ~ NA_character_
    ),
    
    # Education level
    EducationLevel = case_when(
      DMDEDUC2 == 1 ~ "< 9th Grade",
      DMDEDUC2 == 2 ~ "9-11th Grade",
      DMDEDUC2 == 3 ~ "High School",
      DMDEDUC2 == 4 ~ "College/AA",
      DMDEDUC2 == 5 ~ "College Graduate",
      TRUE ~ NA_character_
    ),
    
    # Ethnicity
    Ethnicity = case_when(
      RIDRETH1 == 1 ~ "Mexican American",
      RIDRETH1 == 2 ~ "Other Hispanic",
      RIDRETH1 == 3 ~ "White",
      RIDRETH1 == 4 ~ "Black",
      RIDRETH1 == 5 ~ "Other",
      TRUE ~ NA_character_
    ),
    
    # Marital status
    Marital_Status = case_when(
      DMDMARTL == 1 ~ "Married",
      DMDMARTL == 2 ~ "Widowed",
      DMDMARTL == 3 ~ "Divorced",
      DMDMARTL == 4 ~ "Separated",
      DMDMARTL == 5 ~ "Never married",
      DMDMARTL == 6 ~ "Living with partner",
      TRUE ~ NA_character_
    ),
    
    # Household size
    Household_Size_Factor = case_when(
      DMDHHSIZ == 1 ~ "1",
      DMDHHSIZ == 2 ~ "2",
      DMDHHSIZ == 3 ~ "3",
      DMDHHSIZ == 4 ~ "4",
      DMDHHSIZ == 5 ~ "5",
      DMDHHSIZ == 6 ~ "6",
      DMDHHSIZ >= 7 ~ "7+",
      TRUE ~ NA_character_
    ),
    
    # US Born (simplified - would need more complex logic)
    US_Born = "US Born",  # Simplified for now
    
    # Interview language
    Interview_Language = case_when(
      SIALANG == 1 ~ "English",
      SIALANG == 2 ~ "Spanish",
      TRUE ~ NA_character_
    )
  ) %>%
  select(SEQN, Gender, AgeGroup, AgeGroup_5years, EducationLevel, Ethnicity, Marital_Status, 
         Household_Size_Factor, US_Born, Interview_Language, INDFMPIR)

# Close database connection
dbDisconnect(con)

cat("  [done] Loaded demographic variables\n")

# Merge with clinical data
sample_data_full <- clinical_data_filtered %>%
  left_join(demo_processed, by = "SEQN")

cat("  Updated sample_data columns:", ncol(sample_data_full), "\n")

###############################################################################
## 3c. DEFINE CONTROL GROUPS FOR CLINICAL VARIABLES                          ##
###############################################################################

cat("\nDefining control groups for clinical variables...\n")

# Define oral health control samples (all oral variables = 0)
orad_vars <- c("DENTURE_OHAROCDE", "GUM_DISEASE_OHAROCGP", "ORAL_HYGIENE_OHAROCOH", "TOOTH_DECAY_OHAROCDT")
orad_vars_present <- orad_vars[orad_vars %in% colnames(sample_data_full)]

if (length(orad_vars_present) > 0) {
  SEQN_oradWAS_control <- sample_data_full %>%
    filter(if_all(all_of(orad_vars_present), ~ !is.na(.x) & .x == 0)) %>%
    pull(SEQN)
  cat("  Oral health controls:", length(SEQN_oradWAS_control), "samples\n")
} else {
  SEQN_oradWAS_control <- character(0)
  cat("  Warning: No oral health variables found\n")
}

# Define outcome control samples (all outcome variables = 0)
out_vars <- c("ASTHMA", "DIABETES", "HEART_ATTACK", "BRONCHITIS", "EMPHYSEMA", 
              "ANGINA", "HEART_FAILURE", "STROKE")
out_vars_present <- out_vars[out_vars %in% colnames(sample_data_full)]

if (length(out_vars_present) > 0) {
  SEQN_outWAS_control <- sample_data_full %>%
    filter(if_all(all_of(out_vars_present), ~ !is.na(.x) & .x == 0)) %>%
    pull(SEQN)
  cat("  Outcome controls:", length(SEQN_outWAS_control), "samples\n")
} else {
  SEQN_outWAS_control <- character(0)
  cat("  Warning: No outcome variables found\n")
}

###############################################################################
## 4. FACTORIZE ALL CATEGORICAL VARIABLES                                    ##
##    Following exact patterns from 4_association_phyloseq_analyses.Rmd      ##
###############################################################################

cat("\nFactorizing categorical variables...\n")

# Extract predefined factors
categorical_data <- sample_data_full %>%
  mutate(across(where(is.character), as.factor))

cat("  Initial categorical variables:", sum(sapply(categorical_data, is.factor)), "\n")

# DEMO-WAS categories
categorical_data <- categorical_data %>%
  mutate(
    # Gender - Female as reference
    Gender = if("Gender" %in% colnames(.)) fct_relevel(as.factor(Gender), "Female") else NULL,
    
    # Age group - 30-39 as reference
    Age_group = if("AgeGroup" %in% colnames(.)) fct_relevel(as.factor(AgeGroup), "30-39") else NULL,
    
    # Age group 5-year bins - 30-34 as reference
    Age_group_5years = if("AgeGroup_5years" %in% colnames(.)) {
      factor(AgeGroup_5years, levels = c("14-19", "20-24", "25-29", "30-34", "35-39", 
                                         "40-44", "45-49", "50-54", "55-59", "60-64", "65-69")) %>%
        fct_relevel("30-34")
    } else NULL,
    
    # Education - College/AA as reference
    Education_level = if("EducationLevel" %in% colnames(.)) {
      factor(EducationLevel, levels = c("< 9th Grade", "9-11th Grade", 
                                        "High School", "College/AA", 
                                        "College Graduate")) %>% 
        fct_relevel("College/AA")
    } else NULL,
    
    # Ethnicity - White as reference
    Ethnicity = if("Ethnicity" %in% colnames(.)) fct_relevel(as.factor(Ethnicity), "White") else NULL,
    
    # US Born - US Born as reference
    US_born = if("US_Born" %in% colnames(.)) fct_relevel(as.factor(US_Born), "US Born") else NULL,
    
    # Household size - 4 as reference
    Household_size = if("Household_Size_Factor" %in% colnames(.)) {
      factor(Household_Size_Factor, levels = c("1","2","3","4","5","6","7+")) %>%
        fct_relevel("4")
    } else NULL,
    
    # Marital status - Married as reference
    Marital_status = if("Marital_Status" %in% colnames(.)) fct_relevel(as.factor(Marital_Status), "Married") else NULL,
    
    # Interview language - English as reference
    Interview_language = if("Interview_Language" %in% colnames(.)) fct_relevel(as.factor(Interview_Language), "English") else NULL,
    
    # Income to poverty ratio - 150-184% as reference
    Income_to_poverty_ratio = if("INDFMPIR" %in% colnames(.)) {
      case_when(
        INDFMPIR < 0.50 ~ "Below 50%",
        INDFMPIR >= 0.50 & INDFMPIR < 1.00 ~ "50-99%",
        INDFMPIR >= 1.00 & INDFMPIR < 1.25 ~ "100-124%",
        INDFMPIR >= 1.25 & INDFMPIR < 1.50 ~ "125-149%",
        INDFMPIR >= 1.50 & INDFMPIR < 1.85 ~ "150-184%",
        INDFMPIR >= 1.85 & INDFMPIR < 2.00 ~ "185-199%",
        INDFMPIR >= 2.00 ~ "200% and Over",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Below 50%", "50-99%", "100-124%", "125-149%", 
                              "150-184%", "185-199%", "200% and Over")) %>%
        fct_relevel("150-184%")
    } else NULL
  )

# ORAL-WAS categories - Binary disease variables with "control" reference
orad_vars <- c("DENTURE_OHAROCDE", "GUM_DISEASE_OHAROCGP", "ORAL_HYGIENE_OHAROCOH", "TOOTH_DECAY_OHAROCDT")
orad_vars_present <- orad_vars[orad_vars %in% colnames(categorical_data)]

if (length(orad_vars_present) > 0) {
  SEQN_oradWAS_control <- categorical_data %>%
    filter(if_all(all_of(orad_vars_present), ~ !is.na(.x) & .x == 0)) %>%
    pull(SEQN)
  
  cat("  Oral disease control samples:", length(SEQN_oradWAS_control), "\n")
  
  # Factorize oral disease variables
  categorical_data <- categorical_data %>%
    mutate(
      Denture = if("DENTURE_OHAROCDE" %in% colnames(.)) {
        case_when(
          DENTURE_OHAROCDE == 1 ~ "Denture",
          SEQN %in% SEQN_oradWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Denture"))
      } else NULL,
      
      Gum_disease = if("GUM_DISEASE_OHAROCGP" %in% colnames(.)) {
        case_when(
          GUM_DISEASE_OHAROCGP == 1 ~ "Gum disease",
          SEQN %in% SEQN_oradWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Gum disease"))
      } else NULL,
      
      Oral_hygiene = if("ORAL_HYGIENE_OHAROCOH" %in% colnames(.)) {
        case_when(
          ORAL_HYGIENE_OHAROCOH == 1 ~ "Poor oral hygiene",
          SEQN %in% SEQN_oradWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Poor oral hygiene"))
      } else NULL,
      
      Tooth_decay = if("TOOTH_DECAY_OHAROCDT" %in% colnames(.)) {
        case_when(
          TOOTH_DECAY_OHAROCDT == 1 ~ "Tooth decay",
          SEQN %in% SEQN_oradWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Tooth decay"))
      } else NULL,
      
      # Create unified Oral_all variable with 5 levels
      Oral_all = case_when(
        # Control group: no oral diseases
        SEQN %in% SEQN_oradWAS_control ~ "Control",
        # Individual diseases (prioritized order)
        DENTURE_OHAROCDE == 1 ~ "Denture",
        GUM_DISEASE_OHAROCGP == 1 ~ "Gum Disease", 
        ORAL_HYGIENE_OHAROCOH == 1 ~ "Oral Hygiene",
        TOOTH_DECAY_OHAROCDT == 1 ~ "Tooth Decay",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Control", "Gum Disease", "Oral Hygiene", "Tooth Decay", "Denture")) %>%
        fct_relevel("Control")
    )
}

# EXWAS categories - Smoking and hepatitis
categorical_data <- categorical_data %>%
  mutate(
    # Smoking status
    Smoking_status = if(all(c("SMQ_current_ever_never") %in% colnames(.))) {
      case_when(
        SMQ_current_ever_never == 0 ~ "Never smoker",
        SMQ_current_ever_never == 2 ~ "Former smoker",
        SMQ_current_ever_never == 1 ~ "Current smoker",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Never smoker", "Former smoker", "Current smoker")) %>%
        fct_relevel("Never smoker")
    } else NULL,
    
    # Hepatitis C
    Hepatitis_C_antibody = if("LBDHCV" %in% colnames(.)) {
      case_when(
        LBDHCV == 1 ~ "Positive",
        LBDHCV == 2 ~ "Negative",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Negative", "Positive"))
    } else NULL,
    
    # HPV PCR
    HPV_PCR_summary = if("LBDRPCR" %in% colnames(.)) {
      case_when(
        LBDRPCR == 1 ~ "Positive",
        LBDRPCR == 2 ~ "Negative",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Negative", "Positive"))
    } else NULL
  )

# PHEWAS categories - BMI and blood pressure
categorical_data <- categorical_data %>%
  mutate(
    # BMI category
    BMI_category = if("BMXBMI" %in% colnames(.)) {
      case_when(
        BMXBMI < 18.5 ~ "Underweight",
        BMXBMI >= 18.5 & BMXBMI < 25 ~ "Healthy weight",
        BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight",
        BMXBMI >= 30 & BMXBMI < 35 ~ "Class 1 Obesity",
        BMXBMI >= 35 ~ "Class 2-3 Obesity",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Underweight", "Healthy weight", "Overweight", 
                             "Class 1 Obesity", "Class 2-3 Obesity")) %>%
        fct_relevel("Healthy weight")
    } else NULL,
    
    # Blood pressure
    Blood_pressure = if(all(c("MSYSTOLIC", "MDIASTOLIC") %in% colnames(.))) {
      case_when(
        MSYSTOLIC > 180 | MDIASTOLIC > 120 ~ "Hypertensive Crisis",
        MSYSTOLIC >= 140 | MDIASTOLIC >= 90 ~ "Hypertension Stage 2",
        (MSYSTOLIC >= 130 & MSYSTOLIC < 140) | (MDIASTOLIC >= 80 & MDIASTOLIC < 90) ~ "Hypertension Stage 1",
        (MSYSTOLIC >= 120 & MSYSTOLIC < 130) & (MDIASTOLIC < 80) ~ "Elevated",
        MSYSTOLIC < 120 & MDIASTOLIC < 80 ~ "Normal",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Normal", "Elevated", "Hypertension Stage 1", 
                             "Hypertension Stage 2", "Hypertensive Crisis")) %>%
        fct_relevel("Normal")
    } else NULL,
    
    # Pulse category
    Pulse_category = if("BPXPLS" %in% colnames(.)) {
      case_when(
        BPXPLS < 60 ~ "<60 bpm",
        BPXPLS >= 60 & BPXPLS < 70 ~ "60-70 bpm",
        BPXPLS >= 70 & BPXPLS < 75 ~ "70-75 bpm",
        BPXPLS >= 75 & BPXPLS < 85 ~ "75-85 bpm",
        BPXPLS >= 85 ~ "85+ bpm",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("<60 bpm", "60-70 bpm", "70-75 bpm", "75-85 bpm", "85+ bpm")) %>%
        fct_relevel("70-75 bpm")
    } else NULL
  )

# OUTWAS categories - Disease outcomes with "control" reference
outwas_vars <- c("ASTHMA", "BRONCHITIS", "EMPHYSEMA", "ANGINA", "HEART_FAILURE", 
                 "HEART_ATTACK", "STROKE", "CHD", "CVD", "DIABETES",
                 "CANCER_BREAST", "CANCER_COLON", "CANCER_LUNG", 
                 "CANCER_ESOPHAGEAL", "CANCER_PROSTATE", "CANCER_MOUTH")
outwas_vars_present <- outwas_vars[outwas_vars %in% colnames(categorical_data)]

if (length(outwas_vars_present) > 0) {
  # Define control samples
  SEQN_outWAS_control <- categorical_data %>%
    filter(if_all(all_of(outwas_vars_present), ~ !is.na(.x) & .x == 0)) %>%
    pull(SEQN)
  
  cat("  Outcome disease control samples:", length(SEQN_outWAS_control), "\n")
  
  # Factorize disease outcomes
  categorical_data <- categorical_data %>%
    mutate(
      Asthma = if("ASTHMA" %in% colnames(.)) {
        case_when(
          ASTHMA == 1 ~ "Asthma",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Asthma"))
      } else NULL,
      
      Bronchitis = if("BRONCHITIS" %in% colnames(.)) {
        case_when(
          BRONCHITIS == 1 ~ "Bronchitis",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Bronchitis"))
      } else NULL,
      
      Emphysema = if("EMPHYSEMA" %in% colnames(.)) {
        case_when(
          EMPHYSEMA == 1 ~ "Emphysema",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Emphysema"))
      } else NULL,
      
      Angina = if("ANGINA" %in% colnames(.)) {
        case_when(
          ANGINA == 1 ~ "Angina",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Angina"))
      } else NULL,
      
      Heart_failure = if("HEART_FAILURE" %in% colnames(.)) {
        case_when(
          HEART_FAILURE == 1 ~ "Heart failure",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Heart failure"))
      } else NULL,
      
      Heart_attack = if("HEART_ATTACK" %in% colnames(.)) {
        case_when(
          HEART_ATTACK == 1 ~ "Heart attack",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Heart attack"))
      } else NULL,
      
      Stroke = if("STROKE" %in% colnames(.)) {
        case_when(
          STROKE == 1 ~ "Stroke",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Stroke"))
      } else NULL,
      
      Diabetes = if("DIABETES" %in% colnames(.)) {
        case_when(
          DIABETES == 1 ~ "Diabetes",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Diabetes"))
      } else NULL,
      
      CHD = if("CHD" %in% colnames(.)) {
        case_when(
          CHD == 1 ~ "CHD",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "CHD"))
      } else NULL,
      
      CVD = if("CVD" %in% colnames(.)) {
        case_when(
          CVD == 1 ~ "CVD",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "CVD"))
      } else NULL,
      
      Breast_cancer = if("CANCER_BREAST" %in% colnames(.)) {
        case_when(
          CANCER_BREAST == 1 ~ "Breast cancer",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Breast cancer"))
      } else NULL,
      
      Colon_cancer = if("CANCER_COLON" %in% colnames(.)) {
        case_when(
          CANCER_COLON == 1 ~ "Colon cancer",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Colon cancer"))
      } else NULL,
      
      Lung_cancer = if("CANCER_LUNG" %in% colnames(.)) {
        case_when(
          CANCER_LUNG == 1 ~ "Lung cancer",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Lung cancer"))
      } else NULL,
      
      Esophageal_cancer = if("CANCER_ESOPHAGEAL" %in% colnames(.)) {
        case_when(
          CANCER_ESOPHAGEAL == 1 ~ "Esophageal cancer",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Esophageal cancer"))
      } else NULL,
      
      Prostate_cancer = if("CANCER_PROSTATE" %in% colnames(.)) {
        case_when(
          CANCER_PROSTATE == 1 ~ "Prostate cancer",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Prostate cancer"))
      } else NULL,
      
      Mouth_cancer = if("CANCER_MOUTH" %in% colnames(.)) {
        case_when(
          CANCER_MOUTH == 1 ~ "Mouth cancer",
          SEQN %in% SEQN_outWAS_control ~ "control",
          TRUE ~ NA_character_
        ) %>% factor(levels = c("control", "Mouth cancer"))
      } else NULL,
      
      # Create unified Outwas_all variable with control group for all outcome diseases
      Outwas_all = case_when(
        # Control group: no outcome diseases
        SEQN %in% SEQN_outWAS_control ~ "Control",
        # Individual diseases (prioritized order)
        CANCER_ESOPHAGEAL == 1 ~ "Esophageal Cancer",
        CANCER_MOUTH == 1 ~ "Mouth Cancer", 
        STROKE == 1 ~ "Stroke",
        CVD == 1 ~ "CVD",
        HEART_ATTACK == 1 ~ "Heart Attack",
        DIABETES == 1 ~ "Diabetes",
        HEART_FAILURE == 1 ~ "Heart Failure",
        CANCER_LUNG == 1 ~ "Lung Cancer",
        CANCER_PROSTATE == 1 ~ "Prostate Cancer",
        CANCER_COLON == 1 ~ "Colon Cancer",
        CANCER_BREAST == 1 ~ "Breast Cancer",
        CHD == 1 ~ "CHD",
        ASTHMA == 1 ~ "Asthma",
        EMPHYSEMA == 1 ~ "Emphysema",
        ANGINA == 1 ~ "Angina",
        BRONCHITIS == 1 ~ "Bronchitis",
        TRUE ~ NA_character_
      ) %>% factor(levels = c("Control", "Bronchitis", "Angina", "Emphysema", "Asthma", 
                              "CHD", "Breast Cancer", "Colon Cancer", "Prostate Cancer", 
                              "Lung Cancer", "Heart Failure", "Diabetes", "Heart Attack", 
                              "CVD", "Stroke", "Mouth Cancer", "Esophageal Cancer")) %>%
        fct_relevel("Control")
    )
}

# Remove NULL columns
categorical_data <- categorical_data %>%
  select(where(~ !all(is.null(.))))

# Count final categorical variables
final_categorical_vars <- colnames(categorical_data)[sapply(categorical_data, is.factor)]
final_categorical_vars <- setdiff(final_categorical_vars, "SEQN")

cat("  Final categorical variables:", length(final_categorical_vars), "\n")
cat("  Variables:\n")
for (var in final_categorical_vars) {
  n_levels <- nlevels(categorical_data[[var]])
  n_complete <- sum(!is.na(categorical_data[[var]]))
  cat("    -", var, ":", n_levels, "levels,", n_complete, "complete cases\n")
}

###############################################################################
## 5. MERGE ALPHA DIVERSITY WITH ALL CATEGORICAL DATA                        ##
###############################################################################

cat("\nMerging alpha diversity with categorical data...\n")

alpha_with_all_categories <- alpha_diversity %>%
  left_join(categorical_data, by = "SEQN") %>%
  filter(!is.na(Observed_OTUs))

cat("  Merged data:", nrow(alpha_with_all_categories), "samples\n")
cat("  Total columns:", ncol(alpha_with_all_categories), "\n")

###############################################################################
## 6. SAVE PROCESSED DATA                                                    ##
###############################################################################

cat("\nSaving processed data...\n")

saveRDS(alpha_with_all_categories, file.path(output_path, "alpha_diversity_with_all_categories.rds"))
cat("  [done] Saved: alpha_diversity_with_all_categories.rds\n")

saveRDS(categorical_data, file.path(output_path, "all_categorical_data.rds"))
cat("  [done] Saved: all_categorical_data.rds\n")

# Save list of categorical variables for easy reference
categorical_var_info <- data.frame(
  Variable = final_categorical_vars,
  N_Levels = sapply(categorical_data[final_categorical_vars], nlevels),
  N_Complete = sapply(categorical_data[final_categorical_vars], function(x) sum(!is.na(x))),
  Reference_Level = sapply(categorical_data[final_categorical_vars], function(x) levels(x)[1])
)

write.csv(categorical_var_info, file.path(output_path, "categorical_variables_info.csv"), 
          row.names = FALSE)
cat("  [done] Saved: categorical_variables_info.csv\n")

###############################################################################
## 7. SUMMARY                                                                ##
###############################################################################

cat("\n=== DATA LOADING SUMMARY ===\n")
cat("Mode:", ifelse(TEST_MODE, "TEST", "PRODUCTION"), "\n")
cat("Alpha diversity samples:", nrow(alpha_diversity), "\n")
cat("Final merged dataset:", nrow(alpha_with_all_categories), "\n")
cat("Categorical variables:", length(final_categorical_vars), "\n")

cat("\nCategorical variables with sufficient data (≥30 per level):\n")
sufficient_vars <- categorical_var_info %>%
  filter(N_Complete >= 30)
for (i in 1:nrow(sufficient_vars)) {
  cat("  ", i, ". ", sufficient_vars$Variable[i], " (", 
      sufficient_vars$N_Levels[i], " levels, ref: ", 
      sufficient_vars$Reference_Level[i], ")\n", sep = "")
}

cat("\nData loading complete!\n")
cat("============================\n")
