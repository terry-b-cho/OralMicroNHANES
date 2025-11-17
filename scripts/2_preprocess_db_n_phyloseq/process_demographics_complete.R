#!/usr/bin/env Rscript

# ══════════════════════════════════════════════════════════════════════════════
# COMPREHENSIVE DEMOGRAPHIC PROCESSING PIPELINE
# ══════════════════════════════════════════════════════════════════════════════
# 
# This script consolidates demographic variable derivation, quartile calculation,
# and factor application into a single comprehensive workflow.
# Saves results to a new processed SQLite database.
#
# Author: Assistant
# Date: 2025-06-18
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(fs)
  library(forcats)
})

# Create logs directory
logs_dir <- "scripts/2_preprocess_db_n_phyloseq/logs"
dir_create(logs_dir, recurse = TRUE)

# ── CLI options ─────────────────────────────────────────────
opt_list <- list(
  make_option(c("--input_db"), type = "character", 
              default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite",
              help = "Path to input SQLite database [default: %default]"),
  make_option(c("--output_db"), type = "character",
              default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite",
              help = "Path to output processed SQLite database [default: %default]"),
  make_option(c("--test"), type = "logical", action = "store_true", default = FALSE,
              help = "Test mode - don't create output database [default: %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Set up logging
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(logs_dir, paste0("process_demographics_complete_", timestamp, ".log"))
log_con <- file(log_file, open = "wt")

# Function to log both to console and file
log_message <- function(...) {
  msg <- paste(..., sep = "")
  message(msg)
  writeLines(paste(Sys.time(), ":", msg), log_con)
  flush(log_con)
}

log_message("Starting comprehensive demographic processing")
log_message("Input database: ", opt$input_db)
log_message("Output database: ", opt$output_db)
log_message("Test mode: ", ifelse(opt$test, "YES", "NO"))
log_message("Log file: ", log_file)

# Connect to input database
con_input <- dbConnect(SQLite(), opt$input_db)
log_message("✅ Input database connection established")

# ── DERIVE DEMOGRAPHIC VARIABLES FUNCTION ──────────────────

derive_demographic_variables <- function(demo_data, cycle) {
  log_message("Processing demographic variables for cycle ", cycle, "...")
  
  demo_derived <- demo_data %>%
    mutate(
      # ═══════════════════════════════════════════════════════
      # CORE IDENTIFIERS & BASIC VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Ensure SEQN is character and create sample alias
      SEQN = as.character(SEQN),
      sample = as.character(SEQN),  # Alternative naming for phyloseq
      cycle = cycle,
      
      # Release cycle with full labels (for compatibility)
      Release_Cycle = case_when(
        cycle == "F" ~ "2009-2010",
        cycle == "G" ~ "2011-2012", 
        TRUE ~ NA_character_
      ),
      
      # Binary gender (0/1)
      RIAGENDR_01 = case_when(
        RIAGENDR == 1 ~ 0,  # Male = 0
        RIAGENDR == 2 ~ 1,  # Female = 1
        TRUE ~ NA_real_
      ),
      
      # ═══════════════════════════════════════════════════════
      # CATEGORICAL DEMOGRAPHIC VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Gender (categorical)
      Gender = case_when(
        RIAGENDR == 1 ~ "Male",
        RIAGENDR == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      
      # Age groups
      AgeGroup = case_when(
        RIDAGEYR >= 60 & RIDAGEYR <= 69 ~ "60-69",
        RIDAGEYR >= 50 & RIDAGEYR < 60 ~ "50-59", 
        RIDAGEYR >= 40 & RIDAGEYR < 50 ~ "40-49",
        RIDAGEYR >= 30 & RIDAGEYR < 40 ~ "30-39",
        RIDAGEYR >= 20 & RIDAGEYR < 30 ~ "20-29",
        RIDAGEYR >= 14 & RIDAGEYR < 20 ~ "14-19",
        TRUE ~ NA_character_
      ),
      
      # Alternative age group naming (with different formatting)
      age_group = case_when(
        RIDAGEYR >= 60 & RIDAGEYR <= 69 ~ "60–69",
        RIDAGEYR >= 50 & RIDAGEYR < 60 ~ "50–59", 
        RIDAGEYR >= 40 & RIDAGEYR < 50 ~ "40–49",
        RIDAGEYR >= 30 & RIDAGEYR < 40 ~ "30–39",
        RIDAGEYR >= 20 & RIDAGEYR < 30 ~ "20–29",
        RIDAGEYR >= 14 & RIDAGEYR < 20 ~ "14–19",
        TRUE ~ NA_character_
      ),
      
      # Age (capped at 85)
      Age = case_when(
        RIDAGEYR > 85 ~ 85,
        RIDAGEYR >= 0 ~ RIDAGEYR,
        TRUE ~ NA_real_
      ),
      
      # ═══════════════════════════════════════════════════════
      # EDUCATION VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Education level (from existing dummy variables)
      EducationLevel = case_when(
        EDUCATION_LESS9 == 1 ~ "< 9th Grade",
        EDUCATION_9_11 == 1 ~ "9-11th Grade", 
        EDUCATION_HSGRAD == 1 ~ "High School",
        EDUCATION_AA == 1 ~ "College/AA",
        EDUCATION_COLLEGEGRAD == 1 ~ "College Graduate",
        TRUE ~ NA_character_
      ),
      
      # Alternative education naming (multiple versions)
      Education_Level = EducationLevel,
      
      # Education level with simplified naming for phyloseq
      EducationLevel_Simple = case_when(
        EDUCATION_HSGRAD == 1 ~ "High School",
        EDUCATION_AA == 1 ~ "College/AA", 
        EDUCATION_COLLEGEGRAD == 1 ~ "College Graduate",
        EDUCATION_LESS9 == 1 ~ "< 9th Grade",
        EDUCATION_9_11 == 1 ~ "9-11th Grade",
        TRUE ~ NA_character_
      ),
      
      # ═══════════════════════════════════════════════════════
      # ETHNICITY & CITIZENSHIP VARIABLES  
      # ═══════════════════════════════════════════════════════
      
      # Ethnicity (from existing dummy variables)
      Ethnicity = case_when(
        ETHNICITY_NONHISPANICWHITE == 1 ~ "White",
        ETHNICITY_NONHISPANICBLACK == 1 ~ "Black", 
        ETHNICITY_MEXICAN == 1 ~ "Mexican",
        ETHNICITY_OTHERHISPANIC == 1 ~ "Other Hispanic",
        ETHNICITY_OTHER == 1 ~ "Other",
        TRUE ~ NA_character_
      ),
      
      # US Born status (from existing BORN_INUSA)
      BornInUSA = case_when(
        BORN_INUSA == 1 ~ "US Born",
        BORN_INUSA == 0 ~ "non-US Born", 
        TRUE ~ NA_character_
      ),
      
      # Alternative naming
      US_Born = BornInUSA,
      US_Citizen = BornInUSA,
      
      # ═══════════════════════════════════════════════════════
      # INCOME & POVERTY VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Poverty Income Ratio categories
      PIRcat = case_when(
        INDFMPIR < 1 ~ "<1",
        INDFMPIR >= 1 & INDFMPIR < 2 ~ "1-2",
        INDFMPIR >= 2 & INDFMPIR < 4 ~ "2-4", 
        INDFMPIR >= 4 ~ ">4",
        TRUE ~ NA_character_
      ),
      
      # Clean PIR variable
      Ratio_Family_Income_Poverty = as.numeric(INDFMPIR),
      
      # ═══════════════════════════════════════════════════════
      # HOUSEHOLD & FAMILY VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Household size (capped at 7+)
      Household_Size = case_when(
        DMDHHSIZ %in% 1:6 ~ as.numeric(DMDHHSIZ),
        DMDHHSIZ >= 7 ~ 7,
        DMDHHSIZ %in% c(77, 99) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      
      # Family size (capped at 7+)
      Family_Size = case_when(
        DMDFMSIZ %in% 1:6 ~ as.numeric(DMDFMSIZ), 
        DMDFMSIZ >= 7 ~ 7,
        DMDFMSIZ %in% c(77, 99) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      
      # Household size as factor
      Household_Size_Factor = case_when(
        Household_Size <= 6 ~ as.character(Household_Size),
        Household_Size >= 7 ~ "7+",
        TRUE ~ NA_character_
      ),
      
      # ═══════════════════════════════════════════════════════
      # MARITAL STATUS
      # ═══════════════════════════════════════════════════════
      
      Marital_Status = case_when(
        DMDMARTL == 1 ~ "Married",
        DMDMARTL == 2 ~ "Widowed",
        DMDMARTL == 3 ~ "Divorced", 
        DMDMARTL == 4 ~ "Separated",
        DMDMARTL == 5 ~ "Never Married",
        DMDMARTL == 6 ~ "Living with Partner",
        DMDMARTL %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      
      # ═══════════════════════════════════════════════════════
      # INTERVIEW LANGUAGE & PROXY VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Interview language
      Interview_Language = case_when(
        SIALANG == 1 ~ "English",
        SIALANG == 2 ~ "Spanish", 
        SIALANG %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      
      # Proxy used in interview
      Interview_Proxy = case_when(
        SIAPROXY == 1 ~ "Yes",
        SIAPROXY == 2 ~ "No",
        SIAPROXY %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      
      # Interpreter used 
      Interview_Interpreter = case_when(
        SIAINTRP == 1 ~ "Yes",
        SIAINTRP == 2 ~ "No",
        SIAINTRP %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      
      # ═══════════════════════════════════════════════════════
      # HOUSEHOLD REFERENCE PERSON VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # HH Reference Person's Gender
      HH_Reference_Gender = case_when(
        DMDHRGND == 1 ~ "Male",
        DMDHRGND == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      
      # HH Reference Person's Age (capped at 85)
      HH_Reference_Age = case_when(
        DMDHRAGE >= 18 & DMDHRAGE <= 84 ~ as.numeric(DMDHRAGE),
        DMDHRAGE == 85 ~ 85,
        TRUE ~ NA_real_
      ),
      
      # HH Reference Person's Education
      HH_Reference_Education_Level = case_when(
        DMDHREDU == 1 ~ "Less Than 9th Grade",
        DMDHREDU == 2 ~ "9-11th Grade",
        DMDHREDU == 3 ~ "High School Graduate/GED",
        DMDHREDU == 4 ~ "Some College or AA Degree", 
        DMDHREDU == 5 ~ "College Graduate or Above",
        TRUE ~ NA_character_
      ),
      
      # HH Reference Person's Marital Status
      HH_Reference_Marital_Status = case_when(
        DMDHRMAR == 1 ~ "Married",
        DMDHRMAR == 2 ~ "Widowed",
        DMDHRMAR == 3 ~ "Divorced",
        DMDHRMAR == 4 ~ "Separated", 
        DMDHRMAR == 5 ~ "Never Married",
        DMDHRMAR == 6 ~ "Living with Partner",
        TRUE ~ NA_character_
      ),
      
      # ═══════════════════════════════════════════════════════
      # INCOME CATEGORIES (DETAILED)
      # ═══════════════════════════════════════════════════════
      
      # Annual Household Income (detailed categories)
      Annual_Household_Income = case_when(
        INDHHIN2 == 1 ~ "$0 to $4,999",
        INDHHIN2 == 2 ~ "$5,000 to $9,999", 
        INDHHIN2 == 3 ~ "$10,000 to $14,999",
        INDHHIN2 == 4 ~ "$15,000 to $19,999",
        INDHHIN2 == 5 ~ "$20,000 to $24,999",
        INDHHIN2 == 6 ~ "$25,000 to $34,999",
        INDHHIN2 == 7 ~ "$35,000 to $44,999",
        INDHHIN2 == 8 ~ "$45,000 to $54,999",
        INDHHIN2 == 9 ~ "$55,000 to $64,999",
        INDHHIN2 == 10 ~ "$65,000 to $74,999",
        INDHHIN2 == 12 ~ "$20,000 and Over",
        INDHHIN2 == 13 ~ "Under $20,000", 
        INDHHIN2 == 14 ~ "$75,000 to $99,999",
        INDHHIN2 == 15 ~ "$100,000 and Over",
        TRUE ~ NA_character_
      ),
      
      # Annual Family Income (detailed categories)
      Annual_Family_Income = case_when(
        INDFMIN2 == 1 ~ "$0 to $4,999",
        INDFMIN2 == 2 ~ "$5,000 to $9,999",
        INDFMIN2 == 3 ~ "$10,000 to $14,999", 
        INDFMIN2 == 4 ~ "$15,000 to $19,999",
        INDFMIN2 == 5 ~ "$20,000 to $24,999",
        INDFMIN2 == 6 ~ "$25,000 to $34,999",
        INDFMIN2 == 7 ~ "$35,000 to $44,999",
        INDFMIN2 == 8 ~ "$45,000 to $54,999",
        INDFMIN2 == 9 ~ "$55,000 to $64,999",
        INDFMIN2 == 10 ~ "$65,000 to $74,999",
        INDFMIN2 == 12 ~ "$20,000 and Over",
        INDFMIN2 == 13 ~ "Under $20,000",
        INDFMIN2 == 14 ~ "$75,000 to $99,999", 
        INDFMIN2 == 15 ~ "$100,000 and Over",
        TRUE ~ NA_character_
      )
    )
  
  return(demo_derived)
}

# ── QUARTILE CALCULATION FUNCTIONS ─────────────────────────

calculate_age_quartiles <- function(demo_data) {
  log_message("Calculating age quartiles...")
  
  # Calculate quartiles from actual data
  age_quartiles <- quantile(demo_data$Age, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  
  log_message("  Age quartile breaks: ", paste(round(age_quartiles, 1), collapse = ", "))
  
  # Create quartile variable
  demo_data$AgeQuartile <- with(demo_data, case_when(
    is.na(Age) ~ "Missing",
    TRUE ~ as.character(cut(
      Age,
      breaks = age_quartiles,
      include.lowest = TRUE,
      labels = c("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%")
    ))
  ))
  
  return(demo_data)
}

calculate_pir_quartiles <- function(demo_data) {
  log_message("Calculating PIR quartiles...")
  
  # Calculate quartiles from actual data
  pir_quartiles <- quantile(demo_data$Ratio_Family_Income_Poverty, 
                           probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  
  log_message("  PIR quartile breaks: ", paste(round(pir_quartiles, 2), collapse = ", "))
  
  # Create quartile variable
  demo_data$PIRQuartile <- with(demo_data, case_when(
    is.na(Ratio_Family_Income_Poverty) ~ "Missing",
    TRUE ~ as.character(cut(
      Ratio_Family_Income_Poverty,
      breaks = pir_quartiles,
      include.lowest = TRUE,
      labels = c("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%")
    ))
  ))
  
  return(demo_data)
}

# ── FACTOR APPLICATION FUNCTIONS ───────────────────────────

apply_demographic_factors <- function(demo_data) {
  log_message("🏷️ Applying proper factor levels...")
  
  demo_data %>%
    mutate(
      # ═══════════════════════════════════════════════════════
      # BASIC CATEGORICAL VARIABLES
      # ═══════════════════════════════════════════════════════
      
      # Gender
      Gender = factor(Gender, levels = c("Male", "Female")),
      
      # Age groups (ordered)
      AgeGroup = factor(AgeGroup, 
                        levels = c("14-19", "20-29", "30-39", "40-49", "50-59", "60-69"),
                        ordered = TRUE),
      
      age_group = factor(age_group,
                         levels = c("14–19", "20–29", "30–39", "40–49", "50–59", "60–69"),
                         ordered = TRUE),
      
      # ═══════════════════════════════════════════════════════
      # EDUCATION VARIABLES (ORDERED)
      # ═══════════════════════════════════════════════════════
      
      EducationLevel = factor(EducationLevel,
                              levels = c("< 9th Grade", "9-11th Grade", "High School", 
                                         "College/AA", "College Graduate"),
                              ordered = TRUE),
      
      Education_Level = factor(Education_Level,
                               levels = c("< 9th Grade", "9-11th Grade", "High School", 
                                          "College/AA", "College Graduate"),
                               ordered = TRUE),
      
      EducationLevel_Simple = factor(EducationLevel_Simple,
                                     levels = c("< 9th Grade", "9-11th Grade", "High School", 
                                                "College/AA", "College Graduate"),
                                     ordered = TRUE),
      
      # ═══════════════════════════════════════════════════════
      # ETHNICITY & CITIZENSHIP
      # ═══════════════════════════════════════════════════════
      
      Ethnicity = factor(Ethnicity, 
                         levels = c("White", "Black", "Mexican", "Other Hispanic", "Other")),
      
      BornInUSA = factor(BornInUSA, levels = c("US Born", "non-US Born")),
      US_Born = factor(US_Born, levels = c("US Born", "non-US Born")),
      US_Citizen = factor(US_Citizen, levels = c("US Born", "non-US Born")),
      
      # ═══════════════════════════════════════════════════════
      # INCOME & POVERTY (ORDERED)
      # ═══════════════════════════════════════════════════════
      
      PIRcat = factor(PIRcat, 
                      levels = c("<1", "1-2", "2-4", ">4"),
                      ordered = TRUE),
      
      PIRQuartile = factor(PIRQuartile,
                           levels = c("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%", "Missing"),
                           ordered = TRUE),
      
      AgeQuartile = factor(AgeQuartile,
                           levels = c("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%", "Missing"),
                           ordered = TRUE),
      
      # ═══════════════════════════════════════════════════════
      # HOUSEHOLD VARIABLES
      # ═══════════════════════════════════════════════════════
      
      Household_Size_Factor = factor(Household_Size_Factor,
                                     levels = c("1", "2", "3", "4", "5", "6", "7+"),
                                     ordered = TRUE),
      
      # ═══════════════════════════════════════════════════════
      # MARITAL STATUS
      # ═══════════════════════════════════════════════════════
      
      Marital_Status = factor(Marital_Status,
                              levels = c("Married", "Widowed", "Divorced", "Separated", 
                                         "Never Married", "Living with Partner")),
      
      # ═══════════════════════════════════════════════════════
      # INTERVIEW VARIABLES
      # ═══════════════════════════════════════════════════════
      
      Interview_Language = factor(Interview_Language, levels = c("English", "Spanish")),
      Interview_Proxy = factor(Interview_Proxy, levels = c("No", "Yes")),
      Interview_Interpreter = factor(Interview_Interpreter, levels = c("No", "Yes")),
      
      # ═══════════════════════════════════════════════════════
      # HOUSEHOLD REFERENCE PERSON VARIABLES
      # ═══════════════════════════════════════════════════════
      
      HH_Reference_Gender = factor(HH_Reference_Gender, levels = c("Male", "Female")),
      
      HH_Reference_Education_Level = factor(HH_Reference_Education_Level,
                                            levels = c("Less Than 9th Grade", "9-11th Grade",
                                                       "High School Graduate/GED", 
                                                       "Some College or AA Degree",
                                                       "College Graduate or Above"),
                                            ordered = TRUE),
      
      HH_Reference_Marital_Status = factor(HH_Reference_Marital_Status,
                                           levels = c("Married", "Widowed", "Divorced", 
                                                      "Separated", "Never Married", 
                                                      "Living with Partner")),
      
      # ═══════════════════════════════════════════════════════
      # DETAILED INCOME CATEGORIES (ORDERED)
      # ═══════════════════════════════════════════════════════
      
      Annual_Household_Income = factor(Annual_Household_Income,
                                       levels = c("$0 to $4,999", "$5,000 to $9,999", 
                                                  "$10,000 to $14,999", "$15,000 to $19,999",
                                                  "Under $20,000", "$20,000 to $24,999", 
                                                  "$20,000 and Over", "$25,000 to $34,999",
                                                  "$35,000 to $44,999", "$45,000 to $54,999",
                                                  "$55,000 to $64,999", "$65,000 to $74,999",
                                                  "$75,000 to $99,999", "$100,000 and Over"),
                                       ordered = TRUE),
      
      Annual_Family_Income = factor(Annual_Family_Income,
                                    levels = c("$0 to $4,999", "$5,000 to $9,999",
                                               "$10,000 to $14,999", "$15,000 to $19,999",
                                               "Under $20,000", "$20,000 to $24,999",
                                               "$20,000 and Over", "$25,000 to $34,999",
                                               "$35,000 to $44,999", "$45,000 to $54,999",
                                               "$55,000 to $64,999", "$65,000 to $74,999",
                                               "$75,000 to $99,999", "$100,000 and Over"),
                                    ordered = TRUE),
      
      # Additional factor variables
      cycle = factor(cycle, levels = c("F", "G")),
      Release_Cycle = factor(Release_Cycle, levels = c("2009-2010", "2011-2012"))
    )
}

# ── METADATA DEFINITIONS ───────────────────────────────────

create_demographic_metadata <- function(cycle) {
  begin_year <- if(cycle == "F") 2009 else 2011
  end_year <- if(cycle == "F") 2010 else 2012
  table_name <- paste0("DEMO_", cycle)
  
  # Define all new variables with descriptions
  new_variables <- tribble(
    ~Variable.Name, ~Variable.Description, ~Data.Type,
    
    # Core identifiers
    "sample", "Sample ID (alternative name for SEQN)", "categorical",
    "cycle", "NHANES Cycle (F=2009-2010, G=2011-2012)", "categorical",
    "Release_Cycle", "Release cycle with full year labels", "categorical",
    "RIAGENDR_01", "Gender (0=Male, 1=Female)", "binary",
    "Gender", "Gender (categorical)", "categorical", 
    "AgeGroup", "Age groups (14-19, 20-29, 30-39, 40-49, 50-59, 60-69)", "categorical",
    "age_group", "Age groups with en-dash formatting (14–19, 20–29, etc.)", "categorical",
    "Age", "Age in years (capped at 85)", "continuous",
    "AgeQuartile", "Age quartiles (Q1-Q4)", "categorical",
    
    # Education variables
    "EducationLevel", "Education level (5 categories)", "categorical",
    "Education_Level", "Education level (alternative naming)", "categorical",
    "EducationLevel_Simple", "Education level with simplified phyloseq naming", "categorical",
    
    # Ethnicity and citizenship
    "Ethnicity", "Race/ethnicity (5 categories)", "categorical", 
    "BornInUSA", "Born in USA status", "categorical",
    "US_Born", "US born status (alternative naming)", "categorical",
    "US_Citizen", "US citizenship status (alternative naming)", "categorical",
    
    # Income and poverty
    "PIRcat", "Poverty income ratio categories (<1, 1-2, 2-4, >4)", "categorical",
    "PIRQuartile", "PIR quartiles (Q1-Q4)", "categorical",
    "Ratio_Family_Income_Poverty", "Family income to poverty ratio (continuous)", "continuous",
    
    # Household variables
    "Household_Size", "Household size (capped at 7)", "continuous",
    "Family_Size", "Family size (capped at 7)", "continuous", 
    "Household_Size_Factor", "Household size as factor (1,2,3,4,5,6,7+)", "categorical",
    
    # Marital status
    "Marital_Status", "Marital status (6 categories)", "categorical",
    
    # Interview variables
    "Interview_Language", "Language of interview (English/Spanish)", "categorical",
    "Interview_Proxy", "Proxy used in interview (Yes/No)", "categorical",
    "Interview_Interpreter", "Interpreter used in interview (Yes/No)", "categorical",
    
    # Household reference person
    "HH_Reference_Gender", "Household reference person gender", "categorical",
    "HH_Reference_Age", "Household reference person age", "continuous",
    "HH_Reference_Education_Level", "Household reference person education", "categorical",
    "HH_Reference_Marital_Status", "Household reference person marital status", "categorical",
    
    # Detailed income categories  
    "Annual_Household_Income", "Annual household income (detailed categories)", "categorical",
    "Annual_Family_Income", "Annual family income (detailed categories)", "categorical"
  ) %>%
    mutate(
      Use.Constraints = "None",
      Data.File.Name = table_name,
      Data.File.Description = paste0("Demographics (", cycle, ") - Fully Processed with Factors and Quartiles"),
      Begin.Year = begin_year,
      EndYear = end_year,
      Component = "Demographics"
    )
  
  return(new_variables)
}

# ── COMPLETE PROCESSING PIPELINE ───────────────────────────

process_demographics_complete <- function(demo_data, cycle) {
  log_message("🔄 Running complete processing pipeline for cycle ", cycle, "...")
  
  demo_data %>%
    derive_demographic_variables(cycle) %>%
    calculate_age_quartiles() %>%
    calculate_pir_quartiles() %>%
    apply_demographic_factors()
}

# ── MAIN EXECUTION ─────────────────────────────────────────

log_message("Starting comprehensive demographic processing...")

# Create output database (copy from input first)
if (!opt$test) {
  log_message("📂 Creating output database...")
  dir_create(dirname(opt$output_db), recurse = TRUE)
  file_copy(opt$input_db, opt$output_db, overwrite = TRUE)
  log_message("✅ Output database created: ", opt$output_db)
  
  # Connect to output database
  con_output <- dbConnect(SQLite(), opt$output_db)
  log_message("✅ Output database connection established")
} else {
  log_message("ℹ️ TEST MODE: Skipping output database creation")
  con_output <- con_input  # Use input connection for testing
}

# Process DEMO_F and DEMO_G tables
for (cycle in c("F", "G")) {
  table_name <- paste0("DEMO_", cycle)
  log_message("\nProcessing ", table_name, "...")
  
  # Check if table exists
  if (!table_name %in% dbListTables(con_input)) {
    log_message("⚠️ Table ", table_name, " not found, skipping...")
    next
  }
  
  # Read original data
  original_data <- dbReadTable(con_input, table_name)
  log_message("  Original data: ", nrow(original_data), " rows, ", ncol(original_data), " columns")
  
  # Apply complete processing pipeline
  processed_data <- process_demographics_complete(original_data, cycle)
  log_message("  ✅ Processed data: ", nrow(processed_data), " rows, ", ncol(processed_data), " columns")
  log_message("  Added ", ncol(processed_data) - ncol(original_data), " new variables")
  
  # Update output database (if not in test mode)
  if (!opt$test) {
    dbBegin(con_output)
    tryCatch({
      dbWriteTable(con_output, table_name, processed_data, overwrite = TRUE)
      log_message("  ✅ Updated ", table_name, " in output database")
      dbCommit(con_output)
    }, error = function(e) {
      dbRollback(con_output)
      log_message("❌ Error updating ", table_name, ": ", e$message)
      stop("❌ Error updating ", table_name, ": ", e$message)
    })
  } else {
    log_message("  ℹ️ TEST MODE: Skipping database update")
  }
  
  # Create and add metadata
  new_metadata <- create_demographic_metadata(cycle)
  
  # Filter to only new variables (not in existing metadata)
  if (!opt$test) {
    existing_vars <- dbGetQuery(con_output, 
      sprintf("SELECT \"Variable.Name\" FROM variable_names_epcf WHERE \"Data.File.Name\" = '%s'", table_name)
    )$Variable.Name
    
    vars_to_add <- anti_join(new_metadata, 
                             tibble(Variable.Name = existing_vars),
                             by = "Variable.Name")
    
    if (nrow(vars_to_add) > 0) {
      dbBegin(con_output)
      tryCatch({
        # Read existing metadata
        existing_metadata <- dbReadTable(con_output, "variable_names_epcf")
        # Combine with new metadata
        updated_metadata <- bind_rows(existing_metadata, vars_to_add)
        # Write back to database
        dbWriteTable(con_output, "variable_names_epcf", updated_metadata, overwrite = TRUE)
        log_message("  📝 Added ", nrow(vars_to_add), " variable metadata entries")
        dbCommit(con_output)
      }, error = function(e) {
        dbRollback(con_output)
        log_message("⚠️ Error updating metadata for ", table_name, ": ", e$message)
      })
    } else {
      log_message("  ℹ️ All variable metadata already exists")
    }
  } else {
    log_message("  ℹ️ TEST MODE: Would add ", nrow(new_metadata), " metadata entries")
  }
}

# ── SUMMARY REPORT ─────────────────────────────────────────

log_message("\n", paste(rep("=", 60), collapse=""))
log_message("COMPREHENSIVE DEMOGRAPHIC PROCESSING SUMMARY")
log_message(paste(rep("=", 60), collapse=""))

# Check final status
for (cycle in c("F", "G")) {
  table_name <- paste0("DEMO_", cycle)
  if (table_name %in% dbListTables(if(!opt$test) con_output else con_input)) {
    final_data <- dbGetQuery(if(!opt$test) con_output else con_input, 
                             sprintf("SELECT COUNT(*) as row_count FROM %s", table_name))
    final_cols <- length(dbListFields(if(!opt$test) con_output else con_input, table_name))
    log_message("✅ ", table_name, ": ", final_data$row_count, " rows, ", final_cols, " columns")
  }
}

# Count total variables in metadata
total_vars <- dbGetQuery(if(!opt$test) con_output else con_input, 
                         "SELECT COUNT(*) as count FROM variable_names_epcf")$count
log_message("📝 Total variables in metadata: ", total_vars)

# Create summary report file
summary_file <- file.path(logs_dir, paste0("process_demographics_complete_summary_", timestamp, ".txt"))
writeLines(c(
  "COMPREHENSIVE DEMOGRAPHIC PROCESSING SUMMARY",
  paste("Timestamp:", Sys.time()),
  paste("Input database:", opt$input_db),
  paste("Output database:", opt$output_db),
  paste("Test mode:", ifelse(opt$test, "YES", "NO")),
  "",
  "PROCESSING STEPS COMPLETED:",
  paste("1. Demographic variable derivation (31 variables)"),
  paste("2. Age quartile calculation"),
  paste("3. PIR quartile calculation"),
  paste("4. Factor level application"),
  paste("5. Metadata registration"),
  "",
  "VARIABLES CREATED:",
  paste("- cycle/Release_Cycle (NHANES Cycle identifiers)"),
  paste("- sample/SEQN (Sample identifiers)"),
  paste("- RIAGENDR_01/Gender (Gender variables)"),
  paste("- AgeGroup/age_group/Age/AgeQuartile (Age variables)"),
  paste("- EducationLevel/Education_Level/EducationLevel_Simple (Education)"),
  paste("- Ethnicity/BornInUSA/US_Born/US_Citizen (Ethnicity/Citizenship)"),
  paste("- PIRcat/PIRQuartile/Ratio_Family_Income_Poverty (Income/Poverty)"),
  paste("- Household_Size/Family_Size/Household_Size_Factor (Household)"),
  paste("- Marital_Status (Marital status)"),
  paste("- Interview_Language/Proxy/Interpreter (Interview variables)"),
  paste("- HH_Reference_* variables (Reference person info)"),
  paste("- Annual_Household_Income/Annual_Family_Income (Detailed income)"),
  "",
  "STATUS:",
  if(opt$test) "TEST MODE - No output database created" else paste("APPLIED TO:", opt$output_db)
), summary_file)

log_message("📄 Summary report saved: ", summary_file)

# Close database connections
dbDisconnect(con_input)
if (!opt$test && exists("con_output")) {
  dbDisconnect(con_output)
}

# Final messages before closing log
log_message(paste(rep("=", 60), collapse=""))
log_message("✅ Comprehensive demographic processing completed!")
log_message("📁 Log file saved: ", log_file)
if (!opt$test) {
  log_message("📂 Processed database saved: ", opt$output_db)
}
log_message(paste(rep("=", 60), collapse=""))

close(log_con)

# Final console message
message("✅ All logs saved to: ", logs_dir)
if (!opt$test) {
  message("📂 Processed database available at: ", opt$output_db)
} 