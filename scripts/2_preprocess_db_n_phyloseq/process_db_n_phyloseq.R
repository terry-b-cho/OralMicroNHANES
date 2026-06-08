#!/usr/bin/env Rscript

# =============================================================================
# DATABASE PREPROCESSING + PHYLOSEQ OBJECT CREATION
# =============================================================================
# Stage 1: derive comprehensive demographic variables from DEMO_F/DEMO_G in the
# input NHANES SQLite database and write them back into a processed copy of
# the DB at:
#   data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite
#
# Stage 2: read the processed DB and build six phyloseq objects (counts,
# relative, relative_none, relative_clr, relative_lognorm, relative_hellinger),
# saved as RDS files under:
#   results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/
#
# Environment: R >= 4.5 with DBI, RSQLite, dplyr, tidyr, tibble, purrr,
# stringr, forcats, fs, ape, data.tree, phyloseq.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(forcats)
  library(fs)
  library(ape)
  library(data.tree)
  library(phyloseq)
})

base_path <- PROJECT_ROOT
input_db  <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite")
output_db <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
out_root  <- file.path(base_path, "results/analyses_results/2_preprocess_db_n_phyloseq_out")
intermediate_files_path <- file.path(out_root, "intermediate")
viz_out_path            <- file.path(out_root, "figures_out")
dir_create(intermediate_files_path, recurse = TRUE)
dir_create(viz_out_path,            recurse = TRUE)

# =============================================================================
# STAGE 1 - DEMOGRAPHIC PROCESSING
# =============================================================================

derive_demographic_variables <- function(demo_data, cycle) {
  demo_data %>%
    mutate(
      SEQN = as.character(SEQN),
      sample = as.character(SEQN),
      cycle = cycle,
      Release_Cycle = case_when(
        cycle == "F" ~ "2009-2010",
        cycle == "G" ~ "2011-2012",
        TRUE ~ NA_character_
      ),
      RIAGENDR_01 = case_when(
        RIAGENDR == 1 ~ 0,
        RIAGENDR == 2 ~ 1,
        TRUE ~ NA_real_
      ),
      Gender = case_when(
        RIAGENDR == 1 ~ "Male",
        RIAGENDR == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      AgeGroup = case_when(
        RIDAGEYR >= 60 & RIDAGEYR <= 69 ~ "60-69",
        RIDAGEYR >= 50 & RIDAGEYR < 60 ~ "50-59",
        RIDAGEYR >= 40 & RIDAGEYR < 50 ~ "40-49",
        RIDAGEYR >= 30 & RIDAGEYR < 40 ~ "30-39",
        RIDAGEYR >= 20 & RIDAGEYR < 30 ~ "20-29",
        RIDAGEYR >= 14 & RIDAGEYR < 20 ~ "14-19",
        TRUE ~ NA_character_
      ),
      age_group = case_when(
        RIDAGEYR >= 60 & RIDAGEYR <= 69 ~ "60–69",
        RIDAGEYR >= 50 & RIDAGEYR < 60 ~ "50–59",
        RIDAGEYR >= 40 & RIDAGEYR < 50 ~ "40–49",
        RIDAGEYR >= 30 & RIDAGEYR < 40 ~ "30–39",
        RIDAGEYR >= 20 & RIDAGEYR < 30 ~ "20–29",
        RIDAGEYR >= 14 & RIDAGEYR < 20 ~ "14–19",
        TRUE ~ NA_character_
      ),
      Age = case_when(
        RIDAGEYR > 85 ~ 85,
        RIDAGEYR >= 0 ~ RIDAGEYR,
        TRUE ~ NA_real_
      ),
      EducationLevel = case_when(
        EDUCATION_LESS9 == 1 ~ "< 9th Grade",
        EDUCATION_9_11 == 1 ~ "9-11th Grade",
        EDUCATION_HSGRAD == 1 ~ "High School",
        EDUCATION_AA == 1 ~ "College/AA",
        EDUCATION_COLLEGEGRAD == 1 ~ "College Graduate",
        TRUE ~ NA_character_
      ),
      Education_Level = EducationLevel,
      EducationLevel_Simple = case_when(
        EDUCATION_HSGRAD == 1 ~ "High School",
        EDUCATION_AA == 1 ~ "College/AA",
        EDUCATION_COLLEGEGRAD == 1 ~ "College Graduate",
        EDUCATION_LESS9 == 1 ~ "< 9th Grade",
        EDUCATION_9_11 == 1 ~ "9-11th Grade",
        TRUE ~ NA_character_
      ),
      Ethnicity = case_when(
        ETHNICITY_NONHISPANICWHITE == 1 ~ "White",
        ETHNICITY_NONHISPANICBLACK == 1 ~ "Black",
        ETHNICITY_MEXICAN == 1 ~ "Mexican",
        ETHNICITY_OTHERHISPANIC == 1 ~ "Other Hispanic",
        ETHNICITY_OTHER == 1 ~ "Other",
        TRUE ~ NA_character_
      ),
      BornInUSA = case_when(
        BORN_INUSA == 1 ~ "US Born",
        BORN_INUSA == 0 ~ "non-US Born",
        TRUE ~ NA_character_
      ),
      US_Born = BornInUSA,
      US_Citizen = BornInUSA,
      PIRcat = case_when(
        INDFMPIR < 1 ~ "<1",
        INDFMPIR >= 1 & INDFMPIR < 2 ~ "1-2",
        INDFMPIR >= 2 & INDFMPIR < 4 ~ "2-4",
        INDFMPIR >= 4 ~ ">4",
        TRUE ~ NA_character_
      ),
      Ratio_Family_Income_Poverty = as.numeric(INDFMPIR),
      Household_Size = case_when(
        DMDHHSIZ %in% 1:6 ~ as.numeric(DMDHHSIZ),
        DMDHHSIZ >= 7 ~ 7,
        DMDHHSIZ %in% c(77, 99) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      Family_Size = case_when(
        DMDFMSIZ %in% 1:6 ~ as.numeric(DMDFMSIZ),
        DMDFMSIZ >= 7 ~ 7,
        DMDFMSIZ %in% c(77, 99) ~ NA_real_,
        TRUE ~ NA_real_
      ),
      Household_Size_Factor = case_when(
        Household_Size <= 6 ~ as.character(Household_Size),
        Household_Size >= 7 ~ "7+",
        TRUE ~ NA_character_
      ),
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
      Interview_Language = case_when(
        SIALANG == 1 ~ "English",
        SIALANG == 2 ~ "Spanish",
        SIALANG %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      Interview_Proxy = case_when(
        SIAPROXY == 1 ~ "Yes",
        SIAPROXY == 2 ~ "No",
        SIAPROXY %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      Interview_Interpreter = case_when(
        SIAINTRP == 1 ~ "Yes",
        SIAINTRP == 2 ~ "No",
        SIAINTRP %in% c(77, 99) ~ NA_character_,
        TRUE ~ NA_character_
      ),
      HH_Reference_Gender = case_when(
        DMDHRGND == 1 ~ "Male",
        DMDHRGND == 2 ~ "Female",
        TRUE ~ NA_character_
      ),
      HH_Reference_Age = case_when(
        DMDHRAGE >= 18 & DMDHRAGE <= 84 ~ as.numeric(DMDHRAGE),
        DMDHRAGE == 85 ~ 85,
        TRUE ~ NA_real_
      ),
      HH_Reference_Education_Level = case_when(
        DMDHREDU == 1 ~ "Less Than 9th Grade",
        DMDHREDU == 2 ~ "9-11th Grade",
        DMDHREDU == 3 ~ "High School Graduate/GED",
        DMDHREDU == 4 ~ "Some College or AA Degree",
        DMDHREDU == 5 ~ "College Graduate or Above",
        TRUE ~ NA_character_
      ),
      HH_Reference_Marital_Status = case_when(
        DMDHRMAR == 1 ~ "Married",
        DMDHRMAR == 2 ~ "Widowed",
        DMDHRMAR == 3 ~ "Divorced",
        DMDHRMAR == 4 ~ "Separated",
        DMDHRMAR == 5 ~ "Never Married",
        DMDHRMAR == 6 ~ "Living with Partner",
        TRUE ~ NA_character_
      ),
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
}

create_demographic_dummy_variables <- function(demo_data) {
  cols <- names(demo_data)
  d <- demo_data
  if ("RIAGENDR" %in% cols) {
    d <- d %>% mutate(RIAGENDR2 = case_when(
      RIAGENDR == 2 ~ 1,
      RIAGENDR == 1 ~ 0,
      TRUE ~ NA_real_
    ))
  }
  if ("DMDEDUC2" %in% cols) {
    d <- d %>% mutate(
      DMDEDUC22 = case_when(DMDEDUC2 == 2 ~ 1, DMDEDUC2 %in% c(1,3,4,5) ~ 0, TRUE ~ NA_real_),
      DMDEDUC23 = case_when(DMDEDUC2 == 3 ~ 1, DMDEDUC2 %in% c(1,2,4,5) ~ 0, TRUE ~ NA_real_),
      DMDEDUC24 = case_when(DMDEDUC2 == 4 ~ 1, DMDEDUC2 %in% c(1,2,3,5) ~ 0, TRUE ~ NA_real_),
      DMDEDUC25 = case_when(DMDEDUC2 == 5 ~ 1, DMDEDUC2 %in% c(1,2,3,4) ~ 0, TRUE ~ NA_real_)
    )
  }
  if ("RIDRETH1" %in% cols) {
    d <- d %>% mutate(
      RIDRETH12 = case_when(RIDRETH1 == 2 ~ 1, RIDRETH1 %in% c(1,3,4,5) ~ 0, TRUE ~ NA_real_),
      RIDRETH13 = case_when(RIDRETH1 == 3 ~ 1, RIDRETH1 %in% c(1,2,4,5) ~ 0, TRUE ~ NA_real_),
      RIDRETH14 = case_when(RIDRETH1 == 4 ~ 1, RIDRETH1 %in% c(1,2,3,5) ~ 0, TRUE ~ NA_real_),
      RIDRETH15 = case_when(RIDRETH1 == 5 ~ 1, RIDRETH1 %in% c(1,2,3,4) ~ 0, TRUE ~ NA_real_)
    )
  }
  if ("DMDMARTL" %in% cols) {
    d <- d %>% mutate(
      DMDMARTL2 = case_when(DMDMARTL == 2 ~ 1, DMDMARTL %in% c(1,3,4,5,6) ~ 0, TRUE ~ NA_real_),
      DMDMARTL3 = case_when(DMDMARTL == 3 ~ 1, DMDMARTL %in% c(1,2,4,5,6) ~ 0, TRUE ~ NA_real_),
      DMDMARTL4 = case_when(DMDMARTL == 4 ~ 1, DMDMARTL %in% c(1,2,3,5,6) ~ 0, TRUE ~ NA_real_),
      DMDMARTL5 = case_when(DMDMARTL == 5 ~ 1, DMDMARTL %in% c(1,2,3,4,6) ~ 0, TRUE ~ NA_real_),
      DMDMARTL6 = case_when(DMDMARTL == 6 ~ 1, DMDMARTL %in% c(1,2,3,4,5) ~ 0, TRUE ~ NA_real_)
    )
  }
  if ("DMDBORN2" %in% cols) {
    d <- d %>% mutate(DMDBORN_FOREIGN = case_when(
      DMDBORN2 == 2 ~ 1,
      DMDBORN2 == 1 ~ 0,
      TRUE ~ NA_real_
    ))
  } else if ("DMDBORN4" %in% cols) {
    d <- d %>% mutate(DMDBORN_FOREIGN = case_when(
      DMDBORN4 == 2 ~ 1,
      DMDBORN4 == 1 ~ 0,
      TRUE ~ NA_real_
    ))
  }
  if ("SIALANG" %in% cols) {
    d <- d %>% mutate(SIALANG2 = case_when(
      SIALANG == 2 ~ 1,
      SIALANG == 1 ~ 0,
      TRUE ~ NA_real_
    ))
  }
  if ("DMDHHSIZ" %in% cols) {
    d <- d %>% mutate(
      DMDHHSIZ2 = case_when(DMDHHSIZ == 2 ~ 1, DMDHHSIZ %in% c(1,3,4,5,6) | DMDHHSIZ >= 7 ~ 0, TRUE ~ NA_real_),
      DMDHHSIZ3 = case_when(DMDHHSIZ == 3 ~ 1, DMDHHSIZ %in% c(1,2,4,5,6) | DMDHHSIZ >= 7 ~ 0, TRUE ~ NA_real_),
      DMDHHSIZ4 = case_when(DMDHHSIZ == 4 ~ 1, DMDHHSIZ %in% c(1,2,3,5,6) | DMDHHSIZ >= 7 ~ 0, TRUE ~ NA_real_),
      DMDHHSIZ5 = case_when(DMDHHSIZ == 5 ~ 1, DMDHHSIZ %in% c(1,2,3,4,6) | DMDHHSIZ >= 7 ~ 0, TRUE ~ NA_real_),
      DMDHHSIZ6 = case_when(DMDHHSIZ == 6 ~ 1, DMDHHSIZ %in% c(1,2,3,4,5) | DMDHHSIZ >= 7 ~ 0, TRUE ~ NA_real_),
      DMDHHSIZ7 = case_when(DMDHHSIZ >= 7 ~ 1, DMDHHSIZ %in% c(1,2,3,4,5,6) ~ 0, TRUE ~ NA_real_)
    )
  }
  d
}

calculate_age_quartiles <- function(demo_data) {
  q <- quantile(demo_data$Age, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  demo_data$AgeQuartile <- with(demo_data, case_when(
    is.na(Age) ~ "Missing",
    TRUE ~ as.character(cut(
      Age, breaks = q, include.lowest = TRUE,
      labels = c("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%")
    ))
  ))
  demo_data
}

calculate_pir_quartiles <- function(demo_data) {
  q <- quantile(demo_data$Ratio_Family_Income_Poverty,
                probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  demo_data$PIRQuartile <- with(demo_data, case_when(
    is.na(Ratio_Family_Income_Poverty) ~ "Missing",
    TRUE ~ as.character(cut(
      Ratio_Family_Income_Poverty, breaks = q, include.lowest = TRUE,
      labels = c("Q1: 0–25%", "Q2: 25–50%", "Q3: 50–75%", "Q4: 75–100%")
    ))
  ))
  demo_data
}

apply_demographic_factors <- function(demo_data) {
  demo_data %>%
    mutate(
      Gender = factor(Gender, levels = c("Male", "Female")),
      AgeGroup = factor(AgeGroup, levels = c("14-19","20-29","30-39","40-49","50-59","60-69"), ordered = TRUE),
      age_group = factor(age_group, levels = c("14–19","20–29","30–39","40–49","50–59","60–69"), ordered = TRUE),
      EducationLevel = factor(EducationLevel, levels = c("< 9th Grade","9-11th Grade","High School","College/AA","College Graduate"), ordered = TRUE),
      Education_Level = factor(Education_Level, levels = c("< 9th Grade","9-11th Grade","High School","College/AA","College Graduate"), ordered = TRUE),
      EducationLevel_Simple = factor(EducationLevel_Simple, levels = c("< 9th Grade","9-11th Grade","High School","College/AA","College Graduate"), ordered = TRUE),
      Ethnicity = factor(Ethnicity, levels = c("White","Black","Mexican","Other Hispanic","Other")),
      BornInUSA = factor(BornInUSA, levels = c("US Born","non-US Born")),
      US_Born = factor(US_Born, levels = c("US Born","non-US Born")),
      US_Citizen = factor(US_Citizen, levels = c("US Born","non-US Born")),
      PIRcat = factor(PIRcat, levels = c("<1","1-2","2-4",">4"), ordered = TRUE),
      PIRQuartile = factor(PIRQuartile, levels = c("Q1: 0–25%","Q2: 25–50%","Q3: 50–75%","Q4: 75–100%","Missing"), ordered = TRUE),
      AgeQuartile = factor(AgeQuartile, levels = c("Q1: 0–25%","Q2: 25–50%","Q3: 50–75%","Q4: 75–100%","Missing"), ordered = TRUE),
      Household_Size_Factor = factor(Household_Size_Factor, levels = c("1","2","3","4","5","6","7+"), ordered = TRUE),
      Marital_Status = factor(Marital_Status, levels = c("Married","Widowed","Divorced","Separated","Never Married","Living with Partner")),
      Interview_Language = factor(Interview_Language, levels = c("English","Spanish")),
      Interview_Proxy = factor(Interview_Proxy, levels = c("No","Yes")),
      Interview_Interpreter = factor(Interview_Interpreter, levels = c("No","Yes")),
      HH_Reference_Gender = factor(HH_Reference_Gender, levels = c("Male","Female")),
      HH_Reference_Education_Level = factor(HH_Reference_Education_Level, levels = c("Less Than 9th Grade","9-11th Grade","High School Graduate/GED","Some College or AA Degree","College Graduate or Above"), ordered = TRUE),
      HH_Reference_Marital_Status = factor(HH_Reference_Marital_Status, levels = c("Married","Widowed","Divorced","Separated","Never Married","Living with Partner")),
      Annual_Household_Income = factor(Annual_Household_Income, levels = c("$0 to $4,999","$5,000 to $9,999","$10,000 to $14,999","$15,000 to $19,999","Under $20,000","$20,000 to $24,999","$20,000 and Over","$25,000 to $34,999","$35,000 to $44,999","$45,000 to $54,999","$55,000 to $64,999","$65,000 to $74,999","$75,000 to $99,999","$100,000 and Over"), ordered = TRUE),
      Annual_Family_Income = factor(Annual_Family_Income, levels = c("$0 to $4,999","$5,000 to $9,999","$10,000 to $14,999","$15,000 to $19,999","Under $20,000","$20,000 to $24,999","$20,000 and Over","$25,000 to $34,999","$35,000 to $44,999","$45,000 to $54,999","$55,000 to $64,999","$65,000 to $74,999","$75,000 to $99,999","$100,000 and Over"), ordered = TRUE),
      cycle = factor(cycle, levels = c("F","G")),
      Release_Cycle = factor(Release_Cycle, levels = c("2009-2010","2011-2012"))
    )
}

create_demographic_metadata <- function(cycle) {
  begin_year <- if (cycle == "F") 2009 else 2011
  end_year   <- if (cycle == "F") 2010 else 2012
  table_name <- paste0("DEMO_", cycle)
  tribble(
    ~Variable.Name, ~Variable.Description, ~Data.Type,
    "sample", "Sample ID (alternative name for SEQN)", "categorical",
    "cycle", "NHANES Cycle (F=2009-2010, G=2011-2012)", "categorical",
    "Release_Cycle", "Release cycle with full year labels", "categorical",
    "RIAGENDR_01", "Gender (0=Male, 1=Female)", "binary",
    "Gender", "Gender (categorical)", "categorical",
    "AgeGroup", "Age groups (14-19, 20-29, 30-39, 40-49, 50-59, 60-69)", "categorical",
    "age_group", "Age groups with en-dash formatting (14–19, 20–29, etc.)", "categorical",
    "Age", "Age in years (capped at 85)", "continuous",
    "AgeQuartile", "Age quartiles (Q1-Q4)", "categorical",
    "EducationLevel", "Education level (5 categories)", "categorical",
    "Education_Level", "Education level (alternative naming)", "categorical",
    "EducationLevel_Simple", "Education level with simplified phyloseq naming", "categorical",
    "Ethnicity", "Race/ethnicity (5 categories)", "categorical",
    "BornInUSA", "Born in USA status", "categorical",
    "US_Born", "US born status (alternative naming)", "categorical",
    "US_Citizen", "US citizenship status (alternative naming)", "categorical",
    "PIRcat", "Poverty income ratio categories (<1, 1-2, 2-4, >4)", "categorical",
    "PIRQuartile", "PIR quartiles (Q1-Q4)", "categorical",
    "Ratio_Family_Income_Poverty", "Family income to poverty ratio (continuous)", "continuous",
    "Household_Size", "Household size (capped at 7)", "continuous",
    "Family_Size", "Family size (capped at 7)", "continuous",
    "Household_Size_Factor", "Household size as factor (1,2,3,4,5,6,7+)", "categorical",
    "Marital_Status", "Marital status (6 categories)", "categorical",
    "Interview_Language", "Language of interview (English/Spanish)", "categorical",
    "Interview_Proxy", "Proxy used in interview (Yes/No)", "categorical",
    "Interview_Interpreter", "Interpreter used in interview (Yes/No)", "categorical",
    "HH_Reference_Gender", "Household reference person gender", "categorical",
    "HH_Reference_Age", "Household reference person age", "continuous",
    "HH_Reference_Education_Level", "Household reference person education", "categorical",
    "HH_Reference_Marital_Status", "Household reference person marital status", "categorical",
    "Annual_Household_Income", "Annual household income (detailed categories)", "categorical",
    "Annual_Family_Income", "Annual family income (detailed categories)", "categorical",
    "RIAGENDR2", "Gender: Female (RIAGENDR=2, ref: Male=1)", "binary",
    "DMDEDUC22", "Education: 9-11th grade (DMDEDUC2=2, ref: <9th=1)", "binary",
    "DMDEDUC23", "Education: High school grad (DMDEDUC2=3, ref: <9th=1)", "binary",
    "DMDEDUC24", "Education: Some college/AA (DMDEDUC2=4, ref: <9th=1)", "binary",
    "DMDEDUC25", "Education: College grad+ (DMDEDUC2=5, ref: <9th=1)", "binary",
    "RIDRETH12", "Ethnicity: Other Hispanic (RIDRETH1=2, ref: Mexican=1)", "binary",
    "RIDRETH13", "Ethnicity: Non-Hispanic White (RIDRETH1=3, ref: Mexican=1)", "binary",
    "RIDRETH14", "Ethnicity: Non-Hispanic Black (RIDRETH1=4, ref: Mexican=1)", "binary",
    "RIDRETH15", "Ethnicity: Other race (RIDRETH1=5, ref: Mexican=1)", "binary",
    "DMDMARTL2", "Marital: Widowed (DMDMARTL=2, ref: Married=1)", "binary",
    "DMDMARTL3", "Marital: Divorced (DMDMARTL=3, ref: Married=1)", "binary",
    "DMDMARTL4", "Marital: Separated (DMDMARTL=4, ref: Married=1)", "binary",
    "DMDMARTL5", "Marital: Never married (DMDMARTL=5, ref: Married=1)", "binary",
    "DMDMARTL6", "Marital: Living with partner (DMDMARTL=6, ref: Married=1)", "binary",
    "DMDBORN_FOREIGN", "Birth: Born elsewhere (cycle-specific, ref: US born)", "binary",
    "SIALANG2", "Interview: Spanish (SIALANG=2, ref: English=1)", "binary",
    "DMDHHSIZ2", "Household size: 2 (DMDHHSIZ=2, ref: 1)", "binary",
    "DMDHHSIZ3", "Household size: 3 (DMDHHSIZ=3, ref: 1)", "binary",
    "DMDHHSIZ4", "Household size: 4 (DMDHHSIZ=4, ref: 1)", "binary",
    "DMDHHSIZ5", "Household size: 5 (DMDHHSIZ=5, ref: 1)", "binary",
    "DMDHHSIZ6", "Household size: 6 (DMDHHSIZ=6, ref: 1)", "binary",
    "DMDHHSIZ7", "Household size: 7+ (DMDHHSIZ≥7, ref: 1)", "binary"
  ) %>%
    mutate(
      Use.Constraints = "None",
      Data.File.Name = table_name,
      Data.File.Description = paste0("Demographics (", cycle, ") - Fully Processed with Factors and Quartiles"),
      Begin.Year = begin_year,
      EndYear = end_year,
      Component = "Demographics"
    )
}

process_demographics_complete <- function(demo_data, cycle) {
  demo_data %>%
    derive_demographic_variables(cycle) %>%
    create_demographic_dummy_variables() %>%
    calculate_age_quartiles() %>%
    calculate_pir_quartiles() %>%
    apply_demographic_factors()
}

message("=== STAGE 1: demographic processing ===")
dir_create(dirname(output_db), recurse = TRUE)
file_copy(input_db, output_db, overwrite = TRUE)
con_input  <- dbConnect(SQLite(), input_db)
con_output <- dbConnect(SQLite(), output_db)
on.exit({ try(dbDisconnect(con_input), silent = TRUE)
          try(dbDisconnect(con_output), silent = TRUE) }, add = TRUE)

for (cycle in c("F", "G")) {
  table_name <- paste0("DEMO_", cycle)
  if (!table_name %in% dbListTables(con_input)) next
  original_data  <- dbReadTable(con_input, table_name)
  processed_data <- process_demographics_complete(original_data, cycle)
  dbBegin(con_output)
  tryCatch({
    dbWriteTable(con_output, table_name, processed_data, overwrite = TRUE)
    dbCommit(con_output)
  }, error = function(e) { dbRollback(con_output); stop(e) })
  new_metadata <- create_demographic_metadata(cycle)
  existing_vars <- dbGetQuery(con_output,
    sprintf("SELECT \"Variable.Name\" FROM variable_names_epcf WHERE \"Data.File.Name\" = '%s'", table_name)
  )$Variable.Name
  vars_to_add <- anti_join(new_metadata, tibble(Variable.Name = existing_vars), by = "Variable.Name")
  if (nrow(vars_to_add) > 0) {
    dbBegin(con_output)
    tryCatch({
      existing_metadata <- dbReadTable(con_output, "variable_names_epcf")
      updated_metadata  <- bind_rows(existing_metadata, vars_to_add)
      dbWriteTable(con_output, "variable_names_epcf", updated_metadata, overwrite = TRUE)
      dbCommit(con_output)
    }, error = function(e) { dbRollback(con_output); warning(e) })
  }
}
dbDisconnect(con_input)
dbDisconnect(con_output)

# =============================================================================
# STAGE 2 - PHYLOSEQ OBJECT CREATION
# =============================================================================

message("=== STAGE 2: phyloseq object creation ===")
con <- dbConnect(SQLite(), output_db)
on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)

create_otu_matrix_from_tables <- function(table_f_name, table_g_name, use_base_names = TRUE) {
  genus_f <- tbl(con, table_f_name) %>% collect()
  genus_g <- tbl(con, table_g_name) %>% collect()
  taxa_f <- colnames(genus_f)[colnames(genus_f) != "SEQN"]
  taxa_g <- colnames(genus_g)[colnames(genus_g) != "SEQN"]
  if (!identical(colnames(genus_f), colnames(genus_g))) {
    all_taxa <- union(taxa_f, taxa_g)
    missing_in_f <- setdiff(taxa_g, taxa_f)
    for (taxa in missing_in_f) genus_f[[taxa]] <- 0
    missing_in_g <- setdiff(taxa_f, taxa_g)
    for (taxa in missing_in_g) genus_g[[taxa]] <- 0
    all_cols <- c("SEQN", all_taxa)
    genus_f <- genus_f %>% select(all_of(all_cols))
    genus_g <- genus_g %>% select(all_of(all_cols))
  }
  genus_combined <- bind_rows(genus_f, genus_g)
  if (use_base_names) {
    genus_combined %>%
      mutate(SEQN = as.character(SEQN)) %>%
      tidyr::pivot_longer(cols = matches("^RSV_genus"), names_to = "otu", values_to = "Abundance") %>%
      mutate(otu_base = gsub("_(count|relative)$", "", otu)) %>%
      select(-otu) %>%
      tidyr::pivot_wider(names_from = SEQN, values_from = Abundance) %>%
      column_to_rownames("otu_base")
  } else {
    genus_combined %>%
      mutate(SEQN = as.character(SEQN)) %>%
      tidyr::pivot_longer(cols = matches("^RSV_genus"), names_to = "otu", values_to = "Abundance") %>%
      tidyr::pivot_wider(names_from = SEQN, values_from = Abundance) %>%
      column_to_rownames("otu")
  }
}

otu_mat_counts             <- create_otu_matrix_from_tables("DADA2RSV_GENUS_COUNT_F",          "DADA2RSV_GENUS_COUNT_G")
otu_mat_relative           <- create_otu_matrix_from_tables("DADA2RSV_GENUS_RELATIVE_F",       "DADA2RSV_GENUS_RELATIVE_G")
otu_mat_relative_none      <- create_otu_matrix_from_tables("DADA2RSV_GENUS_RELATIVE_F_none",      "DADA2RSV_GENUS_RELATIVE_G_none")
otu_mat_relative_clr       <- create_otu_matrix_from_tables("DADA2RSV_GENUS_RELATIVE_F_clr",       "DADA2RSV_GENUS_RELATIVE_G_clr")
otu_mat_relative_lognorm   <- create_otu_matrix_from_tables("DADA2RSV_GENUS_RELATIVE_F_lognorm",   "DADA2RSV_GENUS_RELATIVE_G_lognorm")
otu_mat_relative_hellinger <- create_otu_matrix_from_tables("DADA2RSV_GENUS_RELATIVE_F_hellinger", "DADA2RSV_GENUS_RELATIVE_G_hellinger")

otu_taxa_names  <- rownames(otu_mat_counts)
available_vars  <- tbl(con, "variable_names_epcf") %>% collect() %>% filter(Data.File.Name == "DADA2RSV_GENUS_COUNT_G")
otu_with_suffix <- paste0(otu_taxa_names, "_count")

tax_mat_nhanes <- available_vars %>%
  filter(Variable.Name %in% otu_with_suffix) %>%
  mutate(desc_clean = as.character(ifelse(is.na(Variable.Description) | Variable.Description == "", "", Variable.Description))) %>%
  mutate(
    desc_parts = strsplit(desc_clean, ";"),
    Domain = sapply(desc_parts, function(x) ifelse(length(x) >= 1 && trimws(x[1]) != "", trimws(x[1]), NA_character_)),
    Phylum = sapply(desc_parts, function(x) ifelse(length(x) >= 2 && trimws(x[2]) != "", trimws(x[2]), NA_character_)),
    Class  = sapply(desc_parts, function(x) ifelse(length(x) >= 3 && trimws(x[3]) != "", trimws(x[3]), NA_character_)),
    Order  = sapply(desc_parts, function(x) ifelse(length(x) >= 4 && trimws(x[4]) != "", trimws(x[4]), NA_character_)),
    Family = sapply(desc_parts, function(x) ifelse(length(x) >= 5 && trimws(x[5]) != "", trimws(x[5]), NA_character_)),
    Genus  = sapply(desc_parts, function(x) ifelse(length(x) >= 6 && trimws(x[6]) != "", trimws(x[6]), NA_character_))
  ) %>%
  mutate(
    Domain = ifelse(Domain == "NA" | Domain == "", NA_character_, Domain),
    Phylum = ifelse(Phylum == "NA" | Phylum == "", NA_character_, Phylum),
    Class  = ifelse(Class  == "NA" | Class  == "", NA_character_, Class),
    Order  = ifelse(Order  == "NA" | Order  == "", NA_character_, Order),
    Family = ifelse(Family == "NA" | Family == "", NA_character_, Family),
    Genus  = ifelse(Genus  == "NA" | Genus  == "", NA_character_, Genus)
  ) %>%
  mutate(otu_base = gsub("_count$", "", Variable.Name)) %>%
  select(otu_base, Domain, Phylum, Class, Order, Family, Genus) %>%
  column_to_rownames("otu_base") %>%
  as.matrix()

demo_f <- tbl(con, "DEMO_F") %>% collect() %>% mutate(SEQN = as.character(SEQN)) %>% filter(SEQN %in% colnames(otu_mat_counts))
demo_g <- tbl(con, "DEMO_G") %>% collect() %>% mutate(SEQN = as.character(SEQN)) %>% filter(SEQN %in% colnames(otu_mat_counts))
common_cols  <- intersect(names(demo_f), names(demo_g))
demo_combined <- bind_rows(demo_f %>% select(all_of(common_cols)), demo_g %>% select(all_of(common_cols)))
samples_df_nhanes <- demo_combined %>% rename(sample_id = SEQN) %>% column_to_rownames("sample_id")

physeq_tax_df <- tax_mat_nhanes %>%
  as.data.frame() %>%
  rownames_to_column("otu") %>%
  mutate(across(Domain:Genus, ~str_trim(.))) %>%
  mutate(across(Domain:Genus, function(x) {
    x <- gsub("_\\(.*?\\)", "", x)
    x <- gsub("_", " ", x)
    x <- str_trim(x)
    x
  })) %>%
  mutate(across(c(Domain, Phylum, Class, Order, Family, Genus), ~ {
    x <- gsub("(?i)(Unknown|Unclassified)", "unclassified", .)
    ifelse(is.na(x), "unclassified", x)
  }))
physeq_tax_df$FullPath <- with(physeq_tax_df,
  paste(Domain, Phylum, Class, Order, Family, Genus, otu, sep = "/"))

taxonomy_tree <- data.tree::as.Node(physeq_tax_df, pathName = "FullPath", pathDelimiter = "/")
sanitize_node_names <- function(node) {
  node$name <- gsub("[^a-zA-Z0-9._-]+", "_", node$name)
  if (!node$isLeaf) for (child in node$children) sanitize_node_names(child)
}
sanitize_node_names(taxonomy_tree)
convert_to_newick <- function(node) {
  if (node$isLeaf) return(node$name)
  children_newick <- sapply(node$children, convert_to_newick)
  paste0("(", paste(children_newick, collapse = ","), ")", node$name)
}
newick_str <- paste0(convert_to_newick(taxonomy_tree), ";")
phylo_obj  <- ape::read.tree(text = newick_str)
phylo_obj  <- ape::multi2di(phylo_obj)
if (is.null(phylo_obj$edge.length)) phylo_obj$edge.length <- rep(1, nrow(phylo_obj$edge))
phylo_obj$edge.length[is.na(phylo_obj$edge.length) | is.infinite(phylo_obj$edge.length)] <- 1
common_tips <- intersect(phylo_obj$tip.label, rownames(otu_mat_counts))
phylo_obj   <- ape::keep.tip(phylo_obj, common_tips)

OTU_counts             <- otu_table(otu_mat_counts,             taxa_are_rows = TRUE)
OTU_relative           <- otu_table(otu_mat_relative,           taxa_are_rows = TRUE)
OTU_relative_none      <- otu_table(otu_mat_relative_none,      taxa_are_rows = TRUE)
OTU_relative_clr       <- otu_table(otu_mat_relative_clr,       taxa_are_rows = TRUE)
OTU_relative_lognorm   <- otu_table(otu_mat_relative_lognorm,   taxa_are_rows = TRUE)
OTU_relative_hellinger <- otu_table(otu_mat_relative_hellinger, taxa_are_rows = TRUE)
TAX     <- tax_table(tax_mat_nhanes)
SAMPLES <- sample_data(samples_df_nhanes)

ubiome_counts             <- phyloseq(OTU_counts,             TAX, SAMPLES)
ubiome_relative           <- phyloseq(OTU_relative,           TAX, SAMPLES)
ubiome_relative_none      <- phyloseq(OTU_relative_none,      TAX, SAMPLES)
ubiome_relative_clr       <- phyloseq(OTU_relative_clr,       TAX, SAMPLES)
ubiome_relative_lognorm   <- phyloseq(OTU_relative_lognorm,   TAX, SAMPLES)
ubiome_relative_hellinger <- phyloseq(OTU_relative_hellinger, TAX, SAMPLES)

saveRDS(ubiome_counts,             file.path(intermediate_files_path, "ubiome_counts.rds"))
saveRDS(ubiome_relative,           file.path(intermediate_files_path, "ubiome_relative.rds"))
saveRDS(ubiome_relative_none,      file.path(intermediate_files_path, "ubiome_relative_none.rds"))
saveRDS(ubiome_relative_clr,       file.path(intermediate_files_path, "ubiome_relative_clr.rds"))
saveRDS(ubiome_relative_lognorm,   file.path(intermediate_files_path, "ubiome_relative_lognorm.rds"))
saveRDS(ubiome_relative_hellinger, file.path(intermediate_files_path, "ubiome_relative_hellinger.rds"))

phyloseq_objects <- list(
  counts             = ubiome_counts,
  relative           = ubiome_relative,
  relative_none      = ubiome_relative_none,
  relative_clr       = ubiome_relative_clr,
  relative_lognorm   = ubiome_relative_lognorm,
  relative_hellinger = ubiome_relative_hellinger
)
saveRDS(phyloseq_objects, file.path(intermediate_files_path, "nhanes_phyloseq_objects_all.rds"))

dbDisconnect(con)
message("Done. Outputs under ", intermediate_files_path)
