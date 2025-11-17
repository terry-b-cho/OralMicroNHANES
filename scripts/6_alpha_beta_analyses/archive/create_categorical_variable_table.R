#!/usr/bin/env Rscript
################################################################################
##  CREATE COMPREHENSIVE CATEGORICAL VARIABLE DEFINITIONS TABLE              ##
##  Generates supplementary_table_categorical_variable_definitions.csv/md    ##
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(knitr)
  library(readr)
})

cat("=== CREATING CATEGORICAL VARIABLE DEFINITIONS TABLE ===\n\n")

# Load processed data
base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
output_path <- file.path(base_path, "scripts/6_alpha_beta_analyses")

# Load categorical data
categorical_data <- readRDS(file.path(data_path, "all_categorical_data.rds"))

# Define all categorical variable definitions
variable_definitions <- tribble(
  ~categorical_variable_name, ~definition_in_R_code, ~definition_description_in_text, ~reference_level, ~table_source,
  
  # DEMO-WAS Variables
  "Gender", 
  'if("Gender" %in% colnames(.)) fct_relevel(as.factor(Gender), "Female")', 
  "Biological sex of participant from NHANES demographics (RIAGENDR). Male=1, Female=2", 
  "Female",
  "NHANES Demographics (DEMO_F, DEMO_G)",
  
  "Age_group",
  'if("AgeGroup" %in% colnames(.)) fct_relevel(as.factor(AgeGroup), "30-39")',
  "Age categorized into 10-year groups: 18-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+ years from RIDAGEYR",
  "30-39",
  "NHANES Demographics (DEMO_F, DEMO_G) - derived from RIDAGEYR",
  
  "Education_level",
  'factor(EducationLevel, levels = c("< 9th Grade", "9-11th Grade", "High School", "College/AA", "College Graduate")) %>% fct_relevel("College/AA")',
  "Highest education level completed from DMDEDUC2. Categories: <9th Grade, 9-11th Grade, High School, College/AA (Some college or Associate\'s degree), College Graduate (Bachelor\'s degree or higher)",
  "College/AA",
  "NHANES Demographics (DEMO_F, DEMO_G) - DMDEDUC2",
  
  "Ethnicity",
  'if("Ethnicity" %in% colnames(.)) fct_relevel(as.factor(Ethnicity), "White")',
  "Race/ethnicity from RIDRETH3. Categories: White (Non-Hispanic White), Black (Non-Hispanic Black), Hispanic (Mexican American, Other Hispanic), Asian (Non-Hispanic Asian), Other (Other Race including Multi-racial)",
  "White",
  "NHANES Demographics (DEMO_F, DEMO_G) - RIDRETH3",
  
  "US_born",
  'if("US_Born" %in% colnames(.)) fct_relevel(as.factor(US_Born), "US Born")',
  "Country of birth from DMDBORN4. US Born vs Non-US Born (1=Born in US, 2=Not born in US)",
  "US Born",
  "NHANES Demographics (DEMO_F, DEMO_G) - DMDBORN4",
  
  "Household_size",
  'factor(Household_Size_Factor, levels = c("1","2","3","4","5","6","7+")) %>% fct_relevel("4")',
  "Number of people in household from DMDHHSIZ. Categories: 1, 2, 3, 4, 5, 6, 7+ members",
  "4",
  "NHANES Demographics (DEMO_F, DEMO_G) - DMDHHSIZ",
  
  "Marital_status",
  'if("Marital_Status" %in% colnames(.)) fct_relevel(as.factor(Marital_Status), "Married")',
  "Marital status from DMDMARTL. Categories: Married, Widowed, Divorced, Separated, Never married, Living with partner",
  "Married",
  "NHANES Demographics (DEMO_F, DEMO_G) - DMDMARTL",
  
  "Interview_language",
  'if("Interview_Language" %in% colnames(.)) fct_relevel(as.factor(Interview_Language), "English")',
  "Language of interview from SIALANG. English vs Spanish (1=English, 2=Spanish)",
  "English",
  "NHANES Demographics (DEMO_F, DEMO_G) - SIALANG",
  
  "Income_to_poverty_ratio",
  'case_when(INDFMPIR < 0.50 ~ "Below 50%", INDFMPIR >= 0.50 & INDFMPIR < 1.00 ~ "50-99%", INDFMPIR >= 1.00 & INDFMPIR < 1.25 ~ "100-124%", INDFMPIR >= 1.25 & INDFMPIR < 1.50 ~ "125-149%", INDFMPIR >= 1.50 & INDFMPIR < 1.85 ~ "150-184%", INDFMPIR >= 1.85 & INDFMPIR < 2.00 ~ "185-199%", INDFMPIR >= 2.00 ~ "200% and Over", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("150-184%")',
  "Family income to poverty threshold ratio from INDFMPIR. Categorized into 7 levels based on percentage of federal poverty level: Below 50%, 50-99%, 100-124%, 125-149%, 150-184%, 185-199%, 200% and Over. Reference level chosen as median category (150-184%)",
  "150-184%",
  "NHANES Demographics (DEMO_F, DEMO_G) - INDFMPIR",
  
  # ORAL-WAS Variables
  "Denture",
  'case_when(DENTURE_OHAROCDE == 1 ~ "Denture", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Denture"))',
  "Presence of dentures from OHAROCDE. Binary: Denture (OHAROCDE=1) vs control (all oral health indicators=0)",
  "control",
  "NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCDE",
  
  "Gum_disease",
  'case_when(GUM_DISEASE_OHAROCGP == 1 ~ "Gum disease", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Gum disease"))',
  "Self-reported gum disease from OHAROCGP. Binary: Gum disease (OHAROCGP=1) vs control (all oral health indicators=0)",
  "control",
  "NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCGP",
  
  "Oral_hygiene",
  'case_when(ORAL_HYGIENE_OHAROCOH == 1 ~ "Poor oral hygiene", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Poor oral hygiene"))',
  "Poor oral hygiene status from OHAROCOH. Binary: Poor oral hygiene (OHAROCOH=1) vs control (all oral health indicators=0)",
  "control",
  "NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCOH",
  
  "Tooth_decay",
  'case_when(TOOTH_DECAY_OHAROCDT == 1 ~ "Tooth decay", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Tooth decay"))',
  "Presence of tooth decay from OHAROCDT. Binary: Tooth decay (OHAROCDT=1) vs control (all oral health indicators=0)",
  "control",
  "NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCDT",
  
  # EXPOSURE-WAS Variables
  "Smoking_status",
  'case_when(SMQ_current_ever_never == 0 ~ "Never smoker", SMQ_current_ever_never == 2 ~ "Former smoker", SMQ_current_ever_never == 1 ~ "Current smoker", TRUE ~ NA_character_) %>% factor(levels = c("Never smoker", "Former smoker", "Current smoker")) %>% fct_relevel("Never smoker")',
  "Smoking status derived from SMQ020 (ever smoked 100 cigarettes) and SMQ040 (current smoking). Categories: Never smoker (SMQ020=2), Former smoker (SMQ020=1 & SMQ040 in 1,2), Current smoker (SMQ020=1 & SMQ040=3)",
  "Never smoker",
  "NHANES Smoking (SMQ_F, SMQ_G) - derived from SMQ020, SMQ040",
  
  "Hepatitis_C_antibody",
  'case_when(LBDHCV == 1 ~ "Positive", LBDHCV == 2 ~ "Negative", TRUE ~ NA_character_) %>% factor(levels = c("Negative", "Positive"))',
  "Hepatitis C antibody test result from LBDHCV. Binary: Positive (LBDHCV=1) vs Negative (LBDHCV=2). Indeterminate results (LBDHCV=5) treated as NA",
  "Negative",
  "NHANES Laboratory (HEPC_F, HEPC_G) - LBDHCV",
  
  "HPV_PCR_summary",
  'case_when(LBDRPCR == 1 ~ "Positive", LBDRPCR == 2 ~ "Negative", TRUE ~ NA_character_) %>% factor(levels = c("Negative", "Positive"))',
  "Human Papillomavirus (HPV) PCR summary result from LBDRPCR. Binary: Positive (LBDRPCR=1) vs Negative (LBDRPCR=2). Inadequate samples (LBDRPCR=3) treated as NA",
  "Negative",
  "NHANES Laboratory (HPV_F, HPV_G) - LBDRPCR",
  
  # PHEWAS Variables
  "BMI_category",
  'case_when(BMXBMI < 18.5 ~ "Underweight", BMXBMI >= 18.5 & BMXBMI < 25 ~ "Healthy weight", BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight", BMXBMI >= 30 & BMXBMI < 35 ~ "Class 1 Obesity", BMXBMI >= 35 ~ "Class 2-3 Obesity", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("Healthy weight")',
  "Body Mass Index category from BMXBMI using CDC criteria. Categories: Underweight (<18.5), Healthy weight (18.5-24.9), Overweight (25.0-29.9), Class 1 Obesity (30.0-34.9), Class 2-3 Obesity (≥35.0 kg/m²)",
  "Healthy weight",
  "NHANES Examination (BMX_F, BMX_G) - BMXBMI",
  
  "Blood_pressure",
  'case_when(MSYSTOLIC > 180 | MDIASTOLIC > 120 ~ "Hypertensive Crisis", MSYSTOLIC >= 140 | MDIASTOLIC >= 90 ~ "Hypertension Stage 2", (MSYSTOLIC >= 130 & MSYSTOLIC < 140) | (MDIASTOLIC >= 80 & MDIASTOLIC < 90) ~ "Hypertension Stage 1", (MSYSTOLIC >= 120 & MSYSTOLIC < 130) & (MDIASTOLIC < 80) ~ "Elevated", MSYSTOLIC < 120 & MDIASTOLIC < 80 ~ "Normal", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("Normal")',
  "Blood pressure category using AHA/ACC 2017 guidelines. MSYSTOLIC and MDIASTOLIC are means of three BP readings (BPXSY1-3, BPXDI1-3). Categories: Normal (<120/<80), Elevated (120-129/<80), Stage 1 (130-139/80-89), Stage 2 (≥140/≥90), Crisis (>180/>120 mmHg)",
  "Normal",
  "NHANES Examination (BPX_F, BPX_G) - derived from BPXSY1-3, BPXDI1-3",
  
  "Pulse_category",
  'case_when(BPXPLS < 60 ~ "<60 bpm", BPXPLS >= 60 & BPXPLS < 70 ~ "60-70 bpm", BPXPLS >= 70 & BPXPLS < 75 ~ "70-75 bpm", BPXPLS >= 75 & BPXPLS < 85 ~ "75-85 bpm", BPXPLS >= 85 ~ "85+ bpm", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("70-75 bpm")',
  "Pulse rate category from BPXPLS (60-second radial pulse). Categories based on quintile distribution: <60 bpm, 60-70 bpm, 70-75 bpm (reference, median), 75-85 bpm, 85+ bpm",
  "70-75 bpm",
  "NHANES Examination (BPX_F, BPX_G) - BPXPLS",
  
  # OUTCOME-WAS Variables
  "Asthma",
  'case_when(ASTHMA == 1 ~ "Asthma", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Asthma"))',
  "Ever diagnosed with asthma from MCQ010. Binary: Asthma (ASTHMA=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable ASTHMA",
  
  "Bronchitis",
  'case_when(BRONCHITIS == 1 ~ "Bronchitis", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Bronchitis"))',
  "Ever diagnosed with chronic bronchitis from MCQ160K. Binary: Bronchitis (BRONCHITIS=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable BRONCHITIS",
  
  "Emphysema",
  'case_when(EMPHYSEMA == 1 ~ "Emphysema", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Emphysema"))',
  "Ever diagnosed with emphysema from MCQ160L. Binary: Emphysema (EMPHYSEMA=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable EMPHYSEMA",
  
  "Angina",
  'case_when(ANGINA == 1 ~ "Angina", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Angina"))',
  "Ever diagnosed with angina pectoris from MCQ160C. Binary: Angina (ANGINA=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable ANGINA",
  
  "Heart_failure",
  'case_when(HEART_FAILURE == 1 ~ "Heart failure", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Heart failure"))',
  "Ever diagnosed with congestive heart failure from MCQ160B. Binary: Heart failure (HEART_FAILURE=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable HEART_FAILURE",
  
  "Heart_attack",
  'case_when(HEART_ATTACK == 1 ~ "Heart attack", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Heart attack"))',
  "Ever diagnosed with heart attack (myocardial infarction) from MCQ160E. Binary: Heart attack (HEART_ATTACK=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable HEART_ATTACK",
  
  "Stroke",
  'case_when(STROKE == 1 ~ "Stroke", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Stroke"))',
  "Ever diagnosed with stroke from MCQ160F. Binary: Stroke (STROKE=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable STROKE",
  
  "Diabetes",
  'case_when(DIABETES == 1 ~ "Diabetes", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Diabetes"))',
  "Ever diagnosed with diabetes from DIQ010. Binary: Diabetes (DIABETES=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (DIQ_F, DIQ_G) - derived variable DIABETES",
  
  "CHD",
  'case_when(CHD == 1 ~ "CHD", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "CHD"))',
  "Ever diagnosed with coronary heart disease from MCQ160D. Binary: CHD (CHD=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CHD",
  
  "CVD",
  'case_when(CVD == 1 ~ "CVD", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "CVD"))',
  "Any cardiovascular disease (composite of heart attack, CHD, angina, stroke, heart failure). Binary: CVD (CVD=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - composite derived variable",
  
  "Breast_cancer",
  'case_when(CANCER_BREAST == 1 ~ "Breast cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Breast cancer"))',
  "Ever diagnosed with breast cancer from MCQ220. Binary: Breast cancer (CANCER_BREAST=1) vs control (all outcome indicators=0). Female participants only",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_BREAST",
  
  "Colon_cancer",
  'case_when(CANCER_COLON == 1 ~ "Colon cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Colon cancer"))',
  "Ever diagnosed with colon cancer from MCQ220. Binary: Colon cancer (CANCER_COLON=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_COLON",
  
  "Lung_cancer",
  'case_when(CANCER_LUNG == 1 ~ "Lung cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Lung cancer"))',
  "Ever diagnosed with lung cancer from MCQ220. Binary: Lung cancer (CANCER_LUNG=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_LUNG",
  
  "Esophageal_cancer",
  'case_when(CANCER_ESOPHAGEAL == 1 ~ "Esophageal cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Esophageal cancer"))',
  "Ever diagnosed with esophageal cancer from MCQ220. Binary: Esophageal cancer (CANCER_ESOPHAGEAL=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_ESOPHAGEAL",
  
  "Prostate_cancer",
  'case_when(CANCER_PROSTATE == 1 ~ "Prostate cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Prostate cancer"))',
  "Ever diagnosed with prostate cancer from MCQ220. Binary: Prostate cancer (CANCER_PROSTATE=1) vs control (all outcome indicators=0). Male participants only",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_PROSTATE",
  
  "Mouth_cancer",
  'case_when(CANCER_MOUTH == 1 ~ "Mouth cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Mouth cancer"))',
  "Ever diagnosed with oral/pharyngeal cancer from MCQ220. Binary: Mouth cancer (CANCER_MOUTH=1) vs control (all outcome indicators=0)",
  "control",
  "NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_MOUTH"
)

cat("Defined", nrow(variable_definitions), "categorical variables\n")

# Add count statistics from actual data
cat("Calculating statistics from data...\n")

final_table <- variable_definitions %>%
  rowwise() %>%
  mutate(
    NA_count = if(categorical_variable_name %in% colnames(categorical_data)) {
      sum(is.na(categorical_data[[categorical_variable_name]]))
    } else NA_integer_,
    
    count_per_category = if(categorical_variable_name %in% colnames(categorical_data)) {
      var_data <- categorical_data[[categorical_variable_name]]
      if(is.factor(var_data)) {
        counts <- table(var_data, useNA = "no")
        paste(names(counts), "=", counts, collapse = "; ")
      } else "Not factor"
    } else "Variable not found"
  ) %>%
  ungroup() %>%
  select(categorical_variable_name, definition_in_R_code, definition_description_in_text, 
         reference_level, NA_count, count_per_category, table_source)

cat("Statistics calculated\n")

# Save as CSV
csv_file <- file.path(output_path, "supplementary_table_categorical_variable_definitions.csv")
write_csv(final_table, csv_file)
cat("✓ Saved CSV:", csv_file, "\n")

# Create markdown table
md_file <- file.path(output_path, "supplementary_table_categorical_variable_definitions.md")

md_content <- c(
  "# Supplementary Table: Categorical Variable Definitions",
  "",
  "## Complete Definitions for All Categorical Variables Used in Diversity Analysis",
  "",
  "This table provides comprehensive documentation of all categorical variables used in the integrated diversity analysis. Each variable includes:",
  "- **Categorical Variable Name**: Variable name as used in analysis scripts",
  "- **Definition in R Code**: Exact R code used to create the categorical variable",
  "- **Definition Description**: Human-readable description of the variable and its categories",
  "- **Reference Level**: The reference (baseline) category used in statistical comparisons",
  "- **NA Count**: Number of participants with missing data for this variable",
  "- **Count Per Category**: Sample size for each category level",
  "- **Table Source**: Original NHANES table(s) and variable(s) used to derive this categorical variable",
  "",
  "---",
  ""
)

# Create markdown table manually for better formatting
md_content <- c(md_content, 
                "| Categorical Variable Name | Definition Description | Reference Level | NA Count | Count Per Category | Table Source |",
                "|---------------------------|------------------------|-----------------|----------|---------------------|--------------|")

for(i in 1:nrow(final_table)) {
  row <- final_table[i,]
  md_content <- c(md_content, sprintf(
    "| %s | %s | %s | %s | %s | %s |",
    row$categorical_variable_name,
    row$definition_description_in_text,
    row$reference_level,
    ifelse(is.na(row$NA_count), "N/A", as.character(row$NA_count)),
    row$count_per_category,
    row$table_source
  ))
}

md_content <- c(md_content, "", "---", "",
                "## R Code Definitions",
                "",
                "For reproducibility, the exact R code used to create each categorical variable is provided below:",
                "")

for(i in 1:nrow(final_table)) {
  row <- final_table[i,]
  md_content <- c(md_content,
                  sprintf("### %s", row$categorical_variable_name),
                  "```r",
                  row$definition_in_R_code,
                  "```",
                  "")
}

md_content <- c(md_content, "",
                "## Notes",
                "",
                "1. **Control Groups**: For disease outcome variables (Oral-WAS, Outcome-WAS), the 'control' group consists of participants with ALL outcome indicators equal to 0 (no missing values allowed). This ensures a true disease-free comparison group.",
                "",
                "2. **Reference Levels**: Reference levels were chosen based on:",
                "   - Most common category (e.g., 'Married' for marital status)",
                "   - Clinical baseline (e.g., 'Healthy weight' for BMI, 'Normal' for blood pressure)",
                "   - Unexposed group (e.g., 'Never smoker' for smoking status, 'control' for disease outcomes)",
                "",
                "3. **Missing Data**: Variables with 'Indeterminate' or insufficient sample results are coded as NA and excluded from analysis.",
                "",
                "4. **Data Source**: All variables derived from NHANES 2009-2010 (cycles F) and 2011-2012 (cycles G).",
                "",
                "---",
                "",
                sprintf("*Table generated on: %s*", Sys.Date()),
                sprintf("*Total variables documented: %d*", nrow(final_table)),
                sprintf("*Total participants: %d*", nrow(categorical_data)))

writeLines(md_content, md_file)
cat("✓ Saved Markdown:", md_file, "\n")

cat("\n=== TABLE GENERATION COMPLETE ===\n")
cat("Files created:\n")
cat("  1.", csv_file, "\n")
cat("  2.", md_file, "\n")
cat("\nSummary:\n")
cat("  Total variables:", nrow(final_table), "\n")
cat("  Variables with data:", sum(!is.na(final_table$NA_count)), "\n")
cat("  Variables missing from dataset:", sum(is.na(final_table$NA_count)), "\n")

