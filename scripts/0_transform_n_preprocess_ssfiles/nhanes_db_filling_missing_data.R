#!/usr/bin/env Rscript
################################################################################
# nhanes_db_filling_missing_data.R
#
# Step 2: complete the transformed SQLite database. Derives RIAGENDR_01,
# recreates SQL views for the `_none` tables, and fills missing metadata rows
# in `variable_names_epcf` by mapping each new variable to its source table.
#
# Environment: project conda env at envs/.conda/envs/nhanes-analysis (R >= 4.5).
#   module load conda/miniforge3/24.11.3-0
#   eval "$(conda shell.bash hook)"
#   conda activate envs/.conda/envs/nhanes-analysis
#
# Usage:
#   Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
#     --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
#     --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(fs)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-i", "--in_db"), type = "character",
    default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite",
    help    = "input transformed DB"),
  make_option(c("-o", "--out_db"), type = "character",
    default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite",
    help    = "output complete DB")
)))

dir_create(path_dir(opt$out_db), recurse = TRUE)
file_copy(opt$in_db, opt$out_db, overwrite = TRUE)

con <- dbConnect(SQLite(), opt$out_db)

# Recreate SQL views for "_none" tables in the complete DB so downstream
# pipelines querying the complete DB can resolve them.
rel_tables <- c(
  "DADA2RSV_GENUS_RELATIVE_F",  "DADA2RSV_GENUS_RELATIVE_G",
  "DADA2RSV_FAMILY_RELATIVE_F", "DADA2RSV_FAMILY_RELATIVE_G",
  "DADA2RSV_ORDER_RELATIVE_F",  "DADA2RSV_ORDER_RELATIVE_G",
  "DADA2RSV_CLASS_RELATIVE_F",  "DADA2RSV_CLASS_RELATIVE_G",
  "DADA2RSV_PHYLUM_RELATIVE_F", "DADA2RSV_PHYLUM_RELATIVE_G"
)
for (tbl in rel_tables) {
  if (tbl %in% dbListTables(con)) {
    view_name <- paste0(tbl, "_none")
    dbExecute(con, sprintf("CREATE VIEW IF NOT EXISTS %s AS SELECT * FROM %s",
                           view_name, tbl))
  }
}

# Derive RIAGENDR_01 in-place: 0 = Male, 1 = Female.
for (series in c("F", "G")) {
  tbln <- paste0("DEMO_", series)
  if (!tbln %in% dbListTables(con)) next
  cols <- dbListFields(con, tbln)
  if (!"RIAGENDR_01" %in% cols) {
    dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN RIAGENDR_01 INTEGER", tbln))
    dbExecute(con, sprintf("UPDATE %s SET RIAGENDR_01 = RIAGENDR - 1", tbln))
  }
}

# Manual descriptions for derived variables.
annotations <- tribble(
  ~var_name,              ~var_description,
  "AGE_SQUARED",          "Age squared (RIDAGEYR^2)",
  "BORN_INUSA",           "Born in USA (DMDBORN4=1)",
  "EDUCATION_LESS9",      "Edu: <9th grade (DMDEDUC2=1)",
  "EDUCATION_9_11",       "Edu: 9-11th grade (DMDEDUC2=2)",
  "EDUCATION_AA",         "Edu: AA/some college (DMDEDUC2=4)",
  "EDUCATION_COLLEGEGRAD","Edu: college grad./above (DMDEDUC2=5)",
  "ETHNICITY_MEXICAN",    "Mexican American (RIDRETH1=1)",
  "ETHNICITY_OTHERHISPANIC","Other Hispanic (RIDRETH1=2)",
  "ETHNICITY_OTHER",      "Other race/multi-racial (RIDRETH1=5)",
  "ETHNICITY_NONHISPANICBLACK","Non-Hispanic Black (RIDRETH1=4)",
  "RIAGENDR_01",          "Gender (0=Male,1=Female)",
  "LBXP1",           "Total prostate specific antigen (ng/mL)",
  "LBXP2",           "Free prostate specific antigen (ng/mL)",
  "LBDP3",           "Prostate specific antigen ratio (%)",
  "LBDRFO",          "RBC folate (ng/mL)",
  "LBXSF6SI",        "Mefox oxidation product (nmol/L)",
  "LBXSCK",          "Creatine kinase (U/L)",
  "LBXFOLSI",        "Serum folate (nmol/L)",
  "LBXPS4",          "Complex prostate specific antigen (ng/mL)",
  "LBXTBM",          "TB Mitogen control result (Positive)",
  "BMXSUB",          "Subscapular Skinfold (mm)",
  "LBDFOT",          "Serum total folate (ng/mL)",
  "LBXMMASI",        "Mono-methyl arsenic (ng/mL)",
  "BMXTRI",          "Triceps Skinfold (mm)",
  "URX1DC",          "NAC-(1,2-dichlorovinyl)-L-cys (ng/mL)",
  "ASTHMA",          "Asthma (self-reported)",
  "BRONCHITIS",      "Bronchitis (self-reported)",
  "EMPHYSEMA",       "Emphysema (self-reported)",
  "ANGINA",          "Angina (self-reported)",
  "HEART_FAILURE",   "Heart failure (self-reported)",
  "HEART_ATTACK",    "Heart attack (self-reported)",
  "STROKE",          "Stroke (self-reported)",
  "CHD",             "Coronary heart disease (self-reported)",
  "CVD",             "Cardiovascular disease (self-reported)",
  "CANCER_BREAST",   "Breast cancer history",
  "CANCER_COLON",    "Colon cancer history",
  "CANCER_LUNG",     "Lung cancer history",
  "CANCER_ESOPHAGEAL","Esophageal cancer history",
  "CANCER_PROSTATE", "Prostate cancer history",
  "CANCER_MOUTH",    "Oral cancer history",
  "DIABETES",        "Diabetes (self-reported)"
)

# Regression dummy variable annotations -- only those that are actually
# created by universal_was_analysis.R from multi-level categorical sources.
categorical_vars <- DBI::dbGetQuery(con, "
  SELECT \"Variable.Name\", GROUP_CONCAT(\"values\") as all_levels
  FROM e_variable_levels
  WHERE \"values\" IS NOT NULL
  GROUP BY \"Variable.Name\"
  HAVING COUNT(DISTINCT \"values\") > 1
") %>%
  filter(grepl("^[0-9.,\\s]+$", all_levels)) %>%
  mutate(
    levels_list = strsplit(all_levels, ","),
    has_only_integers = map_lgl(levels_list, function(x) {
      nums <- as.numeric(trimws(x))
      all(!is.na(nums)) && all(nums == round(nums))
    })
  ) %>%
  filter(has_only_integers) %>%
  mutate(
    sorted_levels = map(levels_list, function(x) sort(as.numeric(trimws(x)))),
    dummy_vars = map2(`Variable.Name`, sorted_levels, function(var_name, lvls) {
      if (length(lvls) > 1) paste0(var_name, lvls[-1]) else character(0)
    })
  )

actual_dummy_vars <- categorical_vars %>%
  select(`Variable.Name`, dummy_vars) %>%
  unnest(dummy_vars) %>%
  pull(dummy_vars)

dummy_annotations <- tribble(
  ~var_name,              ~var_description,
  "RIAGENDR2",            "Gender: Female (RIAGENDR=2, ref: Male=1)",
  "DS2DS2",               "Any dietary supplements taken: No (DS2DS=2, ref: Yes=1)",
  "DS2AN2",               "Any antacids taken: No (DS2AN=2, ref: Yes=1)",
  "LBXHBC2",              "Hepatitis B core antibody: Negative (LBXHBC=2, ref: Positive=1)",
  "LBXHBC3",              "Hepatitis B core antibody: Indeterminate (LBXHBC=3, ref: Positive=1)",
  "LBDHBG2",              "Hepatitis B surface antigen: Negative (LBDHBG=2, ref: Positive=1)",
  "LBDHCV2",              "Hepatitis C antibody: Negative (LBDHCV=2, ref: Positive=1)",
  "LBDHCV5",              "Hepatitis C antibody: Indeterminate (LBDHCV=5, ref: Positive=1)",
  "LBXHA2",               "Hepatitis A antibody: Negative (LBXHA=2, ref: Positive=1)",
  "LBXHA3",               "Hepatitis A antibody: Indeterminate (LBXHA=3, ref: Positive=1)",
  "LBXHBS2",              "Hepatitis B surface antibody: Negative (LBXHBS=2, ref: Positive=1)",
  "LBXHBS3",              "Hepatitis B surface antibody: Indeterminate (LBXHBS=3, ref: Positive=1)",
  "LBXHE12",              "Herpes Simplex Virus I: Negative (LBXHE1=2, ref: Positive=1)",
  "LBXHE13",              "Herpes Simplex Virus I: Indeterminate (LBXHE1=3, ref: Positive=1)",
  "LBXHE22",              "Herpes Simplex Virus II: Negative (LBXHE2=2, ref: Positive=1)",
  "LBXHE23",              "Herpes Simplex Virus II: Indeterminate (LBXHE2=3, ref: Positive=1)",
  "LBXTTG2",              "Tissue transglutaminase IgA-TTG: Negative (LBXTTG=2, ref: Positive=1)",
  "LBXTTG3",              "Tissue transglutaminase IgA-TTG: Weakly positive (LBXTTG=3, ref: Positive=1)",
  "URXUCL2",              "Chlamydia (Urine): Negative (URXUCL=2, ref: Positive=1)",
  "URXUCL3",              "Chlamydia (Urine): Indeterminate (URXUCL=3, ref: Positive=1)",
  "LBX062",               "HPV type 06: Negative (LBX06=2, ref: Positive=1)",
  "LBX112",               "HPV type 11: Negative (LBX11=2, ref: Positive=1)",
  "LBX162",               "HPV type 16: Negative (LBX16=2, ref: Positive=1)",
  "LBX182",               "HPV type 18: Negative (LBX18=2, ref: Positive=1)",
  "ORXHPV2",              "Oral HPV result: Negative (ORXHPV=2, ref: Positive=1)",
  "ORXHPV3",              "Oral HPV result: Not evaluated (ORXHPV=3, ref: Positive=1)",
  "HOQ0652",              "Home ownership: Rented (HOQ065=2, ref: Owned/buying=1)",
  "HOQ0653",              "Home ownership: Other arrangement (HOQ065=3, ref: Owned/buying=1)",
  "HOQ0657",              "Home ownership: Refused to answer (HOQ065=7, ref: Owned/buying=1)",
  "HOQ0659",              "Home ownership: Don't know (HOQ065=9, ref: Owned/buying=1)",
  "LBDHD2",               "Hepatitis D antibody: Negative (LBDHD=2, ref: Positive=1)",
  "LBDHEG2",              "Hepatitis E IgG (anti-HEV): Negative (LBDHEG=2, ref: Positive=1)",
  "LBDHEM2",              "Hepatitis E IgM (anti-HEV): Negative (LBDHEM=2, ref: Positive=1)",
  "LBDRPCR2",             "Roche HPV linear array (LA): Negative (LBDRPCR=2, ref: Positive=1)",
  "LBDRPCR3",             "Roche HPV linear array (LA): Inadequate (LBDRPCR=3, ref: Positive=1)",
  "LBDVWCGP2",            "von Willebrand factor: Elevated (LBDVWCGP=2, ref: Normal=1)",
  "LBXML13",              "Levofloxacin 1 MRSA: Level 3 (LBXML1=3, ref: Level 1=1)",
  "URXUAS3",              "Total arsenic - urine: Level 3 (URXUAS=3, ref: Level 1=1)",
  "URXUAS5",              "Total arsenic - urine: Level 5 (URXUAS=5, ref: Level 1=1)",
  "LBXEPAH",              "2-(N-Ethyl-perfluorooctane sulfonamido) acetic acid: High level (LBXEPAH=high, ref: low)",
  "DS2ANCNT",             "Total number of antacids taken (0-2)",
  "DS2DSCNT",             "Total number of dietary supplements taken (0-18)",
  "SMQ7102",              "Days smoked past 5 days: Level 2 (SMQ710=2, ref: Level 1=1)",
  "SMQ7103",              "Days smoked past 5 days: Level 3 (SMQ710=3, ref: Level 1=1)",
  "SMQ7104",              "Days smoked past 5 days: Level 4 (SMQ710=4, ref: Level 1=1)",
  "SMQ7105",              "Days smoked past 5 days: Level 5 (SMQ710=5, ref: Level 1=1)",
  "SMQ7252",              "Last time smoked: Yesterday (SMQ725=2, ref: Today=1)",
  "SMQ7253",              "Last time smoked: 3-5 days ago (SMQ725=3, ref: Today=1)",
  "SMQ_current_ever_never1", "Smoking status: Former smoker (SMQ_current_ever_never=1, ref: Never=0)",
  "SMQ_current_ever_never2", "Smoking status: Current smoker (SMQ_current_ever_never=2, ref: Never=0)",
  "SMD4152",              "Total smokers in home: Two (SMD415=2, ref: One=1)",
  "SMD4153",              "Total smokers in home: Three or more (SMD415=3, ref: One=1)",
  "SMD415A2",             "Cigarette smokers in home: Two (SMD415A=2, ref: One=1)",
  "SMD415A3",             "Cigarette smokers in home: Three or more (SMD415A=3, ref: One=1)"
) %>%
  filter(var_name %in% actual_dummy_vars)

quartile_annotations <- tribble(
  ~var_name, ~var_description
)

all_annotations <- bind_rows(annotations, dummy_annotations, quartile_annotations)

template <- tbl(con, "variable_names_epcf") %>%
  select(Data.File.Name, Data.File.Description, Begin.Year, EndYear,
         Component, Use.Constraints) %>%
  distinct() %>%
  collect()

all_vars   <- all_annotations$var_name
all_tables <- dbListTables(con)
fields_map <- setNames(
  lapply(all_tables, function(tbl) dbListFields(con, tbl)),
  all_tables
)

mapping_long <- tibble(Variable = all_vars) %>%
  mutate(Data.File.Name = map(Variable, function(v) {
    names(fields_map)[vapply(fields_map, function(flds) v %in% flds, logical(1))]
  })) %>%
  unnest_longer(Data.File.Name) %>%
  filter(!is.na(Data.File.Name)) %>%
  distinct(Variable, Data.File.Name)

# Map dummy variables (which do not exist as physical columns) to their
# source categorical variable's table.
source_var_mapping <- categorical_vars %>%
  select(source_var = `Variable.Name`, dummy_vars) %>%
  unnest(dummy_vars) %>%
  select(dummy_var = dummy_vars, source_var)

unmapped_dummies <- setdiff(source_var_mapping$dummy_var, mapping_long$Variable)
if (length(unmapped_dummies) > 0) {
  dummy_mapping_expanded <- source_var_mapping %>%
    filter(dummy_var %in% unmapped_dummies) %>%
    mutate(Data.File.Name = map(source_var, function(src) {
      names(fields_map)[vapply(fields_map, function(flds) src %in% flds, logical(1))]
    })) %>%
    unnest_longer(Data.File.Name) %>%
    filter(!is.na(Data.File.Name)) %>%
    select(Variable = dummy_var, Data.File.Name) %>%
    distinct()
  mapping_long <- bind_rows(mapping_long, dummy_mapping_expanded)
}

existing_meta <- tbl(con, "variable_names_epcf") %>%
  select(Variable.Name, Data.File.Name) %>%
  collect()

to_add <- anti_join(mapping_long, existing_meta,
                    by = c("Variable" = "Variable.Name", "Data.File.Name"))

if (nrow(to_add) > 0) {
  new_meta <- to_add %>%
    left_join(all_annotations, by = c("Variable" = "var_name")) %>%
    rename(Variable.Description = var_description) %>%
    left_join(template, by = "Data.File.Name") %>%
    select(
      Variable.Name = Variable,
      Variable.Description,
      Use.Constraints,
      Data.File.Name,
      Data.File.Description,
      Begin.Year,
      EndYear,
      Component
    )

  dbBegin(con)
  tryCatch({
    old <- tbl(con, "variable_names_epcf") %>% collect()
    dbWriteTable(con, "variable_names_epcf",
                 bind_rows(old, new_meta), overwrite = TRUE)
    dbCommit(con)
  }, error = function(e) {
    dbRollback(con)
    stop("Error appending metadata: ", e$message)
  })
}

dbDisconnect(con)
