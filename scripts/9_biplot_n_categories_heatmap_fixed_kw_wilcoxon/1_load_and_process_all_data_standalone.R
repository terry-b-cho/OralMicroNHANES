#!/usr/bin/env Rscript

# =============================================================================
# BIPLOT + CATEGORIES HEATMAP DATA PREPARATION
# =============================================================================
# Builds every intermediate object consumed by the plotting script:
#   - phyloseq objects (with NHANES metadata merged in)
#   - per-WAS sample_data subsets
#   - factorised metadata for the categorical heatmaps
#   - GOLD genus mapping
#
# All outputs are written under:
#   results/analyses_results/9_biplot_n_categories_heatmap_fixed_kw_wilcoxon_out/intermediate/
#
# Environment: R >= 4.5 with phyloseq, DBI, RSQLite, dplyr, tidyr, stringr,
# purrr, tibble, readr, forcats, data.table.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(readr)
  library(forcats)
  library(data.table)
})

utils::globalVariables(c("SEQN"))
SEQN <- NULL

set.seed(42)

message("=== BIPLOT + CATEGORIES HEATMAP DATA PREPARATION (STANDALONE) ===")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
base_path <- PROJECT_ROOT

db_path <- file.path(
  base_path,
  "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"
)
phyloseq_obj_files_path <- file.path(
  base_path,
  "results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate"
)
config_dir_path <- file.path(base_path, "configs")
gold_mapping_path <- file.path(
  base_path,
  "results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv"
)

output_root <- file.path(
  base_path,
  "results/analyses_results/9_biplot_n_categories_heatmap_fixed_kw_wilcoxon_out"
)
data_path <- file.path(output_root, "intermediate")

dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

message("Paths configured:")
message("  Database:            ", db_path)
message("  Phyloseq objects:    ", phyloseq_obj_files_path)
message("  Config directory:    ", config_dir_path)
message("  GOLD mapping:        ", gold_mapping_path)
message("  Intermediate output: ", data_path)

# -----------------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------------
required_files <- c(
  db_path,
  file.path(phyloseq_obj_files_path, "ubiome_relative.rds"),
  file.path(config_dir_path, "1_demoWAS_vars.txt"),
  file.path(config_dir_path, "2_oradWAS_vars.txt"),
  file.path(config_dir_path, "3_exWAS_vars.txt"),
  file.path(config_dir_path, "4_pheWAS_vars.txt"),
  file.path(config_dir_path, "5_outWAS_vars.txt"),
  gold_mapping_path
)

missing_inputs <- required_files[!file.exists(required_files)]
if (length(missing_inputs)) {
  stop("Missing required input(s):\n", paste("  -", missing_inputs, collapse = "\n"))
}

# -----------------------------------------------------------------------------
# Database connection
# -----------------------------------------------------------------------------
message("Connecting to SQLite database...")
con <- dbConnect(SQLite(), dbname = db_path)
on.exit({
  dbDisconnect(con)
  message("Database connection closed.")
}, add = TRUE)

all_tables <- dbListTables(con)
message("  [OK] Found ", length(all_tables), " tables")

# -----------------------------------------------------------------------------
# Load table & variable descriptions (for documentation/parity)
# -----------------------------------------------------------------------------
fg_tables_for_filter <- all_tables[
  grepl("_[FG]$", all_tables, ignore.case = TRUE) |
    all_tables %in% c(
      "d_outcome_mcq_f", "d_outcome_mcq_g",
      "d_outcome_oh_f",  "d_outcome_oh_g",
      "OralDisease_F",   "OralDisease_G",
      "phenotype_vars_f",
      "table_names_epcf", "testosterone_vars_f", "testosterone_vars_g",
      "variable_names_epcf"
    )
]

message("Loading variable metadata...")
table_description <- tbl(con, "table_names_epcf") %>% 
  collect() %>% 
  filter(Data.File.Name %in% fg_tables_for_filter)

variable_description <- tbl(con, "variable_names_epcf") %>% 
  collect() %>% 
  filter(Data.File.Name %in% fg_tables_for_filter)

saveRDS(table_description, file.path(data_path, "table_description.rds"))
saveRDS(variable_description, file.path(data_path, "variable_description.rds"))
message("  [OK] Table descriptions:    ", nrow(table_description))
message("  [OK] Variable descriptions: ", nrow(variable_description))

# -----------------------------------------------------------------------------
# Load phyloseq objects
# -----------------------------------------------------------------------------
message("Loading phyloseq objects...")
phyloseq_files <- c(
  "ubiome_counts"        = "ubiome_counts.rds",
  "ubiome_relative"      = "ubiome_relative.rds",
  "ubiome_relative_none" = "ubiome_relative_none.rds",
  "ubiome_relative_clr"  = "ubiome_relative_clr.rds"
)

phyloseq_objects <- imap(
  phyloseq_files,
  function(fname, key) {
    obj_path <- file.path(phyloseq_obj_files_path, fname)
    ps <- readRDS(obj_path)
    message(
      sprintf(
        "  [OK] %-24s samples: %4d | taxa: %4d",
        key,
        nsamples(ps),
        ntaxa(ps)
      )
    )
    ps
  }
)

# -----------------------------------------------------------------------------
# Helper: load variables from config files
# -----------------------------------------------------------------------------
load_vars <- function(file_path) {
  vars <- readLines(file_path, warn = FALSE)
  vars <- vars[nzchar(vars) & !grepl("^#", vars)]
  vars <- stringr::str_trim(vars)
  unique(vars)
}

var_files <- c(
  demoWAS = "1_demoWAS_vars.txt",
  oradWAS = "2_oradWAS_vars.txt",
  exWAS   = "3_exWAS_vars.txt",
  pheWAS  = "4_pheWAS_vars.txt",
  outWAS  = "5_outWAS_vars.txt"
)

message("Loading analysis variable lists...")
analysis_vars <- map(
  file.path(config_dir_path, var_files),
  load_vars
) %>%
  set_names(names(var_files)) %>%
  map(~ .x[!.x %in% c("", NA_character_)]) %>%
  map(unique)

walk2(
  names(analysis_vars),
  analysis_vars,
  ~ message(sprintf("  [OK] %s variables: %d", .x, length(.y)))
)

all_vars_union <- sort(unique(unlist(analysis_vars)))

# -----------------------------------------------------------------------------
# Retrieve supplemental sample metadata from SQLite
# -----------------------------------------------------------------------------
message("Retrieving NHANES metadata for ", length(all_vars_union), " variables...")

retrieve_sample_data <- function(seqns, variables, tables, con) {
  result <- data.frame(SEQN = seqns, stringsAsFactors = FALSE)
  found_vars <- character(0)

  for (var in variables) {
    var_added <- FALSE

    for (table_name in tables) {
      tryCatch({
        table_cols <- dbListFields(con, table_name)

        if (all(c("SEQN", var) %in% table_cols)) {
          var_data <- tbl(con, table_name) %>%
            dplyr::select(SEQN, !!sym(var)) %>%
            dplyr::filter(SEQN %in% !!seqns) %>%
          collect() %>%
            mutate(SEQN = as.character(SEQN))

          if (nrow(var_data) > 0) {
            overlap <- sum(seqns %in% var_data$SEQN)

            if (overlap > 0) {
              result <- result %>%
                left_join(var_data, by = "SEQN", suffix = c("", ".new"))

      found_vars <- c(found_vars, var)
      var_added <- TRUE

      break
            }
          }
        }
      }, error = function(e) {
        # silently skip tables that error
      })
    }

    if (!var_added) {
      result[[var]] <- NA
    }
  }

  list(
    data = result,
    found = unique(found_vars),
    missing = setdiff(variables, unique(found_vars))
  )
}

current_seqns <- sample_names(phyloseq_objects[["ubiome_relative"]])
nhanes_meta <- retrieve_sample_data(
  seqns = current_seqns,
  variables = all_vars_union,
  tables = fg_tables_for_filter,
  con = con
)

message("  [OK] Retrieved: ", length(nhanes_meta$found), " variables")
if (length(nhanes_meta$missing)) {
  message("  [WARN] Missing variables (kept as NA): ", length(nhanes_meta$missing))
}

sample_data_new <- nhanes_meta$data %>%
  mutate(SEQN = as.character(SEQN)) %>%
  column_to_rownames("SEQN")

# -----------------------------------------------------------------------------
# Merge supplemental metadata into each phyloseq object
# -----------------------------------------------------------------------------
message("Updating phyloseq sample_data with retrieved variables...")
phyloseq_objects <- map(
  phyloseq_objects,
  function(ps) {
    existing <- data.frame(sample_data(ps))
    new_cols <- setdiff(colnames(sample_data_new), colnames(existing))
    merged <- cbind(
      existing,
      sample_data_new[rownames(existing), new_cols, drop = FALSE]
    )
    sample_data(ps) <- sample_data(merged)
    ps
  }
)

saveRDS(phyloseq_objects, file.path(data_path, "phyloseq_objects.rds"))
message("  [OK] Updated phyloseq objects saved")

# Convenience references
ubiome_relative_none <- phyloseq_objects[["ubiome_relative_none"]]

sample_data_df <- data.frame(sample_data(ubiome_relative_none)) %>%
  rownames_to_column("SEQN")

# -----------------------------------------------------------------------------
# Build per-analysis metadata subsets (mirrors Rmd logic)
# -----------------------------------------------------------------------------
nhanes_id_col <- "SEQN"
nhanes_weight_cols <- character(0) # not required for plotting

subset_with_ids <- function(vars) {
  cols <- unique(c(nhanes_id_col, nhanes_weight_cols, vars))
  sample_data_df %>% dplyr::select(any_of(cols))
}

sample_data_demoWAS_vars_subset <- subset_with_ids(analysis_vars$demoWAS)
sample_data_oradWAS_vars_subset <- subset_with_ids(analysis_vars$oradWAS)
sample_data_exWAS_vars_subset   <- subset_with_ids(analysis_vars$exWAS)
sample_data_pheWAS_vars_subset  <- subset_with_ids(analysis_vars$pheWAS)
sample_data_outWAS_vars_subset  <- subset_with_ids(analysis_vars$outWAS)

sample_data_subsets <- list(
  demoWAS = sample_data_demoWAS_vars_subset,
  oradWAS = sample_data_oradWAS_vars_subset,
  exWAS   = sample_data_exWAS_vars_subset,
  pheWAS  = sample_data_pheWAS_vars_subset,
  outWAS  = sample_data_outWAS_vars_subset
)

saveRDS(sample_data_subsets, file.path(data_path, "sample_data_subsets.rds"))
message("  [OK] Sample-data subsets saved")

# -----------------------------------------------------------------------------
# Factorisation logic for heatmaps (verbatim from Rmd)
# -----------------------------------------------------------------------------
message("Constructing factorised metadata for categorical heatmaps...")

# Pre-defined socio-demographic factors
pre_defined_factors <- sample_data_df %>%
  mutate(across(where(is.character), as.factor)) %>%
  dplyr::select(
    SEQN,
    any_of(c(
      "Gender",
      "AgeGroup",
      "EducationLevel",
      "Ethnicity",
      "US_Born",
      "Household_Size_Factor",
      "Marital_Status",
      "Interview_Language"
    ))
  )

demoWAS_vars_subset_factorized_for_heatmapheatmap <- sample_data_demoWAS_vars_subset %>%
  left_join(pre_defined_factors, by = "SEQN") %>%
  mutate(
    Gender = fct_relevel(as.factor(Gender), "Female"),
    Age_group = fct_relevel(as.factor(AgeGroup), "30-39"),
    Education_level = factor(
      EducationLevel,
      levels = c(
        "< 9th Grade",
        "9-11th Grade",
        "High School",
        "College/AA",
        "College Graduate"
      )
    ) %>% fct_relevel("College/AA"),
    Ethnicity = fct_relevel(as.factor(Ethnicity), "White"),
    US_born = fct_relevel(as.factor(US_Born), "US Born"),
    Income_to_poverty_ratio = case_when(
      INDFMPIR < 0.50 ~ "Below 50%",
      INDFMPIR >= 0.50 & INDFMPIR < 1.00 ~ "50–99%",
      INDFMPIR >= 1.00 & INDFMPIR < 1.25 ~ "100–124%",
      INDFMPIR >= 1.25 & INDFMPIR < 1.50 ~ "125–149%",
      INDFMPIR >= 1.50 & INDFMPIR < 1.85 ~ "150–184%",
      INDFMPIR >= 1.85 & INDFMPIR < 2.00 ~ "185–199%",
      INDFMPIR >= 2.00 ~ "200% and Over",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c(
        "Below 50%",
        "50–99%",
        "100–124%",
        "125–149%",
        "150–184%",
        "185–199%",
        "200% and Over"
      )) %>%
      fct_relevel("150–184%"),
    Household_size = factor(
      Household_Size_Factor,
      levels = c("1", "2", "3", "4", "5", "6", "7+")
    ) %>%
      fct_relevel("4"),
    Marital_status = fct_relevel(as.factor(Marital_Status), "Married"),
    Interview_language = fct_relevel(as.factor(Interview_Language), "English")
  ) %>%
  dplyr::select(
    SEQN,
    Gender,
    Age_group,
    Education_level,
    Ethnicity,
    US_born,
    Income_to_poverty_ratio,
    Household_size,
    Marital_status,
    Interview_language
  )

# OUTWAS + ORADWAS helper vectors
outwas_vars_only <- setdiff(names(sample_data_outWAS_vars_subset), "SEQN")
oradwas_vars_only <- setdiff(names(sample_data_oradWAS_vars_subset), "SEQN")

SEQNs_outWAS_control <- sample_data_outWAS_vars_subset %>%
  filter(if_all(all_of(outwas_vars_only), ~ !is.na(.x) & .x == 0)) %>%
  pull(SEQN)

SEQNs_oradWAS_control <- sample_data_oradWAS_vars_subset %>%
  filter(if_all(all_of(oradwas_vars_only), ~ !is.na(.x) & .x == 0)) %>%
  pull(SEQN)

outdWAS_vars_subset_factorized_for_heatmapheatmap <- sample_data_outWAS_vars_subset %>%
  mutate(
    Asthma = case_when(
      ASTHMA == 1 ~ "Asthma",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Asthma")),
    Bronchitis = case_when(
      BRONCHITIS == 1 ~ "Bronchitis",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Bronchitis")),
    Emphysema = case_when(
      EMPHYSEMA == 1 ~ "Emphysema",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Emphysema")),
    Angina = case_when(
      ANGINA == 1 ~ "Angina",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Angina")),
    Heart_failure = case_when(
      HEART_FAILURE == 1 ~ "Heart failure",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Heart failure")),
    Heart_attack = case_when(
      HEART_ATTACK == 1 ~ "Heart attack",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Heart attack")),
    Stroke = case_when(
      STROKE == 1 ~ "Stroke",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Stroke")),
    CHD = case_when(
      CHD == 1 ~ "CHD",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "CHD")),
    CVD = case_when(
      CVD == 1 ~ "CVD",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "CVD")),
    Breast_cancer = case_when(
      CANCER_BREAST == 1 ~ "Breast cancer",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Breast cancer")),
    Colon_cancer = case_when(
      CANCER_COLON == 1 ~ "Colon cancer",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Colon cancer")),
    Lung_cancer = case_when(
      CANCER_LUNG == 1 ~ "Lung cancer",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Lung cancer")),
    Esophageal_cancer = case_when(
      CANCER_ESOPHAGEAL == 1 ~ "Esophageal cancer",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Esophageal cancer")),
    Prostate_cancer = case_when(
      CANCER_PROSTATE == 1 ~ "Prostate cancer",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Prostate cancer")),
    Diabetes = case_when(
      DIABETES == 1 ~ "Diabetes",
      SEQN %in% SEQNs_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Diabetes"))
  ) %>%
  dplyr::select(
    SEQN,
    Asthma,
    Bronchitis,
    Emphysema,
    Angina,
    Heart_failure,
    Heart_attack,
    Stroke,
    Breast_cancer,
    Colon_cancer,
    Lung_cancer,
    Esophageal_cancer,
    Prostate_cancer,
    Diabetes
  )

oradWAS_vars_subset_factorized_for_heatmapheatmap <- sample_data_oradWAS_vars_subset %>%
  mutate(
    Denture = case_when(
      DENTURE_OHAROCDE == 1 ~ "Denture",
      SEQN %in% SEQNs_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Denture")),
    Gum_disease = case_when(
      GUM_DISEASE_OHAROCGP == 1 ~ "Gum disease",
      SEQN %in% SEQNs_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Gum disease")),
    Oral_hygene = case_when(
      ORAL_HYGIENE_OHAROCOH == 1 ~ "Poor oral hygiene",
      SEQN %in% SEQNs_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Poor oral hygiene")),
    Tooth_decay = case_when(
      TOOTH_DECAY_OHAROCDT == 1 ~ "Tooth decay",
      SEQN %in% SEQNs_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Tooth decay"))
  ) %>%
  dplyr::select(SEQN, Denture, Gum_disease, Oral_hygene, Tooth_decay)

# --- EXWAS factors (viral serology + smoking behaviour) ----------------------
exWAS_vars_subset_factorized_for_heatmap <- sample_data_exWAS_vars_subset %>%
  left_join(
    sample_data_demoWAS_vars_subset %>%
      dplyr::select(SEQN, RIDAGEYR, AGE_SQUARED, RIAGENDR),
    by = "SEQN"
  ) %>%
    mutate(
    Hepatitis_B_surface_antigen = case_when(
      LBDHBG == 1 ~ "Positive",
      LBDHBG == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    Hepatitis_C_antibody = case_when(
      LBDHCV == 1 ~ "Positive",
      LBDHCV == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    Hepatitis_D_anti_HDV = case_when(
      LBDHD == 1 ~ "Positive",
      LBDHD == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    Hepatitis_E_IgG = case_when(
      LBDHEG == 1 ~ "Positive",
      LBDHEG == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    Hepatitis_E_IgM = case_when(
      LBDHEM == 1 ~ "Positive",
      LBDHEM == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    HPV_PCR_summary = case_when(
      LBDRPCR == 1 ~ "Positive",
      LBDRPCR == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    Varicella_antibody = case_when(
      LBDVWCGP == 1 ~ "Positive",
      LBDVWCGP == 2 ~ "Negative",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Negative", "Positive")),
    Days_smoked_in_5_days = case_when(
      SMQ_current_ever_never == 0 ~ "Never smoked",
      SMQ710 >= 1 & SMQ710 <= 5 ~ as.character(SMQ710),
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("Never smoked", "1", "2", "3", "4", "5")) %>%
      fct_relevel("Never smoked"),
    Cigarettes_per_day = case_when(
      SMQ_current_ever_never == 0 ~ "Never smoked",
      SMQ720 >= 1 & SMQ720 <= 3 ~ "1–3 per day",
      SMQ720 > 3 & SMQ720 <= 6 ~ "4–6 per day",
      SMQ720 > 6 & SMQ720 <= 10 ~ "7–10 per day",
      SMQ720 > 10 & SMQ720 <= 20 ~ "11–20 per day",
      SMQ720 > 20 & SMQ720 <= 95 ~ "20+ per day",
      SMQ720 == 999 | is.na(SMQ720) ~ NA_character_,
      TRUE ~ NA_character_
    ) %>%
      factor(
        levels = c(
          "Never smoked",
          "1–3 per day",
          "4–6 per day",
          "7–10 per day",
          "11–20 per day",
          "20+ per day"
        )
      ),
    Last_smoked_cigarette = case_when(
      SMQ725 == 1 ~ "Today",
      SMQ725 == 2 ~ "Yesterday",
      SMQ725 == 3 ~ "3–5 days ago",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("Today", "Yesterday", "3–5 days ago")),
    Smoking_status = case_when(
      SMQ_current_ever_never == 0 ~ "Never smoker",
      SMQ_current_ever_never == 2 ~ "Former smoker",
      SMQ_current_ever_never == 1 ~ "Current smoker",
      TRUE ~ NA_character_
    ) %>%
      factor(levels = c("Never smoker", "Former smoker", "Current smoker")) %>%
      fct_relevel("Never smoker")
  )

# --- PHEWAS factors (anthropometry & vitals) ---------------------------------
pheWAS_vars_subset_factorized_for_heatmap <- sample_data_pheWAS_vars_subset %>%
  left_join(
    sample_data_demoWAS_vars_subset %>%
      dplyr::select(SEQN, RIDAGEYR, AGE_SQUARED, RIAGENDR),
    by = "SEQN"
  ) %>%
  mutate(
    BMI_category = case_when(
      BMXBMI < 18.5 ~ "Underweight",
      BMXBMI >= 18.5 & BMXBMI < 25 ~ "Healthy weight",
      BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight",
      BMXBMI >= 30 & BMXBMI < 35 ~ "Class 1 Obesity",
      BMXBMI >= 35 ~ "Class 2-3 Obesity",
      TRUE ~ NA_character_
    ) %>%
      factor(
        levels = c(
          "Underweight",
          "Healthy weight",
          "Overweight",
          "Class 1 Obesity",
          "Class 2-3 Obesity"
        )
      ) %>%
      fct_relevel("Healthy weight%"),
    pulse_category = case_when(
      BPXPLS < 60 ~ "<60 bpm",
      BPXPLS >= 60 & BPXPLS < 70 ~ "60–70 bpm",
      BPXPLS >= 70 & BPXPLS < 75 ~ "70–75 bpm",
      BPXPLS >= 75 & BPXPLS < 85 ~ "75–85 bpm",
      BPXPLS >= 85 ~ "85+ bpm",
      TRUE ~ NA_character_
    ) %>%
      factor(
        levels = c(
          "<60 bpm",
          "60–70 bpm",
          "70–75 bpm",
          "75–85 bpm",
          "85+ bpm"
        )
      ) %>%
      fct_relevel("70–75 bpm"),
    Blood_pressure = case_when(
      MSYSTOLIC > 180 | MDIASTOLIC > 120 ~ "Hypertensive Crisis",
      MSYSTOLIC >= 140 | MDIASTOLIC >= 90 ~ "Hypertension Stage 2",
      (MSYSTOLIC >= 130 & MSYSTOLIC < 140) |
        (MDIASTOLIC >= 80 & MDIASTOLIC < 90) ~ "Hypertension Stage 1",
      (MSYSTOLIC >= 120 & MSYSTOLIC < 130) & (MDIASTOLIC < 80) ~ "Elevated",
      MSYSTOLIC < 120 & MDIASTOLIC < 80 ~ "Normal",
      TRUE ~ NA_character_
    ) %>%
      factor(
        levels = c(
          "Normal",
          "Elevated",
          "Hypertension Stage 1",
          "Hypertension Stage 2",
          "Hypertensive Crisis"
        )
      ) %>%
      fct_relevel("Normal"),
    Waist_circumference_quintiles = case_when(
      (RIAGENDR == 1 & BMXWAIST >= 58 & BMXWAIST <= 82) |
        (RIAGENDR == 2 & BMXWAIST >= 58 & BMXWAIST <= 80) ~ "Male 58-82cm, Female 58-80cm",
      (RIAGENDR == 1 & BMXWAIST > 82 & BMXWAIST <= 92) |
        (RIAGENDR == 2 & BMXWAIST > 80 & BMXWAIST <= 89) ~ "Male 83-92cm, Female 81-89cm",
      (RIAGENDR == 1 & BMXWAIST > 92 & BMXWAIST <= 100) |
        (RIAGENDR == 2 & BMXWAIST > 89 & BMXWAIST <= 97) ~ "Male 93-100cm, Female 90-97cm",
      (RIAGENDR == 1 & BMXWAIST > 100 & BMXWAIST <= 110) |
        (RIAGENDR == 2 & BMXWAIST > 97 & BMXWAIST <= 109) ~ "Male 101-110cm, Female 98-109cm",
      (RIAGENDR == 1 & BMXWAIST > 110 & BMXWAIST <= 179) |
        (RIAGENDR == 2 & BMXWAIST > 109 & BMXWAIST <= 165) ~ "Male 111-179cm, Female 110-165cm",
      TRUE ~ NA_character_
    ) %>%
      factor(
        levels = c(
          "Male 58-82cm, Female 58-80cm",
          "Male 83-92cm, Female 81-89cm",
          "Male 93-100cm, Female 90-97cm",
          "Male 101-110cm, Female 98-109cm",
          "Male 111-179cm, Female 110-165cm"
        )
      )
  )

# -----------------------------------------------------------------------------
# Bundle & persist heatmap factorised datasets
# -----------------------------------------------------------------------------
heatmap_factorized_datasets <- list(
  demoWAS = demoWAS_vars_subset_factorized_for_heatmapheatmap,
  oradWAS = oradWAS_vars_subset_factorized_for_heatmapheatmap,
  exWAS   = exWAS_vars_subset_factorized_for_heatmap,
  pheWAS  = pheWAS_vars_subset_factorized_for_heatmap,
  outWAS  = outdWAS_vars_subset_factorized_for_heatmapheatmap
)

saveRDS(
  heatmap_factorized_datasets,
  file.path(data_path, "heatmap_factorized_datasets.rds")
)
message("  [OK] Heatmap factorised datasets saved")

# -----------------------------------------------------------------------------
# Save additional reference data
# -----------------------------------------------------------------------------
ubiome_genus_mapping_complete <- fread(gold_mapping_path)
saveRDS(
  ubiome_genus_mapping_complete,
  file.path(data_path, "ubiome_genus_mapping_complete.rds")
)
message("  [OK] GOLD genus mapping saved")

# -----------------------------------------------------------------------------
# Final summary
# -----------------------------------------------------------------------------
message("=== DATA PREPARATION COMPLETE ===")
message("Intermediates written to: ", data_path)

