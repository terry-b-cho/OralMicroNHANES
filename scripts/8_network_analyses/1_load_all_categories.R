#!/usr/bin/env Rscript

# =============================================================================
# NETWORK ANALYSIS - DATA LOADING SCRIPT
# =============================================================================
# Loads phyloseq objects, NHANES sample metadata, and aggregated WAS results;
# applies prevalence + significance filters; writes every intermediate object
# consumed by the downstream network scripts.
#
# All outputs are written under:
#   results/analyses_results/8_network_analyses_out/intermediate/
#
# Environment: R >= 4.5 with phyloseq, DBI, RSQLite, dplyr, tidyr, stringr,
# purrr, tibble, readr.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

cat("=== NETWORK ANALYSIS - DATA LOADING ===\n")

# Load required libraries
library(phyloseq)
library(DBI)
library(RSQLite)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(readr)

set.seed(42)

determine_roles <- function(df) {
  dep_matches <- sum(df$dependent_var == df$otu, na.rm = TRUE)
  ind_matches <- sum(df$independent_var == df$otu, na.rm = TRUE)

  if (dep_matches > ind_matches) {
    list(microbe_col = "dependent_var", host_col = "independent_var")
  } else {
    list(microbe_col = "independent_var", host_col = "dependent_var")
  }
}

# =============================================================================
# SETUP PATHS
# =============================================================================

base_path <- PROJECT_ROOT

# Input paths (EXACTLY as original)
db_path <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
phyloseq_obj_files_path <- file.path(base_path, "results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate")
config_dir_path <- file.path(base_path, "configs")
aggregated_association_res_path <- file.path(base_path, "results")

# Output paths
data_path <- file.path(base_path, "results/analyses_results/8_network_analyses_out/intermediate")
dir.create(data_path, recursive = TRUE, showWarnings = FALSE)

cat("Paths configured:\n")
cat("  Database:", db_path, "\n")
cat("  Phyloseq objects:", phyloseq_obj_files_path, "\n")
cat("  Config files:", config_dir_path, "\n")
cat("  Association results:", aggregated_association_res_path, "\n")
cat("  Output data:", data_path, "\n\n")

# =============================================================================
# VERIFY REQUIRED FILES
# =============================================================================

required_files <- c(
  db_path,
  file.path(phyloseq_obj_files_path, "ubiome_relative.rds"),
  file.path(config_dir_path, "1_demoWAS_vars.txt"),
  file.path(config_dir_path, "2_oradWAS_vars.txt"),
  file.path(config_dir_path, "3_exWAS_vars.txt"),
  file.path(config_dir_path, "4_pheWAS_vars.txt"),
  file.path(config_dir_path, "5_outWAS_vars.txt")
)

for (file in required_files) {
  if (!file.exists(file)) {
    stop("Required file missing: ", file)
  }
}

cat("[done] All required files found\n\n")

# =============================================================================
# CONNECT TO DATABASE
# =============================================================================

con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# Get available tables
all_tables <- dbListTables(con)
cat("Database connected. Available tables:", length(all_tables), "\n")

# =============================================================================
# LOAD CONFIG FILES (EXACTLY AS ORIGINAL)
# =============================================================================

cat("Loading variable configuration files...\n")

# Variable configuration files (EXACTLY as original)
var_files <- c(
  "demoWAS" = file.path(config_dir_path, "1_demoWAS_vars.txt"),
  "oradWAS" = file.path(config_dir_path, "2_oradWAS_vars.txt"),
  "exWAS" = file.path(config_dir_path, "3_exWAS_vars.txt"),
  "pheWAS" = file.path(config_dir_path, "4_pheWAS_vars.txt"),
  "outWAS" = file.path(config_dir_path, "5_outWAS_vars.txt")
)

# Load variables from config files (EXACTLY as original)
load_vars <- function(file_path) {
  if (!file.exists(file_path)) return(character(0))
  vars <- readLines(file_path, warn = FALSE)
  vars[nzchar(vars) & !grepl("^#", vars)]  # Remove empty lines and comments
}

all_vars <- purrr::map(var_files, load_vars)
all_vars <- purrr::map(all_vars, unique)
all_vars_union <- unique(unlist(all_vars))
cat("Variables loaded:\n")
for (group in names(all_vars)) {
  cat("  ", group, ":", length(all_vars[[group]]), "variables\n")
}

# =============================================================================
# LOAD DATABASE TABLES (EXACTLY AS ORIGINAL)
# =============================================================================

cat("\nLoading database tables...\n")

# Identify tables for oral microbiome cycles (F and G) - EXACTLY as original
years_F_or_G <- dbListTables(con)[grepl("_F|f", dbListTables(con)) | grepl("_G|g", dbListTables(con))]
years_F_or_G <- c(years_F_or_G, "d_outcome_mcq_f","d_outcome_mcq_g",
                  "d_outcome_oh_f","d_outcome_oh_g",
                  "OralDisease_F", "OralDisease_G",
                  "phenotype_vars_f",
                  "table_names_epcf","testosterone_vars_f","testosterone_vars_g","variable_names_epcf")

# Load table and variable descriptions (EXACTLY as original)
table_description <- tbl(con, "table_names_epcf") %>% 
  collect() %>% 
  filter(Data.File.Name %in% years_F_or_G)

variable_description <- tbl(con, "variable_names_epcf") %>% 
  collect() %>% 
  filter(Data.File.Name %in% years_F_or_G)

cat("  Table descriptions:", nrow(table_description), "tables\n")
cat("  Variable descriptions:", nrow(variable_description), "variables\n")

# =============================================================================
# LOAD PHYLOSEQ OBJECTS (EXACTLY AS ORIGINAL)
# =============================================================================

cat("\nLoading phyloseq objects...\n")

phyloseq_files <- c(
  "ubiome_counts" = "ubiome_counts.rds",
  "ubiome_relative" = "ubiome_relative.rds", 
  "ubiome_relative_none" = "ubiome_relative_none.rds",
  "ubiome_relative_clr" = "ubiome_relative_clr.rds",
  "ubiome_relative_hellinger" = "ubiome_relative_hellinger.rds",
  "ubiome_relative_lognorm" = "ubiome_relative_lognorm.rds"
)

phyloseq_objects <- list()
for (name in names(phyloseq_files)) {
  file_path <- file.path(phyloseq_obj_files_path, phyloseq_files[[name]])
  if (file.exists(file_path)) {
    phyloseq_objects[[name]] <- readRDS(file_path)
    cat("  [done] Loaded", name, "- Samples:", nsamples(phyloseq_objects[[name]]), 
        "Taxa:", ntaxa(phyloseq_objects[[name]]), "\n")
  } else {
    stop("[ERROR] phyloseq object not found:", file_path)
  }
}

# Extract individual objects for convenience (EXACTLY as original)
ubiome_counts <- phyloseq_objects[["ubiome_counts"]]
ubiome_relative <- phyloseq_objects[["ubiome_relative"]]
ubiome_relative_none <- phyloseq_objects[["ubiome_relative_none"]]
ubiome_relative_clr <- phyloseq_objects[["ubiome_relative_clr"]]
ubiome_relative_hellinger <- phyloseq_objects[["ubiome_relative_hellinger"]]
ubiome_relative_lognorm <- phyloseq_objects[["ubiome_relative_lognorm"]]

# =============================================================================
# REBUILD SAMPLE METADATA
# =============================================================================

cat("\nRebuilding sample metadata...\n")

current_seqns <- sample_names(ubiome_relative)
nhanes_id_col <- "SEQN"
nhanes_weight_cols <- character(0)

fg_tables <- all_tables[grepl("_[FG]$", all_tables, ignore.case = TRUE)]
fg_tables <- unique(c(fg_tables,
                      "d_outcome_mcq_f", "d_outcome_mcq_g",
                      "d_outcome_oh_f", "d_outcome_oh_g",
                      "OralDisease_F", "OralDisease_G",
                      "phenotype_vars_f"))

retrieve_sample_data <- function(seqns, variables, tables, con) {
  # Always use a character filter vector: SEQN columns in this DB are stored
  # as TEXT in some F/G tables (e.g. DEMO_F/G) and NUMERIC in others (e.g.
  # BMX_F/G). A numeric filter dbplyr-formats values as "57816.0", which fails
  # a strict TEXT IN-comparison; a character filter works for both column
  # types. Pre-computed outside dplyr to avoid an inline if/else that dbplyr
  # would translate into invalid SQL.
  seqns_chr <- as.character(seqns)
  # Cache per-table column lists so dbListFields() runs once per table, not
  # once per (var, table) pair.
  fields_cache <- new.env(parent = emptyenv())
  get_fields <- function(table_name) {
    if (!exists(table_name, envir = fields_cache, inherits = FALSE)) {
      assign(
        table_name,
        tryCatch(dbListFields(con, table_name), error = function(e) NULL),
        envir = fields_cache
      )
    }
    get(table_name, envir = fields_cache, inherits = FALSE)
  }
  result <- tibble(SEQN = seqns_chr)
  found_vars <- character(0)
  for (var in variables) {
    if (var %in% colnames(result)) {
      next
    }
    var_added <- FALSE
    for (table_name in tables) {
      table_cols <- get_fields(table_name)
      if (is.null(table_cols)) next
      if (!all(c("SEQN", var) %in% table_cols)) next
      var_data <- tryCatch(
        tbl(con, table_name) %>%
          dplyr::select(!!sym("SEQN"), !!sym(var)) %>%
          dplyr::filter(.data$SEQN %in% !!seqns_chr) %>%
          collect() %>%
          mutate(SEQN = as.character(.data$SEQN)),
        error = function(e) NULL
      )
      if (is.null(var_data) || nrow(var_data) == 0) {
        next
      }
      result <- result %>% left_join(var_data, by = "SEQN")
      found_vars <- c(found_vars, var)
      var_added <- TRUE
      break
    }
    if (!var_added) {
      result[[var]] <- NA
    }
  }
  list(
    data = result,
    found_vars = unique(found_vars),
    missing = setdiff(variables, unique(found_vars))
  )
}

cat("  Retrieving variables for", length(all_vars_union), "fields...\n")
metadata_result <- retrieve_sample_data(current_seqns, all_vars_union, fg_tables, con)
sample_data_new <- metadata_result$data

cat("  [done] Retrieved", length(metadata_result$found_vars), "variables\n")
missing_vars <- setdiff(all_vars_union, metadata_result$found_vars)
if (length(missing_vars)) {
  cat("  [warn] Variables unavailable in DB:", length(missing_vars), "(logged for review)\n")
}

sample_data_new <- sample_data_new %>% column_to_rownames("SEQN")

phyloseq_objects_updated <- purrr::map(phyloseq_objects, function(ps_obj) {
  existing_data <- data.frame(sample_data(ps_obj))
  new_columns <- setdiff(colnames(sample_data_new), colnames(existing_data))
  combined_data <- cbind(existing_data,
                         sample_data_new[rownames(existing_data), new_columns, drop = FALSE])
  sample_data(ps_obj) <- sample_data(combined_data)
  ps_obj
})

ubiome_counts <- phyloseq_objects_updated[["ubiome_counts"]]
ubiome_relative <- phyloseq_objects_updated[["ubiome_relative"]]
ubiome_relative_none <- phyloseq_objects_updated[["ubiome_relative_none"]]
ubiome_relative_clr <- phyloseq_objects_updated[["ubiome_relative_clr"]]
ubiome_relative_hellinger <- phyloseq_objects_updated[["ubiome_relative_hellinger"]]
ubiome_relative_lognorm <- phyloseq_objects_updated[["ubiome_relative_lognorm"]]

sample_data_df <- data.frame(sample_data(ubiome_relative_none)) %>%
  rownames_to_column("SEQN")

sample_data_subsets <- imap(all_vars, function(vars, group_name) {
  cols <- unique(c(nhanes_id_col, nhanes_weight_cols, vars))
  sample_data_df %>% dplyr::select(any_of(cols))
})

# Load ubiome genus mapping for downstream node metadata
genus_mapping_path <- file.path(base_path, "results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv")
ubiome_genus_mapping_complete <- read_csv(genus_mapping_path, show_col_types = FALSE)

# =============================================================================
# LOAD ASSOCIATION RESULTS (EXACTLY AS ORIGINAL)
# =============================================================================

cat("\nLoading association results...\n")

# Load WAS Results from New Association Pipeline (EXACTLY as original)
demoWAS_none_tidied  <- if (file.exists(f <- file.path(aggregated_association_res_path, "1_demoWAS_out/result_none/1_demoWAS_none_tidied_complete.rds"))) readRDS(f) else tibble()
demoWAS_none_glanced <- if (file.exists(f <- file.path(aggregated_association_res_path, "1_demoWAS_out/result_none/1_demoWAS_none_glanced_complete.rds"))) readRDS(f) else tibble()
demoWAS_none_rsq     <- if (file.exists(f <- file.path(aggregated_association_res_path, "1_demoWAS_out/result_none/1_demoWAS_none_rsq_complete.rds"))) readRDS(f) else tibble()

demoWAS_clr_tidied   <- if (file.exists(f <- file.path(aggregated_association_res_path, "1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds"))) readRDS(f) else tibble()
demoWAS_clr_glanced  <- if (file.exists(f <- file.path(aggregated_association_res_path, "1_demoWAS_out/result_clr/1_demoWAS_clr_glanced_complete.rds"))) readRDS(f) else tibble()
demoWAS_clr_rsq      <- if (file.exists(f <- file.path(aggregated_association_res_path, "1_demoWAS_out/result_clr/1_demoWAS_clr_rsq_complete.rds"))) readRDS(f) else tibble()

oradWAS_none_tidied  <- if (file.exists(f <- file.path(aggregated_association_res_path, "2_oradWAS_out/result_none/2_oradWAS_none_tidied_complete.rds"))) readRDS(f) else tibble()
oradWAS_none_glanced <- if (file.exists(f <- file.path(aggregated_association_res_path, "2_oradWAS_out/result_none/2_oradWAS_none_glanced_complete.rds"))) readRDS(f) else tibble()
oradWAS_none_rsq     <- if (file.exists(f <- file.path(aggregated_association_res_path, "2_oradWAS_out/result_none/2_oradWAS_none_rsq_complete.rds"))) readRDS(f) else tibble()

oradWAS_clr_tidied   <- if (file.exists(f <- file.path(aggregated_association_res_path, "2_oradWAS_out/result_clr/2_oradWAS_clr_tidied_complete.rds"))) readRDS(f) else tibble()
oradWAS_clr_glanced  <- if (file.exists(f <- file.path(aggregated_association_res_path, "2_oradWAS_out/result_clr/2_oradWAS_clr_glanced_complete.rds"))) readRDS(f) else tibble()
oradWAS_clr_rsq      <- if (file.exists(f <- file.path(aggregated_association_res_path, "2_oradWAS_out/result_clr/2_oradWAS_clr_rsq_complete.rds"))) readRDS(f) else tibble()

exWAS_none_tidied    <- if (file.exists(f <- file.path(aggregated_association_res_path, "3_exWAS_out/result_none/3_exWAS_none_tidied_complete.rds"))) readRDS(f) else tibble()
exWAS_none_glanced   <- if (file.exists(f <- file.path(aggregated_association_res_path, "3_exWAS_out/result_none/3_exWAS_none_glanced_complete.rds"))) readRDS(f) else tibble()
exWAS_none_rsq       <- if (file.exists(f <- file.path(aggregated_association_res_path, "3_exWAS_out/result_none/3_exWAS_none_rsq_complete.rds"))) readRDS(f) else tibble()

exWAS_clr_tidied     <- if (file.exists(f <- file.path(aggregated_association_res_path, "3_exWAS_out/result_clr/3_exWAS_clr_tidied_complete.rds"))) readRDS(f) else tibble()
exWAS_clr_glanced    <- if (file.exists(f <- file.path(aggregated_association_res_path, "3_exWAS_out/result_clr/3_exWAS_clr_glanced_complete.rds"))) readRDS(f) else tibble()
exWAS_clr_rsq        <- if (file.exists(f <- file.path(aggregated_association_res_path, "3_exWAS_out/result_clr/3_exWAS_clr_rsq_complete.rds"))) readRDS(f) else tibble()

pheWAS_none_tidied   <- if (file.exists(f <- file.path(aggregated_association_res_path, "4_pheWAS_out/result_none/4_pheWAS_none_tidied_complete.rds"))) readRDS(f) else tibble()
pheWAS_none_glanced  <- if (file.exists(f <- file.path(aggregated_association_res_path, "4_pheWAS_out/result_none/4_pheWAS_none_glanced_complete.rds"))) readRDS(f) else tibble()
pheWAS_none_rsq      <- if (file.exists(f <- file.path(aggregated_association_res_path, "4_pheWAS_out/result_none/4_pheWAS_none_rsq_complete.rds"))) readRDS(f) else tibble()

pheWAS_clr_tidied    <- if (file.exists(f <- file.path(aggregated_association_res_path, "4_pheWAS_out/result_clr/4_pheWAS_clr_tidied_complete.rds"))) readRDS(f) else tibble()
pheWAS_clr_glanced   <- if (file.exists(f <- file.path(aggregated_association_res_path, "4_pheWAS_out/result_clr/4_pheWAS_clr_glanced_complete.rds"))) readRDS(f) else tibble()
pheWAS_clr_rsq       <- if (file.exists(f <- file.path(aggregated_association_res_path, "4_pheWAS_out/result_clr/4_pheWAS_clr_rsq_complete.rds"))) readRDS(f) else tibble()

outWAS_none_tidied   <- if (file.exists(f <- file.path(aggregated_association_res_path, "5_outWAS_out/result_none/5_outWAS_none_tidied_complete.rds"))) readRDS(f) else tibble()
outWAS_none_glanced  <- if (file.exists(f <- file.path(aggregated_association_res_path, "5_outWAS_out/result_none/5_outWAS_none_glanced_complete.rds"))) readRDS(f) else tibble()
outWAS_none_rsq      <- if (file.exists(f <- file.path(aggregated_association_res_path, "5_outWAS_out/result_none/5_outWAS_none_rsq_complete.rds"))) readRDS(f) else tibble()

outWAS_clr_tidied    <- if (file.exists(f <- file.path(aggregated_association_res_path, "5_outWAS_out/result_clr/5_outWAS_clr_tidied_complete.rds"))) readRDS(f) else tibble()
outWAS_clr_glanced   <- if (file.exists(f <- file.path(aggregated_association_res_path, "5_outWAS_out/result_clr/5_outWAS_clr_glanced_complete.rds"))) readRDS(f) else tibble()
outWAS_clr_rsq       <- if (file.exists(f <- file.path(aggregated_association_res_path, "5_outWAS_out/result_clr/5_outWAS_clr_rsq_complete.rds"))) readRDS(f) else tibble()

cat("Association results loaded:\n")
cat("  demoWAS none:", nrow(demoWAS_none_tidied), "rows\n")
cat("  demoWAS clr:", nrow(demoWAS_clr_tidied), "rows\n")
cat("  oradWAS none:", nrow(oradWAS_none_tidied), "rows\n")
cat("  oradWAS clr:", nrow(oradWAS_clr_tidied), "rows\n")
cat("  exWAS none:", nrow(exWAS_none_tidied), "rows\n")
cat("  exWAS clr:", nrow(exWAS_clr_tidied), "rows\n")
cat("  pheWAS none:", nrow(pheWAS_none_tidied), "rows\n")
cat("  pheWAS clr:", nrow(pheWAS_clr_tidied), "rows\n")
cat("  outWAS none:", nrow(outWAS_none_tidied), "rows\n")
cat("  outWAS clr:", nrow(outWAS_clr_tidied), "rows\n")

# =============================================================================
# BINARY CASE-BALANCE FILTERS (align with original pipeline)
# =============================================================================

cat("\nComputing binary case-balance filters...\n")

get_binary_case_counts <- function(var_list, base_name, con, genus_f, genus_g) {
  table_f <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_f" else paste0(base_name, "_F")
  table_g <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_g" else paste0(base_name, "_G")

  db_tables <- dbListTables(con)
  if (!(table_f %in% db_tables && table_g %in% db_tables)) {
    warning(sprintf("[warn] One or both tables missing: %s, %s", table_f, table_g))
    return(tibble(var_name = var_list, cases_count = NA_integer_))
  }

  df_f <- tbl(con, table_f) %>% collect() %>% mutate(SEQN = as.character(.data$SEQN))
  df_g <- tbl(con, table_g) %>% collect() %>% mutate(SEQN = as.character(.data$SEQN))

  df_all <- bind_rows(
    df_f %>% filter(.data$SEQN %in% genus_f$SEQN),
    df_g %>% filter(.data$SEQN %in% genus_g$SEQN)
  )

  map_dfr(var_list, function(v) {
    if (v %in% names(df_all)) {
      tibble(
        var_name = v,
        cases_count = sum(df_all[[v]] %in% c(1, "1", TRUE), na.rm = TRUE)
      )
    } else {
      tibble(var_name = v, cases_count = NA_integer_)
    }
  })
}

genus_seqn_f <- tbl(con, "DADA2RSV_GENUS_RELATIVE_F") %>% select(SEQN) %>% collect() %>% mutate(SEQN = as.character(.data$SEQN))
genus_seqn_g <- tbl(con, "DADA2RSV_GENUS_RELATIVE_G") %>% select(SEQN) %>% collect() %>% mutate(SEQN = as.character(.data$SEQN))

orad_vars <- readLines(file.path(config_dir_path, "2_oradWAS_vars.txt"), warn = FALSE)
outwas_vars <- readLines(file.path(config_dir_path, "5_outWAS_vars.txt"), warn = FALSE)

oradWAS_case_counts <- get_binary_case_counts(orad_vars, "OralDisease", con, genus_seqn_f, genus_seqn_g)
outWAS_case_counts  <- get_binary_case_counts(outwas_vars, "d_outcome_mcq", con, genus_seqn_f, genus_seqn_g)

binary_case_threshold <- 0.005  # ≥0.5% of 9847 participants
required_cases <- 9847 * binary_case_threshold

oradWAS_case_counts_pass_list <- oradWAS_case_counts %>%
  filter(!is.na(cases_count), cases_count > required_cases) %>%
  pull(var_name)

outWAS_case_counts_pass_list <- outWAS_case_counts %>%
  filter(!is.na(cases_count), cases_count > required_cases) %>%
  pull(var_name)

cat("  oradWAS variables passing case threshold:", length(oradWAS_case_counts_pass_list), "\n")
cat("  outWAS variables passing case threshold:", length(outWAS_case_counts_pass_list), "\n")

# =============================================================================
# PROCESS PREVALENCE THRESHOLD (EXACTLY AS ORIGINAL)
# =============================================================================

cat("\nProcessing prevalence threshold...\n")

# Compute OTUs passing the prevalence threshold (EXACTLY as original)
prevalence_threshold <- 0.01
prevalence_threshold_all <- prevalence_threshold

rsv_genus_relative <- bind_rows(
  tbl(con, 'DADA2RSV_GENUS_RELATIVE_F') %>% collect(),
  tbl(con, 'DADA2RSV_GENUS_RELATIVE_G') %>% collect()
)

otu_non_zero <- rsv_genus_relative %>%
  summarise(across(where(is.numeric), ~ sum(. > 0), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "otu", values_to = "non_zero_count")

otu_pass_prevalance_list_all <- otu_non_zero %>%
  filter(non_zero_count > nrow(rsv_genus_relative) * prevalence_threshold_all,
         otu != "SEQN") %>%
  pull(otu)

cat("  Total samples:", nrow(rsv_genus_relative), "\n")
cat("  Prevalence threshold:", prevalence_threshold_all, "\n")
cat("  OTUs passing threshold:", length(otu_pass_prevalance_list_all), "\n")

# =============================================================================
# PROCESS SIGNIFICANT RESULTS (EXACTLY AS ORIGINAL)
# =============================================================================

cat("\nProcessing significant results...\n")

# fdr_threshold <- 0.05
# fdr_threshold <- 0.001
fdr_threshold <- 0.01
q_threshold <- 0.05

was_specs <- c(
  "demoWAS_none", "demoWAS_clr",
  "oradWAS_none", "oradWAS_clr", 
  "exWAS_none", "exWAS_clr",
  "pheWAS_none", "pheWAS_clr",
  "outWAS_none", "outWAS_clr"
)

clean_otu_id <- function(x) {
  stringr::str_remove(x, "_relative$")
}

significant_results <- list()

for (base in was_specs) {
  tidied_name <- paste0(base, "_tidied")
  sig_res_name <- paste0(base, "_sig_res")
  
  if (!exists(tidied_name, inherits = FALSE)) {
    cat("  ", sig_res_name, ": tidied data not found\n")
    next
  }
  
  dat <- get(tidied_name, inherits = FALSE)
  if (!nrow(dat)) {
    cat("  ", sig_res_name, ": empty table\n")
    next
  }
  
  sig_res <- dat %>%
    mutate(
      otu = case_when(
        str_detect(dependent_var, "^RSV_") ~ dependent_var,
        str_detect(independent_var, "^RSV_") ~ independent_var,
        str_detect(term, "^RSV_") ~ term,
        TRUE ~ phenotype
      ),
      feature_id = otu,
      otu_clean = clean_otu_id(otu)
    ) %>%
    filter(!is.na(independent_var), !is.na(term)) %>%
    filter(str_starts(term, independent_var)) %>%
    filter(!is.na(std.error), std.error != 0) %>%
    filter(!is.na(otu_clean)) %>%
    filter(otu %in% otu_pass_prevalance_list_all | otu_clean %in% otu_pass_prevalance_list_all) %>%
    filter(p.value.fdr <= fdr_threshold, q.value <= q_threshold)

  if (nrow(sig_res)) {
    roles <- determine_roles(sig_res)

    sig_res$host_variable_column <- sig_res[[roles$host_col]]
    sig_res$microbe_variable_column <- sig_res[[roles$microbe_col]]

    if (str_starts(base, "oradWAS") && length(oradWAS_case_counts_pass_list)) {
      sig_res <- sig_res[sig_res$host_variable_column %in% oradWAS_case_counts_pass_list, ]
    }
    if (str_starts(base, "outWAS") && length(outWAS_case_counts_pass_list)) {
      sig_res <- sig_res[sig_res$host_variable_column %in% outWAS_case_counts_pass_list, ]
    }
  }
  
  assign(sig_res_name, sig_res, envir = .GlobalEnv)
  significant_results[[sig_res_name]] <- sig_res
  
  cat("  ", sig_res_name, ":", nrow(sig_res), "significant associations\n")
}

# =============================================================================
# SAVE PROCESSED DATA
# =============================================================================

cat("\nSaving processed data...\n")

# Save phyloseq objects
saveRDS(phyloseq_objects_updated, file.path(data_path, "phyloseq_objects.rds"))

# Save sample data subsets
saveRDS(sample_data_subsets, file.path(data_path, "sample_data_subsets.rds"))

# Save genus mapping
saveRDS(ubiome_genus_mapping_complete, file.path(data_path, "ubiome_genus_mapping_complete.rds"))

# Save significant results
saveRDS(significant_results, file.path(data_path, "significant_results.rds"))

# Save variable descriptions
saveRDS(variable_description, file.path(data_path, "variable_description.rds"))

# Save prevalence list
saveRDS(otu_pass_prevalance_list_all, file.path(data_path, "otu_pass_prevalence_list.rds"))

# Save variable groups
saveRDS(all_vars, file.path(data_path, "variable_groups.rds"))
saveRDS(oradWAS_case_counts_pass_list, file.path(data_path, "oradWAS_case_pass_list.rds"))
saveRDS(outWAS_case_counts_pass_list, file.path(data_path, "outWAS_case_pass_list.rds"))
saveRDS(oradWAS_case_counts, file.path(data_path, "oradWAS_case_counts.rds"))
saveRDS(outWAS_case_counts, file.path(data_path, "outWAS_case_counts.rds"))

cat("[done] Data saved to:", data_path, "\n")

# =============================================================================
# CLEANUP
# =============================================================================

dbDisconnect(con)

cat("\n=== DATA LOADING COMPLETE ===\n")
cat("Ready for network analysis!\n")