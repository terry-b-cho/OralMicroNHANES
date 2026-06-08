#!/usr/bin/env Rscript

# =============================================================================
# Module 10 - prediction pipeline, step 1: build features + targets + weights
# =============================================================================
# Inputs (relative to PROJECT_ROOT):
#   - results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_clr.rds
#   - data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt
#   - data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite
#
# Outputs (under results/analyses_results/10_prediction_analyses_out/intermediate/):
#   - clr_features.rds              prevalence-filtered CLR features (taxa x samples)
#   - alpha_features.rds            three alpha metrics, mean across 10 resamplings
#   - covariates_and_targets.rds    demographics + survey design vars + targets
#   - sample_universe.rds           SEQNs with CLR + alpha + DEMO data
#
# Environment: R >= 4.5 with phyloseq, DBI, RSQLite, dplyr, tidyr, tibble,
# data.table, readr, stringr. Conda spec: envs/nhanes-analysis_for_reviewers.yml.
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
  library(tibble)
  library(data.table)
  library(readr)
  library(stringr)
})

set.seed(42)

base_path     <- PROJECT_ROOT
clr_rds_path  <- file.path(base_path, "results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative_clr.rds")
alpha_path    <- file.path(base_path, "data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt")
db_path       <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")

out_root      <- file.path(base_path, "results/analyses_results/10_prediction_analyses_out")
out_inter     <- file.path(out_root, "intermediate")
dir.create(out_inter, recursive = TRUE, showWarnings = FALSE)

# ---- 1. CLR features (prevalence-filtered from COUNT data) ------------------
# Note: prevalence must be computed on COUNT data (zeros are zeros there). CLR-
# transformed values are almost never exactly zero — using counts gives the
# correct ~219-taxa subset that downstream analyses rely on.
message("== Loading CLR + COUNT phyloseq objects (CLR -> features, counts -> prevalence)")
ps_clr     <- readRDS(clr_rds_path)
counts_path <- file.path(base_path,
  "results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_counts.rds")
ps_counts  <- readRDS(counts_path)

count_mat <- as(otu_table(ps_counts), "matrix")
if (!taxa_are_rows(ps_counts)) count_mat <- t(count_mat)
prev_counts  <- rowMeans(count_mat > 0, na.rm = TRUE)
keep_taxa    <- names(prev_counts)[prev_counts >= 0.01]

clr_mat <- as(otu_table(ps_clr), "matrix")
if (!taxa_are_rows(ps_clr)) clr_mat <- t(clr_mat)
clr_mat <- clr_mat[intersect(keep_taxa, rownames(clr_mat)), , drop = FALSE]
rownames(clr_mat) <- paste0("clr_", rownames(clr_mat))

clr_features <- as.data.frame(t(clr_mat)) |>
  tibble::rownames_to_column("SEQN") |>
  mutate(SEQN = as.character(SEQN))
message(sprintf("   CLR features: %d samples x %d taxa (after >=1%% non-zero count prevalence)",
                nrow(clr_features), ncol(clr_features) - 1))
saveRDS(clr_features, file.path(out_inter, "clr_features.rds"))

# ---- 2. Alpha diversity features (depth 10k, averaged over 10 resamplings) -
message("== Loading alpha diversity")
alpha_raw <- fread(alpha_path)
alpha_raw[, SEQN := as.character(SEQN)]

mean_across_resamples <- function(prefix) {
  cols <- grep(paste0("^", prefix, "_10000_"), names(alpha_raw), value = TRUE)
  if (length(cols) == 0) return(rep(NA_real_, nrow(alpha_raw)))
  m <- as.matrix(alpha_raw[, ..cols])
  m <- apply(m, 2, as.numeric)
  rowMeans(m, na.rm = TRUE)
}

alpha_features <- data.frame(
  SEQN              = alpha_raw$SEQN,
  alpha_observed    = mean_across_resamples("RSV_ObservedOTUs"),
  alpha_shannon     = mean_across_resamples("RSV_ShanWienDiv"),
  alpha_inv_simpson = mean_across_resamples("RSV_InverseSimpson"),
  stringsAsFactors  = FALSE
) |> filter(!is.na(alpha_observed) | !is.na(alpha_shannon) | !is.na(alpha_inv_simpson))
message(sprintf("   Alpha features: %d samples x 3 metrics", nrow(alpha_features)))
saveRDS(alpha_features, file.path(out_inter, "alpha_features.rds"))

# ---- 3. Covariates + targets + survey design (DEMO_F + DEMO_G) -------------
message("== Loading DEMO_F + DEMO_G for covariates + targets + survey design")
con <- dbConnect(SQLite(), db_path)
on.exit(try(dbDisconnect(con), silent = TRUE), add = TRUE)

demo_cols <- c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1", "DMDEDUC2", "INDFMPIR",
               "WTMEC2YR", "SDMVSTRA", "SDMVPSU")

read_demo <- function(tab, cycle_label) {
  fields <- dbListFields(con, tab)
  pick   <- intersect(demo_cols, fields)
  d <- dbGetQuery(con, sprintf("SELECT %s FROM %s",
                               paste(sprintf('"%s"', pick), collapse = ","), tab))
  d$SEQN  <- as.character(d$SEQN)
  d$cycle <- cycle_label
  d
}

demo_f <- read_demo("DEMO_F", "F")
demo_g <- read_demo("DEMO_G", "G")
demo   <- bind_rows(demo_f, demo_g) |>
  filter(!is.na(SEQN)) |>
  distinct(SEQN, .keep_all = TRUE)

age_group_labels <- c("14-19","20-29","30-39","40-49","50-59","60-69")
covariates_and_targets <- demo |>
  transmute(
    SEQN,
    cycle,
    age           = as.numeric(RIDAGEYR),
    age_group     = factor(case_when(
      age >= 14 & age < 20 ~ "14-19",
      age >= 20 & age < 30 ~ "20-29",
      age >= 30 & age < 40 ~ "30-39",
      age >= 40 & age < 50 ~ "40-49",
      age >= 50 & age < 60 ~ "50-59",
      age >= 60 & age <= 69 ~ "60-69",
      TRUE ~ NA_character_
    ), levels = age_group_labels, ordered = TRUE),
    gender        = factor(case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ), levels = c("Male","Female")),
    ethnicity     = factor(case_when(
      RIDRETH1 == 1 ~ "Mexican",
      RIDRETH1 == 2 ~ "OtherHispanic",
      RIDRETH1 == 3 ~ "White",
      RIDRETH1 == 4 ~ "Black",
      RIDRETH1 == 5 ~ "Other",
      TRUE ~ NA_character_
    )),
    education_lt9 = as.integer(DMDEDUC2 == 1),
    education_hs  = as.integer(DMDEDUC2 == 3),
    education_aa  = as.integer(DMDEDUC2 == 4),
    education_cg  = as.integer(DMDEDUC2 == 5),
    pir           = as.numeric(INDFMPIR),
    wtmec2yr      = as.numeric(WTMEC2YR),
    sdmvstra      = as.integer(SDMVSTRA),
    sdmvpsu       = as.integer(SDMVPSU)
  )

dbDisconnect(con)
message(sprintf("   Covariates + targets: %d samples x %d columns",
                nrow(covariates_and_targets), ncol(covariates_and_targets)))
saveRDS(covariates_and_targets, file.path(out_inter, "covariates_and_targets.rds"))

# ---- 4. Sample universe = CLR n alpha n DEMO -------------------------------
sample_universe <- Reduce(intersect, list(
  clr_features$SEQN,
  alpha_features$SEQN,
  covariates_and_targets$SEQN
))
message(sprintf("== Sample universe (CLR n alpha n DEMO): %d SEQNs", length(sample_universe)))
saveRDS(sample_universe, file.path(out_inter, "sample_universe.rds"))

message("Done. Intermediates written to: ", out_inter)
