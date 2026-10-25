#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 2
# Inventory every phenotype held by every ACQUIRED external cohort
#
# Purpose : Establish, from what is actually on disk, which cohorts can support
#           a model we fit ourselves (Refit cohort: needs BOTH a genus x sample table
#           AND a phenotype) versus only a published target.
#           Nothing here is taken from a paper's claims - a phenotype counts
#           only if the column exists and has non-missing values.
#
# Usage   : Rscript 02_inventory_external_phenotypes.R
# Outputs : audit/external_phenotype_dictionary.tsv
#           audit/external_dataset_capability.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 2: external phenotype inventory ===")

DSD <- file.path(M12_OUT, "datasets")
ids <- sort(list.dirs(DSD, full.names = FALSE, recursive = FALSE))
asay("acquired datasets: ", paste(ids, collapse = ", "))

# Values that mean "absent" rather than a real level.
NULLISH <- c("", "NA", "na", "unknown", "Unknown", "not provided", "not collected",
             "missing", "none", "No Sabe/No Contesta", "not applicable")

is_present <- function(x) {
  x <- as.character(x)
  !is.na(x) & !(trimws(x) %in% NULLISH)
}

#' Locate a genus x sample abundance table for a dataset, if one exists.
#' A cohort without one can never be a refit cohort, no matter what its paper reports.
find_abundance <- function(id) {
  pm <- file.path(DSD, id, "05_processed_microbiome")
  vr <- file.path(DSD, id, "09_validation_ready")
  cand <- c(file.path(vr, "analysis_table.tsv"),
            file.path(pm, "genus_table.tsv"), file.path(pm, "genus_table_relative.tsv"))
  cand <- cand[file.exists(cand)]
  if (!length(cand)) return(NA_character_)
  sub(paste0(DSD, "/"), "", cand[1], fixed = TRUE)
}

dict <- list(); cap <- list()

for (id in ids) {
  smf <- file.path(DSD, id, "07_processed_metadata", "sample_metadata.tsv")
  abund <- find_abundance(id)
  targets <- list.files(file.path(DSD, id, "02_supplementary", "derived"),
                        pattern = "^replication_", full.names = FALSE)
  targets <- c(targets, list.files(file.path(DSD, id, "09_validation_ready"),
                                   pattern = "^replication_", full.names = FALSE))

  n_samples <- NA_integer_; ncols <- 0L
  if (file.exists(smf)) {
    m <- suppressWarnings(readr::read_tsv(smf, show_col_types = FALSE, progress = FALSE,
                                          guess_max = 5000))
    n_samples <- nrow(m)
    # Harmonized columns are the analysis-ready ones; *_original is the raw
    # source kept beside them and is never used for modelling.
    keep <- setdiff(names(m), grep("_original$", names(m), value = TRUE))
    ncols <- length(keep)
    for (v in keep) {
      x <- m[[v]]
      pres <- is_present(x)
      lv <- sort(table(as.character(x)[pres]), decreasing = TRUE)
      dict[[length(dict) + 1]] <- tibble(
        external_dataset = id, external_variable = v,
        n_total = nrow(m), n_present = sum(pres),
        pct_present = round(100 * mean(pres), 1),
        n_levels = length(lv),
        is_numeric = is.numeric(x),
        top_levels = paste(utils::head(names(lv), 4), collapse = " | "))
    }
  }
  # Disposition comes from the dataset config, the single source of truth.
  cfgp <- file.path(M12, "config", paste0(id, ".yaml"))
  disp <- if (file.exists(cfgp)) (yaml::read_yaml(cfgp)$validation_disposition %||% NA_character_) else NA_character_
  cap[[length(cap) + 1]] <- tibble(
    external_dataset = id,
    validation_disposition = disp,
    has_sample_metadata = file.exists(smf),
    n_samples = n_samples,
    n_harmonized_columns = ncols,
    abundance_table = abund %||% NA_character_,
    has_abundance_table = !is.na(abund),
    # refit-cohort requires BOTH: a per-sample abundance table and per-sample
    # phenotypes. Everything else can contribute published targets only.
    tier = ifelse(!is.na(abund) & file.exists(smf), "A", "B"),
    n_published_targets = length(targets),
    published_targets = paste(targets, collapse = "; "))
}

dict_all <- bind_rows(dict); cap_all <- bind_rows(cap)
write_tsv(dict_all, file.path(AUD, "external_phenotype_dictionary.tsv"))
write_tsv(cap_all,  file.path(AUD, "external_dataset_capability.tsv"))

asay("")
asay("=== dataset capability ===")
print(as.data.frame(cap_all %>% select(external_dataset, validation_disposition, tier,
                                       n_samples, n_harmonized_columns,
                                       has_abundance_table, n_published_targets)))
asay("")
asay("=== phenotypes with >=50% coverage in a TIER A cohort (fittable ourselves) ===")
tierA <- cap_all$external_dataset[cap_all$tier == "A"]
print(as.data.frame(dict_all %>%
  filter(external_dataset %in% tierA, pct_present >= 50,
         !external_variable %in% c("sample_id", "subject_id", "study_id",
                                   "n_samples_this_subject", "is_subject_representative")) %>%
  select(external_dataset, external_variable, n_present, pct_present, n_levels, top_levels)))
