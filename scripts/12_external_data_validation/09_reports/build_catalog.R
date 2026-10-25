#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Build the cross-dataset catalog and availability matrices
#
# Purpose : One machine-readable row per dataset/cohort, plus the two matrices
#           and the ranked shortlist. Everything is DERIVED from each dataset's
#           config + measured outputs - never hand-written - so the catalog
#           cannot drift from what is actually on disk.
#
# Usage   : Rscript 09_reports/build_catalog.R
# Outputs : dataset_catalog.tsv
#           phenotype_availability_matrix.tsv
#           sequencing_method_matrix.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(yaml); library(tibble); library(tidyr) })

PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
set.seed(42)
`%||%` <- function(a, b) if (is.null(a)) b else a

M12  <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
OUT  <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation")
DSD  <- file.path(OUT, "datasets")
say  <- function(...) cat(..., "\n", sep = "")

cfgs <- list.files(file.path(M12, "config"), pattern = "\\.yaml$", full.names = TRUE)
say("configs found: ", length(cfgs))

rows <- list(); phen <- list(); seqm <- list()

for (cf in cfgs) {
  cfg <- yaml::read_yaml(cf)
  id  <- cfg$dataset_id %||% sub("\\.yaml$", "", basename(cf))
  d   <- file.path(DSD, id)

  # --- measured, not declared ---
  n_samples_pheno <- NA_integer_; n_join <- NA_integer_; n_genera <- NA_integer_
  n_overlap <- NA_integer_; cov <- list()

  mf <- file.path(d, "07_processed_metadata", "sample_metadata.tsv")
  if (file.exists(mf))
    n_samples_pheno <- nrow(suppressWarnings(readr::read_tsv(mf, show_col_types = FALSE,
                                                             progress = FALSE)))

  # Coverage is reported for the DELIVERED table, not the intermediate phenotype
  # layer: supplementary phenotype sources are merged during assembly, so
  # 07_processed_metadata understates what an analyst actually receives.
  af <- file.path(d, "09_validation_ready", "analysis_table.tsv")
  cov_src <- if (file.exists(af)) af else mf
  if (file.exists(cov_src)) {
    m <- suppressWarnings(readr::read_tsv(cov_src, show_col_types = FALSE, progress = FALSE))
    if (identical(cov_src, af)) n_join <- nrow(m)
    for (v in c("age_years", "sex_harmonized", "smoking_status_harmonized", "bmi",
                "ethnicity", "disease_status", "cotinine", "alcohol_g_per_day",
                "gum_disease", "assay_batch")) {
      cov[[v]] <- if (!v %in% names(m)) 0 else {
        x <- m[[v]]
        if (is.character(x) || is.factor(x)) mean(!is.na(x) & x != "unknown") else mean(!is.na(x))
      }
    }
  }
  of <- file.path(d, "05_processed_microbiome", "genus_overlap_nhanes.tsv")
  if (file.exists(of)) {
    o <- readr::read_tsv(of, show_col_types = FALSE, progress = FALSE)
    n_genera <- nrow(o); n_overlap <- sum(o$in_nhanes)
  }

  scs <- cfg$subcohorts %||% list()
  accs <- paste(vapply(scs, function(s) s$accession %||% "", character(1)), collapse = ";")
  sites <- paste(unique(vapply(scs, function(s) s$oral_site %||% "", character(1))), collapse = ";")
  regs  <- paste(unique(vapply(scs, function(s) s$region_16S %||% "", character(1))), collapse = ";")
  plats <- paste(unique(vapply(scs, function(s) s$platform %||% "", character(1))), collapse = ";")

  # A dataset with no individual-level data is not "pending" - it is a finished
  # result. Distinguish "we have not looked" from "we looked and there is
  # nothing", otherwise the catalog understates completed work.
  has_report  <- file.exists(file.path(d, "11_documentation", "DATASET_REPORT.md"))
  has_targets <- length(Sys.glob(file.path(d, "02_supplementary", "derived", "replication_*"))) > 0 ||
                 length(Sys.glob(file.path(d, "09_validation_ready", "replication_*"))) > 0
  status <- if (file.exists(af)) "validation_ready"
            else if (file.exists(mf) && has_targets) "phenotype_and_targets"
            else if (file.exists(mf)) "phenotype_built"
            else if (has_targets && has_report) "targets_only"
            else if (has_report) "documented_no_data"
            else if (dir.exists(file.path(d, "03_repository_metadata")) &&
                     length(list.files(file.path(d, "03_repository_metadata")))) "metadata_downloaded"
            else "pending"

  rows[[id]] <- tibble(
    dataset_id = id,
    citation = paste0(cfg$first_author %||% "", " ", cfg$year %||% "", ", ", cfg$journal %||% ""),
    doi = cfg$doi %||% NA, pmid = cfg$pmid %||% NA, pmcid = cfg$pmcid %||% NA,
    year = cfg$year %||% NA,
    n_subcohorts = length(scs),
    all_accessions = accs,
    oral_site = sites, sequencing_region = regs, sequencing_platform = plats,
    n_samples_phenotype = n_samples_pheno,
    n_samples_validation_ready = n_join,
    n_genera = n_genera, n_genera_shared_nhanes = n_overlap,
    tier = cfg$microbiome$tier %||% NA,
    code_repository = cfg$code_repository %||% NA,
    n_critical_caveats = length(cfg$critical_caveats %||% list()),
    has_replication_targets = has_targets,
    status = status,
    # "retained" means CONTRIBUTED EVIDENCE to the validation - never merely
    # acquired. Declared per dataset in config/<id>.yaml; the controlled
    # vocabulary is asserted below so a typo cannot quietly create a new class.
    validation_disposition = cfg$validation_disposition %||% NA_character_,
    validation_disposition_reason = cfg$validation_disposition_reason %||% NA_character_,
    dataset_report = if (file.exists(file.path(d, "11_documentation", "DATASET_REPORT.md")))
                       file.path("datasets", id, "11_documentation/DATASET_REPORT.md") else NA
  )

  for (v in names(cov)) {
    phen[[length(phen)+1]] <- tibble(dataset_id = id, variable = v,
                                     coverage = round(100 * cov[[v]], 1))
  }
  for (s in scs) {
    seqm[[length(seqm)+1]] <- tibble(
      dataset_id = id, subcohort = s$id %||% NA, accession = s$accession %||% NA,
      design = s$design %||% NA, oral_site = s$oral_site %||% NA,
      region_16S = s$region_16S %||% NA, platform = s$platform %||% NA,
      layout = s$layout %||% NA, n_runs = s$n_runs %||% NA,
      nhanes_comparability =
        (cfg$nhanes_comparability_by_subcohort[[s$id %||% ""]]$rating) %||% NA)
  }
}

cat_tbl <- bind_rows(rows) %>% arrange(dataset_id)

# YAML folded scalars (`>`) carry a trailing newline, and prose fields can carry
# tabs. Writing those straight to TSV splits one logical row across several
# physical lines and silently corrupts the file for every downstream reader.
# Flatten every character cell and assert the result.
flatten_cells <- function(df) {
  dplyr::mutate(df, dplyr::across(dplyr::where(is.character), function(x) {
    x <- gsub("[\r\n\t]+", " ", x); x <- gsub(" {2,}", " ", x); trimws(x)
  }))
}
cat_tbl <- flatten_cells(cat_tbl)
stopifnot(!any(grepl("[\r\n\t]", unlist(cat_tbl[sapply(cat_tbl, is.character)]),
                     useBytes = TRUE)))

DISPOSITIONS <- c("retained", "excluded_context_only", "excluded_reference_only")
if (any(is.na(cat_tbl$validation_disposition)))
  stop("MISSING validation_disposition in config for: ",
       paste(cat_tbl$dataset_id[is.na(cat_tbl$validation_disposition)], collapse = ", "))
if (!all(cat_tbl$validation_disposition %in% DISPOSITIONS))
  stop("UNKNOWN validation_disposition: ",
       paste(setdiff(cat_tbl$validation_disposition, DISPOSITIONS), collapse = ", "))
say("disposition: ", sum(cat_tbl$validation_disposition == "retained"), " retained, ",
    sum(cat_tbl$validation_disposition != "retained"), " excluded (acquired but unused)")
readr::write_tsv(cat_tbl, file.path(OUT, "dataset_catalog.tsv"))
say("wrote dataset_catalog.tsv (", nrow(cat_tbl), " datasets)")

if (length(phen)) {
  pm <- bind_rows(phen) %>% pivot_wider(names_from = variable, values_from = coverage)
  readr::write_tsv(flatten_cells(pm), file.path(OUT, "phenotype_availability_matrix.tsv"))
  say("wrote phenotype_availability_matrix.tsv (coverage %, per dataset)")
}
if (length(seqm)) {
  sm <- bind_rows(seqm)
  readr::write_tsv(flatten_cells(sm), file.path(OUT, "sequencing_method_matrix.tsv"))
  say("wrote sequencing_method_matrix.tsv (", nrow(sm), " subcohorts)")
}

say("")
say("=== catalog ===")
print(as.data.frame(cat_tbl %>% select(dataset_id, validation_disposition, status, n_samples_validation_ready,
                                       n_genera_shared_nhanes, n_critical_caveats)))
