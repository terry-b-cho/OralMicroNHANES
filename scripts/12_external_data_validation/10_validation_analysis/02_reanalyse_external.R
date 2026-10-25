#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 2 — Step 2
# Refit cohorts: re-analyse external datasets with a harmonized model
#
# Purpose : For the three external datasets where we hold a genus x sample
#           table, fit the SAME association NHANES fitted, so effects are
#           directly comparable and a non-replication is genuinely discordant
#           rather than merely unreported.
#
#           Deliberate design choices, each of which changes the answer:
#
#           * AGE is fitted TWICE - unadjusted (matching NHANES demoWAS, which
#             has zero covariates) and adjusted for sex. The unadjusted fit is
#             the like-for-like comparison; the adjusted fit measures how much
#             the NHANES adjustment asymmetry matters.
#           * SMOKING has no refit cohort. The smoking exposure builders below
#             are retained and tested but currently unexercised; smoking
#             evidence is published-target only.
#           * Adults only (>= 18). No adolescent cohort contributes anywhere.
#           * Huang is deduplicated to one sample per subject before fitting
#             (2,550 samples come from 731 people; four contribute 971).
#
# Usage   : Rscript 02_reanalyse_external.R
# Outputs : validation/external_effects_<dataset>_<variable>[_<fit>].tsv
#           validation/external_fit_audit.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(broom)
})
# When another script source()s this one for its fitting functions, commandArgs
# points at THAT script's path, so honour an explicit override first.
HERE <- Sys.getenv("M12_V2_DIR", "")
if (!nzchar(HERE))
  HERE <- dirname(sub("^--file=", "",
                      grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/10_validation_analysis"
source(file.path(HERE, "00_validation_config.R"))

vsay("=== Step 2: refit-cohort external re-analysis ===")

DSD <- file.path(M12_OUT, "datasets")
audit <- list(); outputs <- list()

# --- CLR with a value-type-aware pseudocount ---------------------------------
# Counts: integer eps = 1, matching module 0.
# Relative abundance: eps = half the smallest non-zero value (module 5.5's
# convention). Using 1 on proportions would swamp the data entirely.
clr_for <- function(mat, value_type) {
  m <- as.matrix(mat); storage.mode(m) <- "double"; m[is.na(m)] <- 0
  pc <- if (identical(value_type, "count")) 1 else {
    nz <- m[m > 0]; if (!length(nz)) 1e-6 else min(nz) / 2
  }
  out <- m12_clr(m, pseudocount = pc)
  attr(out, "pseudocount") <- pc
  out
}

#' Fit genus_clr ~ exposure (+ covariates) for every genus and tidy the
#' exposure term only.
fit_genus_models <- function(clr_mat, meta, exposure, covariates, label) {
  covs <- covariates[covariates %in% names(meta)]
  covs <- covs[vapply(covs, function(c) dplyr::n_distinct(meta[[c]][!is.na(meta[[c]])]) > 1,
                      logical(1))]
  rhs <- paste(c(exposure, covs), collapse = " + ")
  res <- lapply(colnames(clr_mat), function(g) {
    df <- meta; df$.y <- clr_mat[, g]
    df <- df[stats::complete.cases(df[, c(".y", exposure, covs), drop = FALSE]), , drop = FALSE]
    if (nrow(df) < 20) return(NULL)
    fit <- try(stats::lm(stats::as.formula(paste(".y ~", rhs)), data = df), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    tt <- broom::tidy(fit)
    # Match the exposure term, allowing for a factor-dummy suffix.
    row <- tt[startsWith(tt$term, exposure) & tt$term != "(Intercept)", , drop = FALSE]
    if (!nrow(row)) return(NULL)
    tibble(genus = g, term = row$term[1], estimate = row$estimate[1],
           std.error = row$std.error[1], statistic = row$statistic[1],
           p.value = row$p.value[1], n_obs = nrow(df))
  })
  out <- bind_rows(res)
  if (!nrow(out)) return(out)
  out %>% mutate(p.value.fdr = stats::p.adjust(p.value, "BH"),
                 direction = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", NA)),
                 sig_fdr = p.value.fdr < EXTERNAL_FDR,
                 model = label, covariates_used = paste(covs, collapse = "+"),
                 n_covariates = length(covs))
}

load_dataset <- function(id) {
  gt <- file.path(DSD, id, "05_processed_microbiome", "genus_table.tsv")
  at <- file.path(DSD, id, "09_validation_ready", "analysis_table.tsv")
  stopifnot(file.exists(gt), file.exists(at))
  g <- suppressWarnings(readr::read_tsv(gt, show_col_types = FALSE, progress = FALSE))
  a <- suppressWarnings(readr::read_tsv(at, show_col_types = FALSE, progress = FALSE))
  key_g <- names(g)[1]
  key_a <- intersect(c(key_g, "sample_id", "run_accession"), names(a))[1]
  stopifnot(!is.na(key_a))
  list(g = g, a = a, key_g = key_g, key_a = key_a)
}

run_dataset <- function(id, variable, value_type, adults_only, subject_dedup,
                        exposure_builder, covariate_sets, extra_filter = NULL) {
  vsay("--- ", id, " / ", variable, " ---")
  D <- load_dataset(id)
  g <- D$g; a <- D$a

  # Align the abundance matrix to the metadata on the join key.
  ids <- as.character(g[[D$key_g]])
  a$.key <- as.character(a[[D$key_a]])
  meta <- a[match(ids, a$.key), , drop = FALSE]
  keep <- !is.na(meta$.key)
  vsay("  samples: ", length(ids), " | matched to metadata: ", sum(keep))

  # Filters, each asserted rather than assumed.
  if (adults_only) {
    ok <- keep & !is.na(meta$age_years) & meta$age_years >= MIN_ADULT_AGE
    vsay("  adults (>=", MIN_ADULT_AGE, "): ", sum(ok, na.rm = TRUE), " of ", sum(keep))
    keep <- ok
  }
  if (subject_dedup) {
    stopifnot("is_subject_representative" %in% names(meta))
    rep_flag <- meta$is_subject_representative
    rep_flag <- if (is.logical(rep_flag)) rep_flag else toupper(as.character(rep_flag)) == "TRUE"
    keep <- keep & !is.na(rep_flag) & rep_flag
    vsay("  after subject dedup: ", sum(keep, na.rm = TRUE))
  }
  if (!is.null(extra_filter)) {
    keep <- keep & extra_filter(meta)
    vsay("  after extra filter: ", sum(keep, na.rm = TRUE))
  }

  gm <- as.matrix(g[keep, setdiff(names(g), D$key_g), drop = FALSE])
  meta <- meta[keep, , drop = FALSE]
  # Drop genera that are all-zero in the retained samples.
  nz <- colSums(gm, na.rm = TRUE) > 0
  gm <- gm[, nz, drop = FALSE]
  vsay("  genera with any signal: ", ncol(gm))

  clr <- clr_for(gm, value_type)
  stopifnot(max(abs(rowSums(clr))) < 1e-8)   # CLR rows must sum to zero
  vsay("  CLR pseudocount: ", signif(attr(clr, "pseudocount"), 3),
       " | row-sum check passed")

  eb <- exposure_builder(meta)
  meta <- eb$meta; exposure <- eb$exposure
  vsay("  exposure '", exposure, "': ", eb$describe)

  out <- list()
  for (nm in names(covariate_sets)) {
    r <- fit_genus_models(clr, meta, exposure, covariate_sets[[nm]], nm)
    if (!nrow(r)) { vsay("  [warn] no models fitted for '", nm, "'"); next }
    r <- r %>% mutate(dataset_id = id, variable = variable, value_type = value_type)
    vsay(sprintf("  fit '%s': %d genera, n=%d, covars={%s}, FDR<%.2f: %d",
                 nm, nrow(r), max(r$n_obs), r$covariates_used[1], EXTERNAL_FDR,
                 sum(r$sig_fdr)))
    f <- file.path(VAL, sprintf("external_effects_%s_%s_%s.tsv", id, variable, nm))
    readr::write_tsv(r, f)
    out[[nm]] <- r
    audit[[length(audit) + 1]] <<- tibble(
      dataset_id = id, variable = variable, fit = nm,
      n_participants = max(r$n_obs), n_genera_fitted = nrow(r),
      n_fdr_sig = sum(r$sig_fdr), covariates = r$covariates_used[1],
      adults_only = adults_only, subject_dedup = subject_dedup,
      value_type = value_type, pseudocount = signif(attr(clr, "pseudocount"), 3))
  }
  out
}

# --- Exposure builders -------------------------------------------------------
age_exposure <- function(meta) {
  list(meta = meta, exposure = "age_years",
       describe = sprintf("continuous, %.0f-%.0f y",
                          min(meta$age_years, na.rm = TRUE), max(meta$age_years, na.rm = TRUE)))
}
# Current vs never: former is DROPPED so the contrast matches datasets that
# cannot separate former from never.
smoking_cn <- function(meta) {
  s <- as.character(meta$smoking_status_harmonized)
  meta$smk_cn <- factor(ifelse(s == "current", "current",
                        ifelse(s == "never", "never", NA)), levels = c("never", "current"))
  list(meta = meta, exposure = "smk_cn",
       describe = sprintf("current=%d vs never=%d (former dropped)",
                          sum(meta$smk_cn == "current", na.rm = TRUE),
                          sum(meta$smk_cn == "never", na.rm = TRUE)))
}
# Sensitivity: current vs everything else, i.e. the binary that cohorts which
# cannot separate former from never are forced to use.
smoking_cnc <- function(meta) {
  s <- as.character(meta$smoking_status_harmonized)
  meta$smk_cnc <- factor(ifelse(s == "current", "current",
                         ifelse(s %in% c("never", "former"), "non_current", NA)),
                         levels = c("non_current", "current"))
  list(meta = meta, exposure = "smk_cnc",
       describe = sprintf("current=%d vs non-current=%d",
                          sum(meta$smk_cnc == "current", na.rm = TRUE),
                          sum(meta$smk_cnc == "non_current", na.rm = TRUE)))
}

# --- Run ---------------------------------------------------------------------
# Phase 3's audit reuses clr_for(), fit_genus_models(), load_dataset() and
# run_dataset() to fit ADDITIONAL phenotypes in the same cohort under the same
# conventions. Sourcing this file with M12_HELPERS_ONLY=1 defines them without
# re-running the age analysis below, so there is exactly one implementation.
if (identical(Sys.getenv("M12_HELPERS_ONLY"), "1")) {
  vsay("[helpers-only] fitting functions defined; skipping the Phase 2 run block")
} else {

# Huang 2020 is the only refit cohort: age only (5 current smokers make smoking
# unusable). Tongue-swab dominated - the site caveat travels with every result.
#
# Two fits, because one adjustment is not like the other:
#   unadjusted - matches NHANES demoWAS exactly (zero covariates). PRIMARY.
#   adjusted   - sex. Measures how much the NHANES adjustment asymmetry matters.
run_dataset("HuangS_2020_AgePrediction", "age", "relative", TRUE, TRUE,
            age_exposure, list(unadjusted = character(0),
                               adjusted = c("sex_harmonized")))

aud <- bind_rows(audit)
readr::write_tsv(aud, file.path(VAL, "external_fit_audit.tsv"))
vsay("")
vsay("=== fit audit ===")
print(as.data.frame(aud))

}  # end of the Phase 2 run block
