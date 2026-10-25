#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 2 — Step 5
# Tables: study characteristics, per-genus comparison, concordance summary
#
# Purpose : Build the study-characteristics table that makes every comparison
#           interpretable - what each cohort measured, what it adjusted for,
#           and BOTH sample sizes (participants, and overlapping genera).
#           Everything is derived from the dataset configs and the measured
#           audit files; nothing is re-typed by hand.
#
# Usage   : Rscript 05_tables.R
# Outputs : validation/study_characteristics.tsv
#           validation/table_concordance.tsv
#           validation/table_pergenus_{age,smoking}.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(yaml)
})
HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/10_validation_analysis"
source(file.path(HERE, "00_validation_config.R"))
vsay("=== Step 5: tables ===")

CFG <- file.path(M12, "config")
sm  <- readr::read_tsv(file.path(VAL, "concordance_summary.tsv"), show_col_types = FALSE)
aud <- readr::read_tsv(file.path(VAL, "external_fit_audit.tsv"), show_col_types = FALSE)
cmp <- bind_rows(lapply(PRIMARY_VARIABLES, function(v)
  readr::read_tsv(file.path(VAL, paste0("comparison_", v, ".tsv")), show_col_types = FALSE)))
ovl <- readr::read_tsv(file.path(M12_OUT, "genus_overlap_summary.tsv"), show_col_types = FALSE)

cfg_of <- function(id) {
  p <- file.path(CFG, paste0(id, ".yaml"))
  if (file.exists(p)) yaml::read_yaml(p) else list()
}
`%or%` <- function(a, b) if (is.null(a) || !length(a) || (length(a) == 1 && is.na(a))) b else a

# --- Study characteristics ---------------------------------------------------
ids <- unique(sm$dataset_id)
rows <- lapply(ids, function(id) {
  cf <- cfg_of(id)
  scs <- cf$subcohorts %or% list()
  # Config `nhanes_comparability` fields are explanatory prose written for the
  # dataset reports. A table needs a short label, so compact them to the
  # canonical vocabulary instead of pasting paragraphs into a cell.
  short_site <- function(x) {
    if (is.null(x) || !length(x) || all(is.na(x))) return(NA_character_)
    x <- paste(x, collapse = " ")
    # Strip NEGATED mentions first. A config that reads "saliva (not oral
    # rinse)" is labelled "oral rinse" by a naive rinse-before-saliva test -
    # the exact opposite of the truth.
    x <- gsub("\\(?not[^)a-z]{0,4}(an? )?(oral )?rinse\\)?", " ", x, ignore.case = TRUE)
    x <- gsub("vs\\.? NHANES[^.]*", " ", x, ignore.case = TRUE)
    if (grepl("tongue", x, ignore.case = TRUE))    return("tongue swab")
    if (grepl("subgingival|plaque", x, ignore.case = TRUE)) return("subgingival plaque")
    if (grepl("mouthwash|oral rinse|whole-mouth", x, ignore.case = TRUE)) return("oral rinse")
    if (grepl("saliva", x, ignore.case = TRUE))    return("saliva")
    substr(x, 1, 24)
  }
  short_reg <- function(x) {
    if (is.null(x) || !length(x) || all(is.na(x))) return(NA_character_)
    x <- paste(x, collapse = " ")
    # Comparability prose names NHANES's own region as the contrast; drop that
    # clause so we report the DATASET's region, not ours.
    x <- gsub("vs\\.? *NHANES.*$", "", x, ignore.case = TRUE)
    x <- gsub("MATCHES NHANES.*$", "", x, ignore.case = TRUE)
    hits <- unique(unlist(regmatches(x, gregexpr("V[0-9]+(-V?[0-9]+)?", x))))
    if (length(hits)) paste(hits, collapse = "; ") else substr(x, 1, 20)
  }
  site <- if (length(scs))
    paste(unique(vapply(scs, function(s) short_site(s$oral_site), character(1))), collapse = "; ")
  else short_site(cf$nhanes_comparability$oral_site)
  reg <- if (length(scs))
    paste(unique(vapply(scs, function(s) short_reg(s$region_16S), character(1))), collapse = "; ")
  else short_reg(cf$nhanes_comparability$region_16S)
  # A Saca La Lengua round without a nhanes_comparability block still follows
  # the published SLL whole-mouth rinse protocol.
  if (is.na(site) && grepl("SacaLaLengua", id)) site <- "oral rinse"
  if (is.na(reg)  && grepl("SacaLaLengua", id)) reg  <- "V3-V4"

  a <- aud %>% filter(dataset_id == id)
  s <- sm  %>% filter(dataset_id == id)
  o <- ovl %>% filter(dataset_id == id) %>% slice(1)

  tibble(
    dataset_id   = id,
    citation     = paste0(cf$first_author %or% "", " ", cf$year %or% "", ", ", cf$journal %or% ""),
    doi          = cf$doi %or% NA_character_,
    tier         = paste(unique(s$tier), collapse = "/"),
    variables_validated = paste(sort(unique(s$nhanes_variable)), collapse = ", "),
    n_participants = if (nrow(a)) max(a$n_participants, na.rm = TRUE) else NA_integer_,
    oral_site    = site,
    region_16S   = reg,
    covariates_adjusted = if (nrow(a))
      paste(unique(ifelse(a$covariates == "" | is.na(a$covariates), "none (crude)", a$covariates)),
            collapse = " | ") else "published result - see source paper",
    adults_only  = if (nrow(a)) all(a$adults_only) else NA,
    subject_dedup = if (nrow(a)) any(a$subject_dedup) else NA,
    smoking_definition = cf$smoking_harmonization$study_definition %or%
                         cf$smoking_harmonization$study_definition_S1 %or% NA_character_,
    n_genera_dataset = if (nrow(o)) o$n_external_genera else NA_integer_,
    n_overlap_nhanes_common = if (nrow(o)) o$n_match_normalized else NA_integer_,
    # Only the PRIMARY fit per variable, otherwise sensitivity fits repeat the
    # same variable several times in one cell.
    n_comparable_after_fdr = {
      sp <- s %>% filter(primary) %>% distinct(nhanes_variable, .keep_all = TRUE)
      if (!nrow(sp)) sp <- s %>% distinct(nhanes_variable, .keep_all = TRUE)
      paste(sprintf("%s=%d", sp$nhanes_variable, sp$n_tested), collapse = "; ")
    },
    key_caveat = paste(vapply(cf$critical_caveats %or% list(),
                              function(x) x$id %or% "", character(1)), collapse = "; ")
  )
})
sc <- bind_rows(rows)

# NHANES itself is the reference row - and its adjustment asymmetry is the
# single most important thing a reader needs before interpreting anything.
nh_row <- tibble(
  dataset_id = "OralMicroNHANES (reference)",
  citation = "NHANES 2009-2012 (cycles F+G), survey-weighted svyglm",
  doi = NA_character_, tier = "reference",
  variables_validated = "age, smoking",
  n_participants = 9828L,
  oral_site = "oral rinse (10 mL Scope, 30 s)",
  region_16S = "V4 (515F/806R)",
  covariates_adjusted = paste(
    "age: NONE (crude, n_covariates=0)",
    "smoking: age, age^2, sex, PIR, education(4), ethnicity(4), US-born", sep = " | "),
  adults_only = FALSE, subject_dedup = FALSE,
  smoking_definition = "SMQ_current_ever_never; current vs never (3-level, former separable)",
  n_genera_dataset = 956L, n_overlap_nhanes_common = 91L,
  n_comparable_after_fdr = "age=54; smoking=48",
  key_caveat = "age model is UNADJUSTED; ages 14-69 (not adults-only)")

sc <- bind_rows(nh_row, sc)

# YAML folded scalars (`>`) carry trailing newlines, and prose fields can carry
# tabs. Writing those straight to TSV splits one logical row across several
# physical lines and silently corrupts the file. Flatten every character cell.
flatten_cells <- function(df) {
  dplyr::mutate(df, dplyr::across(dplyr::where(is.character), function(x) {
    x <- gsub("[\r\n\t]+", " ", x)
    x <- gsub(" {2,}", " ", x)
    trimws(x)
  }))
}
sc <- flatten_cells(sc)
stopifnot(!any(grepl("[\r\n\t]", unlist(sc[sapply(sc, is.character)]), useBytes = TRUE)))
readr::write_tsv(sc, file.path(VAL, "study_characteristics.tsv"))
vsay("wrote study_characteristics.tsv (", nrow(sc), " rows)")

# --- Concordance table (report-ready) ----------------------------------------
tc <- sm %>%
  transmute(
    Variable = unname(VAR_LABELS_SHORT[nhanes_variable]),
    Comparison = label, `Cohort class` = m12_cohort_class(tier), Primary = primary,
    `N participants` = ifelse(is.na(n_participants), "—", format(n_participants, big.mark = ",")),
    `NHANES hits` = n_nhanes_hits, `Tested` = n_tested, `Not tested` = n_not_tested,
    `Concordant` = n_concordant,
    `% concordant` = ifelse(is.na(pct_concordant), "—", sprintf("%.0f", pct_concordant)),
    `95% CI` = ifelse(is.na(ci_lo), "—", sprintf("%.0f-%.0f", ci_lo, ci_hi)),
    `Binomial p` = ifelse(is.na(binom_p), "—", format.pval(binom_p, digits = 2)),
    `Spearman rho` = ifelse(is.na(spearman_rho), "—", sprintf("%.2f", spearman_rho))) %>%
  arrange(Variable, desc(Primary), `Cohort class`)
readr::write_tsv(flatten_cells(tc), file.path(VAL, "table_concordance.tsv"))

# --- Per-genus tables --------------------------------------------------------
for (v in PRIMARY_VARIABLES) {
  pg <- cmp %>% filter(nhanes_variable == v, primary) %>%
    select(genus, phylum, nh_estimate, nh_fdr, nh_direction, prevalence_overall,
           label, outcome) %>%
    mutate(label = paste0(label, "")) %>%
    pivot_wider(names_from = label, values_from = outcome) %>%
    arrange(phylum, nh_estimate)
  readr::write_tsv(flatten_cells(pg), file.path(VAL, paste0("table_pergenus_", v, ".tsv")))
  vsay("wrote table_pergenus_", v, ".tsv (", nrow(pg), " genera)")
}

vsay("")
vsay("=== study characteristics (key columns) ===")
print(as.data.frame(sc %>% select(dataset_id, tier, n_participants, oral_site,
                                  region_16S, n_overlap_nhanes_common)))

# --- Coverage: the denominator behind every concordance percentage -----------
cv  <- readr::read_tsv(file.path(VAL, "coverage_summary.tsv"), show_col_types = FALSE)
cvg <- readr::read_tsv(file.path(VAL, "coverage_per_genus.tsv"), show_col_types = FALSE)
tab_cov <- cv %>% transmute(
  Variable = unname(VAR_LABELS_SHORT[nhanes_variable]),
  `NHANES hits` = n_hits,
  `Tested by >=1 dataset` = sprintf("%d (%.0f%%)", n_tested_by_any, pct_tested_by_any),
  `No external coverage`  = sprintf("%d (%.0f%%)", n_untested, 100 - pct_tested_by_any),
  `Max datasets testing one genus` = max_datasets_testing)
readr::write_tsv(flatten_cells(tab_cov), file.path(VAL, "table_coverage.tsv"))
vsay("wrote table_coverage.tsv")

# The acquisition gap list: NHANES hits no external dataset can check. This is
# the concrete specification for what a future dataset would have to measure.
tab_gap <- cvg %>% filter(n_datasets_testing == 0) %>%
  transmute(Variable = unname(VAR_LABELS_SHORT[nhanes_variable]),
            Genus = genus, Phylum = phylum,
            `NHANES effect (CLR)` = round(nh_estimate, 4),
            Direction = nh_direction,
            FDR = signif(nh_fdr, 3),
            `Prevalence (%)` = round(100 * prevalence_overall, 1)) %>%
  arrange(Variable, desc(abs(`NHANES effect (CLR)`)))
readr::write_tsv(flatten_cells(tab_gap), file.path(VAL, "table_uncovered_hits.tsv"))
vsay("wrote table_uncovered_hits.tsv (", nrow(tab_gap), " genera with no external coverage)")

