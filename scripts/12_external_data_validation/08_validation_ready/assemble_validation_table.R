#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Assemble the validation-ready analysis table
#
# Purpose : Join the genus abundance table to harmonized phenotype via the
#           explicit ID crosswalk, quantify (never discard) unmatched rows, and
#           emit the single table a downstream analyst uses.
#
#           Nothing reaches 09_validation_ready/ unless microbiome and phenotype
#           rows are demonstrably linkable and no join loses rows silently.
#
# Usage   : Rscript assemble_validation_table.R <dataset_id>
# Inputs  : datasets/<id>/05_processed_microbiome/genus_table.tsv
#           datasets/<id>/08_sample_mapping/sample_crosswalk.tsv
#           config/<id>.yaml -> join.microbiome_key
# Outputs : datasets/<id>/09_validation_ready/analysis_table.tsv
#           datasets/<id>/09_validation_ready/analysis_table_README.md
#           datasets/<id>/09_validation_ready/join_diagnostics.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(yaml); library(tibble) })

# === USER CONFIG =============================================================
PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
# =============================================================================
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: assemble_validation_table.R <dataset_id>")
dataset_id <- args[1]

M12    <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
DS_DIR <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation",
                    "datasets", dataset_id)
OUT    <- file.path(DS_DIR, "09_validation_ready")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(grepl("12_external_data_validation", DS_DIR, fixed = TRUE))

source(file.path(M12, "00_shared", "harmonize_helpers.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a
cfg <- yaml::read_yaml(file.path(M12, "config", paste0(dataset_id, ".yaml")))
say <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")

say("=== assemble_validation_table: ", dataset_id, " ===")

gt <- readr::read_tsv(file.path(DS_DIR, "05_processed_microbiome", "genus_table.tsv"),
                      show_col_types = FALSE, progress = FALSE)
xw <- readr::read_tsv(file.path(DS_DIR, "08_sample_mapping", "sample_crosswalk.tsv"),
                      show_col_types = FALSE, progress = FALSE)

mkey <- cfg$join$microbiome_key %||% "run_accession"
say("microbiome rows: ", nrow(gt), " | crosswalk rows: ", nrow(xw), " | join key: ", mkey)
stopifnot(mkey %in% names(xw))

# Two crosswalk shapes exist in this module:
#   - repository datasets: the crosswalk already carries harmonized phenotype
#   - flat-file datasets:  the crosswalk carries only identifiers
# In the second case, attach sample_metadata.tsv here so the delivered table
# always contains phenotype regardless of which upstream builder produced it.
mf <- file.path(DS_DIR, "07_processed_metadata", "sample_metadata.tsv")
harmonized <- c("age_years", "sex_harmonized", "smoking_status_harmonized")
if (file.exists(mf) && !any(harmonized %in% names(xw))) {
  md <- readr::read_tsv(mf, show_col_types = FALSE, progress = FALSE)
  pkey <- cfg$join$phenotype_key %||% mkey
  if (pkey %in% names(md)) {
    dup_cols <- setdiff(intersect(names(md), names(xw)), pkey)
    md <- md %>% select(-any_of(dup_cols))
    n0 <- nrow(xw)
    xw <- xw %>% left_join(md, by = setNames(pkey, mkey))
    m12_assert_no_row_loss(n0, nrow(xw), "crosswalk <- sample_metadata")
    say("attached phenotype from sample_metadata.tsv (", ncol(md) - 1, " columns)")
  } else {
    say("[warn] phenotype key '", pkey, "' not in sample_metadata.tsv; no phenotype attached")
  }
}

genus_cols <- setdiff(names(gt), "sample_id")

# --- Join diagnostics: quantify, never silently drop -------------------------
mic_ids <- as.character(gt$sample_id)
xw_ids  <- as.character(xw[[mkey]])
matched      <- intersect(mic_ids, xw_ids)
mic_orphans  <- setdiff(mic_ids, xw_ids)
pheno_orphan <- setdiff(xw_ids, mic_ids)

say("matched: ", length(matched))
say("microbiome rows with no phenotype: ", length(mic_orphans))
say("phenotype rows with no microbiome: ", length(pheno_orphan))

diag <- tibble(
  metric = c("microbiome_rows", "crosswalk_rows", "matched",
             "microbiome_without_phenotype", "phenotype_without_microbiome",
             "duplicate_microbiome_ids", "duplicate_crosswalk_keys"),
  n = c(length(mic_ids), length(xw_ids), length(matched),
        length(mic_orphans), length(pheno_orphan),
        sum(duplicated(mic_ids)), sum(duplicated(xw_ids)))
)
readr::write_tsv(diag, file.path(OUT, "join_diagnostics.tsv"))

if (diag$n[diag$metric == "duplicate_microbiome_ids"] > 0 ||
    diag$n[diag$metric == "duplicate_crosswalk_keys"] > 0) {
  stop("Duplicate join keys present - a many-to-many join would silently ",
       "inflate rows. Resolve and document before proceeding.")
}
say("[OK] join keys are unique on both sides (no many-to-many risk)")

# --- Join --------------------------------------------------------------------
before <- nrow(gt)
tab <- gt %>%
  rename(!!mkey := sample_id) %>%
  inner_join(xw, by = mkey)
say("joined rows: ", nrow(tab), " (from ", before, " microbiome rows)")
if (nrow(tab) != length(matched))
  stop("Join produced ", nrow(tab), " rows but ", length(matched), " ids matched.")

# --- Supplementary phenotype sources ----------------------------------------
# Repository BioSample attributes are not always the only public phenotype. A
# code repo's own metadata table often covers subcohorts the BioSamples omit.
# Sources are declared in config so this is auditable, and any variable present
# in BOTH sources is cross-checked per sample rather than silently preferred.
for (sp in cfg$supplementary_phenotype %||% list()) {
  f <- file.path(DS_DIR, sp$file)
  if (!file.exists(f)) { say("[warn] supplementary phenotype missing: ", sp$file); next }
  say("")
  say("--- merging supplementary phenotype: ", sp$file, " ---")

  sup <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE,
                         na = c("", "NA", "na", "missing"))
  names(sup) <- sub("^﻿", "", names(sup))          # strip BOM if present
  stopifnot(sp$key %in% names(sup))

  keep <- unlist(sp$variables)
  keep <- keep[keep %in% names(sup)]
  sub_df <- sup %>% select(all_of(c(sp$key, unname(keep))))
  names(sub_df) <- c(sp$join_to, names(keep))

  n0 <- nrow(tab)
  tab <- tab %>% left_join(sub_df, by = sp$join_to)
  m12_assert_no_row_loss(n0, nrow(tab), paste("merge", sp$file))
  say("  merged ", length(keep), " variables; matched ",
      sum(!is.na(tab[[names(keep)[1]]])), "/", nrow(tab), " rows")

  # Cross-check overlapping variables. Disagreements are reported, never hidden.
  for (a in names(sp$crosscheck %||% list())) {
    b <- sp$crosscheck[[a]]
    if (!all(c(a, b) %in% names(tab))) next
    x <- tab[[a]]; y <- tab[[b]]
    both <- !is.na(x) & !is.na(y) & (!is.factor(y) | y != "unknown")
    if (!any(both)) { say("  crosscheck ", a, " vs ", b, ": no overlap"); next }
    if (is.numeric(x) || suppressWarnings(!all(is.na(as.numeric(as.character(x[both])))))) {
      xa <- suppressWarnings(as.numeric(as.character(x))); ya <- suppressWarnings(as.numeric(as.character(y)))
      dis <- sum(both & abs(xa - ya) > 0.5, na.rm = TRUE)
    } else {
      dis <- sum(both & tolower(as.character(x)) != tolower(as.character(y)), na.rm = TRUE)
    }
    say(sprintf("  crosscheck %-14s vs %-26s overlap=%d disagreements=%d",
                a, b, sum(both), dis))
    if (dis > 0)
      say("    [warn] ", dis, " disagreement(s) - recorded, NOT auto-resolved")
  }

  # Fill gaps in the primary harmonized fields from this source where the
  # primary is missing. Both originals are retained; nothing is overwritten.
  for (fill in sp$fill %||% list()) {
    tgt <- fill$target; src <- fill$source
    if (!all(c(tgt, src) %in% names(tab))) next
    before_n <- sum(!is.na(tab[[tgt]]) & (!is.factor(tab[[tgt]]) | tab[[tgt]] != "unknown"))
    if (identical(fill$type, "numeric")) {
      v <- suppressWarnings(as.numeric(as.character(tab[[src]])))
      tab[[tgt]] <- ifelse(is.na(tab[[tgt]]), v, tab[[tgt]])
    } else if (identical(fill$type, "sex")) {
      h <- m12_harmonize_sex(tab[[src]])
      tab[[tgt]] <- as.character(tab[[tgt]])
      tab[[tgt]] <- ifelse(is.na(tab[[tgt]]) | tab[[tgt]] == "unknown", as.character(h), tab[[tgt]])
    } else {
      tab[[tgt]] <- ifelse(is.na(tab[[tgt]]), as.character(tab[[src]]), as.character(tab[[tgt]]))
    }
    after_n <- sum(!is.na(tab[[tgt]]) & (!is.factor(tab[[tgt]]) | tab[[tgt]] != "unknown"))
    say(sprintf("  filled %-24s %d -> %d (+%d)", tgt, before_n, after_n, after_n - before_n))
  }
}

# Recompute the 5-year bin after any age fill.
if ("age_years" %in% names(tab)) tab$age_group_5yr <- m12_age_group_5yr(tab$age_years)

# Order columns: identifiers, phenotype, then genera.
genus_cols <- intersect(genus_cols, names(tab))
pheno_cols <- setdiff(names(tab), c(genus_cols, mkey))
tab <- tab %>% select(all_of(mkey), all_of(pheno_cols), all_of(genus_cols))

f_out <- file.path(OUT, "analysis_table.tsv")
readr::write_tsv(tab, f_out)
say("wrote ", f_out, " (", nrow(tab), " x ", ncol(tab), ")")

# --- Coverage of the validation variables ------------------------------------
cov <- function(v) {
  if (!v %in% names(tab)) return(NA_real_)
  x <- tab[[v]]
  if (is.character(x) || is.factor(x)) mean(!is.na(x) & x != "unknown") else mean(!is.na(x))
}
vars <- c("age_years", "sex_harmonized", "smoking_status_harmonized",
          "bmi", "ethnicity", "disease_status", "cotinine")
say("")
say("=== validation variable coverage in the joined table ===")
for (v in vars) {
  c_ <- cov(v)
  if (!is.na(c_)) say(sprintf("  %-28s %5.1f%% (n=%d)", v, 100*c_, round(c_*nrow(tab))))
}
if ("smoking_status_harmonized" %in% names(tab)) {
  say("  smoking levels: ",
      paste(sprintf("%s=%d", names(table(tab$smoking_status_harmonized)),
                    table(tab$smoking_status_harmonized)), collapse = ", "))
}
if ("subcohort" %in% names(tab)) {
  say("  by subcohort: ",
      paste(sprintf("%s=%d", names(table(tab$subcohort)), table(tab$subcohort)), collapse = ", "))
}

# --- README for the analyst ---------------------------------------------------
n_gen <- length(genus_cols)
readme <- c(
  paste0("# ", dataset_id, " — validation-ready analysis table"),
  "",
  paste0("`analysis_table.tsv`: **", nrow(tab), " samples x ", ncol(tab), " columns** ",
         "(", length(pheno_cols), " identifier/phenotype columns + ", n_gen, " genus abundance columns)."),
  "",
  "## Column layout",
  paste0("1. `", mkey, "` — join key (run accession)."),
  "2. Identifier and phenotype columns. Harmonized fields end in `_years`, `_harmonized`;",
  "   the corresponding `*_original` column is always retained beside them.",
  paste0("3. The final ", n_gen, " columns are genus-level **counts** (not relative abundance)."),
  "   Use `../05_processed_microbiome/genus_table_relative.tsv` for proportions.",
  "",
  "## Join provenance",
  paste0("- microbiome rows: ", length(mic_ids), "; crosswalk rows: ", length(xw_ids),
         "; matched: ", length(matched)),
  paste0("- unmatched microbiome rows: ", length(mic_orphans),
         "; unmatched phenotype rows: ", length(pheno_orphan),
         " (quantified in `join_diagnostics.tsv`, not discarded silently)"),
  "- Join keys verified unique on both sides; no many-to-many join was performed.",
  "",
  "## Read this before analysing",
  "See `../11_documentation/DATASET_REPORT.md` for the full caveat list. The three that",
  "most affect interpretation:",
  "",
  "1. **Subcohorts must not be naively pooled.** Sequencing platform is confounded with",
  "   cohort (454/V3-V4 for ACS+PLCO vs MiSeq/V4 for NYC HANES).",
  "2. **A 43% pre-diagnostic cancer stratum exists in ACS/PLCO** and is not mentioned in the",
  "   source publication. Stratify or restrict on `disease_status`.",
  "3. **The feature space is a pruned cross-study intersection.** A genus absent here was",
  "   filtered upstream; it is not evidence of absence from saliva.",
  ""
)
writeLines(readme, file.path(OUT, "analysis_table_README.md"))
say("wrote analysis_table_README.md")
