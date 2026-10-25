#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Extract published per-taxon association results as REPLICATION TARGETS
#
# Purpose : Some external studies publish taxon-level association statistics we
#           can compare our NHANES results against directly, WITHOUT needing
#           their individual-level data. Those are strictly more valuable than
#           a raw abundance table for validation, because they encode the
#           original authors' own analysis.
#
#           This script normalizes such a table to a common schema, rolls it up
#           to GENUS using the same parser as NHANES, and reports overlap.
#
# Usage   : Rscript build_replication_targets.R <dataset_id>
# Inputs  : config/<id>.yaml -> replication_targets[]
# Outputs : datasets/<id>/09_validation_ready/replication_target_<name>.tsv
#           datasets/<id>/09_validation_ready/replication_targets_summary.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(yaml); library(tibble); library(stringr)
})

PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
set.seed(42)
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: build_replication_targets.R <dataset_id>")
dataset_id <- args[1]

M12    <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
DS_DIR <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation",
                    "datasets", dataset_id)
OUT    <- file.path(DS_DIR, "09_validation_ready")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(grepl("12_external_data_validation", DS_DIR, fixed = TRUE))

source(file.path(M12, "00_shared", "harmonize_helpers.R"))
m12_assert_parse_genus_agrees()
cfg <- yaml::read_yaml(file.path(M12, "config", paste0(dataset_id, ".yaml")))
say <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")

targets <- cfg$replication_targets %||% list()
if (!length(targets)) { say("no replication_targets configured for ", dataset_id); quit(status = 0) }

say("=== build_replication_targets: ", dataset_id, " ===")
summ <- list()

for (tg in targets) {
  f <- file.path(DS_DIR, tg$file)
  if (!file.exists(f)) { say("[warn] missing: ", tg$file); next }
  say("")
  say("--- ", tg$name, " (", tg$file, ") ---")

  d <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE, guess_max = 10000)
  say("  ", nrow(d), " rows x ", ncol(d), " cols")

  cm <- tg$columns
  need <- cm$lineage
  if (!need %in% names(d)) { say("[warn] lineage column '", need, "' absent; skipping"); next }

  res <- tibble(
    dataset_id  = dataset_id,
    target_name = tg$name,
    lineage     = as.character(d[[cm$lineage]]),
    stratum     = if (!is.null(cm$stratum) && cm$stratum %in% names(d)) as.character(d[[cm$stratum]]) else NA_character_,
    effect      = if (!is.null(cm$effect)  && cm$effect  %in% names(d)) suppressWarnings(as.numeric(d[[cm$effect]]))  else NA_real_,
    p_value     = if (!is.null(cm$p_value) && cm$p_value %in% names(d)) suppressWarnings(as.numeric(d[[cm$p_value]])) else NA_real_,
    q_value     = if (!is.null(cm$q_value) && cm$q_value %in% names(d)) suppressWarnings(as.numeric(d[[cm$q_value]])) else NA_real_,
    direction   = if (!is.null(cm$direction) && cm$direction %in% names(d)) as.character(d[[cm$direction]]) else NA_character_,
    significant = if (!is.null(cm$significant) && cm$significant %in% names(d)) as.character(d[[cm$significant]]) else NA_character_
  )

  # Genus via the SAME parser NHANES used, so the vocabularies are comparable.
  g <- m12_genus_from_lineage(res$lineage, style = tg$lineage_style %||% "qiime")
  g[is.na(g) | !nzchar(g)] <- "unclassified"
  res$genus <- vapply(g, parse_genus, character(1), USE.NAMES = FALSE)

  n_unc <- sum(res$genus == "unclassified")
  say(sprintf("  genus assigned for %d/%d rows (%d unclassified)",
              nrow(res) - n_unc, nrow(res), n_unc))

  # Significance filter, if the source provides one
  sig <- res
  if (!is.na(tg$significance_rule %||% NA)) {
    rule <- tg$significance_rule
    if (identical(rule, "q_value") && any(!is.na(res$q_value))) {
      thr <- tg$q_threshold %||% 0.05
      sig <- res %>% filter(!is.na(q_value), q_value < thr)
      say("  significant at q < ", thr, ": ", nrow(sig), " rows")
    } else if (identical(rule, "flag") && any(!is.na(res$significant))) {
      sig <- res %>% filter(tolower(significant) %in% c("true", "yes", "sig", "1"))
      say("  significant by flag: ", nrow(sig), " rows")
    }
  }

  readr::write_tsv(res, file.path(OUT, paste0("replication_target_", tg$name, ".tsv")))

  # Genus-level overlap with the NHANES universe
  gen_all <- setdiff(unique(res$genus), "unclassified")
  gen_sig <- setdiff(unique(sig$genus), "unclassified")
  ov_all <- m12_genus_overlap(gen_all)
  ov_sig <- m12_genus_overlap(gen_sig)
  say(sprintf("  genera: %d total (%d in NHANES, %.1f%%) | %d significant (%d in NHANES, %.1f%%)",
              ov_all$n_external, ov_all$n_shared, ov_all$pct_of_external,
              ov_sig$n_external, ov_sig$n_shared, ov_sig$pct_of_external))
  if (length(ov_sig$shared))
    say("  significant genera present in NHANES: ",
        paste(head(sort(ov_sig$shared), 30), collapse = ", "),
        if (length(ov_sig$shared) > 30) " ..." else "")

  summ[[length(summ)+1]] <- tibble(
    dataset_id = dataset_id, target_name = tg$name, file = tg$file,
    n_rows = nrow(res), n_genera = ov_all$n_external,
    n_genera_in_nhanes = ov_all$n_shared,
    n_significant_rows = nrow(sig), n_significant_genera = ov_sig$n_external,
    n_significant_genera_in_nhanes = ov_sig$n_shared,
    has_effect_size = any(!is.na(res$effect)),
    has_direction = any(!is.na(res$direction)),
    has_qvalue = any(!is.na(res$q_value)),
    comparison_mode = tg$comparison_mode %||% "unspecified",
    caveat = tg$caveat %||% NA_character_)
}

if (length(summ)) {
  s <- bind_rows(summ)
  readr::write_tsv(s, file.path(OUT, "replication_targets_summary.tsv"))
  say("")
  say("wrote replication_targets_summary.tsv (", nrow(s), " targets)")
  print(as.data.frame(s %>% select(target_name, n_rows, n_genera, n_genera_in_nhanes,
                                   n_significant_genera_in_nhanes, has_effect_size, has_direction)))
}
