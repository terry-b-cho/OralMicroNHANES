#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Cross-dataset genus overlap against the NHANES cycle F/G reference
#
# Purpose : Answer, per external dataset: which NHANES genera does it contain,
#           and how much of the NHANES oral community does that actually cover?
#
#           Three things this does that a raw intersection count does not:
#
#           1. USES AN HONEST DENOMINATOR. NHANES catalogues 956 genera but only
#              900 are detected, and 735 of those are rare (<1% of participants).
#              Reporting "40/956" is misleading in both directions. Overlap is
#              reported against detected, core (>=50%) and common (>=10%) sets.
#           2. WEIGHTS BY ABUNDANCE. Covering 40 genera means little if they are
#              rare; it means a lot if they are 80% of the community. The
#              abundance-weighted coverage says which.
#           3. NORMALIZES NOMENCLATURE. NHANES is SILVA (Prevotella_7); external
#              sets are eHOMD/Greengenes (Prevotella). Reported strictly AND
#              normalized so the effect of that choice is visible.
#
# Usage   : Rscript 09_reports/build_genus_overlap_matrix.R
# Inputs  : nhanes_genus_reference.tsv (build_nhanes_genus_reference.R)
#           datasets/<id>/05_processed_microbiome/genus_table.tsv         (abundance sets)
#           datasets/<id>/.../replication_*.tsv                           (target sets)
# Outputs : genus_overlap_summary.tsv          one row per external dataset
#           genus_presence_matrix.tsv          NHANES genus x dataset presence
#           genus_overlap_report.md            human-readable
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(purrr)
})
PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
set.seed(42)
`%||%` <- function(a, b) if (is.null(a)) b else a

M12 <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
OUT <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation")
DSD <- file.path(OUT, "datasets")
source(file.path(M12, "00_shared", "harmonize_helpers.R"))
m12_assert_parse_genus_agrees(verbose = FALSE)
say <- function(...) cat(..., "\n", sep = "")

ref <- m12_nhanes_genus_reference(PROJECT_ROOT) %>%
  filter(genus != "unclassified") %>%
  mutate(key = m12_genus_crossdb_key(genus))

detected <- ref %>% filter(detected_either)
core     <- detected %>% filter(prevalence_overall >= 0.50)
common   <- detected %>% filter(prevalence_overall >= 0.10)
say("NHANES reference: ", nrow(detected), " detected genera (",
    nrow(core), " core >=50%, ", nrow(common), " common >=10%)")

# Total NHANES community mass carried by the detected genera - the denominator
# for abundance-weighted coverage.
total_abund <- sum(detected$mean_abund_overall)

# --- Collect the genus set contributed by each external dataset --------------
# Two kinds of contribution, kept distinct because they mean different things:
#   abundance : we hold an actual genus x sample table (can be re-analysed)
#   target    : we hold a published list of associated genera (comparison only)
sets <- list()
add_set <- function(id, kind, genera, note = "") {
  g <- unique(genera[!is.na(genera) & nzchar(genera) & genera != "unclassified"])
  if (!length(g)) return(invisible(NULL))
  sets[[length(sets)+1]] <<- tibble(dataset_id = id, kind = kind, genus_external = g, note = note)
}

for (d in list.dirs(DSD, recursive = FALSE)) {
  id <- basename(d)
  gt <- file.path(d, "05_processed_microbiome", "genus_table.tsv")
  if (file.exists(gt)) {
    hdr <- names(suppressWarnings(readr::read_tsv(gt, n_max = 0, show_col_types = FALSE)))
    add_set(id, "abundance", setdiff(hdr, c("sample_id", "run_accession")), basename(gt))
  }
  for (f in c(Sys.glob(file.path(d, "02_supplementary", "derived", "replication_*.tsv")),
              Sys.glob(file.path(d, "09_validation_ready", "replication_*.tsv")))) {
    t <- suppressWarnings(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE))
    col <- intersect(c("genus", "taxon"), names(t))
    if (!length(col)) next

    # A target set must be the SELECTED/significant genera, not the whole feature
    # space. Huang's RF-importance file scores all 1,070 ASVs; counting it as a
    # target would imply 141 age-associated genera when it is simply every taxon
    # measured, and would duplicate the abundance row exactly.
    if (grepl("feature_importance|feature_imps", basename(f))) {
      say("  [skip] ", basename(f), " - full feature space, not a selected result")
      next
    }
    if ("q_value" %in% names(t)) {
      q <- suppressWarnings(as.numeric(t$q_value))
      if (any(!is.na(q))) {
        keep <- !is.na(q) & q < 0.05
        say("  [filter] ", basename(f), ": ", sum(keep), "/", nrow(t), " rows at q<0.05")
        t <- t[keep, , drop = FALSE]
      }
    }
    if (!nrow(t)) next
    g <- as.character(t[[col[1]]])
    if ("rank" %in% names(t)) g <- g[t$rank %in% c("genus", NA)]
    add_set(id, "target", vapply(g[!is.na(g)], parse_genus, character(1), USE.NAMES = FALSE),
            basename(f))
  }
}
ext <- bind_rows(sets)
stopifnot(nrow(ext) > 0)

# --- Overlap per dataset -----------------------------------------------------
summ <- list(); pres <- list()
for (id in unique(ext$dataset_id)) {
  for (k in unique(ext$kind[ext$dataset_id == id])) {
    e <- ext %>% filter(dataset_id == id, kind == k)
    ekey <- unique(m12_genus_crossdb_key(e$genus_external))

    strict <- intersect(e$genus_external, detected$genus)
    norm   <- detected %>% filter(key %in% ekey)

    cov_ab <- sum(norm$mean_abund_overall) / total_abund
    summ[[length(summ)+1]] <- tibble(
      dataset_id = id, kind = k,
      n_external_genera = n_distinct(e$genus_external),
      n_match_strict = length(strict),
      n_match_normalized = nrow(norm),
      gained_by_normalization = nrow(norm) - length(strict),
      # An external genus can map to SEVERAL NHANES genera (external
      # "Prevotella" -> NHANES Prevotella, _2, _6, _7, _9). Matched can
      # therefore exceed the external count; record the expansion explicitly.
      n_external_keys_matched = length(intersect(ekey, detected$key)),
      one_to_many_expansion = nrow(norm) - length(intersect(ekey, detected$key)),
      pct_of_external_matched = round(100 * nrow(norm) / n_distinct(e$genus_external), 1),
      n_nhanes_detected = nrow(detected),
      pct_of_nhanes_detected = round(100 * nrow(norm) / nrow(detected), 1),
      n_core_covered = sum(norm$prevalence_overall >= 0.50),
      pct_core_covered = round(100 * sum(norm$prevalence_overall >= 0.50) / nrow(core), 1),
      n_common_covered = sum(norm$prevalence_overall >= 0.10),
      pct_common_covered = round(100 * sum(norm$prevalence_overall >= 0.10) / nrow(common), 1),
      pct_abundance_covered = round(100 * cov_ab, 1),
      n_covered_in_F = sum(norm$detected_F), n_covered_in_G = sum(norm$detected_G),
      n_covered_in_both = sum(norm$detected_both))

    if (k == "abundance" || !any(ext$kind[ext$dataset_id == id] == "abundance"))
      pres[[length(pres)+1]] <- tibble(dataset_id = id, genus = norm$genus, present = TRUE)
  }
}
s <- bind_rows(summ) %>% arrange(desc(pct_abundance_covered))
readr::write_tsv(s, file.path(OUT, "genus_overlap_summary.tsv"))

# --- NHANES genus x dataset presence matrix ----------------------------------
pm <- bind_rows(pres) %>% distinct() %>%
  pivot_wider(names_from = dataset_id, values_from = present, values_fill = FALSE)
mat <- detected %>%
  select(genus, prevalence_overall, prevalence_F, prevalence_G, mean_abund_overall,
         detected_F, detected_G, detected_both, prevalence_tier) %>%
  left_join(pm, by = "genus") %>%
  mutate(across(where(is.logical), ~ tidyr::replace_na(.x, FALSE)))
dcols <- setdiff(names(mat), c("genus", "prevalence_overall", "prevalence_F", "prevalence_G",
                               "mean_abund_overall", "detected_F", "detected_G",
                               "detected_both", "prevalence_tier"))
mat$n_datasets_covering <- rowSums(mat[dcols])
mat <- mat %>% arrange(desc(prevalence_overall))
readr::write_tsv(mat, file.path(OUT, "genus_presence_matrix.tsv"))

say("")
say("=== overlap per dataset (normalized nomenclature) ===")
print(as.data.frame(s %>% select(dataset_id, kind, n_external_genera, n_match_normalized,
                                 pct_core_covered, pct_common_covered, pct_abundance_covered)))

# --- Report ------------------------------------------------------------------
L <- c("# Module 12 — NHANES genus coverage by external dataset", "",
  paste0("_Generated ", format(Sys.time(), tz="UTC", "%Y-%m-%d %H:%M UTC"),
         " by `09_reports/build_genus_overlap_matrix.R`._"), "",
  "## The denominator problem", "",
  sprintf("NHANES catalogues **956** genus names, but only **%d are actually detected** in cycle F or G.",
          nrow(detected)),
  sprintf("Of those, **%d are rare** (<1%% of participants). Only **%d are common (>=10%%)** and **%d are core (>=50%%)**.",
          sum(detected$prevalence_overall < 0.01), nrow(common), nrow(core)),
  "",
  "Quoting overlap against 956 therefore flatters nothing and hides everything. This report uses",
  "three honest denominators plus an abundance weighting.", "",
  "| Definition | N |", "|---|---|",
  sprintf("| catalogued in the OTU map | 956 |"),
  sprintf("| detected in cycle F | %d |", sum(detected$detected_F)),
  sprintf("| detected in cycle G | %d |", sum(detected$detected_G)),
  sprintf("| detected in either cycle | %d |", nrow(detected)),
  sprintf("| detected in **both** cycles | %d |", sum(detected$detected_both)),
  sprintf("| common (>=10%% prevalence) | %d |", nrow(common)),
  sprintf("| core (>=50%% prevalence) | %d |", nrow(core)),
  "",
  "> **The F/G asymmetry is real and large.** Cycle F detects %d genera, cycle G only %d, with",
  "> just %d shared. Any genus-level claim should say which cycle it rests on.", "",
  "## Coverage by dataset", "",
  "`kind = abundance` means we hold a genus x sample table and could re-analyse it.",
  "`kind = target` means we hold a published list of associated genera for comparison only.", "",
  "`Matched` counts NHANES genera, and one external genus can match several",
  "(external *Prevotella* matches NHANES *Prevotella*, *_2*, *_6*, *_7*, *_9*), so it can",
  "exceed the external count.", "",
  "| Dataset | Kind | Ext. genera | Matched | % of NHANES core | % of common | **% of abundance** |",
  "|---|---|---|---|---|---|---|")
L <- sub("Cycle F detects %d genera, cycle G only %d, with",
         sprintf("Cycle F detects %d genera, cycle G only %d, with", sum(detected$detected_F), sum(detected$detected_G)), L)
L <- sub("> just %d shared.", sprintf("> just %d shared.", sum(detected$detected_both)), L)
for (i in seq_len(nrow(s)))
  L <- c(L, sprintf("| `%s` | %s | %d | %d | %.1f%% | %.1f%% | **%.1f%%** |",
                    s$dataset_id[i], s$kind[i], s$n_external_genera[i], s$n_match_normalized[i],
                    s$pct_core_covered[i], s$pct_common_covered[i], s$pct_abundance_covered[i]))
gain <- s %>% filter(gained_by_normalization > 0)
if (nrow(gain)) {
  L <- c(L, "", "### Genera recovered by nomenclature normalization", "",
    "NHANES is SILVA (`Prevotella_7`); external sets are eHOMD/Greengenes (`Prevotella`).",
    "Strict string matching silently drops highly prevalent NHANES genera.", "",
    "| Dataset | Kind | Strict | Normalized | Gained |", "|---|---|---|---|---|")
  for (i in seq_len(nrow(gain)))
    L <- c(L, sprintf("| `%s` | %s | %d | %d | **+%d** |", gain$dataset_id[i], gain$kind[i],
                      gain$n_match_strict[i], gain$n_match_normalized[i], gain$gained_by_normalization[i]))
}
topu <- mat %>% filter(n_datasets_covering == 0, prevalence_overall >= 0.50) %>% head(20)
if (nrow(topu)) {
  L <- c(L, "", "### Core NHANES genera covered by NO external dataset", "",
         "These are the blind spots — prevalent in NHANES, absent from every external set:", "")
  L <- c(L, paste0("- *", topu$genus, "* (", round(100*topu$prevalence_overall, 1), "% of participants)"))
}
L <- c(L, "", "Full detail: `genus_overlap_summary.tsv`, `genus_presence_matrix.tsv`,",
       "`nhanes_genus_reference.tsv`.")
# The narrative .md is a PRIVATE document (team reading), while the TSVs above
# are public. It is written only when the private driver asks for it by setting
# M12_NARRATIVE_MD_DIR; a public run emits no narrative file at all, so the
# public results tree never contains one. Same source of truth either way.
m12_narrative_dir <- Sys.getenv("M12_NARRATIVE_MD_DIR", "")
if (nzchar(m12_narrative_dir)) {
  writeLines(L, file.path(m12_narrative_dir, "genus_overlap_report.md"))
  say("wrote genus_overlap_report.md -> ", m12_narrative_dir)
}
say("")
say("wrote genus_overlap_summary.tsv, genus_presence_matrix.tsv")
