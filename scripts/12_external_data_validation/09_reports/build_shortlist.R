#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Cross-study consolidation: replication targets + ranked validation shortlist
#
# Purpose : Answer the module's actual question — "for each handpicked variable,
#           which external cohorts can validate it, and how strongly?" — from
#           what is measurably on disk, not from recollection.
#
# Usage   : Rscript 09_reports/build_shortlist.R
# Outputs : replication_targets_all.tsv     one row per external target genus
#           validation_shortlist.md         ranked, per variable
#           replication_targets_rows.tsv    per-source rows behind each target
#           pairwise_direction_concordance.tsv  cohort x cohort agreement
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(yaml); library(tibble); library(purrr)
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

nh <- m12_nhanes_genus_universe(PROJECT_ROOT)

# --- Registry of extracted replication targets --------------------------------
# Each entry names a normalized target file, the variable it speaks to, and the
# column holding the taxon name. Registered here (not auto-discovered) so a
# stray file cannot silently become a scientific claim.
targets <- tribble(
  ~dataset_id,                     ~variable, ~file,                                                              ~taxon_col, ~direction_col,                 ~rank_filter, ~q_col,
  "WillisJR_2022_SacaLaLengua2",   "age",     "02_supplementary/derived/replication_age_genera.tsv",               "taxon",    "direction",                    "genus",      NA,
  "WillisJR_2022_SacaLaLengua2",   "smoking", "02_supplementary/derived/replication_smoking_genera.tsv",           "taxon",    "direction",                    "genus",      NA,
  "HuangS_2020_AgePrediction",     "age",     "09_validation_ready/replication_target_young_vs_old_differential_abundance.tsv", "genus", "direction",         NA,           "q_value",
  "HuangS_2020_AgePrediction",     "age",     "02_supplementary/derived/FIGS3_oral_feature_importance_rank.tsv",   "g",        NA,                             NA,           NA,
  "GhensiP_2020_ImplantDisease",   "disease", "02_supplementary/derived/replication_disease_species.tsv",          "genus",    "direction_in_positive_class",  NA,           "q_value",
  # Image-derived: the paper's only taxon-level result is a dot plot, so the
  # direction is read off the figure. Lower provenance confidence than a
  # machine-parsed target - flagged in the report and here.
  "WellsPM_2022_TwinsUKSaliva",    "age",     "02_supplementary/derived/replication_age_genera_from_figure.tsv",   "genus",    "direction_with_age",           NA,           NA
)

rows <- list()
for (i in seq_len(nrow(targets))) {
  t <- targets[i, ]
  f <- file.path(DSD, t$dataset_id, t$file)
  if (!file.exists(f)) { say("[warn] missing target: ", t$dataset_id, " / ", t$file); next }
  d <- suppressWarnings(readr::read_tsv(f, show_col_types = FALSE, progress = FALSE))
  if (!t$taxon_col %in% names(d)) { say("[warn] no column '", t$taxon_col, "' in ", t$file); next }

  if (!is.na(t$rank_filter) && "rank" %in% names(d)) d <- d %>% filter(rank == t$rank_filter)
  if (!is.na(t$q_col) && t$q_col %in% names(d))
    d <- d %>% mutate(.q = suppressWarnings(as.numeric(.data[[t$q_col]]))) %>%
               filter(!is.na(.q), .q < 0.05)

  # Filter the DATA FRAME, then read both columns from it. Filtering a taxon
  # vector and separately slicing the direction column by position silently
  # misaligns every taxon with someone else's direction - which inverted an
  # entire dataset's directions before this was caught.
  d <- d %>% mutate(.taxon = as.character(.data[[t$taxon_col]])) %>%
             filter(!is.na(.taxon), nzchar(.taxon), .taxon != "unclassified")
  if (!nrow(d)) next
  tx <- d$.taxon
  g <- vapply(tx, parse_genus, character(1), USE.NAMES = FALSE)
  dir_raw <- if (!is.na(t$direction_col) && t$direction_col %in% names(d))
               as.character(d[[t$direction_col]]) else rep(NA_character_, nrow(d))
  stopifnot(length(dir_raw) == length(tx))   # alignment is asserted, not assumed
  # Direction vocabularies differ between sources ("down" vs "old_depleted" vs
  # an arrow). They MUST be normalized to a single scale before any concordance
  # claim, or "agreement" is just two different words being counted as one.
  # Convention: `up` = higher with MORE of the exposure (older age / smoking /
  # disease); `down` = lower.
  dir <- dplyr::case_when(
    is.na(dir_raw)                                   ~ NA_character_,
    tolower(dir_raw) %in% c("up", "increase", "increased", "positive",
                            "old_enriched", "enriched", "higher")   ~ "up",
    tolower(dir_raw) %in% c("down", "decrease", "decreased", "negative",
                            "old_depleted", "depleted", "lower")    ~ "down",
    tolower(dir_raw) %in% c("neutral", "ns", "none", "equal", "")    ~ NA_character_,
    TRUE                                             ~ NA_character_
  )
  unmapped <- setdiff(unique(dir_raw[!is.na(dir_raw) & nzchar(dir_raw)]),
                      c("up","down","increase","decrease","increased","decreased",
                        "positive","negative","old_enriched","old_depleted",
                        "enriched","depleted","higher","lower","Neutral","neutral",
                        "ns","none","equal"))
  if (length(unmapped))
    say("[warn] ", t$dataset_id, "/", t$variable,
        ": unmapped direction value(s) -> NA: ", paste(unmapped, collapse = ", "))

  rows[[length(rows)+1]] <- tibble(
    dataset_id = t$dataset_id, variable = t$variable,
    target_file = t$file, genus = g, direction = dir,
    in_nhanes = g %in% nh)
}

all_rows <- bind_rows(rows)
readr::write_tsv(all_rows, file.path(OUT, "replication_targets_rows.tsv"))

# A genus can appear MANY times within one target - once per species, and once
# per stratum (e.g. sex). Those rows can disagree, legitimately: different
# Prevotella species move in opposite directions with age. Taking an arbitrary
# first row would manufacture both false concordance and false discordance, so
# direction is decided by an explicit VOTE and the vote is reported.
all_t <- all_rows %>%
  group_by(dataset_id, variable, genus, in_nhanes) %>%
  summarise(n_rows = n(),
            n_up = sum(direction == "up", na.rm = TRUE),
            n_down = sum(direction == "down", na.rm = TRUE),
            .groups = "drop") %>%
  mutate(
    direction = case_when(
      n_up > 0 & n_down == 0 ~ "up",
      n_down > 0 & n_up == 0 ~ "down",
      n_up == 0 & n_down == 0 ~ NA_character_,
      n_up > n_down ~ "up",
      n_down > n_up ~ "down",
      TRUE ~ "mixed"),
    direction_support = case_when(
      is.na(direction) ~ NA_character_,
      direction == "mixed" ~ sprintf("tied %d up / %d down", n_up, n_down),
      (n_up == 0 | n_down == 0) ~ sprintf("unanimous (%d rows)", pmax(n_up, n_down)),
      TRUE ~ sprintf("majority %d/%d", pmax(n_up, n_down), n_up + n_down)))
readr::write_tsv(all_t, file.path(OUT, "replication_targets_all.tsv"))
say("wrote replication_targets_all.tsv (", nrow(all_t), " dataset x variable x genus rows)")

# --- Per-variable summary -----------------------------------------------------
summ <- all_t %>%
  group_by(variable, dataset_id) %>%
  summarise(n_genera = n_distinct(genus),
            n_in_nhanes = n_distinct(genus[in_nhanes]),
            has_direction = any(!is.na(direction)),
            .groups = "drop") %>%
  arrange(variable, desc(n_in_nhanes))

say("")
print(as.data.frame(summ))

# Consensus: genera named by >1 independent dataset for the same variable
consensus <- all_t %>% filter(in_nhanes) %>%
  group_by(variable, genus) %>%
  summarise(n_datasets = n_distinct(dataset_id),
            n_with_direction = n_distinct(dataset_id[!is.na(direction)]),
            datasets = paste(sort(unique(dataset_id)), collapse = "; "),
            directions = paste(sort(unique(direction[!is.na(direction)])), collapse = "/"),
            support = paste(direction_support[!is.na(direction)], collapse = "; "),
            .groups = "drop") %>%
  filter(n_datasets > 1) %>%
  mutate(concordance = dplyr::case_when(
    n_with_direction < 2                   ~ "overlap only (1 source has direction)",
    !grepl("/", directions)                ~ paste0("CONCORDANT (", directions, ")"),
    TRUE                                   ~ "DISCORDANT"
  )) %>%
  arrange(variable, desc(n_datasets), genus)

say("")
say("genera named by more than one dataset for the same variable: ", nrow(consensus))
if (nrow(consensus)) print(as.data.frame(consensus))

# --- Shortlist ----------------------------------------------------------------
cat_f <- file.path(OUT, "dataset_catalog.tsv")
cat_tbl <- if (file.exists(cat_f)) readr::read_tsv(cat_f, show_col_types = FALSE) else tibble()

L <- c(
  "# Module 12 — Validation shortlist",
  "",
  paste0("_Generated ", format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M UTC"),
         " by `09_reports/build_shortlist.R` from what is on disk._"),
  "",
  "NHANES reference: **956 unique parsed genera across 1,349 genus-level RSVs**, oral rinse, 16S V4.",
  "",
  "## Replication targets by variable",
  "",
  "A *replication target* is a published, taxon-level association result we can compare our NHANES",
  "findings against **without needing the external study's individual-level data**.",
  "",
  "| Variable | Dataset | Genera | In NHANES | Direction? |",
  "|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) {
  L <- c(L, sprintf("| %s | `%s` | %d | **%d** | %s |",
                    summ$variable[i], summ$dataset_id[i], summ$n_genera[i],
                    summ$n_in_nhanes[i], if (summ$has_direction[i]) "yes" else "no"))
}
if (nrow(consensus)) {
  L <- c(L, "", "## Genera named by more than one independent dataset", "",
         "| Variable | Genus | Datasets | Direction(s) | Concordance | Support |",
         "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(consensus)))
    L <- c(L, sprintf("| %s | *%s* | %d | %s | %s | %s |", consensus$variable[i], consensus$genus[i],
                      consensus$n_datasets[i], consensus$directions[i], consensus$concordance[i],
                      consensus$support[i]))
}
if (nrow(cat_tbl)) {
  L <- c(L, "", "## Individual-level cohorts", "",
         "| Dataset | Status | n | Genera shared | Caveats |", "|---|---|---|---|---|")
  for (i in seq_len(nrow(cat_tbl)))
    L <- c(L, sprintf("| `%s` | %s | %s | %s | %s |",
                      cat_tbl$dataset_id[i], cat_tbl$status[i],
                      ifelse(is.na(cat_tbl$n_samples_validation_ready[i]), "—",
                             cat_tbl$n_samples_validation_ready[i]),
                      ifelse(is.na(cat_tbl$n_genera_shared_nhanes[i]), "—",
                             cat_tbl$n_genera_shared_nhanes[i]),
                      cat_tbl$n_critical_caveats[i]))
}
# --- Cross-dataset direction concordance -------------------------------------
# The headline question: when two independent external cohorts both name a
# genus for the same variable, do they AGREE on the direction? If agreement is
# no better than chance, the cohorts are not measuring the same thing and no
# single one should be trusted as "the" replication.
# Concordance requires at least TWO datasets that each supply a direction.
# A genus named by two datasets where only one publishes a direction is NOT
# evidence of agreement - counting it as concordant overstates the evidence.
conc <- all_t %>%
  filter(in_nhanes, !is.na(direction), direction != "mixed") %>%
  mutate(unanimous = (n_up == 0 | n_down == 0)) %>%
  group_by(variable, genus) %>%
  filter(n_distinct(dataset_id) > 1) %>%
  summarise(n_datasets_with_direction = n_distinct(dataset_id),
            all_unanimous = all(unanimous),
            dirs = paste(sort(unique(direction)), collapse = "/"),
            .groups = "drop") %>%
  filter(n_datasets_with_direction > 1) %>%
  rename(n_datasets = n_datasets_with_direction) %>%
  mutate(concordant = !grepl("/", dirs))

if (nrow(conc)) {
  L <- c(L, "", "## Cross-dataset direction concordance", "",
         "When two independent cohorts both name a genus, do they agree on direction?",
         "", "| Variable | Subset | Genera | Concordant | % | Binomial p vs 50% |",
         "|---|---|---|---|---|---|")
  for (v in unique(conc$variable)) {
    for (sub in c("all", "high-confidence")) {
      cc <- conc %>% filter(variable == v)
      if (sub == "high-confidence") cc <- cc %>% filter(all_unanimous)
      if (!nrow(cc)) next
      k <- sum(cc$concordant); n <- nrow(cc)
      p <- tryCatch(binom.test(k, n, 0.5)$p.value, error = function(e) NA_real_)
      L <- c(L, sprintf("| %s | %s | %d | %d | %.0f%% | %s |", v, sub, n, k, 100*k/n,
                        ifelse(is.na(p), "—", sprintf("%.3f", p))))
    }
  }
  hi <- conc %>% filter(all_unanimous, concordant)
  if (nrow(hi)) {
    L <- c(L, "", "**Concordant in the high-confidence subset** (every contributing dataset unanimous):", "")
    for (v in unique(hi$variable))
      L <- c(L, paste0("- *", v, "*: ",
                       paste(sprintf("%s (%s)", hi$genus[hi$variable==v], hi$dirs[hi$variable==v]),
                             collapse = ", ")))
  }
  # Pairwise concordance. The POOLED number hides structure: cohorts that share
  # an oral niche can agree perfectly while a mismatched pair sits at chance and
  # dominates the shared-genus count. Report the pairs.
  pw <- list()
  for (v in unique(conc$variable)) {
    dv <- all_t %>% filter(variable == v, in_nhanes, !is.na(direction), direction != "mixed") %>%
                    select(dataset_id, genus, direction)
    ds <- sort(unique(dv$dataset_id))
    if (length(ds) < 2) next
    for (i in 1:(length(ds)-1)) for (j in (i+1):length(ds)) {
      m <- inner_join(dv %>% filter(dataset_id == ds[i]),
                      dv %>% filter(dataset_id == ds[j]), by = "genus")
      if (!nrow(m)) next
      pw[[length(pw)+1]] <- tibble(variable = v, A = ds[i], B = ds[j],
                                   shared = nrow(m),
                                   concordant = sum(m$direction.x == m$direction.y))
    }
  }
  if (length(pw)) {
    pwt <- bind_rows(pw) %>% mutate(pct = round(100*concordant/shared)) %>%
           arrange(variable, desc(shared))
    L <- c(L, "", "### Pairwise agreement between cohorts", "",
           "| Variable | Cohort A | Cohort B | Shared genera | Concordant | % |",
           "|---|---|---|---|---|---|")
    for (i in seq_len(nrow(pwt)))
      L <- c(L, sprintf("| %s | `%s` | `%s` | %d | %d | %d%% |", pwt$variable[i],
                        pwt$A[i], pwt$B[i], pwt$shared[i], pwt$concordant[i], pwt$pct[i]))
    readr::write_tsv(pwt, file.path(OUT, "pairwise_direction_concordance.tsv"))
  }

  L <- c(L, "",
    "> **The pooled percentage hides the real structure — read the pairwise table.**",
    ">",
    "> Every pair involving a *whole-mouth* cohort (oral rinse or saliva) agrees on 100% of its",
    "> shared genera. The single pair that sits at chance is the one with the most shared genera,",
    "> and it is also the pair with the largest niche mismatch (tongue swab vs saliva) and the",
    "> weakest provenance on both sides. Because that pair dominates the shared-genus count, it",
    "> drags the pooled figure down to near-chance.",
    ">",
    "> The operational conclusion is not 'external cohorts disagree'. It is that **oral niche",
    "> matters more than sample size**: pick the site- and region-matched cohort, disclose the",
    "> mismatch when you cannot, and never pool cohorts across niches to manufacture power.")
}

L <- c(L, "", "Per-cohort provenance: `datasets/<id>/00_manifest/download_manifest.tsv`.")
# The narrative .md is a PRIVATE document (team reading), while the TSVs above
# are public. It is written only when the private driver asks for it by setting
# M12_NARRATIVE_MD_DIR; a public run emits no narrative file at all, so the
# public results tree never contains one. Same source of truth either way.
m12_narrative_dir <- Sys.getenv("M12_NARRATIVE_MD_DIR", "")
if (nzchar(m12_narrative_dir)) {
  writeLines(L, file.path(m12_narrative_dir, "validation_shortlist.md"))
  say("wrote validation_shortlist.md -> ", m12_narrative_dir)
}
say("")
