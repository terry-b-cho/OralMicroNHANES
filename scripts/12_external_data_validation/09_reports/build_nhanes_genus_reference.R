#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Build the NHANES genus reference (per cycle) used for all overlap reporting
#
# Purpose : Every "genus overlap" number in this module needs a defensible
#           denominator. A flat list of 956 parsed genera is the WRONG one:
#
#             - It is the taxonomy of the 1,349 genus-level RSVs, not the set of
#               genera actually OBSERVED in NHANES.
#             - The two survey cycles do not share a taxon set. Both tables are
#               0-filled to the union of 1,349 columns, so "the column exists"
#               says nothing about whether the genus was detected in that cycle.
#             - A genus seen in one participant is not a validation target; a
#               genus in 90% of participants is. Raw counts hide that entirely.
#
#           This computes, per cycle F (2009-10) and G (2011-12), which genera
#           are actually detected, in how many participants, and at what
#           abundance - so overlap can be reported both as a count and weighted
#           by how much of the NHANES community it accounts for.
#
# Usage   : Rscript 09_reports/build_nhanes_genus_reference.R
# Inputs  : data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
#             tables DADA2RSV_GENUS_RELATIVE_F / _G  (columns "<otu>_relative")
#           results/.../03_gold_db_microbial_phenotype_out/intermediate/
#             ubiome_genus_mapping_complete.csv      (otu -> parsed_genus)
# Outputs : results/.../12_external_data_validation/nhanes_genus_reference.tsv
# Env     : R 4.5.1 (conda nhanes-analysis); DBI, RSQLite
# =============================================================================

suppressPackageStartupMessages({
  library(DBI); library(RSQLite); library(dplyr); library(readr); library(tibble)
})

PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
set.seed(42)

M12 <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
OUT <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation")
source(file.path(M12, "00_shared", "harmonize_helpers.R"))
m12_assert_parse_genus_agrees()

say <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")
DB <- file.path(PROJECT_ROOT,
                "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite")
stopifnot(file.exists(DB))

say("=== NHANES genus reference ===")
con <- dbConnect(RSQLite::SQLite(), DB)
on.exit(dbDisconnect(con), add = TRUE)

# --- OTU -> genus map --------------------------------------------------------
# GOTCHA (verified): empty Genus fields must be read as NA, otherwise parse_genus
# passes "" through instead of mapping it to "unclassified".
map <- read.csv(file.path(PROJECT_ROOT,
          "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv"),
          stringsAsFactors = FALSE, na.strings = c("NA", ""))
say("otu->genus map: ", nrow(map), " OTUs, ",
    length(setdiff(unique(map$parsed_genus), c("unclassified", NA))), " named genera")
genus_of <- setNames(map$parsed_genus, map$otu)

per_cycle <- list()
for (cyc in c("F", "G")) {
  tb <- paste0("DADA2RSV_GENUS_RELATIVE_", cyc)
  say("--- cycle ", cyc, " (", tb, ") ---")
  d <- dbGetQuery(con, sprintf("SELECT * FROM %s", tb))
  n_samp <- nrow(d)
  otu_cols <- setdiff(names(d), "SEQN")
  otus <- sub("_relative$", "", otu_cols)
  say("  ", n_samp, " participants x ", length(otu_cols), " genus RSVs")

  m <- as.matrix(d[, otu_cols, drop = FALSE])
  storage.mode(m) <- "double"
  m[is.na(m)] <- 0

  g <- unname(genus_of[otus])
  g[is.na(g) | !nzchar(g)] <- "unclassified"
  n_unmapped <- sum(is.na(genus_of[otus]))
  if (n_unmapped) say("  [warn] ", n_unmapped, " RSV(s) absent from the otu->genus map")

  # Collapse RSVs to genus (sum relative abundances within a genus).
  gm <- t(rowsum(t(m), group = factor(g), na.rm = TRUE))
  stopifnot(isTRUE(all.equal(sum(m), sum(gm))))   # mass preserved

  per_cycle[[cyc]] <- tibble(
    genus       = colnames(gm),
    n_detected  = as.integer(colSums(gm > 0)),
    n_samples   = n_samp,
    prevalence  = as.numeric(colMeans(gm > 0)),
    mean_abund  = as.numeric(colMeans(gm)),
    max_abund   = as.numeric(apply(gm, 2, max))
  )
  det <- sum(per_cycle[[cyc]]$n_detected > 0 & per_cycle[[cyc]]$genus != "unclassified")
  say("  genera detected in >=1 participant: ", det,
      " of ", sum(per_cycle[[cyc]]$genus != "unclassified"), " named")
}

# --- Merge cycles ------------------------------------------------------------
ref <- full_join(
  per_cycle$F %>% rename(n_detected_F = n_detected, n_samples_F = n_samples,
                         prevalence_F = prevalence, mean_abund_F = mean_abund,
                         max_abund_F = max_abund),
  per_cycle$G %>% rename(n_detected_G = n_detected, n_samples_G = n_samples,
                         prevalence_G = prevalence, mean_abund_G = mean_abund,
                         max_abund_G = max_abund),
  by = "genus") %>%
  mutate(across(starts_with("n_detected"), ~ tidyr::replace_na(.x, 0L)),
         across(c(starts_with("prevalence"), starts_with("mean_abund"),
                  starts_with("max_abund")), ~ tidyr::replace_na(.x, 0)),
         n_samples_F = max(n_samples_F, na.rm = TRUE),
         n_samples_G = max(n_samples_G, na.rm = TRUE),
         detected_F = n_detected_F > 0,
         detected_G = n_detected_G > 0,
         detected_either = detected_F | detected_G,
         detected_both   = detected_F & detected_G,
         n_detected_total = n_detected_F + n_detected_G,
         prevalence_overall = n_detected_total / (n_samples_F + n_samples_G),
         mean_abund_overall = (mean_abund_F * n_samples_F + mean_abund_G * n_samples_G) /
                              (n_samples_F + n_samples_G),
         # Tiers make "is this genus worth validating?" explicit rather than
         # leaving every reader to invent their own threshold.
         prevalence_tier = case_when(
           prevalence_overall >= 0.50 ~ "core (>=50%)",
           prevalence_overall >= 0.10 ~ "common (10-50%)",
           prevalence_overall >= 0.01 ~ "uncommon (1-10%)",
           prevalence_overall >  0    ~ "rare (<1%)",
           TRUE                       ~ "not detected")) %>%
  arrange(desc(prevalence_overall))

readr::write_tsv(ref, file.path(OUT, "nhanes_genus_reference.tsv"))

named <- ref %>% filter(genus != "unclassified")
say("")
say("wrote nhanes_genus_reference.tsv (", nrow(ref), " rows)")
say("")
say("=== NHANES genus universe, by definition used ===")
say(sprintf("  taxonomy (all parsed genera in the OTU map) : %d",
            length(setdiff(unique(map$parsed_genus), c("unclassified", NA)))))
say(sprintf("  detected in cycle F                        : %d", sum(named$detected_F)))
say(sprintf("  detected in cycle G                        : %d", sum(named$detected_G)))
say(sprintf("  detected in EITHER cycle                   : %d", sum(named$detected_either)))
say(sprintf("  detected in BOTH cycles                    : %d", sum(named$detected_both)))
say(sprintf("  F-only                                     : %d", sum(named$detected_F & !named$detected_G)))
say(sprintf("  G-only                                     : %d", sum(named$detected_G & !named$detected_F)))
say("")
say("=== prevalence tiers (named, detected) ===")
tt <- named %>% filter(detected_either) %>% count(prevalence_tier)
for (i in seq_len(nrow(tt))) say(sprintf("  %-18s %4d", tt$prevalence_tier[i], tt$n[i]))
say("")
say("=== 15 most prevalent NHANES genera ===")
top <- named %>% filter(detected_either) %>% head(15)
for (i in seq_len(nrow(top)))
  say(sprintf("  %-26s prev %5.1f%%  (F %5.1f%% / G %5.1f%%)  mean abund %6.4f",
              top$genus[i], 100*top$prevalence_overall[i],
              100*top$prevalence_F[i], 100*top$prevalence_G[i], top$mean_abund_overall[i]))
