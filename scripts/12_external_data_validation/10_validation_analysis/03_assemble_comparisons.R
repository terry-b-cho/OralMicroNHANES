#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 2 — Step 3
# Assemble NHANES x external comparisons and score concordance
#
# Purpose : For every NHANES hit (FDR<0.05, prevalence>=10%), decide whether
#           each external dataset agrees on DIRECTION, and where possible how
#           strongly the effects correlate in rank.
#
#           THE CENTRAL DISTINCTION, enforced by assertion:
#
#             Refit cohort - we fitted the model ourselves, so every NHANES hit that
#                      the dataset carries was actually TESTED. A mismatch is
#                      genuinely DISCORDANT.
#             Published-target cohort - we hold only a published list of significant genera. A
#                      NHANES hit absent from that list was NOT TESTED (or was
#                      tested and not reported). Calling it discordant would
#                      manufacture disagreement.
#
#           Every row therefore carries one of three outcomes:
#           concordant / discordant / not tested.
#
# Usage   : Rscript 03_assemble_comparisons.R
# Outputs : validation/comparison_{age,smoking}.tsv     per-genus, three-state
#           validation/concordance_summary.tsv          per dataset x variable
#           validation/reverse_direction_summary.tsv    external hits vs NHANES
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(purrr)
})
HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/10_validation_analysis"
source(file.path(HERE, "00_validation_config.R"))

vsay("=== Step 3: assemble comparisons ===")
DSD <- file.path(M12_OUT, "datasets")

# --- Which external source answers which NHANES variable ---------------------
# `fit` names the refit-cohort model file to use as PRIMARY. For age that is the
# UNADJUSTED fit, because NHANES demoWAS carries zero covariates and a
# like-for-like comparison must match it.
#
# Smoking has NO refit-cohort row: no retained external cohort carries an
# individual-level smoking variable. Smoking evidence is published-target only, and the
# absence is asserted below rather than left to be noticed.
TIERA_MAP <- tibble::tribble(
  ~nhanes_variable, ~dataset_id,                     ~ext_variable,        ~fit,         ~primary, ~label,
  "age",     "HuangS_2020_AgePrediction",     "age",                "unadjusted",         TRUE,  COHORT_LABELS[["huang_unadj"]],
  "age",     "HuangS_2020_AgePrediction",     "age",                "adjusted",           FALSE, COHORT_LABELS[["huang_adj"]]
)

# Published-target cohorts: significant-genus lists only.
TIERB_MAP <- tibble::tribble(
  ~nhanes_variable, ~dataset_id,                    ~file,                                                                   ~taxon_col, ~dir_col,                     ~label,
  "age",     "WillisJR_2022_SacaLaLengua2",  "02_supplementary/derived/replication_age_genera.tsv",                    "taxon", "direction",            COHORT_LABELS[["sll2_age"]],
  "age",     "WellsPM_2022_TwinsUKSaliva",   "02_supplementary/derived/replication_age_genera_from_figure.tsv",        "genus", "direction_with_age",   COHORT_LABELS[["twins_age"]],
  "smoking", "WillisJR_2022_SacaLaLengua2",  "02_supplementary/derived/replication_smoking_genera.tsv",                "taxon", "direction",            COHORT_LABELS[["sll2_smoke"]],
  # Ghensi: genera enriched in peri-implant DISEASE. NHANES oradWAS models
  # gum-disease probability per unit genus CLR, so both sides encode "genus up
  # with disease" and the SIGNS align without transformation. The comparison is
  # cross-niche (subgingival plaque vs whole-mouth rinse; shotgun vs 16S V4;
  # clinically diagnosed peri-implantitis vs self-reported gum disease, n=113
  # with 23 cases) and is labelled as such wherever it appears.
  "gum_disease", "GhensiP_2020_ImplantDisease", "02_supplementary/derived/replication_disease_species.tsv", "genus", "direction_in_positive_class", COHORT_LABELS[["ghensi_gum"]]
)

norm <- function(x) m12_genus_crossdb_key(vapply(as.character(x), parse_genus,
                                                 character(1), USE.NAMES = FALSE))

comparisons <- list(); summaries <- list(); reverse <- list()

for (v in PRIMARY_VARIABLES) {
  vsay("")
  vsay("### variable: ", v)
  nh_hits <- readr::read_tsv(file.path(VAL, paste0("nhanes_effects_", v, ".tsv")),
                             show_col_types = FALSE, progress = FALSE) %>%
    mutate(key = norm(genus))
  nh_all  <- readr::read_tsv(file.path(VAL, paste0("nhanes_effects_", v, "_all.tsv")),
                             show_col_types = FALSE, progress = FALSE) %>%
    mutate(key = norm(genus))
  vsay("  NHANES hits: ", nrow(nh_hits), " genera (",
       sum(nh_hits$direction == "up"), " up / ", sum(nh_hits$direction == "down"), " down)")

  # ---- Refit cohorts ----
  ta <- TIERA_MAP %>% filter(nhanes_variable == v)
  for (i in seq_len(nrow(ta))) {
    r <- ta[i, ]
    f <- file.path(VAL, sprintf("external_effects_%s_%s_%s.tsv", r$dataset_id, r$ext_variable, r$fit))
    if (!file.exists(f)) { vsay("  [warn] missing ", basename(f)); next }
    e <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE) %>%
      mutate(key = norm(genus))
    # One row per key: if several external genera collapse to one key, keep the
    # most significant, and record that it happened.
    e1 <- e %>% group_by(key) %>% arrange(p.value, .by_group = TRUE) %>%
      slice(1) %>% ungroup()

    cmp <- nh_hits %>%
      select(genus, key, phylum, nh_estimate = estimate, nh_fdr = p.value.fdr,
             nh_direction = direction, prevalence_overall, n_obs_nhanes = n_obs) %>%
      left_join(e1 %>% select(key, ext_genus = genus, ext_estimate = estimate,
                              ext_p = p.value, ext_fdr = p.value.fdr,
                              ext_direction = direction, ext_sig = sig_fdr,
                              n_obs_ext = n_obs), by = "key") %>%
      mutate(
        tier = "A", dataset_id = r$dataset_id, label = r$label,
        primary = r$primary, fit = r$fit, nhanes_variable = v,
        # Refit cohort: tested means the genus exists in the external table.
        tested = !is.na(ext_direction),
        outcome = case_when(
          !tested                          ~ "not tested",
          ext_direction == nh_direction     ~ "concordant",
          TRUE                              ~ "discordant"))
    comparisons[[length(comparisons) + 1]] <- cmp

    tested <- cmp %>% filter(tested)
    b <- m12_binom(sum(tested$outcome == "concordant"), nrow(tested))
    rho <- if (nrow(tested) >= 4)
      suppressWarnings(stats::cor(tested$nh_estimate, tested$ext_estimate,
                                  method = "spearman", use = "complete.obs")) else NA_real_
    rho_p <- if (nrow(tested) >= 4)
      suppressWarnings(stats::cor.test(tested$nh_estimate, tested$ext_estimate,
                                       method = "spearman")$p.value) else NA_real_
    summaries[[length(summaries) + 1]] <- tibble(
      nhanes_variable = v, dataset_id = r$dataset_id, label = r$label, tier = "A",
      fit = r$fit, primary = r$primary,
      n_participants = suppressWarnings(max(cmp$n_obs_ext, na.rm = TRUE)),
      n_nhanes_hits = nrow(cmp), n_tested = nrow(tested),
      n_not_tested = sum(!cmp$tested),
      n_concordant = b$k, pct_concordant = b$pct,
      ci_lo = b$ci_lo, ci_hi = b$ci_hi, binom_p = b$p_value,
      spearman_rho = rho, spearman_p = rho_p,
      n_ext_sig_among_tested = sum(tested$ext_sig, na.rm = TRUE))
    vsay(sprintf("  [A] %-34s tested %2d/%2d | concordant %2d (%.0f%%) p=%.3f | rho=%s",
                 r$label, nrow(tested), nrow(cmp), b$k, b$pct, b$p_value,
                 ifelse(is.na(rho), "NA", sprintf("%.2f", rho))))

    # Reverse direction: external hits that NHANES tested but did not call.
    if (r$primary) {
      ext_sig <- e1 %>% filter(sig_fdr)
      rev <- ext_sig %>% select(key, ext_genus = genus, ext_direction = direction) %>%
        left_join(nh_all %>% select(key, nh_direction = direction,
                                    nh_fdr = p.value.fdr), by = "key") %>%
        mutate(in_nhanes = !is.na(nh_direction),
               nhanes_sig = !is.na(nh_fdr) & nh_fdr < FDR_THRESHOLD,
               agree = in_nhanes & ext_direction == nh_direction)
      reverse[[length(reverse) + 1]] <- tibble(
        nhanes_variable = v, dataset_id = r$dataset_id, label = r$label,
        n_ext_sig = nrow(rev),
        n_present_in_nhanes = sum(rev$in_nhanes),
        n_also_sig_in_nhanes = sum(rev$nhanes_sig),
        n_direction_agree = sum(rev$agree, na.rm = TRUE))
    }
  }

  # ---- Published-target cohorts ----
  tb <- TIERB_MAP %>% filter(nhanes_variable == v)
  for (i in seq_len(nrow(tb))) {
    r <- tb[i, ]
    f <- file.path(DSD, r$dataset_id, r$file)
    if (!file.exists(f)) { vsay("  [warn] missing ", r$file); next }
    t <- readr::read_tsv(f, show_col_types = FALSE, progress = FALSE)
    # Targets are published at whichever rank the source used. Ghensi's is
    # species-level and is rolled up to genus by the `taxon_col` below, so
    # filtering to rank == "genus" alone silently discarded all 259 rows and
    # reported 0 of 52 testable. Keep both ranks and assert we kept something.
    if ("rank" %in% names(t)) t <- t %>% filter(rank %in% c("genus", "species"))
    if (!nrow(t))
      stop("TARGET EMPTY after rank filter for ", r$dataset_id, " / ", v,
           " - check the `rank` values in ", r$file)
    t <- t %>% mutate(key = norm(.data[[r$taxon_col]]),
                      ext_direction = as.character(.data[[r$dir_col]])) %>%
      filter(!is.na(ext_direction), ext_direction %in% c("up", "down")) %>%
      group_by(key) %>% slice(1) %>% ungroup()

    cmp <- nh_hits %>%
      select(genus, key, phylum, nh_estimate = estimate, nh_fdr = p.value.fdr,
             nh_direction = direction, prevalence_overall, n_obs_nhanes = n_obs) %>%
      left_join(t %>% select(key, ext_direction), by = "key") %>%
      mutate(
        tier = "B", dataset_id = r$dataset_id, label = r$label,
        primary = TRUE, fit = "published", nhanes_variable = v,
        ext_genus = NA_character_, ext_estimate = NA_real_, ext_p = NA_real_,
        ext_fdr = NA_real_, ext_sig = NA, n_obs_ext = NA_integer_,
        # Published-target cohort: only genera on the published list were "tested" from our side.
        tested = !is.na(ext_direction),
        outcome = case_when(
          !tested                       ~ "not tested",
          ext_direction == nh_direction  ~ "concordant",
          TRUE                           ~ "discordant"))
    comparisons[[length(comparisons) + 1]] <- cmp

    tested <- cmp %>% filter(tested)
    b <- m12_binom(sum(tested$outcome == "concordant"), nrow(tested))
    summaries[[length(summaries) + 1]] <- tibble(
      nhanes_variable = v, dataset_id = r$dataset_id, label = r$label, tier = "B",
      fit = "published", primary = TRUE, n_participants = NA_integer_,
      n_nhanes_hits = nrow(cmp), n_tested = nrow(tested), n_not_tested = sum(!cmp$tested),
      n_concordant = b$k, pct_concordant = b$pct, ci_lo = b$ci_lo, ci_hi = b$ci_hi,
      binom_p = b$p_value, spearman_rho = NA_real_, spearman_p = NA_real_,
      n_ext_sig_among_tested = nrow(tested))
    vsay(sprintf("  [B] %-34s tested %2d/%2d | concordant %2d (%s)",
                 r$label, nrow(tested), nrow(cmp), b$k,
                 ifelse(is.na(b$pct), "n/a", sprintf("%.0f%%", b$pct))))
  }
}

cmp_all <- bind_rows(comparisons)
sm_all  <- bind_rows(summaries)
rv_all  <- bind_rows(reverse)

# --- Assertions --------------------------------------------------------------
# 1. The three states must partition the NHANES hits exactly.
chk <- cmp_all %>% count(nhanes_variable, dataset_id, label, fit, outcome) %>%
  tidyr::pivot_wider(names_from = outcome, values_from = n, values_fill = 0) %>%
  left_join(sm_all %>% select(nhanes_variable, dataset_id, label, fit, n_nhanes_hits),
            by = c("nhanes_variable", "dataset_id", "label", "fit"))
bad <- chk %>% mutate(tot = concordant + discordant + `not tested`) %>%
  filter(tot != n_nhanes_hits)
if (nrow(bad)) { print(as.data.frame(bad)); stop("PARTITION ASSERTION FAILED") }
vsay("")
vsay("[OK] concordant + discordant + not tested == n_nhanes_hits for all ",
     nrow(chk), " dataset x fit combinations")

# 2. No published-target row may be 'discordant' for a genus absent from its list.
badB <- cmp_all %>% filter(tier == "B", outcome == "discordant", is.na(ext_direction))
if (nrow(badB)) stop("TIER-B ASSERTION FAILED: discordant call without a published direction")
vsay("[OK] no published-target discordant call rests on an absent genus")

# 3. Smoking is published-target only. No retained external cohort carries an
#    individual-level smoking variable, so no smoking row may claim a Spearman
#    correlation. If a refit-cohort smoking cohort is ever added, this assertion
#    fires and the report's smoking limitation must be rewritten with it.
smk <- sm_all %>% filter(nhanes_variable == "smoking")
if (any(smk$tier == "A"))
  stop("SMOKING ASSERTION FAILED: a refit-cohort smoking fit appeared; update the ",
       "report's smoking limitation, which states there is none")
if (any(!is.na(smk$spearman_rho)))
  stop("SMOKING ASSERTION FAILED: a smoking correlation exists without a refit-cohort fit")
vsay("[OK] smoking is published-target only (", nrow(smk), " published source",
     ifelse(nrow(smk) == 1, "", "s"), "; no effect-size correlation possible)")

for (v in unique(cmp_all$nhanes_variable))
  readr::write_tsv(cmp_all %>% filter(nhanes_variable == v),
                   file.path(VAL, paste0("comparison_", v, ".tsv")))
readr::write_tsv(sm_all, file.path(VAL, "concordance_summary.tsv"))
readr::write_tsv(rv_all, file.path(VAL, "reverse_direction_summary.tsv"))

vsay("")
vsay("=== concordance summary ===")
print(as.data.frame(sm_all %>%
  transmute(variable = nhanes_variable, label, tier, primary,
            N = n_participants, hits = n_nhanes_hits, tested = n_tested,
            conc = n_concordant, pct = round(pct_concordant),
            p = signif(binom_p, 2), rho = round(spearman_rho, 2))))
vsay("")
vsay("=== reverse direction (external hits vs NHANES) ===")
print(as.data.frame(rv_all))

# =============================================================================
# Coverage: how much of the NHANES result is externally checkable AT ALL
#
# Concordance percentages are conditional on a genus having been tested. That
# denominator hides the more basic question - what fraction of our hits any
# external dataset can even see. Reporting 100% concordance on 9 of 48 genera
# without that context would materially overstate the evidence.
# =============================================================================
vsay("")
vsay("=== coverage: externally checkable fraction ===")
cov_genus <- cmp_all %>%
  filter(primary) %>%
  group_by(nhanes_variable, genus, phylum, nh_estimate, nh_fdr, nh_direction,
           prevalence_overall) %>%
  summarise(n_datasets_testing = sum(tested),
            datasets_testing = paste(sort(label[tested]), collapse = "; "),
            n_concordant = sum(outcome == "concordant"),
            n_discordant = sum(outcome == "discordant"),
            .groups = "drop") %>%
  mutate(coverage = ifelse(n_datasets_testing == 0, "no external coverage",
                           "tested externally"))
readr::write_tsv(cov_genus, file.path(VAL, "coverage_per_genus.tsv"))

cov_sum <- cov_genus %>%
  group_by(nhanes_variable) %>%
  summarise(n_hits = n(),
            n_tested_by_any = sum(n_datasets_testing > 0),
            n_untested      = sum(n_datasets_testing == 0),
            pct_tested_by_any = 100 * mean(n_datasets_testing > 0),
            max_datasets_testing = max(n_datasets_testing),
            .groups = "drop")
readr::write_tsv(cov_sum, file.path(VAL, "coverage_summary.tsv"))
for (i in seq_len(nrow(cov_sum))) {
  r <- cov_sum[i, ]
  vsay(sprintf("  %-8s %2d hits | tested by >=1: %2d (%.0f%%) | NO coverage: %2d (%.0f%%)",
               r$nhanes_variable, r$n_hits, r$n_tested_by_any, r$pct_tested_by_any,
               r$n_untested, 100 - r$pct_tested_by_any))
}
# Every hit is either covered or not - assert the partition so a future join
# bug cannot silently drop genera from the denominator.
stopifnot(all(cov_sum$n_tested_by_any + cov_sum$n_untested == cov_sum$n_hits))

# NOTE: the serum-biomarker measurement-robustness check lives in
# notebooks/03_internal_robustness.qmd, which computes it inline and writes
# no files. It is a within-NHANES check, not external validation, and is
# deliberately kept out of this pipeline and its report.
