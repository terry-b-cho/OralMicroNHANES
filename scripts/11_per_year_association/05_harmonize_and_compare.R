#!/usr/bin/env Rscript

# =============================================================================
# MODULE 11 — STEP 05: harmonized master table + F-vs-G + pooled-vs-cycle
# 05_harmonize_and_compare.R   (run AFTER 03 aggregation and 04 availability)
# =============================================================================
# Builds the harmonized comparison table (one row per main-effect test, with
# pooled/F/G estimates side by side, taxonomy, prevalence, availability), then:
#   * F_vs_G/           Pearson r + direction agreement per analysis x subset
#                       (incl. prevalence>10% primary subset), per-variable
#                       correlations, significance counts
#   * pooled_vs_F(_G)/  replication of pooled-significant associations,
#                       testability-aware (never counts a one-cycle-only
#                       variable as failed replication)
#   * tables/           table3 (significance counts), table4 (per-variable
#                       correlation summary), table5 (pooled-sig row level)
#
# Main effects only (fdr_corrected == TRUE), joined on
# (analysis, dependent_var, exposure, term). Thresholds: FDR/q <= 0.05 primary,
# 0.10 context. Direction > magnitude: betas are compared by sign and Pearson
# correlation, never by equality.
# =============================================================================

source("scripts/11_per_year_association/00_module11_config.R")  # run from the repo root

cat("MODULE 11 STEP 05: harmonize + compare\n")

tax  <- m11_load_taxonomy()
prev <- m11_load_prevalence()
avail <- read_csv(file.path(M11_AVAIL, "variable_availability.csv"), show_col_types = FALSE)

var_col_for <- function(analysis) {
  if (analysis %in% c("1_demoWAS", "3_exWAS")) "exposure" else "dependent_var"
}
otu_col_for <- function(analysis) {
  if (analysis %in% c("1_demoWAS", "3_exWAS")) "dependent_var" else "exposure"
}

# --- Per-cycle p-value CALIBRATION (design degrees of freedom) ---------------
# svyglm sets residual df = (PSUs - strata) - #coefficients + 1. A single NHANES
# cycle has only ~16 design df, so the 13-covariate models are left with ~2
# residual df; a t2 reference has a hard p-value floor (~1e-4 even at |t|=250),
# which makes BH rejection mathematically impossible at family sizes of
# 5,396-638,565 — the reason per-cycle FDR discoveries collapsed to zero
# outside covariate-free 1_demoWAS. Following NCHS / Korn-Graubard practice we
# recalibrate per-cycle p-values on t with df = DESIGN df (PSUs - strata,
# computed from the DEMO tables), then recompute BH-FDR and Storey q within
# each (analysis x cycle) family. Standard errors are unchanged (they are
# design-correct); only the reference distribution changes. The svyglm-default
# quantities are retained as *_svydf columns. Pooled results (design df ~31,
# residual df ~19) are unaffected and remain Module 1's published reference.
m11_design_df <- local({
  con <- DBI::dbConnect(RSQLite::SQLite(), M11_DB)
  on.exit(DBI::dbDisconnect(con))
  vapply(c(F = "DEMO_F", G = "DEMO_G"), function(tb) {
    x <- DBI::dbGetQuery(con, sprintf(
      "SELECT COUNT(DISTINCT SDMVSTRA) AS ns,
              COUNT(DISTINCT SDMVSTRA || '_' || SDMVPSU) AS np FROM %s", tb))
    as.numeric(x$np - x$ns)
  }, numeric(1))
})
names(m11_design_df) <- c("F", "G")
cat("Per-cycle design df (PSUs - strata): F =", m11_design_df[["F"]],
    ", G =", m11_design_df[["G"]], "\n")

load_main_effects_recal <- function(an, cyc) {
  d <- m11_load_main_effects(an, cyc)
  if (is.null(d)) return(NULL)
  dfd <- m11_design_df[[cyc]]
  tstat <- abs(d$estimate / d$std.error)
  p_re <- 2 * pt(-tstat, df = dfd)
  fdr_re <- p.adjust(p_re, "fdr")
  q_re <- fdr_re
  okp <- is.finite(p_re)
  if (sum(okp) >= 10) {
    q_ok <- tryCatch(qvalue::qvalue(p_re[okp])$qvalues, error = function(e) NULL)
    if (!is.null(q_ok)) q_re[okp] <- q_ok
  }
  d %>% mutate(p.value_svydf = p.value,
               p.value.fdr_svydf = p.value.fdr,
               q.value_svydf = q.value,
               p.value = p_re, p.value.fdr = fdr_re, q.value = q_re,
               design_df_used = dfd)
}

# ---------------- build harmonized master table ------------------------------
key <- c("dependent_var", "exposure", "term")
harmonized <- list()

calib_rows <- list()
for (an in M11_ANALYSES) {
  pm <- m11_load_main_effects(an, "pooled")
  fm <- load_main_effects_recal(an, "F")
  gm <- load_main_effects_recal(an, "G")
  if (is.null(fm) || is.null(gm)) stop("Missing per-year aggregation for ", an)
  if (is.null(pm)) stop("Missing POOLED results for ", an, " at ", m11_tidied_path(an, "pooled"))

  for (cc in c("F", "G")) {
    d <- if (cc == "F") fm else gm
    calib_rows[[paste(an, cc)]] <- tibble(
      analysis = an, cycle = cc, design_df = m11_design_df[[cc]],
      n_family = nrow(d),
      n_fdr05_svyglm_default = sum(d$p.value.fdr_svydf < 0.05, na.rm = TRUE),
      n_fdr05_design_df      = sum(d$p.value.fdr < 0.05, na.rm = TRUE),
      min_p_svyglm_default   = min(d$p.value_svydf, na.rm = TRUE),
      min_p_design_df        = min(d$p.value, na.rm = TRUE))
  }

  slim <- function(d, sfx) d %>%
    select(all_of(key), estimate, std.error, p.value, p.value.fdr, q.value, n_obs,
           any_of(c("p.value_svydf", "p.value.fdr_svydf", "q.value_svydf"))) %>%
    rename_with(~ paste0(.x, "_", sfx), -all_of(key))

  h <- slim(pm, "pooled") %>%
    full_join(slim(fm, "F"), by = key) %>%
    full_join(slim(gm, "G"), by = key) %>%
    mutate(analysis = an,
           variable = .data[[var_col_for(an)]],
           otu_variable = .data[[otu_col_for(an)]])
  harmonized[[an]] <- h
}

harmonized <- bind_rows(harmonized) %>%
  left_join(tax %>% select(variable, otu, phylum_std, Phylum, genus_label),
            by = c("otu_variable" = "variable")) %>%
  left_join(prev %>% select(variable, prevalence_pooled, prevalence_F, prevalence_G),
            by = c("otu_variable" = "variable")) %>%
  left_join(avail %>% select(analysis, variable, availability,
                             direct_FG_comparison_valid, reason_not_comparable),
            by = c("analysis", "variable")) %>%
  mutate(
    is_logistic = analysis %in% M11_LOGISTIC_ANALYSES,
    testable_F = is.finite(estimate_F),
    testable_G = is.finite(estimate_G),
    stable_pooled = m11_stable_fit(estimate_pooled, std.error_pooled, is_logistic),
    stable_F = m11_stable_fit(estimate_F, std.error_F, is_logistic),
    stable_G = m11_stable_fit(estimate_G, std.error_G, is_logistic),
    sign_pooled = sign(estimate_pooled),
    sign_F = sign(estimate_F), sign_G = sign(estimate_G),
    sign_agree_FG       = testable_F & testable_G & sign_F == sign_G,
    sign_agree_F_pooled = testable_F & is.finite(estimate_pooled) & sign_F == sign_pooled,
    sign_agree_G_pooled = testable_G & is.finite(estimate_pooled) & sign_G == sign_pooled)

for (th in c(M11_THRESH_PRIMARY, M11_THRESH_CONTEXT)) {
  tag <- gsub("\\.", "", sprintf("%02d", th * 100))
  harmonized[[paste0("sig_fdr", tag, "_pooled")]] <- !is.na(harmonized$p.value.fdr_pooled) & harmonized$p.value.fdr_pooled <= th
  harmonized[[paste0("sig_fdr", tag, "_F")]]      <- !is.na(harmonized$p.value.fdr_F) & harmonized$p.value.fdr_F <= th
  harmonized[[paste0("sig_fdr", tag, "_G")]]      <- !is.na(harmonized$p.value.fdr_G) & harmonized$p.value.fdr_G <= th
  harmonized[[paste0("sig_q", tag, "_pooled")]]   <- !is.na(harmonized$q.value_pooled) & harmonized$q.value_pooled <= th
  harmonized[[paste0("sig_q", tag, "_F")]]        <- !is.na(harmonized$q.value_F) & harmonized$q.value_F <= th
  harmonized[[paste0("sig_q", tag, "_G")]]        <- !is.na(harmonized$q.value_G) & harmonized$q.value_G <= th
}

saveRDS(harmonized, file.path(M11_CMP, "harmonized_master_table.rds"))
write_csv(harmonized, file.path(M11_CMP, "harmonized_master_table.csv.gz"))
cat("Harmonized master table:", nrow(harmonized), "main-effect tests\n")

calib <- bind_rows(calib_rows)
write_csv(calib, file.path(M11_FIGURES_ACCESSORY, "optional_table9_percycle_pvalue_calibration.csv"))
cat("\nPer-cycle p-value calibration (svyglm-default residual df vs design df):\n")
print(as.data.frame(calib), digits = 3)

# Participant N per cycle/analysis (for reporting alongside every correlation)
pN <- harmonized %>% group_by(analysis) %>%
  summarise(participant_n_F = max(n_obs_F, na.rm = TRUE),
            participant_n_G = max(n_obs_G, na.rm = TRUE), .groups = "drop")

# ---------------- F vs G agreement -------------------------------------------
# Comparable set: testable AND numerically stable in both cycles AND the
# variable's coding is directly comparable (availability == F_AND_G).
# One-cycle-only, non-comparable, and unstable-fit tests are excluded here and
# inventoried separately (never silently dropped).
fg <- harmonized %>% filter(testable_F, testable_G, stable_F, stable_G,
                            direct_FG_comparison_valid)

n_unstable <- harmonized %>%
  filter(testable_F, testable_G, direct_FG_comparison_valid, !(stable_F & stable_G)) %>% nrow()
cat("F-vs-G comparable tests:", nrow(fg), "of", nrow(harmonized),
    "(", n_unstable, "excluded as numerically unstable fits )\n")

subset_defs <- function(d, th) {
  tag <- gsub("\\.", "", sprintf("%02d", th * 100))
  list(
    all_shared              = rep(TRUE, nrow(d)),
    prevalence_gt10         = !is.na(d$prevalence_pooled) & d$prevalence_pooled > M11_PREV_THRESHOLD,
    pooled_fdr_sig          = d[[paste0("sig_fdr", tag, "_pooled")]],
    pooled_q_sig            = d[[paste0("sig_q", tag, "_pooled")]],
    pooled_fdr_sig_prev10   = d[[paste0("sig_fdr", tag, "_pooled")]] &
                              !is.na(d$prevalence_pooled) & d$prevalence_pooled > M11_PREV_THRESHOLD,
    F_fdr_sig               = d[[paste0("sig_fdr", tag, "_F")]],
    G_fdr_sig               = d[[paste0("sig_fdr", tag, "_G")]],
    sig_either_cycle_fdr    = d[[paste0("sig_fdr", tag, "_F")]] | d[[paste0("sig_fdr", tag, "_G")]],
    sig_both_cycles_fdr     = d[[paste0("sig_fdr", tag, "_F")]] & d[[paste0("sig_fdr", tag, "_G")]],
    sig_both_cycles_q       = d[[paste0("sig_q", tag, "_F")]] & d[[paste0("sig_q", tag, "_G")]]
  )
}

agree_rows <- list()
for (an in c(M11_ANALYSES, "ALL")) {
  d <- if (an == "ALL") fg else fg %>% filter(analysis == an)
  if (nrow(d) == 0) next
  nF <- if (an == "ALL") max(pN$participant_n_F) else pN$participant_n_F[pN$analysis == an]
  nG <- if (an == "ALL") max(pN$participant_n_G) else pN$participant_n_G[pN$analysis == an]
  for (th in c(M11_THRESH_PRIMARY, M11_THRESH_CONTEXT)) {
    ss <- subset_defs(d, th)
    for (sn in names(ss)) {
      idx <- ss[[sn]]; idx[is.na(idx)] <- FALSE
      pr <- m11_pearson(d$estimate_F[idx], d$estimate_G[idx])
      da <- m11_dir_agree(d$estimate_F[idx], d$estimate_G[idx])
      agree_rows[[length(agree_rows) + 1]] <- tibble(
        analysis = an, threshold = th, subset = sn,
        n_tests = sum(idx),
        pearson_r = pr$r, pearson_p = pr$p, n_in_correlation = pr$n,
        n_direction_compared = da$n, n_direction_agree = da$n_agree,
        direction_agreement_pct = da$pct,
        participant_n_F = nF, participant_n_G = nG)
    }
  }
}
fg_agreement <- bind_rows(agree_rows)
write_csv(fg_agreement, file.path(M11_CMP, "F_vs_G", "fg_agreement_summary.csv"))

# Significance counts (Table 3)
cnt_rows <- list()
for (an in c(M11_ANALYSES, "ALL")) {
  d <- if (an == "ALL") fg else fg %>% filter(analysis == an)
  if (nrow(d) == 0) next
  for (th in c(M11_THRESH_PRIMARY, M11_THRESH_CONTEXT)) {
    tag <- gsub("\\.", "", sprintf("%02d", th * 100))
    for (metric in c("fdr", "q")) {
      sf <- d[[paste0("sig_", metric, tag, "_F")]]
      sg <- d[[paste0("sig_", metric, tag, "_G")]]
      cnt_rows[[length(cnt_rows) + 1]] <- tibble(
        analysis = an, metric = metric, threshold = th,
        n_shared_tests = nrow(d),
        n_sig_F = sum(sf), n_sig_G = sum(sg),
        n_sig_both = sum(sf & sg), n_sig_either = sum(sf | sg),
        n_sig_both_same_direction = sum(sf & sg & d$sign_agree_FG))
    }
  }
}
sig_counts <- bind_rows(cnt_rows)
write_csv(sig_counts, file.path(M11_CMP, "F_vs_G", "significance_counts.csv"))
write_csv(sig_counts, file.path(M11_FIGURES_ACCESSORY, "optional_table3_significance_counts.csv"))

# Per-variable correlations (Table 4 / Figure 4)
pv_rows <- list()
for (an in M11_ANALYSES) {
  d_an <- fg %>% filter(analysis == an)
  for (v in unique(d_an$variable)) {
    d <- d_an %>% filter(variable == v)
    d10 <- d %>% filter(!is.na(prevalence_pooled), prevalence_pooled > M11_PREV_THRESHOLD)
    pr   <- m11_pearson(d$estimate_F, d$estimate_G)
    pr10 <- m11_pearson(d10$estimate_F, d10$estimate_G)
    da   <- m11_dir_agree(d$estimate_F, d$estimate_G)
    tag  <- gsub("\\.", "", sprintf("%02d", M11_THRESH_PRIMARY * 100))
    pv_rows[[length(pv_rows) + 1]] <- tibble(
      analysis = an, variable = v,
      n_overlapping_tests = nrow(d),
      pearson_r = pr$r, pearson_p = pr$p,
      pearson_r_prev10 = pr10$r, n_prev10 = pr10$n,
      direction_agreement_pct = da$pct, n_direction_compared = da$n,
      n_fdr_sig_F = sum(d[[paste0("sig_fdr", tag, "_F")]]),
      n_fdr_sig_G = sum(d[[paste0("sig_fdr", tag, "_G")]]),
      n_fdr_sig_both = sum(d[[paste0("sig_fdr", tag, "_F")]] & d[[paste0("sig_fdr", tag, "_G")]]),
      n_same_direction = sum(d$sign_agree_FG),
      participant_n_F = suppressWarnings(max(d$n_obs_F, na.rm = TRUE)),
      participant_n_G = suppressWarnings(max(d$n_obs_G, na.rm = TRUE)))
  }
}
per_variable <- bind_rows(pv_rows) %>% arrange(analysis, desc(n_overlapping_tests))
write_csv(per_variable, file.path(M11_CMP, "F_vs_G", "per_variable_correlations.csv"))
write_csv(per_variable, file.path(M11_FIGURES_ACCESSORY, "optional_table4_per_variable_correlations.csv"))

# Scale-free overall F-vs-G summary: betas of different variables live in
# different native units (especially 3_exWAS / 4_pheWAS), so the PRIMARY overall
# magnitude summary is the distribution of per-variable Pearson r values, not a
# single pooled-across-units correlation.
pv_overall <- per_variable %>%
  filter(is.finite(pearson_r), n_overlapping_tests >= 10) %>%
  group_by(analysis) %>%
  summarise(n_variables = n(),
            median_per_variable_r = median(pearson_r),
            iqr_lo_r = quantile(pearson_r, 0.25), iqr_hi_r = quantile(pearson_r, 0.75),
            weighted_mean_r = sum(pearson_r * n_overlapping_tests) / sum(n_overlapping_tests),
            median_direction_agreement_pct = median(direction_agreement_pct, na.rm = TRUE),
            .groups = "drop")
write_csv(pv_overall, file.path(M11_CMP, "F_vs_G", "per_variable_r_overall_summary.csv"))

# Inventory of tests excluded from F-vs-G (one-cycle-only, not comparable, or unstable)
excluded <- harmonized %>%
  filter(!(testable_F & testable_G & stable_F & stable_G & direct_FG_comparison_valid)) %>%
  mutate(exclusion_reason = case_when(
    !testable_F & !testable_G ~ "not fitted in either cycle",
    !testable_F ~ "not testable in F (variable/pair absent)",
    !testable_G ~ "not testable in G (variable/pair absent)",
    !direct_FG_comparison_valid ~ coalesce(reason_not_comparable, "coding not comparable"),
    !(stable_F & stable_G) ~ "numerically unstable fit (quasi-separation guard)")) %>%
  select(analysis, variable, dependent_var, exposure, term, availability,
         exclusion_reason, testable_F, testable_G,
         estimate_pooled, p.value.fdr_pooled, q.value_pooled)
write_csv(excluded, file.path(M11_AVAIL, "fg_excluded_tests.csv.gz"))
cat("Excluded from F-vs-G comparison (inventoried):", nrow(excluded), "tests\n")

# ---------------- pooled-vs-cycle replication --------------------------------
tag <- gsub("\\.", "", sprintf("%02d", M11_THRESH_PRIMARY * 100))

classify_replication <- function(d, metric) {
  sp <- d[[paste0("sig_", metric, tag, "_pooled")]]
  ds <- d[sp & !is.na(sp), ]
  if (nrow(ds) == 0) return(NULL)
  sf <- ds[[paste0("sig_fdr", tag, "_F")]]
  sg <- ds[[paste0("sig_fdr", tag, "_G")]]
  ds %>% mutate(
    per_year_sig_F = sf, per_year_sig_G = sg,
    replication_category = case_when(
      !testable_F & !testable_G ~ "not testable in either cycle",
      !testable_F ~ "not testable in F",
      !testable_G ~ "not testable in G",
      !(stable_F & stable_G) ~ "per-year fit unstable (direction not assessable)",
      sf & sg & sign_agree_F_pooled & sign_agree_G_pooled ~ "significant in F & G, same direction",
      sf & sg ~ "significant in F & G, direction discordant",
      sf & sign_agree_F_pooled & sign_agree_G_pooled ~ "significant in F only, same direction in G",
      sg & sign_agree_F_pooled & sign_agree_G_pooled ~ "significant in G only, same direction in F",
      sf | sg ~ "significant in one cycle, direction discordant in the other",
      sign_agree_F_pooled & sign_agree_G_pooled ~ "not significant in either, same direction in both",
      sign_agree_F_pooled | sign_agree_G_pooled ~ "directionally concordant in one cycle only",
      TRUE ~ "directionally discordant in both cycles"),
    pooled_metric = metric)
}

rep_rowlevel <- bind_rows(
  classify_replication(harmonized, "fdr"),
  classify_replication(harmonized, "q"))

# --- Replication tiers (message-level summary) --------------------------------
# For a pooled-significant, pre-specified association, the appropriate per-cycle
# replication test is DIRECTION-INFORMED: one-sided replication p = p/2 when the
# per-cycle beta has the pooled sign, 1 - p/2 otherwise. Per-cycle BH-FDR (which
# re-corrects within each cycle as if hypotheses were new) is reported as the
# strictest tier. Tiers are computed only where the fit is testable AND stable.
rep_rowlevel <- rep_rowlevel %>% mutate(
  dir_concord_F = testable_F & stable_F & is.finite(estimate_pooled) &
                  sign(estimate_F) == sign(estimate_pooled),
  dir_concord_G = testable_G & stable_G & is.finite(estimate_pooled) &
                  sign(estimate_G) == sign(estimate_pooled),
  onesided_p_F = ifelse(dir_concord_F, p.value_F / 2, 1 - p.value_F / 2),
  onesided_p_G = ifelse(dir_concord_G, p.value_G / 2, 1 - p.value_G / 2),
  onesided_rep_F = dir_concord_F & !is.na(onesided_p_F) & onesided_p_F < 0.05,
  onesided_rep_G = dir_concord_G & !is.na(onesided_p_G) & onesided_p_G < 0.05,
  nom_sig_F = dir_concord_F & !is.na(p.value_F) & p.value_F < 0.05,
  nom_sig_G = dir_concord_G & !is.na(p.value_G) & p.value_G < 0.05,
  fdrsig_dir_F = dir_concord_F & per_year_sig_F,
  fdrsig_dir_G = dir_concord_G & per_year_sig_G)

# REPLICATION FDR: for pre-specified pooled-significant hypotheses the correct
# multiplicity family is the pooled-significant set itself (replication-panel
# practice), not the whole per-cycle WAS family. BH over the one-sided
# directional replication p-values within each (pooled metric x analysis),
# among tests testable+stable in both cycles.
bh_in <- function(p, ok) { out <- rep(NA_real_, length(p)); out[ok] <- p.adjust(p[ok], "fdr"); out }
rep_rowlevel <- rep_rowlevel %>%
  group_by(pooled_metric, analysis) %>%
  mutate(ok_ts = testable_F & testable_G & stable_F & stable_G,
         repfdr_p_F = bh_in(onesided_p_F, ok_ts),
         repfdr_p_G = bh_in(onesided_p_G, ok_ts),
         repfdr_F = dir_concord_F & !is.na(repfdr_p_F) & repfdr_p_F < 0.05,
         repfdr_G = dir_concord_G & !is.na(repfdr_p_G) & repfdr_p_G < 0.05) %>%
  ungroup() %>% select(-ok_ts)

# Table 6: tier counts per analysis x pooled metric, on the testable+stable set
tier_rows <- list()
for (met in c("fdr", "q")) {
  for (an in c(M11_ANALYSES, "ALL")) {
    d <- rep_rowlevel %>% filter(pooled_metric == met)
    if (an != "ALL") d <- d %>% filter(analysis == an)
    if (nrow(d) == 0) next
    ds <- d %>% filter(testable_F, testable_G, stable_F, stable_G)
    n0 <- nrow(ds)
    bt <- if (n0 > 0) binom.test(sum(ds$dir_concord_F & ds$dir_concord_G), n0, p = 0.5) else NULL
    tier_rows[[length(tier_rows) + 1]] <- tibble(
      analysis = an, pooled_metric = met,
      n_pooled_sig = nrow(d),
      n_testable_stable_both = n0,
      n_dir_concordant_F = sum(ds$dir_concord_F), n_dir_concordant_G = sum(ds$dir_concord_G),
      n_dir_concordant_both = sum(ds$dir_concord_F & ds$dir_concord_G),
      pct_dir_concordant_both = ifelse(n0 > 0, 100 * mean(ds$dir_concord_F & ds$dir_concord_G), NA),
      dir_concord_ci_lo = if (n0 > 0) 100 * bt$conf.int[1] else NA_real_,
      dir_concord_ci_hi = if (n0 > 0) 100 * bt$conf.int[2] else NA_real_,
      dir_concord_binom_p_vs_50pct = if (n0 > 0) bt$p.value else NA_real_,
      n_onesided_rep_F = sum(ds$onesided_rep_F), n_onesided_rep_G = sum(ds$onesided_rep_G),
      n_onesided_rep_both = sum(ds$onesided_rep_F & ds$onesided_rep_G),
      pct_onesided_rep_both = ifelse(n0 > 0, 100 * mean(ds$onesided_rep_F & ds$onesided_rep_G), NA),
      n_nominal_sig_F = sum(ds$nom_sig_F), n_nominal_sig_G = sum(ds$nom_sig_G),
      n_nominal_sig_both = sum(ds$nom_sig_F & ds$nom_sig_G),
      pct_nominal_sig_both = ifelse(n0 > 0, 100 * mean(ds$nom_sig_F & ds$nom_sig_G), NA),
      n_repfdr_F = sum(ds$repfdr_F), n_repfdr_G = sum(ds$repfdr_G),
      n_repfdr_both = sum(ds$repfdr_F & ds$repfdr_G),
      pct_repfdr_both = ifelse(n0 > 0, 100 * mean(ds$repfdr_F & ds$repfdr_G), NA),
      n_percycle_fdr_F = sum(ds$fdrsig_dir_F), n_percycle_fdr_G = sum(ds$fdrsig_dir_G),
      n_percycle_fdr_both = sum(ds$fdrsig_dir_F & ds$fdrsig_dir_G),
      pct_percycle_fdr_both = ifelse(n0 > 0, 100 * mean(ds$fdrsig_dir_F & ds$fdrsig_dir_G), NA))
  }
}
tiers <- bind_rows(tier_rows)
write_csv(tiers, file.path(M11_TABLES, "table3_replication_tiers.csv"))
cat("\nReplication tiers (pooled FDR set, testable+stable both):\n")
print(tiers %>% filter(pooled_metric == "fdr") %>%
        select(analysis, n_testable_stable_both, pct_dir_concordant_both,
               pct_onesided_rep_both, pct_repfdr_both, pct_percycle_fdr_both))

# Table 7: gold list — pooled FDR-sig, direction-concordant AND one-sided
# replicated in BOTH independent cycles
# GOLD LIST = the VALIDATED core: pooled FDR-significant associations that are
# replication-FDR<=0.05 with the pooled direction in BOTH independent cycles.
gold <- rep_rowlevel %>%
  filter(pooled_metric == "fdr", repfdr_F & repfdr_G) %>%
  arrange(p.value.fdr_pooled) %>%
  transmute(analysis, variable, dependent_var, exposure, term,
            genus = genus_label, phylum = phylum_std,
            beta_pooled = estimate_pooled, se_pooled = std.error_pooled,
            p_pooled = p.value_pooled, fdr_pooled = p.value.fdr_pooled, q_pooled = q.value_pooled,
            beta_F = estimate_F, se_F = std.error_F, p_F = p.value_F,
            onesided_p_F, repfdr_p_F, fdr_F = p.value.fdr_F,
            beta_G = estimate_G, se_G = std.error_G, p_G = p.value_G,
            onesided_p_G, repfdr_p_G, fdr_G = p.value.fdr_G,
            n_obs_pooled, n_obs_F, n_obs_G, prevalence_pooled,
            onesided_rep_both = onesided_rep_F & onesided_rep_G,
            percycle_fdr_both = fdrsig_dir_F & fdrsig_dir_G)
write_csv(gold, file.path(M11_TABLES, "table4_validated_core.csv"))
cat("\nGold list (VALIDATED: replication-FDR<=0.05 + pooled direction in BOTH cycles):",
    nrow(gold), "associations\n")

# --- Table 10: VALIDATION of pooled-significant hits (the headline table) ----
# Denominator: pooled-significant tests (per metric). "Validated" in a cycle =
# replication-FDR<=0.05 (BH within the pooled-significant family) AND the beta
# carries the pooled direction. A pooled hit counts as validated overall if
# validated in EITHER cycle; "both" is the strongest tier.
val_rows <- list()
for (met in c("fdr", "q")) {
  for (an in c(M11_ANALYSES, "ALL")) {
    d <- rep_rowlevel %>% filter(pooled_metric == met)
    if (an != "ALL") d <- d %>% filter(analysis == an)
    if (nrow(d) == 0) next
    ds <- d %>% filter(testable_F, testable_G, stable_F, stable_G)
    n0 <- nrow(ds)
    vF <- ds$repfdr_F; vG <- ds$repfdr_G
    val_rows[[length(val_rows) + 1]] <- tibble(
      analysis = an, pooled_metric = met,
      n_pooled_sig = nrow(d),
      n_testable_stable_both = n0,
      n_not_testable_or_unstable = nrow(d) - n0,
      n_validated_F = sum(vF), n_validated_G = sum(vG),
      n_validated_F_only = sum(vF & !vG), n_validated_G_only = sum(!vF & vG),
      n_validated_both = sum(vF & vG),
      n_validated_either = sum(vF | vG),
      n_not_validated = sum(!vF & !vG),
      pct_validated_either_of_testable = ifelse(n0 > 0, 100 * mean(vF | vG), NA),
      pct_validated_both_of_testable   = ifelse(n0 > 0, 100 * mean(vF & vG), NA),
      pct_validated_either_of_pooled   = 100 * sum(vF | vG) / nrow(d))
  }
}
validation <- bind_rows(val_rows)
write_csv(validation, file.path(M11_TABLES, "table2_validation_summary.csv"))
cat("\nVALIDATION of pooled hits (replication-FDR + pooled direction; fdr metric):\n")
print(as.data.frame(validation %>% filter(pooled_metric == "fdr") %>%
  select(analysis, n_pooled_sig, n_testable_stable_both, n_validated_F_only,
         n_validated_G_only, n_validated_both, n_validated_either,
         pct_validated_either_of_testable, pct_validated_both_of_testable)), digits = 3)

# --- Table 8: significance-overlap counts (F vs G, and pooled -> per-cycle) ---
# Three questions, per analysis and overall:
#  Q1 nominal_p05_all   : two-sided p<0.05 in each cycle independently, over all
#                         F-vs-G-comparable (testable+stable) tests
#  Q2 percycle_fdr05_all: per-cycle BH-FDR<=0.05 in each cycle, same base
#  Q3 pooled_sig_p05    : among pooled BH-FDR<=0.05 tests, two-sided p<0.05 in
#                         each cycle (direction ignored; direction-consistent
#                         'both' subset reported as an extra category)
base8 <- harmonized %>%
  filter(testable_F, testable_G, stable_F, stable_G, direct_FG_comparison_valid)
overlap_counts <- function(d, sF, sG, question) {
  tibble(question = question,
         category = c("cycle F only", "cycle G only", "both cycles", "neither cycle"),
         n = c(sum(sF & !sG), sum(!sF & sG), sum(sF & sG), sum(!sF & !sG)),
         denominator = nrow(d))
}
t8_rows <- list()
for (an in c(M11_ANALYSES, "ALL")) {
  d <- if (an == "ALL") base8 else base8 %>% filter(analysis == an)
  if (nrow(d) == 0) next
  nomF <- !is.na(d$p.value_F) & d$p.value_F < 0.05
  nomG <- !is.na(d$p.value_G) & d$p.value_G < 0.05
  fdrF <- d$sig_fdr05_F; fdrG <- d$sig_fdr05_G
  t8_rows[[length(t8_rows) + 1]] <- overlap_counts(d, nomF, nomG, "nominal_p05_all") %>%
    mutate(analysis = an)
  t8_rows[[length(t8_rows) + 1]] <- overlap_counts(d, fdrF, fdrG, "percycle_fdr05_all") %>%
    mutate(analysis = an)
  dp <- d %>% filter(sig_fdr05_pooled)
  pF <- !is.na(dp$p.value_F) & dp$p.value_F < 0.05
  pG <- !is.na(dp$p.value_G) & dp$p.value_G < 0.05
  dirok <- sign(dp$estimate_F) == sign(dp$estimate_pooled) &
           sign(dp$estimate_G) == sign(dp$estimate_pooled)
  t8_rows[[length(t8_rows) + 1]] <- bind_rows(
    overlap_counts(dp, pF, pG, "pooled_sig_p05"),
    tibble(question = "pooled_sig_p05",
           category = "both cycles, direction-consistent",
           n = sum(pF & pG & dirok), denominator = nrow(dp))) %>%
    mutate(analysis = an)
}
t8 <- bind_rows(t8_rows) %>% select(analysis, question, category, n, denominator)
write_csv(t8, file.path(M11_FIGURES_ACCESSORY, "optional_table8_significance_overlap.csv"))
cat("Table 8 (significance overlap) written:",
    nrow(t8), "rows\n")

REP_LEVELS <- c("significant in F & G, same direction",
                "significant in F only, same direction in G",
                "significant in G only, same direction in F",
                "not significant in either, same direction in both",
                "directionally concordant in one cycle only",
                "significant in one cycle, direction discordant in the other",
                "significant in F & G, direction discordant",
                "directionally discordant in both cycles",
                "per-year fit unstable (direction not assessable)",
                "not testable in F", "not testable in G", "not testable in either cycle")
rep_rowlevel$replication_category <- factor(rep_rowlevel$replication_category, levels = REP_LEVELS)

saveRDS(rep_rowlevel, file.path(M11_CMP, "pooled_vs_F_vs_G", "pooled_replication_rowlevel.rds"))

# Table 5: row-level pooled-significant validation (FDR set; q flag included)
t5 <- rep_rowlevel %>%
  filter(pooled_metric == "fdr") %>%
  transmute(analysis, variable, dependent_var, exposure, term,
            genus = genus_label, phylum = phylum_std,
            beta_pooled = estimate_pooled, se_pooled = std.error_pooled,
            fdr_pooled = p.value.fdr_pooled, q_pooled = q.value_pooled,
            beta_F = estimate_F, se_F = std.error_F, fdr_F = p.value.fdr_F, q_F = q.value_F,
            beta_G = estimate_G, se_G = std.error_G, fdr_G = p.value.fdr_G, q_G = q.value_G,
            n_obs_pooled, n_obs_F, n_obs_G,
            prevalence_pooled, prevalence_F, prevalence_G,
            availability, testable_F, testable_G,
            sign_agree_F_pooled, sign_agree_G_pooled, sign_agree_FG,
            replication_category)
write_csv(t5, file.path(M11_FIGURES_ACCESSORY, "optional_table5_pooled_significant_validation.csv"))

# Summary counts per analysis x metric
rep_summary <- rep_rowlevel %>%
  group_by(analysis, pooled_metric, replication_category, .drop = FALSE) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(analysis, pooled_metric) %>%
  mutate(n_pooled_sig = sum(n), pct = ifelse(n_pooled_sig > 0, 100 * n / n_pooled_sig, NA)) %>%
  ungroup()
write_csv(rep_summary, file.path(M11_CMP, "pooled_vs_F_vs_G", "replication_category_summary.csv"))

# Per-cycle pooled-vs-cycle scatter data + beta correlations among pooled-sig
for (cyc in M11_CYCLES) {
  d <- rep_rowlevel %>% filter(pooled_metric == "fdr",
                               .data[[paste0("testable_", cyc)]],
                               .data[[paste0("stable_", cyc)]], stable_pooled)
  pr <- m11_pearson(d$estimate_pooled, d[[paste0("estimate_", cyc)]])
  da <- m11_dir_agree(d$estimate_pooled, d[[paste0("estimate_", cyc)]])
  out <- tibble(cycle = cyc, n_pooled_sig_testable = nrow(d),
                pearson_r_pooled_vs_cycle = pr$r, pearson_p = pr$p, n_in_correlation = pr$n,
                direction_agreement_pct = da$pct, n_direction_compared = da$n)
  write_csv(out, file.path(M11_CMP, paste0("pooled_vs_", cyc), "pooled_vs_cycle_summary.csv"))
  write_csv(d,   file.path(M11_CMP, paste0("pooled_vs_", cyc), "pooled_sig_rowlevel.csv.gz"))
}

cat("\nReplication category summary (pooled FDR<=", M11_THRESH_PRIMARY, "):\n", sep = "")
print(rep_summary %>% filter(pooled_metric == "fdr", n > 0) %>%
        select(analysis, replication_category, n, pct), n = 60)

cat("\nSTEP 05 complete.\n")
