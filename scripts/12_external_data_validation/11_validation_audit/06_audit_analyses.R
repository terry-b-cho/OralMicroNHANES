#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 6
# The analyses the audit found missing (work order section 10, A-I)
#
# Everything is computed FROM master_validation_table.tsv so the figures,
# tables and report cannot drift from the evidence base.
#
# Every correlation reports its n TAXA. Every percentage reports its
# DENOMINATOR. Participant N and taxon N are kept in separate columns and are
# never interchanged.
#
# Usage   : Rscript 06_audit_analyses.R
# Outputs : audit/concordance_by_phenotype.tsv   (A, D)
#           audit/beta_correlations.tsv          (B)
#           audit/enrichment_tests.tsv           (C)
#           audit/multi_cohort_agreement.tsv     (E)
#           audit/leave_one_study_out.tsv        (G)
#           audit/prevalence_stratified.tsv      (H)
#           audit/evidence_levels.tsv            (I)
#           audit/meta_analysis_assessment.md    (F)
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 6: audit analyses (A-I) ===")

mt <- read_tsv(file.path(AUD, "master_validation_table.tsv"), show_col_types = FALSE)
P  <- mt %>% filter(primary)

# =============================================================================
# A + D  Sign concordance and replication rate, with explicit denominators
# =============================================================================
conc <- P %>%
  group_by(analysis_type, nhanes_variable, phenotype_label, external_dataset,
           label, tier, oral_site, sequencing_type, adjustment_comparability) %>%
  summarise(
    n_nhanes_primary       = n(),                       # denominator 1: eligible
    n_testable             = sum(tested),               # denominator 2: testable
    n_not_testable         = sum(!tested),
    n_same_direction       = sum(same_direction, na.rm = TRUE),
    n_discordant           = sum(!same_direction, na.rm = TRUE),
    n_nominal_significant  = sum(external_nominal_significant, na.rm = TRUE),
    n_fdr_significant      = sum(external_fdr_significant, na.rm = TRUE),
    n_significant_same_dir = sum(significant_same_direction, na.rm = TRUE),
    n_participants_external = suppressWarnings(max(ext_n, na.rm = TRUE)),
    n_participants_nhanes   = suppressWarnings(max(nhanes_n, na.rm = TRUE)),
    .groups = "drop") %>%
  mutate(n_participants_external = ifelse(is.finite(n_participants_external),
                                          n_participants_external, NA_integer_)) %>%
  rowwise() %>%
  mutate(b = list(m12_binom(n_same_direction, n_testable))) %>%
  ungroup() %>%
  mutate(pct_same_direction = vapply(b, function(x) x$pct, numeric(1)),
         ci_lo = vapply(b, function(x) x$ci_lo, numeric(1)),
         ci_hi = vapply(b, function(x) x$ci_hi, numeric(1)),
         binom_p = vapply(b, function(x) x$p_value, numeric(1)),
         # The denominator is carried IN the label so it cannot be quoted alone.
         concordance_statement = sprintf("%d/%d (%.0f%%)", n_same_direction,
                                         n_testable, pct_same_direction),
         replication_rate_fdr = sprintf("%d/%d", n_significant_same_dir, n_testable)) %>%
  select(-b)
write_tsv(conc, file.path(AUD, "concordance_by_phenotype.tsv"))
asay("A/D concordance + replication rate: ", nrow(conc), " phenotype x cohort rows")

# =============================================================================
# B  Beta correlations, refit-cohort only (published-target publishes no coefficient)
# =============================================================================
bc <- P %>% filter(tier == "A", tested, !is.na(ext_beta)) %>%
  group_by(analysis_type, nhanes_variable, phenotype_label, label,
           adjustment_comparability, nhanes_model_role) %>%
  summarise(
    n_taxa = n(),                                     # ALWAYS reported
    n_participants_external = max(ext_n, na.rm = TRUE),
    pearson_r  = ifelse(n() > 2, stats::cor(nhanes_beta, ext_beta), NA_real_),
    pearson_p  = ifelse(n() > 2, stats::cor.test(nhanes_beta, ext_beta)$p.value, NA_real_),
    spearman_rho = ifelse(n() > 2, suppressWarnings(
      stats::cor(nhanes_beta, ext_beta, method = "spearman")), NA_real_),
    pct_same_direction = 100 * mean(sign(nhanes_beta) == sign(ext_beta)),
    .groups = "drop") %>%
  # Scales are only comparable when BOTH sides model genus as the outcome.
  # "INDEPENDENT" CONTAINS the substring "DEPENDENT", so a naive grepl marks
  # pheWAS/oradWAS as scale-comparable when they are the reversed orientation.
  # Key off the scheme, which is unambiguous.
  mutate(scales_comparable = analysis_type %in% c("demoWAS", "exWAS"),
         magnitude_note = ifelse(scales_comparable,
           "both model genus_CLR ~ phenotype: magnitudes on the same scale",
           "NHANES models phenotype ~ genus while the external fit is genus ~ phenotype: RANK only, magnitudes not comparable"))
write_tsv(bc, file.path(AUD, "beta_correlations.tsv"))
asay("B beta correlations: ", nrow(bc), " refit-cohort phenotype rows (n_taxa always reported)")

# =============================================================================
# C  Enrichment: do NHANES-SIGNIFICANT taxa agree more than non-significant?
#
# Needs the non-significant NHANES associations too, so the source WAS is
# re-read for the refit-cohort phenotypes only. Fisher's exact compares directional
# agreement in the significant vs non-significant strata.
# =============================================================================
map <- read.csv(file.path(PROJECT_ROOT,
  "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv"),
  stringsAsFactors = FALSE, na.strings = c("NA", "")) %>% select(otu, parsed_genus)
norm <- function(x) m12_genus_crossdb_key(vapply(as.character(x), parse_genus,
                                                 character(1), USE.NAMES = FALSE))

TIERA <- tibble::tribble(
  ~nhanes_variable, ~rds_rel, ~gcol, ~pcol, ~ext_var,
  "RIDAGEYR", "results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds", "dependent_var", "independent_var", "age",
  "RIAGENDR", "results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds", "dependent_var", "independent_var", "sex",
  "EDUCATION_COLLEGEGRAD", "results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds", "dependent_var", "independent_var", "education",
  "BMXBMI", "results/4_pheWAS_out/result_clr/4_pheWAS_clr_tidied_complete.rds", "independent_var", "dependent_var", "bmi")

enr <- list()
for (i in seq_len(nrow(TIERA))) {
  r <- TIERA[i, ]
  d <- readRDS(file.path(PROJECT_ROOT, r$rds_rel)) %>%
    filter(fdr_corrected, .data[[r$pcol]] == r$nhanes_variable)
  if (!nrow(d)) next
  d <- d %>% mutate(otu = sub("_relative$", "", .data[[r$gcol]])) %>%
    left_join(map, by = "otu") %>%
    filter(!is.na(parsed_genus), !grepl("^unclassified", parsed_genus)) %>%
    group_by(parsed_genus) %>% arrange(p.value.fdr, .by_group = TRUE) %>%
    slice(1) %>% ungroup() %>% mutate(key = norm(parsed_genus))
  f <- file.path(VAL, sprintf("external_effects_HuangS_2020_AgePrediction_%s_unadjusted.tsv", r$ext_var))
  if (!file.exists(f)) next
  ex <- read_tsv(f, show_col_types = FALSE, progress = FALSE) %>%
    mutate(key = norm(genus)) %>% group_by(key) %>%
    arrange(p.value, .by_group = TRUE) %>% slice(1) %>% ungroup()
  j <- d %>% inner_join(ex %>% select(key, e_est = estimate), by = "key") %>%
    mutate(sig = p.value.fdr < FDR_THRESHOLD,
           agree = sign(estimate) == sign(e_est))
  tb <- table(factor(j$sig, c(FALSE, TRUE)), factor(j$agree, c(FALSE, TRUE)))
  ft <- suppressWarnings(stats::fisher.test(tb))
  enr[[length(enr) + 1]] <- tibble(
    nhanes_variable = r$nhanes_variable, external_dataset = "HuangS_2020_AgePrediction",
    n_taxa_total = nrow(j),
    n_taxa_nhanes_sig = sum(j$sig), n_taxa_nhanes_nonsig = sum(!j$sig),
    pct_agree_sig = 100 * mean(j$agree[j$sig]),
    pct_agree_nonsig = 100 * mean(j$agree[!j$sig]),
    fisher_or = unname(ft$estimate), fisher_p = ft$p.value,
    interpretation = ifelse(ft$p.value < 0.05,
      ifelse(unname(ft$estimate) > 1,
             "NHANES-significant taxa agree MORE than non-significant - evidence the signal is real",
             "NHANES-significant taxa agree LESS - no enrichment"),
      "no detectable difference between significant and non-significant taxa"))
  rm(d); gc(verbose = FALSE)
}
enr_all <- bind_rows(enr)
write_tsv(enr_all, file.path(AUD, "enrichment_tests.tsv"))
asay("C enrichment: ", nrow(enr_all), " phenotypes tested")

# =============================================================================
# E  Multi-cohort agreement - only age has >= 2 external cohorts
# =============================================================================
multi <- P %>% filter(tested) %>%
  group_by(analysis_type, nhanes_variable, phenotype_label, genus, phylum,
           nhanes_direction, nhanes_beta, nhanes_fdr, nhanes_prevalence) %>%
  summarise(n_cohorts_testing = n(),
            n_cohorts_agree = sum(same_direction, na.rm = TRUE),
            cohorts = paste(sort(label), collapse = "; "),
            all_agree = n_cohorts_agree == n_cohorts_testing,
            n_fdr_sig = sum(external_fdr_significant, na.rm = TRUE),
            .groups = "drop") %>%
  filter(n_cohorts_testing > 1) %>%
  arrange(desc(n_cohorts_testing), desc(all_agree))
write_tsv(multi, file.path(AUD, "multi_cohort_agreement.tsv"))
asay("E multi-cohort: ", nrow(multi), " associations testable in >1 cohort; ",
     sum(multi$all_agree), " agree in every cohort")

# =============================================================================
# G  Leave-one-study-out - does concordance depend on a single cohort?
# =============================================================================
loo <- list()
for (v in unique(P$nhanes_variable)) {
  sub <- P %>% filter(nhanes_variable == v, tested)
  ds <- unique(sub$label)
  if (length(ds) < 2) next
  for (d in ds) {
    kept <- sub %>% filter(label != d)
    b <- m12_binom(sum(kept$same_direction, na.rm = TRUE), nrow(kept))
    loo[[length(loo) + 1]] <- tibble(
      nhanes_variable = v, phenotype_label = sub$phenotype_label[1],
      dropped_cohort = d, n_cohorts_remaining = length(ds) - 1,
      n_taxa_remaining = nrow(kept), n_same_direction = b$k,
      pct_same_direction = b$pct, binom_p = b$p_value)
  }
}
loo_all <- bind_rows(loo)
if (nrow(loo_all)) write_tsv(loo_all, file.path(AUD, "leave_one_study_out.tsv"))
asay("G leave-one-study-out: ", nrow(loo_all), " drops evaluated")

# =============================================================================
# H  Prevalence-stratified concordance at 10 / 25 / 50%
# =============================================================================
prev <- P %>% filter(tested) %>%
  mutate(prev_stratum = cut(nhanes_prevalence, c(0.10, 0.25, 0.50, 1.0),
                            labels = c("10-25%", "25-50%", ">50%"),
                            include.lowest = TRUE)) %>%
  group_by(phenotype_label, prev_stratum) %>%
  summarise(n_taxa = n(), n_same_direction = sum(same_direction, na.rm = TRUE),
            .groups = "drop") %>%
  rowwise() %>% mutate(b = list(m12_binom(n_same_direction, n_taxa))) %>% ungroup() %>%
  mutate(pct_same_direction = vapply(b, function(x) x$pct, numeric(1)),
         ci_lo = vapply(b, function(x) x$ci_lo, numeric(1)),
         ci_hi = vapply(b, function(x) x$ci_hi, numeric(1))) %>% select(-b)
write_tsv(prev, file.path(AUD, "prevalence_stratified.tsv"))
asay("H prevalence strata: ", nrow(prev), " phenotype x stratum rows")

# =============================================================================
# I  Evidence levels reported SEPARATELY (section 8)
# =============================================================================
lev <- P %>% group_by(phenotype_label, label, tier) %>%
  summarise(
    eligible_nhanes_primary = n(),
    L0_testable            = sum(tested),
    L1_same_direction      = sum(same_direction, na.rm = TRUE),
    L3_nominal_significant = sum(external_nominal_significant, na.rm = TRUE),
    L4_fdr_significant     = sum(external_fdr_significant, na.rm = TRUE),
    L5_strong_replication  = sum(significant_same_direction, na.rm = TRUE),
    discordant             = sum(!same_direction, na.rm = TRUE), .groups = "drop")
write_tsv(lev, file.path(AUD, "evidence_levels.tsv"))

# =============================================================================
# F  Meta-analysis assessment - stated, not forced
# =============================================================================
n_tierA_per_pheno <- P %>% filter(tier == "A") %>% distinct(phenotype_label, label) %>%
  count(phenotype_label, name = "n_tierA_cohorts")
# The narrative .md is a PRIVATE document (team reading), while the TSVs above
# are public. It is written only when the private driver asks for it by setting
# M12_NARRATIVE_MD_DIR; a public run emits no narrative file at all, so the
# public results tree never contains one. Same source of truth either way.
m12_narrative_dir <- Sys.getenv("M12_NARRATIVE_MD_DIR", "")
meta_md <- c(
"# Meta-analysis assessment",
"",
"**A fixed- or random-effects meta-analysis is NOT scientifically defensible here, and none is",
"performed.** The work order asks for one only where at least two genuinely independent cohorts",
"estimate sufficiently comparable effects. That condition fails on three separate grounds:",
"",
"1. **No phenotype has two refit cohorts.** Huang is the only retained cohort with a genus x",
"   sample table, so every phenotype we fit ourselves has exactly one effect estimate:",
paste0("   ", paste(sprintf("%s (%d)", n_tierA_per_pheno$phenotype_label,
                            n_tierA_per_pheno$n_tierA_cohorts), collapse = "; ")), ".",
"2. **published-target contributes no variance.** SLL2, TwinsUK and Ghensi publish significant-genus lists",
"   with directions but no coefficients or standard errors. There is nothing to weight.",
"3. **Coefficient scales differ by design.** oradWAS and pheWAS model phenotype ~ genus while the",
"   external fits model genus ~ phenotype; pooling those magnitudes would be meaningless even if",
"   the standard errors existed.",
"",
"Pooling anyway would manufacture a precision the data do not support. Direction concordance with",
"exact binomial CIs, and rank correlation where both sides yield coefficients, are the strongest",
"defensible summaries and are what this module reports.")
if (nzchar(m12_narrative_dir)) {
  writeLines(meta_md, file.path(m12_narrative_dir, "meta_analysis_assessment.md"))
  asay("F meta-analysis: assessed as NOT defensible; reasons written -> ", m12_narrative_dir)
} else {
  asay("F meta-analysis: assessed as NOT defensible")
}

asay("")
asay("=== A/D concordance summary ===")
print(as.data.frame(conc %>%
  mutate(p = signif(binom_p, 3)) %>%
  select(phenotype_label, label, tier, n_nhanes_primary, n_testable,
         concordance_statement, p, L5 = n_significant_same_dir, disc = n_discordant)))
asay("")
asay("=== B beta correlations ===")
print(as.data.frame(bc %>%
  mutate(across(c(pearson_r, spearman_rho, pct_same_direction), ~ round(.x, 3))) %>%
  select(phenotype_label, n_taxa, pearson_r, spearman_rho, pct_same_direction,
         scales_comparable)))
asay("")
asay("=== C enrichment ===")
print(as.data.frame(enr_all %>%
  mutate(across(c(pct_agree_sig, pct_agree_nonsig, fisher_or), ~ round(.x, 2)),
         fisher_p = signif(fisher_p, 3)) %>%
  select(nhanes_variable, n_taxa_total, n_taxa_nhanes_sig, pct_agree_sig,
         pct_agree_nonsig, fisher_or, fisher_p)))
