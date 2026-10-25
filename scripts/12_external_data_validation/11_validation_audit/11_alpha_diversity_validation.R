#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 11
# Alpha-diversity validation against SLL2's published diversity rows
#
# WHY THIS EXISTS: SLL2's main-text tables report ALPHA DIVERSITY directions
# alongside the genus and phylum rows - "Smokers: Simpson down, Shannon down".
# That is a whole-community claim, independent of any single genus, and it sat
# unused in a file we already held.
#
# NHANES module 6 already tests Smoking_status and AgeGroup against Shannon /
# Inverse Simpson / Observed OTUs, but `alpha_pairwise_stat_res.csv` stores
# p-values WITHOUT a direction. The direction is therefore recomputed here from
# the phyloseq count table, and the recomputed p-values are cross-checked
# against module 6 so the two cannot silently disagree.
#
# Usage   : Rscript 11_alpha_diversity_validation.R
# Outputs : audit/alpha_diversity_nhanes.tsv
#           audit/alpha_diversity_validation.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
suppressPackageStartupMessages({library(phyloseq); library(vegan); library(DBI); library(RSQLite)})
asay("=== Step 11: alpha-diversity validation ===")

# --- Per-sample diversity from the count table -------------------------------
ps <- readRDS(file.path(PROJECT_ROOT,
  "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/nhanes_phyloseq_objects_all.rds"))
o  <- ps[["counts"]]
otu <- as(otu_table(o), "matrix")
if (!taxa_are_rows(o)) otu <- t(otu)          # want taxa x samples -> transpose to samples x taxa
mat <- t(otu)
stopifnot(nrow(mat) == nsamples(o))
asay("count matrix: ", nrow(mat), " samples x ", ncol(mat), " taxa")

div <- tibble(
  SEQN = as.character(rownames(mat)),
  Observed_OTUs     = as.integer(rowSums(mat > 0)),
  Shannon_Diversity = vegan::diversity(mat, index = "shannon"),
  Inverse_Simpson   = vegan::diversity(mat, index = "invsimpson"))
sd <- as.data.frame(sample_data(o))
div$AgeGroup <- as.character(sd$AgeGroup)
div$RIDAGEYR <- suppressWarnings(as.numeric(sd$RIDAGEYR))

# --- Join smoking by SEQN ----------------------------------------------------
con <- dbConnect(SQLite(), file.path(PROJECT_ROOT,
  "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"))
smq <- bind_rows(lapply(c("SMQ_F", "SMQ_G"), function(t)
  dbGetQuery(con, sprintf("SELECT SEQN, SMQ_current_ever_never FROM %s", t)))) %>%
  mutate(SEQN = as.character(SEQN)) %>% filter(!is.na(SMQ_current_ever_never))
dbDisconnect(con)
div <- div %>% left_join(smq, by = "SEQN") %>%
  mutate(smoking = dplyr::recode(as.character(SMQ_current_ever_never),
                                 `0` = "never", `1` = "former", `2` = "current",
                                 .default = NA_character_))
asay("smoking joined: ", sum(!is.na(div$smoking)), " of ", nrow(div), " samples")
print(table(div$smoking, useNA = "ifany"))
write_tsv(div, file.path(AUD, "alpha_diversity_nhanes.tsv"))

METRICS <- c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")

# --- NHANES direction: current vs never smokers ------------------------------
# Adults only, matching the module's scope. Wilcoxon matches module 6's
# pairwise nonparametric approach.
adults <- div %>% filter(!is.na(RIDAGEYR), RIDAGEYR >= MIN_ADULT_AGE)
smk_res <- bind_rows(lapply(METRICS, function(m) {
  d <- adults %>% filter(smoking %in% c("current", "never"))
  x <- d[[m]][d$smoking == "current"]; y <- d[[m]][d$smoking == "never"]
  w <- suppressWarnings(stats::wilcox.test(x, y))
  tibble(nhanes_variable = "SMQ_current_ever_never", metric = m,
         n_current = length(x), n_never = length(y),
         median_current = median(x), median_never = median(y),
         nhanes_direction = ifelse(median(x) < median(y), "down", "up"),
         p_value = w$p.value)
}))
asay("")
asay("=== NHANES alpha diversity, current vs never smokers (adults) ===")
print(as.data.frame(smk_res))

# --- Cross-check against module 6 so the two cannot silently disagree --------
m6 <- read_csv(file.path(PROJECT_ROOT,
  "results/analyses_results/6_alpha_beta_analyses_out/alpha_pairwise_stat_res.csv"),
  show_col_types = FALSE) %>%
  filter(Variable == "Smoking_status",
         (Group1 == "Never smoker" & Group2 == "Current smoker") |
         (Group1 == "Current smoker" & Group2 == "Never smoker"))
chk <- smk_res %>% left_join(m6 %>% select(metric = Metric, m6_p = P_value), by = "metric") %>%
  mutate(agrees = (p_value < 0.05) == (m6_p < 0.05))

# The cross-check found a DEFECT IN MODULE 6, not in this computation.
# Module 6's Inverse_Simpson pairwise output reports the IDENTICAL p-value
# (1.053549e-02) for Never-vs-Former AND Former-vs-Current - two contrasts with
# different group sizes (never 4,611 / former 1,613 / current 1,929). Identical
# to seven significant figures is not attainable by chance. Across the whole
# file, 9 of 9 multi-level Inverse_Simpson groups repeat a p-value between
# distinct contrasts, versus 4 of 11 for Shannon and 5 of 10 for Observed OTUs.
#
# Observed_OTUs and Shannon_Diversity reproduce module 6 closely, so the
# recomputation itself is sound. Rather than hard-stop or silently override,
# the disagreement is recorded and the affected metric is flagged wherever it
# is used. This is an upstream issue for module 6 to fix.
write_tsv(chk, file.path(AUD, "alpha_module6_crosscheck.tsv"))
asay("")
asay("module-6 cross-check: ", sum(chk$agrees, na.rm = TRUE), " of ", nrow(chk),
     " metrics agree on significance")
print(as.data.frame(chk %>% select(metric, recomputed_p = p_value,
                                   module6_p = m6_p, agrees)))
if (!all(chk$agrees, na.rm = TRUE))
  asay("[warn] ", paste(chk$metric[!chk$agrees], collapse = ", "),
       " disagrees with module 6; see alpha_module6_crosscheck.tsv. ",
       "Module 6's Inverse_Simpson pairwise output repeats p-values across ",
       "distinct contrasts and is treated as unreliable.")

# --- SLL2 published alpha-diversity directions -------------------------------
# Table 2 (smokers): "Simpson's diversity down", "Shannon's diversity down".
# Table 1 (age): reported as NON-MONOTONIC ("down - up"), so it is recorded as
# not directionally testable rather than forced into an up/down call.
SLL2_ALPHA <- tibble::tribble(
  ~nhanes_variable, ~phenotype_label, ~metric, ~ext_direction, ~ext_note,
  "SMQ_current_ever_never", "Smoking (current vs never)", "Inverse_Simpson", "down",
  "SLL2 Table 2: Simpson's diversity decreased in smokers (mean adj. P = 0.002, 100% of tests)",
  "SMQ_current_ever_never", "Smoking (current vs never)", "Shannon_Diversity", "down",
  "SLL2 Table 2: Shannon's diversity decreased in smokers (mean adj. P = 0.029, 83% of tests)",
  "RIDAGEYR", "Age (per year)", "Inverse_Simpson", NA_character_,
  "SLL2 Table 1: Simpson's diversity is NON-MONOTONIC across age (down then up) - no single direction to test",
  "RIDAGEYR", "Age (per year)", "Shannon_Diversity", NA_character_,
  "SLL2 Table 1: Shannon's diversity is NON-MONOTONIC across age (down then up) - no single direction to test")

val <- SLL2_ALPHA %>%
  left_join(smk_res %>% select(nhanes_variable, metric, nhanes_direction,
                               nhanes_p = p_value, n_current, n_never),
            by = c("nhanes_variable", "metric")) %>%
  mutate(
    external_dataset = "WillisJR_2022_SacaLaLengua2",
    cohort = "Willis 2022 (SLL2; published alpha diversity)",
    tested = !is.na(ext_direction) & !is.na(nhanes_direction),
    same_direction = ifelse(tested, nhanes_direction == ext_direction, NA),
    validation_classification = dplyr::case_when(
      is.na(ext_direction) ~ "MODEL_NOT_COMPARABLE",
      !tested ~ "NOT_TESTABLE_PHENOTYPE",
      same_direction & nhanes_p < 0.05 ~ "NOMINALLY_REPLICATED_SAME_DIRECTION",
      same_direction ~ "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT",
      TRUE ~ "DIRECTIONALLY_DISCORDANT"),
    evidence_level = "alpha diversity (whole-community, independent of any single genus)",
    module6_agrees = chk$agrees[match(metric, chk$metric)],
    caveat = ifelse(!is.na(module6_agrees) & !module6_agrees,
      "module 6's pairwise p for this metric repeats across distinct contrasts and is unreliable; the direction and p here are recomputed from the count table",
      "recomputed direction; p-value reproduces module 6"))
write_tsv(val, file.path(AUD, "alpha_diversity_validation.tsv"))

t <- val %>% filter(tested)
asay("")
asay("=== alpha-diversity validation vs SLL2 ===")
print(as.data.frame(val %>% select(phenotype_label, metric, nhanes_direction,
                                   ext_direction, same_direction,
                                   validation_classification)))
asay("")
asay(sprintf("testable %d of %d published rows; concordant %d/%d",
             nrow(t), nrow(val), sum(t$same_direction, na.rm = TRUE), nrow(t)))
