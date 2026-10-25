#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 5
# The MASTER association-level validation table
#
# Purpose : One row per (NHANES variable x genus x external cohort) for every
#           phenotype the audit found testable - age, smoking, sex, BMI,
#           education and oral disease - across refit and published-target cohorts. This is the
#           canonical machine-readable evidence base; every figure and table
#           downstream is generated FROM it, never re-derived independently.
#
#           Evidence is graded on five levels (work order section 8) rather
#           than a binary replicated yes/no:
#             L1 direction        sign(beta_NHANES) == sign(beta_external)
#             L2 magnitude        handled at the cohort level in step 6
#             L3 nominal          external p < 0.05
#             L4 FDR              external FDR < 0.05
#             L5 strong           L1 & L4
#
#           SIGN ALIGNMENT is applied BEFORE any concordance is computed, from
#           the machine-readable sign_alignment_log.tsv written in step 4. A
#           comparison whose orientation has not been verified is refused.
#
# Usage   : Rscript 05_master_validation_table.R
# Outputs : audit/master_validation_table.tsv
#           audit/comparison_registry.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 5: master validation table ===")

DSD  <- file.path(M12_OUT, "datasets")
prim <- read_tsv(file.path(AUD, "nhanes_primary_associations.tsv"), show_col_types = FALSE)
sgn  <- read_tsv(file.path(AUD, "sign_alignment_log.tsv"), show_col_types = FALSE)

norm <- function(x) m12_genus_crossdb_key(vapply(as.character(x), parse_genus,
                                                 character(1), USE.NAMES = FALSE))

# =============================================================================
# COMPARISON REGISTRY
#
# Every NHANES variable x external cohort pair that the opportunity matrix
# marked testable. `fit` names the refit-cohort model file; for published-target it is
# "published". `primary` marks the headline model per pair (the LIKE-FOR-LIKE
# fit: unadjusted where NHANES itself is crude).
# =============================================================================
REG <- tibble::tribble(
  ~nhanes_variable, ~analysis_type, ~phenotype_label, ~dataset_id, ~tier, ~ext_variable, ~fit, ~primary, ~label, ~target_file, ~taxon_col, ~dir_col,

  # ---------- AGE ----------
  "RIDAGEYR", "demoWAS", "Age (per year)", "HuangS_2020_AgePrediction", "A", "age", "unadjusted", TRUE,
  COHORT_LABELS[["huang_unadj"]], NA, NA, NA,
  "RIDAGEYR", "demoWAS", "Age (per year)", "HuangS_2020_AgePrediction", "A", "age", "adjusted", FALSE,
  COHORT_LABELS[["huang_adj"]], NA, NA, NA,
  "RIDAGEYR", "demoWAS", "Age (per year)", "WillisJR_2022_SacaLaLengua2", "B", NA, "published", TRUE,
  COHORT_LABELS[["sll2_age"]], "02_supplementary/derived/replication_age_genera.tsv", "taxon", "direction",
  "RIDAGEYR", "demoWAS", "Age (per year)", "WellsPM_2022_TwinsUKSaliva", "B", NA, "published", TRUE,
  COHORT_LABELS[["twins_age"]], "02_supplementary/derived/replication_age_genera_from_figure.tsv", "genus", "direction_with_age",

  # ---------- SMOKING ----------
  "SMQ_current_ever_never", "exWAS", "Smoking (current vs never)", "WillisJR_2022_SacaLaLengua2", "B", NA, "published", TRUE,
  COHORT_LABELS[["sll2_smoke"]], "02_supplementary/derived/replication_smoking_genera.tsv", "taxon", "direction",

  # ---------- SEX ----------
  "RIAGENDR", "demoWAS", "Sex (female vs male)", "HuangS_2020_AgePrediction", "A", "sex", "unadjusted", TRUE,
  COHORT_LABELS[["huang_unadj"]], NA, NA, NA,
  "RIAGENDR", "demoWAS", "Sex (female vs male)", "HuangS_2020_AgePrediction", "A", "sex", "adjusted", FALSE,
  "Huang 2020 (adj. age)", NA, NA, NA,

  # ---------- BMI ----------
  "BMXBMI", "pheWAS", "Body mass index", "HuangS_2020_AgePrediction", "A", "bmi", "unadjusted", TRUE,
  COHORT_LABELS[["huang_unadj"]], NA, NA, NA,
  "BMXBMI", "pheWAS", "Body mass index", "HuangS_2020_AgePrediction", "A", "bmi", "adjusted", FALSE,
  "Huang 2020 (adj. age+sex)", NA, NA, NA,

  # ---------- EDUCATION ----------
  "EDUCATION_COLLEGEGRAD", "demoWAS", "Education (college grad+)", "HuangS_2020_AgePrediction", "A", "education", "unadjusted", TRUE,
  COHORT_LABELS[["huang_unadj"]], NA, NA, NA,
  "EDUCATION_COLLEGEGRAD", "demoWAS", "Education (college grad+)", "HuangS_2020_AgePrediction", "A", "education", "adjusted", FALSE,
  "Huang 2020 (adj. age+sex)", NA, NA, NA,

  # ---------- ORAL DISEASE ----------
  # Heavily caveated: subgingival plaque + shotgun vs whole-mouth rinse + 16S,
  # peri-implantitis vs self-reported gum disease, n = 113 (23 cases).
  "GUM_DISEASE_OHAROCGP", "oradWAS", "Gum disease", "GhensiP_2020_ImplantDisease", "B", NA, "published", TRUE,
  COHORT_LABELS[["ghensi_gum"]], "02_supplementary/derived/replication_disease_species.tsv",
  "genus", "direction_in_positive_class"
)
write_tsv(REG, file.path(AUD, "comparison_registry.tsv"))

# --- Per-cohort technical descriptors (stated on every row, section 18) ------
COHORT_META <- tibble::tribble(
  ~dataset_id, ~oral_site, ~sequencing_type, ~taxonomy_method, ~external_n_note,
  "HuangS_2020_AgePrediction", "tongue-dorsum swab 73.5% / oral rinse / saliva",
  "16S V4 (515F/806R)", "Greengenes/eHOMD genus rollup", "adults, one sample per subject",
  "WillisJR_2022_SacaLaLengua2", "whole-mouth PBS rinse", "16S V3-V4",
  "published genus names", "published table; no participant-level data",
  "WellsPM_2022_TwinsUKSaliva", "saliva", "16S V4", "published genus names",
  "published figure; no participant-level data",
  "GhensiP_2020_ImplantDisease", "subgingival plaque (peri-implant)", "shotgun metagenomics",
  "MetaPhlAn species -> genus", "n=113 samples, 23 peri-implantitis, 70 healthy controls")

# --- Covariate comparability (section 7) ------------------------------------
NHANES_COVS <- c(demoWAS = "NONE (crude; n_covariates = 0)",
                 exWAS = "age, age^2, sex, PIR, education x4, ethnicity x4, US-born (13)",
                 pheWAS = "age, age^2, sex, PIR, education x4, ethnicity x4, US-born (13)",
                 oradWAS = "age, age^2, sex, PIR, education x4, ethnicity x4, US-born (13)",
                 outWAS = "age, age^2, sex, PIR, education x4, ethnicity x4, US-born (13)")

rows <- list()

for (i in seq_len(nrow(REG))) {
  r <- REG[i, ]
  nh <- prim %>% filter(nhanes_variable == r$nhanes_variable) %>% mutate(key = norm(genus))
  if (!nrow(nh)) { asay("[warn] no NHANES primary hits for ", r$nhanes_variable); next }

  # --- Sign alignment must be on record BEFORE anything is compared ---------
  sl <- sgn %>% filter(nhanes_variable == r$nhanes_variable, external_dataset == r$dataset_id)
  if (!nrow(sl))
    stop("SIGN ALIGNMENT MISSING for ", r$nhanes_variable, " x ", r$dataset_id,
         ". Refusing to compute concordance on an unverified orientation.")
  stopifnot(isTRUE(sl$signs_comparable[1]))

  if (r$tier == "A") {
    f <- file.path(VAL, sprintf("external_effects_%s_%s_%s.tsv", r$dataset_id, r$ext_variable, r$fit))
    if (!file.exists(f)) { asay("[warn] missing ", basename(f)); next }
    ex <- read_tsv(f, show_col_types = FALSE, progress = FALSE) %>%
      mutate(key = norm(genus)) %>%
      group_by(key) %>% arrange(p.value, .by_group = TRUE) %>% slice(1) %>% ungroup()
    j <- nh %>%
      left_join(ex %>% select(key, ext_genus = genus, ext_estimate = estimate,
                              ext_se = std.error, ext_p = p.value, ext_fdr = p.value.fdr,
                              ext_direction = direction, ext_n = n_obs,
                              ext_covariates = covariates_used), by = "key")
    ext_cov <- ifelse(all(is.na(j$ext_covariates)) || all(j$ext_covariates == "", na.rm = TRUE),
                      "none (crude)", j$ext_covariates[!is.na(j$ext_covariates)][1])
  } else {
    tf <- file.path(DSD, r$dataset_id, r$target_file)
    if (!file.exists(tf)) { asay("[warn] missing target ", r$target_file); next }
    t <- read_tsv(tf, show_col_types = FALSE, progress = FALSE)
    if ("rank" %in% names(t)) t <- t %>% filter(rank == "genus" | rank == "species")
    t <- t %>%
      mutate(key = norm(.data[[r$taxon_col]]),
             ext_direction = as.character(.data[[r$dir_col]])) %>%
      filter(!is.na(ext_direction), ext_direction %in% c("up", "down")) %>%
      group_by(key) %>% slice(1) %>% ungroup()
    j <- nh %>%
      left_join(t %>% select(key, ext_direction), by = "key") %>%
      mutate(ext_genus = NA_character_, ext_estimate = NA_real_, ext_se = NA_real_,
             ext_p = NA_real_, ext_fdr = NA_real_, ext_n = NA_integer_)
    ext_cov <- "published result - covariates as reported by the source paper"
  }

  cm <- COHORT_META %>% filter(dataset_id == r$dataset_id)
  nh_cov <- unname(NHANES_COVS[r$analysis_type])

  # Adjustment comparability (section 7). NHANES demoWAS is crude, so an
  # unadjusted external fit is MATCHED, not a deficiency.
  adj <- if (r$tier == "B") "NOT_RECONSTRUCTABLE" else {
    nh_crude <- grepl("^NONE", nh_cov)
    ex_crude <- identical(ext_cov, "none (crude)")
    if (nh_crude && ex_crude) "MATCHED"
    else if (ex_crude) "UNADJUSTED_EXTERNAL"
    else "PARTIALLY_MATCHED"
  }

  rows[[length(rows) + 1]] <- j %>%
    transmute(
      analysis_type = r$analysis_type,
      nhanes_variable = r$nhanes_variable,
      phenotype_label = r$phenotype_label,
      genus, phylum,
      nhanes_prevalence = prevalence_overall,
      nhanes_beta = estimate, nhanes_se = std.error, nhanes_p = p.value,
      nhanes_fdr = p.value.fdr, nhanes_q = q.value, nhanes_n = n_obs,
      nhanes_direction = direction,
      nhanes_model_role = ifelse(r$analysis_type %in% c("demoWAS", "exWAS"),
                                 "genus is DEPENDENT (genus_CLR ~ phenotype)",
                                 "genus is INDEPENDENT (phenotype ~ genus_CLR)"),
      external_dataset = r$dataset_id, tier = r$tier, fit = r$fit,
      label = r$label, primary = r$primary,
      ext_genus, ext_beta = ext_estimate, ext_se, ext_p, ext_fdr, ext_n,
      ext_direction,
      nhanes_covariates = nh_cov, external_covariates = ext_cov,
      covariates_identical = adj == "MATCHED",
      covariates_partially_comparable = adj == "PARTIALLY_MATCHED",
      important_missing_external_covariates = ifelse(
        r$tier == "B", "unknown - published result",
        "PIR, education, ethnicity, US-born (not measured externally)"),
      adjustment_comparability = adj,
      oral_site = cm$oral_site, sequencing_type = cm$sequencing_type,
      taxonomy_method = cm$taxonomy_method,
      orientation_transform = sl$orientation_transform[1],
      signs_comparable = sl$signs_comparable[1],
      # ---- evidence levels ------------------------------------------------
      tested = !is.na(ext_direction),
      same_direction = ifelse(is.na(ext_direction), NA, ext_direction == direction),
      external_nominal_significant = ifelse(is.na(ext_p), NA, ext_p < 0.05),
      external_fdr_significant = ifelse(is.na(ext_fdr), NA, ext_fdr < EXTERNAL_FDR),
      significant_same_direction = ifelse(is.na(ext_fdr) | is.na(ext_direction), NA,
                                          ext_fdr < EXTERNAL_FDR & ext_direction == direction))
}

mt <- bind_rows(rows)

# --- Validation classification (section 12, controlled vocabulary) -----------
mt <- mt %>%
  mutate(validation_classification = case_when(
    !tested & tier == "B" ~ "NOT_TESTABLE_TAXON",     # absent from the published list
    !tested               ~ "NOT_TESTABLE_TAXON",     # genus not measured in the cohort
    same_direction & !is.na(external_fdr_significant) & external_fdr_significant ~
      "FDR_REPLICATED_SAME_DIRECTION",
    same_direction & !is.na(external_nominal_significant) & external_nominal_significant ~
      "NOMINALLY_REPLICATED_SAME_DIRECTION",
    same_direction ~ "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT",
    !same_direction ~ "DIRECTIONALLY_DISCORDANT",
    TRUE ~ "MODEL_NOT_COMPARABLE"),
    limitations = paste0(
      "site: ", oral_site, "; seq: ", sequencing_type,
      "; adjustment: ", adjustment_comparability,
      ifelse(tier == "B", "; published direction only (no beta/SE)", "")))

stopifnot(all(mt$validation_classification %in% VALIDATION_CLASSES))

# --- Assertions --------------------------------------------------------------
# 1. Three-state partition: every NHANES hit is tested-and-classified or
#    explicitly not testable, and the parts sum to the hit count.
chk <- mt %>% group_by(nhanes_variable, external_dataset, fit) %>%
  summarise(n_rows = n(), n_hits = n_distinct(genus), .groups = "drop")
stopifnot(all(chk$n_rows == chk$n_hits))

# 2. No published-target "discordant" without a published direction.
badB <- mt %>% filter(tier == "B", validation_classification == "DIRECTIONALLY_DISCORDANT",
                      is.na(ext_direction))
stopifnot(nrow(badB) == 0)

# 3. Participant N and taxon N must never be conflated: ext_n is participants.
stopifnot(!"n_taxa" %in% names(mt))

write_tsv(mt, file.path(AUD, "master_validation_table.tsv"))
asay("wrote master_validation_table.tsv (", format(nrow(mt), big.mark = ","), " rows)")
asay("[OK] partition, published-target and N-conflation assertions passed")

asay("")
asay("=== classification by phenotype x cohort (primary fits) ===")
print(as.data.frame(mt %>% filter(primary) %>%
  count(phenotype_label, label, validation_classification) %>%
  tidyr::pivot_wider(names_from = validation_classification, values_from = n, values_fill = 0)))
