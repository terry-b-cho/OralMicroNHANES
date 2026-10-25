#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 8
# Tables 1-6, generated FROM master_validation_table.tsv
#
# Every percentage carries its denominator; participant N and taxon N are kept
# in separate, differently-named columns throughout.
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 8: tables 1-6 ===")

mt   <- read_tsv(file.path(AUD, "master_validation_table.tsv"), show_col_types = FALSE)
conc <- read_tsv(file.path(AUD, "concordance_by_phenotype.tsv"), show_col_types = FALSE)
opp  <- read_tsv(file.path(AUD, "validation_opportunity_matrix.tsv"), show_col_types = FALSE)
P <- mt %>% filter(primary)

flat <- function(df) dplyr::mutate(df, dplyr::across(dplyr::where(is.character),
  ~ trimws(gsub(" {2,}", " ", gsub("[\r\n\t]+", " ", .x)))))
pc <- function(k, n) ifelse(n == 0, "-", sprintf("%d/%d (%.0f%%)", k, n, 100 * k / n))

# --- Table 1: headline summary ----------------------------------------------
t1 <- conc %>% transmute(
  `External cohort` = label, `Cohort class` = m12_cohort_class(tier), Phenotype = phenotype_label,
  `NHANES sig+prevalent assoc (n taxa)` = n_nhanes_primary,
  `Externally testable (n taxa)` = n_testable,
  `Same direction` = pc(n_same_direction, n_testable),
  `95% CI` = ifelse(is.na(ci_lo), "-", sprintf("%.0f-%.0f", ci_lo, ci_hi)),
  `Binomial p` = ifelse(is.na(binom_p), "-", format.pval(binom_p, digits = 2)),
  `External nominal p<0.05` = pc(n_nominal_significant, n_testable),
  `External FDR<0.05` = pc(n_fdr_significant, n_testable),
  `FDR sig + same direction` = pc(n_significant_same_dir, n_testable),
  `Discordant` = pc(n_discordant, n_testable),
  `Not testable (n taxa)` = n_not_testable,
  `NHANES participants` = n_participants_nhanes,
  `External participants` = ifelse(is.na(n_participants_external), "n/a (published table)",
                                   format(n_participants_external, big.mark = ",")),
  `Oral site` = oral_site, `Adjustment comparability` = adjustment_comparability) %>%
  arrange(Phenotype, desc(`Externally testable (n taxa)`))
write_tsv(flat(t1), file.path(AUD, "table1_headline_summary.tsv"))

# --- Table 2: strongest replicated findings ---------------------------------
t2 <- P %>% filter(validation_classification %in%
                     c("FDR_REPLICATED_SAME_DIRECTION", "NOMINALLY_REPLICATED_SAME_DIRECTION")) %>%
  arrange(validation_classification, nhanes_fdr) %>%
  transmute(Phenotype = phenotype_label, Genus = genus, Phylum = phylum,
            `NHANES prevalence (%)` = round(100 * nhanes_prevalence, 1),
            `NHANES beta` = signif(nhanes_beta, 3), `NHANES FDR` = signif(nhanes_fdr, 3),
            `NHANES q` = signif(nhanes_q, 3),
            `External beta` = signif(ext_beta, 3),
            `External p` = signif(ext_p, 3), `External FDR` = signif(ext_fdr, 3),
            `External cohort` = label,
            `NHANES participants` = nhanes_n, `External participants` = ext_n,
            `NHANES covariates` = nhanes_covariates, `External covariates` = external_covariates,
            `Adjustment comparability` = adjustment_comparability,
            `Validation class` = validation_classification)
write_tsv(flat(t2), file.path(AUD, "table2_strongest_replicated.tsv"))

# --- Table 3: directionally concordant but underpowered ---------------------
t3 <- P %>% filter(validation_classification == "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT") %>%
  arrange(phenotype_label, nhanes_fdr) %>%
  transmute(Phenotype = phenotype_label, Genus = genus, Phylum = phylum,
            `NHANES prevalence (%)` = round(100 * nhanes_prevalence, 1),
            `NHANES beta` = signif(nhanes_beta, 3), `NHANES FDR` = signif(nhanes_fdr, 3),
            `External beta` = signif(ext_beta, 3),
            `External 95% CI` = ifelse(is.na(ext_se), "n/a (published direction only)",
              sprintf("%.3f to %.3f", ext_beta - 1.96 * ext_se, ext_beta + 1.96 * ext_se)),
            `External p` = signif(ext_p, 3),
            `External participants` = ifelse(is.na(ext_n), NA_integer_, ext_n),
            `External cohort` = label,
            Interpretation = "same direction; external cohort is smaller and does not reach significance - NOT evidence against the NHANES finding")
write_tsv(flat(t3), file.path(AUD, "table3_concordant_underpowered.tsv"))

# --- Table 4: discordant findings, with their comparability differences -----
t4 <- P %>% filter(validation_classification == "DIRECTIONALLY_DISCORDANT") %>%
  arrange(phenotype_label, nhanes_fdr) %>%
  transmute(Phenotype = phenotype_label, Genus = genus, Phylum = phylum,
            `NHANES prevalence (%)` = round(100 * nhanes_prevalence, 1),
            `NHANES beta` = signif(nhanes_beta, 3), `NHANES direction` = nhanes_direction,
            `NHANES FDR` = signif(nhanes_fdr, 3),
            `External beta` = signif(ext_beta, 3), `External direction` = ext_direction,
            `External p` = signif(ext_p, 3), `External FDR` = signif(ext_fdr, 3),
            `External cohort` = label, `External participants` = ext_n,
            `Oral site (external)` = oral_site, `Sequencing (external)` = sequencing_type,
            `Adjustment comparability` = adjustment_comparability,
            `Phenotype definition gap` = dplyr::case_when(
              nhanes_variable == "GUM_DISEASE_OHAROCGP" ~
                "self-reported gum disease vs clinically diagnosed peri-implantitis",
              nhanes_variable == "EDUCATION_COLLEGEGRAD" ~
                "US census education categories vs a volunteer cohort skewed to graduate degrees",
              nhanes_variable == "BMXBMI" ~ "measured BMI vs partly self-reported BMI",
              nhanes_variable == "RIAGENDR" ~ "same construct; the gap is oral site, not definition",
              TRUE ~ "same construct"))
write_tsv(flat(t4), file.path(AUD, "table4_discordant.tsv"))

# --- Table 5: cross-study validation opportunities ---------------------------
t5 <- opp %>%
  filter(n_nhanes_fdr_significant_prevalence_gt_10pct > 0) %>%
  transmute(`Analysis` = analysis_type, `NHANES variable` = nhanes_variable,
            `Label` = nhanes_variable_label,
            `Sig+prevalent assoc (n taxa)` = n_nhanes_fdr_significant_prevalence_gt_10pct,
            `External cohort` = external_dataset, `Cohort class` = m12_cohort_class(tier),
            `External variable` = external_variable_harmonized,
            `Comparison status` = comparison_status,
            `Already performed` = validation_analysis_already_performed,
            `Additional analysis possible` = additional_analysis_possible,
            Recommendation = recommendation,
            `Reason if not possible` = reason_if_not_possible) %>%
  arrange(desc(`Sig+prevalent assoc (n taxa)`), `NHANES variable`, `External cohort`)
write_tsv(flat(t5), file.path(AUD, "table5_validation_opportunities.tsv"))

# --- Table 6: covariate comparability ---------------------------------------
t6 <- P %>% distinct(phenotype_label, analysis_type, label, tier, .keep_all = TRUE) %>%
  transmute(Phenotype = phenotype_label, `NHANES scheme` = analysis_type,
            `NHANES model` = nhanes_model_role, `NHANES covariates` = nhanes_covariates,
            `External cohort` = label, `External model` = ifelse(tier == "A",
              "genus_CLR ~ exposure (+ covariates), OLS, unweighted", "published model (see source)"),
            `External covariates` = external_covariates,
            `Missing external covariates` = important_missing_external_covariates,
            `Covariates identical` = covariates_identical,
            `Adjustment comparability` = adjustment_comparability,
            `Orientation transform` = orientation_transform,
            `Signs comparable` = signs_comparable) %>%
  arrange(Phenotype, `External cohort`)
write_tsv(flat(t6), file.path(AUD, "table6_covariate_comparability.tsv"))

for (n in c("table1_headline_summary", "table2_strongest_replicated",
            "table3_concordant_underpowered", "table4_discordant",
            "table5_validation_opportunities", "table6_covariate_comparability")) {
  f <- file.path(AUD, paste0(n, ".tsv"))
  asay("wrote ", n, ".tsv (", nrow(read_tsv(f, show_col_types = FALSE, progress = FALSE)), " rows)")
}
asay("")
asay("=== Table 1 ===")
print(as.data.frame(t1 %>% select(Phenotype, `External cohort`, `Cohort class`,
  testable = `Externally testable (n taxa)`, `Same direction`,
  `FDR sig + same direction`, Discordant)))
