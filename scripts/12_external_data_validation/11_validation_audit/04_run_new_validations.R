#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 4
# Run the validation analyses the audit identified as possible but unperformed
#
# Four opportunities survived the opportunity matrix:
#   RUN                  sex        Huang 2020 (refit) - 55 NHANES primary hits
#   RUN_WITH_LIMITATIONS BMI        Huang 2020 (refit) -  4 hits, orientation differs
#   RUN_WITH_LIMITATIONS education  Huang 2020 (refit) - 47 hits, 36% coverage
#   RUN_WITH_LIMITATIONS gum disease Ghensi 2020 (published-target) - 56 hits, published target
#
# SIGN ALIGNMENT (work order section 17) is handled here, once, and recorded in
# a machine-readable field. It is NOT assumed:
#   * NHANES demoWAS/exWAS: genus_CLR ~ phenotype -> beta is "genus per unit
#     phenotype". Huang is fitted the same way, so no transformation is needed.
#   * NHANES oradWAS/pheWAS: phenotype ~ genus_CLR -> beta is "phenotype per
#     unit genus". The external comparators (Huang BMI fit, Ghensi's published
#     disease direction) express "genus higher in/with the phenotype". Both
#     encode the same co-movement, so SIGNS are directly comparable while
#     MAGNITUDES are on different scales and are never compared.
#   * RIAGENDR2 is FEMALE vs male reference (verified against the project's own
#     recode in scripts/0_transform_n_preprocess_ssfiles). Huang's sex factor is
#     releveled to male so the two betas mean the same thing.
#
# Usage   : Rscript 04_run_new_validations.R
# Outputs : validation/external_effects_<dataset>_<variable>_<fit>.tsv
#           audit/new_fit_audit.tsv
#           audit/sign_alignment_log.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))

# Reuse Phase 2's fitting machinery verbatim - one implementation, not two.
Sys.setenv(M12_HELPERS_ONLY = "1", M12_V2_DIR = V2)
audit <- list()
source(file.path(V2, "02_reanalyse_external.R"))
Sys.unsetenv("M12_HELPERS_ONLY"); Sys.unsetenv("M12_V2_DIR")
asay("=== Step 4: run the newly identified validations ===")

# --- Exposure builders for the new phenotypes --------------------------------

# Sex. RIAGENDR2 = Female vs Male reference, so Huang is releveled to male.
# A positive beta then means "genus CLR higher in females" in BOTH cohorts.
sex_exposure <- function(meta) {
  s <- tolower(trimws(as.character(meta$sex_harmonized)))
  s[!s %in% c("male", "female")] <- NA
  meta$sex_f <- factor(s, levels = c("male", "female"))   # male = reference
  list(meta = meta, exposure = "sex_f",
       describe = sprintf("female=%d vs male=%d (male is the reference, matching RIAGENDR2)",
                          sum(meta$sex_f == "female", na.rm = TRUE),
                          sum(meta$sex_f == "male", na.rm = TRUE)))
}

# BMI, continuous. NHANES pheWAS fits BMI ~ genus; we fit genus ~ BMI. Signs
# are comparable, magnitudes are not.
bmi_exposure <- function(meta) {
  b <- suppressWarnings(as.numeric(meta$bmi))
  b[!is.na(b) & (b < 10 | b > 80)] <- NA          # implausible values
  meta$bmi_num <- b
  list(meta = meta, exposure = "bmi_num",
       describe = sprintf("continuous, n=%d, median %.1f kg/m2 (range %.1f-%.1f)",
                          sum(!is.na(b)), median(b, na.rm = TRUE),
                          min(b, na.rm = TRUE), max(b, na.rm = TRUE)))
}

# Education collapsed to the NHANES EDUCATION_COLLEGEGRAD contrast:
# college graduate or above vs everything below.
edu_exposure <- function(meta) {
  e <- trimws(as.character(meta$education))
  grad <- c("Bachelor's degree", "Graduate or Professional degree",
            "Some graduate school or professional")
  below <- c("Did not complete high school", "High School or GED equilivant",
             "Some college or technical school", "Associate's degree")
  v <- ifelse(e %in% grad, "collegegrad", ifelse(e %in% below, "below", NA))
  meta$edu_f <- factor(v, levels = c("below", "collegegrad"))
  list(meta = meta, exposure = "edu_f",
       describe = sprintf("collegegrad=%d vs below=%d",
                          sum(meta$edu_f == "collegegrad", na.rm = TRUE),
                          sum(meta$edu_f == "below", na.rm = TRUE)))
}

HU <- "HuangS_2020_AgePrediction"

# Each fitted twice where NHANES itself is crude, so the like-for-like model is
# available AND the adjustment sensitivity is measurable.
run_dataset(HU, "sex", "relative", TRUE, TRUE, sex_exposure,
            list(unadjusted = character(0), adjusted = c("age_years")))

run_dataset(HU, "bmi", "relative", TRUE, TRUE, bmi_exposure,
            list(unadjusted = character(0), adjusted = c("age_years", "sex_harmonized")))

run_dataset(HU, "education", "relative", TRUE, TRUE, edu_exposure,
            list(unadjusted = character(0), adjusted = c("age_years", "sex_harmonized")))

new_aud <- bind_rows(audit)
write_tsv(new_aud, file.path(AUD, "new_fit_audit.tsv"))
asay("")
asay("=== new refit-cohort fits ===")
print(as.data.frame(new_aud))

# --- Sign-alignment log (machine readable, work order section 17) ------------
sign_log <- tibble::tribble(
  ~nhanes_variable, ~analysis_type, ~external_dataset, ~nhanes_positive_beta_means, ~external_positive_beta_means, ~orientation_transform, ~signs_comparable,
  "RIDAGEYR", "demoWAS", HU,
  "genus CLR increases per year of age", "genus CLR increases per year of age",
  "none - both model genus ~ phenotype", TRUE,
  "RIAGENDR", "demoWAS", HU,
  "genus CLR higher in FEMALES (RIAGENDR2 = female, ref male)",
  "genus CLR higher in FEMALES (factor releveled so male is the reference)",
  "external factor releveled to male reference; no sign flip applied", TRUE,
  "BMXBMI", "pheWAS", HU,
  "BMI increases per unit genus CLR (genus is the PREDICTOR)",
  "genus CLR increases per unit BMI (genus is the OUTCOME)",
  "orientation differs; both encode positive co-movement so signs are comparable, magnitudes are NOT",
  TRUE,
  "EDUCATION_COLLEGEGRAD", "demoWAS", HU,
  "genus CLR higher in college graduates",
  "genus CLR higher in college graduates (below = reference)",
  "external factor releveled so 'below' is the reference", TRUE,
  "GUM_DISEASE_OHAROCGP", "oradWAS", "GhensiP_2020_ImplantDisease",
  "gum-disease probability increases per unit genus CLR (genus is the PREDICTOR)",
  "genus is enriched in the peri-implant DISEASE class",
  "orientation differs; both encode 'genus up with disease' so signs are comparable, magnitudes are NOT",
  TRUE,
  "SMQ_current_ever_never", "exWAS", "WillisJR_2022_SacaLaLengua2",
  "genus CLR higher in current vs never smokers",
  "genus reported as increased in smokers", "none", TRUE,
  "RIDAGEYR", "demoWAS", "WillisJR_2022_SacaLaLengua2",
  "genus CLR increases per year of age",
  "genus reported as increasing with age in the published table",
  "none - published direction already expressed as 'with age'", TRUE,
  "RIDAGEYR", "demoWAS", "WellsPM_2022_TwinsUKSaliva",
  "genus CLR increases per year of age",
  "genus reported as increasing with age in the published figure",
  "none - published direction already expressed as 'with age'", TRUE)
write_tsv(sign_log, file.path(AUD, "sign_alignment_log.tsv"))
asay("")
asay("wrote sign_alignment_log.tsv (", nrow(sign_log), " comparisons, all sign-aligned before use)")
