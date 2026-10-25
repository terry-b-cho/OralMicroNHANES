#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 3
# The NHANES -> external validation-opportunity matrix
#
# Purpose : One row per (OralMicroNHANES variable x ACQUIRED external dataset)
#           for EVERY variable carrying at least one primary association, so an
#           opportunity cannot be missed merely because nobody thought to look.
#           Curated mappings are declared explicitly and justified; every other
#           cell gets a computed, specific reason for non-testability.
#
# Usage   : Rscript 03_build_opportunity_matrix.R
# Outputs : audit/validation_opportunity_matrix.tsv
#           audit/additional_validation_opportunities.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
suppressPackageStartupMessages({library(DBI); library(RSQLite)})
asay("=== Step 3: validation-opportunity matrix ===")

inv <- read_tsv(file.path(AUD, "nhanes_variable_inventory.tsv"), show_col_types = FALSE)
cap <- read_tsv(file.path(AUD, "external_dataset_capability.tsv"), show_col_types = FALSE)
dict<- read_tsv(file.path(AUD, "external_phenotype_dictionary.tsv"), show_col_types = FALSE)

# --- NHANES variable labels --------------------------------------------------
# Raw NHANES names resolve from the shipped codebook; project-derived variables
# (EDUCATION_*, ETHNICITY_*, SMQ_current_ever_never, the *_OHAROC* oral-health
# recodes) do not exist there and are labelled from their construction.
con <- dbConnect(SQLite(), file.path(PROJECT_ROOT,
  "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"))
cb <- dbGetQuery(con, "SELECT `Variable.Name` v, `Variable.Description` d, Component c
                       FROM variable_names_epcf") %>%
  filter(!is.na(d), d != "") %>% distinct(v, .keep_all = TRUE)
dbDisconnect(con)

DERIVED_LABELS <- c(
  RIAGENDR = "Sex (RIAGENDR2 = Female; reference Male)",
  AGE_SQUARED = "Age squared (years^2)",
  BORN_INUSA = "Born in the USA (yes vs no)",
  INDFMPIR = "Family income-to-poverty ratio",
  SMQ_current_ever_never = "Smoking status, 3-level (term ...2 = current vs never)",
  EDUCATION_LESS9 = "Education: less than 9th grade",
  EDUCATION_9_11 = "Education: 9th-11th grade",
  EDUCATION_AA = "Education: associate degree / some college",
  EDUCATION_COLLEGEGRAD = "Education: college graduate or above",
  ETHNICITY_MEXICAN = "Ethnicity: Mexican American",
  ETHNICITY_OTHERHISPANIC = "Ethnicity: Other Hispanic",
  ETHNICITY_NONHISPANICBLACK = "Ethnicity: Non-Hispanic Black",
  ETHNICITY_OTHER = "Ethnicity: Other/multi-racial",
  GUM_DISEASE_OHAROCGP = "Self-reported gum/periodontal disease",
  ORAL_HYGIENE_OHAROCOH = "Oral hygiene status recode",
  TOOTH_DECAY_OHAROCDT = "Tooth decay recode",
  DENTURE_OHAROCDE = "Denture wear recode")

lab_of <- function(v) {
  out <- DERIVED_LABELS[v]
  hit <- cb$d[match(v, cb$v)]
  out[is.na(out)] <- hit[is.na(out)]
  ifelse(is.na(out), v, out)
}

# =============================================================================
# CURATED MAPPINGS
#
# Each row is a scientific judgement about whether an external variable is a
# defensible analogue of a NHANES variable. Declared explicitly, with the
# reason, so a reader can disagree with a specific call rather than with an
# opaque automatic match. Anything absent from this table is NOT claimed to be
# comparable.
# =============================================================================
MAPPING <- tibble::tribble(
  ~nhanes_variable, ~analysis_type, ~external_dataset, ~external_variable_original, ~external_variable_harmonized, ~external_definition, ~external_units, ~external_coding, ~comparability, ~already_done, ~recommendation, ~rationale,

  # ---- Refit cohort: Huang is the only cohort with a genus x sample table ---------
  "RIDAGEYR", "demoWAS", "HuangS_2020_AgePrediction", "age", "age_years",
  "Chronological age at sampling", "years", "continuous",
  "DIRECTLY_COMPARABLE", TRUE, "RUN",
  "Same construct, continuous in both. Validated in Phase 2.",

  "RIAGENDR", "demoWAS", "HuangS_2020_AgePrediction", "sex", "sex_harmonized",
  "Reported sex", "binary", "female vs male (male = reference, matching RIAGENDR2)",
  "DIRECTLY_COMPARABLE", FALSE, "RUN",
  "Same binary construct and the SAME reference level, so betas are sign-aligned without transformation. 55 NHANES primary hits; Huang has sex at 98.7%. The single largest missed opportunity.",

  "BMXBMI", "pheWAS", "HuangS_2020_AgePrediction", "bmi", "bmi",
  "Body mass index", "kg/m2", "continuous",
  "COMPARABLE_WITH_LIMITATIONS", FALSE, "RUN_WITH_LIMITATIONS",
  "Same construct, but NHANES pheWAS models BMI ~ genus while the external fit is genus ~ BMI. Signs are comparable (both mean 'genus and BMI move together'); magnitudes are NOT and are never compared. Huang BMI is 66% complete and self-reported in several contributing studies.",

  "EDUCATION_COLLEGEGRAD", "demoWAS", "HuangS_2020_AgePrediction", "education", "education",
  "Educational attainment", "categorical", "college graduate or above vs below",
  "COMPARABLE_WITH_LIMITATIONS", FALSE, "RUN_WITH_LIMITATIONS",
  "Collapsible to the same binary contrast, but Huang education is only 36% complete (264 adults) and severely skewed to graduate degrees (128/264) because the contributing cohorts are research volunteers. Low power and non-representative.",

  "ETHNICITY_NONHISPANICBLACK", "demoWAS", "HuangS_2020_AgePrediction", "race", "race",
  "Race/ethnicity", "categorical", "African American vs Caucasian",
  "NOT_COMPARABLE", FALSE, "DO_NOT_RUN",
  "Only 8 African American adults against 439 Caucasian. NHANES ethnicity dummies are US-census categories; Huang pools 22 countries. Cell sizes cannot support the contrast.",

  "SMQ_current_ever_never", "exWAS", "HuangS_2020_AgePrediction", "smoking", "smoking_status_harmonized",
  "Smoking status", "categorical", "current vs never",
  "NOT_COMPARABLE", FALSE, "DO_NOT_RUN",
  "Only 5 current smokers among adults. Documented in Phase 1 as unusable.",

  # ---- Published-target cohort: published significant-genus targets --------------------------
  "RIDAGEYR", "demoWAS", "WillisJR_2022_SacaLaLengua2", "Age", "published Table 1",
  "Age-associated genera, published directional table", "n/a", "direction only",
  "COMPARABLE_WITH_LIMITATIONS", TRUE, "RUN",
  "Published significant-genus list only; ~80% adolescent cohort. Validated in Phase 2.",

  "RIDAGEYR", "demoWAS", "WellsPM_2022_TwinsUKSaliva", "Age", "published Fig 5",
  "Age-associated genera read from a published figure", "n/a", "direction only",
  "COMPARABLE_WITH_LIMITATIONS", TRUE, "RUN",
  "Direction transcribed from a bitmap figure; adults. Validated in Phase 2.",

  "SMQ_current_ever_never", "exWAS", "WillisJR_2022_SacaLaLengua2", "Tobacco", "published Table 2",
  "Smoking-associated genera, published directional table", "n/a", "direction only",
  "COMPARABLE_WITH_LIMITATIONS", TRUE, "RUN",
  "Binary smoking ('No' conflates never+former) vs NHANES 3-level. Validated in Phase 2.",

  "GUM_DISEASE_OHAROCGP", "oradWAS", "GhensiP_2020_ImplantDisease", "peri-implantitis vs healthy", "published disease target",
  "Genera enriched in peri-implant disease vs healthy implants", "n/a", "direction in positive (disease) class",
  "COMPARABLE_WITH_LIMITATIONS", FALSE, "RUN_WITH_LIMITATIONS",
  "Both are oral inflammatory dysbiosis with the genus as predictor of disease, so signs align without transformation. But the site differs sharply: Ghensi is SUBGINGIVAL PLAQUE around implants (n=113, 23 peri-implantitis), NHANES is whole-mouth oral rinse with self-reported gum disease. A different ecological niche and a different case definition.",

  "TOOTH_DECAY_OHAROCDT", "oradWAS", "GhensiP_2020_ImplantDisease", "peri-implantitis vs healthy", "published disease target",
  "Genera enriched in peri-implant disease", "n/a", "direction in positive class",
  "NOT_COMPARABLE", FALSE, "DO_NOT_RUN",
  "Caries and peri-implant disease are distinct pathologies with different keystone taxa; equating them would manufacture agreement.",

  "ORAL_HYGIENE_OHAROCOH", "oradWAS", "GhensiP_2020_ImplantDisease", "peri-implantitis vs healthy", "published disease target",
  "Genera enriched in peri-implant disease", "n/a", "direction in positive class",
  "NOT_COMPARABLE", FALSE, "DO_NOT_RUN",
  "Oral-hygiene status is a behavioural recode, not a disease diagnosis. Not a defensible analogue."
)
stopifnot(all(MAPPING$comparability %in% COMPARABILITY_LEVELS),
          all(MAPPING$recommendation %in% RECOMMENDATIONS))

# --- Full matrix: every primary variable x every ACQUIRED dataset ------------
vars <- inv %>% filter(n_fdr_sig_prev_gt10 > 0)
asay("variables with >=1 primary association: ", nrow(vars))
asay("acquired datasets: ", nrow(cap))

grid <- tidyr::crossing(vars %>% select(analysis_type, genus_role, nhanes_variable,
                                        n_genus_tested, n_fdr_significant,
                                        n_fdr_sig_prev_gt10, n_covariates, n_obs_median),
                        cap %>% select(external_dataset, tier, n_samples,
                                       has_abundance_table, n_published_targets))

mm <- grid %>%
  left_join(MAPPING, by = c("nhanes_variable", "analysis_type", "external_dataset")) %>%
  mutate(
    nhanes_variable_label = lab_of(nhanes_variable),
    nhanes_variable_definition = sprintf("%s [%s; genus %s; %s]",
      lab_of(nhanes_variable), analysis_type, genus_role,
      ifelse(n_covariates == 0, "crude", paste0(n_covariates, " covariates"))),
    nhanes_model_role = ifelse(genus_role == "dependent",
                               "genus is DEPENDENT (genus_CLR ~ phenotype)",
                               "genus is INDEPENDENT (phenotype ~ genus_CLR)"),
    nhanes_transformation = "CLR (ln-ratio, centered; pseudocount 1)",
    nhanes_covariate_schema = ifelse(n_covariates == 0, "NONE (crude)",
      sprintf("%d covariates (age, age^2, sex, PIR, education x4, ethnicity x4, US-born)", n_covariates)),
    covariates_available = case_when(
      is.na(comparability) ~ NA_character_,
      external_dataset == "HuangS_2020_AgePrediction" ~ "age, sex, BMI (partial), country, body site",
      external_dataset == "GhensiP_2020_ImplantDisease" ~ "published model only: age, sex, smoking per source paper",
      TRUE ~ "published result - covariates as reported by the source"),
    oral_site_comparable = case_when(
      external_dataset == "HuangS_2020_AgePrediction" ~ "NO - 73.5% tongue-dorsum swab vs NHANES oral rinse",
      external_dataset == "GhensiP_2020_ImplantDisease" ~ "NO - subgingival plaque vs NHANES oral rinse",
      external_dataset == "WillisJR_2022_SacaLaLengua2" ~ "YES - whole-mouth rinse",
      external_dataset == "WellsPM_2022_TwinsUKSaliva" ~ "PARTIAL - saliva vs oral rinse",
      TRUE ~ NA_character_),
    sequencing_type = case_when(
      external_dataset == "HuangS_2020_AgePrediction" ~ "16S V4 (515F/806R) - matches NHANES",
      external_dataset == "GhensiP_2020_ImplantDisease" ~ "shotgun metagenomics",
      external_dataset == "WillisJR_2022_SacaLaLengua2" ~ "16S V3-V4",
      external_dataset == "WellsPM_2022_TwinsUKSaliva" ~ "16S V4",
      TRUE ~ NA_character_),
    taxonomy_level = "genus",
    microbiome_available = has_abundance_table | n_published_targets > 0,
    phenotype_available = !is.na(external_variable_harmonized),
    sample_mapping_available = has_abundance_table,
    validation_analysis_already_performed = tidyr::replace_na(already_done, FALSE),
    # Anything not curated gets a SPECIFIC computed reason, never a blank.
    comparison_status = case_when(
      !is.na(comparability) ~ comparability,
      !has_abundance_table & n_published_targets == 0 ~ "METADATA_UNAVAILABLE",
      !has_abundance_table ~ "METADATA_UNAVAILABLE",
      TRUE ~ "NOT_COMPARABLE"),
    additional_analysis_possible = !validation_analysis_already_performed &
      comparison_status %in% c("DIRECTLY_COMPARABLE", "COMPARABLE_WITH_LIMITATIONS"),
    reason_if_not_possible = case_when(
      comparison_status %in% c("DIRECTLY_COMPARABLE", "COMPARABLE_WITH_LIMITATIONS") ~ NA_character_,
      !is.na(rationale) ~ rationale,
      !has_abundance_table & n_published_targets == 0 ~
        "No genus x sample abundance table and no published per-taxon target: nothing to compare against.",
      !has_abundance_table ~
        sprintf("No abundance table; the %d published target(s) cover other phenotypes, so this variable cannot be tested here.",
                n_published_targets),
      TRUE ~ "No comparable phenotype measured in this cohort."),
    major_limitations = dplyr::coalesce(rationale, reason_if_not_possible),
    recommendation = dplyr::coalesce(recommendation, "DO_NOT_RUN")) %>%
  # label and definition are BOTH required by the spec and must be distinct
  # columns; selecting one column twice silently keeps only the last.
  select(analysis_type, nhanes_variable, nhanes_variable_label,
         nhanes_variable_definition,
         nhanes_model_role, nhanes_transformation, nhanes_covariate_schema,
         n_nhanes_genus_associations = n_genus_tested,
         n_nhanes_fdr_significant = n_fdr_significant,
         n_nhanes_fdr_significant_prevalence_gt_10pct = n_fdr_sig_prev_gt10,
         external_dataset, tier, external_variable_original,
         external_variable_harmonized, external_definition, external_units,
         external_coding, external_sample_size = n_samples,
         microbiome_available, phenotype_available, sample_mapping_available,
         covariates_available, oral_site_comparable, sequencing_type, taxonomy_level,
         comparison_status, validation_analysis_already_performed,
         additional_analysis_possible, recommendation,
         reason_if_not_possible, major_limitations)

write_tsv(mm, file.path(AUD, "validation_opportunity_matrix.tsv"))
asay("wrote validation_opportunity_matrix.tsv (", format(nrow(mm), big.mark = ","), " rows)")

# --- The actionable subset ---------------------------------------------------
opp <- mm %>%
  filter(comparison_status %in% c("DIRECTLY_COMPARABLE", "COMPARABLE_WITH_LIMITATIONS") |
           recommendation != "DO_NOT_RUN" | !is.na(external_variable_harmonized)) %>%
  arrange(validation_analysis_already_performed,
          desc(n_nhanes_fdr_significant_prevalence_gt_10pct))
write_tsv(opp, file.path(AUD, "additional_validation_opportunities.tsv"))

asay("")
asay("=== every curated NHANES <-> external mapping ===")
print(as.data.frame(opp %>% select(analysis_type, nhanes_variable, external_dataset,
                                   external_variable_harmonized,
                                   n = n_nhanes_fdr_significant_prevalence_gt_10pct,
                                   comparison_status, recommendation,
                                   done = validation_analysis_already_performed)))
asay("")
asay("NEW opportunities not yet performed: ",
     sum(opp$additional_analysis_possible))
asay("")
asay("=== why the other ", format(nrow(mm) - nrow(opp), big.mark = ","),
     " variable x dataset cells are not testable ===")
print(as.data.frame(mm %>% filter(!nhanes_variable %in% opp$nhanes_variable) %>%
  count(comparison_status, reason_if_not_possible, name = "cells") %>% arrange(desc(cells))))
