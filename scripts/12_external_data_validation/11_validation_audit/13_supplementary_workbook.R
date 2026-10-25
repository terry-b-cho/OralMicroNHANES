#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 13
# Supplementary workbook: one file, three sheets, replacing six tables
#
# Purpose : Build S10_summary_of_external_concordance_validation_data_and_results.xlsx
#           with sheets `cohorts`, `summary` and `per genus`.
#
#           A reviewer asked us to cut the supplementary table count and widen
#           what survives. Old S3/S4/S5 (table2/table3/table4) were the same
#           column set over three DISJOINT row subsets of one per-genus space -
#           25 + 79 + 61 rows, pairwise key overlap exactly zero. They merge
#           losslessly into one sheet keyed by phenotype x genus x cohort, and
#           the merged sheet additionally carries the rows those three never
#           showed: the NHANES hits no cohort could test.
#
#           Old S2/S6 (table1/table6) share a grain - one row per phenotype x
#           cohort - so they join rather than stack.
#
#           Cohort-level attributes (site, assay, covariates) were repeated on
#           every per-genus row. They are lifted onto their own sheets, but only
#           after asserting they really are constant within a comparison; that
#           assertion is what makes the lift lossless rather than merely
#           plausible.
#
# Invariant: every count on `summary` is reproducible by filtering `per genus`.
#           That is the point of the workbook, so it is asserted, not assumed.
#
# Usage   : Rscript 13_supplementary_workbook.R
# Env     : R 4.5.1 (conda nhanes-analysis) + system python3 with openpyxl
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 13: supplementary workbook ===")

OUT_XLSX <- file.path(M12_OUT,
  "S10_summary_of_external_concordance_validation_data_and_results.xlsx")
# Output guard (module convention): never write outside the module tree.
stopifnot(grepl("12_external_data_validation", OUT_XLSX, fixed = TRUE))

# --- Sources -----------------------------------------------------------------
mt    <- read_tsv(file.path(AUD, "master_validation_table.tsv"),   show_col_types = FALSE)
conc  <- read_tsv(file.path(AUD, "concordance_by_phenotype.tsv"),  show_col_types = FALSE)
bcor  <- read_tsv(file.path(AUD, "beta_correlations.tsv"),         show_col_types = FALSE)
enr   <- read_tsv(file.path(AUD, "enrichment_tests.tsv"),          show_col_types = FALSE)
phy   <- read_tsv(file.path(AUD, "phylum_level_validation.tsv"),   show_col_types = FALSE)
alp   <- read_tsv(file.path(AUD, "alpha_diversity_validation.tsv"),show_col_types = FALSE)
sc    <- read_tsv(file.path(VAL, "study_characteristics.tsv"),     show_col_types = FALSE)
ovl   <- read_tsv(file.path(M12_OUT, "genus_overlap_summary.tsv"), show_col_types = FALSE)

# Every downstream aggregate in this module is computed on the primary fits;
# the sensitivity fits (Huang adj. sex / adj. age) are excluded so the workbook
# cannot disagree with the reports.
P <- mt %>% filter(primary)
asay("master: ", nrow(mt), " rows -> ", nrow(P), " primary")

# Excel cells must not contain tabs or newlines, and the TSV handoff would
# corrupt on either.
flat <- function(df) dplyr::mutate(df, dplyr::across(dplyr::where(is.character),
  ~ trimws(gsub(" {2,}", " ", gsub("[\r\n\t]+", " ", .x)))))
clip <- function(x, n = 220) ifelse(!is.na(x) & nchar(x) > n,
                                    paste0(substr(x, 1, n - 3), "..."), x)

# =============================================================================
# ASSERTION: the columns about to be lifted off `per genus` are constant within
# a comparison. If any varied, moving it to another sheet would destroy
# information rather than de-duplicate it.
# =============================================================================
LIFTED <- c("analysis_type", "phenotype_label", "tier", "external_dataset",
            "nhanes_model_role", "oral_site", "sequencing_type", "taxonomy_method",
            "nhanes_covariates", "external_covariates", "covariates_identical",
            "covariates_partially_comparable", "important_missing_external_covariates",
            "adjustment_comparability", "orientation_transform", "signs_comparable",
            "limitations")
varies <- P %>% group_by(nhanes_variable, label) %>%
  summarise(across(all_of(LIFTED), ~ dplyr::n_distinct(.x)), .groups = "drop") %>%
  select(all_of(LIFTED)) %>% summarise(across(everything(), max)) %>% unlist()
if (any(varies > 1))
  stop("ASSERTION FAILED: not constant within a comparison, cannot lift: ",
       paste(names(varies)[varies > 1], collapse = ", "))
asay("lift check: all ", length(LIFTED), " cohort-level columns constant per comparison")

# =============================================================================
# SHEET 1 — cohorts
# One row per data source, including the NHANES reference row.
# =============================================================================
# Abundance coverage is the number that matters and the one most often misread:
# Huang carries 16% of the NHANES genera by count but 88% of the community by
# abundance. Take the abundance table where a cohort has one, else its target list.
ov <- ovl %>%
  mutate(rank = ifelse(kind == "abundance", 1L, 2L)) %>%
  group_by(dataset_id) %>% slice_min(rank, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(dataset_id, ov_kind = kind, ov_n_genera = n_external_genera,
         ov_matched = n_match_normalized, ov_pct_common = pct_common_covered,
         ov_pct_abund = pct_abundance_covered)

# study_characteristics leaves specimen/assay blank for cohorts with no
# per-sample data (Ghensi is shotgun metagenomics, so it has no 16S region
# field), but master records both for every comparison. Coalesce so the sheet
# states what is known rather than what one source happened to carry.
cohort_attrs <- P %>%
  distinct(dataset_id = external_dataset, m_site = oral_site,
           m_seq = sequencing_type, m_tax = taxonomy_method)
stopifnot(!anyDuplicated(cohort_attrs$dataset_id))

cohorts <- sc %>%
  left_join(ov, by = "dataset_id") %>%
  left_join(cohort_attrs, by = "dataset_id") %>%
  transmute(
    `Cohort`                          = dataset_id,
    `Role`                            = ifelse(tier == "reference", "Reference (this study)",
                                               m12_cohort_class(tier)),
    `Citation`                        = citation,
    `DOI`                             = doi,
    `Participants (n)`                = n_participants,
    `Specimen`                        = clip(coalesce(oral_site, m_site)),
    `Assay / 16S region`              = clip(coalesce(region_16S, m_seq)),
    `Taxonomy`                        = clip(m_tax),
    `Phenotypes validated`            = variables_validated,
    `Covariates in model`             = clip(covariates_adjusted),
    `Adults only`                     = adults_only,
    `Subject de-duplicated`           = subject_dedup,
    `Genera in dataset (n)`           = coalesce(ov_n_genera, n_genera_dataset),
    `Genera matched to NHANES (n)`    = coalesce(ov_matched, n_overlap_nhanes_common),
    `NHANES common genera covered (%)`= round(ov_pct_common, 1),
    `NHANES abundance covered (%)`    = round(ov_pct_abund, 1),
    `Comparable genera after FDR`     = n_comparable_after_fdr,
    `Smoking definition`              = clip(smoking_definition),
    `Key caveat`                      = clip(key_caveat)) %>%
  arrange(match(Role, "Reference (this study)") %>% is.na(), Cohort)
asay("sheet 'cohorts': ", nrow(cohorts), " rows")

# =============================================================================
# SHEET 2 — summary
# The 8 genus-level comparisons, plus the coarser resolutions at which the
# validation was also run. `resolution` is the leading column so the two are
# never conflated.
# =============================================================================
cmp_meta <- P %>%
  distinct(nhanes_variable, label, nhanes_covariates, external_covariates,
           important_missing_external_covariates, covariates_identical,
           orientation_transform, signs_comparable, limitations)

genus_rows <- conc %>%
  left_join(bcor %>% select(nhanes_variable, label, pearson_r, spearman_rho),
            by = c("nhanes_variable", "label")) %>%
  left_join(enr %>% select(nhanes_variable, external_dataset, fisher_or, fisher_p),
            by = c("nhanes_variable", "external_dataset")) %>%
  left_join(cmp_meta, by = c("nhanes_variable", "label")) %>%
  transmute(
    `Resolution`                       = "genus",
    `Phenotype`                        = phenotype_label,
    `NHANES variable`                  = nhanes_variable,
    `WAS scheme`                       = analysis_type,
    `Cohort`                           = label,
    `Cohort class`                     = m12_cohort_class(tier),
    `NHANES associations (n taxa)`     = n_nhanes_primary,
    `Externally testable (n taxa)`     = n_testable,
    `Not testable (n taxa)`            = n_not_testable,
    `Same direction (n taxa)`          = n_same_direction,
    `Discordant (n taxa)`              = n_discordant,
    `Same direction (%)`               = round(pct_same_direction, 1),
    `95% CI lower`                     = round(ci_lo, 1),
    `95% CI upper`                     = round(ci_hi, 1),
    `Binomial p`                       = binom_p,
    `Spearman rho (effect sizes)`      = spearman_rho,
    `Pearson r (effect sizes)`         = pearson_r,
    `External nominal p<0.05 (n taxa)` = n_nominal_significant,
    `External FDR<0.05 (n taxa)`       = n_fdr_significant,
    `External FDR + same direction (n taxa)` = n_significant_same_dir,
    `Enrichment odds ratio`            = fisher_or,
    `Enrichment p`                     = fisher_p,
    `NHANES participants (n)`          = n_participants_nhanes,
    `External participants (n)`        = n_participants_external,
    `NHANES covariates`                = clip(nhanes_covariates),
    `External covariates`              = clip(external_covariates),
    `Missing external covariates`      = clip(important_missing_external_covariates),
    `Covariates identical`             = covariates_identical,
    `Adjustment comparability`         = adjustment_comparability,
    `Orientation transform`            = clip(orientation_transform),
    `Signs comparable`                 = signs_comparable,
    `Limitations`                      = clip(limitations, 300))

# Coarser resolutions. Both are published-target comparisons against SLL2, so
# there is no effect size on either side and the numeric columns stay empty.
coarse <- bind_rows(
  phy %>% group_by(Phenotype = phenotype_label, Cohort = cohort,
                   `NHANES variable` = nhanes_variable) %>%
    summarise(Resolution = "phylum",
              n_elig = dplyr::n(), n_test = sum(tested),
              n_same = sum(same_direction, na.rm = TRUE), .groups = "drop"),
  alp %>% group_by(Phenotype = phenotype_label, Cohort = cohort,
                   `NHANES variable` = nhanes_variable) %>%
    summarise(Resolution = "alpha diversity",
              n_elig = dplyr::n(), n_test = sum(tested),
              n_same = sum(same_direction, na.rm = TRUE), .groups = "drop")) %>%
  transmute(Resolution, Phenotype, `NHANES variable`, `WAS scheme` = NA_character_,
            Cohort, `Cohort class` = unname(COHORT_CLASS[["B"]]),
            `NHANES associations (n taxa)` = n_elig,
            `Externally testable (n taxa)` = n_test,
            `Not testable (n taxa)`        = n_elig - n_test,
            `Same direction (n taxa)`      = n_same,
            `Discordant (n taxa)`          = n_test - n_same,
            `Same direction (%)` = ifelse(n_test > 0, round(100 * n_same / n_test, 1), NA_real_),
            Limitations = paste0("Published ", Resolution,
                                 " rows only; no effect size on either side."))

summary_tbl <- bind_rows(genus_rows, coarse) %>%
  arrange(factor(Resolution, levels = c("genus", "phylum", "alpha diversity")),
          Phenotype, desc(`Externally testable (n taxa)`))
asay("sheet 'summary': ", nrow(summary_tbl), " rows (",
     sum(summary_tbl$Resolution == "genus"), " genus + ",
     sum(summary_tbl$Resolution != "genus"), " coarser)")

# =============================================================================
# SHEET 3 — per genus
# A column projection of master (primary fits), not a reshape. Long in cohort,
# so a genus tested by three cohorts appears three times; that keeps the columns
# identical across all six phenotypes instead of leaving empty cohort blocks.
# =============================================================================
per_genus <- P %>%
  transmute(
    `Phenotype`                = phenotype_label,
    `NHANES variable`          = nhanes_variable,
    `WAS scheme`               = analysis_type,
    `Genus`                    = genus,
    `Phylum`                   = phylum,
    `NHANES prevalence (%)`    = round(100 * nhanes_prevalence, 1),
    `NHANES beta (CLR)`        = nhanes_beta,
    `NHANES SE`                = nhanes_se,
    `NHANES p-value`           = nhanes_p,
    `NHANES FDR`               = nhanes_fdr,
    `NHANES q-value`           = nhanes_q,
    `NHANES direction`         = nhanes_direction,
    `NHANES participants (n)`  = nhanes_n,
    `Cohort`                   = label,
    `Cohort class`             = m12_cohort_class(tier),
    `External taxon as published` = ext_genus,
    `External beta`            = ext_beta,
    `External SE`              = ext_se,
    `External p-value`         = ext_p,
    `External FDR`             = ext_fdr,
    `External direction`       = ext_direction,
    `External participants (n)`= ext_n,
    `Tested`                   = tested,
    `Same direction`           = same_direction,
    `External nominal significant` = external_nominal_significant,
    `External FDR significant` = external_fdr_significant,
    `Validation class`         = validation_classification) %>%
  arrange(Phenotype, Cohort, desc(`Tested`), Phylum, Genus)
asay("sheet 'per genus': ", nrow(per_genus), " rows (",
     sum(per_genus$Tested), " tested / ", sum(!per_genus$Tested), " not tested)")

# =============================================================================
# INTEGRITY — the summary must be reproducible from the raw rows
# =============================================================================
chk <- per_genus %>% group_by(Phenotype, Cohort) %>%
  summarise(k_assoc = dplyr::n(), k_test = sum(Tested),
            k_nottest = sum(!Tested), k_same = sum(`Same direction`, na.rm = TRUE),
            k_disc = sum(!`Same direction`, na.rm = TRUE),
            k_nom  = sum(`External nominal significant`, na.rm = TRUE),
            k_fdr  = sum(`External FDR significant`, na.rm = TRUE),
            .groups = "drop")
cross <- summary_tbl %>% filter(Resolution == "genus") %>%
  select(Phenotype, Cohort, s_assoc = `NHANES associations (n taxa)`,
         s_test = `Externally testable (n taxa)`, s_nottest = `Not testable (n taxa)`,
         s_same = `Same direction (n taxa)`, s_disc = `Discordant (n taxa)`,
         s_nom = `External nominal p<0.05 (n taxa)`,
         s_fdr = `External FDR<0.05 (n taxa)`) %>%
  full_join(chk, by = c("Phenotype", "Cohort"))
if (anyNA(cross$s_assoc) || anyNA(cross$k_assoc))
  stop("ASSERTION FAILED: phenotype x cohort sets differ between 'summary' and 'per genus'")
bad <- cross %>% filter(s_assoc != k_assoc | s_test != k_test | s_nottest != k_nottest |
                        s_same != k_same | s_disc != k_disc | s_nom != k_nom | s_fdr != k_fdr)
if (nrow(bad)) {
  print(as.data.frame(bad)); stop("ASSERTION FAILED: 'summary' counts not reproducible from 'per genus'")
}
asay("integrity: all ", nrow(cross), " comparisons - summary counts reproduce from per-genus rows")

# The old S3/S4/S5 rows must all survive in the merged sheet.
old_keys <- bind_rows(lapply(
  c("table2_strongest_replicated.tsv", "table3_concordant_underpowered.tsv",
    "table4_discordant.tsv"),
  function(f) read_tsv(file.path(AUD, f), show_col_types = FALSE) %>%
    select(Phenotype, Genus, Cohort = `External cohort`))) %>% distinct()
new_keys <- per_genus %>% filter(Tested) %>% distinct(Phenotype, Genus, Cohort)
missing <- dplyr::anti_join(old_keys, new_keys, by = c("Phenotype", "Genus", "Cohort"))
if (nrow(missing)) {
  print(as.data.frame(head(missing, 20)))
  stop("ASSERTION FAILED: rows from tables 2-4 absent from 'per genus'")
}
asay("supersession: all ", nrow(old_keys), " rows of old tables 2-4 present among the ",
     nrow(new_keys), " tested rows")

# Excel hard limits, checked on exactly the text that will be written. Cell
# content is flattened first: a tab or newline surviving into the TSV handoff
# would split one logical row across several physical lines.
SHEETS <- list(cohorts = flat(cohorts), summary = flat(summary_tbl),
               `per genus` = flat(per_genus))
for (nm in names(SHEETS)) {
  ch <- unlist(SHEETS[[nm]][vapply(SHEETS[[nm]], is.character, logical(1))])
  if (length(ch) && any(nchar(ch) > 32767, na.rm = TRUE))
    stop("ASSERTION FAILED: cell over Excel's 32,767-character limit in '", nm, "'")
  if (length(ch) && any(grepl("[\r\n\t]", ch, useBytes = TRUE)))
    stop("ASSERTION FAILED: tab or newline inside a cell in '", nm, "'")
}
asay("cell limits: no tabs, newlines or over-length cells in any sheet")

# =============================================================================
# WRITE — R hands three flat TSVs to the packager; nothing but the .xlsx lands
# in the results tree.
# =============================================================================
tmp <- file.path(tempdir(), "m12_supp"); dir.create(tmp, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
paths <- c(cohorts = file.path(tmp, "cohorts.tsv"),
           summary = file.path(tmp, "summary.tsv"),
           `per genus` = file.path(tmp, "per_genus.tsv"))
for (nm in names(SHEETS)) write_tsv(SHEETS[[nm]], paths[[nm]], na = "")

PY <- Sys.getenv("M12_PYTHON", "/usr/bin/python3")
if (!nzchar(Sys.which(PY)) && !file.exists(PY)) stop("FATAL: python3 not found: ", PY)
rc <- system2(PY, c(file.path(HERE, "tsv_to_xlsx.py"), shQuote(OUT_XLSX),
                    shQuote(paste0("cohorts=",   paths[["cohorts"]])),
                    shQuote(paste0("summary=",   paths[["summary"]])),
                    shQuote(paste0("per genus=", paths[["per genus"]]))),
              stdout = TRUE, stderr = TRUE)
cat(paste0("  ", rc, collapse = "\n"), "\n")
if (!is.null(attr(rc, "status")) && attr(rc, "status") != 0)
  stop("FATAL: workbook packaging failed")
if (!file.exists(OUT_XLSX)) stop("FATAL: workbook not written: ", OUT_XLSX)

asay("wrote ", basename(OUT_XLSX), " (", round(file.size(OUT_XLSX) / 1024), " KB, 3 sheets)")
asay("")
asay("=== summary sheet ===")
print(as.data.frame(summary_tbl %>% select(Resolution, Phenotype, Cohort, `Cohort class`,
        `NHANES associations (n taxa)`, `Externally testable (n taxa)`,
        `Same direction (n taxa)`, `Same direction (%)`, `Binomial p`)))
