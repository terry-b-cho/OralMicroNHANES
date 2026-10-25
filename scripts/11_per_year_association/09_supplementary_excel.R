#!/usr/bin/env Rscript

# =============================================================================
# MODULE 11 — STEP 09: publication supplementary workbook
# 09_supplementary_excel.R                          (run AFTER steps 05-06)
# =============================================================================
# Builds OUT_ROOT/Table_S11_summary_of_cross_cycle_validation_data_and_results.xlsx
# — the supplementary workbook restricted to the tables NECESSARY for the
# Figure S16 panels (A: cycle-pattern partitions; B/C: the row-level pooled +
# per-cycle estimates): cycle characteristics, validation summary, replication
# tiers, the row-level pooled-significant results, and the validated core.
# Accessory tables live as CSVs in figures/additional_accessory_figures/ and
# are deliberately NOT part of this workbook. The workbook itself NEVER uses
# internal repo terminology ("module 11" etc.).
#
# Mechanics: R stages each sheet as CSV + a manifest; the module-local Python
# venv (.pyxlsx_venv, openpyxl) assembles the styled workbook via
# 09_make_xlsx.py (R package installation is broken in this conda env's
# build configuration, so Excel writing is delegated to Python).
# =============================================================================

source("scripts/11_per_year_association/00_module11_config.R")  # run from the repo root

suppressPackageStartupMessages(library(jsonlite))

OUT_XLSX <- file.path(M11_OUT, "Table_S11_summary_of_cross_cycle_validation_data_and_results.xlsx")
STAGE <- file.path(tempdir(), "s11_stage")
dir.create(STAGE, recursive = TRUE, showWarnings = FALSE)

th <- M11_THRESH_PRIMARY
rd <- function(f) read_csv(file.path(M11_TABLES, f), show_col_types = FALSE)
t1  <- rd("table1_cycle_characteristics.csv")
t6  <- rd("table3_replication_tiers.csv")
t7  <- rd("table4_validated_core.csv")
t10 <- rd("table2_validation_summary.csv")
# design df joined into sheet 3 comes from the (accessory) calibration table
t9  <- read_csv(file.path(M11_FIGURES_ACCESSORY, "optional_table9_percycle_pvalue_calibration.csv"),
                show_col_types = FALSE)
rep_rowlvl <- readRDS(file.path(M11_CMP, "pooled_vs_F_vs_G", "pooled_replication_rowlevel.rds"))

# ---- Sheet 6: row-level pooled-significant results (the key data sheet) -----
rowlevel <- rep_rowlvl %>%
  filter(pooled_metric == "fdr") %>%
  mutate(validation_category = case_when(
    !testable_F | !testable_G ~ "not testable per cycle (variable absent)",
    !(stable_F & stable_G) ~ "per-cycle fit unstable (excluded)",
    repfdr_F & repfdr_G ~ "validated in both cycles",
    repfdr_F ~ "validated in cycle F only",
    repfdr_G ~ "validated in cycle G only",
    dir_concord_F & dir_concord_G ~ "same direction in both, not validated",
    TRUE ~ "direction inconsistent in >=1 cycle")) %>%
  transmute(analysis_type = analysis, variable, dependent_variable = dependent_var,
            independent_variable = exposure, model_term = term,
            genus = genus_label, phylum = phylum_std,
            beta_pooled = estimate_pooled, se_pooled = std.error_pooled,
            p_pooled = p.value_pooled, fdr_pooled = p.value.fdr_pooled,
            q_pooled = q.value_pooled, n_pooled = n_obs_pooled,
            beta_cycleF = estimate_F, se_cycleF = std.error_F,
            p_cycleF = p.value_F, replication_fdr_cycleF = repfdr_p_F,
            n_cycleF = n_obs_F,
            beta_cycleG = estimate_G, se_cycleG = std.error_G,
            p_cycleG = p.value_G, replication_fdr_cycleG = repfdr_p_G,
            n_cycleG = n_obs_G,
            genus_prevalence_pooled = prevalence_pooled,
            genus_prevalence_cycleF = prevalence_F,
            genus_prevalence_cycleG = prevalence_G,
            variable_availability = availability,
            validated_cycleF = repfdr_F, validated_cycleG = repfdr_G,
            direction_concordant_cycleF = dir_concord_F,
            direction_concordant_cycleG = dir_concord_G,
            validation_category)

# ---- Sheet 1: contents ------------------------------------------------------
contents <- tibble(Item = c(
  "Title",
  "Description",
  "Survey design (per cycle)",
  "Per-cycle p-value calibration",
  "Definition: validated",
  "Definition: replication FDR",
  "Significance thresholds",
  "Stability guard",
  "Variable availability",
  "Genus prevalence",
  "",
  "Sheet 3", "Sheet 4", "Sheet 5", "Sheet 6", "Sheet 7"
), Description = c(
  "Cross-cycle validation of the pooled NHANES 2009-2012 oral microbiome association results",
  "The survey-weighted CLR association analyses originally estimated on the pooled NHANES 2009-2012 sample (cycles F 2009-10 and G 2011-12) were re-estimated independently within each cycle using the identical regression framework (same schemas, covariates, reference levels, and models). This workbook contains the row-level results for every pooled FDR<0.05 association with its independent cycle-F and cycle-G estimates, the pooled-vs-per-cycle validation summaries, and the aggregate cross-cycle comparison tables underlying the accompanying figures.",
  "svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE) with unscaled 2-year MEC weights and survey.lonely.psu='certainty' (NCHS single-cycle design); the pooled analysis used 4-year reweighting (WTMEC2YR/2) with cycle-unique PSU/strata identifiers.",
  "A single cycle provides only 16 (F) / 17 (G) design degrees of freedom (PSUs minus strata); with 13 covariates the default svyglm residual df (~2) imposes a hard p-value floor. Per-cycle p-values therefore use t with df = the design df (Korn-Graubard / NCHS practice); BH-FDR and Storey q were recomputed within each analysis-by-cycle family. Standard errors are unchanged; default-reference quantities are retained in the released row-level data.",
  "A pooled FDR<0.05 association is VALIDATED in a cycle when its per-cycle replication FDR is <0.05 AND its per-cycle beta carries the pooled direction. Overall validation = validated in either cycle; the strongest tier = validated in both cycles.",
  "Benjamini-Hochberg FDR applied to the one-sided directional replication p-values (per-cycle p/2 when the beta carries the pooled direction, 1 - p/2 otherwise) WITHIN the pre-specified pooled FDR<0.05 family of each analysis type (replication-panel practice).",
  "Pooled significance: BH-FDR<0.05 (Storey q<0.05 reported in parallel). Per-cycle validation: replication FDR<0.05. Primary prevalence subset: pooled genus prevalence >10%.",
  "Quasi-separated per-cycle survey-logistic fits (|log-odds| or SE > 15) are numerically meaningless and are excluded from beta comparisons ('per-cycle fit unstable'), never counted as failed validation.",
  "Variables measured in only one cycle are classified F_ONLY/G_ONLY, analysed where measured, and marked 'not testable' in cross-cycle comparisons - never counted as failed validation. AVAILABLE_BUT_NOT_COMPARABLE marks variables whose exposure coding differed between cycles.",
  "Share of oral-rinse microbiome samples with relative abundance > 0, computed from the untransformed genus-level relative-abundance tables (pooled and per cycle).",
  "",
  "Cycle characteristics: analytic samples, design variables, design df, test counts per cycle",
  "Validation summary: per analysis type, counts and percentages of pooled FDR<0.05 associations validated in cycle F / G / both / either",
  "Replication tiers: direction concordance (with exact binomial 95% CI) and replication tiers of pooled-significant associations",
  "Pooled-significant results: row-level data for EVERY pooled FDR<0.05 association (pooled, cycle-F and cycle-G estimates, validation status)",
  "Validated core: associations validated in BOTH cycles"
))

# ---- Sheet 2: column dictionary --------------------------------------------
dict <- tribble(
  ~Sheet, ~Column, ~Description,
  "3-7", "analysis / analysis_type", "Association analysis type: 1_demoWAS (demographics), 2_oradWAS (oral disease), 3_exWAS (exposures), 4_pheWAS (clinical phenotypes), 5_outWAS (medical outcomes)",
  "3-7", "cycle", "NHANES cycle: F = 2009-2010, G = 2011-2012",
  "6,7", "variable", "The non-microbiome association variable (exposure, phenotype or outcome)",
  "6,7", "dependent_variable / independent_variable", "Model dependent and independent variables (the microbiome genus CLR abundance is the dependent variable in 1_demoWAS/3_exWAS and the independent variable otherwise)",
  "6,7", "model_term", "Regression term the estimate belongs to (equals the variable for continuous coding; a contrast label for categorical coding)",
  "6,7", "genus / phylum", "Genus label and phylum of the microbial feature (SILVA taxonomy)",
  "6,7", "beta_* / se_*", "Survey-weighted regression coefficient and its design-based standard error (pooled, cycle F, cycle G)",
  "6,7", "p_*", "Two-sided p-value (per-cycle values use the design-df t reference; see Contents)",
  "6,7", "fdr_pooled / q_pooled", "Benjamini-Hochberg FDR and Storey q of the pooled test within its analysis-wide family",
  "6,7", "replication_fdr_cycle*", "BH-FDR of the one-sided directional replication p-value within the pooled FDR<0.05 family",
  "6,7", "n_*", "Complete-case participant count entering the model",
  "6,7", "genus_prevalence_*", "Share of samples with relative abundance > 0 for the genus",
  "6", "validated_cycle*", "TRUE when replication FDR<0.05 AND the per-cycle beta carries the pooled direction",
  "6", "validation_category", "Mutually exclusive validation outcome (validated both / F only / G only / same direction / inconsistent / unstable / not testable)",
  "4", "n_validated_* / pct_validated_*", "Counts and percentages of pooled FDR<0.05 associations validated per cycle; denominators as named (of_testable = testable with stable fits in both cycles)",
  "5", "pct_dir_concordant_both, dir_concord_ci_lo/hi", "Percent of pooled-significant tests whose per-cycle betas both carry the pooled direction, with exact (Clopper-Pearson) binomial 95% CI",
  "5", "n_onesided_rep_* / n_repfdr_* / n_percycle_fdr_*", "Tier counts: one-sided directional p<0.05; replication FDR<0.05; per-cycle analysis-wide FDR<0.05 (each with pooled direction)")

# ---- Sheet 3: cycle characteristics (+ design df) ---------------------------
t3 <- t1 %>% left_join(
  t9 %>% distinct(cycle, design_df), by = "cycle") %>%
  relocate(design_df, .after = nhanes_years)

# ---- Stage all sheets -------------------------------------------------------
sheets <- list(
  list(sheet = "1. Contents",                    df = contents, long_text = TRUE),
  list(sheet = "2. Column dictionary",           df = dict,     long_text = TRUE),
  list(sheet = "3. Cycle characteristics",       df = t3),
  list(sheet = "4. Validation summary",          df = t10),
  list(sheet = "5. Replication tiers",           df = t6),
  list(sheet = "6. Pooled-significant results",  df = rowlevel),
  list(sheet = "7. Validated core",              df = t7))

manifest <- lapply(seq_along(sheets), function(i) {
  s <- sheets[[i]]
  csv_path <- file.path(STAGE, sprintf("sheet_%02d.csv", i))
  write_csv(s$df, csv_path, na = "")
  list(sheet = s$sheet, csv = csv_path, long_text = isTRUE(s$long_text))
})
manifest_path <- file.path(STAGE, "manifest.json")
write_json(manifest, manifest_path, auto_unbox = TRUE)

# ---- Assemble via the module-local Python venv ------------------------------
py <- file.path(M11_BASE, "scripts/11_per_year_association/.pyxlsx_venv/bin/python")
mk <- file.path(M11_BASE, "scripts/11_per_year_association/09_make_xlsx.py")
st <- system2(py, c(shQuote(mk), shQuote(manifest_path), shQuote(OUT_XLSX)),
              stdout = TRUE, stderr = TRUE)
cat(st, sep = "\n")
stopifnot(file.exists(OUT_XLSX), file.size(OUT_XLSX) > 1e5)

# Guard: internal repo terminology must not appear anywhere in the workbook
for (s in sheets) {
  txt <- paste(capture.output(write.csv(s$df, row.names = FALSE)), collapse = " ")
  if (grepl("module\\s*11|module11", txt, ignore.case = TRUE))
    stop("Internal terminology leaked into sheet: ", s$sheet)
}
cat(sprintf("STEP 09 complete: %s (%.1f MB, %d sheets, %s rows in the row-level sheet)\n",
            OUT_XLSX, file.size(OUT_XLSX) / 1e6, length(sheets),
            format(nrow(rowlevel), big.mark = ",")))