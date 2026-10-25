#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Validation audit: shared configuration
#
# Purpose : Single source of truth for the audit that asks what OralMicroNHANES
#           findings CAN be externally validated, beyond the age/smoking pair
#           that Phase 2 implemented.
#
#           Inherits every Phase 2 convention (thresholds, Arial 6 pt, the
#           canonical phylum palette, figure writers) by sourcing its config,
#           so the audit's figures are visually identical to Phase 2's.
#
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(purrr)
})
PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
V2 <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation/10_validation_analysis")
source(file.path(V2, "00_validation_config.R"))

AUD <- file.path(M12_OUT, "audit")
AUD_FIG_PDF <- file.path(AUD, "figures", "pdf")
AUD_FIG_PNG <- file.path(AUD, "figures", "png")
for (d in c(AUD, AUD_FIG_PDF, AUD_FIG_PNG, file.path(AUD, "logs")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
stopifnot(grepl("12_external_data_validation", AUD, fixed = TRUE))

# =============================================================================
# The five WAS schemes.
#
# ORIENTATION IS NOT UNIFORM AND THIS IS THE SINGLE MOST DANGEROUS FACT IN THE
# AUDIT (see section 17 of the work order):
#   * demoWAS / exWAS model  genus_CLR ~ phenotype
#       -> the GENUS is the dependent variable; a positive beta means the genus
#          increases with the phenotype.
#   * oradWAS / pheWAS / outWAS model  phenotype ~ genus_CLR
#       -> the GENUS is the INDEPENDENT variable; a positive beta means the
#          phenotype increases with the genus.
# Verified by inspection, not assumed: an assertion in step 1 fails the build
# if the genus is not where this table says it is.
#
# Both encode "genus and phenotype move together", so the SIGN is comparable
# across schemes, but the effect SCALE is not (CLR units per phenotype unit vs
# phenotype units per CLR unit). Magnitudes are therefore never pooled across
# orientations; only signs and ranks are.
# =============================================================================
WAS_SPEC <- tibble::tribble(
  ~analysis_type, ~rds_rel,                                                            ~genus_role,   ~note,
  "demoWAS", "results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds", "dependent",   "demographics; UNADJUSTED (n_covariates = 0)",
  "oradWAS", "results/2_oradWAS_out/result_clr/2_oradWAS_clr_tidied_complete.rds", "independent", "oral/dental health; genus is the PREDICTOR",
  "exWAS",   "results/3_exWAS_out/result_clr/3_exWAS_clr_tidied_complete.rds",     "dependent",   "exposures; adjusted",
  "pheWAS",  "results/4_pheWAS_out/result_clr/4_pheWAS_clr_tidied_complete.rds",   "independent", "clinical/anthropometric phenotypes; genus is the PREDICTOR",
  "outWAS",  "results/5_outWAS_out/result_clr/5_outWAS_clr_tidied_complete.rds",   "independent", "disease outcomes; genus is the PREDICTOR"
)

# Controlled vocabularies (work order sections 4, 7, 12)
COMPARABILITY_LEVELS <- c("DIRECTLY_COMPARABLE", "COMPARABLE_WITH_LIMITATIONS",
                          "NOT_COMPARABLE", "METADATA_UNAVAILABLE")
ADJUSTMENT_LEVELS    <- c("MATCHED", "PARTIALLY_MATCHED", "UNADJUSTED_EXTERNAL",
                          "NOT_RECONSTRUCTABLE")
VALIDATION_CLASSES   <- c("FDR_REPLICATED_SAME_DIRECTION",
                          "NOMINALLY_REPLICATED_SAME_DIRECTION",
                          "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT",
                          "DIRECTIONALLY_DISCORDANT",
                          "NOT_TESTABLE_PHENOTYPE", "NOT_TESTABLE_TAXON",
                          "NOT_TESTABLE_MAPPING", "MODEL_NOT_COMPARABLE")
RECOMMENDATIONS      <- c("RUN", "RUN_WITH_LIMITATIONS", "DO_NOT_RUN")

# Evidence-level colours: reuse the project's direction palette so "supported"
# is always the warm colour and "discordant" always the cool one, matching
# association_colors in every other module.
# validation_colors now lives in the shared Phase 2 config, because both
# phases draw this scale and they must stay colour-identical.

asay <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")

#' Save an audit figure using the Phase 2 writer (Arial, 6 pt, <= 7 in, PNG
#' rasterised from the PDF) but into the audit figure tree.
save_audit_figure <- function(name, plot_obj, width_mm, height_mm, dpi = 300) {
  old_pdf <- VAL_FIG_PDF; old_png <- VAL_FIG_PNG
  assign("VAL_FIG_PDF", AUD_FIG_PDF, envir = .GlobalEnv)
  assign("VAL_FIG_PNG", AUD_FIG_PNG, envir = .GlobalEnv)
  on.exit({assign("VAL_FIG_PDF", old_pdf, envir = .GlobalEnv)
           assign("VAL_FIG_PNG", old_png, envir = .GlobalEnv)}, add = TRUE)
  save_validation_figure(name, plot_obj, width_mm, height_mm, dpi)
}

#' Heat map with module-wide fixed cell geometry, written to the audit tree.
save_audit_heatmap <- function(name, p, n_cols, n_rows, ...) {
  old_pdf <- VAL_FIG_PDF; old_png <- VAL_FIG_PNG
  assign("VAL_FIG_PDF", AUD_FIG_PDF, envir = .GlobalEnv)
  assign("VAL_FIG_PNG", AUD_FIG_PNG, envir = .GlobalEnv)
  on.exit({assign("VAL_FIG_PDF", old_pdf, envir = .GlobalEnv)
           assign("VAL_FIG_PNG", old_png, envir = .GlobalEnv)}, add = TRUE)
  m12_save_heatmap(name, p, n_cols = n_cols, n_rows = n_rows, ...)
}

asay("audit config loaded (FDR<", FDR_THRESHOLD, ", prevalence>=", PREVALENCE_THRESHOLD, ")")
