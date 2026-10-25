#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 2 — External Validation Analysis
# Shared configuration: thresholds, palettes, theme, fonts, figure writers
#
# Purpose : Single source of truth for the validation analysis. Sourced by every
#           0N_*.R script. Defines the comparison thresholds, the project's
#           canonical phylum palette, the compact 6pt Arial theme, and the PDF
#           writer + PDF->PNG converter used for the report.
#
# Usage   : source(".../10_validation_analysis/00_validation_config.R")
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(ggplot2); library(grid)
})

# === USER CONFIG =============================================================
PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
# =============================================================================
set.seed(42)

M12      <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
M12_OUT  <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation")
VAL      <- file.path(M12_OUT, "validation")
VAL_FIG_PDF <- file.path(VAL, "figures", "pdf")
VAL_FIG_PNG <- file.path(VAL, "figures", "png")
for (d in c(VAL, VAL_FIG_PDF, VAL_FIG_PNG, file.path(VAL, "logs")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Output guard: never write outside the validation tree (module 5.5 convention).
stopifnot(grepl("12_external_data_validation", VAL, fixed = TRUE),
          grepl("validation", VAL, fixed = TRUE))

source(file.path(M12, "00_shared", "harmonize_helpers.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a
vsay <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")

# --- Comparison thresholds ---------------------------------------------------
# Set once, reported in the methods section of the report.
FDR_THRESHOLD        <- 0.05   # NHANES significance (BH, scheme-wide family)
PREVALENCE_THRESHOLD <- 0.10   # NHANES prevalence floor -> 91 genera
MIN_ADULT_AGE        <- 18     # "adults only"
EXTERNAL_FDR         <- 0.05   # BH within external dataset x variable

# --- NHANES source specification --------------------------------------------
# `term` holds the literal model-matrix name; there is NO "indep_var" value.
# Main effects are the rows where `fdr_corrected == TRUE` (a "was-corrected"
# flag, NOT a significance flag - covariate rows carry an uncorrected
# p.value.fdr equal to the raw p).
#
# `genus_role` is NOT cosmetic. demoWAS and exWAS model genus_CLR ~ phenotype,
# so the genus is the DEPENDENT variable and the effect is "CLR per unit
# exposure". oradWAS models phenotype ~ genus_CLR, so the genus is the
# INDEPENDENT variable and the effect is "phenotype units per CLR unit". Both
# encode co-movement, so DIRECTIONS are comparable across schemes; MAGNITUDES
# are not, and no figure places them on a shared axis.
NHANES_SPEC <- tibble::tribble(
  ~variable,            ~role,      ~genus_role,   ~rds_rel,                                                    ~independent_var,          ~term,                        ~expected_n_main, ~expected_fdr_sig, ~n_covariates_expected,
  "age",                "primary",  "dependent",   "results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds", "RIDAGEYR",               "RIDAGEYR",                   1349L,            1132L,             0L,
  "smoking",            "primary",  "dependent",   "results/3_exWAS_out/result_clr/3_exWAS_clr_tidied_complete.rds",     "SMQ_current_ever_never", "SMQ_current_ever_never2",    1349L,            79L,               13L,
  "gum_disease",        "primary",  "independent", "results/2_oradWAS_out/result_clr/2_oradWAS_clr_tidied_complete.rds", "GUM_DISEASE_OHAROCGP",   "GUM_DISEASE_OHAROCGP",       1349L,            1067L,             13L,
  # Secondary: extracted so notebooks/03_internal_robustness.qmd can run its
  # measurement-robustness check. It is NOT part of the external validation and
  # appears in no pipeline figure, table or report.
  "smoking_cotinine",   "secondary","dependent",   "results/3_exWAS_out/result_clr/3_exWAS_clr_tidied_complete.rds",     "LBXCOT",                 "LBXCOT",                     1349L,            72L,               13L
)

# The variables carried through comparison, figures, tables and the report.
# Declared once so a fourth phenotype cannot be added to some steps and not
# others (it previously appeared hardcoded in four separate files).
PRIMARY_VARIABLES <- NHANES_SPEC$variable[NHANES_SPEC$role == "primary"]

# Display labels, used by every figure and table.
VAR_LABELS <- c(age = "Age (per year)",
                smoking = "Smoking (current vs never)",
                gum_disease = "Gum disease")
VAR_LABELS_SHORT <- c(age = "Age", smoking = "Smoking", gum_disease = "Gum disease")
stopifnot(all(PRIMARY_VARIABLES %in% names(VAR_LABELS)))

# Direction semantics, stated once and reused in every caption.
DIRECTION_NOTE <- paste(
  "The taxon is the DEPENDENT variable, so a positive estimate means the genus CLR abundance",
  "INCREASES with the exposure (per year of age; in current vs never smokers).")

# `SMQ_current_ever_never2` = CURRENT vs never. Modules 6 and 9 label these
# levels backwards; the SQLite SMQ_F/SMQ_G tables and the exWAS schema are
# authoritative. 01_extract asserts the level counts rather than trusting either.
SMQ_LEVELS <- c("0" = "never", "1" = "former", "2" = "current")

# --- Cohort class registry ---------------------------------------------------
# Two distinctions govern every validation method downstream, so the cohort
# classes are named for what they let us do rather than ranked A/B:
#
# REFIT COHORT ("refit"): supplies a genus-by-sample abundance matrix AND
#   linkable per-sample phenotypes, so the model is fitted here. A disagreeing
#   direction is therefore genuinely DISCORDANT. Limitation: sampling-method
#   differences between the cohort and NHANES remain.
#
# PUBLISHED-TARGET COHORT ("published-target"): supplies only the taxa its
#   authors reported significant, with the published direction and magnitude
#   relative to the host variable of interest. A NHANES finding absent from
#   that list is NOT TESTED, never discordant, because no complementary null
#   set was published.
#
# The code-level key stays "A"/"B" because it is a join key in every result table
# and in config/<id>.yaml; COHORT_CLASS is the only spelling shown to a reader.
COHORT_CLASS       <- c(A = "refit cohort", B = "published-target cohort")
COHORT_CLASS_SHORT <- c(A = "refit", B = "published-target")
m12_cohort_class <- function(tier, short = FALSE)
  unname((if (short) COHORT_CLASS_SHORT else COHORT_CLASS)[as.character(tier)])

# Dataset labels carry the PUBLICATION (first author + year) so a reader of any
# figure knows the source without consulting a key. Defined once here and
# referenced by name wherever a label is written or matched, so a rename cannot
# leave a stale hardcoded string behind.
COHORT_LABELS <- c(
  huang_unadj = "Huang 2020 (adults, subject-dedup)",
  huang_adj   = "Huang 2020 (adj. sex)",
  sll2_age    = "Willis 2022 (SLL2; published, Table 1)",
  sll2_smoke  = "Willis 2022 (SLL2; published, Table 2)",
  twins_age   = "Wells 2022 (TwinsUK; published, Fig 5)",
  ghensi_gum  = "Ghensi 2020 (published, peri-implantitis)")

# Short forms for narrow axes. Still publication-anchored: the cohort acronym
# alone ("SLL2") does not tell a reader whose data it is.
COHORT_SHORT <- c(HuangS_2020_AgePrediction   = "Huang 2020",
                  WillisJR_2022_SacaLaLengua2 = "Willis 2022 (SLL2)",
                  WellsPM_2022_TwinsUKSaliva  = "Wells 2022 (TwinsUK)",
                  GhensiP_2020_ImplantDisease = "Ghensi 2020")

# --- Supplementary Figure S15 panel registry ---------------------------------
# The nine Phase 2 figures are NOT standalone figures; they are sub-panels of a
# single manuscript figure, Supplementary Figure S15. File stems, captions and
# the order of the report's Results sections all derive from this one table, so
# "Panel D-2" names exactly one file, one caption and one section.
#
# The stems sort alphabetically into panel order, so `ls figures/pdf` reads in
# the order the panels are assembled.
S15_FIG <- c(
  `A-1` = "S15_A1_sample_sizes",
  `A-2` = "S15_A2_concordance_summary",
  `A-3` = "S15_A3_external_coverage",
  `B`   = "S15_B_effect_scatter",
  `C`   = "S15_C_concordance_strata",
  `D-1` = "S15_D1_concordance_grid_age",
  `D-2` = "S15_D2_concordance_grid_gum_disease",
  `D-3` = "S15_D3_concordance_grid_smoking",
  `E`   = "S15_E_effect_vs_prevalence",
  # Panel D version 2: the same grids graded by strength of evidence rather
  # than by direction alone. Additional to D-1..D-3, which are unchanged.
  # Version 2 keeps version 1's D-1 / D-2 / D-3 sub-panel numbering - same
  # phenotype in the same slot - and marks the version in the stem, so the two
  # versions of a phenotype sort next to each other.
  `D-1 v2` = "S15_D1_v2_classified_grid_age",
  `D-2 v2` = "S15_D2_v2_classified_grid_gum_disease",
  `D-3 v2` = "S15_D3_v2_classified_grid_smoking")

# The D panels run age, gum disease, smoking - NOT the PRIMARY_VARIABLES order
# (age, smoking, gum_disease) that the generating loop iterates. The mapping is
# therefore explicit; deriving it from the loop index would silently mislabel
# two of the three grids.
S15_GRID_PANEL    <- c(age = "D-1",    gum_disease = "D-2",    smoking = "D-3")
S15_GRID_PANEL_V2 <- c(age = "D-1 v2", gum_disease = "D-2 v2", smoking = "D-3 v2")

#' Prefix a caption with its S15 panel, so no caption can drift from its file.
m12_s15_cap <- function(panel, text) {
  stopifnot(panel %in% names(S15_FIG))
  sprintf("Figure S15, Panel %s. %s", panel, text)
}

# NOTE: no refit cohort carries a usable smoking variable. Smoking therefore
# has published-target evidence only, and no external effect-size correlation
# is possible for it. This is a stated limitation of the analysis, asserted in
# 03 so it cannot be forgotten, not an oversight.
#
# TERMINOLOGY: "retained" means a dataset CONTRIBUTED EVIDENCE to the
# validation. A dataset that was acquired and reconstructed but yielded no
# comparable result is NOT retained; it is excluded, with the reason recorded.
# The per-dataset disposition lives in config/<id>.yaml.
REFIT_COHORTS <- tibble::tribble(
  ~dataset_id,                      ~variables,          ~adults_only, ~subject_dedup, ~note,
  "HuangS_2020_AgePrediction",      "age",               TRUE,         TRUE,           "2,550 samples / 731 subjects; tongue-swab dominated"
)
PUBLISHED_TARGET_COHORTS <- c("WillisJR_2022_SacaLaLengua2",
                              "WellsPM_2022_TwinsUKSaliva", "GhensiP_2020_ImplantDisease")
# Back-compatible aliases: `tier` remains the join key in every result table.
TIER_A <- REFIT_COHORTS
TIER_B <- PUBLISHED_TARGET_COHORTS

RETAINED_IDS <- c(REFIT_COHORTS$dataset_id, PUBLISHED_TARGET_COHORTS)
stopifnot(setequal(names(COHORT_SHORT), RETAINED_IDS))
# Datasets that were assessed and excluded are no longer carried in the module;
# the evidence that disqualified each is in notebooks/02_qc_and_exclusions.html.
# Every cohort present here therefore contributed evidence: assert that.
stopifnot(setequal(RETAINED_IDS, list.files(file.path(M12, "config"),
                                            pattern = "[.]yaml$") |>
                     sub(pattern = "[.]yaml$", replacement = "")))

# Adding a phenotype must not leave its grid without a panel id.
stopifnot(setequal(names(S15_GRID_PANEL), PRIMARY_VARIABLES),
          setequal(names(S15_GRID_PANEL_V2), PRIMARY_VARIABLES),
          all(c(S15_GRID_PANEL, S15_GRID_PANEL_V2) %in% names(S15_FIG)))

OUTCOME_LEVELS <- c("concordant", "discordant", "not tested")

# --- Fonts -------------------------------------------------------------------
# The user requires Arial. `ArialMT` is a built-in pdfFonts() entry and works
# for the PDF device. We assert it rather than let R silently substitute.
#
# NOTE: `systemfonts` resolves "Arial" to DejaVuSans on this system, so a
# raster device would NOT produce Arial. Report PNGs are therefore rendered
# from the PDFs with pdftoppm (see save_validation_png), never re-drawn.
PDF_FAMILY      <- "ArialMT"
THEME_BASE_SIZE <- 6

m12_assert_arial <- function() {
  if (!PDF_FAMILY %in% names(grDevices::pdfFonts()))
    stop("FONT ASSERTION FAILED: '", PDF_FAMILY, "' is not in pdfFonts(). ",
         "Refusing to emit figures in a substituted font.")
  tf <- tempfile(fileext = ".pdf")
  # Draw with grid, not base graphics: base graphics enforce default margins
  # and error with "figure margins too large" on a small canvas, which would
  # look like a font failure and mask the thing we are actually testing.
  opened <- FALSE
  err <- NULL
  tryCatch({
    grDevices::pdf(tf, width = 2, height = 2, family = PDF_FAMILY)
    opened <- TRUE
    grid::grid.text("Ag", gp = grid::gpar(fontfamily = PDF_FAMILY, fontsize = 6))
  }, error = function(e) err <<- conditionMessage(e))
  if (opened) try(grDevices::dev.off(), silent = TRUE)   # never leak the device
  if (!is.null(err))
    stop("FONT ASSERTION FAILED: pdf(family='", PDF_FAMILY, "') errored: ", err)
  if (!file.exists(tf) || file.size(tf) < 500)
    stop("FONT ASSERTION FAILED: pdf(family='", PDF_FAMILY, "') produced no output.")
  unlink(tf)
  invisible(TRUE)
}

# --- Phylum palette ----------------------------------------------------------
# Canonical 17-entry form, copied VERBATIM from
# scripts/8_network_analyses/2_network.R:112 (identical in modules 4, 5, 9 apart
# from their trailing catch-all entries). Do not re-order or re-code: figures in
# this module must be colour-comparable with modules 4/5/7/8/9.
phylum_colors <- c(
  "Firmicutes"             = "#F38400",  # kelly_3   - orange
  "Bacteroidetes"          = "#0067A5",  # kelly_10  - blue
  "Actinobacteria"         = "#8DB600",  # kelly_17  - yellow green
  "Proteobacteria"         = "#E68FAC",  # kelly_9   - pink
  "Fusobacteria"           = "#BE0032",  # kelly_5   - red
  "Spirochaetae"           = "#F3C300",  # kelly_1   - yellow
  "Cyanobacteria"          = "#875692",  # kelly_2   - purple
  "Acidobacteria"          = "#F6A600",  # kelly_13  - orange yellow
  "Candidate division SR1" = "#2B3D26",  # kelly_20  - forest green
  "Planctomycetes"         = "#332288",  # safe_violet
  "Saccharibacteria"       = "#B3446C",  # kelly_14  - purple red
  "Synergistetes"          = "#A1CAF1",  # kelly_4   - light blue
  "Tenericutes"            = "#654522",  # kelly_18  - dark brown
  "Verrucomicrobia"        = "#C2B280",  # kelly_6   - tan
  "unclassified"           = "#DDDDDD",  # pale_grey
  "NA"                     = "#DDDDDD",  # very light grey
  "TM7"                    = "#000000"   # Contrast Black
)

# Module 8's coercion, verbatim in behaviour: anything not in the palette
# becomes "unclassified" rather than silently dropping out of the scale.
m12_phylum_std <- function(v) {
  s <- trimws(as.character(v))
  s[is.na(s)] <- "NA"
  s[s == ""]  <- "unclassified"
  s[!(s %in% names(phylum_colors))] <- "unclassified"
  s
}

# Effect-direction colours (module 4, verbatim)
association_colors <- c("Negative" = "#0072B2", "Positive" = "#CC6677")

# Three-state outcome colours: concordant / discordant reuse the direction
# palette; "not tested" is deliberately the palette's neutral grey so it never
# reads as a result.
outcome_colors <- c("concordant" = "#CC6677", "discordant" = "#0072B2",
                    "not tested" = "#DDDDDD")

# --- Validation classification, graded by strength of evidence ---------------
# The three-state outcome (concordant / discordant / not tested) says WHETHER a
# direction agreed; this grades HOW WELL it agreed. Defined here rather than in
# the audit config because both phases draw it: Phase 2 for Panel D version 2,
# Phase 3 for its evidence matrices, and the two must be colour-identical.
#
# A published-target cohort publishes no p-value, so a concordant tile from one
# can only ever be "same direction, not significant" - that is a property of the
# source, not a weak result.
validation_colors <- c(
  "FDR_REPLICATED_SAME_DIRECTION"            = "#8B2F42",
  "NOMINALLY_REPLICATED_SAME_DIRECTION"      = "#CC6677",
  "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT" = "#E8B4BC",
  "DIRECTIONALLY_DISCORDANT"                 = "#0072B2",
  "NOT_TESTABLE_PHENOTYPE"                   = "#EDEDED",
  "NOT_TESTABLE_TAXON"                       = "#DDDDDD",
  "NOT_TESTABLE_MAPPING"                     = "#CCCCCC",
  "MODEL_NOT_COMPARABLE"                     = "#BBBBBB")

# Compact forms for a figure legend, where the full sentences below do not fit
# a 50 mm-wide panel. Same keys, same order.
VALIDATION_CLASS_LABELS_SHORT <- c(
  FDR_REPLICATED_SAME_DIRECTION            = "FDR replicated",
  NOMINALLY_REPLICATED_SAME_DIRECTION      = "nominal p<0.05",
  DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT = "same direction, no p",
  DIRECTIONALLY_DISCORDANT                 = "opposite direction",
  NOT_TESTABLE_TAXON                       = "not testable",
  NOT_TESTABLE_PHENOTYPE                   = "not testable (phenotype)",
  NOT_TESTABLE_MAPPING                     = "not testable (mapping)",
  MODEL_NOT_COMPARABLE                     = "model not comparable")

# The five grades that can arise in Phase 2. The other three describe failures
# only the audit's wider phenotype set can produce, and listing them here would
# put three permanently empty rows in every legend.
PHASE2_CLASSES <- c("FDR_REPLICATED_SAME_DIRECTION",
                    "NOMINALLY_REPLICATED_SAME_DIRECTION",
                    "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT",
                    "DIRECTIONALLY_DISCORDANT",
                    "NOT_TESTABLE_TAXON")

VALIDATION_CLASS_LABELS <- c(
  FDR_REPLICATED_SAME_DIRECTION            = "FDR replicated, same direction",
  NOMINALLY_REPLICATED_SAME_DIRECTION      = "nominal p<0.05, same direction",
  DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT = "same direction, no p-value",
  DIRECTIONALLY_DISCORDANT                 = "opposite direction",
  NOT_TESTABLE_TAXON                       = "not testable (taxon absent)",
  NOT_TESTABLE_PHENOTYPE                   = "not testable (phenotype)",
  NOT_TESTABLE_MAPPING                     = "not testable (mapping)",
  MODEL_NOT_COMPARABLE                     = "model not comparable")
stopifnot(setequal(names(validation_colors), names(VALIDATION_CLASS_LABELS)),
          setequal(names(validation_colors), names(VALIDATION_CLASS_LABELS_SHORT)),
          all(PHASE2_CLASSES %in% names(validation_colors)))

#' Grade a Phase 2 comparison row by strength of external evidence.
#'
#' Reproduces the audit's `validation_classification` from the columns Phase 2
#' already holds, so Panel D version 2 needs no Phase 3 output. Verified to
#' match master_validation_table.tsv exactly on the three Phase 2 phenotypes.
m12_classify_evidence <- function(tested, outcome, ext_p, ext_fdr) {
  dplyr::case_when(
    !tested                                        ~ "NOT_TESTABLE_TAXON",
    outcome == "discordant"                        ~ "DIRECTIONALLY_DISCORDANT",
    !is.na(ext_fdr) & ext_fdr < EXTERNAL_FDR       ~ "FDR_REPLICATED_SAME_DIRECTION",
    !is.na(ext_p)   & ext_p   < 0.05               ~ "NOMINALLY_REPLICATED_SAME_DIRECTION",
    TRUE                                           ~ "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT")
}

# FDR -> dot size, module 4's calc_dot_size / sig_bin_from_dot_size, verbatim
calc_dot_size <- function(p) {
  dplyr::case_when(
    p < 1e-6 ~ 3.5, p < 1e-5 ~ 3.0, p < 1e-4 ~ 2.5,
    p < 1e-3 ~ 2.0, p < 1e-2 ~ 1.5, TRUE ~ 1.0)
}
SIG_LEVELS <- c("FDR < 0.05", "FDR < 1e-2", "FDR < 1e-3",
                "FDR < 1e-4", "FDR < 1e-5", "FDR < 1e-6")
sig_bin_from_dot_size <- function(d) {
  dplyr::case_when(
    d >= 3.5 ~ "FDR < 1e-6", d >= 3.0 ~ "FDR < 1e-5", d >= 2.5 ~ "FDR < 1e-4",
    d >= 2.0 ~ "FDR < 1e-3", d >= 1.5 ~ "FDR < 1e-2", TRUE ~ "FDR < 0.05")
}
SIG_SIZE_VALUES <- setNames(c(0.8, 1.1, 1.4, 1.7, 2.0, 2.4), SIG_LEVELS)

# --- Theme -------------------------------------------------------------------
# Module 5.5's theme_compact(), built on egg::theme_article so the look matches
# the rest of the repo. Every figure in this module uses it - there is no base
# graphics and no other theme.
#
# EVERY text element is EXACTLY `base` (6 pt). The earlier form shrank axis
# text, captions and legend text to 4.2-5 pt and inflated titles to 7 pt, which
# violated the uniform-6pt requirement. Nothing here may deviate from `base`:
# m12_assert_font_sizes() fails the build if it does.
theme_compact <- function(base = THEME_BASE_SIZE, family = PDF_FAMILY) {
  txt <- function(...) ggplot2::element_text(size = base, family = family, ...)
  egg::theme_article(base_size = base, base_family = family) +
    ggplot2::theme(
      text             = txt(),
      plot.title       = txt(face = "bold"),
      plot.subtitle    = txt(color = "grey30"),
      plot.caption     = txt(color = "grey45", hjust = 0),
      axis.title       = txt(),
      axis.text        = txt(),
      legend.title     = txt(),
      legend.text      = txt(),
      strip.text       = txt(face = "bold"),
      # ggplot defaults plot.tag to rel(1.2) = 7.2 pt. Pin it: patchwork can
      # draw tags, and the assertion rightly refuses any non-6pt element.
      plot.tag         = txt(face = "bold"),
      # Align titles to the whole canvas, not to the panel. With ggplot's
      # default ("panel") a title starts after the y-axis label gutter, so on a
      # narrow heat map it has only the tile block to live in and is clipped.
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      # Keys and margins are the only levers left for compactness once the type
      # size is fixed, so they are held tight.
      legend.key.size  = grid::unit(0.32, "lines"),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt"),
      legend.spacing   = grid::unit(2, "pt"),
      panel.grid.minor = ggplot2::element_blank(),
      # 1 pt at the top clipped the ascenders of the 6 pt bold title.
      plot.margin      = ggplot2::margin(3, 2, 1, 1, unit = "pt")
    )
}

# geom_text/geom_label/ggrepel take `size` in MILLIMETRES, not points. Passing
# the point size directly (size = 6) would draw ~17 pt text. This constant is
# the only correct way to write 6 pt inside a layer.
GEOM_TEXT_SIZE <- THEME_BASE_SIZE / ggplot2::.pt

#' Fail the build if any text in a figure is not exactly THEME_BASE_SIZE.
#'
#' Walks a ggplot or patchwork object: theme element_text sizes, and the `size`
#' aesthetic of every text-drawing layer (converted mm -> pt). Guards the
#' uniform-6pt requirement against future edits that reintroduce a 4.5 pt axis.
m12_assert_font_sizes <- function(p, name = "figure", base = THEME_BASE_SIZE) {
  bad <- character(0)
  walk_plot <- function(g) {
    if (inherits(g, "patchwork")) {
      for (q in g$patches$plots) walk_plot(q)
      g$patches <- NULL
    }
    if (!inherits(g, "gg")) return(invisible(NULL))
    th <- g$theme
    for (nm in names(th)) {
      el <- th[[nm]]
      if (inherits(el, "element_text") && !is.null(el$size) && is.numeric(el$size)) {
        # element_text(size = rel(x)) multiplies the inherited size; resolve it
        # so a rel() cannot smuggle in a non-6pt element.
        pt <- if (inherits(el$size, "rel")) base * as.numeric(el$size) else as.numeric(el$size)
        if (!isTRUE(all.equal(pt, base)))
          bad <<- c(bad, sprintf("theme$%s = %.2f pt%s", nm, pt,
                                 if (inherits(el$size, "rel")) " (via rel)" else ""))
      }
    }
    for (ly in g$layers) {
      gm <- class(ly$geom)[1]
      if (grepl("GeomText|GeomLabel", gm)) {
        sz <- ly$aes_params$size
        if (!is.null(sz) && is.numeric(sz)) {
          pt <- sz * ggplot2::.pt
          if (!isTRUE(all.equal(pt, base, tolerance = 0.02)))
            bad <<- c(bad, sprintf("%s size = %.2f mm (%.2f pt)", gm, sz, pt))
        }
      }
    }
    invisible(NULL)
  }
  walk_plot(p)
  if (length(bad))
    stop("FONT SIZE ASSERTION FAILED for ", name, ": every text element must be ",
         base, " pt. Offenders: ", paste(unique(bad), collapse = "; "))
  invisible(TRUE)
}

# --- Figure writers ----------------------------------------------------------
MAX_FIG_INCHES <- 7   # user requirement: no plot exceeds 7 in on either axis

#' Write a ggplot to PDF in Arial at an asserted size, then render a matching
#' PNG from that PDF. Sizes are given in mm (project convention).
save_validation_figure <- function(name, plot_obj, width_mm, height_mm, dpi = 300) {
  m12_assert_arial()
  # A gtable has no theme to walk; m12_save_heatmap asserts on the ggplot first
  # and only then converts, so nothing reaches the device unchecked.
  if (inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork"))
    m12_assert_font_sizes(plot_obj, name)
  else if (!inherits(plot_obj, "gtable"))
    stop("save_validation_figure: unsupported object for ", name)
  w_in <- width_mm / 25.4; h_in <- height_mm / 25.4
  if (w_in > MAX_FIG_INCHES || h_in > MAX_FIG_INCHES)
    stop(sprintf("FIGURE TOO LARGE: %s is %.2f x %.2f in (max %d). Reduce the size.",
                 name, w_in, h_in, MAX_FIG_INCHES))

  pdf_path <- file.path(VAL_FIG_PDF, paste0(name, ".pdf"))
  grDevices::pdf(pdf_path, width = w_in, height = h_in, family = PDF_FAMILY,
                 colormodel = "rgb", pointsize = THEME_BASE_SIZE, useDingbats = FALSE)
  if (inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) print(plot_obj)
  else grid::grid.draw(plot_obj)
  grDevices::dev.off()

  png_path <- save_validation_png(pdf_path, dpi = dpi)
  vsay("  figure: ", basename(pdf_path), sprintf("  (%.2f x %.2f in)", w_in, h_in))
  invisible(list(pdf = pdf_path, png = png_path))
}

# --- Title wrapping ----------------------------------------------------------
# ggplot does not wrap titles; it clips them. Now that captions carry the unit
# and both denominators they are long, and the figures are deliberately narrow,
# so every title must be wrapped to the canvas it will actually be drawn on.
#
# Mean Arial glyph is ~0.62 em, and 6 pt = 2.117 mm, so one character occupies
# ~1.31 mm. Bold titles run slightly wider, hence the 0.94 safety factor.
MM_PER_CHAR_6PT <- 0.62 * THEME_BASE_SIZE * (25.4 / 72)

m12_wrap <- function(x, width_mm, bold = FALSE) {
  if (is.null(x) || !nzchar(x)) return(x)
  n <- max(12L, floor(width_mm * (if (bold) 0.94 else 1) / MM_PER_CHAR_6PT))
  paste(strwrap(x, width = n), collapse = "\n")
}

#' Wrap a plot's own title and subtitle to a known canvas width.
m12_wrap_labs <- function(p, width_mm) {
  if (!inherits(p, "gg")) return(p)
  if (!is.null(p$labels$title))
    p$labels$title <- m12_wrap(p$labels$title, width_mm, bold = TRUE)
  if (!is.null(p$labels$subtitle))
    p$labels$subtitle <- m12_wrap(p$labels$subtitle, width_mm)
  p
}

# --- Heat maps: one cell is the same size in every figure ---------------------
# A tile in the smoking grid must be exactly as wide as a tile in the 3-column
# age grid, or a reader comparing panels silently mis-reads the scale. ggplot
# sizes panels to fill whatever canvas it is given, which is precisely what
# breaks that, so the panel is pinned with egg::set_panel_size() and the DEVICE
# is then sized to the resulting gtable instead of the other way round.
#
# Cell width is set by the widest thing that has to fit inside one column: the
# rotated cohort label. Cell height is bounded below by the type size - 6 pt is
# 2.12 mm, so 2.35 mm is the tightest row pitch that does not collide.
HEATMAP_CELL_W_MM <- 6.5
HEATMAP_CELL_H_MM <- 2.35

# Rows that will not fit on one page at the fixed pitch are wrapped into
# side-by-side blocks, like a phone book, rather than being squashed. Cell size
# is the invariant; the page adapts around it.
HEATMAP_MAX_ROWS_PER_BLOCK <- 56L

#' Assign each row of a heat map to a vertical block so a tall matrix becomes
#' several columns of the SAME cell height instead of one over-long page.
#'
#' @param row_levels ordered row identities (the y factor levels)
#' @return a data frame of row_id -> block label, balanced across blocks
m12_row_blocks <- function(row_levels, max_rows = HEATMAP_MAX_ROWS_PER_BLOCK) {
  n <- length(row_levels)
  k <- max(1L, ceiling(n / max_rows))
  per <- ceiling(n / k)
  tibble::tibble(row_id = row_levels,
                 block = factor(sprintf("%d of %d", ((seq_len(n) - 1L) %/% per) + 1L, k)))
}

#' Save a tile heat map with an exactly specified cell size.
#'
#' @param n_cols,n_rows tile counts; the panel becomes n * cell, exactly.
#' Fonts are asserted on the ggplot BEFORE conversion, because the assertion
#' walks theme elements and a gtable no longer carries them.
m12_save_heatmap <- function(name, p, n_cols, n_rows,
                             cell_w_mm = HEATMAP_CELL_W_MM,
                             cell_h_mm = HEATMAP_CELL_H_MM, dpi = 300) {
  m12_assert_arial()
  m12_assert_font_sizes(p, name)
  fix <- function(q) egg::set_panel_size(
    q,
    width  = grid::unit(n_cols * cell_w_mm, "mm"),
    height = grid::unit(n_rows * cell_h_mm, "mm"))
  # Two passes: the first learns the canvas width the fixed panel implies, the
  # second wraps the title to it. A title never widens the gtable - it is
  # clipped instead - so the width from pass one is the width that is drawn.
  w_mm <- grid::convertWidth(sum(fix(p)$widths), "mm", valueOnly = TRUE)
  gt   <- fix(m12_wrap_labs(p, w_mm))
  h_mm <- grid::convertHeight(sum(gt$heights), "mm", valueOnly = TRUE)
  stopifnot(is.finite(w_mm), is.finite(h_mm), w_mm > 0, h_mm > 0)
  # The whole point of this writer: one tile is the same size in every figure.
  # Faceted gtables name their panels "panel-1-1", "panel-1-2", ... so match on
  # the prefix; every panel must be identical, and each is checked.
  panels <- grep("^panel", gt$layout$name)
  stopifnot(length(panels) > 0)
  pw <- vapply(unique(gt$layout$l[panels]),
               function(i) grid::convertWidth(gt$widths[i], "mm", valueOnly = TRUE),
               numeric(1))
  if (any(abs(pw - n_cols * cell_w_mm) > 0.05))
    stop(sprintf("HEATMAP CELL ASSERTION FAILED: %s panel width %s mm, expected %.2f mm",
                 name, paste(sprintf("%.2f", pw), collapse = "/"), n_cols * cell_w_mm))
  res <- save_validation_figure(name, gt, width_mm = w_mm, height_mm = h_mm, dpi = dpi)
  vsay(sprintf("    cells: %d x %d at %.2f x %.2f mm", n_cols, n_rows, cell_w_mm, cell_h_mm))
  invisible(res)
}

#' Save a figure whose PANEL is exactly square, sizing the page around it.
#'
#' For a scatter of NHANES effect against external effect the two axes carry the
#' SAME quantity in the SAME units, so a non-square panel silently exaggerates
#' one of them: the identity line stops being the 45-degree diagonal a reader
#' assumes it is. Pinning the panel with egg::set_panel_size() and sizing the
#' device to the resulting gtable is the only way to guarantee it, since ggplot
#' otherwise stretches the panel to fill whatever canvas it is handed.
m12_save_square_panel <- function(name, p, side_mm = 48, dpi = 300) {
  m12_assert_arial()
  m12_assert_font_sizes(p, name)
  side <- grid::unit(side_mm, "mm")
  fix  <- function(q) egg::set_panel_size(q, width = side, height = side)
  # Pass one learns the canvas width, pass two wraps the title into it.
  w_mm <- grid::convertWidth(sum(fix(p)$widths), "mm", valueOnly = TRUE)
  gt   <- fix(m12_wrap_labs(p, w_mm))
  h_mm <- grid::convertHeight(sum(gt$heights), "mm", valueOnly = TRUE)
  stopifnot(is.finite(w_mm), is.finite(h_mm), w_mm > 0, h_mm > 0)
  # Assert the invariant rather than trusting it.
  panels <- grep("^panel", gt$layout$name)
  stopifnot(length(panels) > 0)
  pw <- grid::convertWidth(gt$widths[unique(gt$layout$l[panels])],   "mm", valueOnly = TRUE)
  ph <- grid::convertHeight(gt$heights[unique(gt$layout$t[panels])], "mm", valueOnly = TRUE)
  if (any(abs(pw - side_mm) > 0.05) || any(abs(ph - side_mm) > 0.05))
    stop(sprintf("SQUARE PANEL ASSERTION FAILED: %s panel is %s x %s mm, expected %.2f square",
                 name, paste(sprintf("%.2f", pw), collapse = "/"),
                 paste(sprintf("%.2f", ph), collapse = "/"), side_mm))
  res <- save_validation_figure(name, gt, width_mm = w_mm, height_mm = h_mm, dpi = dpi)
  vsay(sprintf("    panel: %.2f x %.2f mm (square)", pw[1], ph[1]))
  invisible(res)
}

#' Rasterize a PDF to PNG with pdftoppm.
#'
#' Deliberately NOT a fresh ggsave to PNG: `systemfonts` resolves "Arial" to
#' DejaVuSans on this system, so re-drawing would embed a different typeface
#' than the PDF. Rasterizing the PDF guarantees the report image and the
#' publication figure are typographically identical.
save_validation_png <- function(pdf_path, dpi = 300) {
  stem <- tools::file_path_sans_ext(basename(pdf_path))
  out  <- file.path(VAL_FIG_PNG, stem)
  tool <- Sys.which("pdftoppm")
  if (!nzchar(tool)) {
    warning("pdftoppm not found; report PNG not generated for ", stem)
    return(NA_character_)
  }
  system2(tool, c("-png", "-r", dpi, "-singlefile", shQuote(pdf_path), shQuote(out)),
          stdout = FALSE, stderr = FALSE)
  png <- paste0(out, ".png")
  if (!file.exists(png)) { warning("PNG conversion failed for ", stem); return(NA_character_) }
  png
}

# --- Small helpers -----------------------------------------------------------
#' Wilson-free exact binomial test with CI, returned tidily. NA-safe for n = 0.
m12_binom <- function(k, n, p = 0.5) {
  if (is.na(n) || n < 1) return(tibble(k = k, n = n, pct = NA_real_,
                                       ci_lo = NA_real_, ci_hi = NA_real_, p_value = NA_real_))
  bt <- stats::binom.test(k, n, p)
  tibble(k = k, n = n, pct = 100 * k / n,
         ci_lo = 100 * bt$conf.int[1], ci_hi = 100 * bt$conf.int[2],
         p_value = bt$p.value)
}

#' CLR with the SAME rule as module 0: integer pseudocount eps = 1 applied to
#' COUNTS, then ln(x / geometric mean). Applied row-wise (per sample).
m12_clr <- function(mat, pseudocount = 1) {
  m <- as.matrix(mat); storage.mode(m) <- "double"
  m[is.na(m)] <- 0
  m <- m + pseudocount
  gm <- exp(rowMeans(log(m)))
  out <- log(m / gm)
  stopifnot(all(is.finite(out)))
  out
}

vsay("validation config loaded (FDR<", FDR_THRESHOLD, ", prevalence>=", PREVALENCE_THRESHOLD,
     ", adults>=", MIN_ADULT_AGE, ")")
