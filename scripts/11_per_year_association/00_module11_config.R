# =============================================================================
# MODULE 11 — shared configuration
# 00_module11_config.R
# =============================================================================
# source()d by every module-11 R analysis script (04-07, 09). Follows the
# copy-of-record pattern of module 12's 00_validation_config.R: paths,
# significance thresholds, the project phylum palette, the ArialMT 6-pt figure
# conventions, and loaders for taxonomy / prevalence / main-effect results.
#
# Everything module 11 writes lives under OUT_ROOT; the pooled Module 1 results
# and the SQLite databases are opened READ-ONLY.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

set.seed(42)

# --- Paths -------------------------------------------------------------------
# PROJECT_ROOT: absolute path of your local clone of this repository.
# Export PROJECT_ROOT before running, or run every script from the repository
# root (the default below then resolves to the working directory).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", unset = getwd())
M11_BASE     <- PROJECT_ROOT
M11_OUT      <- file.path(M11_BASE, "results/analyses_results/11_per_year_association")
M11_POOLED   <- file.path(M11_BASE, "results")   # Module 1 pooled results (READ-ONLY)
M11_DB       <- file.path(M11_BASE, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite")

M11_CMP      <- file.path(M11_OUT, "comparisons")
M11_TABLES   <- file.path(M11_OUT, "tables")
M11_FIGURES  <- file.path(M11_OUT, "figures")
# Figure tiers: the four MAIN figures (manuscript designation Figure S16A-C)
# live directly in figures/; every other figure is an additional accessory
# figure, saved to figures/additional_accessory_figures/ with an "optional_"
# filename prefix. Routing and naming are by figure name (single source of
# truth here), so no per-call-site flags are needed.
M11_MAIN_FIGURES <- c("figS16A_left_direction_reproduced_by_cycle",
                      "figS16A_right_validated_by_cycle",
                      "figS16B_pooled_vs_per_cycle_beta",
                      "figS16C_beta_F_vs_G_pooled_sig")
M11_FIGURES_ACCESSORY <- file.path(M11_FIGURES, "additional_accessory_figures")
m11_figure_dir <- function(name)
  if (name %in% M11_MAIN_FIGURES) M11_FIGURES else M11_FIGURES_ACCESSORY
m11_figure_filename <- function(name) {
  if (name %in% M11_MAIN_FIGURES || startsWith(name, "optional_")) name
  else paste0("optional_", name)
}
M11_REPORTS  <- file.path(M11_OUT, "reports")
M11_AVAIL    <- file.path(M11_CMP, "variable_availability")

for (d in c(M11_CMP, M11_TABLES, M11_FIGURES, M11_FIGURES_ACCESSORY, M11_REPORTS, M11_AVAIL,
            file.path(M11_CMP, "F_vs_G"), file.path(M11_CMP, "pooled_vs_F"),
            file.path(M11_CMP, "pooled_vs_G"), file.path(M11_CMP, "pooled_vs_F_vs_G"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

M11_ANALYSES <- c("1_demoWAS", "2_oradWAS", "3_exWAS", "4_pheWAS", "5_outWAS")
M11_NORM     <- "clr"
M11_CYCLES   <- c("F", "G")

# --- Significance thresholds (user decision 2026-08-30) ----------------------
M11_THRESH_PRIMARY <- 0.05   # headline cutoff, displayed as "FDR<0.05" / "q<0.05"; code compares
                             # <= but no FDR/q value equals the cutoff exactly, so the sets match
M11_THRESH_CONTEXT <- 0.10   # also tabulated for context
M11_PREV_THRESHOLD <- 0.10   # primary correlation subset: pooled prevalence > 10%
                             # (matches module 12's PREVALENCE_THRESHOLD)

# --- Fit-stability guard ------------------------------------------------------
# Survey-weighted quasibinomial fits can quasi-separate on sparse taxa, producing
# astronomically large log-odds coefficients (|beta| up to ~1e14 observed) whose
# sign and magnitude are numerical noise. Such fits are flagged UNSTABLE and
# excluded from beta correlations / direction assessment (always inventoried,
# never silently dropped). |log-odds| > 15 means OR > 3e6 — beyond any plausible
# biology. Linear fits get only a numerical-explosion guard (betas are in native
# outcome units, which legitimately vary by orders of magnitude across variables).
M11_LOGISTIC_ANALYSES <- c("2_oradWAS", "5_outWAS")
M11_LOGIT_ABS_MAX     <- 15      # |beta| and SE bound for logistic fits
M11_LINEAR_ABS_MAX    <- 1e6    # numerical-explosion bound for linear fits

m11_stable_fit <- function(beta, se, is_logistic) {
  ok <- is.finite(beta) & is.finite(se) & se > 0
  lim_b <- ifelse(is_logistic, M11_LOGIT_ABS_MAX, M11_LINEAR_ABS_MAX)
  lim_s <- ifelse(is_logistic, M11_LOGIT_ABS_MAX, M11_LINEAR_ABS_MAX)
  ok & abs(beta) <= lim_b & se <= lim_s
}

# --- Phylum palette ----------------------------------------------------------
# Canonical 17-entry form, copied from
# scripts/8_network_analyses/2_network.R:112 via the module 12 copy-of-record
# (scripts/12_external_data_validation/10_validation_analysis/00_validation_config.R).
# Do not re-order or re-code: module-11 figures must be colour-comparable with
# modules 4/5/8/9/12.
# ONE deliberate deviation (colour-vision-accessibility update, 2026-09-17):
# Acidobacteria moved
# from Kelly "orange yellow" #F6A600 (too close to Spirochaetae's yellow
# #F3C300 under colour-vision deficiency) to Kelly "olive green" #2B3D26 —
# a manuscript-wide phylum recolouring. (#2B3D26 also codes Candidate division
# SR1 below, which never appears in module-11 figures, so no collision here.)
phylum_colors <- c(
  "Firmicutes"             = "#F38400",  # kelly_3   - orange
  "Bacteroidetes"          = "#0067A5",  # kelly_10  - blue
  "Actinobacteria"         = "#8DB600",  # kelly_17  - yellow green
  "Proteobacteria"         = "#E68FAC",  # kelly_9   - pink
  "Fusobacteria"           = "#BE0032",  # kelly_5   - red
  "Spirochaetae"           = "#F3C300",  # kelly_1   - yellow
  "Cyanobacteria"          = "#875692",  # kelly_2   - purple
  "Acidobacteria"          = "#2B3D26",  # kelly_20  - olive green (was #F6A600; see note)
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
m11_phylum_std <- function(v) {
  s <- trimws(as.character(v))
  s[is.na(s)] <- "NA"
  s[s == ""]  <- "unclassified"
  s[!(s %in% names(phylum_colors))] <- "unclassified"
  s
}

# Effect-direction colours (module 4, verbatim)
association_colors <- c("Negative" = "#0072B2", "Positive" = "#CC6677")

# --- Colour-blind-friendly palette (all NON-phylum colours) -------------------
# Policy: phylum colours keep the project's predefined hex palette above;
# every other discrete colour comes from the Okabe-Ito colour-blind-safe
# palette (as shipped in grafify's `okabe_ito`); sequential scales use viridis;
# diverging scales use the Okabe-Ito blue <-> white <-> vermillion anchors.
M11_CBF <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
             yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
             purple = "#CC79A7", grey = "#999999", black = "#000000")
# Shared validation-class colours (used identically in fig3 / fig8 / fig10)
M11_VAL_COLS <- c(both = unname(M11_CBF["vermillion"]),
                  F_only = unname(M11_CBF["blue"]),
                  G_only = unname(M11_CBF["green"]),
                  samedir = unname(M11_CBF["orange"]),
                  inconsistent = unname(M11_CBF["grey"]))
M11_DIVERGING <- c(low = unname(M11_CBF["blue"]), mid = "white",
                   high = unname(M11_CBF["vermillion"]))

# --- Fonts / figure conventions ----------------------------------------------
# ArialMT at base size 6, PDF-first, PNG rasterized from the PDF (module 12
# guardrails: never let a raster device silently substitute DejaVu for Arial).
PDF_FAMILY      <- "ArialMT"
THEME_BASE_SIZE <- 6
MAX_FIG_INCHES  <- 7

m11_register_fonts <- function() {
  suppressMessages(try(extrafont::loadfonts(device = "pdf", quiet = TRUE), silent = TRUE))
  if (!PDF_FAMILY %in% names(grDevices::pdfFonts())) {
    try(extrafont::font_import(path = file.path(M11_BASE, "data/fonts"),
                               pattern = "arial.*\\.ttf$", prompt = FALSE), silent = TRUE)
    suppressMessages(try(extrafont::loadfonts(device = "pdf", quiet = TRUE), silent = TRUE))
  }
  invisible(TRUE)
}

m11_assert_arial <- function() {
  if (!PDF_FAMILY %in% names(grDevices::pdfFonts()))
    stop("FONT ASSERTION FAILED: '", PDF_FAMILY, "' is not in pdfFonts(). ",
         "Refusing to emit figures in a substituted font.")
  invisible(TRUE)
}

# Based on module 5.5's theme_compact(), with EVERY text element pinned to
# exactly `base` (6 pt) per user requirement — titles, axis ticks, legend,
# strips, and captions all render at the same 6 pt ArialMT.
theme_compact <- function(base = THEME_BASE_SIZE, family = PDF_FAMILY) {
  egg::theme_article(base_size = base, base_family = family) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = base, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = base, color = "grey30"),
      plot.caption     = ggplot2::element_text(size = base, color = "grey45", hjust = 0),
      axis.title       = ggplot2::element_text(size = base),
      axis.text        = ggplot2::element_text(size = base),
      legend.title     = ggplot2::element_text(size = base),
      legend.text      = ggplot2::element_text(size = base),
      legend.key.size  = grid::unit(0.3, "lines"),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = base, face = "bold"),
      plot.margin      = ggplot2::margin(2, 4, 2, 2, unit = "pt")
    )
}

#' Write a ggplot to PDF (ArialMT, mm sizing, <= 7x7 in enforced), then render
#' a matching PNG from that PDF with pdftoppm (module 12 idiom).
m11_save_figure <- function(name, plot_obj, width_mm, height_mm, dpi = 300,
                            allow_oversize = FALSE) {
  m11_assert_arial()
  w_in <- width_mm / 25.4; h_in <- height_mm / 25.4
  if (!allow_oversize && (w_in > MAX_FIG_INCHES || h_in > MAX_FIG_INCHES))
    stop(sprintf("FIGURE TOO LARGE: %s is %.2f x %.2f in (max %d).",
                 name, w_in, h_in, MAX_FIG_INCHES))
  fname <- m11_figure_filename(name)
  pdf_path <- file.path(m11_figure_dir(name), paste0(fname, ".pdf"))
  grDevices::pdf(pdf_path, width = w_in, height = h_in, family = PDF_FAMILY,
                 colormodel = "rgb", pointsize = 5, useDingbats = FALSE)
  if (inherits(plot_obj, "gg")) print(plot_obj) else grid::grid.draw(plot_obj)
  grDevices::dev.off()
  # PNG rasterized from the PDF so the typeface is identical
  png_out <- file.path(m11_figure_dir(name), fname)
  tool <- Sys.which("pdftoppm")
  if (nzchar(tool)) {
    system2(tool, c("-png", "-r", dpi, "-singlefile", shQuote(pdf_path), shQuote(png_out)),
            stdout = FALSE, stderr = FALSE)
  }
  cat(sprintf("  figure: %s (%.2f x %.2f in)\n", basename(pdf_path), w_in, h_in))
  invisible(pdf_path)
}

#' Multi-page variant: one PDF with one plot per page (catalog figures).
#' PNG for report embedding is rasterized from PAGE 1 only.
m11_save_figure_multipage <- function(name, plot_list, width_mm, height_mm, dpi = 300) {
  m11_assert_arial()
  w_in <- width_mm / 25.4; h_in <- height_mm / 25.4
  if (w_in > MAX_FIG_INCHES || h_in > MAX_FIG_INCHES)
    stop(sprintf("FIGURE TOO LARGE: %s is %.2f x %.2f in (max %d).",
                 name, w_in, h_in, MAX_FIG_INCHES))
  fname <- m11_figure_filename(name)
  pdf_path <- file.path(m11_figure_dir(name), paste0(fname, ".pdf"))
  grDevices::pdf(pdf_path, width = w_in, height = h_in, family = PDF_FAMILY,
                 colormodel = "rgb", pointsize = 5, useDingbats = FALSE, onefile = TRUE)
  for (p in plot_list) print(p)
  grDevices::dev.off()
  png_out <- file.path(m11_figure_dir(name), fname)
  tool <- Sys.which("pdftoppm")
  if (nzchar(tool)) {
    system2(tool, c("-png", "-r", dpi, "-f", "1", "-l", "1", "-singlefile",
                    shQuote(pdf_path), shQuote(png_out)),
            stdout = FALSE, stderr = FALSE)
  }
  cat(sprintf("  figure: %s (%d pages, %.2f x %.2f in)\n",
              basename(pdf_path), length(plot_list), w_in, h_in))
  invisible(pdf_path)
}

# --- Taxonomy loader (Route B: variable_names_epcf, self-contained) ----------
# Variable.Description is "Domain;Phylum;Class;Order;Family;Genus" with literal
# "NA" for missing ranks. Keys on the association variable name
# ("RSV_genusN_relative"). Genus label preference follows module 8:
# coalesce(GOLD parsed_genus, taxonomy Genus, otu id); unnamed -> "NA; <otu>".
m11_load_taxonomy <- function() {
  cache <- file.path(M11_AVAIL, "taxonomy_annotations.rds")
  if (file.exists(cache)) return(readRDS(cache))
  con <- DBI::dbConnect(RSQLite::SQLite(), M11_DB)
  on.exit(DBI::dbDisconnect(con))
  tax <- DBI::dbGetQuery(con, "
    SELECT DISTINCT \"Variable.Name\" AS variable, \"Variable.Description\" AS lineage
    FROM variable_names_epcf
    WHERE \"Data.File.Name\" LIKE 'DADA2RSV_GENUS_RELATIVE_%'
      AND \"Variable.Name\" LIKE 'RSV_genus%_relative'") %>%
    distinct(variable, .keep_all = TRUE) %>%
    mutate(parts  = strsplit(lineage, ";"),
           Phylum = vapply(parts, function(x) trimws(x[2]), character(1)),
           Genus  = vapply(parts, function(x) if (length(x) >= 6) trimws(x[6]) else NA_character_, character(1)),
           across(c(Phylum, Genus), ~ ifelse(.x %in% c("NA", ""), NA_character_, .x)),
           otu = str_remove(variable, "_relative$")) %>%
    select(variable, otu, Phylum, Genus)
  gold_path <- file.path(M11_BASE,
    "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv")
  if (file.exists(gold_path)) {
    gold <- read_csv(gold_path, show_col_types = FALSE) %>%
      transmute(otu, parsed_genus)
    tax <- tax %>% left_join(gold, by = "otu")
  } else {
    tax$parsed_genus <- NA_character_
  }
  tax <- tax %>%
    mutate(phylum_std  = m11_phylum_std(Phylum),
           genus_label = coalesce(parsed_genus, Genus),
           genus_label = if_else(is.na(genus_label) | genus_label == "unclassified",
                                 paste0("NA; ", otu), genus_label))
  saveRDS(tax, cache)
  tax
}

# --- Prevalence loader --------------------------------------------------------
# Pooled per-RSV prevalence from module 10's ready-made QC CSV (bare RSV_genusN);
# per-cycle prevalence computed once from the untransformed RELATIVE tables
# (mean(x > 0) per column) and cached. All inputs READ-ONLY.
m11_load_prevalence <- function() {
  cache <- file.path(M11_AVAIL, "otu_prevalence.rds")
  if (file.exists(cache)) return(readRDS(cache))
  con <- DBI::dbConnect(RSQLite::SQLite(), M11_DB)
  on.exit(DBI::dbDisconnect(con))
  prev_cycle <- function(tbl_name, cyc) {
    d <- DBI::dbGetQuery(con, sprintf('SELECT * FROM "%s"', tbl_name))
    d$SEQN <- NULL
    tibble(variable = names(d),
           !!paste0("prevalence_", cyc) := vapply(d, function(x) mean(x > 0, na.rm = TRUE), numeric(1)),
           !!paste0("n_samples_", cyc)  := nrow(d))
  }
  pf <- prev_cycle("DADA2RSV_GENUS_RELATIVE_F", "F")
  pg <- prev_cycle("DADA2RSV_GENUS_RELATIVE_G", "G")
  prev <- full_join(pf, pg, by = "variable") %>%
    mutate(otu = str_remove(variable, "_relative$"))
  m10 <- file.path(M11_BASE,
    "results/analyses_results/10_prediction_analyses_out/intermediate/qc/qc_microbiome_prevalence_summary.csv")
  if (file.exists(m10)) {
    pooled <- read_csv(m10, show_col_types = FALSE) %>%
      transmute(otu = taxon_clean, prevalence_pooled = prevalence)
    prev <- prev %>% left_join(pooled, by = "otu")
  }
  # Fallback / fill: pooled prevalence from the per-cycle counts if CSV absent
  prev <- prev %>%
    mutate(prevalence_pooled = coalesce(
      prevalence_pooled,
      (prevalence_F * n_samples_F + prevalence_G * n_samples_G) / (n_samples_F + n_samples_G)))
  saveRDS(prev, cache)
  prev
}

# --- Result loaders -----------------------------------------------------------
m11_tidied_path <- function(analysis, where = c("pooled", "F", "G")) {
  where <- match.arg(where)
  root <- switch(where,
    pooled = file.path(M11_POOLED, paste0(analysis, "_out"), paste0("result_", M11_NORM)),
    F      = file.path(M11_OUT, "per_year_F", paste0(analysis, "_out"), paste0("result_", M11_NORM)),
    G      = file.path(M11_OUT, "per_year_G", paste0(analysis, "_out"), paste0("result_", M11_NORM)))
  file.path(root, paste0(analysis, "_", M11_NORM, "_tidied_complete.rds"))
}

#' Load MAIN-EFFECT rows (fdr_corrected == TRUE) of one aggregated tidied file.
m11_load_main_effects <- function(analysis, where) {
  p <- m11_tidied_path(analysis, where)
  if (!file.exists(p)) return(NULL)
  x <- readRDS(p)
  if (!"fdr_corrected" %in% names(x)) return(NULL)
  x <- x %>% filter(fdr_corrected == TRUE)
  if (!"q.value" %in% names(x)) x$q.value <- NA_real_
  x %>% transmute(dependent_var, exposure, term,
                  estimate, std.error, p.value, p.value.fdr, q.value,
                  n_obs, cycle_mode, formula_used, n_covariates)
}

# Direction / correlation helpers (r never reported without n by construction:
# both are computed from the same frame in the calling scripts)
m11_pearson <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  n <- sum(ok)
  if (n < 3) return(list(r = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(a[ok], b[ok], method = "pearson"))
  list(r = unname(ct$estimate), p = ct$p.value, n = n)
}

m11_dir_agree <- function(a, b) {
  ok <- is.finite(a) & is.finite(b) & a != 0 & b != 0
  n <- sum(ok)
  if (n == 0) return(list(n = 0, n_agree = 0, pct = NA_real_))
  k <- sum(sign(a[ok]) == sign(b[ok]))
  list(n = n, n_agree = k, pct = 100 * k / n)
}
