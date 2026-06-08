#!/usr/bin/env Rscript
# =============================================================================
# 5.5_smoking_analyses_additional.R
#
# Smoking-exposure multicollinearity diagnostics for the NHANES oral microbiome
# analysis. Produces figures and tables documenting:
#   (1) the 24 NHANES smoking/combustion-related variables span 4 substantively
#       different measurement modalities;
#   (2) mutual adjustment is impractical (joint complete-case N collapses;
#       VIF inflates);
#   (3) the per-variable microbial signatures are partially non-overlapping in
#       a way that tracks modality.
#
# All input reads are read-only. Every artifact is written under
#   <PROJECT_ROOT>/results/analyses_results/5.5_smoking_analyses_out_additional/
#
# Environment: R >= 4.5 with phyloseq, dplyr, tidyr, ggplot2, egg, patchwork,
# scales, grid, gridExtra, ComplexHeatmap, circlize, ggrepel, viridis,
# RColorBrewer, car.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
#
# Variable transformation policy:
#   For each continuous biomarker, zero values are handled by adding an
#   analyte-specific pseudocount = min(x[x>0])/2, followed by natural-log
#   transform. PCA and VIF additionally z-score. Correlations are scale-
#   invariant and use only log(x + pc) without standardization.
#   SMQ_current_ever_never is treated as an ordinal 0/1/2 integer (never =
#   reference) - no log transform.
#
# Run:
#   Rscript scripts/5.5_smoking_analyses/5.5_smoking_analyses_additional.R
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

# ---- SECTION 0 - preamble ----
set.seed(42)
SCRIPT_VERSION <- "0.1.0"
RUN_TIMESTAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Redirect R's default graphics device to a null PDF - anything that implicitly
# opens a screen device on this headless Linux node (e.g. some
# ComplexHeatmap/grid internals) will go to /dev/null instead of leaving a
# stray Rplots.pdf at the cwd.
options(device = function(...) grDevices::pdf(NULL))

base_path <- PROJECT_ROOT

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(egg)
  library(patchwork)
  library(scales)
  library(grid)
  library(gridExtra)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(viridis)
  library(RColorBrewer)
  library(car)
})

# ---- paths ----
viz_out_root <- file.path(base_path, "results/analyses_results")
viz_out_path <- file.path(viz_out_root,
                          "5.5_smoking_analyses_out_additional")

# Output-root guard - never write outside this tree.
stopifnot(
  grepl("5.5_smoking_analyses_out_additional",
        viz_out_path, fixed = TRUE)
)

PLOTS_DIR   <- file.path(viz_out_path, "plots")
TABLES_DIR  <- file.path(viz_out_path, "tables")
INPUTS_DIR  <- file.path(viz_out_path, "inputs")
for (d in c(viz_out_path, PLOTS_DIR, TABLES_DIR, INPUTS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(viz_out_path, "run.log")
log_con  <- file(LOG_FILE, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(),                 silent = TRUE)
  try(close(log_con),         silent = TRUE)
}, add = TRUE)

cat("=== 5.5_smoking_analyses_additional.R ===\n")
cat("script version:", SCRIPT_VERSION, "\n")
cat("run timestamp:",  RUN_TIMESTAMP, "\n")
cat("base_path:",      base_path, "\n")
cat("viz_out_path:",   viz_out_path, "\n")
cat("R version:",      R.version.string, "\n\n")

# ---- input paths (read-only) ----
INPUTS <- list(
  var_mapping_csv   = file.path(base_path,
                                "results/intermediate/ubiome_variable_mapping.csv"),
  host_meta_rds     = file.path(base_path,
                                "results/intermediate/ubiome_relative_none_updated.rds"),
  clr_phyloseq_rds  = file.path(base_path,
                                "results/analyses_results/2_preprocess_db_n_phyloseq_out",
                                "intermediate/ubiome_relative_clr.rds"),
  tidied_rds        = file.path(base_path,
                                "results/3_exWAS_out/result_clr/3_exWAS_clr_tidied_complete.rds"),
  glanced_rds       = file.path(base_path,
                                "results/3_exWAS_out/result_clr/3_exWAS_clr_glanced_complete.rds"),
  rsq_rds           = file.path(base_path,
                                "results/3_exWAS_out/result_clr/3_exWAS_clr_rsq_complete.rds"),
  provenance_txt    = file.path(base_path,
                                "results/3_exWAS_out/result_clr/3_exWAS_clr_aggregation_summary.txt"),
  # Precomputed microbial-signature host-host correlation matrix (weighted
  # inverse-variance Pearson on per-taxon beta across shared microbes; FDR via
  # BH). Produced by scripts/7_microbial_signature_heatmap/. Read-only.
  sig_corr_rds      = file.path(base_path,
                                "results/analyses_results/7_microbial_signature_heatmap_out",
                                "intermediate/microbial_signature/exWAS_clr_all_results",
                                "correlation_matrix.rds")
)

# Hard guard: no input path may match the disallowed CLR file or the
# 4_association tree.
for (p in INPUTS) {
  stopifnot(
    !grepl("ubiome_relative_clr_updated\\.rds$", p),
    !grepl("/4_association_phyloseq_analyses_out/intermediate/ubiome_relative_clr_updated",
           p, fixed = TRUE)
  )
}

# ---- provenance: md5 + mtime of every input read ----
prov <- do.call(rbind, lapply(names(INPUTS), function(nm) {
  fp <- INPUTS[[nm]]
  exists_ok <- file.exists(fp)
  data.frame(
    role      = nm,
    path      = fp,
    exists    = exists_ok,
    bytes     = if (exists_ok) unname(file.info(fp)$size) else NA_real_,
    mtime     = if (exists_ok) format(file.info(fp)$mtime) else NA_character_,
    md5       = if (exists_ok) unname(tools::md5sum(fp)) else NA_character_,
    stringsAsFactors = FALSE
  )
}))
write.table(prov, file.path(INPUTS_DIR, "input_provenance.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
cat("[provenance] wrote input_provenance.txt; missing inputs:\n")
print(prov[!prov$exists, c("role", "path"), drop = FALSE])
stopifnot(all(prov$exists))   # abort if any input missing

# ---- SECTION 2 - build the 24-var modality table ----
TARGET_VARS <- c(
  # Questionnaire (Q)
  "SMQ_current_ever_never",
  # Urinary PAH metabolites (U-PAH)
  "URXP01", "URXP03", "URXP04", "URXP05",
  # Urinary metals & mercapturates (U-MM)
  "URXUCD", "URX1DC",
  # Blood VOCs + combustion analytes (B-VOC)
  "LBXCOT", "LBX2DF", "LBXVFN", "LBXVBZ",
  "LBXV06", "LBXVEB", "LBXVXY", "LBXVCB",
  "LBXV1D", "LBXV3B", "LBXVMC", "LBXVDM",
  "LBXVTC", "LBXV4C", "LBXV2A", "LBXV2T", "LBXV4E"
)
stopifnot(length(TARGET_VARS) == 24L, !anyDuplicated(TARGET_VARS))

modality_for_var <- c(
  SMQ_current_ever_never = "Q",
  URXP01 = "U-PAH",  URXP03 = "U-PAH",  URXP04 = "U-PAH",  URXP05 = "U-PAH",
  URXUCD = "U-MM",   URX1DC = "U-MM",
  LBXCOT = "B-VOC",  LBX2DF = "B-VOC",  LBXVFN = "B-VOC",  LBXVBZ = "B-VOC",
  LBXV06 = "B-VOC",  LBXVEB = "B-VOC",  LBXVXY = "B-VOC",  LBXVCB = "B-VOC",
  LBXV1D = "B-VOC",  LBXV3B = "B-VOC",  LBXVMC = "B-VOC",  LBXVDM = "B-VOC",
  LBXVTC = "B-VOC",  LBXV4C = "B-VOC",  LBXV2A = "B-VOC",  LBXV2T = "B-VOC",
  LBXV4E = "B-VOC"
)
stopifnot(setequal(names(modality_for_var), TARGET_VARS))

MODALITY_LEVELS <- c("Q", "U-PAH", "U-MM", "B-VOC")
MODALITY_LABEL  <- c(Q = "Questionnaire",
                     `U-PAH` = "Urinary PAH metabolites",
                     `U-MM`  = "Urinary metals & mercapturates",
                     `B-VOC` = "Blood VOCs & combustion analytes")
# Okabe-Ito-derived; visually distinct on B/W print
MODALITY_COLOR  <- c(Q     = "#0072B2",
                     `U-PAH` = "#009E73",
                     `U-MM`  = "#E69F00",
                     `B-VOC` = "#CC79A7")

# Read canonical var_name -> var_description mapping (640 rows).
var_map <- utils::read.csv(INPUTS$var_mapping_csv, stringsAsFactors = FALSE)
stopifnot(all(c("var_name", "var_description") %in% names(var_map)))
miss_desc <- setdiff(TARGET_VARS, var_map$var_name)
if (length(miss_desc)) stop("missing var_descriptions for: ",
                            paste(miss_desc, collapse = ", "))
desc_lookup <- setNames(var_map$var_description, var_map$var_name)

modality_tbl <- data.frame(
  var_name        = TARGET_VARS,
  var_description = unname(desc_lookup[TARGET_VARS]),
  modality_code   = unname(modality_for_var[TARGET_VARS]),
  stringsAsFactors = FALSE
)
modality_tbl$modality_label <- unname(MODALITY_LABEL[modality_tbl$modality_code])
modality_tbl$plot_order     <- seq_len(nrow(modality_tbl))
write.csv(modality_tbl,
          file.path(INPUTS_DIR, "smoking_variable_modality_mapping.csv"),
          row.names = FALSE)
cat("\n[section 2] modality mapping written (",
    nrow(modality_tbl), " rows):\n", sep = "")
print(modality_tbl[, c("var_name", "modality_code", "var_description")])
cat("\n")

# =============================================================================
# Helpers used across plotting sections.
# =============================================================================

THEME_BASE_SIZE   <- 6
THEME_BASE_FAMILY <- "ArialMT"
PDF_FAMILY        <- "ArialMT"  # base pdf() device has ArialMT in pdfFonts()

theme_compact <- function(base = THEME_BASE_SIZE,
                          family = THEME_BASE_FAMILY) {
  egg::theme_article(base_size = base, base_family = family) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = base + 1, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = base, color = "grey30"),
      plot.caption     = ggplot2::element_text(size = base - 1, color = "grey45",
                                               hjust = 0),
      axis.title       = ggplot2::element_text(size = base),
      axis.text        = ggplot2::element_text(size = base - 1),
      legend.title     = ggplot2::element_text(size = base),
      legend.text      = ggplot2::element_text(size = base - 1),
      legend.key.size  = grid::unit(0.3, "lines"),
      legend.margin    = ggplot2::margin(0, 0, 0, 0),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = base, face = "bold"),
      plot.margin      = ggplot2::margin(2, 4, 2, 2, unit = "pt")
    )
}

open_pdf <- function(file_name, width_in, height_in) {
  fp <- file.path(PLOTS_DIR, file_name)
  pdf(file = fp,
      width = width_in, height = height_in,
      family = PDF_FAMILY, colormodel = "rgb",
      pointsize = THEME_BASE_SIZE, useDingbats = FALSE)
  fp
}

# Map a value vector to an analyte-specific log transform per the variable-
# transformation policy. Zeros become log(pc); negatives drop to NA.
log_transform_with_pc <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  pos <- x[is.finite(x) & x > 0]
  if (!length(pos)) {
    return(list(x_log = rep(NA_real_, length(x)), pc = NA_real_,
                n_pos = 0L, min_pos = NA_real_))
  }
  pc <- min(pos) / 2
  x_log <- ifelse(is.finite(x) & x >= 0, log(x + pc), NA_real_)
  list(x_log = x_log, pc = pc, n_pos = length(pos), min_pos = min(pos))
}

# =============================================================================
# SECTION 3 - host-level summaries
#   - pairwise_complete_case_N.csv + correlation/N heatmap PDF
# =============================================================================
cat("=== SECTION 3 - host-level summaries ===\n")

host_ps <- readRDS(INPUTS$host_meta_rds)
host_sd <- as(sample_data(host_ps), "data.frame")
cat("host_ps nsamples=", nsamples(host_ps),
    "  ntaxa=", ntaxa(host_ps),
    "  taxa_are_rows=", taxa_are_rows(host_ps), "\n", sep = "")
cat("host sample_data dim: ", paste(dim(host_sd), collapse = " x "), "\n")
host_sd$sample_id <- rownames(host_sd)

# ---- pull the 24 columns; build raw, log-transformed, and z matrices ----
miss_cols <- setdiff(TARGET_VARS, names(host_sd))
if (length(miss_cols)) stop("missing target columns in host sample_data: ",
                            paste(miss_cols, collapse = ", "))
host_raw <- host_sd[, TARGET_VARS, drop = FALSE]
rownames(host_raw) <- host_sd$sample_id

# SMQ_current_ever_never -> ordinal numeric (0/1/2); other 23 -> numeric raw
host_raw$SMQ_current_ever_never <- suppressWarnings(
  as.numeric(as.character(host_raw$SMQ_current_ever_never)))
for (v in setdiff(TARGET_VARS, "SMQ_current_ever_never")) {
  host_raw[[v]] <- suppressWarnings(as.numeric(as.character(host_raw[[v]])))
}

# Per-variable transformation table (pseudocount, n_pos, min_pos)
transform_meta <- do.call(rbind, lapply(TARGET_VARS, function(v) {
  if (v == "SMQ_current_ever_never") {
    return(data.frame(var_name = v, pc = NA_real_,
                      n_pos = sum(is.finite(host_raw[[v]])),
                      min_pos = NA_real_,
                      transform = "ordinal 0/1/2",
                      stringsAsFactors = FALSE))
  }
  z <- log_transform_with_pc(host_raw[[v]])
  data.frame(var_name = v, pc = z$pc, n_pos = z$n_pos,
             min_pos = z$min_pos,
             transform = "log(x + min(x>0)/2)",
             stringsAsFactors = FALSE)
}))
cat("\n[transform_meta]\n"); print(transform_meta)

# Build x_log matrix (rows = samples, cols = 24 vars)
x_log <- matrix(NA_real_,
                nrow = nrow(host_raw), ncol = length(TARGET_VARS),
                dimnames = list(rownames(host_raw), TARGET_VARS))
for (v in TARGET_VARS) {
  if (v == "SMQ_current_ever_never") {
    x_log[, v] <- as.numeric(host_raw[[v]])   # already ordinal 0/1/2
  } else {
    x_log[, v] <- log_transform_with_pc(host_raw[[v]])$x_log
  }
}

# ---- pairwise complete-case N matrix ----
cat("\n[pairwise_N] computing 24x24 complete-case N matrix...\n")
pairwise_N <- matrix(0L, nrow = 24, ncol = 24,
                     dimnames = list(TARGET_VARS, TARGET_VARS))
for (i in seq_along(TARGET_VARS)) {
  for (j in seq_along(TARGET_VARS)) {
    pairwise_N[i, j] <- sum(is.finite(host_raw[[i]]) &
                             is.finite(host_raw[[j]]))
  }
}
write.csv(as.data.frame(pairwise_N),
          file.path(TABLES_DIR, "pairwise_complete_case_N.csv"))
joint_N_all <- sum(stats::complete.cases(host_raw))
cat("[pairwise_N] joint complete-case across all 24 vars: N =",
    joint_N_all, "\n")

# ---- pairwise N ggplot (will be combined with Spearman heatmap below) ----
{
  pn_long <- as.data.frame(as.table(pairwise_N))
  names(pn_long) <- c("var_i", "var_j", "n")
  pn_long$var_i <- factor(as.character(pn_long$var_i), levels = TARGET_VARS)
  pn_long$var_j <- factor(as.character(pn_long$var_j), levels = TARGET_VARS)

  # Modality bar separators (between modality blocks)
  blk_breaks <- cumsum(table(factor(modality_tbl$modality_code,
                                    levels = MODALITY_LEVELS)))
  blk_breaks <- blk_breaks[-length(blk_breaks)] + 0.5

  p_pn <- ggplot(pn_long, aes(x = var_i, y = var_j, fill = n)) +
    geom_tile(color = "white", linewidth = 0.05) +
    geom_text(aes(label = scales::comma(n)),
              size = 1.4, color = ifelse(pn_long$n < 200, "red3", "grey25")) +
    geom_vline(xintercept = blk_breaks, color = "black", linewidth = 0.4) +
    geom_hline(yintercept = blk_breaks, color = "black", linewidth = 0.4) +
    scale_fill_viridis_c(option = "mako", direction = -1,
                         name = "pairwise N",
                         labels = scales::comma) +
    scale_y_discrete(limits = rev(TARGET_VARS)) +
    coord_fixed() +
    labs(
      title = "Pairwise complete-case N across 24 smoking variables",
      subtitle = paste0(
        "Cells with N < 200 (red text) cannot support stable joint regression. ",
        "Joint complete-case N (all 24 vars) = ", joint_N_all, "."),
      x = NULL, y = NULL,
      caption = "Black grid lines: modality-block boundaries (Q | U-PAH | U-MM | B-VOC)."
    ) +
    theme_compact() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# The defense of mutual-adjustment infeasibility rests on
# pairwise_complete_case_N.csv directly (joint complete-case N across all 24
# vars = 0; see the pairwise-N panel of the correlation heatmap PDF).
saturated_24_N <- sum(stats::complete.cases(x_log))
cat("\n[note] saturated 24-var joint complete-case N =", saturated_24_N, "\n")

# ---- Spearman correlation heatmap (24x24) | Pairwise N heatmap ----
# Side-by-side composition: Spearman ComplexHeatmap (left) + Pairwise N ggplot
# (right) rendered into one PDF.
{
  spear <- suppressWarnings(stats::cor(x_log, use = "pairwise.complete.obs",
                                       method = "spearman"))
  # Display strings (2 decimals; "" for NA)
  disp_s <- matrix("", nrow = nrow(spear), ncol = ncol(spear),
                   dimnames = dimnames(spear))
  for (i in seq_len(nrow(spear))) for (j in seq_len(ncol(spear))) {
    v <- spear[i, j]
    disp_s[i, j] <- if (is.na(v)) "" else sprintf("%.2f", v)
  }

  # Modality column annotation
  mod_codes <- modality_tbl$modality_code[match(TARGET_VARS,
                                                modality_tbl$var_name)]
  col_ann <- ComplexHeatmap::HeatmapAnnotation(
    Modality = mod_codes,
    col = list(Modality = MODALITY_COLOR),
    annotation_name_gp = gpar(fontsize = THEME_BASE_SIZE,
                              fontfamily = THEME_BASE_FAMILY),
    annotation_legend_param = list(
      Modality = list(title_gp = gpar(fontsize = THEME_BASE_SIZE,
                                      fontfamily = THEME_BASE_FAMILY),
                      labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                                       fontfamily = THEME_BASE_FAMILY))),
    show_legend = TRUE,
    simple_anno_size = grid::unit(2, "mm")
  )
  pretty_lab <- modality_tbl$var_description
  spear_safe <- spear; spear_safe[is.na(spear_safe)] <- 0
  col_fun <- circlize::colorRamp2(c(-1, 0, 1),
                                  c("#3B4992", "white", "#A20056"))

  ht_s <- ComplexHeatmap::Heatmap(
    spear_safe,
    name = "Spearman rho",
    col = col_fun,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = TARGET_VARS, column_order = TARGET_VARS,
    row_labels = pretty_lab, column_labels = TARGET_VARS,
    top_annotation = col_ann,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_side = "left",
    column_names_rot = 90,
    row_names_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                        fontfamily = THEME_BASE_FAMILY),
    column_names_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                           fontfamily = THEME_BASE_FAMILY),
    column_title = paste0(
      "Spearman rho (pairwise.complete.obs on log-transformed values; ",
      "median pairwise N = ", median(pairwise_N[upper.tri(pairwise_N)]),
      ", min = ", min(pairwise_N[upper.tri(pairwise_N)]),
      ", max = ", max(pairwise_N[upper.tri(pairwise_N)]), ")"),
    column_title_gp = gpar(fontsize = THEME_BASE_SIZE + 1, fontface = "bold",
                           fontfamily = THEME_BASE_FAMILY),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(disp_s[i, j], x, y,
                gp = gpar(fontsize = 4.0,
                          fontfamily = THEME_BASE_FAMILY))
    },
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = THEME_BASE_SIZE,
                      fontfamily = THEME_BASE_FAMILY),
      labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                       fontfamily = THEME_BASE_FAMILY),
      legend_height = grid::unit(2, "cm")
    ),
    width  = grid::unit(120, "mm"),
    height = grid::unit(120, "mm"),
    border = TRUE
  )

  # Combine: Spearman ComplexHeatmap (grabExpr -> grob) | Pairwise N ggplot.
  # Use a null pdf device to capture the heatmap grob without R opening the
  # default Rplots.pdf device on a headless Linux node.
  grDevices::pdf(NULL)
  ht_grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(ht_s, merge_legend = TRUE,
                         heatmap_legend_side = "right",
                         annotation_legend_side = "right"),
    warn = FALSE
  )
  grDevices::dev.off()
  pn_grob <- ggplot2::ggplotGrob(p_pn)
  side_by_side <- gridExtra::arrangeGrob(ht_grob, pn_grob, ncol = 2,
                                          widths = grid::unit(c(1, 1),
                                                              "null"))

  fp <- open_pdf("b2_correlation_heatmap_24vars.pdf",
                 14.0 * 0.85, 7.5 * 0.85)
  grid::grid.draw(side_by_side)
  dev.off()
  cat("[plot] b2_correlation_heatmap_24vars.pdf saved (side-by-side) -> ",
      fp, "\n")
}

cat("\n=== SECTION 3 DONE ===\n\n")

# =============================================================================
# SECTION 4 - Per-modality PCAs
#   d_pca_modality_marginal_PCs.pdf
#
# Four independent PCAs, one per modality block (Q is skipped - single var).
# Each panel uses the standard prcomp route on the within-modality complete-
# case subset (NHANES sub-samples DO overlap enough within a modality to
# give usable N), not the correlation-matrix workaround used for the 24-var
# global summary (which is undefined because joint complete-case N = 0).
# =============================================================================
cat("=== SECTION 4 - Per-modality PCAs ===\n")

{
  panels <- list()
  pca_meta <- list()  # for the figure-level caption
  for (b in MODALITY_LEVELS) {
    vars_b <- modality_tbl$var_name[modality_tbl$modality_code == b]
    if (length(vars_b) < 2) {
      cat("[pca][", b, "] skipped (single-var modality)\n", sep = "")
      next
    }
    x_b <- x_log[, vars_b, drop = FALSE]
    cc <- stats::complete.cases(x_b)
    n_cc <- sum(cc)
    if (n_cc < 50) {
      cat("[pca][", b, "] skipped (N=", n_cc, " < 50)\n", sep = "")
      next
    }
    x_bcc <- x_b[cc, , drop = FALSE]
    # Drop zero-variance LOD-pinned vars
    sd_cols <- apply(x_bcc, 2, stats::sd, na.rm = TRUE)
    keep    <- is.finite(sd_cols) & sd_cols > 0
    dropped <- vars_b[!keep]
    if (length(dropped)) {
      cat("[pca][", b, "] dropped zero-variance vars: ",
          paste(dropped, collapse = ", "), "\n", sep = "")
    }
    x_bcc <- x_bcc[, keep, drop = FALSE]
    if (ncol(x_bcc) < 2) next

    x_bz <- scale(x_bcc)
    set.seed(42)
    pcr  <- stats::prcomp(x_bz, center = FALSE, scale. = FALSE)
    vex  <- (pcr$sdev^2) / sum(pcr$sdev^2)
    loadings_b <- pcr$rotation[, 1:2] * matrix(pcr$sdev[1:2],
                                               nrow = nrow(pcr$rotation),
                                               ncol = 2, byrow = TRUE)
    load_b <- data.frame(
      var_name = rownames(loadings_b),
      PC1      = loadings_b[, 1],
      PC2      = loadings_b[, 2],
      stringsAsFactors = FALSE
    )
    lim_b <- max(abs(c(load_b$PC1, load_b$PC2))) * 1.25
    pca_meta[[b]] <- list(
      n_cc    = n_cc,
      n_vars  = ncol(x_bcc),
      vars    = colnames(x_bcc),
      dropped = dropped,
      vex     = vex
    )
    p <- ggplot(load_b, aes(x = PC1, y = PC2)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
      geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
      geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                   color = MODALITY_COLOR[b], linewidth = 0.35,
                   arrow = grid::arrow(length = grid::unit(1, "mm"),
                                       type = "closed")) +
      geom_point(color = MODALITY_COLOR[b], size = 1.2) +
      ggrepel::geom_text_repel(aes(label = var_name),
                               size = 1.8, max.overlaps = 30,
                               segment.size = 0.2,
                               min.segment.length = 0) +
      coord_equal(xlim = c(-lim_b, lim_b), ylim = c(-lim_b, lim_b)) +
      labs(title = paste0(MODALITY_LABEL[b], "  (", b, ")"),
           subtitle = paste0(
             "N = ", scales::comma(n_cc),
             ";  vars used = ", ncol(x_bcc),
             if (length(dropped)) paste0(
               " (dropped ", length(dropped),
               " zero-variance after subset)") else "",
             ";  PC1 = ", scales::percent(vex[1], 0.1),
             ";  PC2 = ", scales::percent(vex[2], 0.1)),
           x = paste0("PC1 (", scales::percent(vex[1], 0.1), ")"),
           y = paste0("PC2 (", scales::percent(vex[2], 0.1), ")")) +
      theme_compact()
    panels[[b]] <- p
  }

  if (length(panels) >= 1) {
    # Build the methods annotation text from the actual run values
    methods_lines <- c(
      "Four independent PCAs, one per modality block.",
      paste0(
        "Q (single var) skipped — PCA on a single column is undefined. ",
        "Each remaining panel uses the standard prcomp route on the within-",
        "modality complete-case subset (NHANES sub-samples DO overlap enough ",
        "within a modality to give usable N), NOT the correlation-matrix ",
        "workaround used for the cross-modality summary."),
      paste0(
        "U-PAH: N = ", scales::comma(pca_meta$`U-PAH`$n_cc),
        " participants with all ", pca_meta$`U-PAH`$n_vars,
        " PAH metabolites measured; prcomp on standardized ",
        scales::comma(pca_meta$`U-PAH`$n_cc),
        " × ", pca_meta$`U-PAH`$n_vars, " matrix."),
      paste0(
        "U-MM: N = ", scales::comma(pca_meta$`U-MM`$n_cc),
        " with all ", pca_meta$`U-MM`$n_vars,
        " urinary metals/mercapturates; prcomp on standardized ",
        scales::comma(pca_meta$`U-MM`$n_cc),
        " × ", pca_meta$`U-MM`$n_vars, " matrix."),
      paste0(
        "B-VOC: N = ", scales::comma(pca_meta$`B-VOC`$n_cc),
        " with all 16 blood VOCs / combustion analytes ",
        if (length(pca_meta$`B-VOC`$dropped))
          paste0("(", length(pca_meta$`B-VOC`$dropped),
                 " zero-variance LOD-pinned vars dropped from the subset: ",
                 paste(pca_meta$`B-VOC`$dropped, collapse = ", "), "); ")
        else "; ",
        "prcomp on standardized ", scales::comma(pca_meta$`B-VOC`$n_cc),
        " × ", pca_meta$`B-VOC`$n_vars, " matrix."),
      paste0(
        "Each panel plots PC1 vs PC2 variable loadings ",
        "(eigenvector × sqrt(eigenvalue)). ",
        "Joint complete-case N across ALL 24 vars = ",
        sum(stats::complete.cases(x_log)),
        " — saturated 24-var PCA is infeasible.")
    )
    methods_text <- paste(methods_lines, collapse = "\n")

    # 1x3 grid: 3 PCAs side-by-side at the same PDF width
    p_marg <- patchwork::wrap_plots(panels, ncol = 3) +
      patchwork::plot_annotation(
        title    = "Per-modality PCAs",
        subtitle = "PC1/PC2 variable loadings, one panel per modality block",
        caption  = methods_text,
        theme    = theme_compact()
      ) &
      ggplot2::theme(plot.caption = element_text(size = THEME_BASE_SIZE - 1,
                                                 color = "grey25",
                                                 hjust = 0,
                                                 lineheight = 1.1))
    fp <- open_pdf("d_pca_modality_marginal_PCs.pdf",
                   8.0 * 0.75, 5.5 * 0.75)
    print(p_marg)
    dev.off()
    cat("[plot] d_pca_modality_marginal_PCs.pdf saved -> ", fp, "\n")
  }
}

cat("\n=== SECTION 4 DONE ===\n\n")

# =============================================================================
# SECTION 5 - extract canonical CLR association rows for the 24 vars
#   - per-variable significant-taxa bars PDF
# =============================================================================
cat("=== SECTION 5 - canonical CLR associations ===\n")

# Read the canonical tidied table (~525 MB, ~9.5M rows). We filter immediately
# to the rows of interest (independent_var %in% TARGET_VARS AND term %in%
# {non-SMQ + SMQ dummies}) and discard everything else from memory.
cat("[load] reading tidied (this may take ~30s)...\n")
tidied_all <- readRDS(INPUTS$tidied_rds)
cat("[load] tidied dim = ",
    paste(dim(tidied_all), collapse = " x "), "\n", sep = "")

NON_SMQ_VARS <- setdiff(TARGET_VARS, "SMQ_current_ever_never")
SMQ_DUMMY_TERMS <- c("SMQ_current_ever_never1", "SMQ_current_ever_never2")
EXPOSURE_TERMS <- c(NON_SMQ_VARS, SMQ_DUMMY_TERMS)

tidied_smoke <- tidied_all %>%
  dplyr::filter(.data$independent_var %in% TARGET_VARS,
                .data$term %in% EXPOSURE_TERMS) %>%
  dplyr::mutate(
    # Strip the "_relative" suffix the phyloseq naming uses
    taxon = sub("_relative$", "", .data$phenotype),
    # Collapse SMQ dummies into one parent var_name
    var_name = ifelse(grepl("^SMQ_current_ever_never", .data$term),
                      "SMQ_current_ever_never", .data$independent_var)
  ) %>%
  dplyr::select(var_name, term, taxon, estimate, std.error, statistic,
                p.value, p.value.fdr, q.value, n_obs, formula_used,
                cycle_mode)
rm(tidied_all); invisible(gc())
cat("[filter] rows after smoking filter = ", nrow(tidied_smoke), "\n", sep = "")

# Per-taxon SMQ collapse: keep the dummy with the smaller q (i.e. stronger
# evidence). This becomes the "SMQ" row used in downstream plots.
tidied_smq <- tidied_smoke %>%
  dplyr::filter(.data$var_name == "SMQ_current_ever_never") %>%
  dplyr::group_by(.data$taxon) %>%
  dplyr::slice(which.min(.data$q.value)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(term = "SMQ_current_ever_never")
tidied_non_smq <- tidied_smoke %>%
  dplyr::filter(.data$var_name != "SMQ_current_ever_never")
tidied_smoke_collapsed <- dplyr::bind_rows(tidied_non_smq, tidied_smq)
cat("[collapse] after SMQ-dummy collapse rows = ",
    nrow(tidied_smoke_collapsed), "  unique vars = ",
    length(unique(tidied_smoke_collapsed$var_name)),
    "  unique taxa = ",
    length(unique(tidied_smoke_collapsed$taxon)), "\n", sep = "")

# Stash the filtered association subset for downstream sections (small, ~MBs).
saveRDS(tidied_smoke_collapsed,
        file.path(INPUTS_DIR, "tidied_smoking_subset.rds"))
cat("[input] tidied_smoking_subset.rds written for downstream reuse\n")

# Per-variable significant-taxa counts using a joint rule: both
# p.value.fdr <= 0.01 (BH-corrected) AND q.value <= 0.01 (Storey).
sig_counts <- tidied_smoke_collapsed %>%
  dplyr::group_by(.data$var_name) %>%
  dplyr::summarise(
    n_taxa             = dplyr::n(),
    n_joint_fdr_q_0_01 = sum(.data$p.value.fdr <= 0.01 &
                             .data$q.value     <= 0.01, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(modality_tbl[, c("var_name", "modality_code",
                                    "modality_label", "plot_order")],
                   by = "var_name") %>%
  dplyr::arrange(.data$plot_order)
cat("[note] sig_counts kept in memory for the per-variable bar plot.\n")

# Per-variable significant-taxa bar plot
{
  sc_p <- sig_counts %>%
    dplyr::mutate(var_name = factor(.data$var_name,
                                    levels = sig_counts$var_name),
                  mod_color = MODALITY_COLOR[.data$modality_code])

  p_sig <- ggplot(sc_p,
                  aes(x = var_name, y = n_joint_fdr_q_0_01,
                      fill = .data$modality_code)) +
    geom_col(color = "grey25", linewidth = 0.1, width = 0.7) +
    geom_text(aes(label = n_joint_fdr_q_0_01),
              vjust = -0.3, size = 1.6, color = "grey25") +
    scale_fill_manual(values = MODALITY_COLOR, name = "Modality") +
    labs(title = "Per-variable significant-taxa counts (CLR associations, 1349 taxa)",
         subtitle = paste0(
           "Joint rule: BH-FDR ≤ 0.01 AND Storey q ≤ 0.01. ",
           "Vars ordered by modality block."),
         x = NULL, y = "Number of taxa",
         caption = paste0(
           "Storey q-value applied scheme-wise upstream by ",
           "results/3_exWAS_out/result_clr/.")) +
    theme_compact() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     face = "bold"))

  # Modality-block separator lines
  block_x <- which(diff(as.integer(factor(
    modality_tbl$modality_code[match(levels(sc_p$var_name),
                                     modality_tbl$var_name)],
    levels = MODALITY_LEVELS))) != 0) + 0.5
  p_sig <- p_sig +
    geom_vline(xintercept = block_x, color = "black", linewidth = 0.3,
               linetype = "dashed")

  fp <- open_pdf("b1_per_variable_sig_taxa_bars.pdf",
                 8.5 * 0.75 * 0.70, 4.5 * 0.75 * 0.70)
  print(p_sig)
  dev.off()
  cat("[plot] b1_per_variable_sig_taxa_bars.pdf saved -> ", fp, "\n")
}

# Read the CLR phyloseq once to stash tax_tab + genus column used by the
# taxon-by-variable heatmap row labels.
{
  clr_ps <- readRDS(INPUTS$clr_phyloseq_rds)
  tax_tab <- as(phyloseq::tax_table(clr_ps), "matrix")
  genus_col <- if ("Genus" %in% colnames(tax_tab)) "Genus" else
               if ("genus" %in% colnames(tax_tab)) "genus" else NA_character_
  saveRDS(list(tax_tab = tax_tab, genus_col = genus_col,
               sample_names = sample_names(clr_ps),
               nsamples = nsamples(clr_ps),
               ntaxa = ntaxa(clr_ps)),
          file.path(INPUTS_DIR, "clr_tax_metadata.rds"))
  cat("[note] clr_tax_metadata.rds stashed for taxon-by-variable heatmap labels.\n")
  rm(clr_ps); invisible(gc())
}

cat("\n=== SECTION 5 DONE ===\n\n")

# =============================================================================
# SECTION 6 - signature analyses
#   e_taxon_by_variable_significance_heatmap.pdf
#   a_signature_similarity_heatmap.pdf  (precomputed; size = FDR, color = rho)
# =============================================================================
cat("=== SECTION 6 - signature analyses ===\n")

# Build per-var FDR-significant taxa lists (q <= 0.01) - used by both heatmaps
sig_list <- split(tidied_smoke_collapsed$taxon[
  tidied_smoke_collapsed$q.value <= 0.01],
  factor(tidied_smoke_collapsed$var_name[
    tidied_smoke_collapsed$q.value <= 0.01],
    levels = modality_tbl$var_name))
sig_list <- lapply(sig_list, function(z) unique(z[!is.na(z)]))
cat("[sig_list] sizes per var:\n")
print(sapply(sig_list, length))

# ---- taxon x variable significance heatmap ----
{
  union_taxa <- unique(unlist(sig_list))
  cat("[heatmap] union sig taxa (q<=0.01): ",
      length(union_taxa), "\n", sep = "")
  if (length(union_taxa) >= 2) {
    # Cell value = sign(beta) * log(FDR), where FDR is the Storey q-value
    # (labeled "FDR" in the figure for clarity). log() is natural log.
    # Values are NEGATIVE wherever the association is significant
    # (FDR < 1 -> log(FDR) < 0), so the color ramp is read as:
    #   strongly NEGATIVE cell -> highly significant POSITIVE association
    #   strongly POSITIVE cell -> highly significant NEGATIVE association
    mat <- matrix(NA_real_, nrow = length(union_taxa), ncol = 24,
                  dimnames = list(union_taxa, modality_tbl$var_name))
    for (v in modality_tbl$var_name) {
      sub <- tidied_smoke_collapsed[
        tidied_smoke_collapsed$var_name == v &
        tidied_smoke_collapsed$taxon %in% union_taxa &
        !is.na(tidied_smoke_collapsed$q.value) &
        tidied_smoke_collapsed$q.value <= 0.01, ]
      if (!nrow(sub)) next
      idx <- match(sub$taxon, union_taxa)
      mat[idx, v] <- sign(sub$estimate) *
        log(pmax(sub$q.value, 1e-300))
    }
    # Cap absolute scale at the 99th percentile for color stability
    cap <- stats::quantile(abs(mat), 0.99, na.rm = TRUE)
    if (!is.finite(cap) || cap <= 0) cap <- 5
    # NA -> 0 so clustering works (0 = white in the diverging ramp = not sig)
    mat[is.na(mat)] <- 0

    # Reuse the lead-taxa genus labeler (clr_ps tax_tab in metadata RDS)
    meta <- readRDS(file.path(INPUTS_DIR, "clr_tax_metadata.rds"))
    tax_tab <- meta$tax_tab; genus_col <- meta$genus_col
    row_labs <- vapply(union_taxa, function(tid) {
      if (is.na(genus_col)) return(tid)
      i <- match(tid, rownames(tax_tab))
      if (is.na(i)) return(tid)
      g <- trimws(as.character(tax_tab[i, genus_col]))
      if (!nzchar(g) || is.na(g) || g %in% c("NA", "Unassigned")) return(tid)
      g <- sub("^g__", "", g); g <- sub("^G__", "", g)
      paste0(g, " (", tid, ")")
    }, character(1))

    mod_codes <- modality_tbl$modality_code
    col_ann <- ComplexHeatmap::HeatmapAnnotation(
      Modality = mod_codes,
      col = list(Modality = MODALITY_COLOR),
      simple_anno_size = grid::unit(2, "mm"),
      annotation_name_gp = gpar(fontsize = THEME_BASE_SIZE,
                                fontfamily = THEME_BASE_FAMILY,
                                fontface = "bold"),
      annotation_legend_param = list(
        Modality = list(title_gp = gpar(fontsize = THEME_BASE_SIZE,
                                        fontfamily = THEME_BASE_FAMILY),
                        labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                                         fontfamily = THEME_BASE_FAMILY)))
    )
    col_fun <- circlize::colorRamp2(
      c(-cap, -cap/2, 0, cap/2, cap),
      c("#3B4992", "#9FB4D8", "white", "#E0A9B6", "#A20056"))

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "signed\nlog(FDR)",
      col = col_fun,
      cluster_rows = TRUE, cluster_columns = FALSE,
      column_order = modality_tbl$var_name,
      row_dend_width = grid::unit(8, "mm"),
      show_row_dend = TRUE,
      row_labels = row_labs,
      column_labels = modality_tbl$var_description,
      top_annotation = col_ann,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = THEME_BASE_SIZE - 3,
                          fontfamily = THEME_BASE_FAMILY),
      column_names_gp = gpar(fontsize = THEME_BASE_SIZE,
                             fontfamily = THEME_BASE_FAMILY,
                             fontface = "bold"),
      column_names_rot = 90,
      column_title = paste0(
        "Taxon x variable signed-significance map (FDR ≤ 0.01)\n",
        "Cell = sign(beta) * log(FDR)  [natural log; FDR = Storey q].\n",
        "Rows: union of FDR-sig taxa across 24 vars (n=", length(union_taxa),
        "). Columns: 24 vars, modality-blocked."),
      column_title_gp = gpar(fontsize = THEME_BASE_SIZE + 1,
                             fontface = "bold",
                             fontfamily = THEME_BASE_FAMILY),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = THEME_BASE_SIZE,
                        fontfamily = THEME_BASE_FAMILY),
        labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                         fontfamily = THEME_BASE_FAMILY),
        legend_height = grid::unit(2, "cm")),
      width  = grid::unit(110, "mm"),
      height = grid::unit(220 * 1.9 * 1.9 * 1.25, "mm"),
      border = TRUE
    )

    fp <- open_pdf("e_taxon_by_variable_significance_heatmap.pdf",
                   9.0 * 0.85, 12.0 * 1.7 * 1.9 * 1.9 * 1.25)
    ComplexHeatmap::draw(ht, merge_legend = TRUE)
    dev.off()
    cat("[plot] e_taxon_by_variable_significance_heatmap.pdf saved -> ",
        fp, "\n")
  } else {
    cat("[plot] taxon-by-variable heatmap skipped - union sig set too small\n")
  }
}

# ---- Signature similarity (sized cells, color = rho, area = FDR) -----------
# Uses the precomputed canonical host-host weighted signature correlation
# matrix produced by scripts/7_microbial_signature_heatmap/. We read it,
# subset to the 24 smoking vars, and render with the "size encodes FDR"
# convention (inspired by 2_microbial_signature_heatmap_sized.R):
#   color = rho   (diverging blue-white-red)
#   area  = FDR tier  (FDR<0.001 -> 1.00; <0.01 -> 0.80; <0.05 -> 0.55;
#                      >=0.05 -> 0.25)
# We never modify the source file; this is a read-only consumer.
{
  cat("[similarity] reading precomputed signature correlation matrix...\n")
  sig_corr <- readRDS(INPUTS$sig_corr_rds)
  miss_in_corr <- setdiff(TARGET_VARS, rownames(sig_corr$correlation_matrix))
  if (length(miss_in_corr)) {
    stop("signature-similarity plot: precomputed matrix missing vars: ",
         paste(miss_in_corr, collapse = ", "))
  }
  rho_mat <- sig_corr$correlation_matrix[TARGET_VARS, TARGET_VARS]
  fdr_mat <- sig_corr$fdr_matrix[TARGET_VARS, TARGET_VARS]
  diag(fdr_mat) <- 0   # diagonals should be perfect-correlation

  # FDR tier -> cell-area fraction (inspired by 7_microbial_signature)
  fdr_size <- function(q) {
    s <- rep(0.25, length(q))                       # >= 0.05  (smallest)
    s[!is.na(q) & q < 0.05]  <- 0.55                # [0.01, 0.05)
    s[!is.na(q) & q < 0.01]  <- 0.80                # [0.001, 0.01)
    s[!is.na(q) & q < 0.001] <- 1.00                # < 0.001  (largest)
    s
  }
  size_mat <- matrix(fdr_size(as.vector(fdr_mat)),
                     nrow = nrow(fdr_mat), ncol = ncol(fdr_mat),
                     dimnames = dimnames(fdr_mat))

  rho_safe <- rho_mat; rho_safe[is.na(rho_safe)] <- 0
  col_fun_rho <- circlize::colorRamp2(
    seq(-1, 1, length.out = 11),
    c("#2B5B8A", "#4071A0", "#5789B6", "#6FA3CB", "#A2BCCF",
      "#D8D4C9", "#F0AC72", "#EF8530", "#DA6524", "#BE4E21", "#9E3D21"))

  # Modality annotation (top + side)
  col_ann2 <- ComplexHeatmap::HeatmapAnnotation(
    Modality = modality_tbl$modality_code,
    col = list(Modality = MODALITY_COLOR),
    simple_anno_size = grid::unit(2, "mm"),
    annotation_name_gp = gpar(fontsize = THEME_BASE_SIZE,
                              fontfamily = THEME_BASE_FAMILY,
                              fontface = "bold"),
    annotation_legend_param = list(
      Modality = list(title_gp = gpar(fontsize = THEME_BASE_SIZE,
                                      fontfamily = THEME_BASE_FAMILY),
                      labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                                       fontfamily = THEME_BASE_FAMILY)))
  )

  # Heatmap with a white background and a cell_fun that draws sized colored
  # tiles. We pass a "white" color matrix (zeros) as the heatmap data so the
  # background remains blank; sizing/coloring happens entirely inside cell_fun.
  white_mat <- rho_safe; white_mat[] <- 0
  ht_sim <- ComplexHeatmap::Heatmap(
    white_mat,
    name = "Weighted\nrho",
    col = circlize::colorRamp2(c(-1, 0, 1), c("white", "white", "white")),
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_order = TARGET_VARS, column_order = TARGET_VARS,
    row_labels = modality_tbl$var_description,
    column_labels = TARGET_VARS,
    top_annotation = col_ann2,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                        fontfamily = THEME_BASE_FAMILY),
    column_names_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                           fontfamily = THEME_BASE_FAMILY,
                           fontface = "bold"),
    column_names_rot = 90,
    column_title = paste0(
      "Signature similarity (precomputed weighted rho)\n",
      "Color = rho   |   Cell area = FDR tier   ",
      "(<0.001 largest; <0.01; <0.05; ≥0.05 smallest)"),
    column_title_gp = gpar(fontsize = THEME_BASE_SIZE + 1,
                           fontface = "bold",
                           fontfamily = THEME_BASE_FAMILY),
    rect_gp = gpar(col = "white", lwd = 0.4),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (i == j) {
        grid::grid.rect(x, y, w, h,
                        gp = gpar(fill = "grey90", col = "white", lwd = 0.4))
        return(invisible())
      }
      rho_v <- rho_mat[i, j]
      sz    <- size_mat[i, j]
      if (is.na(rho_v)) return(invisible())
      grid::grid.rect(
        x, y, width = w * sz, height = h * sz,
        gp = gpar(fill = col_fun_rho(rho_v), col = "white", lwd = 0.4)
      )
    },
    width  = grid::unit(120, "mm"),
    height = grid::unit(120, "mm"),
    border = TRUE
  )

  # Build a 2-component legend: rho color ramp + FDR-tier size legend.
  lgd_rho <- ComplexHeatmap::Legend(
    col_fun = col_fun_rho,
    title   = "Weighted rho",
    at      = c(-1, -0.5, 0, 0.5, 1),
    title_gp = gpar(fontsize = THEME_BASE_SIZE,
                    fontfamily = THEME_BASE_FAMILY,
                    fontface = "bold"),
    labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                     fontfamily = THEME_BASE_FAMILY),
    legend_height = grid::unit(3, "cm")
  )

  # FDR-tier "size" legend: a column of squares of decreasing size
  fdr_tiers <- c("FDR < 0.001" = 1.0,
                 "FDR < 0.01"  = 0.8,
                 "FDR < 0.05"  = 0.55,
                 "FDR ≥ 0.05"  = 0.25)
  lgd_fdr <- ComplexHeatmap::Legend(
    title = "FDR tier (area)",
    labels = names(fdr_tiers),
    type   = "points",
    pch    = 22,
    size   = grid::unit(fdr_tiers * 5.0, "mm"),
    legend_gp = gpar(fill = "grey45", col = "grey25"),
    background = "white",
    title_gp = gpar(fontsize = THEME_BASE_SIZE,
                    fontfamily = THEME_BASE_FAMILY,
                    fontface = "bold"),
    labels_gp = gpar(fontsize = THEME_BASE_SIZE - 1,
                     fontfamily = THEME_BASE_FAMILY)
  )

  fp <- open_pdf("a_signature_similarity_heatmap.pdf",
                 10.5 * 0.85 * 0.85 * 0.80,
                 9.5  * 0.85 * 0.85 * 0.80)
  ComplexHeatmap::draw(ht_sim,
                       annotation_legend_list = list(lgd_rho, lgd_fdr),
                       merge_legend = TRUE,
                       annotation_legend_side = "right",
                       heatmap_legend_side = "right")
  dev.off()
  cat("[plot] a_signature_similarity_heatmap.pdf saved -> ", fp, "\n")
}

cat("\n=== SECTION 6 DONE ===\n\n")

# =============================================================================
# SECTION 7 - modality-level summary
#   c_modality_modality_jaccard.pdf
#   modality_stratified_jaccard.csv
# =============================================================================
cat("=== SECTION 7 - modality summary ===\n")

# Modality-level sig taxa sets (union of per-var sig taxa within modality)
mod_sig <- lapply(MODALITY_LEVELS, function(b) {
  vs <- modality_tbl$var_name[modality_tbl$modality_code == b]
  unique(unlist(sig_list[vs]))
})
names(mod_sig) <- MODALITY_LEVELS
cat("[modality] sig taxa counts per modality:\n")
print(sapply(mod_sig, length))

# 4x4 Jaccard
jac_mat <- matrix(NA_real_, 4, 4,
                  dimnames = list(MODALITY_LEVELS, MODALITY_LEVELS))
for (i in seq_along(MODALITY_LEVELS)) {
  for (j in seq_along(MODALITY_LEVELS)) {
    a <- mod_sig[[i]]; b <- mod_sig[[j]]
    if (length(a) == 0 && length(b) == 0) { jac_mat[i, j] <- NA; next }
    jac_mat[i, j] <- length(intersect(a, b)) / max(length(union(a, b)), 1)
  }
}
write.csv(jac_mat,
          file.path(TABLES_DIR, "modality_stratified_jaccard.csv"))
cat("[table] modality_stratified_jaccard.csv (4x4 Jaccard):\n")
print(round(jac_mat, 3))

# Modality-modality Jaccard similarity plot (single panel).
{
  jac_long <- as.data.frame(as.table(jac_mat))
  names(jac_long) <- c("m_i", "m_j", "jaccard")
  jac_long$m_i <- factor(as.character(jac_long$m_i), levels = MODALITY_LEVELS)
  jac_long$m_j <- factor(as.character(jac_long$m_j), levels = MODALITY_LEVELS)
  p_b <- ggplot(jac_long, aes(x = m_i, y = m_j, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(is.na(.data$jaccard), "",
                                 sprintf("%.2f", .data$jaccard))),
              size = 2.4, color = "grey15") +
    scale_fill_viridis_c(option = "mako", direction = -1,
                         name = "Jaccard",
                         limits = c(0, 1)) +
    scale_y_discrete(limits = rev(MODALITY_LEVELS)) +
    coord_equal() +
    labs(title = "Modality-modality Jaccard similarity",
         subtitle = paste0(
           "Jaccard overlap among the per-modality FDR-sig taxa sets ",
           "(q ≤ 0.01).  Diagonals = 1 by construction."),
         x = NULL, y = NULL,
         caption = "Signatures differ substantially across modalities.") +
    theme_compact() +
    theme(axis.text.x = element_text(angle = 0))

  fp <- open_pdf("c_modality_modality_jaccard.pdf",
                 4.5 * 0.40, 4.0 * 0.40)
  print(p_b)
  dev.off()
  cat("[plot] c_modality_modality_jaccard.pdf saved -> ", fp, "\n")
}

cat("\n=== SECTION 7 DONE ===\n\n")

# =============================================================================
# SECTION 8 - manifest + assertions
# =============================================================================
cat("=== SECTION 8 - manifest + assertions ===\n")

all_files <- list.files(viz_out_path, recursive = TRUE, full.names = TRUE)
manifest <- data.frame(
  path  = sub(paste0("^", viz_out_path, "/"), "", all_files),
  bytes = file.info(all_files)$size,
  md5   = unname(tools::md5sum(all_files)),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(viz_out_path, "manifest.csv"), row.names = FALSE)
cat("[manifest] wrote manifest.csv (", nrow(manifest), " entries)\n", sep = "")

# Assertion suite
expected_pdfs <- c(
  "a_signature_similarity_heatmap.pdf",
  "b1_per_variable_sig_taxa_bars.pdf",
  "b2_correlation_heatmap_24vars.pdf",
  "c_modality_modality_jaccard.pdf",
  "d_pca_modality_marginal_PCs.pdf",
  "e_taxon_by_variable_significance_heatmap.pdf"
)
expected_tables <- c(
  "modality_stratified_jaccard.csv",
  "pairwise_complete_case_N.csv"
)
assertions <- list(
  `all expected PDFs exist` =
      all(file.exists(file.path(PLOTS_DIR, expected_pdfs))),
  `exactly the published tables exist` =
      setequal(list.files(TABLES_DIR, pattern = "\\.csv$"),
               expected_tables),
  `modality mapping has 24 rows` =
      nrow(modality_tbl) == 24,
  `input provenance exists` =
      file.exists(file.path(INPUTS_DIR, "input_provenance.txt")),
  `disallowed CLR file not in inputs` =
      !any(grepl("ubiome_relative_clr_updated\\.rds$", prov$path)),
  `4_association tree not in inputs` =
      !any(grepl("/4_association_phyloseq_analyses_out/", prov$path,
                 fixed = TRUE)),
  `precomputed signature corr matrix was read` =
      any(grepl("sig_corr_rds", prov$role))
)
for (nm in names(assertions)) {
  if (!isTRUE(assertions[[nm]])) {
    stop("ASSERTION FAILED: ", nm)
  } else {
    cat("[ok] ", nm, "\n", sep = "")
  }
}

cat("\n=== run complete: ", viz_out_path, " ===\n", sep = "")
