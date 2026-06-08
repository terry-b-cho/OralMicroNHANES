#!/usr/bin/env Rscript

## -----------------------------------------------------------------------------
## Terry Plot PDF Generation (Host / OTU / *_indices)
## -----------------------------------------------------------------------------
## Reads intermediates produced by
## 1_otu_host_and_indicies_intermediate_processing_additional.R and writes all
## PDFs under
## results/analyses_results/4_otu_host_plot_out_additional/figures_out/{compact,comprehensive}/...
##
## Environment: R >= 4.5 with data.table, dplyr, tidyr, tibble, purrr, stringr,
## forcats, readr, ggplot2, ggrepel, egg, grid, gridExtra, extrafont.
## Conda spec: envs/nhanes-analysis_for_reviewers.yml.
## -----------------------------------------------------------------------------

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

set.seed(42)

resolve_base_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
  if (length(script_path) == 0) {
    return(normalizePath("."))
  }
  normalizePath(file.path(dirname(script_path), "..", ".."))
}

`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

base_path <- PROJECT_ROOT

## Output tree before library() so paths exist even if ggplot/phyloseq stack fails
output_base_path <- file.path(base_path, "results/analyses_results/4_otu_host_plot_out_additional")
viz_out_path <- file.path(output_base_path, "figures_out")
intermediate_files_path <- file.path(output_base_path, "intermediate")
supplementary_files_path <- file.path(output_base_path, "supplementary")
dir.create(supplementary_files_path, recursive = TRUE, showWarnings = FALSE)
dir.create(viz_out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_files_path, recursive = TRUE, showWarnings = FALSE)

## Terry plots: separate compact (main) vs comprehensive (supplement) export trees
terry_modes <- c("compact", "comprehensive")
transformation_dirs <- c("none", "hellinger", "clr")
for (m in terry_modes) {
  for (dir_name in transformation_dirs) {
    dir.create(file.path(viz_out_path, m, dir_name), recursive = TRUE, showWarnings = FALSE)
    dir.create(
      file.path(viz_out_path, m, "indices", dir_name),
      recursive = TRUE, showWarnings = FALSE
    )
  }
}

message("### Terry Plot Standalone Plotting ###")
message("Base path: ", base_path)
message("Output root (created if missing): ", output_base_path)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(forcats)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(egg)
  library(grid)
  library(gridExtra)
  library(extrafont)
})

ensure_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required intermediate not found: ", path)
  }
  path
}

message("Loading intermediates...")

datasets <- readRDS(ensure_exists(file.path(intermediate_files_path, "datasets_for_plots.rds")))
taxonomy_annotations <- readRDS(ensure_exists(file.path(intermediate_files_path, "taxonomy_annotations.rds")))
ubiome_variable_mapping <- readRDS(ensure_exists(file.path(intermediate_files_path, "ubiome_variable_mapping.rds")))
bucket_definitions <- readRDS(ensure_exists(file.path(intermediate_files_path, "bucket_definitions.rds")))
ubiome_genus_mapping_complete <- readRDS(ensure_exists(file.path(intermediate_files_path, "ubiome_genus_mapping_complete.rds")))

plot_configs_none <- readRDS(ensure_exists(file.path(intermediate_files_path, "plot_configs_none.rds")))
plot_configs_hellinger <- readRDS(ensure_exists(file.path(intermediate_files_path, "plot_configs_hellinger.rds")))
plot_configs_clr <- readRDS(ensure_exists(file.path(intermediate_files_path, "plot_configs_clr.rds")))

indices_plot_configs_none <- readRDS(ensure_exists(file.path(intermediate_files_path, "indices_plot_configs_none.rds")))
indices_plot_configs_hellinger <- readRDS(ensure_exists(file.path(intermediate_files_path, "indices_plot_configs_hellinger.rds")))
indices_plot_configs_clr <- readRDS(ensure_exists(file.path(intermediate_files_path, "indices_plot_configs_clr.rds")))

list2env(bucket_definitions, envir = environment())

phylum_info <- taxonomy_annotations$phylum_info
genus_mapping <- taxonomy_annotations$genus_mapping

message("Ensuring Arial fonts are available...")
try({
  if (!"ArialMT" %in% extrafont::fonts()) {
    suppressWarnings(extrafont::font_import(path = file.path(base_path, "data/fonts"),
                                           pattern = "arial.*\\.ttf$",
                                           prompt = FALSE))
  }
  extrafont::loadfonts(device = "pdf", quiet = TRUE)
}, silent = TRUE)

## Default sans avoids extrafont/ArialMT metric issues on some O2/PDF paths. Override: NHANES_TERRY_PDF_FONT=ArialMT
pdf_base_family <- Sys.getenv("NHANES_TERRY_PDF_FONT", unset = "sans")
message("PDF theme base family: ", pdf_base_family)

phylum_colors <- c(
  "Firmicutes"             = "#F38400",
  "Bacteroidetes"          = "#0067A5",
  "Actinobacteria"         = "#8DB600",
  "Proteobacteria"         = "#E68FAC",
  "Fusobacteria"           = "#BE0032",
  "Spirochaetae"           = "#F3C300",
  "Cyanobacteria"          = "#875692",
  "Acidobacteria"          = "#F6A600",
  "Candidate division SR1" = "#2B3D26",
  "Planctomycetes"         = "#332288",
  "Saccharibacteria"       = "#B3446C",
  "Synergistetes"          = "#A1CAF1",
  "Tenericutes"            = "#654522",
  "Verrucomicrobia"        = "#C2B280",
  "unclassified"           = "#DDDDDD",
  "TM7"                    = "#000000"
)

host_variable_colors <- c(
  "DUM" = "#FFFFFF",
  "RID" = "#E69F00", "AGE" = "#56B4E9", "BOR" = "#009E73", "ETH" = "#F0E442", "RIA" = "#CC6677",
  "IND" = "#D55E00", "EDU" = "#DDA0DD", "RSV" = "#999999", "LBX" = "#332288", "URX" = "#44AA99",
  "SMQ" = "#88CCEE", "LBD" = "#117733", "DS1" = "#DDDDDD", "DSQ" = "#0072B2", "HOQ" = "#000000",
  "DR1" = "#AD7700", "DS2" = "#882255", "DR2" = "#661100", "PAQ" = "#6699CC", "GUM" = "#AA4499",
  "ORA" = "#DDCC77", "TOO" = "#117733", "DEN" = "#88CCEE", "CAN" = "#332288", "HEA" = "#E69F00",
  "EMP" = "#56B4E9", "STR" = "#009E73", "ANG" = "#F0E442", "CHD" = "#CC6677", "CVD" = "#D55E00",
  "DIA" = "#DDA0DD", "AST" = "#999999", "BRO" = "#AD7700", "BMX" = "#44AA99", "BPX" = "#88CCEE",
  "CANCER" = "#FF6B6B", "DIABETES" = "#45B7D1", "HEART" = "#4ECDC4", "ASTHMA" = "#96CEB4",
  "BRONCHITIS" = "#FFEAA7", "EMPHYSEMA" = "#DDA0DD", "ANGINA" = "#74B9FF", "STROKE" = "#A0E7E5"
)

grafify_all_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#332288", "#88ccee", "#44aa99", "#117733", "#999933",
  "#ddcc77", "#cc6677", "#882255", "#aa4499", "#77aadd", "#99ddff",
  "#44bb99", "#bbcc33", "#aaaa00", "#eedd88", "#ee8866", "#ffaabb",
  "#bbccee", "#cceeff", "#ccddaa", "#eeeebb", "#ffcccc", "#dddddd",
  "#555555", "#bbbbbb", "#999999", "#ffffff", "#ddaa33", "#bb5566",
  "#004488", "#000000", "#0077bb", "#33bbee", "#009988", "#ee7733",
  "#cc3311", "#ee3377", "#222255", "#225555", "#225522", "#666633",
  "#663333", "#4477aa", "#66ccee", "#228833", "#ccbb44", "#ee6677",
  "#aa3377"
)

FONT_SIZE_PT <- 6
ROW_HEIGHT_PT <- 5
POINTS_TO_MM <- 25.4 / 72
TEXT_MARGIN_PT <- 1

DUMMY_LABEL <- "{DUMMY_PADDING_ROW_FOR_CONSISTENT_MARGINS_ACROSS_ALL_PLOTS_DUMMY_}"

association_colors <- c(
  "Negative" = "#0072B2",
  "Positive" = "#CC6677"
)

## Compact layer: top 5 positive + top 5 negative by |estimate| (same filtered results as full plot)
COMPACT_MAX_POS <- 5L
COMPACT_MAX_NEG <- 5L
## ggrepel annotation labels: top 2 positive + top 2 negative per row (|estimate|, then dot_size, then otu_clean if present)
ANNOTATION_MAX_POS <- 2L
ANNOTATION_MAX_NEG <- 2L
## Base font size (pt) for axis/strip text
AXIS_TEXT_PT <- 6L

## Styling knobs
## - jitter is vertical only (width=0) so effect sizes are unchanged
JITTER_H_COMPACT <- 0.18 * 1.1
JITTER_H_COMPREHENSIVE <- 0.09
## Compact Terry: top 5 pos + top 5 neg per row at full opacity; other sig hits faded
COMPACT_HIGHLIGHT_ALPHA <- 0.65
COMPACT_FADE_ALPHA <- 0.1
## Significance dot-size scale factor (1 = use mapped sizes directly)
POINT_SIZE_SCALE <- 1
## Frozen Terry repel label size (mm); compact host/OTU repel slightly smaller only
LABEL_SIZE_MM <- 1.8
PT_TO_MM <- 25.4 / 72
LABEL_SIZE_MM_COMPACT <- LABEL_SIZE_MM - PT_TO_MM

## Legend bins for significance (maps dot_size -> readable FDR thresholds)
sig_bin_from_dot_size <- function(dot_size) {
  dplyr::case_when(
    dot_size >= 3.5 ~ "FDR < 1e-6",
    dot_size >= 3.0 ~ "FDR < 1e-5",
    dot_size >= 2.5 ~ "FDR < 1e-4",
    dot_size >= 2.0 ~ "FDR < 1e-3",
    dot_size >= 1.5 ~ "FDR < 1e-2",
    TRUE ~ "FDR < 0.05"
  )
}

sig_levels <- c("FDR < 0.05", "FDR < 1e-2", "FDR < 1e-3", "FDR < 1e-4", "FDR < 1e-5", "FDR < 1e-6")
sig_size_values <- setNames(c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5) * POINT_SIZE_SCALE, sig_levels)

clean_taxon_label <- function(x) {
  # Strip square brackets, then keep prefix before underscore.
  # Examples:
  #  "[Eubacterium]_brachy_group" -> "Eubacterium"
  #  "Lachnospiraceae_AC2044_group" -> "Lachnospiraceae"
  #  "Ruminococcaceae_UCG−003" -> "Ruminococcaceae"
  x2 <- gsub("\\[|\\]", "", x)
  sub("_.*$", "", x2)
}

annotation_pos_neg_slice <- function(d, k_pos = ANNOTATION_MAX_POS, k_neg = ANNOTATION_MAX_NEG) {
  if (nrow(d) == 0) {
    return(d[0, , drop = FALSE])
  }
  arrange_ann <- function(dx, positive) {
    dx <- if (positive) {
      dplyr::filter(dx, .data$estimate > 0)
    } else {
      dplyr::filter(dx, .data$estimate < 0)
    }
    if (nrow(dx) == 0) {
      return(dx)
    }
    if ("otu_clean" %in% names(dx)) {
      dx %>% dplyr::arrange(dplyr::desc(abs(.data$estimate)), dplyr::desc(.data$dot_size), .data$otu_clean)
    } else {
      dx %>% dplyr::arrange(dplyr::desc(abs(.data$estimate)), dplyr::desc(.data$dot_size))
    }
  }
  pos <- arrange_ann(d, TRUE) %>% dplyr::slice_head(n = k_pos)
  neg <- arrange_ann(d, FALSE) %>% dplyr::slice_head(n = k_neg)
  dplyr::bind_rows(neg, pos)
}

top_pos_neg_slice <- function(d, k_pos = COMPACT_MAX_POS, k_neg = COMPACT_MAX_NEG) {
  if (nrow(d) == 0) {
    return(d)
  }
  pos <- d %>%
    filter(estimate > 0) %>%
    arrange(desc(abs(estimate)), p.value.fdr, otu_clean) %>%
    slice_head(n = k_pos)
  neg <- d %>%
    filter(estimate < 0) %>%
    arrange(desc(abs(estimate)), p.value.fdr, otu_clean) %>%
    slice_head(n = k_neg)
  bind_rows(neg, pos)
}

## Mark rows that belong to the compact "emphasis" set (top COMPACT_MAX_POS / COMPACT_MAX_NEG per group).
## Same grouping as subset_assoc_rows_for_compact(); used for alpha instead of dropping rows.
mark_assoc_compact_highlight <- function(df2, plot_type, split_by_prefix, host_prefix_ready = FALSE) {
  if (nrow(df2) == 0) {
    return(df2)
  }
  df2 <- df2 %>% dplyr::mutate(rid = dplyr::row_number())

  out <- if (plot_type == "otu") {
    row_key <- if ("Genus" %in% names(df2)) {
      dplyr::if_else(df2$Genus == "unclassified" | is.na(df2$Genus),
                       paste0("NA; ", df2$otu_clean), df2$Genus)
    } else {
      df2$otu_clean
    }
    df2 %>%
      dplyr::mutate(grp_key = row_key) %>%
      dplyr::group_by(.data$grp_key) %>%
      dplyr::group_modify(~ {
        sl <- top_pos_neg_slice(.x)
        .x$compact_highlight <- .x$rid %in% sl$rid
        .x[, names(.x) != "grp_key", drop = FALSE]
      }) %>%
      dplyr::ungroup()
  } else if (plot_type == "host") {
    if (isTRUE(split_by_prefix) && isTRUE(host_prefix_ready) &&
        all(c("host_prefix", "row_label") %in% names(df2))) {
      df2 %>%
        dplyr::group_by(.data$host_prefix, .data$row_label) %>%
        dplyr::group_modify(~ {
          sl <- top_pos_neg_slice(.x)
          .x$compact_highlight <- .x$rid %in% sl$rid
          .x
        }) %>%
        dplyr::ungroup()
    } else if ("host_variable_column" %in% names(df2)) {
      df2 %>%
        dplyr::group_by(.data$host_variable_column) %>%
        dplyr::group_modify(~ {
          sl <- top_pos_neg_slice(.x)
          .x$compact_highlight <- .x$rid %in% sl$rid
          .x
        }) %>%
        dplyr::ungroup()
    } else if ("row_label" %in% names(df2)) {
      df2 %>%
        dplyr::group_by(.data$row_label) %>%
        dplyr::group_modify(~ {
          sl <- top_pos_neg_slice(.x)
          .x$compact_highlight <- .x$rid %in% sl$rid
          .x
        }) %>%
        dplyr::ungroup()
    } else {
      sl <- top_pos_neg_slice(df2)
      df2$compact_highlight <- df2$rid %in% sl$rid
      df2
    }
  } else {
    sl <- top_pos_neg_slice(df2)
    df2$compact_highlight <- df2$rid %in% sl$rid
    df2
  }

  dplyr::select(out, -"rid")
}

subset_assoc_rows_for_compact <- function(df2, plot_type, split_by_prefix, host_prefix_ready = FALSE) {
  if (nrow(df2) == 0) {
    return(df2)
  }

  ## Compact requirement: for each rendered "row", keep up to 5 positive + 5 negative points.
  ## Note: this function is called *before* row_label is constructed for some plot types.
  if (plot_type == "otu") {
    # eventual row_label is based on Genus (or NA/unclassified); build a stable grouping key
    row_key <- if ("Genus" %in% names(df2)) {
      if_else(df2$Genus == "unclassified" | is.na(df2$Genus),
              paste0("NA; ", df2$otu_clean),
              df2$Genus)
    } else {
      df2$otu_clean
    }
    df2 %>%
      mutate(.row_key = row_key) %>%
      group_by(.row_key) %>%
      group_modify(~ top_pos_neg_slice(.x)) %>%
      ungroup() %>%
      select(-.row_key)
  } else if (plot_type == "host") {
    if (isTRUE(split_by_prefix) && isTRUE(host_prefix_ready) &&
        all(c("host_prefix", "row_label") %in% names(df2))) {
      df2 %>%
        group_by(host_prefix, row_label) %>%
        group_modify(~ top_pos_neg_slice(.x)) %>%
        ungroup()
    } else if ("host_variable_column" %in% names(df2)) {
      df2 %>%
        group_by(host_variable_column) %>%
        group_modify(~ top_pos_neg_slice(.x)) %>%
        ungroup()
    } else if ("row_label" %in% names(df2)) {
      df2 %>%
        group_by(row_label) %>%
        group_modify(~ top_pos_neg_slice(.x)) %>%
        ungroup()
    } else {
      top_pos_neg_slice(df2)
    }
  } else {
    top_pos_neg_slice(df2)
  }
}

subset_indices_df_compact <- function(df) {
  top_pos_neg_slice(df)
}

## Compact indices PDF height matches comprehensive: n_rows from full `df` before row subset.
indices_host_compact_dim_row_count <- function(df_dt, mapping_dt, variable_mapping_dt,
                                               dummy_top_bottom, DUMMY_LABEL) {
  plot_data_host <- merge(
    df_dt[, .(otu_clean, estimate, host_var = host_variable_column)],
    mapping_dt[, .(otu, STAIN_INDEX, OXYGEN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX)],
    by.x = "otu_clean", by.y = "otu", all.x = TRUE
  )

  dummy_rows <- data.table(
    otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
    estimate = c(0, 0),
    host_var = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
    STAIN_INDEX = c(0, 1),
    OXYGEN_INDEX = c(0, 1),
    MOTILITY_INDEX = c(0, 1),
    SPORULATION_INDEX = c(0, 1)
  )
  plot_data_host <- rbind(plot_data_host, dummy_rows)
  plot_data_host[, dummy_flag := grepl("^DUMMY_OTU_", otu_clean)]

  value_cols <- c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
  valid_mask <- plot_data_host[, rowSums(!is.na(.SD)) > 0, .SDcols = value_cols]
  valid_assoc_host <- plot_data_host[valid_mask & !dummy_flag, .(host_var, estimate)]

  plot_data_host <- merge(
    plot_data_host,
    variable_mapping_dt,
    by.x = "host_var", by.y = "var_name", all.x = TRUE
  )
  plot_data_host[host_var %in% c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")), var_description := host_var]

  assoc_data_host <- if (nrow(valid_assoc_host)) {
    valid_assoc_host[, .(
      Neg_Assoc = sum(estimate < 0, na.rm = TRUE),
      Pos_Assoc = sum(estimate > 0, na.rm = TRUE)
    ), by = host_var]
  } else {
    data.table(host_var = character(), Neg_Assoc = numeric(), Pos_Assoc = numeric())
  }

  assoc_data_host <- merge(
    assoc_data_host,
    variable_mapping_dt[, .(host_var = var_name, var_description)],
    by = "host_var",
    all.x = TRUE
  )
  assoc_data_host[is.na(var_description), var_description := host_var]
  assoc_data_host <- assoc_data_host[, .(var_description, Neg_Assoc, Pos_Assoc)]
  if (nrow(assoc_data_host)) {
    assoc_data_host <- assoc_data_host[Neg_Assoc != 0 | Pos_Assoc != 0]
  }

  max_abs_count <- if (nrow(assoc_data_host)) {
    max(abs(c(assoc_data_host$Neg_Assoc, assoc_data_host$Pos_Assoc)), na.rm = TRUE)
  } else {
    0
  }
  dummy_assoc <- data.table(
    var_description = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
    Neg_Assoc = c(max_abs_count + 8.5, 0),
    Pos_Assoc = c(0, max_abs_count + 8.5)
  )
  assoc_data_host <- rbind(assoc_data_host, dummy_assoc, fill = TRUE)
  assoc_data_host <- assoc_data_host[, .(
    Neg_Assoc = sum(Neg_Assoc, na.rm = TRUE),
    Pos_Assoc = sum(Pos_Assoc, na.rm = TRUE)
  ), by = var_description]
  assoc_data_host[, Total_Assoc := Neg_Assoc + Pos_Assoc]

  assoc_counts <- assoc_data_host[, .(var_description, Total_Assoc)]
  dummy_rows_order <- assoc_counts[str_detect(var_description, fixed(DUMMY_LABEL))]
  real_rows_order <- assoc_counts[!str_detect(var_description, fixed(DUMMY_LABEL))][order(Total_Assoc, dplyr::desc(var_description))]

  final_order <- if (dummy_top_bottom == "top") {
    c(dummy_rows_order$var_description, real_rows_order$var_description)
  } else if (dummy_top_bottom == "split") {
    dummy_1 <- dummy_rows_order$var_description[1]
    dummy_2 <- dummy_rows_order$var_description[2]
    c(dummy_1, real_rows_order$var_description, dummy_2)
  } else {
    c(real_rows_order$var_description, dummy_rows_order$var_description)
  }

  length(final_order)
}

indices_otu_compact_dim_row_count <- function(df_dt, mapping_dt, dummy_top_bottom, DUMMY_LABEL) {
  unique_otus <- unique(df_dt$otu_clean)
  plot_data_otu <- merge(
    df_dt[, .(otu_clean, estimate, host_var = host_variable_column)],
    mapping_dt[otu %in% unique_otus, .(otu, parsed_genus, STAIN_INDEX, OXYGEN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX)],
    by.x = "otu_clean", by.y = "otu", all.x = TRUE
  )

  dummy_rows <- data.table(
    otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
    estimate = c(0, 0),
    host_var = c("DUMMY_HOST_1", "DUMMY_HOST_2"),
    parsed_genus = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
    STAIN_INDEX = c(0, 1),
    OXYGEN_INDEX = c(0, 1),
    MOTILITY_INDEX = c(0, 1),
    SPORULATION_INDEX = c(0, 1)
  )
  plot_data_otu <- rbind(plot_data_otu, dummy_rows)
  plot_data_otu[, dummy_flag := grepl("^DUMMY_OTU_", otu_clean)]

  value_cols <- c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
  valid_mask <- plot_data_otu[, rowSums(!is.na(.SD)) > 0, .SDcols = value_cols]
  valid_assoc_otu <- plot_data_otu[valid_mask & !dummy_flag, .(parsed_genus, estimate)]

  assoc_data_otu <- if (nrow(valid_assoc_otu)) {
    valid_assoc_otu[, .(
      Neg_Assoc = sum(estimate < 0, na.rm = TRUE),
      Pos_Assoc = sum(estimate > 0, na.rm = TRUE)
    ), by = parsed_genus]
  } else {
    data.table(parsed_genus = character(), Neg_Assoc = numeric(), Pos_Assoc = numeric())
  }

  if (nrow(assoc_data_otu)) {
    assoc_data_otu <- assoc_data_otu[Neg_Assoc != 0 | Pos_Assoc != 0]
  }

  max_abs_count <- if (nrow(assoc_data_otu)) {
    max(abs(c(assoc_data_otu$Neg_Assoc, assoc_data_otu$Pos_Assoc)), na.rm = TRUE)
  } else {
    0
  }
  dummy_assoc <- data.table(
    parsed_genus = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
    Neg_Assoc = c(max_abs_count + 15.5, 0),
    Pos_Assoc = c(0, max_abs_count + 15.5)
  )
  assoc_data_otu <- rbind(assoc_data_otu, dummy_assoc, fill = TRUE)
  assoc_data_otu <- assoc_data_otu[, .(
    Neg_Assoc = sum(Neg_Assoc, na.rm = TRUE),
    Pos_Assoc = sum(Pos_Assoc, na.rm = TRUE)
  ), by = parsed_genus]
  assoc_data_otu[, Total_Assoc := Neg_Assoc + Pos_Assoc]

  assoc_counts <- assoc_data_otu[, .(parsed_genus, Total_Assoc)]
  dummy_rows_order <- assoc_counts[str_detect(parsed_genus, fixed(DUMMY_LABEL))]
  real_rows_order <- assoc_counts[!str_detect(parsed_genus, fixed(DUMMY_LABEL))][order(Total_Assoc, dplyr::desc(parsed_genus))]

  final_order <- if (dummy_top_bottom == "top") {
    c(dummy_rows_order$parsed_genus, real_rows_order$parsed_genus)
  } else if (dummy_top_bottom == "split") {
    dummy_1 <- dummy_rows_order$parsed_genus[1]
    dummy_2 <- dummy_rows_order$parsed_genus[2]
    c(dummy_1, real_rows_order$parsed_genus, dummy_2)
  } else {
    c(real_rows_order$parsed_genus, dummy_rows_order$parsed_genus)
  }

  length(final_order)
}

calc_exact_plot_dimensions <- function(n_rows, plot_width_mm = 300, calibrate_y_mm = 0) {
  total_plot_height_pt <- n_rows * ROW_HEIGHT_PT
  plot_height_mm <- total_plot_height_pt * POINTS_TO_MM
  margin_mm <- 10
  calibrated_height_mm <- plot_height_mm + margin_mm + calibrate_y_mm
  list(width_mm = plot_width_mm,
       height_mm = calibrated_height_mm,
       plot_height_pt = total_plot_height_pt,
       n_rows = n_rows,
       row_height_pt = ROW_HEIGHT_PT,
       calibration_mm = calibrate_y_mm)
}

save_pdf_figure <- function(relative_path, plot_obj, width_mm = 180, height_mm = 170, pointsize = 5) {
  target_path <- file.path(viz_out_path, relative_path)
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  message("  ↳ Saving PDF: ", target_path)
  width_mm <- max(width_mm, 5, na.rm = TRUE)
  height_mm <- max(height_mm, 5, na.rm = TRUE)
  ggsave(
    filename = target_path,
    plot = plot_obj,
    device = "pdf",
    width = width_mm,
    height = height_mm,
    units = "mm",
    family = "ArialMT",
    colormodel = "rgb",
    pointsize = pointsize,
    useDingbats = FALSE
  )
}

filter_by_flexible_prefixes <- function(df, prefix_subsets, variable_column = "host_variable_column") {
  if (is.null(prefix_subsets)) {
    return(df)
  }

  matches <- logical(nrow(df))

  for (prefix in prefix_subsets) {
    prefix_matches <- startsWith(toupper(df[[variable_column]]), toupper(prefix))
    matches <- matches | prefix_matches
  }

  cat("  Prefix filtering results:\n")
  for (prefix in prefix_subsets) {
    prefix_count <- sum(startsWith(toupper(df[[variable_column]]), toupper(prefix)))
    cat("    '", prefix, "': ", prefix_count, " variables significant \t", sep = "")
  }
  cat("  Total significant variables after filtering:", sum(matches), "out of", nrow(df), "\n")

  df[matches, ]
}

extract_variable_prefix <- function(variable_names, target_prefixes = NULL, default_length = 3) {
  if (is.null(target_prefixes)) {
    return(substr(variable_names, 1, default_length))
  }
  prefixes <- character(length(variable_names))
  ordered_prefixes <- target_prefixes[order(nchar(target_prefixes), decreasing = TRUE)]
  for (i in seq_along(variable_names)) {
    var_name <- variable_names[i]
    matched <- FALSE
    for (prefix in ordered_prefixes) {
      if (startsWith(toupper(var_name), toupper(prefix))) {
        prefixes[i] <- prefix
        matched <- TRUE
        break
      }
    }
    if (!matched) prefixes[i] <- substr(var_name, 1, default_length)
  }
  prefixes
}

calc_dot_size <- function(p) {
  case_when(
    p < 1e-6 ~ 3.5,
    p < 1e-5 ~ 3.0,
    p < 1e-4 ~ 2.5,
    p < 1e-3 ~ 2.0,
    p < 1e-2 ~ 1.5,
    TRUE ~ 1.0
  )
}

get_variable_colors <- function(variable_names, base_colors = host_variable_colors) {
  unique_vars <- unique(variable_names)
  color_map <- base_colors
  unmapped_vars <- setdiff(unique_vars, names(color_map))
  if (length(unmapped_vars) > 0) {
    n_needed <- length(unmapped_vars)
    assigned_colors <- grafify_all_colors[seq_len(n_needed)]
    names(assigned_colors) <- unmapped_vars
    color_map <- c(color_map, assigned_colors)
  }
  color_map
}

prepare_dataset <- function(df) {
  dep_matches_otu <- sum(df$dependent_var == df$otu, na.rm = TRUE)
  ind_matches_otu <- sum(df$independent_var == df$otu, na.rm = TRUE)
  df_prepared <- if (dep_matches_otu > ind_matches_otu) {
    df %>% mutate(otu_column = dependent_var, host_variable_column = independent_var)
  } else {
    df %>% mutate(otu_column = independent_var, host_variable_column = dependent_var)
  }
  df_prepared %>% mutate(otu_clean = str_remove(otu, "_relative$"))
}

create_annotations <- function(df_plot, plot_type,
                               k_pos = ANNOTATION_MAX_POS, k_neg = ANNOTATION_MAX_NEG,
                               annotate_na = FALSE) {
  dummy_label <- "{DUMMY_PADDING_ROW_FOR_CONSISTENT_MARGINS_ACROSS_ALL_PLOTS_DUMMY_}"
  ann <- df_plot %>%
    dplyr::filter(.data$row_label != dummy_label) %>%
    dplyr::group_by(.data$row_label) %>%
    dplyr::group_modify(~ annotation_pos_neg_slice(.x, k_pos, k_neg)) %>%
    dplyr::ungroup()

  if (nrow(ann) == 0) {
    return(tibble())
  }

  if (plot_type == "otu") {
    ann %>%
      dplyr::mutate(annotation_label = paste0(.data$host_prefix, " (", round(.data$estimate, 2), ")"))
  } else {
    annotations <- ann %>%
      dplyr::mutate(
        annotation_label_raw = dplyr::if_else(
          .data$Genus == "unclassified" | is.na(.data$Genus), "NA", .data$Genus
        ),
        annotation_label = clean_taxon_label(.data$annotation_label_raw)
      )

    # Compact plots do not annotate NA / Incertae / sedis labels.
    annotations <- annotations %>%
      dplyr::filter(.data$annotation_label != "NA") %>%
      dplyr::filter(!str_detect(tolower(.data$annotation_label_raw), "incertae")) %>%
      dplyr::filter(!str_detect(tolower(.data$annotation_label_raw), "sedis"))

    if (!annotate_na) {
      annotations <- annotations %>%
        dplyr::filter(.data$Genus != "unclassified" & !is.na(.data$Genus))
    }

    annotations %>% dplyr::select(-.data$annotation_label_raw)
  }
}

## --- Terry y-axis order: compact matches comprehensive (order from full association set) ---

terry_dummy_label_string <- function() {
  "{DUMMY_PADDING_ROW_FOR_CONSISTENT_MARGINS_ACROSS_ALL_PLOTS_DUMMY_}"
}

## OTU Terry: add row_label + host_prefix + dummy rows (same as plot_and_save body).
terry_otu_add_labels_and_dummy <- function(df2, prefix_subsets, color_by_variable, dummy_label = terry_dummy_label_string()) {
  carry_highlight <- "compact_highlight" %in% names(df2)
  df2 <- df2 %>%
    mutate(
      row_label = if_else(Genus == "unclassified" | is.na(Genus), paste0("NA; ", otu_clean), Genus),
      host_prefix = if (color_by_variable) host_variable_column else extract_variable_prefix(host_variable_column, prefix_subsets)
    )
  dummy_prefix <- if (!is.null(prefix_subsets) && length(prefix_subsets) > 0) prefix_subsets[1] else "DUM"
  dummy_rows <- tibble(
    otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
    estimate = c(0, 0),
    p.value.fdr = c(1e-12, 1e-12),
    Genus = c(dummy_label, dummy_label),
    row_label = c(dummy_label, paste0(dummy_label, "_2")),
    dot_size = c(3.5, 3.5),
    host_prefix = c(dummy_prefix, dummy_prefix)
  )
  if (carry_highlight) {
    dummy_rows$compact_highlight <- c(TRUE, TRUE)
  }
  bind_rows(df2, dummy_rows)
}

terry_otu_y_order_levels <- function(df2_labeled, dummy_top_bottom) {
  assoc_counts <- df2_labeled %>%
    dplyr::group_by(.data$row_label) %>%
    dplyr::summarise(Total_Assoc = dplyr::n(), .groups = "drop")
  dummy_rows_df <- assoc_counts %>% dplyr::filter(str_detect(.data$row_label, "DUMMY_PADDING"))
  real_rows_df <- assoc_counts %>%
    dplyr::filter(!str_detect(.data$row_label, "DUMMY_PADDING")) %>%
    dplyr::arrange(.data$Total_Assoc, dplyr::desc(.data$row_label))
  switch(dummy_top_bottom,
    top = c(dummy_rows_df$row_label, real_rows_df$row_label),
    split = c(dummy_rows_df$row_label[1], real_rows_df$row_label, dummy_rows_df$row_label[2]),
    bottom = c(real_rows_df$row_label, dummy_rows_df$row_label)
  )
}

terry_factor_row_order <- function(level_order, present_labels) {
  keep <- level_order[level_order %in% present_labels]
  keep
}

terry_host_join_meta <- function(df2, split_by_prefix, prefix_subsets, phylum_df, ubiome_variable_mapping) {
  df2 %>%
    left_join(phylum_df %>% distinct(otu, .keep_all = TRUE), by = c("otu_clean" = "otu")) %>%
    left_join(ubiome_variable_mapping %>% select(var_name, var_description), by = c("host_variable_column" = "var_name")) %>%
    mutate(
      host_prefix = if (split_by_prefix) extract_variable_prefix(host_variable_column, prefix_subsets) else "All",
      row_label = if_else(is.na(var_description), host_variable_column, var_description)
    )
}

terry_host_bind_dummy_rows <- function(df2, split_by_prefix, dummy_label = terry_dummy_label_string()) {
  carry_highlight <- "compact_highlight" %in% names(df2)
  dummy_rows <- if (split_by_prefix) {
    map_dfr(unique(df2$host_prefix), ~ {
      dr <- tibble(
        otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
        estimate = c(0, 0),
        p.value.fdr = c(1e-12, 1e-12),
        Phylum = c("Unclassified", "Unclassified"),
        Genus = c(dummy_label, dummy_label),
        row_label = c(dummy_label, paste0(dummy_label, "_2")),
        dot_size = c(3.5, 3.5),
        host_prefix = c(.x, .x),
        host_variable_column = c(dummy_label, dummy_label)
      )
      if (carry_highlight) dr$compact_highlight <- c(TRUE, TRUE)
      dr
    })
  } else {
    dr <- tibble(
      otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
      estimate = c(0, 0),
      p.value.fdr = c(1e-12, 1e-12),
      Phylum = c("Unclassified", "Unclassified"),
      Genus = c(dummy_label, dummy_label),
      row_label = c(dummy_label, paste0(dummy_label, "_2")),
      dot_size = c(3.5, 3.5),
      host_prefix = c("All", "All"),
      host_variable_column = c(dummy_label, dummy_label)
    )
    if (carry_highlight) dr$compact_highlight <- c(TRUE, TRUE)
    dr
  }
  bind_rows(df2, dummy_rows)
}

terry_host_y_order_levels_nosplit <- function(df2_with_dummy, dummy_top_bottom) {
  assoc_counts <- df2_with_dummy %>%
    dplyr::group_by(.data$row_label) %>%
    dplyr::summarise(Total_Assoc = dplyr::n(), .groups = "drop")
  dummy_rows_df <- assoc_counts %>% dplyr::filter(str_detect(.data$row_label, "DUMMY_PADDING"))
  real_rows_df <- assoc_counts %>%
    dplyr::filter(!str_detect(.data$row_label, "DUMMY_PADDING")) %>%
    dplyr::arrange(.data$Total_Assoc, dplyr::desc(.data$row_label))
  switch(dummy_top_bottom,
    top = c(dummy_rows_df$row_label, real_rows_df$row_label),
    split = c(dummy_rows_df$row_label[1], real_rows_df$row_label, dummy_rows_df$row_label[2]),
    bottom = c(real_rows_df$row_label, dummy_rows_df$row_label)
  )
}

## Compact split-prefix host plots: same facet row-label order as comprehensive (assoc_counts from full data).
terry_host_apply_split_y_order <- function(df2_plot_with_dummy, assoc_counts_full) {
  df2_plot_with_dummy %>%
    dplyr::group_by(.data$host_prefix) %>%
    dplyr::group_modify(function(rows, key) {
      hp <- key$host_prefix[[1]]
      present <- unique(as.character(rows$row_label))
      full_ord <- unique(assoc_counts_full$row_label[assoc_counts_full$host_prefix == hp])
      lv <- terry_factor_row_order(full_ord, present)
      rows$row_label <- factor(rows$row_label, levels = lv)
      rows
    }) %>%
    dplyr::ungroup()
}

save_plot_object <- function(plot_list, transformation, genus_mapping, ubiome_variable_mapping,
                             phylum_info, phylum_colors) {
  for (config in plot_list) {
    if (is.null(config$df) || nrow(config$df) == 0) {
      next
    }
    for (fig_mode in c("comprehensive", "compact")) {
      message("  Plot: ", config$file_prefix, " [", fig_mode, "]")
      plot_and_save(
        df = config$df,
        plot_type = config$plot_type,
        file_prefix = file.path(fig_mode, transformation, config$file_prefix),
        width_mm = config$width_mm,
        calibrate_y_mm = config$calibrate_y_mm %||% 0,
        genus_mapping = genus_mapping,
        ubiome_variable_mapping = ubiome_variable_mapping,
        phylum_df = phylum_info,
        phylum_colors = phylum_colors,
        split_by_prefix = config$split_by_prefix %||% FALSE,
        prefix_subsets = config$prefix_subsets,
        remove_na = config$remove_na %||% FALSE,
        remove_na_unclassified = config$remove_na_unclassified %||% FALSE,
        annotate_na = config$annotate_na %||% FALSE,
        xlim = config$xlim,
        add_annotations = config$add_annotations %||% TRUE,
        show_legend = config$show_legend %||% FALSE,
        color_by_variable = config$color_by_variable %||% FALSE,
        dummy_top_bottom = config$dummy_top_bottom %||% "bottom",
        fig_mode = fig_mode
      )
    }
  }
}

plot_and_save <- function(df, plot_type = "otu", file_prefix, width_mm = 115, calibrate_y_mm = 0,
                          genus_mapping, ubiome_variable_mapping, phylum_df, phylum_colors,
                          split_by_prefix = FALSE, prefix_subsets = NULL, remove_na = FALSE,
                          remove_na_unclassified = FALSE, annotate_na = FALSE, xlim = NULL,
                          add_annotations = TRUE, show_legend = FALSE,
                          color_by_variable = FALSE, dummy_top_bottom = "bottom",
                          fig_mode = "comprehensive") {

  df <- prepare_dataset(df)
  if (!is.null(prefix_subsets)) {
    df <- filter_by_flexible_prefixes(df, prefix_subsets)
    if (nrow(df) == 0) {
      warning("No variables found matching prefixes for ", file_prefix)
      return(NULL)
    }
  }

  df2_core <- df %>%
    distinct(otu_clean, estimate, .keep_all = TRUE) %>%
    left_join(genus_mapping %>% distinct(otu, .keep_all = TRUE), by = c("otu_clean" = "otu")) %>%
    mutate(dot_size = calc_dot_size(p.value.fdr)) %>%
    filter(!is.na(dot_size))

  if (remove_na_unclassified) {
    df2_core <- df2_core %>% filter(Genus != "unclassified" & !is.na(Genus))
  }

  if (nrow(df2_core) == 0) {
    warning("No significant results for ", plot_type, " plot: ", file_prefix)
    return(NULL)
  }

  ## Compact Terry: keep all points; mark top 5 + / 5 − per group for full opacity (see mark_assoc_compact_highlight).
  if (fig_mode == "compact") {
    if (plot_type == "otu") {
      df2_core <- mark_assoc_compact_highlight(df2_core, "otu", isTRUE(split_by_prefix), FALSE)
    } else if (plot_type == "host" && !isTRUE(split_by_prefix)) {
      df2_core <- mark_assoc_compact_highlight(df2_core, "host", FALSE, FALSE)
    }
  }

  dummy_label <- terry_dummy_label_string()

  label_pt <- AXIS_TEXT_PT
  genus_y_axis_pt <- label_pt
  repel_genus_mm <- if (fig_mode == "compact") LABEL_SIZE_MM_COMPACT else LABEL_SIZE_MM

  if (plot_type == "otu") {
    df_full_labeled <- terry_otu_add_labels_and_dummy(df2_core, prefix_subsets, color_by_variable, dummy_label)
    comp_order <- terry_otu_y_order_levels(df_full_labeled, dummy_top_bottom)
    df2 <- terry_otu_add_labels_and_dummy(df2_core, prefix_subsets, color_by_variable, dummy_label)
    df2$row_label <- factor(
      df2$row_label,
      levels = terry_factor_row_order(comp_order, unique(as.character(df2$row_label)))
    )

    n_rows <- length(levels(df2$row_label))
    ## Jitter (vertical only): separate overlapping points without changing effect sizes (x).
    ## Apply the same position to points and annotation anchors so leader lines still target dots.
    pos_jitter <- position_jitter(
      width = 0,
      height = if (fig_mode == "compact") JITTER_H_COMPACT else JITTER_H_COMPREHENSIVE
    )

    ## Labels: top ANNOTATION_MAX_POS positive + ANNOTATION_MAX_NEG negative per row (see annotation_pos_neg_slice).
    df_annotate <- if (add_annotations) create_annotations(df2, "otu", annotate_na = TRUE) else tibble()

    if (is.null(xlim)) {
      est_excl <- abs(df2$estimate[as.character(df2$row_label) != dummy_label])
      max_abs_estimate <- if (length(est_excl)) max(est_excl, na.rm = TRUE) else 0
      if (!is.finite(max_abs_estimate) || max_abs_estimate < 1e-12) {
        xlim <- c(-1, 1)
      } else {
        xlim <- c(-max_abs_estimate * 1.15, max_abs_estimate * 1.15)
      }
    }

    color_map <- get_variable_colors(df2$host_prefix)

    df2 <- df2 %>% mutate(sig_bin = factor(sig_bin_from_dot_size(dot_size), levels = sig_levels))
    df2 <- df2 %>%
      dplyr::mutate(
        point_alpha = if (fig_mode == "compact") {
          dplyr::case_when(
            str_detect(.data$row_label, "DUMMY_PADDING") ~ COMPACT_HIGHLIGHT_ALPHA,
            .data$compact_highlight %in% TRUE ~ COMPACT_HIGHLIGHT_ALPHA,
            TRUE ~ COMPACT_FADE_ALPHA
          )
        } else {
          COMPACT_HIGHLIGHT_ALPHA
        }
      )
    p <- ggplot(df2, aes(x = estimate, y = row_label, fill = host_prefix, size = sig_bin)) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.3, alpha = 0.8) +
      geom_point(aes(alpha = .data$point_alpha), shape = 21, stroke = 0.2, position = pos_jitter) +
      scale_alpha_identity(guide = "none") +
      {if (add_annotations && nrow(df_annotate) > 0) geom_text_repel(
        data = df_annotate,
        aes(label = annotation_label),
        size = if (fig_mode == "compact") LABEL_SIZE_MM_COMPACT else LABEL_SIZE_MM,
        box.padding = 0.5,
        point.padding = 0,
        segment.color = "gray60",
        segment.size = 0.15,
        max.overlaps = 25,
        min.segment.length = 0.05,
        position = pos_jitter
      ) else NULL} +
      {if (show_legend) scale_size_manual(
         name = "Significance (BH FDR)",
         values = sig_size_values,
         drop = FALSE,
         guide = guide_legend(override.aes = list(fill = "gray50"))
       ) else scale_size_manual(values = sig_size_values, guide = "none")} +
      scale_fill_manual(values = color_map) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(limits = xlim, expand = c(0.01, 0.01)) +
      egg::theme_article(base_size = label_pt, base_family = pdf_base_family) +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = genus_y_axis_pt, margin = margin(r = 1, unit = "pt"), hjust = 1, lineheight = 0.9, family = pdf_base_family),
        axis.text.x = element_text(size = label_pt, margin = margin(t = 1, unit = "pt"), family = pdf_base_family),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL
      )
  } else {
    df_core_joined <- terry_host_join_meta(df2_core, split_by_prefix, prefix_subsets, phylum_df, ubiome_variable_mapping)

    if (fig_mode == "compact" && isTRUE(split_by_prefix)) {
      df_core_joined <- mark_assoc_compact_highlight(df_core_joined, "host", TRUE, TRUE)
    }

    df_plot_joined <- df_core_joined

    df_full_dummy <- terry_host_bind_dummy_rows(df_core_joined, split_by_prefix, dummy_label)

    if (!split_by_prefix) {
      comp_order_host <- terry_host_y_order_levels_nosplit(df_full_dummy, dummy_top_bottom)
      df2 <- terry_host_bind_dummy_rows(df_plot_joined, split_by_prefix, dummy_label)
      df2$row_label <- factor(
        df2$row_label,
        levels = terry_factor_row_order(comp_order_host, unique(as.character(df2$row_label)))
      )
    } else {
      assoc_counts_full <- df_full_dummy %>%
        dplyr::group_by(.data$host_prefix, .data$row_label) %>%
        dplyr::summarise(Total_Assoc = dplyr::n(), .groups = "drop")
      df2 <- terry_host_bind_dummy_rows(df_plot_joined, split_by_prefix, dummy_label)
      df2 <- terry_host_apply_split_y_order(df2, assoc_counts_full)
    }

    n_rows <- length(levels(df2$row_label))
    pos_jitter <- position_jitter(
      width = 0,
      height = if (fig_mode == "compact") JITTER_H_COMPACT else JITTER_H_COMPREHENSIVE
    )

    df_annotate <- if (add_annotations) create_annotations(df2, "host", annotate_na = TRUE) else tibble()

    if (is.null(xlim)) {
      est_excl <- abs(df2$estimate[as.character(df2$row_label) != dummy_label])
      max_abs_estimate <- if (length(est_excl)) max(est_excl, na.rm = TRUE) else 0
      if (!is.finite(max_abs_estimate) || max_abs_estimate < 1e-12) {
        xlim <- c(-1, 1)
      } else {
        xlim <- c(-max_abs_estimate * 1.15, max_abs_estimate * 1.15)
      }
    }

    df2 <- df2 %>% mutate(sig_bin = factor(sig_bin_from_dot_size(dot_size), levels = sig_levels))
    df2 <- df2 %>%
      dplyr::mutate(
        point_alpha = if (fig_mode == "compact") {
          dplyr::case_when(
            str_detect(.data$row_label, "DUMMY_PADDING") ~ COMPACT_HIGHLIGHT_ALPHA,
            .data$compact_highlight %in% TRUE ~ COMPACT_HIGHLIGHT_ALPHA,
            TRUE ~ COMPACT_FADE_ALPHA
          )
        } else {
          COMPACT_HIGHLIGHT_ALPHA
        }
      )
    p <- ggplot(df2, aes(x = estimate, y = row_label, fill = Phylum, size = sig_bin)) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.3, alpha = 0.8) +
      geom_point(aes(alpha = .data$point_alpha), shape = 21, stroke = 0.15, position = pos_jitter) +
      scale_alpha_identity(guide = "none") +
      {if (add_annotations && nrow(df_annotate) > 0) geom_text_repel(
        data = df_annotate,
        aes(label = annotation_label),
        size = repel_genus_mm,
        box.padding = 0.5,
        point.padding = 0,
        segment.color = "gray60",
        segment.size = 0.15,
        max.overlaps = 25,
        min.segment.length = 0.05,
        position = pos_jitter
      ) else NULL} +
      {if (show_legend) scale_size_manual(
         name = "Significance (BH FDR)",
         values = sig_size_values,
         drop = FALSE,
         guide = guide_legend(override.aes = list(fill = "gray50"))
       ) else scale_size_manual(values = sig_size_values, guide = "none")} +
      scale_fill_manual(values = phylum_colors, na.value = "grey90") +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(limits = xlim, expand = c(0.0005, 0.0005)) +
      egg::theme_article(base_size = label_pt, base_family = pdf_base_family) +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = label_pt, margin = margin(r = 1, unit = "pt"), hjust = 1, lineheight = 0.9, family = pdf_base_family),
        axis.text.x = element_text(size = label_pt, margin = margin(t = 1, unit = "pt"), family = pdf_base_family),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL
      )
    if (split_by_prefix) {
      p <- p + facet_wrap(~host_prefix, scales = "free_y")
    }
  }

  dims <- calc_exact_plot_dimensions(n_rows, width_mm, calibrate_y_mm)
  out_w <- dims$width_mm
  out_h <- dims$height_mm
  output_file <- paste0(file_prefix, ".pdf")
  save_pdf_figure(
    output_file, p, width_mm = out_w, height_mm = out_h,
    pointsize = 5
  )
  invisible(p)
}

plot_indices <- function(df, plot_type = "host", file_prefix = "indices", width_mm = 165, calibrate_y_mm = 0,
                         mapping_data, genus_mapping, ubiome_variable_mapping,
                         prefix_subsets = NULL, remove_na = FALSE, remove_na_unclassified = FALSE,
                         show_legend = TRUE, dummy_top_bottom = "split",
                         hide_axis_title = FALSE, geom_point_shape = 21, geom_point_size = 3.5,
                         fig_mode = "comprehensive") {

  idx_label_pt <- AXIS_TEXT_PT
  idx_y_genus_pt <- idx_label_pt

  df <- prepare_dataset(df)

  if (!is.null(prefix_subsets)) {
    cat("  Original variables:", length(unique(df$host_variable_column)), "\n")
    df <- filter_by_flexible_prefixes(df, prefix_subsets, "host_variable_column")

    if (nrow(df) == 0) {
      warning(paste("No variables found matching prefixes:", paste(prefix_subsets, collapse = ", ")))
      return(NULL)
    }
  }

  df <- df %>%
    left_join(genus_mapping %>% distinct(otu, .keep_all = TRUE), by = c("otu_clean" = "otu"))

  if (remove_na_unclassified) {
    df <- df %>% filter(Genus != "unclassified" & !is.na(Genus))
  }

  if (nrow(df) == 0) {
    warning(paste("No data remaining after filtering for", file_prefix))
    return(NULL)
  }

  df_for_compact_dims <- df

  if (fig_mode == "compact") {
    df <- as.data.frame(subset_indices_df_compact(df))
    if (nrow(df) == 0) {
      message("  [compact] Skipping indices (no rows after subset): ", file_prefix)
      return(invisible(NULL))
    }
  }

  setDT(df)
  mapping_dt <- as.data.table(mapping_data)
  variable_mapping_dt <- as.data.table(ubiome_variable_mapping)

  if (plot_type == "host") {
    plot_data_host <- merge(
      df[, .(otu_clean, estimate, host_var = host_variable_column)],
      mapping_dt[, .(otu, STAIN_INDEX, OXYGEN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX)],
      by.x = "otu_clean", by.y = "otu", all.x = TRUE
    )

    dummy_rows <- data.table(
      otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
      estimate = c(0, 0),
      host_var = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
      STAIN_INDEX = c(0, 1),
      OXYGEN_INDEX = c(0, 1),
      MOTILITY_INDEX = c(0, 1),
      SPORULATION_INDEX = c(0, 1)
    )
    plot_data_host <- rbind(plot_data_host, dummy_rows)
    plot_data_host[, dummy_flag := grepl("^DUMMY_OTU_", otu_clean)]

    value_cols <- c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
    valid_mask <- plot_data_host[, rowSums(!is.na(.SD)) > 0, .SDcols = value_cols]
    valid_assoc_host <- plot_data_host[valid_mask & !dummy_flag, .(host_var, estimate)]

    plot_data_host <- merge(
      plot_data_host,
      variable_mapping_dt,
      by.x = "host_var", by.y = "var_name", all.x = TRUE
    )
    plot_data_host[host_var %in% c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")), var_description := host_var]

    plot_data_long_host <- melt(
      plot_data_host,
      id.vars = c("otu_clean", "host_var", "var_description", "estimate"),
      measure.vars = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
      variable.name = "Index_Type",
      value.name = "Index_Value"
    )
    plot_data_long_host <- plot_data_long_host[!is.na(Index_Value)]
    plot_data_long_host[, dummy_flag := grepl("^DUMMY_OTU_", otu_clean)]
    plot_data_long_host <- plot_data_long_host[
      ,
      if (var_description %in% c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")) || any(!dummy_flag)) .SD else NULL,
      by = var_description
    ]
    plot_data_long_host[, dummy_flag := NULL]

    if (hide_axis_title) {
      plot_data_long_host[, Index_Type := factor(
        Index_Type,
        levels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
        labels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
      )]
    } else {
      plot_data_long_host[, Index_Type := factor(
        Index_Type,
        levels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
        labels = c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)")
      )]
    }
    plot_data_long_host[, Index_Type := as.character(Index_Type)]
    plot_data_long_host[, Index_Value := as.numeric(Index_Value)]
    plot_data_long_host[, Assoc_Type := factor(ifelse(estimate < 0, "Negative", "Positive"),
                                               levels = c("Negative", "Positive"))]

    assoc_data_host <- if (nrow(valid_assoc_host)) {
      valid_assoc_host[, .(
        Neg_Assoc = sum(estimate < 0, na.rm = TRUE),
        Pos_Assoc = sum(estimate > 0, na.rm = TRUE)
      ), by = host_var]
    } else {
      data.table(host_var = character(), Neg_Assoc = numeric(), Pos_Assoc = numeric())
    }

    assoc_data_host <- merge(
      assoc_data_host,
      variable_mapping_dt[, .(host_var = var_name, var_description)],
      by = "host_var",
      all.x = TRUE
    )
    assoc_data_host[is.na(var_description), var_description := host_var]
    assoc_data_host <- assoc_data_host[, .(var_description, Neg_Assoc, Pos_Assoc)]
    if (nrow(assoc_data_host)) {
      assoc_data_host <- assoc_data_host[Neg_Assoc != 0 | Pos_Assoc != 0]
    }

    max_abs_count <- if (nrow(assoc_data_host)) {
      max(abs(c(assoc_data_host$Neg_Assoc, assoc_data_host$Pos_Assoc)), na.rm = TRUE)
    } else {
      0
    }
    dummy_assoc <- data.table(
      var_description = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
      Neg_Assoc = c(max_abs_count + 8.5, 0),
      Pos_Assoc = c(0, max_abs_count + 8.5)
    )
    assoc_data_host <- rbind(assoc_data_host, dummy_assoc, fill = TRUE)
    assoc_data_host <- assoc_data_host[, .(
      Neg_Assoc = sum(Neg_Assoc, na.rm = TRUE),
      Pos_Assoc = sum(Pos_Assoc, na.rm = TRUE)
    ), by = var_description]
    assoc_data_host[, Total_Assoc := Neg_Assoc + Pos_Assoc]

    assoc_counts <- assoc_data_host[, .(var_description, Total_Assoc)]
    dummy_rows_order <- assoc_counts[str_detect(var_description, fixed(DUMMY_LABEL))]
    real_rows_order <- assoc_counts[!str_detect(var_description, fixed(DUMMY_LABEL))][order(Total_Assoc, dplyr::desc(var_description))]

    final_order <- if (dummy_top_bottom == "top") {
      c(dummy_rows_order$var_description, real_rows_order$var_description)
    } else if (dummy_top_bottom == "split") {
      dummy_1 <- dummy_rows_order$var_description[1]
      dummy_2 <- dummy_rows_order$var_description[2]
      c(dummy_1, real_rows_order$var_description, dummy_2)
    } else {
      c(real_rows_order$var_description, dummy_rows_order$var_description)
    }

    assoc_data_long_host <- melt(
      assoc_data_host,
      id.vars = "var_description",
      measure.vars = c("Neg_Assoc", "Pos_Assoc"),
      variable.name = "Assoc_Type",
      value.name = "Count"
    )
    assoc_data_long_host[Assoc_Type == "Neg_Assoc", Count := -Count]
    assoc_data_long_host[, Assoc_Type := factor(
      Assoc_Type,
      levels = c("Neg_Assoc", "Pos_Assoc"),
      labels = c("Negative", "Positive")
    )]
    assoc_data_long_host[, Index_Type := "Association Direction Count"]
    assoc_data_long_host[, Index_Type := as.character(Index_Type)]

    plot_data_combined <- rbind(
      plot_data_long_host[, .(var_description, Index_Type, Value = Index_Value, Assoc_Type)],
      assoc_data_long_host[, .(var_description, Index_Type, Value = Count, Assoc_Type)]
    )

    if (hide_axis_title) {
      plot_data_combined[, Index_Type := factor(
        Index_Type,
        levels = c("Association Direction Count", "STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
      )]
    } else {
      plot_data_combined[, Index_Type := factor(
        Index_Type,
        levels = c("Association Direction Count", "Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)")
      )]
    }

    plot_data_combined[, var_description := factor(var_description, levels = final_order)]
    n_rows_plot <- length(levels(plot_data_combined$var_description))
    n_rows <- if (fig_mode == "compact") {
      indices_host_compact_dim_row_count(
        data.table::as.data.table(as.data.frame(df_for_compact_dims)),
        mapping_dt,
        variable_mapping_dt,
        dummy_top_bottom,
        DUMMY_LABEL
      )
    } else {
      n_rows_plot
    }

    mean_data <- plot_data_combined[Index_Type != "Association Direction Count",
                                   .(Mean_Value = mean(Value, na.rm = TRUE)),
                                   by = .(var_description, Index_Type, Assoc_Type)]

    assoc_max <- ceiling(max(abs(assoc_data_long_host$Count), na.rm = TRUE) / 5) * 5
    assoc_breaks <- seq(-assoc_max, assoc_max, by = 5)

    p <- ggplot(plot_data_combined, aes(x = Value, y = var_description)) +
      geom_violin(
        data = plot_data_combined[Index_Type != "Association Direction Count"],
        aes(fill = Assoc_Type),
        alpha = 0.65, color = "black", linewidth = 0.1, na.rm = TRUE,
        width = 0.8, scale = "width", position = "identity"
      ) +
      geom_point(
        data = mean_data,
        aes(x = Mean_Value, y = var_description, colour = Assoc_Type, fill = Assoc_Type),
        shape = geom_point_shape, size = geom_point_size / 2, alpha = 0.8, stroke = 0.1, na.rm = TRUE
      ) +
      geom_bar(
        data = plot_data_combined[Index_Type == "Association Direction Count"],
        aes(fill = Assoc_Type), stat = "identity", width = 0.8, alpha = 0.65
      ) +
      geom_text(
        data = plot_data_combined[Index_Type == "Association Direction Count" & Value != 0],
        aes(
          x = ifelse(Value < 0, Value - max(abs(Value)) * 0.02, Value + max(abs(Value)) * 0.02),
          label = abs(Value)
        ),
        size = 1.5, color = "#555555",
        hjust = ifelse(plot_data_combined[Index_Type == "Association Direction Count" & Value != 0]$Value < 0, 1, 0)
      ) +
      scale_fill_manual(values = association_colors, na.translate = FALSE) +
      scale_color_manual(values = association_colors, na.translate = FALSE) +
      geom_vline(
        aes(xintercept = ifelse(Index_Type == "Association Direction Count", 0, 0.5)),
        linetype = "dashed", color = "black", linewidth = 0.2
      ) +
      geom_blank(data = data.frame(
        Index_Type = if (hide_axis_title) {
          factor(
            c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
            levels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
          )
        } else {
          factor(
            c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)"),
            levels = c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)")
          )
        },
        Value = rep(c(0, 1), 2),
        y_var = rep(levels(plot_data_combined$var_description)[1], 8)
      ), aes(x = Value, y = y_var), inherit.aes = FALSE) +
      facet_wrap(~Index_Type, nrow = 1, scales = "free_x") +
      scale_x_continuous(
        breaks = function(x) {
          if (length(x) == 0) return(NULL)
          if (max(x, na.rm = TRUE) <= 1) {
            c(0, 0.5, 1)
          } else {
            assoc_breaks
          }
        },
        labels = function(x) {
          if (length(x) == 0) return(NULL)
          x
        }
      ) +
      coord_cartesian(xlim = c(NA, NA)) +
      scale_y_discrete(expand = c(0, 0)) +
      egg::theme_article(base_size = idx_label_pt, base_family = pdf_base_family) +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = idx_label_pt, margin = margin(r = TEXT_MARGIN_PT, unit = "pt"), hjust = 1, lineheight = 0.9),
        axis.text.x = element_text(size = idx_label_pt, margin = margin(t = TEXT_MARGIN_PT, unit = "pt")),
        axis.ticks.y = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.x = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL,
        legend.title = element_blank(),
        legend.text = element_text(size = idx_label_pt),
        strip.text = if (hide_axis_title) element_blank() else element_text(size = idx_label_pt, angle = 90, hjust = 0)
      )

    file_suffix <- "_host_indices.pdf"

  } else if (plot_type == "otu") {
    unique_otus <- unique(df$otu_clean)
    plot_data_otu <- merge(
      df[, .(otu_clean, estimate, host_var = host_variable_column)],
      mapping_dt[otu %in% unique_otus, .(otu, parsed_genus, STAIN_INDEX, OXYGEN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX)],
      by.x = "otu_clean", by.y = "otu", all.x = TRUE
    )

    dummy_rows <- data.table(
      otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
      estimate = c(0, 0),
      host_var = c("DUMMY_HOST_1", "DUMMY_HOST_2"),
      parsed_genus = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
      STAIN_INDEX = c(0, 1),
      OXYGEN_INDEX = c(0, 1),
      MOTILITY_INDEX = c(0, 1),
      SPORULATION_INDEX = c(0, 1)
    )
    plot_data_otu <- rbind(plot_data_otu, dummy_rows)
    plot_data_otu[, dummy_flag := grepl("^DUMMY_OTU_", otu_clean)]

    value_cols <- c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
    valid_mask <- plot_data_otu[, rowSums(!is.na(.SD)) > 0, .SDcols = value_cols]
    valid_assoc_otu <- plot_data_otu[valid_mask & !dummy_flag, .(parsed_genus, estimate)]

    plot_data_long_otu <- melt(
      plot_data_otu,
      id.vars = c("otu_clean", "parsed_genus", "estimate", "host_var"),
      measure.vars = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
      variable.name = "Index_Type",
      value.name = "Index_Value"
    )
    plot_data_long_otu <- plot_data_long_otu[!is.na(Index_Value)]
    plot_data_long_otu[, dummy_flag := grepl("^DUMMY_OTU_", otu_clean)]
    plot_data_long_otu <- plot_data_long_otu[
      ,
      if (parsed_genus %in% c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")) || any(!dummy_flag)) .SD else NULL,
      by = parsed_genus
    ]
    plot_data_long_otu[, dummy_flag := NULL]

    if (hide_axis_title) {
      plot_data_long_otu[, Index_Type := factor(
        Index_Type,
        levels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
        labels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
      )]
    } else {
      plot_data_long_otu[, Index_Type := factor(
        Index_Type,
        levels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
        labels = c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)")
      )]
    }
    plot_data_long_otu[, Index_Type := as.character(Index_Type)]
    plot_data_long_otu[, Index_Value := as.numeric(Index_Value)]
    plot_data_long_otu[, Assoc_Type := factor(ifelse(estimate < 0, "Negative", "Positive"),
                                             levels = c("Negative", "Positive"))]

    assoc_data_otu <- if (nrow(valid_assoc_otu)) {
      valid_assoc_otu[, .(
        Neg_Assoc = sum(estimate < 0, na.rm = TRUE),
        Pos_Assoc = sum(estimate > 0, na.rm = TRUE)
      ), by = parsed_genus]
    } else {
      data.table(parsed_genus = character(), Neg_Assoc = numeric(), Pos_Assoc = numeric())
    }

    if (nrow(assoc_data_otu)) {
      assoc_data_otu <- assoc_data_otu[Neg_Assoc != 0 | Pos_Assoc != 0]
    }

    max_abs_count <- if (nrow(assoc_data_otu)) {
      max(abs(c(assoc_data_otu$Neg_Assoc, assoc_data_otu$Pos_Assoc)), na.rm = TRUE)
    } else {
      0
    }
    dummy_assoc <- data.table(
      parsed_genus = c(DUMMY_LABEL, paste0(DUMMY_LABEL, "_2")),
      Neg_Assoc = c(max_abs_count + 15.5, 0),
      Pos_Assoc = c(0, max_abs_count + 15.5)
    )
    assoc_data_otu <- rbind(assoc_data_otu, dummy_assoc, fill = TRUE)
    assoc_data_otu <- assoc_data_otu[, .(
      Neg_Assoc = sum(Neg_Assoc, na.rm = TRUE),
      Pos_Assoc = sum(Pos_Assoc, na.rm = TRUE)
    ), by = parsed_genus]
    assoc_data_otu[, Total_Assoc := Neg_Assoc + Pos_Assoc]

    assoc_counts <- assoc_data_otu[, .(parsed_genus, Total_Assoc)]
    dummy_rows_order <- assoc_counts[str_detect(parsed_genus, fixed(DUMMY_LABEL))]
    real_rows_order <- assoc_counts[!str_detect(parsed_genus, fixed(DUMMY_LABEL))][order(Total_Assoc, dplyr::desc(parsed_genus))]

    final_order <- if (dummy_top_bottom == "top") {
      c(dummy_rows_order$parsed_genus, real_rows_order$parsed_genus)
    } else if (dummy_top_bottom == "split") {
      dummy_1 <- dummy_rows_order$parsed_genus[1]
      dummy_2 <- dummy_rows_order$parsed_genus[2]
      c(dummy_1, real_rows_order$parsed_genus, dummy_2)
    } else {
      c(real_rows_order$parsed_genus, dummy_rows_order$parsed_genus)
    }

    assoc_data_long_otu <- melt(
      assoc_data_otu,
      id.vars = "parsed_genus",
      measure.vars = c("Neg_Assoc", "Pos_Assoc"),
      variable.name = "Assoc_Type",
      value.name = "Count"
    )
    assoc_data_long_otu[Assoc_Type == "Neg_Assoc", Count := -Count]
    assoc_data_long_otu[, Assoc_Type := factor(
      Assoc_Type,
      levels = c("Neg_Assoc", "Pos_Assoc"),
      labels = c("Negative", "Positive")
    )]
    assoc_data_long_otu[, Index_Type := "Association Direction Count"]
    assoc_data_long_otu[, Index_Type := as.character(Index_Type)]

    plot_data_combined <- rbind(
      plot_data_long_otu[, .(parsed_genus, Index_Type, Value = Index_Value, Assoc_Type)],
      assoc_data_long_otu[, .(parsed_genus, Index_Type, Value = Count, Assoc_Type)]
    )

    if (hide_axis_title) {
      plot_data_combined[, Index_Type := factor(
        Index_Type,
        levels = c("Association Direction Count", "STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
      )]
    } else {
      plot_data_combined[, Index_Type := factor(
        Index_Type,
        levels = c("Association Direction Count", "Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)")
      )]
    }

    plot_data_combined[, parsed_genus := factor(parsed_genus, levels = final_order)]
    n_rows_plot <- length(levels(plot_data_combined$parsed_genus))
    n_rows <- if (fig_mode == "compact") {
      indices_otu_compact_dim_row_count(
        data.table::as.data.table(as.data.frame(df_for_compact_dims)),
        mapping_dt,
        dummy_top_bottom,
        DUMMY_LABEL
      )
    } else {
      n_rows_plot
    }

    mean_data <- plot_data_combined[Index_Type != "Association Direction Count",
                                   .(Mean_Value = mean(Value, na.rm = TRUE)),
                                   by = .(parsed_genus, Index_Type, Assoc_Type)]

    assoc_max <- ceiling(max(abs(assoc_data_long_otu$Count), na.rm = TRUE) / 5) * 5
    assoc_breaks <- seq(-assoc_max, assoc_max, by = 5)

    p <- ggplot(plot_data_combined, aes(x = Value, y = parsed_genus)) +
      geom_violin(
        data = plot_data_combined[Index_Type != "Association Direction Count"],
        aes(fill = Assoc_Type),
        alpha = 0.65, color = "black", linewidth = 0.1, na.rm = TRUE,
        width = 0.8, scale = "width", position = "identity"
      ) +
      geom_point(
        data = mean_data,
        aes(x = Mean_Value, y = parsed_genus, colour = Assoc_Type, fill = Assoc_Type),
        shape = geom_point_shape, size = geom_point_size / 2, alpha = 0.8, stroke = 0.1, na.rm = TRUE
      ) +
      geom_bar(
        data = plot_data_combined[Index_Type == "Association Direction Count"],
        aes(fill = Assoc_Type), stat = "identity", width = 0.8, alpha = 0.65
      ) +
      geom_text(
        data = plot_data_combined[Index_Type == "Association Direction Count" & Value != 0],
        aes(
          x = ifelse(Value < 0, Value - max(abs(Value)) * 0.02, Value + max(abs(Value)) * 0.02),
          label = abs(Value)
        ),
        size = 1.5, color = "#555555",
        hjust = ifelse(plot_data_combined[Index_Type == "Association Direction Count" & Value != 0]$Value < 0, 1, 0)
      ) +
      scale_fill_manual(values = association_colors, na.translate = FALSE) +
      scale_color_manual(values = association_colors, na.translate = FALSE) +
      geom_vline(
        aes(xintercept = ifelse(Index_Type == "Association Direction Count", 0, 0.5)),
        linetype = "dashed", color = "black", linewidth = 0.2
      ) +
      geom_blank(data = data.frame(
        Index_Type = if (hide_axis_title) {
          factor(
            c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
            levels = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
          )
        } else {
          factor(
            c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)"),
            levels = c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)")
          )
        },
        Value = rep(c(0, 1), 2),
        y_var = rep(levels(plot_data_combined$parsed_genus)[1], 8)
      ), aes(x = Value, y = y_var), inherit.aes = FALSE) +
      facet_wrap(~Index_Type, nrow = 1, scales = "free_x") +
      scale_x_continuous(
        breaks = function(x) {
          if (length(x) == 0) return(NULL)
          if (max(x, na.rm = TRUE) <= 1) {
            c(0, 0.5, 1)
          } else {
            assoc_breaks
          }
        },
        labels = function(x) {
          if (length(x) == 0) return(NULL)
          x
        }
      ) +
      coord_cartesian(xlim = c(NA, NA)) +
      scale_y_discrete(expand = c(0, 0)) +
      egg::theme_article(base_size = idx_label_pt, base_family = pdf_base_family) +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = idx_y_genus_pt, margin = margin(r = TEXT_MARGIN_PT, unit = "pt"), hjust = 1, lineheight = 0.9),
        axis.text.x = element_text(size = idx_label_pt, margin = margin(t = TEXT_MARGIN_PT, unit = "pt")),
        axis.ticks.y = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.x = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL,
        legend.title = element_blank(),
        legend.text = element_text(size = idx_label_pt),
        strip.text = if (hide_axis_title) element_blank() else element_text(size = idx_label_pt, angle = 90, hjust = 0)
      )

    file_suffix <- "_otu_indices.pdf"

  } else {
    stop("Unsupported plot_type provided to plot_indices")
  }

  dims <- calc_exact_plot_dimensions(n_rows, width_mm, calibrate_y_mm)
  out_w <- dims$width_mm
  out_h <- dims$height_mm
  save_pdf_figure(
    paste0(file_prefix, file_suffix), p,
    width_mm = out_w, height_mm = out_h,
    pointsize = 5
  )
  invisible(p)
}

message("Generating Terry plots for none transformation...")
save_plot_object(plot_configs_none, "none", genus_mapping, ubiome_variable_mapping, phylum_info, phylum_colors)

message("Generating Terry plots for hellinger transformation...")
save_plot_object(plot_configs_hellinger, "hellinger", genus_mapping, ubiome_variable_mapping, phylum_info, phylum_colors)

message("Generating Terry plots for clr transformation...")
save_plot_object(plot_configs_clr, "clr", genus_mapping, ubiome_variable_mapping, phylum_info, phylum_colors)

message("Generating indices plots for none transformation...")
for (config in indices_plot_configs_none) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  for (fig_mode in c("comprehensive", "compact")) {
    message("  Plot: ", config$file_prefix, " [", fig_mode, "]")
    plot_indices(df = config$df,
                 plot_type = config$plot_type,
                 file_prefix = file.path(fig_mode, "indices", "none", config$file_prefix),
                 width_mm = config$width_mm,
                 calibrate_y_mm = config$calibrate_y_mm %||% 0,
                 mapping_data = ubiome_genus_mapping_complete,
                 genus_mapping = genus_mapping,
                 ubiome_variable_mapping = ubiome_variable_mapping,
                 prefix_subsets = config$prefix_subsets,
                 remove_na = config$remove_na %||% FALSE,
                 remove_na_unclassified = config$remove_na_unclassified %||% FALSE,
                 show_legend = config$show_legend %||% TRUE,
                 dummy_top_bottom = config$dummy_top_bottom %||% "split",
                 hide_axis_title = config$hide_axis_title %||% FALSE,
                 geom_point_shape = config$geom_point_shape %||% 21,
                 geom_point_size = config$geom_point_size %||% 3.5,
                 fig_mode = fig_mode)
  }
}

message("Generating indices plots for hellinger transformation...")
for (config in indices_plot_configs_hellinger) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  for (fig_mode in c("comprehensive", "compact")) {
    message("  Plot: ", config$file_prefix, " [", fig_mode, "]")
    plot_indices(df = config$df,
                 plot_type = config$plot_type,
                 file_prefix = file.path(fig_mode, "indices", "hellinger", config$file_prefix),
                 width_mm = config$width_mm,
                 calibrate_y_mm = config$calibrate_y_mm %||% 0,
                 mapping_data = ubiome_genus_mapping_complete,
                 genus_mapping = genus_mapping,
                 ubiome_variable_mapping = ubiome_variable_mapping,
                 prefix_subsets = config$prefix_subsets,
                 remove_na = config$remove_na %||% FALSE,
                 remove_na_unclassified = config$remove_na_unclassified %||% FALSE,
                 show_legend = config$show_legend %||% TRUE,
                 dummy_top_bottom = config$dummy_top_bottom %||% "split",
                 hide_axis_title = config$hide_axis_title %||% FALSE,
                 geom_point_shape = config$geom_point_shape %||% 21,
                 geom_point_size = config$geom_point_size %||% 3.5,
                 fig_mode = fig_mode)
  }
}

message("Generating indices plots for clr transformation...")
for (config in indices_plot_configs_clr) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  for (fig_mode in c("comprehensive", "compact")) {
    message("  Plot: ", config$file_prefix, " [", fig_mode, "]")
    plot_indices(df = config$df,
                 plot_type = config$plot_type,
                 file_prefix = file.path(fig_mode, "indices", "clr", config$file_prefix),
                 width_mm = config$width_mm,
                 calibrate_y_mm = config$calibrate_y_mm %||% 0,
                 mapping_data = ubiome_genus_mapping_complete,
                 genus_mapping = genus_mapping,
                 ubiome_variable_mapping = ubiome_variable_mapping,
                 prefix_subsets = config$prefix_subsets,
                 remove_na = config$remove_na %||% FALSE,
                 remove_na_unclassified = config$remove_na_unclassified %||% FALSE,
                 show_legend = config$show_legend %||% TRUE,
                 dummy_top_bottom = config$dummy_top_bottom %||% "split",
                 hide_axis_title = config$hide_axis_title %||% FALSE,
                 geom_point_shape = config$geom_point_shape %||% 21,
                 geom_point_size = config$geom_point_size %||% 3.5,
                 fig_mode = fig_mode)
  }
}

message("All Terry plots and indices plots have been generated in ", viz_out_path)

