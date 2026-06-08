#!/usr/bin/env Rscript
################################################################################
##  MICROBIAL SIGNATURE CORRELATION HEATMAPS - OPTIMAL CLUSTER SPLITTING      ##
##  - Renders host-host correlation heatmaps with silhouette-selected k       ##
##    splits for both rows and columns                                        ##
##  - Reads pre-computed correlation matrices + significant-results bundle    ##
##  - CLI flags: otu_threshold, fdr_threshold, prevalence_threshold, k_min,   ##
##    k_max                                                                   ##
################################################################################

# Environment: R >= 4.5 with ggplot2, dplyr, tidyr, data.table, ComplexHeatmap,
# circlize, RColorBrewer, grid, DBI, RSQLite, extrafont, cluster.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(DBI)
  library(RSQLite)
  library(extrafont)
  library(cluster)
})

# Load fonts
loadfonts(device = "pdf", quiet = TRUE)

cat("=== MICROBIAL SIGNATURE HEATMAPS (OPTIMAL k SPLITS) ===\n\n")

base_path <- PROJECT_ROOT

# Load variable descriptions from database (optional)
cat("Loading variable descriptions from database...\n")
db_path_check <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")

if (file.exists(db_path_check)) {
  tryCatch({
    con <- dbConnect(RSQLite::SQLite(), db_path_check)
    variable_description <- tbl(con, "variable_names_epcf") %>%
      collect() %>%
      data.table::as.data.table()
    dbDisconnect(con)
    
    # Create variable mapping
    ubiome_variable_mapping <- variable_description %>%
      dplyr::select(Variable.Name, Variable.Description) %>%
      distinct(Variable.Name, .keep_all = TRUE) %>%
      dplyr::rename(var_name = Variable.Name, var_description = Variable.Description) %>%
      dplyr::select(var_name, var_description) %>%
      distinct(var_name, .keep_all = TRUE)
    setDT(ubiome_variable_mapping)
    
    cat("Variable descriptions loaded\n")
  }, error = function(e) {
    cat("Warning: Could not load variable descriptions from database:", e$message, "\n")
    cat("Will use variable names directly.\n")
    ubiome_variable_mapping <- NULL
  })
} else {
  cat("Warning: Database not found at", db_path_check, "\n")
  cat("Will use variable names directly.\n")
  ubiome_variable_mapping <- NULL
}

###############################################################################
## 1. CONFIGURATION                                                          ##
###############################################################################

# Input paths
db_path <- file.path(
  base_path,
  "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"
)

# Output paths
output_path <- file.path(base_path, "results/analyses_results/7_microbial_signature_heatmap_out")
intermediate_path <- file.path(output_path, "intermediate")

# Create output directories
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_path, recursive = TRUE, showWarnings = FALSE)

# Input intermediates (written by 1_correlation_signature_matrix_intermediate_processing.R)
input_root <- intermediate_path
microbial_signature_input <- file.path(input_root, "microbial_signature")
heatmap_data_input <- file.path(microbial_signature_input, "heatmap_data")

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

extract_significant_count <- function(cor_result) {
  `%||%`(
    cor_result$n_significant_correlations,
    `%||%`(
      cor_result$n_significant_correlations_005,
      `%||%`(
        cor_result$n_significant,
        NA_integer_
      )
    )
  )
}

# Analysis configuration
analysis_types <- c("demoWAS", "oradWAS", "exWAS", "pheWAS", "outWAS")
transformations <- c("clr")

# Color scheme for heatmaps - EXACT match to original
orbl_colors <- c(
  "#2B5B8A", "#4071A0", "#5789B6", "#6FA3CB", "#A2BCCF",
  "#D8D4C9", "#F0AC72", "#EF8530", "#DA6524", "#BE4E21", "#9E3D21"
)

# Candidate k range (silhouette selection)
k_min <- 2L
k_max <- 20L
otu_threshold <- 2L
fdr_threshold <- 0.01
prevalence_threshold <- 0.01

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  for (arg in args) {
    if (grepl("^k_min=", arg)) {
      value <- sub("^k_min=", "", arg)
      suppressWarnings(parsed <- as.integer(value))
      if (!is.na(parsed) && parsed >= 2) {
        k_min <- parsed
      } else {
        warning("Ignoring k_min argument '", arg, "'; expected integer ≥ 2.")
      }
    } else if (grepl("^k_max=", arg)) {
      value <- sub("^k_max=", "", arg)
      suppressWarnings(parsed <- as.integer(value))
      if (!is.na(parsed) && parsed >= 2) {
        k_max <- parsed
      } else {
        warning("Ignoring k_max argument '", arg, "'; expected integer ≥ 2.")
      }
    } else if (grepl("^otu_threshold=", arg)) {
      value <- sub("^otu_threshold=", "", arg)
      suppressWarnings(parsed <- as.integer(value))
      if (!is.na(parsed) && parsed >= 2) {
        otu_threshold <- parsed
      } else {
        warning("Ignoring otu_threshold argument '", arg, "'; expected integer ≥ 2.")
      }
    } else if (grepl("^fdr_threshold=", arg)) {
      value <- sub("^fdr_threshold=", "", arg)
      suppressWarnings(parsed <- as.numeric(value))
      if (!is.na(parsed) && parsed > 0 && parsed < 1) {
        fdr_threshold <- parsed
      } else {
        warning("Ignoring fdr_threshold argument '", arg, "'; expected numeric between 0 and 1.")
      }
    } else if (grepl("^prevalence_threshold=", arg)) {
      value <- sub("^prevalence_threshold=", "", arg)
      suppressWarnings(parsed <- as.numeric(value))
      if (!is.na(parsed) && parsed > 0 && parsed < 1) {
        prevalence_threshold <- parsed
      } else {
        warning("Ignoring prevalence_threshold argument '", arg, "'; expected numeric between 0 and 1.")
      }
    }
  }
}

if (k_max < k_min) {
  warning("k_max < k_min; resetting to k_min = 2, k_max = 20.")
  k_min <- 2L
  k_max <- 20L
}

candidate_k_values <- seq.int(k_min, k_max)
message("Candidate k range: ", k_min, ":", k_max)
message("OTU threshold (distinct OTUs per host variable) set to: ", otu_threshold)
message("FDR threshold set to: ", fdr_threshold)
message("Prevalence threshold (descriptive) set to: ", prevalence_threshold)

heatmap_targets <- expand.grid(
  analysis = c("demoWAS", "oradWAS", "exWAS", "pheWAS", "outWAS"),
  transformation = c("clr"),
  stringsAsFactors = FALSE
) %>%
  mutate(dataset_name = paste0(analysis, "_", transformation, "_all_results"))
heatmap_targets <- bind_rows(
  heatmap_targets,
  tibble(analysis = "allWAS", transformation = "clr", dataset_name = "allWAS_clr_all_results")
)

# ---------------------------------------------------------------------------
# Load significant association results to enforce ≥ OTU rule with FDR filter
# ---------------------------------------------------------------------------
sigres_path <- file.path(input_root, "heatmap_significant_results.rds")

allowed_terms_map <- list()

if (file.exists(sigres_path)) {
  message("Loading significant association results (for ≥2 OTU filter)...")
  significant_results <- readRDS(sigres_path)

  get_host_variable_column <- function(analysis) {
    if (analysis %in% c("demoWAS", "exWAS")) {
      "independent_var"
    } else if (analysis %in% c("oradWAS", "pheWAS", "outWAS")) {
      "dependent_var"
    } else {
      "term"
    }
  }

  for (name in names(significant_results)) {
    df <- significant_results[[name]]
    if (is.null(df) || !nrow(df)) {
      next
    }
    analysis <- sub("_.*$", "", name)
    df_filtered <- df %>%
      filter(!is.na(p.value.fdr), p.value.fdr < fdr_threshold)

    if (!nrow(df_filtered)) {
      next
    }

    host_col <- get_host_variable_column(analysis)
    if (!host_col %in% names(df_filtered)) {
      host_col <- if ("term" %in% names(df_filtered)) "term" else host_col
    }
    df_filtered$host_variable <- df_filtered[[host_col]]

    keep_terms <- df_filtered %>%
      filter(!is.na(host_variable), host_variable != "") %>%
      group_by(host_variable) %>%
      summarise(n_otus = n_distinct(otu), .groups = "drop") %>%
      filter(n_otus >= otu_threshold) %>%
      pull(host_variable)

    if (length(keep_terms) == 0) {
      next
    }

    dataset_name <- paste0(sub("_sig_res$", "", name), "_all_results")
    if (dataset_name %in% heatmap_targets$dataset_name) {
      allowed_terms_map[[dataset_name]] <- unique(keep_terms)
    }
  }

  # Build aggregated allWAS set from component analyses
  allwas_sources <- c(
    "oradWAS_clr_all_results",
    "outWAS_clr_all_results",
    "pheWAS_clr_all_results",
    "exWAS_clr_all_results"
  )
  available_sources <- intersect(allwas_sources, names(allowed_terms_map))
  if (length(available_sources) > 0) {
    all_terms <- unique(unlist(allowed_terms_map[available_sources], use.names = FALSE))
    if (length(all_terms) >= 2) {
      allowed_terms_map[["allWAS_clr_all_results"]] <- all_terms
    }
  }

  message("  Filterable datasets: ", length(allowed_terms_map))
} else {
  warning(
    "Significant results file not found at ",
    sigres_path,
    ". No heatmaps will be generated."
  )
}

format_otu_label <- function(x) as.character(as.integer(round(x)))
format_fdr_label <- function(x) as.character(as.integer(round(x * 100)))
format_prev_label <- function(x) as.character(as.integer(round(x * 1000)))

combo_prefix <- paste0(
  "otu", format_otu_label(otu_threshold),
  "_fdr", format_fdr_label(fdr_threshold),
  "_prev", format_prev_label(prevalence_threshold)
)

cat("Configuration loaded\n")

###############################################################################
## 2. DATA LOADING FUNCTIONS                                                 ##
###############################################################################

compute_heatmap_dimensions <- function(n_rows,
                                       min_dim = 300,
                                       base_dim = 80,
                                       per_row = 2.4,
                                       max_dim = 1250) {
  if (is.na(n_rows) || n_rows <= 0) n_rows <- 1
  dim_estimate <- base_dim + n_rows * per_row
  dim_limited <- max(min_dim, min(max_dim, dim_estimate))
  list(width_mm = dim_limited, height_mm = dim_limited)
}

# Load pre-computed correlation matrices from RDS files
load_correlation_matrix <- function(analysis_type, transformation) {
  dataset_name <- paste0(analysis_type, "_", transformation, "_all_results")

  possible_paths <- if (dataset_name == "allWAS_clr_all_results") {
    c(file.path(microbial_signature_input, "merged_clr_data", "correlation_matrix.rds"))
  } else {
    c(file.path(microbial_signature_input, dataset_name, "correlation_matrix.rds"))
  }

  for (file_path in possible_paths) {
    if (file.exists(file_path)) {
      cat("  Loading correlation matrix from:", file_path, "\n")
      return(readRDS(file_path))
    }
  }

  cat(
    "  WARNING: Could not find correlation matrix for",
    dataset_name,
    "under",
    microbial_signature_input,
    "\n"
  )
  return(NULL)
}

# Helper functions from original
create_orbl_color_function <- function() {
  colorRamp2(
    seq(-1, 1, length.out = length(orbl_colors)),
    orbl_colors
  )
}

get_variable_description <- function(var_names) {
  if (exists("ubiome_variable_mapping", envir = .GlobalEnv)) {
    variable_mapping <- ubiome_variable_mapping
    descriptions <- variable_mapping$var_description[match(var_names, variable_mapping$var_name)]
    descriptions[is.na(descriptions)] <- var_names[is.na(descriptions)]
    return(descriptions)
  }
  var_names
}

save_pdf_figure <- function(filename, plot_object, width_mm, height_mm) {
  full_path <- file.path(output_path, filename)
  dir.create(dirname(full_path), showWarnings = FALSE, recursive = TRUE)

  tryCatch({
    pdf(
      full_path,
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      family = "ArialMT",
      useDingbats = FALSE
    )
    if (inherits(plot_object, "Heatmap")) {
      draw(plot_object)
    } else if (inherits(plot_object, "ggplot")) {
      grid::grid.newpage()
      grid::grid.draw(ggplot2::ggplotGrob(plot_object))
    } else {
      grid::grid.newpage()
      grid::grid.draw(plot_object)
    }
    dev.off()
    message("  Saved: ", filename, " (", width_mm, "×", height_mm, "mm)")
    TRUE
  }, error = function(e) {
    message("  Error saving PDF: ", e$message)
    FALSE
  })
}

save_placeholder_heatmap <- function(analysis_type, transformation, combo_prefix,
                                     reason, n_rows = 1) {
  placeholder_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      geom = "text",
      x = 0.5,
      y = 0.5,
      label = reason,
      family = "ArialMT",
      size = 3,
      lineheight = 1.1,
      color = "#333333"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()

  file_name <- compose_heatmap_filename(
    row_k = NA,
    column_k = NA,
    analysis_type = analysis_type,
    transformation = transformation,
    prefix = combo_prefix
  )

  dims <- compute_heatmap_dimensions(max(1, n_rows))
  save_pdf_figure(file_name, placeholder_plot, dims$width_mm, dims$height_mm)
}

###############################################################################
## 3. OPTIMAL CLUSTER SELECTION                                              ##
###############################################################################

format_k_label <- function(row_k, column_k) {
  if (is.na(row_k) && is.na(column_k)) {
    "na"
  } else if (!is.na(row_k) && !is.na(column_k) && row_k != column_k) {
    paste0(row_k, "c", column_k)
  } else {
    as.character(ifelse(is.na(row_k), column_k, row_k))
  }
}

compose_heatmap_filename <- function(row_k, column_k, analysis_type,
                                     transformation, prefix = combo_prefix) {
  paste0(
    prefix,
    "_k", format_k_label(row_k, column_k), "_",
    analysis_type, "_", transformation,
    "_siges_heatmap.pdf"
  )
}

select_optimal_k <- function(hclust_obj, dist_mat, k_values) {
  n <- length(hclust_obj$order)
  k_candidates <- k_values[k_values >= 2 & k_values < n]
  if (length(k_candidates) == 0) {
    return(list(k = NA_integer_, silhouette = NA_real_))
  }

  scores <- sapply(k_candidates, function(k) {
    clusters <- cutree(hclust_obj, k = k)
    if (length(unique(clusters)) < 2) {
      return(NA_real_)
    }
    sil <- tryCatch(cluster::silhouette(clusters, dist_mat), error = function(e) NULL)
    if (is.null(sil)) {
      return(NA_real_)
    }
    mean(sil[, "sil_width"], na.rm = TRUE)
  })

  if (all(is.na(scores))) {
    return(list(k = NA_integer_, silhouette = NA_real_))
  }

  best_idx <- which.max(scores)
  list(k = k_candidates[best_idx], silhouette = scores[best_idx])
}

###############################################################################
## 4. HEATMAP CREATION FUNCTIONS                                            ##
###############################################################################

create_conventional_heatmap_sized <- function(cor_result, dataset_name,
                                               analysis_type, transformation,
                                               combo_prefix = combo_prefix) {
  title <- paste("Host Variable Correlations:", analysis_type, transformation)
  cor_matrix <- cor_result$correlation_matrix
  fdr_matrix <- cor_result$fdr_matrix
  n_shared_matrix <- cor_result$n_shared_matrix
  n_eff_matrix <- cor_result$n_effective_matrix

  valid_vars <- apply(cor_matrix, 1, function(x) sum(!is.na(x)) > 0)
  if (sum(valid_vars) < 2) {
    placeholder_reason <- sprintf(
      "%s\nCorrelation matrix has fewer than %d non-NA variables\nunder current thresholds.",
      dataset_name,
      2
    )
    save_placeholder_heatmap(
      analysis_type,
      transformation,
      combo_prefix,
      placeholder_reason,
      n_rows = sum(valid_vars)
    )
    return(list(
      row_k = NA_integer_,
      row_silhouette = NA_real_,
      column_k = NA_integer_,
      column_silhouette = NA_real_,
      retained_variables = sum(valid_vars)
    ))
  }

  cor_clean <- cor_matrix[valid_vars, valid_vars]
  fdr_clean <- fdr_matrix[valid_vars, valid_vars]
  n_shared_clean <- n_shared_matrix[valid_vars, valid_vars]
  n_eff_clean <- n_eff_matrix[valid_vars, valid_vars]

  allowed_terms <- allowed_terms_map[[dataset_name]]
  if (is.null(allowed_terms)) {
    message("  No host variables pass ≥", otu_threshold, " OTU (FDR < ", fdr_threshold, ") filter; skipping.")
    placeholder_reason <- sprintf(
      "%s\nNo host variables pass FDR ≤ %.3f\nwith OTU ≥ %d and prevalence ≥ %.2f%%.",
      dataset_name,
      fdr_threshold,
      otu_threshold,
      prevalence_threshold * 100
    )
    save_placeholder_heatmap(
      analysis_type,
      transformation,
      combo_prefix,
      placeholder_reason,
      n_rows = length(allowed_terms)
    )
    return(list(
      row_k = NA_integer_,
      row_silhouette = NA_real_,
      column_k = NA_integer_,
      column_silhouette = NA_real_,
      retained_variables = 0
    ))
  }

  missing_terms <- setdiff(allowed_terms, rownames(cor_clean))
  if (length(missing_terms) > 0) {
    message("  Host variables absent from correlation matrix: ", length(missing_terms),
            " (appending NA-structured rows)")
    desired_terms <- unique(c(rownames(cor_clean), missing_terms))
    expand_square <- function(mat) {
      full <- matrix(NA_real_, nrow = length(desired_terms), ncol = length(desired_terms),
                     dimnames = list(desired_terms, desired_terms))
      if (!is.null(mat) && length(rownames(mat))) {
        existing <- intersect(rownames(mat), desired_terms)
        full[existing, existing] <- mat[existing, existing, drop = FALSE]
      }
      full
    }
    cor_clean <- expand_square(cor_clean)
    fdr_clean <- expand_square(fdr_clean)
    n_shared_clean <- expand_square(n_shared_clean)
    n_eff_clean <- expand_square(n_eff_clean)
    diag(cor_clean)[is.na(diag(cor_clean))] <- 1
    fdr_clean[is.na(fdr_clean)] <- 1
    diag(fdr_clean) <- 0
    n_shared_clean[is.na(n_shared_clean)] <- 0
    n_eff_clean[is.na(n_eff_clean)] <- 0
  }

  message("  Allowed host variables: ", length(allowed_terms))

  keep_idx <- rownames(cor_clean) %in% allowed_terms
  if (sum(keep_idx) < 2) {
    message("  Fewer than two host variables meet ≥", otu_threshold, " OTU (FDR < ", fdr_threshold, ") filter; skipping.")
    placeholder_reason <- sprintf(
      "%s\nFewer than %d host variables remain after applying\nFDR ≤ %.3f, prevalence ≥ %.2f%%, OTU ≥ %d.",
      dataset_name,
      2,
      fdr_threshold,
      prevalence_threshold * 100,
      otu_threshold
    )
    save_placeholder_heatmap(
      analysis_type,
      transformation,
      combo_prefix,
      placeholder_reason,
      n_rows = sum(keep_idx)
    )
    return(list(
      row_k = NA_integer_,
      row_silhouette = NA_real_,
      column_k = NA_integer_,
      column_silhouette = NA_real_,
      retained_variables = sum(keep_idx)
    ))
  }

  cor_clean <- cor_clean[keep_idx, keep_idx, drop = FALSE]
  fdr_clean <- fdr_clean[keep_idx, keep_idx, drop = FALSE]
  n_shared_clean <- n_shared_clean[keep_idx, keep_idx, drop = FALSE]
  n_eff_clean <- n_eff_clean[keep_idx, keep_idx, drop = FALSE]

  message("  Retained host variables after filtering: ", nrow(cor_clean))

  # Remove exact duplicate host-variable profiles (identical correlation rows/columns)
  if (analysis_type != "allWAS" && nrow(cor_clean) > 1) {
    row_signature <- apply(cor_clean, 1, function(x) paste(signif(x, 12), collapse = "|"))
    unique_mask <- !duplicated(row_signature)
    if (any(!unique_mask)) {
      removed <- sum(!unique_mask)
      message(sprintf("  Removing %d duplicated host-variable profile(s).", removed))
      cor_clean <- cor_clean[unique_mask, unique_mask, drop = FALSE]
      fdr_clean <- fdr_clean[unique_mask, unique_mask, drop = FALSE]
      n_shared_clean <- n_shared_clean[unique_mask, unique_mask, drop = FALSE]
      n_eff_clean <- n_eff_clean[unique_mask, unique_mask, drop = FALSE]
    }
  }

  if (nrow(cor_clean) < 2) {
    message("  Less than two unique host-variable profiles remain after duplicate removal; skipping.")
    placeholder_reason <- sprintf(
      "%s\nDuplicate removal leaves fewer than %d unique host profiles\nunder FDR ≤ %.3f, prevalence ≥ %.2f%%, OTU ≥ %d.",
      dataset_name,
      2,
      fdr_threshold,
      prevalence_threshold * 100,
      otu_threshold
    )
    save_placeholder_heatmap(
      analysis_type,
      transformation,
      combo_prefix,
      placeholder_reason,
      n_rows = nrow(cor_clean)
    )
    return(list(
      row_k = NA_integer_,
      row_silhouette = NA_real_,
      column_k = NA_integer_,
      column_silhouette = NA_real_,
      retained_variables = nrow(cor_clean)
    ))
  }

  var_names <- rownames(cor_clean)
  var_descriptions <- get_variable_description(var_names)
  rownames(cor_clean) <- var_descriptions
  colnames(cor_clean) <- var_descriptions
  rownames(fdr_clean) <- var_descriptions
  colnames(fdr_clean) <- var_descriptions

  cor_for_dist <- cor_clean
  cor_for_dist[is.na(cor_for_dist)] <- 0
  diag(cor_for_dist) <- 1

  row_dist <- as.dist(1 - abs(cor_for_dist))
  row_hclust <- hclust(row_dist, method = "complete")

  col_dist <- as.dist(1 - abs(cor_for_dist))
  col_hclust <- hclust(col_dist, method = "complete")

  row_sel <- select_optimal_k(row_hclust, row_dist, candidate_k_values)
  col_sel <- select_optimal_k(col_hclust, col_dist, candidate_k_values)

  row_split_arg <- if (!is.na(row_sel$k)) row_sel$k else NULL
  column_split_arg <- if (!is.na(col_sel$k)) col_sel$k else NULL

  fdr_sizes <- fdr_clean
  fdr_sizes[is.na(fdr_sizes) | fdr_sizes > 0.05] <- 0.25
  fdr_sizes[fdr_sizes <= 0.001 & fdr_sizes > 0] <- 1.0
  fdr_sizes[fdr_sizes <= 0.01 & fdr_sizes > 0.001] <- 0.85
  fdr_sizes[fdr_sizes <= 0.1 & fdr_sizes > 0.01] <- 0.7
  fdr_sizes[fdr_sizes <= 0.05 & fdr_sizes > 0.1] <- 0.55

  col_fun <- create_orbl_color_function()

  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(cor_clean[i, j]) && i != j) {
      size_factor <- fdr_sizes[i, j]
      rect_width <- width * size_factor
      rect_height <- height * size_factor
      corr_val <- cor_clean[i, j]
      cell_color <- col_fun(corr_val)
      grid.rect(
        x, y,
        width = rect_width,
        height = rect_height,
        gp = gpar(fill = cell_color, col = "white", lwd = 0.5)
      )
    }
  }

  white_matrix <- cor_clean
  white_matrix[!is.na(white_matrix)] <- 0

  ht <- Heatmap(
    white_matrix,
    name = "Correlation",
    col = colorRamp2(c(-1, 0, 1), c("white", "white", "white")),
    na_col = "grey90",
    cluster_rows = row_hclust,
    cluster_columns = col_hclust,
    clustering_distance_rows = function(x) as.dist(1 - abs(cor_clean)),
    clustering_distance_columns = function(x) as.dist(1 - abs(cor_clean)),
    row_split = row_split_arg,
    column_split = column_split_arg,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    row_dend_gp = gpar(lwd = 1),
    column_dend_gp = gpar(lwd = 1),
    width = unit(5 * ncol(cor_clean), "pt"),
    height = unit(5 * nrow(cor_clean), "pt"),
    cell_fun = cell_fun,
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_gap = unit(2, "mm"),
    column_gap = unit(2, "mm"),
    row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title = "Weighted\nCorrelation",
      title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "ArialMT"),
      labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
      legend_height = unit(3, "cm"),
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
      col = col_fun
    ),
    column_title = title,
    column_title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "ArialMT")
  )

  file_name <- compose_heatmap_filename(
    row_k = row_sel$k,
    column_k = col_sel$k,
    analysis_type = analysis_type,
    transformation = transformation,
    prefix = combo_prefix
  )

  dims <- compute_heatmap_dimensions(nrow(cor_clean))
  save_pdf_figure(file_name, ht, dims$width_mm, dims$height_mm)

  message(sprintf("    Selected splits -> rows: %s, columns: %s",
                  ifelse(is.na(row_sel$k), "none", row_sel$k),
                  ifelse(is.na(col_sel$k), "none", col_sel$k)))

  list(
    row_k = row_sel$k,
    row_silhouette = row_sel$silhouette,
    column_k = col_sel$k,
    column_silhouette = col_sel$silhouette,
    retained_variables = nrow(cor_clean)
  )
}

create_legend_pdf <- function(prefix_base = combo_prefix) {
  col_fun <- create_orbl_color_function()

  legend_matrix <- matrix(seq(-1, 1, length.out = 100), nrow = 1, ncol = 100)
  rownames(legend_matrix) <- "Legend"
  colnames(legend_matrix) <- paste0("V", 1:100)

  ht_legend <- Heatmap(
    legend_matrix,
    name = "Weighted\nCorrelation",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = "Weighted\nCorrelation",
      title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "ArialMT"),
      labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
      legend_height = unit(3, "cm"),
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
      col = col_fun
    ),
    width = unit(10, "cm"),
    height = unit(2, "cm")
  )

  save_pdf_figure(paste0(prefix_base, "_weighted_correlation_legend.pdf"), ht_legend, 200, 100)
  TRUE
}

cat("Heatmap creation functions defined\n")

###############################################################################
## 5. MAIN PROCESSING FUNCTION                                              ##
###############################################################################

process_single_dataset <- function(analysis_type, transformation, dataset_name,
                                   prefix_base = combo_prefix) {
  message("Processing: ", dataset_name)

  cor_result <- load_correlation_matrix(analysis_type, transformation)
  if (is.null(cor_result)) {
    message("  Skipping: No correlation matrix found")
    return(list(cor_result = NULL))
  }

  message("  Creating heatmap...")
  split_info <- tryCatch({
    create_conventional_heatmap_sized(
      cor_result = cor_result,
      dataset_name = dataset_name,
      analysis_type = analysis_type,
      transformation = transformation,
      combo_prefix = prefix_base
    )
  }, error = function(e) {
    message("  Error creating heatmap: ", e$message)
    list(row_k = NA_integer_, row_silhouette = NA_real_,
         column_k = NA_integer_, column_silhouette = NA_real_,
         retained_variables = 0)
  })

  list(
    cor_result = cor_result,
    split_info = split_info
  )
}

cat("Main processing function defined\n")

###############################################################################
## 6. MAIN EXECUTION                                                        ##
###############################################################################

cat("=== STARTING OPTIMAL-k MICROBIAL SIGNATURE CORRELATION ANALYSIS ===\n")

message("Creating weighted correlation legend...")
create_legend_pdf()

all_correlation_results <- list()
dataset_summary <- data.frame(
  Dataset = character(),
  Analysis_Type = character(),
  Transformation = character(),
  N_Host_Variables = integer(),
  N_Correlations = integer(),
  N_Significant = integer(),
  Percent_Significant = numeric(),
  Row_K = integer(),
  Row_Silhouette = numeric(),
  Column_K = integer(),
  Column_Silhouette = numeric(),
  Retained_Variables = integer(),
  FDR_Threshold = numeric(),
  OTU_Threshold = integer(),
  Prevalence_Threshold = numeric(),
  stringsAsFactors = FALSE
)

message("Processing ", length(analysis_types) * length(transformations), " individual datasets...")

for (analysis_type in analysis_types) {
  for (transformation in transformations) {

    dataset_name <- paste0(analysis_type, "_", transformation, "_all_results")
    dataset_key <- paste0(analysis_type, "_", transformation)

    terms_available <- allowed_terms_map[[dataset_name]]
    if (is.null(terms_available) || length(terms_available) < 2) {
      message("Skipping ", dataset_name, " (fewer than 2 host variables pass ≥", otu_threshold, " OTU filter).")
      placeholder_reason <- sprintf(
        "%s\nFewer than %d host variables pass FDR ≤ %.3f,\nprevalence ≥ %.2f%%, OTU ≥ %d.",
        dataset_name,
        2,
        fdr_threshold,
        prevalence_threshold * 100,
        otu_threshold
      )
      save_placeholder_heatmap(
        analysis_type,
        transformation,
        combo_prefix,
        placeholder_reason,
        n_rows = length(terms_available)
      )
      dataset_summary <- rbind(dataset_summary, data.frame(
        Dataset = dataset_key,
        Analysis_Type = analysis_type,
        Transformation = transformation,
        N_Host_Variables = 0,
        N_Correlations = 0,
        N_Significant = 0,
        Percent_Significant = NA_real_,
        Row_K = NA_integer_,
        Row_Silhouette = NA_real_,
        Column_K = NA_integer_,
        Column_Silhouette = NA_real_,
        Retained_Variables = 0,
        FDR_Threshold = fdr_threshold,
        OTU_Threshold = otu_threshold,
        Prevalence_Threshold = prevalence_threshold,
        stringsAsFactors = FALSE
      ))
      next
    }

    res <- process_single_dataset(analysis_type, transformation, dataset_name)
    if (!is.null(res$cor_result)) {
      all_correlation_results[[dataset_key]] <- res$cor_result
      n_valid <- `%||%`(res$cor_result$n_valid_correlations, NA_integer_)
      n_sig <- extract_significant_count(res$cor_result)
      percent_sig <- if (!is.na(n_sig) && !is.na(n_valid) && n_valid > 0) {
        round(n_sig / n_valid * 100, 1)
      } else {
        NA_real_
      }

      dataset_summary <- rbind(dataset_summary, data.frame(
        Dataset = dataset_key,
        Analysis_Type = analysis_type,
        Transformation = transformation,
        N_Host_Variables = `%||%`(res$cor_result$n_variables, NA_integer_),
        N_Correlations = n_valid,
        N_Significant = n_sig,
        Percent_Significant = percent_sig,
        Row_K = res$split_info$row_k,
        Row_Silhouette = res$split_info$row_silhouette,
        Column_K = res$split_info$column_k,
        Column_Silhouette = res$split_info$column_silhouette,
        Retained_Variables = ifelse(
          is.null(res$split_info$retained_variables),
          NA_integer_,
          res$split_info$retained_variables
        ),
        FDR_Threshold = fdr_threshold,
        OTU_Threshold = otu_threshold,
        Prevalence_Threshold = prevalence_threshold,
        stringsAsFactors = FALSE
      ))
    }

    gc()
  }
}

# Process aggregated allWAS CLR dataset if available
if ("allWAS_clr_all_results" %in% names(allowed_terms_map)) {
  analysis_type <- "allWAS"
  transformation <- "clr"
  dataset_name <- "allWAS_clr_all_results"
  dataset_key <- paste0(analysis_type, "_", transformation)

  terms_allwas <- allowed_terms_map[[dataset_name]]
  if (!is.null(terms_allwas) && length(terms_allwas) >= 2) {
    message("Processing: ", dataset_name)
    res <- process_single_dataset(analysis_type, transformation, dataset_name)
    if (!is.null(res$cor_result)) {
      all_correlation_results[[dataset_key]] <- res$cor_result
      n_valid <- `%||%`(res$cor_result$n_valid_correlations, NA_integer_)
      n_sig <- extract_significant_count(res$cor_result)
      percent_sig <- if (!is.na(n_sig) && !is.na(n_valid) && n_valid > 0) {
        round(n_sig / n_valid * 100, 1)
      } else {
        NA_real_
      }

      dataset_summary <- rbind(dataset_summary, data.frame(
        Dataset = dataset_key,
        Analysis_Type = analysis_type,
        Transformation = transformation,
        N_Host_Variables = `%||%`(res$cor_result$n_variables, NA_integer_),
        N_Correlations = n_valid,
        N_Significant = n_sig,
        Percent_Significant = percent_sig,
        Row_K = res$split_info$row_k,
        Row_Silhouette = res$split_info$row_silhouette,
        Column_K = res$split_info$column_k,
        Column_Silhouette = res$split_info$column_silhouette,
        Retained_Variables = ifelse(
          is.null(res$split_info$retained_variables),
          NA_integer_,
          res$split_info$retained_variables
        ),
        FDR_Threshold = fdr_threshold,
        OTU_Threshold = otu_threshold,
        Prevalence_Threshold = prevalence_threshold,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    message("Skipping ", dataset_name, " (fewer than 2 host variables aggregated across WAS analyses).")
    placeholder_reason <- sprintf(
      "%s\nCombined WAS analyses yield fewer than %d host variables\nwith FDR ≤ %.3f, prevalence ≥ %.2f%%, OTU ≥ %d.",
      dataset_name,
      2,
      fdr_threshold,
      prevalence_threshold * 100,
      otu_threshold
    )
    n_placeholder <- if (exists("terms_allwas") && length(terms_allwas)) length(terms_allwas) else 1
    save_placeholder_heatmap(
      analysis_type,
      transformation,
      combo_prefix,
      placeholder_reason,
      n_rows = n_placeholder
    )
    dataset_summary <- rbind(dataset_summary, data.frame(
      Dataset = dataset_key,
      Analysis_Type = analysis_type,
      Transformation = transformation,
      N_Host_Variables = 0,
      N_Correlations = 0,
      N_Significant = 0,
      Percent_Significant = NA_real_,
      Row_K = NA_integer_,
      Row_Silhouette = NA_real_,
      Column_K = NA_integer_,
      Column_Silhouette = NA_real_,
      Retained_Variables = 0,
      FDR_Threshold = fdr_threshold,
      OTU_Threshold = otu_threshold,
      Prevalence_Threshold = prevalence_threshold,
      stringsAsFactors = FALSE
    ))
  }
}

summary_filename <- paste0(combo_prefix, "_dataset_summary.csv")

write.csv(
  dataset_summary,
  file.path(output_path, summary_filename),
  row.names = FALSE
)

message("\n=== ANALYSIS COMPLETE ===")
message("Total datasets processed: ", nrow(dataset_summary))
message("Results saved to: ", output_path)
message("Summary file: ", summary_filename)

cat("\n[done] Complete\n")


