#!/usr/bin/env Rscript
################################################################################
##  MICROBIAL SIGNATURE CORRELATION HEATMAPS - LEGACY SIZED PANELS           ##
##  - Reads pre-computed host-host weighted correlation matrices             ##
##  - Renders one sized heatmap per (analysis, transformation) combo         ##
##    where cell area scales with FDR tier and clustering uses 1-|r|         ##
################################################################################

# Environment: R >= 4.5 with ggplot2, dplyr, tidyr, data.table, ComplexHeatmap,
# circlize, RColorBrewer, grid, DBI, RSQLite, extrafont.
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
})

# Load fonts
loadfonts(device = "pdf", quiet = TRUE)

cat("=== MICROBIAL SIGNATURE CORRELATION HEATMAPS ===\n\n")

base_path <- PROJECT_ROOT

# Load variable descriptions from database
cat("Loading variable descriptions from database...\n")
con <- dbConnect(
  RSQLite::SQLite(),
  file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
)

variable_description <- tbl(con, "variable_names_epcf") %>% 
  collect() %>% 
  data.table::as.data.table()

# Create variable mapping
ubiome_variable_mapping <- variable_description %>%
  dplyr::select(Variable.Name, Variable.Description) %>%
  distinct(Variable.Name, .keep_all = TRUE) %>%
  dplyr::rename(var_name = Variable.Name, var_description = Variable.Description) %>%
  dplyr::select(var_name, var_description) %>%
  distinct(var_name, .keep_all = TRUE)
setDT(ubiome_variable_mapping)

cat("Variable descriptions loaded\n")

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

microbial_signature_input <- file.path(intermediate_path, "microbial_signature")

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
transformations <- c("clr", "lognorm")

# Color scheme for heatmaps - EXACT match to original
orbl_colors <- c("#2B5B8A", "#4071A0", "#5789B6", "#6FA3CB", "#A2BCCF", 
                 "#D8D4C9", "#F0AC72", "#EF8530", "#DA6524", "#BE4E21", "#9E3D21")

cat("Configuration loaded\n")

###############################################################################
## 2. DATA LOADING FUNCTIONS                                                 ##
###############################################################################

# Load pre-computed correlation matrices from RDS files
load_correlation_matrix <- function(analysis_type, transformation) {
  dataset_name <- paste0(analysis_type, "_", transformation, "_all_results")
  
  # Locate correlation matrices strictly within the new intermediate outputs
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

# Function to save PDF figures with proper settings
save_pdf_figure <- function(filename, plot_object, width_mm, height_mm) {
  full_path <- file.path(output_path, filename)
  
  # Force create directory if it doesn't exist
  dir.create(dirname(full_path), showWarnings = FALSE, recursive = TRUE)
  
  tryCatch({
    pdf(full_path, 
        width = width_mm / 25.4, 
        height = height_mm / 25.4,
        family = "ArialMT", 
        useDingbats = FALSE)
    draw(plot_object)
    dev.off()
    message("  Saved: ", filename, " (", width_mm, "×", height_mm, "mm)")
    return(TRUE)
  }, error = function(e) {
    message("  Error saving PDF: ", e$message)
    return(FALSE)
  })
}

# Function to get variable description
get_variable_description <- function(var_names) {
  if (exists("ubiome_variable_mapping", envir = .GlobalEnv)) {
    variable_mapping <- ubiome_variable_mapping
    descriptions <- variable_mapping$var_description[match(var_names, variable_mapping$var_name)]
    descriptions[is.na(descriptions)] <- var_names[is.na(descriptions)]
    return(descriptions)
  } else {
    return(var_names)
  }
}

cat("Data loading functions defined\n")

###############################################################################
## 3. HEATMAP CREATION FUNCTIONS                                            ##
###############################################################################

# Create conventional heatmap with sized squares based on FDR - EXACT match to original
create_conventional_heatmap_sized <- function(cor_result, title, filename, width_mm = 1250, height_mm = 1250) {
  
  cor_matrix <- cor_result$correlation_matrix
  fdr_matrix <- cor_result$fdr_matrix
  n_shared_matrix <- cor_result$n_shared_matrix
  n_eff_matrix <- cor_result$n_effective_matrix
  
  # Clean matrices
  valid_vars <- apply(cor_matrix, 1, function(x) sum(!is.na(x)) > 1)
  if (sum(valid_vars) < 2) return(NULL)
  
  cor_clean <- cor_matrix[valid_vars, valid_vars]
  fdr_clean <- fdr_matrix[valid_vars, valid_vars]
  n_shared_clean <- n_shared_matrix[valid_vars, valid_vars]
  n_eff_clean <- n_eff_matrix[valid_vars, valid_vars]
  
  # Get variable descriptions for labels
  var_names <- rownames(cor_clean)
  var_descriptions <- get_variable_description(var_names)
  
  # Update row and column names with descriptions
  rownames(cor_clean) <- var_descriptions
  colnames(cor_clean) <- var_descriptions
  rownames(fdr_clean) <- var_descriptions
  colnames(fdr_clean) <- var_descriptions
  
  # FDR rectangle sizes as specified
  fdr_sizes <- fdr_clean
  fdr_sizes[is.na(fdr_sizes) | fdr_sizes > 0.05] <- 0.25  # Not significant
  fdr_sizes[fdr_sizes <= 0.001 & fdr_sizes > 0] <- 1.0   # Most significant
  fdr_sizes[fdr_sizes <= 0.01 & fdr_sizes > 0.001] <- 0.85
  fdr_sizes[fdr_sizes <= 0.1 & fdr_sizes > 0.01] <- 0.7
  fdr_sizes[fdr_sizes <= 0.05 & fdr_sizes > 0.1] <- 0.55
  
  # Use corrected color function
  col_fun <- create_orbl_color_function()
  
  # Custom cell function to draw sized rectangles with proper coloring
  cell_fun <- function(j, i, x, y, width, height, fill) {
    if (!is.na(cor_clean[i, j]) && i != j) {
      # Get FDR size factor for this cell
      size_factor <- fdr_sizes[i, j]
      
      # Calculate rectangle size (centered in cell)
      rect_width <- width * size_factor
      rect_height <- height * size_factor
      
      # Get the correlation value and map to color
      corr_val <- cor_clean[i, j]
      cell_color <- col_fun(corr_val)
      
      # Draw sized rectangle with correlation-based color
      grid.rect(x, y, 
               width = rect_width, 
               height = rect_height,
               gp = gpar(fill = cell_color, col = "white", lwd = 0.5))
    }
  }
  
  # Create white background matrix for original heatmap
  white_matrix <- cor_clean
  white_matrix[!is.na(white_matrix)] <- 0  # Set all values to 0 for white background
  
  # Create heatmap with white background and custom cell function for colored rectangles
  ht <- Heatmap(
    white_matrix,  # Use white background matrix
    name = "Correlation",
    col = colorRamp2(c(-1, 0, 1), c("white", "white", "white")),  # Force white background
    na_col = "grey90",
    
    # Clustering - always enabled with enhanced dendrograms
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = function(x) as.dist(1 - abs(cor_clean)),  # Use actual correlation matrix for clustering
    clustering_distance_columns = function(x) as.dist(1 - abs(cor_clean)),
    
    # Enhanced dendrogram display and sizing
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    row_dend_gp = gpar(lwd = 1),
    column_dend_gp = gpar(lwd = 1),
    
    # Fixed cell size - 5pt
    width = unit(5 * ncol(cor_clean), "pt"),
    height = unit(5 * nrow(cor_clean), "pt"),
    
    # Custom cell function for sized squares with proper colors
    cell_fun = cell_fun,
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_gap = unit(0.5, "pt"),
    column_gap = unit(0.5, "pt"),
    
    # Labels - all 6pt font with Arial/Helvetica
    row_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    column_names_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
    column_names_rot = 90,
    
    # Legend - use the actual correlation color scale
    heatmap_legend_param = list(
      title = "Weighted\nCorrelation",
      title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "ArialMT"),
      labels_gp = gpar(fontsize = 6, fontfamily = "ArialMT"),
      legend_height = unit(3, "cm"),
      at = c(-1, -0.5, 0, 0.5, 1),
      labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
      col = col_fun  # Use the corrected color function for legend
    ),
    
    # Title
    column_title = title,
    column_title_gp = gpar(fontsize = 6, fontface = "bold", fontfamily = "ArialMT")
  )
  
  # Save heatmap using the original function
  save_pdf_figure(filename, ht, width_mm, height_mm)
  
  return(ht)
}

# Create standalone legend PDF
create_legend_pdf <- function() {
  # Create color function
  col_fun <- create_orbl_color_function()
  
  # Create a simple heatmap just for the legend
  legend_matrix <- matrix(seq(-1, 1, length.out = 100), nrow = 1, ncol = 100)
  rownames(legend_matrix) <- "Legend"
  colnames(legend_matrix) <- paste0("V", 1:100)
  
  # Create heatmap with legend
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
  
  # Save legend PDF
  save_pdf_figure("0_weighted_correlation_legend.pdf", ht_legend, 200, 100)
  
  return(TRUE)
}

cat("Heatmap creation functions defined\n")

###############################################################################
## 5. MAIN PROCESSING FUNCTION                                              ##
###############################################################################

# Process a single dataset
process_single_dataset <- function(analysis_type, transformation, dataset_name) {
  
  message("Processing: ", dataset_name)
  
  # Load pre-computed correlation matrix
  cor_result <- load_correlation_matrix(analysis_type, transformation)
  if (is.null(cor_result)) {
    message("  Skipping: No correlation matrix found")
    return(NULL)
  }
  
  # Create heatmap
  message("  Creating heatmap...")
  tryCatch({
    create_conventional_heatmap_sized(cor_result,
                                     paste("Host Variable Correlations:", analysis_type, transformation),
                                     paste0(dataset_name, "_conventional_heatmap_sized.pdf"))
  }, error = function(e) {
    message("  Error creating heatmap: ", e$message)
  })
  
  return(cor_result)
}

cat("Main processing function defined\n")

###############################################################################
## 6. MAIN EXECUTION                                                        ##
###############################################################################

cat("=== STARTING MICROBIAL SIGNATURE CORRELATION ANALYSIS ===\n")

# Create standalone legend PDF
message("Creating weighted correlation legend...")
create_legend_pdf()

# Process individual datasets
all_correlation_results <- list()
dataset_summary <- data.frame(
  Dataset = character(),
  Analysis_Type = character(),
  Transformation = character(),
  N_Host_Variables = integer(),
  N_Correlations = integer(),
  N_Significant = integer(),
  Percent_Significant = numeric(),
  stringsAsFactors = FALSE
)

message("Processing ", length(analysis_types) * length(transformations), " individual datasets...")

# Process all individual datasets
for (analysis_type in analysis_types) {
  for (transformation in transformations) {
    
    dataset_name <- paste0(analysis_type, "_", transformation, "_all_results")
    dataset_key <- paste0(analysis_type, "_", transformation)
    
    tryCatch({
      cor_result <- process_single_dataset(analysis_type, transformation, dataset_name)
      if (!is.null(cor_result)) {
        all_correlation_results[[dataset_key]] <- cor_result

        n_valid <- `%||%`(cor_result$n_valid_correlations, NA_integer_)
        n_sig <- extract_significant_count(cor_result)
        percent_sig <- if (!is.na(n_sig) && !is.na(n_valid) && n_valid > 0) {
          round(n_sig / n_valid * 100, 1)
        } else {
          NA_real_
        }

        dataset_summary <- rbind(dataset_summary, data.frame(
          Dataset = dataset_key,
          Analysis_Type = analysis_type,
          Transformation = transformation,
          N_Host_Variables = `%||%`(cor_result$n_variables, NA_integer_),
          N_Correlations = n_valid,
          N_Significant = n_sig,
          Percent_Significant = percent_sig,
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      message("ERROR processing ", dataset_name, ": ", e$message)
    })
    
    # Force garbage collection
    gc()
  }
}

# Save summary
write.csv(dataset_summary, file.path(output_path, "dataset_summary.csv"), row.names = FALSE)

# Final summary
message("\n=== ANALYSIS COMPLETE ===")
message("Total datasets processed: ", nrow(dataset_summary))
message("Results saved to: ", output_path)

cat("\n[done] Complete\n")
