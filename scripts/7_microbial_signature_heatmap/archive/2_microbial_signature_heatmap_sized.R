#!/usr/bin/env Rscript
################################################################################
##  MICROBIAL SIGNATURE CORRELATION HEATMAPS - STANDALONE SCRIPT            ##
##  - Creates weighted correlation matrices between host variables           ##
##  - Uses microbial association signatures (β values) as weights           ##
##  - Generates sized heatmaps with FDR-based cell sizing                    ##
##  - Processes all 28 schemes: 24 individual + 4 merged datasets          ##
################################################################################

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

# Load variable descriptions from database
cat("Loading variable descriptions from database...\n")
con <- dbConnect(RSQLite::SQLite(), "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho/data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")

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

base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"

# Input paths
results_path <- file.path(base_path, "results")
db_path <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
config_dir_path <- file.path(base_path, "configs")

# Output paths
output_path <- file.path(base_path, "results/analyses_results/7_microbial_signature_heatmap_out")
intermediate_path <- file.path(output_path, "intermediate")

# Create output directories
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_path, recursive = TRUE, showWarnings = FALSE)

# Analysis configuration
analysis_types <- c("demoWAS", "oradWAS", "exWAS", "pheWAS", "outWAS", "zimWAS")
# transformations <- c("none", "clr", "hellinger", "lognorm")  # COMMENTED OUT - all transformations
transformations <- c("clr", "lognorm")  # ACTIVE - only CLR and lognorm
# transformations <- c("none", "hellinger")  # COMMENTED OUT - can uncomment when needed

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
  
  # Try different possible file locations for RDS files
  possible_paths <- c(
    file.path(results_path, "analyses_results", "4_association_phyloseq_analyses_out", "intermediate", "arcive", "New Folder With Items", dataset_name, "correlation_matrix.rds"),
    file.path(results_path, "analyses_results", "4_association_phyloseq_analyses_out", "intermediate", dataset_name, "correlation_matrix.rds"),
    file.path(results_path, "analyses_results", "4_association_phyloseq_analyses_out", "microbial_signature", dataset_name, "correlation_matrix.rds")
  )
  
  for (file_path in possible_paths) {
    if (file.exists(file_path)) {
      cat("  Loading correlation matrix from:", file_path, "\n")
      return(readRDS(file_path))
    }
  }
  
  cat("  WARNING: Could not find correlation matrix for", dataset_name, "\n")
  return(NULL)
}

# Load variable descriptions
load_variable_descriptions <- function() {
  # Try to load from database
  con <- dbConnect(RSQLite::SQLite(), db_path)
  
  tryCatch({
    # Get variable descriptions from database
    var_desc <- tbl(con, 'VARIABLE_DESCRIPTION') %>% collect()
    dbDisconnect(con)
    return(var_desc)
  }, error = function(e) {
    if (exists("con")) dbDisconnect(con)
    cat("  Could not load variable descriptions from database\n")
    return(NULL)
  })
}

# Helper functions from original
create_orbl_color_function <- function() {
  colorRamp2(
    seq(-1, 1, length.out = length(orbl_colors)),
    orbl_colors
  )
}

# Function to get variable description
get_variable_description <- function(var_names) {
  # Simplified version - just return formatted names
  descriptions <- gsub("_", " ", var_names)
  descriptions <- gsub("\\b(\\w)", "\\U\\1", descriptions, perl = TRUE)
  return(descriptions)
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
## 3. CORE ANALYSIS FUNCTIONS                                               ##
###############################################################################

# Create signature matrix from association results
create_signature_matrix <- function(results_data, min_prevalence = 0.00) {
  
  if (is.null(results_data) || nrow(results_data) == 0) {
    return(NULL)
  }
  
  # Filter and prepare data
  signature_matrix <- results_data %>%
    filter(!is.na(.data$estimate), !is.na(.data$statistic), !is.na(.data$std.error), .data$std.error > 0) %>%
    mutate(
      weight = 1 / (.data$std.error^2),
      abs_estimate = abs(.data$estimate),
      abs_statistic = abs(.data$statistic)
    )
  
  # Add p-value processing if available
  if ("p.value.fdr" %in% colnames(signature_matrix)) {
    signature_matrix <- signature_matrix %>%
      mutate(log_p_fdr = -log10(pmax(.data$p.value.fdr, 1e-300)))
  } else if ("p.value" %in% colnames(signature_matrix)) {
    signature_matrix <- signature_matrix %>%
      mutate(log_p_value = -log10(pmax(.data$p.value, 1e-300)))
  } else {
    # Compute p-values from t-statistics if missing
    signature_matrix <- signature_matrix %>%
      mutate(
        computed_p_value = 2 * pt(abs(.data$statistic), df = Inf, lower.tail = FALSE),
        log_p_value = -log10(pmax(.data$computed_p_value, 1e-300))
      )
  }
  
  # Rename columns for consistency
  if ("term" %in% colnames(signature_matrix)) {
    signature_matrix <- signature_matrix %>%
      rename(host_var = .data$term)
  }
  
  if ("otu" %in% colnames(signature_matrix)) {
    signature_matrix <- signature_matrix %>%
      rename(microbe_clean = .data$otu)
  }
  
  message("  Signature matrix: ", nrow(signature_matrix), " associations")
  message("  Host variables: ", n_distinct(signature_matrix$host_var))
  message("  Microbes: ", n_distinct(signature_matrix$microbe_clean))
  
  return(signature_matrix)
}

# Compute weighted correlation between two host variables
compute_weighted_correlation <- function(sig_matrix_dt, host_var1, host_var2) {
  
  # Use data.table for efficient joins
  data1 <- sig_matrix_dt[host_var == host_var1, .(microbe_clean, est1 = estimate, se1 = std.error)]
  data2 <- sig_matrix_dt[host_var == host_var2, .(microbe_clean, est2 = estimate, se2 = std.error)]
  
  # Set keys for efficient merging
  setkey(data1, microbe_clean)
  setkey(data2, microbe_clean)
  
  # Join on shared microbes
  shared_data <- merge(data1, data2, by = "microbe_clean", all = FALSE, allow.cartesian = FALSE)
  n_shared <- nrow(shared_data)
  
  if (n_shared < 3) {
    return(list(correlation = NA, p_value = NA, n_shared = n_shared, n_effective = NA))
  }
  
  # Handle duplicates if present
  if (any(duplicated(shared_data$microbe_clean))) {
    shared_data <- shared_data[!duplicated(microbe_clean)]
    n_shared <- nrow(shared_data)
  }
  
  if (n_shared < 3) {
    return(list(correlation = NA, p_value = NA, n_shared = n_shared, n_effective = NA))
  }
  
  # Meta-analytic weighting: inverse-variance weighting
  shared_data[, combined_weight := 1 / (se1^2 + se2^2)]
  
  # Remove zero or infinite weights
  shared_data <- shared_data[is.finite(combined_weight) & combined_weight > 0]
  
  if (nrow(shared_data) < 3) {
    return(list(correlation = NA, p_value = NA, n_shared = n_shared, n_effective = NA))
  }
  
  # Extract vectors
  x <- shared_data$est1
  y <- shared_data$est2
  w <- shared_data$combined_weight
  
  # Compute weighted correlation
  sum_w <- sum(w)
  x_mean <- sum(w * x) / sum_w
  y_mean <- sum(w * y) / sum_w
  
  numerator <- sum(w * (x - x_mean) * (y - y_mean))
  x_var <- sum(w * (x - x_mean)^2)
  y_var <- sum(w * (y - y_mean)^2)
  denominator <- sqrt(x_var * y_var)
  
  if (denominator == 0 || is.na(denominator)) {
    return(list(correlation = NA, p_value = NA, n_shared = n_shared, n_effective = NA))
  }
  
  correlation <- numerator / denominator
  
  # Effective sample size: (Σw)² / Σw²
  n_effective <- sum_w^2 / sum(w^2)
  
  # Compute p-value
  if (n_effective > 2 && abs(correlation) < 0.999) {
    t_stat <- correlation * sqrt((n_effective - 2) / (1 - correlation^2))
    p_value <- 2 * pt(abs(t_stat), df = n_effective - 2, lower.tail = FALSE)
  } else {
    p_value <- ifelse(abs(correlation) > 0.999, 0, 1)
  }
  
  return(list(
    correlation = correlation,
    p_value = p_value,
    n_shared = n_shared,
    n_effective = n_effective
  ))
}

# Compute correlation matrix for all host variables
compute_host_correlation_matrix <- function(sig_matrix, min_shared_microbes = 3) {
  
  # Convert to data.table for efficiency
  sig_matrix_dt <- data.table::as.data.table(sig_matrix)
  setkey(sig_matrix_dt, host_var, microbe_clean)
  
  host_vars <- unique(sig_matrix$host_var)
  n_vars <- length(host_vars)
  
  message("  Computing ", n_vars, " × ", n_vars, " correlation matrix...")
  
  # Generate all pairs
  var_pairs <- combn(host_vars, 2, simplify = FALSE)
  message("  Processing ", length(var_pairs), " host variable pairs...")
  
  # Process all pairs
  pair_results <- lapply(var_pairs, function(pair) {
    result <- compute_weighted_correlation(sig_matrix_dt, pair[1], pair[2])
    result$var1 <- pair[1]
    result$var2 <- pair[2]
    return(result)
  })
  
  # Initialize matrices
  cor_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  pval_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  n_shared_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  n_eff_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  
  # Fill diagonal (perfect correlation)
  diag(cor_matrix) <- 1.0
  diag(pval_matrix) <- 0
  for (i in seq_len(n_vars)) {
    n_shared_matrix[i, i] <- sum(sig_matrix_dt$host_var == host_vars[i])
    n_eff_matrix[i, i] <- n_shared_matrix[i, i]
  }
  
  # Fill matrices from pair results
  for (result in pair_results) {
    if (result$n_shared >= min_shared_microbes && !is.na(result$correlation)) {
      i <- which(host_vars == result$var1)
      j <- which(host_vars == result$var2)
      
      cor_matrix[i, j] <- cor_matrix[j, i] <- result$correlation
      pval_matrix[i, j] <- pval_matrix[j, i] <- result$p_value
      n_eff_matrix[i, j] <- n_eff_matrix[j, i] <- result$n_effective
    }
    
    # Always fill n_shared
    i <- which(host_vars == result$var1)
    j <- which(host_vars == result$var2)
    n_shared_matrix[i, j] <- n_shared_matrix[j, i] <- result$n_shared
  }
  
  # Apply FDR correction
  message("  Applying FDR correction...")
  
  # Extract p-values from pair_results
  p_vec <- vapply(pair_results, function(x) {
    if (!is.na(x$correlation) && x$n_shared >= min_shared_microbes) x$p_value else NA_real_
  }, numeric(1))
  
  # Apply FDR correction
  p_clean <- p_vec[!is.na(p_vec)]
  if (length(p_clean) > 0) {
    fdr_corrected <- p.adjust(p_clean, method = "fdr")
    
    # Create FDR matrix
    fdr_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
    diag(fdr_matrix) <- 0
    
    # Fill FDR matrix
    fdr_idx <- 1
    for (k in seq_along(pair_results)) {
      result <- pair_results[[k]]
      if (!is.na(p_vec[k])) {
        i <- which(host_vars == result$var1)
        j <- which(host_vars == result$var2)
        fdr_matrix[i, j] <- fdr_matrix[j, i] <- fdr_corrected[fdr_idx]
        fdr_idx <- fdr_idx + 1
      }
    }
  } else {
    fdr_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  }
  
  # Calculate summary statistics
  n_valid_correlations <- sum(!is.na(cor_matrix) & cor_matrix != 1)
  n_significant_correlations <- sum(fdr_matrix < 0.05, na.rm = TRUE)
  
  return(list(
    correlation_matrix = cor_matrix,
    pvalue_matrix = pval_matrix,
    fdr_matrix = fdr_matrix,
    n_shared_matrix = n_shared_matrix,
    n_effective_matrix = n_eff_matrix,
    n_variables = n_vars,
    n_valid_correlations = n_valid_correlations,
    n_significant_correlations = n_significant_correlations
  ))
}

cat("Core analysis functions defined\n")

###############################################################################
## 4. HEATMAP CREATION FUNCTIONS                                            ##
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

# Load variable descriptions
var_descriptions <- load_variable_descriptions()

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
        
        # Add to summary
        dataset_summary <- rbind(dataset_summary, data.frame(
          Dataset = dataset_key,
          Analysis_Type = analysis_type,
          Transformation = transformation,
          N_Host_Variables = cor_result$n_variables,
          N_Correlations = cor_result$n_valid_correlations,
          N_Significant = cor_result$n_significant_correlations,
          Percent_Significant = round(cor_result$n_significant_correlations / cor_result$n_valid_correlations * 100, 1),
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

cat("\n✓ Complete!\n")
