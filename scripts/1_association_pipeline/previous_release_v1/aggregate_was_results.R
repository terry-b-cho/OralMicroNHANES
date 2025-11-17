#!/usr/bin/env Rscript

# =====================================================================================
# WAS Results Aggregation Script with Multiple Comparisons Correction
# aggregate_was_results.R
# =====================================================================================
# This script aggregates individual dependent variable results into analysis-type
# specific files and applies multiple comparisons correction (FDR and Bonferroni)
#
# Usage: Rscript aggregate_was_results.R <analysis_type> <normalization> <results_dir>
# Example: Rscript aggregate_was_results.R 1_demoWAS clr /path/to/results
# =====================================================================================

library(dplyr)
library(readr)
library(glue)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript aggregate_was_results.R <analysis_type> <normalization> <results_dir>\n")
  cat("Example: Rscript aggregate_was_results.R 1_demoWAS clr /path/to/results\n")
  quit(status = 1)
}

analysis_type <- args[1]
normalization <- args[2]
results_dir <- args[3]

cat("=============================================================================\n")
cat("WAS Results Aggregation with Multiple Comparisons Correction\n")
cat("=============================================================================\n")
cat("Analysis Type:", analysis_type, "\n")
cat("Normalization:", normalization, "\n")
cat("Results Directory:", results_dir, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")

# Define paths
input_dir <- file.path(results_dir, paste0(analysis_type, "_out"), paste0("result_", normalization))
output_dir <- input_dir  # Save aggregated results in the same directory

if (!dir.exists(input_dir)) {
  cat("ERROR: Input directory does not exist:", input_dir, "\n")
  quit(status = 1)
}

# Find all .rds result files
rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
# Exclude model files and aggregated files
rds_files <- rds_files[!grepl("_models\\.rds$", rds_files)]
rds_files <- rds_files[!grepl("_complete\\.rds$", rds_files)]

if (length(rds_files) == 0) {
  cat("ERROR: No individual result files found in:", input_dir, "\n")
  quit(status = 1)
}

cat("Found", length(rds_files), "individual result files\n")

# Function to safely read RDS files
safe_read_rds <- function(file_path) {
  tryCatch({
    result <- readRDS(file_path)
    if (!is.list(result) || !all(c("pe_tidied", "pe_glanced", "rsq") %in% names(result))) {
      cat("WARNING: Invalid result structure in", basename(file_path), "\n")
      return(NULL)
    }
    return(result)
  }, error = function(e) {
    cat("ERROR reading", basename(file_path), ":", e$message, "\n")
    return(NULL)
  })
}

# Function to apply multiple comparisons correction
apply_multiple_comparisons_correction <- function(tidied_data) {
  if (nrow(tidied_data) == 0) {
    tidied_data$p.value.fdr <- numeric(0)
    tidied_data$p.value.bonferroni <- numeric(0)
    return(tidied_data)
  }
  
  # Only adjust p-values for the main effect terms (indep_var), not intercepts or covariates
  main_effect_rows <- tidied_data$term == "indep_var"
  
  if (sum(main_effect_rows) == 0) {
    cat("WARNING: No main effect terms found for p-value adjustment\n")
    tidied_data$p.value.fdr <- tidied_data$p.value
    tidied_data$p.value.bonferroni <- tidied_data$p.value
    return(tidied_data)
  }
  
  # Extract p-values for main effects only
  main_effect_pvals <- tidied_data$p.value[main_effect_rows]
  
  # Apply corrections
  fdr_adjusted <- p.adjust(main_effect_pvals, method = "fdr")
  bonferroni_adjusted <- p.adjust(main_effect_pvals, method = "bonferroni")
  
  # Initialize adjustment columns
  tidied_data$p.value.fdr <- tidied_data$p.value
  tidied_data$p.value.bonferroni <- tidied_data$p.value
  
  # Apply adjustments only to main effect terms
  tidied_data$p.value.fdr[main_effect_rows] <- fdr_adjusted
  tidied_data$p.value.bonferroni[main_effect_rows] <- bonferroni_adjusted
  
  return(tidied_data)
}

# Initialize aggregated data containers
all_tidied <- list()
all_glanced <- list()
all_rsq <- list()

# Process each file
cat("Processing individual result files...\n")
successful_files <- 0

for (i in seq_along(rds_files)) {
  file_path <- rds_files[i]
  file_name <- basename(file_path)
  
  cat(sprintf("[%d/%d] Processing %s...\n", i, length(rds_files), file_name))
  
  result <- safe_read_rds(file_path)
  
  if (is.null(result)) {
    next
  }
  
  # Extract and store each component
  if (!is.null(result$pe_tidied) && nrow(result$pe_tidied) > 0) {
    all_tidied[[file_name]] <- result$pe_tidied
  }
  
  if (!is.null(result$pe_glanced) && nrow(result$pe_glanced) > 0) {
    all_glanced[[file_name]] <- result$pe_glanced
  }
  
  if (!is.null(result$rsq) && nrow(result$rsq) > 0) {
    all_rsq[[file_name]] <- result$rsq
  }
  
  successful_files <- successful_files + 1
}

cat("Successfully processed", successful_files, "files\n\n")

# Combine all results
cat("Aggregating results...\n")

# Combine tidied results
if (length(all_tidied) > 0) {
  combined_tidied <- bind_rows(all_tidied)
  cat("Combined tidied results:", nrow(combined_tidied), "rows\n")
  
  # Apply multiple comparisons correction
  cat("Applying multiple comparisons correction...\n")
  combined_tidied <- apply_multiple_comparisons_correction(combined_tidied)
  
  # Report correction results
  main_effects <- combined_tidied[combined_tidied$term == "indep_var", ]
  if (nrow(main_effects) > 0) {
    n_tests <- nrow(main_effects)
    n_sig_raw <- sum(main_effects$p.value < 0.05, na.rm = TRUE)
    n_sig_fdr <- sum(main_effects$p.value.fdr < 0.05, na.rm = TRUE)
    n_sig_bonf <- sum(main_effects$p.value.bonferroni < 0.05, na.rm = TRUE)
    
    cat("Multiple comparisons correction summary:\n")
    cat("  Total tests:", n_tests, "\n")
    cat("  Significant (raw p < 0.05):", n_sig_raw, sprintf("(%.2f%%)", 100*n_sig_raw/n_tests), "\n")
    cat("  Significant (FDR p < 0.05):", n_sig_fdr, sprintf("(%.2f%%)", 100*n_sig_fdr/n_tests), "\n")
    cat("  Significant (Bonferroni p < 0.05):", n_sig_bonf, sprintf("(%.2f%%)", 100*n_sig_bonf/n_tests), "\n")
  }
} else {
  combined_tidied <- tibble()
  cat("No tidied results to combine\n")
}

# Combine glanced results
if (length(all_glanced) > 0) {
  combined_glanced <- bind_rows(all_glanced)
  cat("Combined glanced results:", nrow(combined_glanced), "rows\n")
} else {
  combined_glanced <- tibble()
  cat("No glanced results to combine\n")
}

# Combine R-squared results
if (length(all_rsq) > 0) {
  combined_rsq <- bind_rows(all_rsq)
  cat("Combined R-squared results:", nrow(combined_rsq), "rows\n")
} else {
  combined_rsq <- tibble()
  cat("No R-squared results to combine\n")
}

# Save aggregated results
cat("\nSaving aggregated results...\n")

# Define output file names
tidied_file <- file.path(output_dir, paste0(analysis_type, "_", normalization, "_tidied_complete.rds"))
glanced_file <- file.path(output_dir, paste0(analysis_type, "_", normalization, "_glanced_complete.rds"))
rsq_file <- file.path(output_dir, paste0(analysis_type, "_", normalization, "_rsq_complete.rds"))

# Save each component
if (nrow(combined_tidied) > 0) {
  saveRDS(combined_tidied, tidied_file)
  cat("Saved tidied results to:", basename(tidied_file), "\n")
  cat("  Dimensions:", nrow(combined_tidied), "x", ncol(combined_tidied), "\n")
}

if (nrow(combined_glanced) > 0) {
  saveRDS(combined_glanced, glanced_file)
  cat("Saved glanced results to:", basename(glanced_file), "\n")
  cat("  Dimensions:", nrow(combined_glanced), "x", ncol(combined_glanced), "\n")
}

if (nrow(combined_rsq) > 0) {
  saveRDS(combined_rsq, rsq_file)
  cat("Saved R-squared results to:", basename(rsq_file), "\n")
  cat("  Dimensions:", nrow(combined_rsq), "x", ncol(combined_rsq), "\n")
}

# Create summary report
summary_file <- file.path(output_dir, paste0(analysis_type, "_", normalization, "_aggregation_summary.txt"))

cat("\nCreating summary report...\n")
sink(summary_file)

cat("=============================================================================\n")
cat("WAS Results Aggregation Summary\n")
cat("=============================================================================\n")
cat("Analysis Type:", analysis_type, "\n")
cat("Normalization:", normalization, "\n")
cat("Aggregation Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input Directory:", input_dir, "\n")
cat("\n")

cat("INPUT FILES:\n")
cat("Total .rds files found:", length(rds_files), "\n")
cat("Successfully processed:", successful_files, "\n")
cat("Failed to process:", length(rds_files) - successful_files, "\n")
cat("\n")

cat("AGGREGATED RESULTS:\n")
cat("Tidied results:", nrow(combined_tidied), "rows x", ncol(combined_tidied), "columns\n")
cat("Glanced results:", nrow(combined_glanced), "rows x", ncol(combined_glanced), "columns\n")
cat("R-squared results:", nrow(combined_rsq), "rows x", ncol(combined_rsq), "columns\n")
cat("\n")

if (nrow(combined_tidied) > 0) {
  main_effects <- combined_tidied[combined_tidied$term == "indep_var", ]
  if (nrow(main_effects) > 0) {
    cat("MULTIPLE COMPARISONS CORRECTION:\n")
    n_tests <- nrow(main_effects)
    n_sig_raw <- sum(main_effects$p.value < 0.05, na.rm = TRUE)
    n_sig_fdr <- sum(main_effects$p.value.fdr < 0.05, na.rm = TRUE)
    n_sig_bonf <- sum(main_effects$p.value.bonferroni < 0.05, na.rm = TRUE)
    
    cat("Total hypothesis tests:", n_tests, "\n")
    cat("Raw p-value < 0.05:", n_sig_raw, sprintf("(%.2f%%)\n", 100*n_sig_raw/n_tests))
    cat("FDR adjusted p-value < 0.05:", n_sig_fdr, sprintf("(%.2f%%)\n", 100*n_sig_fdr/n_tests))
    cat("Bonferroni adjusted p-value < 0.05:", n_sig_bonf, sprintf("(%.2f%%)\n", 100*n_sig_bonf/n_tests))
    cat("\n")
    
    # Top significant results
    if (n_sig_fdr > 0) {
      cat("TOP 10 SIGNIFICANT RESULTS (FDR < 0.05):\n")
      top_results <- main_effects %>%
        filter(p.value.fdr < 0.05) %>%
        arrange(p.value.fdr) %>%
        head(10) %>%
        select(phenotype, exposure, estimate, p.value, p.value.fdr, p.value.bonferroni)
      
      print(top_results)
      cat("\n")
    }
  }
}

cat("OUTPUT FILES:\n")
if (nrow(combined_tidied) > 0) cat("Tidied:", basename(tidied_file), "\n")
if (nrow(combined_glanced) > 0) cat("Glanced:", basename(glanced_file), "\n")
if (nrow(combined_rsq) > 0) cat("R-squared:", basename(rsq_file), "\n")
cat("Summary:", basename(summary_file), "\n")

sink()

cat("Summary report saved to:", basename(summary_file), "\n")
cat("\n=============================================================================\n")
cat("Aggregation completed successfully!\n")
cat("=============================================================================\n") 