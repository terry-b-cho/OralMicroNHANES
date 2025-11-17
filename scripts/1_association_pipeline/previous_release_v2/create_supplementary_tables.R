#!/usr/bin/env Rscript

# =====================================================================================
#  Create Supplementary Tables from WAS Aggregation Summaries
# create_supplementary_tables.R
# =====================================================================================
# This script creates publication-ready supplementary tables from existing 
# aggregation summary files
#
# Usage: Rscript create_supplementary_tables.R <results_dir>
# Example: Rscript create_supplementary_tables.R $(pwd)/results
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript create_supplementary_tables.R <results_dir>\n")
  cat("Example: Rscript create_supplementary_tables.R $(pwd)/results\n")
  quit(status = 1)
}

results_dir <- args[1]

# Define analysis configurations
ANALYSIS_TYPES <- c("1_demoWAS", "2_oradWAS", "3_exWAS", "4_pheWAS", "5_outWAS", "6_zimWAS")
NORMALIZATIONS <- c("clr", "lognorm", "none")

cat("=============================================================================\n")
cat(" Creating Supplementary Tables from WAS Aggregation Summaries\n")
cat("=============================================================================\n")
cat("Results Directory:", results_dir, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n")

# Function to parse individual summary files
parse_summary_file <- function(file_path, analysis_type, normalization) {
  lines <- readLines(file_path)
  
  # Initialize data
  data <- list(
    analysis_type = analysis_type,
    normalization = normalization,
    total_files = NA,
    successful_files = NA,
    failed_files = NA,
    tidied_rows = NA,
    tidied_cols = NA,
    glanced_rows = NA,
    glanced_cols = NA,
    rsq_rows = NA,
    rsq_cols = NA,
    total_tests = NA,
    dependent_vars = NA,
    avg_tests_per_depvar = NA,
    raw_significant = NA,
    raw_significant_pct = NA,
    fdr_significant = NA,
    fdr_significant_pct = NA,
    bonf_significant = NA,
    bonf_significant_pct = NA
  )
  
  # Parse lines
  for (line in lines) {
    if (grepl("Total .rds files found:", line)) {
      data$total_files <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Successfully processed:", line)) {
      data$successful_files <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Failed to process:", line)) {
      data$failed_files <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Tidied results:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      if (length(nums) >= 2) {
        data$tidied_rows <- as.numeric(nums[1])
        data$tidied_cols <- as.numeric(nums[2])
      }
    } else if (grepl("Glanced results:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      if (length(nums) >= 2) {
        data$glanced_rows <- as.numeric(nums[1])
        data$glanced_cols <- as.numeric(nums[2])
      }
    } else if (grepl("R-squared results:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      if (length(nums) >= 2) {
        data$rsq_rows <- as.numeric(nums[1])
        data$rsq_cols <- as.numeric(nums[2])
      }
    } else if (grepl("Total hypothesis tests:", line)) {
      data$total_tests <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Dependent variables:", line)) {
      data$dependent_vars <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Average tests per dependent variable:", line)) {
      data$avg_tests_per_depvar <- as.numeric(str_extract(line, "[0-9.]+"))
    } else if (grepl("Raw p-value < 0.05:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$raw_significant <- as.numeric(nums[1])
        data$raw_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    } else if (grepl("FDR adjusted p-value < 0.05:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$fdr_significant <- as.numeric(nums[1])
        data$fdr_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    } else if (grepl("Bonferroni adjusted p-value < 0.05:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$bonf_significant <- as.numeric(nums[1])
        data$bonf_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    }
  }
  
  return(data)
}

# Initialize data collection
all_summaries <- list()

# Collect all summary files
cat("Parsing aggregation summary files...\n")
for (analysis_type in ANALYSIS_TYPES) {
  for (normalization in NORMALIZATIONS) {
    summary_file <- file.path(results_dir, paste0(analysis_type, "_out"), 
                             paste0("result_", normalization),
                             paste0(analysis_type, "_", normalization, "_aggregation_summary.txt"))
    
    if (file.exists(summary_file)) {
      cat("Reading:", basename(summary_file), "\n")
      
      # Parse summary file
      summary_data <- parse_summary_file(summary_file, analysis_type, normalization)
      all_summaries[[paste(analysis_type, normalization, sep = "_")]] <- summary_data
    } else {
      cat("WARNING: Missing summary file:", summary_file, "\n")
    }
  }
}

# Combine all summaries
combined_summaries <- bind_rows(all_summaries, .id = "pipeline")

# Create output directory for supplementary tables
supp_dir <- file.path(results_dir, "supplementary_tables")
if (!dir.exists(supp_dir)) {
  dir.create(supp_dir, recursive = TRUE)
}

cat("\nCreating supplementary tables...\n")

# 1. Overall Summary Table
cat("Creating overall summary table...\n")
summary_table <- combined_summaries %>%
  select(analysis_type, normalization, total_files, successful_files, failed_files,
         total_tests, dependent_vars, fdr_significant, fdr_significant_pct) %>%
  arrange(analysis_type, normalization)

output_file <- file.path(supp_dir, "s_table_aggregation_summary_overall.csv")
write_csv(summary_table, output_file)
cat("Created:", basename(output_file), "\n")

# 2. Significance Summary Table
cat("Creating significance summary table...\n")
significance_table <- combined_summaries %>%
  select(analysis_type, normalization, total_tests, 
         raw_significant, raw_significant_pct,
         fdr_significant, fdr_significant_pct,
         bonf_significant, bonf_significant_pct) %>%
  arrange(analysis_type, normalization)

output_file <- file.path(supp_dir, "s_table_aggregation_summary_significance.csv")
write_csv(significance_table, output_file)
cat("Created:", basename(output_file), "\n")

# 3. Detailed Statistics Table
cat("Creating detailed statistics table...\n")
detailed_table <- combined_summaries %>%
  select(analysis_type, normalization, total_files, successful_files, 
         tidied_rows, glanced_rows, rsq_rows, 
         total_tests, dependent_vars, avg_tests_per_depvar) %>%
  arrange(analysis_type, normalization)

output_file <- file.path(supp_dir, "s_table_aggregation_summary_detailed_statistics.csv")
write_csv(detailed_table, output_file)
cat("Created:", basename(output_file), "\n")

cat("\n=============================================================================\n")
cat("✅ Supplementary tables created successfully!\n")
cat("Output directory:", supp_dir, "\n")
cat("=============================================================================\n")

# Print summary of what was created
cat("\nSUMMARY:\n")
cat("Processed", nrow(combined_summaries), "pipeline summaries\n")
cat("Created 3 supplementary tables:\n")
cat("  1. s_table_aggregation_summary_overall.csv\n")
cat("  2. s_table_aggregation_summary_significance.csv\n")
cat("  3. s_table_aggregation_summary_detailed_statistics.csv\n")
cat("\nThese tables are ready for manuscript submission.\n")
