#!/usr/bin/env Rscript

# =====================================================================================
#  WAS Results Aggregation Script with Multiple Comparisons Correction
# aggregate_was_results.R
# =====================================================================================
# This script aggregates pipeline results and applies multiple comparisons correction
#
# USAGE MODES:
# 1. Individual aggregation:
#    Rscript aggregate_was_results.R <analysis_type> <normalization> <results_dir>
#    Example: Rscript aggregate_was_results.R 1_demoWAS clr /path/to/results
#
# 2. Run ALL 18 aggregations + create supplementary tables:
#    Rscript aggregate_was_results.R --run-all <results_dir>
#    Example: Rscript aggregate_was_results.R --run-all $(pwd)/results
#
# 3. Create supplementary tables only (after aggregations are done):
#    Rscript aggregate_was_results.R --create-tables <results_dir>
#    Example: Rscript aggregate_was_results.R --create-tables $(pwd)/results
#
# KEY FEATURES: 
# - FDR correction applied WITHIN each dependent variable (microbiome OTU)
# - Maintains FDR < 0.05 as standard threshold
# - Automatic batch processing of all 18 pipelines
# - Creates publication-ready supplementary tables
# =====================================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(DBI)
  library(RSQLite)
  library(getopt)
})

# Define analysis configurations
ANALYSIS_TYPES <- c("1_demoWAS", "2_oradWAS", "3_exWAS", "4_pheWAS", "5_outWAS", "6_zimWAS")
NORMALIZATIONS <- c("clr", "lognorm", "none")

# =====================================================================================
# FUNCTION DEFINITIONS (must come before main execution)
# =====================================================================================

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
  
  # Initialize adjustment columns with original p-values
  tidied_data$p.value.fdr <- tidied_data$p.value
  tidied_data$p.value.bonferroni <- tidied_data$p.value
  
  # Only adjust p-values for the main effect terms (indep_var), not intercepts or covariates
  main_effect_mask <- tidied_data$term == "indep_var"
  
  if (sum(main_effect_mask) == 0) {
    # Try alternative term names from nhanespewas results
    main_effect_mask <- tidied_data$term == "expo"
    if (sum(main_effect_mask) == 0) {
      # Try other possible main effect term names
      main_effect_mask <- tidied_data$term %in% c("exposure", "independent_var", "xvar")
      if (sum(main_effect_mask) == 0) {
        cat("WARNING: No main effect terms found for p-value adjustment\n")
        cat("Available terms:", paste(unique(tidied_data$term), collapse = ", "), "\n")
        return(tidied_data)
      }
    }
  }
  
  cat("Applying FDR correction by dependent variable...\n")
  
  # Determine which column contains the dependent variable
  dependent_var_col <- NULL
  if ("phenotype" %in% names(tidied_data)) {
    dependent_var_col <- "phenotype"
  } else if ("yvar" %in% names(tidied_data)) {
    dependent_var_col <- "yvar"
  } else if ("dependent_var" %in% names(tidied_data)) {
    dependent_var_col <- "dependent_var"
  } else {
    cat("ERROR: Cannot identify dependent variable column\n")
    return(tidied_data)
  }
  
  # Extract main effects data
  main_effects_indices <- which(main_effect_mask)
  main_effects_data <- tidied_data[main_effects_indices, ]
  
  cat("Found", length(unique(main_effects_data[[dependent_var_col]])), "unique dependent variables\n")
  
  # EFFICIENT VECTORIZED APPROACH: Use dplyr group operations
  cat("Performing correction (vectorized approach)...\n")
  
  # Apply FDR correction within each dependent variable using group_by
  corrected_data <- main_effects_data %>%
    group_by(!!sym(dependent_var_col)) %>%
    mutate(
      p.value.fdr = p.adjust(p.value, method = "fdr"),
      p.value.bonferroni = p.adjust(p.value, method = "bonferroni")
    ) %>%
    ungroup()
  
  # Update the original data using vectorized indexing
  tidied_data$p.value.fdr[main_effects_indices] <- corrected_data$p.value.fdr
  tidied_data$p.value.bonferroni[main_effects_indices] <- corrected_data$p.value.bonferroni
  
  return(tidied_data)
}

# =====================================================================================
# FUNCTION: Individual Aggregation
# =====================================================================================
run_individual_aggregation <- function(analysis_type, normalization, results_dir) {
  cat("=============================================================================\n")
  cat(" WAS Results Aggregation with Multiple Comparisons Correction\n")
  cat("=============================================================================\n")
  cat("Analysis Type:", analysis_type, "\n")
  cat("Normalization:", normalization, "\n")
  cat("Results Directory:", results_dir, "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("METHOD: FDR correction by dependent variable (FDR < 0.05)\n")
  cat("PIPELINE: Restored nhanespewas logic, no pre-filtering\n")
  cat("\n")

  # Define paths
  input_dir <- file.path(results_dir, paste0(analysis_type, "_out"), paste0("result_", normalization))
  output_dir <- input_dir  # Save aggregated results in the same directory

  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }

  # Find all .rds result files
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  # Exclude model files and aggregated files
  rds_files <- rds_files[!grepl("_models\\.rds$", rds_files)]
  rds_files <- rds_files[!grepl("_complete\\.rds$", rds_files)]
  rds_files <- rds_files[!grepl("_stratified\\.rds$", rds_files)]

  if (length(rds_files) == 0) {
    stop("No individual result files found in: ", input_dir)
  }

  cat("Found", length(rds_files), "individual result files\n")

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
    cat("Applying multiple comparisons correction (FDR < 0.05 standard)...\n")
    combined_tidied <- apply_multiple_comparisons_correction(combined_tidied)
    
    # Report correction results
    main_effects <- combined_tidied[combined_tidied$term %in% c("indep_var", "expo"), ]
    if (nrow(main_effects) > 0) {
      n_tests <- nrow(main_effects)
      n_sig_raw <- sum(main_effects$p.value < 0.05, na.rm = TRUE)
      n_sig_fdr <- sum(main_effects$p.value.fdr < 0.05, na.rm = TRUE)
      n_sig_bonf <- sum(main_effects$p.value.bonferroni < 0.05, na.rm = TRUE)
      
      # Determine dependent variable column
      dependent_var_col <- if ("phenotype" %in% names(main_effects)) "phenotype" else if ("yvar" %in% names(main_effects)) "yvar" else "dependent_var"
      n_dependent_vars <- length(unique(main_effects[[dependent_var_col]]))
      
      cat(" multiple comparisons correction summary:\n")
      cat("  Total tests:", n_tests, "\n")
      cat("  Dependent variables:", n_dependent_vars, "\n")
      cat("  Average tests per dependent variable:", round(n_tests / n_dependent_vars, 1), "\n")
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

  # Define output file names with prefix
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
  cat(" WAS Results Aggregation Summary - CORRECTION\n")
  cat("=============================================================================\n")
  cat("Analysis Type:", analysis_type, "\n")
  cat("Normalization:", normalization, "\n")
  cat("Aggregation Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Input Directory:", input_dir, "\n")
  cat("METHOD: FDR correction by dependent variable (FDR < 0.05)\n")
  cat("PIPELINE: Restored nhanespewas logic, no pre-filtering\n")
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
    main_effects <- combined_tidied[combined_tidied$term %in% c("indep_var", "expo"), ]
    if (nrow(main_effects) > 0) {
      cat(" MULTIPLE COMPARISONS CORRECTION:\n")
      n_tests <- nrow(main_effects)
      n_sig_raw <- sum(main_effects$p.value < 0.05, na.rm = TRUE)
      n_sig_fdr <- sum(main_effects$p.value.fdr < 0.05, na.rm = TRUE)
      n_sig_bonf <- sum(main_effects$p.value.bonferroni < 0.05, na.rm = TRUE)
      
      # Determine dependent variable column
      dependent_var_col <- if ("phenotype" %in% names(main_effects)) "phenotype" else if ("yvar" %in% names(main_effects)) "yvar" else "dependent_var"
      n_dependent_vars <- length(unique(main_effects[[dependent_var_col]]))
      
      cat("Total hypothesis tests:", n_tests, "\n")
      cat("Dependent variables:", n_dependent_vars, "\n")
      cat("Average tests per dependent variable:", round(n_tests / n_dependent_vars, 1), "\n")
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
          select(any_of(c("phenotype", "yvar", "dependent_var", "exposure", "xvar", "estimate", "p.value", "p.value.fdr", "p.value.bonferroni")))
        
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
  cat(" aggregation completed successfully!\n")
  cat("Expected: Hundreds of significant associations with restored nhanespewas logic\n")
  cat("Standard: FDR < 0.05 maintained as proper threshold\n")
  cat("=============================================================================\n")
}

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

# Create overall summary table
create_overall_summary_table <- function(combined_summaries, output_dir) {
  summary_table <- combined_summaries %>%
    select(analysis_type, normalization, total_files, successful_files, failed_files,
           total_tests, dependent_vars, fdr_significant, fdr_significant_pct) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_aggregation_summary_overall.csv")
  write_csv(summary_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# Create significance summary table
create_significance_summary_table <- function(combined_summaries, output_dir) {
  significance_table <- combined_summaries %>%
    select(analysis_type, normalization, total_tests, 
           raw_significant, raw_significant_pct,
           fdr_significant, fdr_significant_pct,
           bonf_significant, bonf_significant_pct) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_aggregation_summary_significance.csv")
  write_csv(significance_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# Create detailed statistics table
create_detailed_statistics_table <- function(combined_summaries, output_dir) {
  detailed_table <- combined_summaries %>%
    select(analysis_type, normalization, total_files, successful_files, 
           tidied_rows, glanced_rows, rsq_rows, 
           total_tests, dependent_vars, avg_tests_per_depvar) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_aggregation_summary_detailed_statistics.csv")
  write_csv(detailed_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# =====================================================================================
# FUNCTION: Create Supplementary Tables
# =====================================================================================
create_supplementary_tables <- function(results_dir) {
  cat("Parsing all aggregation summary files...\n")
  
  # Initialize data collection
  all_summaries <- list()
  
  # Collect all summary files
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
  if (length(all_summaries) == 0) {
    cat("ERROR: No summary files found. Cannot create supplementary tables.\n")
    cat("This usually means the aggregation step failed for all analyses.\n")
    return()
  }
  
  combined_summaries <- bind_rows(all_summaries, .id = "pipeline")
  
  if (nrow(combined_summaries) == 0) {
    cat("ERROR: No data in combined summaries. Cannot create supplementary tables.\n")
    return()
  }
  
  # Create output directory for supplementary tables
  supp_dir <- file.path(results_dir, "supplementary_tables")
  if (!dir.exists(supp_dir)) {
    dir.create(supp_dir, recursive = TRUE)
  }
  
  # Create different supplementary tables
  
  # 1. Overall Summary Table
  create_overall_summary_table(combined_summaries, supp_dir)
  
  # 2. Significance Summary Table
  create_significance_summary_table(combined_summaries, supp_dir)
  
  # 3. Detailed Statistics Table
  create_detailed_statistics_table(combined_summaries, supp_dir)
  
  cat("\n✅ Supplementary tables created in:", supp_dir, "\n")
}

# =====================================================================================
# MAIN EXECUTION LOGIC
# =====================================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for special modes
if (length(args) >= 1 && args[1] == "--run-all") {
  if (length(args) < 2) {
    cat("Usage: Rscript aggregate_was_results.R --run-all <results_dir>\n")
    quit(status = 1)
  }
  
  results_dir <- args[2]
  cat("=============================================================================\n")
  cat(" BATCH MODE: Running ALL 18 WAS Aggregations + Creating Supplementary Tables\n")
  cat("=============================================================================\n")
  cat("Results Directory:", results_dir, "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("\n")
  
  # Run all 18 aggregations
  total_analyses <- length(ANALYSIS_TYPES) * length(NORMALIZATIONS)
  current_analysis <- 0
  
  for (analysis_type in ANALYSIS_TYPES) {
    for (normalization in NORMALIZATIONS) {
      current_analysis <- current_analysis + 1
      cat(sprintf("[%d/%d] Processing %s %s...\n", current_analysis, total_analyses, analysis_type, normalization))
      
      # Call the individual aggregation function
      tryCatch({
        run_individual_aggregation(analysis_type, normalization, results_dir)
        cat("✅ Completed successfully\n\n")
      }, error = function(e) {
        cat("❌ Error:", e$message, "\n\n")
      })
    }
  }
  
  cat("=============================================================================\n")
  cat(" Creating Supplementary Tables from All Aggregation Summaries\n")
  cat("=============================================================================\n")
  
  # Create supplementary tables
  create_supplementary_tables(results_dir)
  
  cat("=============================================================================\n")
  cat(" BATCH PROCESSING COMPLETED!\n")
  cat("=============================================================================\n")
  quit(status = 0)
  
} else if (length(args) >= 1 && args[1] == "--create-tables") {
  if (length(args) < 2) {
    cat("Usage: Rscript aggregate_was_results.R --create-tables <results_dir>\n")
    quit(status = 1)
  }
  
  results_dir <- args[2]
  cat("=============================================================================\n")
  cat(" Creating Supplementary Tables from Existing Aggregation Summaries\n")
  cat("=============================================================================\n")
  cat("Results Directory:", results_dir, "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("\n")
  
  create_supplementary_tables(results_dir)
  quit(status = 0)
  
} else {
  # Individual aggregation mode
  if (length(args) < 3) {
    cat("Usage Options:\n")
    cat("1. Individual: Rscript aggregate_was_results.R <analysis_type> <normalization> <results_dir>\n")
    cat("2. Batch All:  Rscript aggregate_was_results.R --run-all <results_dir>\n")
    cat("3. Tables Only: Rscript aggregate_was_results.R --create-tables <results_dir>\n")
    quit(status = 1)
  }
  
  analysis_type <- args[1]
  normalization <- args[2]
  results_dir <- args[3]
  
  # Run individual aggregation
  run_individual_aggregation(analysis_type, normalization, results_dir)
}

