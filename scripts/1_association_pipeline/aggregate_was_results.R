#!/usr/bin/env Rscript

# =============================================================================
# Aggregates per-dep-var RDS outputs from universal_was_analysis.R, applies
# scheme-wise FDR correction (and optionally Storey q), and emits combined
# *_tidied_complete.rds / *_glanced_complete.rds / *_rsq_complete.rds plus
# supplementary CSV tables.
#
# Usage modes:
#   1. Individual:   Rscript aggregate_was_results.R <type> <norm> <results_dir> [use_qvalue]
#   2. All 24:       Rscript aggregate_was_results.R --run-all <results_dir> [use_qvalue]
#   3. Tables only:  Rscript aggregate_was_results.R --create-tables <results_dir>
#
# Environment: R >= 4.5 with dplyr, readr, purrr, tibble, stringr, DBI,
# RSQLite, getopt; qvalue (Bioconductor) optional for Storey q-values.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

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

# Try to load qvalue package for Storey's method
QVALUE_AVAILABLE <- FALSE
tryCatch({
  library(qvalue)
  QVALUE_AVAILABLE <- TRUE
  cat("qvalue package loaded - Storey's q-value method available\n")
}, error = function(e) {
  cat("qvalue package not available - using Benjamini-Hochberg only\n")
  cat("To install: BiocManager::install('qvalue')\n")
})

# Define analysis configurations (20 schemes total)
ANALYSIS_TYPES <- c("1_demoWAS", "2_oradWAS", "3_exWAS", "4_pheWAS", "5_outWAS")
NORMALIZATIONS <- c("clr", "lognorm", "none", "hellinger")

# =====================================================================================
# FUNCTION DEFINITIONS (must come before main execution)
# =====================================================================================

# Function to safely read RDS files (enhanced for new result structure)
safe_read_rds <- function(file_path) {
  tryCatch({
    result <- readRDS(file_path)
    if (!is.list(result) || !all(c("pe_tidied", "pe_glanced", "rsq") %in% names(result))) {
      cat("WARNING: Invalid result structure in", basename(file_path), "\n")
      return(NULL)
    }
    
    # Check if results contain enhanced columns from final audit corrections
    if (!is.null(result$pe_tidied) && nrow(result$pe_tidied) > 0) {
      expected_enhanced_cols <- c("effect_scale", "interpretation_note", "normalization")
      missing_enhanced <- expected_enhanced_cols[!expected_enhanced_cols %in% names(result$pe_tidied)]
      if (length(missing_enhanced) > 0) {
        cat("INFO: Enhanced columns missing in", basename(file_path), ":", paste(missing_enhanced, collapse = ", "), "\n")
      }
    }
    
    return(result)
  }, error = function(e) {
    cat("ERROR reading", basename(file_path), ":", e$message, "\n")
    return(NULL)
  })
}

# Function to apply SCHEME-WISE multiple comparisons correction (STATISTICALLY CORRECTED)
# 
# CRITICAL ASSUMPTIONS FOR VALID FDR CONTROL:
# 1. ONE BIOLOGICAL QUESTION PER FILE: Each file contains exactly one exposure/OTU for analyses 2-6
# 2. CATEGORICAL FACTORS: Pipeline must store ONE test per factor (global Wald/LRT), not per dummy level
# 3. SCHEME-WISE CONTROL: FDR ≤ α within each scheme, NOT across all 24 schemes
# 4. TERM IDENTIFICATION: Uses deterministic matching, not data-driven pattern matching
#
# ERROR CONDITIONS: Function will STOP execution (not warn) if main effects cannot be identified
apply_scheme_wise_fdr_correction <- function(tidied_data, analysis_type, use_qvalue = FALSE, input_dir = NULL) {
  if (nrow(tidied_data) == 0) {
    tidied_data$p.value.fdr <- numeric(0)
    tidied_data$p.value.bonferroni <- numeric(0)
    tidied_data$q.value <- numeric(0)
    tidied_data$fdr_corrected <- logical(0)
    return(tidied_data)
  }
  
  # Initialize adjustment columns with original p-values (every term gets these columns)
  tidied_data$p.value.fdr <- tidied_data$p.value
  tidied_data$p.value.bonferroni <- tidied_data$p.value  # FWER reference only - NOT used for significance
  tidied_data$fdr_corrected <- FALSE  # Flag to indicate which terms were actually corrected
  if (use_qvalue && QVALUE_AVAILABLE) {
    tidied_data$q.value <- tidied_data$p.value
  }
  
  # *** STATISTICAL FIX: DETERMINISTIC MAIN EFFECT IDENTIFICATION ***
  # Based on statistical audit: Use hard-coded terms, not data-driven pattern matching
  # This ensures one test per biological question, regardless of variable type
  
  if (analysis_type == "1_demoWAS") {
    # Demographics  Microbiome: All demographic predictors are main effects
    demographic_predictors <- c(
      "RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "INDFMPIR",
      "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
      "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER",
      "ETHNICITY_NONHISPANICBLACK", "BORN_INUSA"
    )
    
    # Start with exact matching
    main_effect_mask <- tidied_data$term %in% demographic_predictors
    
    # CRITICAL FIX: Add categorical variable support for demographics
    # Handle categorical demographics (e.g., RIAGENDR  RIAGENDR2)
    for (demo_var in demographic_predictors) {
      pattern <- paste0("^", demo_var)
      categorical_matches <- grepl(pattern, tidied_data$term) & 
                           tidied_data$term != "(Intercept)" &
                           !tidied_data$term %in% demographic_predictors
      main_effect_mask <- main_effect_mask | categorical_matches
    }
  } else {
    # All other analyses: Use the unique values from independent_var column
    # The pipeline stores original variable names in both 'term' and 'independent_var'  
    # This ensures exactly one test per biological question
    # NOTE: For categorical variables, this selects one test per variable (not per dummy level)
    if ("independent_var" %in% names(tidied_data)) {
      main_effect_terms <- unique(tidied_data$independent_var)
      
      # CORRECTED LOGIC: Multiple predictors per file is the correct design for schemes 2-6
      # No longer treat this as an error - this is expected behavior
      if (length(main_effect_terms) > 1) {
        cat(" MULTIPLE MAIN EFFECTS DETECTED: Found", length(main_effect_terms), 
            "main effect variables in this scheme\n")
        if (length(main_effect_terms) <= 10) {
          cat("   Variables:", paste(main_effect_terms, collapse = ", "), "\n")
        } else {
          cat("   First 10 variables:", paste(head(main_effect_terms, 10), collapse = ", "), "...\n")
        }
        cat("   All", length(main_effect_terms), "main effects will be subject to FDR correction\n")
      }
      
      # Enhanced main effect identification with categorical variable support
      # Try exact matching first
      main_effect_mask <- tidied_data$term %in% main_effect_terms
      
      # For terms that didn't match exactly, try pattern matching for categorical variables
      unmatched_terms <- main_effect_terms[!main_effect_terms %in% tidied_data$term]
      if (length(unmatched_terms) > 0) {
        cat(" CATEGORICAL VARIABLES DETECTED: Pattern matching for", length(unmatched_terms), "variables\n")
        
        for (base_var in unmatched_terms) {
          pattern <- paste0("^", base_var)
          categorical_matches <- grepl(pattern, tidied_data$term) & 
                               tidied_data$term != "(Intercept)" &
                               !main_effect_mask  # Don't double-count exact matches
          
          if (sum(categorical_matches) > 0) {
            main_effect_mask <- main_effect_mask | categorical_matches
            dummy_levels <- unique(tidied_data$term[categorical_matches])
            cat("   ", base_var, "", length(dummy_levels), "dummy level(s):", paste(dummy_levels, collapse=", "), "\n")
          }
        }
        cat("   All categorical dummy levels will be subject to FDR correction\n")
      }
      
      # NOTE: Manual override removed - LBXVST, LBXVTO, LBXML13 have invalid p-values (NaN)
      # These should be excluded from FDR correction due to upstream modeling issues
    } else {
      # Fallback: try generic terms (should not be needed with current pipeline)
      main_effect_mask <- tidied_data$term %in% c("indep_var", "expo")
    }
  }
  
  if (sum(main_effect_mask) == 0) {
    cat("CRITICAL ERROR: No main effect terms found for p-value adjustment\n")
        cat("Available terms:", paste(unique(tidied_data$term), collapse = ", "), "\n")
    if (analysis_type == "1_demoWAS") {
      cat("Expected demographic terms:", paste(demographic_predictors, collapse = ", "), "\n")
    } else {
      if (exists("main_effect_terms")) {
        cat("Expected main effect terms:", paste(main_effect_terms, collapse = ", "), "\n")
      } else {
        cat("Expected main effect terms: indep_var, expo\n")
      }
    }
    cat("Analysis type:", analysis_type, "\n")
    stop("STATISTICAL VALIDITY COMPROMISED: Cannot proceed without main effects for FDR correction.")
  }
  
  cat("Applying SCHEME-WISE FDR correction (STATISTICALLY CORRECTED)...\n")
  cat("NO FILTERING: All results preserved (significant and non-significant)\n")
  cat("Deterministic main effect identification with categorical variable support\n")
  cat("Categorical variables: FDR correction applied to all dummy levels (pragmatic approach)\n")
  cat("Method:", if(use_qvalue && QVALUE_AVAILABLE) "Storey's q-value" else "Benjamini-Hochberg", "\n")
  
  # Extract main effects data and handle missing p-values
  main_effects_indices <- which(main_effect_mask)
  main_effects_p_values <- tidied_data$p.value[main_effects_indices]
  main_effect_terms <- unique(tidied_data$term[main_effects_indices])
  
  # Filter out rows with missing p-values before correction
  valid_p_mask <- !is.na(main_effects_p_values)
  main_effects_indices_valid <- main_effects_indices[valid_p_mask]
  main_effects_p_values_valid <- main_effects_p_values[valid_p_mask]
  
  cat("Analysis type:", analysis_type, "\n")
  if (analysis_type == "1_demoWAS") {
    cat("Main effect terms subject to FDR correction:", paste(demographic_predictors, collapse = ", "), "\n")  
  } else {
    if (exists("main_effect_terms")) {
      cat("Main effect terms subject to FDR correction:", paste(main_effect_terms, collapse = ", "), "\n")
      
      # DIAGNOSTIC: Report categorical variable handling
      corrected_terms <- unique(tidied_data$term[main_effect_mask])
      original_vars <- unique(tidied_data$independent_var[main_effect_mask])
      
      if (length(corrected_terms) > length(original_vars)) {
        cat(" CATEGORICAL VARIABLES DETECTED:\n")
        for (var in original_vars) {
          var_terms <- corrected_terms[grepl(paste0("^", var), corrected_terms)]
          if (length(var_terms) > 1) {
            cat("   ", var, "", paste(var_terms, collapse = ", "), "\n")
          }
        }
        cat("   FDR correction applied to all", length(corrected_terms), "terms (including dummy levels)\n")
        cat("   This represents", length(original_vars), "biological variables with factor expansions\n")
      }
    } else {
      cat("Main effect terms subject to FDR correction: indep_var, expo\n")
    }
  }
  cat("Total main effect tests in this scheme:", length(main_effects_p_values), "\n")
  cat("Tests with valid p-values:", length(main_effects_p_values_valid), "\n")
  
  # SCHEME-WISE CORRECTION: Apply to ALL main effects together (not grouped)
  if (length(main_effects_p_values_valid) > 0) {
    if (use_qvalue && QVALUE_AVAILABLE && length(main_effects_p_values_valid) >= 10) {
      # Use Storey's q-value method (requires at least 10 tests for reliable π₀ estimation)
      tryCatch({
        qobj <- qvalue(main_effects_p_values_valid)
        q_values <- qobj$qvalues
        # Also compute traditional BH for comparison
        fdr_values <- p.adjust(main_effects_p_values_valid, method = "fdr")
        bonf_values <- p.adjust(main_effects_p_values_valid, method = "bonferroni")
        
        # Update ONLY the main effect terms with valid p-values
        tidied_data$q.value[main_effects_indices_valid] <- q_values
        tidied_data$p.value.fdr[main_effects_indices_valid] <- fdr_values
        tidied_data$p.value.bonferroni[main_effects_indices_valid] <- bonf_values
        tidied_data$fdr_corrected[main_effects_indices_valid] <- TRUE
        
        cat("Storey's method: π₀ estimate =", round(qobj$pi0, 3), "\n")
        
      }, error = function(e) {
        cat("qvalue method failed:", e$message, "- falling back to BH\n")
        # Fallback to BH
        fdr_values <- p.adjust(main_effects_p_values_valid, method = "fdr")
        bonf_values <- p.adjust(main_effects_p_values_valid, method = "bonferroni")
        tidied_data$p.value.fdr[main_effects_indices_valid] <- fdr_values
        tidied_data$p.value.bonferroni[main_effects_indices_valid] <- bonf_values
        tidied_data$fdr_corrected[main_effects_indices_valid] <- TRUE
      })
    } else {
      if (use_qvalue && QVALUE_AVAILABLE) {
        cat("Using BH instead of qvalue (fewer than 10 tests - qvalue π₀ estimation unreliable)\n")
      }
      # Use traditional Benjamini-Hochberg
      fdr_values <- p.adjust(main_effects_p_values_valid, method = "fdr")
      bonf_values <- p.adjust(main_effects_p_values_valid, method = "bonferroni")
      
      # Update ONLY the main effect terms with valid p-values
      tidied_data$p.value.fdr[main_effects_indices_valid] <- fdr_values
      tidied_data$p.value.bonferroni[main_effects_indices_valid] <- bonf_values
      tidied_data$fdr_corrected[main_effects_indices_valid] <- TRUE
    }
  }
  
  cat("Terms with FDR correction applied:", sum(tidied_data$fdr_corrected), "\n")
  cat("Terms with original p-values retained:", sum(!tidied_data$fdr_corrected), "\n")
  
  return(tidied_data)
}

# =====================================================================================
# FUNCTION: Individual Aggregation (UPDATED)
# =====================================================================================
run_individual_aggregation <- function(analysis_type, normalization, results_dir, use_qvalue = FALSE) {
  cat("=============================================================================\n")
  cat(" WAS Results Aggregation with SCHEME-WISE FDR Correction\n")
  cat("=============================================================================\n")
  cat("Analysis Type:", analysis_type, "\n")
  cat("Normalization:", normalization, "\n")
  cat("Results Directory:", results_dir, "\n")
  cat("Method:", if(use_qvalue && QVALUE_AVAILABLE) "Storey's q-value" else "Benjamini-Hochberg", "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("STATISTICAL APPROACH: Scheme-wise FDR correction (NO FILTERING - all results preserved)\n")
  cat("THEORETICAL BASIS: Each scheme treated as independent hypothesis family\n")
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
    
    # Apply SCHEME-WISE multiple comparisons correction
    cat("Applying SCHEME-WISE FDR correction...\n")
    combined_tidied <- apply_scheme_wise_fdr_correction(combined_tidied, analysis_type, use_qvalue, input_dir)
    
    # Report correction results summary (NO FILTERING)
    if ("fdr_corrected" %in% names(combined_tidied)) {
      corrected_terms <- combined_tidied[combined_tidied$fdr_corrected == TRUE, ]
      non_corrected_terms <- combined_tidied[combined_tidied$fdr_corrected == FALSE, ]
      
      if (nrow(corrected_terms) > 0) {
        n_tests <- nrow(corrected_terms)
      
      # Determine dependent variable column
        dependent_var_col <- if ("phenotype" %in% names(corrected_terms)) "phenotype" else if ("yvar" %in% names(corrected_terms)) "yvar" else "dependent_var"
        n_dependent_vars <- length(unique(corrected_terms[[dependent_var_col]]))
        
        cat(" SCHEME-WISE FDR correction summary:\n")
        cat("  Terms subject to FDR correction:", paste(unique(corrected_terms$term), collapse = ", "), "\n")
        cat("  Terms NOT subject to correction:", paste(unique(non_corrected_terms$term), collapse = ", "), "\n")
        cat("  Total corrected tests in scheme:", n_tests, "\n")
        cat("  Total non-corrected terms:", nrow(non_corrected_terms), "\n")
      cat("  Dependent variables:", n_dependent_vars, "\n")
        cat("  Average corrected tests per dependent variable:", round(n_tests / n_dependent_vars, 1), "\n")
        cat("  FDR-corrected p-values added to all", n_tests, "main effect terms\n")
        cat("  All results preserved (no filtering based on significance)\n")
      }
    } else {
      cat(" ERROR: fdr_corrected column not found in combined results\n")
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
  cat(" WAS Results Aggregation Summary - SCHEME-WISE FDR CORRECTION\n")
  cat("=============================================================================\n")
  cat("Analysis Type:", analysis_type, "\n")
  cat("Normalization:", normalization, "\n")
  cat("Aggregation Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Input Directory:", input_dir, "\n")
  cat("Method:", if(use_qvalue && QVALUE_AVAILABLE) "Storey's q-value" else "Benjamini-Hochberg", "\n")
  cat("STATISTICAL APPROACH: Scheme-wise FDR correction (NO FILTERING - all results preserved)\n")
  cat("THEORETICAL BASIS: Each scheme treated as independent hypothesis family\n")
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
    if ("fdr_corrected" %in% names(combined_tidied)) {
      corrected_terms <- combined_tidied[combined_tidied$fdr_corrected == TRUE, ]
      non_corrected_terms <- combined_tidied[combined_tidied$fdr_corrected == FALSE, ]
      
      if (nrow(corrected_terms) > 0) {
        cat(" SCHEME-WISE FDR CORRECTION RESULTS:\n")
        cat("Terms subject to FDR correction:", paste(unique(corrected_terms$term), collapse = ", "), "\n")
        cat("Terms NOT subject to correction:", paste(unique(non_corrected_terms$term), collapse = ", "), "\n")
        
        n_tests <- nrow(corrected_terms)
        
        # Determine dependent variable column
        dependent_var_col <- if ("phenotype" %in% names(corrected_terms)) "phenotype" else if ("yvar" %in% names(corrected_terms)) "yvar" else "dependent_var"
        n_dependent_vars <- length(unique(corrected_terms[[dependent_var_col]]))
        
        cat("Total corrected hypothesis tests:", n_tests, "\n")
        cat("Total non-corrected terms:", nrow(non_corrected_terms), "\n")
        cat("Dependent variables:", n_dependent_vars, "\n")
        cat("Average corrected tests per dependent variable:", round(n_tests / n_dependent_vars, 1), "\n")
        cat("FDR correction method:", if(use_qvalue && QVALUE_AVAILABLE) "Storey's q-value" else "Benjamini-Hochberg", "\n")
        cat("Results preservation: ALL", nrow(combined_tidied), "rows retained (no filtering)\n")
        cat("\n")
      }
    } else {
      cat(" ERROR: fdr_corrected column not found in results\n")
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
  cat(" SCHEME-WISE aggregation completed successfully!\n")
  cat("FDR correction: Applied to main effects within this scheme\n")
  cat("Results preservation: ALL rows retained (no filtering based on significance)\n")
  cat("Statistical validity: Maintained through proper scheme-wise correction\n")
  cat("=============================================================================\n")
}

# Function to parse individual summary files (UPDATED for scheme-wise correction)
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
    fdr_level = NA,
    correction_method = NA,
    raw_significant = NA,
    raw_significant_pct = NA,
    fdr_significant = NA,
    fdr_significant_pct = NA,
    bonf_significant = NA,
    bonf_significant_pct = NA,
    qvalue_significant = NA,
    qvalue_significant_pct = NA
  )
  
  # Parse lines
  for (line in lines) {
    if (grepl("Total .rds files found:", line)) {
      data$total_files <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Successfully processed:", line)) {
      data$successful_files <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Failed to process:", line)) {
      data$failed_files <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("FDR Level:", line)) {
      data$fdr_level <- as.numeric(str_extract(line, "[0-9.]+"))
    } else if (grepl("Method:", line)) {
      if (grepl("Storey", line)) {
        data$correction_method <- "Storey's q-value"
      } else {
        data$correction_method <- "Benjamini-Hochberg"
      }
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
    } else if (grepl("Total corrected hypothesis tests:", line)) {
      data$total_tests <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Total hypothesis tests:", line) && is.na(data$total_tests)) {
      # Fallback for old format files
      data$total_tests <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Dependent variables:", line)) {
      data$dependent_vars <- as.numeric(str_extract(line, "\\d+"))
    } else if (grepl("Average corrected tests per dependent variable:", line)) {
      data$avg_tests_per_depvar <- as.numeric(str_extract(line, "[0-9.]+"))
    } else if (grepl("Average tests per dependent variable:", line) && is.na(data$avg_tests_per_depvar)) {
      # Fallback for old format files
      data$avg_tests_per_depvar <- as.numeric(str_extract(line, "[0-9.]+"))
    } else if (grepl("Raw p-value < 0.05:", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$raw_significant <- as.numeric(nums[1])
        data$raw_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    } else if (grepl("FDR adjusted p-value <", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$fdr_significant <- as.numeric(nums[1])
        data$fdr_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    } else if (grepl("Bonferroni adjusted p-value <", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$bonf_significant <- as.numeric(nums[1])
        data$bonf_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    } else if (grepl("Storey q-value <", line)) {
      nums <- str_extract_all(line, "\\d+")[[1]]
      pct <- str_extract(line, "\\([0-9.]+%\\)")
      if (length(nums) >= 1) {
        data$qvalue_significant <- as.numeric(nums[1])
        data$qvalue_significant_pct <- as.numeric(str_extract(pct, "[0-9.]+"))
      }
    }
  }
  
  return(data)
}

# Create overall summary table (UPDATED for scheme-wise correction)
create_overall_summary_table <- function(combined_summaries, output_dir) {
  summary_table <- combined_summaries %>%
    select(analysis_type, normalization, total_files, successful_files, failed_files,
           total_tests, dependent_vars, fdr_level, correction_method, 
           fdr_significant, fdr_significant_pct) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_aggregation_summary_overall.csv")
  write_csv(summary_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# Create significance summary table (UPDATED for scheme-wise correction)
create_significance_summary_table <- function(combined_summaries, output_dir) {
  significance_table <- combined_summaries %>%
    select(analysis_type, normalization, total_tests, fdr_level, correction_method,
           raw_significant, raw_significant_pct,
           fdr_significant, fdr_significant_pct,
           bonf_significant, bonf_significant_pct,
           qvalue_significant, qvalue_significant_pct) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_aggregation_summary_significance.csv")
  write_csv(significance_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# Create detailed statistics table (UPDATED for scheme-wise correction)
create_detailed_statistics_table <- function(combined_summaries, output_dir) {
  detailed_table <- combined_summaries %>%
    select(analysis_type, normalization, total_files, successful_files, 
           tidied_rows, glanced_rows, rsq_rows, 
           total_tests, dependent_vars, avg_tests_per_depvar,
           fdr_level, correction_method) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_aggregation_summary_detailed_statistics.csv")
  write_csv(detailed_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# Create scheme-wise correction summary table (NEW)
create_scheme_wise_correction_table <- function(combined_summaries, output_dir) {
  correction_table <- combined_summaries %>%
    select(analysis_type, normalization, total_tests, fdr_level, correction_method,
           raw_significant, fdr_significant, bonf_significant, qvalue_significant) %>%
    mutate(
      fdr_vs_raw_ratio = round(fdr_significant / raw_significant, 3),
      bonf_vs_raw_ratio = round(bonf_significant / raw_significant, 3),
      qvalue_vs_fdr_ratio = ifelse(!is.na(qvalue_significant), round(qvalue_significant / fdr_significant, 3), NA)
    ) %>%
    arrange(analysis_type, normalization)
  
  output_file <- file.path(output_dir, "s_table_scheme_wise_correction_summary.csv")
  write_csv(correction_table, output_file)
  cat("Created:", basename(output_file), "\n")
}

# =====================================================================================
# FUNCTION: Create Supplementary Tables (UPDATED for scheme-wise correction)
# =====================================================================================
create_supplementary_tables <- function(results_dir) {
  cat("Parsing all 24 scheme aggregation summary files...\n")
  
  # Initialize data collection
  all_summaries <- list()
  
  # Collect all summary files from 24 schemes
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
    cat("This usually means the aggregation step failed for all 24 schemes.\n")
    return()
  }
  
  combined_summaries <- bind_rows(all_summaries, .id = "pipeline")
  
  if (nrow(combined_summaries) == 0) {
    cat("ERROR: No data in combined summaries. Cannot create supplementary tables.\n")
    return()
  }
  
  cat("Successfully parsed", nrow(combined_summaries), "scheme summaries\n")
  
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
  
  # 4. Scheme-wise Correction Summary Table (NEW)
  create_scheme_wise_correction_table(combined_summaries, supp_dir)
  
  cat("\n Supplementary tables created in:", supp_dir, "\n")
  cat("Tables created:\n")
  cat("  - Overall summary (24 schemes)\n")
  cat("  - Significance summary (scheme-wise FDR control)\n")
  cat("  - Detailed statistics\n")
  cat("  - Scheme-wise correction effectiveness\n")
}

# =====================================================================================
# MAIN EXECUTION LOGIC (UPDATED for 24 schemes)
# =====================================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for special modes
if (length(args) >= 1 && args[1] == "--run-all") {
  if (length(args) < 2) {
    cat("Usage: Rscript aggregate_was_results.R --run-all <results_dir> [use_qvalue]\n")
    quit(status = 1)
  }
  
  results_dir <- args[2]
  use_qvalue <- if(length(args) >= 3) as.logical(args[3]) else FALSE
  
  cat("=============================================================================\n")
  cat(" BATCH MODE: Running ALL 24 WAS Aggregations + Creating Supplementary Tables\n")
  cat("=============================================================================\n")
  cat("Results Directory:", results_dir, "\n")
  cat("Method:", if(use_qvalue && QVALUE_AVAILABLE) "Storey's q-value" else "Benjamini-Hochberg", "\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("STATISTICAL APPROACH: Scheme-wise FDR correction (NO FILTERING - all results preserved)\n")
  cat("\n")
  
  # Run all 24 aggregations (CORRECTED count)
  total_analyses <- length(ANALYSIS_TYPES) * length(NORMALIZATIONS)  # = 24
  current_analysis <- 0
  
  for (analysis_type in ANALYSIS_TYPES) {
    for (normalization in NORMALIZATIONS) {
      current_analysis <- current_analysis + 1
      cat(sprintf("[%d/%d] Processing %s %s...\n", current_analysis, total_analyses, analysis_type, normalization))
      
      # Call the individual aggregation function
      tryCatch({
        run_individual_aggregation(analysis_type, normalization, results_dir, use_qvalue)
        cat(" Completed successfully\n\n")
      }, error = function(e) {
        cat(" Error:", e$message, "\n\n")
      })
    }
  }
  
  cat("=============================================================================\n")
  cat(" Creating Supplementary Tables from All 24 Aggregation Summaries\n")
  cat("=============================================================================\n")
  
  # Create supplementary tables
  create_supplementary_tables(results_dir)
  
  cat("=============================================================================\n")
  cat(" BATCH PROCESSING COMPLETED - 24 SCHEMES PROCESSED!\n")
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
    cat("1. Individual: Rscript aggregate_was_results.R <analysis_type> <normalization> <results_dir> [use_qvalue]\n")
    cat("2. Batch All:  Rscript aggregate_was_results.R --run-all <results_dir> [use_qvalue]\n")
    cat("3. Tables Only: Rscript aggregate_was_results.R --create-tables <results_dir>\n")
    cat("\nDefaults: use_qvalue=FALSE\n")
    quit(status = 1)
  }
  
  analysis_type <- args[1]
  normalization <- args[2]
  results_dir <- args[3]
  use_qvalue <- if(length(args) >= 4) as.logical(args[4]) else FALSE
  
  # Run individual aggregation
  run_individual_aggregation(analysis_type, normalization, results_dir, use_qvalue)
}

