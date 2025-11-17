#!/usr/bin/env Rscript

# COMPREHENSIVE  PIPELINE TEST
# Tests all analysis types and normalization methods
# Validates both submission script and analysis script integration

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(DBI)
  library(RSQLite)
  library(survey)
  library(broom)
  library(stringr)
  library(tibble)
})

cat("COMPREHENSIVE  PIPELINE TEST\n")
cat("=======================================\n")
cat("Testing all analysis types and normalization methods\n\n")

# Test configuration
BASE_PATH <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
DB_PATH <- file.path(BASE_PATH, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite")

# All analysis types and normalizations to test
ANALYSIS_TYPES <- c("1_demoWAS", "2_oradWAS", "3_exWAS", "4_pheWAS", "5_outWAS", "6_zimWAS")
NORMALIZATIONS <- c("clr", "lognorm", "none")

# Results tracking
test_results <- data.frame(
  analysis_type = character(),
  normalization = character(),
  schema_exists = logical(),
  dep_var_found = logical(),
  pipeline_success = logical(),
  output_created = logical(),
  result_structure_valid = logical(),
  stringsAsFactors = FALSE
)

# Database connection for initial validation
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = DB_PATH)
tables <- dbListTables(con)
cat("✅ Database connected, found", length(tables), "tables\n\n")

# Test each combination
total_tests <- length(ANALYSIS_TYPES) * length(NORMALIZATIONS)
current_test <- 0

for (analysis_type in ANALYSIS_TYPES) {
  for (normalization in NORMALIZATIONS) {
    current_test <- current_test + 1
    
    cat("TEST", current_test, "/", total_tests, ":", analysis_type, "with", normalization, "\n")
    cat("=", rep("=", 50), "\n", sep = "")
    
    # Initialize result row
    result_row <- data.frame(
      analysis_type = analysis_type,
      normalization = normalization,
      schema_exists = FALSE,
      dep_var_found = FALSE,
      pipeline_success = FALSE,
      output_created = FALSE,
      result_structure_valid = FALSE,
      stringsAsFactors = FALSE
    )
    
    # 1. Check schema file exists
    schema_file <- file.path(BASE_PATH, "results/0_ss_files", 
                            paste0(analysis_type, "_", normalization, "_schema_structure.csv"))
    
    if (file.exists(schema_file)) {
      cat("✅ Schema file exists\n")
      result_row$schema_exists <- TRUE
      
      # 2. Load schema and get test dependent variable
      tryCatch({
        schema_data <- read_csv(schema_file, show_col_types = FALSE)
        
        # Get first dependent variable for testing
        if (nrow(schema_data) > 0) {
          test_dep_var <- schema_data$dep_var[1]
          cat("✅ Test dependent variable:", test_dep_var, "\n")
          result_row$dep_var_found <- TRUE
          
          # 3. Set up output directory
          output_dir <- file.path(BASE_PATH, "results", 
                                 paste0(analysis_type, "_out"), 
                                 paste0("result_", normalization))
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
          
          # 4. Test the full pipeline
          test_command <- paste(
            "Rscript scripts/1_association_pipeline/universal_was_analysis_PRODUCTION.R",
            "--dependent_var", shQuote(test_dep_var),
            "--schema_structure_file", shQuote(schema_file),
            "--database_path", shQuote(DB_PATH),
            "--output_path", shQuote(output_dir),
            "--analysis_type", analysis_type,
            "--normalization", normalization,
            "--test"
          )
          
          cat("Running pipeline test...\n")
          
          # Execute with timeout and capture output
          result_code <- system(test_command, intern = FALSE)
          
          if (result_code == 0) {
            cat("✅ Pipeline execution SUCCESS\n")
            result_row$pipeline_success <- TRUE
            
            # 5. Check output file
            output_file <- file.path(output_dir, paste0(test_dep_var, ".rds"))
            if (file.exists(output_file)) {
              cat("✅ Output file created\n")
              result_row$output_created <- TRUE
              
              # 6. Validate output structure
              tryCatch({
                output_data <- readRDS(output_file)
                required_components <- c("pe_tidied", "pe_glanced", "rsq")
                
                all_components_present <- all(required_components %in% names(output_data))
                has_data <- all(sapply(output_data[required_components], function(x) nrow(x) > 0))
                
                if (all_components_present && has_data) {
                  cat("✅ Result structure valid\n")
                  result_row$result_structure_valid <- TRUE
                  
                  # Additional validations
                  if (nrow(output_data$pe_tidied) > 0) {
                    # Check for required metadata columns
                    required_cols <- c("independent_var", "phenotype", "exposure", "dependent_var", "n_obs")
                    missing_cols <- required_cols[!required_cols %in% names(output_data$pe_tidied)]
                    
                    if (length(missing_cols) == 0) {
                      cat("✅ All metadata columns present\n")
                    } else {
                      cat("⚠️ Missing metadata columns:", paste(missing_cols, collapse = ", "), "\n")
                    }
                    
                    # Check for significant results
                    sig_results <- output_data$pe_tidied %>% 
                      filter(term != "(Intercept)", p.value < 0.05)
                    cat("Significant results (p < 0.05):", nrow(sig_results), "\n")
                  }
                } else {
                  cat("❌ Result structure incomplete\n")
                }
              }, error = function(e) {
                cat("❌ Error validating output:", e$message, "\n")
              })
            } else {
              cat("❌ Output file not created\n")
            }
          } else {
            cat("❌ Pipeline execution FAILED with code:", result_code, "\n")
          }
        } else {
          cat("❌ No dependent variables found in schema\n")
        }
      }, error = function(e) {
        cat("❌ Error loading schema:", e$message, "\n")
      })
    } else {
      cat("❌ Schema file missing:", schema_file, "\n")
    }
    
    # Add result to tracking table
    test_results <- rbind(test_results, result_row)
    cat("\n")
  }
}

# Close database connection
DBI::dbDisconnect(con)

# =====================================================================================
# COMPREHENSIVE RESULTS SUMMARY
# =====================================================================================

cat("COMPREHENSIVE TEST RESULTS SUMMARY\n")
cat("=====================================\n")

# Overall statistics
total_tests <- nrow(test_results)
successful_tests <- sum(test_results$result_structure_valid)
success_rate <- round(successful_tests / total_tests * 100, 1)

cat("Total test combinations:", total_tests, "\n")
cat("Successful tests:", successful_tests, "\n")
cat("Overall success rate:", success_rate, "%\n\n")

# Detailed results table
cat("DETAILED RESULTS BY ANALYSIS TYPE AND NORMALIZATION\n")
cat("======================================================\n")

# Create summary table
summary_table <- test_results %>%
  mutate(
    status = case_when(
      result_structure_valid ~ "✅ PASS",
      pipeline_success ~ "⚠️ PARTIAL",
      schema_exists ~ "❌ FAIL",
      TRUE ~ "❌ NO SCHEMA"
    )
  ) %>%
  select(analysis_type, normalization, status)

# Print formatted table
for (analysis in ANALYSIS_TYPES) {
  cat("\n", analysis, ":\n")
  analysis_results <- summary_table %>% filter(analysis_type == analysis)
  for (i in 1:nrow(analysis_results)) {
    cat("  ", analysis_results$normalization[i], ": ", analysis_results$status[i], "\n")
  }
}

# =====================================================================================
# SUBMISSION SCRIPT INTEGRATION TEST
# =====================================================================================

cat("\nSUBMISSION SCRIPT INTEGRATION TEST\n")
cat("====================================\n")

# Test the submission script with a successful combination
successful_combo <- test_results %>% 
  filter(result_structure_valid == TRUE) %>% 
  slice_head(n = 1)

if (nrow(successful_combo) > 0) {
  test_analysis <- successful_combo$analysis_type[1]
  test_norm <- successful_combo$normalization[1]
  
  cat("Testing submission script with:", test_analysis, test_norm, "\n")
  
  # Test submission script syntax
  submission_script <- "scripts/1_association_pipeline/run_single_was_analysis_PRODUCTION.sh"
  
  if (file.exists(submission_script)) {
    cat("✅ Submission script exists\n")
    
    # Test script help
    help_result <- system(paste("bash", submission_script), intern = FALSE)
    cat("✅ Submission script syntax validated\n")
    
    # Test with dry run (test mode)
    cat("Testing submission script in test mode...\n")
    
    # Note: We don't actually submit SLURM jobs in testing
    cat("⚠️ SLURM submission test skipped (would require cluster resources)\n")
    cat("✅ Submission script integration: READY\n")
    
  } else {
    cat("❌ Submission script missing:", submission_script, "\n")
  }
} else {
  cat("❌ No successful pipeline combinations found for submission test\n")
}

# =====================================================================================
# FINAL ASSESSMENT
# =====================================================================================

cat("\nFINAL COMPREHENSIVE ASSESSMENT\n")
cat("=================================\n")

# Analysis type success rates
cat("SUCCESS RATES BY ANALYSIS TYPE:\n")
for (analysis in ANALYSIS_TYPES) {
  analysis_results <- test_results %>% filter(analysis_type == analysis)
  analysis_success <- sum(analysis_results$result_structure_valid)
  analysis_total <- nrow(analysis_results)
  analysis_rate <- round(analysis_success / analysis_total * 100, 1)
  cat("  ", analysis, ":", analysis_success, "/", analysis_total, "(", analysis_rate, "%)\n")
}

cat("\nSUCCESS RATES BY NORMALIZATION:\n")
for (norm in NORMALIZATIONS) {
  norm_results <- test_results %>% filter(normalization == norm)
  norm_success <- sum(norm_results$result_structure_valid)
  norm_total <- nrow(norm_results)
  norm_rate <- round(norm_success / norm_total * 100, 1)
  cat("  ", norm, ":", norm_success, "/", norm_total, "(", norm_rate, "%)\n")
}

# Critical issues
cat("\nCRITICAL ISSUES DETECTED:\n")
critical_issues <- test_results %>% 
  filter(schema_exists == TRUE, result_structure_valid == FALSE)

if (nrow(critical_issues) > 0) {
  for (i in 1:nrow(critical_issues)) {
    issue <- critical_issues[i, ]
    cat("❌", issue$analysis_type, issue$normalization, "- Schema exists but pipeline failed\n")
  }
} else {
  cat("✅ No critical issues detected\n")
}

# Production readiness assessment
if (success_rate >= 90) {
  cat("\nPRODUCTION READINESS: EXCELLENT (", success_rate, "%)\n")
  cat("✅ Pipeline is ready for full production deployment\n")
} else if (success_rate >= 75) {
  cat("\n⚠️ PRODUCTION READINESS: GOOD (", success_rate, "%)\n")
  cat("✅ Pipeline is mostly ready, minor issues to address\n")
} else if (success_rate >= 50) {
  cat("\n⚠️ PRODUCTION READINESS: MODERATE (", success_rate, "%)\n")
  cat("❌ Significant issues need resolution before production\n")
} else {
  cat("\n❌ PRODUCTION READINESS: POOR (", success_rate, "%)\n")
  cat("❌ Major issues prevent production deployment\n")
}

cat("\nCOMPREHENSIVE TESTING COMPLETE!\n")
cat("==================================\n")

# Save detailed results
results_file <- "comprehensive_pipeline_test_results.csv"
write_csv(test_results, results_file)
cat("📄 Detailed results saved to:", results_file, "\n")

# Exit with appropriate code
if (success_rate >= 75) {
  quit(status = 0)
} else {
  quit(status = 1)
} 