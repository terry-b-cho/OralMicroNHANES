#!/usr/bin/env Rscript
# =============================================================================
# NHANES Oral Microbiome Data Summary Generator - Enhanced Version
# nhanes_oral_microbiome_data_summary.R
# =============================================================================
# This script generates comprehensive publication-ready supplementary tables 
# for NHANES oral microbiome data, including:
# 1. All WAS analysis variables with NHANES codes and database descriptions
# 2. Comprehensive variable type classification and statistics
# 3. Detailed binary variable distribution analysis
# 4. Extended variable coverage beyond WAS analyses
#
# Usage:
#   Rscript scripts/2.1_general_summary_generation/nhanes_oral_microbiome_data_summary.R
#
# Output: Creates comprehensive supplementary tables in results/supplementary_tables/
# =============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(readr)
  library(stringr)
  library(survey)
  library(purrr)
  library(tidyr)
  library(knitr)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

DB_PATH <- "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"
OUTPUT_DIR <- "results/supplementary_tables"
CYCLES <- c("F", "G")  # Focus on cycles F and G only

# WAS analysis configuration files
WAS_CONFIG_FILES <- list(
  "1_demoWAS" = "configs/1_demoWAS_vars.txt",
  "2_oradWAS" = "configs/2_oradWAS_vars.txt", 
  "3_exWAS" = "configs/3_exWAS_vars.txt",
  "4_pheWAS" = "configs/4_pheWAS_vars.txt",
  "5_outWAS" = "configs/5_outWAS_vars.txt",
  "6_zimWAS" = "configs/6_zimWAS_vars.txt"
)

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Load WAS variables from config files
load_was_variables <- function() {
  cat("Loading WAS analysis variables from config files...\n")
  
  was_vars <- list()
  
  for (analysis_name in names(WAS_CONFIG_FILES)) {
    config_file <- WAS_CONFIG_FILES[[analysis_name]]
    if (file.exists(config_file)) {
      vars <- readLines(config_file, warn = FALSE)
      vars <- vars[vars != "" & !startsWith(vars, "#")]  # Remove empty lines and comments
      was_vars[[analysis_name]] <- vars
      cat(sprintf("  %s: %d variables\n", analysis_name, length(vars)))
    } else {
      cat(sprintf("  Warning: %s not found\n", config_file))
    }
  }
  
  # Get all unique variables
  all_was_vars <- unique(unlist(was_vars))
  cat(sprintf("Total unique WAS variables: %d\n", length(all_was_vars)))
  
  return(list(by_analysis = was_vars, all_variables = all_was_vars))
}

# Get participants with oral microbiome data
get_microbiome_participants <- function(con) {
  cat("Identifying participants with oral microbiome data...\n")
  
  microbiome_seqn <- list()
  
  for (cycle in CYCLES) {
    # Check genus relative abundance table (main microbiome data)
    genus_table <- paste0("DADA2RSV_GENUS_RELATIVE_", cycle, "_none")
    
    if (genus_table %in% dbListTables(con)) {
      participants <- tbl(con, genus_table) %>%
        select(SEQN) %>%
        distinct() %>%
        collect() %>%
        mutate(
          SEQN = as.character(SEQN),
          cycle = cycle
        )
      
      microbiome_seqn[[cycle]] <- participants
      cat(sprintf("  Cycle %s: %d participants with microbiome data\n", 
                  cycle, nrow(participants)))
    }
  }
  
  if (length(microbiome_seqn) > 0) {
    combined <- bind_rows(microbiome_seqn)
    cat(sprintf("Total unique participants with microbiome data: %d\n", 
                length(unique(combined$SEQN))))
    return(combined)
  } else {
    stop("No microbiome data found in database")
  }
}

# Get variable metadata from database
get_database_metadata <- function(con, variables) {
  cat("Extracting variable metadata from database...\n")
  
  # Get all available tables
  all_tables <- dbListTables(con)
  relevant_tables <- all_tables[str_detect(all_tables, paste0("(", paste(CYCLES, collapse = "|"), ")$"))]
  
  metadata_list <- list()
  
  for (table_name in relevant_tables) {
    tryCatch({
      # Get column info for this table
      table_info <- dbGetQuery(con, sprintf("PRAGMA table_info(%s)", table_name))
      
      if (nrow(table_info) > 0) {
        table_metadata <- tibble(
          variable = table_info$name,
          table_source = table_name,
          sql_type = table_info$type
        ) %>%
          filter(variable %in% variables, variable != "SEQN")
        
        if (nrow(table_metadata) > 0) {
          metadata_list[[table_name]] <- table_metadata
        }
      }
    }, error = function(e) {
      cat(sprintf("Warning: Could not get metadata for table %s: %s\n", table_name, e$message))
    })
  }
  
  if (length(metadata_list) > 0) {
    combined_metadata <- bind_rows(metadata_list)
    
    # Remove duplicates, keeping first occurrence
    unique_metadata <- combined_metadata %>%
      group_by(variable) %>%
      slice_head(n = 1) %>%
      ungroup()
    
    cat(sprintf("Found metadata for %d/%d variables\n", 
                nrow(unique_metadata), length(variables)))
    
    return(unique_metadata)
  } else {
    return(tibble(variable = character(), table_source = character(), sql_type = character()))
  }
}

# Classify variable types based on data and metadata
classify_variable_type <- function(data, variable, sql_type = "REAL") {
  if (!variable %in% names(data)) {
    return("unknown")
  }
  
  values <- data[[variable]][!is.na(data[[variable]])]
  
  if (length(values) == 0) {
    return("unknown")
  }
  
  # Check if binary (only 0/1 or 1/2 values)
  unique_vals <- unique(values)
  if (length(unique_vals) <= 2) {
    if (all(unique_vals %in% c(0, 1)) || all(unique_vals %in% c(1, 2))) {
      return("binary")
    }
  }
  
  # Check if categorical (limited unique values or character data)
  if (is.character(values) || length(unique_vals) <= 10) {
    return("categorical")
  }
  
  # Check if ordinal (sequential integers)
  if (all(values == round(values)) && length(unique_vals) <= 20) {
    return("ordinal")
  }
  
  # Otherwise continuous
  return("continuous")
}

# Collect comprehensive data for all variables
collect_comprehensive_data <- function(con, microbiome_participants, variables) {
  cat("Collecting comprehensive data for all variables...\n")
  
  # Start with demographic data (includes survey weights)
  demo_data <- list()
  for (cycle in CYCLES) {
    demo_table <- paste0("DEMO_", cycle)
    if (demo_table %in% dbListTables(con)) {
      demo <- tbl(con, demo_table) %>%
        collect() %>%
        mutate(SEQN = as.character(SEQN), cycle = cycle)
      
      demo_data[[cycle]] <- demo
      cat(sprintf("  %s: %d participants, %d total variables\n", 
                  demo_table, nrow(demo), ncol(demo)))
    }
  }
  
  combined_data <- bind_rows(demo_data)
  
  # Collect data from all relevant tables
  other_tables <- dbListTables(con)
  relevant_patterns <- paste0("(", paste(CYCLES, collapse = "|"), ")$")
  relevant_tables <- other_tables[str_detect(other_tables, relevant_patterns)]
  
  for (table_name in relevant_tables) {
    if (!str_detect(table_name, "DEMO_")) {  # Skip demo tables (already processed)
      tryCatch({
        # Get all data from table
        table_data <- tbl(con, table_name) %>%
          collect() %>%
          mutate(SEQN = as.character(SEQN))
        
        # Find variables of interest in this table
        available_vars <- intersect(variables, names(table_data))
        if (length(available_vars) > 0) {
          keep_vars <- c("SEQN", available_vars)
          table_data <- table_data[, names(table_data) %in% keep_vars, drop = FALSE]
          
          cat(sprintf("  %s: %d variables found (%s)\n", 
                      table_name, length(available_vars),
                      paste(head(available_vars, 3), collapse = ", ")))
          
          # Merge with combined data
          combined_data <- combined_data %>%
            left_join(table_data, by = "SEQN", suffix = c("", "_new"))
          
          # Handle duplicate columns
          for (var in available_vars) {
            new_col <- paste0(var, "_new")
            if (new_col %in% names(combined_data)) {
              combined_data[[var]] <- coalesce(combined_data[[var]], 
                                               combined_data[[new_col]])
              combined_data[[new_col]] <- NULL
            }
          }
        }
      }, error = function(e) {
        cat(sprintf("Warning: Could not process table %s: %s\n", 
                    table_name, e$message))
      })
    }
  }
  
  # Filter to only microbiome participants
  microbiome_seqn <- unique(microbiome_participants$SEQN)
  final_data <- combined_data %>%
    filter(SEQN %in% microbiome_seqn)
  
  cat(sprintf("Final dataset: %d participants with %d total variables\n", 
              nrow(final_data), ncol(final_data)))
  
  return(final_data)
}

# Generate comprehensive variable statistics
generate_comprehensive_stats <- function(data, variable, db_metadata, was_analysis = NULL) {
  # Get metadata for this variable
  var_metadata <- db_metadata[db_metadata$variable == variable, ]
  
  n_total <- nrow(data)
  n_available <- sum(!is.na(data[[variable]]))
  n_missing <- n_total - n_available
  percent_available <- round(100 * n_available / n_total, 1)
  
  # Determine variable type
  sql_type <- ifelse(nrow(var_metadata) > 0, var_metadata$sql_type[1], "REAL")
  var_type <- classify_variable_type(data, variable, sql_type)
  
  # Generate appropriate statistics
  if (n_available == 0) {
    statistic <- "No data available"
    additional_stats <- ""
  } else {
    values <- data[[variable]][!is.na(data[[variable]])]
    
    if (var_type == "continuous") {
      mean_val <- mean(values)
      sd_val <- sd(values)
      se_val <- sd_val / sqrt(length(values))
      median_val <- median(values)
      q25 <- quantile(values, 0.25)
      q75 <- quantile(values, 0.75)
      
      statistic <- sprintf("%.2f ± %.2f", mean_val, se_val)
      additional_stats <- sprintf("Median: %.2f (IQR: %.2f-%.2f)", median_val, q25, q75)
      
    } else if (var_type == "binary") {
      freq_table <- table(values, useNA = "no")
      unique_vals <- names(freq_table)
      
      # Determine which value represents "positive" case
      if (all(unique_vals %in% c(0, 1))) {
        positive_val <- "1"
      } else if (all(unique_vals %in% c(1, 2))) {
        positive_val <- "1"  # Usually "1" = Yes, "2" = No in NHANES
      } else {
        positive_val <- unique_vals[1]  # Default to first value
      }
      
      positive_n <- ifelse(positive_val %in% names(freq_table), freq_table[positive_val], 0)
      positive_pct <- round(100 * positive_n / n_available, 1)
      
      statistic <- sprintf("Cases: %d (%.1f%%)", positive_n, positive_pct)
      additional_stats <- sprintf("Controls: %d (%.1f%%)", n_available - positive_n, 100 - positive_pct)
      
    } else if (var_type == "categorical" || var_type == "ordinal") {
      freq_table <- table(values, useNA = "no")
      n_categories <- length(freq_table)
      most_common <- names(freq_table)[which.max(freq_table)]
      most_common_n <- max(freq_table)
      most_common_pct <- round(100 * most_common_n / n_available, 1)
      
      statistic <- sprintf("Mode: %s (%d, %.1f%%)", most_common, most_common_n, most_common_pct)
      additional_stats <- sprintf("%d categories total", n_categories)
      
    } else {
      statistic <- "Unknown data type"
      additional_stats <- ""
    }
  }
  
  # Get table source
  table_source <- ifelse(nrow(var_metadata) > 0, var_metadata$table_source[1], "Unknown")
  
  tibble(
    nhanes_code = variable,
    description = ifelse(variable %in% names(data), "Available in dataset", "Not found"),
    was_analysis = ifelse(is.null(was_analysis), "Extended", was_analysis),
    table_source = table_source,
    variable_type = var_type,
    n_total = n_total,
    n_available = n_available,
    n_missing = n_missing,
    percent_available = percent_available,
    statistic = statistic,
    additional_stats = additional_stats
  )
}

# Generate detailed binary variable analysis
generate_binary_analysis <- function(data, binary_vars, db_metadata) {
  cat("Generating detailed binary variable analysis...\n")
  
  binary_results <- list()
  
  for (var in binary_vars) {
    if (var %in% names(data)) {
      values <- data[[var]][!is.na(data[[var]])]
      
      if (length(values) > 0) {
        freq_table <- table(values, useNA = "no")
        unique_vals <- sort(unique(values))
        
        # Determine case/control coding
        if (all(unique_vals %in% c(0, 1))) {
          cases <- sum(values == 1, na.rm = TRUE)
          controls <- sum(values == 0, na.rm = TRUE)
          case_label <- "1 (Yes/Positive)"
          control_label <- "0 (No/Negative)"
        } else if (all(unique_vals %in% c(1, 2))) {
          cases <- sum(values == 1, na.rm = TRUE)
          controls <- sum(values == 2, na.rm = TRUE)
          case_label <- "1 (Yes)"
          control_label <- "2 (No)"
        } else {
          cases <- freq_table[1]
          controls <- sum(freq_table[-1])
          case_label <- paste0(names(freq_table)[1], " (Reference)")
          control_label <- "Other"
        }
        
        total_valid <- cases + controls
        case_pct <- round(100 * cases / total_valid, 1)
        control_pct <- round(100 - case_pct, 1)
        
        # Calculate prevalence and 95% CI
        prevalence <- cases / total_valid
        se_prev <- sqrt(prevalence * (1 - prevalence) / total_valid)
        ci_lower <- max(0, prevalence - 1.96 * se_prev)
        ci_upper <- min(1, prevalence + 1.96 * se_prev)
        
        binary_results[[var]] <- tibble(
          nhanes_code = var,
          total_n = nrow(data),
          valid_n = total_valid,
          missing_n = nrow(data) - total_valid,
          cases_n = cases,
          cases_pct = case_pct,
          controls_n = controls,
          controls_pct = control_pct,
          prevalence = round(prevalence, 3),
          prevalence_ci_lower = round(ci_lower, 3),
          prevalence_ci_upper = round(ci_upper, 3),
          case_label = case_label,
          control_label = control_label
        )
      }
    }
  }
  
  if (length(binary_results) > 0) {
    return(bind_rows(binary_results))
  } else {
    return(tibble())
  }
}

# =============================================================================
# MAIN ANALYSIS
# =============================================================================

cat("=== NHANES Oral Microbiome Data Summary Generator - Enhanced Version ===\n\n")

# Connect to database
cat("Connecting to database...\n")
con <- dbConnect(SQLite(), DB_PATH)

# Load WAS variables
was_variables <- load_was_variables()

# Get microbiome participants
microbiome_participants <- get_microbiome_participants(con)

# Get database metadata for all variables
db_metadata <- get_database_metadata(con, was_variables$all_variables)

# Collect comprehensive data
participant_data <- collect_comprehensive_data(con, microbiome_participants, was_variables$all_variables)

# Close database connection
dbDisconnect(con)

# =============================================================================
# GENERATE COMPREHENSIVE SUPPLEMENTARY TABLES
# =============================================================================

cat("\nGenerating comprehensive supplementary tables...\n")

# Table S1: Complete WAS Variable Summary
cat("Creating Table S1: Complete WAS Variable Summary...\n")

was_stats_list <- list()
for (analysis_name in names(was_variables$by_analysis)) {
  analysis_vars <- was_variables$by_analysis[[analysis_name]]
  
  for (var in analysis_vars) {
    stats <- generate_comprehensive_stats(participant_data, var, db_metadata, analysis_name)
    was_stats_list[[paste(analysis_name, var, sep = "_")]] <- stats
  }
}

table_s1 <- bind_rows(was_stats_list) %>%
  arrange(was_analysis, nhanes_code) %>%
  select(
    `WAS Analysis` = was_analysis,
    `NHANES Code` = nhanes_code,
    `Variable Description` = description,
    `Table Source` = table_source,
    `Variable Type` = variable_type,
    `N Total` = n_total,
    `N Available` = n_available,
    `N Missing` = n_missing,
    `% Available` = percent_available,
    `Primary Statistic` = statistic,
    `Additional Statistics` = additional_stats
  )

write_csv(table_s1, file.path(OUTPUT_DIR, "table_s1_complete_was_variable_summary.csv"))

# Table S2: WAS Analysis Summary by Type
cat("Creating Table S2: WAS Analysis Summary by Type...\n")

analysis_summary <- bind_rows(was_stats_list) %>%
  group_by(was_analysis) %>%
  summarise(
    `Total Variables` = n(),
    `Variables Available (≥50%)` = sum(percent_available >= 50),
    `Variables Available (≥75%)` = sum(percent_available >= 75),
    `Variables Available (≥90%)` = sum(percent_available >= 90),
    `Mean Availability (%)` = round(mean(percent_available), 1),
    `Median Availability (%)` = round(median(percent_available), 1),
    `Continuous Variables` = sum(variable_type == "continuous"),
    `Binary Variables` = sum(variable_type == "binary"),
    `Categorical Variables` = sum(variable_type == "categorical"),
    `Ordinal Variables` = sum(variable_type == "ordinal"),
    .groups = "drop"
  ) %>%
  arrange(desc(`Mean Availability (%)`))

write_csv(analysis_summary, file.path(OUTPUT_DIR, "table_s2_was_analysis_summary.csv"))

# Table S3: Variable Type Distribution
cat("Creating Table S3: Variable Type Distribution...\n")

type_summary <- bind_rows(was_stats_list) %>%
  group_by(variable_type) %>%
  summarise(
    `Count` = n(),
    `Percentage` = round(100 * n() / nrow(bind_rows(was_stats_list)), 1),
    `Mean Availability (%)` = round(mean(percent_available), 1),
    `Variables ≥90% Available` = sum(percent_available >= 90),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

write_csv(type_summary, file.path(OUTPUT_DIR, "table_s3_variable_type_distribution.csv"))

# Table S4: Detailed Binary Variable Analysis
cat("Creating Table S4: Detailed Binary Variable Analysis...\n")

binary_vars <- bind_rows(was_stats_list) %>%
  filter(variable_type == "binary") %>%
  pull(nhanes_code)

if (length(binary_vars) > 0) {
  binary_analysis <- generate_binary_analysis(participant_data, binary_vars, db_metadata)
  
  if (nrow(binary_analysis) > 0) {
    table_s4 <- binary_analysis %>%
      select(
        `NHANES Code` = nhanes_code,
        `Total N` = total_n,
        `Valid N` = valid_n,
        `Missing N` = missing_n,
        `Cases N` = cases_n,
        `Cases %` = cases_pct,
        `Controls N` = controls_n,
        `Controls %` = controls_pct,
        `Prevalence` = prevalence,
        `95% CI Lower` = prevalence_ci_lower,
        `95% CI Upper` = prevalence_ci_upper,
        `Case Definition` = case_label,
        `Control Definition` = control_label
      )
    
    write_csv(table_s4, file.path(OUTPUT_DIR, "table_s4_binary_variable_analysis.csv"))
  }
}

# Summary report
total_participants <- nrow(participant_data)
total_was_vars <- length(was_variables$all_variables)
available_vars <- sum(bind_rows(was_stats_list)$percent_available >= 50)

cat("\n=== COMPREHENSIVE SUMMARY REPORT ===\n")
cat(sprintf("Study Population: %d participants with oral microbiome data\n", total_participants))
cat(sprintf("NHANES Cycles: %s\n", paste(CYCLES, collapse = ", ")))
cat(sprintf("Total WAS variables assessed: %d\n", total_was_vars))
cat(sprintf("Variables with ≥50%% availability: %d (%.1f%%)\n", 
            available_vars, 100 * available_vars / total_was_vars))

for (analysis in names(was_variables$by_analysis)) {
  n_vars <- length(was_variables$by_analysis[[analysis]])
  cat(sprintf("  %s: %d variables\n", analysis, n_vars))
}

cat(sprintf("\nAll comprehensive tables saved to: %s\n", OUTPUT_DIR))
cat("\nComprehensive supplementary tables generated successfully!\n")
cat("\nTables generated:\n")
cat("- Table S1: Complete WAS Variable Summary (with NHANES codes, types, statistics)\n")
cat("- Table S2: WAS Analysis Summary by Type\n") 
cat("- Table S3: Variable Type Distribution\n")
cat("- Table S4: Detailed Binary Variable Analysis\n")
cat("\nFootnotes:\n")
cat("- Statistics shown as mean ± SE for continuous variables\n")
cat("- Binary variables show case/control distributions with 95% CI\n")
cat("- Variable types classified based on data distribution patterns\n")