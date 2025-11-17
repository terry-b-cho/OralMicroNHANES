#!/usr/bin/env Rscript

# =====================================================================================
# Universal WAS Analysis Script - NHANES Survey-Aware Implementation
# universal_was_analysis_debugging_report_information.R
# =====================================================================================
# For each dependent variable this script produces a single .rds file
# containing three broom tables:
#   • pe_tidied   – per-term coefficients (one row per indep_var + covariates)  
#   • pe_glanced  – model-level statistics (one row per indep_var)
#   • rsq         – R^2 / pseudo-R^2 (one row per indep_var)
#
# CRITICAL: This script implements NHANES complex survey design with proper
# stratification, clustering, and weight scaling across multiple cycles.
# Standard lm/glm results are INVALID for NHANES data.
#
# Usage: Rscript universal_was_analysis.R --dependent_var=<VAR> [options]
# =====================================================================================

# Configure parallel processing to respect SLURM allocation
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if(is.na(num_cores)) num_cores <- 1
Sys.setenv(OMP_NUM_THREADS = num_cores)

# =====================================================================================
# Package Requirements Hardening
# =====================================================================================

check_required_packages <- function() {
  required_packages <- c("survey", "broom", "dplyr", "glue", "DBI", "RSQLite", "getopt", "logger", "readr")

  missing <- sapply(required_packages, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  })
  
  if (any(missing)) {
    cat("Missing required packages:", paste(names(missing)[missing], collapse = ", "), "\n")
    cat("Install via conda environment or environment.yml\n")
    cat("Current .libPaths():\n")
    print(.libPaths())
    quit(status=1)
}

  # Load packages
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE, quietly = TRUE)
  }
  
  invisible(TRUE)
}

# Check and load all required packages
check_required_packages()

options(dplyr.threads = num_cores)

# =====================================================================================
# Helper functions for direct database operations
# =====================================================================================

# Check if a variable exists in a table
check_variable_in_table <- function(conn, table_name, var_name) {
  tryCatch({
    col_info <- dbGetQuery(conn, sprintf("PRAGMA table_info(%s);", table_name))
    var_name %in% col_info$name
  }, error = function(e) {
    log_error("Error checking variable {var_name} in table {table_name}: {e$message}")
    FALSE
  })
}

# List all columns in a table
list_table_columns <- function(conn, table_name) {
  tryCatch({
    col_info <- dbGetQuery(conn, sprintf("PRAGMA table_info(%s);", table_name))
    col_info$name
  }, error = function(e) {
    log_error("Error listing columns in table {table_name}: {e$message}")
    character(0)
  })
}

# Check if table exists in the database
check_table_exists <- function(conn, table_name) {
  tryCatch({
    table_name %in% dbListTables(conn)
  }, error = function(e) {
    log_error("Error checking if table {table_name} exists: {e$message}")
    FALSE
  })
}

# Get the series letter (F, G, etc.) from a table name
get_series_from_table_name <- function(table_name) {
  series_match <- regexpr("_([A-Z])($|_)", table_name)
  if(series_match > 0) {
    substr(table_name, series_match + 1, series_match + 1)
  } else {
    NA_character_
  }
}

# Extract the demo table name for a given series
get_demo_table_for_series <- function(conn, series) {
  tryCatch({
    demo_table <- paste0("DEMO_", series)
    if(check_table_exists(conn, demo_table)) {
      return(demo_table)
    } else {
      log_warn("No demographic table found for series {series}")
      return(NULL)
    }
  }, error = function(e) {
    log_error("Error finding demo table for series {series}: {e$message}")
    return(NULL)
  })
}

# Function to clean up memory during processing
clean_memory <- function() {
  gc(full = TRUE)
  invisible(NULL)
}

# Function to estimate memory usage
report_memory <- function(message = "Current memory usage") {
  mem_used <- gc(full = FALSE)
  mem_mb <- sum(mem_used[,2]) / 1024  # Convert to MB
  log_info("{message}: {round(mem_mb, 2)} MB")
  invisible(NULL)
}

# =====================================================================================
# Canonical SQL Template (Drop-in Replacement)
# =====================================================================================

build_sql_select <- function(dependent_var, independent_var, use_covariates, 
                           dependent_table, independent_table, demo_table) {
  
  core_cols <- glue("
    t1.SEQN,
    t1.{dependent_var}  AS dep_var,
    t2.{independent_var} AS indep_var,
    d.WTINT2YR,
    d.WTMEC2YR,
    d.SDMVSTRA,
    d.SDMVPSU,
    d.SDDSRVYR")

  covar_cols <- "
    ,d.RIDAGEYR
    ,d.AGE_SQUARED
    ,d.RIAGENDR
    ,d.INDFMPIR
    ,d.EDUCATION_LESS9
    ,d.EDUCATION_9_11
    ,d.EDUCATION_AA
    ,d.EDUCATION_COLLEGEGRAD
    ,d.ETHNICITY_MEXICAN
    ,d.ETHNICITY_OTHERHISPANIC
    ,d.ETHNICITY_OTHER
    ,d.ETHNICITY_NONHISPANICBLACK
    ,d.BORN_INUSA"

  select_clause <- if (use_covariates) paste0(core_cols, covar_cols) else core_cols

  # Assemble final query
  sql_query <- glue("
    SELECT {select_clause}
    FROM  {dependent_table}  AS t1
    INNER JOIN {independent_table} AS t2 ON t1.SEQN = t2.SEQN
    INNER JOIN {demo_table} AS d      ON t1.SEQN = d.SEQN")
  
  return(sql_query)
}

# =====================================================================================
# Pre-Modelling Preprocessing Pipeline
# =====================================================================================

prepare_model_data <- function(df, covars, use_mec = FALSE, regression_type = "linear") {
  # 1. Keep complete cases across *all* modelling and design variables
  must_have <- c("dep_var", "indep_var",
                 if (length(covars)) covars,
                 "SDMVSTRA", "SDMVPSU",        # design
                 if (use_mec) "WTMEC2YR" else "WTINT2YR",
                 "SDDSRVYR")
  
  df <- df[complete.cases(df[, must_have]), , drop = FALSE]
  
  if (nrow(df) == 0) {
    stop("No complete cases found after filtering")
  }
  
  # 2. Choose correct weight and scale by #cycles
  wt_col <- if (use_mec) "WTMEC2YR" else "WTINT2YR"
  n_cycles <- length(unique(df$SDDSRVYR))
  df$weight <- df[[wt_col]] / n_cycles     # new unified column 'weight'
  
  # 3. Validate binary outcomes for logistic regression
  if (regression_type == "logistic" && 
      !all(df$dep_var %in% c(0, 1, NA))) {
    stop(sprintf("Logistic regression requested but dep_var contains non-binary values: %s", 
                paste(unique(df$dep_var), collapse = ", ")))
  }
  
  df
}

# =====================================================================================
# Survey-Aware Model Fitting
# =====================================================================================

fit_survey_model <- function(model_data, formula_str, regression_type) {
  
  # Build survey design
  design <- svydesign(
    id      = ~SDMVPSU,
    strata  = ~SDMVSTRA,
    weights = ~weight,      # unified weight column
    nest    = TRUE,
    data    = model_data)
  
  # Choose family
  family_choice <- if (regression_type == "linear") gaussian() else quasibinomial()
  
  # Fit model
  fit <- svyglm(as.formula(formula_str), design = design, family = family_choice)
  
  return(fit)
}

# =====================================================================================
# Results Extraction with Broom
# =====================================================================================

extract_survey_results <- function(fit, regression_type) {
  
  # Extract coefficient table
  tidied <- broom::tidy(fit)
  
  # Extract model-level statistics
  glanced <- broom::glance(fit)
  
  # Calculate R-squared manually for survey models
  if (regression_type == "linear") {
    # For linear survey models, calculate R-squared manually
    tryCatch({
      # Method 1: Use survey-weighted R-squared calculation
      design <- fit$survey.design
      
      # Get fitted values and residuals
      fitted_vals <- fitted(fit)
      residuals_vals <- residuals(fit)
      observed_vals <- fitted_vals + residuals_vals
      
      # Calculate weighted means
      weights <- weights(design, "sampling")
      y_mean <- weighted.mean(observed_vals, weights, na.rm = TRUE)
      
      # Total sum of squares (weighted)
      tss <- sum(weights * (observed_vals - y_mean)^2, na.rm = TRUE)
      
      # Residual sum of squares (weighted)
      rss <- sum(weights * residuals_vals^2, na.rm = TRUE)
      
      # R-squared calculation
      r_squared <- ifelse(tss > 0, 1 - (rss / tss), 0)
      
      # Adjusted R-squared (approximate for survey data)
      n <- nobs(fit)
      p <- length(coef(fit)) - 1  # number of predictors (excluding intercept)
      adj_r_squared <- ifelse(n > p + 1, 1 - ((1 - r_squared) * (n - 1) / (n - p - 1)), NA_real_)
      
      rsq <- tibble(
        r.squared = as.numeric(pmax(0, pmin(1, r_squared))),  # Bound between 0 and 1
        adj.r.squared = as.numeric(adj_r_squared)
      )
    }, error = function(e) {
      # Fallback: simple R-squared approximation
      rsq <- tibble(
        r.squared = NA_real_,
        adj.r.squared = NA_real_
      )
    })
  } else {
    # Pseudo R-squared for logistic survey models
    rsq <- tibble(
      r.squared = 1 - fit$deviance / fit$null.deviance,
      adj.r.squared = NA_real_
    )
  }
  
  return(list(tidied = tidied, glanced = glanced, rsq = rsq))
}

# =====================================================================================
# Command line arguments
# =====================================================================================

# Define command line options
spec <- matrix(c(
  # Required arguments
  'dependent_var',          'd', 1, "character", "Dependent variable name",
  'schema_structure_file',  's', 1, "character", "File with schema structure information",
  'database_path',          'b', 1, "character", "Path to SQLite database",
  'output_path',            'o', 1, "character", "Path to save output files",
  
  # Optional arguments
  'analysis_type',          'a', 2, "character", "Type of analysis (e.g., 1_demoWAS)",
  'regression_type',        'r', 2, "character", "Regression type (linear or logistic)",
  'use_covariates',         'c', 2, "character", "Whether to include covariates (TRUE/FALSE)",
  'covariate_model',        'm', 2, "character", "Covariate model type (base, full, etc.)",
  'quantile_analysis',      'q', 2, "character", "Use quantile transformation for independent vars",
  'sample_size_threshold',  't', 2, "integer",   "Minimum sample size for analysis",
  'scale_variables',        'v', 2, "character", "Scale variables in regression",
  'transform_dependent',    'p', 2, "character", "Log transform the dependent variable",
  'transform_independent',  'e', 2, "character", "Log transform the independent variable",
  'debug_mode',             'g', 2, "character", "Run in debug mode with extra logging",
  'test_mode',              'x', 2, "character", "Run in test mode with limited variables",
  'batch_size',             'k', 2, "integer",   "Number of variables to process in each batch"
), byrow=TRUE, ncol=5)

# Parse command line arguments
opt <- getopt(spec)

# Check for required arguments
required_args <- c('dependent_var', 'schema_structure_file', 'database_path', 'output_path')
missing_args <- required_args[!required_args %in% names(opt)]

if(length(missing_args) > 0) {
  cat("Missing required arguments:", paste(missing_args, collapse=", "), "\n")
  cat("Usage: Rscript universal_was_analysis.R -d dependent_var -s schema_file -b db_path -o output_path [options]\n")
  quit(status=1)
}

# Convert string parameters to appropriate types
parse_bool <- function(x, default = FALSE) {
  if(is.null(x)) return(default)
  if(is.logical(x)) return(x)
  if(is.numeric(x)) return(as.logical(x))
  
  x_upper <- toupper(as.character(x))
  if(x_upper %in% c("TRUE", "T", "YES", "Y", "1", "TEST")) return(TRUE)
  if(x_upper %in% c("FALSE", "F", "NO", "N", "0")) return(FALSE)
  
  # Default for unrecognized values
  log_warn("Unrecognized boolean value: {x}, using default: {default}")
  return(default)
}

# Set default values for optional parameters
if(is.null(opt$analysis_type))          opt$analysis_type <- 'exWAS'
if(is.null(opt$regression_type))        opt$regression_type <- 'linear'
if(is.null(opt$use_covariates))         opt$use_covariates <- TRUE else opt$use_covariates <- parse_bool(opt$use_covariates, TRUE)
if(is.null(opt$covariate_model))        opt$covariate_model <- 'age_sex_ethnicity_income_education'
if(is.null(opt$quantile_analysis))      opt$quantile_analysis <- FALSE else opt$quantile_analysis <- parse_bool(opt$quantile_analysis, FALSE)
if(is.null(opt$sample_size_threshold))  opt$sample_size_threshold <- 0
if(is.null(opt$scale_variables))        opt$scale_variables <- TRUE else opt$scale_variables <- parse_bool(opt$scale_variables, TRUE)
if(is.null(opt$transform_dependent))    opt$transform_dependent <- FALSE else opt$transform_dependent <- parse_bool(opt$transform_dependent, FALSE)
if(is.null(opt$transform_independent))  opt$transform_independent <- FALSE else opt$transform_independent <- parse_bool(opt$transform_independent, FALSE)
if(is.null(opt$debug_mode))             opt$debug_mode <- FALSE else opt$debug_mode <- parse_bool(opt$debug_mode, FALSE)
if(is.null(opt$test_mode))              opt$test_mode <- FALSE else opt$test_mode <- parse_bool(opt$test_mode, FALSE)
if(is.null(opt$batch_size))             opt$batch_size <- 50

# Override covariate settings based on analysis type
if (startsWith(opt$analysis_type, "1_demoWAS")) {
  opt$use_covariates <- FALSE
}

# Set regression type based on analysis type if not explicitly set
if (startsWith(opt$analysis_type, "2_oradWAS") || startsWith(opt$analysis_type, "5_outWAS")) {
  opt$regression_type <- "logistic"
} else {
  opt$regression_type <- "linear"
}

# Set up logging
log_threshold(INFO)
if(opt$debug_mode) {
  log_threshold(DEBUG)
  log_debug("Debug mode enabled")
}

# Initialize the analysis environment
log_info("Initializing analysis for dependent variable: {opt$dependent_var}")
log_info("Analysis type: {opt$analysis_type}")
log_info("Regression type: {opt$regression_type}")
log_info("Process ID: {Sys.getpid()}")
log_info("Using {num_cores} CPU cores")

# =====================================================================================
# Database connection and data loading
# =====================================================================================

# Connect to the database
log_info("Connecting to database at {opt$database_path}")
con <- tryCatch({
  DBI::dbConnect(RSQLite::SQLite(), dbname = opt$database_path)
}, error = function(e) {
  log_error("Failed to connect to database: {e$message}")
  quit(status=1)
})

# Read schema structure file for this dependent variable
log_info("Reading schema structure file: {opt$schema_structure_file}")
schema_data <- tryCatch({
  read_csv(opt$schema_structure_file, show_col_types = FALSE)
}, error = function(e) {
  log_error("Failed to read schema structure file: {e$message}")
  DBI::dbDisconnect(con)
  quit(status=1)
})

# Check if there are any variable pairs to analyze
if(nrow(schema_data) == 0) {
  log_warn("No variable pairs found for analysis of {opt$dependent_var}, quitting")
  DBI::dbDisconnect(con)
  quit(status=0)
}

# Log the column names for debugging
log_debug("Schema columns: {paste(names(schema_data), collapse=', ')}")

# Determine dependent and independent variable columns in schema
dep_col <- NULL
indep_col <- NULL
if("dep_var" %in% names(schema_data)) dep_col <- "dep_var"
if("indep_var" %in% names(schema_data)) indep_col <- "indep_var"

if(is.null(dep_col) || is.null(indep_col)) {
  log_error("Schema must contain 'dep_var' and 'indep_var' columns")
  DBI::dbDisconnect(con)
  quit(status=1)
}

# Filter for the specific dependent variable
schema_data_filtered <- schema_data %>%
  filter(.data[[dep_col]] == opt$dependent_var)

# Filter by sample size threshold
if (opt$sample_size_threshold > 0 && "n" %in% colnames(schema_data_filtered)) {
  schema_data_filtered <- schema_data_filtered %>% filter(n >= opt$sample_size_threshold)
}

# Limit number of pairs in test mode
if(opt$test_mode) {
  log_info("Running in test mode - limiting to 5 variable pairs")
  set.seed(123)
  schema_data_filtered <- schema_data_filtered %>% sample_n(min(5, nrow(schema_data_filtered)))
}

# Log the number of variable pairs
num_independent_vars <- nrow(schema_data_filtered)
log_info("Found {num_independent_vars} independent variables for analysis with {opt$dependent_var}")

# Report memory usage after loading data
report_memory("Memory usage after loading data")

# =====================================================================================
# Run the analysis - Main Loop Refactoring
# =====================================================================================

log_info("Starting analysis for {opt$dependent_var} with {num_independent_vars} independent variables")

# Create a list to store models
models_list <- vector("list", length = num_independent_vars)

# Process in batches to manage memory better
batch_size <- min(opt$batch_size, num_independent_vars)
num_batches <- ceiling(num_independent_vars / batch_size)

log_info("Processing in {num_batches} batches of up to {batch_size} variables each")

for(batch in 1:num_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, num_independent_vars)
  
  log_info("Processing batch {batch}/{num_batches} - variables {start_idx} to {end_idx}")
  
  # Loop through variables in this batch
  for(i in start_idx:end_idx) {
    tryCatch({
      # Get variable information
      curr_pair <- schema_data_filtered[i, ]
      dependent_var <- curr_pair[[dep_col]]
      independent_var <- curr_pair[[indep_col]]
      
      # Get table names
      dependent_table <- curr_pair$dep_table
      independent_table <- curr_pair$indep_table
      
      # Log progress
      log_info("[{i}/{num_independent_vars}] Analyzing {dependent_var} ~ {independent_var}")
      log_debug("Tables: {dependent_table} ~ {independent_table}")
      
      # Get series information
      series1 <- get_series_from_table_name(dependent_table)
      series2 <- get_series_from_table_name(independent_table)
      
      # Use the first valid series
      series <- if(!is.na(series1)) series1 else if(!is.na(series2)) series2 else "F"
      log_info("Series determined from table names: {series}")
      
      # Get demo table for this series
      demo_table <- get_demo_table_for_series(con, series)
      if(is.null(demo_table)) {
        log_warn("No demographic table found for series {series}. Skipping.")
        models_list[[i]] <- list(error = "No demographic table found")
        next
      }
      
      # First, check if tables and variables exist
      if(!check_table_exists(con, dependent_table)) {
        log_warn("Dependent variable table {dependent_table} not found. Skipping.")
        models_list[[i]] <- list(error = "Dependent variable table not found")
        next
      }
      
      if(!check_table_exists(con, independent_table)) {
        log_warn("Independent variable table {independent_table} not found. Skipping.")
        models_list[[i]] <- list(error = "Independent variable table not found")
        next
      }
      
      # 1. Build and execute SQL query using new template
      sql_query <- build_sql_select(
        dependent_var = dependent_var,
        independent_var = independent_var,
        use_covariates = opt$use_covariates,
        dependent_table = dependent_table,
        independent_table = independent_table,
        demo_table = demo_table
      )
      
      # Execute the query
      raw_data <- tryCatch({
        dbGetQuery(con, sql_query)
      }, error = function(e) {
        log_error("SQL query failed: {e$message}")
        log_debug("Query: {sql_query}")
        return(NULL)
      })
      
      if(is.null(raw_data) || nrow(raw_data) == 0) {
        log_warn("No data retrieved for this variable pair. Skipping.")
        models_list[[i]] <- list(error = "No data retrieved")
        next
      }
      
      log_info("Retrieved {nrow(raw_data)} rows of data")
      
      # Define available covariates
      covar_cols <- c(
        "RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "INDFMPIR", 
        "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
        "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", 
        "ETHNICITY_NONHISPANICBLACK", "BORN_INUSA"
      )
      
      available_covars <- covar_cols[covar_cols %in% names(raw_data)]
      
      # 2. Preprocess data with survey design
      use_mec <- if("needs_mec_weight" %in% names(curr_pair)) isTRUE(curr_pair$needs_mec_weight) else FALSE
      model_data <- tryCatch({
        prepare_model_data(raw_data, available_covars, use_mec, opt$regression_type)
      }, error = function(e) {
        log_error("Data preprocessing failed: {e$message}")
        return(NULL)
      })
      
      if(is.null(model_data) || nrow(model_data) == 0) {
        log_warn("No complete cases for {dependent_var} vs {independent_var}")
        models_list[[i]] <- list(error = "No complete cases after preprocessing")
        next
      }
      
      # Store original values before transforms
      orig_dep_mean <- mean(model_data$dep_var, na.rm = TRUE)
      orig_dep_sd <- sd(model_data$dep_var, na.rm = TRUE)
      orig_indep_mean <- mean(model_data$indep_var, na.rm = TRUE)
      orig_indep_sd <- sd(model_data$indep_var, na.rm = TRUE)
      
      # Apply transformations if requested
      if(opt$transform_dependent && opt$regression_type != "logistic") {
        # For log transform, ensure data is positive
        if(all(model_data$dep_var >= 0, na.rm = TRUE)) {
          model_data$dep_var <- log10(model_data$dep_var + 1)
          log_debug("Applied log10(x+1) transformation to dependent variable")
        }
      }
      
      if(opt$transform_independent) {
        # For log transform, ensure data is positive
        if(all(model_data$indep_var >= 0, na.rm = TRUE)) {
          model_data$indep_var <- log10(model_data$indep_var + 1)
          log_debug("Applied log10(x+1) transformation to independent variable")
        }
      }
      
      # Scale variables if requested
      if(opt$scale_variables) {
        # Scale only if we have variation
        if(sd(model_data$indep_var, na.rm = TRUE) > 0) {
          model_data$indep_var <- as.vector(scale(model_data$indep_var))
          log_debug("Scaled independent variable")
        }
        
        if(opt$regression_type == "linear" && sd(model_data$dep_var, na.rm = TRUE) > 0) {
          model_data$dep_var <- as.vector(scale(model_data$dep_var))
          log_debug("Scaled dependent variable")
        }
      }
      
      # 3. Build formula
      formula_str <- "dep_var ~ indep_var"
      
      if(opt$use_covariates && length(available_covars) > 0) {
        log_info("Using {length(available_covars)} available covariates")
        formula_str <- paste0(formula_str, " + ", paste(available_covars, collapse = " + "))
      }
      
      log_debug("Formula: {formula_str}")
      
      # ------------------------------------------------------------------
      # CONTROLLED DEBUG DUMP - Only if debug mode is enabled
      # ------------------------------------------------------------------
      if(opt$debug_mode) {
      tryCatch({
        dbg_dir  <- file.path(opt$output_path, "debug_dumps")
        dir.create(dbg_dir, showWarnings = FALSE, recursive = TRUE)

        dbg_file <- file.path(
          dbg_dir,
          sprintf("%s__%s__%s.txt",
                  gsub("[^[:alnum:]_]", "", dependent_var),
                  gsub("[^[:alnum:]_]", "", independent_var),
                  format(Sys.time(), "%Y%m%d%H%M%S"))
        )

        # Create file and write initial content
        file.create(dbg_file)
        
        ## helper that appends nicely with error handling
        append_ln <- function(..., sep = " ") {
          tryCatch({
            cat(..., "\n", file = dbg_file, sep = sep, append = TRUE)
          }, error = function(e) {
            log_error("Error writing to debug file: {e$message}")
          })
        }

        ## helper to append dput objects
        write_dput <- function(obj) {
            con_file <- file(dbg_file, open = "a")
            on.exit(close(con_file))
            dput(obj, con_file)
            cat("\n", file = con_file)      # tidy newline
        }

          # Debug information collection
          append_ln("===== ANALYSIS METADATA =====")
          append_ln("Analysis Type: ", opt$analysis_type)
          append_ln("Regression Type: ", opt$regression_type)
          append_ln("Dependent Variable: ", dependent_var)
          append_ln("Independent Variable: ", independent_var)
          append_ln("Use Covariates: ", opt$use_covariates)
          append_ln("Scale Variables: ", opt$scale_variables)
          append_ln("Transform Dependent: ", opt$transform_dependent)
          append_ln("Transform Independent: ", opt$transform_independent)
          append_ln("Series: ", series)
          append_ln("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
          append_ln()

        append_ln("===== SQL QUERY =====")
        append_ln(sql_query, "\n")

          # Schema of every table touched
        tryCatch({
          tb_list <- unique(c(dependent_table, independent_table, demo_table))
          for (tb in tb_list[!is.na(tb_list)]) {
            append_ln("===== PRAGMA table_info(", tb, ") =====")
            capture.output(
              dbGetQuery(con, sprintf("PRAGMA table_info(%s);", tb)),
              file = dbg_file, append = TRUE
            )
            append_ln()
          }
        }, error = function(e) {
          append_ln("Error getting table info: ", e$message)
        })

          # Raw data diagnostics
          append_ln("===== RAW DATA DIAGNOSTICS =====")
          append_ln("Raw data rows: ", nrow(raw_data))
          append_ln("Raw data columns: ", ncol(raw_data))
          append_ln("Raw data column names: ", paste(names(raw_data), collapse = ", "))
          
          # Check for missing values in key columns
          key_cols <- c("dep_var", "indep_var", "WTINT2YR", "WTMEC2YR", "SDMVSTRA", "SDMVPSU", "SDDSRVYR")
          for(col in key_cols) {
            if(col %in% names(raw_data)) {
              na_count <- sum(is.na(raw_data[[col]]))
              append_ln(sprintf("Missing values in %s: %d (%.2f%%)", col, na_count, 100*na_count/nrow(raw_data)))
            }
          }
          append_ln()

          # Data-frame overview before modelling
          append_ln("===== PROCESSED MODEL DATA =====")
          append_ln("Processed data rows: ", nrow(model_data))
          append_ln("Processed data columns: ", ncol(model_data))
          
        tryCatch({
            append_ln("===== str(model_data) =====")
          capture.output(str(model_data), file = dbg_file, append = TRUE)
        }, error = function(e) {
          append_ln("Error in str(model_data): ", e$message)
        })

        append_ln("\n===== Summary(model_data) =====")
        tryCatch({
          capture.output(summary(model_data), file = dbg_file, append = TRUE)
        }, error = function(e) {
          append_ln("Error in summary(model_data): ", e$message)
        })

        append_ln("\n===== NA COUNTS =====")
        tryCatch({
          capture.output(colSums(is.na(model_data)), file = dbg_file, append = TRUE)
        }, error = function(e) {
          append_ln("Error in NA counts: ", e$message)
        })

          # Weight diagnostics
          append_ln("\n===== WEIGHT DIAGNOSTICS =====")
          tryCatch({
            wt_cols <- grep("^WT", names(model_data), value = TRUE)
            append_ln("Weight columns detected: ", if (length(wt_cols)) paste(wt_cols, collapse = ", ") else "NONE")
            append_ln("Number of cycles detected: ", n_cycles)
            append_ln("Weight scaling factor: ", 1/n_cycles)
            append_ln("Unified weight column range: ", paste(range(model_data$weight, na.rm = TRUE), collapse = " - "))
            append_ln("Unified weight column mean: ", round(mean(model_data$weight, na.rm = TRUE), 2))
            append_ln("Unified weight column sum: ", round(sum(model_data$weight, na.rm = TRUE), 2))
            
            # Check for zero or negative weights
            zero_weights <- sum(model_data$weight <= 0, na.rm = TRUE)
            if(zero_weights > 0) {
              append_ln("WARNING: ", zero_weights, " zero or negative weights detected")
            }
          }, error = function(e) {
            append_ln("Error in weight diagnostics: ", e$message)
          })

          # Survey design diagnostics
          append_ln("\n===== SURVEY DESIGN DIAGNOSTICS =====")
        tryCatch({
            # Check strata and PSU distributions
            strata_counts <- table(model_data$SDMVSTRA)
            psu_counts <- table(model_data$SDMVPSU)
            append_ln("Number of strata: ", length(strata_counts))
            append_ln("Strata sizes (min, max, mean): ", 
                     paste(c(min(strata_counts), max(strata_counts), round(mean(strata_counts), 1)), collapse = ", "))
            append_ln("Number of PSUs: ", length(psu_counts))
            append_ln("PSU sizes (min, max, mean): ", 
                     paste(c(min(psu_counts), max(psu_counts), round(mean(psu_counts), 1)), collapse = ", "))
            
            # Check for singleton PSUs (problematic for variance estimation)
            singleton_strata <- names(strata_counts)[strata_counts == 1]
            if(length(singleton_strata) > 0) {
              append_ln("WARNING: ", length(singleton_strata), " singleton strata detected: ", paste(singleton_strata, collapse = ", "))
          }
        }, error = function(e) {
            append_ln("Error in survey design diagnostics: ", e$message)
        })

          # Variable distribution diagnostics
          append_ln("\n===== VARIABLE DISTRIBUTIONS =====")
        tryCatch({
            # Dependent variable
            append_ln("Dependent variable (", dependent_var, ") statistics:")
            append_ln("  Original mean: ", round(orig_dep_mean, 4))
            append_ln("  Original SD: ", round(orig_dep_sd, 4))
            append_ln("  Processed mean: ", round(mean(model_data$dep_var, na.rm = TRUE), 4))
            append_ln("  Processed SD: ", round(sd(model_data$dep_var, na.rm = TRUE), 4))
            append_ln("  Range: ", paste(round(range(model_data$dep_var, na.rm = TRUE), 4), collapse = " to "))
            
            # For binary outcomes, check distribution
            if(opt$regression_type == "logistic") {
              dep_table <- table(model_data$dep_var)
              append_ln("  Binary distribution: ", paste(names(dep_table), "=", dep_table, collapse = ", "))
              if(length(dep_table) == 2) {
                prevalence <- dep_table["1"] / sum(dep_table)
                append_ln("  Prevalence: ", round(prevalence * 100, 2), "%")
              }
            }
            
            # Independent variable
            append_ln("Independent variable (", independent_var, ") statistics:")
            append_ln("  Original mean: ", round(orig_indep_mean, 4))
            append_ln("  Original SD: ", round(orig_indep_sd, 4))
            append_ln("  Processed mean: ", round(mean(model_data$indep_var, na.rm = TRUE), 4))
            append_ln("  Processed SD: ", round(sd(model_data$indep_var, na.rm = TRUE), 4))
            append_ln("  Range: ", paste(round(range(model_data$indep_var, na.rm = TRUE), 4), collapse = " to "))
        }, error = function(e) {
            append_ln("Error in variable distribution diagnostics: ", e$message)
        })

          # Covariate information
          append_ln("\n===== COVARIATE INFORMATION =====")
          append_ln("Available covariates: ", paste(available_covars, collapse = ", "))
          append_ln("Number of covariates used: ", length(available_covars))
          append_ln("Formula: ", formula_str)
          append_ln()

          # Model fitting diagnostics (after model is fitted)
          if(!is.null(fit)) {
            append_ln("===== MODEL FITTING DIAGNOSTICS =====")
        tryCatch({
              append_ln("Model class: ", class(fit)[1])
              append_ln("Model family: ", fit$family$family)
              append_ln("Model link: ", fit$family$link)
              append_ln("Number of observations: ", nobs(fit))
              append_ln("Number of parameters: ", length(coef(fit)))
              append_ln("Degrees of freedom: ", fit$df.residual)
              
              # Check for convergence issues
              if(!is.null(fit$converged)) {
                append_ln("Model converged: ", fit$converged)
              }
              
              # Check for fitted values issues
              fitted_vals <- fitted(fit)
              append_ln("Fitted values range: ", paste(round(range(fitted_vals, na.rm = TRUE), 4), collapse = " to "))
              
              # For logistic models, check for separation issues
              if(opt$regression_type == "logistic") {
                extreme_fitted <- sum(fitted_vals < 0.01 | fitted_vals > 0.99, na.rm = TRUE)
                if(extreme_fitted > 0) {
                  append_ln("WARNING: ", extreme_fitted, " extreme fitted values (< 0.01 or > 0.99) - possible separation")
                }
              }
              
              # Check coefficient magnitudes
              coefs <- coef(fit)
              large_coefs <- sum(abs(coefs) > 10, na.rm = TRUE)
              if(large_coefs > 0) {
                append_ln("WARNING: ", large_coefs, " large coefficients (|coef| > 10) detected")
              }
              
        }, error = function(e) {
              append_ln("Error in model diagnostics: ", e$message)
            })
          }

          # Results validation
          if(!is.null(results)) {
            append_ln("\n===== RESULTS VALIDATION =====")
            tryCatch({
              append_ln("Tidied results rows: ", nrow(results$tidied))
              append_ln("Glanced results rows: ", nrow(results$glanced))
              append_ln("R-squared results rows: ", nrow(results$rsq))
              
              # Check for missing p-values or extreme values
              if(nrow(results$tidied) > 0) {
                na_pvals <- sum(is.na(results$tidied$p.value))
                extreme_pvals <- sum(results$tidied$p.value < 1e-10 | results$tidied$p.value > 1, na.rm = TRUE)
                append_ln("Missing p-values: ", na_pvals)
                append_ln("Extreme p-values (< 1e-10 or > 1): ", extreme_pvals)
                
                # Check standard errors
                na_se <- sum(is.na(results$tidied$std.error))
                zero_se <- sum(results$tidied$std.error <= 0, na.rm = TRUE)
                append_ln("Missing standard errors: ", na_se)
                append_ln("Zero or negative standard errors: ", zero_se)
              }
              
              # R-squared validation
              if(nrow(results$rsq) > 0) {
                rsq_val <- results$rsq$r.squared[1]
                append_ln("R-squared value: ", round(rsq_val, 6))
                if(rsq_val < 0 || rsq_val > 1) {
                  append_ln("WARNING: R-squared outside [0,1] range")
                }
              }
              
            }, error = function(e) {
              append_ln("Error in results validation: ", e$message)
            })
          }

          # Session & memory
          append_ln("\n===== SESSION INFORMATION =====")
        tryCatch({
          capture.output(sessionInfo(), file = dbg_file, append = TRUE)
        }, error = function(e) {
          append_ln("Error in sessionInfo: ", e$message)
        })

        mem_mb <- sum(gc(full = FALSE)[, 2]) / 1024
          append_ln("\n===== MEMORY USAGE =====")
          append_ln("Memory used (MB): ", round(mem_mb, 2))
          append_ln("R max memory (MB): ", round(as.numeric(object.size(model_data)) / 1024^2, 2))

          # Final summary
          append_ln("\n===== DEBUG SUMMARY =====")
          append_ln("Analysis completed successfully: ", !is.null(results))
          append_ln("Model fitted successfully: ", !is.null(fit))
          append_ln("Results extracted successfully: ", !is.null(results))
          append_ln("Debug dump completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

      }, error = function(e) {
        log_error("Error in debug dump: {e$message}")
      })
      }
      # ------------------------------------------------------------------
      
      # 4. Fit survey-aware model
      fit <- tryCatch({
        fit_survey_model(model_data, formula_str, opt$regression_type)
      }, error = function(e) {
        log_error("Model fitting failed: {e$message}")
        return(NULL)
      })
      
      if(is.null(fit)) {
        models_list[[i]] <- list(error = "Model fitting failed")
        next
      }
      
      # 5. Extract results
      results <- tryCatch({
        extract_survey_results(fit, opt$regression_type)
      }, error = function(e) {
        log_error("Results extraction failed: {e$message}")
        return(NULL)
      })
      
      if(is.null(results)) {
        models_list[[i]] <- list(error = "Results extraction failed")
        next
      }
      
      # 6. Store in models_list
      models_list[[i]] <- list(
        pair_info = curr_pair,
        model = fit,
        results = results,
        n_obs = nrow(model_data),
          series = series,
          orig_dep_mean = orig_dep_mean,
          orig_dep_sd = orig_dep_sd,
          orig_indep_mean = orig_indep_mean,
          orig_indep_sd = orig_indep_sd
      )
      
      log_info("Successfully fit model for {dependent_var} ~ {independent_var}")
      
    }, error = function(e) {
      log_error("Error in analysis for variables at index {i}: {e$message}")
      models_list[[i]] <- list(error = e$message)
    })
    
    # Periodically clean up memory
    if(i %% 10 == 0) {
      clean_memory()
    }
  }
  
  # Clean up memory after each batch
  report_memory(sprintf("Memory usage after batch %d", batch))
  clean_memory()
  report_memory("Memory usage after cleanup")
}

# =====================================================================================
# Final Results Collation
# =====================================================================================

extract_final_results <- function(models_list) {
  
  # Initialize result containers
  pe_tidied_list <- list()
  pe_glanced_list <- list()
  rsq_list <- list()
  
  for (i in seq_along(models_list)) {
    if (is.null(models_list[[i]]) || !is.null(models_list[[i]]$error)) next
    
    model_info <- models_list[[i]]
    results <- model_info$results
    
    # Add identifying information
    results$tidied$independent_var <- model_info$pair_info[[indep_col]]
    results$tidied$n_obs <- model_info$n_obs
    results$tidied$series <- model_info$series
    results$tidied$phenotype <- opt$dependent_var
    results$tidied$exposure <- model_info$pair_info[[indep_col]]
    results$tidied$log_p <- opt$transform_dependent
    results$tidied$log_e <- opt$transform_independent
    results$tidied$scaled_p <- opt$scale_variables
    results$tidied$scaled_e <- opt$scale_variables
    results$tidied$orig_dep_mean <- model_info$orig_dep_mean
    results$tidied$orig_dep_sd <- model_info$orig_dep_sd
    results$tidied$orig_indep_mean <- model_info$orig_indep_mean
    results$tidied$orig_indep_sd <- model_info$orig_indep_sd
    
    results$glanced$independent_var <- model_info$pair_info[[indep_col]]
    results$glanced$n_obs <- model_info$n_obs
    results$glanced$series <- model_info$series
    results$glanced$phenotype <- opt$dependent_var
    results$glanced$exposure <- model_info$pair_info[[indep_col]]
    results$glanced$log_p <- opt$transform_dependent
    results$glanced$log_e <- opt$transform_independent
    results$glanced$scaled_p <- opt$scale_variables
    results$glanced$scaled_e <- opt$scale_variables
    
    results$rsq$independent_var <- model_info$pair_info[[indep_col]]
    results$rsq$n_obs <- model_info$n_obs
    results$rsq$series <- model_info$series
    results$rsq$phenotype <- opt$dependent_var
    results$rsq$exposure <- model_info$pair_info[[indep_col]]
    results$rsq$log_p <- opt$transform_dependent
    results$rsq$log_e <- opt$transform_independent
    results$rsq$scaled_p <- opt$scale_variables
    results$rsq$scaled_e <- opt$scale_variables
    
    # Collect results
    pe_tidied_list[[i]] <- results$tidied
    pe_glanced_list[[i]] <- results$glanced
    rsq_list[[i]] <- results$rsq
  }
  
  # Combine into final tables
  pe_tidied <- if(length(pe_tidied_list) > 0) do.call(rbind, pe_tidied_list) else data.frame()
  pe_glanced <- if(length(pe_glanced_list) > 0) do.call(rbind, pe_glanced_list) else data.frame()
  rsq <- if(length(rsq_list) > 0) do.call(rbind, rsq_list) else data.frame()
  
  return(list(
    pe_tidied = pe_tidied,
    pe_glanced = pe_glanced,
    rsq = rsq
  ))
}

# Extract and format results
log_info("Extracting and formatting results")
output_results <- extract_final_results(models_list)

# Add aggregate_base_model and dependent_var columns for compatibility
if(nrow(output_results$pe_tidied) > 0) {
  output_results$pe_tidied$aggregate_base_model <- FALSE
  output_results$pe_tidied$dependent_var <- opt$dependent_var
}

if(nrow(output_results$pe_glanced) > 0) {
  output_results$pe_glanced$aggregate_base_model <- FALSE
  output_results$pe_glanced$dependent_var <- opt$dependent_var
}

if(nrow(output_results$rsq) > 0) {
  output_results$rsq$aggregate_base_model <- FALSE
  output_results$rsq$dependent_var <- opt$dependent_var
}

# Save the results
output_filename <- sprintf('%s.rds', opt$dependent_var)
log_info("Saving results to {file.path(opt$output_path, output_filename)}")
saveRDS(output_results, file = file.path(opt$output_path, output_filename))

# Optionally save model objects separately if in debug mode
if(opt$debug_mode) {
  models_filename <- sprintf('%s_models.rds', opt$dependent_var)
  log_info("Saving model objects to {file.path(opt$output_path, models_filename)}")
  saveRDS(models_list, file = file.path(opt$output_path, models_filename))
}

# Summary stats
successful_models <- sum(sapply(models_list, function(x) !is.null(x$results)))
log_info("Summary: {successful_models}/{length(models_list)} models successfully fitted")

# Close the database connection
DBI::dbDisconnect(con)

log_info("Analysis complete for {opt$dependent_var}")