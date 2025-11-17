#!/usr/bin/env Rscript

#  UNIVERSAL WAS ANALYSIS - PRODUCTION VERSION (FULLY CORRECTED)
# Matches the original pipeline structure: per-dependent-variable processing
# Outputs: pe_tidied, pe_glanced, rsq (like universal_was_analysis_debugging_report_information.R)
# INCLUDES: Proper survey result extraction, log transform function, robust R-squared calculation

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(readr)
  library(survey)
  library(broom)
  library(stringr)
  library(tibble)
})

# =====================================================================================
# HELPER FUNCTIONS (RESTORED FROM ORIGINAL NHANESPEWAS)
# =====================================================================================

# Log transform function with proper zero handling
log10_xform_variable <- function(x) {
  if(any(x==0, na.rm = TRUE)) {
    x[x==0 & !is.na(x)] <- sqrt(min(x[x>0 & !is.na(x)]))
  }
  log10(x)
}

# Variable type detection function
check_e_data_type <- function(varname, con = NULL) {
  ret <- list(vartype="continuous", varlevels=NULL)
  
  if(grepl('CNT$', varname)) {
    return(list(vartype="continuous-rank", varlevels=NULL))
  }
  
  if(grepl("^PAQ", varname)) {
    return(list(vartype="continuous", varlevels=NULL))
  }
  
  if (!is.null(con) && "e_variable_levels" %in% dbListTables(con)) {
    elvl <- tbl(con, 'e_variable_levels') %>%  
              filter(`Variable.Name` == !!varname, !is.na(values)) %>% 
              collect() %>%
              pull(values) %>% 
              unique()
  } else {
    elvl <- numeric(0)
  }
  
  if(length(elvl) == 0) {
    return(ret)
  } else if(length(elvl) == 1) {
    return(list(vartype="continuous", varlevels=NULL))
  } else if(any(elvl < 1 & elvl > 0) | any(round(elvl) != elvl)) {
    return(list(vartype="continuous-rank", varlevels=sort(elvl)))
  } else if(all(round(elvl) == elvl)) {
    return(list(vartype="categorical", varlevels=sort(elvl)))
  }
  
  return(ret)
}

# Survey model fitting function with enhanced error handling
fit_survey_model <- function(model_data, formula_str, regression_type) {
  # Create survey design
  dsn <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA, 
    weights = ~WTMEC2YR,
    nest = TRUE,
    data = model_data
  )
  
  # Choose family
  if (regression_type == "logistic") {
    family_choice <- quasibinomial()
  } else {
    family_choice <- gaussian()
  }
  
  # Fit model with enhanced error handling
  fit <- tryCatch({
    if (regression_type == "logistic") {
      # For logistic models, use more robust fitting options
      survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice,
                     control = glm.control(maxit = 100, epsilon = 1e-8))
    } else {
      # For linear models, standard fitting
      survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice)
    }
  }, warning = function(w) {
    # Handle warnings but continue
    if (grepl("glm.fit: algorithm did not converge", w$message)) {
      # Try with reduced tolerance for convergence issues
      tryCatch({
        survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice,
                       control = glm.control(maxit = 200, epsilon = 1e-6))
      }, error = function(e2) {
        stop("Model convergence failed: ", e2$message)
      })
    } else {
      # For other warnings, proceed normally
      survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice)
    }
  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })
  
  return(fit)
}

# ROBUST survey results extraction function with enhanced error handling
extract_survey_results <- function(fit, regression_type) {
  
  # Check for computational singularity issues
  if (!is.null(fit$qr)) {
    qr_rank <- fit$qr$rank
    n_coef <- length(coef(fit))
    if (qr_rank < n_coef) {
      stop("Model is computationally singular: rank deficient design matrix")
    }
  }
  
  # Extract coefficient table with error handling
  tidied <- tryCatch({
    broom::tidy(fit)
  }, error = function(e) {
    # If broom fails, create minimal output
    data.frame(
      term = "(Intercept)",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    )
  })
  
  # Extract model-level statistics with error handling
  glanced <- tryCatch({
    broom::glance(fit)
  }, error = function(e) {
    # If broom fails, create minimal output
    data.frame(
      null.deviance = NA_real_,
      df.null = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      deviance = NA_real_,
      df.residual = NA_real_,
      nobs = nobs(fit)
    )
  })
  
  # Calculate R-squared manually for survey models with enhanced error handling
  if (regression_type == "linear") {
    rsq <- tryCatch({
      # Enhanced method with better error handling
      design <- fit$survey.design
      
      # Check if model converged properly
      if (is.null(fitted(fit)) || any(is.na(fitted(fit)))) {
        stop("Model fitting produced invalid fitted values")
      }
      
      # Get fitted values and residuals
      fitted_vals <- fitted(fit)
      residuals_vals <- residuals(fit)
      observed_vals <- fitted_vals + residuals_vals
      
      # Check for valid values
      if (any(is.na(observed_vals)) || any(is.infinite(observed_vals))) {
        stop("Invalid observed values detected")
      }
      
      # Calculate weighted means
      weights <- weights(design, "sampling")
      if (any(is.na(weights)) || any(weights <= 0)) {
        stop("Invalid survey weights detected")
      }
      
      y_mean <- weighted.mean(observed_vals, weights, na.rm = TRUE)
      
      # Total sum of squares (weighted)
      tss <- sum(weights * (observed_vals - y_mean)^2, na.rm = TRUE)
      
      # Residual sum of squares (weighted)
      rss <- sum(weights * residuals_vals^2, na.rm = TRUE)
      
      # Check for valid sums of squares
      if (is.na(tss) || is.na(rss) || tss <= 0) {
        stop("Invalid sums of squares calculation")
      }
      
      # R-squared calculation
      r_squared <- ifelse(tss > 0, 1 - (rss / tss), 0)
      
      # Bound R-squared between 0 and 1
      r_squared <- pmax(0, pmin(1, r_squared))
      
      # Adjusted R-squared (approximate for survey data)
      n <- nobs(fit)
      p <- length(coef(fit)) - 1  # number of predictors (excluding intercept)
      adj_r_squared <- ifelse(n > p + 1, 1 - ((1 - r_squared) * (n - 1) / (n - p - 1)), NA_real_)
      
      tibble(
        r.squared = as.numeric(r_squared),
        adj.r.squared = as.numeric(adj_r_squared)
      )
    }, error = function(e) {
      # Enhanced fallback for computational issues
      tibble(
        r.squared = NA_real_,
        adj.r.squared = NA_real_
      )
    })
  } else {
    # Pseudo R-squared for logistic survey models with error handling
    rsq <- tryCatch({
      if (is.null(fit$deviance) || is.null(fit$null.deviance) || 
          is.na(fit$deviance) || is.na(fit$null.deviance) ||
          fit$null.deviance <= 0) {
        stop("Invalid deviance values")
      }
      
      pseudo_r2 <- 1 - fit$deviance / fit$null.deviance
      pseudo_r2 <- pmax(0, pmin(1, pseudo_r2))  # Bound between 0 and 1
      
      tibble(
        r.squared = as.numeric(pseudo_r2),
        adj.r.squared = NA_real_
      )
    }, error = function(e) {
      tibble(
        r.squared = NA_real_,
        adj.r.squared = NA_real_
      )
    })
  }
  
  return(list(tidied = tidied, glanced = glanced, rsq = rsq))
}

# =====================================================================================
# COMMAND LINE INTERFACE
# =====================================================================================

# CLI options - FIXED to match old pipeline
opt_list <- list(
  make_option(c("-d", "--dependent_var"), type = "character", 
              help = "Dependent variable name (single variable)"),
  make_option(c("-s", "--schema_structure_file"), type = "character", 
              help = "Schema structure file path"),
  make_option(c("-b", "--database_path"), type = "character", 
              help = "Database file path"),
  make_option(c("-o", "--output_path"), type = "character", 
              help = "Output directory path"),
  make_option(c("-a", "--analysis_type"), type = "character", 
              help = "Analysis type (1_demoWAS, 2_oradWAS, etc.)"),
  make_option(c("-n", "--normalization"), type = "character", 
              help = "Normalization method (clr, lognorm, none)"),
  make_option(c("-t", "--test"), action = "store_true", default = FALSE,
              help = "Test mode - process only first 5 rows")
)

opt <- parse_args(OptionParser(option_list = opt_list))

# Check required arguments
required_args <- c('dependent_var', 'schema_structure_file', 'database_path', 'output_path')
missing_args <- required_args[!required_args %in% names(opt)]

if (length(missing_args) > 0) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

# Set defaults
if (is.null(opt$analysis_type)) opt$analysis_type <- "exWAS"
if (is.null(opt$normalization)) opt$normalization <- "clr"
if (is.null(opt$test)) opt$test <- FALSE

cat(" UNIVERSAL WAS ANALYSIS - PRODUCTION (FULLY CORRECTED)\n")
cat("================================================================\n")
cat("Dependent variable:", opt$dependent_var, "\n")
cat("Analysis type:", opt$analysis_type, "\n")
cat("Normalization:", opt$normalization, "\n")
cat("Test mode:", opt$test, "\n\n")

# =====================================================================================
# DATABASE CONNECTION AND SCHEMA LOADING
# =====================================================================================

# Database connection
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = opt$database_path)
cat("✅ Database connected\n")

# Load schema for this specific dependent variable
schema_data <- read_csv(opt$schema_structure_file, show_col_types = FALSE) %>%
  filter(dep_var == opt$dependent_var)

if (opt$test) {
  schema_data <- schema_data %>% slice_head(n = 5)
  cat("Test mode: Processing first", nrow(schema_data), "rows\n")
}

cat("✅ Schema loaded:", nrow(schema_data), "variable pairs for", opt$dependent_var, "\n")

if (nrow(schema_data) == 0) {
  cat("No variable pairs found for", opt$dependent_var, "\n")
  DBI::dbDisconnect(con)
  quit(status = 0)
}

# =====================================================================================
# ANALYSIS CONFIGURATION
# =====================================================================================

# Determine regression type and covariate usage
regression_type <- if (opt$analysis_type %in% c("2_oradWAS", "5_outWAS")) {
  "logistic"
} else {
  "linear"
}

use_covariates <- if (opt$analysis_type == "1_demoWAS") {
  FALSE
} else {
  TRUE
}

# Determine if transformations should be applied
apply_transformations <- opt$normalization %in% c("clr", "lognorm")

cat("Regression type:", regression_type, "\n")
cat("Use covariates:", use_covariates, "\n")
cat("Apply transformations:", apply_transformations, "\n")
cat("Normalization method:", opt$normalization, "\n\n")

# =====================================================================================
# MAIN ANALYSIS LOOP
# =====================================================================================

# Process each variable pair
models_list <- vector("list", nrow(schema_data))
successful_models <- 0
failed_models <- 0

for (i in 1:nrow(schema_data)) {
  rw <- schema_data[i, ]
  
  if (i %% 10 == 0) {
    cat("Processing", i, "/", nrow(schema_data), "...\n")
  }
  
  tryCatch({
    # 1. Variable type detection
    e_levels <- check_e_data_type(rw$indep_var, con)
    
    # 2. Check tables exist
    if (!rw$dep_table %in% dbListTables(con) || !rw$indep_table %in% dbListTables(con)) {
      failed_models <- failed_models + 1
      next
    }
    
    # 3. Data loading and merging
    if (opt$analysis_type == "1_demoWAS") {
      # For 1_demoWAS: indep_table IS the demo table
      demo_data <- dbReadTable(con, rw$indep_table)
      dep_data <- dbReadTable(con, rw$dep_table) %>% select(SEQN, !!sym(rw$dep_var))
      
      merged_data <- demo_data %>%
        inner_join(dep_data, by = "SEQN") %>%
        filter(!is.na(!!sym(rw$dep_var)), !is.na(!!sym(rw$indep_var)),
               !is.na(WTMEC2YR), WTMEC2YR > 0,
               !is.na(SDMVSTRA), !is.na(SDMVPSU))
      
    } else {
      # For other analyses: extract cycle and find demo table
      dep_data <- dbReadTable(con, rw$dep_table) %>% select(SEQN, !!sym(rw$dep_var))
      indep_data <- dbReadTable(con, rw$indep_table) %>% select(SEQN, !!sym(rw$indep_var))
      
      # Extract cycle letter with better error handling
      if (opt$analysis_type == "3_exWAS") {
        target_table <- rw$indep_table
      } else {
        target_table <- rw$dep_table
      }
      
      cycle_letter <- str_extract(target_table, "_([FG])_", group = 1)
      if (is.na(cycle_letter)) {
        cycle_letter <- str_extract(target_table, "_([FGfg])$", group = 1) %>% str_to_upper()
      }
      
      if (is.na(cycle_letter)) {
        cat("WARNING: Could not extract cycle letter from table:", target_table, "\n")
        failed_models <- failed_models + 1
        next
      }
      
      demo_table <- paste0("DEMO_", cycle_letter)
      if (!demo_table %in% dbListTables(con)) {
        cat("WARNING: Demo table not found:", demo_table, "\n")
        failed_models <- failed_models + 1
        next
      }
      
      demo_data <- dbReadTable(con, demo_table)
      
      # Merge all datasets
      merged_data <- demo_data %>%
        inner_join(dep_data, by = "SEQN") %>%
        inner_join(indep_data, by = "SEQN") %>%
        filter(!is.na(!!sym(rw$dep_var)), !is.na(!!sym(rw$indep_var)),
               !is.na(WTMEC2YR), WTMEC2YR > 0,
               !is.na(SDMVSTRA), !is.na(SDMVPSU))
    }
    
    # 4. Check minimum sample size
    if (nrow(merged_data) < 10) {
      failed_models <- failed_models + 1
      next
    }
    
    # 5. Rename variables for model
    merged_data$dep_var <- merged_data[[rw$dep_var]]
    merged_data$indep_var <- merged_data[[rw$indep_var]]
    
    # 6. Apply normalization transformations ONLY to non-microbiome dependent variables
    # Microbiome data (RSV_*_relative) is already transformed in the database tables
    is_microbiome_data <- grepl("^RSV_.*_relative$", rw$dep_var)
    is_microbiome_indep <- grepl("^RSV_.*_relative$", rw$indep_var)
    
    if (apply_transformations && !is_microbiome_data) {
      if (opt$normalization == "lognorm") {
        # Apply log transformation to non-microbiome dependent variables only
        if (any(merged_data$dep_var <= 0, na.rm = TRUE)) {
          merged_data$dep_var <- log10_xform_variable(merged_data$dep_var)
        } else {
          merged_data$dep_var <- log10(merged_data$dep_var)
        }
      }
      # Note: CLR transformation not applicable to non-microbiome data
    }
    
    # ENHANCED ZERO HANDLING for microbiome data
    # Add pseudocount to ALL microbiome dependent variables for 'none' normalization
    if (is_microbiome_data && opt$normalization == "none") {
      pseudo <- 1e-6
      merged_data$dep_var[!is.na(merged_data$dep_var)] <- merged_data$dep_var[!is.na(merged_data$dep_var)] + pseudo
    }
    
    # CRITICAL FIX: Add pseudocount to ALL microbiome independent variables for 'none' normalization
    # This is essential for 2_oradWAS, 4_pheWAS, 6_zimWAS with 'none' normalization
    if (is_microbiome_indep && opt$normalization == "none") {
      pseudo <- 1e-6
      merged_data$indep_var[!is.na(merged_data$indep_var)] <- merged_data$indep_var[!is.na(merged_data$indep_var)] + pseudo
    }
    
    # 7. Apply variable type transformations
    if (e_levels$vartype == "continuous-rank") {
      merged_data$indep_var <- as.numeric(cut(merged_data$indep_var, 
                                            breaks = quantile(merged_data$indep_var, 
                                                            probs = seq(0, 1, 0.25), na.rm = TRUE),
                                            include.lowest = TRUE))
    } else if (e_levels$vartype == "categorical" && !is.null(e_levels$varlevels)) {
      merged_data$indep_var <- factor(merged_data$indep_var, levels = e_levels$varlevels)
    }
    
    # 8. Build formula with FLEXIBLE MAXIMUM covariate selection
    if (use_covariates) {
      # Define the FULL set of covariates that should be used for all non-demoWAS analyses
      full_covariates <- c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "INDFMPIR", 
                          "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
                          "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", 
                          "ETHNICITY_NONHISPANICBLACK", "BORN_INUSA")
      
      # FLEXIBLE ASSESSMENT: Find maximum usable covariates for this specific context
      usable_covariates <- character()
      
      for (covar in full_covariates) {
        if (covar %in% names(merged_data)) {
          # Check data quality for this covariate
          covar_data <- merged_data[[covar]]
          non_missing_count <- sum(!is.na(covar_data))
          missing_rate <- 1 - (non_missing_count / nrow(merged_data))
          
          # Quality checks for covariate inclusion:
          # 1. Must exist in data
          # 2. Must have reasonable non-missing data (≤30% missing for better coverage)
          # 3. Must have variation (not all same value)
          # 4. Must not cause numerical issues
          
          if (non_missing_count > 0 && missing_rate <= 0.3) {
            # Check for variation
            unique_vals <- length(unique(covar_data[!is.na(covar_data)]))
            if (unique_vals > 1) {
              # Check for potential numerical issues
              if (is.numeric(covar_data)) {
                # For numeric variables, check for extreme values or perfect correlations
                finite_vals <- covar_data[is.finite(covar_data)]
                if (length(finite_vals) > 0) {
                  var_covar <- var(finite_vals, na.rm = TRUE)
                  if (!is.na(var_covar) && var_covar > 0) {
                    usable_covariates <- c(usable_covariates, covar)
                  }
                }
              } else {
                # For categorical variables, ensure sufficient representation in each category
                table_covar <- table(covar_data, useNA = "no")
                min_category_size <- min(table_covar)
                if (min_category_size >= 5) {  # At least 5 observations per category
                  usable_covariates <- c(usable_covariates, covar)
                }
              }
            }
          }
        }
      }
      
      cat("   Assessed", length(full_covariates), "potential covariates,", 
          length(usable_covariates), "are usable\n")
      
      # Build formula with maximum usable covariates
      if (length(usable_covariates) > 0) {
        formula_str <- paste("dep_var ~ indep_var +", paste(usable_covariates, collapse = " + "))
        cat("   Using MAXIMUM", length(usable_covariates), "covariates:", 
            paste(usable_covariates, collapse = ", "), "\n")
      } else {
        # Fallback to minimal model only if NO covariates are usable
        formula_str <- "dep_var ~ indep_var"
        cat("   No usable covariates found, using minimal model\n")
      }
      
      # Store for later use
      available_covariates <- usable_covariates
    } else {
      # For 1_demoWAS, use no covariates
      formula_str <- "dep_var ~ indep_var"
      available_covariates <- character()
    }
    
    # 9. Validate dependent variable for logistic regression
    if (regression_type == "logistic") {
      unique_vals <- unique(merged_data$dep_var)
      unique_vals <- unique_vals[!is.na(unique_vals)]
      if (!all(unique_vals %in% c(0, 1))) {
        failed_models <- failed_models + 1
        next
      }
    }
    
    # 10. Fit survey-aware model with FLEXIBLE fallback if needed
    fit <- NULL
    final_formula <- formula_str
    covariates_used <- available_covariates
    
    # First attempt: Use optimally selected covariates
    fit <- tryCatch({
      fit_survey_model(merged_data, formula_str, regression_type)
    }, error = function(e) {
      if (length(available_covariates) > 0) {
        cat("   WARNING: Optimal covariate set failed (", e$message, "), trying reduced sets...\n")
        
        # FLEXIBLE fallback: Try reducing covariates systematically
        # Priority order: keep demographic basics, then education, then ethnicity
        priority_order <- c("RIDAGEYR", "RIAGENDR", "AGE_SQUARED", "INDFMPIR",
                           "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA", "EDUCATION_COLLEGEGRAD",
                           "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", 
                           "ETHNICITY_NONHISPANICBLACK", "BORN_INUSA")
        
        # Try with progressively fewer covariates based on priority
        available_priority <- priority_order[priority_order %in% available_covariates]
        
        # Try with 75% of covariates first
        n_to_try <- max(1, floor(length(available_priority) * 0.75))
        reduced_covars <- available_priority[1:n_to_try]
        
        if (length(reduced_covars) > 0) {
          reduced_formula <- paste("dep_var ~ indep_var +", paste(reduced_covars, collapse = " + "))
          cat("   Trying", length(reduced_covars), "priority covariates:", paste(reduced_covars, collapse = ", "), "\n")
          
          tryCatch({
            fit_result <- fit_survey_model(merged_data, reduced_formula, regression_type)
            final_formula <<- reduced_formula
            covariates_used <<- reduced_covars
            return(fit_result)
          }, error = function(e2) {
            # If still failing, try with just demographic basics
            basic_covars <- c("RIDAGEYR", "RIAGENDR")[c("RIDAGEYR", "RIAGENDR") %in% available_covariates]
            if (length(basic_covars) > 0) {
              cat("   Trying basic demographics only:", paste(basic_covars, collapse = ", "), "\n")
              basic_formula <- paste("dep_var ~ indep_var +", paste(basic_covars, collapse = " + "))
              tryCatch({
                fit_result <- fit_survey_model(merged_data, basic_formula, regression_type)
                final_formula <<- basic_formula
                covariates_used <<- basic_covars
                return(fit_result)
              }, error = function(e3) {
                cat("   All covariate combinations failed, using minimal model\n")
                minimal_formula <- "dep_var ~ indep_var"
                final_formula <<- minimal_formula
                covariates_used <<- character()
                return(fit_survey_model(merged_data, minimal_formula, regression_type))
              })
            } else {
              cat("   No basic demographics available, using minimal model\n")
              minimal_formula <- "dep_var ~ indep_var"
              final_formula <<- minimal_formula
              covariates_used <<- character()
              return(fit_survey_model(merged_data, minimal_formula, regression_type))
            }
          })
        } else {
          cat("   No priority covariates available, using minimal model\n")
          minimal_formula <- "dep_var ~ indep_var"
          final_formula <<- minimal_formula
          covariates_used <<- character()
          return(fit_survey_model(merged_data, minimal_formula, regression_type))
        }
      } else {
        cat("WARNING: Model fitting failed for", rw$dep_var, "~", rw$indep_var, ":", e$message, "\n")
        return(NULL)
      }
    })
    
    if(is.null(fit)) {
      failed_models <- failed_models + 1
      next
    }
    
    # 11. Extract results using ROBUST function
    results <- tryCatch({
      extract_survey_results(fit, regression_type)
    }, error = function(e) {
      cat("WARNING: Result extraction failed for", rw$dep_var, "~", rw$indep_var, ":", e$message, "\n")
      return(NULL)
    })
    
    if(is.null(results)) {
      failed_models <- failed_models + 1
      next
    }
    
    # 12. Add metadata (like original pipeline) with formula information
    results$tidied$independent_var <- rw$indep_var
    results$tidied$n_obs <- nrow(merged_data)
    results$tidied$phenotype <- opt$dependent_var
    results$tidied$exposure <- rw$indep_var
    results$tidied$dependent_var <- opt$dependent_var
    results$tidied$formula_used <- final_formula
    results$tidied$n_covariates <- length(covariates_used)
    
    results$glanced$independent_var <- rw$indep_var
    results$glanced$n_obs <- nrow(merged_data)
    results$glanced$phenotype <- opt$dependent_var
    results$glanced$exposure <- rw$indep_var
    results$glanced$dependent_var <- opt$dependent_var
    results$glanced$formula_used <- final_formula
    results$glanced$n_covariates <- length(covariates_used)
    
    results$rsq$independent_var <- rw$indep_var
    results$rsq$n_obs <- nrow(merged_data)
    results$rsq$phenotype <- opt$dependent_var
    results$rsq$exposure <- rw$indep_var
    results$rsq$dependent_var <- opt$dependent_var
    results$rsq$formula_used <- final_formula
    results$rsq$n_covariates <- length(covariates_used)
    
    # 13. Store results
    models_list[[i]] <- results
    successful_models <- successful_models + 1
    
  }, error = function(e) {
    cat("ERROR: Analysis failed for", rw$dep_var, "~", rw$indep_var, ":", e$message, "\n")
    failed_models <- failed_models + 1
  })
}

# =====================================================================================
# RESULTS COLLATION AND OUTPUT
# =====================================================================================

# ROBUST rbind function to handle column mismatches
robust_rbind <- function(df_list) {
  if (length(df_list) == 0) {
    return(data.frame())
  }
  
  # Remove NULL elements
  df_list <- df_list[!sapply(df_list, is.null)]
  
  if (length(df_list) == 0) {
    return(data.frame())
  }
  
  # Remove empty data frames (0 rows) to avoid column addition issues
  df_list <- df_list[sapply(df_list, function(df) nrow(df) > 0)]
  
  if (length(df_list) == 0) {
    return(data.frame())
  }
  
  # If only one data frame, return it
  if (length(df_list) == 1) {
    return(df_list[[1]])
  }
  
  # Get all unique column names across all data frames
  all_cols <- unique(unlist(lapply(df_list, names)))
  
  # Standardize all data frames to have the same columns
  standardized_list <- lapply(df_list, function(df) {
    # Skip empty data frames (should already be filtered out, but double-check)
    if (nrow(df) == 0) {
      return(df)
    }
    
    missing_cols <- setdiff(all_cols, names(df))
    
    # Add missing columns with appropriate NA values
    for (col in missing_cols) {
      # Determine appropriate NA type based on existing similar columns in other dfs
      na_val <- NA_real_  # default to numeric NA
      
      # Try to infer type from other data frames
      for (other_df in df_list) {
        if (col %in% names(other_df) && nrow(other_df) > 0) {
          if (is.character(other_df[[col]])) {
            na_val <- NA_character_
          } else if (is.logical(other_df[[col]])) {
            na_val <- NA
          } else if (is.integer(other_df[[col]])) {
            na_val <- NA_integer_
          } else {
            na_val <- NA_real_
          }
          break
        }
      }
      
      # Create a vector of appropriate length with the correct NA type
      df[[col]] <- rep(na_val, nrow(df))
    }
    
    # Reorder columns to match the standard order
    df[all_cols]
  })
  
  # Filter out any remaining empty data frames
  standardized_list <- standardized_list[sapply(standardized_list, function(df) nrow(df) > 0)]
  
  if (length(standardized_list) == 0) {
    return(data.frame())
  }
  
  # Now safely rbind
  tryCatch({
    do.call(rbind, standardized_list)
  }, error = function(e) {
    cat("WARNING: rbind still failed after standardization:", e$message, "\n")
    # Return the first data frame as fallback
    if (length(standardized_list) > 0) {
      standardized_list[[1]]
    } else {
      data.frame()
    }
  })
}

# Combine results (like original pipeline)
pe_tidied_list <- list()
pe_glanced_list <- list()
rsq_list <- list()

for (i in seq_along(models_list)) {
  if (!is.null(models_list[[i]])) {
    pe_tidied_list[[i]] <- models_list[[i]]$tidied
    pe_glanced_list[[i]] <- models_list[[i]]$glanced
    rsq_list[[i]] <- models_list[[i]]$rsq
  }
}

# Create final output using ROBUST rbind function
output_results <- list(
  pe_tidied = robust_rbind(pe_tidied_list),
  pe_glanced = robust_rbind(pe_glanced_list),
  rsq = robust_rbind(rsq_list)
)

# Add compatibility columns
if(nrow(output_results$pe_tidied) > 0) {
  output_results$pe_tidied$aggregate_base_model <- FALSE
}

if(nrow(output_results$pe_glanced) > 0) {
  output_results$pe_glanced$aggregate_base_model <- FALSE
}

if(nrow(output_results$rsq) > 0) {
  output_results$rsq$aggregate_base_model <- FALSE
}

# Save results (like original pipeline)
output_filename <- paste0(opt$dependent_var, '.rds')
output_filepath <- file.path(opt$output_path, output_filename)

# Create output directory
dir.create(opt$output_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(output_results, file = output_filepath)

# =====================================================================================
# SUMMARY REPORTING
# =====================================================================================

# Summary
cat("\nANALYSIS COMPLETE\n")
cat("===================\n")
cat("Dependent variable:", opt$dependent_var, "\n")
cat("Total variable pairs:", nrow(schema_data), "\n")
cat("Successful models:", successful_models, "\n")
cat("Failed models:", failed_models, "\n")
cat("Success rate:", round(successful_models / nrow(schema_data) * 100, 1), "%\n")
cat("Results saved to:", output_filepath, "\n")

# Summary of significant results
if (nrow(output_results$pe_tidied) > 0) {
  sig_results <- output_results$pe_tidied %>% 
    filter(term != "(Intercept)", p.value < 0.05)
  cat("Significant associations (p < 0.05):", nrow(sig_results), "\n")
}

DBI::dbDisconnect(con)
cat("\n✅ Analysis complete with ROBUST survey result extraction!\n") 
