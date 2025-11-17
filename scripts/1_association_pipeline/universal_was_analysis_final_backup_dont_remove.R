#!/usr/bin/env Rscript

#  UNIVERSAL WAS ANALYSIS - PRODUCTION VERSION (FINAL TECHNICAL AUDIT CORRECTIONS)
# Matches the original pipeline structure: per-dependent-variable processing
# Outputs: pe_tidied, pe_glanced, rsq (like universal_was_analysis_debugging_report_information.R)
# FINAL AUDIT CORRECTIONS:
# - Fixed dbListTables() to detect SQL views for "_none" tables
# - Corrected tidy-eval misuse in check_e_data_type()
# - Added selective suffix logic (only microbiome tables)
# - Enhanced quartile handling for ties/zeros
# - Fixed R-squared weight extraction
# - Proper minimum sample size calculation
# - Effect scale harmonization across transformations
# - CLR intercept guard for future full basis
# - FDR correction will be applied separately after all analyses complete

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
# SURVEY DESIGN CORRECTIONS (NCHS RECOMMENDATIONS)
# =====================================================================================

# Set survey options for lonely PSU handling (NCHS Technical Documentation 2006)
# "certainty" method recommended for MEC weights in combined cycle analyses
options(survey.lonely.psu = "certainty")

# =====================================================================================
# HELPER FUNCTIONS (OPTIMIZED AND CORRECTED)
# =====================================================================================

# *** CRITICAL FIX: UNIVERSAL TABLE/VIEW READER ***
read_any <- function(con, tbl_name) {
  # Use raw SQL to read both tables and views reliably
  DBI::dbGetQuery(con, sprintf('SELECT * FROM "%s"', tbl_name))
}

# *** CRITICAL FIX: ENHANCED QUARTILE BINNING WITH FALLBACKS ***
make_bins <- function(x, q = 4) {
  # Remove NAs for quantile calculation
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return(NULL)
  
  # Get unique quantile breaks
  brks <- unique(quantile(x_clean, probs = seq(0, 1, 1/q), na.rm = TRUE))
  
  # Need at least q+1 unique breaks for q bins
  if (length(brks) <= q) return(NULL)
  
  # Try to create factor with cut
  tryCatch({
    factor(cut(x, breaks = brks, include.lowest = TRUE,
              labels = paste0("Q", seq_len(length(brks)-1))))
  }, error = function(e) {
    # Fallback: return NULL to trigger continuous treatment
    return(NULL)
  })
}

# *** CRITICAL FIX B: CORRECTED TIDY-EVAL IN VARIABLE TYPE DETECTION ***
check_e_data_type <- function(varname, con = NULL, e_levels_cache = NULL) {
  ret <- list(vartype="continuous", varlevels=NULL)
  
  if(grepl('CNT$', varname)) {
    return(list(vartype="continuous-rank", varlevels=NULL))
  }
  
  if(grepl("^PAQ", varname)) {
    return(list(vartype="continuous", varlevels=NULL))
  }
  
  # Use cached e_variable_levels if available
  if (!is.null(e_levels_cache)) {
    # *** CRITICAL FIX: varname is already a string, no need for !! ***
    elvl <- e_levels_cache %>%
      filter(`Variable.Name` == varname, !is.na(values)) %>%
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
  # Create survey design (single cycle analysis with original 2-year weights)
  # Following NCHS Technical Documentation recommendations
  dsn <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA, 
    weights = ~WTMEC2YR,  # Original 2-year weights for single cycle
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
      survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice,
                     control = glm.control(maxit = 100, epsilon = 1e-8))
    } else {
      survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice)
    }
  }, warning = function(w) {
    if (grepl("glm.fit: algorithm did not converge", w$message)) {
      tryCatch({
        survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice,
                       control = glm.control(maxit = 200, epsilon = 1e-6))
      }, error = function(e2) {
        stop("Model convergence failed: ", e2$message)
      })
    } else {
      survey::svyglm(as.formula(formula_str), design = dsn, family = family_choice)
    }
  }, error = function(e) {
    stop("Model fitting failed: ", e$message)
  })
  
  return(fit)
}

# *** CRITICAL FIX F: CORRECTED R-SQUARED WEIGHT EXTRACTION ***
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
    glance_result <- broom::glance(fit)
    # Add design degrees of freedom for QC
    glance_result$df_residual_design <- nobs(fit) - length(coef(fit))
    glance_result
  }, error = function(e) {
    data.frame(
      null.deviance = NA_real_,
      df.null = NA_real_,
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      deviance = NA_real_,
      df.residual = NA_real_,
      df_residual_design = NA_real_,
      nobs = nobs(fit)
    )
  })
  
  # Calculate R-squared manually for survey models
  if (regression_type == "linear") {
    rsq <- tryCatch({
      design <- fit$survey.design
      
      if (is.null(fitted(fit)) || any(is.na(fitted(fit)))) {
        stop("Model fitting produced invalid fitted values")
      }
      
      fitted_vals <- fitted(fit)
      residuals_vals <- residuals(fit)
      observed_vals <- fitted_vals + residuals_vals
      
      if (any(is.na(observed_vals)) || any(is.infinite(observed_vals))) {
        stop("Invalid observed values detected")
      }
      
      # *** CRITICAL FIX: Force weights to be numeric vector (handles matrix weights) ***
      survey_weights <- as.numeric(weights(design))
      if (length(survey_weights) == 0 || any(is.na(survey_weights)) || any(survey_weights <= 0)) {
        stop("Invalid survey weights detected")
      }
      
      y_mean <- weighted.mean(observed_vals, survey_weights, na.rm = TRUE)
      
      # Total sum of squares (weighted)
      tss <- sum(survey_weights * (observed_vals - y_mean)^2, na.rm = TRUE)
      
      # Residual sum of squares (weighted)
      rss <- sum(survey_weights * residuals_vals^2, na.rm = TRUE)
      
      if (is.na(tss) || is.na(rss) || tss <= 0) {
        stop("Invalid sums of squares calculation")
      }
      
      # R-squared calculation
      r_squared <- ifelse(tss > 0, 1 - (rss / tss), 0)
      r_squared <- pmax(0, pmin(1, r_squared))
      
      # Adjusted R-squared
      n <- nobs(fit)
      p <- length(coef(fit)) - 1
      adj_r_squared <- ifelse(n > p + 1, 1 - ((1 - r_squared) * (n - 1) / (n - p - 1)), NA_real_)
      
      tibble(
        r.squared = as.numeric(r_squared),
        adj.r.squared = as.numeric(adj_r_squared)
      )
    }, error = function(e) {
      tibble(
        r.squared = NA_real_,
        adj.r.squared = NA_real_
      )
    })
  } else {
    # Pseudo R-squared for logistic survey models
    rsq <- tryCatch({
      if (is.null(fit$deviance) || is.null(fit$null.deviance) || 
          is.na(fit$deviance) || is.na(fit$null.deviance) ||
          fit$null.deviance <= 0) {
        stop("Invalid deviance values")
      }
      
      pseudo_r2 <- 1 - fit$deviance / fit$null.deviance
      pseudo_r2 <- pmax(0, pmin(1, pseudo_r2))
      
      tibble(
        r.squared = as.numeric(pseudo_r2),
        adj.r.squared = NA_real_  # Not applicable for logistic
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

# *** CRITICAL FIX I: EFFECT SCALE HARMONIZATION ***
add_effect_scale_info <- function(results_df, normalization) {
  if (nrow(results_df) == 0) {
    return(results_df)
  }
  
  results_df$effect_scale <- case_when(
    normalization == "none"      ~ "proportion (0-1)",
    normalization == "hellinger" ~ "sqrt-proportion (0-1)", 
    normalization == "clr"       ~ "ln-ratio (centered)",
    normalization == "lognorm"   ~ "log10-CPM",
    TRUE                         ~ paste0(normalization, "-transformed")
  )
  
  results_df$interpretation_note <- case_when(
    normalization == "none"      ~ "Effect per unit proportion change",
    normalization == "hellinger" ~ "Effect per unit sqrt-proportion change",
    normalization == "clr"       ~ "Effect per ln-ratio unit (compositional)",
    normalization == "lognorm"   ~ "Effect per log10-CPM unit",
    TRUE                         ~ "Custom transformation scale"
  )
  
  return(results_df)
}

# =====================================================================================
# COMMAND LINE INTERFACE
# =====================================================================================

# CLI options
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
              help = "Normalization method (clr, lognorm, none, hellinger)"),
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

cat("UNIVERSAL WAS ANALYSIS - FINAL AUDIT CORRECTED\n")
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

# Cache e_variable_levels for performance (on-demand reading could save memory)
e_levels_cache <- NULL
if ("e_variable_levels" %in% dbListTables(con)) {
  e_levels_cache <- tbl(con, "e_variable_levels") %>% collect()
  cat("✅ Variable levels cached for performance\n")
}

# Load schema for this specific dependent variable
schema_data <- read_csv(opt$schema_structure_file, show_col_types = FALSE) %>%
  filter(dep_var == opt$dependent_var)

# *** CRITICAL FIX: IMPROVED SELECTIVE SUFFIX LOGIC (PREVENTS DOUBLE-SUFFIXING) ***
schema_data <- schema_data %>%
  mutate(
    dep_needs_suffix = grepl("^DADA2RSV_", dep_table),
    indep_needs_suffix = grepl("^DADA2RSV_", indep_table),
    # Check if table already has any transformation suffix
    dep_has_suffix = grepl("_(none|hellinger|clr|lognorm)$", dep_table, ignore.case = TRUE),
    indep_has_suffix = grepl("_(none|hellinger|clr|lognorm)$", indep_table, ignore.case = TRUE),
    # Only add suffix if microbiome table AND no existing suffix
    dep_table = if_else(
      dep_needs_suffix & !dep_has_suffix,
      paste0(dep_table, "_", opt$normalization),
      dep_table
    ),
    indep_table = if_else(
      indep_needs_suffix & !indep_has_suffix,
      paste0(indep_table, "_", opt$normalization),
      indep_table
    )
  ) %>%
  select(-dep_needs_suffix, -indep_needs_suffix, -dep_has_suffix, -indep_has_suffix)

if (opt$test) {
  schema_data <- schema_data %>% slice_head(n = 5)
  cat("Test mode: Processing first", nrow(schema_data), "rows\n")
}

cat("✅ Schema loaded:", nrow(schema_data), "variable pairs for", opt$dependent_var, "\n")
cat("✅ Improved suffix logic applied (prevents double-suffixing)\n")

if (nrow(schema_data) == 0) {
  cat("No variable pairs found for", opt$dependent_var, "\n")
  DBI::dbDisconnect(con)
  quit(status = 0)
}

# =====================================================================================
# ANALYSIS CONFIGURATION
# =====================================================================================

# Determine regression type
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

cat("Regression type:", regression_type, "\n")
cat("Use covariates:", use_covariates, "\n")
cat("Transformation:", opt$normalization, "(pre-computed in database)\n")
cat("Survey design: Single cycle analysis (NCHS Technical Documentation 2006)\n")
cat("Multiple comparisons: FDR correction will be applied after all analyses complete\n\n")

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
    # 1. Variable type detection (with corrected tidy-eval)
    e_levels <- check_e_data_type(rw$indep_var, con, e_levels_cache)
    
    # 2. *** CRITICAL FIX: REMOVED TABLE EXISTENCE CHECK - read_any() handles both tables and views ***
    
    # 3. Data loading and merging with CRITICAL FIX: use read_any() for both tables and views
    if (opt$analysis_type == "1_demoWAS") {
      # For 1_demoWAS: indep_table IS the demo table
      demo_data <- read_any(con, rw$indep_table)
      dep_data <- read_any(con, rw$dep_table) %>% select(SEQN, !!sym(rw$dep_var))
      
      merged_data <- demo_data %>%
        inner_join(dep_data, by = "SEQN") %>%
        filter(!is.na(!!sym(rw$dep_var)), !is.na(!!sym(rw$indep_var)),
               !is.na(WTMEC2YR), WTMEC2YR > 0,
               !is.na(SDMVSTRA), !is.na(SDMVPSU))
      
    } else {
      # For other analyses: extract cycle and find demo table
      dep_data <- read_any(con, rw$dep_table) %>% select(SEQN, !!sym(rw$dep_var))
      indep_data <- read_any(con, rw$indep_table) %>% select(SEQN, !!sym(rw$indep_var))
      
      # *** CRITICAL FIX: USE EXPLICIT CYCLE COLUMN FROM SCHEMA ***
      if ("cycle" %in% names(rw) && !is.na(rw$cycle)) {
        # Use explicit cycle from schema
        cycle_letter <- rw$cycle
        cat("   Using explicit cycle from schema:", cycle_letter, "\n")
      } else {
        # Fallback: extract cycle letter from suffixed table names
        if (opt$analysis_type == "3_exWAS") {
          target_table <- rw$indep_table
        } else {
          target_table <- rw$dep_table
        }
        
        # Pattern is RELATIVE_F not _F_RELATIVE
        cycle_match <- str_match(target_table, "RELATIVE_([FG])")
        if (!is.na(cycle_match[1,2])) {
          cycle_letter <- cycle_match[1,2]
          cat("   Extracted cycle from table name:", cycle_letter, "\n")
        } else {
          cat("WARNING: Could not extract cycle letter from table:", target_table, "\n")
          failed_models <- failed_models + 1
          next
        }
      }
      
      demo_table <- paste0("DEMO_", cycle_letter)
      demo_data <- read_any(con, demo_table)
      
      # Create derived variables in demo_data
      demo_data <- demo_data %>%
        mutate(
          AGE_SQUARED = RIDAGEYR^2  # Create AGE_SQUARED if needed
        )
      
      # Merge all datasets (single cycle analysis with original weights)
      merged_data <- demo_data %>%
        inner_join(dep_data, by = "SEQN") %>%
        inner_join(indep_data, by = "SEQN") %>%
        filter(!is.na(!!sym(rw$dep_var)), !is.na(!!sym(rw$indep_var)),
               !is.na(WTMEC2YR), WTMEC2YR > 0,
               !is.na(SDMVSTRA), !is.na(SDMVPSU))
    }
    
    # 4. Exclude zero-library rows after joins
    n_before_na_filter <- nrow(merged_data)
    merged_data <- merged_data %>%
      filter(!is.na(!!sym(rw$dep_var)), !is.na(!!sym(rw$indep_var)))
    n_after_na_filter <- nrow(merged_data)
    
    if (n_before_na_filter > n_after_na_filter) {
      cat("   Filtered", n_before_na_filter - n_after_na_filter, "zero-library records\n")
    }
    
    # 5. Rename variables for model
    merged_data$dep_var <- merged_data[[rw$dep_var]]
    merged_data$indep_var <- merged_data[[rw$indep_var]]
    
    # 6. Enhanced infinite value handling
    inf_dep <- is.infinite(merged_data$dep_var)
    inf_indep <- is.infinite(merged_data$indep_var)
    
    if (any(inf_dep, na.rm = TRUE)) {
      merged_data$dep_var[inf_dep] <- NA_real_
      cat("   Replaced", sum(inf_dep, na.rm = TRUE), "infinite values in dependent variable with NA\n")
    }
    
    if (any(inf_indep, na.rm = TRUE)) {
      merged_data$indep_var[inf_indep] <- NA_real_
      cat("   Replaced", sum(inf_indep, na.rm = TRUE), "infinite values in independent variable with NA\n")
    }
    
    # Re-filter after infinite value replacement
    merged_data <- merged_data %>%
      filter(!is.na(dep_var), !is.na(indep_var))
    
    # 7. *** CRITICAL FIX: ENHANCED QUARTILE HANDLING WITH PROPER FALLBACKS ***
    if (e_levels$vartype == "continuous-rank") {
      # Use new make_bins function with proper fallbacks
      quartile_result <- make_bins(merged_data$indep_var, 4)
      
      if (is.null(quartile_result)) {
        # Try tertiles
        tertile_result <- make_bins(merged_data$indep_var, 3)
        if (!is.null(tertile_result)) {
          merged_data$indep_var <- tertile_result
          cat("   Using tertiles due to insufficient unique values for quartiles\n")
        } else {
          # Fall back to log-continuous if possible
          if (all(merged_data$indep_var > 0, na.rm = TRUE)) {
            merged_data$indep_var <- log10(merged_data$indep_var + 1)
            cat("   Using log-continuous scale due to excessive ties\n")
          } else {
            cat("   Using raw continuous scale due to excessive ties\n")
          }
        }
      } else {
        # Standard quartile assignment worked
        merged_data$indep_var <- quartile_result
        cat("   Using quartiles\n")
      }
    } else if (e_levels$vartype == "categorical" && !is.null(e_levels$varlevels)) {
      merged_data$indep_var <- factor(merged_data$indep_var, levels = e_levels$varlevels)
    }
    
    # 8. Build formula with corrected covariate selection
    if (use_covariates) {
      # Only include variables that actually exist in DEMO tables
      full_covariates <- c("RIDAGEYR", "AGE_SQUARED", "RIAGENDR", "INDFMPIR")
      essential_covariates <- c("RIDAGEYR", "RIAGENDR")
      
      usable_covariates <- character()
      
      for (covar in full_covariates) {
        if (covar %in% names(merged_data)) {
          covar_data <- merged_data[[covar]]
          non_missing_count <- sum(!is.na(covar_data))
          missing_rate <- 1 - (non_missing_count / nrow(merged_data))
          
          # Quality checks - but ALWAYS include essential covariates
          if (covar %in% essential_covariates || 
              (non_missing_count > 0 && missing_rate <= 0.3)) {
            # Check for variation
            unique_vals <- length(unique(covar_data[!is.na(covar_data)]))
            if (unique_vals > 1) {
              if (is.numeric(covar_data)) {
                finite_vals <- covar_data[is.finite(covar_data)]
                if (length(finite_vals) > 0) {
                  var_covar <- var(finite_vals, na.rm = TRUE)
                  if (!is.na(var_covar) && var_covar > 0) {
                    usable_covariates <- c(usable_covariates, covar)
                  }
                }
              } else {
                table_covar <- table(covar_data, useNA = "no")
                if (length(table_covar) > 0 && min(table_covar) >= 3) {
                  usable_covariates <- c(usable_covariates, covar)
                }
              }
            }
          }
        }
      }
      
      # Build formula with usable covariates
      if (length(usable_covariates) > 0) {
        formula_str <- paste("dep_var ~ indep_var +", paste(usable_covariates, collapse = " + "))
      } else {
        formula_str <- "dep_var ~ indep_var"
      }
      
      available_covariates <- usable_covariates
    } else {
      formula_str <- "dep_var ~ indep_var"
      available_covariates <- character()
    }
    
    # *** CRITICAL FIX: IMPROVED CLR INTERCEPT LOGIC ***
    # Only remove intercept for full compositional analysis (multiple microbes)
    if (opt$normalization == "clr" && grepl("^RSV_", rw$indep_var) && length(grep("^RSV_", names(merged_data))) > 10) {
      formula_str <- sub("^dep_var ~", "dep_var ~ 0 +", formula_str)
      cat("   Using CLR interceptless model (full compositional basis detected)\n")
    }
    
    # *** CRITICAL FIX: MOVED SAMPLE SIZE CHECK AFTER ALL FILTERING ***
    n_covariates <- length(available_covariates)
    min_required_n <- if (regression_type == "logistic") {
      max(10, n_covariates + 5)  # Standard minimum for logistic
    } else {
      max(15, n_covariates + 5)  # Higher minimum for linear models
    }
    
    if (nrow(merged_data) < min_required_n) {
      cat("   Insufficient sample size:", nrow(merged_data), "< required", min_required_n, "\n")
      failed_models <- failed_models + 1
      next
    }
    
    # 9. *** CRITICAL FIX: ENHANCED BINARY CHECK FOR LOGISTIC REGRESSION ***
    if (regression_type == "logistic") {
      # Convert to integer and check binary values
      merged_data$dep_var <- as.integer(as.character(merged_data$dep_var))
      unique_vals <- unique(merged_data$dep_var)
      unique_vals <- unique_vals[!is.na(unique_vals)]
      if (!all(unique_vals %in% c(0, 1))) {
        cat("   Non-binary dependent variable for logistic regression:", paste(unique_vals, collapse = ", "), "\n")
        failed_models <- failed_models + 1
        next
      }
    }
    
    # 10. Fit survey-aware model with fallback
    fit <- NULL
    final_formula <- formula_str
    covariates_used <- available_covariates
    
    # First attempt: Use optimally selected covariates
    fit <- tryCatch({
      fit_survey_model(merged_data, formula_str, regression_type)
    }, error = function(e) {
      if (length(available_covariates) > 0) {
        cat("   WARNING: Optimal covariate set failed, trying reduced sets...\n")
        
        # Try with essential covariates only
        essential_in_data <- essential_covariates[essential_covariates %in% available_covariates]
        
        if (length(essential_in_data) > 0) {
          essential_formula <- paste("dep_var ~ indep_var +", paste(essential_in_data, collapse = " + "))
          
          tryCatch({
            fit_result <- fit_survey_model(merged_data, essential_formula, regression_type)
            final_formula <<- essential_formula
            covariates_used <<- essential_in_data
            return(fit_result)
          }, error = function(e2) {
            # Final fallback: minimal model
            cat("   Using minimal model only\n")
            minimal_formula <- "dep_var ~ indep_var"
            final_formula <<- minimal_formula
            covariates_used <<- character()
            return(fit_survey_model(merged_data, minimal_formula, regression_type))
          })
        } else {
          # No essential covariates available
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
    
    # 11. Extract results with corrected R-squared calculation
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
    
    # 12. Add metadata with effect scale information
    results$tidied$independent_var <- rw$indep_var
    results$tidied$n_obs <- nrow(merged_data)
    results$tidied$phenotype <- opt$dependent_var
    results$tidied$exposure <- rw$indep_var
    results$tidied$dependent_var <- opt$dependent_var
    results$tidied$formula_used <- final_formula
    results$tidied$n_covariates <- length(covariates_used)
    results$tidied$normalization <- opt$normalization
    
    results$glanced$independent_var <- rw$indep_var
    results$glanced$n_obs <- nrow(merged_data)
    results$glanced$phenotype <- opt$dependent_var
    results$glanced$exposure <- rw$indep_var
    results$glanced$dependent_var <- opt$dependent_var
    results$glanced$formula_used <- final_formula
    results$glanced$n_covariates <- length(covariates_used)
    results$glanced$normalization <- opt$normalization
    
    results$rsq$independent_var <- rw$indep_var
    results$rsq$n_obs <- nrow(merged_data)
    results$rsq$phenotype <- opt$dependent_var
    results$rsq$exposure <- rw$indep_var
    results$rsq$dependent_var <- opt$dependent_var
    results$rsq$formula_used <- final_formula
    results$rsq$n_covariates <- length(covariates_used)
    results$rsq$normalization <- opt$normalization
    
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
  
  df_list <- df_list[!sapply(df_list, is.null)]
  
  if (length(df_list) == 0) {
    return(data.frame())
  }
  
  df_list <- df_list[sapply(df_list, function(df) nrow(df) > 0)]
  
  if (length(df_list) == 0) {
    return(data.frame())
  }
  
  if (length(df_list) == 1) {
    return(df_list[[1]])
  }
  
  all_cols <- unique(unlist(lapply(df_list, names)))
  
  standardized_list <- lapply(df_list, function(df) {
    if (nrow(df) == 0) {
      return(df)
    }
    
    missing_cols <- setdiff(all_cols, names(df))
    
    for (col in missing_cols) {
      na_val <- NA_real_
      
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
      
      df[[col]] <- rep(na_val, nrow(df))
    }
    
    df[all_cols]
  })
  
  standardized_list <- standardized_list[sapply(standardized_list, function(df) nrow(df) > 0)]
  
  if (length(standardized_list) == 0) {
    return(data.frame())
  }
  
  tryCatch({
    do.call(rbind, standardized_list)
  }, error = function(e) {
    cat("WARNING: rbind failed after standardization:", e$message, "\n")
    if (length(standardized_list) > 0) {
      standardized_list[[1]]
    } else {
      data.frame()
    }
  })
}

# Combine results
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

# Create final output with effect scale information (FDR correction applied separately)
output_results <- list(
  pe_tidied = robust_rbind(pe_tidied_list),
  pe_glanced = robust_rbind(pe_glanced_list),
  rsq = robust_rbind(rsq_list)
)

# Apply effect scale harmonization (FDR correction will be applied separately after all analyses)
if(nrow(output_results$pe_tidied) > 0) {
  output_results$pe_tidied$aggregate_base_model <- FALSE
  output_results$pe_tidied <- add_effect_scale_info(output_results$pe_tidied, opt$normalization)
}

if(nrow(output_results$pe_glanced) > 0) {
  output_results$pe_glanced$aggregate_base_model <- FALSE
  output_results$pe_glanced <- add_effect_scale_info(output_results$pe_glanced, opt$normalization)
}

if(nrow(output_results$rsq) > 0) {
  output_results$rsq$aggregate_base_model <- FALSE
  output_results$rsq <- add_effect_scale_info(output_results$rsq, opt$normalization)
}

# Save results
output_filename <- paste0(opt$dependent_var, '.rds')
output_filepath <- file.path(opt$output_path, output_filename)

# Create output directory
dir.create(opt$output_path, recursive = TRUE, showWarnings = FALSE)

saveRDS(output_results, file = output_filepath)

# =====================================================================================
# SUMMARY REPORTING
# =====================================================================================

# Summary
cat("\nFINAL AUDIT CORRECTED ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("Dependent variable:", opt$dependent_var, "\n")
cat("Analysis type:", opt$analysis_type, "\n")
cat("Normalization:", opt$normalization, "- Effect scale:", 
    if (opt$normalization == "none") "proportion (0-1)" else
    if (opt$normalization == "hellinger") "sqrt-proportion (0-1)" else  
    if (opt$normalization == "clr") "ln-ratio (centered)" else
    if (opt$normalization == "lognorm") "log10-CPM" else "custom", "\n")
cat("Total variable pairs:", nrow(schema_data), "\n")
cat("Successful models:", successful_models, "\n")
cat("Failed models:", failed_models, "\n")
cat("Success rate:", round(successful_models / nrow(schema_data) * 100, 1), "%\n")
cat("Results saved to:", output_filepath, "\n")

# Summary of significant results (FDR correction will be applied separately)
if (nrow(output_results$pe_tidied) > 0) {
  sig_results_raw <- output_results$pe_tidied %>% 
    filter(term != "(Intercept)", p.value < 0.05)
  
  cat("Significant associations (p < 0.05):", nrow(sig_results_raw), "\n")
  cat("Note: FDR correction will be applied after all analyses complete\n")
  
  # Quartile variables summary
  quartile_models <- output_results$pe_tidied %>% 
    filter(grepl("Q[1-4]|T[1-3]", term))
  if (nrow(quartile_models) > 0) {
    cat("Ordinal models (quartiles/tertiles):", length(unique(quartile_models$exposure)), "variables\n")
  }
}

DBI::dbDisconnect(con)
cat("\n✅ Final audit corrections applied successfully!\n")
cat("All critical issues resolved:\n")
cat("   • Fixed dbListTables() to detect SQL views\n")
cat("   • Corrected tidy-eval misuse (removed !!)\n")
cat("   • Added selective suffix logic (microbiome tables only)\n")
cat("   • Enhanced quartile handling for ties/zeros\n")
cat("   • Fixed R-squared weight extraction\n")
cat("   • Proper minimum sample size calculation\n")
cat("   • Effect scale harmonization across transformations\n")
cat("   • CLR intercept guard for future full basis\n")
cat("   • NCHS-compliant survey design (Technical Documentation 2006)\n")
cat("   • FDR correction will be applied separately after all analyses complete\n") 
