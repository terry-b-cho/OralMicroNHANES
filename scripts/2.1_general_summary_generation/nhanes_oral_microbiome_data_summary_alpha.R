# ═══════════════════════════════════════════════════════════════════════════════
# NHANES ORAL MICROBIOME DATA SUMMARY WITH ALPHA DIVERSITY STRATIFICATION
# ═══════════════════════════════════════════════════════════════════════════════
# Author: Assistant
# Purpose: Generate comprehensive data summary tables for NHANES oral microbiome study
# Database: nhanes_oral_transformed_complete_processed.sqlite
# Output: Publication-ready supplementary tables with alpha diversity stratification
# ═══════════════════════════════════════════════════════════════════════════════

# ── LIBRARY LOADING ────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# ── CONFIGURATION ──────────────────────────────────────────────────────────────
DB_PATH <- "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"
ALPHA_DIV_PATH <- "data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt"
ALPHA_VAR_PATH <- "data/00_nhanes_omp_diversity_db/dada2rsv-alpha-variablelist.txt"
OUTPUT_DIR <- "results/supplementary_tables"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── HELPER FUNCTIONS ───────────────────────────────────────────────────────────

# Logging function
log_message <- function(...) {
  cat(paste0("[", Sys.time(), "] ", ..., "\n"))
}

# Calculate Wilson confidence interval for proportions
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(lower = 0, upper = 0))
  
  z <- qnorm((1 + conf.level) / 2)
  p <- x / n
  
  # Wilson score interval
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  
  return(c(lower = max(0, center - width), upper = min(1, center + width)))
}

# Variable type classification
classify_variable_type <- function(x) {
  if (all(is.na(x))) return("Missing")
  
  unique_vals <- unique(x[!is.na(x)])
  n_unique <- length(unique_vals)
  
  # Binary check
  if (n_unique == 2 && all(unique_vals %in% c(0, 1))) {
    return("Binary")
  }
  
  # Categorical check (non-numeric or few unique values)
  if (is.character(x) || is.factor(x) || (is.numeric(x) && n_unique <= 10)) {
    return("Categorical")
  }
  
  # Ordinal check (ordered factor or specific patterns)
  if (is.ordered(x)) {
    return("Ordinal")
  }
  
  # Continuous (many unique numeric values)
  if (is.numeric(x) && n_unique > 10) {
    return("Continuous")
  }
  
  return("Other")
}

# ── MAIN FUNCTIONS ─────────────────────────────────────────────────────────────

# Load WAS configuration files
load_was_variables <- function() {
  log_message("Loading WAS configuration files...")
  
  config_files <- c(
    "1_demoWAS" = "configs/1_demoWAS_vars.txt",
    "2_oradWAS" = "configs/2_oradWAS_vars.txt", 
    "3_exWAS" = "configs/3_exWAS_vars.txt",
    "4_pheWAS" = "configs/4_pheWAS_vars.txt",
    "5_outWAS" = "configs/5_outWAS_vars.txt",
    "6_zimWAS" = "configs/6_zimWAS_vars.txt"
  )
  
  was_vars <- map(config_files, ~ {
    if (file.exists(.x)) {
      vars <- read_lines(.x) %>% str_trim() %>% .[. != ""]
      log_message("  ", basename(.x), ": ", length(vars), " variables")
      return(vars)
    } else {
      log_message("  ⚠️ File not found: ", .x)
      return(character(0))
    }
  })
  
  return(was_vars)
}

# Get microbiome participants
get_microbiome_participants <- function(con) {
  log_message("🦠 Identifying participants with microbiome data...")
  
  microbiome_participants <- dbGetQuery(con, "
    SELECT DISTINCT SEQN, 'F' as cycle
    FROM DEMO_F 
    WHERE SEQN IS NOT NULL
    UNION
    SELECT DISTINCT SEQN, 'G' as cycle
    FROM DEMO_G
    WHERE SEQN IS NOT NULL
  ")
  
  microbiome_participants$SEQN <- as.character(microbiome_participants$SEQN)
  
  log_message("  Total participants with microbiome data: ", nrow(microbiome_participants))
  log_message("    Cycle F: ", sum(microbiome_participants$cycle == "F"))
  log_message("    Cycle G: ", sum(microbiome_participants$cycle == "G"))
  
  return(microbiome_participants)
}

# Load and analyze alpha diversity data
load_alpha_diversity_data <- function() {
  log_message("Loading alpha diversity data...")
  
  alpha_div_file <- "data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt"
  if (!file.exists(alpha_div_file)) {
    log_message("❌ Alpha diversity file not found: ", alpha_div_file)
    stop("Alpha diversity file not found")
  }
  
  alpha_div_raw <- read.table(alpha_div_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  log_message(sprintf("    Loaded %d rows, %d columns", nrow(alpha_div_raw), ncol(alpha_div_raw)))
  
  # Convert SEQN to numeric if it's character
  if (is.character(alpha_div_raw$SEQN)) {
    alpha_div_raw$SEQN <- as.numeric(alpha_div_raw$SEQN)
  }
  
  # Look at actual column names to understand the structure
  col_names <- names(alpha_div_raw)
  log_message("    Sample column names: ", paste(head(col_names, 10), collapse = ", "))
  
  # Find columns that contain diversity metrics at 10000 read depth using the actual format
  observed_cols <- grep("RSV_ObservedOTUs_10000_9", col_names, value = TRUE)
  faith_cols <- grep("RSV_FaPhyloDiv_10000_9", col_names, value = TRUE)
  shannon_cols <- grep("RSV_ShanWienDiv_10000_9", col_names, value = TRUE)
  simpson_cols <- grep("RSV_InverseSimpson_10000_9", col_names, value = TRUE)
  
  log_message(sprintf("    Found %d observed, %d faith, %d shannon, %d simpson columns at 10000 depth", 
                     length(observed_cols), length(faith_cols), length(shannon_cols), length(simpson_cols)))
  
  # Convert alpha diversity columns to numeric
  for (col in c(observed_cols, faith_cols, shannon_cols, simpson_cols)) {
    if (col %in% names(alpha_div_raw)) {
      alpha_div_raw[[col]] <- as.numeric(alpha_div_raw[[col]])
    }
  }
  
  # Create summary metrics by averaging across resampling iterations
  log_message("📈 Creating summary alpha diversity metrics...")
  alpha_summary <- alpha_div_raw %>%
    filter(!is.na(SEQN)) %>%
    group_by(SEQN) %>%
    summarise(
      # Use the specific columns for 10000 read depth, 9th resampling
      Observed_ASVs_10000 = if(length(observed_cols) > 0) {
        rowMeans(select(pick(everything()), all_of(observed_cols)), na.rm = TRUE)
      } else { NA_real_ },
      
      Faith_PD_10000 = if(length(faith_cols) > 0) {
        rowMeans(select(pick(everything()), all_of(faith_cols)), na.rm = TRUE)
      } else { NA_real_ },
      
      Shannon_10000 = if(length(shannon_cols) > 0) {
        rowMeans(select(pick(everything()), all_of(shannon_cols)), na.rm = TRUE)
      } else { NA_real_ },
      
      InvSimpson_10000 = if(length(simpson_cols) > 0) {
        rowMeans(select(pick(everything()), all_of(simpson_cols)), na.rm = TRUE)
      } else { NA_real_ },
      
      .groups = 'drop'
    )
  
  log_message(sprintf("    Created summary metrics for %d participants", nrow(alpha_summary)))
  return(alpha_summary)
}

# Analyze alpha diversity stratification
analyze_alpha_diversity_stratification <- function(demographics_data, alpha_summary) {
  log_message("Analyzing alpha diversity stratification...")
  
  # Merge alpha diversity with demographics
  log_message("🔗 Merging alpha diversity with demographics...")
  alpha_with_demo <- merge(alpha_summary, demographics_data, by = "SEQN", all.x = TRUE)
  
  log_message(sprintf("    Merged data for %d participants", nrow(alpha_with_demo)))
  
  # Define basic stratification variables that we can create from raw NHANES data
  strat_variables <- c(
    "Gender" = "Gender",
    "Age Quartile" = "AgeQuartile", 
    "Education Level" = "EducationLevel",
    "Ethnicity" = "Ethnicity",
    "PIR Quartile" = "PIRQuartile"
  )
  
  # Alpha diversity metrics
  alpha_metrics <- c("Observed_ASVs_10000", "Faith_PD_10000", "Shannon_10000", "InvSimpson_10000")
  metric_names <- c("Observed ASVs", "Faith's PD", "Shannon Diversity", "Inverse Simpson")
  
  results <- list()
  
  for (i in seq_along(strat_variables)) {
    strat_var <- strat_variables[i]
    strat_name <- names(strat_variables)[i]
    
    log_message("    Analyzing ", strat_name, "...")
    
    # Check if the stratification variable exists in the data
    if (!strat_var %in% names(alpha_with_demo)) {
      log_message("      ⚠️ Variable not found: ", strat_var)
      next
    }
    
    for (j in seq_along(alpha_metrics)) {
      metric <- alpha_metrics[j]
      metric_name <- metric_names[j]
      
      # Get data for this stratification variable and metric
      strat_data <- alpha_with_demo %>%
        filter(!is.na(.data[[metric]]) & !is.na(.data[[strat_var]])) %>%
        group_by(.data[[strat_var]]) %>%
        summarise(
          N = n(),
          Mean = mean(.data[[metric]], na.rm = TRUE),
          SE = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
          Median = median(.data[[metric]], na.rm = TRUE),
          Q1 = quantile(.data[[metric]], 0.25, na.rm = TRUE),
          Q3 = quantile(.data[[metric]], 0.75, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        mutate(
          Stratification_Variable = strat_name,
          Alpha_Diversity_Metric = metric_name,
          Mean_SE = sprintf("%.3f ± %.3f", Mean, SE),
          Median_IQR = sprintf("%.3f (%.3f-%.3f)", Median, Q1, Q3),
          Category = .data[[strat_var]]
        ) %>%
        select(Stratification_Variable, Category, Alpha_Diversity_Metric, N, Mean_SE, Median_IQR)
      
      results[[paste(strat_name, metric_name, sep = "_")]] <- strat_data
    }
  }
  
  # Combine all results
  combined_results <- bind_rows(results)
  
  # Debug: Check column names
  log_message("    Debug: Column names in combined results: ", paste(names(combined_results), collapse = ", "))
  
  log_message("  Alpha diversity stratification analysis complete")
  log_message(sprintf("    Generated %d stratification results", nrow(combined_results)))
  return(combined_results)
}

# Get variable metadata
get_variable_metadata <- function(con, variables) {
  log_message("📝 Extracting variable metadata...")
  
  var_metadata <- dbGetQuery(con, "
    SELECT `Data.File.Name` as nhanes_code, 
           `Data.File.Description` as description
    FROM variable_names_epcf
  ")
  
  metadata_df <- tibble(
    variable = variables,
    nhanes_code = variables,
    description = map_chr(variables, ~ {
      desc <- var_metadata$description[var_metadata$nhanes_code == .x]
      if (length(desc) > 0 && !is.na(desc[1]) && desc[1] != "") {
        return(desc[1])
      } else {
        return(paste("Variable:", .x))
      }
    })
  )
  
  log_message("  Metadata extracted for ", nrow(metadata_df), " variables")
  return(metadata_df)
}

# Collect variable data from database
collect_variable_data <- function(con, variables, microbiome_participants) {
  log_message("Collecting variable data from database...")
  
  tables <- dbListTables(con)
  log_message("  Total tables in database: ", length(tables))
  
  all_data <- microbiome_participants %>% select(SEQN, cycle)
  variables_found <- character(0)
  
  for (var in variables) {
    found_in_table <- FALSE
    
    for (table in tables) {
      if (grepl("^(sqlite_|variable_names)", table)) next
      
      tryCatch({
        columns <- dbListFields(con, table)
        if (var %in% columns && "SEQN" %in% columns) {
          
          query <- sprintf("SELECT SEQN, `%s` as `%s` FROM `%s` WHERE `%s` IS NOT NULL", 
                          var, var, table, var)
          var_data <- dbGetQuery(con, query)
          var_data$SEQN <- as.character(var_data$SEQN)
          
          if (nrow(var_data) > 0) {
            all_data <- all_data %>%
              left_join(var_data, by = "SEQN", suffix = c("", paste0(".", table)))
            
            variables_found <- c(variables_found, var)
            found_in_table <- TRUE
            log_message("    ✓ ", var, " found in ", table, " (", nrow(var_data), " records)")
            break
          }
        }
      }, error = function(e) {
        # Skip problematic tables silently
      })
    }
    
    if (!found_in_table) {
      log_message("    ✗ ", var, " not found in any table")
    }
  }
  
  log_message("  Variables successfully collected: ", length(variables_found), "/", length(variables))
  return(list(data = all_data, variables_found = variables_found))
}

# Load demographics data from processed database
load_demographics_data <- function(con) {
  log_message("Loading demographics data...")
  
  # Try to get demographics from the processed database
  demographics_data <- tryCatch({
    # Get from individual DEMO tables with only common columns
    demo_f_query <- "SELECT SEQN, RIDAGEYR, RIAGENDR, DMDEDUC2, RIDRETH1, INDFMPIR FROM DEMO_F"
    demo_g_query <- "SELECT SEQN, RIDAGEYR, RIAGENDR, DMDEDUC2, RIDRETH1, INDFMPIR FROM DEMO_G"
    
    demo_f <- dbGetQuery(con, demo_f_query)
    demo_g <- dbGetQuery(con, demo_g_query)
    
    # Combine demo data
    demo_combined <- rbind(demo_f, demo_g)
    
    # Create basic stratification variables from raw NHANES codes
    demo_processed <- demo_combined %>%
      mutate(
        # Age quartiles
        AgeQuartile = case_when(
          is.na(RIDAGEYR) ~ "Missing",
          RIDAGEYR <= quantile(RIDAGEYR, 0.25, na.rm = TRUE) ~ "Q1 (14-27)",
          RIDAGEYR <= quantile(RIDAGEYR, 0.50, na.rm = TRUE) ~ "Q2 (28-40)", 
          RIDAGEYR <= quantile(RIDAGEYR, 0.75, na.rm = TRUE) ~ "Q3 (41-55)",
          TRUE ~ "Q4 (56-69)"
        ),
        
        # Gender
        Gender = case_when(
          RIAGENDR == 1 ~ "Male",
          RIAGENDR == 2 ~ "Female",
          TRUE ~ "Missing"
        ),
        
        # Education level
        EducationLevel = case_when(
          DMDEDUC2 == 1 ~ "< 9th Grade",
          DMDEDUC2 == 2 ~ "9-11th Grade",
          DMDEDUC2 == 3 ~ "High School",
          DMDEDUC2 == 4 ~ "College/AA",
          DMDEDUC2 == 5 ~ "College Graduate",
          TRUE ~ "Missing"
        ),
        
        # Ethnicity
        Ethnicity = case_when(
          RIDRETH1 == 1 ~ "Mexican American",
          RIDRETH1 == 2 ~ "Other Hispanic",
          RIDRETH1 == 3 ~ "Non-Hispanic White",
          RIDRETH1 == 4 ~ "Non-Hispanic Black",
          RIDRETH1 == 5 ~ "Other/Multi-racial",
          TRUE ~ "Missing"
        ),
        
        # PIR quartiles
        PIRQuartile = case_when(
          is.na(INDFMPIR) ~ "Missing",
          INDFMPIR <= quantile(INDFMPIR, 0.25, na.rm = TRUE) ~ "Q1 (Lowest)",
          INDFMPIR <= quantile(INDFMPIR, 0.50, na.rm = TRUE) ~ "Q2",
          INDFMPIR <= quantile(INDFMPIR, 0.75, na.rm = TRUE) ~ "Q3",
          TRUE ~ "Q4 (Highest)"
        )
      )
    
    log_message(sprintf("    Loaded demographics for %d participants", nrow(demo_processed)))
    return(demo_processed)
    
  }, error = function(e) {
    log_message("❌ Error loading demographics: ", e$message)
    return(data.frame())
  })
  
  return(demographics_data)
}

# ── TABLE GENERATION FUNCTIONS ─────────────────────────────────────────────────

# Generate Table S1: Complete WAS Variable Summary
generate_table_s1 <- function(was_vars, metadata, all_data) {
  log_message("Generating Table S1: Complete WAS Variable Summary...")
  
  all_vars <- unlist(was_vars, use.names = FALSE) %>% unique()
  
  var_stats <- map_dfr(all_vars, ~ {
    var_name <- .x
    
    if (var_name %in% names(all_data)) {
      var_data <- all_data[[var_name]]
      n_total <- length(var_data)
      n_available <- sum(!is.na(var_data))
      availability_pct <- round((n_available / n_total) * 100, 1)
      
      var_type <- classify_variable_type(var_data)
      
      desc <- metadata$description[metadata$variable == var_name]
      if (length(desc) == 0) desc <- paste("Variable:", var_name)
      
      if (var_type == "Continuous") {
        stats <- sprintf("Mean: %.2f ± %.2f, Median: %.2f (%.2f-%.2f)",
                        mean(var_data, na.rm = TRUE),
                        sd(var_data, na.rm = TRUE) / sqrt(sum(!is.na(var_data))),
                        median(var_data, na.rm = TRUE),
                        quantile(var_data, 0.25, na.rm = TRUE),
                        quantile(var_data, 0.75, na.rm = TRUE))
      } else if (var_type == "Binary") {
        n_cases <- sum(var_data == 1, na.rm = TRUE)
        n_controls <- sum(var_data == 0, na.rm = TRUE)
        stats <- sprintf("Cases: %d (%.1f%%), Controls: %d (%.1f%%)",
                        n_cases, (n_cases/n_available)*100,
                        n_controls, (n_controls/n_available)*100)
      } else {
        mode_val <- names(sort(table(var_data, useNA = "no"), decreasing = TRUE))[1]
        mode_freq <- max(table(var_data, useNA = "no"))
        stats <- sprintf("Mode: %s (%d, %.1f%%)", mode_val, mode_freq, (mode_freq/n_available)*100)
      }
      
      was_analysis <- names(was_vars)[map_lgl(was_vars, ~ var_name %in% .x)]
      if (length(was_analysis) == 0) was_analysis <- "Unknown"
      
      tibble(
        Variable = var_name,
        NHANES_Code = var_name,
        Description = desc[1],
        WAS_Analysis = was_analysis[1],
        Variable_Type = var_type,
        N_Available = n_available,
        N_Total = n_total,
        Availability_Percent = availability_pct,
        Statistics = stats
      )
    } else {
      tibble(
        Variable = var_name,
        NHANES_Code = var_name,
        Description = paste("Variable:", var_name),
        WAS_Analysis = names(was_vars)[map_lgl(was_vars, ~ var_name %in% .x)][1],
        Variable_Type = "Not Found",
        N_Available = 0,
        N_Total = nrow(all_data),
        Availability_Percent = 0.0,
        Statistics = "Variable not found in database"
      )
    }
  })
  
  table_s1 <- var_stats %>%
    arrange(WAS_Analysis, desc(Availability_Percent)) %>%
    mutate(
      WAS_Analysis = case_when(
        WAS_Analysis == "1_demoWAS" ~ "Demographics & Socioeconomic",
        WAS_Analysis == "2_oradWAS" ~ "Oral Health Behaviors", 
        WAS_Analysis == "3_exWAS" ~ "Environmental Exposures & Diet",
        WAS_Analysis == "4_pheWAS" ~ "Clinical Phenotypes & Biomarkers",
        WAS_Analysis == "5_outWAS" ~ "Health Outcomes & Disease History",
        WAS_Analysis == "6_zimWAS" ~ "Laboratory Measurements",
        TRUE ~ WAS_Analysis
      )
    )
  
  log_message("  Table S1 generated with ", nrow(table_s1), " variables")
  return(table_s1)
}

# Generate Table S2: Detailed Binary Variable Analysis
generate_table_s2 <- function(all_data, metadata) {
  log_message("Generating Table S2: Detailed Binary Variable Analysis...")
  
  binary_vars <- names(all_data)[map_lgl(names(all_data), ~ {
    if (.x %in% c("SEQN", "cycle")) return(FALSE)
    var_data <- all_data[[.x]]
    if (all(is.na(var_data))) return(FALSE)
    unique_vals <- unique(var_data[!is.na(var_data)])
    length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
  })]
  
  log_message("  Binary variables identified: ", length(binary_vars))
  
  binary_stats <- map_dfr(binary_vars, ~ {
    var_name <- .x
    var_data <- all_data[[var_name]]
    
    n_total <- sum(!is.na(var_data))
    n_cases <- sum(var_data == 1, na.rm = TRUE)
    n_controls <- sum(var_data == 0, na.rm = TRUE)
    
    prevalence <- n_cases / n_total
    ci <- wilson_ci(n_cases, n_total)
    
    desc <- metadata$description[metadata$variable == var_name]
    if (length(desc) == 0) desc <- paste("Variable:", var_name)
    
    tibble(
      Variable = var_name,
      NHANES_Code = var_name,
      Description = desc[1],
      Cases_N = n_cases,
      Cases_Percent = sprintf("%.1f%%", prevalence * 100),
      Controls_N = n_controls,
      Controls_Percent = sprintf("%.1f%%", (n_controls/n_total) * 100),
      Total_N = n_total,
      Prevalence = sprintf("%.1f%%", prevalence * 100),
      Prevalence_95CI = sprintf("(%.1f%%-%.1f%%)", ci["lower"] * 100, ci["upper"] * 100)
    )
  })
  
  table_s2 <- binary_stats %>%
    arrange(desc(Cases_N))
  
  log_message("  Table S2 generated with ", nrow(table_s2), " binary variables")
  return(table_s2)
}

# Generate Table S3: Alpha Diversity Stratification Analysis
generate_table_s3 <- function(alpha_strat_results) {
  log_message("Generating Table S3: Alpha Diversity Stratification Analysis...")
  
  # Debug: Check what columns are available
  log_message("    Available columns: ", paste(names(alpha_strat_results), collapse = ", "))
  
  # Use the actual column names from the results
  table_s3 <- alpha_strat_results %>%
    arrange(Stratification_Variable, Alpha_Diversity_Metric)
  
  log_message("  Table S3 generated with ", nrow(table_s3), " stratification results")
  return(table_s3)
}

# Generate Table S4: Variable Type Distribution by Analysis
generate_table_s4 <- function(table_s1_data) {
  log_message("Generating Table S4: Variable Type Distribution by Analysis...")
  
  table_s4 <- table_s1_data %>%
    group_by(WAS_Analysis, Variable_Type) %>%
    summarise(
      Count = n(),
      Mean_Availability = round(mean(Availability_Percent, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = Variable_Type, values_from = Count, values_fill = 0) %>%
    mutate(
      Total_Variables = rowSums(select(., -WAS_Analysis, -Mean_Availability), na.rm = TRUE)
    ) %>%
    arrange(desc(Mean_Availability))
  
  log_message("  Table S4 generated with ", nrow(table_s4), " analysis types")
  return(table_s4)
}

# ── MAIN EXECUTION ─────────────────────────────────────────────────────────────

log_message("Starting NHANES Oral Microbiome Data Summary Generation...")

# Connect to database
log_message("🔗 Connecting to database...")
con <- dbConnect(SQLite(), DB_PATH)

tryCatch({
  # Load configurations
  was_vars <- load_was_variables()
  
  # Get microbiome participants
  microbiome_participants <- get_microbiome_participants(con)
  
  # Load demographics data from processed database
  demographics_data <- load_demographics_data(con)
  
  # Load alpha diversity data
  alpha_summary <- load_alpha_diversity_data()
  
  # Collect all variables
  all_vars <- unlist(was_vars, use.names = FALSE) %>% unique()
  log_message("Total unique variables across all WAS analyses: ", length(all_vars))
  
  # Get variable metadata
  metadata <- get_variable_metadata(con, all_vars)
  
  # Collect variable data
  collection_result <- collect_variable_data(con, all_vars, microbiome_participants)
  all_data <- collection_result$data
  variables_found <- collection_result$variables_found
  
  # Analyze alpha diversity stratification
  alpha_strat_results <- analyze_alpha_diversity_stratification(demographics_data, alpha_summary)
  
  # Generate tables
  log_message("Generating supplementary tables...")
  
  table_s1 <- generate_table_s1(was_vars, metadata, all_data)
  table_s2 <- generate_table_s2(all_data, metadata)
  table_s3 <- generate_table_s3(alpha_strat_results)
  table_s4 <- generate_table_s4(table_s1)
  
  # Save tables
  log_message("💾 Saving supplementary tables...")
  
  write_csv(table_s1, file.path(OUTPUT_DIR, "Table_S1_Complete_WAS_Variable_Summary.csv"))
  write_csv(table_s2, file.path(OUTPUT_DIR, "Table_S2_Binary_Variables_Analysis.csv"))
  write_csv(table_s3, file.path(OUTPUT_DIR, "Table_S3_Alpha_Diversity_Stratification.csv"))
  write_csv(table_s4, file.path(OUTPUT_DIR, "Table_S4_Variable_Type_Distribution.csv"))
  
  # Print summary
  log_message("✅ Analysis complete!")
  log_message("Summary Statistics:")
  log_message("  • Total participants with microbiome data: ", nrow(microbiome_participants))
  log_message("  • Total variables analyzed: ", length(all_vars))
  log_message("  • Variables successfully found: ", length(variables_found))
  log_message("  • Binary variables identified: ", nrow(table_s2))
  log_message("  • Alpha diversity stratifications: ", nrow(alpha_strat_results))
  
  # WAS analysis breakdown
  was_summary <- table_s1 %>%
    group_by(WAS_Analysis) %>%
    summarise(
      Variables = n(),
      Mean_Availability = round(mean(Availability_Percent, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Availability))
  
  log_message("📈 Data Availability by WAS Analysis:")
  for (i in seq_len(nrow(was_summary))) {
    log_message("  • ", was_summary$WAS_Analysis[i], ": ", 
                was_summary$Mean_Availability[i], "% mean availability (", 
                was_summary$Variables[i], " variables)")
  }
  
  # Summary report
  log_message("Tables Generated:")
  log_message("- Table S1: Complete WAS Variable Summary")
  log_message("- Table S2: Detailed Binary Variable Analysis") 
  log_message("- Table S3: Alpha Diversity Stratification Analysis")
  log_message("- Table S4: Variable Type Distribution by Analysis")
  log_message("\nFootnotes:")
  log_message("- Variable descriptions from NHANES database metadata")
  log_message("- Statistics shown as mean ± SE for continuous variables")
  log_message("- Binary variables show case/control distributions with 95% CI")
  log_message("- Alpha diversity metrics by demographic groups")
  log_message("- All analyses focus on participants with oral microbiome data")

}, finally = {
  if (exists("con") && dbIsValid(con)) {
    dbDisconnect(con)
    log_message("🔌 Database connection closed")
  }
}) 