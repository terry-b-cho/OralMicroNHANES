#!/usr/bin/env Rscript
# -----------------------------------------------------------
# missing_tables_n_variables_diagnostic.R
# Comprehensive diagnostic script to check what tables, variables, 
# and derived data already exist in the NHANES SQLite database.
# This prevents redundant data creation and processing.
# -----------------------------------------------------------
# Usage:
#   Rscript scripts/2_preprocess_db_n_phyloseq/missing_tables_n_variables_diagnostic.R \
#     --db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
#     --output_dir results/diagnostics/
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(knitr)
  library(tibble)
  library(fs)
  library(readr)
})

# Create logs directory
logs_dir <- "scripts/2_preprocess_db_n_phyloseq/logs"
dir_create(logs_dir, recurse = TRUE)

# ── CLI options ─────────────────────────────────────────────
opt_list <- list(
  make_option(c("-d","--db"), type = "character",
              default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite",
              help = "Path to SQLite database"),
  make_option(c("-o","--output_dir"), type = "character", 
              default = "results/analyses_results/02_preprocess_db_n_phyloseq_out/diagnostics/",
              help = "Output directory for diagnostic reports"),
  make_option(c("-f","--format"), type = "character",
              default = "html",
              help = "Output format: html, csv, or both"),
  make_option(c("--processed"), action = "store_true", default = FALSE,
              help = "Check processed database instead of original [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Override database path if processed flag is used
if (opt$processed) {
  opt$db <- "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite"
  # Change output directory to processed subdirectory
  opt$output_dir <- gsub("/diagnostics/$", "/diagnostics_processed/", opt$output_dir)
}

# Create output directory
dir_create(opt$output_dir, recurse = TRUE)

# Set up logging
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(logs_dir, paste0("diagnostic_", timestamp, ".log"))
log_con <- file(log_file, open = "wt")

# Function to log both to console and file
log_message <- function(...) {
  msg <- paste(..., sep = "")
  message(msg)
  writeLines(paste(Sys.time(), ":", msg), log_con)
  flush(log_con)
}

# Connect to database
log_message("Connecting to database: ", opt$db)
con <- dbConnect(SQLite(), opt$db)

# ── DEFINE EXPECTED TABLES AND VARIABLES ──────────────────────

# Define all variables that might be derived/created
expected_derived_vars <- list(
  # Demographics original derived variables
  demographics_original = c("AGE_SQUARED", "BORN_INUSA", "EDUCATION_LESS9", "EDUCATION_9_11", 
                           "EDUCATION_AA", "EDUCATION_COLLEGEGRAD", "ETHNICITY_MEXICAN", 
                           "ETHNICITY_OTHERHISPANIC", "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", 
                           "RIAGENDR_01"),
  
  # Comprehensive demographic variables (new processing pipeline)
  demographics_comprehensive = c(
    # Core identifiers
    "sample", "cycle", "Release_Cycle", "RIAGENDR_01", "Gender",
    # Age variables
    "AgeGroup", "age_group", "Age", "AgeQuartile",
    # Education variables  
    "EducationLevel", "Education_Level", "EducationLevel_Simple",
    # Ethnicity and citizenship
    "Ethnicity", "BornInUSA", "US_Born", "US_Citizen",
    # Income and poverty
    "PIRcat", "PIRQuartile", "Ratio_Family_Income_Poverty",
    # Household variables
    "Household_Size", "Family_Size", "Household_Size_Factor",
    # Marital status
    "Marital_Status",
    # Interview variables
    "Interview_Language", "Interview_Proxy", "Interview_Interpreter",
    # Household reference person
    "HH_Reference_Gender", "HH_Reference_Age", "HH_Reference_Education_Level", "HH_Reference_Marital_Status",
    # Detailed income categories
    "Annual_Household_Income", "Annual_Family_Income"
  ),
  
  # Oral health derived variables  
  oral_health = c("TOOTH_DECAY_OHAROCDT", "GUM_DISEASE_OHAROCGP", 
                  "ORAL_HYGIENE_OHAROCOH", "DENTURE_OHAROCDE"),
  
  # Disease outcome derived variables
  cardiovascular = c("CHD", "STROKE", "HEART_ATTACK", "HEART_FAILURE", "ANGINA", "CVD"),
  respiratory = c("ASTHMA", "BRONCHITIS", "EMPHYSEMA"),
  cancer = c("CANCER_BREAST", "CANCER_COLON", "CANCER_LUNG", "CANCER_ESOPHAGEAL", 
             "CANCER_PROSTATE", "CANCER_MOUTH"),
  metabolic = c("DIABETES"),
  
  # Lab measurements
  lab_measures = c("LBXP1", "LBXP2", "LBDP3", "LBDRFO", "LBXSF6SI", "LBXSCK",
                   "LBXFOLSI", "LBXPS4", "LBXTBM", "BMXSUB", "LBDFOT", 
                   "LBXMMASI", "BMXTRI", "URX1DC")
)

# Define expected derived tables
expected_derived_tables <- list(
  oral_health = c("OralDisease_F", "OralDisease_G", "d_outcome_oh_f", "d_outcome_oh_g"),
  disease_outcomes = c("d_outcome_mcq_f", "d_outcome_mcq_g"),
  microbiome_transformed = c(
    # CLR transformed
    "DADA2RSV_GENUS_RELATIVE_F_clr", "DADA2RSV_GENUS_RELATIVE_G_clr",
    "DADA2RSV_FAMILY_RELATIVE_F_clr", "DADA2RSV_FAMILY_RELATIVE_G_clr",
    "DADA2RSV_ORDER_RELATIVE_F_clr", "DADA2RSV_ORDER_RELATIVE_G_clr",
    "DADA2RSV_CLASS_RELATIVE_F_clr", "DADA2RSV_CLASS_RELATIVE_G_clr",
    "DADA2RSV_PHYLUM_RELATIVE_F_clr", "DADA2RSV_PHYLUM_RELATIVE_G_clr",
    # Log-normal transformed
    "DADA2RSV_GENUS_RELATIVE_F_lognorm", "DADA2RSV_GENUS_RELATIVE_G_lognorm",
    "DADA2RSV_FAMILY_RELATIVE_F_lognorm", "DADA2RSV_FAMILY_RELATIVE_G_lognorm",
    "DADA2RSV_ORDER_RELATIVE_F_lognorm", "DADA2RSV_ORDER_RELATIVE_G_lognorm",
    "DADA2RSV_CLASS_RELATIVE_F_lognorm", "DADA2RSV_CLASS_RELATIVE_G_lognorm",
    "DADA2RSV_PHYLUM_RELATIVE_F_lognorm", "DADA2RSV_PHYLUM_RELATIVE_G_lognorm",
    # None (raw) transformed
    "DADA2RSV_GENUS_RELATIVE_F_none", "DADA2RSV_GENUS_RELATIVE_G_none",
    "DADA2RSV_FAMILY_RELATIVE_F_none", "DADA2RSV_FAMILY_RELATIVE_G_none",
    "DADA2RSV_ORDER_RELATIVE_F_none", "DADA2RSV_ORDER_RELATIVE_G_none",
    "DADA2RSV_CLASS_RELATIVE_F_none", "DADA2RSV_CLASS_RELATIVE_G_none",
    "DADA2RSV_PHYLUM_RELATIVE_F_none", "DADA2RSV_PHYLUM_RELATIVE_G_none"
  )
)

# ── DIAGNOSTIC FUNCTIONS ──────────────────────────────────────

check_database_structure <- function(con) {
  log_message("Checking database structure...")
  
  # Get all tables
  all_tables <- dbListTables(con)
  
  # Get metadata tables info
  metadata_tables <- c("variable_names_epcf", "table_names_epcf")
  metadata_exists <- metadata_tables %in% all_tables
  
  # Get table sizes
  table_info <- map_dfr(all_tables, ~{
    tryCatch({
      count_query <- sprintf("SELECT COUNT(*) as row_count FROM %s", .x)
      row_count <- dbGetQuery(con, count_query)$row_count
      
      cols_query <- sprintf("PRAGMA table_info(%s)", .x)
      col_count <- nrow(dbGetQuery(con, cols_query))
      
      tibble(
        table_name = .x,
        row_count = row_count,
        column_count = col_count,
        table_type = case_when(
          str_detect(.x, "_clr$|_lognorm$|_none$") ~ "transformed_microbiome",
          str_detect(.x, "^d_outcome_") ~ "derived_outcome",
          str_detect(.x, "OralDisease_") ~ "derived_oral",
          .x %in% metadata_tables ~ "metadata",
          TRUE ~ "original"
        )
      )
    }, error = function(e) {
      tibble(table_name = .x, row_count = NA, column_count = NA, table_type = "error")
    })
  })
  
  return(list(
    total_tables = length(all_tables),
    metadata_exists = metadata_exists,
    table_info = table_info
  ))
}

check_expected_tables <- function(con, expected_tables) {
  log_message("Checking expected derived tables...")
  
  all_tables <- dbListTables(con)
  
  table_status <- map_dfr(names(expected_tables), ~{
    category <- .x
    tables <- expected_tables[[category]]
    
    map_dfr(tables, function(tbl) {
      exists <- tbl %in% all_tables
      row_count <- if(exists) {
        tryCatch(dbGetQuery(con, sprintf("SELECT COUNT(*) as cnt FROM %s", tbl))$cnt,
                 error = function(e) NA)
      } else NA
      
      tibble(
        category = category,
        table_name = tbl,
        exists = exists,
        row_count = row_count,
        status = ifelse(exists, "✅ EXISTS", "❌ MISSING")
      )
    })
  })
  
  return(table_status)
}

check_expected_variables <- function(con, expected_vars) {
  log_message("Checking expected derived variables...")
  
  # Get all tables and their fields
  all_tables <- dbListTables(con)
  
  fields_by_table <- setNames(
    map(all_tables, ~dbListFields(con, .x)),
    all_tables
  )
  
  # Check metadata table for variable info
  metadata_vars <- if("variable_names_epcf" %in% all_tables) {
    dbReadTable(con, "variable_names_epcf") %>%
      select(Variable.Name, Data.File.Name, Variable.Description) %>%
      rename(variable_name = Variable.Name, table_name = Data.File.Name)
  } else {
    tibble(variable_name = character(), table_name = character(), Variable.Description = character())
  }
  
  # Check each expected variable
  var_status <- map_dfr(names(expected_vars), ~{
    category <- .x
    vars <- expected_vars[[category]]
    
    map_dfr(vars, function(var) {
      # Find which tables contain this variable
      tables_with_var <- names(fields_by_table)[
        map_lgl(fields_by_table, ~var %in% .x)
      ]
      
      # Check if in metadata
      in_metadata <- var %in% metadata_vars$variable_name
      metadata_info <- if(in_metadata) {
        filter(metadata_vars, variable_name == var)
      } else {
        tibble(table_name = NA, Variable.Description = NA)
      }
      
      tibble(
        category = category,
        variable_name = var,
        exists_in_tables = length(tables_with_var) > 0,
        table_locations = paste(tables_with_var, collapse = ", "),
        in_metadata = in_metadata,
        metadata_description = if(nrow(metadata_info) > 0) metadata_info$Variable.Description[1] else NA,
        status = case_when(
          length(tables_with_var) > 0 & in_metadata ~ "✅ COMPLETE",
          length(tables_with_var) > 0 & !in_metadata ~ "⚠️ EXISTS BUT NO METADATA",
          length(tables_with_var) == 0 ~ "❌ MISSING"
        )
      )
    })
  })
  
  return(var_status)
}

check_microbiome_transformations <- function(con) {
  log_message("Checking microbiome transformation completeness...")
  
  # Expected transformations
  taxa_levels <- c("GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM")
  cycles <- c("F", "G")
  transforms <- c("clr", "lognorm", "none")
  
  all_tables <- dbListTables(con)
  
  transformation_status <- expand_grid(
    taxa_level = taxa_levels,
    cycle = cycles,
    transform = transforms
  ) %>%
    mutate(
      table_name = sprintf("DADA2RSV_%s_RELATIVE_%s_%s", taxa_level, cycle, transform),
      exists = table_name %in% all_tables,
      row_count = map_dbl(table_name, ~{
        if(.x %in% all_tables) {
          tryCatch(dbGetQuery(con, sprintf("SELECT COUNT(*) as cnt FROM %s", .x))$cnt,
                   error = function(e) NA_real_)
        } else NA_real_
      }),
      column_count = map_dbl(table_name, ~{
        if(.x %in% all_tables) {
          length(dbListFields(con, .x))
        } else NA_real_
      }),
      status = ifelse(exists, "✅ EXISTS", "❌ MISSING")
    )
  
  return(transformation_status)
}

check_demographic_processing <- function(con) {
  log_message("👥 Checking demographic variable processing in DEMO_F and DEMO_G...")
  
  # Expected comprehensive demographic variables
  comprehensive_vars <- c(
    # Core identifiers
    "sample", "cycle", "Release_Cycle", "RIAGENDR_01", "Gender",
    # Age variables
    "AgeGroup", "age_group", "Age", "AgeQuartile",
    # Education variables  
    "EducationLevel", "Education_Level", "EducationLevel_Simple",
    # Ethnicity and citizenship
    "Ethnicity", "BornInUSA", "US_Born", "US_Citizen",
    # Income and poverty
    "PIRcat", "PIRQuartile", "Ratio_Family_Income_Poverty",
    # Household variables
    "Household_Size", "Family_Size", "Household_Size_Factor",
    # Marital status
    "Marital_Status",
    # Interview variables
    "Interview_Language", "Interview_Proxy", "Interview_Interpreter",
    # Household reference person
    "HH_Reference_Gender", "HH_Reference_Age", "HH_Reference_Education_Level", "HH_Reference_Marital_Status",
    # Detailed income categories
    "Annual_Household_Income", "Annual_Family_Income"
  )
  
  all_tables <- dbListTables(con)
  
  # Check metadata for demographic variables
  metadata_vars <- if("variable_names_epcf" %in% all_tables) {
    dbReadTable(con, "variable_names_epcf") %>%
      filter(Data.File.Name %in% c("DEMO_F", "DEMO_G")) %>%
      select(Variable.Name, Data.File.Name, Variable.Description) %>%
      rename(variable_name = Variable.Name, table_name = Data.File.Name)
  } else {
    tibble(variable_name = character(), table_name = character(), Variable.Description = character())
  }
  
  # Check each table and variable combination
  demo_status <- expand_grid(
    table_name = c("DEMO_F", "DEMO_G"),
    variable_name = comprehensive_vars
  ) %>%
    mutate(
      table_exists = table_name %in% all_tables,
      variable_exists = case_when(
        !table_exists ~ FALSE,
        TRUE ~ map2_lgl(table_name, variable_name, ~{
          if(.x %in% all_tables) {
            .y %in% dbListFields(con, .x)
          } else FALSE
        })
      ),
      in_metadata = paste(variable_name, table_name) %in% paste(metadata_vars$variable_name, metadata_vars$table_name),
      status = case_when(
        !table_exists ~ "❌ TABLE MISSING",
        variable_exists & in_metadata ~ "✅ COMPLETE",
        variable_exists & !in_metadata ~ "⚠️ EXISTS BUT NO METADATA",
        !variable_exists ~ "❌ VARIABLE MISSING"
      )
    ) %>%
    arrange(table_name, variable_name)
  
  return(demo_status)
}

generate_summary_report <- function(db_structure, table_status, var_status, transformation_status, demographic_status = NULL) {
  log_message("Generating summary report...")
  
  # Overall summary
  summary_stats <- list(
    total_tables = db_structure$total_tables,
    metadata_tables_exist = all(db_structure$metadata_exists),
    
    # Table summaries
    expected_tables_total = nrow(table_status),
    expected_tables_exist = sum(table_status$exists),
    expected_tables_missing = sum(!table_status$exists),
    
    # Variable summaries  
    expected_vars_total = nrow(var_status),
    expected_vars_complete = sum(var_status$status == "✅ COMPLETE"),
    expected_vars_exist_no_meta = sum(var_status$status == "⚠️ EXISTS BUT NO METADATA"),
    expected_vars_missing = sum(var_status$status == "❌ MISSING"),
    
    # Transformation summaries
    transformations_total = nrow(transformation_status),
    transformations_exist = sum(transformation_status$exists),
    transformations_missing = sum(!transformation_status$exists)
  )
  
  # Add demographic summaries if available
  if (!is.null(demographic_status)) {
    summary_stats$demographic_vars_total <- nrow(demographic_status)
    summary_stats$demographic_vars_complete <- sum(demographic_status$status == "✅ COMPLETE")
    summary_stats$demographic_vars_exist_no_meta <- sum(demographic_status$status == "⚠️ EXISTS BUT NO METADATA")
    summary_stats$demographic_vars_missing <- sum(demographic_status$status == "❌ VARIABLE MISSING")
    summary_stats$demographic_table_missing <- sum(demographic_status$status == "❌ TABLE MISSING")
  }
  
  return(summary_stats)
}

create_action_plan <- function(table_status, var_status, transformation_status, demographic_status = NULL) {
  log_message("📝 Creating action plan...")
  
  action_items <- list()
  
  # Missing tables
  missing_tables <- filter(table_status, !exists)
  if(nrow(missing_tables) > 0) {
    action_items$missing_tables <- missing_tables %>%
      group_by(category) %>%
      summarise(
        tables_needed = paste(table_name, collapse = ", "),
        .groups = "drop"
      )
  }
  
  # Missing variables
  missing_vars <- filter(var_status, status == "❌ MISSING")
  if(nrow(missing_vars) > 0) {
    action_items$missing_variables <- missing_vars %>%
      group_by(category) %>%
      summarise(
        variables_needed = paste(variable_name, collapse = ", "),
        .groups = "drop"
      )
  }
  
  # Variables without metadata
  vars_no_meta <- filter(var_status, status == "⚠️ EXISTS BUT NO METADATA")
  if(nrow(vars_no_meta) > 0) {
    action_items$variables_need_metadata <- vars_no_meta %>%
      select(variable_name, table_locations)
  }
  
  # Missing transformations
  missing_transforms <- filter(transformation_status, !exists)
  if(nrow(missing_transforms) > 0) {
    action_items$missing_transformations <- missing_transforms %>%
      select(taxa_level, cycle, transform, table_name)
  }
  
  # Add demographic action items if available
  if (!is.null(demographic_status)) {
    # Missing demographic variables
    missing_demo_vars <- filter(demographic_status, status == "❌ VARIABLE MISSING")
    if(nrow(missing_demo_vars) > 0) {
      action_items$missing_demographic_variables <- missing_demo_vars %>%
        select(table_name, variable_name, status)
    }
    
    # Demographic variables without metadata
    demo_vars_no_meta <- filter(demographic_status, status == "⚠️ EXISTS BUT NO METADATA")
    if(nrow(demo_vars_no_meta) > 0) {
      action_items$demographic_variables_need_metadata <- demo_vars_no_meta %>%
        select(table_name, variable_name, status)
    }
  }
  
  return(action_items)
}

# ── MAIN EXECUTION ─────────────────────────────────────────

log_message("Starting comprehensive database diagnostic...")
log_message("Database: ", opt$db)
log_message("Output directory: ", opt$output_dir)

# Run all diagnostics
db_structure <- check_database_structure(con)
table_status <- check_expected_tables(con, expected_derived_tables)
var_status <- check_expected_variables(con, expected_derived_vars)
transformation_status <- check_microbiome_transformations(con)
demographic_status <- check_demographic_processing(con)

# Generate reports
summary_stats <- generate_summary_report(db_structure, table_status, var_status, transformation_status, demographic_status)
action_plan <- create_action_plan(table_status, var_status, transformation_status, demographic_status)

# ── OUTPUT RESULTS ─────────────────────────────────────────

# Print summary to console
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("DATABASE DIAGNOSTIC SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat(sprintf("Total tables in database: %d\n", summary_stats$total_tables))
cat(sprintf("Metadata tables present: %s\n", ifelse(summary_stats$metadata_tables_exist, "✅ YES", "❌ NO")))
cat("\n")
cat("EXPECTED TABLES:\n")
cat(sprintf("  Total expected: %d\n", summary_stats$expected_tables_total))
cat(sprintf("  ✅ Exist: %d\n", summary_stats$expected_tables_exist))
cat(sprintf("  ❌ Missing: %d\n", summary_stats$expected_tables_missing))
cat("\n")
cat("EXPECTED VARIABLES:\n")
cat(sprintf("  Total expected: %d\n", summary_stats$expected_vars_total))
cat(sprintf("  ✅ Complete (exist + metadata): %d\n", summary_stats$expected_vars_complete))
cat(sprintf("  ⚠️ Exist but no metadata: %d\n", summary_stats$expected_vars_exist_no_meta))
cat(sprintf("  ❌ Missing: %d\n", summary_stats$expected_vars_missing))
cat("\n")
cat("MICROBIOME TRANSFORMATIONS:\n")
cat(sprintf("  Total expected: %d\n", summary_stats$transformations_total))
cat(sprintf("  ✅ Exist: %d\n", summary_stats$transformations_exist))
cat(sprintf("  ❌ Missing: %d\n", summary_stats$transformations_missing))
cat("\n")

# Add demographic processing summary if available
if (!is.null(summary_stats$demographic_vars_total)) {
  cat("👥 DEMOGRAPHIC PROCESSING (DEMO_F & DEMO_G):\n")
  cat(sprintf("  Total expected variables: %d\n", summary_stats$demographic_vars_total))
  cat(sprintf("  ✅ Complete (exist + metadata): %d\n", summary_stats$demographic_vars_complete))
  cat(sprintf("  ⚠️ Exist but no metadata: %d\n", summary_stats$demographic_vars_exist_no_meta))
  cat(sprintf("  ❌ Variables missing: %d\n", summary_stats$demographic_vars_missing))
  if (summary_stats$demographic_table_missing > 0) {
    cat(sprintf("  ❌ Table missing: %d\n", summary_stats$demographic_table_missing))
  }
  cat("\n")
}

# Save detailed results
if(opt$format %in% c("csv", "both")) {
  write_csv(db_structure$table_info, file.path(opt$output_dir, "database_table_info.csv"))
  write_csv(table_status, file.path(opt$output_dir, "expected_tables_status.csv"))
  write_csv(var_status, file.path(opt$output_dir, "expected_variables_status.csv"))
  write_csv(transformation_status, file.path(opt$output_dir, "transformation_status.csv"))
  
  # Save demographic status if available
  if (!is.null(demographic_status)) {
    write_csv(demographic_status, file.path(opt$output_dir, "demographic_processing_status.csv"))
  }
  
  # Save action plan
  if(length(action_plan) > 0) {
    iwalk(action_plan, ~{
      write_csv(.x, file.path(opt$output_dir, paste0("action_", .y, ".csv")))
    })
  }
  
  log_message("✅ CSV reports saved to: ", opt$output_dir)
}

if(opt$format %in% c("html", "both")) {
  # Create HTML report (simplified version)
  html_content <- paste0(
    "<h1>NHANES Database Diagnostic Report</h1>",
    "<h2>Summary</h2>",
    "<ul>",
    "<li>Total tables: ", summary_stats$total_tables, "</li>",
    "<li>Expected tables exist: ", summary_stats$expected_tables_exist, "/", summary_stats$expected_tables_total, "</li>",
    "<li>Expected variables complete: ", summary_stats$expected_vars_complete, "/", summary_stats$expected_vars_total, "</li>",
    "<li>Transformations exist: ", summary_stats$transformations_exist, "/", summary_stats$transformations_total, "</li>",
    "</ul>",
    "<p>Generated: ", Sys.time(), "</p>"
  )
  
  writeLines(html_content, file.path(opt$output_dir, "diagnostic_report.html"))
  log_message("✅ HTML report saved to: ", file.path(opt$output_dir, "diagnostic_report.html"))
}

# Close database connection
dbDisconnect(con)

# Save diagnostic log summary
log_summary_file <- file.path(logs_dir, paste0("diagnostic_summary_", timestamp, ".txt"))
writeLines(c(
  "DATABASE DIAGNOSTIC SUMMARY",
  paste("Timestamp:", Sys.time()),
  paste("Database:", opt$db),
  paste("Output directory:", opt$output_dir),
  paste("Format:", opt$format),
  "",
  "RESULTS:",
  paste("Total tables:", summary_stats$total_tables),
  paste("Expected tables exist:", summary_stats$expected_tables_exist, "/", summary_stats$expected_tables_total),
  paste("Expected variables complete:", summary_stats$expected_vars_complete, "/", summary_stats$expected_vars_total),
  paste("Microbiome transformations exist:", summary_stats$transformations_exist, "/", summary_stats$transformations_total),
  "",
  "FILES GENERATED:",
  if(opt$format %in% c("csv", "both")) paste("- CSV reports in", opt$output_dir) else "- No CSV reports",
  if(opt$format %in% c("html", "both")) paste("- HTML report in", opt$output_dir) else "- No HTML reports"
), log_summary_file)

log_message("📄 Diagnostic summary saved: ", log_summary_file)
log_message(paste(rep("=", 60), collapse=""))
log_message("✅ Diagnostic complete! Check output directory for detailed reports.")
log_message("📁 Log file saved: ", log_file)
log_message(paste(rep("=", 60), collapse=""))

close(log_con)

cat(paste(rep("=", 60), collapse=""), "\n")
cat("✅ Diagnostic complete! Check output directory for detailed reports.\n")
cat("✅ All logs saved to: ", logs_dir, "\n")
cat(paste(rep("=", 60), collapse=""), "\n") 