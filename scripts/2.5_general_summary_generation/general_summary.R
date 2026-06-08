################################################################################
# general_summary.R
#
# Generates supplementary summary tables for the NHANES oral microbiome study
# in two phases:
#   Phase A -- WAS variable summary, binary variable analysis, alpha diversity
#              stratification, variable type distribution by analysis.
#   Phase B -- comprehensive WAS variable summary with table-source metadata,
#              WAS analysis summary by type, variable type distribution,
#              detailed binary variable analysis.
#
# Inputs (paths relative to PROJECT_ROOT):
#   - data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite
#   - data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt
#   - configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt
#
# Outputs (under results/analyses_results/2.5_general_summary_generation/out/):
#   Phase A:
#     - Table_S1_Complete_WAS_Variable_Summary.csv
#     - Table_S2_Binary_Variables_Analysis.csv
#     - Table_S3_Alpha_Diversity_Stratification.csv
#     - Table_S4_Variable_Type_Distribution.csv
#   Phase B:
#     - table_s1_complete_was_variable_summary.csv
#     - table_s2_was_analysis_summary.csv
#     - table_s3_variable_type_distribution.csv
#     - table_s4_binary_variable_analysis.csv
#
# Environment: R >= 4.5 with DBI, RSQLite, dplyr, readr, stringr, purrr,
# tidyr, survey, knitr. Exact versions in module2.5_tool_version_list.txt.
#
# Run:
#   Rscript scripts/2.5_general_summary_generation/general_summary.R
################################################################################

# === USER SETTING =============================================================
# Set PROJECT_ROOT to the absolute path of your local clone of this repository.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ==============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(survey)
  library(knitr)
})

# ---- paths ----
DB_PATH        <- file.path(PROJECT_ROOT,
                            "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
ALPHA_DIV_PATH <- file.path(PROJECT_ROOT,
                            "data/00_nhanes_omp_diversity_db/dada2rsv-alpha.txt")
CONFIG_DIR     <- file.path(PROJECT_ROOT, "configs")
OUTPUT_DIR     <- file.path(PROJECT_ROOT,
                            "results/analyses_results/2.5_general_summary_generation/out")
CYCLES         <- c("F", "G")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

WAS_CONFIG_FILES <- list(
  "1_demoWAS" = file.path(CONFIG_DIR, "1_demoWAS_vars.txt"),
  "2_oradWAS" = file.path(CONFIG_DIR, "2_oradWAS_vars.txt"),
  "3_exWAS"   = file.path(CONFIG_DIR, "3_exWAS_vars.txt"),
  "4_pheWAS"  = file.path(CONFIG_DIR, "4_pheWAS_vars.txt"),
  "5_outWAS"  = file.path(CONFIG_DIR, "5_outWAS_vars.txt")
)

# ---- shared helpers ----
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(lower = 0, upper = 0))
  z <- qnorm((1 + conf.level) / 2)
  p <- x / n
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  width  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(lower = max(0, center - width), upper = min(1, center + width))
}

# Phase A variable-type classifier (CamelCase labels).
classify_variable_type_A <- function(x) {
  if (all(is.na(x))) return("Missing")
  unique_vals <- unique(x[!is.na(x)])
  n_unique <- length(unique_vals)
  if (n_unique == 2 && all(unique_vals %in% c(0, 1))) return("Binary")
  if (is.character(x) || is.factor(x) || (is.numeric(x) && n_unique <= 10)) return("Categorical")
  if (is.ordered(x))                                                          return("Ordinal")
  if (is.numeric(x) && n_unique > 10)                                         return("Continuous")
  "Other"
}

# Phase B variable-type classifier (lowercase labels).
classify_variable_type_B <- function(data, variable, sql_type = "REAL") {
  if (!variable %in% names(data)) return("unknown")
  values <- data[[variable]][!is.na(data[[variable]])]
  if (length(values) == 0) return("unknown")
  unique_vals <- unique(values)
  if (length(unique_vals) <= 2) {
    if (all(unique_vals %in% c(0, 1)) || all(unique_vals %in% c(1, 2))) return("binary")
  }
  if (is.character(values) || length(unique_vals) <= 10) return("categorical")
  if (all(values == round(values)) && length(unique_vals) <= 20) return("ordinal")
  "continuous"
}

# ---- connect once ----
con <- dbConnect(SQLite(), DB_PATH)
on.exit({ if (exists("con") && dbIsValid(con)) dbDisconnect(con) }, add = TRUE)

# ==============================================================================
# PHASE A
# ==============================================================================

load_was_variables_A <- function() {
  was_vars <- map(WAS_CONFIG_FILES, ~ {
    if (file.exists(.x)) read_lines(.x) %>% str_trim() %>% .[. != ""]
    else                 character(0)
  })
  was_vars
}

get_microbiome_participants_A <- function(con) {
  q <- "
    SELECT DISTINCT SEQN, 'F' as cycle FROM DEMO_F WHERE SEQN IS NOT NULL
    UNION
    SELECT DISTINCT SEQN, 'G' as cycle FROM DEMO_G WHERE SEQN IS NOT NULL
  "
  d <- dbGetQuery(con, q)
  d$SEQN <- as.character(d$SEQN)
  d
}

load_demographics_data <- function(con) {
  demo_f <- dbGetQuery(con, "SELECT SEQN, RIDAGEYR, RIAGENDR, DMDEDUC2, RIDRETH1, INDFMPIR FROM DEMO_F")
  demo_g <- dbGetQuery(con, "SELECT SEQN, RIDAGEYR, RIAGENDR, DMDEDUC2, RIDRETH1, INDFMPIR FROM DEMO_G")
  demo_combined <- rbind(demo_f, demo_g)
  demo_combined %>%
    mutate(
      AgeQuartile = case_when(
        is.na(RIDAGEYR) ~ "Missing",
        RIDAGEYR <= quantile(RIDAGEYR, 0.25, na.rm = TRUE) ~ "Q1 (14-27)",
        RIDAGEYR <= quantile(RIDAGEYR, 0.50, na.rm = TRUE) ~ "Q2 (28-40)",
        RIDAGEYR <= quantile(RIDAGEYR, 0.75, na.rm = TRUE) ~ "Q3 (41-55)",
        TRUE ~ "Q4 (56-69)"
      ),
      Gender = case_when(RIAGENDR == 1 ~ "Male", RIAGENDR == 2 ~ "Female", TRUE ~ "Missing"),
      EducationLevel = case_when(
        DMDEDUC2 == 1 ~ "< 9th Grade",
        DMDEDUC2 == 2 ~ "9-11th Grade",
        DMDEDUC2 == 3 ~ "High School",
        DMDEDUC2 == 4 ~ "College/AA",
        DMDEDUC2 == 5 ~ "College Graduate",
        TRUE ~ "Missing"
      ),
      Ethnicity = case_when(
        RIDRETH1 == 1 ~ "Mexican American",
        RIDRETH1 == 2 ~ "Other Hispanic",
        RIDRETH1 == 3 ~ "Non-Hispanic White",
        RIDRETH1 == 4 ~ "Non-Hispanic Black",
        RIDRETH1 == 5 ~ "Other/Multi-racial",
        TRUE ~ "Missing"
      ),
      PIRQuartile = case_when(
        is.na(INDFMPIR) ~ "Missing",
        INDFMPIR <= quantile(INDFMPIR, 0.25, na.rm = TRUE) ~ "Q1 (Lowest)",
        INDFMPIR <= quantile(INDFMPIR, 0.50, na.rm = TRUE) ~ "Q2",
        INDFMPIR <= quantile(INDFMPIR, 0.75, na.rm = TRUE) ~ "Q3",
        TRUE ~ "Q4 (Highest)"
      )
    )
}

load_alpha_diversity_data <- function() {
  alpha_div_raw <- read.table(ALPHA_DIV_PATH, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE)
  if (is.character(alpha_div_raw$SEQN)) {
    alpha_div_raw$SEQN <- as.numeric(alpha_div_raw$SEQN)
  }
  col_names    <- names(alpha_div_raw)
  observed_cols<- grep("RSV_ObservedOTUs_10000_9",   col_names, value = TRUE)
  faith_cols   <- grep("RSV_FaPhyloDiv_10000_9",     col_names, value = TRUE)
  shannon_cols <- grep("RSV_ShanWienDiv_10000_9",    col_names, value = TRUE)
  simpson_cols <- grep("RSV_InverseSimpson_10000_9", col_names, value = TRUE)
  for (cl in c(observed_cols, faith_cols, shannon_cols, simpson_cols)) {
    if (cl %in% names(alpha_div_raw)) alpha_div_raw[[cl]] <- as.numeric(alpha_div_raw[[cl]])
  }
  alpha_div_raw %>%
    filter(!is.na(SEQN)) %>%
    group_by(SEQN) %>%
    summarise(
      Observed_ASVs_10000 = if (length(observed_cols) > 0)
        rowMeans(select(pick(everything()), all_of(observed_cols)), na.rm = TRUE) else NA_real_,
      Faith_PD_10000      = if (length(faith_cols) > 0)
        rowMeans(select(pick(everything()), all_of(faith_cols)),    na.rm = TRUE) else NA_real_,
      Shannon_10000       = if (length(shannon_cols) > 0)
        rowMeans(select(pick(everything()), all_of(shannon_cols)),  na.rm = TRUE) else NA_real_,
      InvSimpson_10000    = if (length(simpson_cols) > 0)
        rowMeans(select(pick(everything()), all_of(simpson_cols)),  na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}

analyze_alpha_diversity_stratification <- function(demographics_data, alpha_summary) {
  alpha_with_demo <- merge(alpha_summary, demographics_data, by = "SEQN", all.x = TRUE)
  strat_variables <- c("Gender" = "Gender", "Age Quartile" = "AgeQuartile",
                       "Education Level" = "EducationLevel", "Ethnicity" = "Ethnicity",
                       "PIR Quartile" = "PIRQuartile")
  alpha_metrics <- c("Observed_ASVs_10000", "Faith_PD_10000", "Shannon_10000", "InvSimpson_10000")
  metric_names  <- c("Observed ASVs", "Faith's PD", "Shannon Diversity", "Inverse Simpson")

  results <- list()
  for (i in seq_along(strat_variables)) {
    strat_var  <- strat_variables[i]
    strat_name <- names(strat_variables)[i]
    if (!strat_var %in% names(alpha_with_demo)) next
    for (j in seq_along(alpha_metrics)) {
      metric      <- alpha_metrics[j]
      metric_name <- metric_names[j]
      strat_data <- alpha_with_demo %>%
        filter(!is.na(.data[[metric]]) & !is.na(.data[[strat_var]])) %>%
        group_by(.data[[strat_var]]) %>%
        summarise(
          N      = n(),
          Mean   = mean(.data[[metric]], na.rm = TRUE),
          SE     = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
          Median = median(.data[[metric]], na.rm = TRUE),
          Q1     = quantile(.data[[metric]], 0.25, na.rm = TRUE),
          Q3     = quantile(.data[[metric]], 0.75, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          Stratification_Variable = strat_name,
          Alpha_Diversity_Metric  = metric_name,
          Mean_SE    = sprintf("%.3f ± %.3f", Mean, SE),
          Median_IQR = sprintf("%.3f (%.3f-%.3f)", Median, Q1, Q3),
          Category   = .data[[strat_var]]
        ) %>%
        select(Stratification_Variable, Category, Alpha_Diversity_Metric, N, Mean_SE, Median_IQR)
      results[[paste(strat_name, metric_name, sep = "_")]] <- strat_data
    }
  }
  bind_rows(results)
}

get_variable_metadata_A <- function(con, variables) {
  var_metadata <- dbGetQuery(con, "
    SELECT `Data.File.Name` as nhanes_code,
           `Data.File.Description` as description
    FROM variable_names_epcf
  ")
  tibble(
    variable    = variables,
    nhanes_code = variables,
    description = map_chr(variables, ~ {
      desc <- var_metadata$description[var_metadata$nhanes_code == .x]
      if (length(desc) > 0 && !is.na(desc[1]) && desc[1] != "") desc[1]
      else                                                     paste("Variable:", .x)
    })
  )
}

collect_variable_data_A <- function(con, variables, microbiome_participants) {
  tables   <- dbListTables(con)
  tables   <- tables[!grepl("^(sqlite_|variable_names)", tables)]
  # Memoize per-table column lists so the inner loop is a hash lookup, not a DB roundtrip.
  table_fields <- setNames(
    lapply(tables, function(t) tryCatch(dbListFields(con, t), error = function(e) character(0))),
    tables
  )
  all_data <- microbiome_participants %>% select(SEQN, cycle)
  variables_found <- character(0)
  for (var in variables) {
    for (table in tables) {
      columns <- table_fields[[table]]
      if (var %in% columns && "SEQN" %in% columns) {
        tryCatch({
          query <- sprintf("SELECT SEQN, `%s` as `%s` FROM `%s` WHERE `%s` IS NOT NULL",
                           var, var, table, var)
          var_data <- dbGetQuery(con, query)
          var_data$SEQN <- as.character(var_data$SEQN)
          if (nrow(var_data) > 0) {
            all_data <- all_data %>%
              left_join(var_data, by = "SEQN", suffix = c("", paste0(".", table)))
            variables_found <- c(variables_found, var)
            break
          }
        }, error = function(e) NULL)
      }
    }
  }
  list(data = all_data, variables_found = variables_found)
}

generate_table_s1_A <- function(was_vars, metadata, all_data) {
  all_vars <- unlist(was_vars, use.names = FALSE) %>% unique()
  var_stats <- map_dfr(all_vars, ~ {
    var_name <- .x
    if (var_name %in% names(all_data)) {
      var_data <- all_data[[var_name]]
      n_total      <- length(var_data)
      n_available  <- sum(!is.na(var_data))
      availability <- round((n_available / n_total) * 100, 1)
      var_type <- classify_variable_type_A(var_data)
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
        n_cases    <- sum(var_data == 1, na.rm = TRUE)
        n_controls <- sum(var_data == 0, na.rm = TRUE)
        stats <- sprintf("Cases: %d (%.1f%%), Controls: %d (%.1f%%)",
                         n_cases, (n_cases / n_available) * 100,
                         n_controls, (n_controls / n_available) * 100)
      } else {
        mode_val  <- names(sort(table(var_data, useNA = "no"), decreasing = TRUE))[1]
        mode_freq <- max(table(var_data, useNA = "no"))
        stats <- sprintf("Mode: %s (%d, %.1f%%)", mode_val, mode_freq,
                         (mode_freq / n_available) * 100)
      }
      was_analysis <- names(was_vars)[map_lgl(was_vars, ~ var_name %in% .x)]
      if (length(was_analysis) == 0) was_analysis <- "Unknown"
      tibble(Variable = var_name, NHANES_Code = var_name, Description = desc[1],
             WAS_Analysis = was_analysis[1], Variable_Type = var_type,
             N_Available = n_available, N_Total = n_total,
             Availability_Percent = availability, Statistics = stats)
    } else {
      tibble(Variable = var_name, NHANES_Code = var_name,
             Description = paste("Variable:", var_name),
             WAS_Analysis = names(was_vars)[map_lgl(was_vars, ~ var_name %in% .x)][1],
             Variable_Type = "Not Found", N_Available = 0L, N_Total = nrow(all_data),
             Availability_Percent = 0.0, Statistics = "Variable not found in database")
    }
  })
  var_stats %>%
    arrange(WAS_Analysis, desc(Availability_Percent)) %>%
    mutate(
      WAS_Analysis = case_when(
        WAS_Analysis == "1_demoWAS" ~ "Demographics & Socioeconomic",
        WAS_Analysis == "2_oradWAS" ~ "Oral Health Behaviors",
        WAS_Analysis == "3_exWAS"   ~ "Environmental Exposures & Diet",
        WAS_Analysis == "4_pheWAS"  ~ "Clinical Phenotypes & Biomarkers",
        WAS_Analysis == "5_outWAS"  ~ "Health Outcomes & Disease History",
        TRUE                        ~ WAS_Analysis
      )
    )
}

generate_table_s2_A <- function(all_data, metadata) {
  binary_vars <- names(all_data)[map_lgl(names(all_data), ~ {
    if (.x %in% c("SEQN", "cycle")) return(FALSE)
    var_data <- all_data[[.x]]
    if (all(is.na(var_data))) return(FALSE)
    unique_vals <- unique(var_data[!is.na(var_data)])
    length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
  })]
  binary_stats <- map_dfr(binary_vars, ~ {
    var_name <- .x
    var_data <- all_data[[var_name]]
    n_total      <- sum(!is.na(var_data))
    n_cases      <- sum(var_data == 1, na.rm = TRUE)
    n_controls   <- sum(var_data == 0, na.rm = TRUE)
    prevalence   <- n_cases / n_total
    ci <- wilson_ci(n_cases, n_total)
    desc <- metadata$description[metadata$variable == var_name]
    if (length(desc) == 0) desc <- paste("Variable:", var_name)
    tibble(
      Variable        = var_name, NHANES_Code = var_name, Description = desc[1],
      Cases_N         = n_cases, Cases_Percent     = sprintf("%.1f%%", prevalence * 100),
      Controls_N      = n_controls, Controls_Percent = sprintf("%.1f%%", (n_controls / n_total) * 100),
      Total_N         = n_total,    Prevalence       = sprintf("%.1f%%", prevalence * 100),
      Prevalence_95CI = sprintf("(%.1f%%-%.1f%%)", ci["lower"] * 100, ci["upper"] * 100)
    )
  })
  binary_stats %>% arrange(desc(Cases_N))
}

generate_table_s3_A <- function(alpha_strat_results) {
  alpha_strat_results %>% arrange(Stratification_Variable, Alpha_Diversity_Metric)
}

generate_table_s4_A <- function(table_s1_data) {
  table_s1_data %>%
    group_by(WAS_Analysis, Variable_Type) %>%
    summarise(Count = n(),
              Mean_Availability = round(mean(Availability_Percent, na.rm = TRUE), 1),
              .groups = "drop") %>%
    pivot_wider(names_from = Variable_Type, values_from = Count, values_fill = 0) %>%
    mutate(Total_Variables = rowSums(select(., -WAS_Analysis, -Mean_Availability),
                                     na.rm = TRUE)) %>%
    arrange(desc(Mean_Availability))
}

# ---- Phase A run ----
cat("[Phase A] start\n")
was_vars_A             <- load_was_variables_A()
microbiome_participants <- get_microbiome_participants_A(con)
demographics_data      <- load_demographics_data(con)
alpha_summary          <- load_alpha_diversity_data()
all_vars_A             <- unlist(was_vars_A, use.names = FALSE) %>% unique()
metadata_A             <- get_variable_metadata_A(con, all_vars_A)
cat("[Phase A] collecting", length(all_vars_A), "variables\n")
collection_A           <- collect_variable_data_A(con, all_vars_A, microbiome_participants)
all_data_A             <- collection_A$data
alpha_strat_results    <- analyze_alpha_diversity_stratification(demographics_data, alpha_summary)

table_s1_A <- generate_table_s1_A(was_vars_A, metadata_A, all_data_A)
table_s2_A <- generate_table_s2_A(all_data_A, metadata_A)
table_s3_A <- generate_table_s3_A(alpha_strat_results)
table_s4_A <- generate_table_s4_A(table_s1_A)

write_csv(table_s1_A, file.path(OUTPUT_DIR, "Table_S1_Complete_WAS_Variable_Summary.csv"))
write_csv(table_s2_A, file.path(OUTPUT_DIR, "Table_S2_Binary_Variables_Analysis.csv"))
write_csv(table_s3_A, file.path(OUTPUT_DIR, "Table_S3_Alpha_Diversity_Stratification.csv"))
write_csv(table_s4_A, file.path(OUTPUT_DIR, "Table_S4_Variable_Type_Distribution.csv"))

# ==============================================================================
# PHASE B
# ==============================================================================

load_was_variables_B <- function() {
  was_vars <- list()
  for (analysis_name in names(WAS_CONFIG_FILES)) {
    cfg <- WAS_CONFIG_FILES[[analysis_name]]
    if (file.exists(cfg)) {
      vars <- readLines(cfg, warn = FALSE)
      vars <- vars[vars != "" & !startsWith(vars, "#")]
      was_vars[[analysis_name]] <- vars
    }
  }
  list(by_analysis = was_vars, all_variables = unique(unlist(was_vars)))
}

get_microbiome_participants_B <- function(con) {
  microbiome_seqn <- list()
  for (cycle in CYCLES) {
    genus_table <- paste0("DADA2RSV_GENUS_RELATIVE_", cycle, "_none")
    if (genus_table %in% dbListTables(con)) {
      microbiome_seqn[[cycle]] <- tbl(con, genus_table) %>%
        select(SEQN) %>% distinct() %>% collect() %>%
        mutate(SEQN = as.character(SEQN), cycle = cycle)
    }
  }
  if (length(microbiome_seqn) == 0) stop("No microbiome data found in database")
  bind_rows(microbiome_seqn)
}

get_database_metadata_B <- function(con, variables) {
  all_tables      <- dbListTables(con)
  relevant_tables <- all_tables[str_detect(all_tables,
                                           paste0("(", paste(CYCLES, collapse = "|"), ")$"))]
  metadata_list <- list()
  for (table_name in relevant_tables) {
    tryCatch({
      table_info <- dbGetQuery(con, sprintf("PRAGMA table_info(%s)", table_name))
      if (nrow(table_info) > 0) {
        table_metadata <- tibble(variable = table_info$name,
                                 table_source = table_name,
                                 sql_type = table_info$type) %>%
          filter(variable %in% variables, variable != "SEQN")
        if (nrow(table_metadata) > 0) metadata_list[[table_name]] <- table_metadata
      }
    }, error = function(e) NULL)
  }
  if (length(metadata_list) > 0) {
    bind_rows(metadata_list) %>% group_by(variable) %>% slice_head(n = 1) %>% ungroup()
  } else {
    tibble(variable = character(), table_source = character(), sql_type = character())
  }
}

collect_comprehensive_data_B <- function(con, microbiome_participants, variables) {
  demo_data <- list()
  for (cycle in CYCLES) {
    demo_table <- paste0("DEMO_", cycle)
    if (demo_table %in% dbListTables(con)) {
      demo_data[[cycle]] <- tbl(con, demo_table) %>% collect() %>%
        mutate(SEQN = as.character(SEQN), cycle = cycle)
    }
  }
  combined_data    <- bind_rows(demo_data)
  other_tables     <- dbListTables(con)
  relevant_pattern <- paste0("(", paste(CYCLES, collapse = "|"), ")$")
  relevant_tables  <- other_tables[str_detect(other_tables, relevant_pattern)]
  for (table_name in relevant_tables) {
    if (!str_detect(table_name, "DEMO_")) {
      tryCatch({
        table_data <- tbl(con, table_name) %>% collect() %>%
          mutate(SEQN = as.character(SEQN))
        available_vars <- intersect(variables, names(table_data))
        if (length(available_vars) > 0) {
          keep_vars <- c("SEQN", available_vars)
          table_data <- table_data[, names(table_data) %in% keep_vars, drop = FALSE]
          combined_data <- combined_data %>%
            left_join(table_data, by = "SEQN", suffix = c("", "_new"))
          for (var in available_vars) {
            new_col <- paste0(var, "_new")
            if (new_col %in% names(combined_data)) {
              combined_data[[var]] <- coalesce(combined_data[[var]], combined_data[[new_col]])
              combined_data[[new_col]] <- NULL
            }
          }
        }
      }, error = function(e) NULL)
    }
  }
  microbiome_seqn <- unique(microbiome_participants$SEQN)
  combined_data %>% filter(SEQN %in% microbiome_seqn)
}

generate_comprehensive_stats_B <- function(data, variable, db_metadata, was_analysis = NULL) {
  var_metadata <- db_metadata[db_metadata$variable == variable, ]
  n_total       <- nrow(data)
  n_available   <- sum(!is.na(data[[variable]]))
  n_missing     <- n_total - n_available
  pct_available <- round(100 * n_available / n_total, 1)
  sql_type <- ifelse(nrow(var_metadata) > 0, var_metadata$sql_type[1], "REAL")
  var_type <- classify_variable_type_B(data, variable, sql_type)

  if (n_available == 0) {
    statistic <- "No data available"; additional_stats <- ""
  } else {
    values <- data[[variable]][!is.na(data[[variable]])]
    if (var_type == "continuous") {
      mean_val   <- mean(values); sd_val <- sd(values); se_val <- sd_val / sqrt(length(values))
      median_val <- median(values); q25 <- quantile(values, 0.25); q75 <- quantile(values, 0.75)
      statistic        <- sprintf("%.2f ± %.2f", mean_val, se_val)
      additional_stats <- sprintf("Median: %.2f (IQR: %.2f-%.2f)", median_val, q25, q75)
    } else if (var_type == "binary") {
      freq_table  <- table(values, useNA = "no")
      unique_vals <- names(freq_table)
      positive_val <- if (all(unique_vals %in% c(0, 1))) "1"
                      else if (all(unique_vals %in% c(1, 2))) "1"
                      else unique_vals[1]
      positive_n  <- ifelse(positive_val %in% names(freq_table), freq_table[positive_val], 0)
      positive_pct <- round(100 * positive_n / n_available, 1)
      statistic        <- sprintf("Cases: %d (%.1f%%)", positive_n, positive_pct)
      additional_stats <- sprintf("Controls: %d (%.1f%%)",
                                  n_available - positive_n, 100 - positive_pct)
    } else if (var_type %in% c("categorical", "ordinal")) {
      freq_table   <- table(values, useNA = "no")
      n_categories <- length(freq_table)
      most_common  <- names(freq_table)[which.max(freq_table)]
      most_common_n   <- max(freq_table)
      most_common_pct <- round(100 * most_common_n / n_available, 1)
      statistic        <- sprintf("Mode: %s (%d, %.1f%%)",
                                  most_common, most_common_n, most_common_pct)
      additional_stats <- sprintf("%d categories total", n_categories)
    } else {
      statistic <- "Unknown data type"; additional_stats <- ""
    }
  }
  table_source <- ifelse(nrow(var_metadata) > 0, var_metadata$table_source[1], "Unknown")
  tibble(
    nhanes_code   = variable,
    description   = ifelse(variable %in% names(data), "Available in dataset", "Not found"),
    was_analysis  = ifelse(is.null(was_analysis), "Extended", was_analysis),
    table_source  = table_source,
    variable_type = var_type,
    n_total       = n_total,
    n_available   = n_available,
    n_missing     = n_missing,
    percent_available = pct_available,
    statistic     = statistic,
    additional_stats = additional_stats
  )
}

generate_binary_analysis_B <- function(data, binary_vars, db_metadata) {
  binary_results <- list()
  for (var in binary_vars) {
    if (var %in% names(data)) {
      values <- data[[var]][!is.na(data[[var]])]
      if (length(values) > 0) {
        freq_table  <- table(values, useNA = "no")
        unique_vals <- sort(unique(values))
        if (all(unique_vals %in% c(0, 1))) {
          cases <- sum(values == 1, na.rm = TRUE); controls <- sum(values == 0, na.rm = TRUE)
          case_label <- "1 (Yes/Positive)";        control_label <- "0 (No/Negative)"
        } else if (all(unique_vals %in% c(1, 2))) {
          cases <- sum(values == 1, na.rm = TRUE); controls <- sum(values == 2, na.rm = TRUE)
          case_label <- "1 (Yes)";                  control_label <- "2 (No)"
        } else {
          cases <- freq_table[1]; controls <- sum(freq_table[-1])
          case_label <- paste0(names(freq_table)[1], " (Reference)"); control_label <- "Other"
        }
        total_valid <- cases + controls
        case_pct    <- round(100 * cases / total_valid, 1)
        control_pct <- round(100 - case_pct, 1)
        prevalence  <- cases / total_valid
        se_prev     <- sqrt(prevalence * (1 - prevalence) / total_valid)
        ci_lower    <- max(0, prevalence - 1.96 * se_prev)
        ci_upper    <- min(1, prevalence + 1.96 * se_prev)
        binary_results[[var]] <- tibble(
          nhanes_code = var,
          total_n     = nrow(data),
          valid_n     = total_valid,
          missing_n   = nrow(data) - total_valid,
          cases_n     = cases,    cases_pct           = case_pct,
          controls_n  = controls, controls_pct        = control_pct,
          prevalence  = round(prevalence, 3),
          prevalence_ci_lower = round(ci_lower, 3),
          prevalence_ci_upper = round(ci_upper, 3),
          case_label  = case_label, control_label    = control_label
        )
      }
    }
  }
  if (length(binary_results) > 0) bind_rows(binary_results) else tibble()
}

cat("[Phase A] done\n")

# ---- Phase B run ----
cat("[Phase B] start\n")
was_variables_B  <- load_was_variables_B()
microbiome_participants_B <- get_microbiome_participants_B(con)
db_metadata_B    <- get_database_metadata_B(con, was_variables_B$all_variables)
cat("[Phase B] collecting", length(was_variables_B$all_variables), "variables\n")
participant_data <- collect_comprehensive_data_B(con, microbiome_participants_B,
                                                 was_variables_B$all_variables)

was_stats_list <- list()
for (analysis_name in names(was_variables_B$by_analysis)) {
  for (var in was_variables_B$by_analysis[[analysis_name]]) {
    was_stats_list[[paste(analysis_name, var, sep = "_")]] <-
      generate_comprehensive_stats_B(participant_data, var, db_metadata_B, analysis_name)
  }
}

table_s1_B <- bind_rows(was_stats_list) %>%
  arrange(was_analysis, nhanes_code) %>%
  select(
    `WAS Analysis`         = was_analysis,
    `NHANES Code`          = nhanes_code,
    `Variable Description` = description,
    `Table Source`         = table_source,
    `Variable Type`        = variable_type,
    `N Total`              = n_total,
    `N Available`          = n_available,
    `N Missing`            = n_missing,
    `% Available`          = percent_available,
    `Primary Statistic`    = statistic,
    `Additional Statistics`= additional_stats
  )

table_s2_B <- bind_rows(was_stats_list) %>%
  group_by(was_analysis) %>%
  summarise(
    `Total Variables`            = n(),
    `Variables Available (>=50%)` = sum(percent_available >= 50),
    `Variables Available (>=75%)` = sum(percent_available >= 75),
    `Variables Available (>=90%)` = sum(percent_available >= 90),
    `Mean Availability (%)`       = round(mean(percent_available), 1),
    `Median Availability (%)`     = round(median(percent_available), 1),
    `Continuous Variables`        = sum(variable_type == "continuous"),
    `Binary Variables`            = sum(variable_type == "binary"),
    `Categorical Variables`       = sum(variable_type == "categorical"),
    `Ordinal Variables`           = sum(variable_type == "ordinal"),
    .groups = "drop"
  ) %>%
  arrange(desc(`Mean Availability (%)`))

table_s3_B <- bind_rows(was_stats_list) %>%
  group_by(variable_type) %>%
  summarise(
    Count                    = n(),
    Percentage               = round(100 * n() / nrow(bind_rows(was_stats_list)), 1),
    `Mean Availability (%)`  = round(mean(percent_available), 1),
    `Variables >=90% Available` = sum(percent_available >= 90),
    .groups = "drop"
  ) %>%
  arrange(desc(Count))

binary_vars_B <- bind_rows(was_stats_list) %>%
  filter(variable_type == "binary") %>%
  pull(nhanes_code)

write_csv(table_s1_B, file.path(OUTPUT_DIR, "table_s1_complete_was_variable_summary.csv"))
write_csv(table_s2_B, file.path(OUTPUT_DIR, "table_s2_was_analysis_summary.csv"))
write_csv(table_s3_B, file.path(OUTPUT_DIR, "table_s3_variable_type_distribution.csv"))

if (length(binary_vars_B) > 0) {
  binary_analysis <- generate_binary_analysis_B(participant_data, binary_vars_B, db_metadata_B)
  if (nrow(binary_analysis) > 0) {
    table_s4_B <- binary_analysis %>%
      select(
        `NHANES Code`        = nhanes_code,
        `Total N`            = total_n,
        `Valid N`            = valid_n,
        `Missing N`          = missing_n,
        `Cases N`            = cases_n,
        `Cases %`            = cases_pct,
        `Controls N`         = controls_n,
        `Controls %`         = controls_pct,
        Prevalence           = prevalence,
        `95% CI Lower`       = prevalence_ci_lower,
        `95% CI Upper`       = prevalence_ci_upper,
        `Case Definition`    = case_label,
        `Control Definition` = control_label
      )
    write_csv(table_s4_B, file.path(OUTPUT_DIR, "table_s4_binary_variable_analysis.csv"))
  }
}

dbDisconnect(con)
cat("[Phase B] done\n")
