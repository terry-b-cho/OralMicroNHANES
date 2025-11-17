#!/usr/bin/env Rscript
################################################################################
##  1. LOAD PRE-COMPUTED DIVERSITY DATA & CATEGORICAL VARIABLES               ##
##  - Loads alpha diversity from dada2rsv-alpha.txt                           ##
##  - Loads beta diversity matrices from dada2rsv-beta-*.txt                  ##
##  - Loads categorical variables from phyloseq sample_data                   ##
##  - Uses SMALLER SUBSET for testing (first 1000 samples)                    ##
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(tibble)
  library(DBI)
})

cat("=== LOADING PRE-COMPUTED DIVERSITY DATA ===\n\n")

###############################################################################
## 1. PATH CONFIGURATION                                                      ##
###############################################################################

base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
phyloseq_obj_path <- file.path(base_path, "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate")
db_path <- file.path(base_path, "data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite")
output_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")

# Create output directory
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

cat("Data path:", data_path, "\n")
cat("Phyloseq path:", phyloseq_obj_path, "\n")
cat("Database path:", db_path, "\n")
cat("Output path:", output_path, "\n\n")

###############################################################################
## 2. LOAD ALPHA DIVERSITY DATA                                              ##
###############################################################################

cat("Loading alpha diversity data...\n")

alpha_file <- file.path(data_path, "dada2rsv-alpha.txt")
if (file.exists(alpha_file)) {
  # Read with explicit column types
  alpha_data_raw <- read_delim(alpha_file, delim = "\t", show_col_types = FALSE,
                               col_types = cols(.default = col_double(), SEQN = col_double()))
  
  cat("  Alpha diversity loaded:", nrow(alpha_data_raw), "samples\n")
  cat("  Alpha diversity columns:", ncol(alpha_data_raw), "\n")
  
  # Extract the rarefaction depth 10000 metrics (most complete)
  # SEQN is numeric in this file
  alpha_diversity <- alpha_data_raw %>%
    select(
      SEQN,
      Observed_OTUs = RSV_ObservedOTUs_10000_0,
      Shannon_Diversity = RSV_ShanWienDiv_10000_0,
      Inverse_Simpson = RSV_InverseSimpson_10000_0
    ) %>%
    mutate(SEQN = as.character(SEQN))  # Convert to character for matching
  
  cat("  Extracted metrics at rarefaction depth 10000\n")
  cat("  SEQN type:", class(alpha_diversity$SEQN), "\n")
  cat("  Final alpha diversity:", nrow(alpha_diversity), "samples,", 
      ncol(alpha_diversity) - 1, "metrics\n")
  
  # USE SUBSET FOR TESTING
  cat("  Using subset of first 1000 samples for testing...\n")
  alpha_diversity <- alpha_diversity[1:min(1000, nrow(alpha_diversity)), ]
  cat("  Test subset:", nrow(alpha_diversity), "samples\n")
} else {
  stop("Alpha diversity file not found: ", alpha_file)
}

###############################################################################
## 3. LOAD CATEGORICAL VARIABLES FROM DATABASE                               ##
###############################################################################

cat("\nLoading categorical variables from database...\n")

# Connect to database
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# Get sample SEQNs from alpha diversity (as character)
test_seqns <- alpha_diversity$SEQN

# Load demographic data from database
demo_query <- sprintf("SELECT SEQN, RIDAGEYR, RIAGENDR FROM demo_f 
                       WHERE SEQN IN (%s)", 
                      paste(test_seqns, collapse = ", "))

demo_data <- dbGetQuery(con, demo_query) %>%
  mutate(SEQN = as.character(SEQN))

cat("  Demographic data loaded:", nrow(demo_data), "samples\n")
cat("  SEQN type:", class(demo_data$SEQN), "\n")

# Create age groups (following 5_age_analysis.Rmd pattern)
age_breaks <- c(14, seq(20, 85, by = 5))
age_labels <- c("14-19", paste(seq(20, 75, by = 5), seq(24, 79, by = 5), sep = "-"), "80-85")

demo_data <- demo_data %>%
  mutate(
    Age_group = cut(RIDAGEYR, breaks = age_breaks, labels = age_labels, right = FALSE),
    Gender = factor(ifelse(RIAGENDR == 1, "Male", "Female"), levels = c("Female", "Male"))
  )

cat("  Age groups created\n")
cat("  Age group distribution:\n")
print(table(demo_data$Age_group, useNA = "ifany"))

# Close database connection
dbDisconnect(con)

###############################################################################
## 4. MERGE ALPHA DIVERSITY WITH CATEGORICAL DATA                            ##
###############################################################################

cat("\nMerging alpha diversity with categorical data...\n")
cat("  Alpha diversity SEQN type:", class(alpha_diversity$SEQN), "\n")
cat("  Demo data SEQN type:", class(demo_data$SEQN), "\n")

alpha_with_categories <- alpha_diversity %>%
  left_join(demo_data, by = "SEQN") %>%
  filter(!is.na(Observed_OTUs), !is.na(Age_group))

cat("  Merged data:", nrow(alpha_with_categories), "samples\n")
cat("  Complete cases:", sum(complete.cases(alpha_with_categories)), "\n")
cat("  Columns:", paste(colnames(alpha_with_categories), collapse = ", "), "\n")

###############################################################################
## 5. LOAD BETA DIVERSITY DATA (SUBSET)                                      ##
###############################################################################

cat("\nLoading beta diversity matrices (SUBSET FOR TESTING)...\n")

# Function to load distance matrix with subset
load_distance_matrix_subset <- function(file_path, metric_name, sample_ids) {
  cat("  Loading", metric_name, "(subset)...")
  
  # Read only first few lines to get structure
  con <- file(file_path, "r")
  header <- readLines(con, n = 1)
  close(con)
  
  # For large files, we'll use data.table fread for speed
  library(data.table)
  dist_matrix <- fread(file_path, nrows = 1001)  # Read subset
  
  # First column should be sample IDs
  all_sample_ids <- as.character(dist_matrix[[1]])
  
  # Find indices of our test samples
  test_indices <- which(all_sample_ids %in% sample_ids)
  
  if (length(test_indices) == 0) {
    cat(" No matching samples found!\n")
    return(NULL)
  }
  
  cat(" Found", length(test_indices), "matching samples\n")
  
  # For testing, just create a small distance matrix
  # In production, you'd load the full subset
  dist_mat <- as.matrix(dist_matrix[test_indices, (test_indices + 1), with = FALSE])
  rownames(dist_mat) <- all_sample_ids[test_indices]
  colnames(dist_mat) <- all_sample_ids[test_indices]
  
  # Convert to dist object
  dist_obj <- as.dist(dist_mat)
  
  return(dist_obj)
}

# Load beta diversity for test samples
beta_diversity <- list()

# Note: Beta diversity files are very large (>500MB each)
# For testing, we create placeholder objects
cat("  Note: Beta diversity files are very large (>500MB)\n")
cat("  Creating placeholder for testing - REPLACE WITH FULL LOAD IN PRODUCTION\n")

# For now, create empty placeholders
beta_diversity$braycurtis <- NULL
beta_diversity$unwunifrac <- NULL  
beta_diversity$wunifrac <- NULL

cat("  Beta diversity: Using NULL placeholders for testing\n")

###############################################################################
## 6. SAVE PROCESSED DATA                                                    ##
###############################################################################

cat("\nSaving processed data...\n")

# Save alpha diversity with categories
saveRDS(alpha_with_categories, file.path(output_path, "alpha_diversity_with_categories.rds"))
cat("  ✓ Saved: alpha_diversity_with_categories.rds\n")

# Save beta diversity matrices
saveRDS(beta_diversity, file.path(output_path, "beta_diversity_matrices.rds"))
cat("  ✓ Saved: beta_diversity_matrices.rds (placeholders)\n")

# Save demographic data for reference
saveRDS(demo_data, file.path(output_path, "demographic_data.rds"))
cat("  ✓ Saved: demographic_data.rds\n")

###############################################################################
## 7. SUMMARY                                                                ##
###############################################################################

cat("\n=== DATA LOADING SUMMARY ===\n")
cat("Alpha diversity samples:", nrow(alpha_diversity), "\n")
cat("Merged dataset:", nrow(alpha_with_categories), "samples\n")
cat("Categorical variables: Age_group, Gender\n")
cat("\nAge group distribution:\n")
print(table(alpha_with_categories$Age_group, useNA = "ifany"))
cat("\nGender distribution:\n")
print(table(alpha_with_categories$Gender, useNA = "ifany"))
cat("\nData loading complete!\n")
cat("============================\n")
