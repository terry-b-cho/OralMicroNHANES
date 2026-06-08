#!/usr/bin/env Rscript

# =============================================================================
# ADDITIONAL SUPPLEMENTARY DISTRIBUTION GENERATION
# =============================================================================
# Four end-to-end phases, all writing under one supplementary dir:
#   Phase 1 — full variable distribution table for all WAS-config variables
#   Phase 2 — alpha + beta diversity group means per categorical level
#   Phase 3 — same, plus Wilcoxon vs reference level (effect size + BH-FDR)
#   Phase 4 — host x OTU regression results supplementary table (+ CLR subset)
#
# All outputs are written under:
#   results/analyses_results/2.5_general_summary_generation/supplementary/
#
# Environment: R >= 4.5 with data.table, dplyr, tidyr, DBI, RSQLite, dbplyr,
# purrr, stringr.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(DBI)
  library(RSQLite)
  library(dbplyr)
  library(purrr)
  library(stringr)
  library(parallel)
})

# Parallel CPU: cap workers so we stay within node limits.
N_CORES <- min(12L, max(1L, parallel::detectCores(logical = TRUE) - 2L))

base_dir         <- PROJECT_ROOT
config_dir       <- file.path(base_dir, "configs")
sqlite_db_path   <- file.path(base_dir, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
diversity_dir   <- file.path(base_dir, "data/00_nhanes_omp_diversity_db")
categorical_dir <- file.path(base_dir, "archive_final_before_submission/archive_result_final/6_archive/archive/data")
module4_inter   <- file.path(base_dir, "results/analyses_results/4_otu_host_plot_out/intermediate")

OUTPUT_DIR <- file.path(base_dir, "results/analyses_results/2.5_general_summary_generation/supplementary")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

ensure_exists <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  path
}

# =============================================================================
# PHASE 1 — full variable distribution table
# =============================================================================

message("=== PHASE 1: full variable distribution table ===")

determine_variable_type <- function(x) {
  if (is.null(x) || length(x) == 0) return("unknown")
  if (all(is.na(x))) return("all_na")
  if (is.factor(x))  return(if (nlevels(x) == 2) "binary" else "categorical")
  if (is.logical(x)) return("binary")
  if (is.numeric(x)) {
    non_na <- unique(x[!is.na(x)])
    if (length(non_na) <= 2 && all(non_na %in% c(0, 1))) "binary" else "continuous"
  } else if (is.character(x)) "categorical" else "unknown"
}

compute_continuous_stats <- function(x) {
  n_non_na <- sum(!is.na(x))
  missing_prop <- sum(is.na(x)) / length(x)
  if (n_non_na == 0) {
    return(list(n_non_na = 0L, missing_prop = missing_prop,
                mean = NA_real_, sd = NA_real_, median = NA_real_,
                q1 = NA_real_, q3 = NA_real_,
                n_levels = NA_integer_, most_common_level = NA_character_, most_common_prop = NA_real_))
  }
  list(
    n_non_na          = as.integer(n_non_na),
    missing_prop      = missing_prop,
    mean              = mean(x, na.rm = TRUE),
    sd                = sd(x, na.rm = TRUE),
    median            = median(x, na.rm = TRUE),
    q1                = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    q3                = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    n_levels          = NA_integer_,
    most_common_level = NA_character_,
    most_common_prop  = NA_real_
  )
}

compute_categorical_stats <- function(x) {
  n_non_na <- sum(!is.na(x))
  missing_prop <- sum(is.na(x)) / length(x)
  if (n_non_na == 0) {
    return(list(n_non_na = 0L, missing_prop = missing_prop,
                mean = NA_real_, sd = NA_real_, median = NA_real_,
                q1 = NA_real_, q3 = NA_real_,
                n_levels = NA_integer_, most_common_level = NA_character_, most_common_prop = NA_real_))
  }
  x <- factor(x, exclude = NA)
  tab <- table(x, useNA = "no")
  if (length(tab) == 0) {
    return(list(n_non_na = as.integer(n_non_na), missing_prop = missing_prop,
                mean = NA_real_, sd = NA_real_, median = NA_real_,
                q1 = NA_real_, q3 = NA_real_,
                n_levels = 0L, most_common_level = NA_character_, most_common_prop = NA_real_))
  }
  idx <- which.max(tab)
  list(
    n_non_na          = as.integer(n_non_na),
    missing_prop      = missing_prop,
    mean = NA_real_, sd = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_,
    n_levels          = as.integer(length(tab)),
    most_common_level = as.character(names(tab)[idx]),
    most_common_prop  = as.numeric(tab[idx] / n_non_na)
  )
}

build_var_table_mapping <- function(con, var_names, cycles = c("F", "G")) {
  var2table <- list()
  for (cyc in cycles) {
    lc <- tolower(cyc)
    var_list <- paste0("'", var_names, "'", collapse = ",")
    query <- sprintf(
      "SELECT \"Variable.Name\", \"Data.File.Name\"
       FROM variable_names_epcf
       WHERE \"Variable.Name\" IN (%s)
         AND lower(\"Data.File.Name\") LIKE '%%%s'",
      var_list, paste0("_", lc))
    mapping_df <- dbGetQuery(con, query)
    if (nrow(mapping_df) > 0) {
      v <- setNames(mapping_df$`Data.File.Name`, mapping_df$`Variable.Name`)
      missing <- setdiff(var_names, names(v))
      if (length(missing) > 0) v[missing] <- NA_character_
      var2table[[cyc]] <- v[var_names]
    } else {
      var2table[[cyc]] <- setNames(rep(NA_character_, length(var_names)), var_names)
    }
  }
  var2table
}

load_vars_from_mapping <- function(con, var_names, var2table, cycles = c("F", "G"), filter_seqns = NULL) {
  all_tables <- dbListTables(con)
  result_list <- list()
  for (cyc in cycles) {
    vars_by_table <- split(var_names, var2table[[cyc]][var_names])
    vars_by_table <- vars_by_table[!is.na(names(vars_by_table))]
    cycle_data_list <- list()
    for (table_name in names(vars_by_table)) {
      if (!table_name %in% all_tables) next
      vars_in_table <- vars_by_table[[table_name]]
      vars_in_table <- vars_in_table[!is.na(vars_in_table)]
      if (length(vars_in_table) == 0) next
      tryCatch({
        q <- tbl(con, table_name) %>% select(SEQN, any_of(vars_in_table))
        if (!is.null(filter_seqns)) q <- q %>% filter(SEQN %in% filter_seqns)
        td <- q %>% collect() %>% mutate(SEQN = as.character(SEQN))
        if (nrow(td) > 0 && length(intersect(vars_in_table, colnames(td))) > 0) {
          cycle_data_list[[table_name]] <- td
        }
      }, error = function(e) NULL)
    }
    if (length(cycle_data_list) > 0) {
      result_list[[cyc]] <- reduce(cycle_data_list, function(x, y) full_join(x, y, by = "SEQN"))
    }
  }
  if (length(result_list) == 0) return(NULL)
  combined <- bind_rows(result_list)
  combined <- combined[!duplicated(combined$SEQN), ]
  if (!is.null(filter_seqns)) combined <- combined %>% filter(SEQN %in% filter_seqns)
  combined
}

config_files <- list(
  demoWAS = file.path(config_dir, "1_demoWAS_vars.txt"),
  oradWAS = file.path(config_dir, "2_oradWAS_vars.txt"),
  exWAS   = file.path(config_dir, "3_exWAS_vars.txt"),
  pheWAS  = file.path(config_dir, "4_pheWAS_vars.txt"),
  outWAS  = file.path(config_dir, "5_outWAS_vars.txt")
)
all_vars_list <- list()
for (group_name in names(config_files)) {
  if (file.exists(config_files[[group_name]])) {
    vars <- readLines(config_files[[group_name]], warn = FALSE)
    vars <- vars[nzchar(vars) & !grepl("^#", vars)]
    all_vars_list[[group_name]] <- vars
  } else {
    all_vars_list[[group_name]] <- character(0)
  }
}
all_vars <- unique(unlist(all_vars_list))

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = ensure_exists(sqlite_db_path))

microbiome_seqns_f <- tbl(con, "DADA2RSV_GENUS_RELATIVE_F_clr") %>%
  select(SEQN) %>% distinct() %>% collect() %>%
  mutate(SEQN = as.character(SEQN)) %>% pull(SEQN)
microbiome_seqns_g <- tbl(con, "DADA2RSV_GENUS_RELATIVE_G_clr") %>%
  select(SEQN) %>% distinct() %>% collect() %>%
  mutate(SEQN = as.character(SEQN)) %>% pull(SEQN)
microbiome_seqns <- unique(c(microbiome_seqns_f, microbiome_seqns_g))

var_descriptions_db <- tbl(con, "variable_names_epcf") %>%
  select("Variable.Name", "Variable.Description", "Data.File.Name") %>%
  collect()

var2table <- build_var_table_mapping(con, all_vars, cycles = c("F", "G"))
all_data  <- load_vars_from_mapping(con, all_vars, var2table, cycles = c("F", "G"),
                                    filter_seqns = microbiome_seqns)
if (is.null(all_data) || nrow(all_data) == 0) stop("Failed to load any variables from database")

n_unique_loaded <- length(unique(all_data$SEQN))
if (nrow(all_data) > n_unique_loaded) {
  all_data <- all_data[!duplicated(all_data$SEQN), ]
}

derived_vars <- c("AGE_SQUARED", "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA",
                  "EDUCATION_COLLEGEGRAD", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC",
                  "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "BORN_INUSA")
for (dv in derived_vars) {
  if (dv %in% all_vars && !dv %in% colnames(all_data)) {
    if (dv == "AGE_SQUARED" && "RIDAGEYR" %in% colnames(all_data)) {
      all_data[[dv]] <- all_data$RIDAGEYR^2
    } else if (dv %in% c("EDUCATION_LESS9","EDUCATION_9_11","EDUCATION_AA","EDUCATION_COLLEGEGRAD") &&
               "DMDEDUC2" %in% colnames(all_data)) {
      val <- switch(dv,
        EDUCATION_LESS9       = 1, EDUCATION_9_11       = 2,
        EDUCATION_AA          = 4, EDUCATION_COLLEGEGRAD = 5)
      all_data[[dv]] <- as.numeric(all_data$DMDEDUC2 == val)
    } else if (dv %in% c("ETHNICITY_MEXICAN","ETHNICITY_OTHERHISPANIC","ETHNICITY_OTHER","ETHNICITY_NONHISPANICBLACK") &&
               "RIDRETH1" %in% colnames(all_data)) {
      val <- switch(dv,
        ETHNICITY_MEXICAN          = 1, ETHNICITY_OTHERHISPANIC    = 2,
        ETHNICITY_OTHER            = 5, ETHNICITY_NONHISPANICBLACK = 4)
      all_data[[dv]] <- as.numeric(all_data$RIDRETH1 == val)
    } else if (dv == "BORN_INUSA" && "DMDCITZN" %in% colnames(all_data)) {
      all_data[[dv]] <- as.numeric(all_data$DMDCITZN == 1)
    }
  }
}

output <- data.frame(
  variable_name = all_vars, variable_group = NA_character_,
  variable_description = NA_character_, source = NA_character_,
  levels = NA_character_, reference_level = NA_character_,
  n_non_na = NA_integer_, missing_prop = NA_real_,
  mean = NA_real_, sd = NA_real_, median = NA_real_,
  q1 = NA_real_, q3 = NA_real_,
  n_levels = NA_integer_, most_common_level = NA_character_, most_common_prop = NA_real_,
  stringsAsFactors = FALSE
)
for (group_name in names(all_vars_list)) {
  output$variable_group[output$variable_name %in% all_vars_list[[group_name]]] <- group_name
}

for (i in seq_along(all_vars)) {
  var_name <- all_vars[i]
  rows <- var_descriptions_db[var_descriptions_db$`Variable.Name` == var_name, ]
  if (nrow(rows) > 0) {
    output$variable_description[i] <- rows$`Variable.Description`[1]
    src <- rows$`Data.File.Name`[1]
    if (!is.na(src)) output$source[i] <- sub("_[FG]$", "", src, ignore.case = TRUE)
  }
  if (!var_name %in% colnames(all_data)) next
  x <- all_data[[var_name]]
  var_type <- determine_variable_type(x)
  if (var_type %in% c("binary", "categorical")) {
    x_fac <- if (is.factor(x)) factor(x, exclude = NA) else factor(x, exclude = NA)
    if (length(levels(x_fac)) > 0) {
      output$levels[i]          <- paste(levels(x_fac), collapse = "; ")
      output$reference_level[i] <- levels(x_fac)[1]
    }
  }
  tryCatch({
    stats <- if (var_type == "continuous") compute_continuous_stats(x)
             else if (var_type %in% c("binary","categorical")) compute_categorical_stats(x)
             else if (var_type == "all_na")
               list(n_non_na = 0L, missing_prop = 1.0,
                    mean = NA_real_, sd = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_,
                    n_levels = NA_integer_, most_common_level = NA_character_, most_common_prop = NA_real_)
             else list(n_non_na = NA_integer_, missing_prop = NA_real_,
                       mean = NA_real_, sd = NA_real_, median = NA_real_, q1 = NA_real_, q3 = NA_real_,
                       n_levels = NA_integer_, most_common_level = NA_character_, most_common_prop = NA_real_)
    output$n_non_na[i]          <- stats$n_non_na
    output$missing_prop[i]      <- stats$missing_prop
    output$mean[i]              <- stats$mean
    output$sd[i]                <- stats$sd
    output$median[i]            <- stats$median
    output$q1[i]                <- stats$q1
    output$q3[i]                <- stats$q3
    output$n_levels[i]          <- stats$n_levels
    output$most_common_level[i] <- stats$most_common_level
    output$most_common_prop[i]  <- stats$most_common_prop
  }, error = function(e) NULL)
}
DBI::dbDisconnect(con)

write.csv(output,
          file.path(OUTPUT_DIR, "full_list_supplementary_variable_distribution.csv"),
          row.names = FALSE)
message("Phase 1 done. ", nrow(output), " variables written.")

# =============================================================================
# PHASE 2 + 3 — alpha / beta diversity group means (+ optional ref comparisons)
# =============================================================================

message("=== PHASE 2/3: alpha + beta diversity group means and ref comparisons ===")

# Reuse microbiome_seqns from Phase 1
alpha_file <- ensure_exists(file.path(diversity_dir, "dada2rsv-alpha.txt"))
alpha_data_raw <- fread(alpha_file)

obs_cols  <- grep("RSV_ObservedOTUs_10000_",  colnames(alpha_data_raw), value = TRUE)
shan_cols <- grep("RSV_ShanWienDiv_10000_",   colnames(alpha_data_raw), value = TRUE)
simp_cols <- grep("RSV_InverseSimpson_10000_", colnames(alpha_data_raw), value = TRUE)

alpha_diversity <- alpha_data_raw %>% mutate(SEQN = as.character(SEQN))
mean_cols <- function(cols) {
  if (length(cols) == 0) return(rep(NA_real_, nrow(alpha_diversity)))
  m <- as.matrix(alpha_diversity[, cols, with = FALSE]); m <- apply(m, 2, as.numeric)
  rowMeans(m, na.rm = TRUE)
}
alpha_diversity$Observed_OTUs     <- mean_cols(obs_cols)
alpha_diversity$Shannon_Diversity <- mean_cols(shan_cols)
alpha_diversity$Inverse_Simpson   <- mean_cols(simp_cols)
alpha_diversity <- alpha_diversity %>%
  select(SEQN, Observed_OTUs, Shannon_Diversity, Inverse_Simpson) %>%
  filter(SEQN %in% microbiome_seqns)

categorical_data <- readRDS(ensure_exists(file.path(categorical_dir, "all_categorical_data.rds")))
categorical_data$SEQN <- as.character(categorical_data$SEQN)
categorical_data <- categorical_data %>% filter(SEQN %in% microbiome_seqns)

alpha_with_categories <- alpha_diversity %>% inner_join(categorical_data, by = "SEQN")

categorical_vars <- names(alpha_with_categories)[
  sapply(alpha_with_categories, function(x) is.factor(x) && nlevels(x) >= 1)
]
categorical_vars <- setdiff(categorical_vars, c("SEQN", "sample", "cycle", "Release_Cycle"))

load_beta_matrix <- function(file_path, sample_ids) {
  if (!file.exists(file_path)) stop("Beta diversity file not found: ", file_path)
  dist_data <- fread(file_path)
  all_seqns <- as.character(dist_data[[1]])
  keep_idx <- which(all_seqns %in% sample_ids)
  if (length(keep_idx) == 0) return(NULL)
  dist_mat <- as.matrix(dist_data[keep_idx, (keep_idx + 1), with = FALSE])
  rownames(dist_mat) <- all_seqns[keep_idx]; colnames(dist_mat) <- all_seqns[keep_idx]
  dist_mat
}

sample_ids   <- alpha_with_categories$SEQN
beta_metrics <- list(
  Braycurtis = load_beta_matrix(file.path(diversity_dir, "dada2rsv-beta-braycurtis.txt"), sample_ids),
  Unwunifrac = load_beta_matrix(file.path(diversity_dir, "dada2rsv-beta-unwunifrac.txt"), sample_ids),
  Wunifrac   = load_beta_matrix(file.path(diversity_dir, "dada2rsv-beta-wunifrac.txt"),   sample_ids)
)

alpha_metric_names <- c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")
beta_metric_labels <- c("Beta Diversity ( Braycurtis )",
                        "Beta Diversity ( Unwunifrac )",
                        "Beta Diversity ( Wunifrac )")
beta_matrix_names  <- c("Braycurtis", "Unwunifrac", "Wunifrac")

centroid_dist <- function(level_samples, dist_mat) {
  # Fully vectorized: mean distance from each sample to the OTHER n-1 samples
  # in the same group, using BLAS rowSums on the submatrix. The diagonal of a
  # distance matrix is 0, so rowSums equals the sum over all "other" samples;
  # divide by (n - 1) to get the mean. ~100x faster than the sapply variant.
  n <- length(level_samples)
  if (n < 2) return(rep(NA_real_, n))
  sub <- dist_mat[level_samples, level_samples, drop = FALSE]
  out <- rowSums(sub, na.rm = TRUE) / (n - 1)
  names(out) <- level_samples
  out
}

# -- Phase 2: alpha + beta group means (no reference comparison) -------------

phase2_one_var <- function(cv) {
  vd <- alpha_with_categories %>%
    filter(!is.na(!!sym(cv))) %>%
    select(SEQN, all_of(alpha_metric_names), all_of(cv))
  if (nrow(vd) == 0) return(list(alpha = list(), beta = list()))
  ulev <- unique(vd[[cv]]); ulev <- ulev[!is.na(ulev)]
  if (length(ulev) == 0) return(list(alpha = list(), beta = list()))
  alpha_local <- list(); beta_local <- list()
  for (lv in ulev) {
    ld <- vd %>% filter(!!sym(cv) == lv)
    if (nrow(ld) == 0) next
    gN <- nrow(ld)
    for (m in alpha_metric_names) {
      mv <- ld[[m]]; mv <- mv[!is.na(mv)]
      if (length(mv) > 0) {
        alpha_local[[length(alpha_local) + 1L]] <- data.frame(
          VARNAME = cv, Metric = m, LEVEL = as.character(lv),
          group_mean = mean(mv, na.rm = TRUE), group_N = gN,
          stringsAsFactors = FALSE)
      }
    }
  }
  for (k in seq_along(beta_matrix_names)) {
    metric_name <- beta_metric_labels[k]
    dist_mat <- beta_metrics[[beta_matrix_names[k]]]
    if (is.null(dist_mat)) next
    for (lv in ulev) {
      ls <- vd %>% filter(!!sym(cv) == lv) %>% pull(SEQN)
      ls <- intersect(ls, rownames(dist_mat))
      if (length(ls) < 2) next
      cd <- centroid_dist(ls, dist_mat); cd <- cd[!is.na(cd)]
      if (length(cd) > 0) {
        beta_local[[length(beta_local) + 1L]] <- data.frame(
          VARNAME = cv, Metric = metric_name, LEVEL = as.character(lv),
          group_mean = mean(cd, na.rm = TRUE), group_N = length(ls),
          stringsAsFactors = FALSE)
      }
    }
  }
  list(alpha = alpha_local, beta = beta_local)
}

phase2_results <- mclapply(categorical_vars, phase2_one_var, mc.cores = N_CORES)
alpha_results_p2 <- unlist(lapply(phase2_results, `[[`, "alpha"), recursive = FALSE)
beta_results_p2  <- unlist(lapply(phase2_results, `[[`, "beta"),  recursive = FALSE)

all_results_p2 <- bind_rows(c(alpha_results_p2, beta_results_p2))
write.csv(all_results_p2,
          file.path(OUTPUT_DIR, "full_list_supplementary_alpha_beta_variable_distribution_per_group_per_level.csv"),
          row.names = FALSE)
message("Phase 2 done. ", nrow(all_results_p2), " rows written.")

# -- Phase 3: alpha + beta group means with reference comparisons ------------

ref_levels_df  <- read.csv(ensure_exists(file.path(categorical_dir, "categorical_variables_info.csv")),
                           stringsAsFactors = FALSE)
ref_levels_map <- setNames(ref_levels_df$Reference_Level, ref_levels_df$Variable)

calculate_wilcox_effsize <- function(x, y) {
  test_result <- wilcox.test(x, y, exact = FALSE)
  n1 <- length(x); n2 <- length(y); n <- n1 + n2
  W   <- as.numeric(test_result$statistic)
  E_W <- n1 * (n + 1) / 2
  Var_W <- n1 * n2 * (n + 1) / 12
  Z <- (W - E_W) / sqrt(Var_W)
  r <- max(-1, min(1, Z / sqrt(n)))
  list(p_value = test_result$p.value, effect_size = r, Z = Z)
}

resolve_ref_level <- function(cv) {
  rl <- ref_levels_map[[cv]]
  if (is.null(rl) || is.na(rl)) {
    vd_tmp <- alpha_with_categories %>% filter(!is.na(!!sym(cv))) %>% select(all_of(cv))
    if (nrow(vd_tmp) > 0 && is.factor(vd_tmp[[cv]])) levels(vd_tmp[[cv]])[1] else NA_character_
  } else rl
}

# Per-variable worker: returns list(alpha = data.frame, beta = data.frame)
# already with wilcox_p_FDR filled in (FDR adjustment is local to each
# (variable, metric) pair). Each variable is independent, so we farm out via
# mclapply.
phase3_one_var <- function(cv) {
  ref_level <- resolve_ref_level(cv)
  if (is.na(ref_level)) return(list(alpha = NULL, beta = NULL))
  vd <- alpha_with_categories %>%
    filter(!is.na(!!sym(cv))) %>%
    select(SEQN, all_of(alpha_metric_names), all_of(cv))
  if (nrow(vd) == 0) return(list(alpha = NULL, beta = NULL))
  ulev <- unique(vd[[cv]]); ulev <- ulev[!is.na(ulev)]
  if (length(ulev) == 0) return(list(alpha = NULL, beta = NULL))
  ref_data <- vd %>% filter(!!sym(cv) == ref_level)
  if (nrow(ref_data) == 0) return(list(alpha = NULL, beta = NULL))

  alpha_block <- vector("list", length(alpha_metric_names))
  for (mi in seq_along(alpha_metric_names)) {
    m <- alpha_metric_names[mi]
    ref_vals <- ref_data[[m]]; ref_vals <- ref_vals[!is.na(ref_vals)]
    if (length(ref_vals) == 0) next
    ref_mean <- mean(ref_vals, na.rm = TRUE)
    rows <- list(); pvals <- numeric(0); pidx <- integer(0)
    for (lv in ulev) {
      ld <- vd %>% filter(!!sym(cv) == lv); if (nrow(ld) == 0) next
      mv <- ld[[m]]; mv <- mv[!is.na(mv)]; if (length(mv) == 0) next
      gN <- nrow(ld); gM <- mean(mv, na.rm = TRUE)
      wp <- NA_real_; we <- NA_real_
      if (as.character(lv) != ref_level && length(ref_vals) > 0 && length(mv) > 0) {
        tryCatch({
          wr <- calculate_wilcox_effsize(mv, ref_vals)
          wp <- wr$p_value; we <- wr$effect_size
        }, error = function(e) NULL)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        VARNAME = cv, Metric = m, LEVEL = as.character(lv),
        group_mean = gM, group_N = gN, REF_LEVEL = ref_level,
        group_mean.REF_group_mean = gM - ref_mean,
        wilcox_effsize = we, wilcox_p_val = wp, wilcox_p_FDR = NA_real_,
        stringsAsFactors = FALSE)
      if (!is.na(wp)) { pvals <- c(pvals, wp); pidx <- c(pidx, length(rows)) }
    }
    if (length(pvals) > 0) {
      fdr <- p.adjust(pvals, method = "fdr")
      for (j in seq_along(pidx)) rows[[pidx[j]]]$wilcox_p_FDR <- fdr[j]
    }
    alpha_block[[mi]] <- bind_rows(rows)
  }
  alpha_df <- bind_rows(alpha_block)

  beta_block <- vector("list", length(beta_matrix_names))
  for (k in seq_along(beta_matrix_names)) {
    metric_name <- beta_metric_labels[k]
    dist_mat <- beta_metrics[[beta_matrix_names[k]]]; if (is.null(dist_mat)) next
    ref_samples <- alpha_with_categories %>%
      filter(!is.na(!!sym(cv)), !!sym(cv) == ref_level) %>% pull(SEQN)
    ref_samples <- intersect(ref_samples, rownames(dist_mat))
    if (length(ref_samples) < 2) next
    ref_cd <- centroid_dist(ref_samples, dist_mat); ref_cd <- ref_cd[!is.na(ref_cd)]
    if (length(ref_cd) == 0) next
    ref_mean <- mean(ref_cd, na.rm = TRUE)
    rows <- list(); pvals <- numeric(0); pidx <- integer(0)
    for (lv in ulev) {
      ls <- alpha_with_categories %>%
        filter(!is.na(!!sym(cv)), !!sym(cv) == lv) %>% pull(SEQN)
      ls <- intersect(ls, rownames(dist_mat))
      if (length(ls) < 2) next
      cd <- centroid_dist(ls, dist_mat); cd <- cd[!is.na(cd)]
      if (length(cd) == 0) next
      gM <- mean(cd, na.rm = TRUE)
      wp <- NA_real_; we <- NA_real_
      if (as.character(lv) != ref_level && length(ref_cd) > 0 && length(cd) > 0) {
        tryCatch({
          wr <- calculate_wilcox_effsize(cd, ref_cd)
          wp <- wr$p_value; we <- wr$effect_size
        }, error = function(e) NULL)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        VARNAME = cv, Metric = metric_name, LEVEL = as.character(lv),
        group_mean = gM, group_N = length(ls), REF_LEVEL = ref_level,
        group_mean.REF_group_mean = gM - ref_mean,
        wilcox_effsize = we, wilcox_p_val = wp, wilcox_p_FDR = NA_real_,
        stringsAsFactors = FALSE)
      if (!is.na(wp)) { pvals <- c(pvals, wp); pidx <- c(pidx, length(rows)) }
    }
    if (length(pvals) > 0) {
      fdr <- p.adjust(pvals, method = "fdr")
      for (j in seq_along(pidx)) rows[[pidx[j]]]$wilcox_p_FDR <- fdr[j]
    }
    beta_block[[k]] <- bind_rows(rows)
  }
  beta_df <- bind_rows(beta_block)
  list(alpha = alpha_df, beta = beta_df)
}

phase3_results <- mclapply(categorical_vars, phase3_one_var, mc.cores = N_CORES)
alpha_results_p3 <- bind_rows(lapply(phase3_results, `[[`, "alpha"))
beta_results_p3  <- bind_rows(lapply(phase3_results, `[[`, "beta"))
all_results_p3   <- bind_rows(alpha_results_p3, beta_results_p3)
colnames(all_results_p3)[colnames(all_results_p3) == "group_mean.REF_group_mean"] <- "group_mean-REF_group_mean"
write.csv(all_results_p3,
          file.path(OUTPUT_DIR, "full_list_supplementary_alpha_beta_variable_distribution_per_group_per_level_compared_to_ref.csv"),
          row.names = FALSE)
message("Phase 3 done. ", nrow(all_results_p3), " rows written.")

# =============================================================================
# PHASE 4 — host x OTU regression results supplementary table (+ CLR subset)
# =============================================================================

message("=== PHASE 4: host x OTU regression results table ===")

tidied_with_otu_file <- file.path(module4_inter, "tidied_results_with_otu.rds")
if (!file.exists(tidied_with_otu_file)) {
  tidied_raw_file <- file.path(module4_inter, "tidied_results_raw.rds")
  if (!file.exists(tidied_raw_file)) {
    stop("Neither tidied_results_with_otu.rds nor tidied_results_raw.rds found in ", module4_inter)
  }
  tidied_results <- readRDS(tidied_raw_file)
  use_dependent  <- c("demoWAS", "exWAS")
  use_independent <- c("oradWAS", "pheWAS", "outWAS")
  for (nm in names(tidied_results)) {
    dat <- tidied_results[[nm]]; if (nrow(dat) == 0) next
    if (any(startsWith(nm, use_dependent))) {
      dat <- dat %>% mutate(otu = dependent_var,
                            otu_generic = sub("_relative$", "", dependent_var))
    } else if (any(startsWith(nm, use_independent))) {
      dat <- dat %>% mutate(otu = independent_var,
                            otu_generic = sub("_relative$", "", independent_var))
    } else {
      dat <- dat %>% mutate(otu = NA_character_, otu_generic = NA_character_)
    }
    tidied_results[[nm]] <- dat
  }
} else {
  tidied_results <- readRDS(tidied_with_otu_file)
}

taxonomy_annotations    <- readRDS(ensure_exists(file.path(module4_inter, "taxonomy_annotations.rds")))
phylum_info             <- taxonomy_annotations$phylum_info
genus_mapping           <- taxonomy_annotations$genus_mapping
ubiome_variable_mapping <- readRDS(ensure_exists(file.path(module4_inter, "ubiome_variable_mapping.rds")))
prevalence_filters      <- readRDS(ensure_exists(file.path(module4_inter, "prevalence_filters.rds")))
otu_pass_prevalance_list <- prevalence_filters$otu_pass_prevalance_list
case_filters            <- readRDS(ensure_exists(file.path(module4_inter, "binary_case_filters.rds")))
oradWAS_pass_list       <- case_filters$oradWAS_case_counts_pass_list
outWAS_pass_list        <- case_filters$outWAS_case_counts_pass_list

all_results_list <- list()
for (key in names(tidied_results)) {
  dat <- tidied_results[[key]]
  if (is.null(dat) || nrow(dat) == 0) next
  parts <- strsplit(key, "_")[[1]]
  at <- if (length(parts) >= 2) parts[1] else key
  tr <- if (length(parts) >= 2) paste(parts[-1], collapse = "_") else "unknown"
  all_results_list[[key]] <- dat %>%
    mutate(analysis_type = at, transformation = tr, result_key = key)
}
all_results <- bind_rows(all_results_list) %>%
  mutate(
    host_variable_name = case_when(
      analysis_type %in% c("demoWAS","exWAS") ~ independent_var,
      analysis_type %in% c("oradWAS","pheWAS","outWAS") ~ dependent_var,
      TRUE ~ NA_character_
    ),
    otu_name = otu_generic
  ) %>%
  filter(str_starts(term, independent_var))

taxonomy_combined <- phylum_info %>%
  left_join(genus_mapping, by = "otu") %>%
  mutate(otu_taxonomy_annotation = case_when(
    !is.na(Genus) & Genus != "unclassified" & !is.na(Phylum) ~ paste(Phylum, Genus, sep = "; "),
    !is.na(Phylum) ~ paste(Phylum, "unclassified", sep = "; "),
    TRUE ~ "unclassified"
  )) %>%
  select(otu, otu_taxonomy_annotation, Phylum, Genus)

all_results <- all_results %>% left_join(taxonomy_combined, by = c("otu_name" = "otu"))

all_results_dt <- as.data.table(all_results)
ubiome_dt      <- as.data.table(ubiome_variable_mapping)
all_results_dt <- merge(all_results_dt, ubiome_dt[, .(var_name, var_description)],
                        by.x = "host_variable_name", by.y = "var_name", all.x = TRUE)
all_results <- as_tibble(all_results_dt)

all_results <- all_results %>%
  mutate(
    otu_has_prevalence_0.01 = (otu %in% otu_pass_prevalance_list) |
                              (paste0(otu_generic, "_relative") %in% otu_pass_prevalance_list),
    passes_term_alignment = str_starts(term, independent_var),
    passes_nonzero_se     = !is.na(std.error) & std.error != 0,
    passes_prevalence     = otu_has_prevalence_0.01,
    passes_fdr            = !is.na(p.value.fdr) & p.value.fdr < 0.05,
    passes_q              = !is.na(q.value)     & q.value     < 0.05,
    passes_case_balance   = case_when(
      analysis_type == "outWAS"  ~ dependent_var %in% outWAS_pass_list,
      analysis_type == "oradWAS" ~ dependent_var %in% oradWAS_pass_list,
      TRUE ~ TRUE
    ),
    passes_all_filters = passes_term_alignment & passes_nonzero_se & passes_prevalence &
                         passes_fdr & passes_q & passes_case_balance
  )

summary_stats <- all_results %>%
  group_by(host_variable_name, analysis_type, transformation) %>%
  summarise(
    n_fdr_sig_0.05            = sum(!is.na(p.value.fdr) & p.value.fdr < 0.05, na.rm = TRUE),
    n_fdr_sig_0.05_prev_0.01  = sum(!is.na(p.value.fdr) & p.value.fdr < 0.05 & otu_has_prevalence_0.01, na.rm = TRUE),
    n_pass_all_filters        = sum(passes_all_filters, na.rm = TRUE),
    .groups = "drop"
  )
all_results <- all_results %>%
  left_join(summary_stats, by = c("host_variable_name", "analysis_type", "transformation"))

final_table <- all_results %>%
  select(
    analysis_type, transformation,
    host_variable_name, host_variable_description = var_description,
    otu_name, otu_taxonomy_annotation, Phylum, Genus,
    estimate, std.error, p.value, p.value.fdr, q.value,
    n_fdr_sig_0.05_per_host_variable           = n_fdr_sig_0.05,
    n_fdr_sig_0.05_prev_0.01_per_host_variable = n_fdr_sig_0.05_prev_0.01,
    n_pass_all_filters_per_host_variable       = n_pass_all_filters,
    statistic, term, dependent_var, independent_var, otu, otu_generic,
    available_cycles, effect_scale, interpretation_note, fdr_corrected
  ) %>%
  arrange(analysis_type, transformation, host_variable_name, p.value.fdr)

write.csv(final_table,
          file.path(OUTPUT_DIR, "full_list_supplementary_host_otu_regression_res.csv"),
          row.names = FALSE)

subset_clr <- final_table %>%
  filter(transformation == "clr") %>%
  mutate(analysis_type = case_when(
    analysis_type == "demoWAS" ~ "Demographics -> Microbiome (Linear)",
    analysis_type == "oradWAS" ~ "Microbiome -> Oral Conditions (Logistic)",
    analysis_type == "exWAS"   ~ "Blood/Urine Markers -> Microbiome (Linear)",
    analysis_type == "pheWAS"  ~ "Microbiome -> Measured Phenotypes (Linear)",
    analysis_type == "outWAS"  ~ "Microbiome -> Disease Incidents (Logistic)",
    TRUE ~ analysis_type
  )) %>%
  select(-any_of(c("term","dependent_var","independent_var","otu","otu_generic",
                   "available_cycles","effect_scale","interpretation_note","fdr_corrected")))

write.csv(subset_clr,
          file.path(OUTPUT_DIR, "full_list_supplementary_host_otu_regression_res_subset_clr.csv"),
          row.names = FALSE)

message("Phase 4 done. Full table: ", nrow(final_table), " rows. CLR subset: ", nrow(subset_clr), " rows.")
message("All outputs under: ", OUTPUT_DIR)
