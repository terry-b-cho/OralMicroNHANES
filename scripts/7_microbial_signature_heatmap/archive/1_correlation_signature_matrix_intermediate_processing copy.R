#!/usr/bin/env Rscript

## -----------------------------------------------------------------------------
## Microbial Signature Heatmap – Intermediate Data Builder
## -----------------------------------------------------------------------------
## This script recreates every intermediate object needed by the standalone
## microbial signature heatmap plotting scripts.  It replaces the portion of
## `4_association_phyloseq_analyses.Rmd` that generated:
##   * per-analysis microbial signature matrices (`signature_matrix.rds`)
##   * per-analysis host–host correlation results   (`correlation_matrix.rds`)
##   * CSV exports used for Python or QA workflows (`heatmap_data/*.csv`)
##   * merged (allWAS) signature/correlation outputs and summaries
##   * `heatmap_significant_results.rds` for downstream filtering
##
## All artefacts are written underneath:
##   `results/analyses_results/7_microbial_signature_heatmap_out/intermediate/`
##
## Running this script alone is sufficient to feed both
## `microbial_signature_heatmap_sized.R` and
## `microbial_signature_heatmap_sized_cutree_optimal_k_sigres.R`.
## -----------------------------------------------------------------------------

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(DBI)
  library(RSQLite)
})

resolve_base_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
  if (length(script_path) == 0) {
    return(normalizePath("."))
  }
  normalizePath(file.path(dirname(script_path), "..", ".."))
}

base_path <- resolve_base_path()
message("Base path: ", base_path)

## -----------------------------------------------------------------------------
## Directory layout
## -----------------------------------------------------------------------------
output_root <- file.path(base_path, "results/analyses_results/7_microbial_signature_heatmap_out")
intermediate_root <- file.path(output_root, "intermediate")
microbial_signature_dir <- file.path(intermediate_root, "microbial_signature")
heatmap_data_dir <- file.path(microbial_signature_dir, "heatmap_data")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_root, recursive = TRUE, showWarnings = FALSE)
dir.create(microbial_signature_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_data_dir, recursive = TRUE, showWarnings = FALSE)

clean_dir <- function(path) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

## -----------------------------------------------------------------------------
## Configuration
## -----------------------------------------------------------------------------
analysis_specs <- list(
  demoWAS = list(prefix = "1_demoWAS"),
  oradWAS = list(prefix = "2_oradWAS"),
  exWAS   = list(prefix = "3_exWAS"),
  pheWAS  = list(prefix = "4_pheWAS"),
  outWAS  = list(prefix = "5_outWAS")
)
analysis_types <- names(analysis_specs)
transformation <- "clr"
target_analyses_for_merge <- c("oradWAS", "exWAS", "pheWAS", "outWAS")
use_dependent <- c("demoWAS", "exWAS")
use_independent <- c("oradWAS", "pheWAS", "outWAS")

prevalence_threshold <- 0.01
significance_fdr <- 0.05
min_shared_microbes <- 3
max_variables <- 5000

## -----------------------------------------------------------------------------
## Utility helpers
## -----------------------------------------------------------------------------
safe_save <- function(object, filepath, type = "rds") {
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
  tryCatch({
    if (type == "rds") {
      saveRDS(object, filepath)
    } else if (type == "csv") {
      write.csv(object, filepath, row.names = FALSE)
    } else {
      stop("Unsupported save type: ", type)
    }
    message("  Saved: ", filepath)
    TRUE
  }, error = function(e) {
    warning("  Failed to save ", filepath, ": ", e$message)
    FALSE
  })
}

compute_heatmap_dimensions <- function(n_rows,
                                       min_dim = 300,
                                       base_dim = 80,
                                       per_row = 2.4,
                                       max_dim = 1250) {
  if (is.na(n_rows) || n_rows <= 0) n_rows <- 1
  dim_estimate <- base_dim + n_rows * per_row
  dim_limited <- max(min_dim, min(max_dim, dim_estimate))
  list(width_mm = dim_limited, height_mm = dim_limited)
}

load_tidied_result <- function(analysis, transformation, results_root) {
  prefix <- analysis_specs[[analysis]]$prefix
  dir_path <- file.path(results_root, paste0(prefix, "_out"), paste0("result_", transformation))
  file_path <- file.path(dir_path, paste0(prefix, "_", transformation, "_tidied_complete.rds"))
  if (!file.exists(file_path)) {
    stop("Missing tidied result: ", file_path)
  }
  readRDS(file_path)
}

attach_otu_columns <- function(dataset, analysis) {
  dataset %>%
    mutate(
      dependent_var = as.character(dependent_var),
      independent_var = as.character(independent_var),
      otu = if (analysis %in% use_dependent) dependent_var else independent_var,
      host_var = if (analysis %in% use_dependent) independent_var else dependent_var,
      otu_generic = sub("_relative$", "", otu)
    )
}

apply_all_results_filters <- function(dataset, analysis, prevalence_list) {
  attach_otu_columns(dataset, analysis) %>%
    filter(str_starts(term, independent_var)) %>%
    filter(std.error != 0) %>%
    filter(otu %in% prevalence_list)
}

apply_significant_filters <- function(dataset,
                                      analysis,
                                      prevalence_list,
                                      case_pass_lists) {
  case_vec <- switch(
    analysis,
    oradWAS = case_pass_lists$oradWAS,
    outWAS = case_pass_lists$outWAS,
    character(0)
  )

  attach_otu_columns(dataset, analysis) %>%
    filter(str_starts(term, independent_var)) %>%
    filter(std.error != 0) %>%
    { if (length(case_vec)) filter(., host_var %in% case_vec) else . } %>%
    filter(otu %in% prevalence_list) %>%
    filter(!is.na(p.value.fdr), p.value.fdr < significance_fdr) %>%
    filter(!is.na(q.value), q.value < significance_fdr) %>%
    arrange(p.value.fdr) %>%
    select(
      otu, term, dependent_var, independent_var,
      estimate, statistic, std.error, p.value, p.value.fdr, q.value,
      otu_generic, available_cycles, effect_scale, interpretation_note, fdr_corrected
    )
}

build_signature_matrix <- function(filtered_dataset) {
  filtered_dataset %>%
    select(host_var, microbe_clean = otu_generic, estimate, statistic, std.error) %>%
    filter(!is.na(estimate), !is.na(std.error), std.error > 0)
}

compute_weighted_correlation <- function(sig_dt, var1, var2) {
  data1 <- sig_dt[host_var == var1, .(microbe_clean, est1 = estimate, se1 = std.error)]
  data2 <- sig_dt[host_var == var2, .(microbe_clean, est2 = estimate, se2 = std.error)]

  if (nrow(data1) == 0 || nrow(data2) == 0) {
    return(list(
      correlation = NA_real_,
      p_value = NA_real_,
      n_shared = 0,
      n_effective = NA_real_
    ))
  }

  data1 <- unique(data1, by = "microbe_clean")
  data2 <- unique(data2, by = "microbe_clean")

  setkey(data1, microbe_clean)
  setkey(data2, microbe_clean)

  shared <- merge(data1, data2, by = "microbe_clean", all = FALSE)
  shared <- shared[!duplicated(microbe_clean)]

  n_shared <- nrow(shared)
  if (n_shared < min_shared_microbes) {
    return(list(
      correlation = NA_real_,
      p_value = NA_real_,
      n_shared = n_shared,
      n_effective = NA_real_
    ))
  }

  shared[, combined_weight := 1 / (se1^2 + se2^2)]
  shared <- shared[is.finite(combined_weight) & combined_weight > 0]

  if (nrow(shared) < min_shared_microbes) {
    return(list(
      correlation = NA_real_,
      p_value = NA_real_,
      n_shared = nrow(shared),
      n_effective = NA_real_
    ))
  }

  x <- shared$est1
  y <- shared$est2
  w <- shared$combined_weight

  sum_w <- sum(w)
  x_mean <- sum(w * x) / sum_w
  y_mean <- sum(w * y) / sum_w

  numerator <- sum(w * (x - x_mean) * (y - y_mean))
  x_var <- sum(w * (x - x_mean)^2)
  y_var <- sum(w * (y - y_mean)^2)

  denominator <- sqrt(x_var * y_var)
  if (denominator == 0 || is.na(denominator)) {
    return(list(
      correlation = NA_real_,
      p_value = NA_real_,
      n_shared = nrow(shared),
      n_effective = NA_real_
    ))
  }

  correlation <- numerator / denominator
  n_effective <- sum(w > 0)

  if (n_effective > 3 && abs(correlation) < 0.999) {
    t_stat <- correlation * sqrt((n_effective - 2) / (1 - correlation^2))
    p_value <- 2 * pt(abs(t_stat), df = n_effective - 2, lower.tail = FALSE)
  } else {
    p_value <- ifelse(abs(correlation) > 0.999, 0, 1)
  }

  list(
    correlation = correlation,
    p_value = p_value,
    n_shared = nrow(shared),
    n_effective = n_effective
  )
}

compute_sparse_correlation_matrix <- function(sig_matrix,
                                              min_shared_microbes = 3,
                                              max_vars = 5000) {
  sig_dt <- as.data.table(sig_matrix)
  setkey(sig_dt, host_var, microbe_clean)

  host_vars <- unique(sig_dt$host_var)
  n_vars <- length(host_vars)
  message("    Host variables available: ", n_vars)

  if (n_vars == 0) {
    return(list(
      correlation_matrix = matrix(numeric(0)),
      pvalue_matrix = matrix(numeric(0)),
      fdr_matrix = matrix(numeric(0)),
      n_shared_matrix = matrix(numeric(0)),
      n_effective_matrix = matrix(numeric(0)),
      n_variables = 0,
      n_valid_correlations = 0,
      n_significant_correlations_005 = 0,
      n_significant_correlations_001 = 0,
      n_shared_distribution = numeric(0),
      n_effective_distribution = numeric(0),
      correlation_distribution = numeric(0),
      correlation_results_table = data.table()
    ))
  }

  if (n_vars > max_vars) {
    message("    Limiting to top ", max_vars, " variables by microbe coverage")
    var_counts <- sig_dt[, .N, by = host_var][order(-N)]
    keep_vars <- var_counts[1:max_vars, host_var]
    sig_dt <- sig_dt[host_var %in% keep_vars]
    host_vars <- keep_vars
    n_vars <- length(host_vars)
  }

  results <- data.table()
  total_pairs <- n_vars * (n_vars - 1) / 2
  progress_interval <- max(1, floor(total_pairs / 20))
  pair_counter <- 0L

  for (i in seq_len(n_vars - 1)) {
    for (j in (i + 1):n_vars) {
      res <- compute_weighted_correlation(sig_dt, host_vars[i], host_vars[j])
      if (!is.na(res$correlation) && res$n_shared >= min_shared_microbes) {
        results <- rbind(
          results,
          data.table(
            var1 = host_vars[i],
            var2 = host_vars[j],
            correlation = res$correlation,
            p_value = res$p_value,
            n_shared = res$n_shared,
            n_effective = res$n_effective
          ),
          fill = TRUE
        )
      }
      pair_counter <- pair_counter + 1
      if (pair_counter %% progress_interval == 0) {
        pct <- round(100 * pair_counter / total_pairs, 1)
        message("      Progress: ", pct, "% (", pair_counter, "/", total_pairs, ")")
      }
    }
  }

  results[, fdr_pvalue := p.adjust(p_value, method = "fdr")]

  cor_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  pvalue_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  fdr_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  n_shared_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))
  n_eff_matrix <- matrix(NA_real_, nrow = n_vars, ncol = n_vars, dimnames = list(host_vars, host_vars))

  diag(cor_matrix) <- 1
  diag(pvalue_matrix) <- 0
  diag(fdr_matrix) <- 0
  host_counts <- sig_dt[, .N, by = host_var][match(host_vars, host_var)]$N
  diag(n_shared_matrix) <- host_counts
  diag(n_eff_matrix) <- host_counts

  if (nrow(results)) {
    for (idx in seq_len(nrow(results))) {
      i <- match(results$var1[idx], host_vars)
      j <- match(results$var2[idx], host_vars)
      cor_matrix[i, j] <- cor_matrix[j, i] <- results$correlation[idx]
      pvalue_matrix[i, j] <- pvalue_matrix[j, i] <- results$p_value[idx]
      fdr_matrix[i, j] <- fdr_matrix[j, i] <- results$fdr_pvalue[idx]
      n_shared_matrix[i, j] <- n_shared_matrix[j, i] <- results$n_shared[idx]
      n_eff_matrix[i, j] <- n_eff_matrix[j, i] <- results$n_effective[idx]
    }
  }

  list(
    correlation_matrix = cor_matrix,
    pvalue_matrix = pvalue_matrix,
    fdr_matrix = fdr_matrix,
    n_shared_matrix = n_shared_matrix,
    n_effective_matrix = n_eff_matrix,
    n_variables = n_vars,
    n_valid_correlations = nrow(results),
    n_significant_correlations_005 = sum(results$fdr_pvalue <= 0.05, na.rm = TRUE),
    n_significant_correlations_001 = sum(results$fdr_pvalue <= 0.01, na.rm = TRUE),
    n_shared_distribution = if (nrow(results)) results$n_shared else numeric(0),
    n_effective_distribution = if (nrow(results)) results$n_effective else numeric(0),
    correlation_distribution = if (nrow(results)) results$correlation else numeric(0),
    correlation_results_table = results
  )
}

save_heatmap_data <- function(correlation_result, dataset_key, variable_mapping) {
  dir.create(heatmap_data_dir, showWarnings = FALSE, recursive = TRUE)

  cor_matrix <- correlation_result$correlation_matrix
  fdr_matrix <- correlation_result$fdr_matrix
  n_shared_matrix <- correlation_result$n_shared_matrix
  n_eff_matrix <- correlation_result$n_effective_matrix

  vars <- rownames(cor_matrix)
  descriptions <- variable_mapping$var_description[match(vars, variable_mapping$var_name)]
  descriptions[is.na(descriptions)] <- vars[is.na(descriptions)]

  var_map <- data.frame(
    var_name = vars,
    var_description = descriptions,
    stringsAsFactors = FALSE
  )

  write.csv(
    as.data.frame(cor_matrix) %>%
      mutate(var_name = vars),
    file.path(heatmap_data_dir, paste0(dataset_key, "_correlation_matrix.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(fdr_matrix) %>%
      mutate(var_name = vars),
    file.path(heatmap_data_dir, paste0(dataset_key, "_fdr_matrix.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(n_eff_matrix) %>%
      mutate(var_name = vars),
    file.path(heatmap_data_dir, paste0(dataset_key, "_n_effective_matrix.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(n_shared_matrix) %>%
      mutate(var_name = vars),
    file.path(heatmap_data_dir, paste0(dataset_key, "_n_shared_matrix.csv")),
    row.names = FALSE
  )

  write.csv(
    var_map,
    file.path(heatmap_data_dir, paste0(dataset_key, "_variable_mapping.csv")),
    row.names = FALSE
  )

  summary_stats <- data.frame(
    Dataset = dataset_key,
    N_Variables = correlation_result$n_variables,
    N_Valid_Correlations = correlation_result$n_valid_correlations,
    N_Significant_005 = correlation_result$n_significant_correlations_005,
    N_Significant_001 = correlation_result$n_significant_correlations_001,
    Max_Abs_Correlation = if (length(correlation_result$correlation_distribution)) {
      round(max(abs(correlation_result$correlation_distribution), na.rm = TRUE), 3)
    } else NA_real_,
    Median_N_Shared = if (length(correlation_result$n_shared_distribution)) {
      round(median(correlation_result$n_shared_distribution, na.rm = TRUE), 1)
    } else NA_real_,
    Median_N_Effective = if (length(correlation_result$n_effective_distribution)) {
      round(median(correlation_result$n_effective_distribution, na.rm = TRUE), 1)
    } else NA_real_,
    stringsAsFactors = FALSE
  )

  write.csv(
    summary_stats,
    file.path(heatmap_data_dir, paste0(dataset_key, "_summary_stats.csv")),
    row.names = FALSE
  )
}

## -----------------------------------------------------------------------------
## Data sources
## -----------------------------------------------------------------------------
aggregated_results_root <- file.path(base_path, "results")
config_dir_path <- file.path(base_path, "configs")
db_path <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")

if (!file.exists(db_path)) {
  stop("SQLite database not found at ", db_path)
}

message("Connecting to SQLite database...")
con <- dbConnect(SQLite(), db_path)
on.exit({
  try(dbDisconnect(con), silent = TRUE)
}, add = TRUE)

message("Loading variable metadata...")
variable_description <- tbl(con, "variable_names_epcf") %>%
  collect() %>%
  transmute(
    var_name = Variable.Name,
    var_description = Variable.Description
  )

lookup_description <- function(var_names) {
  desc <- variable_description$var_description[match(var_names, variable_description$var_name)]
  desc[is.na(desc) | desc == ""] <- var_names[is.na(desc) | desc == ""]
  tibble(var_name = var_names, var_description = desc)
}

message("Computing OTU prevalence (threshold ", prevalence_threshold * 100, "%)...")
rsv_genus_relative <- bind_rows(
  tbl(con, "DADA2RSV_GENUS_RELATIVE_F") %>% collect(),
  tbl(con, "DADA2RSV_GENUS_RELATIVE_G") %>% collect()
)

total_samples <- nrow(rsv_genus_relative)
otu_non_zero <- rsv_genus_relative %>%
  summarise(across(where(is.numeric), ~ sum(. > 0), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "otu", values_to = "non_zero_count")

otu_pass_prevalence_list <- otu_non_zero %>%
  filter(non_zero_count > total_samples * prevalence_threshold,
         otu != "SEQN") %>%
  pull(otu)

message("Loaded ", length(otu_pass_prevalence_list), " OTUs passing prevalence filter")

## -----------------------------------------------------------------------------
## Case-count filters for binary outcomes
## -----------------------------------------------------------------------------
message("Computing binary case-count filters...")
genus_f <- tbl(con, "DADA2RSV_GENUS_RELATIVE_F") %>%
  select(SEQN) %>%
  collect() %>%
  mutate(SEQN = as.character(SEQN))
genus_g <- tbl(con, "DADA2RSV_GENUS_RELATIVE_G") %>%
  select(SEQN) %>%
  collect() %>%
  mutate(SEQN = as.character(SEQN))

load_vars <- function(path) {
  if (!file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  lines[nzchar(lines) & !grepl("^#", lines)]
}

var_files <- list(
  demoWAS = file.path(config_dir_path, "1_demoWAS_vars.txt"),
  oradWAS = file.path(config_dir_path, "2_oradWAS_vars.txt"),
  exWAS   = file.path(config_dir_path, "3_exWAS_vars.txt"),
  pheWAS  = file.path(config_dir_path, "4_pheWAS_vars.txt"),
  outWAS  = file.path(config_dir_path, "5_outWAS_vars.txt"),
  zimWAS  = file.path(config_dir_path, "6_zimWAS_vars.txt")
)

var_source_mapping <- map(var_files, load_vars)

get_binary_case_counts <- function(var_list, base_name, con, genus_f, genus_g) {
  if (!length(var_list)) return(tibble(var_name = character(), cases_count = integer()))

  table_f <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_f" else paste0(base_name, "_F")
  table_g <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_g" else paste0(base_name, "_G")
  db_tables <- dbListTables(con)

  if (!(table_f %in% db_tables && table_g %in% db_tables)) {
    return(tibble(var_name = character(), cases_count = integer()))
  }

  df_f <- tbl(con, table_f) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  df_g <- tbl(con, table_g) %>% collect() %>% mutate(SEQN = as.character(SEQN))

  df_all <- bind_rows(
    df_f %>% filter(SEQN %in% genus_f$SEQN),
    df_g %>% filter(SEQN %in% genus_g$SEQN)
  )

  map_dfr(var_list, function(v) {
    tibble(
      var_name = v,
      cases_count = if (v %in% names(df_all)) {
        sum(df_all[[v]] %in% c(1, "1", TRUE), na.rm = TRUE)
      } else {
        NA_integer_
      }
    )
  })
}

orad_case_counts <- get_binary_case_counts(var_source_mapping$oradWAS, "OralDisease", con, genus_f, genus_g)
out_case_counts <- get_binary_case_counts(var_source_mapping$outWAS, "d_outcome_mcq", con, genus_f, genus_g)

total_genus_samples <- length(unique(c(genus_f$SEQN, genus_g$SEQN)))
case_cutoff <- total_genus_samples * 0.005

case_pass_lists <- list(
  oradWAS = orad_case_counts %>% filter(cases_count > case_cutoff) %>% pull(var_name),
  outWAS = out_case_counts %>% filter(cases_count > case_cutoff) %>% pull(var_name)
)

## -----------------------------------------------------------------------------
## Variable mapping table (for CSV exports)
## -----------------------------------------------------------------------------
variable_mapping <- lookup_description(unique(variable_description$var_name))

## -----------------------------------------------------------------------------
## Per-analysis processing
## -----------------------------------------------------------------------------
tidied_results <- list()
significant_results <- list()
signature_matrices <- list()
summary_rows <- list()

message("\n=== Processing individual analyses ===")

for (analysis in analysis_types) {
  dataset_name <- paste0(analysis, "_", transformation, "_all_results")
  dataset_dir <- file.path(microbial_signature_dir, dataset_name)
  clean_dir(dataset_dir)

  message("\n--- Dataset: ", dataset_name, " ---")
  tidied <- load_tidied_result(analysis, transformation, aggregated_results_root)
  tidied_results[[dataset_name]] <- tidied

  all_filtered <- apply_all_results_filters(tidied, analysis, otu_pass_prevalence_list)
  if (!nrow(all_filtered)) {
    message("  No rows after baseline filtering; skipping correlation")
    significant_results[[paste0(analysis, "_", transformation, "_sig_res")]] <- tibble()
    next
  }

  sig_filtered <- apply_significant_filters(
    tidied,
    analysis,
    otu_pass_prevalence_list,
    case_pass_lists
  )
  significant_results[[paste0(analysis, "_", transformation, "_sig_res")]] <- sig_filtered

  sig_matrix <- build_signature_matrix(all_filtered)
  signature_matrices[[dataset_name]] <- sig_matrix
  safe_save(sig_matrix, file.path(dataset_dir, "signature_matrix.rds"))

  cor_result <- compute_sparse_correlation_matrix(
    sig_matrix,
    min_shared_microbes = min_shared_microbes,
    max_vars = max_variables
  )
  safe_save(cor_result, file.path(dataset_dir, "correlation_matrix.rds"))
  save_heatmap_data(cor_result, dataset_name, variable_mapping)

  dims <- compute_heatmap_dimensions(cor_result$n_variables)
  message("  Heatmap dimension estimate: ", round(dims$width_mm, 1), "mm")

  summary_rows[[dataset_name]] <- tibble(
    Dataset = dataset_name,
    Analysis_Type = analysis,
    Transformation = transformation,
    N_Signature_Rows = nrow(sig_matrix),
    N_Host_Variables = cor_result$n_variables,
    N_Valid_Correlations = cor_result$n_valid_correlations,
    N_Significant_FDR_005 = cor_result$n_significant_correlations_005,
    N_Significant_FDR_001 = cor_result$n_significant_correlations_001,
    Output_Dir = dataset_dir
  )
}

## -----------------------------------------------------------------------------
## Save tidied results bundle (optional but useful for QA)
## -----------------------------------------------------------------------------
safe_save(tidied_results, file.path(intermediate_root, "tidied_results_clr.rds"))

## -----------------------------------------------------------------------------
## Save significant results bundle for plotting scripts
## -----------------------------------------------------------------------------
significant_results_path <- file.path(intermediate_root, "heatmap_significant_results.rds")
safe_save(significant_results, significant_results_path)

## -----------------------------------------------------------------------------
## Build merged (allWAS) dataset
## -----------------------------------------------------------------------------
message("\n=== Building merged CLR dataset (allWAS) ===")

create_merged_dataset <- function(signature_matrices,
                                  transformation,
                                  target_analyses,
                                  var_source_mapping,
                                  variable_mapping) {
  merged_dir <- file.path(microbial_signature_dir, paste0("merged_", transformation, "_data"))
  clean_dir(merged_dir)

  combined_sig <- data.table()
  analysis_contributions <- tibble()

  for (analysis in target_analyses) {
    dataset_name <- paste0(analysis, "_", transformation, "_all_results")
    sig_matrix <- signature_matrices[[dataset_name]]
    if (is.null(sig_matrix)) {
      message("  WARNING: missing signature matrix for ", dataset_name)
      next
    }

    source_vars <- var_source_mapping[[analysis]]
    if (!length(source_vars)) {
      message("  WARNING: no configured variables for ", analysis)
      next
    }

    sig_dt <- as.data.table(sig_matrix)
    filtered_sig <- sig_dt[host_var %in% source_vars]

    analysis_contributions <- bind_rows(
      analysis_contributions,
      tibble(
        Analysis = analysis,
        Original_Rows = nrow(sig_dt),
        Original_Variables = dplyr::n_distinct(sig_dt$host_var),
        Filtered_Rows = nrow(filtered_sig),
        Filtered_Variables = dplyr::n_distinct(filtered_sig$host_var)
      )
    )

    combined_sig <- rbind(combined_sig, filtered_sig, fill = TRUE)
  }

  if (!nrow(combined_sig)) {
    stop("Merged dataset is empty.")
  }

  unique_vars <- unique(combined_sig$host_var)
  var_sources <- map_chr(unique_vars, function(v) {
    found <- names(var_source_mapping)[map_lgl(var_source_mapping, ~ v %in% .x)]
    if (!length(found)) "Unknown" else if (found[1] == "demoWAS") "oradWAS" else found[1]
  })

  unwanted_sources <- c("demoWAS", "zimWAS", "Unknown")
  unwanted_vars <- unique_vars[var_sources %in% unwanted_sources]
  if (length(unwanted_vars)) {
    message("  Removing ", length(unwanted_vars), " unwanted host variables from merged dataset")
    combined_sig <- combined_sig[!host_var %in% unwanted_vars]
    unique_vars <- unique(combined_sig$host_var)
    var_sources <- map_chr(unique_vars, function(v) {
      found <- names(var_source_mapping)[map_lgl(var_source_mapping, ~ v %in% .x)]
      if (!length(found)) "Unknown" else if (found[1] == "demoWAS") "oradWAS" else found[1]
    })
  }

  composition_summary <- tibble(
    Transformation = transformation,
    Total_Rows = nrow(combined_sig),
    Unique_Host_Variables = length(unique_vars),
    Unique_Microbes = dplyr::n_distinct(combined_sig$microbe_clean),
    oradWAS_Variables = sum(var_sources == "oradWAS"),
    exWAS_Variables = sum(var_sources == "exWAS"),
    pheWAS_Variables = sum(var_sources == "pheWAS"),
    outWAS_Variables = sum(var_sources == "outWAS"),
    Unknown_Variables = sum(var_sources == "Unknown")
  )

  safe_save(analysis_contributions, file.path(merged_dir, "analysis_contributions.csv"), "csv")
  safe_save(composition_summary, file.path(merged_dir, "composition_summary.csv"), "csv")
  safe_save(combined_sig, file.path(merged_dir, "signature_matrix.rds"))

  cor_result <- compute_sparse_correlation_matrix(
    combined_sig,
    min_shared_microbes = min_shared_microbes,
    max_vars = max_variables
  )
  safe_save(cor_result, file.path(merged_dir, "correlation_matrix.rds"))
  save_heatmap_data(cor_result, "allWAS_clr_all_results", variable_mapping)

  detailed_summary <- data.frame(
    Dataset = "allWAS_clr_all_results",
    Analysis_Type = "merged",
    Transformation = transformation,
    Total_Variables = cor_result$n_variables,
    Total_Signature_Rows = nrow(combined_sig),
    Unique_Microbes = dplyr::n_distinct(combined_sig$microbe_clean),
    N_Valid_Correlations = cor_result$n_valid_correlations,
    N_Significant_005 = cor_result$n_significant_correlations_005,
    N_Significant_001 = cor_result$n_significant_correlations_001,
    Max_Abs_Correlation = if (length(cor_result$correlation_distribution)) {
      round(max(abs(cor_result$correlation_distribution), na.rm = TRUE), 3)
    } else NA_real_,
    Min_Correlation = if (length(cor_result$correlation_distribution)) {
      round(min(cor_result$correlation_distribution, na.rm = TRUE), 3)
    } else NA_real_,
    Max_Correlation = if (length(cor_result$correlation_distribution)) {
      round(max(cor_result$correlation_distribution, na.rm = TRUE), 3)
    } else NA_real_,
    Median_N_Shared = if (length(cor_result$n_shared_distribution)) {
      round(median(cor_result$n_shared_distribution, na.rm = TRUE), 1)
    } else NA_real_,
    Median_N_Effective = if (length(cor_result$n_effective_distribution)) {
      round(median(cor_result$n_effective_distribution, na.rm = TRUE), 1)
    } else NA_real_,
    stringsAsFactors = FALSE
  )

  safe_save(detailed_summary, file.path(merged_dir, "detailed_correlation_summary.csv"), "csv")

  list(
    dataset_dir = merged_dir,
    signature_matrix = combined_sig,
    correlation_result = cor_result,
    summary_row = tibble(
      Dataset = "allWAS_clr_all_results",
      Analysis_Type = "allWAS",
      Transformation = transformation,
      N_Signature_Rows = nrow(combined_sig),
      N_Host_Variables = cor_result$n_variables,
      N_Valid_Correlations = cor_result$n_valid_correlations,
      N_Significant_FDR_005 = cor_result$n_significant_correlations_005,
      N_Significant_FDR_001 = cor_result$n_significant_correlations_001,
      Output_Dir = merged_dir
    )
  )
}

merged_info <- create_merged_dataset(
  signature_matrices = signature_matrices,
  transformation = transformation,
  target_analyses = target_analyses_for_merge,
  var_source_mapping = var_source_mapping,
  variable_mapping = variable_mapping
)

summary_rows[["allWAS_clr_all_results"]] <- merged_info$summary_row

## -----------------------------------------------------------------------------
## Verification summary
## -----------------------------------------------------------------------------
dataset_summary <- bind_rows(summary_rows)
safe_save(dataset_summary, file.path(microbial_signature_dir, "dataset_summary.csv"), "csv")

message("\n=== Verification summary ===")
if (nrow(dataset_summary)) {
  print(dataset_summary)
}

verification <- dataset_summary %>%
  mutate(
    Signature_Matrix = file.exists(file.path(Output_Dir, "signature_matrix.rds")),
    Correlation_Matrix = file.exists(file.path(Output_Dir, "correlation_matrix.rds"))
  ) %>%
  select(Dataset, Signature_Matrix, Correlation_Matrix)

safe_save(verification, file.path(microbial_signature_dir, "complete_verification_results.csv"), "csv")

message("\nIntermediate processing complete.")
message("Outputs written to: ", intermediate_root)

