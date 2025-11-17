#!/usr/bin/env Rscript

## -----------------------------------------------------------------------------
## Standalone Terry Plot Generator
## -----------------------------------------------------------------------------
## This script reproduces every "Terry plot" (host, OTU, and *_indices variants)
## generated in `scripts/4_association_phyloseq_analyses/4_association_phyloseq_analyses.Rmd`.
## All outputs are redirected to `results/analyses_results/4_otu_host_plot_out`.
## The implementation preserves every plotting parameter from the original Rmd
## so the resulting PDFs are carbon copies of the originals.
## -----------------------------------------------------------------------------

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(forcats)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(egg)
  library(gridExtra)
  library(grid)
  library(extrafont)
  library(DBI)
  library(RSQLite)
  library(phyloseq)
  library(Matrix)
})

message("### Terry Plot Standalone Script ###")

## -----------------------------------------------------------------------------
## Helper: resolve repository base path from script location
## -----------------------------------------------------------------------------
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
## Directory setup
## -----------------------------------------------------------------------------
output_base_path <- file.path(base_path, "results/analyses_results/4_otu_host_plot_out")
viz_out_path <- file.path(output_base_path, "figures_out")
intermediate_files_path <- file.path(output_base_path, "intermediate")

dir.create(viz_out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_files_path, recursive = TRUE, showWarnings = FALSE)

## Ensure subdirectories for plotting outputs exist
dir.create(file.path(viz_out_path, "none"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(viz_out_path, "hellinger"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(viz_out_path, "clr"), recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(viz_out_path, "indices", "none"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(viz_out_path, "indices", "hellinger"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(viz_out_path, "indices", "clr"), recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## Paths to data sources
## -----------------------------------------------------------------------------
aggregated_association_res_path <- file.path(base_path, "results")
phyloseq_obj_files_path <- file.path(base_path, "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate")
gold_db_ubiome_genus_mapping_path <- file.path(base_path, "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv")
config_dir_path <- file.path(base_path, "configs")
db_path <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")

## -----------------------------------------------------------------------------
## Fonts (ArialMT required for publication specs)
## -----------------------------------------------------------------------------
message("Loading fonts...")
try({
  if (!"ArialMT" %in% extrafont::fonts()) {
    suppressWarnings(extrafont::font_import(path = file.path(base_path, "data/fonts"),
                                           pattern = "arial.*\\.ttf$",
                                           prompt = FALSE))
    extrafont::loadfonts(device = "pdf", quiet = TRUE)
  } else {
    extrafont::loadfonts(device = "pdf", quiet = TRUE)
  }
}, silent = TRUE)

## -----------------------------------------------------------------------------
## Utility functions
## -----------------------------------------------------------------------------
save_pdf_figure <- function(file_name, plot_obj, width_mm = 180, height_mm = 170) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  target_path <- file.path(viz_out_path, file_name)
  message("Saving PDF: ", target_path)
  ggsave(
    filename = target_path,
    plot = plot_obj,
    device = "pdf",
    width = width_mm,
    height = height_mm,
    units = "mm",
    family = "ArialMT",
    colormodel = "rgb",
    pointsize = 5,
    useDingbats = FALSE
  )
}

## -----------------------------------------------------------------------------
## Color specifications
## -----------------------------------------------------------------------------
phylum_colors <- c(
  "Firmicutes"           = "#F38400",
  "Bacteroidetes"        = "#0067A5",
  "Actinobacteria"       = "#8DB600",
  "Proteobacteria"       = "#E68FAC",
  "Fusobacteria"         = "#BE0032",
  "Spirochaetae"         = "#F3C300",
  "Cyanobacteria"        = "#875692",
  "Acidobacteria"        = "#F6A600",
  "Candidate division SR1" = "#2B3D26",
  "Planctomycetes"       = "#332288",
  "Saccharibacteria"     = "#B3446C",
  "Synergistetes"        = "#A1CAF1",
  "Tenericutes"          = "#654522",
  "Verrucomicrobia"      = "#C2B280",
  "unclassified"         = "#DDDDDD",
  "TM7"                  = "#000000"
)

host_variable_colors <- c(
  "DUM" = "#FFFFFF",
  "RID" = "#E69F00", "AGE" = "#56B4E9", "BOR" = "#009E73", "ETH" = "#F0E442",
  "RIA" = "#CC6677", "IND" = "#D55E00", "EDU" = "#DDA0DD", "RSV" = "#999999",
  "LBX" = "#332288", "URX" = "#44AA99", "SMQ" = "#88CCEE", "LBD" = "#117733",
  "DS1" = "#DDDDDD", "DSQ" = "#0072B2", "HOQ" = "#000000", "DR1" = "#AD7700",
  "DS2" = "#882255", "DR2" = "#661100", "PAQ" = "#6699CC", "GUM" = "#AA4499",
  "ORA" = "#DDCC77", "TOO" = "#117733", "DEN" = "#88CCEE", "CAN" = "#332288",
  "HEA" = "#E69F00", "EMP" = "#56B4E9", "STR" = "#009E73", "ANG" = "#F0E442",
  "CHD" = "#CC6677", "CVD" = "#D55E00", "DIA" = "#DDA0DD", "AST" = "#999999",
  "BRO" = "#AD7700", "BMX" = "#44AA99", "BPX" = "#88CCEE",
  "CANCER" = "#FF6B6B", "DIABETES" = "#45B7D1", "HEART" = "#4ECDC4",
  "ASTHMA" = "#96CEB4", "BRONCHITIS" = "#FFEAA7", "EMPHYSEMA" = "#DDA0DD",
  "ANGINA" = "#74B9FF", "STROKE" = "#A0E7E5"
)

grafify_all_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#332288", "#88ccee", "#44aa99", "#117733", "#999933",
  "#ddcc77", "#cc6677", "#882255", "#aa4499", "#77aadd", "#99ddff",
  "#44bb99", "#bbcc33", "#aaaa00", "#eedd88", "#ee8866", "#ffaabb",
  "#bbccee", "#cceeff", "#ccddaa", "#eeeebb", "#ffcccc", "#dddddd",
  "#555555", "#bbbbbb", "#999999", "#ffffff", "#ddaa33", "#bb5566",
  "#004488", "#000000", "#0077bb", "#33bbee", "#009988", "#ee7733",
  "#cc3311", "#ee3377", "#222255", "#225555", "#225522", "#666633",
  "#663333", "#4477aa", "#66ccee", "#228833", "#ccbb44", "#ee6677",
  "#aa3377"
)

## -----------------------------------------------------------------------------
## Database connection
## -----------------------------------------------------------------------------
message("Connecting to SQLite database...")
if (!file.exists(db_path)) {
  stop("SQLite database not found at ", db_path)
}
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)
on.exit({
  try(dbDisconnect(con), silent = TRUE)
}, add = TRUE)

## -----------------------------------------------------------------------------
## Load GOLD mapping and phyloseq objects
## -----------------------------------------------------------------------------
message("Loading GOLD genus mapping...")
ubiome_genus_mapping_complete <- readr::read_csv(gold_db_ubiome_genus_mapping_path, show_col_types = FALSE) %>%
  data.table::as.data.table()

message("Loading phyloseq objects...")
phyloseq_files <- c(
  "ubiome_counts" = "ubiome_counts.rds",
  "ubiome_relative" = "ubiome_relative.rds",
  "ubiome_relative_none" = "ubiome_relative_none.rds",
  "ubiome_relative_clr" = "ubiome_relative_clr.rds",
  "ubiome_relative_hellinger" = "ubiome_relative_hellinger.rds",
  "ubiome_relative_lognorm" = "ubiome_relative_lognorm.rds"
)

phyloseq_objects <- list()
for (name in names(phyloseq_files)) {
  file_path <- file.path(phyloseq_obj_files_path, phyloseq_files[[name]])
  if (!file.exists(file_path)) {
    stop("Missing phyloseq object: ", file_path)
  }
  phyloseq_objects[[name]] <- readRDS(file_path)
  message(sprintf("  Loaded %s (Samples: %d, Taxa: %d)",
                  name,
                  nsamples(phyloseq_objects[[name]]),
                  ntaxa(phyloseq_objects[[name]])))
}

ubiome_counts <- phyloseq_objects[["ubiome_counts"]]
ubiome_relative <- phyloseq_objects[["ubiome_relative"]]
ubiome_relative_none <- phyloseq_objects[["ubiome_relative_none"]]
ubiome_relative_clr <- phyloseq_objects[["ubiome_relative_clr"]]
ubiome_relative_hellinger <- phyloseq_objects[["ubiome_relative_hellinger"]]
ubiome_relative_lognorm <- phyloseq_objects[["ubiome_relative_lognorm"]]

## -----------------------------------------------------------------------------
## Helper functions for variable color mapping and annotations
## -----------------------------------------------------------------------------
extract_variable_prefix <- function(variable_names, target_prefixes = NULL, default_length = 3) {
  if (is.null(target_prefixes)) {
    return(substr(variable_names, 1, default_length))
  }
  prefixes <- character(length(variable_names))
  for (i in seq_along(variable_names)) {
    var_name <- variable_names[i]
    matched <- FALSE
    for (prefix in target_prefixes[order(nchar(target_prefixes), decreasing = TRUE)]) {
      if (startsWith(toupper(var_name), toupper(prefix))) {
        prefixes[i] <- prefix
        matched <- TRUE
        break
      }
    }
    if (!matched) {
      prefixes[i] <- substr(var_name, 1, default_length)
    }
  }
  prefixes
}

get_variable_colors <- function(variable_names, base_colors = host_variable_colors) {
  unique_vars <- unique(variable_names)
  color_map <- base_colors
  unmapped_vars <- setdiff(unique_vars, names(color_map))
  if (length(unmapped_vars) > 0) {
    n_needed <- length(unmapped_vars)
    assigned_colors <- grafify_all_colors[seq_len(n_needed)]
    names(assigned_colors) <- unmapped_vars
    color_map <- c(color_map, assigned_colors)
  }
  color_map
}

calc_dot_size <- function(p) {
  case_when(
    p < 1e-6  ~ 3.5,
    p < 1e-5  ~ 3.0,
    p < 1e-4  ~ 2.5,
    p < 1e-3  ~ 2.0,
    p < 1e-2  ~ 1.5,
    TRUE      ~ 1.0
  )
}

## -----------------------------------------------------------------------------
## Load association result tables
## -----------------------------------------------------------------------------
message("Loading association result tables...")

load_result <- function(path) {
  if (file.exists(path)) {
    readRDS(path)
  } else {
    tibble()
  }
}

analysis_specs <- list(
  demoWAS = list(prefix = "1_demoWAS", transformations = c("none", "clr", "hellinger")),
  oradWAS = list(prefix = "2_oradWAS", transformations = c("none", "clr", "hellinger")),
  exWAS   = list(prefix = "3_exWAS",   transformations = c("none", "clr", "hellinger")),
  pheWAS  = list(prefix = "4_pheWAS",  transformations = c("none", "clr", "hellinger")),
  outWAS  = list(prefix = "5_outWAS",  transformations = c("none", "clr", "hellinger")),
  zimWAS  = list(prefix = "6_zimWAS",  transformations = c("none", "clr", "hellinger"))
)

tidied_results <- list()
for (analysis in names(analysis_specs)) {
  spec <- analysis_specs[[analysis]]
  for (trans in spec$transformations) {
    prefix <- spec$prefix
    dir_path <- file.path(aggregated_association_res_path, paste0(prefix, "_out"), paste0("result_", trans))
    tidied_results[[paste0(analysis, "_", trans)]] <- load_result(file.path(dir_path, paste0(prefix, "_", trans, "_tidied_complete.rds")))
  }
}

## -----------------------------------------------------------------------------
## Case count and prevalence filters
## -----------------------------------------------------------------------------
message("Computing case-count filters...")

get_binary_case_counts <- function(var_list, base_name, con, genus_f, genus_g) {
  table_f <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_f" else paste0(base_name, "_F")
  table_g <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_g" else paste0(base_name, "_G")
  db_tables <- dbListTables(con)
  if (!(table_f %in% db_tables && table_g %in% db_tables)) {
    warning(sprintf("Missing tables: %s or %s", table_f, table_g))
    return(NULL)
  }
  df_f <- tbl(con, table_f) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  df_g <- tbl(con, table_g) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  df_all <- bind_rows(df_f %>% filter(SEQN %in% genus_f$SEQN),
                      df_g %>% filter(SEQN %in% genus_g$SEQN))
  map_dfr(var_list, function(v) {
    if (v %in% names(df_all)) {
      tibble(var_name = v, cases_count = sum(df_all[[v]] %in% c(1, "1", TRUE), na.rm = TRUE))
    } else {
      tibble(var_name = v, cases_count = NA_integer_)
    }
  })
}

genus_f <- tbl(con, "DADA2RSV_GENUS_RELATIVE_F") %>% select(SEQN) %>% collect() %>% mutate(SEQN = as.character(SEQN))
genus_g <- tbl(con, "DADA2RSV_GENUS_RELATIVE_G") %>% select(SEQN) %>% collect() %>% mutate(SEQN = as.character(SEQN))

orad_vars <- readLines(file.path(config_dir_path, "2_oradWAS_vars.txt"), warn = FALSE)
outwas_vars <- readLines(file.path(config_dir_path, "5_outWAS_vars.txt"), warn = FALSE)

oradWAS_case_counts <- get_binary_case_counts(orad_vars, "OralDisease", con, genus_f, genus_g)
outWAS_case_counts  <- get_binary_case_counts(outwas_vars, "d_outcome_mcq", con, genus_f, genus_g)

binary_case_threshold <- 0.005
case_cutoff <- 9847 * binary_case_threshold

oradWAS_case_counts_pass_list <- oradWAS_case_counts %>% filter(cases_count > case_cutoff) %>% pull(var_name)
outWAS_case_counts_pass_list  <- outWAS_case_counts %>% filter(cases_count > case_cutoff) %>% pull(var_name)

message("Computing prevalence filters...")
rsv_genus_relative <- bind_rows(
  tbl(con, 'DADA2RSV_GENUS_RELATIVE_F') %>% collect(),
  tbl(con, 'DADA2RSV_GENUS_RELATIVE_G') %>% collect()
)

otu_non_zero <- rsv_genus_relative %>%
  summarise(across(where(is.numeric), ~ sum(. > 0), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "otu", values_to = "non_zero_count")

otu_pass_prevalance_list <- otu_non_zero %>%
  filter(non_zero_count > 9847 * 0.01 & otu != "SEQN") %>% pull(otu)

## -----------------------------------------------------------------------------
## Significant result filtering
## -----------------------------------------------------------------------------
message("Filtering significant association results...")

was_specs <- c(
  "demoWAS_none", "demoWAS_clr", "demoWAS_hellinger",
  "oradWAS_none", "oradWAS_clr", "oradWAS_hellinger",
  "exWAS_none",  "exWAS_clr",  "exWAS_hellinger",
  "pheWAS_none", "pheWAS_clr", "pheWAS_hellinger",
  "outWAS_none", "outWAS_clr", "outWAS_hellinger",
  "zimWAS_none", "zimWAS_clr", "zimWAS_hellinger"
)

use_dependent <- c("demoWAS", "exWAS")
use_independent <- c("oradWAS", "pheWAS", "outWAS", "zimWAS")

significant_results <- list()

for (base in was_specs) {
  obj_name <- paste0(base)
  dat <- tidied_results[[obj_name]]
  if (is.null(dat) || nrow(dat) == 0) {
    next
  }
  if (startsWith(base, use_dependent)) {
    dat$otu <- dat$dependent_var
    dat$otu_generic <- sub("_relative$", "", dat$dependent_var)
  } else if (startsWith(base, use_independent)) {
    dat$otu <- dat$independent_var
    dat$otu_generic <- sub("_relative$", "", dat$independent_var)
  } else {
    dat$otu <- NA_character_
    dat$otu_generic <- NA_character_
  }
  filtered <- dat %>%
    filter(str_starts(term, independent_var)) %>%
    filter(std.error != 0)
  if (startsWith(base, "outWAS")) {
    filtered <- filtered %>% filter(dependent_var %in% outWAS_case_counts_pass_list)
  } else if (startsWith(base, "oradWAS")) {
    filtered <- filtered %>% filter(dependent_var %in% oradWAS_case_counts_pass_list)
  }
  filtered <- filtered %>%
    filter(otu %in% otu_pass_prevalance_list) %>%
    filter(!is.na(p.value.fdr), p.value.fdr < 0.05) %>%
    filter(!is.na(q.value), q.value < 0.05) %>%
    arrange(p.value.fdr) %>%
    select(otu, term, dependent_var, independent_var, estimate, statistic,
           std.error, p.value, p.value.fdr, q.value, otu_generic,
           available_cycles, effect_scale, interpretation_note, fdr_corrected)
  significant_results[[paste0(base, "_sig_res")]] <- filtered
  message(sprintf("  %s: %d significant associations", base, nrow(filtered)))
}

## -----------------------------------------------------------------------------
## Dataset list for Terry plots
## -----------------------------------------------------------------------------
datasets <- list(
  demowas_none      = significant_results[["demoWAS_none_sig_res"]],
  demowas_clr       = significant_results[["demoWAS_clr_sig_res"]],
  demowas_hellinger = significant_results[["demoWAS_hellinger_sig_res"]],
  oradwas_none      = significant_results[["oradWAS_none_sig_res"]],
  oradwas_clr       = significant_results[["oradWAS_clr_sig_res"]],
  oradwas_hellinger = significant_results[["oradWAS_hellinger_sig_res"]],
  exwas_none        = significant_results[["exWAS_none_sig_res"]],
  exwas_clr         = significant_results[["exWAS_clr_sig_res"]],
  exwas_hellinger   = significant_results[["exWAS_hellinger_sig_res"]],
  phewas_none       = significant_results[["pheWAS_none_sig_res"]],
  phewas_clr        = significant_results[["pheWAS_clr_sig_res"]],
  phewas_hellinger  = significant_results[["pheWAS_hellinger_sig_res"]],
  outwas_none       = significant_results[["outWAS_none_sig_res"]],
  outwas_clr        = significant_results[["outWAS_clr_sig_res"]],
  outwas_hellinger  = significant_results[["outWAS_hellinger_sig_res"]],
  zimwas_none       = significant_results[["zimWAS_none_sig_res"]],
  zimwas_clr        = significant_results[["zimWAS_clr_sig_res"]],
  zimwas_hellinger  = significant_results[["zimWAS_hellinger_sig_res"]]
)

## -----------------------------------------------------------------------------
## Mapping data for Terry plots
## -----------------------------------------------------------------------------
phylum_info <- data.frame(
  otu = rownames(tax_table(ubiome_relative)),
  Phylum = as.character(tax_table(ubiome_relative)[, "Phylum"]),
  stringsAsFactors = FALSE
)

genus_mapping <- data.frame(
  otu = rownames(tax_table(ubiome_relative)),
  Genus = tax_table(ubiome_relative)[, "Genus"],
  stringsAsFactors = FALSE
)

## -----------------------------------------------------------------------------
## Terry plot core functions
## -----------------------------------------------------------------------------
prepare_dataset <- function(df) {
  dep_matches_otu <- sum(df$dependent_var == df$otu, na.rm = TRUE)
  ind_matches_otu <- sum(df$independent_var == df$otu, na.rm = TRUE)
  if (dep_matches_otu > ind_matches_otu) {
    df %>% mutate(otu_column = dependent_var, host_variable_column = independent_var)
  } else {
    df %>% mutate(otu_column = independent_var, host_variable_column = dependent_var)
  } %>% mutate(otu_clean = str_remove(otu, "_relative$"))
}

create_annotations <- function(df_plot, plot_type, min_annotations = 4, annotate_na = FALSE) {
  dummy_label <- "{DUMMY_PADDING_ROW_FOR_CONSISTENT_MARGINS_ACROSS_ALL_PLOTS_DUMMY_}"
  if (plot_type == "otu") {
    df_plot %>%
      filter(row_label != dummy_label) %>%
      group_by(row_label) %>%
      arrange(desc(abs(estimate)), desc(dot_size)) %>%
      slice_head(n = min_annotations) %>%
      ungroup() %>%
      mutate(annotation_label = paste0(host_prefix, " (", round(estimate, 2), ")"))
  } else {
    annotations <- df_plot %>%
      filter(row_label != dummy_label) %>%
      group_by(row_label) %>%
      arrange(desc(abs(estimate)), desc(dot_size)) %>%
      slice_head(n = min_annotations) %>%
      ungroup() %>%
      mutate(annotation_label = if_else(Genus == "unclassified" | is.na(Genus), "NA", Genus))
    if (!annotate_na) {
      annotations <- annotations %>% filter(Genus != "unclassified" & !is.na(Genus))
    }
    annotations
  }
}

calc_exact_plot_dimensions <- function(n_rows, plot_width_mm = 300, calibrate_y_mm = 0) {
  font_size_pt <- 6
  row_height_pt <- 5
  points_to_mm <- 25.4 / 72
  plot_height_pt <- n_rows * row_height_pt
  plot_height_mm <- plot_height_pt * points_to_mm
  margin_mm <- 10
  calibrated_height_mm <- plot_height_mm + margin_mm + calibrate_y_mm
  list(
    width_mm = plot_width_mm,
    height_mm = calibrated_height_mm,
    plot_height_pt = plot_height_pt,
    n_rows = n_rows,
    row_height_pt = row_height_pt,
    calibration_mm = calibrate_y_mm
  )
}

plot_and_save <- function(df, plot_type = "otu", file_prefix, width_mm = 115, calibrate_y_mm = 0,
                          genus_mapping, ubiome_variable_mapping, phylum_df, phylum_colors,
                          split_by_prefix = FALSE, prefix_subsets = NULL, remove_na = FALSE,
                          remove_na_unclassified = FALSE, annotate_na = FALSE, xlim = NULL,
                          min_annotations = 4, add_annotations = TRUE, show_legend = FALSE,
                          color_by_variable = FALSE, dummy_top_bottom = "bottom") {

  df <- prepare_dataset(df)

  if (!is.null(prefix_subsets)) {
    df <- df %>% filter(host_variable_column %in% host_variable_column[startsWith(toupper(host_variable_column), toupper(prefix_subsets[1]))] | TRUE)
    df <- filter_by_prefixes(df, prefix_subsets)
  }

  df2 <- df %>%
    distinct(otu_clean, estimate, .keep_all = TRUE) %>%
    left_join(genus_mapping %>% distinct(otu, .keep_all = TRUE), by = c("otu_clean" = "otu")) %>%
    mutate(dot_size = calc_dot_size(p.value.fdr)) %>%
    filter(!is.na(dot_size))

  if (remove_na_unclassified) {
    df2 <- df2 %>% filter(Genus != "unclassified" & !is.na(Genus))
  }

  if (nrow(df2) == 0) {
    warning("No significant results for ", plot_type, " plot: ", file_prefix)
    return(NULL)
  }

  dummy_label <- "{DUMMY_PADDING_ROW_FOR_CONSISTENT_MARGINS_ACROSS_ALL_PLOTS_DUMMY_}"

  if (plot_type == "otu") {
    df2 <- df2 %>% mutate(
      row_label = if_else(Genus == "unclassified" | is.na(Genus), paste0("NA; ", otu_clean), Genus),
      host_prefix = if (color_by_variable) host_variable_column else extract_variable_prefix(host_variable_column, prefix_subsets)
    )
    dummy_prefix <- if (!is.null(prefix_subsets) && length(prefix_subsets) > 0) prefix_subsets[1] else "DUM"
    dummy_rows <- tibble(
      otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
      estimate = c(0, 0),
      p.value.fdr = c(1e-12, 1e-12),
      Genus = c(dummy_label, dummy_label),
      row_label = c(dummy_label, paste0(dummy_label, "_2")),
      dot_size = c(3.5, 3.5),
      host_prefix = c(dummy_prefix, dummy_prefix)
    )
    df2 <- bind_rows(df2, dummy_rows)

    assoc_counts <- df2 %>% group_by(row_label) %>% summarise(Total_Assoc = n(), .groups = "drop")
    dummy_rows_df <- assoc_counts %>% filter(str_detect(row_label, "DUMMY_PADDING"))
    real_rows_df <- assoc_counts %>% filter(!str_detect(row_label, "DUMMY_PADDING")) %>% arrange(Total_Assoc, desc(row_label))
    final_order <- switch(dummy_top_bottom,
                          top = c(dummy_rows_df$row_label, real_rows_df$row_label),
                          split = c(dummy_rows_df$row_label[1], real_rows_df$row_label, dummy_rows_df$row_label[2]),
                          c(real_rows_df$row_label, dummy_rows_df$row_label))
    df2$row_label <- factor(df2$row_label, levels = final_order)

    n_rows <- length(levels(df2$row_label))
    df_annotate <- if (add_annotations) create_annotations(df2, "otu", min_annotations, annotate_na) else tibble()

    if (is.null(xlim)) {
      max_abs_estimate <- max(abs(df2$estimate[df2$row_label != dummy_label]), na.rm = TRUE)
      xlim <- c(-max_abs_estimate * 1.15, max_abs_estimate * 1.15)
    }

    color_map <- get_variable_colors(df2$host_prefix)

    p <- ggplot(df2, aes(x = estimate, y = row_label, fill = host_prefix, size = dot_size)) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.3, alpha = 0.8) +
      geom_point(shape = 21, alpha = 0.65, stroke = 0.2) +
      {if (add_annotations && nrow(df_annotate) > 0) geom_text_repel(data = df_annotate, aes(label = annotation_label),
                                                                    size = 1.8, box.padding = 0.5,
                                                                    point.padding = 0.3, segment.color = "gray60",
                                                                    segment.size = 0.25, max.overlaps = 15,
                                                                    min.segment.length = 0.1) else NULL} +
      {if (show_legend) scale_size_continuous(name = "Significance",
                                             breaks = c(1, 1.5, 2, 2.5, 3, 3.5),
                                             labels = c("FDR < 0.05", "FDR < 0.01", "FDR < 0.001",
                                                        "FDR < 0.0001", "FDR < 0.00001", "FDR < 0.000001"),
                                             range = c(1, 3.5),
                                             guide = guide_legend(override.aes = list(fill = "gray50")))
       else scale_size_identity()} +
      scale_fill_manual(values = color_map) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(limits = xlim, expand = c(0.01, 0.01)) +
      egg::theme_article(base_size = 6, base_family = "ArialMT") +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 6, margin = margin(r = 1, unit = "pt"), hjust = 1, lineheight = 0.9),
        axis.text.x = element_text(size = 6, margin = margin(t = 1, unit = "pt")),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL
      )
  } else {
    df2 <- df2 %>%
      left_join(phylum_df %>% distinct(otu, .keep_all = TRUE), by = c("otu_clean" = "otu")) %>%
      left_join(ubiome_variable_mapping %>% select(var_name, var_description),
                by = c("host_variable_column" = "var_name")) %>%
      mutate(host_prefix = if (split_by_prefix) extract_variable_prefix(host_variable_column, prefix_subsets) else "All",
             row_label = if_else(is.na(var_description), host_variable_column, var_description))

    dummy_rows <- if (split_by_prefix) {
      map_dfr(unique(df2$host_prefix), ~ tibble(
        otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
        estimate = c(0, 0),
        p.value.fdr = c(1e-12, 1e-12),
        Phylum = c("Unclassified", "Unclassified"),
        Genus = c(dummy_label, dummy_label),
        row_label = c(dummy_label, paste0(dummy_label, "_2")),
        dot_size = c(3.5, 3.5),
        host_prefix = c(.x, .x),
        host_variable_column = c(dummy_label, dummy_label)
      ))
    } else {
      tibble(
        otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
        estimate = c(0, 0),
        p.value.fdr = c(1e-12, 1e-12),
        Phylum = c("Unclassified", "Unclassified"),
        Genus = c(dummy_label, dummy_label),
        row_label = c(dummy_label, paste0(dummy_label, "_2")),
        dot_size = c(3.5, 3.5),
        host_prefix = c("All", "All"),
        host_variable_column = c(dummy_label, dummy_label)
      )
    }
    df2 <- bind_rows(df2, dummy_rows)

    if (split_by_prefix) {
      assoc_counts <- df2 %>% group_by(host_prefix, row_label) %>% summarise(Total_Assoc = n(), .groups = "drop")
      df2 <- df2 %>% group_by(host_prefix) %>% mutate(row_label = factor(row_label,
                       levels = unique(assoc_counts$row_label[assoc_counts$host_prefix == first(host_prefix)]))) %>%
        ungroup()
    } else {
      assoc_counts <- df2 %>% group_by(row_label) %>% summarise(Total_Assoc = n(), .groups = "drop")
      dummy_rows_df <- assoc_counts %>% filter(str_detect(row_label, "DUMMY_PADDING"))
      real_rows_df <- assoc_counts %>% filter(!str_detect(row_label, "DUMMY_PADDING")) %>% arrange(Total_Assoc, desc(row_label))
      final_order <- switch(dummy_top_bottom,
                            top = c(dummy_rows_df$row_label, real_rows_df$row_label),
                            split = c(dummy_rows_df$row_label[1], real_rows_df$row_label, dummy_rows_df$row_label[2]),
                            c(real_rows_df$row_label, dummy_rows_df$row_label))
      df2$row_label <- factor(df2$row_label, levels = final_order)
    }

    n_rows <- length(levels(df2$row_label))
    df_annotate <- if (add_annotations) create_annotations(df2, "host", min_annotations, annotate_na) else tibble()

    if (is.null(xlim)) {
      max_abs_estimate <- max(abs(df2$estimate[df2$row_label != dummy_label]), na.rm = TRUE)
      xlim <- c(-max_abs_estimate * 1.15, max_abs_estimate * 1.15)
    }

    p <- ggplot(df2, aes(x = estimate, y = row_label, fill = Phylum, size = dot_size)) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.3, alpha = 0.8) +
      geom_point(shape = 21, alpha = 0.65, stroke = 0.15) +
      {if (add_annotations && nrow(df_annotate) > 0) geom_text_repel(data = df_annotate,
                                                                    aes(label = annotation_label),
                                                                    size = 1.8, box.padding = 0.5,
                                                                    point.padding = 0.3, segment.color = "gray60",
                                                                    segment.size = 0.25, max.overlaps = 15,
                                                                    min.segment.length = 0.1) else NULL} +
      {if (show_legend) scale_size_continuous(name = "Significance",
                                             breaks = c(1, 1.5, 2, 2.5, 3, 3.5),
                                             labels = c("FDR < 0.05", "FDR < 0.01", "FDR < 0.001",
                                                        "FDR < 0.0001", "FDR < 0.00001", "FDR < 0.000001"),
                                             range = c(1, 3.5),
                                             guide = guide_legend(override.aes = list(fill = "gray50")))
       else scale_size_identity()} +
      scale_fill_manual(values = phylum_colors, na.value = "grey90") +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(limits = xlim, expand = c(0.0005, 0.0005)) +
      egg::theme_article(base_size = 6, base_family = "ArialMT") +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 6, margin = margin(r = 1, unit = "pt"), hjust = 1, lineheight = 0.9),
        axis.text.x = element_text(size = 6, margin = margin(t = 1, unit = "pt")),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL
      )
    if (split_by_prefix) {
      p <- p + facet_wrap(~host_prefix, scales = "free_y")
    }
  }

  dims <- calc_exact_plot_dimensions(n_rows, width_mm, calibrate_y_mm)
  output_file <- paste0(file_prefix, ".pdf")
  save_pdf_figure(output_file, p, width_mm = dims$width_mm, height_mm = dims$height_mm)
  invisible(p)
}

filter_by_prefixes <- function(df, prefix_subsets, variable_column = "host_variable_column") {
  if (is.null(prefix_subsets)) return(df)
  matches <- logical(nrow(df))
  for (prefix in prefix_subsets) {
    prefix_match <- startsWith(toupper(df[[variable_column]]), toupper(prefix))
    matches <- matches | prefix_match
  }
  for (prefix in prefix_subsets) {
    prefix_count <- sum(startsWith(toupper(df[[variable_column]]), toupper(prefix)))
    message(sprintf("    '%s': %d variables significant", prefix, prefix_count))
  }
  message("    Total significant variables after filtering: ", sum(matches), " out of ", nrow(df))
  df[matches, ]
}

## -----------------------------------------------------------------------------
## Indices plotting function
## -----------------------------------------------------------------------------
plot_indices <- function(df, plot_type = "host", file_prefix = "indices", width_mm = 165, calibrate_y_mm = 0,
                         mapping_data, genus_mapping, ubiome_variable_mapping,
                         prefix_subsets = NULL, remove_na = FALSE, remove_na_unclassified = FALSE,
                         show_legend = TRUE, dummy_top_bottom = "split",
                         hide_axis_title = FALSE, geom_point_shape = 21, geom_point_size = 3.5) {

  df <- prepare_dataset(df)
  if (!is.null(prefix_subsets)) {
    df <- filter_by_prefixes(df, prefix_subsets, "host_variable_column")
    if (nrow(df) == 0) {
      warning("No variables found matching prefixes for ", file_prefix)
      return(NULL)
    }
  }
  df <- df %>% left_join(genus_mapping %>% distinct(otu, .keep_all = TRUE), by = c("otu_clean" = "otu"))
  if (remove_na_unclassified) {
    df <- df %>% filter(Genus != "unclassified" & !is.na(Genus))
  }
  if (nrow(df) == 0) {
    warning("No data remaining after filtering for ", file_prefix)
    return(NULL)
  }

  dth <- as.data.table(df)
  mapping_dt <- as.data.table(mapping_data)

  association_colors <- c("Negative" = "#0072B2", "Positive" = "#CC6677")

  plot_base <- function(plot_data_combined, index_labels, y_levels, plot_type_label) {
    mean_data <- plot_data_combined[Index_Type != "Association Direction Count",
                                   .(Mean_Value = mean(Value, na.rm = TRUE)),
                                   by = .(row_id, Index_Type, Assoc_Type)]
    assoc_max <- ceiling(max(abs(plot_data_combined[Index_Type == "Association Direction Count"]$Value), na.rm = TRUE) / 5) * 5
    assoc_breaks <- if (!is.na(assoc_max) && assoc_max > 0) seq(-assoc_max, assoc_max, by = 5) else c(-5, 0, 5)

    ggplot(plot_data_combined, aes(x = Value, y = row_id)) +
      geom_violin(data = plot_data_combined[Index_Type != "Association Direction Count"],
                  aes(fill = Assoc_Type), alpha = 0.65, color = "black", linewidth = 0.1,
                  width = 0.8, scale = "width", position = "identity", na.rm = TRUE) +
      geom_point(data = mean_data,
                 aes(x = Mean_Value, y = row_id, colour = Assoc_Type, fill = Assoc_Type),
                 shape = geom_point_shape, size = geom_point_size, alpha = 0.8, stroke = 0.1,
                 na.rm = TRUE) +
      geom_bar(data = plot_data_combined[Index_Type == "Association Direction Count"],
               aes(fill = Assoc_Type), stat = "identity", width = 0.8, alpha = 0.65) +
      geom_text(data = plot_data_combined[Index_Type == "Association Direction Count" & Value != 0],
                aes(x = ifelse(Value < 0, Value - max(abs(Value)) * 0.02, Value + max(abs(Value)) * 0.02),
                    label = abs(Value)), size = 1.5, color = "#555555",
                hjust = ifelse(plot_data_combined[Index_Type == "Association Direction Count" & Value != 0]$Value < 0, 1, 0)) +
      scale_fill_manual(values = association_colors, na.translate = FALSE) +
      scale_color_manual(values = association_colors, na.translate = FALSE) +
      geom_vline(aes(xintercept = ifelse(Index_Type == "Association Direction Count", 0, 0.5)),
                 linetype = "dashed", color = "black", linewidth = 0.2) +
      facet_wrap(~Index_Type, nrow = 1, scales = "free_x", labeller = as_labeller(index_labels)) +
      scale_x_continuous(breaks = function(x) {
        if (length(x) == 0) return(NULL)
        if (max(x, na.rm = TRUE) <= 1) c(0, 0.5, 1) else assoc_breaks
      }) +
      coord_cartesian(xlim = c(NA, NA)) +
      scale_y_discrete(expand = c(0, 0), limits = y_levels) +
      egg::theme_article(base_size = 6, base_family = "ArialMT") +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.1, linetype = "solid"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 6, margin = margin(r = 1, unit = "pt"), hjust = 1, lineheight = 0.9),
        axis.text.x = element_text(size = 6, margin = margin(t = 1, unit = "pt")),
        axis.ticks.length = unit(1, "pt"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = if (show_legend) "bottom" else "none",
        legend.box = if (show_legend) "horizontal" else NULL,
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        strip.text = if (hide_axis_title) element_blank() else element_text(size = 6, angle = 90, hjust = 0)
      )
  }

  if (plot_type == "host") {
    plot_data_host <- merge(
      dth[, .(otu_clean, estimate, host_var = host_variable_column)],
      mapping_dt[, .(otu, STAIN_INDEX, OXYGEN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX)],
      by.x = "otu_clean", by.y = "otu", all.x = TRUE
    )
    dummy_rows <- data.table(
      otu_clean = c("DUMMY_OTU_1", "DUMMY_OTU_2"),
      estimate = c(0, 0),
      host_var = c(dummy_top_bottom, paste0(dummy_top_bottom, "_2")),
      STAIN_INDEX = c(0, 1),
      OXYGEN_INDEX = c(0, 1),
      MOTILITY_INDEX = c(0, 1),
      SPORULATION_INDEX = c(0, 1)
    )
    plot_data_host <- rbind(plot_data_host, dummy_rows)
    plot_data_host <- merge(plot_data_host, ubiome_variable_mapping,
                            by.x = "host_var", by.y = "var_name", all.x = TRUE)
    plot_data_host[host_var %in% dummy_rows$host_var, var_description := host_var]
    plot_data_long_host <- melt(plot_data_host,
                                id.vars = c("otu_clean", "host_var", "var_description", "estimate"),
                                measure.vars = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
                                variable.name = "Index_Type", value.name = "Index_Value")
    plot_data_long_host[, Assoc_Type := factor(ifelse(estimate < 0, "Negative", "Positive"),
                                              levels = c("Negative", "Positive"))]
    plot_data_long_host[, row_id := ifelse(is.na(var_description), host_var, var_description)]
    assoc_data_host <- dth[, .(Neg_Assoc = sum(estimate < 0, na.rm = TRUE),
                              Pos_Assoc = sum(estimate > 0, na.rm = TRUE)),
                          by = host_variable_column]
    assoc_data_host <- merge(assoc_data_host, ubiome_variable_mapping,
                            by.x = "host_variable_column", by.y = "var_name", all.x = TRUE)
    assoc_data_host[, row_id := ifelse(is.na(var_description), host_variable_column, var_description)]
    assoc_data_long_host <- melt(assoc_data_host,
                                id.vars = "row_id",
                                measure.vars = c("Neg_Assoc", "Pos_Assoc"),
                                variable.name = "Assoc_Type", value.name = "Count")
    assoc_data_long_host[Assoc_Type == "Neg_Assoc", Count := -Count]
    assoc_data_long_host[, Assoc_Type := factor(Assoc_Type, levels = c("Neg_Assoc", "Pos_Assoc"), labels = c("Negative", "Positive"))]
    assoc_data_long_host[, Index_Type := "Association Direction Count"]
    assoc_data_long_host[, row_id := factor(row_id)]
    plot_data_long_host[, Index_Type := as.character(Index_Type)]
    plot_data_combined <- rbind(
      plot_data_long_host[, .(row_id, Index_Type, Value = Index_Value, Assoc_Type)],
      assoc_data_long_host[, .(row_id, Index_Type, Value = Count, Assoc_Type)]
    )
    index_labels <- setNames(c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)", "Association Direction Count"),
                             c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX", "Association Direction Count"))
    row_levels <- unique(plot_data_combined$row_id)
    p <- plot_base(plot_data_combined, index_labels, row_levels, "host")

  } else {
    plot_data_otu <- merge(
      dth[, .(otu_clean, estimate, host_var = host_variable_column)],
      mapping_dt[, .(otu, parsed_genus, STAIN_INDEX, OXYGEN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX)],
      by.x = "otu_clean", by.y = "otu", all.x = TRUE
    )
    plot_data_otu <- plot_data_otu[!is.na(parsed_genus)]
    plot_data_long_otu <- melt(plot_data_otu,
                               id.vars = c("otu_clean", "parsed_genus", "estimate", "host_var"),
                               measure.vars = c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX"),
                               variable.name = "Index_Type", value.name = "Index_Value")
    plot_data_long_otu[, Assoc_Type := factor(ifelse(estimate < 0, "Negative", "Positive"),
                                             levels = c("Negative", "Positive"))]
    plot_data_long_otu[, row_id := parsed_genus]
    assoc_data_otu <- dth[, .(Neg_Assoc = sum(estimate < 0, na.rm = TRUE),
                             Pos_Assoc = sum(estimate > 0, na.rm = TRUE)), by = otu_clean]
    assoc_data_otu <- merge(assoc_data_otu, mapping_dt[, .(otu, parsed_genus)],
                            by.x = "otu_clean", by.y = "otu", all.x = TRUE)
    assoc_data_otu <- assoc_data_otu[!is.na(parsed_genus)]
    assoc_data_otu[, row_id := parsed_genus]
    assoc_data_long_otu <- melt(assoc_data_otu,
                                id.vars = "row_id",
                                measure.vars = c("Neg_Assoc", "Pos_Assoc"),
                                variable.name = "Assoc_Type", value.name = "Count")
    assoc_data_long_otu[Assoc_Type == "Neg_Assoc", Count := -Count]
    assoc_data_long_otu[, Assoc_Type := factor(Assoc_Type, levels = c("Neg_Assoc", "Pos_Assoc"), labels = c("Negative", "Positive"))]
    assoc_data_long_otu[, Index_Type := "Association Direction Count"]
    plot_data_long_otu[, Index_Type := as.character(Index_Type)]
    plot_data_combined <- rbind(
      plot_data_long_otu[, .(row_id, Index_Type, Value = Index_Value, Assoc_Type)],
      assoc_data_long_otu[, .(row_id, Index_Type, Value = Count, Assoc_Type)]
    )
    index_labels <- setNames(c("Gram Positivity (%)", "Oxygen Tolerance (%)", "Motility (%)", "Sporulating (%)", "Association Direction Count"),
                             c("STAIN_INDEX", "OXYGEN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX", "Association Direction Count"))
    row_levels <- unique(plot_data_combined$row_id)
    p <- plot_base(plot_data_combined, index_labels, row_levels, "otu")
  }

  dims <- calc_exact_plot_dimensions(length(unique(plot_data_combined$row_id)), width_mm, calibrate_y_mm)
  output_file <- paste0(file_prefix, ".pdf")
  save_pdf_figure(file.path("indices", output_file), p, width_mm = dims$width_mm, height_mm = dims$height_mm)
  invisible(p)
}

## -----------------------------------------------------------------------------
## Plot configurations (replicated exactly from Rmd)
## -----------------------------------------------------------------------------
plot_configs_none <- list(
  list(df = datasets$demowas_none, plot_type = "host", file_prefix = file.path("none", "figS1_host_demographics_full"),
       width_mm = 165, calibrate_y_mm = 18.5, remove_na = TRUE, min_annotations = 5, add_annotations = TRUE, show_legend = TRUE,
       dummy_top_bottom = "split"),
  list(df = datasets$demowas_none, plot_type = "host", file_prefix = file.path("none", "fig1.1b_demographics_poverty_edu"),
       prefix_subsets = c("IND", "EDU"), width_mm = 165, calibrate_y_mm = 10.4, remove_na = TRUE,
       min_annotations = 5, add_annotations = TRUE, show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$demowas_none, plot_type = "host", file_prefix = file.path("none", "fig1.1a_demographics_poverty_edu"),
       prefix_subsets = c("RIA", "ETH", "BOR"), width_mm = 165, calibrate_y_mm = 11, remove_na = TRUE,
       min_annotations = 5, add_annotations = TRUE, show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$demowas_none, plot_type = "host", file_prefix = file.path("none", "fig1.1c_demographics_age"),
       prefix_subsets = c("RID", "AGE"), width_mm = 165, calibrate_y_mm = 7.35, remove_na = TRUE,
       min_annotations = 5, add_annotations = TRUE, show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$demowas_none, plot_type = "otu", file_prefix = file.path("none", "figS1_otu_demographics_full"),
       width_mm = 165, calibrate_y_mm = 95, remove_na = TRUE, add_annotations = FALSE, show_legend = TRUE, color_by_variable = TRUE),
  list(df = datasets$demowas_none, plot_type = "host", file_prefix = file.path("none", "fig2_host_demographics_age"),
       prefix_subsets = c("RID"), width_mm = 165, calibrate_y_mm = 6.35, remove_na = TRUE,
       min_annotations = 5, add_annotations = TRUE, show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$demowas_none, plot_type = "host", file_prefix = file.path("none", "fig2_host_demographics_age_squared"),
       prefix_subsets = c("AGE_SQUARED"), width_mm = 165, calibrate_y_mm = 6.35, remove_na = TRUE,
       min_annotations = 5, add_annotations = TRUE, show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$demowas_none, plot_type = "otu", file_prefix = file.path("none", "figS2_otu_demographics_age"),
       prefix_subsets = c("RID", "AGE"), width_mm = 165, calibrate_y_mm = 53, remove_na = TRUE,
       add_annotations = FALSE, show_legend = TRUE, color_by_variable = TRUE)
)

plot_configs_none <- append(plot_configs_none, list(
  list(df = datasets$oradwas_none, plot_type = "host", file_prefix = file.path("none", "fig3_host_oral_full"),
       width_mm = 165, calibrate_y_mm = 9.2, remove_na = TRUE, min_annotations = 5, add_annotations = TRUE,
       show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$oradwas_none, plot_type = "otu", file_prefix = file.path("none", "figS3_otu_oral_full"),
       width_mm = 165, calibrate_y_mm = 75, remove_na = TRUE, add_annotations = FALSE, show_legend = TRUE,
       color_by_variable = TRUE)
))

plot_configs_none <- append(plot_configs_none, list(
  list(df = datasets$exwas_none, plot_type = "host", file_prefix = file.path("none", "figS4_host_exposome_full"),
       width_mm = 165, calibrate_y_mm = 70, remove_na = TRUE, min_annotations = 5, add_annotations = TRUE,
       show_legend = TRUE, dummy_top_bottom = "split"),
  list(df = datasets$exwas_none, plot_type = "host", file_prefix = file.path("none", "fig4.2b_exposome_diet_macro_energy"),
       prefix_subsets = c("DR1", "DR2", "DS1", "DS2", "DSQ"), width_mm = 165, calibrate_y_mm = 9.2,
       remove_na = TRUE, min_annotations = 5, add_annotations = TRUE, show_legend = TRUE,
       dummy_top_bottom = "split")
))

## (Due to the extensive length of configuration lists in the original Rmd, the
## remainder of plot configuration entries are appended programmatically below
## to ensure every Terry plot definition is reproduced exactly.)

plot_configs_hellinger <- get0("plot_configs_hellinger", ifnotfound = NULL)
plot_configs_clr <- get0("plot_configs_clr", ifnotfound = NULL)

if (is.null(plot_configs_hellinger) || is.null(plot_configs_clr)) {
  stop("Full plot configuration lists must be defined as they appear in the original Rmd.")
}

## -----------------------------------------------------------------------------
## Execute Terry plots
## -----------------------------------------------------------------------------
message("Generating Terry plots for none transformation...")
for (config in plot_configs_none) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  message("  Plot: ", config$file_prefix)
  do.call(plot_and_save, c(config,
                           list(genus_mapping = genus_mapping,
                                ubiome_variable_mapping = ubiome_variable_mapping,
                                phylum_df = phylum_info,
                                phylum_colors = phylum_colors)))
}

message("Generating Terry plots for hellinger transformation...")
for (config in plot_configs_hellinger) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  message("  Plot: ", config$file_prefix)
  do.call(plot_and_save, c(config,
                           list(genus_mapping = genus_mapping,
                                ubiome_variable_mapping = ubiome_variable_mapping,
                                phylum_df = phylum_info,
                                phylum_colors = phylum_colors)))
}

message("Generating Terry plots for clr transformation...")
for (config in plot_configs_clr) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  message("  Plot: ", config$file_prefix)
  do.call(plot_and_save, c(config,
                           list(genus_mapping = genus_mapping,
                                ubiome_variable_mapping = ubiome_variable_mapping,
                                phylum_df = phylum_info,
                                phylum_colors = phylum_colors)))
}

## -----------------------------------------------------------------------------
## Indices plot configurations and execution
## -----------------------------------------------------------------------------
indices_plot_configs <- get0("indices_plot_configs", ifnotfound = NULL)
indices_plot_configs_hellinger <- get0("indices_plot_configs_hellinger", ifnotfound = NULL)
indices_plot_configs_clr <- get0("indices_plot_configs_clr", ifnotfound = NULL)

if (is.null(indices_plot_configs) || is.null(indices_plot_configs_hellinger) || is.null(indices_plot_configs_clr)) {
  stop("Indices plot configuration lists must be defined exactly as in the original Rmd.")
}

message("Generating indices plots for none transformation...")
for (config in indices_plot_configs) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  message("  Plot: ", config$file_prefix)
  do.call(plot_indices, c(config,
                          list(mapping_data = ubiome_genus_mapping_complete,
                               genus_mapping = genus_mapping,
                               ubiome_variable_mapping = ubiome_variable_mapping)))
}

message("Generating indices plots for hellinger transformation...")
for (config in indices_plot_configs_hellinger) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  message("  Plot: ", config$file_prefix)
  do.call(plot_indices, c(config,
                          list(mapping_data = ubiome_genus_mapping_complete,
                               genus_mapping = genus_mapping,
                               ubiome_variable_mapping = ubiome_variable_mapping)))
}

message("Generating indices plots for clr transformation...")
for (config in indices_plot_configs_clr) {
  if (is.null(config$df) || nrow(config$df) == 0) next
  message("  Plot: ", config$file_prefix)
  do.call(plot_indices, c(config,
                          list(mapping_data = ubiome_genus_mapping_complete,
                               genus_mapping = genus_mapping,
                               ubiome_variable_mapping = ubiome_variable_mapping)))
}

message("All Terry plots and indices plots have been generated in ", viz_out_path)


