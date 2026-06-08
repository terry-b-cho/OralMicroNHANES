#!/usr/bin/env Rscript

# =============================================================================
# NETWORK ANALYSIS SCRIPT
# =============================================================================
# Builds host-microbe association networks from the intermediates produced by
# 1_load_all_categories.R.
#
# All outputs are written under:
#   results/analyses_results/8_network_analyses_out/
#
# Environment: R >= 4.5 with phyloseq, igraph, ggraph, ggplot2, dplyr, tidyr,
# stringr, purrr, gridExtra, extrafont, ggrepel, scales, rlang.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

cat("=== MICROBIAL ASSOCIATION NETWORK ANALYSIS ===\n")

# Load required libraries
library(phyloseq)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(gridExtra)
library(extrafont)
library(ggrepel)
library(scales)
library(rlang)

set.seed(42)

utils::globalVariables(c(
  "dependent_var", "independent_var", "feature_id", "phenotype", "otu_clean",
  "estimate", "otu", "otu1", "otu2", "est1", "est2", "pct_positive",
  "same_direction", "both_positive", "both_negative", "agreement_score",
  "Phylum", "Genus", "parsed_genus", "layout_weight", "visual_weight",
  "edge_color", "edge_colour_value", "edge_width_value", "final_phylum",
  "final_genus", "node_size", "genus_label", "node_phylum", "node_size_value",
  "node_label_value"
))

otu1 <- otu2 <- otu <- edge_color <- edge_colour_value <- edge_width_value <- visual_weight <- final_phylum <- node_size <- genus_label <- node_phylum <- node_size_value <- node_label_value <- NULL

# Load fonts
loadfonts(device = "pdf", quiet = TRUE)

# =============================================================================
# CONFIGURATION
# =============================================================================

base_path <- PROJECT_ROOT
data_path <- file.path(base_path, "results/analyses_results/8_network_analyses_out/intermediate")
output_path <- file.path(base_path, "results/analyses_results/8_network_analyses_out")

# Create output directories
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters (EXACTLY as original, prevalence updated to 5%)
transformations <- c("none", "clr")
alpha_threshold <- 0.05
min_shared_traits <- 2
# min_prevalence <- 0.01
# min_prevalence <- 0.05
min_prevalence <- 0.10

cat("Configuration:\n")
cat("  Transformations:", paste(transformations, collapse = ", "), "\n")
cat("  Alpha threshold:", alpha_threshold, "\n")
cat("  Min shared traits:", min_shared_traits, "\n")
cat("  Min prevalence:", min_prevalence, "\n")
cat("  Output path:", output_path, "\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

# Load processed data
phyloseq_objects <- readRDS(file.path(data_path, "phyloseq_objects.rds"))
significant_results <- readRDS(file.path(data_path, "significant_results.rds"))
variable_description <- readRDS(file.path(data_path, "variable_description.rds"))
otu_pass_prevalence_list <- readRDS(file.path(data_path, "otu_pass_prevalence_list.rds"))
variable_groups <- readRDS(file.path(data_path, "variable_groups.rds"))
sample_data_subsets <- readRDS(file.path(data_path, "sample_data_subsets.rds"))
ubiome_genus_mapping_complete <- readRDS(file.path(data_path, "ubiome_genus_mapping_complete.rds"))
oradWAS_case_pass_list <- readRDS(file.path(data_path, "oradWAS_case_pass_list.rds"))
outWAS_case_pass_list <- readRDS(file.path(data_path, "outWAS_case_pass_list.rds"))

cat("  [done] Loaded phyloseq objects:", length(phyloseq_objects), "\n")
cat("  [done] Loaded significant results:", length(significant_results), "\n")
cat("  [done] Loaded variable descriptions:", nrow(variable_description), "\n")
cat("  [done] Loaded prevalence list:", length(otu_pass_prevalence_list), "OTUs\n")
cat("  [done] Loaded variable groups:", length(variable_groups), "groups\n")
cat("  [done] Loaded sample-data subsets:", length(sample_data_subsets), "datasets\n")
cat("  [done] Loaded genus mapping rows:", nrow(ubiome_genus_mapping_complete), "\n\n")

# =============================================================================
# DEFINE COLOR SCHEMES (EXACTLY AS ORIGINAL)
# =============================================================================

# Exact Kelly color scheme for phylum colors
phylum_colors <- c(
  "Firmicutes"           = "#F38400",  # kelly_3   - orange
  "Bacteroidetes"        = "#0067A5",  # kelly_10  - blue
  "Actinobacteria"       = "#8DB600",  # kelly_17  - yellow green
  "Proteobacteria"       = "#E68FAC",  # kelly_9   - pink
  "Fusobacteria"         = "#BE0032",  # kelly_5   - red
  "Spirochaetae"         = "#F3C300",  # kelly_1   - yellow
  "Cyanobacteria"        = "#875692",  # kelly_2   - purple
  "Acidobacteria"        = "#F6A600",  # kelly_13  - orange yellow
  "Candidate division SR1" = "#2B3D26",  # kelly_20  - forest green
  "Planctomycetes"       = "#332288",  # safe_violet
  "Saccharibacteria"     = "#B3446C",  # kelly_14  - purple red
  "Synergistetes"        = "#A1CAF1",  # kelly_4   - light blue
  "Tenericutes"          = "#654522",  # kelly_18  - dark brown
  "Verrucomicrobia"      = "#C2B280",  # kelly_6   - tan
  "unclassified"         = "#DDDDDD", # pale_grey
  "NA"                  = "#DDDDDD", # very light grey
  "TM7"                 = "#000000"  # Contrast Black
)

# Color scheme for edges (OrBl_div colors)
orbl_colors <- c("#2B5B8A", "#4071A0", "#5789B6", "#6FA3CB", "#A2BCCF", 
                 "#D8D4C9", "#F0AC72", "#EF8530", "#DA6524", "#BE4E21", "#9E3D21")

cat("Color schemes defined:\n")
cat("  Phylum colors:", length(phylum_colors), "colors\n")
cat("  Edge colors:", length(orbl_colors), "colors\n\n")

variable_description_map <- variable_description %>%
  filter(!is.na(.data$Variable.Name)) %>%
  mutate(
    Variable.Name = as.character(.data$Variable.Name),
    Variable.Description = stringr::str_squish(as.character(.data$Variable.Description))
  ) %>%
  filter(!is.na(.data$Variable.Description), .data$Variable.Description != "") %>%
  arrange(.data$Variable.Name) %>%
  distinct(.data$Variable.Name, .keep_all = TRUE)

variable_description_map <- setNames(variable_description_map$Variable.Description,
                                     variable_description_map$Variable.Name)

get_variable_description <- function(var_name) {
  desc <- variable_description_map[[var_name]]
  if (!is.null(desc) && !is.na(desc) && nzchar(desc)) {
    return(desc)
  }
  var_name
}

sanitize_for_filename <- function(text) {
  if (is.null(text) || is.na(text) || !nzchar(text)) {
    return("no_description")
  }
  cleaned <- stringr::str_replace_all(text, "[^A-Za-z0-9]+", "_")
  cleaned <- stringr::str_replace_all(cleaned, "^_+|_+$", "")
  if (!nzchar(cleaned)) cleaned <- "no_description"
  stringr::str_trunc(cleaned, 120, ellipsis = "")
}

# =============================================================================
# HELPER FUNCTIONS (EXACTLY AS ORIGINAL)
# =============================================================================

# Map agreement percentage to OrBl palette (11 buckets)
map_to_orbl_color <- function(pct_positive) {
  idx <- pmax(1, pmin(length(orbl_colors), floor(pct_positive * (length(orbl_colors) - 1)) + 1))
  orbl_colors[idx]
}

# Function to create nice variable names
nice <- function(var_name) {
  var_name %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

clean_feature_id <- function(x) stringr::str_remove(x, "_relative$")

`%||%` <- function(x, y) if (!is.null(x)) x else y

# Function to get taxonomy table (EXACTLY as original)
tax_tbl <- function(ps) {
  df <- as.data.frame(unclass(tax_table(ps)), stringsAsFactors = FALSE)
  df$OTU <- rownames(df)
  df
}

# Function to filter taxa based on prevalence (EXACTLY as original)
# keep_taxa <- function(ps, min_prev = 0.01, trans = "none"){
# keep_taxa <- function(ps, min_prev = 0.05, trans = "none"){
keep_taxa <- function(ps, min_prev = 0.10, trans = "none"){
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  
  # EXACTLY as original: For CLR/lognorm use variance > 0, otherwise use prevalence
  ok1 <- if (trans %in% c("clr", "lognorm"))
            apply(m, 1, var, na.rm = TRUE) > 0
          else
            rowSums(m > 0) / ncol(m) >= min_prev
  
  ok2 <- !is.na(tax_tbl(ps)$Genus)
  which(ok1 & ok2)
}

# Function to create phylum colors with alpha (EXACTLY as original)
make_phylum_col <- function(v, alpha = 0.65){
  standard <- trimws(v)
  standard[standard == ""] <- "unclassified"
  standard[is.na(standard)] <- "NA"
  standard[!(standard %in% names(phylum_colors))] <- "unclassified"
  unique_standard <- unique(standard)
  color_vec <- phylum_colors[unique_standard]
  color_vec <- grDevices::adjustcolor(color_vec, alpha.f = alpha)
  names(color_vec) <- unique_standard
  color_vec
}

format_genus_label <- function(genus_vec) {
  vapply(genus_vec, function(g) {
    if (is.na(g) || g == "") {
      return(g)
    }
    parts <- strsplit(g, "_", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      paste(parts[-length(parts)], collapse = "_")
    } else {
      g
    }
  }, character(1), USE.NAMES = FALSE)
}

cat("Helper functions defined\n")

# =============================================================================
# NETWORK ANALYSIS FUNCTIONS (EXACTLY AS ORIGINAL)
# =============================================================================

analyze_network_for_variable <- function(sig_res_table, host_var_name, ps, trans,
                                         min_shared = 2, title_prefix = "") {

  if (is.null(sig_res_table) || !nrow(sig_res_table)) {
    cat("    Analyzing:", host_var_name, "... no associations\n")
    return(NULL)
  }

  var_associations <- sig_res_table %>%
    filter((.data$dependent_var == host_var_name | .data$independent_var == host_var_name)) %>%
    mutate(
      feature_id = if_else(is.na(.data$feature_id), .data$phenotype, .data$feature_id),
      otu = if_else(is.na(.data$otu_clean) | .data$otu_clean == "", clean_feature_id(.data$feature_id), .data$otu_clean),
      effect = .data$estimate
    ) %>%
    filter(!is.na(.data$otu))

  n_assoc <- nrow(var_associations)
  cat("    Analyzing:", host_var_name, "...", n_assoc, "associations")

  if (n_assoc < min_shared) {
    cat(" - insufficient associations\n")
    return(NULL)
  }

  otus <- unique(var_associations$otu)
  n_otus <- length(otus)
  cat(",", n_otus, "unique OTUs")
  if (n_otus < 2) {
    cat(" - insufficient OTUs\n")
    return(NULL)
  }

  otu_pairs <- expand.grid(otu1 = otus, otu2 = otus, stringsAsFactors = FALSE) %>%
    filter(.data$otu1 < .data$otu2)

  cat(",", nrow(otu_pairs), "potential pairs")

  edge_list <- otu_pairs %>%
    left_join(var_associations %>% transmute(otu = .data$otu, est1 = .data$effect), by = c("otu1" = "otu")) %>%
    left_join(var_associations %>% transmute(otu = .data$otu, est2 = .data$effect), by = c("otu2" = "otu")) %>%
    filter(!is.na(.data$est1), !is.na(.data$est2)) %>%
    mutate(
      same_direction = sign(.data$est1) == sign(.data$est2),
      both_positive = .data$est1 > 0 & .data$est2 > 0,
      both_negative = .data$est1 < 0 & .data$est2 < 0,
      agreement_score = sqrt(abs(.data$est1 * .data$est2)),
      pct_positive = case_when(
        .data$both_positive ~ 1,
        .data$both_negative ~ 0,
        TRUE ~ 0.5
      ),
      edge_color = map_to_orbl_color(.data$pct_positive),
      layout_weight = ifelse(.data$same_direction, .data$agreement_score * 2, .data$agreement_score * 0.2),
      visual_weight = pmax(0.1, .data$agreement_score)
    )

  cat(",", nrow(edge_list), "valid edges")

  if (nrow(edge_list) < min_shared) {
    cat(" - insufficient edges\n")
    return(NULL)
  }

  node_stats <- var_associations %>%
    group_by(.data$otu) %>%
    summarise(
      association_strength = sum(abs(.data$effect)),
      mean_effect = mean(.data$effect),
      .groups = "drop"
    )

  taxonomy_df <- tax_tbl(ps)
  nodes_df <- node_stats %>%
    left_join(taxonomy_df, by = c("otu" = "OTU")) %>%
    left_join(ubiome_genus_mapping_complete %>% transmute(otu = .data$otu, parsed_genus = .data$parsed_genus), by = "otu") %>%
    mutate(
      final_phylum = coalesce(.data$Phylum, "unknown"),
      final_genus = coalesce(.data$parsed_genus, .data$Genus, .data$otu),
      node_size = .data$association_strength
    )

  if (nrow(nodes_df) > 0) {
    if (length(unique(nodes_df$node_size)) > 1) {
      nodes_df$node_size <- scales::rescale(nodes_df$node_size, to = c(3, 10))
    } else {
      nodes_df$node_size <- 5
    }
  }

  node_metadata <- nodes_df %>%
    mutate(
      genus_label = .data$final_genus
    )
  node_metadata$name <- node_metadata$otu
  node_metadata$node_phylum <- trimws(node_metadata$final_phylum)
  node_metadata$node_phylum[node_metadata$node_phylum == ""] <- "unclassified"
  node_metadata$node_phylum[is.na(node_metadata$node_phylum)] <- "NA"
  node_metadata$node_phylum[!(node_metadata$node_phylum %in% names(phylum_colors))] <- "unclassified"
  node_metadata$node_size_value <- node_metadata$node_size
  node_metadata$node_label_value <- dplyr::case_when(
    is.na(node_metadata$genus_label) ~ NA_character_,
    tolower(node_metadata$genus_label) %in% c("unclassified", "na") ~ NA_character_,
    TRUE ~ format_genus_label(node_metadata$genus_label)
  )

  phylum_palette <- make_phylum_col(node_metadata$node_phylum)
  node_metadata$node_color <- phylum_palette[node_metadata$node_phylum]

  graph_df <- edge_list %>%
    transmute(from = .data$otu1, to = .data$otu2, layout_weight = .data$layout_weight,
              edge_width_value = .data$visual_weight, edge_colour_value = .data$edge_color)

  set.seed(1)
  network_graph <- graph_from_data_frame(
    d = graph_df,
    vertices = node_metadata,
    directed = FALSE
  )

  E(network_graph)$layout_weight <- graph_df$layout_weight
  E(network_graph)$visual_weight <- graph_df$edge_width_value
  E(network_graph)$edge_color <- graph_df$edge_colour_value

  V(network_graph)$node_phylum <- node_metadata$node_phylum
  V(network_graph)$genus_label <- node_metadata$node_label_value
  V(network_graph)$node_size <- node_metadata$node_size_value

  set.seed(42)
  layout_coords <- layout_with_fr(network_graph, weights = E(network_graph)$layout_weight)

  description_text <- get_variable_description(host_var_name)
  plot_title <- sprintf("%s | %s — %s", title_prefix, host_var_name, description_text)
  plot_subtitle <- sprintf("%d OTUs • %d edges", vcount(network_graph), ecount(network_graph))

  network_plot <- ggraph(network_graph, layout = layout_coords) +
    geom_edge_link(aes(color = !!sym("edge_colour_value"), edge_width = !!sym("edge_width_value")), alpha = 0.325, show.legend = FALSE) +
    scale_edge_color_identity() +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_point(aes(fill = !!sym("node_phylum"), size = !!sym("node_size_value")), shape = 21, colour = "black", alpha = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = phylum_palette) +
    geom_node_text(
      aes(label = !!sym("node_label_value")),
      size = 5.0 / ggplot2::.pt,
      family = "Helvetica",
      repel = TRUE,
      box.padding = grid::unit(0.2, "lines"),
      point.padding = grid::unit(0, "lines"),
      segment.size = 0.2,
      segment.color = "#4D4D4D",
      min.segment.length = 0,
      max.overlaps = Inf,
      force = 0.7,
      force_pull = 1.8,
      max.iter = 5000,
      na.rm = TRUE
    ) +
    coord_equal(clip = "off") +
    theme_void(base_size = 5.0) +
    theme(
      text = element_text(size = 5.0, family = "Helvetica"),
      plot.title = element_text(size = 5.0, family = "Helvetica", hjust = 0.5),
      plot.subtitle = element_text(size = 5.0, family = "Helvetica"),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    labs(title = plot_title, subtitle = plot_subtitle)

  cat(" [done] network created\n")

  list(
    plot = network_plot,
    description = description_text
  )
}

# Function to run network analysis for one dataset (EXACTLY as original)
run_network_dataset <- function(ds_name, sig_res_table, ps, trans, out_dir) {

  message("== Network analysis: ", nice(ds_name), " (", trans, ")")

  keep <- keep_taxa(ps, min_prevalence, trans)
  if (!length(keep)) {
    message("   no taxa left")
    return(list())
  }

  keep_names <- taxa_names(ps)[keep]
  ps_filtered <- prune_taxa(keep_names, ps)
  taxa_keep_ids <- taxa_names(ps_filtered)

  vars_config <- variable_groups[[ds_name]] %||% character(0)
  metadata_subset <- sample_data_subsets[[ds_name]]
  available_columns <- if (!is.null(metadata_subset)) setdiff(colnames(metadata_subset), "SEQN") else character(0)
  vars_to_analyze <- intersect(vars_config, available_columns)

  if (identical(ds_name, "oradWAS")) {
    vars_to_analyze <- intersect(vars_to_analyze, oradWAS_case_pass_list)
  }
  if (identical(ds_name, "outWAS")) {
    vars_to_analyze <- intersect(vars_to_analyze, outWAS_case_pass_list)
  }

  missing_vars <- setdiff(vars_config, available_columns)
  if (length(missing_vars)) {
    cat("  [warn] Missing metadata for:", paste(missing_vars, collapse = ", "), "\n")
  }

  if (!length(vars_to_analyze)) {
    message("   no variables to analyze")
    return(list())
  }

  cat("  Variables to analyze:", paste(vars_to_analyze, collapse = ", "), "\n")

  if (!is.null(sig_res_table) && nrow(sig_res_table)) {
    sig_res_table <- sig_res_table %>%
      filter(is.na(.data$otu_clean) | .data$otu_clean %in% taxa_keep_ids | clean_feature_id(.data$feature_id) %in% taxa_keep_ids)
  }

  plot_list <- list()

  for (v in vars_to_analyze) {
    network_result <- analyze_network_for_variable(
      sig_res_table = sig_res_table,
      host_var_name = v,
      ps = ps_filtered,
      trans = trans,
      min_shared = min_shared_traits,
      title_prefix = paste(nice(ds_name), trans)
    )

    if (!is.null(network_result)) {
      plot_list[[v]] <- network_result
      cat("       -> Plot stored for", v, "\n")
    }
  }

  cat("  -> Total plots created:", length(plot_list), "\n")
  plot_list
}

cat("Network analysis functions defined\n")

# =============================================================================
# MAIN EXECUTION (EXACTLY AS ORIGINAL)
# =============================================================================

cat("=== STARTING NETWORK ANALYSIS ===\n")

# Use original datasets for correct variable name matching (EXACTLY as original)
datasets_corrected <- list(
  demoWAS = "demoWAS",
  oradWAS = "oradWAS", 
  exWAS = "exWAS",
  pheWAS = "pheWAS",
  outWAS = "outWAS"
)

# Organize significant results (EXACTLY as original)
sig_results <- list(
  demoWAS = list(none = significant_results[["demoWAS_none_sig_res"]], 
                 clr = significant_results[["demoWAS_clr_sig_res"]]),
  oradWAS = list(none = significant_results[["oradWAS_none_sig_res"]], 
                 clr = significant_results[["oradWAS_clr_sig_res"]]),
  exWAS = list(none = significant_results[["exWAS_none_sig_res"]], 
               clr = significant_results[["exWAS_clr_sig_res"]]),
  pheWAS = list(none = significant_results[["pheWAS_none_sig_res"]], 
                clr = significant_results[["pheWAS_clr_sig_res"]]),
  outWAS = list(none = significant_results[["outWAS_none_sig_res"]], 
                clr = significant_results[["outWAS_clr_sig_res"]])
)

# Organize phyloseq objects (EXACTLY as original)
phylo_sets <- list(
  none = phyloseq_objects[["ubiome_relative_none"]], 
  clr = phyloseq_objects[["ubiome_relative_clr"]]
)

grand_total <- 0

# Main execution loop (EXACTLY as original)
for (tr in names(phylo_sets)) {
  ps <- phylo_sets[[tr]]
  tr_dir <- file.path(output_path, "network_plots_final", tr)
  dir.create(tr_dir, recursive = TRUE, showWarnings = FALSE)
  dataset_plot_lists <- list()
  
  cat("\nProcessing", tr, "transformation...\n")
  
  for (ds in names(datasets_corrected)) {
    if (ds %in% names(sig_results) && tr %in% names(sig_results[[ds]])) {
      
      plot_list <- run_network_dataset(
        ds_name = ds,
        sig_res_table = sig_results[[ds]][[tr]],
        ps = ps,
        trans = tr,
        out_dir = tr_dir
      )
      
      dataset_plot_lists[[ds]] <- plot_list

      if (length(plot_list) > 0) {
        cat("  Saving", length(plot_list), "plots for", ds, "\n")
        
        # Save individual plots (EXACTLY as original)
        for (var_name in names(plot_list)) {
          desc_text <- plot_list[[var_name]]$description
          desc_slug <- sanitize_for_filename(desc_text)
          filename <- file.path(tr_dir, paste0(ds, "_", var_name, "_", desc_slug, "_", tr, "_network.pdf"))
          cat("    Saving:", filename, "\n")
          
          pdf(filename, 2.5, 2.5, family = "Helvetica")
          print(plot_list[[var_name]]$plot)
          dev.off()
        }
      } else {
        cat("  No plots created for", ds, "(", tr, ")\n")
      }

      grand_total <- grand_total + length(plot_list)
    }
  }
  
  # Create booklet for this transformation (new five-page layout)
  dataset_order <- c("demoWAS", "exWAS", "outWAS", "oradWAS", "pheWAS")
  dataset_order <- dataset_order[dataset_order %in% names(dataset_plot_lists)]

  if (length(dataset_order) > 0) {
    booklet_file <- file.path(tr_dir, paste0("ALL_", tr, "_network_booklet.pdf"))
    cat("Creating booklet:", booklet_file, "with", length(dataset_order), "pages (demoWAS, exWAS, outWAS, oradWAS, pheWAS order)\n")

    max_grid_side <- max(1, sapply(dataset_order, function(ds) {
      pl <- dataset_plot_lists[[ds]]
      if (length(pl) == 0) return(1)
      ceiling(sqrt(length(pl)))
    }))

    pdf(booklet_file, width = max_grid_side * 2.5, height = max_grid_side * 2.5, family = "Helvetica")
    blank_cell <- ggplot() + theme_void()

    for (ds in dataset_order) {
      plot_list_ds <- dataset_plot_lists[[ds]]
      grid_side <- if (length(plot_list_ds)) ceiling(sqrt(length(plot_list_ds))) else 1

      if (length(plot_list_ds)) {
        grobs <- lapply(plot_list_ds, `[[`, "plot")
      } else {
        grobs <- list(
          ggplot() +
            theme_void() +
            annotate("text", x = 0.5, y = 0.5,
                     label = paste0("No significant networks for ", ds),
                     family = "Helvetica", size = 5.5 / ggplot2::.pt)
        )
      }

      num_cells <- grid_side^2
      if (length(grobs) < num_cells) {
        grobs <- c(grobs, replicate(num_cells - length(grobs), blank_cell, simplify = FALSE))
      }

      arranged <- gridExtra::arrangeGrob(grobs = grobs, ncol = grid_side)

      grid::grid.newpage()
      grid::grid.text(paste0(ds, " networks"),
                      y = grid::unit(0.98, "npc"),
                      gp = grid::gpar(fontfamily = "Helvetica", fontsize = 5.5, fontface = "bold"))

      content_vp <- grid::viewport(width = grid_side / max_grid_side,
                                   height = grid_side / max_grid_side,
                                   y = 0.45)
      grid::pushViewport(content_vp)
      grid::grid.draw(arranged)
      grid::popViewport()
    }
    dev.off()
  }
}

# =============================================================================
# FINAL SUMMARY (EXACTLY AS ORIGINAL)
# =============================================================================

message("\n════════ FINAL NETWORK SUMMARY ════════")
message("Total network plots created: ", grand_total)
message("Output directory: ", file.path(output_path, "network_plots_final"))
message("Transformations processed: ", paste(transformations, collapse = ", "))
message("Variable groups processed: ", length(datasets_corrected))
message("")
message("Files created:")
output_dir <- file.path(output_path, "network_plots_final")
if (dir.exists(output_dir)) {
  all_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  for (f in all_files) {
    message("  ", f)
  }
} else {
  message("  No output directory found")
}
message("════════════════════════════════════════")

cat("\n=== NETWORK ANALYSIS COMPLETE ===\n")