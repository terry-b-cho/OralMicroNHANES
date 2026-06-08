#!/usr/bin/env Rscript

# =============================================================================
# GROUPED NETWORK ANALYSIS SCRIPT
# =============================================================================
# Variant of 2_network_additional.R that groups host variables before
# constructing the network. Reads the intermediates produced by
# 1_load_all_categories_additional.R. Emits comprehensive/ and compact/
# layouts using the same K-per-direction logic as 2_network_additional.R.
#
# All outputs are written under:
#   results/analyses_results/8_network_analyses_out_additional/
#
# Environment: R >= 4.5 with the dependencies declared inside the script.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

cat("=== GROUPED NETWORK ANALYSIS ===\n")

library(phyloseq)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(grid)
library(gridExtra)
library(extrafont)
library(scales)
library(rlang)

set.seed(42)

loadfonts(device = "pdf", quiet = TRUE)

# =============================================================================
# PATHS / DATA
# =============================================================================

base_path <- PROJECT_ROOT
data_path <- file.path(base_path, "results/analyses_results/8_network_analyses_out_additional/intermediate")
output_path <- file.path(base_path, "results/analyses_results/8_network_analyses_out_additional")
pdf_font_family <- Sys.getenv("NHANES_NETWORK_PDF_FONT", unset = "Helvetica")
COMPACT_TOP_K_PER_DIRECTION <- 15L

if (!dir.exists(data_path)) {
  stop("Intermediate data not found. Run 1_load_all_categories_additional.R first.\nMissing: ", data_path)
}

cat("Data directory:", data_path, "\n")
cat("Output directory:", output_path, "\n\n")

phyloseq_objects <- readRDS(file.path(data_path, "phyloseq_objects.rds"))
significant_results <- readRDS(file.path(data_path, "significant_results.rds"))
variable_description <- readRDS(file.path(data_path, "variable_description.rds"))
sample_data_subsets <- readRDS(file.path(data_path, "sample_data_subsets.rds"))
ubiome_genus_mapping_complete <- readRDS(file.path(data_path, "ubiome_genus_mapping_complete.rds"))
otu_pass_prevalence_list <- readRDS(file.path(data_path, "otu_pass_prevalence_list.rds"))

# Optional case counts (used for outWAS grouping summaries)
outwas_case_counts <- NULL
case_counts_path <- file.path(data_path, "outWAS_case_counts.rds")
if (file.exists(case_counts_path)) {
  outwas_case_counts <- readRDS(case_counts_path)
}

# Map transformation to phyloseq object
phylo_sets <- list(
  none = phyloseq_objects[["ubiome_relative_none"]],
  clr = phyloseq_objects[["ubiome_relative_clr"]]
)

min_shared_traits <- 2
min_prevalence <- 0.10

# =============================================================================
# COLOR PALETTES / HELPERS (copied from 2_network.R to keep styling identical)
# =============================================================================

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
  "unclassified"         = "#DDDDDD",  # pale_grey
  "NA"                   = "#DDDDDD",  # very light grey
  "TM7"                  = "#000000"   # Contrast Black
)

orbl_colors <- c("#2B5B8A", "#4071A0", "#5789B6", "#6FA3CB", "#A2BCCF",
                 "#D8D4C9", "#F0AC72", "#EF8530", "#DA6524", "#BE4E21", "#9E3D21")

map_to_orbl_color <- function(pct_positive) {
  idx <- pmax(1, pmin(length(orbl_colors), floor(pct_positive * (length(orbl_colors) - 1)) + 1))
  orbl_colors[idx]
}

tax_tbl <- function(ps) {
  df <- as.data.frame(unclass(tax_table(ps)), stringsAsFactors = FALSE)
  df$OTU <- rownames(df)
  df
}

make_phylum_col <- function(v, alpha = 0.65) {
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

clean_feature_id <- function(x) stringr::str_remove(x, "_relative$")

`%||%` <- function(x, y) if (!is.null(x)) x else y

suppress_compact_annotation <- function(raw) {
  if (length(raw) != 1L) {
    return(TRUE)
  }
  if (is.na(raw)) {
    return(TRUE)
  }
  raw <- trimws(as.character(raw))
  if (!nzchar(raw)) {
    return(TRUE)
  }
  gl <- tolower(raw)
  if (gl %in% c("na", "unclassified")) {
    return(TRUE)
  }
  grepl("incertae", gl, fixed = TRUE) || grepl("ncertae", gl, fixed = TRUE) || grepl("sedis", gl, fixed = TRUE)
}

clean_taxon_label_compact <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(NA_character_)
  }
  x2 <- gsub("\\[|\\]", "", as.character(x))
  sub("_.*$", "", x2)
}

prepare_compact_sig_subset <- function(sig_res_table, host_var_name, taxa_keep_ids, k_per_dir) {
  if (is.null(sig_res_table) || !nrow(sig_res_table)) {
    return(list(tbl = sig_res_table, summary = tibble::tibble()))
  }
  d <- sig_res_table %>%
    filter(.data$dependent_var == host_var_name | .data$independent_var == host_var_name) %>%
    mutate(
      feature_id = if_else(is.na(.data$feature_id), .data$phenotype, .data$feature_id),
      otu = if_else(is.na(.data$otu_clean) | .data$otu_clean == "", clean_feature_id(.data$feature_id), .data$otu_clean)
    ) %>%
    filter(!is.na(.data$otu)) %>%
    filter(is.na(.data$otu_clean) | .data$otu_clean %in% taxa_keep_ids | clean_feature_id(.data$feature_id) %in% taxa_keep_ids)

  if (!nrow(d)) {
    return(list(tbl = d, summary = tibble::tibble()))
  }

  d_dedup <- d %>%
    group_by(.data$otu) %>%
    slice_max(order_by = abs(.data$estimate), n = 1L, with_ties = FALSE) %>%
    ungroup()

  pos <- d_dedup %>%
    filter(.data$estimate > 0) %>%
    arrange(desc(abs(.data$estimate))) %>%
    slice_head(n = k_per_dir)
  neg <- d_dedup %>%
    filter(.data$estimate < 0) %>%
    arrange(desc(abs(.data$estimate))) %>%
    slice_head(n = k_per_dir)

  keep_otus <- unique(c(pos$otu, neg$otu))
  if (!length(keep_otus)) {
    return(list(tbl = d[0L, , drop = FALSE], summary = tibble::tibble()))
  }

  sub_tbl <- d %>% filter(.data$otu %in% keep_otus)
  summ <- dplyr::bind_rows(
    pos %>% mutate(direction = "positive", rank_within_direction = dplyr::row_number()),
    neg %>% mutate(direction = "negative", rank_within_direction = dplyr::row_number())
  )
  list(tbl = sub_tbl, summary = summ)
}

keep_taxa <- function(ps, min_prev = 0.10, trans = "none") {
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)

  ok1 <- if (trans %in% c("clr", "lognorm"))
    apply(m, 1, var, na.rm = TRUE) > 0
  else
    rowSums(m > 0) / ncol(m) >= min_prev

  ok2 <- !is.na(tax_tbl(ps)$Genus)
  which(ok1 & ok2)
}

analyze_network_for_variable <- function(sig_res_table, host_var_name, ps, trans,
                                         min_shared = 2, title_prefix = "",
                                         annotation_mode = "comprehensive") {

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
    filter(!is.na(.data$otu)) %>%
    group_by(.data$otu) %>%
    slice_max(order_by = abs(.data$effect), n = 1L, with_ties = FALSE) %>%
    ungroup()

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
  if (identical(annotation_mode, "compact")) {
    node_metadata$node_label_value <- vapply(node_metadata$genus_label, function(g) {
      rg <- as.character(g)
      if (suppress_compact_annotation(rg)) {
        return(NA_character_)
      }
      clean_taxon_label_compact(rg)
    }, character(1), USE.NAMES = FALSE)
  } else {
    node_metadata$node_label_value <- dplyr::case_when(
      is.na(node_metadata$genus_label) ~ NA_character_,
      tolower(node_metadata$genus_label) %in% c("unclassified", "na") ~ NA_character_,
      TRUE ~ format_genus_label(node_metadata$genus_label)
    )
  }

  phylum_palette <- make_phylum_col(node_metadata$node_phylum)
  node_metadata$node_color <- phylum_palette[node_metadata$node_phylum]

  graph_df <- edge_list %>%
    transmute(from = .data$otu1, to = .data$otu2, layout_weight = .data$layout_weight,
              edge_width_value = .data$visual_weight, edge_colour_value = .data$edge_color)

  # set.seed(1)
  set.seed(42)
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

  # set.seed(1)
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
      family = pdf_font_family,
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
      text = element_text(size = 5.0, family = pdf_font_family),
      plot.title = element_text(size = 5.0, family = pdf_font_family, hjust = 0.5),
      plot.subtitle = element_text(size = 5.0, family = pdf_font_family),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    labs(title = plot_title, subtitle = plot_subtitle)

  cat(" [done] network created\n")

  list(
    plot = network_plot,
    description = description_text
  )
}

variable_description_map <- variable_description %>%
  filter(!is.na(.data$Variable.Name)) %>%
  mutate(
    Variable.Name = as.character(.data$Variable.Name),
    Variable.Description = str_squish(as.character(.data$Variable.Description))
  ) %>%
  filter(!is.na(.data$Variable.Description), .data$Variable.Description != "") %>%
  arrange(.data$Variable.Name) %>%
  distinct(.data$Variable.Name, .keep_all = TRUE) %>%
  {
    setNames(.$Variable.Description, .$Variable.Name)
  }

get_variable_description <- function(var_name) {
  desc <- variable_description_map[[var_name]]
  if (!is.null(desc) && !is.na(desc) && nzchar(desc)) {
    desc
  } else {
    var_name
  }
}

resolve_host_variable_name <- function(entry, host_vars_available, dataset) {
  if (length(host_vars_available) == 0) {
    return(entry)
  }

  entry_clean <- entry
  prefix <- paste0(dataset, "_")
  if (startsWith(entry_clean, prefix)) {
    entry_clean <- substr(entry_clean, nchar(prefix) + 1, nchar(entry_clean))
  }

  if (entry_clean %in% host_vars_available) {
    return(entry_clean)
  }

  for (hv in host_vars_available) {
    if (startsWith(entry_clean, hv)) {
      return(hv)
    }
  }

  tokens <- strsplit(entry_clean, "_")[[1]]
  if (length(tokens) > 1) {
    combos <- vapply(seq_along(tokens), function(k) paste(tokens[seq_len(k)], collapse = "_"), character(1))
    for (combo in combos) {
      if (combo %in% host_vars_available) {
        return(combo)
      }
    }
    combos_rev <- vapply(seq_along(tokens), function(k) paste(tokens[k:length(tokens)], collapse = "_"), character(1))
    for (combo in combos_rev) {
      if (combo %in% host_vars_available) {
        return(combo)
      }
    }
  }

  entry_clean
}

# =============================================================================
# GROUP DEFINITIONS (extendable)
# =============================================================================

# To add another umbrella network, append a new list element matching this
# structure (dataset, variables, transformations, optional label and notes).
group_specs <- list(
  Respiratory_diseases = list(
    dataset = "outWAS",
    variables = c(
      "ASTHMA_ASTHMA_MCQ",
      "EMPHYSEMA_EMPHYSEMA_MCQ",
      "BRONCHITIS_BRONCHITIS_MCQ"
    ),
    transformations = c("none", "clr"),
    label = "Respiratory diseases",
    subtitle_note = "Aggregated ASTHMA, EMPHYSEMA, BRONCHITIS"
  )
)

# =============================================================================
# CORE LOGIC
# =============================================================================

build_group_network <- function(group_name, spec, transformation, network_mode = "comprehensive") {
  dataset <- spec$dataset
  variables <- spec$variables
  group_label <- spec$label %||% group_name

  sig_key <- paste0(dataset, "_", transformation, "_sig_res")
  if (!sig_key %in% names(significant_results)) {
    message("[FAIL] No significant results found for ", sig_key)
    return(NULL)
  }

  sig_res_table <- significant_results[[sig_key]]
  if (is.null(sig_res_table) || !nrow(sig_res_table)) {
    message("[FAIL] Empty significant results for ", sig_key)
    return(NULL)
  }

  host_vars_available <- unique(sig_res_table$host_variable_column)
  resolved_vars <- unique(vapply(variables, resolve_host_variable_name, character(1), host_vars_available = host_vars_available, dataset = dataset))
  matched_vars <- intersect(resolved_vars, host_vars_available)
  unmatched_vars <- setdiff(resolved_vars, matched_vars)

  if (length(unmatched_vars)) {
    message("  [warn] Host variables not found for ", group_name, " (", transformation, "): ", paste(unmatched_vars, collapse = ", "))
  }

  filtered <- sig_res_table %>%
    filter(.data$host_variable_column %in% matched_vars)

  if (!nrow(filtered)) {
    message("[FAIL] No overlapping associations for ", group_name, " (", transformation, ")")
    return(NULL)
  }

  agg_id <- paste0(dataset, "_", group_name)

  filtered_aug <- filtered %>%
    mutate(host_var_orig = .data$host_variable_column) %>%
    mutate(
      dependent_var = if_else(.data$host_var_orig == .data$dependent_var, agg_id, .data$dependent_var),
      independent_var = if_else(.data$host_var_orig == .data$independent_var, agg_id, .data$independent_var),
      host_variable_column = agg_id
    ) %>%
    select(-.data$host_var_orig)

  if (!nrow(filtered_aug)) {
    message("[FAIL] Aggregated associations empty for ", group_name, " (", transformation, ")")
    return(NULL)
  }

  ps <- phylo_sets[[transformation]]
  keep <- keep_taxa(ps, min_prevalence, transformation)
  if (!length(keep)) {
    message("[FAIL] No taxa after prevalence filter for ", group_name, " (", transformation, ")")
    return(NULL)
  }
  keep_names <- taxa_names(ps)[keep]
  ps_filtered <- prune_taxa(keep_names, ps)

  filtered_aug <- filtered_aug %>%
    filter(is.na(.data$otu_clean) | .data$otu_clean %in% keep_names | clean_feature_id(.data$feature_id) %in% keep_names)

  if (!nrow(filtered_aug)) {
    message("[FAIL] No associations remain after taxa filter for ", group_name, " (", transformation, ")")
    return(NULL)
  }

  if (identical(network_mode, "compact")) {
    prep <- prepare_compact_sig_subset(filtered_aug, agg_id, keep_names, COMPACT_TOP_K_PER_DIRECTION)
    filtered_aug <- prep$tbl
    if (!nrow(filtered_aug)) {
      message("[FAIL] No associations after compact filter for ", group_name, " (", transformation, ")")
      return(NULL)
    }
  }

  variable_description_map[agg_id] <<- group_label

  ann_mode <- if (identical(network_mode, "compact")) "compact" else "comprehensive"

  result <- analyze_network_for_variable(
    sig_res_table = filtered_aug,
    host_var_name = agg_id,
    ps = ps_filtered,
    trans = transformation,
    min_shared = min_shared_traits,
    title_prefix = paste(group_label, transformation),
    annotation_mode = ann_mode
  )

  if (is.null(result)) {
    return(NULL)
  }

  var_descs <- purrr::map_chr(matched_vars, get_variable_description)
  vars_label <- paste(matched_vars, collapse = ", ")
  desc_label <- paste(var_descs, collapse = "; ")

  total_cases <- NA
  if (!is.null(outwas_case_counts) && dataset == "outWAS") {
    total_cases <- outwas_case_counts %>%
      filter(.data$var_name %in% matched_vars) %>%
      summarise(total_cases = sum(.data$cases_count, na.rm = TRUE)) %>%
      pull(.data$total_cases)
    if (length(total_cases) == 0) total_cases <- NA
  }

  subtitle_suffix <- if (!is.na(total_cases)) sprintf(" cases≈%s", format(total_cases, big.mark = ",")) else ""

  adjusted_plot <- result$plot +
    labs(
      title = sprintf("%s — %s (%s, %s)", group_label, dataset, transformation, network_mode),
      subtitle = paste0(result$plot$labels$subtitle, subtitle_suffix)
    )

  list(
    plot = adjusted_plot,
    variables = vars_label,
    descriptions = desc_label
  )
}

# =============================================================================
# EXECUTION
# =============================================================================

saved_files <- list()

for (network_mode in c("comprehensive", "compact")) {
  group_dir_base <- file.path(output_path, "group_networks", network_mode)
  dir.create(group_dir_base, recursive = TRUE, showWarnings = FALSE)

  for (group_name in names(group_specs)) {
    spec <- group_specs[[group_name]]
    transformations <- spec$transformations %||% names(phylo_sets)

    for (trans in transformations) {
      message("Processing group: ", group_name, " (", trans, ", ", network_mode, ")")
      result <- build_group_network(group_name, spec, trans, network_mode = network_mode)
      if (is.null(result)) {
        next
      }

      group_dir <- file.path(group_dir_base, trans)
      dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

      filename <- file.path(
        group_dir,
        paste0(group_name, "_", trans, "_network.pdf")
      )

      grDevices::pdf(filename, width = 2.5, height = 2.5, family = pdf_font_family)
      print(result$plot)
      grDevices::dev.off()

      saved_files[[length(saved_files) + 1]] <- filename
      message("  [done] Saved ", filename)
    }
  }
}

cat("\n════════ GROUP NETWORK SUMMARY ════════\n")
if (length(saved_files)) {
  cat("Files created:\n")
  purrr::walk(saved_files, ~ cat("  ", ., "\n", sep = ""))
} else {
  cat("No group networks were generated.\n")
}
cat("════════════════════════════════════\n")

cat("\n=== GROUPED NETWORK ANALYSIS COMPLETE ===\n")

