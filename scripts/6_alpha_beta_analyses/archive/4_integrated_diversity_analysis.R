#!/usr/bin/env Rscript
################################################################################
##  4. INTEGRATED ALPHA & BETA DIVERSITY ANALYSIS                             ##
##  - Combined 2×3 grid: Alpha diversity (top) + Beta centroid (bottom)       ##
##  - Violin + boxplot overlay                                                ##
##  - One PDF per categorical variable showing all 6 metrics                  ##
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ape)
  library(egg)
  library(cowplot)
  library(data.table)
})

cat("=== INTEGRATED ALPHA & BETA DIVERSITY ANALYSIS ===\n\n")

###############################################################################
## 1. CONFIGURATION                                                          ##
###############################################################################

TEST_MODE <- TRUE
TEST_SAMPLES <- 500

base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
beta_data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
output_path <- file.path(base_path, "results/analyses_results/6_alpha_beta_analyses_out")
plots_path <- file.path(output_path, "integrated_diversity_plots")

if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

cat("TEST_MODE:", TEST_MODE, "\n")
cat("Output:", plots_path, "\n\n")

grafify_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 
  "#332288", "#88ccee", "#44aa99", "#117733"
)

###############################################################################
## 2. LOAD DATA                                                              ##
###############################################################################

cat("Loading data...\n")
alpha_with_categories <- readRDS(file.path(data_path, "alpha_diversity_with_all_categories.rds"))

if (TEST_MODE) {
  alpha_with_categories <- alpha_with_categories[1:min(TEST_SAMPLES, nrow(alpha_with_categories)), ]
}

categorical_vars <- names(alpha_with_categories)[sapply(alpha_with_categories, function(x) {
  is.factor(x) && nlevels(x) >= 2 && nlevels(x) <= 10 && sum(!is.na(x)) >= 30
})]

categorical_vars <- setdiff(categorical_vars, c("sample", "cycle", "Release_Cycle"))

###############################################################################
## 3. LOAD BETA & CALCULATE CENTROIDS                                        ##
###############################################################################

cat("Loading beta diversity...\n")

load_beta_subset <- function(file_path, metric_name, sample_ids) {
  cat("  ", metric_name, "...")
  if (!file.exists(file_path)) { cat(" NOT FOUND\n"); return(NULL) }
  
  dist_data <- fread(file_path, nrows = min(length(sample_ids) + 1, 1000))
  all_seqns <- as.character(dist_data[[1]])
  keep_idx <- which(all_seqns %in% sample_ids)
  
  if (length(keep_idx) < 10) { cat(" Too few\n"); return(NULL) }
  
  dist_mat <- as.matrix(dist_data[keep_idx, (keep_idx + 1), with = FALSE])
  rownames(dist_mat) <- all_seqns[keep_idx]
  colnames(dist_mat) <- all_seqns[keep_idx]
  
  cat(" Done\n")
  return(list(dist = as.dist(dist_mat), matrix = dist_mat))
}

sample_ids <- alpha_with_categories$SEQN

beta_metrics <- list()
beta_metrics$Braycurtis <- load_beta_subset(
  file.path(beta_data_path, "dada2rsv-beta-braycurtis.txt"), "Bray-Curtis", sample_ids)
beta_metrics$Unwunifrac <- load_beta_subset(
  file.path(beta_data_path, "dada2rsv-beta-unwunifrac.txt"), "Unweighted UniFrac", sample_ids)
beta_metrics$Wunifrac <- load_beta_subset(
  file.path(beta_data_path, "dada2rsv-beta-wunifrac.txt"), "Weighted UniFrac", sample_ids)

beta_metrics <- beta_metrics[!sapply(beta_metrics, is.null)]

# Calculate centroid distances
cat("\nCalculating centroid distances...\n")

centroid_data_all <- data.frame()

for (metric_name in names(beta_metrics)) {
  dist_mat <- beta_metrics[[metric_name]]$matrix
  
  for (cat_var in categorical_vars) {
    if (!cat_var %in% colnames(alpha_with_categories)) next
    
    for (group_level in unique(alpha_with_categories[[cat_var]])) {
      if (is.na(group_level)) next
      
      group_samples <- alpha_with_categories %>%
        filter(!!sym(cat_var) == group_level) %>%
        pull(SEQN)
      
      group_samples <- intersect(group_samples, rownames(dist_mat))
      if (length(group_samples) < 2) next
      
      group_dist <- dist_mat[group_samples, group_samples]
      
      for (sample in group_samples) {
        dist_to_centroid <- mean(group_dist[sample, group_samples[group_samples != sample]])
        
        centroid_data_all <- rbind(centroid_data_all, data.frame(
          Metric = metric_name,
          Variable = cat_var,
          Group = as.character(group_level),
          Distance_to_Centroid = dist_to_centroid
        ))
      }
    }
  }
}

###############################################################################
## 4. RUN PERMANOVA                                                          ##
###############################################################################

cat("Running PERMANOVA...\n")

permanova_results <- data.frame()

for (metric_name in names(beta_metrics)) {
  dist_obj <- beta_metrics[[metric_name]]$dist
  
  common_samples <- intersect(labels(dist_obj), alpha_with_categories$SEQN)
  dist_mat_full <- beta_metrics[[metric_name]]$matrix
  dist_subset <- as.dist(dist_mat_full[common_samples, common_samples])
  
  meta_subset <- alpha_with_categories %>%
    filter(SEQN %in% common_samples) %>%
    arrange(match(SEQN, common_samples))
  
  for (cat_var in categorical_vars) {
    if (!cat_var %in% colnames(meta_subset)) next
    
    var_data <- meta_subset[[cat_var]]
    if (sum(!is.na(var_data)) < 30 || nlevels(droplevels(var_data[!is.na(var_data)])) < 2) next
    
    tryCatch({
      perm_result <- adonis2(dist_subset ~ get(cat_var), data = meta_subset, permutations = 199)
      
      permanova_results <- rbind(permanova_results, data.frame(
        Metric = metric_name,
        Variable = cat_var,
        R2 = perm_result$R2[1],
        P_value = perm_result$`Pr(>F)`[1]
      ))
    }, error = function(e) {})
  }
}

cat("  ✓ PERMANOVA:", nrow(permanova_results), "tests\n")

###############################################################################
## 5. CREATE INTEGRATED 2×3 PLOTS                                            ##
###############################################################################

cat("\nCreating integrated 2×3 plots...\n")

# Get top 10 variables
top_vars <- permanova_results %>%
  filter(P_value < 0.05) %>%
  group_by(Variable) %>%
  summarise(Avg_R2 = mean(R2), N_sig = n(), .groups = "drop") %>%
  filter(N_sig >= 2) %>%
  arrange(desc(Avg_R2)) %>%
  slice_head(n = 10) %>%
  pull(Variable)

cat("  Creating plots for", length(top_vars), "variables\n\n")

for (cat_var in top_vars) {
  
  plot_data_alpha <- alpha_with_categories %>% filter(!is.na(!!sym(cat_var)))
  plot_data_beta <- centroid_data_all %>% filter(Variable == cat_var)
  
  if (nrow(plot_data_alpha) < 30) next
  
  n_levels <- nlevels(plot_data_alpha[[cat_var]])
  plot_colors <- grafify_colors[1:min(n_levels, length(grafify_colors))]
  
  plot_list <- list()
  
  # ALPHA PLOTS (top row)
  for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
    metric_label <- gsub("_", " ", metric_col)
    
    p <- ggplot(plot_data_alpha, aes(x = !!sym(cat_var), y = !!sym(metric_col), fill = !!sym(cat_var))) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
      geom_jitter(alpha = 0.3, size = 0.2, width = 0.15) +
      scale_fill_manual(values = plot_colors) +
      labs(title = metric_label, x = "", y = "") +
      egg::theme_article(base_size = 5) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5),
        axis.text.y = element_text(size = 3.5)
      )
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # BETA PLOTS (bottom row)
  for (metric_name in names(beta_metrics)) {
    beta_plot_data <- plot_data_beta %>% filter(Metric == metric_name)
    
    if (nrow(beta_plot_data) < 30) {
      plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
      next
    }
    
    perm_row <- permanova_results %>% filter(Metric == metric_name, Variable == cat_var)
    subtitle <- if(nrow(perm_row) > 0) {
      paste0("R²=", round(perm_row$R2, 3), ", P=", 
             if(perm_row$P_value < 0.001) "<0.001" else round(perm_row$P_value, 3))
    } else ""
    
    p <- ggplot(beta_plot_data, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
      geom_jitter(alpha = 0.3, size = 0.2, width = 0.15) +
      scale_fill_manual(values = plot_colors) +
      labs(title = paste(metric_name, "Centroid"), subtitle = subtitle, x = "", y = "") +
      egg::theme_article(base_size = 5) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 3.5, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5),
        axis.text.y = element_text(size = 3.5)
      )
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # Combine
  combined <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 2,
                       labels = c("A", "B", "C", "D", "E", "F"), label_size = 6)
  
  filename <- file.path(plots_path, paste0("Integrated_diversity_", cat_var, ".pdf"))
  ggsave(filename, combined, width = 9, height = 6, dpi = 300)
  
  cat("  ✓ Integrated_diversity_", cat_var, ".pdf\n", sep = "")
}

cat("\n✓ Complete! Integrated 2×3 plots created\n")
cat("  Output:", plots_path, "\n")
cat("========================================\n")
