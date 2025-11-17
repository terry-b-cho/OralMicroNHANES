#!/usr/bin/env Rscript
################################################################################
##  3. COMPREHENSIVE BETA DIVERSITY ANALYSIS - ALL CATEGORICAL VARIABLES      ##
##  - ONE PDF per categorical variable with 3 metrics side-by-side            ##
##  - Plus sign (+) centroids, violin plots for centroid distances            ##
##  - Combined scree plots and combined centroid plots                        ##
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

cat("=== BETA DIVERSITY ANALYSIS - ALL CATEGORICAL VARIABLES ===\n\n")

###############################################################################
## 1. CONFIGURATION                                                          ##
###############################################################################

TEST_MODE <- TRUE
TEST_SAMPLES <- 500

base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
beta_data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
output_path <- file.path(base_path, "results/analyses_results/6_alpha_beta_analyses_out")
plots_path <- file.path(output_path, "beta_diversity_plots")

if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

cat("TEST_MODE:", TEST_MODE, "\n")
if (TEST_MODE) cat("TEST_SAMPLES:", TEST_SAMPLES, "\n\n")

grafify_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 
  "#332288", "#88ccee", "#44aa99", "#117733", "#999933", "#ddcc77", "#cc6677"
)

###############################################################################
## 2. LOAD DATA                                                              ##
###############################################################################

cat("Loading data...\n")
alpha_with_categories <- readRDS(file.path(data_path, "alpha_diversity_with_all_categories.rds"))

if (TEST_MODE) {
  alpha_with_categories <- alpha_with_categories[1:min(TEST_SAMPLES, nrow(alpha_with_categories)), ]
}

cat("  Data loaded:", nrow(alpha_with_categories), "samples\n")

categorical_vars <- names(alpha_with_categories)[sapply(alpha_with_categories, function(x) {
  is.factor(x) && nlevels(x) >= 2 && nlevels(x) <= 15 && sum(!is.na(x)) >= 30
})]

categorical_vars <- setdiff(categorical_vars, c("sample", "cycle", "Release_Cycle"))

cat("  Categorical variables:", length(categorical_vars), "\n\n")

###############################################################################
## 3. LOAD BETA DIVERSITY MATRICES                                           ##
###############################################################################

cat("Loading beta diversity matrices...\n")

load_beta_subset <- function(file_path, metric_name, sample_ids) {
  cat("  Loading", metric_name, "...")
  if (!file.exists(file_path)) { cat(" NOT FOUND\n"); return(NULL) }
  
  dist_data <- fread(file_path, nrows = min(length(sample_ids) + 1, 1000))
  all_seqns <- as.character(dist_data[[1]])
  keep_idx <- which(all_seqns %in% sample_ids)
  
  if (length(keep_idx) < 10) { cat(" Too few\n"); return(NULL) }
  
  dist_mat <- as.matrix(dist_data[keep_idx, (keep_idx + 1), with = FALSE])
  rownames(dist_mat) <- all_seqns[keep_idx]
  colnames(dist_mat) <- all_seqns[keep_idx]
  
  cat(" Done (", length(keep_idx), ")\n")
  return(as.dist(dist_mat))
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
cat("\n")

###############################################################################
## 4. PERFORM PCoA                                                           ##
###############################################################################

cat("Performing PCoA...\n")

pcoa_results <- list()

for (metric_name in names(beta_metrics)) {
  pcoa_res <- pcoa(beta_metrics[[metric_name]])
  
  pcoa_df <- data.frame(
    SEQN = rownames(pcoa_res$vectors),
    PC1 = pcoa_res$vectors[, 1],
    PC2 = pcoa_res$vectors[, 2]
  ) %>% left_join(alpha_with_categories, by = "SEQN")
  
  var_exp <- pcoa_res$values$Eigenvalues / sum(pcoa_res$values$Eigenvalues) * 100
  
  pcoa_results[[metric_name]] <- list(data = pcoa_df, var_explained = var_exp)
  cat("  ", metric_name, ": PC1=", round(var_exp[1], 1), "%, PC2=", round(var_exp[2], 1), "%\n")
}

###############################################################################
## 5. PERMANOVA FOR ALL VARIABLES                                            ##
###############################################################################

cat("\nRunning PERMANOVA...\n")

permanova_results <- data.frame()

for (metric_name in names(beta_metrics)) {
  dist_obj <- beta_metrics[[metric_name]]
  pcoa_data <- pcoa_results[[metric_name]]$data
  
  common_samples <- intersect(labels(dist_obj), pcoa_data$SEQN)
  dist_mat <- as.matrix(dist_obj)
  dist_subset <- as.dist(dist_mat[common_samples, common_samples])
  
  meta_subset <- pcoa_data %>%
    filter(SEQN %in% common_samples) %>%
    arrange(match(SEQN, common_samples))
  
  for (cat_var in categorical_vars) {
    if (!cat_var %in% colnames(meta_subset)) next
    
    var_data <- meta_subset[[cat_var]]
    if (sum(!is.na(var_data)) < 30) next
    if (nlevels(droplevels(var_data[!is.na(var_data)])) < 2) next
    
    tryCatch({
      perm_result <- adonis2(dist_subset ~ get(cat_var), data = meta_subset, permutations = 199)
      
      permanova_results <- rbind(permanova_results, data.frame(
        Metric = metric_name,
        Variable = cat_var,
        R2 = perm_result$R2[1],
        P_value = perm_result$`Pr(>F)`[1],
        N_samples = nrow(meta_subset)
      ))
    }, error = function(e) {})
  }
}

write.csv(permanova_results, file.path(output_path, "permanova_all_variables.csv"), row.names = FALSE)
cat("  ✓ Completed:", nrow(permanova_results), "tests,", sum(permanova_results$P_value < 0.05), "significant\n")

###############################################################################
## 6. CALCULATE CENTROID DISTANCES                                           ##
###############################################################################

cat("\nCalculating centroid distances...\n")

centroid_distances <- data.frame()

for (metric_name in names(beta_metrics)) {
  dist_mat <- as.matrix(beta_metrics[[metric_name]])
  pcoa_data <- pcoa_results[[metric_name]]$data
  
  for (cat_var in categorical_vars) {
    if (!cat_var %in% colnames(pcoa_data)) next
    
    for (group_level in unique(pcoa_data[[cat_var]])) {
      if (is.na(group_level)) next
      
      group_samples <- pcoa_data %>%
        filter(!!sym(cat_var) == group_level) %>%
        pull(SEQN)
      
      group_samples <- intersect(group_samples, rownames(dist_mat))
      if (length(group_samples) < 2) next
      
      group_dist <- dist_mat[group_samples, group_samples]
      
      for (sample in group_samples) {
        dist_to_centroid <- mean(group_dist[sample, group_samples[group_samples != sample]])
        
        centroid_distances <- rbind(centroid_distances, data.frame(
          Metric = metric_name,
          Variable = cat_var,
          Group = as.character(group_level),
          Distance_to_Centroid = dist_to_centroid
        ))
      }
    }
  }
}

write.csv(centroid_distances, file.path(output_path, "beta_diversity_centroid_distances_all.csv"), row.names = FALSE)
cat("  ✓ Calculated:", nrow(centroid_distances), "centroid distances\n")

###############################################################################
## 7. CREATE FACETED PCoA PLOTS (ONE PDF PER CATEGORICAL VARIABLE)           ##
###############################################################################

cat("\nCreating faceted PCoA plots (one PDF per categorical variable)...\n")

# Get variables with significant effects in at least one metric
sig_variables <- permanova_results %>%
  filter(P_value < 0.05) %>%
  pull(Variable) %>%
  unique()

cat("  Creating plots for", length(sig_variables), "significant variables\n\n")

plot_counter <- 0

for (cat_var in sig_variables) {
  
  plot_list <- list()
  
  for (metric_name in names(beta_metrics)) {
    
    pcoa_data <- pcoa_results[[metric_name]]$data %>% filter(!is.na(!!sym(cat_var)))
    if (nrow(pcoa_data) < 30) next
    
    var_exp <- pcoa_results[[metric_name]]$var_explained
    
    # Get PERMANOVA results for annotation
    perm_row <- permanova_results %>% 
      filter(Metric == metric_name, Variable == cat_var)
    
    if (nrow(perm_row) == 0) next
    
    anno_text <- paste0("R² = ", round(perm_row$R2, 3), ", P = ", 
                       if(perm_row$P_value < 0.001) "< 0.001" else round(perm_row$P_value, 3))
    
    # Calculate centroids
    centroids <- pcoa_data %>%
      group_by(!!sym(cat_var)) %>%
      summarise(PC1_center = mean(PC1), PC2_center = mean(PC2), .groups = "drop")
    
    n_levels <- nlevels(pcoa_data[[cat_var]])
    plot_colors <- grafify_colors[1:min(n_levels, length(grafify_colors))]
    
    # Create plot with PLUS SIGN centroids (shape = 3)
    p <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = !!sym(cat_var), fill = !!sym(cat_var))) +
      geom_point(alpha = 0.5, size = 0.6) +
      stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.4, level = 0.95) +
      geom_point(data = centroids, aes(x = PC1_center, y = PC2_center),
                 size = 3, shape = 3, stroke = 1.5, show.legend = FALSE) +  # PLUS SIGN
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      labs(
        title = metric_name,
        subtitle = anno_text,
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
        color = gsub("_", " ", cat_var)
      ) +
      egg::theme_article(base_size = 5) +
      theme(
        plot.title = element_text(size = 5, face = "bold"),
        plot.subtitle = element_text(size = 4),
        axis.title = element_text(size = 4),
        axis.text = element_text(size = 3.5),
        legend.title = element_text(size = 4),
        legend.text = element_text(size = 3.5),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm")
      )
    
    plot_list[[metric_name]] <- p
  }
  
  if (length(plot_list) >= 2) {
    combined <- plot_grid(plotlist = plot_list, ncol = 3, labels = LETTERS[1:length(plot_list)], label_size = 6)
    
    filename <- file.path(plots_path, paste0("PCoA_", cat_var, "_all_metrics.pdf"))
    ggsave(filename, combined, width = 4.5, height = 1.5, dpi = 300)  # 50% smaller
    
    plot_counter <- plot_counter + 1
    cat("  ✓", plot_counter, ". PCoA_", cat_var, "_all_metrics.pdf\n", sep = "")
  }
}

cat("\n  Total faceted PCoA plots:", plot_counter, "\n")

###############################################################################
## 8. COMBINED SCREE PLOTS (ALL 3 METRICS SIDE-BY-SIDE)                      ##
###############################################################################

cat("\nCreating combined scree plot...\n")

scree_data_all <- data.frame()

for (metric_name in names(beta_metrics)) {
  var_exp <- pcoa_results[[metric_name]]$var_explained[1:10]
  
  scree_data <- data.frame(
    Metric = metric_name,
    Axis = paste0("PC", 1:length(var_exp)),
    Variance = var_exp
  )
  
  scree_data_all <- rbind(scree_data_all, scree_data)
}

scree_data_all$Axis <- factor(scree_data_all$Axis, levels = unique(scree_data_all$Axis))
scree_data_all$Metric <- factor(scree_data_all$Metric, levels = names(beta_metrics))

p_scree_combined <- ggplot(scree_data_all, aes(x = Axis, y = Variance)) +
  geom_bar(stat = "identity", fill = "#56B4E9", alpha = 0.7) +
  geom_line(aes(group = 1), color = "#E69F00", linewidth = 0.6) +
  geom_point(color = "#E69F00", size = 1.5) +
  facet_wrap(~ Metric, ncol = 3) +
  labs(
    title = "Variance Explained by Principal Coordinates",
    x = "Principal Coordinate",
    y = "Variance Explained (%)"
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    plot.title = element_text(size = 5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5),
    strip.text = element_text(size = 5, face = "bold")
  )

ggsave(file.path(plots_path, "Scree_all_metrics_combined.pdf"), p_scree_combined, 
       width = 4.5, height = 1.5, dpi = 300)  # 50% smaller
cat("  ✓ Saved: Scree_all_metrics_combined.pdf\n")

###############################################################################
## 9. COMBINED CENTROID DISTANCE PLOTS (VIOLIN + POINTS)                     ##
###############################################################################

cat("\nCreating combined centroid distance plots...\n")

# Get top 5 variables by average R² across metrics
top_centroid_vars <- permanova_results %>%
  filter(P_value < 0.05) %>%
  group_by(Variable) %>%
  summarise(Avg_R2 = mean(R2), .groups = "drop") %>%
  arrange(desc(Avg_R2)) %>%
  slice_head(n = 5) %>%
  pull(Variable)

cat("  Creating centroid plots for top", length(top_centroid_vars), "variables\n")

for (cat_var in top_centroid_vars) {
  
  plot_data <- centroid_distances %>%
    filter(Variable == cat_var)
  
  if (nrow(plot_data) == 0) next
  
  plot_data$Metric <- factor(plot_data$Metric, levels = names(beta_metrics))
  
  # VIOLIN PLOT with BOXPLOT overlay and individual points (alpha = 0.3)
  p_centroid <- ggplot(plot_data, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
    geom_violin(alpha = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA, 
                 position = position_dodge(width = 0.9)) +  # Boxplot overlay showing IQR
    geom_jitter(alpha = 0.3, size = 0.3, width = 0.15) +  # Individual points with alpha 0.3
    facet_wrap(~ Metric, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = grafify_colors) +
    labs(
      title = paste("Distance to Centroid −", gsub("_", " ", cat_var)),
      x = gsub("_", " ", cat_var),
      y = "Distance to Centroid"
    ) +
    egg::theme_article(base_size = 5) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5),
      axis.text.y = element_text(size = 3.5),
      strip.text = element_text(size = 4, face = "bold")
    )
  
  filename <- file.path(plots_path, paste0("Centroid_distance_", cat_var, "_all_metrics.pdf"))
  ggsave(filename, p_centroid, width = 4.5, height = 1.5, dpi = 300)  # 50% smaller
  cat("  ✓ Saved: Centroid_distance_", cat_var, "_all_metrics.pdf\n", sep = "")
}

###############################################################################
## 10. SUMMARY                                                               ##
###############################################################################

cat("\n=== FINAL SUMMARY ===\n")
cat("Mode:", ifelse(TEST_MODE, "TEST", "PRODUCTION"), "\n")
cat("Samples:", nrow(alpha_with_categories), "\n")
cat("Categorical variables:", length(categorical_vars), "\n")
cat("Beta diversity metrics:", length(beta_metrics), "\n")
cat("PERMANOVA tests:", nrow(permanova_results), "\n")
cat("Significant (P<0.05):", sum(permanova_results$P_value < 0.05), "\n")
cat("\nPlots created:\n")
cat("  - Faceted PCoA plots:", plot_counter, "(one per variable, 3 metrics side-by-side)\n")
cat("  - Combined scree plot: 1\n")
cat("  - Combined centroid plots:", length(top_centroid_vars), "\n")

cat("\nTop 5 associations:\n")
top5 <- permanova_results %>% arrange(desc(R2)) %>% head(5)
print(top5[, c("Metric", "Variable", "R2", "P_value")])

cat("\n✓ All outputs saved to:", output_path, "\n")
cat("========================================\n")
