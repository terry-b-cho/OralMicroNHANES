#!/usr/bin/env Rscript
################################################################################
##  3. COMPREHENSIVE BETA DIVERSITY ANALYSIS WITH PCoA                        ##
##  - Principal Coordinate Analysis (PCoA) for all beta diversity metrics     ##
##  - Plots with 95% confidence ellipses and centroids                        ##
##  - Marginal distributions on both axes                                     ##
##  - PERMANOVA tests (adonis2) with R² and p-values                          ##
##  - Distance to centroid analysis                                           ##
##  - TEST MODE: Use smaller subset for testing                               ##
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ape)
  library(egg)
  library(cowplot)
  library(readr)
  library(data.table)
})

cat("=== BETA DIVERSITY ANALYSIS WITH PCoA ===\n\n")

###############################################################################
## 1. CONFIGURATION AND PATHS                                                ##
###############################################################################

base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
beta_data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
output_path <- file.path(base_path, "results/analyses_results/6_alpha_beta_analyses_out")
plots_path <- file.path(output_path, "beta_diversity_plots")

# Create directories
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

# TEST MODE CONFIGURATION
TEST_MODE <- TRUE  # Set to FALSE for full dataset
TEST_SAMPLES <- 500  # Number of samples to use in test mode

cat("Configuration:\n")
cat("  TEST_MODE:", TEST_MODE, "\n")
if (TEST_MODE) cat("  TEST_SAMPLES:", TEST_SAMPLES, "\n")
cat("  Output path:", output_path, "\n")
cat("  Plots path:", plots_path, "\n\n")

# Grafify colors
grafify_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 
  "#332288", "#88ccee", "#44aa99", "#117733"
)

###############################################################################
## 2. LOAD DEMOGRAPHIC DATA                                                  ##
###############################################################################

cat("Loading demographic data...\n")
demo_data <- readRDS(file.path(data_path, "demographic_data.rds"))
cat("  Demographic data loaded:", nrow(demo_data), "samples\n")

# Apply test mode subset if needed
if (TEST_MODE) {
  demo_data <- demo_data[1:min(TEST_SAMPLES, nrow(demo_data)), ]
  cat("  TEST MODE: Using", nrow(demo_data), "samples\n")
}

###############################################################################
## 3. LOAD BETA DIVERSITY DISTANCE MATRICES                                  ##
###############################################################################

cat("\nLoading beta diversity distance matrices...\n")

# Function to load and subset distance matrix
load_beta_distance <- function(file_path, metric_name, sample_ids, test_mode = FALSE) {
  cat("  Loading", metric_name, "...")
  
  if (!file.exists(file_path)) {
    cat(" FILE NOT FOUND\n")
    return(NULL)
  }
  
  # Get file size
  file_size_mb <- file.size(file_path) / 1024^2
  cat(" (", round(file_size_mb, 1), "MB)...")
  
  if (test_mode) {
    # In test mode, load only the subset we need
    cat(" TEST MODE - loading subset...")
    
    # Use fread for faster loading
    dist_data <- fread(file_path, nrows = length(sample_ids) + 1, showProgress = FALSE)
    
    # First column is SEQN
    all_seqns <- as.character(dist_data[[1]])
    
    # Find indices of our samples
    keep_idx <- which(all_seqns %in% sample_ids)
    
    if (length(keep_idx) == 0) {
      cat(" No matching samples\n")
      return(NULL)
    }
    
    # Subset the matrix
    dist_mat <- as.matrix(dist_data[keep_idx, (keep_idx + 1), with = FALSE])
    rownames(dist_mat) <- all_seqns[keep_idx]
    colnames(dist_mat) <- all_seqns[keep_idx]
    
  } else {
    # Full mode - load entire matrix
    cat(" FULL MODE - loading complete matrix...")
    dist_data <- fread(file_path, showProgress = TRUE)
    
    all_seqns <- as.character(dist_data[[1]])
    keep_idx <- which(all_seqns %in% sample_ids)
    
    dist_mat <- as.matrix(dist_data[keep_idx, (keep_idx + 1), with = FALSE])
    rownames(dist_mat) <- all_seqns[keep_idx]
    colnames(dist_mat) <- all_seqns[keep_idx]
  }
  
  # Convert to dist object
  dist_obj <- as.dist(dist_mat)
  
  cat(" Done (", length(keep_idx), "samples)\n")
  return(dist_obj)
}

# Load all three beta diversity metrics
beta_metrics <- list()

sample_ids_char <- demo_data$SEQN

bray_file <- file.path(beta_data_path, "dada2rsv-beta-braycurtis.txt")
beta_metrics$braycurtis <- load_beta_distance(bray_file, "Bray-Curtis", sample_ids_char, TEST_MODE)

unwunifrac_file <- file.path(beta_data_path, "dada2rsv-beta-unwunifrac.txt")
beta_metrics$unwunifrac <- load_beta_distance(unwunifrac_file, "Unweighted UniFrac", sample_ids_char, TEST_MODE)

wunifrac_file <- file.path(beta_data_path, "dada2rsv-beta-wunifrac.txt")
beta_metrics$wunifrac <- load_beta_distance(wunifrac_file, "Weighted UniFrac", sample_ids_char, TEST_MODE)

# Remove NULL entries
beta_metrics <- beta_metrics[!sapply(beta_metrics, is.null)]

cat("\nBeta diversity metrics loaded:", length(beta_metrics), "\n")
for (name in names(beta_metrics)) {
  cat("  ", name, ":", attr(beta_metrics[[name]], "Size"), "samples\n")
}

###############################################################################
## 4. PRINCIPAL COORDINATE ANALYSIS (PCoA)                                   ##
###############################################################################

cat("\nPerforming Principal Coordinate Analysis...\n")

# Function to perform PCoA
perform_pcoa_analysis <- function(dist_obj, metric_name) {
  cat("  Computing PCoA for", metric_name, "...\n")
  
  # Perform PCoA using ape::pcoa
  pcoa_result <- pcoa(dist_obj)
  
  # Extract coordinates (first 3 axes)
  pcoa_coords <- pcoa_result$vectors[, 1:3]
  colnames(pcoa_coords) <- c("PC1", "PC2", "PC3")
  
  # Calculate variance explained
  eigenvalues <- pcoa_result$values$Eigenvalues
  var_explained <- eigenvalues / sum(eigenvalues) * 100
  
  # Create dataframe with SEQN
  pcoa_df <- data.frame(
    SEQN = rownames(pcoa_coords),
    PC1 = pcoa_coords[, 1],
    PC2 = pcoa_coords[, 2],
    PC3 = pcoa_coords[, 3],
    stringsAsFactors = FALSE
  )
  
  cat("    Variance explained - PC1:", round(var_explained[1], 2), "%,",
      "PC2:", round(var_explained[2], 2), "%\n")
  
  return(list(
    coords = pcoa_df,
    var_explained = var_explained,
    eigenvalues = eigenvalues
  ))
}

# Perform PCoA for each metric
pcoa_results <- list()

for (metric_name in names(beta_metrics)) {
  pcoa_results[[metric_name]] <- perform_pcoa_analysis(beta_metrics[[metric_name]], metric_name)
}

###############################################################################
## 5. MERGE PCoA WITH CATEGORICAL DATA                                       ##
###############################################################################

cat("\nMerging PCoA coordinates with categorical data...\n")

for (metric_name in names(pcoa_results)) {
  pcoa_results[[metric_name]]$data <- pcoa_results[[metric_name]]$coords %>%
    left_join(demo_data, by = "SEQN") %>%
    filter(!is.na(Age_group))
  
  cat("  ", metric_name, ": Merged", nrow(pcoa_results[[metric_name]]$data), "samples\n")
}

###############################################################################
## 6. CREATE PCoA PLOTS WITH ELLIPSES AND CENTROIDS                          ##
###############################################################################

cat("\nCreating PCoA plots with 95% confidence ellipses...\n")

# Function to create PCoA plot with ellipses and centroids
create_pcoa_plot <- function(pcoa_data, var_explained, metric_name, color_var, 
                            var_label, colors) {
  
  # Calculate centroids for each group
  centroids <- pcoa_data %>%
    group_by(!!sym(color_var)) %>%
    summarise(
      PC1_centroid = mean(PC1, na.rm = TRUE),
      PC2_centroid = mean(PC2, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # Create the plot
  p <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = !!sym(color_var), fill = !!sym(color_var))) +
    # Add points
    geom_point(alpha = 0.6, size = 1) +
    
    # Add 95% confidence ellipses (dashed)
    stat_ellipse(geom = "polygon", alpha = 0.1, linetype = "dashed", size = 0.5, level = 0.95) +
    
    # Add centroids as larger points
    geom_point(data = centroids, 
               aes(x = PC1_centroid, y = PC2_centroid, color = !!sym(color_var)),
               size = 3, shape = 18, stroke = 1.5, show.legend = FALSE) +
    
    # Color scales
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    
    # Labels
    labs(
      title = paste(metric_name, "PCoA -", var_label),
      x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
      color = var_label,
      fill = var_label
    ) +
    
    # Theme
    egg::theme_article(base_size = 5) +
    theme(
      plot.title = element_text(size = 5, face = "bold"),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 4),
      legend.title = element_text(size = 5, face = "bold"),
      legend.text = element_text(size = 4),
      legend.position = "right"
    )
  
  return(p)
}

# Create plots for each metric and categorical variable
plot_counter <- 0

for (metric_name in names(pcoa_results)) {
  pcoa_data <- pcoa_results[[metric_name]]$data
  var_explained <- pcoa_results[[metric_name]]$var_explained
  
  metric_label <- gsub("_", " ", metric_name)
  metric_label <- paste0(toupper(substring(metric_label, 1, 1)), substring(metric_label, 2))
  
  # Plot by Age Group
  p_age <- create_pcoa_plot(pcoa_data, var_explained, metric_label, "Age_group", 
                           "Age Group", grafify_colors)
  
  filename <- file.path(plots_path, paste0("pcoa_", metric_name, "_by_age_group.pdf"))
  ggsave(filename, p_age, width = 4, height = 3, dpi = 300)
  cat("  ✓ Saved:", basename(filename), "\n")
  plot_counter <- plot_counter + 1
  
  # Plot by Gender
  p_gender <- create_pcoa_plot(pcoa_data, var_explained, metric_label, "Gender", 
                               "Gender", grafify_colors[c(2, 1)])
  
  filename <- file.path(plots_path, paste0("pcoa_", metric_name, "_by_gender.pdf"))
  ggsave(filename, p_gender, width = 3, height = 3, dpi = 300)
  cat("  ✓ Saved:", basename(filename), "\n")
  plot_counter <- plot_counter + 1
}

###############################################################################
## 7. CREATE SCREE PLOTS                                                     ##
###############################################################################

cat("\nCreating scree plots...\n")

for (metric_name in names(pcoa_results)) {
  var_explained <- pcoa_results[[metric_name]]$var_explained[1:10]  # First 10 axes
  
  scree_data <- data.frame(
    Axis = paste0("PC", 1:length(var_explained)),
    Variance = var_explained
  )
  
  scree_data$Axis <- factor(scree_data$Axis, levels = scree_data$Axis)
  
  metric_label <- gsub("_", " ", metric_name)
  metric_label <- paste0(toupper(substring(metric_label, 1, 1)), substring(metric_label, 2))
  
  p_scree <- ggplot(scree_data, aes(x = Axis, y = Variance)) +
    geom_bar(stat = "identity", fill = "#56B4E9", alpha = 0.7) +
    geom_line(aes(group = 1), color = "#E69F00", size = 0.8) +
    geom_point(color = "#E69F00", size = 2) +
    labs(
      title = paste(metric_label, "- Variance Explained"),
      x = "Principal Coordinate",
      y = "Variance Explained (%)"
    ) +
    egg::theme_article(base_size = 5) +
    theme(
      plot.title = element_text(size = 5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 4)
    )
  
  filename <- file.path(plots_path, paste0("scree_", metric_name, ".pdf"))
  ggsave(filename, p_scree, width = 3, height = 2.5, dpi = 300)
  cat("  ✓ Saved:", basename(filename), "\n")
}

###############################################################################
## 8. PERMANOVA TESTS (ADONIS2)                                              ##
###############################################################################

cat("\nPerforming PERMANOVA tests...\n")

permanova_results <- data.frame()

for (metric_name in names(beta_metrics)) {
  dist_obj <- beta_metrics[[metric_name]]
  pcoa_data <- pcoa_results[[metric_name]]$data
  
  # Ensure sample order matches
  common_samples <- intersect(labels(dist_obj), pcoa_data$SEQN)
  
  if (length(common_samples) < 10) {
    cat("  ", metric_name, ": Too few samples for PERMANOVA\n")
    next
  }
  
  # Subset distance matrix to common samples
  dist_mat <- as.matrix(dist_obj)
  dist_subset <- as.dist(dist_mat[common_samples, common_samples])
  
  # Subset metadata
  meta_subset <- pcoa_data %>%
    filter(SEQN %in% common_samples) %>%
    arrange(match(SEQN, common_samples))
  
  # PERMANOVA for Age Group
  cat("  ", metric_name, "- Age Group...")
  perm_age <- adonis2(dist_subset ~ Age_group, data = meta_subset, permutations = 999)
  
  permanova_results <- rbind(permanova_results, data.frame(
    Metric = metric_name,
    Variable = "Age_group",
    R2 = perm_age$R2[1],
    F_statistic = perm_age$F[1],
    P_value = perm_age$`Pr(>F)`[1],
    DF = perm_age$Df[1],
    N_samples = nrow(meta_subset)
  ))
  cat(" R² =", round(perm_age$R2[1], 4), "P =", round(perm_age$`Pr(>F)`[1], 4), "\n")
  
  # PERMANOVA for Gender
  cat("  ", metric_name, "- Gender...")
  perm_gender <- adonis2(dist_subset ~ Gender, data = meta_subset, permutations = 999)
  
  permanova_results <- rbind(permanova_results, data.frame(
    Metric = metric_name,
    Variable = "Gender",
    R2 = perm_gender$R2[1],
    F_statistic = perm_gender$F[1],
    P_value = perm_gender$`Pr(>F)`[1],
    DF = perm_gender$Df[1],
    N_samples = nrow(meta_subset)
  ))
  cat(" R² =", round(perm_gender$R2[1], 4), "P =", round(perm_gender$`Pr(>F)`[1], 4), "\n")
}

# Save PERMANOVA results
write.csv(permanova_results, file.path(output_path, "beta_diversity_permanova_results.csv"), 
          row.names = FALSE)

cat("\n  ✓ PERMANOVA results saved\n")

###############################################################################
## 9. DISTANCE TO CENTROID ANALYSIS                                          ##
###############################################################################

cat("\nCalculating distances to group centroids...\n")

centroid_distances <- data.frame()

for (metric_name in names(beta_metrics)) {
  dist_mat <- as.matrix(beta_metrics[[metric_name]])
  pcoa_data <- pcoa_results[[metric_name]]$data
  
  # Calculate distance to centroid for each group
  for (group_var in c("Age_group", "Gender")) {
    
    for (group_level in unique(pcoa_data[[group_var]])) {
      # Get samples in this group
      group_samples <- pcoa_data %>%
        filter(!!sym(group_var) == group_level) %>%
        pull(SEQN)
      
      group_samples <- intersect(group_samples, rownames(dist_mat))
      
      if (length(group_samples) < 2) next
      
      # Calculate pairwise distances within group
      group_dist <- dist_mat[group_samples, group_samples]
      
      # Distance to centroid = mean distance to all other samples in group
      for (sample in group_samples) {
        dist_to_centroid <- mean(group_dist[sample, group_samples[group_samples != sample]])
        
        centroid_distances <- rbind(centroid_distances, data.frame(
          Metric = metric_name,
          Variable = group_var,
          Group = as.character(group_level),
          SEQN = sample,
          Distance_to_Centroid = dist_to_centroid
        ))
      }
    }
  }
}

# Save centroid distances
write.csv(centroid_distances, file.path(output_path, "beta_diversity_centroid_distances.csv"),
          row.names = FALSE)

cat("  ✓ Centroid distances calculated and saved\n")

###############################################################################
## 10. DISTANCE TO CENTROID PLOTS                                            ##
###############################################################################

cat("\nCreating distance-to-centroid plots...\n")

for (metric_name in names(beta_metrics)) {
  
  metric_label <- gsub("_", " ", metric_name)
  metric_label <- paste0(toupper(substring(metric_label, 1, 1)), substring(metric_label, 2))
  
  # Plot for Age Group
  plot_data_age <- centroid_distances %>%
    filter(Metric == metric_name, Variable == "Age_group")
  
  if (nrow(plot_data_age) > 0) {
    p_centroid_age <- ggplot(plot_data_age, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(alpha = 0.3, size = 0.5, width = 0.2) +
      scale_fill_manual(values = grafify_colors) +
      labs(
        title = paste(metric_label, "- Distance to Centroid by Age"),
        x = "Age Group",
        y = "Distance to Centroid"
      ) +
      egg::theme_article(base_size = 5) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
        plot.title = element_text(size = 5, face = "bold")
      )
    
    filename <- file.path(plots_path, paste0("centroid_", metric_name, "_by_age.pdf"))
    ggsave(filename, p_centroid_age, width = 4, height = 3, dpi = 300)
    cat("  ✓ Saved:", basename(filename), "\n")
  }
  
  # Plot for Gender
  plot_data_gender <- centroid_distances %>%
    filter(Metric == metric_name, Variable == "Gender")
  
  if (nrow(plot_data_gender) > 0) {
    p_centroid_gender <- ggplot(plot_data_gender, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(alpha = 0.4, size = 0.8, width = 0.1) +
      scale_fill_manual(values = grafify_colors[c(2, 1)]) +
      labs(
        title = paste(metric_label, "- Distance to Centroid by Gender"),
        x = "Gender",
        y = "Distance to Centroid"
      ) +
      egg::theme_article(base_size = 5) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 5, face = "bold")
      )
    
    filename <- file.path(plots_path, paste0("centroid_", metric_name, "_by_gender.pdf"))
    ggsave(filename, p_centroid_gender, width = 3, height = 3, dpi = 300)
    cat("  ✓ Saved:", basename(filename), "\n")
  }
}

###############################################################################
## 11. COMBINED COMPARISON PLOTS                                             ##
###############################################################################

cat("\nCreating combined comparison plots...\n")

# Combine all PCoA plots for Age Group
if (length(pcoa_results) >= 3) {
  pcoa_plots_age <- list()
  
  for (metric_name in names(pcoa_results)) {
    pcoa_data <- pcoa_results[[metric_name]]$data
    var_explained <- pcoa_results[[metric_name]]$var_explained
    metric_label <- gsub("_", " ", metric_name)
    metric_label <- paste0(toupper(substring(metric_label, 1, 1)), substring(metric_label, 2))
    
    pcoa_plots_age[[metric_name]] <- create_pcoa_plot(
      pcoa_data, var_explained, metric_label, "Age_group", "Age Group", grafify_colors
    ) + theme(legend.position = "none")
  }
  
  # Combine plots
  combined_pcoa <- plot_grid(plotlist = pcoa_plots_age, ncol = 3, labels = c("A", "B", "C"))
  
  filename <- file.path(plots_path, "pcoa_all_metrics_by_age_combined.pdf")
  ggsave(filename, combined_pcoa, width = 9, height = 3, dpi = 300)
  cat("  ✓ Saved: pcoa_all_metrics_by_age_combined.pdf\n")
}

###############################################################################
## 12. FINAL SUMMARY                                                         ##
###############################################################################

cat("\n=== BETA DIVERSITY ANALYSIS SUMMARY ===\n")
cat("Mode:", ifelse(TEST_MODE, "TEST", "PRODUCTION"), "\n")
cat("Samples analyzed:", nrow(demo_data), "\n")
cat("Beta diversity metrics:", length(beta_metrics), "\n")
cat("PCoA plots created:", plot_counter, "\n")
cat("PERMANOVA tests:", nrow(permanova_results), "\n")

cat("\nPERMANOVA Results:\n")
print(permanova_results)

cat("\nSignificant PERMANOVA results (P < 0.05):\n")
sig_perm <- permanova_results %>% filter(P_value < 0.05)
if (nrow(sig_perm) > 0) {
  print(sig_perm)
} else {
  cat("  No significant results\n")
}

cat("\nOutputs saved to:", output_path, "\n")
cat("Plots saved to:", plots_path, "\n")
cat("========================================\n")

# Save workspace for potential further analysis
save.image(file.path(output_path, "beta_diversity_workspace.RData"))
cat("\n✓ Workspace saved to: beta_diversity_workspace.RData\n")

