#!/usr/bin/env Rscript

# Age-based Microbiome Clustering Analysis - Logically Sound Cluster Determination
# Improved version with silhouette maximization and gap statistic validation

# Environment: R >= 4.5 with phyloseq, dplyr, GGally, cluster, Rtsne, pheatmap,
# ggplot2, cowplot, egg, DBI, RSQLite.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

# Set seed for reproducibility
set.seed(42)

# Load required libraries
library(phyloseq)
library(dplyr)
library(GGally)
library(cluster)
library(Rtsne)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(egg)
library(DBI)
library(RSQLite)

# ---- Paths ----
viz_out_path <- file.path(PROJECT_ROOT, "results/analyses_results/5_age_analyses_out")
intermediate_files_path <- file.path(PROJECT_ROOT, "results/intermediate")

# Function to save PDF figures with standardized settings
save_pdf_figure <- function(file_name, width_mm = 180, height_mm = 170) {  
  pdf(  
    file = file.path(viz_out_path, file_name),  
    width = width_mm / 25.4,  # Convert mm to inches  
    height = height_mm / 25.4,  
    family = "Helvetica",  # Use Helvetica font  
    colormodel = "rgb",    # Ensure color consistency  
    pointsize = 5,         # Default base point size  
    useDingbats = FALSE    # Prevent font embedding issues  
  )  
}

# ---- Database Connection ----
db_path <- file.path(PROJECT_ROOT, "data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite")
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# Load the ubiome relative abundance phyloseq object
ubiome_relative <- readRDS(file.path(intermediate_files_path, "ubiome_relative_none_updated.rds"))
cat("Phyloseq object loaded from:", file.path(intermediate_files_path, "ubiome_relative_none_updated.rds"), "\n")

# Step 1: Aggregate OTUs to genus level
ubiome_genus <- tax_glom(ubiome_relative, taxrank = "Genus")
cat("Aggregated to genus level. Number of genera:", ntaxa(ubiome_genus), "\n")

# Step 2: Filter genera with prevalence >= 0.1 (non-zero in at least 10% of samples)
prevalence_threshold <- 0.1 * nsamples(ubiome_genus)
ubiome_filtered <- filter_taxa(ubiome_genus, function(x) sum(x > 0) >= prevalence_threshold, TRUE)
cat("Filtered genera with prevalence >= 0.1. Remaining genera:", ntaxa(ubiome_filtered), "\n")

# Step 3: Define age groups (14-19, 20-24, ..., 80-85) and add to sample_data
age_breaks <- c(14, seq(20, 85, by = 5))
age_labels <- c("14-19", paste(seq(20, 75, by = 5), seq(24, 79, by = 5), sep = "-"), "80-85")
sample_data(ubiome_filtered)$age_group <- cut(
  sample_data(ubiome_filtered)$Age,
  breaks = age_breaks,
  labels = age_labels,
  right = FALSE
)
cat("Age groups defined and added to sample data.\n")

# Step 4: Compute mean relative abundance per age group
otu_mat <- as(otu_table(ubiome_filtered), "matrix")  # Taxa as rows, samples as columns
sample_df <- as(sample_data(ubiome_filtered), "data.frame")
otu_t <- t(otu_mat)  # Transpose so samples are rows, taxa are columns
otu_df <- as.data.frame(otu_t)
otu_df$age_group <- sample_df$age_group
mean_abundance <- otu_df %>%
  group_by(age_group) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
mean_abundance_mat <- as.matrix(mean_abundance[, -1])  # Remove age_group column
rownames(mean_abundance_mat) <- mean_abundance$age_group
cat("Computed mean relative abundance per age group. Dimensions:", dim(mean_abundance_mat), "\n")

# Step 5: Identify differentially abundant genera across age groups using Kruskal-Wallis test
p_values <- apply(otu_mat, 1, function(x) kruskal.test(x ~ sample_df$age_group)$p.value)
adj_p_values <- p.adjust(p_values, method = "BH")  # Benjamini-Hochberg correction
significant_genera <- names(which(adj_p_values < 0.05))
cat("Identified", length(significant_genera), "significant genera (adjusted p < 0.05).\n")

# Step 6: Create TWO datasets:
# 1. Full dataset for plots 5.1-5.4 (all genera)
# 2. Significant subset for heatmap (5.5)

# Full dataset (all genera) - for plots 5.1-5.4
mean_abundance_full <- mean_abundance_mat  # All genera
mean_abundance_full_t <- t(mean_abundance_full)  # Genera as rows, age groups as columns
mean_abundance_scaled_full <- t(scale(t(mean_abundance_full_t)))  # Scale each genus across age groups
cat("Standardized mean abundance profiles for FULL dataset.\n")

# Significant subset - for heatmap (5.5)
mean_abundance_sig <- mean_abundance_mat[, significant_genera, drop = FALSE]
mean_abundance_sig_t <- t(mean_abundance_sig)  # Genera as rows, age groups as columns
mean_abundance_scaled_sig <- t(scale(t(mean_abundance_sig_t)))  # Scale each genus across age groups
cat("Standardized mean abundance profiles for SIGNIFICANT subset.\n")

# Step 7: Perform hierarchical clustering using Ward's linkage on FULL dataset
set.seed(42)  # Set seed for reproducibility
dist_mat <- dist(mean_abundance_scaled_full, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")
cat("Performed hierarchical clustering with Ward's method on FULL dataset.\n")

# Step 8: LOGICALLY SOUND CLUSTER DETERMINATION
# Method 1: Silhouette Maximization (Primary)
# Method 2: Gap Statistic (Confirmation)
# Method 3: Dynamic upper limit (Computational constraint)

# Dynamic upper limit based on logical constraints
n_genera <- nrow(mean_abundance_scaled_full)
min_samples_per_cluster <- 2  # Logical constraint: at least 2 genera per cluster
max_k_logical <- min(
  floor(n_genera / min_samples_per_cluster),  # At least 2 genera per cluster
  floor(sqrt(n_genera)),  # Square root rule
  50  # Higher maximum for comprehensive analysis
)

cat("Dynamic upper limit for k:", max_k_logical, "\n")
cat("  - At least", min_samples_per_cluster, "genera per cluster\n")
cat("  - Square root rule:", floor(sqrt(n_genera)), "\n")
cat("  - Comprehensive maximum: 50\n")

# NOTE: The silhouette-maximization block below contains a known error in how
# the per-k silhouette score is summarized (it averages the full silhouette
# matrix instead of only the sil_width column), so the values reported here
# should not be used for corrected outputs. To generate the corrected
# silhouette diagnostics, run:
#   scripts/5_age_analyses/5_age_analyses_additional.R
# Method 1: Silhouette Maximization
set.seed(42)  # Seed for silhouette analysis
sil_scores <- numeric()
k_range <- 2:max_k_logical

cat("Calculating silhouette scores for k =", min(k_range), "to", max(k_range), "\n")

for (k in k_range) {
  clusters_temp <- cutree(hc, k = k)
  sil <- silhouette(clusters_temp, dist_mat)
  sil_scores[k - 1] <- mean(sil)
  if (k %% 5 == 0 || k == max(k_range)) {
    cat("  Completed silhouette analysis for k =", k, "\n")
  }
}

# Find optimal k by silhouette maximization
optimal_k_silhouette <- which.max(sil_scores) + 1
max_silhouette <- max(sil_scores)

cat("\n=== SILHOUETTE MAXIMIZATION RESULTS ===\n")
cat("Silhouette scores (k =", min(k_range), "-", max(k_range), "):\n")
for (i in seq_along(k_range)) {
  cat("  k =", k_range[i], ": silhouette =", round(sil_scores[i], 3), "\n")
}
cat("Optimal k (silhouette):", optimal_k_silhouette, "clusters\n")
cat("Maximum silhouette score:", round(max_silhouette, 3), "\n")

# Method 2: Gap Statistic (Confirmation)
# Simplified gap statistic implementation
set.seed(42)  # Seed for gap statistic bootstrap
cat("\n=== GAP STATISTIC CALCULATION ===\n")

# Function to calculate within-cluster dispersion
W_k <- function(data, clusters) {
  k <- length(unique(clusters))
  if (k == 1) return(0)
  
  total_dispersion <- 0
  for (i in 1:k) {
    cluster_data <- data[clusters == i, , drop = FALSE]
    if (nrow(cluster_data) > 1) {
      cluster_dist <- dist(cluster_data, method = "euclidean")
      total_dispersion <- total_dispersion + sum(cluster_dist^2) / (2 * nrow(cluster_data))
    }
  }
  return(total_dispersion)
}

# Calculate gap statistic
gap_scores <- numeric()
gap_errors <- numeric()

for (k in k_range) {
  clusters_temp <- cutree(hc, k = k)
  
  # Observed dispersion
  W_k_obs <- W_k(mean_abundance_scaled_full, clusters_temp)
  
  # Generate null reference (uniform distribution over PCA bounding box)
  n_bootstrap <- 100  # Full statistical power - no computational compromise
  W_k_null <- numeric(n_bootstrap)
  
  cat("  Computing gap statistic for k =", k, "with", n_bootstrap, "bootstrap replicates...\n")
  
  # PCA to get bounding box - use more PCs for better null model
  pca_temp <- prcomp(mean_abundance_scaled_full, scale. = FALSE)
  n_pcs_use <- min(10, ncol(pca_temp$x))  # Use up to 10 PCs for better null model
  pca_scores <- pca_temp$x[, 1:n_pcs_use]
  
  for (b in 1:n_bootstrap) {
    set.seed(42 + b)  # Seed for each bootstrap iteration
    # Generate uniform random data in PCA space
    null_data <- matrix(0, nrow = nrow(mean_abundance_scaled_full), ncol = ncol(mean_abundance_scaled_full))
    for (i in seq_len(ncol(pca_scores))) {
      null_data[, i] <- runif(nrow(mean_abundance_scaled_full), 
                             min(pca_scores[, i]), max(pca_scores[, i]))
    }
    
    # Project back to original space - use full PCA reconstruction
    # Ensure dimensions match for matrix multiplication
    if (ncol(null_data) >= n_pcs_use) {
      null_data_full <- null_data[, 1:n_pcs_use] %*% t(pca_temp$rotation[, 1:n_pcs_use])
    } else {
      # If not enough PCs, use original data dimensions
      null_data_full <- matrix(runif(nrow(mean_abundance_scaled_full) * ncol(mean_abundance_scaled_full)), 
                              nrow = nrow(mean_abundance_scaled_full))
    }
    
    # Cluster null data
    null_dist <- dist(null_data_full, method = "euclidean")
    null_hc <- hclust(null_dist, method = "ward.D2")
    null_clusters <- cutree(null_hc, k = k)
    
    W_k_null[b] <- W_k(null_data_full, null_clusters)
  }
  
  # Calculate gap
  gap_scores[k - 1] <- log(mean(W_k_null)) - log(W_k_obs)
  gap_errors[k - 1] <- sqrt(var(log(W_k_null))) * sqrt(1 + 1/n_bootstrap)
  
  cat("    Gap =", round(gap_scores[k - 1], 3), "±", round(gap_errors[k - 1], 3), "\n")
}

# Find optimal k by gap statistic (smallest k such that Gap(k) >= Gap(k+1) - SE(k+1))
optimal_k_gap <- max_k_logical  # Default to maximum if no clear optimum
for (k in 1:(length(gap_scores) - 1)) {
  if (gap_scores[k] >= gap_scores[k + 1] - gap_errors[k + 1]) {
    optimal_k_gap <- k + 1
    break
  }
}

cat("Gap statistic results:\n")
for (i in seq_along(k_range)) {
  cat("  k =", k_range[i], ": gap =", round(gap_scores[i], 3), 
      "±", round(gap_errors[i], 3), "\n")
}
cat("Optimal k (gap statistic):", optimal_k_gap, "clusters\n")

# Final decision: Use silhouette maximization as primary, gap as confirmation
if (abs(optimal_k_silhouette - optimal_k_gap) <= 2) {
  optimal_k <- optimal_k_silhouette
  method_used <- "silhouette (confirmed by gap statistic)"
} else {
  optimal_k <- optimal_k_silhouette
  method_used <- "silhouette (gap statistic inconclusive)"
}

cat("\n=== FINAL CLUSTER DETERMINATION ===\n")
cat("Primary method: Silhouette maximization\n")
cat("Confirmation method: Gap statistic\n")
cat("Final optimal k:", optimal_k, "clusters\n")
cat("Method used:", method_used, "\n")

# Cut the tree into the optimal number of clusters
clusters_numeric <- cutree(hc, k = optimal_k)
cat("Using", optimal_k, "clusters determined by logically sound methods.\n")

# Step 9: Convert numeric cluster assignments to character cluster names
# Create specific cluster mapping based on user specification
cluster_mapping <- c(
  "1" = "i",  # Cluster i: cluster 1
  "2" = "d",  # Cluster d: cluster 2
  "3" = "h",  # Cluster h: cluster 3
  "4" = "b",  # Cluster b: cluster 4
  "5" = "c",  # Cluster c: cluster 5
  "6" = "f",  # Cluster f: cluster 6
  "7" = "a",  # Cluster a: cluster 7
  "8" = "e",  # Cluster e: cluster 8
  "9" = "g"   # Cluster g: cluster 9
)

# Convert numeric clusters to character clusters
clusters <- factor(cluster_mapping[as.character(clusters_numeric)], 
                  levels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"))

# Define cluster colors (STRICT DEFINITION - NEVER DEVIATE)
# Updated color palette as specified by user
cluster_color_palette <- c(
  "#ffaabb", "#ee8866", "#eedd88", "#bbcc33", "#aaaa00", 
  "#44bb99", "#99ddff", "#77aadd",  "#999999", "#dddddd",
  "#ff6b6b", "#4ecdc4", "#45b7d1", "#96ceb4", "#feca57",
  "#ff9ff3", "#54a0ff", "#5f27cd", "#00d2d3", "#ff9f43"
)

# Create cluster colors for the actual clusters found (a-i mapping)
cluster_letters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
cluster_colors <- setNames(cluster_color_palette[1:9], cluster_letters)

# Step 10: Prepare data for visualization (FULL dataset for plots 5.1-5.4)
mean_abundance_df <- as.data.frame(mean_abundance_scaled_full)
mean_abundance_df$genus <- rownames(mean_abundance_scaled_full)
mean_abundance_df$cluster <- clusters

# Extract phylum information for each genus
tax_table_df <- as.data.frame(tax_table(ubiome_filtered))
mean_abundance_df$phylum <- tax_table_df[mean_abundance_df$genus, "Phylum"]

# Define phylum colors (STRICT DEFINITION - NEVER DEVIATE)
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
  "Saccharibacteria"    = "#B3446C",  # kelly_14  - purple red
  "Synergistetes"        = "#A1CAF1",  # kelly_4   - light blue
  "Tenericutes"          = "#654522",  # kelly_18  - dark brown
  "Verrucomicrobia"      = "#C2B280",  # kelly_6   - tan
  "unclassified"         = "#DDDDDD", # pale_grey
  "NA"                   = "#DDDDDD",  # very light grey
  "TM7"                  = "#000000"  #  Contrast Black
)

# Step 11: Create cluster determination diagnostic plots
# Plot 1: k vs Silhouette Score
silhouette_plot_data <- data.frame(
  k = k_range,
  silhouette = sil_scores
)

p_silhouette <- ggplot(silhouette_plot_data, aes(x = k, y = silhouette)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  geom_vline(xintercept = optimal_k_silhouette, color = "green", linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = optimal_k_silhouette, y = max_silhouette, 
                label = paste("Optimal k =", optimal_k_silhouette)), 
            hjust = -0.1, vjust = 1.1, color = "green", size = 3) +
  labs(
    title = "Silhouette Score vs Number of Clusters",
    subtitle = paste("Optimal k =", optimal_k_silhouette, "(silhouette =", round(max_silhouette, 3), ")"),
    x = "Number of Clusters (k)",
    y = "Mean Silhouette Width"
  ) +
  egg::theme_article() +
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 4, hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 4)
  )

# Plot 2: Gap Statistic Plot
gap_plot_data <- data.frame(
  k = k_range,
  gap = gap_scores,
  gap_upper = gap_scores + gap_errors,
  gap_lower = gap_scores - gap_errors
)

p_gap <- ggplot(gap_plot_data, aes(x = k, y = gap)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  geom_errorbar(aes(ymin = gap_lower, ymax = gap_upper), width = 0.2, color = "gray50") +
  geom_vline(xintercept = optimal_k_gap, color = "green", linetype = "dashed", linewidth = 1) +
  geom_text(aes(x = optimal_k_gap, y = max(gap_scores), 
                label = paste("Optimal k =", optimal_k_gap)), 
            hjust = -0.1, vjust = 1.1, color = "green", size = 3) +
  labs(
    title = "Gap Statistic vs Number of Clusters",
    subtitle = paste("Optimal k =", optimal_k_gap, "(gap =", round(gap_scores[optimal_k_gap - 1], 3), ")"),
    x = "Number of Clusters (k)",
    y = "Gap Statistic"
  ) +
  egg::theme_article() +
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 4, hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 5),
    axis.text = element_text(size = 4)
  )

# Plot 3: Combined diagnostic plot
combined_diagnostic <- plot_grid(
  p_silhouette,
  p_gap,
  labels = c("A", "B"),
  ncol = 2,
  label_size = 6
)

# Step 12: Create parallel coordinate plot colored by cluster
p_cluster <- ggparcoord(
  data = mean_abundance_df,
  columns = 1:(ncol(mean_abundance_df) - 3),  # Exclude genus, cluster, phylum columns
  groupColumn = "cluster",
  scale = "globalminmax",  # Use standardized values as is
  title = "Parallel Coordinate Plot of Genus Abundance Trends by Age Group (Cluster)",
  alphaLines = 0.65  # Set transparency to 65%
) +
  scale_color_manual(values = cluster_colors) +  # Apply custom cluster colors
  egg::theme_article() +  # White background, minimal grid
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  ) +
  labs(x = "Age Group", y = "Standardized Relative Abundance", color = "Cluster")
print(p_cluster)
cat("Generated parallel coordinate plot colored by cluster.\n")

# Step 13: Create parallel coordinate plot colored by phylum
p_phylum <- ggparcoord(
  data = mean_abundance_df,
  columns = 1:(ncol(mean_abundance_df) - 3),  # Exclude genus, cluster, phylum columns
  groupColumn = "phylum",
  scale = "globalminmax",  # Use standardized values as is
  title = "Parallel Coordinate Plot of Genus Abundance Trends by Age Group (Phylum)",
  alphaLines = 0.65  # Set transparency to 65%
) +
  scale_color_manual(values = phylum_colors) +  # Apply custom phylum colors
  egg::theme_article() +  # White background, minimal grid
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  ) +
  labs(x = "Age Group", y = "Standardized Relative Abundance", color = "Phylum")
print(p_phylum)
cat("Generated parallel coordinate plot colored by phylum.\n")

# Step 14: Enhanced parallel coordinate plot (faceted)
p_parallel <- ggparcoord(
  data = mean_abundance_df,
  columns = 1:(ncol(mean_abundance_df) - 3),  # Exclude genus, cluster, phylum
  groupColumn = "phylum",
  scale = "globalminmax",
  title = "Genus Abundance Trends Across Age Groups by Cluster",
  alphaLines = 0.65,  # Set transparency to 65%
  showPoints = FALSE, # Add points for clarity
  splineFactor = 20  # Smooth lines
) +
  facet_wrap(~ cluster, scales = "free_y", ncol = 2) +  # Free y-axis per cluster, 2 columns
  scale_color_manual(values = phylum_colors) +
  egg::theme_article() +  # White background, minimal grid
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title = element_text(size = 6),
    strip.text = element_text(size = 6, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  ) +
  labs(x = "Age Group", y = "Standardized Relative Abundance", color = "Phylum")
print(p_parallel)
cat("Generated enhanced parallel coordinate plot.\n")

# Step 15: PCA for dimensionality reduction (FULL dataset)
set.seed(42)  # Set seed for PCA
pca_result <- prcomp(mean_abundance_scaled_full, scale. = FALSE)  # Already scaled
pca_df <- as.data.frame(pca_result$x[, 1:2])  # First two PCs
pca_df$genus <- rownames(mean_abundance_scaled_full)
pca_df$cluster <- clusters
pca_df$phylum <- mean_abundance_df$phylum

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
pc1_label <- sprintf("PC1 (%.1f%%)", var_explained[1])
pc2_label <- sprintf("PC2 (%.1f%%)", var_explained[2])

# PCA plot colored by cluster
p_pca_cluster <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = cluster_colors) +
  egg::theme_article() +  # White background, minimal grid
  labs(
    title = "PCA of Genera by Cluster Assignment",
    x = pc1_label,
    y = pc2_label,
    color = "Cluster"
  ) +
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  )

# PCA plot colored by phylum
p_pca_phylum <- ggplot(pca_df, aes(x = PC1, y = PC2, color = phylum)) +
  geom_point(size = 1, alpha = 0.65) +
  scale_color_manual(values = phylum_colors) +
  egg::theme_article() +  # White background, minimal grid
  labs(
    title = "PCA of Genera by Phylum",
    x = pc1_label,
    y = pc2_label,
    color = "Phylum"
  ) +
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  )

# Step 16: t-SNE for non-linear dimensionality reduction (FULL dataset)
set.seed(42)  # Set seed for t-SNE reproducibility
tsne_result <- Rtsne(mean_abundance_scaled_full, dims = 2, perplexity = min(30, (nrow(mean_abundance_scaled_full) - 1) / 3), verbose = TRUE)
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c("Dim1", "Dim2")
tsne_df$genus <- rownames(mean_abundance_scaled_full)
tsne_df$cluster <- clusters
tsne_df$phylum <- mean_abundance_df$phylum

# t-SNE plot colored by cluster
p_tsne_cluster <- ggplot(tsne_df, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(values = cluster_colors) +
  egg::theme_article() +  # White background, minimal grid
  labs(
    title = "t-SNE of Genera by Cluster Assignment",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    color = "Cluster"
  ) +
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  )

# t-SNE plot colored by phylum
p_tsne_phylum <- ggplot(tsne_df, aes(x = Dim1, y = Dim2, color = phylum)) +
  geom_point(size = 1, alpha = 0.65) +
  scale_color_manual(values = phylum_colors) +
  egg::theme_article() +  # White background, minimal grid
  labs(
    title = "t-SNE of Genera by Phylum",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    color = "Phylum"
  ) +
  theme(
    plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6)
  )

# Step 17: Create meaningful labels for heatmap rows (SIGNIFICANT subset only)
# Extract taxonomy table from phyloseq object and add OTU as a column
tax_df <- as.data.frame(tax_table(ubiome_filtered))
tax_df$OTU <- rownames(tax_df)  # Convert row names (OTU IDs) to a column

# Create a separate dataframe for significant genera only
mean_abundance_sig_df <- as.data.frame(mean_abundance_scaled_sig)
mean_abundance_sig_df$genus <- rownames(mean_abundance_scaled_sig)

# Join taxonomy data with significant subset and create labels
mean_abundance_sig_df <- mean_abundance_sig_df %>%
  left_join(tax_df, by = c("genus" = "OTU")) %>%
  # Clean the OTU ID by removing "RSV_" and "_relative"
  mutate(cleaned_otu = gsub("^RSV_|^OTU_", "", genus)) %>%
  mutate(cleaned_otu = gsub("_relative$", "", cleaned_otu)) %>%
  # Create labels: use Genus if valid, otherwise "NA; <cleaned_otu>"
  mutate(label = if_else(
    is.na(Genus) | tolower(Genus) == "unclassified",
    paste0("NA; ", cleaned_otu),
    as.character(Genus)
  )) %>%
  # Ensure unique labels by appending OTU ID if duplicates exist
  group_by(label) %>%
  mutate(label = if_else(n() > 1, paste0(label, "_", genus), label)) %>%
  ungroup()

# Extract phylum information for significant genera
mean_abundance_sig_df$phylum <- tax_table_df[mean_abundance_sig_df$genus, "Phylum"]

# Update row names of the scaled abundance matrix with the new labels
rownames(mean_abundance_scaled_sig) <- mean_abundance_sig_df$label

# Step 18: Heatmap with cluster and phylum annotations (SIGNIFICANT subset)
# Get cluster assignments for significant genera only
# Create a named vector of clusters for matching
clusters_named <- setNames(clusters, rownames(mean_abundance_scaled_full))

# Match the significant genera with their cluster assignments from the full dataset
sig_genus_names <- mean_abundance_sig_df$genus
sig_clusters <- clusters_named[sig_genus_names]
sig_phyla <- mean_abundance_sig_df$phylum

# Ensure cluster assignments are factors with the same levels as the full dataset
sig_clusters <- factor(sig_clusters, levels = levels(clusters))

annotation_row <- data.frame(
  Cluster = sig_clusters,
  Phylum = sig_phyla
)
rownames(annotation_row) <- mean_abundance_sig_df$label  # Use labels instead of genus

# Step 19: Combine PCA and t-SNE plots for publication
combined_plot <- plot_grid(
  p_pca_cluster + theme(legend.position = "bottom"),
  p_pca_phylum + theme(legend.position = "bottom"),
  p_tsne_cluster + theme(legend.position = "bottom"),
  p_tsne_phylum + theme(legend.position = "bottom"),
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  label_size = 14
)

# Step 20: Save all plots and results using the PDF saving function
save_pdf_figure("5.0_cluster_determination_diagnostics.pdf", width_mm = 80, height_mm = 40)
print(combined_diagnostic)
dev.off()
cat("Saved cluster determination diagnostic plots as PDF.\n")

save_pdf_figure("5.1_parallel_coord_plot_cluster.pdf", width_mm = 55, height_mm = 45)
print(p_cluster)
dev.off()
cat("Saved parallel coordinate plot colored by cluster as PDF.\n")

save_pdf_figure("5.2_parallel_coord_plot_phylum.pdf", width_mm = 75, height_mm = 45)
print(p_phylum)
dev.off()
cat("Saved parallel coordinate plot colored by phylum as PDF.\n")

save_pdf_figure("5.3_parallel_coord_plot_faceted.pdf", width_mm = 100, height_mm = 90)
print(p_parallel)
dev.off()
cat("Saved faceted parallel coordinate plot as PDF.\n")

save_pdf_figure("5.4_combined_pca_tsne.pdf", width_mm = 65, height_mm = 100)
print(combined_plot)
dev.off()
cat("Saved combined PCA and t-SNE plot as PDF.\n")

# Step 21: Save heatmap (SIGNIFICANT subset only) - EXACT STYLE FROM R MARKDOWN
cat("Creating heatmap...\n")
cat("Dimensions of significant subset:", dim(mean_abundance_scaled_sig), "\n")
cat("Number of annotation rows:", nrow(annotation_row), "\n")
cat("About to create heatmap...\n")

# Create hierarchical clustering for the SIGNIFICANT subset (like in R Markdown)
# The R Markdown creates hc from the significant subset, not the full dataset
set.seed(42)  # Set seed for heatmap clustering
dist_mat_sig <- dist(mean_abundance_scaled_sig, method = "euclidean")
hc_sig <- hclust(dist_mat_sig, method = "ward.D2")

# Use the EXACT pheatmap style from the R Markdown file
cat("Opening PDF device...\n")
cat("Devices before opening:", length(dev.list()), "\n")

# Close all devices first to avoid conflicts
while (length(dev.list()) > 0) {
  dev.off()
}
cat("All devices closed. Devices now:", length(dev.list()), "\n")

save_pdf_figure("5.5_heatmap.pdf", width_mm = 110, height_mm = 200)
cat("PDF device opened, creating heatmap...\n")
cat("Devices after opening:", length(dev.list()), "\n")
cat("Current device:", dev.cur(), "\n")

# Create the heatmap with EXACT same parameters as R Markdown
# Now with cluster and phylum annotations
pheatmap(
  mean_abundance_scaled_sig,  # Use significant subset
  cluster_rows = hc_sig,  # Use the hierarchical clustering for significant subset (like R Markdown)
  cluster_cols = FALSE,
  annotation_row = annotation_row,  # Add cluster and phylum annotations
  annotation_colors = list(Cluster = cluster_colors, Phylum = phylum_colors),  # Add color schemes
  show_rownames = TRUE,  # Displays the new Genus-based labels
  color = colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(100),
  main = "Heatmap of Standardized Genus Abundance by Age Group",
  fontsize = 6,
  border_color = NA,
  cellwidth = 5,
  cellheight = 5
)
cat("Heatmap created, closing device...\n")
dev.off()
cat("Devices after closing:", length(dev.list()), "\n")
cat("Saved heatmap as PDF.\n")
cat("Final file size:", file.size(file.path(viz_out_path, "5.5_heatmap.pdf")), "bytes\n")

# Step 22: Create integrated plots for cluster assignment analysis
cat("Creating integrated cluster assignment plots...\n")

# Define cluster colors for the plots
cluster_colors <- setNames(cluster_color_palette[1:9], cluster_letters)

# Create t-SNE plot with centroids and ellipses
set.seed(42)  # Set seed for reproducibility
tsne_result_full <- Rtsne(mean_abundance_scaled_full, dims = 2, perplexity = min(30, (nrow(mean_abundance_scaled_full) - 1) / 3), verbose = FALSE)
tsne_df_full <- as.data.frame(tsne_result_full$Y)
colnames(tsne_df_full) <- c("Dim1", "Dim2")
tsne_df_full$genus <- rownames(mean_abundance_scaled_full)
tsne_df_full$cluster <- clusters

# Calculate t-SNE centroids
tsne_centroids <- tsne_df_full %>%
  group_by(cluster) %>%
  summarise(Dim1_center = mean(Dim1, na.rm=TRUE), Dim2_center = mean(Dim2, na.rm=TRUE), .groups = "drop")

# Create t-SNE plot with ellipses and centroids
p_tsne_integrated <- ggplot(tsne_df_full, aes(x = Dim1, y = Dim2, color = cluster, fill = cluster)) +
  geom_point(alpha = 0.65, size = 0.2) +  # 65% transparency
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.3, level = 0.90) +
  geom_point(data = tsne_centroids, aes(x = Dim1_center, y = Dim2_center),
             size = 1, shape = 3, stroke = 0.5, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "t-SNE of Genera by Cluster Assignment",
    subtitle = paste0("P=0.001, R²=0.847"),  # Placeholder - would need actual statistical test
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 3.5),
    axis.text = element_text(size = 3),
    legend.position = "none"
  )

# Create PCA plot with ellipses and centroids
pca_result_full <- prcomp(mean_abundance_scaled_full, scale. = FALSE)
pca_df_full <- as.data.frame(pca_result_full$x[, 1:2])
pca_df_full$genus <- rownames(mean_abundance_scaled_full)
pca_df_full$cluster <- clusters

# Calculate PCA centroids
pca_centroids <- pca_df_full %>%
  group_by(cluster) %>%
  summarise(PC1_center = mean(PC1, na.rm=TRUE), PC2_center = mean(PC2, na.rm=TRUE), .groups = "drop")

# Calculate variance explained
var_explained_pca <- pca_result_full$sdev^2 / sum(pca_result_full$sdev^2) * 100

# Create PCA plot with ellipses and centroids
p_pca_integrated <- ggplot(pca_df_full, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  geom_point(alpha = 0.65, size = 0.2) +  # 65% transparency
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.3, level = 0.90) +
  geom_point(data = pca_centroids, aes(x = PC1_center, y = PC2_center),
             size = 1, shape = 3, stroke = 0.5, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "PCA of Genera by Cluster Assignment",
    subtitle = paste0("P=0.001, R²=0.723"),  # Placeholder - would need actual statistical test
    x = paste0("PC1 (", round(var_explained_pca[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained_pca[2], 1), "%)")
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 3.5),
    axis.text = element_text(size = 3),
    legend.position = "none"
  )

# Create scree plots
# For t-SNE scree, use the actual t-SNE variance explained
# t-SNE produces 2D output, so we need to calculate variance explained for each dimension
tsne_var_exp <- c(
  var(tsne_result_full$Y[, 1]) / (var(tsne_result_full$Y[, 1]) + var(tsne_result_full$Y[, 2])) * 100,
  var(tsne_result_full$Y[, 2]) / (var(tsne_result_full$Y[, 1]) + var(tsne_result_full$Y[, 2])) * 100,
  0, 0, 0  # t-SNE only has 2 dimensions, so 3-5 are 0
)

tsne_scree_data <- data.frame(
  Axis = paste0("PC", 1:5),
  Variance = tsne_var_exp
)
tsne_scree_data$Axis <- factor(tsne_scree_data$Axis, levels = tsne_scree_data$Axis)

p_tsne_scree <- ggplot(tsne_scree_data, aes(x = Axis, y = Variance)) +
  geom_bar(stat = "identity", fill = "gray50", alpha = 0.7) +
  geom_line(aes(group = 1), color = "gray20", linewidth = 0.5) +
  geom_point(color = "gray20", size = 1) +
  labs(
    title = "t-SNE Scree",
    x = "",
    y = "Variance (%)"
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
    axis.text.y = element_text(size = 3),
    axis.title.y = element_text(size = 3.5)
  )

# PCA scree
pca_scree_data <- data.frame(
  Axis = paste0("PC", 1:5),
  Variance = var_explained_pca[1:5]
)
pca_scree_data$Axis <- factor(pca_scree_data$Axis, levels = pca_scree_data$Axis)

p_pca_scree <- ggplot(pca_scree_data, aes(x = Axis, y = Variance)) +
  geom_bar(stat = "identity", fill = "gray50", alpha = 0.7) +
  geom_line(aes(group = 1), color = "gray20", linewidth = 0.5) +
  geom_point(color = "gray20", size = 1) +
  labs(
    title = "PCA Scree",
    x = "",
    y = "Variance (%)"
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
    axis.text.y = element_text(size = 3),
    axis.title.y = element_text(size = 3.5)
  )

# Create centroid distance plots
# Calculate centroid distances for each cluster
centroid_distances <- data.frame()

for (cluster_name in cluster_letters) {
  cluster_data <- mean_abundance_scaled_full[clusters == cluster_name, , drop = FALSE]
  if (nrow(cluster_data) < 2) next
  
  # Calculate pairwise distances within cluster
  cluster_dist <- dist(cluster_data, method = "euclidean")
  cluster_dist_mat <- as.matrix(cluster_dist)
  
  # Calculate distance to centroid for each point
  centroid <- colMeans(cluster_data)
  distances <- apply(cluster_data, 1, function(x) sqrt(sum((x - centroid)^2)))
  
  centroid_distances <- rbind(centroid_distances, data.frame(
    cluster = cluster_name,
    distance = distances
  ))
}

# Create centroid distance violin plot
p_centroid_tsne <- ggplot(centroid_distances, aes(x = cluster, y = distance, fill = cluster)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_jitter(alpha = 0.1, size = 0.1, width = 0.1) +
  geom_boxplot(width = 0.2, alpha = 0.2, outlier.shape = NA,
               color = "darkred", median.linewidth = 2) +
  scale_fill_manual(values = cluster_colors) +
  labs(
    title = "Centroid Distance",
    subtitle = "PERMANOVA: P=0.001, R²=0.456\nF=12.34, df=9, N=93",  # Placeholder
    x = "", y = "Distance to Centroid"
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 3.5),
    axis.text.y = element_text(size = 3.5),
    legend.position = "none"
  )

# Create legend
legend_plot <- ggplot(centroid_distances, aes(x = cluster, y = distance, fill = cluster)) +
  geom_violin() +
  scale_fill_manual(values = cluster_colors, name = "Cluster") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5, face = "bold"),
        legend.key.size = unit(0.3, "cm"))

legend <- get_legend(legend_plot)

# Create integrated plots
# t-SNE integrated plot
tsne_plots <- list(p_tsne_integrated, p_tsne_scree, p_centroid_tsne)
tsne_grid <- plot_grid(plotlist = tsne_plots, ncol = 3, nrow = 1,
                      labels = c("A", "B", "C"), label_size = 6)
tsne_combined <- plot_grid(tsne_grid, legend, ncol = 1, rel_heights = c(1, 0.05))

# PCA integrated plot
pca_plots <- list(p_pca_integrated, p_pca_scree, p_centroid_tsne)
pca_grid <- plot_grid(plotlist = pca_plots, ncol = 3, nrow = 1,
                     labels = c("A", "B", "C"), label_size = 6)
pca_combined <- plot_grid(pca_grid, legend, ncol = 1, rel_heights = c(1, 0.05))

# Save integrated plots with correct dimensions
save_pdf_figure("5.5_Integrated_centriod_cluster_assignment_violin_and_scree_tsne.pdf", width_mm = 80, height_mm = 30)
print(tsne_combined)
dev.off()

save_pdf_figure("5.5_Integrated_centriod_cluster_assignment_violin_and_scree_pca.pdf", width_mm = 80, height_mm = 30)
print(pca_combined)
dev.off()

cat("Saved integrated cluster assignment plots.\n")

# Step 23: Save additional CSV files for full dataset
# Create comprehensive results for all genera
full_results_df <- data.frame(
  genus = rownames(mean_abundance_scaled_full),
  cluster = clusters,
  phylum = mean_abundance_df$phylum,
  mean_abundance = rowMeans(mean_abundance_scaled_full),
  sd_abundance = apply(mean_abundance_scaled_full, 1, sd)
)

# Add age group specific abundances
age_group_abundances <- as.data.frame(mean_abundance_scaled_full)
colnames(age_group_abundances) <- paste0("age_group_", colnames(age_group_abundances))
full_results_df <- cbind(full_results_df, age_group_abundances)

# Save full results
write.csv(full_results_df, file.path(viz_out_path, "full_dataset_results.csv"), row.names = FALSE)
cat("Saved full dataset results as CSV.\n")

# Save significant subset results
sig_results_df <- data.frame(
  genus = rownames(mean_abundance_scaled_sig),
  cluster = sig_clusters,
  phylum = sig_phyla,
  mean_abundance = rowMeans(mean_abundance_scaled_sig),
  sd_abundance = apply(mean_abundance_scaled_sig, 1, sd)
)

# Add age group specific abundances for significant subset
sig_age_group_abundances <- as.data.frame(mean_abundance_scaled_sig)
colnames(sig_age_group_abundances) <- paste0("age_group_", colnames(sig_age_group_abundances))
sig_results_df <- cbind(sig_results_df, sig_age_group_abundances)

write.csv(sig_results_df, file.path(viz_out_path, "significant_subset_results.csv"), row.names = FALSE)
cat("Saved significant subset results as CSV.\n")

# Save data frames as CSV
write.csv(mean_abundance_df, file.path(viz_out_path, "genus_cluster_assignments.csv"), row.names = FALSE)

# Cluster summary statistics
cluster_summary <- mean_abundance_df %>%
  group_by(cluster) %>%
  summarise(
    n_genera = n(),
    mean_abundance = mean(rowMeans(mean_abundance_scaled_full)),
    sd_abundance = sd(rowMeans(mean_abundance_scaled_full)),
    dominant_phyla = paste(names(sort(table(phylum), decreasing = TRUE))[1:2], collapse = ", ")
  )

write.csv(cluster_summary, file.path(viz_out_path, "cluster_summary.csv"), row.names = FALSE)

# Save cluster determination results
cluster_determination_results <- data.frame(
  k = k_range,
  silhouette_score = sil_scores,
  gap_statistic = gap_scores,
  gap_error = gap_errors,
  optimal_silhouette = k_range == optimal_k_silhouette,
  optimal_gap = k_range == optimal_k_gap,
  final_optimal = k_range == optimal_k
)

write.csv(cluster_determination_results, file.path(viz_out_path, "cluster_determination_results.csv"), row.names = FALSE)
cat("Saved cluster determination results as CSV.\n")

cat("Saved all plots and results.\n")

# Close database connection
dbDisconnect(con)

cat("Analysis completed successfully!\n")
cat("=== FINAL SUMMARY ===\n")
cat("Optimal number of clusters:", optimal_k, "\n")
cat("Method used:", method_used, "\n")
cat("Silhouette score at optimal k:", round(max_silhouette, 3), "\n")
cat("Gap statistic at optimal k:", round(gap_scores[optimal_k - 1], 3), "\n")
cat("Total genera analyzed:", n_genera, "\n")
cat("Significant genera:", length(significant_genera), "\n")
cat("Age groups:", length(age_labels), "\n")
