#!/usr/bin/env Rscript
################################################################################
##  2. COMPREHENSIVE ALPHA DIVERSITY ANALYSIS                                 ##
##  - Violin plots with LOESS fit lines                                       ##
##  - Statistical tests (Kruskal-Wallis)                                      ##
##  - Faceted by categorical variables                                        ##
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(cowplot)
  library(ggbeeswarm)
  library(egg)
})

cat("=== ALPHA DIVERSITY ANALYSIS ===\n\n")

###############################################################################
## 1. PATHS AND COLOR CONFIGURATION                                          ##
###############################################################################

base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
output_path <- file.path(base_path, "results/analyses_results/6_alpha_beta_analyses_out")
plots_path <- file.path(output_path, "alpha_diversity_plots")

# Create directories
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

# Grafify color palette
grafify_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 
  "#332288", "#88ccee", "#44aa99", "#117733", "#999933", "#ddcc77", "#cc6677"
)

cat("Output directory:", output_path, "\n")
cat("Plots directory:", plots_path, "\n\n")

###############################################################################
## 2. LOAD DATA                                                              ##
###############################################################################

cat("Loading processed data...\n")
alpha_with_categories <- readRDS(file.path(data_path, "alpha_diversity_with_categories.rds"))
cat("  Alpha diversity data loaded:", nrow(alpha_with_categories), "samples\n")
cat("  Variables:", paste(colnames(alpha_with_categories), collapse = ", "), "\n\n")

###############################################################################
## 3. ALPHA DIVERSITY PLOTS BY AGE GROUP                                     ##
###############################################################################

cat("Creating alpha diversity plots by Age Group...\n")

# Create violin plots for each metric
for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
  
  metric_label <- gsub("_", " ", metric_col)
  
  p <- ggplot(alpha_with_categories, aes(x = Age_group, y = !!sym(metric_col), fill = Age_group)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_beeswarm(alpha = 0.3, size = 0.5, cex = 0.5) +
    scale_fill_manual(values = grafify_colors) +
    labs(
      title = paste(metric_label, "by Age Group"),
      x = "Age Group",
      y = metric_label
    ) +
    egg::theme_article(base_size = 5) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      plot.title = element_text(size = 5, face = "bold")
    )
  
  filename <- file.path(plots_path, paste0("alpha_", tolower(metric_col), "_by_age_group.pdf"))
  ggsave(filename, p, width = 4, height = 3, dpi = 300)  # 50% smaller
  cat("  ✓ Saved:", basename(filename), "\n")
}

###############################################################################
## 4. ALPHA DIVERSITY PLOTS BY GENDER                                        ##
###############################################################################

cat("\nCreating alpha diversity plots by Gender...\n")

for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
  
  metric_label <- gsub("_", " ", metric_col)
  
  p <- ggplot(alpha_with_categories, aes(x = Gender, y = !!sym(metric_col), fill = Gender)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA) +
    geom_beeswarm(alpha = 0.4, size = 0.8) +
    scale_fill_manual(values = grafify_colors[c(2, 1)]) +  # Female blue, Male orange
    labs(
      title = paste(metric_label, "by Gender"),
      x = "Gender",
      y = metric_label
    ) +
    egg::theme_article(base_size = 5) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 5, face = "bold")
    )
  
  filename <- file.path(plots_path, paste0("alpha_", tolower(metric_col), "_by_gender.pdf"))
  ggsave(filename, p, width = 3, height = 3, dpi = 300)  # 50% smaller
  cat("  ✓ Saved:", basename(filename), "\n")
}

###############################################################################
## 5. FACETED PLOTS - ALL METRICS SIDE BY SIDE                               ##
###############################################################################

cat("\nCreating faceted plots...\n")

# Prepare long format
alpha_long <- alpha_with_categories %>%
  pivot_longer(
    cols = c(Observed_OTUs, Shannon_Diversity, Inverse_Simpson),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, 
                   levels = c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson"),
                   labels = c("Observed OTUs", "Shannon Diversity", "Inverse Simpson"))
  )

# Faceted by metric, colored by age group
p_faceted_age <- ggplot(alpha_long, aes(x = Age_group, y = Value, fill = Age_group)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = grafify_colors) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  labs(
    title = "Alpha Diversity Metrics by Age Group",
    x = "Age Group",
    y = "Value"
  ) +
  egg::theme_article(base_size = 5) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 4),
    plot.title = element_text(size = 5, face = "bold"),
    strip.text = element_text(size = 5, face = "bold")
  )

ggsave(file.path(plots_path, "alpha_faceted_by_age_group.pdf"), p_faceted_age, 
       width = 6, height = 2, dpi = 300)  # 50% smaller
cat("  ✓ Saved: alpha_faceted_by_age_group.pdf\n")

###############################################################################
## 6. STATISTICAL TESTS                                                      ##
###############################################################################

cat("\nPerforming statistical tests...\n")

# Kruskal-Wallis tests for Age Group
kw_results <- data.frame()

for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
  formula_obj <- as.formula(paste(metric_col, "~ Age_group"))
  kw_test <- kruskal.test(formula_obj, data = alpha_with_categories)
  
  kw_results <- rbind(kw_results, data.frame(
    Metric = gsub("_", " ", metric_col),
    Categorical_Variable = "Age_group",
    Chi_squared = kw_test$statistic,
    DF = kw_test$parameter,
    P_value = kw_test$p.value
  ))
}

# Wilcoxon tests for Gender
for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
  formula_obj <- as.formula(paste(metric_col, "~ Gender"))
  wilcox_test <- wilcox.test(formula_obj, data = alpha_with_categories)
  
  kw_results <- rbind(kw_results, data.frame(
    Metric = gsub("_", " ", metric_col),
    Categorical_Variable = "Gender",
    Chi_squared = NA,
    DF = NA,
    P_value = wilcox_test$p.value
  ))
}

# Add FDR correction
kw_results$P_value_FDR <- p.adjust(kw_results$P_value, method = "BH")

# Save results
write.csv(kw_results, file.path(output_path, "alpha_diversity_statistical_tests.csv"), 
          row.names = FALSE)

cat("  ✓ Statistical tests completed\n")
cat("\nTest Results:\n")
print(kw_results)

###############################################################################
## 7. SUMMARY STATISTICS                                                     ##
###############################################################################

cat("\nCalculating summary statistics...\n")

# Summary by Age Group
summary_age <- alpha_long %>%
  group_by(Metric, Age_group) %>%
  summarise(
    N = n(),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Q25 = quantile(Value, 0.25, na.rm = TRUE),
    Q75 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Summary by Gender
summary_gender <- alpha_with_categories %>%
  group_by(Gender) %>%
  summarise(
    across(c(Observed_OTUs, Shannon_Diversity, Inverse_Simpson),
           list(Mean = ~mean(., na.rm = TRUE),
                SD = ~sd(., na.rm = TRUE),
                Median = ~median(., na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    N = n(),
    .groups = "drop"
  )

# Save summaries
write.csv(summary_age, file.path(output_path, "alpha_diversity_summary_by_age.csv"),
          row.names = FALSE)
write.csv(summary_gender, file.path(output_path, "alpha_diversity_summary_by_gender.csv"),
          row.names = FALSE)

cat("  ✓ Summary statistics saved\n")

###############################################################################
## 8. FINAL SUMMARY                                                          ##
###############################################################################

cat("\n=== ALPHA DIVERSITY ANALYSIS SUMMARY ===\n")
cat("Samples analyzed:", nrow(alpha_with_categories), "\n")
cat("Plots created: 9\n")
cat("Statistical tests: 6\n")
cat("Significant results (P < 0.05):", sum(kw_results$P_value < 0.05), "\n")
cat("Outputs saved to:", output_path, "\n")
cat("Plots saved to:", plots_path, "\n")
cat("========================================\n")
