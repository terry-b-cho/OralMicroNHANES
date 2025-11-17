#!/usr/bin/env Rscript
################################################################################
##  4. INTEGRATED ALPHA & BETA DIVERSITY ANALYSIS - HPC VERSION               ##
##  - 4×3 grid: Alpha + Beta centroid + PCoA + Scree                          ##
##  - Safe Grafify colors, dark red median, grayscale scree                   ##
##  - ALL 40+ categorical variables (when TEST_MODE = FALSE)                  ##
##  - HPC version without vegan dependency                                     ##
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ape)
  library(egg)
  library(cowplot)
  library(data.table)
  library(ggdist)
  library(ggsignif)
})

cat("=== INTEGRATED DIVERSITY ANALYSIS - SIMPLIFIED ===\n\n")

###############################################################################
## 1. CONFIGURATION                                                          ##
###############################################################################

TEST_MODE <- FALSE   # ← FULL DATASET (N=9,349) for all modes
TEST_SAMPLES <- 500  # Only used when TEST_MODE = TRUE

# Plot styling options
NO_TICK_AND_BRACKET <- TRUE  # Set to TRUE to remove x-tick labels and brackets, FALSE for normal plots
NO_TICK_AND_BRACKET_ASTERISK <- FALSE  # Set to TRUE to show asterisks (*, **, ***) instead of exact P-values
HALF_PLOTS <- TRUE  # Set to TRUE to use half violin plots with box+jitter on left side

base_path <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
beta_data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
output_path <- file.path(base_path, "results/analyses_results/6_alpha_beta_analyses_out")
plots_path <- file.path(output_path, "integrated_diversity_plots")

if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

# Safe Grafify palette
safe_colors <- c(
  "#88CCEE",  # safe_blue
  "#CC6677",  # safe_red
  "#DDCC77",  # safe_yellow
  "#117733",  # safe_green
  "#332288",  # safe_violet
  "#AA4499",  # safe_purple
  "#44AA99",  # safe_bluegreen
  "#999933",  # safe_bush
  "#882255",  # safe_reddish
  "#661100",  # safe_wine
  "#6699CC"   # safe_skyblue
)

cat("TEST_MODE:", TEST_MODE, "\n\n")

###############################################################################
## 2. LOAD DATA                                                              ##
###############################################################################

alpha_with_categories <- readRDS(file.path(data_path, "alpha_diversity_with_all_categories.rds"))

if (TEST_MODE) {
  alpha_with_categories <- alpha_with_categories[1:min(TEST_SAMPLES, nrow(alpha_with_categories)), ]
}

categorical_vars <- names(alpha_with_categories)[sapply(alpha_with_categories, function(x) {
  is.factor(x) && nlevels(x) >= 1 && nlevels(x) <= 10 && sum(!is.na(x)) >= 30
})]

categorical_vars <- setdiff(categorical_vars, c("sample", "cycle", "Release_Cycle"))

cat("Found", length(categorical_vars), "categorical variables for analysis\n")

###############################################################################
## 3. SIMPLIFIED BETA DIVERSITY ANALYSIS                                    ##
###############################################################################

cat("Loading beta diversity data...\n")

# Load beta diversity matrices
load_beta_subset <- function(file_path, sample_ids) {
  if (TEST_MODE) {
    dist_data <- fread(file_path, nrows = min(length(sample_ids) + 1, 1000))
  } else {
    dist_data <- fread(file_path)
  }
  all_seqns <- as.character(dist_data[[1]])
  keep_idx <- which(all_seqns %in% sample_ids)
  
  dist_mat <- as.matrix(dist_data[keep_idx, (keep_idx + 1), with = FALSE])
  rownames(dist_mat) <- all_seqns[keep_idx]
  colnames(dist_mat) <- all_seqns[keep_idx]
  
  return(list(dist = as.dist(dist_mat), matrix = dist_mat))
}

sample_ids <- alpha_with_categories$SEQN

beta_metrics <- list()
beta_metrics$Braycurtis <- load_beta_subset(
  file.path(beta_data_path, "dada2rsv-beta-braycurtis.txt"), sample_ids)
beta_metrics$Unwunifrac <- load_beta_subset(
  file.path(beta_data_path, "dada2rsv-beta-unwunifrac.txt"), sample_ids)
beta_metrics$Wunifrac <- load_beta_subset(
  file.path(beta_data_path, "dada2rsv-beta-wunifrac.txt"), sample_ids)

cat("  ✓ Loaded beta diversity matrices\n")

# Perform PCoA for all metrics using ape::pcoa
cat("Performing PCoA...\n")

pcoa_results <- list()

for (metric_name in names(beta_metrics)) {
  pcoa_res <- pcoa(beta_metrics[[metric_name]]$dist)
  
  pcoa_df <- data.frame(
    SEQN = rownames(pcoa_res$vectors),
    PC1 = pcoa_res$vectors[, 1],
    PC2 = pcoa_res$vectors[, 2]
  ) %>% 
    mutate(SEQN = as.character(SEQN)) %>%
    left_join(alpha_with_categories %>% mutate(SEQN = as.character(SEQN)), by = "SEQN")
  
  var_exp <- pcoa_res$values$Eigenvalues / sum(pcoa_res$values$Eigenvalues) * 100
  
  pcoa_results[[metric_name]] <- list(data = pcoa_df, var_explained = var_exp)
  cat("  ", metric_name, ": PC1=", round(var_exp[1], 1), "%, PC2=", round(var_exp[2], 1), "%\n", sep="")
}

cat("  ✓ PCoA complete\n")

# Calculate centroids (optimized for large datasets)
cat("Calculating centroids...\n")

centroid_data_all <- data.frame()

for (metric_name in names(beta_metrics)) {
  cat("  Processing", metric_name, "...\n")
  dist_mat <- beta_metrics[[metric_name]]$matrix
  
  for (cat_var in categorical_vars) {
    if (!cat_var %in% colnames(alpha_with_categories)) next
    
    cat("    Variable:", cat_var, "\n")
    unique_levels <- unique(alpha_with_categories[[cat_var]])
    unique_levels <- unique_levels[!is.na(unique_levels)]
    
    for (group_level in unique_levels) {
      group_samples <- alpha_with_categories %>%
        filter(!!sym(cat_var) == group_level) %>%
        pull(SEQN)
      
      group_samples <- intersect(group_samples, rownames(dist_mat))
      if (length(group_samples) < 2) next
      
      # Optimized centroid calculation
      group_dist <- dist_mat[group_samples, group_samples]
      
      # Calculate centroid distances more efficiently
      centroid_distances <- sapply(group_samples, function(sample) {
        other_samples <- group_samples[group_samples != sample]
        mean(group_dist[sample, other_samples])
      })
      
      # Add to dataframe in batch
      batch_data <- data.frame(
        Metric = metric_name,
        Variable = cat_var,
        Group = as.character(group_level),
        Distance_to_Centroid = centroid_distances,
        stringsAsFactors = FALSE
      )
      
      centroid_data_all <- rbind(centroid_data_all, batch_data)
    }
  }
}

cat("  ✓ Done\n")

###############################################################################
## 4. SIMPLIFIED STATISTICAL TESTS                                          ##
###############################################################################

cat("Running statistical tests...\n")

# Initialize data frames to store results
alpha_pairwise_results <- data.frame()
beta_pairwise_results <- data.frame()

# Simple statistical tests without vegan
for (cat_var in categorical_vars) {
  
  plot_data_alpha <- alpha_with_categories %>% filter(!is.na(!!sym(cat_var)))
  if (nrow(plot_data_alpha) < 30) next
  
  n_levels <- nlevels(plot_data_alpha[[cat_var]])
  plot_colors <- safe_colors[1:min(n_levels, length(safe_colors))]
  
  # ALPHA DIVERSITY TESTS
  for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
    
    # Check if we have enough groups for statistical test
    if (!is.factor(plot_data_alpha[[cat_var]])) {
      plot_data_alpha[[cat_var]] <- as.factor(plot_data_alpha[[cat_var]])
    }
    if (nlevels(droplevels(plot_data_alpha[[cat_var]])) < 2) {
      cat("  Skipping", metric_col, "for", cat_var, "- insufficient groups\n")
      next
    }
    
    # Perform Kruskal-Wallis test with error handling
    tryCatch({
      kw_test <- kruskal.test(as.formula(paste(metric_col, "~", cat_var)), data = plot_data_alpha)
    }, error = function(e) {
      cat("  Error in Kruskal-Wallis test for", metric_col, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(kw_test)) next
    
    # Calculate effect size (epsilon-squared for Kruskal-Wallis)
    n <- nrow(plot_data_alpha)
    k <- nlevels(plot_data_alpha[[cat_var]])
    epsilon_sq <- (kw_test$statistic - k + 1) / (n - k)
    
    # Store results
    if (kw_test$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
      # Perform pairwise Wilcoxon tests
      pwc <- pairwise.wilcox.test(plot_data_alpha[[metric_col]], 
                                 plot_data_alpha[[cat_var]], 
                                 p.adjust.method = "fdr", exact = FALSE)
      
      # Extract significant comparisons
      group_levels <- levels(plot_data_alpha[[cat_var]])
      
      for (i in 1:(length(group_levels)-1)) {
        for (j in (i+1):length(group_levels)) {
          p_val <- pwc$p.value[j-1, i]
          
          if (!is.na(p_val)) {
            alpha_pairwise_results <- rbind(alpha_pairwise_results, data.frame(
              Variable = cat_var,
              Metric = metric_col,
              Group1 = group_levels[i],
              Group2 = group_levels[j],
              P_value = p_val,
              P_value_FDR = p_val,
              Significant = p_val < 0.05,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # BETA DIVERSITY TESTS
  for (metric_name in names(beta_metrics)) {
    beta_plot_data <- centroid_data_all %>% filter(Metric == metric_name, Variable == cat_var)
    
    if (nrow(beta_plot_data) < 30) next
    
    # Check if we have enough groups for statistical test
    if (!is.factor(beta_plot_data$Group)) {
      beta_plot_data$Group <- as.factor(beta_plot_data$Group)
    }
    if (nlevels(droplevels(beta_plot_data$Group)) < 2) {
      next
    }
    
    # Perform Kruskal-Wallis test on centroid distances with error handling
    group_factor <- factor(beta_plot_data$Group, levels = levels(plot_data_alpha[[cat_var]]))
    tryCatch({
      kw_beta <- kruskal.test(Distance_to_Centroid ~ group_factor, data = beta_plot_data)
    }, error = function(e) {
      cat("  Error in beta Kruskal-Wallis test for", metric_name, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(kw_beta)) next
    
    if (kw_beta$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
      # Perform pairwise Wilcoxon tests
      pwc_beta <- pairwise.wilcox.test(beta_plot_data$Distance_to_Centroid, 
                                       group_factor, 
                                       p.adjust.method = "fdr", exact = FALSE)
      
      # Extract significant comparisons
      group_levels <- levels(group_factor)
      
      for (i in 1:(length(group_levels)-1)) {
        for (j in (i+1):length(group_levels)) {
          p_val <- pwc_beta$p.value[j-1, i]
          
          if (!is.na(p_val)) {
            beta_pairwise_results <- rbind(beta_pairwise_results, data.frame(
              Variable = cat_var,
              Metric = metric_name,
              Group1 = group_levels[i],
              Group2 = group_levels[j],
              P_value = p_val,
              P_value_FDR = p_val,
              Significant = p_val < 0.05,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
}

cat("  ✓ Statistical tests complete\n")

###############################################################################
## 5. CREATE INTEGRATED 4×3 PLOTS                                           ##
###############################################################################

cat("\nCreating integrated 4×3 plots...\n")

for (cat_var in categorical_vars) {
  
  plot_data_alpha <- alpha_with_categories %>% filter(!is.na(!!sym(cat_var)))
  if (nrow(plot_data_alpha) < 30) next
  
  n_levels <- nlevels(plot_data_alpha[[cat_var]])
  plot_colors <- safe_colors[1:min(n_levels, length(safe_colors))]
  
  plot_list <- list()
  
  # ALPHA PLOTS (top row) - with comprehensive statistics
  for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
    
    # Perform Kruskal-Wallis test with error handling
    tryCatch({
      kw_test <- kruskal.test(as.formula(paste(metric_col, "~", cat_var)), data = plot_data_alpha)
    }, error = function(e) {
      cat("  Error in plotting Kruskal-Wallis test for", metric_col, ":", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(kw_test)) {
      # Create empty plot if test fails
      p <- ggplot() + theme_void() + labs(title = paste("Error:", metric_col))
      plot_list[[length(plot_list) + 1]] <- p
      next
    }
    
    # Calculate effect size (epsilon-squared for Kruskal-Wallis)
    n <- nrow(plot_data_alpha)
    k <- nlevels(plot_data_alpha[[cat_var]])
    epsilon_sq <- (kw_test$statistic - k + 1) / (n - k)
    
    # Create statistical subtitle
    stat_subtitle <- paste0(
      "K-W: P=", format.pval(kw_test$p.value, digits=3, eps=0.001),
      ", R²=", format(epsilon_sq, digits=3, nsmall=3),
      "\nH=", format(kw_test$statistic, digits=3, nsmall=2),
      ", df=", kw_test$parameter,
      ", N=", n
    )
    
    # Base plot - Half violin with box+jitter on left
    if (HALF_PLOTS) {
      p <- ggplot(plot_data_alpha, aes(x = !!sym(cat_var), y = !!sym(metric_col), fill = !!sym(cat_var))) +
        ggdist::stat_halfeye(
          adjust = .5, 
          width = .6, 
          .width = 0, 
          justification = -.3, 
          point_colour = NA,
          alpha = 0.5
        ) +
        geom_boxplot(
          width = .25, 
          outlier.shape = NA,
          color = "darkred", 
          fatten = 2,
          alpha = 0.2
        ) +
        geom_point(
          size = 0.1,
          alpha = 0.05,  # Twice more transparent than original 0.1
          position = position_jitter(
            seed = 1, width = .1
          )
        ) +
        scale_fill_manual(values = plot_colors) +
        labs(title = gsub("_", " ", metric_col), 
             subtitle = stat_subtitle,
             x = "", y = "") +
        egg::theme_article(base_size = 5) +
        theme(
          plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 3.5),
          axis.text.y = element_text(size = 3.5),
          legend.position = "none"
        )
      
      # Add pairwise comparisons if significant
      if (kw_test$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
        tryCatch({
          # Perform pairwise Wilcoxon tests
          pwc <- pairwise.wilcox.test(plot_data_alpha[[metric_col]], 
                                     plot_data_alpha[[cat_var]], 
                                     p.adjust.method = "fdr", exact = FALSE)
          
          # Add significance brackets using ggsignif
          group_levels <- levels(plot_data_alpha[[cat_var]])
          if (length(group_levels) == 2) {
            # For 2 groups, add single comparison
            p_val <- pwc$p.value[1, 1]
            if (!is.na(p_val) && p_val < 0.05) {
              p <- p + geom_signif(
                comparisons = list(c(group_levels[1], group_levels[2])),
                annotations = if(NO_TICK_AND_BRACKET_ASTERISK) {
                  if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
                } else {
                  format.pval(p_val, digits = 3, eps = 0.001)
                },
                y_position = max(plot_data_alpha[[metric_col]], na.rm = TRUE) * 1.1,
                size = 2,
                textsize = 3
              )
            }
          }
        }, error = function(e) {
          cat("  Error adding pairwise comparisons for", metric_col, ":", e$message, "\n")
        })
      }
    } else {
      # Original full violin plot
      p <- ggplot(plot_data_alpha, aes(x = !!sym(cat_var), y = !!sym(metric_col), fill = !!sym(cat_var))) +
        geom_violin(alpha = 0.5, trim = FALSE) +
        geom_jitter(alpha = 0.1, size = 0.1, width = 0.1) +
        geom_boxplot(width = 0.2, alpha = 0.2, outlier.shape = NA,
                     color = "darkred", fatten = 2) +
        scale_fill_manual(values = plot_colors) +
        labs(title = gsub("_", " ", metric_col), 
             subtitle = stat_subtitle,
             x = "", y = "") +
        egg::theme_article(base_size = 5) +
        theme(
          plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 3.5),
          axis.text.y = element_text(size = 3.5),
          legend.position = "none"
        )
      
      # Add pairwise comparisons if significant
      if (kw_test$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
        tryCatch({
          # Perform pairwise Wilcoxon tests
          pwc <- pairwise.wilcox.test(plot_data_alpha[[metric_col]], 
                                     plot_data_alpha[[cat_var]], 
                                     p.adjust.method = "fdr", exact = FALSE)
          
          # Add significance brackets using ggsignif
          group_levels <- levels(plot_data_alpha[[cat_var]])
          if (length(group_levels) == 2) {
            # For 2 groups, add single comparison
            p_val <- pwc$p.value[1, 1]
            if (!is.na(p_val) && p_val < 0.05) {
              p <- p + geom_signif(
                comparisons = list(c(group_levels[1], group_levels[2])),
                annotations = if(NO_TICK_AND_BRACKET_ASTERISK) {
                  if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
                } else {
                  format.pval(p_val, digits = 3, eps = 0.001)
                },
                y_position = max(plot_data_alpha[[metric_col]], na.rm = TRUE) * 1.1,
                size = 2,
                textsize = 3
              )
            }
          }
        }, error = function(e) {
          cat("  Error adding pairwise comparisons for", metric_col, ":", e$message, "\n")
        })
      }
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # BETA PLOTS (second row) - MATCHING colors with alpha
  for (metric_name in names(beta_metrics)) {
    beta_plot_data <- centroid_data_all %>% filter(Metric == metric_name, Variable == cat_var)
    
    if (nrow(beta_plot_data) < 30) {
      plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
      next
    }
    
    # Ensure Group is a factor with same levels as categorical variable
    beta_plot_data$Group <- factor(beta_plot_data$Group, 
                                    levels = levels(plot_data_alpha[[cat_var]]))
    
    # Perform Kruskal-Wallis test
    group_factor <- factor(beta_plot_data$Group, levels = levels(plot_data_alpha[[cat_var]]))
    kw_beta <- kruskal.test(Distance_to_Centroid ~ group_factor, data = beta_plot_data)
    
    # Calculate effect size
    n <- nrow(beta_plot_data)
    k <- nlevels(group_factor)
    epsilon_sq <- (kw_beta$statistic - k + 1) / (n - k)
    
    stat_subtitle <- paste0(
      "K-W: P=", format.pval(kw_beta$p.value, digits=3, eps=0.001),
      ", R²=", format(epsilon_sq, digits=3, nsmall=3),
      "\nH=", format(kw_beta$statistic, digits=3, nsmall=2),
      ", df=", kw_beta$parameter,
      ", N=", n
    )
    
    # Beta diversity plot - Half violin with box+jitter on left
    if (HALF_PLOTS) {
      p <- ggplot(beta_plot_data, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
        ggdist::stat_halfeye(
          adjust = .5, 
          width = .6, 
          .width = 0, 
          justification = -.3, 
          point_colour = NA,
          alpha = 0.5
        ) +
        geom_boxplot(
          width = .25, 
          outlier.shape = NA,
          color = "darkred", 
          fatten = 2,
          alpha = 0.2
        ) +
        geom_point(
          size = 0.1,
          alpha = 0.05,  # Twice more transparent than original 0.1
          position = position_jitter(
            seed = 1, width = .1
          )
        ) +
        scale_fill_manual(values = plot_colors) +
        labs(title = paste(metric_name, "Centroid"), 
             subtitle = stat_subtitle,
             x = "", y = "") +
        egg::theme_article(base_size = 5) +
        theme(
          plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 3.5),
          axis.text.y = element_text(size = 3.5),
          legend.position = "none"
        )
      
      # Add pairwise comparisons if significant
      if (kw_beta$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
        tryCatch({
          # Perform pairwise Wilcoxon tests
          pwc <- pairwise.wilcox.test(beta_plot_data$Distance_to_Centroid, 
                                     beta_plot_data$Group, 
                                     p.adjust.method = "fdr", exact = FALSE)
          
          # Add significance brackets using ggsignif
          group_levels <- levels(beta_plot_data$Group)
          if (length(group_levels) == 2) {
            # For 2 groups, add single comparison
            p_val <- pwc$p.value[1, 1]
            if (!is.na(p_val) && p_val < 0.05) {
              p <- p + geom_signif(
                comparisons = list(c(group_levels[1], group_levels[2])),
                annotations = if(NO_TICK_AND_BRACKET_ASTERISK) {
                  if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
                } else {
                  format.pval(p_val, digits = 3, eps = 0.001)
                },
                y_position = max(beta_plot_data$Distance_to_Centroid, na.rm = TRUE) * 1.1,
                size = 2,
                textsize = 3
              )
            }
          }
        }, error = function(e) {
          cat("  Error adding pairwise comparisons for", metric_name, ":", e$message, "\n")
        })
      }
    } else {
      # Original full violin plot
      p <- ggplot(beta_plot_data, aes(x = Group, y = Distance_to_Centroid, fill = Group)) +
        geom_violin(alpha = 0.5, trim = FALSE) +
        geom_jitter(alpha = 0.1, size = 0.1, width = 0.1) +
        geom_boxplot(width = 0.2, alpha = 0.2, outlier.shape = NA,
                     color = "darkred", fatten = 2) +
        scale_fill_manual(values = plot_colors) +
        labs(title = paste(metric_name, "Centroid"), 
             subtitle = stat_subtitle,
             x = "", y = "") +
        egg::theme_article(base_size = 5) +
        theme(
          plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 3, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 3.5),
          axis.text.y = element_text(size = 3.5),
          legend.position = "none"
        )
      
      # Add pairwise comparisons if significant
      if (kw_beta$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
        tryCatch({
          # Perform pairwise Wilcoxon tests
          pwc <- pairwise.wilcox.test(beta_plot_data$Distance_to_Centroid, 
                                     beta_plot_data$Group, 
                                     p.adjust.method = "fdr", exact = FALSE)
          
          # Add significance brackets using ggsignif
          group_levels <- levels(beta_plot_data$Group)
          if (length(group_levels) == 2) {
            # For 2 groups, add single comparison
            p_val <- pwc$p.value[1, 1]
            if (!is.na(p_val) && p_val < 0.05) {
              p <- p + geom_signif(
                comparisons = list(c(group_levels[1], group_levels[2])),
                annotations = if(NO_TICK_AND_BRACKET_ASTERISK) {
                  if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns"
                } else {
                  format.pval(p_val, digits = 3, eps = 0.001)
                },
                y_position = max(beta_plot_data$Distance_to_Centroid, na.rm = TRUE) * 1.1,
                size = 2,
                textsize = 3
              )
            }
          }
        }, error = function(e) {
          cat("  Error adding pairwise comparisons for", metric_name, ":", e$message, "\n")
        })
      }
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # PCoA PLOTS (third row)
  for (metric_name in names(beta_metrics)) {
    pcoa_data <- pcoa_results[[metric_name]]$data %>% filter(!is.na(!!sym(cat_var)))
    
    if (nrow(pcoa_data) < 30) {
      plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
      next
    }
    
    var_exp <- pcoa_results[[metric_name]]$var_explained
    
    # Calculate centroids
    centroids <- pcoa_data %>%
      group_by(!!sym(cat_var)) %>%
      summarise(PC1_center = mean(PC1, na.rm=TRUE), PC2_center = mean(PC2, na.rm=TRUE), .groups = "drop")
    
    # Create PCoA plot with proper 95% confidence ellipses and 50% more transparent points
    p_pcoa <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = !!sym(cat_var), fill = !!sym(cat_var))) +
      geom_point(alpha = 0.15, size = 0.2) +  # 50% more transparent (0.3 * 0.5 = 0.15)
      stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.3, level = 0.90) +
      geom_point(data = centroids, aes(x = PC1_center, y = PC2_center),
                 size = 1, shape = 3, stroke = 0.5, show.legend = FALSE) +
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      labs(
        title = paste(metric_name, "PCoA"),
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      ) +
      egg::theme_article(base_size = 5) +
      theme(
        plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 3.5),
        axis.text = element_text(size = 3),
        legend.position = "none"
      )
    
    plot_list[[length(plot_list) + 1]] <- p_pcoa
  }
  
  # SCREE PLOTS (fourth row) - Grayscale
  for (metric_name in names(beta_metrics)) {
    var_exp <- pcoa_results[[metric_name]]$var_explained[1:5]
    
    scree_data <- data.frame(
      Axis = paste0("PC", 1:length(var_exp)),
      Variance = var_exp
    )
    scree_data$Axis <- factor(scree_data$Axis, levels = scree_data$Axis)
    
    p_scree <- ggplot(scree_data, aes(x = Axis, y = Variance)) +
      geom_bar(stat = "identity", fill = "gray50", alpha = 0.7) +
      geom_line(aes(group = 1), color = "gray20", linewidth = 0.5) +
      geom_point(color = "gray20", size = 1) +
      labs(
        title = paste(metric_name, "Scree"),
        x = "",
        y = "Variance (%)"
      ) +
      egg::theme_article(base_size = 5) +
      theme(
        plot.title = element_text(size = 5, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 3),
        axis.text.y = element_text(size = 3),
        axis.title.y = element_text(size = 3.5)
      )
    
    plot_list[[length(plot_list) + 1]] <- p_scree
  }
  
  # Create legend
  legend_plot <- ggplot(plot_data_alpha, aes(x = !!sym(cat_var), y = Observed_OTUs, fill = !!sym(cat_var))) +
    geom_violin() +
    scale_fill_manual(values = plot_colors, name = gsub("_", " ", cat_var)) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 4),
          legend.title = element_text(size = 5, face = "bold"),
          legend.key.size = unit(0.3, "cm"))
  
  legend <- get_legend(legend_plot)
  
  # Combine plots in 4×3 grid
  plots_grid <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 4,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), 
                         label_size = 6)
  
  # Add legend at bottom
  combined <- plot_grid(plots_grid, legend, ncol = 1, rel_heights = c(1, 0.05))
  
  # Save with exact dimensions
  pdf_width <- 80 / 25.4    # Convert mm to inches
  pdf_height <- if(NO_TICK_AND_BRACKET) (135 * 0.9) / 25.4 else 135 / 25.4
  
  # Add suffix based on mode
  if (NO_TICK_AND_BRACKET && !NO_TICK_AND_BRACKET_ASTERISK) {
    filename <- paste0("Integrated_", cat_var, "_no_tick_and_bracket.pdf")
  } else if (NO_TICK_AND_BRACKET_ASTERISK && !NO_TICK_AND_BRACKET) {
    filename <- paste0("Integrated_", cat_var, "_no_tick_and_bracket_asterisk.pdf")
  } else if (NO_TICK_AND_BRACKET && NO_TICK_AND_BRACKET_ASTERISK) {
    filename <- paste0("Integrated_", cat_var, "_no_tick_and_bracket_asterisk.pdf")
  } else {
    filename <- paste0("Integrated_", cat_var, ".pdf")
  }
  
  ggsave(file.path(plots_path, filename), 
         combined, width = pdf_width, height = pdf_height, dpi = 300)
  
  dimensions <- if(NO_TICK_AND_BRACKET) "80×121.5mm" else "80×135mm"
  cat("  ✓ ", filename, " (4×3 layout, ", dimensions, ")\n", sep = "")
}

###############################################################################
## 6. SAVE STATISTICAL RESULTS                                              ##
###############################################################################

cat("\nSaving statistical results...\n")

# Save alpha diversity pairwise results
if (nrow(alpha_pairwise_results) > 0) {
  write.csv(alpha_pairwise_results, 
            file.path(output_path, "alpha_pairwise_stat_res.csv"), 
            row.names = FALSE)
  cat("  ✓ Saved alpha_pairwise_stat_res.csv (", nrow(alpha_pairwise_results), " comparisons)\n", sep = "")
  cat("    Significant comparisons:", sum(alpha_pairwise_results$Significant, na.rm = TRUE), "\n")
}

# Save beta diversity centroid pairwise results
if (nrow(beta_pairwise_results) > 0) {
  write.csv(beta_pairwise_results, 
            file.path(output_path, "beta_centroid_pairwise_stat_res.csv"), 
            row.names = FALSE)
  cat("  ✓ Saved beta_centroid_pairwise_stat_res.csv (", nrow(beta_pairwise_results), " comparisons)\n", sep = "")
  cat("    Significant comparisons:", sum(beta_pairwise_results$Significant, na.rm = TRUE), "\n")
}

cat("\n✓ Complete!\n")
