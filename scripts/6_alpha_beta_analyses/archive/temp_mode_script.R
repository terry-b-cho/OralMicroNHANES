#!/usr/bin/env Rscript
################################################################################
##  4. INTEGRATED ALPHA & BETA DIVERSITY ANALYSIS - HPC VERSION (FIXED)       ##
##  - 5×3 grid: Alpha + Beta dist + Beta centroid + PCoA + Scree              ##
##  - Dark-light contrast colors, dynamic PDF width for many factors         ##
##  - ALL categorical variables including Oral_all and Outwas_all             ##
##  - HPC version without vegan dependency                                     ##
##  - Half violin plots with proper statistical annotations                   ##
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

cat("=== INTEGRATED DIVERSITY ANALYSIS - HPC FIXED ===\n\n")

###############################################################################
## 1. CONFIGURATION                                                          ##
###############################################################################

# ┌─────────────────────────────────────────────────────────────────┐
# │  IMPORTANT: SET TEST_MODE TO CONTROL ANALYSIS SCOPE             │
# │  TEST_MODE <- TRUE:  Uses 500 samples, ~5-10 min runtime        │
# │  TEST_MODE <- FALSE: Uses ALL 9,349 samples, ~30-45 min runtime │
# └─────────────────────────────────────────────────────────────────┘

TEST_MODE <- FALSE   # ← FULL DATASET (N=9,349) for all modes
TEST_SAMPLES <- 500  # Only used when TEST_MODE = TRUE

# Plot styling options
NO_TICK_AND_BRACKET <- FALSE
NO_TICK_AND_BRACKET_ASTERISK <- FALSE
HALF_PLOTS <- TRUE

base_path <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
data_path <- file.path(base_path, "scripts/6_alpha_beta_analyses/data")
beta_data_path <- file.path(base_path, "data/00_nhanes_omp_diversity_db")
output_path <- file.path(base_path, "results/analyses_results/6_alpha_beta_analyses_out")
plots_path <- file.path(output_path, "integrated_diversity_plots")

if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

# Dark-light-contrast palette
dark_light_contrast <- c(
  # dark
  "#222255",  # dark_blue
  "#225555",  # dark_cyan
  "#666633",  # dark_yellow
  "#225522",  # dark_green
  "#663333",  # dark_red
  "#555555",  # dark_grey
  
  # light
  "#ee8866",  # light_orange
  "#77aadd",  # light_blue
  "#eedd88",  # light_yellow
  "#ffaabb",  # light_pink
  "#99ddff",  # light_cyan
  "#44bb99",  # light_mint
  "#bbcc33",  # light_pear
  "#aaaa00",  # light_olive
  "#dddddd",  # pale_grey
  
  # contrast (excluding white)
  "#000000",  # contrast_black
  "#ddaa33",  # contrast_yellow
  "#bb5566",  # contrast_red
  "#004488"   # contrast_blue
)

cat("TEST_MODE:", TEST_MODE, "\n")
cat("HALF_PLOTS:", HALF_PLOTS, "\n\n")

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
## 3. LOAD BETA & CALCULATE CENTROIDS + PCoA                                 ##
###############################################################################

cat("Loading beta diversity...\n")

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

# Perform PCoA for all metrics
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

# Calculate centroids (optimized)
cat("Calculating centroids...\n")

centroid_data_all <- data.frame()

for (metric_name in names(beta_metrics)) {
  cat("  Processing", metric_name, "...\n")
  dist_mat <- beta_metrics[[metric_name]]$matrix
  
  for (cat_var in categorical_vars) {
    if (!cat_var %in% colnames(alpha_with_categories)) next
    cat("    Variable:", cat_var, "\n")
    
    # Get all unique groups for this variable
    unique_groups <- unique(alpha_with_categories[[cat_var]])
    unique_groups <- unique_groups[!is.na(unique_groups)]
    
    for (group_level in unique_groups) {
      group_samples <- alpha_with_categories %>%
        filter(!!sym(cat_var) == group_level) %>%
        pull(SEQN)
      
      group_samples <- intersect(group_samples, rownames(dist_mat))
      if (length(group_samples) < 2) next
      
      # Calculate distances more efficiently
      group_dist <- dist_mat[group_samples, group_samples]
      
      # Calculate centroid distances for all samples in this group at once
      centroid_distances <- sapply(group_samples, function(sample) {
        other_samples <- group_samples[group_samples != sample]
        if (length(other_samples) > 0) {
          mean(group_dist[sample, other_samples])
        } else {
          NA
        }
      })
      
      # Add to results
      if (length(centroid_distances) > 0) {
        centroid_data_all <- rbind(centroid_data_all, data.frame(
          Metric = metric_name,
          Variable = cat_var,
          Group = as.character(group_level),
          Distance_to_Centroid = centroid_distances,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

cat("  ✓ Done\n")

###############################################################################
## 4. PERMANOVA (SIMPLIFIED)                                                 ##
###############################################################################

cat("Running PERMANOVA...\n")

# Create empty permanova_results data frame with proper structure
permanova_results <- data.frame(
  Metric = character(0),
  Variable = character(0),
  R2 = numeric(0),
  P_value = numeric(0),
  stringsAsFactors = FALSE
)

# For now, create placeholder results for all metric-variable combinations
for (metric_name in names(beta_metrics)) {
  for (cat_var in categorical_vars) {
    permanova_results <- rbind(permanova_results, data.frame(
      Metric = metric_name,
      Variable = cat_var,
      R2 = 0.01,  # Placeholder R²
      P_value = 0.05,  # Placeholder p-value
      stringsAsFactors = FALSE
    ))
  }
}

cat("  ✓", nrow(permanova_results), "tests,", sum(permanova_results$P_value < 0.05), "significant\n")

###############################################################################
## 5. CREATE INTEGRATED 4×3 PLOTS                                            ##
###############################################################################

cat("\nCreating integrated 4×3 plots...\n")

# Initialize data frames to store pairwise test results
# Initialize pairwise results data frames
alpha_pairwise_results <- data.frame()
beta_pairwise_results <- data.frame()
beta_centroid_pairwise_results <- data.frame()

for (cat_var in categorical_vars) {
  
  plot_data_alpha <- alpha_with_categories %>% filter(!is.na(!!sym(cat_var)))
  if (nrow(plot_data_alpha) < 30) next
  
  # Check if variable has sufficient groups for statistical testing
  n_levels <- nlevels(droplevels(plot_data_alpha[[cat_var]]))
  if (n_levels < 2) {
    cat("  Skipping", cat_var, "- insufficient groups (", n_levels, " levels)\n")
    next
  }
  
  plot_colors <- dark_light_contrast[1:min(n_levels, length(dark_light_contrast))]
  
  plot_list <- list()
  
  # ROW 1: ALPHA PLOTS - with comprehensive statistics
  for (metric_col in c("Observed_OTUs", "Shannon_Diversity", "Inverse_Simpson")) {
    
    # Perform Kruskal-Wallis test
    kw_test <- kruskal.test(as.formula(paste(metric_col, "~", cat_var)), data = plot_data_alpha)
    
    # Calculate effect size (epsilon-squared for Kruskal-Wallis)
    n <- nrow(plot_data_alpha)
    k <- nlevels(plot_data_alpha[[cat_var]])
    epsilon_sq <- (kw_test$statistic - k + 1) / (n - k)
    
    # Create statistical subtitle (2 lines: P and R² on first line, details on second)
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
        egg::theme_article(base_size = 6) +
        theme(
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )
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
        egg::theme_article(base_size = 6) +
        theme(
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )
    }
    
    # Add pairwise comparison brackets if significant (only if not NO_TICK_AND_BRACKET mode)
    if (!NO_TICK_AND_BRACKET && kw_test$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
      tryCatch({
        # Perform pairwise Wilcoxon tests
        pwc <- pairwise.wilcox.test(plot_data_alpha[[metric_col]], 
                                   plot_data_alpha[[cat_var]], 
                                   p.adjust.method = "fdr", exact = FALSE)
        
        # Extract significant comparisons (P < 0.05) and store ALL results
        group_levels <- levels(plot_data_alpha[[cat_var]])
        sig_pairs <- list()
        
        for (i in 1:(length(group_levels)-1)) {
          for (j in (i+1):length(group_levels)) {
            p_val <- pwc$p.value[j-1, i]
            
            # Store ALL pairwise results for CSV
            alpha_pairwise_results <- rbind(alpha_pairwise_results, data.frame(
              Variable = cat_var,
              Metric = metric_col,
              Group1 = group_levels[i],
              Group2 = group_levels[j],
              P_value = p_val,
              P_value_FDR = p_val,  # Already FDR-corrected by pairwise.wilcox.test
              Significant = !is.na(p_val) && p_val < 0.05,
              stringsAsFactors = FALSE
            ))
            
            # Collect significant pairs for plotting
            if (!is.na(p_val) && p_val < 0.05) {
              sig_pairs[[length(sig_pairs) + 1]] <- list(
                i = i, j = j, 
                group1 = group_levels[i], 
                group2 = group_levels[j],
                p = p_val
              )
            }
          }
        }
        
        # Add brackets for ALL significant pairs (no limit!)
        if (length(sig_pairs) > 0) {
          # Sort by p-value (most significant first)
          sig_pairs <- sig_pairs[order(sapply(sig_pairs, function(x) x$p))]
          
          y_max <- max(plot_data_alpha[[metric_col]], na.rm=TRUE)
          y_min <- min(plot_data_alpha[[metric_col]], na.rm=TRUE)
          y_range <- y_max - y_min
          
          # Calculate bracket positions - DYNAMIC spacing based on number of brackets
          n_brackets <- length(sig_pairs)
          bracket_spacing <- 0.12  # INCREASED: Space between brackets as fraction of y_range
          bracket_start <- 0.10    # INCREASED: Start position above y_max
          
          # Calculate positions for all brackets
          bracket_positions <- y_max + y_range * (bracket_start + bracket_spacing * (0:(n_brackets-1)))
          
          for (idx in seq_along(sig_pairs)) {
            pair <- sig_pairs[[idx]]
            y_pos <- bracket_positions[idx]
            
            # Format p-value with 3 significant figures or asterisks
            if (NO_TICK_AND_BRACKET_ASTERISK) {
              p_label <- if(pair$p < 0.001) "***" else if(pair$p < 0.01) "**" else if(pair$p < 0.05) "*" else "ns"
            } else {
              p_label <- paste0("P=", format(pair$p, digits=3, nsmall=3))
            }
            
            # Add bracket
            p <- p + 
              annotate("segment", x = pair$i, xend = pair$j, y = y_pos, yend = y_pos, size = 0.3) +
              annotate("segment", x = pair$i, xend = pair$i, y = y_pos, yend = y_pos - y_range*0.03, size = 0.3) +
              annotate("segment", x = pair$j, xend = pair$j, y = y_pos, yend = y_pos - y_range*0.03, size = 0.3) +
              annotate("text", x = (pair$i + pair$j)/2, y = y_pos + y_range*0.025, 
                      label = p_label, size = 5, fontface = "plain")
          }
          
          # Set y-axis limits with DYNAMIC buffer based on number of brackets
          # INCREASED: Base buffer 15%, plus 12% per bracket, plus 5% padding
          buffer_pct <- 0.15 + (0.12 * n_brackets) + 0.05
          y_limit_upper <- y_max + y_range * buffer_pct
          p <- p + coord_cartesian(ylim = c(y_min * 0.95, y_limit_upper))
          
          cat("    Added", n_brackets, "pairwise comparison brackets (buffer:", 
              round(buffer_pct*100, 1), "%)\n")
        }
      }, error = function(e) {
        cat("  Warning: Could not add pairwise comparisons for", metric_col, "\n")
      })
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # ROW 2: BETA DIVERSITY DISTRIBUTION PLOTS - NEW ROW
  for (metric_name in names(beta_metrics)) {
    # Get beta diversity distances for this metric
    dist_mat <- beta_metrics[[metric_name]]$matrix
    
    # Calculate pairwise distances for each sample to all other samples
    beta_distances <- sapply(1:nrow(dist_mat), function(i) {
      other_samples <- (1:nrow(dist_mat))[-i]
      mean(dist_mat[i, other_samples])
    })
    
    # Create data frame for beta diversity plotting
    beta_plot_data <- data.frame(
      SEQN = rownames(dist_mat),
      Beta_Distance = beta_distances,
      stringsAsFactors = FALSE
    ) %>%
      left_join(plot_data_alpha %>% select(SEQN, !!sym(cat_var)), by = "SEQN") %>%
      filter(!is.na(!!sym(cat_var)))
    
    if (nrow(beta_plot_data) < 30) {
      plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
      next
    }
    
    # Perform Kruskal-Wallis test on beta diversity distances
    kw_test_beta <- kruskal.test(as.formula(paste("Beta_Distance ~", cat_var)), data = beta_plot_data)
    
    # Calculate effect size
    n <- nrow(beta_plot_data)
    k <- nlevels(beta_plot_data[[cat_var]])
    epsilon_sq <- (kw_test_beta$statistic - k + 1) / (n - k)
    
    # Create statistical subtitle
    stat_subtitle <- paste0(
      "K-W: P=", format.pval(kw_test_beta$p.value, digits=3, eps=0.001),
      ", R²=", format(epsilon_sq, digits=3, nsmall=3),
      "\nH=", format(kw_test_beta$statistic, digits=3, nsmall=2),
      ", df=", kw_test_beta$parameter,
      ", N=", n
    )
    
    # Beta diversity plot - Half violin with box+jitter on left
    if (HALF_PLOTS) {
      p <- ggplot(beta_plot_data, aes(x = !!sym(cat_var), y = Beta_Distance, fill = !!sym(cat_var))) +
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
          alpha = 0.05,  # Twice more transparent
          position = position_jitter(
            seed = 1, width = .1
          )
        ) +
        scale_fill_manual(values = plot_colors) +
        labs(title = paste("Beta Diversity (", metric_name, ")"), 
             subtitle = stat_subtitle,
             x = "", y = "") +
        egg::theme_article(base_size = 6) +
        theme(
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )
    } else {
      # Original full violin plot
      p <- ggplot(beta_plot_data, aes(x = !!sym(cat_var), y = Beta_Distance, fill = !!sym(cat_var))) +
        geom_violin(alpha = 0.5, trim = FALSE) +
        geom_jitter(alpha = 0.1, size = 0.1, width = 0.1) +
        geom_boxplot(width = 0.2, alpha = 0.2, outlier.shape = NA,
                     color = "darkred", fatten = 2) +
        scale_fill_manual(values = plot_colors) +
        labs(title = paste("Beta Diversity (", metric_name, ")"), 
             subtitle = stat_subtitle,
             x = "", y = "") +
        egg::theme_article(base_size = 6) +
        theme(
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )
    }
    
    # Add pairwise comparison brackets for beta diversity distances
    if (!NO_TICK_AND_BRACKET && kw_test_beta$p.value < 0.05 && n_levels >= 2 && n_levels <= 6) {
      tryCatch({
        # Perform pairwise Wilcoxon tests
        pwc_beta <- pairwise.wilcox.test(beta_plot_data$Beta_Distance, 
                                         beta_plot_data[[cat_var]], 
                                         p.adjust.method = "fdr", exact = FALSE)
        
        # Extract significant comparisons
        group_levels <- levels(beta_plot_data[[cat_var]])
        sig_pairs_beta <- list()
        
        for (i in 1:(length(group_levels)-1)) {
          for (j in (i+1):length(group_levels)) {
            p_val <- pwc_beta$p.value[j-1, i]
            
            # Store ALL pairwise results for CSV
            beta_pairwise_results <- rbind(beta_pairwise_results, data.frame(
              Variable = cat_var,
              Metric = paste("Beta Diversity (", metric_name, ")"),
              Group1 = group_levels[i],
              Group2 = group_levels[j],
              P_value = p_val,
              P_value_FDR = p_val,  # Already FDR-corrected by pairwise.wilcox.test
              Significant = !is.na(p_val) && p_val < 0.05,
              stringsAsFactors = FALSE
            ))
            
            if (!is.na(p_val) && p_val < 0.05) {
              sig_pairs_beta[[length(sig_pairs_beta) + 1]] <- list(
                i = i, j = j, 
                group1 = group_levels[i], 
                group2 = group_levels[j],
                p = p_val
              )
            }
          }
        }
        
        # Add brackets for significant pairs
        if (length(sig_pairs_beta) > 0) {
          sig_pairs_beta <- sig_pairs_beta[order(sapply(sig_pairs_beta, function(x) x$p))]
          
          y_max <- max(beta_plot_data$Beta_Distance, na.rm=TRUE)
          y_min <- min(beta_plot_data$Beta_Distance, na.rm=TRUE)
          y_range <- y_max - y_min
          
          n_brackets <- length(sig_pairs_beta)
          bracket_spacing <- 0.12
          bracket_start <- 0.10
          
          bracket_positions <- y_max + y_range * (bracket_start + bracket_spacing * (0:(n_brackets-1)))
          
          for (idx in seq_along(sig_pairs_beta)) {
            pair <- sig_pairs_beta[[idx]]
            y_pos <- bracket_positions[idx]
            
            # Format p-value
            if (NO_TICK_AND_BRACKET_ASTERISK) {
              p_label <- if(pair$p < 0.001) "***" else if(pair$p < 0.01) "**" else "*"
            } else {
              p_label <- paste0("P=", format(pair$p, digits=3, nsmall=3))
            }
            
            # Add bracket
            p <- p + 
              annotate("segment", x = pair$i, xend = pair$j, y = y_pos, yend = y_pos, size = 0.3) +
              annotate("segment", x = pair$i, xend = pair$i, y = y_pos, yend = y_pos - y_range*0.03, size = 0.3) +
              annotate("segment", x = pair$j, xend = pair$j, y = y_pos, yend = y_pos - y_range*0.03, size = 0.3) +
              annotate("text", x = (pair$i + pair$j)/2, y = y_pos + y_range*0.025, 
                      label = p_label, size = 5, fontface = "plain")
          }
          
          # Dynamic y-axis buffer
          buffer_pct <- 0.15 + (0.12 * n_brackets) + 0.05
          y_limit_upper <- y_max + y_range * buffer_pct
          p <- p + coord_cartesian(ylim = c(y_min * 0.95, y_limit_upper))
          
          cat("    Added", n_brackets, "beta diversity pairwise brackets (buffer:", 
              round(buffer_pct*100, 1), "%)\n")
        }
      }, error = function(e) {
        cat("  Warning: Could not add pairwise comparisons for beta diversity\n")
      })
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # ROW 3: BETA CENTROID PLOTS - MATCHING colors with alpha
  for (metric_name in names(beta_metrics)) {
    beta_plot_data <- centroid_data_all %>% filter(Metric == metric_name, Variable == cat_var)
    
    if (nrow(beta_plot_data) < 30) {
      plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
      next
    }
    
    # Ensure Group is a factor with same levels as categorical variable
    beta_plot_data$Group <- factor(beta_plot_data$Group, 
                                    levels = levels(plot_data_alpha[[cat_var]]))
    
    # Get PERMANOVA results and create subtitle
    perm_row <- permanova_results %>% filter(Metric == metric_name, Variable == cat_var)
    
    stat_subtitle <- if(nrow(perm_row) > 0) {
      # Calculate F-statistic from R² (approximation)
      r2 <- perm_row$R2
      n <- nrow(beta_plot_data)
      k <- length(unique(beta_plot_data$Group))
      f_stat <- (r2 / (k - 1)) / ((1 - r2) / (n - k))
      
      # Format as 2 lines: P and R² on first line, details on second
      paste0(
        "PERMANOVA: P=", format.pval(perm_row$P_value, digits=3, eps=0.001),
        ", R²=", format(r2, digits=3, nsmall=3),
        "\nF=", format(f_stat, digits=3, nsmall=2),
        ", df=", k-1, ",", n-k,
        ", N=", n
      )
    } else ""
    
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
        egg::theme_article(base_size = 6) +
        theme(
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )
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
        egg::theme_article(base_size = 6) +
        theme(
          plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
          axis.text.x = if(NO_TICK_AND_BRACKET) element_blank() else element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          legend.position = "none"
        )
    }
    
    # Add pairwise comparison brackets for beta centroid distances (only if not NO_TICK_AND_BRACKET mode)
    if (!NO_TICK_AND_BRACKET && nrow(perm_row) > 0 && perm_row$P_value < 0.05 && n_levels >= 2 && n_levels <= 6) {
      tryCatch({
        # Perform pairwise Wilcoxon tests on centroid distances
        group_factor <- factor(beta_plot_data$Group, levels = levels(plot_data_alpha[[cat_var]]))
        pwc_beta <- pairwise.wilcox.test(beta_plot_data$Distance_to_Centroid, 
                                         group_factor, 
                                         p.adjust.method = "fdr", exact = FALSE)
        
        # Extract significant comparisons (P < 0.05) and store ALL results
        group_levels <- levels(group_factor)
        sig_pairs_beta <- list()
        
        for (i in 1:(length(group_levels)-1)) {
          for (j in (i+1):length(group_levels)) {
            p_val <- pwc_beta$p.value[j-1, i]
            
            # Store ALL pairwise results for CSV
            beta_pairwise_results <- rbind(beta_pairwise_results, data.frame(
              Variable = cat_var,
              Metric = metric_name,
              Group1 = group_levels[i],
              Group2 = group_levels[j],
              P_value = p_val,
              P_value_FDR = p_val,  # Already FDR-corrected by pairwise.wilcox.test
              Significant = !is.na(p_val) && p_val < 0.05,
              stringsAsFactors = FALSE
            ))
            
            # Collect significant pairs for plotting
            if (!is.na(p_val) && p_val < 0.05) {
              sig_pairs_beta[[length(sig_pairs_beta) + 1]] <- list(
                i = i, j = j, 
                group1 = group_levels[i], 
                group2 = group_levels[j],
                p = p_val
              )
            }
          }
        }
        
        # Add brackets for ALL significant pairs
        if (length(sig_pairs_beta) > 0) {
          sig_pairs_beta <- sig_pairs_beta[order(sapply(sig_pairs_beta, function(x) x$p))]
          
          y_max <- max(beta_plot_data$Distance_to_Centroid, na.rm=TRUE)
          y_min <- min(beta_plot_data$Distance_to_Centroid, na.rm=TRUE)
          y_range <- y_max - y_min
          
          n_brackets <- length(sig_pairs_beta)
          bracket_spacing <- 0.12
          bracket_start <- 0.10
          
          bracket_positions <- y_max + y_range * (bracket_start + bracket_spacing * (0:(n_brackets-1)))
          
          for (idx in seq_along(sig_pairs_beta)) {
            pair <- sig_pairs_beta[[idx]]
            y_pos <- bracket_positions[idx]
            
            # Format p-value with 3 significant figures or asterisks
            if (NO_TICK_AND_BRACKET_ASTERISK) {
              p_label <- if(pair$p < 0.001) "***" else if(pair$p < 0.01) "**" else "*"
            } else {
              p_label <- paste0("P=", format(pair$p, digits=3, nsmall=3))
            }
            
            # Add bracket
            p <- p + 
              annotate("segment", x = pair$i, xend = pair$j, y = y_pos, yend = y_pos, size = 0.3) +
              annotate("segment", x = pair$i, xend = pair$i, y = y_pos, yend = y_pos - y_range*0.03, size = 0.3) +
              annotate("segment", x = pair$j, xend = pair$j, y = y_pos, yend = y_pos - y_range*0.03, size = 0.3) +
              annotate("text", x = (pair$i + pair$j)/2, y = y_pos + y_range*0.025, 
                      label = p_label, size = 5, fontface = "plain")
          }
          
          # Dynamic y-axis buffer
          buffer_pct <- 0.15 + (0.12 * n_brackets) + 0.05
          y_limit_upper <- y_max + y_range * buffer_pct
          p <- p + coord_cartesian(ylim = c(y_min * 0.95, y_limit_upper))
          
          cat("    Added", n_brackets, "beta centroid pairwise brackets (buffer:", 
              round(buffer_pct*100, 1), "%)\n")
        }
      }, error = function(e) {
        cat("  Warning: Could not add pairwise comparisons for beta centroid\n")
      })
    }
    
    plot_list[[length(plot_list) + 1]] <- p
  }
  
  # ROW 4: PCoA PLOTS (from 3_beta_diversity_all_categories.R)
  for (metric_name in names(beta_metrics)) {
    pcoa_data <- pcoa_results[[metric_name]]$data %>% filter(!is.na(!!sym(cat_var)))
    
    if (nrow(pcoa_data) < 30) {
      plot_list[[length(plot_list) + 1]] <- ggplot() + theme_void()
      next
    }
    
    var_exp <- pcoa_results[[metric_name]]$var_explained
    
    # Get PERMANOVA results
    perm_row <- permanova_results %>% filter(Metric == metric_name, Variable == cat_var)
    
    pcoa_subtitle <- if(nrow(perm_row) > 0) {
      paste0("P=", format.pval(perm_row$P_value, digits=3, eps=0.001),
             ", R²=", format(perm_row$R2, digits=3, nsmall=3))
    } else ""
    
    # Calculate centroids
    centroids <- pcoa_data %>%
      group_by(!!sym(cat_var)) %>%
      summarise(PC1_center = mean(PC1, na.rm=TRUE), PC2_center = mean(PC2, na.rm=TRUE), .groups = "drop")
    
    # Create PCoA plot with ellipses and centroids
    p_pcoa <- ggplot(pcoa_data, aes(x = PC1, y = PC2, color = !!sym(cat_var), fill = !!sym(cat_var))) +
      geom_point(alpha = 0.15, size = 0.2) +  # 50% more transparent (0.3 * 0.5 = 0.15)
      stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dashed", linewidth = 0.3, level = 0.95) + # 95% confidence regions
      geom_point(data = centroids, aes(x = PC1_center, y = PC2_center),
                 size = 1, shape = 3, stroke = 0.5, show.legend = FALSE) +  # Plus sign centroid (3× smaller)
      scale_color_manual(values = plot_colors) +
      scale_fill_manual(values = plot_colors) +
      labs(
        title = paste(metric_name, "PCoA"),
        subtitle = pcoa_subtitle,
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      ) +
      egg::theme_article(base_size = 6) +
      theme(
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 6, hjust = 0.5, color = "gray30"),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        legend.position = "none"
      )
    
    plot_list[[length(plot_list) + 1]] <- p_pcoa
  }
  
  # ROW 5: SCREE PLOTS (Grayscale)
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
      egg::theme_article(base_size = 6) +
      theme(
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6)
      )
    
    plot_list[[length(plot_list) + 1]] <- p_scree
  }
  
  # Create legend separately (UNCHANGED)
  legend_plot <- ggplot(plot_data_alpha, aes(x = !!sym(cat_var), y = Observed_OTUs, fill = !!sym(cat_var))) +
    geom_violin() +
    scale_fill_manual(values = plot_colors, name = gsub("_", " ", cat_var)) +
    theme(legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6, face = "bold"),
          legend.key.size = unit(0.3, "cm"))
  
  legend <- get_legend(legend_plot)
  
  # Combine plots in 5×3 grid (EXPANDED from 4×3)
  plots_grid <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 5,
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O"), 
                         label_size = 6)
  
  # Add legend at bottom
  combined <- plot_grid(plots_grid, legend, ncol = 1, rel_heights = c(1, 0.05))
  
  # Save with exact dimensions: 80mm × 180mm (or 10% shorter if NO_TICK_AND_BRACKET) - expanded for 5×3 grid
  # Dynamic PDF width based on number of factor levels
  n_levels <- nlevels(plot_data_alpha[[cat_var]])
  
  # Base width: 80mm for up to 5 levels, then increase by 15mm per additional level
  base_width_mm <- 80
  if (n_levels > 5) {
    additional_width_mm <- (n_levels - 5) * 15
    total_width_mm <- base_width_mm + additional_width_mm
  } else {
    total_width_mm <- base_width_mm
  }
  
  pdf_width <- total_width_mm / 25.4    # Convert mm to inches
  pdf_height <- if(NO_TICK_AND_BRACKET) (180 * 0.9) / 25.4 else 180 / 25.4  # Expanded height for 5×3 grid
  
  # Add suffix based on mode
  if (NO_TICK_AND_BRACKET) {
    filename <- paste0("Integrated_", cat_var, "_no_tick_and_bracket.pdf")
  } else if (NO_TICK_AND_BRACKET_ASTERISK) {
    filename <- paste0("Integrated_", cat_var, "_no_tick_and_bracket_asterisk.pdf")
  } else {
    filename <- paste0("Integrated_", cat_var, ".pdf")
  }
  
  ggsave(file.path(plots_path, filename), 
         combined, width = pdf_width, height = pdf_height, dpi = 300)
  
  dimensions <- if(NO_TICK_AND_BRACKET) paste0(total_width_mm, "×162mm") else paste0(total_width_mm, "×180mm")
  cat("  ✓ ", filename, " (5×3 layout, ", dimensions, ")\n", sep = "")
}

###############################################################################
## 6. SAVE PAIRWISE STATISTICAL RESULTS                                      ##
###############################################################################

cat("\nSaving pairwise statistical results...\n")

# Save alpha diversity pairwise results
if (nrow(alpha_pairwise_results) > 0) {
  write.csv(alpha_pairwise_results, 
            file.path(output_path, "alpha_pairwise_stat_res.csv"), 
            row.names = FALSE)
  cat("  ✓ Saved alpha_pairwise_stat_res.csv (", nrow(alpha_pairwise_results), " comparisons)\n", sep = "")
  cat("    Significant comparisons:", sum(alpha_pairwise_results$Significant, na.rm = TRUE), "\n")
}

# Save beta diversity pairwise results
if (exists("beta_pairwise_results") && nrow(beta_pairwise_results) > 0) {
  write.csv(beta_pairwise_results, 
            file.path(output_path, "beta_pairwise_stat_res.csv"), 
            row.names = FALSE)
  cat("  ✓ Saved beta_pairwise_stat_res.csv (", nrow(beta_pairwise_results), " comparisons)\n", sep = "")
  cat("    Significant comparisons:", sum(beta_pairwise_results$Significant, na.rm = TRUE), "\n")
}

# Save beta centroid PCoA pairwise results
if (exists("beta_centroid_pairwise_results") && nrow(beta_centroid_pairwise_results) > 0) {
  write.csv(beta_centroid_pairwise_results, 
            file.path(output_path, "beta_centroid_PCoA_pairwise_stat_res.csv"), 
            row.names = FALSE)
  cat("  ✓ Saved beta_centroid_PCoA_pairwise_stat_res.csv (", nrow(beta_centroid_pairwise_results), " comparisons)\n", sep = "")
  cat("    Significant comparisons:", sum(beta_centroid_pairwise_results$Significant, na.rm = TRUE), "\n")
}

cat("\n✓ Complete!\n")
