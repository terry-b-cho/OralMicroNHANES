#!/usr/bin/env Rscript

# =============================================================================
# Module 10 - prediction pipeline, step 3: performance visualizations
# =============================================================================
# Reads outputs from step 2 (metrics_summary.csv and per-task prediction RDS)
# and writes one PDF per visualization under:
#   results/analyses_results/10_prediction_analyses_out/figures/
#
#   metrics_overview.pdf           grouped bar chart, all metrics x task x feature_set
#   age_pred_vs_true.pdf           regression task scatter + 1:1 line
#   gender_roc.pdf                 binary task ROC curve, both feature sets overlaid
#   age_group_confusion_matrix.pdf multiclass task heat-map confusion matrix
#
# Environment: R >= 4.5 with dplyr, tibble, tidyr, readr, ggplot2, pROC.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pROC)
})

base_path <- PROJECT_ROOT
out_root  <- file.path(base_path, "results/analyses_results/10_prediction_analyses_out")
out_pred  <- file.path(out_root, "predictions")
out_fig   <- file.path(out_root, "figures")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

metrics <- readr::read_csv(file.path(out_root, "metrics_summary.csv"), show_col_types = FALSE)

# ---- 1. Metrics overview: all metrics x task x feature_set x weighting -----
metrics_long <- metrics |>
  pivot_longer(
    cols      = -c(task, feature_set, model, weighting, n),
    names_to  = "metric",
    values_to = "value"
  ) |>
  filter(!is.na(value))

p_overview <- ggplot(metrics_long,
                     aes(x = feature_set, y = value,
                         fill = weighting)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", value)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.5) +
  facet_grid(metric ~ task, scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(unweighted = "#56B4E9", survey_weighted = "#E69F00")) +
  labs(title = "Module 10 - test-set metrics by task and feature set",
       x = NULL, y = NULL, fill = "Weighting") +
  theme_bw(base_size = 9) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(out_fig, "metrics_overview.pdf"),
       p_overview, width = 9, height = 7)
message("Wrote ", file.path(out_fig, "metrics_overview.pdf"))

# ---- 2. Age regression — predicted vs true scatter -------------------------
age_files <- list.files(out_pred, pattern = "^age__.*\\.rds$", full.names = TRUE)
if (length(age_files) > 0) {
  age_preds <- bind_rows(lapply(age_files, function(f) {
    d <- readRDS(f)
    d$feature_set <- sub("^age__", "", tools::file_path_sans_ext(basename(f)))
    d
  }))
  ax_min <- floor(min(age_preds$y_true, age_preds$y_pred,   na.rm = TRUE))
  ax_max <- ceiling(max(age_preds$y_true, age_preds$y_pred, na.rm = TRUE))
  p_age <- ggplot(age_preds, aes(x = y_true, y = y_pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.25, size = 0.6, colour = "#0072B2") +
    geom_smooth(method = "lm", se = FALSE, colour = "#D55E00", linewidth = 0.6) +
    facet_wrap(~ feature_set, ncol = 2) +
    coord_equal(xlim = c(ax_min, ax_max), ylim = c(ax_min, ax_max)) +
    labs(title = "Age regression - predicted vs observed (test set)",
         x = "Observed age (years)", y = "Predicted age (years)") +
    theme_bw(base_size = 10)
  ggsave(file.path(out_fig, "age_pred_vs_true.pdf"),
         p_age, width = 8, height = 4.5)
  message("Wrote ", file.path(out_fig, "age_pred_vs_true.pdf"))
}

# ---- 3. Gender binary classification - ROC curve ---------------------------
gender_files <- list.files(out_pred, pattern = "^gender__.*\\.rds$", full.names = TRUE)
if (length(gender_files) > 0) {
  pdf(file.path(out_fig, "gender_roc.pdf"), width = 6, height = 6)
  on.exit(dev.off(), add = TRUE)
  cols <- c(microbiome_only = "#0072B2", microbiome_plus_basic_demo = "#D55E00")
  first <- TRUE
  for (f in gender_files) {
    fs <- sub("^gender__", "", tools::file_path_sans_ext(basename(f)))
    d  <- readRDS(f)
    r  <- pROC::roc(d$y_true, d$y_pred_prob, quiet = TRUE)
    a  <- as.numeric(pROC::auc(r))
    plot_call <- if (first) {
      plot(r, col = cols[fs], lwd = 2, main = "Gender classification - ROC")
    } else {
      lines(r, col = cols[fs], lwd = 2)
    }
    if (first) {
      legend_labels <- c(sprintf("%s (AUC = %.3f)", fs, a))
      legend_colors <- cols[fs]
    } else {
      legend_labels <- c(legend_labels, sprintf("%s (AUC = %.3f)", fs, a))
      legend_colors <- c(legend_colors, cols[fs])
    }
    first <- FALSE
  }
  legend("bottomright", legend = legend_labels, col = legend_colors,
         lwd = 2, bty = "n", cex = 0.9)
  dev.off()
  message("Wrote ", file.path(out_fig, "gender_roc.pdf"))
}

# ---- 4. Age-group multiclass - confusion matrix heatmap --------------------
ag_files <- list.files(out_pred, pattern = "^age_group__.*\\.rds$", full.names = TRUE)
if (length(ag_files) > 0) {
  ag <- bind_rows(lapply(ag_files, function(f) {
    d <- readRDS(f)
    d$feature_set <- sub("^age_group__", "", tools::file_path_sans_ext(basename(f)))
    d
  }))
  age_levels <- c("14-19","20-29","30-39","40-49","50-59","60-69")
  cm <- ag |>
    count(feature_set, y_true, y_pred_class) |>
    group_by(feature_set, y_true) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    mutate(y_true       = factor(y_true,       levels = age_levels),
           y_pred_class = factor(y_pred_class, levels = age_levels))
  p_cm <- ggplot(cm, aes(x = y_pred_class, y = y_true, fill = prop)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", prop)), size = 2.6) +
    scale_fill_gradient(low = "#FFFFFF", high = "#0072B2", limits = c(0, 1)) +
    facet_wrap(~ feature_set, ncol = 2) +
    labs(title = "Age-group classification - row-normalized confusion matrix",
         x = "Predicted age group", y = "Observed age group", fill = "Row %") +
    coord_equal() +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(out_fig, "age_group_confusion_matrix.pdf"),
         p_cm, width = 9, height = 5)
  message("Wrote ", file.path(out_fig, "age_group_confusion_matrix.pdf"))
}

message("Done. Figures written to: ", out_fig)
