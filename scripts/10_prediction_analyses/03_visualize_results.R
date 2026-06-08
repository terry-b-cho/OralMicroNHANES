#!/usr/bin/env Rscript

# =============================================================================
# Module 10 - prediction pipeline, step 3: publication-grade visualizations
# =============================================================================
# Reads step-2 outputs (metrics_summary.csv + per-task prediction RDS files +
# fitted cv.glmnet model RDS files) and writes one faceted PDF per family
# under: results/analyses_results/10_prediction_analyses_out/figures/
#
#   metrics_overview.pdf            per-task tabular heat-map of test metrics
#   roc_curves_classification.pdf   ROC, gender binary + age_group one-vs-rest
#   precision_recall_curves.pdf     same as above but PR (better when imbalanced)
#   calibration_curves.pdf          decile reliability curves for all classes
#   confusion_matrices.pdf          row-normalised, per (task x feature_set)
#   feature_importance.pdf          top-20 taxa per task by |glmnet coefficient|
#   regression_diagnostics.pdf      observed-vs-predicted + residuals + Q-Q
#
# Environment: R >= 4.5 with dplyr, tibble, tidyr, readr, ggplot2, pROC, glmnet,
# patchwork. Conda spec: envs/nhanes-analysis_for_reviewers.yml.
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
  library(glmnet)
  library(patchwork)
})

base_path  <- PROJECT_ROOT
out_root   <- file.path(base_path, "results/analyses_results/10_prediction_analyses_out")
out_pred   <- file.path(out_root, "predictions")
out_models <- file.path(out_root, "models")
out_fig    <- file.path(out_root, "figures")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

metrics    <- readr::read_csv(file.path(out_root, "metrics_summary.csv"),
                              show_col_types = FALSE)
pred_files <- list.files(out_pred,   pattern = "\\.rds$", full.names = TRUE)
model_files <- list.files(out_models, pattern = "\\.rds$", full.names = TRUE)

# Feature set -> linetype (more legible than colour when many classes share lines)
fs_linetype <- c(microbiome_only = "dashed", microbiome_plus_basic_demo = "solid")

# Class -> grafify (Okabe-Ito + Tol-vibrant) so 6+ classes stay distinguishable
grafify_class <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7",
                   "#0072B2", "#D55E00", "#332288", "#117733",
                   "#88CCEE", "#999933")

parse_task <- function(fname) {
  stem  <- tools::file_path_sans_ext(basename(fname))
  parts <- strsplit(stem, "__", fixed = TRUE)[[1]]
  list(task = parts[1], feature_set = parts[2])
}

# =============================================================================
# 1. metrics_overview.pdf  -- per-task tabular heat-map
# =============================================================================
# One panel per task; rows are (feature_set, weighting), columns are the
# metrics that actually apply to that task type. Cell colour ramps from light
# to deep within each column so the relative ranking is readable at a glance.

task_metrics <- list(
  age       = c("rmse", "mae", "r2", "cor_pearson"),
  age_group = c("accuracy", "macro_auc", "n_classes"),
  gender    = c("auc", "accuracy", "brier")
)
# Direction: "lower" -> smaller is better (RMSE, MAE, Brier);
# default ("higher") -> larger is better
metric_dir <- c(rmse = "lower", mae = "lower", brier = "lower")

make_metric_panel <- function(task_name) {
  m <- metrics |>
    filter(task == task_name) |>
    select(feature_set, weighting, all_of(task_metrics[[task_name]])) |>
    pivot_longer(-c(feature_set, weighting),
                 names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric, levels = task_metrics[[task_name]]),
           row    = paste(feature_set, weighting, sep = " · "))
  # Per-metric rank for the fill colour, oriented so 1.0 = best
  m <- m |>
    group_by(metric) |>
    mutate(rank01 = {
      r <- rank(value, ties.method = "average") / n()
      if (metric_dir[as.character(metric[1])] %in% "lower") 1 - r else r
    }) |> ungroup()

  ggplot(m, aes(x = metric, y = row, fill = rank01)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.3f", value)), size = 3) +
    scale_fill_gradient(low = "#FFFFFF", high = "#2C7FB8", limits = c(0, 1),
                        guide = "none") +
    scale_x_discrete(position = "top") +
    labs(title = task_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.x.top = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}

p_overview <- make_metric_panel("age") +
              make_metric_panel("age_group") +
              make_metric_panel("gender") +
  plot_layout(ncol = 1, heights = c(1, 1, 1)) +
  plot_annotation(title = "Module 10 -- test-set metrics by task and feature set",
                  subtitle = "Cell colour is per-metric rank (deeper = better, accounting for direction)",
                  theme = theme(plot.title = element_text(face = "bold")))
ggsave(file.path(out_fig, "metrics_overview.pdf"), p_overview,
       width = 8, height = 8)
message("Wrote ", file.path(out_fig, "metrics_overview.pdf"))

# =============================================================================
# 2. ROC curves  (carry-over with grafify colour-by-class + linetype-by-fs)
# =============================================================================
roc_long <- list()
for (f in pred_files[grepl("^gender__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  r <- pROC::roc(d$y_true, d$y_pred_prob, quiet = TRUE)
  roc_long[[length(roc_long) + 1L]] <- tibble(
    task_panel = "gender (binary)",
    feature_set = meta$feature_set, class = "Female (positive)",
    fpr = 1 - r$specificities, tpr = r$sensitivities,
    auc = as.numeric(pROC::auc(r)))
}
for (f in pred_files[grepl("^age_group__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  for (cl in sub("^prob_", "", grep("^prob_", names(d), value = TRUE))) {
    yt <- as.integer(as.character(d$y_true) == cl)
    if (length(unique(yt)) < 2) next
    r <- pROC::roc(yt, d[[paste0("prob_", cl)]], quiet = TRUE)
    roc_long[[length(roc_long) + 1L]] <- tibble(
      task_panel = "age_group (one-vs-rest)",
      feature_set = meta$feature_set, class = cl,
      fpr = 1 - r$specificities, tpr = r$sensitivities,
      auc = as.numeric(pROC::auc(r)))
  }
}
roc_tbl <- bind_rows(roc_long)
class_levels <- c("14-19","20-29","30-39","40-49","50-59","60-69",
                  setdiff(unique(roc_tbl$class),
                          c("14-19","20-29","30-39","40-49","50-59","60-69")))
roc_tbl$class <- factor(roc_tbl$class, levels = class_levels)
class_palette <- setNames(grafify_class[seq_along(class_levels)], class_levels)
auc_lab <- roc_tbl |> distinct(task_panel, feature_set, class, auc) |>
  arrange(task_panel, feature_set, class) |>
  group_by(task_panel) |>
  summarise(label = paste0(sprintf("%s [%s]: AUC=%.3f",
                                   class, substr(feature_set, 1, 16), auc),
                           collapse = "\n"), .groups = "drop")

p_roc <- ggplot(roc_tbl, aes(x = fpr, y = tpr,
                             colour = class, linetype = feature_set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_path(linewidth = 0.7) +
  geom_text(data = auc_lab, aes(x = 0.55, y = 0.18, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.3) +
  facet_wrap(~ task_panel, ncol = 2) +
  scale_colour_manual(values = class_palette, drop = FALSE) +
  scale_linetype_manual(values = fs_linetype) +
  coord_equal() +
  labs(title = "ROC curves -- all classification tasks (test set)",
       x = "False positive rate", y = "True positive rate",
       colour = "Class", linetype = "Feature set") +
  guides(colour   = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1,
                                 override.aes = list(colour = "grey20"))) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out_fig, "roc_curves_classification.pdf"), p_roc,
       width = 11, height = 6.5)
message("Wrote ", file.path(out_fig, "roc_curves_classification.pdf"))

# =============================================================================
# 3. Precision-Recall curves  (better than ROC when classes imbalanced)
# =============================================================================
pr_curve <- function(y_true_int, scores) {
  ord <- order(-scores)
  yt  <- y_true_int[ord]
  tp  <- cumsum(yt == 1)
  fp  <- cumsum(yt == 0)
  recall    <- tp / sum(yt == 1)
  precision <- tp / (tp + fp)
  baseline  <- mean(yt == 1)
  list(recall = recall, precision = precision, baseline = baseline,
       auprc = sum(diff(c(0, recall)) * precision))
}

pr_long <- list()
for (f in pred_files[grepl("^gender__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  pr <- pr_curve(d$y_true, d$y_pred_prob)
  pr_long[[length(pr_long) + 1L]] <- tibble(
    task_panel = "gender (binary)",
    feature_set = meta$feature_set, class = "Female (positive)",
    recall = pr$recall, precision = pr$precision,
    baseline = pr$baseline, auprc = pr$auprc)
}
for (f in pred_files[grepl("^age_group__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  for (cl in sub("^prob_", "", grep("^prob_", names(d), value = TRUE))) {
    yt <- as.integer(as.character(d$y_true) == cl)
    if (length(unique(yt)) < 2) next
    pr <- pr_curve(yt, d[[paste0("prob_", cl)]])
    pr_long[[length(pr_long) + 1L]] <- tibble(
      task_panel = "age_group (one-vs-rest)",
      feature_set = meta$feature_set, class = cl,
      recall = pr$recall, precision = pr$precision,
      baseline = pr$baseline, auprc = pr$auprc)
  }
}
pr_tbl <- bind_rows(pr_long)
pr_tbl$class <- factor(pr_tbl$class, levels = class_levels)
baselines <- pr_tbl |> distinct(task_panel, class, baseline)

p_pr <- ggplot(pr_tbl, aes(x = recall, y = precision,
                           colour = class, linetype = feature_set)) +
  geom_hline(data = baselines, aes(yintercept = baseline, colour = class),
             linetype = "dotted", linewidth = 0.3, alpha = 0.6) +
  geom_path(linewidth = 0.7) +
  facet_wrap(~ task_panel, ncol = 2) +
  scale_colour_manual(values = class_palette, drop = FALSE) +
  scale_linetype_manual(values = fs_linetype) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Precision-Recall curves -- all classification tasks (test set)",
       subtitle = "Dotted horizontal = class-prevalence baseline (chance precision)",
       x = "Recall (= TPR)", y = "Precision",
       colour = "Class", linetype = "Feature set") +
  guides(colour   = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1,
                                 override.aes = list(colour = "grey20"))) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out_fig, "precision_recall_curves.pdf"), p_pr,
       width = 11, height = 6.5)
message("Wrote ", file.path(out_fig, "precision_recall_curves.pdf"))

# =============================================================================
# 4. Calibration curves  (decile reliability)
# =============================================================================
decile_cal <- function(p, y, n_bins = 10) {
  if (length(p) == 0) return(tibble(bin = integer(), mean_pred = numeric(),
                                    obs_freq = numeric(), n = integer()))
  bks <- unique(quantile(p, probs = seq(0, 1, length.out = n_bins + 1),
                         na.rm = TRUE))
  if (length(bks) < 3) return(tibble(bin = integer(), mean_pred = numeric(),
                                     obs_freq = numeric(), n = integer()))
  bin <- cut(p, breaks = bks, include.lowest = TRUE, labels = FALSE)
  tibble(p, y, bin) |>
    group_by(bin) |>
    summarise(mean_pred = mean(p, na.rm = TRUE),
              obs_freq  = mean(y, na.rm = TRUE),
              n         = n(), .groups = "drop")
}

cal_long <- list()
for (f in pred_files[grepl("^gender__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  cal <- decile_cal(d$y_pred_prob, d$y_true)
  if (nrow(cal)) cal_long[[length(cal_long) + 1L]] <- cal |>
    mutate(task_panel = "gender (binary)",
           feature_set = meta$feature_set, class = "Female (positive)")
}
for (f in pred_files[grepl("^age_group__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  for (cl in sub("^prob_", "", grep("^prob_", names(d), value = TRUE))) {
    yt <- as.integer(as.character(d$y_true) == cl)
    cal <- decile_cal(d[[paste0("prob_", cl)]], yt)
    if (nrow(cal)) cal_long[[length(cal_long) + 1L]] <- cal |>
      mutate(task_panel = "age_group (one-vs-rest)",
             feature_set = meta$feature_set, class = cl)
  }
}
cal_tbl <- bind_rows(cal_long)
cal_tbl$class <- factor(cal_tbl$class, levels = class_levels)

p_cal <- ggplot(cal_tbl,
                aes(x = mean_pred, y = obs_freq,
                    colour = class, linetype = feature_set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.6) +
  geom_point(aes(size = n), alpha = 0.7, shape = 16) +
  scale_size_continuous(range = c(0.6, 3), guide = "none") +
  facet_wrap(~ task_panel, ncol = 2) +
  scale_colour_manual(values = class_palette, drop = FALSE) +
  scale_linetype_manual(values = fs_linetype) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Calibration (reliability) curves -- decile binning",
       subtitle = "Closer to the diagonal = better calibrated; point size = bin n",
       x = "Mean predicted probability", y = "Observed positive frequency",
       colour = "Class", linetype = "Feature set") +
  guides(colour   = guide_legend(order = 1, nrow = 1),
         linetype = guide_legend(order = 2, nrow = 1,
                                 override.aes = list(colour = "grey20"))) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(out_fig, "calibration_curves.pdf"), p_cal,
       width = 11, height = 6.5)
message("Wrote ", file.path(out_fig, "calibration_curves.pdf"))

# =============================================================================
# 5. Confusion matrices  (carry-over)
# =============================================================================
cm_long <- list()
for (f in pred_files[grepl("^gender__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  pred_class <- ifelse(d$y_pred_prob >= 0.5, "Female", "Male")
  obs_class  <- ifelse(d$y_true == 1, "Female", "Male")
  cm <- tibble(task = "gender", feature_set = meta$feature_set,
               y_true = obs_class, y_pred = pred_class) |>
    count(task, feature_set, y_true, y_pred) |>
    group_by(task, feature_set, y_true) |> mutate(prop = n / sum(n)) |> ungroup()
  cm$y_true <- factor(cm$y_true, levels = c("Male","Female"))
  cm$y_pred <- factor(cm$y_pred, levels = c("Male","Female"))
  cm_long[[length(cm_long) + 1L]] <- cm
}
for (f in pred_files[grepl("^age_group__", basename(pred_files))]) {
  meta <- parse_task(f); d <- readRDS(f)
  cm <- tibble(task = "age_group", feature_set = meta$feature_set,
               y_true = as.character(d$y_true),
               y_pred = as.character(d$y_pred_class)) |>
    count(task, feature_set, y_true, y_pred) |>
    group_by(task, feature_set, y_true) |> mutate(prop = n / sum(n)) |> ungroup()
  age_levels <- c("14-19","20-29","30-39","40-49","50-59","60-69")
  present <- intersect(age_levels, c(cm$y_true, cm$y_pred))
  cm$y_true <- factor(cm$y_true, levels = present)
  cm$y_pred <- factor(cm$y_pred, levels = present)
  cm_long[[length(cm_long) + 1L]] <- cm
}
cm_tbl <- bind_rows(cm_long)
p_cm <- ggplot(cm_tbl, aes(x = y_pred, y = y_true, fill = prop)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", prop)), size = 2.6) +
  scale_fill_gradient(low = "#FFFFFF", high = "#0072B2", limits = c(0, 1)) +
  facet_grid(task ~ feature_set, scales = "free", space = "free") +
  labs(title = "Confusion matrices (row-normalised) -- all classification tasks",
       x = "Predicted class", y = "Observed class", fill = "Row prop.") +
  coord_equal() +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_fig, "confusion_matrices.pdf"), p_cm,
       width = 10, height = 7)
message("Wrote ", file.path(out_fig, "confusion_matrices.pdf"))

# =============================================================================
# 6. Feature importance  -- top-20 taxa per task by |glmnet coefficient|
# =============================================================================
# For multinomial fits, importance = max |coef| across classes (pool by feature).
# For binomial / gaussian fits, importance = |coef|.
get_top_coefs <- function(fit, top_n = 20) {
  co <- coef(fit, s = "lambda.min")
  if (is.list(co)) {
    # multinomial: list of class coefficient matrices
    mat <- do.call(cbind, lapply(co, function(m) as.numeric(m)))
    rownames(mat) <- rownames(co[[1]])
    imp <- apply(abs(mat), 1, max)
  } else {
    imp <- abs(as.numeric(co))
    names(imp) <- rownames(co)
  }
  imp <- imp[names(imp) != "(Intercept)"]
  tibble(feature = names(imp), importance = imp) |>
    filter(importance > 0) |>
    arrange(desc(importance)) |>
    slice_head(n = top_n)
}

fi_long <- list()
for (f in model_files) {
  meta <- parse_task(f)
  fit  <- readRDS(f)
  ti   <- get_top_coefs(fit, top_n = 20)
  if (nrow(ti)) fi_long[[length(fi_long) + 1L]] <- ti |>
    mutate(task = meta$task, feature_set = meta$feature_set)
}
fi_tbl <- bind_rows(fi_long)

if (nrow(fi_tbl) > 0) {
  # Trim long "clr_RSV_genus..." labels to last 24 chars for readability
  fi_tbl$feature_short <- ifelse(nchar(fi_tbl$feature) > 24,
                                 paste0("...", substr(fi_tbl$feature,
                                        nchar(fi_tbl$feature) - 23,
                                        nchar(fi_tbl$feature))),
                                 fi_tbl$feature)
  fi_tbl <- fi_tbl |>
    group_by(task, feature_set) |>
    arrange(importance, .by_group = TRUE) |>
    mutate(rank = row_number()) |>
    ungroup()

  p_fi <- ggplot(fi_tbl,
                 aes(x = importance, y = reorder(feature_short, importance),
                     fill = task)) +
    geom_col() +
    facet_grid(task ~ feature_set, scales = "free", space = "free_y") +
    scale_fill_manual(values = c(age = "#0072B2", age_group = "#009E73",
                                 gender = "#D55E00"), guide = "none") +
    labs(title = "Top-20 features per task and feature set",
         subtitle = "Importance = absolute glmnet coefficient at lambda.min (max over classes for multinomial)",
         x = "|coefficient|", y = NULL) +
    theme_bw(base_size = 8) +
    theme(axis.text.y = element_text(size = 6))
  ggsave(file.path(out_fig, "feature_importance.pdf"), p_fi,
         width = 12, height = 10)
  message("Wrote ", file.path(out_fig, "feature_importance.pdf"))
}

# =============================================================================
# 7. Regression diagnostics  -- pred-vs-true + residuals + Q-Q
# =============================================================================
reg_long <- list()
for (f in pred_files) {
  d <- readRDS(f)
  if (!"y_pred" %in% names(d)) next  # only regression preds carry a y_pred col
  meta <- parse_task(f)
  reg_long[[length(reg_long) + 1L]] <- d |>
    select(SEQN, y_true, y_pred, weight) |>
    mutate(task = meta$task, feature_set = meta$feature_set,
           residual = y_true - y_pred)
}
reg_tbl <- bind_rows(reg_long)

if (nrow(reg_tbl) > 0) {
  p_pvt <- ggplot(reg_tbl, aes(x = y_true, y = y_pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.25, size = 0.55, colour = "#0072B2") +
    geom_smooth(method = "lm", se = FALSE, colour = "#D55E00", linewidth = 0.6) +
    facet_grid(task ~ feature_set, scales = "free") +
    labs(subtitle = "(a) Observed vs predicted (dashed = 1:1, orange = LM fit)",
         x = "Observed", y = "Predicted") +
    theme_bw(base_size = 9) + theme(aspect.ratio = 1)

  p_res <- ggplot(reg_tbl, aes(x = y_pred, y = residual)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.25, size = 0.55, colour = "#009E73") +
    geom_smooth(method = "loess", se = FALSE, colour = "#D55E00", linewidth = 0.6) +
    facet_grid(task ~ feature_set, scales = "free") +
    labs(subtitle = "(b) Residuals vs predicted (loess = bias trend)",
         x = "Predicted", y = "Residual (observed - predicted)") +
    theme_bw(base_size = 9) + theme(aspect.ratio = 1)

  qq_dat <- reg_tbl |> group_by(task, feature_set) |>
    mutate(theoretical = qnorm(ppoints(n()))[rank(residual)]) |> ungroup()
  p_qq <- ggplot(qq_dat, aes(x = theoretical, y = residual)) +
    geom_abline(slope = sd(qq_dat$residual, na.rm = TRUE),
                intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.25, size = 0.55, colour = "#CC79A7") +
    facet_grid(task ~ feature_set, scales = "free") +
    labs(subtitle = "(c) Residual Q-Q plot vs normal",
         x = "Theoretical normal quantile", y = "Residual") +
    theme_bw(base_size = 9) + theme(aspect.ratio = 1)

  p_reg <- (p_pvt + p_res + p_qq) +
    plot_layout(ncol = 3) +
    plot_annotation(title = "Regression diagnostics (test set)",
                    theme = theme(plot.title = element_text(face = "bold")))
  ggsave(file.path(out_fig, "regression_diagnostics.pdf"), p_reg,
         width = 14, height = 6)
  message("Wrote ", file.path(out_fig, "regression_diagnostics.pdf"))
}

message("Done. Figures written to: ", out_fig)
