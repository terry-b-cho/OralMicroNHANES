#!/usr/bin/env Rscript

# =============================================================================
# Module 10 - prediction pipeline, step 2: train + evaluate
# =============================================================================
# Three representative prediction tasks on the oral-microbiome cohort:
#   (1) age            - regression                  (age ~ CLR features [+ covariates])
#   (2) age_group      - multi-class classification  (6 NHANES age strata)
#   (3) gender         - binary classification       (Male vs Female)
#
# For each task we compare two feature sets:
#   A. CLR microbiome only       (clr_*)
#   B. CLR microbiome + basic demographics (ethnicity, education, PIR) when
#      the demographic is NOT the target itself
#
# Model: regularized regression via glmnet (cv.glmnet, alpha = 0.5 elastic net).
# A simple, defensible, interpretable baseline that handles 200+ correlated
# microbiome features well.
#
# Split: 80/20 stratified train/test on a random hold-out (seed 42).
# Survey design (WTMEC2YR + SDMVSTRA + SDMVPSU) is loaded so that
# weighted metrics can be computed on the test set without re-running.
#
# Outputs (under results/analyses_results/10_prediction_analyses_out/):
#   - predictions/<task>__<feature_set>.rds      per-sample y_true + y_pred + weights
#   - metrics_summary.csv                         one row per task x feature_set
#   - models/<task>__<feature_set>.rds            fitted cv.glmnet object
#
# Environment: R >= 4.5 with glmnet, dplyr, tibble, tidyr, readr, pROC.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(pROC)
})

set.seed(42)

base_path  <- PROJECT_ROOT
out_root   <- file.path(base_path, "results/analyses_results/10_prediction_analyses_out")
out_inter  <- file.path(out_root, "intermediate")
out_pred   <- file.path(out_root, "predictions")
out_models <- file.path(out_root, "models")
dir.create(out_pred,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_models, recursive = TRUE, showWarnings = FALSE)

# ---- Load intermediates from step 1 ----------------------------------------
clr_features            <- readRDS(file.path(out_inter, "clr_features.rds"))
covariates_and_targets  <- readRDS(file.path(out_inter, "covariates_and_targets.rds"))
sample_universe         <- readRDS(file.path(out_inter, "sample_universe.rds"))

dataset <- clr_features |>
  inner_join(covariates_and_targets, by = "SEQN") |>
  filter(SEQN %in% sample_universe)
message(sprintf("== Working dataset: %d samples", nrow(dataset)))

clr_cols <- grep("^clr_", names(dataset), value = TRUE)
basic_covariate_cols <- c("ethnicity", "education_lt9", "education_hs",
                          "education_aa", "education_cg", "pir")

# ---- Helpers ----------------------------------------------------------------

build_x_matrix <- function(df, feature_cols) {
  d <- df[, feature_cols, drop = FALSE]
  # One-hot encode any factor / character columns; impute NA in PIR with median
  if ("pir" %in% names(d)) d$pir[is.na(d$pir)] <- median(d$pir, na.rm = TRUE)
  mm <- model.matrix(~ . - 1, data = d)
  mm
}

fit_glmnet_cv <- function(x, y, family) {
  cv.glmnet(x = x, y = y, family = family, alpha = 0.5, nfolds = 5,
            standardize = TRUE, type.measure = if (family == "gaussian") "mse" else "deviance")
}

eval_regression <- function(y_true, y_pred, w = NULL) {
  err <- y_true - y_pred
  if (is.null(w)) w <- rep(1, length(y_true))
  ss_res <- sum(w * err^2)
  ss_tot <- sum(w * (y_true - weighted.mean(y_true, w))^2)
  tibble(
    n       = length(y_true),
    rmse    = sqrt(weighted.mean(err^2, w)),
    mae     = weighted.mean(abs(err), w),
    r2      = 1 - ss_res / ss_tot,
    cor_pearson = stats::cor(y_true, y_pred)
  )
}

eval_binary <- function(y_true, prob, w = NULL) {
  if (is.null(w)) w <- rep(1, length(y_true))
  auc <- as.numeric(pROC::auc(pROC::roc(y_true, prob, quiet = TRUE)))
  thr <- 0.5
  pred_class <- as.integer(prob >= thr)
  tibble(
    n        = length(y_true),
    auc      = auc,
    accuracy = weighted.mean(pred_class == y_true, w),
    brier    = weighted.mean((prob - y_true)^2, w)
  )
}

eval_multiclass <- function(y_true, prob_mat, w = NULL) {
  if (is.null(w)) w <- rep(1, length(y_true))
  pred_class <- colnames(prob_mat)[apply(prob_mat, 1, which.max)]
  acc <- weighted.mean(pred_class == as.character(y_true), w)
  # Macro-AUC (one-vs-rest)
  classes <- colnames(prob_mat)
  aucs <- sapply(classes, function(cl) {
    yt <- as.integer(as.character(y_true) == cl)
    if (length(unique(yt)) < 2) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(yt, prob_mat[, cl], quiet = TRUE)))
  })
  tibble(
    n              = length(y_true),
    accuracy       = acc,
    macro_auc      = mean(aucs, na.rm = TRUE),
    n_classes      = length(classes)
  )
}

stratified_split <- function(y, p_train = 0.8) {
  idx <- seq_along(y)
  # safe per-stratum sampler: empty / size-1 strata go entirely into train
  pick <- function(i) {
    n <- length(i)
    if (n == 0L) return(integer(0))
    if (n == 1L) return(i)
    k <- max(1L, floor(n * p_train))
    if (k >= n) return(i)
    sample(i, k)
  }
  if (is.factor(y) || is.character(y)) {
    train_idx <- unlist(lapply(split(idx, y), pick))
  } else {
    bins <- cut(y, breaks = quantile(y, probs = seq(0, 1, 0.1), na.rm = TRUE),
                include.lowest = TRUE, labels = FALSE)
    train_idx <- unlist(lapply(split(idx, bins), pick))
  }
  list(train = sort(unique(train_idx)), test = setdiff(idx, train_idx))
}

# ---- Define tasks -----------------------------------------------------------

tasks <- list(
  list(name = "age",       target = "age",       family = "gaussian",    type = "regression"),
  list(name = "age_group", target = "age_group", family = "multinomial", type = "multiclass"),
  list(name = "gender",    target = "gender",    family = "binomial",    type = "binary")
)

# Two feature sets per task; covariates dropped if they trivially leak the target
feature_sets_for <- function(target) {
  fs <- list(microbiome_only = clr_cols)
  cov_cols <- setdiff(basic_covariate_cols, target)
  fs$microbiome_plus_basic_demo <- c(clr_cols, cov_cols)
  fs
}

# ---- Run --------------------------------------------------------------------

all_metrics <- list()

for (task in tasks) {
  fs_list <- feature_sets_for(task$target)
  d <- dataset |> filter(!is.na(.data[[task$target]]))
  for (fs_name in names(fs_list)) {
    feat_cols <- fs_list[[fs_name]]
    d_full <- d |> select(all_of(c("SEQN", "wtmec2yr", "sdmvstra", "sdmvpsu",
                                   task$target, feat_cols)))
    d_full <- d_full[complete.cases(d_full[, c(task$target, feat_cols)]), ]
    if (nrow(d_full) < 100) {
      message(sprintf("  [skip] task=%s fs=%s — too few complete cases (%d)",
                      task$name, fs_name, nrow(d_full)))
      next
    }
    # For classification tasks, drop classes with < 10 observations to avoid
    # cv.glmnet 'one class has 1 or 0 observations' on the train-fold split
    if (task$family %in% c("binomial", "multinomial")) {
      yv <- d_full[[task$target]]
      tab <- table(yv, useNA = "no")
      rare <- names(tab)[tab < 10]
      if (length(rare) > 0) {
        message(sprintf("  [info] task=%s fs=%s — dropping rare classes (n<10): %s",
                        task$name, fs_name, paste(rare, collapse = ", ")))
        d_full <- d_full[!as.character(yv) %in% rare, , drop = FALSE]
        if (is.factor(d_full[[task$target]])) d_full[[task$target]] <- droplevels(d_full[[task$target]])
      }
    }
    y <- d_full[[task$target]]
    split <- stratified_split(y)
    x_train <- build_x_matrix(d_full[split$train, , drop = FALSE], feat_cols)
    x_test  <- build_x_matrix(d_full[split$test,  , drop = FALSE], feat_cols)
    # Align columns of x_test to x_train (drop unseen levels)
    common_cols <- intersect(colnames(x_train), colnames(x_test))
    x_train <- x_train[, common_cols, drop = FALSE]
    x_test  <- x_test[,  common_cols, drop = FALSE]
    y_train <- y[split$train]
    y_test  <- y[split$test]
    w_test  <- d_full$wtmec2yr[split$test]

    message(sprintf("== task=%-10s fs=%-26s n_train=%d n_test=%d p=%d",
                    task$name, fs_name, length(y_train), length(y_test), ncol(x_train)))

    if (task$family == "binomial") {
      y_train_int <- as.integer(y_train == levels(y_train)[2])
      y_test_int  <- as.integer(y_test  == levels(y_test)[2])
      fit  <- fit_glmnet_cv(x_train, y_train_int, family = "binomial")
      prob <- as.numeric(predict(fit, newx = x_test, s = "lambda.min", type = "response"))
      m_unw <- eval_binary(y_test_int, prob)
      m_w   <- eval_binary(y_test_int, prob, w = w_test)
      m_unw$weighting <- "unweighted"
      m_w$weighting   <- "survey_weighted"
      m <- bind_rows(m_unw, m_w)
      preds <- tibble(SEQN = d_full$SEQN[split$test], y_true = y_test_int,
                      y_pred_prob = prob, weight = w_test)
    } else if (task$family == "multinomial") {
      fit  <- fit_glmnet_cv(x_train, y_train, family = "multinomial")
      prob <- predict(fit, newx = x_test, s = "lambda.min", type = "response")[, , 1]
      m_unw <- eval_multiclass(y_test, prob)
      m_w   <- eval_multiclass(y_test, prob, w = w_test)
      m_unw$weighting <- "unweighted"
      m_w$weighting   <- "survey_weighted"
      m <- bind_rows(m_unw, m_w)
      preds <- tibble(SEQN = d_full$SEQN[split$test], y_true = as.character(y_test),
                      y_pred_class = colnames(prob)[apply(prob, 1, which.max)],
                      weight = w_test)
    } else {
      fit  <- fit_glmnet_cv(x_train, y_train, family = "gaussian")
      pred <- as.numeric(predict(fit, newx = x_test, s = "lambda.min"))
      m_unw <- eval_regression(y_test, pred)
      m_w   <- eval_regression(y_test, pred, w = w_test)
      m_unw$weighting <- "unweighted"
      m_w$weighting   <- "survey_weighted"
      m <- bind_rows(m_unw, m_w)
      preds <- tibble(SEQN = d_full$SEQN[split$test], y_true = y_test,
                      y_pred = pred, weight = w_test)
    }

    m$task <- task$name; m$feature_set <- fs_name; m$model <- "glmnet_enet"
    all_metrics[[length(all_metrics) + 1]] <- m
    saveRDS(preds, file.path(out_pred,   paste0(task$name, "__", fs_name, ".rds")))
    saveRDS(fit,   file.path(out_models, paste0(task$name, "__", fs_name, ".rds")))
  }
}

metrics_tbl <- bind_rows(all_metrics) |>
  relocate(task, feature_set, model, weighting, n)
write_csv(metrics_tbl, file.path(out_root, "metrics_summary.csv"))

message("\n== Metrics summary ==")
print(metrics_tbl, n = Inf)
message("\nDone. Outputs under: ", out_root)
