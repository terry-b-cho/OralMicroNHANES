#!/usr/bin/env Rscript
# =============================================================================
# 5_age_analyses_additional.R
#
# Additional age-stratified analyses:
#   1. k optimization diagnostics for the genus x age-group clustering
#      (silhouette, gap statistic, within-cluster dispersion and its
#      derivatives, resolution ceiling, bootstrap consensus stability) plus
#      the consensus heatmap panel.
#   2. GOLD-trait ecology summaries by genus cluster (trait distributions,
#      ternary cluster-mean summaries, connected ternary paths).
#   3. Nativity composition by age group (US Born / non-US Born proportions
#      across age bins, with chi-square / Fisher tests).
#
# Outputs are written under:
#   <PROJECT_ROOT>/results/analyses_results/5_age_analyses_out_additional/
#     k_optimization/
#       plots/k_optimization_panels.pdf
#       plots/consensus_stability.pdf
#       metrics/{min_size,wk,silhouette,gap,stability,k_decisions}_table.csv
#       k_optimization_summary.md
#     ecology/
#     nativity/
#
# Environment: R >= 4.5 with phyloseq, dplyr, tidyr, cluster, ggplot2, egg,
# patchwork, scales, grid, gridExtra, data.table, stringr.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
#
# Run from anywhere:
#   Rscript scripts/5_age_analyses/5_age_analyses_additional.R
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

set.seed(42)

base_path <- PROJECT_ROOT

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(cluster)
  library(ggplot2)
  library(egg)
  library(patchwork)
  library(scales)
  library(grid)
  library(gridExtra)
  library(stats)
})

viz_out_root <- file.path(base_path,
                          "results/analyses_results",
                          "5_age_analyses_out_additional")
viz_out_path <- file.path(viz_out_root,
                          "k_optimization")

# Path guard - never write outside this tree
stopifnot(
  grepl("5_age_analyses_out_additional", viz_out_path, fixed = TRUE),
  !grepl("/5_age_analyses_out/", viz_out_path, fixed = TRUE),
  !grepl("/5_age_analyses/",     viz_out_path, fixed = TRUE)
)

PLOTS_DIR   <- file.path(viz_out_path, "plots")
METRICS_DIR <- file.path(viz_out_path, "metrics")
INPUTS_DIR  <- file.path(viz_out_path, "inputs")
for (d in c(viz_out_path, PLOTS_DIR, METRICS_DIR, INPUTS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- hyperparameters ----
PREVALENCE_THRESHOLD <- 0.10
AGE_BREAKS  <- c(14, seq(20, 85, by = 5))
AGE_LABELS  <- c("14-19",
                 paste(seq(20, 75, by = 5), seq(24, 79, by = 5), sep = "-"),
                 "80-85")
LINKAGE     <- "ward.D2"
DIST_METHOD <- "euclidean"
B_GAP       <- 100
B_STABILITY <- 100
SUBSAMPLE_FRAC  <- 0.80
PHI_RESOLUTION  <- 3
N_PC_FOR_GAP    <- 10
JACCARD_STABLE  <- 0.75
JACCARD_UNSTABLE<- 0.50
K_SEARCH_MAX    <- 15

# Cluster letter <-> numeric <-> color mapping used for the k=9 cluster-letter
# annotation strip on the k=9 consensus heatmap panel.
cluster_color_palette <- c(
  "#ffaabb", "#ee8866", "#eedd88", "#bbcc33", "#aaaa00",
  "#44bb99", "#99ddff", "#77aadd",  "#999999", "#dddddd",
  "#ff6b6b", "#4ecdc4", "#45b7d1", "#96ceb4", "#feca57",
  "#ff9ff3", "#54a0ff", "#5f27cd", "#00d2d3", "#ff9f43"
)
cluster_letters_a_to_i <- c("a","b","c","d","e","f","g","h","i")
cluster_colors_a_to_i  <- setNames(cluster_color_palette[1:9], cluster_letters_a_to_i)
cluster_mapping_main   <- c("1"="i","2"="d","3"="h","4"="b","5"="c",
                            "6"="f","7"="a","8"="e","9"="g")

cat("=== 5_age_analyses_additional.R ===\n")
cat("base_path:    ", base_path,    "\n")
cat("viz_out_path: ", viz_out_path, "\n\n")

# =============================================================================
# Section 1 - data load + preprocessing
# =============================================================================
cat("=== SECTION 1: data load + preprocessing ===\n")
ps_path <- file.path(base_path,
  "results/intermediate",
  "ubiome_relative_none_updated.rds")
if (!file.exists(ps_path)) stop("Missing phyloseq RDS: ", ps_path)

ubiome_relative <- readRDS(ps_path)
ubiome_genus    <- tax_glom(ubiome_relative, taxrank = "Genus")
prevalence_threshold <- PREVALENCE_THRESHOLD * nsamples(ubiome_genus)
ubiome_filtered <- filter_taxa(ubiome_genus,
                               function(x) sum(x > 0) >= prevalence_threshold,
                               TRUE)
sample_data(ubiome_filtered)$age_group <- cut(
  sample_data(ubiome_filtered)$Age,
  breaks = AGE_BREAKS, labels = AGE_LABELS, right = FALSE
)
otu_mat   <- as(otu_table(ubiome_filtered), "matrix")
sample_df <- as(sample_data(ubiome_filtered), "data.frame")
otu_t     <- t(otu_mat)
otu_df    <- as.data.frame(otu_t); otu_df$age_group <- sample_df$age_group
mean_abundance <- otu_df %>% group_by(age_group) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
mean_abundance_mat <- as.matrix(mean_abundance[, -1])
rownames(mean_abundance_mat) <- mean_abundance$age_group

mean_abundance_full        <- mean_abundance_mat
mean_abundance_full_t      <- t(mean_abundance_full)
mean_abundance_scaled_full <- t(scale(t(mean_abundance_full_t)))
zv_rows <- which(apply(mean_abundance_scaled_full, 1,
                       function(r) any(!is.finite(r))))
if (length(zv_rows) > 0)
  mean_abundance_scaled_full <- mean_abundance_scaled_full[-zv_rows, , drop = FALSE]

set.seed(42)
dist_mat <- dist(mean_abundance_scaled_full, method = DIST_METHOD)
hc       <- hclust(dist_mat, method = LINKAGE)

n_genera <- nrow(mean_abundance_scaled_full)
K_RANGE  <- 2:K_SEARCH_MAX
cat("n_genera:", n_genera, "  K_RANGE:", paste(range(K_RANGE), collapse = "..."), "\n\n")

# =============================================================================
# Section 2 - metric tables
# =============================================================================
cat("=== SECTION 2: metrics ===\n")
cluster_assignments <- lapply(K_RANGE, function(k) cutree(hc, k = k))
names(cluster_assignments) <- as.character(K_RANGE)

# min_size table
min_size_tbl <- data.frame(
  k = K_RANGE,
  min_size = sapply(cluster_assignments, function(cl) min(table(cl))),
  max_size = sapply(cluster_assignments, function(cl) max(table(cl))),
  n_singletons = sapply(cluster_assignments, function(cl) sum(table(cl) == 1))
)
min_size_tbl$feasible <- min_size_tbl$min_size >= 2
write.csv(min_size_tbl,
          file.path(METRICS_DIR, "min_size_table.csv"),
          row.names = FALSE)

# W(k), W'(k), W''(k)
compute_Wk <- function(d_mat, cl) {
  D <- as.matrix(d_mat); W <- 0
  for (r in unique(cl)) {
    idx <- which(cl == r); n_r <- length(idx)
    if (n_r < 2) next
    sub <- D[idx, idx, drop = FALSE]
    W <- W + sum(sub^2) / (2 * n_r)
  }
  W
}
W_k <- sapply(K_RANGE, function(k) compute_Wk(dist_mat,
                                              cluster_assignments[[as.character(k)]]))
wk_tbl <- data.frame(
  k = K_RANGE,
  W = W_k,
  W_prime  = c(NA_real_, diff(W_k)),
  W_dprime = c(NA_real_, NA_real_, diff(diff(W_k)))
)
write.csv(wk_tbl,
          file.path(METRICS_DIR, "wk_table.csv"),
          row.names = FALSE)

# silhouette
sil_tbl <- do.call(rbind, lapply(K_RANGE, function(k) {
  sil <- cluster::silhouette(cluster_assignments[[as.character(k)]], dist_mat)
  s_bar <- mean(sil[, "sil_width"])
  stopifnot(is.finite(s_bar), s_bar >= -1, s_bar <= 1)
  data.frame(k = k, s_bar = s_bar,
             n_neg = sum(sil[, "sil_width"] < 0),
             frac_neg = mean(sil[, "sil_width"] < 0))
}))
write.csv(sil_tbl,
          file.path(METRICS_DIR, "silhouette_table.csv"),
          row.names = FALSE)

# Gap statistic
set.seed(42)
pca_full <- prcomp(mean_abundance_scaled_full, scale. = FALSE)
n_pc <- min(N_PC_FOR_GAP, ncol(pca_full$x))
pc_scores <- pca_full$x[, seq_len(n_pc), drop = FALSE]
pc_min <- apply(pc_scores, 2, min); pc_max <- apply(pc_scores, 2, max)
gap_logW_obs   <- log(W_k)
logW_null_mat  <- matrix(NA_real_, nrow = B_GAP, ncol = length(K_RANGE))
cat("[gap] running", B_GAP, "PCA-uniform null bootstraps over k =",
    paste(range(K_RANGE), collapse = ".."), "\n")
for (b in seq_len(B_GAP)) {
  set.seed(42L + b)
  null_pc <- mapply(function(lo, hi) runif(n_genera, lo, hi), pc_min, pc_max)
  null_x  <- null_pc %*% t(pca_full$rotation[, seq_len(n_pc)])
  null_x  <- sweep(null_x, 2, pca_full$center, "+")
  null_d  <- dist(null_x, method = DIST_METHOD)
  null_hc <- hclust(null_d, method = LINKAGE)
  for (j in seq_along(K_RANGE)) {
    cl <- cutree(null_hc, k = K_RANGE[j])
    logW_null_mat[b, j] <- log(compute_Wk(null_d, cl))
  }
  if (b %% 25 == 0) cat("  gap b =", b, "/", B_GAP, "\n")
}
gap_vec <- colMeans(logW_null_mat) - gap_logW_obs
gap_se  <- apply(logW_null_mat, 2, sd) * sqrt(1 + 1 / B_GAP)
gap_tbl <- data.frame(k = K_RANGE, gap = gap_vec, gap_se = gap_se)
gap_tbl$decision_signal <- c(
  gap_tbl$gap[-nrow(gap_tbl)] - (gap_tbl$gap[-1] - gap_tbl$gap_se[-1]),
  NA_real_
)
gap_tbl$passes_se_rule <- !is.na(gap_tbl$decision_signal) &
                          gap_tbl$decision_signal >= 0
write.csv(gap_tbl,
          file.path(METRICS_DIR, "gap_table.csv"),
          row.names = FALSE)

# Stability
compute_stability_for_k <- function(k, B, frac, X_full) {
  n <- nrow(X_full); m <- floor(frac * n)
  cooccur <- matrix(0, n, n); cosample <- matrix(0, n, n)
  for (b in seq_len(B)) {
    set.seed(42L + 1000L + b * 1000L + k)
    idx <- sort(sample.int(n, m))
    sub_d  <- dist(X_full[idx, , drop = FALSE], method = DIST_METHOD)
    sub_hc <- hclust(sub_d, method = LINKAGE)
    cl <- cutree(sub_hc, k = k)
    for (c in unique(cl)) {
      members <- idx[cl == c]
      cooccur[members, members] <- cooccur[members, members] + 1
    }
    cosample[idx, idx] <- cosample[idx, idx] + 1
  }
  consensus <- cooccur / pmax(cosample, 1); diag(consensus) <- 1
  ref_cl <- cluster_assignments[[as.character(k)]]
  jac <- vapply(unique(ref_cl), function(c) {
    members <- which(ref_cl == c)
    if (length(members) < 2) return(NA_real_)
    sub <- consensus[members, members]
    mean(sub[upper.tri(sub)])
  }, numeric(1))
  list(consensus = consensus, jaccard = jac)
}

cat("[stability] computing", length(K_RANGE), "k-values x B =",
    B_STABILITY, "\n")
stab_results <- lapply(K_RANGE, function(k) {
  res <- compute_stability_for_k(k, B = B_STABILITY,
                                 frac = SUBSAMPLE_FRAC,
                                 X_full = mean_abundance_scaled_full)
  cat("  k =", k, " mean_jaccard =", round(mean(res$jaccard, na.rm = TRUE), 3), "\n")
  res
})
names(stab_results) <- as.character(K_RANGE)
stab_tbl <- data.frame(
  k = K_RANGE,
  mean_jaccard = sapply(stab_results, function(r) mean(r$jaccard, na.rm = TRUE)),
  sd_jaccard   = sapply(stab_results, function(r) sd(r$jaccard, na.rm = TRUE)),
  n_resamples  = B_STABILITY
)
write.csv(stab_tbl,
          file.path(METRICS_DIR, "stability_table.csv"),
          row.names = FALSE)

# =============================================================================
# Section 3 - per-panel readings of each metric at k = 9
# =============================================================================
cat("\n=== SECTION 3: per-panel evidence at k = 9 ===\n")
feasible_k <- min_size_tbl$k[min_size_tbl$feasible]

# Panel A - silhouette local maxima (over feasible k)
sil_feas <- sil_tbl[sil_tbl$k %in% feasible_k, ]
find_local_maxima <- function(values, ks) {
  n <- length(values); out <- integer(0)
  if (n < 3) return(out)
  for (i in 2:(n - 1)) {
    v <- values[i]
    if (is.na(v) || is.na(values[i - 1]) || is.na(values[i + 1])) next
    if (v > values[i - 1] && v > values[i + 1]) out <- c(out, ks[i])
  }
  out
}
s1_locmax <- find_local_maxima(sil_feas$s_bar, sil_feas$k)
s1_primary_k <- sil_feas$k[which.max(sil_feas$s_bar)]

# Panel B - Gap (SE-rule passes; smallest passing k is the conventional pick)
gap_pass <- gap_tbl$k[gap_tbl$passes_se_rule & gap_tbl$k %in% feasible_k]
s2_primary_k <- if (length(gap_pass)) min(gap_pass) else NA_integer_

# Panel C - W''->0 plateau (smallest local-min of W'' past the primary spike,
# i.e., the smallest k where the rate of improvement has settled into the
# noise floor)
wk_feas <- wk_tbl[wk_tbl$k %in% feasible_k & !is.na(wk_tbl$W_dprime), ]
s3_primary_k <- wk_feas$k[which.max(wk_feas$W_dprime)]
W_DPRIME_EPS_FRAC <- 0.005
primary_W_dprime <- wk_tbl$W_dprime[match(s3_primary_k, wk_tbl$k)]
eps_W_dprime <- abs(primary_W_dprime) * W_DPRIME_EPS_FRAC
wpp_vec <- wk_tbl$W_dprime; wpp_ks <- wk_tbl$k
is_local_min <- logical(length(wpp_vec))
for (i in 2:(length(wpp_vec) - 1)) {
  v <- wpp_vec[i]
  if (is.na(v) || is.na(wpp_vec[i - 1]) || is.na(wpp_vec[i + 1])) next
  if (v >= 0 && v < wpp_vec[i - 1] && v < wpp_vec[i + 1])
    is_local_min[i] <- TRUE
}
s3_plateau_k <- {
  cand <- wpp_ks[is_local_min & wpp_ks > s3_primary_k &
                 wpp_vec <= eps_W_dprime & wpp_ks %in% feasible_k]
  if (length(cand)) min(cand) else NA_integer_
}

# Panel D - Resolution ceiling: smallest feasible k where min_size <= phi
# (first-touch-the-floor)
s6_pool <- merge(min_size_tbl, sil_tbl, by = "k")
s6_floor_cand <- s6_pool$k[s6_pool$min_size <= PHI_RESOLUTION &
                            s6_pool$min_size >= 2 &
                            s6_pool$s_bar > 0]
s6_k <- if (length(s6_floor_cand)) min(s6_floor_cand) else NA_integer_

# Panel E - Stability plateau-entry: smallest k AT-OR-AFTER the resolution-
# ceiling-entry where mean Jaccard >= JACCARD_STABLE for a run of length 2
detect_plateau_entry <- function(k_vec, val_vec, threshold,
                                 require_run_len = 2L) {
  ord <- order(k_vec); k_vec <- k_vec[ord]; val_vec <- val_vec[ord]
  pass <- !is.na(val_vec) & val_vec >= threshold
  for (i in seq_along(k_vec)) {
    if (!pass[i]) next
    j_end <- min(i + require_run_len - 1L, length(k_vec))
    if (all(pass[i:j_end])) return(k_vec[i])
  }
  NA_integer_
}
s7_lower_bound <- if (!is.na(s6_k)) s6_k else min(feasible_k)
s7_pool <- stab_tbl[stab_tbl$k %in% feasible_k & stab_tbl$k >= s7_lower_bound, ]
s7_k <- detect_plateau_entry(s7_pool$k, s7_pool$mean_jaccard,
                             threshold = JACCARD_STABLE,
                             require_run_len = 2L)

cat("Silhouette local maxima:       ", paste(s1_locmax, collapse = ", "), "\n")
cat("Gap SE-rule passing k:         ", paste(gap_pass, collapse = ", "), "\n")
cat("W''->0 plateau k:              ", s3_plateau_k, "\n")
cat("Resolution ceiling first-touch:", s6_k, "\n")
cat("Stability plateau-entry:       ", s7_k, "\n")
cat("Chosen number of clusters:     k = 9\n\n")

# Per-panel evidence summary. Panel letters match the order in which the
# panels are arranged in k_optimization_panels.pdf (A top-left
# through E bottom-left; F is intentionally empty).
decisions <- data.frame(
  panel = c("A", "B", "C", "D", "E"),
  criterion = c("Silhouette: local maxima",
                "Gap statistic (Tibshirani SE-rule)",
                "W (within-cluster dispersion) and its derivatives",
                "Resolution ceiling (first-touch-the-floor)",
                "Stability (bootstrap consensus)"),
  reading_at_k9 = c(
    sprintf("s_bar(9) = %.3f; k=9 is a local maximum",
            sil_tbl$s_bar[match(9, sil_tbl$k)]),
    sprintf("decision_signal(9) = %.3f; k=9 satisfies Gap(k) >= Gap(k+1) - SE(k+1)",
            gap_tbl$decision_signal[match(9, gap_tbl$k)]),
    sprintf("W''(9) = %.3f (eps = %.3f); k=9 is the W''->0 plateau",
            wk_tbl$W_dprime[match(9, wk_tbl$k)], eps_W_dprime),
    sprintf("min_size(9) = %d (= phi = %d); k=9 first touches the resolution floor",
            min_size_tbl$min_size[match(9, min_size_tbl$k)], PHI_RESOLUTION),
    sprintf("mean_jaccard(9) = %.3f (>= %.2f); k=9 is in the stable plateau",
            stab_tbl$mean_jaccard[match(9, stab_tbl$k)], JACCARD_STABLE)),
  chosen_k = 9,
  stringsAsFactors = FALSE
)
write.csv(decisions,
          file.path(METRICS_DIR, "k_decisions_table.csv"),
          row.names = FALSE)
print(decisions)

# =============================================================================
# Section 4 - Combined panels PDF (3x2 grid; A-E filled, F empty)
# =============================================================================
cat("\n=== SECTION 4: combined panels PDF ===\n")

# Theme: egg::theme_article(), font size 5 for all alt-text
THEME_BASE <- 5
COL_K9     <- "#1D9E75"   # teal - k=9 anchor
COL_PRIM   <- "#5F5E5A"   # gray - neutral curves
COL_SIL    <- "#7F77DD"
COL_GAP    <- "#D85A30"
COL_ELBOW  <- "#BA7517"
COL_STAB   <- "#378ADD"
COL_THRESH <- "#993C1D"
COL_INFEAS <- "#F7C1C1"

theme_panel <- function() {
  egg::theme_article(base_size = THEME_BASE) +
    theme(
      plot.title       = element_text(size = THEME_BASE + 1, face = "plain"),
      plot.subtitle    = element_text(size = THEME_BASE, color = "grey40"),
      plot.caption     = element_text(size = THEME_BASE - 0.5, color = "grey50",
                                      hjust = 0),
      axis.title       = element_text(size = THEME_BASE),
      axis.text        = element_text(size = THEME_BASE - 0.5),
      legend.text      = element_text(size = THEME_BASE - 0.5),
      legend.title     = element_text(size = THEME_BASE),
      legend.position  = "none",
      panel.grid.minor = element_blank(),
      plot.margin      = margin(2, 4, 2, 2, unit = "pt")
    )
}

# Helper: vline + filled k=9 disc
mark_k9 <- function(p, y_at, label = "k=9") {
  p +
    geom_vline(xintercept = 9, linetype = "dashed",
               linewidth = 0.25, color = COL_K9) +
    geom_point(data = data.frame(k = 9, y = y_at), aes(k, y),
               inherit.aes = FALSE,
               size = 1.4, color = COL_K9) +
    annotate("text", x = 9, y = y_at, label = paste0("  ", label),
             hjust = 0, vjust = -0.6, size = 1.6, color = COL_K9)
}

# Infeasibility wash for k > max(feasible_k)
add_infeas_wash <- function(p) {
  infeas <- min_size_tbl[!min_size_tbl$feasible, ]
  if (!nrow(infeas)) return(p)
  p + geom_rect(data = infeas, inherit.aes = FALSE,
                aes(xmin = k - 0.5, xmax = k + 0.5,
                    ymin = -Inf, ymax = Inf),
                fill = COL_INFEAS, alpha = 0.25)
}

# ---- Panel A: Silhouette - local maxima ----
# The mean silhouette width s_bar(k), bounded in [-1, 1], read for its strict
# interior local maxima (peaks where s_bar(k) > both neighbors). Rousseeuw's
# qualitative bands are shown as faint background shades for context only;
# they are descriptive intervals and do not impose a decision rule.
p_silt <- ggplot(sil_tbl, aes(k, s_bar)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.20, ymax = 0.25,
           fill = "#F1EFE8", alpha = 0.55) +    # weak band
  annotate("rect", xmin = -Inf, xmax = Inf, ymin =  0.25, ymax = 0.50,
           fill = "#FAEEDA", alpha = 0.55)      # reasonable band
p_silt <- add_infeas_wash(p_silt)
p_silt <- p_silt +
  geom_line(linewidth = 0.3, color = COL_SIL) +
  geom_point(size = 0.7, color = COL_SIL) +
  geom_point(data = subset(sil_tbl, k %in% s1_locmax),
             shape = 21, fill = "white", color = COL_SIL,
             size = 1.3, stroke = 0.45) +
  scale_y_continuous(limits = c(-0.2, 1.0),
                     breaks = c(-0.2, 0, 0.25, 0.5, 0.75, 1.0)) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "k", y = expression(bar(s)(k)),
       title = "Silhouette: local maxima",
       subtitle = paste0(
         "s_bar(k) in [-1, 1]. Open circles = local maxima at k = ",
         paste(s1_locmax, collapse = ", "),
         ". Background shades = Rousseeuw qualitative bands (context only). ",
         sprintf("At k=9: s_bar = %.3f.", sil_tbl$s_bar[match(9, sil_tbl$k)]))) +
  theme_panel()
p_silt <- mark_k9(p_silt, sil_tbl$s_bar[match(9, sil_tbl$k)])

# ---- Panel B: Gap statistic (Tibshirani SE-rule) ----
p_s2 <- ggplot(gap_tbl, aes(k, decision_signal)) +
  geom_hline(yintercept = 0, linewidth = 0.2)
p_s2 <- add_infeas_wash(p_s2)
p_s2 <- p_s2 +
  geom_col(aes(fill = passes_se_rule), width = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#1D9E75", `FALSE` = "#F0997B"),
                    na.value = "grey80", guide = "none") +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "k",
       y = "Gap(k) - [Gap(k+1) - SE(k+1)]",
       title = "Gap statistic (Tibshirani SE-rule)",
       subtitle = paste0(
         "Bars >= 0 (green) satisfy Gap(k) >= Gap(k+1) - SE(k+1). ",
         "Passing k = ", paste(gap_pass, collapse = ", "),
         sprintf(". At k=9: decision_signal = %.3f.",
                 gap_tbl$decision_signal[match(9, gap_tbl$k)]))) +
  theme_panel()
p_s2 <- mark_k9(p_s2, gap_tbl$decision_signal[match(9, gap_tbl$k)])

# ---- Panel C: W (within-cluster dispersion) and its derivatives ----
# Three curves overlaid in a single ggplot: the level (W), the slope
# magnitude (|W'|), and the curvature (W'').
#
# W(k)  = within-cluster dispersion: for each cluster r, sum of squared pairwise
#         distances divided by 2*n_r, summed over r. Monotone non-increasing in
#         k. SAME quantity that the GAP STATISTIC compares against a uniform
#         null (Gap(k) = log(W_null) - log(W_obs)) - so W, W', W'' are upstream
#         of the gap statistic, not separate ideas.
# W'(k) = first difference W(k) - W(k-1); always <= 0. |W'(k)| = magnitude of
#         improvement at one more split. We plot |W'(k)| so it's on the same
#         positive side as W and W''.
# W''(k)= second difference of W; positive = improvement decelerating; -> 0 =
#         improvement has plateaued.
#
# All three curves are min-max normalized independently to a common [0, 1]
# axis so their shapes are directly comparable on one panel. The original-unit
# range of each curve is shown in the legend at the top right.
W_min  <- min(wk_tbl$W,                   na.rm = TRUE)
W_max  <- max(wk_tbl$W,                   na.rm = TRUE)
Wp_min <- min(abs(wk_tbl$W_prime),        na.rm = TRUE)
Wp_max <- max(abs(wk_tbl$W_prime),        na.rm = TRUE)
Wpp_min <- min(wk_tbl$W_dprime,           na.rm = TRUE)
Wpp_max <- max(wk_tbl$W_dprime,           na.rm = TRUE)
nrm <- function(x, lo, hi) (x - lo) / (hi - lo)

wk_tbl_norm <- wk_tbl %>%
  mutate(
    W_norm   = nrm(W,                W_min,   W_max),
    Wp_norm  = nrm(abs(W_prime),     Wp_min,  Wp_max),
    Wpp_norm = nrm(W_dprime,         Wpp_min, Wpp_max)
  )

# Color tokens for the three overlaid curves
COL_W   <- COL_PRIM     # gray - level
COL_WP  <- "#7F77DD"    # purple - slope magnitude
COL_WPP <- COL_ELBOW    # amber - curvature

eps_norm  <- nrm(eps_W_dprime, Wpp_min, Wpp_max)
zero_norm <- nrm(0,            Wpp_min, Wpp_max)  # baseline for W'' bars

p_s3 <- ggplot(wk_tbl_norm, aes(x = k))
p_s3 <- add_infeas_wash(p_s3)
p_s3 <- p_s3 +
  # W''(k) as filled bars (amber for positive, light gray for negative)
  geom_col(data = subset(wk_tbl_norm, !is.na(W_dprime)),
           aes(y = Wpp_norm, fill = W_dprime > 0),
           alpha = 0.55, color = NA) +
  scale_fill_manual(values = c(`TRUE` = COL_WPP, `FALSE` = "#D3D1C7"),
                    guide = "none") +
  # |W'(k)| as a line+points
  geom_line(data = subset(wk_tbl_norm, !is.na(W_prime)),
            aes(y = Wp_norm), color = COL_WP, linewidth = 0.35) +
  geom_point(data = subset(wk_tbl_norm, !is.na(W_prime)),
             aes(y = Wp_norm), color = COL_WP, size = 0.7) +
  # W(k) as a line+points
  geom_line(aes(y = W_norm),  color = COL_W, linewidth = 0.4) +
  geom_point(aes(y = W_norm), color = COL_W, size = 0.7) +
  # Reference lines: zero baseline for W'' (in normalized coords) and the
  # eps band that defines the W'' -> 0 plateau
  geom_hline(yintercept = zero_norm, linewidth = 0.2, color = "grey50") +
  geom_hline(yintercept = eps_norm,  linetype = "dotted",
             color = COL_K9, linewidth = 0.2) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(breaks = K_RANGE) +
  # Per-curve range legend, anchored top-right of the panel
  annotate("text", x = max(K_RANGE), y = 1.02, hjust = 1, vjust = 1,
           color = COL_W, size = 1.5,
           label = sprintf("W(k): %.0f -> %.0f (line)", W_max, W_min)) +
  annotate("text", x = max(K_RANGE), y = 0.92, hjust = 1, vjust = 1,
           color = COL_WP, size = 1.5,
           label = sprintf("|W'(k)|: %.0f -> %.0f (line)", Wp_max, Wp_min)) +
  annotate("text", x = max(K_RANGE), y = 0.82, hjust = 1, vjust = 1,
           color = COL_WPP, size = 1.5,
           label = sprintf("W''(k): peak %.0f (bars)", Wpp_max)) +
  labs(x = "k", y = "min-max normalized to [0, 1]",
       title = "W (within-cluster dispersion) and its derivatives",
       subtitle = sprintf(paste0(
         "W(k) = within-cluster dispersion (same quantity the gap statistic ",
         "compares against a uniform null). Three curves min-max normalized ",
         "to a shared axis; legend gives original ranges. W''(k) returns to ",
         "the noise floor (eps=%.3f) at k=%d, where the rate of improvement ",
         "has plateaued."),
         eps_W_dprime, s3_plateau_k)) +
  theme_panel()
p_s3 <- mark_k9(p_s3, wk_tbl_norm$Wpp_norm[match(9, wk_tbl_norm$k)])

# ---- Panel D: Resolution ceiling (first-touch-the-floor) ----
p_s6 <- ggplot(min_size_tbl, aes(k, min_size,
                                 fill = min_size >= PHI_RESOLUTION))
p_s6 <- add_infeas_wash(p_s6)
p_s6 <- p_s6 +
  geom_col() +
  geom_hline(yintercept = PHI_RESOLUTION, color = COL_THRESH,
             linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 2, color = "#A32D2D", linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#9FE1CB", `FALSE` = COL_INFEAS),
                    guide = "none") +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "k", y = "min cluster size",
       title = "Resolution ceiling (first-touch-the-floor)",
       subtitle = sprintf(paste0(
         "Smallest cluster vs k. Dashed line: phi=%d (resolution floor). ",
         "min_size first reaches the floor at k=%d. Beyond k=%d, the floor is ",
         "broken (singletons appear)."),
         PHI_RESOLUTION, s6_k, max(min_size_tbl$k[min_size_tbl$feasible]))) +
  theme_panel()
p_s6 <- mark_k9(p_s6, min_size_tbl$min_size[match(9, min_size_tbl$k)])

# ---- Panel E: Stability (bootstrap consensus, gated at the resolution-ceiling-entry) ----
p_s7 <- ggplot(stab_tbl, aes(k, mean_jaccard)) +
  geom_hline(yintercept = JACCARD_STABLE, linetype = "dashed",
             color = "#185FA5", linewidth = 0.25)
p_s7 <- add_infeas_wash(p_s7)
p_s7 <- p_s7 +
  geom_errorbar(aes(ymin = mean_jaccard - sd_jaccard,
                    ymax = mean_jaccard + sd_jaccard),
                width = 0.15, linewidth = 0.2, color = COL_STAB) +
  geom_line(color = COL_STAB, linewidth = 0.25) +
  geom_point(color = COL_STAB, size = 0.6) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "k", y = "mean Jaccard",
       title = "Stability (bootstrap consensus)",
       subtitle = sprintf(paste0(
         "Bootstrap mean Jaccard +/- 1 sd vs k. Dashed line: stable >= %.2f. ",
         "At k=9 stability holds (mean Jaccard = %.3f) and persists for k=10."),
         JACCARD_STABLE, stab_tbl$mean_jaccard[match(9, stab_tbl$k)])) +
  theme_panel()
p_s7 <- mark_k9(p_s7, stab_tbl$mean_jaccard[match(9, stab_tbl$k)])

# Compose 3 rows x 2 cols
# Layout: 3x2 grid (six equal-sized cells); the sixth cell is intentionally
# empty so that all five content panels share the same dimensions.
#   row 1 - silhouette                | gap statistic
#   row 2 - W and its derivatives     | resolution ceiling
#   row 3 - stability                 | (empty)
combined <- p_silt + p_s2 + p_s3 + p_s6 + p_s7 +
  plot_layout(design = "AB\nCD\nE#") +
  plot_annotation(
    tag_levels = "A",
    title = "k (number of cluster) optimization (93 genera x 11 age groups)",
    subtitle = "k = 9 is the chosen number of clusters; teal disc + dashed teal vline mark k=9 on each panel.",
    caption = paste0("Gray-washed columns = infeasible (min_size < 2). ",
                     "egg::theme_article(base_size = 5).")
  ) &
  theme(plot.title    = element_text(size = THEME_BASE + 2, face = "bold"),
        plot.subtitle = element_text(size = THEME_BASE, color = "grey30"),
        plot.caption  = element_text(size = THEME_BASE - 0.5, color = "grey45",
                                     hjust = 0))

# PDF dimensions: 7.5 x 8.5 in scaled by 0.7x in both axes.
PDF_W <- 7.5 * 0.7
PDF_H <- 8.5 * 0.7
out_pdf <- file.path(PLOTS_DIR, "k_optimization_panels.pdf")
if (capabilities("cairo")) {
  ggsave(out_pdf, plot = combined,
         width = PDF_W, height = PDF_H, units = "in", device = cairo_pdf)
} else {
  ggsave(out_pdf, plot = combined,
         width = PDF_W, height = PDF_H, units = "in",
         device = grDevices::pdf, useDingbats = FALSE)
}
cat(sprintf("[saved] %s (%.2f x %.2f in)\n", out_pdf, PDF_W, PDF_H))

# =============================================================================
# Section 4b - Stability + bootstrap consensus heatmaps.
# 7.0 x 5.0 in landscape; saved as consensus_stability.pdf.
# =============================================================================
cat("\n=== SECTION 4b: stability + consensus heatmaps ===\n")

cm_pal <- scales::col_numeric(palette = c("white", "#0C447C"),
                              domain = c(0, 1))

build_heatmap <- function(K) {
  cm <- stab_results[[as.character(K)]]$consensus
  ord <- hclust(as.dist(1 - cm), method = "average")$order
  cm <- cm[ord, ord]
  n <- nrow(cm)
  df <- data.frame(
    Var1 = rep(seq_len(n), times = n),
    Var2 = rep(seq_len(n), each  = n),
    fill_hex = cm_pal(as.vector(cm))
  )

  add_strip <- (K == 9L)
  ylim_low  <- if (add_strip) -3.5 else 0.5

  p <- ggplot(df, aes(Var1, Var2, fill = fill_hex)) +
    geom_raster() +
    scale_fill_identity(guide = "none") +
    coord_fixed(ylim = c(ylim_low, n + 0.5), clip = "off") +
    labs(title = sprintf("k = %d", K)) +
    egg::theme_article(base_size = 4) +
    theme(axis.text  = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 4, hjust = 0.5))

  if (!is.na(s7_k) && K == s7_k) {
    p <- p + annotate("rect",
                      xmin = 0.5, xmax = n + 0.5,
                      ymin = 0.5, ymax = n + 0.5,
                      fill = NA, color = COL_K9, linewidth = 0.6)
  }

  if (add_strip) {
    ref_cl <- cluster_assignments[["9"]]
    letters_ord <- unname(cluster_mapping_main[as.character(ref_cl)][ord])
    strip_df <- data.frame(
      x = seq_len(n),
      y = -1.8,
      fill_hex = unname(cluster_colors_a_to_i[letters_ord]),
      stringsAsFactors = FALSE
    )
    p <- p + geom_tile(data = strip_df,
                       aes(x = x, y = y, fill = fill_hex),
                       height = 1.5, color = NA, inherit.aes = FALSE) +
      annotate("text", x = n / 2 + 0.5, y = -3.2,
               label = "cluster a..i",
               size = 1.2, color = "grey25")
  }
  p
}

hm_list   <- lapply(K_RANGE, build_heatmap)
left_grid <- patchwork::wrap_plots(hm_list, ncol = 4)

sub_extra <- if (any(diff(stab_tbl$mean_jaccard) > 0)) "" else
  "  |  stability declines monotonically - naive plateau-entry would trivially favor small k here."

p_right <- ggplot(stab_tbl, aes(k, mean_jaccard)) +
  geom_hline(yintercept = JACCARD_STABLE, linetype = "dashed",
             color = "#185FA5", linewidth = 0.25) +
  geom_hline(yintercept = JACCARD_UNSTABLE, linetype = "dotted",
             color = "#A32D2D", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean_jaccard - sd_jaccard,
                    ymax = mean_jaccard + sd_jaccard),
                width = 0.15, linewidth = 0.2, color = COL_STAB) +
  geom_line(color = COL_STAB, linewidth = 0.25) +
  geom_point(color = COL_STAB, size = 0.6) +
  annotate("text", x = max(K_RANGE), y = JACCARD_STABLE,
           label = sprintf("stable >= %.2f", JACCARD_STABLE),
           size = 1.4, hjust = 1, vjust = -0.4, color = "#185FA5") +
  annotate("text", x = max(K_RANGE), y = JACCARD_UNSTABLE,
           label = sprintf("unstable < %.2f", JACCARD_UNSTABLE),
           size = 1.4, hjust = 1, vjust = -0.4, color = "#A32D2D") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = K_RANGE,
                     limits = range(K_RANGE) + c(-0.4, 0.6)) +
  labs(x = "k", y = "mean Jaccard (per-cluster, averaged)") +
  egg::theme_article(base_size = THEME_BASE) +
  theme(axis.title       = element_text(size = THEME_BASE),
        axis.text        = element_text(size = THEME_BASE - 0.5),
        plot.title       = element_text(size = THEME_BASE + 1, face = "plain"),
        plot.subtitle    = element_text(size = THEME_BASE, color = "grey40"))

if (!is.na(s7_k)) {
  s7_y <- stab_tbl$mean_jaccard[match(s7_k, stab_tbl$k)]
  p_right <- p_right +
    geom_vline(xintercept = s7_k, linetype = "dashed",
               linewidth = 0.25, color = COL_K9) +
    geom_point(data = data.frame(k = s7_k, y = s7_y),
               aes(k, y), inherit.aes = FALSE,
               size = 1.4, color = COL_K9) +
    annotate("text", x = s7_k, y = s7_y,
             label = sprintf("  k* = %d", s7_k),
             hjust = 0, vjust = -0.6, size = 1.6, color = COL_K9)
}

p_s7_combined <- patchwork::wrap_plots(left_grid, p_right,
                                       ncol = 2, widths = c(1, 1)) +
  patchwork::plot_annotation(
    title = "Bootstrap consensus and stability",
    subtitle = paste0(
      sprintf(paste0(
        "Per-pair co-clustering rate across %d bootstrap subsamples ",
        "(left, k=2..%d) and mean Jaccard +/- 1 sd vs k (right). ",
        "k=9 panel boxed in teal."),
        B_STABILITY, max(K_RANGE)),
      sub_extra),
    caption = sprintf(paste0(
      "B = %d resamples; subsample fraction = %.2f. ",
      "k=9 panel annotated with cluster letters."),
      B_STABILITY, SUBSAMPLE_FRAC))

s7_pdf <- file.path(PLOTS_DIR, "consensus_stability.pdf")
if (capabilities("cairo")) {
  ggsave(s7_pdf, plot = p_s7_combined,
         width = 7.0, height = 5.0, units = "in", device = cairo_pdf)
} else {
  ggsave(s7_pdf, plot = p_s7_combined,
         width = 7.0, height = 5.0, units = "in",
         device = grDevices::pdf, useDingbats = FALSE)
}
cat(sprintf("[saved] %s (7.0 x 5.0 in)\n", s7_pdf))

# =============================================================================
# Section 5 - markdown summary alongside the PDF
# =============================================================================
md <- c(
  "# k optimization summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  paste0("Search range: k = ", min(K_RANGE), "..", max(K_RANGE),
         "; n_genera = ", n_genera,
         "; B_GAP = ", B_GAP, "; B_STABILITY = ", B_STABILITY, "."),
  "",
  "**Chosen number of clusters: k = 9.**",
  "",
  "## Evidence by panel",
  "",
  "Each panel of `k_optimization_panels.pdf` shows the data behind a complementary criterion. The reading at k = 9 is summarized below.",
  "",
  "| Panel | Criterion | Reading at k = 9 |",
  "|---|---|---|",
  paste(sprintf("| %s | %s | %s |",
                decisions$panel,
                decisions$criterion,
                decisions$reading_at_k9),
        collapse = "\n"),
  "",
  paste0("Bootstrap consensus heatmaps for k = 2..", max(K_RANGE),
         " are in `consensus_stability.pdf`; ",
         "the k = 9 panel is boxed in teal."),
  ""
)
writeLines(md, file.path(viz_out_path, "k_optimization_summary.md"))
cat("[saved] k_optimization_summary.md\n")

cat("\n=== k-optimization section complete ===\n")

# =============================================================================
# Section 6 - GOLD-trait ecology and nativity composition
# =============================================================================

set.seed(42)

base_path <- PROJECT_ROOT

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(egg)
  library(data.table)
  library(stringr)
  library(grid)
})

viz_out_root <- file.path(base_path, "results/analyses_results/5_age_analyses_out_additional")
dir.create(viz_out_root, recursive = TRUE, showWarnings = FALSE)

derive_genus_cluster_assignments_k9 <- function(repo_root, viz_root) {
  # Derive k=9 genus cluster assignments using the standard preprocessing
  # (genus-level relative abundances, prevalence filter,
  # age-binned means, z-scoring across age bins, Ward.D2 on euclidean distances).
  intermediate_files_path <- file.path(
    repo_root,
    "results/intermediate"
  )
  ps_path <- file.path(intermediate_files_path, "ubiome_relative_none_updated.rds")
  if (!file.exists(ps_path)) stop("Missing required input phyloseq RDS: ", ps_path)
  ps <- readRDS(ps_path)

  ps_genus <- tax_glom(ps, taxrank = "Genus")
  prevalence_threshold <- 0.1 * nsamples(ps_genus)
  ps_filt <- filter_taxa(ps_genus, function(x) sum(x > 0) >= prevalence_threshold, TRUE)

  # Age bins
  age_breaks <- c(14, seq(20, 85, by = 5))
  age_labels <- c("14-19", paste(seq(20, 75, by = 5), seq(24, 79, by = 5), sep = "-"), "80-85")
  sample_data(ps_filt)$age_group <- cut(sample_data(ps_filt)$Age, breaks = age_breaks, labels = age_labels, right = FALSE)

  otu_mat <- as(otu_table(ps_filt), "matrix")
  sample_df <- as(sample_data(ps_filt), "data.frame")
  otu_df <- as.data.frame(t(otu_mat))
  otu_df$age_group <- sample_df$age_group

  mean_abundance <- otu_df %>%
    group_by(age_group) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  mean_abundance_mat <- as.matrix(mean_abundance[, -1])
  rownames(mean_abundance_mat) <- mean_abundance$age_group

  mean_abundance_scaled_full <- t(scale(t(t(mean_abundance_mat))))
  set.seed(42)
  dist_mat <- dist(mean_abundance_scaled_full, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  cl <- cutree(hc, k = 9)

  # Map numeric cluster ids -> letters using the canonical relabeling so that
  # cluster letters refer to the same genera as in 5_age_analyses_full.R and
  # cluster_mapping_main used by the consensus-heatmap annotation strip.
  cl_ids <- sort(unique(as.integer(cl)))
  if (length(cl_ids) != 9) stop("Expected 9 clusters from cutree(k=9); got: ", length(cl_ids))
  cl_map <- c("1" = "i", "2" = "d", "3" = "h", "4" = "b", "5" = "c",
              "6" = "f", "7" = "a", "8" = "e", "9" = "g")
  cl_letters <- unname(cl_map[as.character(as.integer(cl))])

  tax <- as(tax_table(ps_filt), "matrix")
  phylum <- if ("Phylum" %in% colnames(tax)) as.character(tax[, "Phylum"]) else rep(NA_character_, nrow(tax))

  out <- data.frame(
    genus = taxa_names(ps_filt),
    cluster = cl_letters,
    phylum = phylum,
    stringsAsFactors = FALSE
  )

  dir.create(viz_root, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(viz_root, "genus_cluster_assignments.csv")
  write.csv(out, out_csv, row.names = FALSE)
  out_csv
}

run_ecology <- function(repo_root = base_path, viz_root = viz_out_root) {
  viz_out_path <- file.path(viz_root, "ecology")
  dir.create(viz_out_path, recursive = TRUE, showWarnings = FALSE)

  save_pdf_figure <- function(file_name, width_mm = 180, height_mm = 170) {
    pdf(
      file = file.path(viz_out_path, file_name),
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      family = "Helvetica",
      colormodel = "rgb",
      pointsize = 5,
      useDingbats = FALSE
    )
  }

  cluster_assign_path <- file.path(viz_root, "genus_cluster_assignments.csv")
  if (!file.exists(cluster_assign_path)) {
    message("Missing genus_cluster_assignments.csv; deriving k=9 assignments now...")
    cluster_assign_path <- derive_genus_cluster_assignments_k9(repo_root, viz_root)
  }
  gold_map_path <- file.path(
    repo_root,
    "results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv"
  )
  if (!file.exists(cluster_assign_path)) stop("Missing required input: ", cluster_assign_path)
  if (!file.exists(gold_map_path)) stop("Missing required input: ", gold_map_path)

  clusters_df <- read.csv(cluster_assign_path, check.names = FALSE) %>%
    mutate(
      genus = as.character(genus),
      cluster = as.character(cluster),
      phylum = as.character(phylum)
    )

  gold_df <- read.csv(gold_map_path) %>%
    mutate(
      otu = as.character(otu),
      OXYGEN_INDEX = suppressWarnings(as.numeric(OXYGEN_INDEX)),
      STAIN_INDEX = suppressWarnings(as.numeric(STAIN_INDEX)),
      MOTILITY_INDEX = suppressWarnings(as.numeric(MOTILITY_INDEX)),
      SPORULATION_INDEX = suppressWarnings(as.numeric(SPORULATION_INDEX))
    )

  joined <- clusters_df %>%
    left_join(
      gold_df %>% select(
        otu, parsed_genus, Genus,
        OXYGEN_REQUIREMENT, GRAM_STAIN, MOTILITY, SPORULATION,
        OXYGEN_INDEX, STAIN_INDEX, MOTILITY_INDEX, SPORULATION_INDEX
      ),
      by = c("genus" = "otu")
    )

  trait_cols <- c("OXYGEN_INDEX", "STAIN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")

  cluster_color_palette <- c(
    "#ffaabb", "#ee8866", "#eedd88", "#bbcc33", "#aaaa00",
    "#44bb99", "#99ddff", "#77aadd",  "#999999", "#dddddd",
    "#ff6b6b", "#4ecdc4", "#45b7d1", "#96ceb4", "#feca57",
    "#ff9ff3", "#54a0ff", "#5f27cd", "#00d2d3", "#ff9f43"
  )
  cluster_letters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
  cluster_colors <- setNames(cluster_color_palette[1:9], cluster_letters)

  phylum_colors <- c(
    "Firmicutes"             = "#F38400",
    "Bacteroidetes"          = "#0067A5",
    "Actinobacteria"         = "#8DB600",
    "Proteobacteria"         = "#E68FAC",
    "Fusobacteria"           = "#BE0032",
    "Spirochaetae"           = "#F3C300",
    "Cyanobacteria"          = "#875692",
    "Acidobacteria"          = "#F6A600",
    "Candidate division SR1" = "#2B3D26",
    "Planctomycetes"         = "#332288",
    "Saccharibacteria"       = "#B3446C",
    "Synergistetes"          = "#A1CAF1",
    "Tenericutes"            = "#654522",
    "Verrucomicrobia"        = "#C2B280",
    "unclassified"           = "#DDDDDD",
    "NA"                     = "#DDDDDD",
    "TM7"                    = "#000000"
  )

  missing_phyla <- setdiff(na.omit(unique(joined$phylum)), names(phylum_colors))
  if (length(missing_phyla)) {
    stop("Phylum(s) missing from master palette: ", paste(missing_phyla, collapse = ", "))
  }

  trait_specs <- list(
    OXYGEN_INDEX = "Oxygen tolerance index (0=anaerobe, 1=aerobe)",
    STAIN_INDEX = "Gram stain index (0=Gram-, 1=Gram+)",
    MOTILITY_INDEX = "Motility index (0=nonmotile, 1=motile)",
    SPORULATION_INDEX = "Sporulation index (0=nonsporulating, 1=sporulating)"
  )

  joined$cluster <- factor(joined$cluster, levels = cluster_letters)

  cluster_means_trait <- joined %>%
    dplyr::filter(!is.na(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(trait_cols), \(x) mean(x, na.rm = TRUE)),
      .groups = "drop"
    )
  between_cluster_var_rank <- sapply(trait_cols, function(tr) {
    v <- cluster_means_trait[[tr]]
    stats::var(v, na.rm = TRUE)
  })
  between_var_tbl <- data.frame(
    trait = trait_cols,
    between_cluster_variance_of_cluster_means = between_cluster_var_rank,
    stringsAsFactors = FALSE
  )
  total_between_var_sum <- sum(between_var_tbl$between_cluster_variance_of_cluster_means)
  between_var_tbl$pct_of_total_between_cluster_variance <- if (total_between_var_sum > 0) {
    100 * between_var_tbl$between_cluster_variance_of_cluster_means / total_between_var_sum
  } else {
    rep(NA_real_, nrow(between_var_tbl))
  }
  between_var_tbl <- between_var_tbl %>%
    dplyr::arrange(dplyr::desc(.data$between_cluster_variance_of_cluster_means)) %>%
    dplyr::mutate(rank_across_axes = dplyr::row_number())

  bv_lookup <- stats::setNames(between_var_tbl$between_cluster_variance_of_cluster_means, between_var_tbl$trait)
  pct_lookup <- stats::setNames(between_var_tbl$pct_of_total_between_cluster_variance, between_var_tbl$trait)

  subtitle_append_between_cluster_pct <- function(sub) {
    pct_line <- paste(sprintf("%s %.1f%%", trait_cols, unname(pct_lookup[trait_cols])), collapse = "; ")
    paste0(sub, "\nBetween-cluster variance (% of total across four GOLD indices): ", pct_line)
  }

  bubble_legend_title_mean <- function(trait_key) {
    paste0(
      "Mean ", trait_key,
      " (", sprintf("%.1f%%", pct_lookup[[trait_key]]), " of total between-cluster variance)"
    )
  }

  bubble_legend_title_genus <- function(trait_key) {
    paste0(
      "Genus-level ", trait_key,
      " (", sprintf("%.1f%%", pct_lookup[[trait_key]]), " of total between-cluster variance)"
    )
  }

  plot_list <- list()
  for (trait in names(trait_specs)) {
    label <- trait_specs[[trait]]
    df <- joined %>%
      select(genus, cluster, phylum, !!trait) %>%
      rename(value = !!trait) %>%
      filter(!is.na(cluster))

    df_non_na <- df %>% filter(!is.na(value))
    if (nrow(df_non_na) < 3) next

    bv_val <- unname(bv_lookup[[trait]])

    p <- ggplot(df_non_na, aes(x = cluster, y = value)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.25, fill = "gray90") +
      geom_jitter(aes(color = phylum), width = 0.12, height = 0, alpha = 0.65, size = 0.6) +
      scale_color_manual(values = phylum_colors) +
      egg::theme_article() +
      theme(
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 4, hjust = 0.5, color = "gray30"),
        axis.title = element_text(size = 5),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 5),
        legend.position = "none"
      ) +
      labs(
        title = trait,
        subtitle = paste0(
          "Between-cluster var(cluster means)=", signif(bv_val, 4)
        ),
        x = "Cluster (A-I)",
        y = label
      ) +
      scale_y_continuous(limits = c(0, 1))

    plot_list[[trait]] <- p
  }

  if (length(plot_list)) {
    ordered_traits <- intersect(names(trait_specs), names(plot_list))
    grobs <- lapply(ordered_traits, function(t) plot_list[[t]])

    trait_distributions_input <- joined %>%
      dplyr::select(genus, cluster, phylum, dplyr::all_of(trait_cols)) %>%
      dplyr::filter(!is.na(cluster)) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(trait_cols),
        names_to = "trait",
        values_to = "value"
      )
    write.csv(
      trait_distributions_input,
      file.path(viz_out_path, "raw_input_trait_distributions_four_panel.csv"),
      row.names = FALSE
    )

    save_pdf_figure("trait_distributions_four_panel_by_cluster.pdf", width_mm = 170, height_mm = 140)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
    for (i in seq_along(grobs)) {
      row <- ((i - 1) %/% 2) + 1
      col <- ((i - 1) %% 2) + 1
      print(grobs[[i]], vp = grid::viewport(layout.pos.row = row, layout.pos.col = col))
    }
    dev.off()
  }

  # ---- Ternary cluster-mean summaries ----

  make_ternary_projection <- function(df, a_raw, b_raw, c_raw) {
    total <- df[[a_raw]] + df[[b_raw]] + df[[c_raw]]
    total[total == 0] <- 1e-6
    df$a <- df[[a_raw]] / total
    df$b <- df[[b_raw]] / total
    df$c <- df[[c_raw]] / total
    df$x <- df$b + 0.5 * df$c
    df$y <- (sqrt(3) / 2) * df$c
    df
  }

  ternary_triangle <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(0, 0, sqrt(3) / 2, 0)
  )

  make_ternary_grid_segments <- function(tick_values, level = c("major", "minor")) {
    level <- match.arg(level)
    h <- sqrt(3) / 2
    tick_values <- tick_values[is.finite(tick_values) & tick_values > 0 & tick_values < 1]
    if (length(tick_values) == 0) {
      return(data.frame(
        x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0),
        level = character(0), stringsAsFactors = FALSE
      ))
    }
    rows <- lapply(tick_values, function(t) {
      dplyr::bind_rows(
        data.frame(
          x = 0.5 * (1 - t), y = h * (1 - t), xend = 1 - t, yend = 0,
          stringsAsFactors = FALSE
        ),
        data.frame(
          x = t, y = 0, xend = 0.5 * t + 0.5, yend = h * (1 - t),
          stringsAsFactors = FALSE
        ),
        data.frame(
          x = 0.5 * t, y = h * t, xend = 1 - 0.5 * t, yend = h * t,
          stringsAsFactors = FALSE
        )
      )
    })
    out <- dplyr::bind_rows(rows)
    out$level <- level
    out
  }

  ternary_grid_df <- dplyr::bind_rows(
    make_ternary_grid_segments(seq(0.2, 0.8, by = 0.2), "major"),
    make_ternary_grid_segments(
      setdiff(seq(0.1, 0.9, by = 0.1), seq(0, 1, by = 0.2)),
      "minor"
    )
  )

  make_ternary_axes_guides <- function(axis_titles, tick_by = 0.2) {
    stopifnot(length(axis_titles) == 3)

    ticks <- seq(0, 1, by = tick_by)
    h <- sqrt(3) / 2
    ticks_mid <- ticks[ticks > 0 & ticks < 1]

    df_bottom <- data.frame(
      x = ticks_mid,
      y = 0,
      lab = format(ticks_mid, nsmall = 1, trim = TRUE),
      edge = "bottom",
      stringsAsFactors = FALSE
    )
    df_left <- data.frame(
      x = 0.5 * ticks_mid,
      y = h * ticks_mid,
      lab = format(ticks_mid, nsmall = 1, trim = TRUE),
      edge = "left",
      stringsAsFactors = FALSE
    )
    df_right <- data.frame(
      x = 1 - 0.5 * ticks_mid,
      y = h * ticks_mid,
      lab = format(ticks_mid, nsmall = 1, trim = TRUE),
      edge = "right",
      stringsAsFactors = FALSE
    )
    tick_df <- dplyr::bind_rows(df_bottom, df_left, df_right)

    tick_df <- tick_df %>%
      mutate(
        x_off = dplyr::case_when(
          edge == "bottom" ~ x,
          edge == "left" ~ x - 0.035,
          edge == "right" ~ x + 0.035,
          TRUE ~ x
        ),
        y_off = dplyr::case_when(
          edge == "bottom" ~ y - 0.045,
          edge == "left" ~ y,
          edge == "right" ~ y,
          TRUE ~ y
        ),
        angle = dplyr::case_when(
          edge == "bottom" ~ 0,
          edge == "left" ~ 60,
          edge == "right" ~ -60,
          TRUE ~ 0
        )
      )

    title_df <- data.frame(
      x = c(0.5, 0.25, 0.75),
      y = c(-0.10, h / 2, h / 2),
      x_off = c(0.5, 0.25 - 0.06, 0.75 + 0.06),
      y_off = c(-0.10, h / 2, h / 2),
      lab = axis_titles,
      angle = c(0, 60, -60),
      stringsAsFactors = FALSE
    )

    list(ticks = tick_df, titles = title_df)
  }

  ternary_vertex_frame <- function(vertex_labels) {
    data.frame(
      x = c(0, 1, 0.5),
      y = c(-0.07, -0.07, sqrt(3) / 2 + 0.07),
      lab = vertex_labels,
      stringsAsFactors = FALSE
    )
  }

  wrap_vertex_label <- function(trait_key) {
    lab <- if (trait_key %in% names(trait_specs)) trait_specs[[trait_key]] else trait_key
    paste(strwrap(lab, width = 30), collapse = "\n")
  }

  wrap_vertex_label_var <- function(trait_key) {
    pct <- pct_lookup[[trait_key]]
    if (length(pct) == 0 || is.na(pct)) pct <- NA_real_
    paste0(
      wrap_vertex_label(trait_key),
      "\n(",
      sprintf("%.1f%%", pct),
      " of total between-cluster var.)"
    )
  }

  render_ternary_cluster_mean_summary <- function(
      cluster_means,
      title,
      subtitle,
      file_name,
      color_mode = c("cluster", "phylum"),
      vertex_labels,
      bubble_trait_name
  ) {
    color_mode <- match.arg(color_mode)
    vdf <- ternary_vertex_frame(vertex_labels)
    axes_guides <- make_ternary_axes_guides(axis_titles = vertex_labels, tick_by = 0.2)
    leg_title_color <- if (color_mode == "cluster") "Cluster" else "Dominant phylum"

    p <- ggplot() +
      geom_path(data = ternary_triangle, aes(x = x, y = y), linewidth = 0.3, color = "gray30") +
      geom_segment(
        data = dplyr::filter(ternary_grid_df, .data$level == "minor"),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        color = "gray93",
        linewidth = 0.12,
        lineend = "round"
      ) +
      geom_segment(
        data = dplyr::filter(ternary_grid_df, .data$level == "major"),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        color = "gray86",
        linewidth = 0.18,
        lineend = "round"
      ) +
      geom_text(
        data = axes_guides$ticks,
        aes(x = x_off, y = y_off, label = lab, angle = angle),
        size = 2.4,
        color = "gray35"
      ) +
      geom_text(
        data = axes_guides$titles,
        aes(x = x_off, y = y_off, label = lab, angle = angle),
        size = 2.6,
        fontface = "bold",
        color = "gray20"
      ) +
      geom_text(
        data = vdf,
        aes(x = x, y = y, label = lab),
        size = 2.8,
        fontface = "bold",
        lineheight = 0.95
      )

    if (color_mode == "cluster") {
      p <- p +
        geom_point(
          data = cluster_means,
          aes(x = x, y = y, color = cluster, size = bubble),
          alpha = 0.85
        ) +
        scale_color_manual(values = cluster_colors, drop = FALSE, name = leg_title_color)
    } else {
      p <- p +
        geom_point(
          data = cluster_means,
          aes(x = x, y = y, color = dominant_phylum, size = bubble),
          alpha = 0.85
        ) +
        scale_color_manual(values = phylum_colors, drop = FALSE, name = leg_title_color)
    }

    p <- p +
      scale_size_continuous(name = bubble_legend_title_mean(bubble_trait_name), range = c(1.5, 5.0)) +
      coord_equal(xlim = c(-0.08, 1.08), ylim = c(-0.17, sqrt(3) / 2 + 0.17)) +
      egg::theme_article() +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 4, hjust = 0.5, color = "gray30"),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5)
      ) +
      labs(title = title, subtitle = subtitle)

    save_pdf_figure(file_name, width_mm = 170, height_mm = 150)
    print(p)
    dev.off()
  }

  # --- Additional ternary families ---
  # Variance-ranked centroid summary + connected means,
  # plus systematic families where each index is the bubble once.

  render_ternary_connected_cluster_means <- function(
      cluster_means,
      title,
      subtitle,
      file_name,
      vertex_labels,
      bubble_trait_name
  ) {
    cluster_means <- cluster_means %>%
      mutate(cluster_chr = as.character(cluster_chr)) %>%
      arrange(cluster_chr)

    vdf <- ternary_vertex_frame(vertex_labels)
    axes_guides <- make_ternary_axes_guides(axis_titles = vertex_labels, tick_by = 0.2)

    p <- ggplot() +
      geom_path(data = ternary_triangle, aes(x = x, y = y), linewidth = 0.3, color = "gray30") +
      geom_segment(
        data = dplyr::filter(ternary_grid_df, .data$level == "minor"),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        color = "gray93",
        linewidth = 0.12,
        lineend = "round"
      ) +
      geom_segment(
        data = dplyr::filter(ternary_grid_df, .data$level == "major"),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        color = "gray88",
        linewidth = 0.2,
        lineend = "round"
      ) +
      geom_text(
        data = axes_guides$titles,
        aes(x = x_off, y = y_off, label = lab, angle = angle),
        size = 2.2,
        lineheight = 0.9,
        fontface = "bold",
        inherit.aes = FALSE
      ) +
      geom_text(
        data = axes_guides$ticks,
        aes(x = x_off, y = y_off, label = lab, angle = angle),
        size = 2.0,
        color = "gray25",
        inherit.aes = FALSE
      ) +
      geom_text(
        data = vdf,
        aes(x = x, y = y, label = lab),
        size = 2.4,
        lineheight = 0.9,
        inherit.aes = FALSE
      ) +
      geom_path(data = cluster_means, aes(x = x, y = y), linewidth = 0.4, color = "gray30") +
      geom_point(data = cluster_means, aes(x = x, y = y, size = bubble, color = cluster_chr), alpha = 0.9) +
      geom_text(data = cluster_means, aes(x = x, y = y, label = cluster_chr), size = 2.6, vjust = -0.8) +
      scale_color_manual(values = cluster_colors, name = "Cluster") +
      scale_size_continuous(range = c(1.5, 8), name = bubble_legend_title_mean(bubble_trait_name)) +
      coord_equal(xlim = c(-0.12, 1.12), ylim = c(-0.14, sqrt(3) / 2 + 0.16), expand = FALSE) +
      egg::theme_article() +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 4, hjust = 0.5, color = "gray30"),
        plot.caption = element_text(size = 3.5, hjust = 0, color = "gray35"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      labs(
        title = title,
        subtitle = subtitle_append_between_cluster_pct(subtitle),
        caption = paste0(
          "Cluster means of GOLD indices; simplex axes; path connects a..i. ",
          "Grid: constant-share lines (10% minor / 20% major) in the simplex."
        )
      ) +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 3, alpha = 1)))

    save_pdf_figure(file_name, width_mm = 150, height_mm = 125)
    print(p)
    dev.off()
  }

  # Variance-ranked family: bubble = highest between-cluster var; axes = next three
  if (nrow(between_var_tbl) >= 4) {
    axis1_bv <- between_var_tbl$trait[1]
    axes_bv <- between_var_tbl$trait[2:4]

    cluster_means_bv <- joined %>%
      group_by(cluster) %>%
      summarise(
        n_genera = n(),
        bubble = mean(.data[[axis1_bv]], na.rm = TRUE),
        a_raw = mean(.data[[axes_bv[1]]], na.rm = TRUE),
        b_raw = mean(.data[[axes_bv[2]]], na.rm = TRUE),
        c_raw = mean(.data[[axes_bv[3]]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(!is.na(cluster))

    cluster_means_bv <- make_ternary_projection(cluster_means_bv, "a_raw", "b_raw", "c_raw")
    cluster_means_bv$cluster_chr <- as.character(cluster_means_bv$cluster)

    write.csv(
      cluster_means_bv,
      file.path(viz_out_path, "raw_input_ternary_summary_variance_rank.csv"),
      row.names = FALSE
    )
    write.csv(
      cluster_means_bv,
      file.path(viz_out_path, "raw_input_connected_cluster_means_variance_rank.csv"),
      row.names = FALSE
    )

    tern_subtitle_bv <- paste0(
      "Bubble size: ", axis1_bv, "; axes: ", paste(axes_bv, collapse = ", "),
      " (bubble chosen by highest between-cluster variance of cluster means)"
    )

    vertex_labs_bv <- c(
      wrap_vertex_label_var(axes_bv[1]),
      wrap_vertex_label_var(axes_bv[2]),
      wrap_vertex_label_var(axes_bv[3])
    )

    render_ternary_cluster_mean_summary(
      cluster_means = cluster_means_bv,
      title = "Cluster mean profile (ternary): bubble by between-cluster variance rank",
      subtitle = tern_subtitle_bv,
      file_name = "ternary_centroid_summary_cluster_colors_variance_rank.pdf",
      color_mode = "cluster",
      vertex_labels = vertex_labs_bv,
      bubble_trait_name = axis1_bv
    )

    render_ternary_connected_cluster_means(
      cluster_means = cluster_means_bv,
      title = "Connected cluster-mean profiles (variance-ranked bubble trait)",
      subtitle = tern_subtitle_bv,
      file_name = "ternary_connected_cluster_means_variance_rank.pdf",
      vertex_labels = vertex_labs_bv,
      bubble_trait_name = axis1_bv
    )
  }

  # Systematic family: each index is bubble once
  ecology_trait_order <- c("OXYGEN_INDEX", "STAIN_INDEX", "MOTILITY_INDEX", "SPORULATION_INDEX")
  trait_size_tag <- function(tr) sub("_INDEX$", "", tr)

  if (all(ecology_trait_order %in% names(joined))) {
    for (bubble_trait in ecology_trait_order) {
      axis_traits <- ecology_trait_order[ecology_trait_order != bubble_trait]
      size_tag <- trait_size_tag(bubble_trait)

      vertex_labs_sys <- c(
        wrap_vertex_label_var(axis_traits[1]),
        wrap_vertex_label_var(axis_traits[2]),
        wrap_vertex_label_var(axis_traits[3])
      )
      tern_sub_sys <- paste0(
        "Systematic encoding: bubble = ", bubble_trait, "; ternary axes = ",
        paste(axis_traits, collapse = ", ")
      )

      cluster_means_sys <- joined %>%
        group_by(cluster) %>%
        summarise(
          n_genera = n(),
          bubble = mean(.data[[bubble_trait]], na.rm = TRUE),
          a_raw = mean(.data[[axis_traits[1]]], na.rm = TRUE),
          b_raw = mean(.data[[axis_traits[2]]], na.rm = TRUE),
          c_raw = mean(.data[[axis_traits[3]]], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(!is.na(cluster))

      cluster_means_sys <- make_ternary_projection(cluster_means_sys, "a_raw", "b_raw", "c_raw")
      cluster_means_sys$cluster_chr <- as.character(cluster_means_sys$cluster)

      write.csv(
        cluster_means_sys,
        file.path(viz_out_path, paste0("raw_input_ternary_summary_size_", size_tag, ".csv")),
        row.names = FALSE
      )
      write.csv(
        cluster_means_sys,
        file.path(viz_out_path, paste0("raw_input_connected_cluster_means_size_", size_tag, ".csv")),
        row.names = FALSE
      )

      fn_cent_cl <- paste0("ternary_centroid_summary_cluster_colors_size_", size_tag, ".pdf")
      fn_conn <- paste0("ternary_connected_cluster_means_size_", size_tag, ".pdf")

      render_ternary_cluster_mean_summary(
        cluster_means = cluster_means_sys,
        title = paste0("Cluster mean summary (systematic): bubble = ", bubble_trait),
        subtitle = tern_sub_sys,
        file_name = fn_cent_cl,
        color_mode = "cluster",
        vertex_labels = vertex_labs_sys,
        bubble_trait_name = bubble_trait
      )

      render_ternary_connected_cluster_means(
        cluster_means = cluster_means_sys,
        title = paste0("Connected cluster means (systematic): bubble = ", bubble_trait),
        subtitle = tern_sub_sys,
        file_name = fn_conn,
        vertex_labels = vertex_labs_sys,
        bubble_trait_name = bubble_trait
      )
    }
  } else {
    warning("Systematic ternary block skipped: not all GOLD index columns present in joined table.")
  }

  # Manifest (output file -> raw plotting input CSV)
  manifest_rows <- list(
    data.frame(
      output_file = "trait_distributions_four_panel_by_cluster.pdf",
      raw_input_csv = "raw_input_trait_distributions_four_panel.csv",
      plot_family = "trait_distributions_four_panel",
      stringsAsFactors = FALSE
    ),
    data.frame(
      output_file = "ternary_centroid_summary_cluster_colors_variance_rank.pdf",
      raw_input_csv = "raw_input_ternary_summary_variance_rank.csv",
      plot_family = "ternary_centroid_summary_variance_rank",
      stringsAsFactors = FALSE
    ),
    data.frame(
      output_file = "ternary_connected_cluster_means_variance_rank.pdf",
      raw_input_csv = "raw_input_connected_cluster_means_variance_rank.csv",
      plot_family = "connected_cluster_means_variance_rank",
      stringsAsFactors = FALSE
    )
  )

  for (tg in c("OXYGEN", "STAIN", "MOTILITY", "SPORULATION")) {
    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      output_file = paste0("ternary_centroid_summary_cluster_colors_size_", tg, ".pdf"),
      raw_input_csv = paste0("raw_input_ternary_summary_size_", tg, ".csv"),
      plot_family = paste0("ternary_centroid_summary_size_", tg),
      stringsAsFactors = FALSE
    )
    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      output_file = paste0("ternary_connected_cluster_means_size_", tg, ".pdf"),
      raw_input_csv = paste0("raw_input_connected_cluster_means_size_", tg, ".csv"),
      plot_family = paste0("connected_cluster_means_size_", tg),
      stringsAsFactors = FALSE
    )
  }

  write.csv(
    dplyr::bind_rows(manifest_rows),
    file.path(viz_out_path, "retained_figure_input_manifest.csv"),
    row.names = FALSE
  )

  # Save tables
  write.csv(
    between_var_tbl,
    file.path(viz_out_path, "between_cluster_variance_rank_table.csv"),
    row.names = FALSE
  )
  write.csv(
    cluster_means_trait,
    file.path(viz_out_path, "cluster_means_by_trait.csv"),
    row.names = FALSE
  )

  cat("Saved ecology outputs to:", viz_out_path, "\n")
  invisible(TRUE)
}

run_nativity <- function(repo_root = base_path, viz_root = viz_out_root) {
  viz_out_path <- file.path(viz_root, "nativity")
  dir.create(viz_out_path, recursive = TRUE, showWarnings = FALSE)

  intermediate_files_path <- file.path(repo_root, "results/intermediate")

  save_pdf_figure <- function(file_name, width_mm = 180, height_mm = 120) {
    pdf(
      file = file.path(viz_out_path, file_name),
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      family = "Helvetica",
      colormodel = "rgb",
      pointsize = 5,
      useDingbats = FALSE
    )
  }

  ps_path <- file.path(intermediate_files_path, "ubiome_relative_none_updated.rds")
  if (!file.exists(ps_path)) stop("Missing required input phyloseq RDS: ", ps_path)
  ps <- readRDS(ps_path)

  sample_df <- as(sample_data(ps), "data.frame")
  if (!"Age" %in% colnames(sample_df)) stop("Missing `Age` in phyloseq sample_data.")
  if (!"BornInUSA" %in% colnames(sample_df)) stop("Missing `BornInUSA` in phyloseq sample_data.")

  age_breaks <- c(14, seq(20, 85, by = 5))
  age_labels <- c("14-19", paste(seq(20, 75, by = 5), seq(24, 79, by = 5), sep = "-"), "80-85")

  sample_df <- sample_df %>%
    mutate(
      age_group = cut(Age, breaks = age_breaks, labels = age_labels, right = FALSE),
      BornInUSA = as.character(BornInUSA),
      nativity = case_when(
        BornInUSA %in% c("US Born", "US_Born", "US-born", "US born") ~ "US Born",
        BornInUSA %in% c("non-US Born", "Non-US Born", "non US Born", "non-US-born", "non-US born") ~ "non-US Born",
        TRUE ~ BornInUSA
      )
    )

  nativity_levels <- c("US Born", "non-US Born")
  sample_df$nativity <- factor(sample_df$nativity, levels = nativity_levels)

  analysis_df <- sample_df %>%
    filter(!is.na(age_group)) %>%
    mutate(
      nativity_missing = is.na(nativity) | nativity == "",
      nativity_known = !nativity_missing
    )

  tab <- table(analysis_df$age_group, analysis_df$nativity, useNA = "ifany")
  cat("Nativity x age_group table:\n")
  print(tab)

  tab_known <- table(
    analysis_df$age_group[analysis_df$nativity_known],
    analysis_df$nativity[analysis_df$nativity_known]
  )
  tab_known <- tab_known[rowSums(tab_known) > 0, , drop = FALSE]
  chisq_global <- suppressWarnings(chisq.test(tab_known))

  analysis_df <- analysis_df %>% mutate(is_20_24 = age_group == "20-24")
  tab_20 <- table(
    analysis_df$is_20_24[analysis_df$nativity_known],
    analysis_df$nativity[analysis_df$nativity_known]
  )
  fisher_20 <- fisher.test(tab_20)

  df_known <- analysis_df %>% filter(nativity_known)
  p_20 <- mean(df_known$nativity[df_known$is_20_24] == "non-US Born", na.rm = TRUE)
  p_rest <- mean(df_known$nativity[!df_known$is_20_24] == "non-US Born", na.rm = TRUE)
  n_20 <- sum(df_known$is_20_24)
  n_rest <- sum(!df_known$is_20_24)
  diff_prop <- p_20 - p_rest
  se_diff <- sqrt((p_20 * (1 - p_20)) / n_20 + (p_rest * (1 - p_rest)) / n_rest)
  ci_low <- diff_prop - 1.96 * se_diff
  ci_high <- diff_prop + 1.96 * se_diff

  summary_df <- analysis_df %>%
    group_by(age_group) %>%
    summarise(
      n_total = n(),
      n_known = sum(nativity_known),
      n_missing = sum(nativity_missing),
      n_non_us = sum(nativity == "non-US Born", na.rm = TRUE),
      prop_non_us = if_else(n_known > 0, n_non_us / n_known, NA_real_),
      .groups = "drop"
    ) %>%
    arrange(age_group)

  summary_df <- summary_df %>%
    mutate(
      bar_label = if_else(
        n_known > 0,
        paste0(n_non_us, " / ", n_known),
        ""
      )
    )

  overall_prop_non_us <- mean(df_known$nativity == "non-US Born", na.rm = TRUE)

  plot_input_df <- summary_df %>%
    mutate(overall_prop_non_us = overall_prop_non_us)
  write.csv(
    plot_input_df,
    file.path(viz_out_path, "nativity_by_age_group_with_counts_plot_input.csv"),
    row.names = FALSE
  )

  p <- ggplot(summary_df, aes(x = age_group, y = prop_non_us)) +
    geom_col(fill = "gray70", color = "gray30", linewidth = 0.2) +
    geom_text(
      aes(label = bar_label),
      vjust = -0.35,
      size = 2.2,
      lineheight = 0.9,
      color = "gray15"
    ) +
    geom_hline(
      yintercept = overall_prop_non_us,
      linewidth = 0.35,
      linetype = "dashed",
      color = "steelblue4"
    ) +
    geom_point(
      data = subset(summary_df, age_group == "20-24"),
      aes(x = age_group, y = prop_non_us),
      color = "red3",
      size = 2
    ) +
    egg::theme_article() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 5),
      plot.title = element_text(size = 6, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 4, hjust = 0.5, color = "gray30")
    ) +
    labs(
      title = "Nativity composition by age group (unweighted)",
      subtitle = paste0(
        "Bar labels: non-US Born count / known nativity (n); global chisq p=", signif(chisq_global$p.value, 3),
        "; 20-24 vs others Fisher p=", signif(fisher_20$p.value, 3),
        "; overall non-US Born proportion=", signif(overall_prop_non_us, 3),
        "; 20-24 vs others diff=", signif(diff_prop, 3),
        " (Wald 95% CI ", signif(ci_low, 3), ", ", signif(ci_high, 3), ")"
      ),
      x = "Age group",
      y = "Proportion non-US Born (among known)"
    ) +
    scale_y_continuous(limits = c(0, 1.12), expand = c(0, 0))

  save_pdf_figure("nativity_by_age_group_with_counts.pdf", width_mm = 100, height_mm = 62)
  print(p)
  dev.off()

  cat("Saved nativity outputs to:", viz_out_path, "\n")
  invisible(TRUE)
}

run_ecology()
run_nativity()
