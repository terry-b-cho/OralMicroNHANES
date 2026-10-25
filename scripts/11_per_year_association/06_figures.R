#!/usr/bin/env Rscript

# =============================================================================
# MODULE 11 — STEP 06: publication figures
# 06_figures.R                       (run AFTER 05_harmonize_and_compare.R)
# =============================================================================
# Conventions (user-mandated):
#   * ALL text exactly 6 pt ArialMT; every figure <= 7 x 7 in except *_full
#     digital-only supplements; PDFs are pure vector; strict "FDR<0.05" wording.
#   * Colour-blind-safe palettes: phylum colours keep the project's predefined
#     hex palette; every other colour is Okabe-Ito (M11_CBF / M11_VAL_COLS in
#     00_module11_config.R); diverging scales use blue<->white<->vermillion.
#   * "Validated" in a cycle = replication-FDR<0.05 (BH over the one-sided
#     directional replication p-values within the pooled FDR<0.05 family,
#     design-df t reference) AND the per-cycle beta carries the pooled direction.
#
# Fig 1b  2_oradWAS outcomes: F-vs-G betas of pooled-significant tests,
#         phylum-coloured, SQUARE per-outcome panels.
# Fig 2   Pooled-significant betas vs per-cycle betas; 10 identical SQUARE
#         panels, zero centred on both axes (pure vector PDF).
# Fig 3   VALIDATION shares over the EVALUABLE pooled-significant set (bars sum
#         to 100%; %s annotated in-bar; untestable/unstable counts annotated
#         separately, not stacked).
# Fig 5   PER-SCHEMA MULTI-PAGE CATALOGS of ALL pooled FDR<0.05 associations:
#         pooled/F/G beta +- 95% CI per association (error-bar overlap view).
# Fig 7   Direction-reproducibility bars (pooled direction reproduced in BOTH
#         cycles; single metric), all pooled-significant vs prevalence>10%,
#         exact binomial 95% CIs, 50% chance line visible.
# Fig 8   THE F-vs-G scatter (unit-free z), validation-class colouring with
#         per-panel class counts, OLS fit line + r + pooled-sig n.
# Fig 9   Cross-cycle beta-profile correlation heatmaps; ONLY variables with
#         >=1 pooled FDR<0.05 association (and >=10 common genera) enter;
#         exWAS capped at top 25 (full uncapped version saved as *_full.pdf);
#         inclusion/exclusion inventory in fig9_variable_inclusion.csv.
# Fig 10  Validation counts as 100%-stacked horizontal bars over the single
#         pooled FDR<0.05 denominator (validated / direction / nominal p<0.05).
#
# Retired (code in archive/retired_figure_code_20260913.R, files in
# OUT_ROOT/archive/superseded_figures/): fig1, fig4, fig6, fig3_detailed.
# =============================================================================

source("scripts/11_per_year_association/00_module11_config.R")  # run from the repo root

m11_register_fonts()
m11_assert_arial()

harmonized  <- readRDS(file.path(M11_CMP, "harmonized_master_table.rds"))
rep_rowlvl  <- readRDS(file.path(M11_CMP, "pooled_vs_F_vs_G", "pooled_replication_rowlevel.rds"))

th  <- M11_THRESH_PRIMARY
tag <- gsub("\\.", "", sprintf("%02d", th * 100))
TXT <- 6 / .pt   # ALL in-plot text annotations at 6 pt
# Shared FIXED panel edge for the square main scatter figures (figS16B and
# figS16C): both use egg::set_panel_size with this value, so their x- and
# y-axis lengths match EXACTLY. 24mm keeps figS16B (5 facet rows) within 7in.
SQ_PANEL_MM <- 24
# percentage label precise to 1-2 decimals; never a false 0%/100% and never a
# decimal on an exact 0/100
m11_pct_lab <- function(p) ifelse(p %in% c(0, 100), sprintf("%.0f%%", p),
                           ifelse(p > 99.94 | p < 0.05, sprintf("%.2f%%", p),
                                  sprintf("%.1f%%", p)))

# ---------------- Figure 1b: 2_oradWAS outcomes (pooled-significant) ---------
d1b <- harmonized %>%
  filter(analysis == "2_oradWAS", testable_F, testable_G, stable_F, stable_G,
         direct_FG_comparison_valid, .data[[paste0("sig_fdr", tag, "_pooled")]])
if (nrow(d1b) > 0) {
  # three short lines: one long line would overrun the shrunken square panels
  ann1b <- d1b %>% group_by(variable) %>%
    summarise(lab = { pr <- m11_pearson(estimate_F, estimate_G)
                      da <- m11_dir_agree(estimate_F, estimate_G)
                      sprintf("r = %.2f\nconcordant %.0f%%\nn = %d", pr$r, da$pct, pr$n) },
              .groups = "drop")
  p1b <- ggplot(d1b, aes(estimate_F, estimate_G, color = phylum_std)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.25, color = "grey40") +
    geom_point(size = 0.45, alpha = 0.65, stroke = 0) +
    scale_color_manual(values = phylum_colors, name = NULL,
                       guide = guide_legend(override.aes = list(size = 1.3, alpha = 1), ncol = 1)) +
    # two-line strips (outcome name / code): one line is wider than a 28mm panel
    facet_wrap(~variable, scales = "free", nrow = 2,
               labeller = labeller(variable = function(v)
                 sub("_(OHAROC[A-Z]+)$", "\n\\1", v))) +
    geom_label(data = ann1b, aes(x = -Inf, y = Inf, label = lab),
               hjust = -0.02, vjust = 1.02, size = TXT, family = PDF_FAMILY,
               inherit.aes = FALSE, fill = scales::alpha("white", 0.7),
               linewidth = 0, label.padding = grid::unit(0.08, "lines")) +
    theme_compact() +
    labs(x = "beta, cycle F (2009-10)", y = "beta, cycle G (2011-12)",
         caption = sprintf("2_oradWAS tests with pooled BH-FDR<%.2f, comparable in both cycles.", th)) +
    theme(legend.key.height = grid::unit(0.5, "lines"))
  # EXACTLY square panels at a fixed size, device sized from the gtable's
  # natural dimensions (no aspect.ratio padding, no dead space, no clipping)
  p1b_fixed <- egg::set_panel_size(p1b, width = grid::unit(28, "mm"),
                                   height = grid::unit(28, "mm"))
  m11_save_figure("fig1b_beta_F_vs_G_oradWAS_outcomes", p1b_fixed,
                  width_mm = grid::convertWidth(sum(p1b_fixed$widths), "mm",
                                                valueOnly = TRUE) + 3,
                  height_mm = grid::convertHeight(sum(p1b_fixed$heights), "mm",
                                                  valueOnly = TRUE) + 3)
}

# ---------------- Figure 2: pooled vs per-cycle — square, centred panels -----
d2 <- rep_rowlvl %>% filter(pooled_metric == "fdr", stable_pooled) %>%
  mutate(estimate_F = ifelse(stable_F, estimate_F, NA_real_),
         estimate_G = ifelse(stable_G, estimate_G, NA_real_)) %>%
  select(analysis, estimate_pooled, estimate_F, estimate_G, phylum_std) %>%
  pivot_longer(c(estimate_F, estimate_G), names_to = "cycle",
               values_to = "estimate_cycle") %>%
  mutate(cycle = ifelse(cycle == "estimate_F", "cycle F (2009-10)", "cycle G (2011-12)")) %>%
  filter(is.finite(estimate_cycle)) %>%
  mutate(panel_lab = factor(sprintf("%s\n%s", analysis, cycle)))
ann2 <- d2 %>% group_by(analysis, cycle, panel_lab) %>%
  # three short lines: one long line overran the narrow square panels
  summarise(lab = { pr <- m11_pearson(estimate_pooled, estimate_cycle)
                    da <- m11_dir_agree(estimate_pooled, estimate_cycle)
                    sprintf("r = %.2f\nconcordant %s\nn = %s",
                            pr$r, m11_pct_lab(da$pct), format(pr$n, big.mark = ",")) },
            .groups = "drop")
lims2 <- d2 %>% group_by(panel_lab) %>%
  summarise(Mx = max(abs(estimate_pooled)) * 1.05,
            My = max(abs(estimate_cycle)) * 1.05, .groups = "drop")
blank2 <- bind_rows(lims2 %>% transmute(panel_lab, estimate_pooled = -Mx, estimate_cycle = -My),
                    lims2 %>% transmute(panel_lab, estimate_pooled =  Mx, estimate_cycle =  My))
p2 <- ggplot(d2, aes(estimate_pooled, estimate_cycle, color = phylum_std)) +
  geom_blank(data = blank2, inherit.aes = FALSE,
             aes(estimate_pooled, estimate_cycle)) +
  geom_hline(yintercept = 0, linewidth = 0.15, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.15, color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "grey40") +
  geom_point(size = 0.7, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = phylum_colors, name = NULL,
                     guide = guide_legend(override.aes = list(size = 1.2, alpha = 1), ncol = 1)) +
  facet_wrap(~panel_lab, ncol = 2, scales = "free") +
  geom_text(data = ann2, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.2, size = TXT, family = PDF_FAMILY,
            inherit.aes = FALSE, color = "black") +
  theme_compact() +
  labs(x = sprintf("pooled beta (BH-FDR<%.2f associations only)", th),
       y = "independent per-cycle beta") +
  theme(legend.key.height = grid::unit(0.5, "lines"))
# EXACTLY square panels at the shared edge length (axes match figS16C)
p2_fixed <- egg::set_panel_size(p2, width = grid::unit(SQ_PANEL_MM, "mm"),
                                height = grid::unit(SQ_PANEL_MM, "mm"))
m11_save_figure("figS16B_pooled_vs_per_cycle_beta", p2_fixed,
                width_mm = grid::convertWidth(sum(p2_fixed$widths), "mm",
                                              valueOnly = TRUE) + 3,
                height_mm = grid::convertHeight(sum(p2_fixed$heights), "mm",
                                                valueOnly = TRUE) + 3)

# ---------------- Figure 3: validation shares over the EVALUABLE set ---------
# Bars sum to 100% over pooled-significant tests testable with stable fits in
# both cycles; untestable/unstable counts are annotated, not stacked.
CAT5_LEVELS <- c("validated in both cycles",
                 "validated in cycle F only",
                 "validated in cycle G only",
                 "same direction in both, not validated",
                 "direction inconsistent in >=1 cycle")
CAT5_COLS <- setNames(unname(M11_VAL_COLS[c("both", "F_only", "G_only", "samedir", "inconsistent")]),
                      CAT5_LEVELS)
d3_all <- rep_rowlvl %>%
  mutate(metric_lab = ifelse(pooled_metric == "fdr",
                             sprintf("pooled BH-FDR<%.2f", th),
                             sprintf("pooled Storey q<%.2f", th)),
         evaluable = testable_F & testable_G & stable_F & stable_G)
d3 <- d3_all %>% filter(evaluable) %>%
  mutate(cat5 = factor(case_when(
    repfdr_F & repfdr_G ~ CAT5_LEVELS[1],
    repfdr_F ~ CAT5_LEVELS[2],
    repfdr_G ~ CAT5_LEVELS[3],
    dir_concord_F & dir_concord_G ~ CAT5_LEVELS[4],
    TRUE ~ CAT5_LEVELS[5]), levels = CAT5_LEVELS)) %>%
  count(analysis, metric_lab, cat5) %>%
  group_by(analysis, metric_lab) %>%
  mutate(pct = 100 * n / sum(n)) %>% ungroup()
d3_ann <- d3_all %>% group_by(analysis, metric_lab) %>%
  summarise(n_eval = sum(evaluable),
            n_absent = sum(!testable_F | !testable_G),
            n_unstable = sum(testable_F & testable_G & !(stable_F & stable_G)),
            .groups = "drop") %>%
  mutate(top_lab = sprintf("n=%s", format(n_eval, big.mark = ",", trim = TRUE)),
         excl_lab = sprintf("%s / %s",
                            format(n_absent, big.mark = ",", trim = TRUE),
                            format(n_unstable, big.mark = ",", trim = TRUE)))
p3 <- ggplot(d3, aes(analysis, n, fill = cat5)) +
  geom_col(position = "fill", width = 0.72) +
  geom_text(data = d3 %>% filter(pct >= 2.5),
            aes(label = sprintf("%.0f%%", pct)),
            position = position_fill(vjust = 0.5),
            size = TXT, family = PDF_FAMILY, color = "white") +
  geom_text(data = d3_ann, aes(analysis, 1.03, label = top_lab),
            inherit.aes = FALSE, size = TXT, family = PDF_FAMILY, vjust = 0) +
  geom_text(data = d3_ann, aes(analysis, -0.035, label = excl_lab),
            inherit.aes = FALSE, size = TXT, family = PDF_FAMILY, vjust = 1,
            color = "grey35") +
  scale_fill_manual(values = CAT5_COLS, name = NULL, guide = guide_legend(ncol = 1)) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.12, 0.08))) +
  facet_wrap(~metric_lab) +
  theme_compact() +
  labs(x = NULL, y = "share of evaluable pooled-significant associations",
       caption = "Denominator (n above each bar): pooled-significant tests testable with stable fits in BOTH cycles.\nCounts under each bar: excluded as variable-absent / fit-unstable (not stacked).\nValidated: per-cycle replication-FDR<0.05 (BH within the pooled-significant family, design-df t reference) AND pooled direction.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.key.size = grid::unit(0.55, "lines"))
m11_save_figure("fig3_replication_outcomes", p3, width_mm = 160, height_mm = 82)

# ---------------- Figure 3-direction: DIRECTION-ONLY criterion version -------
# Companion to fig3 over the SAME evaluable denominator, but classifying each
# association by the cycle pattern of the pooled DIRECTION alone (no
# significance requirement). Cycle-pattern colours (both/F/G/neither) are
# shared with the validation figures; the criterion is stated in the title,
# filename and caption.
DIR4_LEVELS <- c("pooled direction in both cycles",
                 "pooled direction in cycle F only",
                 "pooled direction in cycle G only",
                 "pooled direction in neither cycle")
DIR4_COLS <- setNames(unname(M11_VAL_COLS[c("both", "F_only", "G_only", "inconsistent")]),
                      DIR4_LEVELS)
d3d <- d3_all %>% filter(evaluable) %>%
  mutate(cat4 = factor(case_when(
    dir_concord_F & dir_concord_G ~ DIR4_LEVELS[1],
    dir_concord_F ~ DIR4_LEVELS[2],
    dir_concord_G ~ DIR4_LEVELS[3],
    TRUE ~ DIR4_LEVELS[4]), levels = DIR4_LEVELS)) %>%
  count(analysis, metric_lab, cat4) %>%
  group_by(analysis, metric_lab) %>%
  mutate(pct = 100 * n / sum(n)) %>% ungroup()
p3d <- ggplot(d3d, aes(analysis, n, fill = cat4)) +
  geom_col(position = "fill", width = 0.72) +
  # one decimal near the boundaries so 99.8% never displays as "100%"
  geom_text(data = d3d %>% filter(pct >= 2.5),
            aes(label = ifelse(pct > 99 & pct < 100, sprintf("%.1f%%", pct),
                               sprintf("%.0f%%", pct))),
            position = position_fill(vjust = 0.5),
            size = TXT, family = PDF_FAMILY, color = "white") +
  geom_text(data = d3_ann, aes(analysis, 1.03, label = top_lab),
            inherit.aes = FALSE, size = TXT, family = PDF_FAMILY, vjust = 0) +
  geom_text(data = d3_ann, aes(analysis, -0.035, label = excl_lab),
            inherit.aes = FALSE, size = TXT, family = PDF_FAMILY, vjust = 1,
            color = "grey35") +
  scale_fill_manual(values = DIR4_COLS, name = NULL, guide = guide_legend(ncol = 1)) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.12, 0.08))) +
  facet_wrap(~metric_lab) +
  theme_compact() +
  labs(x = NULL, y = "share of evaluable pooled-significant associations",
       caption = "DIRECTION-ONLY criterion: each association is classed by which cycles' betas carry the pooled\nsign, with NO per-cycle significance requirement (validation-criterion version: fig3).\nDenominator (n above each bar): pooled-significant tests testable with stable fits in BOTH\ncycles. Counts under each bar: excluded as variable-absent / fit-unstable (not stacked).") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.key.size = grid::unit(0.55, "lines"))
m11_save_figure("fig3_direction_outcomes", p3d, width_mm = 160, height_mm = 82)

# ---------------- Figure 5: catalogs of ALL pooled FDR<0.05 associations -----
ROWS_PER_PAGE <- 44
cat5_src <- rep_rowlvl %>%
  filter(pooled_metric == "fdr") %>%
  mutate(assoc_label = sprintf("%s ~ %s", variable, genus_label),
         phylum = m11_phylum_std(phylum_std),
         estimate_F = ifelse(testable_F & stable_F, estimate_F, NA_real_),
         std.error_F = ifelse(testable_F & stable_F, std.error_F, NA_real_),
         estimate_G = ifelse(testable_G & stable_G, estimate_G, NA_real_),
         std.error_G = ifelse(testable_G & stable_G, std.error_G, NA_real_)) %>%
  arrange(analysis, variable, desc(estimate_pooled))
catalog_counts <- list()
for (an in unique(cat5_src$analysis)) {
  ga <- cat5_src %>% filter(analysis == an) %>%
    mutate(page = ceiling(row_number() / ROWS_PER_PAGE))
  n_pages <- max(ga$page)
  plots <- vector("list", n_pages)
  for (pg in seq_len(n_pages)) {
    dpg <- ga %>% filter(page == pg) %>%
      mutate(assoc_label = factor(make.unique(assoc_label), levels = rev(make.unique(assoc_label))))
    dl <- dpg %>%
      select(assoc_label, phylum, estimate_pooled, std.error_pooled,
             estimate_F, std.error_F, estimate_G, std.error_G) %>%
      pivot_longer(-c(assoc_label, phylum),
                   names_to = c(".value", "which"),
                   names_pattern = "(estimate|std\\.error)_(.*)") %>%
      filter(is.finite(estimate)) %>%
      mutate(which = factor(which, levels = c("pooled", "F", "G"),
                            labels = c("pooled F+G", "cycle F", "cycle G")),
             lo = estimate - 1.96 * std.error, hi = estimate + 1.96 * std.error)
    plots[[pg]] <- ggplot(dl, aes(estimate, assoc_label, color = phylum, shape = which)) +
      geom_vline(xintercept = 0, linewidth = 0.2, color = "grey60") +
      geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                    linewidth = 0.25, position = position_dodge(width = 0.62)) +
      geom_point(size = 0.65, position = position_dodge(width = 0.62)) +
      scale_color_manual(values = phylum_colors, name = NULL, drop = TRUE,
                         guide = guide_legend(override.aes = list(size = 1.2), ncol = 1)) +
      scale_shape_manual(values = c("pooled F+G" = 16, "cycle F" = 17, "cycle G" = 15),
                         name = NULL) +
      theme_compact() +
      labs(title = sprintf("%s — pooled FDR<0.05 associations, page %d of %d", an, pg, n_pages),
           x = "beta +- 95% CI", y = NULL,
           caption = "ALL pooled BH-FDR<0.05 associations, with the independent cycle-F and cycle-G estimates for error-bar comparison.\nMissing cycle points: variable absent in that cycle or fit excluded by the stability guard.") +
      theme(legend.key.height = grid::unit(0.5, "lines"))
  }
  m11_save_figure_multipage(sprintf("fig5_%s_replicated_forest", an), plots,
                            width_mm = 175, height_mm = 175)
  catalog_counts[[an]] <- tibble(analysis = an, n_associations = nrow(ga), n_pages = n_pages)
}
write_csv(bind_rows(catalog_counts), file.path(M11_FIGURES_ACCESSORY, "fig5_catalog_index.csv"))

# ---------------- Figure 7: direction reproducibility, bar design ------------
# ONE metric only — both per-cycle betas carry the POOLED direction (the
# validation criterion the text uses). The earlier second metric (F=G with each
# other) differed by only 8 of 32,725 rows and made every point/label look
# duplicated. Bars run 0-100% so the dashed 50% chance line is visible and the
# gap to it IS the message; labels are plain "reproduced of total" counts.
ts <- rep_rowlvl %>%
  filter(pooled_metric == "fdr", testable_F, testable_G, stable_F, stable_G)
conc_ci <- function(d, subset_lab) {
  d %>% group_by(analysis) %>%
    summarise(n = n(), k = sum(dir_concord_F & dir_concord_G), .groups = "drop") %>%
    rowwise() %>%
    mutate(bt = list(binom.test(k, n, 0.5)),
           pct = 100 * k / n, lo = 100 * bt$conf.int[1], hi = 100 * bt$conf.int[2]) %>%
    ungroup() %>% select(-bt) %>% mutate(subset = subset_lab)
}
ts10 <- ts %>% filter(!is.na(prevalence_pooled), prevalence_pooled > M11_PREV_THRESHOLD)
SUBSET_LABS <- c(sprintf("all pooled FDR<%.2f associations", th),
                 sprintf("pooled FDR<%.2f & genus prevalence>%.0f%%",
                         th, 100 * M11_PREV_THRESHOLD))
d7 <- bind_rows(conc_ci(ts, SUBSET_LABS[1]), conc_ci(ts10, SUBSET_LABS[2])) %>%
  mutate(subset = factor(subset, levels = SUBSET_LABS),
         analysis = factor(analysis, levels = rev(sort(unique(analysis)))))
p7 <- ggplot(d7, aes(pct, analysis, fill = subset)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.22,
                linewidth = 0.3, color = "grey25",
                position = position_dodge(width = 0.78)) +
  # one decimal so 99.8% never displays as a false "100%"
  geom_text(aes(x = hi + 1.5, label = sprintf("%.1f%% (%s of %s)", pct,
                                              format(k, big.mark = ",", trim = TRUE),
                                              format(n, big.mark = ",", trim = TRUE))),
            position = position_dodge(width = 0.78), hjust = 0,
            size = TXT, family = PDF_FAMILY, show.legend = FALSE) +
  geom_vline(xintercept = 50, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  annotate("text", x = 48.5, y = 5.45, label = "chance (50%)", hjust = 1, vjust = 1,
           size = TXT, family = PDF_FAMILY, color = "grey40") +
  scale_fill_manual(values = setNames(unname(M11_CBF[c("blue", "vermillion")]), SUBSET_LABS),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 142), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_compact() +
  labs(x = "pooled direction reproduced in BOTH independent cycles", y = NULL,
       caption = "Bars: percentage of pooled BH-FDR<0.05 associations (testable with stable fits in\nboth cycles) whose cycle-F and cycle-G betas BOTH reproduce the pooled effect\ndirection; exact binomial 95% CI; labels: reproduced of total. Dashed: 50% chance.") +
  theme(legend.position = "top",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
m11_save_figure("fig7_direction_concordance_ci", p7, width_mm = 104, height_mm = 64)

# ---------------- Figure 7-v2: direction reproduced, by cycle pattern --------
# The SAME evaluable denominator split into its four mutually exclusive cycle
# patterns (both / F only / G only / neither reproduce the pooled direction);
# the four counts partition n exactly. Colours echo the fig8/fig10 validation
# classes for cross-figure consistency.
CAT7V <- c("reproduced in BOTH cycles (F & G)",
           "reproduced in cycle F only",
           "reproduced in cycle G only",
           "reproduced in NEITHER cycle")
d7v <- ts %>% group_by(analysis) %>%
  summarise(n = n(),
            `reproduced in BOTH cycles (F & G)` = sum(dir_concord_F & dir_concord_G),
            `reproduced in cycle F only`        = sum(dir_concord_F & !dir_concord_G),
            `reproduced in cycle G only`        = sum(!dir_concord_F & dir_concord_G),
            `reproduced in NEITHER cycle`       = sum(!dir_concord_F & !dir_concord_G),
            .groups = "drop") %>%
  pivot_longer(all_of(CAT7V), names_to = "cat", values_to = "k") %>%
  rowwise() %>%
  mutate(bt = list(binom.test(k, n)),
         pct = 100 * k / n, lo = 100 * bt$conf.int[1], hi = 100 * bt$conf.int[2]) %>%
  ungroup() %>% select(-bt) %>%
  mutate(cat = factor(cat, levels = rev(CAT7V)),
         analysis = factor(analysis, levels = rev(sort(unique(analysis)))))
CAT7V_COLS <- setNames(unname(M11_VAL_COLS[c("both", "F_only", "G_only", "inconsistent")]),
                       CAT7V)
# labels of LONG bars (>50%) go INSIDE the bar in white (just left of the CI
# whisker), short bars keep outside labels — this removes the wide right
# margin the outside labels forced, allowing a much narrower figure
d7v <- d7v %>%
  mutate(lab = sprintf("%s (%s of %s)",
                       ifelse(k == 0, "0%",
                              ifelse(pct < 0.1, sprintf("%.2f%%", pct),
                                     sprintf("%.1f%%", pct))),
                       format(k, big.mark = ",", trim = TRUE),
                       format(n, big.mark = ",", trim = TRUE)),
         inside = pct > 50)
p7v <- ggplot(d7v, aes(pct, analysis, fill = cat)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2,
                linewidth = 0.25, color = "grey25",
                position = position_dodge(width = 0.85)) +
  # ONE layer with per-row x/hjust/colour: subsetting rows per layer breaks
  # position_dodge alignment (each layer re-dodges among its own rows only)
  geom_text(aes(x = ifelse(inside, lo - 1.5, hi + 1.5),
                hjust = ifelse(inside, 1, 0), label = lab, colour = inside,
                group = cat),
            position = position_dodge(width = 0.85),
            size = TXT, family = PDF_FAMILY, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                      guide = "none") +
  scale_fill_manual(values = CAT7V_COLS, name = NULL,
                    breaks = CAT7V, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(limits = c(0, 107), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  theme_compact() +
  labs(x = "share of pooled FDR<0.05 associations,\nby cycle pattern of direction reproduction",
       y = NULL,
       caption = "The four mutually exclusive patterns partition each analysis' pooled BH-FDR<0.05\nassociations testable with stable fits in both cycles (n in the labels; the four counts\nsum to n exactly). Error bars: exact binomial 95% CI; labels: count of total.") +
  theme(legend.position = "top",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
m11_save_figure("figS16A_left_direction_reproduced_by_cycle", p7v,
                width_mm = 82, height_mm = 90)

# ---------------- Figure 7-validated: VALIDATION criterion version -----------
# Same bar design and subsets as fig7, but the criterion is full VALIDATION in
# both cycles: per-cycle replication-FDR<0.05 AND the pooled direction. The
# 50% chance line is omitted (chance level applies to sign agreement, not to
# a significance criterion).
val_ci <- function(d, subset_lab) {
  d %>% group_by(analysis) %>%
    summarise(n = n(), k = sum(repfdr_F & repfdr_G), .groups = "drop") %>%
    rowwise() %>%
    mutate(bt = list(binom.test(k, n)),
           pct = 100 * k / n, lo = 100 * bt$conf.int[1], hi = 100 * bt$conf.int[2]) %>%
    ungroup() %>% select(-bt) %>% mutate(subset = subset_lab)
}
d7val <- bind_rows(val_ci(ts, SUBSET_LABS[1]), val_ci(ts10, SUBSET_LABS[2])) %>%
  mutate(subset = factor(subset, levels = SUBSET_LABS),
         analysis = factor(analysis, levels = rev(sort(unique(analysis)))))
p7val <- ggplot(d7val, aes(pct, analysis, fill = subset)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.22,
                linewidth = 0.3, color = "grey25",
                position = position_dodge(width = 0.78)) +
  geom_text(aes(x = hi + 1.5, label = sprintf("%.1f%% (%s of %s)", pct,
                                              format(k, big.mark = ",", trim = TRUE),
                                              format(n, big.mark = ",", trim = TRUE))),
            position = position_dodge(width = 0.78), hjust = 0,
            size = TXT, family = PDF_FAMILY, show.legend = FALSE) +
  scale_fill_manual(values = setNames(unname(M11_CBF[c("blue", "vermillion")]), SUBSET_LABS),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 142), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_compact() +
  labs(x = "VALIDATED (replication-FDR<0.05 + pooled direction) in BOTH independent cycles",
       y = NULL,
       caption = "VALIDATION criterion: bars give the percentage of pooled BH-FDR<0.05 associations (testable\nwith stable fits in both cycles) that are replication-FDR<0.05 AND carry the pooled direction\nin BOTH cycles; exact binomial 95% CI; labels: validated of total. Direction-only criterion\nversion of the same bars: fig7.") +
  theme(legend.position = "top",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
m11_save_figure("fig7_validated_both_ci", p7val, width_mm = 104, height_mm = 64)

# ---------------- Figure 7v2-validated: VALIDATION partition by cycle --------
# Same partition design as fig7v2, but by which cycles reach full validation
# (replication-FDR<0.05 + pooled direction). "Neither" here includes
# associations directionally concordant in both cycles without per-cycle
# significance — direction-only partition: fig7v2.
CAT7VV <- c("validated in BOTH cycles (F & G)",
            "validated in cycle F only",
            "validated in cycle G only",
            "validated in NEITHER cycle")
d7vv <- ts %>% group_by(analysis) %>%
  summarise(n = n(),
            `validated in BOTH cycles (F & G)` = sum(repfdr_F & repfdr_G),
            `validated in cycle F only`        = sum(repfdr_F & !repfdr_G),
            `validated in cycle G only`        = sum(!repfdr_F & repfdr_G),
            `validated in NEITHER cycle`       = sum(!repfdr_F & !repfdr_G),
            .groups = "drop") %>%
  pivot_longer(all_of(CAT7VV), names_to = "cat", values_to = "k") %>%
  rowwise() %>%
  mutate(bt = list(binom.test(k, n)),
         pct = 100 * k / n, lo = 100 * bt$conf.int[1], hi = 100 * bt$conf.int[2]) %>%
  ungroup() %>% select(-bt) %>%
  mutate(cat = factor(cat, levels = rev(CAT7VV)),
         analysis = factor(analysis, levels = rev(sort(unique(analysis)))))
CAT7VV_COLS <- setNames(unname(M11_VAL_COLS[c("both", "F_only", "G_only", "inconsistent")]),
                        CAT7VV)
# long-bar labels inside / short-bar labels outside, as in fig7v2
d7vv <- d7vv %>%
  mutate(lab = sprintf("%s (%s of %s)",
                       ifelse(k == 0, "0%",
                              ifelse(pct < 0.1, sprintf("%.2f%%", pct),
                                     sprintf("%.1f%%", pct))),
                       format(k, big.mark = ",", trim = TRUE),
                       format(n, big.mark = ",", trim = TRUE)),
         inside = pct > 50)
p7vv <- ggplot(d7vv, aes(pct, analysis, fill = cat)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2,
                linewidth = 0.25, color = "grey25",
                position = position_dodge(width = 0.85)) +
  geom_text(aes(x = ifelse(inside, lo - 1.5, hi + 1.5),
                hjust = ifelse(inside, 1, 0), label = lab, colour = inside,
                group = cat),
            position = position_dodge(width = 0.85),
            size = TXT, family = PDF_FAMILY, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                      guide = "none") +
  scale_fill_manual(values = CAT7VV_COLS, name = NULL,
                    breaks = CAT7VV, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(limits = c(0, 107), breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  theme_compact() +
  labs(x = "share of pooled FDR<0.05 associations, by cycle pattern of VALIDATION",
       y = NULL,
       caption = "VALIDATION criterion: the four mutually exclusive patterns of which cycles reach\nreplication-FDR<0.05 with the pooled direction partition each analysis' evaluable\npooled BH-FDR<0.05 associations (counts sum to n exactly). 'Neither' includes\nassociations directionally concordant in both cycles without per-cycle significance.\nError bars: exact binomial 95% CI. Direction-only partition: fig7v2.") +
  theme(legend.position = "top",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
m11_save_figure("figS16A_right_validated_by_cycle", p7vv, width_mm = 82, height_mm = 90)

# ---------------- Figure 8: THE F-vs-G scatter (validation classes) ----------
zlim <- 12
# short legend labels so the legend fits the empty facet cell at the reduced
# figure size; the full criterion (replication-FDR<0.05 + pooled direction)
# stays in the caption
CLS8_LEVELS <- c("validated in both cycles",
                 "validated in cycle F only",
                 "validated in cycle G only",
                 "same direction, not validated",
                 "direction inconsistent")
CLS8_COLS <- setNames(unname(M11_VAL_COLS[c("both", "F_only", "G_only", "samedir", "inconsistent")]),
                      CLS8_LEVELS)
d8 <- ts %>% mutate(
  zr_F = estimate_F / std.error_F,
  zr_G = estimate_G / std.error_G,
  z_F = pmax(pmin(zr_F, zlim), -zlim),
  z_G = pmax(pmin(zr_G, zlim), -zlim),
  cls = factor(case_when(
    repfdr_F & repfdr_G           ~ CLS8_LEVELS[1],
    repfdr_F                      ~ CLS8_LEVELS[2],
    repfdr_G                      ~ CLS8_LEVELS[3],
    dir_concord_F & dir_concord_G ~ CLS8_LEVELS[4],
    TRUE                          ~ CLS8_LEVELS[5]), levels = CLS8_LEVELS))
fit8 <- d8 %>% group_by(analysis) %>%
  summarise(slope = coef(lm(zr_G ~ zr_F))[2],
            intercept = coef(lm(zr_G ~ zr_F))[1],
            r = suppressWarnings(cor(zr_F, zr_G)),
            n = n(),
            n_both = sum(repfdr_F & repfdr_G),
            n_Fonly = sum(repfdr_F & !repfdr_G),
            n_Gonly = sum(!repfdr_F & repfdr_G), .groups = "drop")
ann8 <- fit8 %>%
  mutate(lab = sprintf("fit r = %.2f\npooled-sig n = %s\nvalidated both: %s\nF only: %s | G only: %s",
                       r, format(n, big.mark = ",", trim = TRUE),
                       format(n_both, big.mark = ",", trim = TRUE),
                       format(n_Fonly, big.mark = ",", trim = TRUE),
                       format(n_Gonly, big.mark = ",", trim = TRUE)))
shade <- expand.grid(xmin = c(0, -zlim), analysis = unique(d8$analysis)) %>%
  mutate(xmax = ifelse(xmin == 0, zlim, 0), ymin = xmin, ymax = xmax)
p8 <- ggplot(d8, aes(z_F, z_G)) +
  geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = M11_CBF[["vermillion"]], alpha = 0.05) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey55") +
  geom_hline(yintercept = c(-1.96, 1.96), linewidth = 0.15, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = c(-1.96, 1.96), linewidth = 0.15, linetype = "dotted", color = "grey55") +
  geom_point(aes(color = cls), size = 0.35, alpha = 0.6, stroke = 0) +
  geom_abline(data = fit8, aes(slope = slope, intercept = intercept),
              linewidth = 0.3, color = "grey20") +
  scale_color_manual(values = CLS8_COLS, name = NULL, drop = FALSE,
                     guide = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1)) +
  facet_wrap(~analysis, nrow = 2) +
  # geom_label with a translucent white backing keeps the annotation readable
  # where it crosses dense point clouds (3_exWAS)
  geom_label(data = ann8, aes(x = -zlim + 0.4, y = zlim - 0.4, label = lab),
             hjust = 0, vjust = 1, size = TXT, family = PDF_FAMILY,
             inherit.aes = FALSE, fill = scales::alpha("white", 0.75),
             linewidth = 0, label.padding = grid::unit(0.08, "lines")) +
  coord_fixed(xlim = c(-zlim, zlim), ylim = c(-zlim, zlim)) +
  theme_compact() +
  labs(x = "z statistic, cycle F (beta/SE)", y = "z statistic, cycle G (beta/SE)",
       caption = sprintf("Pooled BH-FDR<%.2f tests, testable+stable in both cycles. Validated = per-cycle\nreplication-FDR<0.05 (BH over the one-sided directional replication p-values within the\npooled FDR<0.05 family, design-df t reference) with the pooled direction. Solid line: OLS\nfit of z_G on z_F; r = Pearson correlation (unclipped z). Shaded: sign-concordant quadrants.\nDotted: |z|=1.96. Clipped at |z|=%d.", th, zlim)) +
  theme(legend.position = c(0.84, 0.26),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
# EXACTLY square panels, device sized from the gtable's natural dimensions
p8_fixed <- egg::set_panel_size(p8, width = grid::unit(36, "mm"),
                                height = grid::unit(36, "mm"))
m11_save_figure("fig8_zstat_F_vs_G_pooled_sig", p8_fixed,
                width_mm = grid::convertWidth(sum(p8_fixed$widths), "mm",
                                              valueOnly = TRUE) + 3,
                height_mm = grid::convertHeight(sum(p8_fixed$heights), "mm",
                                                valueOnly = TRUE) + 3)

# ---------------- Figure 8-v2 (figS16C): F-vs-G map, F=G concordance x FDR ---
# Cycle-F vs cycle-G betas (per-panel free scales, symmetric about 0). Classes
# combine DIRECTIONAL CONCORDANCE BETWEEN THE TWO CYCLES (sign(beta_F) =
# sign(beta_G)) with per-cycle replication-FDR<0.05, so colour membership and
# quadrant membership align exactly: the four concordant classes live only in
# the shaded quadrants, the discordant class only outside them. (The z map in
# fig8 keeps the pooled-direction validation classes.)
# two-line legend labels: single lines are wider than the empty facet cell
CLS16C <- c("F=G concordant &\nFDR<0.05 in both cycles",
            "F=G concordant &\nFDR<0.05 in cycle F only",
            "F=G concordant &\nFDR<0.05 in cycle G only",
            "F=G concordant,\nFDR<0.05 in neither",
            "F-G direction discordant")
CLS16C_COLS <- setNames(unname(M11_VAL_COLS[c("both", "F_only", "G_only", "samedir", "inconsistent")]),
                        CLS16C)
d8c <- d8 %>% mutate(
  fg_conc = sign(estimate_F) == sign(estimate_G),
  ccls = factor(case_when(
    fg_conc & repfdr_F & repfdr_G   ~ CLS16C[1],
    fg_conc & repfdr_F              ~ CLS16C[2],
    fg_conc & repfdr_G              ~ CLS16C[3],
    fg_conc                         ~ CLS16C[4],
    TRUE                            ~ CLS16C[5]), levels = CLS16C))
lims8v <- d8 %>% group_by(analysis) %>%
  summarise(M = max(abs(c(estimate_F, estimate_G))) * 1.04, .groups = "drop")
blank8v <- bind_rows(lims8v %>% transmute(analysis, estimate_F = -M, estimate_G = -M),
                     lims8v %>% transmute(analysis, estimate_F =  M, estimate_G =  M))
fit8v <- d8c %>% group_by(analysis) %>%
  summarise(slope = coef(lm(estimate_G ~ estimate_F))[2],
            intercept = coef(lm(estimate_G ~ estimate_F))[1],
            r = suppressWarnings(cor(estimate_F, estimate_G)),
            n = n(),
            n_both  = sum(ccls == CLS16C[1]),
            n_Fonly = sum(ccls == CLS16C[2]),
            n_Gonly = sum(ccls == CLS16C[3]),
            n_nosig = sum(ccls == CLS16C[4]),
            n_disc  = sum(ccls == CLS16C[5]), .groups = "drop")
# annotations split across the corners: top-left r and n; bottom-right the
# class counts with 1-2-decimal percentages ("neither" class = remainder,
# defined in the caption); short lines fit the SQ_PANEL_MM square panels
# top-left (a DISCORDANT quadrant, sparse): r, n, and the discordant count;
# bottom-right: the four F=G-concordant classes together — keeps each label
# block <= 4 lines so neither covers the concordant data clusters
ann8v_tl <- fit8v %>%
  mutate(lab = sprintf("r = %.2f\nn = %s\ndiscord.: %s (%s)",
                       r, format(n, big.mark = ",", trim = TRUE),
                       format(n_disc, big.mark = ",", trim = TRUE), m11_pct_lab(100 * n_disc / n)))
ann8v_br <- fit8v %>%
  mutate(lab = sprintf("both: %s (%s)\nF: %s (%s)\nG: %s (%s)\nneither: %s (%s)",
                       format(n_both, big.mark = ",", trim = TRUE), m11_pct_lab(100 * n_both / n),
                       format(n_Fonly, big.mark = ",", trim = TRUE), m11_pct_lab(100 * n_Fonly / n),
                       format(n_Gonly, big.mark = ",", trim = TRUE), m11_pct_lab(100 * n_Gonly / n),
                       format(n_nosig, big.mark = ",", trim = TRUE), m11_pct_lab(100 * n_nosig / n)))
shade8v <- expand.grid(q = c("pos", "neg"), analysis = unique(d8$analysis)) %>%
  mutate(xmin = ifelse(q == "pos", 0, -Inf), xmax = ifelse(q == "pos", Inf, 0),
         ymin = xmin, ymax = xmax)
p8v <- ggplot(d8c, aes(estimate_F, estimate_G)) +
  geom_rect(data = shade8v, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = M11_CBF[["vermillion"]], alpha = 0.05) +
  geom_blank(data = blank8v, inherit.aes = FALSE, aes(estimate_F, estimate_G)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "grey40") +
  geom_point(aes(color = ccls), size = 0.7, alpha = 0.6, stroke = 0) +
  geom_abline(data = fit8v, aes(slope = slope, intercept = intercept),
              linewidth = 0.3, color = "grey20") +
  scale_color_manual(values = CLS16C_COLS, name = NULL, drop = FALSE,
                     guide = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1)) +
  facet_wrap(~analysis, nrow = 2, scales = "free") +
  geom_label(data = ann8v_tl, aes(x = -Inf, y = Inf, label = lab),
             hjust = -0.02, vjust = 1.02, size = TXT, family = PDF_FAMILY,
             inherit.aes = FALSE, fill = scales::alpha("white", 0.75),
             linewidth = 0, label.padding = grid::unit(0.08, "lines")) +
  geom_label(data = ann8v_br, aes(x = Inf, y = -Inf, label = lab),
             hjust = 1.02, vjust = -0.02, size = TXT, family = PDF_FAMILY,
             inherit.aes = FALSE, fill = scales::alpha("white", 0.75),
             linewidth = 0, label.padding = grid::unit(0.08, "lines")) +
  theme_compact() +
  labs(x = "beta, cycle F (2009-10)", y = "beta, cycle G (2011-12)",
       caption = sprintf("Pooled BH-FDR<%.2f associations, testable+stable in both cycles; axes\nare per-panel free scales, symmetric about zero. Dashed: identity.\nSolid: OLS fit of beta_G on beta_F; r = Pearson correlation of betas.\nClasses: F=G concordant = sign(beta_F)=sign(beta_G) (exactly the shaded\nquadrants); FDR<0.05 = per-cycle replication-FDR<0.05. The five classes\npartition n; annotations give all five as count (percent of n) —\ntop-left: discord. = F-G direction discordant; bottom-right:\nboth/F/G/neither = F=G concordant with FDR<0.05 in both cycles/\ncycle F only/cycle G only/neither cycle.", th)) +
  theme(legend.position = c(0.845, 0.235),
        legend.key.height = grid::unit(1.05, "lines"),
        legend.key.spacing.y = grid::unit(1.5, "pt"),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
# panel edge shared with figS16B so the two figures' axis lengths match exactly
p8v_fixed <- egg::set_panel_size(p8v, width = grid::unit(SQ_PANEL_MM, "mm"),
                                 height = grid::unit(SQ_PANEL_MM, "mm"))
m11_save_figure("figS16C_beta_F_vs_G_pooled_sig", p8v_fixed,
                width_mm = grid::convertWidth(sum(p8v_fixed$widths), "mm",
                                              valueOnly = TRUE) + 3,
                height_mm = grid::convertHeight(sum(p8v_fixed$heights), "mm",
                                                valueOnly = TRUE) + 3)

# ---------------- Figure 8v2-direction: DIRECTION-ONLY colouring -------------
# Identical beta scatter to fig8v2, but points classed by the cycle pattern of
# the pooled DIRECTION alone (no significance requirement); colours = the
# shared cycle-pattern palette. Validation-class colouring: fig8v2.
d8d <- d8 %>% mutate(dircls = factor(case_when(
  dir_concord_F & dir_concord_G ~ DIR4_LEVELS[1],
  dir_concord_F ~ DIR4_LEVELS[2],
  dir_concord_G ~ DIR4_LEVELS[3],
  TRUE ~ DIR4_LEVELS[4]), levels = DIR4_LEVELS))
ann8d <- d8d %>% group_by(analysis) %>%
  summarise(r = suppressWarnings(cor(estimate_F, estimate_G)), n = n(),
            n_both = sum(dircls == DIR4_LEVELS[1]),
            n_F = sum(dircls == DIR4_LEVELS[2]),
            n_G = sum(dircls == DIR4_LEVELS[3]), .groups = "drop") %>%
  mutate(lab = sprintf("beta r = %.2f\npooled-sig n = %s\ndirection both: %s\nF only: %s | G only: %s",
                       r, format(n, big.mark = ",", trim = TRUE),
                       format(n_both, big.mark = ",", trim = TRUE),
                       format(n_F, big.mark = ",", trim = TRUE),
                       format(n_G, big.mark = ",", trim = TRUE)))
p8vd <- ggplot(d8d, aes(estimate_F, estimate_G)) +
  geom_rect(data = shade8v, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = M11_CBF[["vermillion"]], alpha = 0.05) +
  geom_blank(data = blank8v, inherit.aes = FALSE, aes(estimate_F, estimate_G)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey55") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.25, color = "grey40") +
  geom_point(aes(color = dircls), size = 0.35, alpha = 0.6, stroke = 0) +
  geom_abline(data = fit8v, aes(slope = slope, intercept = intercept),
              linewidth = 0.3, color = "grey20") +
  scale_color_manual(values = DIR4_COLS, name = NULL, drop = FALSE,
                     guide = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1)) +
  facet_wrap(~analysis, nrow = 2, scales = "free") +
  geom_label(data = ann8d, aes(x = -Inf, y = Inf, label = lab),
             hjust = -0.02, vjust = 1.02, size = TXT, family = PDF_FAMILY,
             inherit.aes = FALSE, fill = scales::alpha("white", 0.75),
             linewidth = 0, label.padding = grid::unit(0.08, "lines")) +
  theme_compact() +
  labs(x = "beta, cycle F (2009-10)", y = "beta, cycle G (2011-12)",
       caption = sprintf("DIRECTION-ONLY criterion: the same pooled BH-FDR<%.2f associations as fig8v2 (testable+stable\nin both cycles), coloured by which cycles' betas carry the pooled sign — no per-cycle\nsignificance requirement (validation-class colouring: fig8v2). Axes: per-panel free scales,\nsymmetric about zero. Dashed: identity. Solid: OLS fit of beta_G on beta_F.", th)) +
  theme(legend.position = c(0.84, 0.26),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
p8vd_fixed <- egg::set_panel_size(p8vd, width = grid::unit(36, "mm"),
                                  height = grid::unit(36, "mm"))
m11_save_figure("fig8v2_beta_F_vs_G_direction", p8vd_fixed,
                width_mm = grid::convertWidth(sum(p8vd_fixed$widths), "mm",
                                              valueOnly = TRUE) + 3,
                height_mm = grid::convertHeight(sum(p8vd_fixed$heights), "mm",
                                                valueOnly = TRUE) + 3)

# ---------------- Figure 9: cross-cycle beta-profile correlation heatmaps ----
# INCLUSION RULE (defensible, stated in captions and inventoried in
# fig9_variable_inclusion.csv): a variable enters iff it has >=1 pooled
# BH-FDR<0.05 association AND >=10 common genera (pooled prevalence>10%,
# comparable+stable in both cycles). Only 3_exWAS exceeds 25 qualifying
# variables and is capped at the 25 with most pooled-significant hits for the
# main figure; the FULL uncapped exWAS matrix is saved as *_full.pdf
# (digital-only supplement, exceeds 7in by necessity).
HM_CELL_MM  <- 4.6
HM_VAR_CAP  <- 25
HM_MARGIN_W <- 62
HM_MARGIN_H <- 50
hm_base <- harmonized %>%
  filter(testable_F, testable_G, stable_F, stable_G, direct_FG_comparison_valid,
         !is.na(prevalence_pooled), prevalence_pooled > M11_PREV_THRESHOLD)
build_C <- function(da, vars) {
  wide <- da %>% filter(variable %in% vars) %>%
    select(otu_variable, variable, estimate_pooled, estimate_F, estimate_G)
  ids <- sort(unique(wide$otu_variable))
  mk <- function(col) {
    m <- wide %>% select(otu_variable, variable, all_of(col)) %>%
      pivot_wider(names_from = variable, values_from = all_of(col)) %>%
      as.data.frame()
    rownames(m) <- m$otu_variable; m$otu_variable <- NULL
    as.matrix(m[ids, vars, drop = FALSE])
  }
  MP <- mk("estimate_pooled"); MF <- mk("estimate_F"); MG <- mk("estimate_G")
  C <- matrix(NA_real_, length(vars), length(vars), dimnames = list(vars, vars))
  for (i in seq_along(vars)) for (j in seq_along(vars)) {
    ok <- is.finite(MF[, i]) & is.finite(MG[, j])
    if (sum(ok) >= 10) C[i, j] <- suppressWarnings(cor(MF[ok, i], MG[ok, j]))
  }
  CP <- suppressWarnings(cor(MP, use = "pairwise.complete.obs"))
  CP[!is.finite(CP)] <- 0
  ord <- if (length(vars) > 2) vars[hclust(as.dist(1 - CP))$order] else vars
  list(C = C, ord = ord)
}
plot_hm <- function(C, ord, an, annotate_diag, extra_caption = "",
                    caption_wrap = 96, legend_bar_mm = 30) {
  dl <- as.data.frame(as.table(C)) %>%
    rename(var_F = Var1, var_G = Var2, r = Freq) %>%
    mutate(var_F = factor(var_F, levels = rev(ord)),
           var_G = factor(var_G, levels = ord),
           on_diag = as.character(var_F) == as.character(var_G))
  diag_r <- dl %>% filter(on_diag) %>% pull(r)
  p <- ggplot(dl, aes(var_G, var_F, fill = r)) +
    geom_tile() +
    geom_tile(data = dl %>% filter(on_diag), fill = NA, color = "black",
              linewidth = ifelse(annotate_diag, 0.35, 0.2)) +
    { if (annotate_diag) geom_text(data = dl %>% filter(on_diag),
        aes(label = ifelse(is.na(r), "", sprintf("%.2f", r))),
        size = TXT, family = PDF_FAMILY) } +
    scale_fill_gradient2(low = M11_DIVERGING[["low"]], mid = M11_DIVERGING[["mid"]],
                         high = M11_DIVERGING[["high"]],
                         limits = c(-1, 1), na.value = "grey88",
                         name = "Pearson r\n(beta profiles)",
                         # display="rectangles": legend gradient drawn as vector
                         # rects, not an embedded raster (pure-vector PDF proof);
                         # barheight capped so the legend never outgrows a small panel
                         guide = guide_colorbar(display = "rectangles", nbin = 100,
                                                barheight = grid::unit(legend_bar_mm, "mm"))) +
    scale_x_discrete(expand = c(0, 0), position = "top") +
    scale_y_discrete(expand = c(0, 0)) +
    theme_compact() +
    # axis titles as a subtitle: the vertical y-title / long x-title grobs are
    # taller/wider than a small panel row and collide with caption and legend
    labs(title = sprintf("%s: cross-cycle beta-profile correlation", an),
         subtitle = "rows: cycle F (2009-10) betas  |  columns: cycle G (2011-12) betas",
         x = NULL, y = NULL,
         # dynamic wrap so the inserted cap-note never pushes a line past the
         # panel edge (fixed \n breaks clipped the exWAS caption)
         caption = paste(strwrap(sprintf(
           "Variables shown: >=1 pooled BH-FDR<0.05 association and >=10 common genera (prevalence>10%%, comparable in both cycles).%s Cell (i,j): Pearson r of variable i's cycle-F betas vs variable j's cycle-G betas across those genera (>=10 pairwise). Median diagonal r = %.2f. Order: pooled-profile clustering.",
           extra_caption, median(diag_r, na.rm = TRUE)), width = caption_wrap),
           collapse = "\n")) +
    # title/caption anchored to the FULL figure ("plot"), not the panel —
    # with egg::set_panel_size, panel-anchored text clips on small heatmaps
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          legend.key.height = grid::unit(0.9, "lines"),
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0))
  list(p = p, diag_r = diag_r)
}
diag_rows <- list(); incl_rows <- list()
for (an in M11_ANALYSES) {
  da <- hm_base %>% filter(analysis == an) %>%
    group_by(variable, otu_variable) %>%
    slice_min(p.value_pooled, n = 1, with_ties = FALSE) %>%
    ungroup()
  if (nrow(da) == 0) next
  vstats <- da %>% group_by(variable) %>%
    summarise(n_taxa = n(),
              n_pooled_sig = sum(.data[[paste0("sig_fdr", tag, "_pooled")]]),
              .groups = "drop")
  qual <- vstats %>% filter(n_pooled_sig >= 1, n_taxa >= 10) %>%
    arrange(desc(n_pooled_sig), desc(n_taxa))
  capped <- nrow(qual) > HM_VAR_CAP
  shown <- if (capped) qual %>% slice_head(n = HM_VAR_CAP) else qual
  incl_rows[[an]] <- vstats %>%
    mutate(analysis = an,
           included_main = variable %in% shown$variable,
           reason = case_when(
             variable %in% shown$variable ~ "included",
             n_pooled_sig < 1 ~ "no pooled BH-FDR<0.05 association",
             n_taxa < 10 ~ "fewer than 10 common genera",
             TRUE ~ sprintf("qualifying but below top-%d cap (exWAS main figure only)", HM_VAR_CAP))) %>%
    select(analysis, variable, n_taxa, n_pooled_sig, included_main, reason)
  vars <- shown$variable
  if (length(vars) < 2) next
  bc <- build_C(da, vars)
  cap_note <- if (capped)
    sprintf(" Showing the %d qualifying variables with most pooled-significant hits (of %d; full matrix: *_full.pdf).",
            HM_VAR_CAP, nrow(qual)) else ""
  # estimate the FINAL figure width from panel + actual row-label length +
  # legend, so the caption wrap never exceeds it (HM_MARGIN_W overestimates
  # for schemas with short labels, which clipped the pheWAS caption)
  fig_w <- length(vars) * HM_CELL_MM + max(nchar(vars)) * 1.15 + 24
  hm <- plot_hm(bc$C, bc$ord, an, annotate_diag = TRUE, extra_caption = cap_note,
                caption_wrap = max(46, floor(fig_w / 1.35)),
                legend_bar_mm = max(12, min(30, length(vars) * HM_CELL_MM - 10)))
  diag_rows[[an]] <- tibble(analysis = an, n_variables = length(vars),
                            n_qualifying = nrow(qual),
                            capped_at = ifelse(capped, HM_VAR_CAP, NA_integer_),
                            median_diagonal_r = median(hm$diag_r, na.rm = TRUE),
                            min_diagonal_r = min(hm$diag_r, na.rm = TRUE),
                            max_diagonal_r = max(hm$diag_r, na.rm = TRUE))
  panel_mm <- length(vars) * HM_CELL_MM
  p_fixed <- egg::set_panel_size(hm$p, width = grid::unit(panel_mm, "mm"),
                                 height = grid::unit(panel_mm, "mm"))
  # size the device from the gtable's NATURAL dimensions (fixed panel + label,
  # title, caption rows) — fixed margins clipped small heatmaps
  m11_save_figure(sprintf("fig9_%s_beta_correlation_heatmap", an), p_fixed,
                  width_mm = grid::convertWidth(sum(p_fixed$widths), "mm",
                                                valueOnly = TRUE) + 3,
                  height_mm = grid::convertHeight(sum(p_fixed$heights), "mm",
                                                  valueOnly = TRUE) + 3)
  # FULL uncapped version (digital-only supplement) for capped schemas
  if (capped) {
    bcf <- build_C(da, qual$variable)
    hmf <- plot_hm(bcf$C, bcf$ord, an, annotate_diag = FALSE,
                   extra_caption = sprintf(" FULL matrix: all %d qualifying variables (digital supplement).", nrow(qual)),
                   caption_wrap = 160)
    cell_full <- 2.6
    panel_full <- nrow(qual) * cell_full
    pf_fixed <- egg::set_panel_size(hmf$p, width = grid::unit(panel_full, "mm"),
                                    height = grid::unit(panel_full, "mm"))
    m11_save_figure(sprintf("fig9_%s_beta_correlation_heatmap_full", an), pf_fixed,
                    width_mm = grid::convertWidth(sum(pf_fixed$widths), "mm",
                                                  valueOnly = TRUE) + 3,
                    height_mm = grid::convertHeight(sum(pf_fixed$heights), "mm",
                                                    valueOnly = TRUE) + 3,
                    allow_oversize = TRUE)
  }
}
write_csv(bind_rows(diag_rows), file.path(M11_FIGURES_ACCESSORY, "fig9_diagonal_summary.csv"))
write_csv(bind_rows(incl_rows), file.path(M11_FIGURES_ACCESSORY, "fig9_variable_inclusion.csv"))

# ---------------- Figure 10: validation counts, 100%-stacked horizontal ------
Q_LEVELS <- c("validated: replication-FDR<0.05 + pooled direction",
              "direction concordant with pooled",
              "nominal p<0.05 (two-sided)")
C10_LEVELS <- c("both cycles", "cycle F only", "cycle G only", "neither cycle")
C10_COLS <- setNames(c(unname(M11_VAL_COLS[c("both", "F_only", "G_only")]), "grey80"),
                     C10_LEVELS)
mk_overlap <- function(d, sF, sG, q) {
  d %>% mutate(.f = sF, .g = sG) %>% group_by(analysis) %>%
    summarise(`both cycles` = sum(.f & .g), `cycle F only` = sum(.f & !.g),
              `cycle G only` = sum(!.f & .g), `neither cycle` = sum(!.f & !.g),
              denominator = n(), .groups = "drop") %>%
    pivot_longer(all_of(C10_LEVELS), names_to = "category", values_to = "n") %>%
    mutate(question = q)
}
d10 <- bind_rows(
  mk_overlap(ts, ts$repfdr_F, ts$repfdr_G, Q_LEVELS[1]),
  mk_overlap(ts, ts$dir_concord_F, ts$dir_concord_G, Q_LEVELS[2]),
  mk_overlap(ts, !is.na(ts$p.value_F) & ts$p.value_F < 0.05,
                 !is.na(ts$p.value_G) & ts$p.value_G < 0.05, Q_LEVELS[3])) %>%
  group_by(question, analysis) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  mutate(question = factor(question, levels = Q_LEVELS),
         category = factor(category, levels = C10_LEVELS),
         axis_lab = sprintf("%s (n=%s)", analysis,
                            format(denominator, big.mark = ",", trim = TRUE)))
d10 <- d10 %>% mutate(axis_lab = factor(axis_lab, levels = rev(sort(unique(axis_lab)))))
p10 <- ggplot(d10, aes(n, axis_lab, fill = category)) +
  geom_col(position = "fill", orientation = "y", width = 0.72) +
  geom_text(data = d10 %>% filter(pct >= 4),
            aes(label = format(n, big.mark = ",", trim = TRUE)),
            position = position_fill(vjust = 0.5), orientation = "y",
            size = TXT, family = PDF_FAMILY, color = "white") +
  scale_fill_manual(values = C10_COLS, name = NULL, drop = FALSE,
                    guide = guide_legend(nrow = 1)) +
  scale_x_continuous(labels = function(x) paste0(100 * x, "%"),
                     expand = c(0, 0)) +
  facet_wrap(~question, ncol = 1) +
  theme_compact() +
  labs(x = "share of pooled FDR<0.05 associations", y = NULL,
       caption = "Single denominator throughout (n in each row label): pooled BH-FDR<0.05 tests, testable with stable fits in both cycles.\nSegment labels: counts. Validated = per-cycle replication-FDR<0.05 (BH within the pooled-significant family) AND pooled direction.") +
  theme(legend.position = "top",
        panel.spacing.y = grid::unit(5, "pt"))
m11_save_figure("fig10_significance_overlap", p10, width_mm = 150, height_mm = 105)

cat("STEP 06 complete. Figures in", M11_FIGURES, "\n")