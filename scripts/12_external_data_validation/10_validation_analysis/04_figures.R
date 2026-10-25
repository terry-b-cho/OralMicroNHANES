#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 2 — Step 4
# Figures: the nine panels of Supplementary Figure S15
#
# Purpose : Render the validation figures in the project's house style so they
#           sit alongside modules 4/5/7/8/9 without re-explaining the colours.
#
#           Conventions enforced by 00_validation_config.R:
#             * Arial (ArialMT), asserted before any device opens
#             * base font size 6, no plot exceeding 7 inches on either axis
#             * the canonical 17-entry phylum palette used repo-wide
#             * PNG rendered FROM the PDF (systemfonts would substitute
#               DejaVuSans, so a re-drawn raster would not be Arial)
#
# Usage   : Rscript 04_figures.R
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
  library(ggplot2); library(ggrepel); library(patchwork)
})
HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/10_validation_analysis"
source(file.path(HERE, "00_validation_config.R"))
m12_assert_arial()
vsay("=== Step 4: figures ===")

cmp <- bind_rows(lapply(PRIMARY_VARIABLES, function(v)
  readr::read_tsv(file.path(VAL, paste0("comparison_", v, ".tsv")),
                  show_col_types = FALSE, progress = FALSE)))
sm  <- readr::read_tsv(file.path(VAL, "concordance_summary.tsv"),
                       show_col_types = FALSE, progress = FALSE)
nh  <- bind_rows(lapply(PRIMARY_VARIABLES, function(v)
  readr::read_tsv(file.path(VAL, paste0("nhanes_effects_", v, ".tsv")),
                  show_col_types = FALSE, progress = FALSE)))

# Labels come from the config so a new phenotype cannot be added to the
# pipeline and then render as NA in a figure.
VAR_LABEL <- VAR_LABELS
VAR_SHORT <- VAR_LABELS_SHORT
stopifnot(all(PRIMARY_VARIABLES %in% names(VAR_LABEL)),
          all(PRIMARY_VARIABLES %in% names(VAR_SHORT)))
cmp$phylum <- m12_phylum_std(cmp$phylum)
nh$phylum  <- m12_phylum_std(nh$phylum)

# Only show phyla actually present, so the legend stays small at 6 pt.
present_phyla <- sort(unique(c(cmp$phylum, nh$phylum)))
pal <- phylum_colors[present_phyla]

# =============================================================================
# Panel B - Effect scatter: NHANES vs external, refit-cohort primary fits
# =============================================================================
f1_data <- cmp %>% filter(tier == "A", primary, tested)

panel_scatter <- function(d) {
  s <- sm %>% filter(dataset_id == d$dataset_id[1], fit == d$fit[1],
                     nhanes_variable == d$nhanes_variable[1])
  # Calibration slope of external ~ NHANES. Conditioning on NHANES significance
  # inflates our estimates, so the winner's-curse prediction is slope < 1.
  # Report it with its SE rather than asserting attenuation qualitatively.
  fit <- stats::lm(ext_estimate ~ nh_estimate, data = d)
  cf  <- stats::coef(summary(fit))
  # Two short lines in the SUBTITLE, not an in-panel annotation: at this panel
  # size a floating stats box sits on top of real points.
  # Commas, not pipes, and no hard line break: m12_wrap() re-flows this to the
  # canvas width, and a "|" that lands at the start of a wrapped line reads as
  # a stray glyph. "Spearman" already implies rank, so it is not repeated.
  lab <- sprintf("Spearman rho = %s, calibration slope = %.2f (SE %.2f), n = %d genera, N = %s participants",
                 ifelse(is.na(s$spearman_rho[1]), "NA", sprintf("%.2f", s$spearman_rho[1])),
                 cf["nh_estimate", "Estimate"], cf["nh_estimate", "Std. Error"],
                 nrow(d), format(s$n_participants[1], big.mark = ","))
  top <- d %>% mutate(m = abs(nh_estimate)) %>% arrange(desc(m)) %>% head(6)
  ggplot(d, aes(nh_estimate, ext_estimate)) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, colour = "grey80", linewidth = 0.2,
                linetype = "dotted") +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey35", linewidth = 0.3) +
    geom_point(aes(fill = phylum, size = sig_bin_from_dot_size(calc_dot_size(nh_fdr))),
               shape = 21, stroke = 0.15, alpha = 0.85) +
    geom_text_repel(data = top, aes(label = genus), size = GEOM_TEXT_SIZE, family = PDF_FAMILY,
                    box.padding = 0.25, point.padding = 0, segment.size = 0.12,
                    segment.colour = "grey60", min.segment.length = 0.05,
                    max.overlaps = 20) +
    scale_fill_manual(values = pal, drop = FALSE, name = "Phylum") +
    scale_size_manual(values = SIG_SIZE_VALUES, drop = FALSE, name = "NHANES FDR") +
    guides(fill = guide_legend(override.aes = list(size = 1.5), ncol = 1, order = 1),
           size = guide_legend(override.aes = list(fill = "grey55"), ncol = 1, order = 2)) +
    labs(title = m12_wrap(d$label[1], 40, bold = TRUE),
         subtitle = m12_wrap(paste0("1 point = 1 genus tested in both. ", lab), 40),
         x = "NHANES effect (CLR per unit exposure)",
         y = "External refitted effect (CLR per unit exposure)") +
    theme_compact()
}

f1_panels <- f1_data %>% group_split(nhanes_variable, dataset_id, fit) %>% lapply(panel_scatter)
n_p <- length(f1_panels)
# Both axes are the same quantity in the same units (CLR per unit exposure), so
# the panel MUST be square: on a panel wider than it is tall the dotted identity
# line is not the visual 45-degree diagonal, and a reader judging agreement
# against it is misled. There is one refit comparison, so the single panel is
# pinned square directly; the patchwork path exists only for the day a second
# refit cohort appears, and cannot pin panel geometry.
if (n_p == 1) {
  m12_save_square_panel(S15_FIG[["B"]],
                        f1_panels[[1]] + ggplot2::theme(legend.position = "right"),
                        side_mm = 48)
} else {
  n_col <- 2L; n_row <- ceiling(n_p / n_col)
  f1 <- patchwork::wrap_plots(f1_panels, ncol = n_col) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "right")
  save_validation_figure(S15_FIG[["B"]], f1,
                         width_mm  = min(170, 48 * n_col + 32),
                         height_mm = min(175, 48 * n_row + 18))
}

# =============================================================================
# Panels D-1 / D-2 / D-3 - Concordance grid: every NHANES hit x every cohort
# =============================================================================
# Row order is needed both to build the plot and to budget the label gutter,
# so it is computed once, outside, and passed in.
row_order <- function(v)
  nh %>% filter(variable == v) %>% arrange(phylum, estimate) %>% pull(genus)

grid_one <- function(v, ord) {
  d <- cmp %>% filter(nhanes_variable == v, primary)
  # Levels must come from the PLOTTED rows. Taking them from `sm` pulled in the
  # non-primary sensitivity fit, adding a phantom 4th factor level that the
  # panel was then sized for - which silently made every drawn tile wider here
  # than in the single-column grids.
  d <- d %>% mutate(genus = factor(genus, levels = ord),
                    outcome = factor(outcome, levels = OUTCOME_LEVELS),
                    label = factor(label, levels = sort(unique(label))))
  ycol <- nh %>% filter(variable == v) %>% arrange(phylum, estimate) %>%
    mutate(col = phylum_colors[phylum]) %>% pull(col)
  n_hits <- sum(nh$variable == v)
  n_test <- sum(d$outcome != "not tested")
  n_coh  <- dplyr::n_distinct(d$label)
  ggplot(d, aes(label, genus, fill = outcome)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = outcome_colors, drop = FALSE, name = NULL) +
    # State the unit and both denominators on the figure. A reader should not
    # have to infer whether a tile is a genus, a participant or a comparison.
    labs(title = sprintf("%s: direction of each NHANES genus vs each external cohort",
                         VAR_LABEL[v]),
         # Name the axes. "N pairs" alone reads as the row count, which it is
         # only for the single-cohort grids; for age the 162 pairs are 54 rows
         # x 3 columns, and that arithmetic has to be visible on the figure.
         subtitle = sprintf(
           "Rows = %d NHANES genera (FDR<0.05, prev>=10%%); columns = %d cohort%s; %d tiles, %d testable",
           n_hits, n_coh, ifelse(n_coh > 1, "s", ""), nrow(d), n_test),
         x = NULL, y = NULL) +
    theme_compact() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(colour = ycol),
          legend.position = "bottom", panel.grid = element_blank())
}
for (v in PRIMARY_VARIABLES) {
  # Cell geometry is FIXED across every heat map in the module, so one tile in
  # the single-column smoking grid is exactly as wide as one tile in the
  # three-column age grid. The device is then sized to the result, which also
  # removes the dead width the old hand-budgeted formula left behind.
  ord    <- row_order(v)
  n_rows <- sum(nh$variable == v)
  n_cols <- dplyr::n_distinct(cmp$label[cmp$nhanes_variable == v & cmp$primary])
  m12_save_heatmap(S15_FIG[[S15_GRID_PANEL[[v]]]], grid_one(v, ord),
                   n_cols = n_cols, n_rows = n_rows)
}

# =============================================================================
# Panels D2-i / D2-ii / D2-iii - Panel D version 2
#
# The same tiles as D-1/D-2/D-3, graded by STRENGTH of external evidence rather
# than by direction alone. D version 1 answers "did the direction agree?"; this
# answers "how well did it agree?", separating a genus that reproduced at
# external FDR from one that merely pointed the same way.
#
# Cell geometry is identical to version 1 (6.50 x 2.35 mm), so a v1 and a v2
# grid can be laid side by side and read tile-for-tile.
#
# The classification is recomputed from the columns Phase 2 already holds, NOT
# read from the Phase 3 audit output: the audit runs after this script, so
# depending on it would invert the pipeline order. m12_classify_evidence() is
# verified to reproduce master_validation_table.tsv exactly.
# =============================================================================
grid_one_v2 <- function(v, ord) {
  d <- cmp %>% filter(nhanes_variable == v, primary) %>%
    mutate(cls = m12_classify_evidence(tested, outcome, ext_p, ext_fdr),
           cls = factor(cls, levels = PHASE2_CLASSES),
           genus = factor(genus, levels = ord),
           label = factor(label, levels = sort(unique(label))))
  ycol <- nh %>% filter(variable == v) %>% arrange(phylum, estimate) %>%
    mutate(col = phylum_colors[phylum]) %>% pull(col)
  n_hits <- sum(nh$variable == v)
  n_coh  <- dplyr::n_distinct(d$label)
  n_fdr  <- sum(d$cls == "FDR_REPLICATED_SAME_DIRECTION")
  # ggplot2 4.0 paints NO legend key for a level absent from the data, even with
  # drop = FALSE, an explicit `limits`, or override.aes - all three were tried
  # and all leave the label beside an empty box. A zero-area layer carrying
  # every level does work, and draws nothing on the panel. This matters because
  # "nominal p<0.05" occurs for none of these three phenotypes, and the reader
  # must see that the grade was assessed and never reached, not omitted.
  legend_seed <- tibble::tibble(label = d$label[1], genus = d$genus[1],
                                cls = factor(PHASE2_CLASSES, levels = PHASE2_CLASSES))
  ggplot(d, aes(label, genus, fill = cls)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_tile(data = legend_seed, aes(fill = cls), width = 0, height = 0) +
    scale_fill_manual(values = validation_colors[PHASE2_CLASSES],
                      labels = VALIDATION_CLASS_LABELS_SHORT[PHASE2_CLASSES],
                      drop = FALSE, name = NULL) +
    labs(title = sprintf("%s: strength of external evidence per genus", VAR_LABEL[v]),
         subtitle = sprintf(
           "Rows = %d NHANES genera; columns = %d cohort%s; %d tiles, %d replicated at external FDR",
           n_hits, n_coh, ifelse(n_coh > 1, "s", ""), nrow(d), n_fdr),
         x = NULL, y = NULL) +
    theme_compact() +
    # Legend on the RIGHT, not the bottom. These grids are already near the
    # 7-inch height cap, and five legend rows underneath pushed the age grid to
    # 7.85 in; moving the legend sideways spends width, which is not scarce.
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          axis.text.y = element_text(colour = ycol),
          legend.position = "right", panel.grid = element_blank()) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE))
}
for (v in PRIMARY_VARIABLES) {
  ord    <- row_order(v)
  n_rows <- sum(nh$variable == v)
  n_cols <- dplyr::n_distinct(cmp$label[cmp$nhanes_variable == v & cmp$primary])
  m12_save_heatmap(S15_FIG[[S15_GRID_PANEL_V2[[v]]]], grid_one_v2(v, ord),
                   n_cols = n_cols, n_rows = n_rows)
}

# =============================================================================
# Panel A-2 - Concordance summary with 95% binomial CIs
# =============================================================================
f3d <- sm %>% filter(!is.na(pct_concordant)) %>%
  mutate(lab2 = sprintf("%s  (%d/%d)", label, n_concordant, n_tested),
         variable = VAR_SHORT[nhanes_variable]) %>%
  arrange(nhanes_variable, pct_concordant) %>%
  mutate(lab2 = factor(lab2, levels = unique(lab2)))
f3 <- ggplot(f3d, aes(pct_concordant, lab2, colour = tier)) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey55", linewidth = 0.25) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.3) +
  geom_point(aes(shape = primary), size = 1.2) +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(A = "#CC6677", B = "#0072B2"),
                      labels = COHORT_CLASS_SHORT, name = "Cohort class") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Primary") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(title = m12_wrap("Share of NHANES genera whose direction agrees, per cohort", 104, bold = TRUE),
       subtitle = m12_wrap("Unit = genus. (concordant / testable genera) shown after each cohort; dashed line = chance", 104),
       x = "Genera agreeing in direction with NHANES (%, exact 95% CI)", y = NULL) +
  theme_compact() +
  # A vertical strip cannot hold "Gum disease" at 6 pt on a one-row facet;
  # rotate to horizontal and buy the width instead of clipping the label.
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
save_validation_figure(S15_FIG[["A-2"]], f3, width_mm = 104, height_mm = 60)

# =============================================================================
# Panel A-1 - The two sample sizes: participants and comparable genera
# =============================================================================
f4d <- sm %>% filter(primary) %>%
  mutate(variable = VAR_SHORT[nhanes_variable],
         lab2 = paste0(label, " [", m12_cohort_class(tier, short = TRUE), "]")) %>%
  select(variable, lab2, tier, `Participants (N)` = n_participants,
         `Comparable genera (n)` = n_tested) %>%
  pivot_longer(c(`Participants (N)`, `Comparable genera (n)`),
               names_to = "measure", values_to = "value")
f4 <- ggplot(f4d, aes(value, lab2, fill = tier)) +
  geom_col(width = 0.65, alpha = 0.85, na.rm = TRUE) +
  geom_text(data = ~ dplyr::filter(.x, !is.na(value)),
            aes(label = format(value, big.mark = ",")),
            hjust = -0.1, size = GEOM_TEXT_SIZE, family = PDF_FAMILY, colour = "grey20") +
  # A published-target cohort holds only a published table, so there is no
  # participant count. Say so in the panel rather than leaving an empty row
  # that reads as missing data.
  geom_text(data = ~ dplyr::filter(.x, is.na(value)), x = 0,
            aes(label = "no participant-level data"),
            hjust = -0.03, size = GEOM_TEXT_SIZE, family = PDF_FAMILY, colour = "grey45") +
  facet_grid(variable ~ measure, scales = "free", space = "free_y") +
  scale_fill_manual(values = c(A = "#CC6677", B = "#0072B2"),
                    labels = COHORT_CLASS_SHORT, name = "Cohort class") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = m12_wrap("Two different denominators: participants contributing, and genera comparable", 112, bold = TRUE),
       subtitle = m12_wrap("Left column counts PEOPLE; right column counts GENERA. They are not interchangeable", 112),
       x = NULL, y = NULL) +
  theme_compact() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0))
save_validation_figure(S15_FIG[["A-1"]], f4, width_mm = 112, height_mm = 58)

# =============================================================================
# Panel E - Which NHANES hits are even testable externally
# =============================================================================
cover <- cmp %>% filter(primary, tested) %>%
  count(nhanes_variable, genus, name = "n_datasets_testing")
f5d <- nh %>%
  left_join(cover, by = c("variable" = "nhanes_variable", "genus")) %>%
  mutate(n_datasets_testing = tidyr::replace_na(n_datasets_testing, 0L),
         variable = VAR_LABEL[variable])
f5 <- ggplot(f5d, aes(100 * prevalence_overall, estimate)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
  geom_point(aes(fill = phylum, size = factor(n_datasets_testing)),
             shape = 21, stroke = 0.15, alpha = 0.85) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values = pal, name = "Phylum",
                    guide = guide_legend(ncol = 2, override.aes = list(size = 1.4))) +
  scale_size_manual(values = c(`0` = 0.6, `1` = 1.1, `2` = 1.6, `3` = 2.1, `4` = 2.6),
                    name = "Number of external cohorts testing this genus") +
  labs(title = m12_wrap("Which NHANES genera are testable externally, by prevalence and effect size", 96, bold = TRUE),
       subtitle = m12_wrap("1 point = 1 NHANES genus; point size = how many cohorts carry it, not a count of samples", 96),
       x = "NHANES prevalence (% of participants carrying the genus)",
       y = "NHANES effect (per-facet scale; see caption)") +
  theme_compact() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(1, "pt"))
save_validation_figure(S15_FIG[["E"]], f5, width_mm = 96, height_mm = 64)

vsay("")
vsay("figures written to ", VAL_FIG_PDF)
vsay("PNGs (rasterized from the PDFs) in ", VAL_FIG_PNG)
print(list.files(VAL_FIG_PDF, pattern = "\\.pdf$"))

# =============================================================================
# Panel A-3 - How much of the NHANES result is externally checkable at all
#
# Every concordance percentage in this report is conditional on a genus having
# been tested. This figure supplies that missing denominator, and checks
# whether the untested genera are systematically the weaker associations
# (which would make the tested subset an optimistic sample).
# =============================================================================
cov_g <- readr::read_tsv(file.path(VAL, "coverage_per_genus.tsv"),
                         show_col_types = FALSE, progress = FALSE)
cov_g$phylum <- m12_phylum_std(cov_g$phylum)
cov_g$variable <- VAR_LABEL[cov_g$nhanes_variable]

# Sub-plot (i) is built from the SAME source as Panel A-1 - the per-cohort tested counts
# in concordance_summary - so every bar here is a number that also appears in
# Panel A-1. The previous form counted DISTINCT GENERA by how many cohorts test them,
# which is a different quantity: for age it reported 41 covered genera where Panel A-1
# reports 35 + 10 + 4 = 49 genus-cohort comparisons, because 7 age genera are
# tested by more than one cohort. Both were correct; only one can be compared
# with Panel A-1, so sub-plot (i) now counts per cohort.
#
# Bars are DODGED, not stacked: a genus tested by two cohorts is counted in each
# cohort's bar, so a stacked total would exceed the genus count and imply a sum
# that does not exist.
f6a_cov <- sm %>% filter(primary) %>%
  transmute(variable = VAR_LABEL[nhanes_variable],
            cohort = unname(COHORT_SHORT[dataset_id]),
            n = n_tested)
stopifnot(!any(is.na(f6a_cov$cohort)))
f6a_none <- cov_g %>% group_by(variable) %>%
  summarise(n = sum(n_datasets_testing == 0), .groups = "drop") %>%
  mutate(cohort = "no external cohort")
f6a_hits <- cov_g %>% count(variable, name = "n_hits")
f6a_d <- bind_rows(f6a_cov, f6a_none) %>%
  mutate(variable = factor(variable, levels = unname(VAR_LABEL[PRIMARY_VARIABLES])),
         cohort = factor(cohort, levels = rev(c(unname(COHORT_SHORT), "no external cohort")))) %>%
  filter(!is.na(cohort)) %>%
  droplevels()
# Every per-cohort bar must equal the count Panel A-1 draws, or the two figures
# disagree. Assert it here rather than discovering it in review.
stopifnot(nrow(dplyr::anti_join(
  f6a_cov, sm %>% filter(primary) %>%
    transmute(variable = VAR_LABEL[nhanes_variable],
              cohort = unname(COHORT_SHORT[dataset_id]), n = n_tested),
  by = c("variable", "cohort", "n"))) == 0)

# Grey for the uncovered genera so "no evidence" never reads as a result.
cov_cols <- setNames(c("#CC6677", "#0072B2", "#4A8CBF", "#9EC5DE"), unname(COHORT_SHORT))
cov_cols <- c(cov_cols, `no external cohort` = "#DDDDDD")
# One row per cohort, faceted by phenotype - deliberately the same layout as
# Panel A-1, so the two figures can be read against each other without re-orienting.
# The cohort names are on the axis, so no colour legend is needed.
f6a <- ggplot(f6a_d, aes(n, cohort, fill = cohort)) +
  geom_col(width = 0.66, colour = "white", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.25, size = GEOM_TEXT_SIZE,
            family = PDF_FAMILY, colour = "grey15") +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(variable = function(x) unname(VAR_SHORT[
               match(x, unname(VAR_LABEL))]))) +
  scale_fill_manual(values = cov_cols, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "NHANES genera that cohort can test (count)", y = NULL,
       title = m12_wrap("i  How many NHANES genera each external cohort can test", 104, bold = TRUE),
       subtitle = m12_wrap(paste0(
         "Counts are GENERA, not participants, and are the same numbers as Panel A-1. ",
         "Grey = genera no cohort carries. A genus tested by several cohorts is counted in each ",
         "bar, so bars do not sum to the ",
         paste(sprintf("%d %s", f6a_hits$n_hits, tolower(sub(" .*", "", f6a_hits$variable))),
               collapse = " / "), " NHANES genera."), 104)) +
  theme_compact() +
  theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0))

# =============================================================================
# Panel (b): the multiplicity view. Panel (a) counts genus-cohort comparisons,
# so a genus carried by three cohorts is counted three times; this panel counts
# each genus ONCE and splits it by how many cohorts carry it. The two answer
# different questions - "how much can each cohort see?" versus "how many
# independent looks does each genus get?" - and only (i) is comparable with Panel A-1.
# =============================================================================
f6mult_d <- cov_g %>%
  count(variable, n_datasets_testing) %>%
  group_by(variable) %>% mutate(frac = n / sum(n), tot = sum(n)) %>% ungroup() %>%
  mutate(variable = factor(variable, levels = unname(VAR_LABEL[PRIMARY_VARIABLES])),
         nd = factor(n_datasets_testing, levels = rev(sort(unique(n_datasets_testing)))))
# Grey for zero so "no evidence" never reads as a result; saturation increases
# with the number of INDEPENDENT cohorts carrying the genus.
mult_cols <- c(`0` = "#DDDDDD", `1` = "#9EC5DE", `2` = "#4A8CBF", `3` = "#0067A5")
# Dodged, not stacked, and on a COUNT axis: as a share-stack the "3 cohorts"
# class for age is 1 of 54 genera, a 2%-wide sliver whose label is unreadable.
# One bar per class keeps every count visible, which is the point of the panel.
f6mult <- ggplot(f6mult_d, aes(n, nd, fill = nd)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25, show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.3, size = GEOM_TEXT_SIZE,
            family = PDF_FAMILY, colour = "grey15") +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(variable = function(x) unname(VAR_SHORT[
               match(x, unname(VAR_LABEL))]))) +
  scale_fill_manual(values = mult_cols, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  # No y-axis title: rotated and centred on a short panel it overruns the top
  # and collides with the subtitle, which can carry the same meaning in-line.
  labs(x = "NHANES genera (count)", y = NULL,
       title = m12_wrap("ii  How many cohorts test each NHANES genus", 104, bold = TRUE),
       subtitle = m12_wrap(paste0(
         "Row = number of cohorts testing the genus (0 = none). Each genus counted ONCE, ",
         "unlike (i)."), 104)) +
  theme_compact() +
  # Both the subtitle and the rotated y-axis title live in the leftmost column
  # (plot.title.position = "plot"), so without this gap they touch.
  theme(panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0))

# Are the untested genera the marginal ones? If so, the tested subset is a
# biased, optimistic sample and concordance is overstated.
f6b_d <- cov_g %>% mutate(cov = factor(ifelse(n_datasets_testing == 0, "no", "yes"),
                                       levels = c("yes", "no")))
wil <- f6b_d %>% group_by(variable) %>%
  summarise(p = tryCatch(stats::wilcox.test(abs(nh_estimate) ~ cov)$p.value,
                         error = function(e) NA_real_), .groups = "drop") %>%
  mutate(lab = ifelse(is.na(p), "n/a", paste0("p = ", format.pval(p, digits = 2))))
f6b <- ggplot(f6b_d, aes(cov, abs(nh_estimate))) +
  geom_boxplot(aes(fill = cov), outlier.shape = NA, linewidth = 0.25, width = 0.5,
               alpha = 0.35, show.legend = FALSE) +
  geom_jitter(aes(fill = phylum), shape = 21, stroke = 0.15, size = 1.0,
              width = 0.16, height = 0, alpha = 0.9) +
  geom_text(data = wil, aes(x = 1.5, y = Inf, label = lab), vjust = 1.5,
            size = GEOM_TEXT_SIZE, family = PDF_FAMILY, colour = "grey25", inherit.aes = FALSE) +
  # Short phenotype names: the full "Smoking (current vs never)" does not fit a
  # one-third-width strip at 6 pt and was being clipped mid-word.
  facet_wrap(~ variable, scales = "free_y",
             labeller = labeller(variable = function(x) unname(VAR_SHORT[
               match(x, unname(VAR_LABEL))]))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  scale_fill_manual(values = c(pal, no = "grey75", yes = "grey55"),
                    breaks = names(pal), name = "Phylum",
                    guide = guide_legend(ncol = 2, override.aes = list(size = 1.4))) +
  labs(x = "Genus carried by at least one external cohort", y = "|NHANES effect| (CLR)",
       title = m12_wrap("iii  Are the genera no cohort tests the weaker associations?", 104, bold = TRUE),
       subtitle = m12_wrap("1 point = 1 NHANES genus; Wilcoxon rank-sum on |effect|, covered vs uncovered", 104)) +
  theme_compact() + theme(legend.position = "bottom")

f6 <- f6a / f6mult / f6b + patchwork::plot_layout(heights = c(1.15, 1.15, 2.0))
save_validation_figure(S15_FIG[["A-3"]], f6, width_mm = 104, height_mm = 142)


# =============================================================================
# Panel C - Is non-replication structured, or is it noise?
#
# If the genera that fail to replicate were concentrated in low-prevalence taxa
# or in one phylum, that would point to a detection artefact or a real
# taxon-specific difference. Both strata are shown with exact binomial CIs so
# the reader can see whether any apparent pattern survives the small n.
# =============================================================================
strat <- cmp %>% filter(tier == "A", primary, tested, nhanes_variable == "age")
strat$phylum <- m12_phylum_std(strat$phylum)
qs <- unname(stats::quantile(strat$prevalence_overall, 0:3 / 3))
# "low / mid / high" is meaningless without the cut points, and the tertiles are
# data-dependent so they cannot be written by hand. Carry the actual prevalence
# range in the label itself.
stopifnot(length(unique(qs)) == 4L)
# Comma, not parentheses: bin_by() appends "(k/n)", and two bracketed groups
# in one label read as a nested expression.
TERT_LAB <- sprintf("%s, %.0f-%.0f%%", c("low", "mid", "high"),
                    100 * qs[1:3], 100 * qs[2:4])
strat <- strat %>%
  mutate(prev_tertile = cut(prevalence_overall, qs, include.lowest = TRUE,
                            labels = TERT_LAB))
bin_by <- function(d, key) {
  d %>% group_by(.grp = .data[[key]]) %>%
    summarise(n = n(), k = sum(outcome == "concordant"), .groups = "drop") %>%
    bind_cols(purrr::map2_dfr(.$k, .$n, ~ m12_binom(.x, .y) %>% select(pct, ci_lo, ci_hi))) %>%
    mutate(lab = sprintf("%s (%d/%d)", .grp, k, n))
}
suppressPackageStartupMessages(library(purrr))
# Short facet labels: the strip is vertical and a long title is truncated.
s_prev <- bin_by(strat, "prev_tertile") %>% mutate(facet = "Prevalence", col = "grey35")
s_phy  <- bin_by(strat, "phylum")       %>%
  mutate(facet = "Phylum", col = unname(phylum_colors[as.character(.grp)]))
s_all  <- bind_rows(s_prev, s_phy) %>%
  mutate(lab = factor(lab, levels = rev(unique(lab))),
         facet = factor(facet, levels = c("Prevalence", "Phylum")))
# Colour the phylum rows with the canonical palette so this figure is
# colour-comparable with Panels B/E/A-3 and modules 4/5/7/8/9.
f8 <- ggplot(s_all, aes(pct, lab)) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey55", linewidth = 0.25) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.3,
                 colour = "grey45") +
  geom_point(aes(fill = I(col)), shape = 21, size = 1.4, stroke = 0.2) +
  facet_grid(facet ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Genera agreeing in direction with NHANES (%, exact 95% CI)", y = NULL,
       title = m12_wrap("Age concordance by stratum, Huang 2020 (the only refit cohort)", 84, bold = TRUE),
       subtitle = m12_wrap("Unit = genus; (concordant / tested genera) after each stratum. Intervals are wide by construction", 84)) +
  theme_compact()
save_validation_figure(S15_FIG[["C"]], f8, width_mm = 84, height_mm = 54)

vsay("")
vsay("all figures: ", paste(list.files(VAL_FIG_PDF, pattern = "\\.pdf$"), collapse = ", "))
