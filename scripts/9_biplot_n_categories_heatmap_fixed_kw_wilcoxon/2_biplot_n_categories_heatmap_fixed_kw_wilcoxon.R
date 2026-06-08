#!/usr/bin/env Rscript

# =============================================================================
# BIPLOT + CATEGORIES HEATMAP PLOTTING
# =============================================================================
# Reads intermediates produced by 1_load_and_process_all_data_standalone.R and
# emits:
#   - fig1a abundance/prevalence biplot
#   - per-factor categorical heatmaps (Kruskal-Wallis / Wilcoxon, BH-FDR)
#     under categories_heatmap_fixed/{none,clr}/
#
# Environment: R >= 4.5 with phyloseq, dplyr, tidyr, stringr, purrr, tibble,
# readr, ggplot2, egg, gridExtra, grid, pheatmap, RColorBrewer, forcats,
# extrafont.
# Conda spec: envs/nhanes-analysis_for_reviewers.yml.
# =============================================================================

# === USER CONFIG ============================================================
# Update PROJECT_ROOT to the absolute path of your local clone of this repo.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(egg)
  library(gridExtra)
  library(grid)
  library(pheatmap)
  library(RColorBrewer)
  library(forcats)
  library(extrafont)
})

set.seed(42)

message("=== BIPLOT + CATEGORIES HEATMAP PLOTTING (STANDALONE) ===")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
base_path <- PROJECT_ROOT
output_root <- file.path(
  base_path,
  "results/analyses_results/9_biplot_n_categories_heatmap_fixed_kw_wilcoxon_out"
)
data_path <- file.path(output_root, "intermediate")

if (!dir.exists(data_path)) {
  stop("Intermediate directory not found: ", data_path,
       "\nRun 1_load_and_process_all_data_standalone.R first.")
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

required_intermediates <- c(
  "phyloseq_objects.rds",
  "heatmap_factorized_datasets.rds",
  "sample_data_subsets.rds"
)

missing_intermediates <- required_intermediates[
  !file.exists(file.path(data_path, required_intermediates))
]
if (length(missing_intermediates)) {
  stop(
    "Missing intermediate file(s):\n",
    paste("  -", missing_intermediates, collapse = "\n"),
    "\nRegenerate them with 1_load_and_process_all_data_standalone.R."
  )
}

# -----------------------------------------------------------------------------
# Load intermediates
# -----------------------------------------------------------------------------
message("Loading intermediates...")
phyloseq_objects <- readRDS(file.path(data_path, "phyloseq_objects.rds"))
heatmap_factorized_datasets <- readRDS(
  file.path(data_path, "heatmap_factorized_datasets.rds")
)
sample_data_subsets <- readRDS(file.path(data_path, "sample_data_subsets.rds"))

message("  [OK] phyloseq objects:             ", length(phyloseq_objects))
message("  [OK] heatmap factorised datasets:  ", length(heatmap_factorized_datasets))
message("  [OK] sample_data subsets:          ", length(sample_data_subsets))

# -----------------------------------------------------------------------------
# Fonts (ArialMT as in Chapter 4)
# -----------------------------------------------------------------------------
loadfonts(device = "pdf", quiet = TRUE)

# -----------------------------------------------------------------------------
# Shared colour palette
# -----------------------------------------------------------------------------
phylum_colors <- c(
  "Firmicutes"           = "#F38400",
  "Bacteroidetes"        = "#0067A5",
  "Actinobacteria"       = "#8DB600",
  "Proteobacteria"       = "#E68FAC",
  "Fusobacteria"         = "#BE0032",
  "Spirochaetae"         = "#F3C300",
  "Cyanobacteria"        = "#875692",
  "Acidobacteria"        = "#F6A600",
  "Candidate division SR1" = "#2B3D26",
  "Planctomycetes"       = "#332288",
  "Saccharibacteria"     = "#B3446C",
  "Synergistetes"        = "#A1CAF1",
  "Tenericutes"          = "#654522",
  "Verrucomicrobia"      = "#C2B280",
  "unclassified"         = "#DDDDDD",
  "unknown"              = "#999999",
  "TM7"                  = "#000000"
)

# -----------------------------------------------------------------------------
# Utility: PDF writer (mirrors Chapter 4)
# -----------------------------------------------------------------------------
save_pdf_figure <- function(file_name, grob_obj, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  file_path <- file.path(output_root, file_name)
  message("  - Saving ", file_path)
  pdf(
    file = file_path,
    width = width_in,
    height = height_in,
    family = "ArialMT",
    pointsize = 5,
    colormodel = "rgb",
    useDingbats = FALSE
  )
  grid.newpage()
  grid.draw(grob_obj)
  dev.off()
}

# -----------------------------------------------------------------------------
# Section 1: Taxonomic abundance/prevalence biplot
# -----------------------------------------------------------------------------
plot_top_taxa_abundance_with_prevalence_horizontal <- function(phy, rank_name, top_n = 15) {
  phy_agg <- tax_glom(phy, taxrank = rank_name)
  tax_df <- as.data.frame(tax_table(phy_agg), stringsAsFactors = FALSE)
  tax_df$OTU <- rownames(tax_df)

  abundance_mat <- as(otu_table(phy_agg), "matrix")
  if (!taxa_are_rows(phy_agg)) {
    abundance_mat <- t(abundance_mat)
  }

  summary_df <- data.frame(
    OTU = rownames(abundance_mat),
    mean_abundance = rowMeans(abundance_mat, na.rm = TRUE),
    prevalence = apply(
      abundance_mat,
      1,
      function(x) sum(x > 0, na.rm = TRUE)
    ) / ncol(abundance_mat),
    stringsAsFactors = FALSE
  ) %>%
    left_join(tax_df, by = "OTU") %>%
    mutate(
      x_label = if_else(
        rank_name == "Genus" & tolower(Genus) == "unclassified",
        paste0("NA; ", OTU),
        as.character(.data[[rank_name]])
      )
    )

  if (!rank_name %in% names(summary_df)) {
    summary_df <- summary_df %>% mutate(x_label = OTU)
  }

  top_df <- summary_df %>%
    arrange(desc(mean_abundance)) %>%
    head(top_n)

  ggplot(top_df) +
    geom_bar(
      aes(x = mean_abundance, y = reorder(x_label, mean_abundance), fill = Phylum),
      stat = "identity",
      alpha = 0.65
    ) +
    geom_bar(
      aes(x = -prevalence * 0.15, y = reorder(x_label, mean_abundance)),
      stat = "identity",
      fill = "grey70",
      alpha = 0.65
    ) +
    geom_text(
      aes(
        x = mean_abundance - 0.015,
        y = x_label,
        label = signif(mean_abundance, digits = 3)
      ),
      hjust = -0.2,
      size = 1.8
    ) +
    geom_text(
      aes(
        x = -prevalence * 0.15 + 0.015,
        y = x_label,
        label = signif(prevalence, digits = 3)
      ),
      hjust = 1.2,
      size = 1.8
    ) +
    geom_vline(xintercept = 0, colour = "black") +
    scale_fill_manual(values = phylum_colors) +
    scale_x_continuous(
      limits = c(-0.3, 0.6),
      breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
    ) +
    labs(
      title = paste("Top", top_n, rank_name, "by Mean Abundance and Prevalence"),
      x = "Prevalence (left), Mean Abundance (right)",
      y = rank_name
    ) +
    egg::theme_article(base_size = 6, base_family = "ArialMT") +
    theme(
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 6),
      axis.title.y = element_text(size = 6),
      plot.title = element_text(size = 6, hjust = 0.5),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 6)
    )
}

message("Generating abundance/prevalence biplot...")
tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus")
biplot_plots <- map(
  tax_ranks,
  ~ plot_top_taxa_abundance_with_prevalence_horizontal(
    phyloseq_objects[["ubiome_relative_none"]],
    rank_name = .x,
    top_n = 15
  )
)

biplot_grob <- gridExtra::arrangeGrob(grobs = biplot_plots, nrow = 1)
save_pdf_figure(
  file_name = "fig1a_Top_Taxa_Abundance_With_Prevalence_biplot.pdf",
  grob_obj = biplot_grob,
  width_mm = 410,
  height_mm = 43
)
message("  [OK] Biplot saved")

# -----------------------------------------------------------------------------
# Section 2: Categories heatmap (fixed thresholds)
# -----------------------------------------------------------------------------
message("Generating categorical heatmaps...")

# Helpers copied from Chapter 4 Rmd
nice <- function(x) gsub("_", " ", x)

tax_tbl <- function(ps) {
  df <- as.data.frame(unclass(tax_table(ps)), stringsAsFactors = FALSE)
  df$OTU <- rownames(df)
  df
}

keep_taxa <- function(ps, min_prev = 0.01, trans = "none") {
  m <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) m <- t(m)
  ok1 <- if (trans %in% c("clr", "lognorm")) {
    apply(m, 1, var, na.rm = TRUE) > 0
  } else {
    rowSums(m > 0) / ncol(m) >= min_prev
  }
  ok2 <- !is.na(tax_tbl(ps)$Genus)
  which(ok1 & ok2)
}

p_one <- function(x, g) {
  tryCatch({
    ok <- !is.na(x) & !is.na(g)
    if (sum(ok) < 4) return(1)
    x <- x[ok]
    g <- droplevels(g[ok])
    if (nlevels(g) < 2 || var(x) == 0) return(1)
    if (nlevels(g) == 2) {
      suppressWarnings(wilcox.test(x ~ g)$p.value)
    } else {
      kruskal.test(x, g)$p.value
    }
  }, error = function(e) 1)
}

make_phylum_col <- function(v, alpha = 0.65) {
  v[is.na(v) | v == ""] <- "unknown"
  need <- unique(v)
  miss <- setdiff(need, names(phylum_colors))
  extra <- if (length(miss)) {
    qual <- brewer.pal(max(3, length(miss)), "Set3")
    setNames(qual[seq_along(miss)], miss)
  } else character(0)
  cols <- c(phylum_colors, extra)[need]
  setNames(adjustcolor(cols, alpha.f = alpha), need)
}

colors_fixed <- colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(100)
breaks_fixed <- list(
  clr = seq(-2.5, 2.5, length.out = 101),
  none = seq(-1.5, 1.5, length.out = 101),
  lognorm = seq(-1, 1, length.out = 101)
)

draw_hm <- function(mat, title, key, anno_row, anno_cols) {
  h <- max(5, ceiling(nrow(mat) * 5 / 16))
  brks <- breaks_fixed[[key]]
  lgd <- c(min(brks), 0, max(brks))

  args <- list(
    mat = mat,
    cluster_rows = nrow(mat) > 1,
    clustering_distance_rows = "euclidean",
    cluster_cols = FALSE,
    treeheight_row = h,
    color = colors_fixed,
    breaks = brks,
    legend_breaks = lgd,
    legend_labels = format(lgd, trim = TRUE),
    annotation_row = anno_row,
    annotation_colors = anno_cols,
    cellwidth = 5.5,
    cellheight = 5.5,
    fontsize = 6,
    fontsize_row = 6,
    fontsize_col = 6,
    main = title,
    silent = TRUE
  )

  if (nrow(mat) > 3) {
    args$cutree_rows <- 4
  }

  do.call(pheatmap::pheatmap, args)
}

analyze_factor <- function(otu, g, ps, trans, top_n = 30, alpha = 0.05) {
  p_vec <- apply(otu, 1, p_one, g = g)
  q_vec <- p.adjust(p_vec, "BH")

  levels_g <- levels(g)
  means_all <- sapply(levels_g, function(lv) {
    rowMeans(otu[, g == lv, drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(means_all))) {
    means_all <- matrix(means_all, nrow = 1, dimnames = list(rownames(otu), levels_g))
  }

  tx <- tax_tbl(ps)[match(rownames(otu), tax_tbl(ps)$OTU), ]
  genus <- ifelse(
    is.na(tx$Genus) | tx$Genus %in% c("", "unclassified"),
    paste0("Uncl ", tx$OTU),
    nice(tx$Genus)
  )
  phylum <- ifelse(is.na(tx$Phylum) | tx$Phylum == "", "unknown", tx$Phylum)

  res_tab <- data.frame(
    OTU = tx$OTU,
    Genus = genus,
    Phylum = phylum,
    p_value = p_vec,
    p_adj_BH = q_vec,
    means_all,
    in_top_panel = FALSE,
    stringsAsFactors = FALSE
  )

  sig_idx <- which(q_vec <= alpha)
  if (!length(sig_idx)) {
    return(list(plots = list(), table = res_tab))
  }

  means_sig <- means_all[sig_idx, , drop = FALSE]
  rng <- apply(means_sig, 1, function(z) diff(range(z, na.rm = TRUE)))
  keep <- head(order(rng, decreasing = TRUE), min(top_n, length(rng)))
  res_tab$in_top_panel[sig_idx[keep]] <- TRUE
  means_top <- means_sig[keep, , drop = FALSE]

  row_lab <- make.unique(genus[sig_idx][keep])
  rownames(means_top) <- row_lab
  phylum_top <- phylum[sig_idx][keep]
  anno_r <- data.frame(Phylum = phylum_top, row.names = row_lab)
  anno_c <- list(Phylum = make_phylum_col(phylum_top))

  if (trans == "clr") {
    panel <- sweep(means_top, 1, means_top[, 1], "-") / log(2)
    plt <- draw_hm(panel, paste(nice(attr(g, "varname")), "log₂ fold change"), "clr", anno_r, anno_c)
    return(list(plots = list(log2FC = plt), table = res_tab))
  }

  if (trans == "lognorm") {
    panel <- sweep(means_top, 1, means_top[, 1], "-")
    plt <- draw_hm(panel, paste(nice(attr(g, "varname")), "log₁₀ diff"), "lognorm", anno_r, anno_c)
    return(list(plots = list(log10Diff = plt), table = res_tab))
  }

  panel <- t(scale(t(means_top)))
  plt <- draw_hm(panel, paste(nice(attr(g, "varname")), "z-score"), "none", anno_r, anno_c)
  list(plots = list(Abundance = plt), table = res_tab)
}

run_dataset <- function(ds_name, meta, ps, trans, out_dir,
                        min_prev = 0.01, top_n = 30, alpha = 0.05) {
  message("== ", nice(ds_name), " (", trans, ")")
  keep <- keep_taxa(ps, min_prev, trans)
  if (!length(keep)) {
    message("   no taxa left")
    return(list())
  }

  ps_filtered <- prune_taxa(taxa_names(ps)[keep], ps)
  otu <- as(otu_table(ps_filtered), "matrix")
  if (!taxa_are_rows(ps_filtered)) otu <- t(otu)

  common <- intersect(colnames(otu), meta$SEQN)
  if (length(common) < 30) {
    message("   too few samples")
    return(list())
  }

  otu <- otu[, common]
  meta <- meta %>%
    filter(SEQN %in% common) %>%
    slice(match(common, SEQN))

  fac_vars <- setdiff(names(Filter(is.factor, meta)), "SEQN")
  if (!length(fac_vars)) {
    message("   no factor vars")
    return(list())
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pl_all <- list()
  for (v in fac_vars) {
    g <- droplevels(meta[[v]])
    if (nlevels(g) < 2 || sum(!is.na(g)) < 30) {
      message("   ", nice(v), " skipped")
      next
    }
    attr(g, "varname") <- v
    ana <- analyze_factor(otu, g, ps_filtered, trans, top_n, alpha)

    csv_path <- file.path(out_dir, paste0(ds_name, "_", v, "_", trans, "_table.csv"))
    write.csv(ana$table, csv_path, row.names = FALSE)

    if (length(ana$plots)) {
      nm <- names(ana$plots)[1]
      plot_obj <- ana$plots[[1]]
      if (inherits(plot_obj, "list") && !is.null(plot_obj$gtable)) {
        pl_all[[paste0(v, "_", nm)]] <- plot_obj$gtable
      } else if (inherits(plot_obj, "gtable")) {
        pl_all[[paste0(v, "_", nm)]] <- plot_obj
      } else {
        pl_all[[paste0(v, "_", nm)]] <- plot_obj
      }
      message("   ", nice(v), " [OK]")
    } else {
      message("   ", nice(v), " – none")
    }
  }

  pl_all
}

make_dir <- function(x) if (!dir.exists(x)) dir.create(x, TRUE)

out_root <- file.path(output_root, "categories_heatmap_fixed")
make_dir(out_root)

# Process NONE + CLR transformations only.
phylo_sets <- list(
  none = phyloseq_objects[["ubiome_relative_none"]],
  clr  = phyloseq_objects[["ubiome_relative_clr"]]
)

grand_total <- 0

for (tr in names(phylo_sets)) {
  ps <- phylo_sets[[tr]]
  tr_dir <- file.path(out_root, tr)
  supp_dir <- file.path(tr_dir, "supplementary")
  make_dir(tr_dir)
  make_dir(supp_dir)
  booklet <- list()

  for (ds in names(heatmap_factorized_datasets)) {
    meta <- heatmap_factorized_datasets[[ds]]
    pl <- run_dataset(ds, meta, ps, tr, supp_dir)
    if (length(pl)) {
      for (nm in names(pl)) {
        pdf_path <- file.path(tr_dir, paste0(ds, "_", nm, "_", tr, ".pdf"))
        message("      - Saving ", pdf_path)
        pdf(pdf_path, width = 10, height = 7, family = "Helvetica")
        grid::grid.newpage()
        grid::grid.draw(pl[[nm]])
        dev.off()
      }
      booklet <- c(booklet, pl)
      grand_total <- grand_total + length(pl)
    }
  }

  if (length(booklet)) {
    flat <- Filter(function(z) inherits(z, "gtable"), booklet)
    if (length(flat)) {
      booklet_file <- file.path(tr_dir, paste0("ALL_", tr, "_booklet.pdf"))
      message("      - Saving ", booklet_file)
      pdf(booklet_file, width = 12, height = 18, family = "Helvetica")
      print(gridExtra::marrangeGrob(grobs = flat, ncol = 2, nrow = 3))
      dev.off()
    }
  }
}

message("\n======== HEATMAP SUMMARY ========")
message("Total individual heatmap panels: ", grand_total)
message("Outputs written to: ", out_root)
message("=================================")
message("=== PLOTTING COMPLETE ===")

