#!/usr/bin/env Rscript
################################################################################
##   STANDALONE DIFFERENTIAL ABUNDANCE HEATMAP PLOTTER                        ##
##   - Uses pre-computed phyloseq objects and factorized categorical datasets ##
##   - Minimal computation, maximum reuse of saved data                       ##
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(dplyr)
  library(forcats)
  library(stringr)
  library(tibble)
})

###############################################################################
## 1. PATH CONFIGURATION                                                      ##
###############################################################################

# Base paths (adjust these to your system)
base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
phyloseq_obj_path <- file.path(base_path, "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate")
viz_out_path <- file.path(base_path, "results/analyses_results/04_association_phyloseq_out/visualizations")
config_dir_path <- file.path(base_path, "config")

# Output directory
out_root <- file.path(viz_out_path, "categories_heatmap_standalone")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

cat("=== DIFFERENTIAL ABUNDANCE HEATMAP PLOTTER ===\n")
cat("Base path:", base_path, "\n")
cat("Phyloseq objects:", phyloseq_obj_path, "\n")
cat("Output directory:", out_root, "\n\n")

###############################################################################
## 2. LOAD PRE-COMPUTED PHYLOSEQ OBJECTS                                     ##
###############################################################################

cat("Loading phyloseq objects...\n")

phyloseq_files <- c(
  "none"    = "ubiome_relative_none.rds",
  "clr"     = "ubiome_relative_clr.rds",
  "lognorm" = "ubiome_relative_lognorm.rds"
)

phylo_sets <- list()
for (trans in names(phyloseq_files)) {
  file_path <- file.path(phyloseq_obj_path, phyloseq_files[[trans]])
  if (file.exists(file_path)) {
    phylo_sets[[trans]] <- readRDS(file_path)
    cat("  ✓", trans, "- Samples:", nsamples(phylo_sets[[trans]]), 
        "Taxa:", ntaxa(phylo_sets[[trans]]), "\n")
  } else {
    cat("  ✗", trans, "- File not found:", file_path, "\n")
  }
}

###############################################################################
## 3. LOAD AND FACTORIZE CATEGORICAL DATASETS                                ##
###############################################################################

cat("\nLoading variable configuration files...\n")

# Load variable lists from config files
load_vars <- function(file_path) {
  if (!file.exists(file_path)) return(character(0))
  vars <- readLines(file_path, warn = FALSE)
  vars[nzchar(vars) & !grepl("^#", vars)]
}

var_files <- c(
  "demoWAS" = file.path(config_dir_path, "1_demoWAS_vars.txt"),
  "oradWAS" = file.path(config_dir_path, "2_oradWAS_vars.txt"),
  "exWAS"   = file.path(config_dir_path, "3_exWAS_vars.txt"),
  "pheWAS"  = file.path(config_dir_path, "4_pheWAS_vars.txt"),
  "outWAS"  = file.path(config_dir_path, "5_outWAS_vars.txt")
)

# Extract sample data from phyloseq objects
# Use the first available phyloseq object to get sample data
ps_ref <- phylo_sets[[1]]
sample_data_full <- data.frame(sample_data(ps_ref)) %>%
  rownames_to_column("SEQN")

cat("Sample data loaded. Dimensions:", dim(sample_data_full), "\n")

# Get variable subsets for each analysis type
demoWAS_vars <- purrr::map(var_files, load_vars)$demoWAS %>% unique()
oradWAS_vars <- purrr::map(var_files, load_vars)$oradWAS %>% unique()
exWAS_vars   <- purrr::map(var_files, load_vars)$exWAS %>% unique()
pheWAS_vars  <- purrr::map(var_files, load_vars)$pheWAS %>% unique()
outWAS_vars  <- purrr::map(var_files, load_vars)$outWAS %>% unique()

# Create subsets
sample_data_demoWAS_vars_subset <- sample_data_full %>% 
  select(SEQN, any_of(demoWAS_vars))
sample_data_oradWAS_vars_subset <- sample_data_full %>% 
  select(SEQN, any_of(oradWAS_vars))
sample_data_exWAS_vars_subset <- sample_data_full %>% 
  select(SEQN, any_of(exWAS_vars))
sample_data_pheWAS_vars_subset <- sample_data_full %>% 
  select(SEQN, any_of(pheWAS_vars))
sample_data_outWAS_vars_subset <- sample_data_full %>% 
  select(SEQN, any_of(outWAS_vars))

cat("\nVariable subsets created:\n")
cat("  demoWAS:", ncol(sample_data_demoWAS_vars_subset) - 1, "vars\n")
cat("  oradWAS:", ncol(sample_data_oradWAS_vars_subset) - 1, "vars\n")
cat("  exWAS:", ncol(sample_data_exWAS_vars_subset) - 1, "vars\n")
cat("  pheWAS:", ncol(sample_data_pheWAS_vars_subset) - 1, "vars\n")
cat("  outWAS:", ncol(sample_data_outWAS_vars_subset) - 1, "vars\n")

###############################################################################
## 4. FACTORIZE CATEGORICAL VARIABLES WITH REFERENCE LEVELS                  ##
###############################################################################

cat("\nFactorizing categorical variables...\n")

# Pre-defined factors from phyloseq sample data
pre_defined_factors <- sample_data_full %>%
  mutate(across(where(is.character), as.factor)) %>%
  select(SEQN, any_of(c("Gender", "AgeGroup", "EducationLevel", "Ethnicity", 
                        "US_Born", "Household_Size_Factor", "Marital_Status", 
                        "Interview_Language")))

# DEMO-WAS factorization
demoWAS_vars_subset_factorized_for_heatmapheatmap <- sample_data_demoWAS_vars_subset %>%
  left_join(pre_defined_factors, by = "SEQN") %>%
  mutate(
    Gender = fct_relevel(as.factor(Gender), "Female"),
    Age_group = fct_relevel(as.factor(AgeGroup), "30-39"),
    Education_level = factor(EducationLevel, levels = c("< 9th Grade", "9-11th Grade", 
                                                         "High School", "College/AA", 
                                                         "College Graduate")) %>% fct_relevel("College/AA"),
    Ethnicity = fct_relevel(as.factor(Ethnicity), "White"),
    US_born = fct_relevel(as.factor(US_Born), "US Born")
  ) %>%
  select(SEQN, any_of(c("Gender", "Age_group", "Education_level", "Ethnicity", "US_born")))

# ORAD-WAS factorization (binary disease variables with control)
SEQN_oradWAS_control <- sample_data_oradWAS_vars_subset %>%
  filter(if_all(all_of(setdiff(oradWAS_vars, "SEQN")), ~ !is.na(.x) & .x == 0)) %>%
  pull(SEQN)

oradWAS_vars_subset_factorized_for_heatmapheatmap <- sample_data_oradWAS_vars_subset %>%
  mutate(
    Denture = case_when(
      get("DENTURE_OHAROCDE", .) == 1 ~ "Denture",
      SEQN %in% SEQN_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Denture")),
    
    Gum_disease = case_when(
      get("GUM_DISEASE_OHAROCGP", .) == 1 ~ "Gum disease",
      SEQN %in% SEQN_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Gum disease")),
    
    Oral_hygiene = case_when(
      get("ORAL_HYGIENE_OHAROCOH", .) == 1 ~ "Poor oral hygiene",
      SEQN %in% SEQN_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Poor oral hygiene")),
    
    Tooth_decay = case_when(
      get("TOOTH_DECAY_OHAROCDT", .) == 1 ~ "Tooth decay",
      SEQN %in% SEQN_oradWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Tooth decay"))
  ) %>%
  select(SEQN, any_of(c("Denture", "Gum_disease", "Oral_hygiene", "Tooth_decay")))

# PHEWAS factorization (BMI, blood pressure categories)
pheWAS_vars_subset_factorized_for_heatmap <- sample_data_pheWAS_vars_subset %>%
  mutate(
    BMI_category = case_when(
      get("BMXBMI", .) < 18.5 ~ "Underweight",
      get("BMXBMI", .) >= 18.5 & get("BMXBMI", .) < 25 ~ "Healthy weight",
      get("BMXBMI", .) >= 25 & get("BMXBMI", .) < 30 ~ "Overweight",
      get("BMXBMI", .) >= 30 & get("BMXBMI", .) < 35 ~ "Class 1 Obesity",
      get("BMXBMI", .) >= 35 ~ "Class 2-3 Obesity",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Underweight", "Healthy weight", "Overweight", 
                           "Class 1 Obesity", "Class 2-3 Obesity")) %>% fct_relevel("Healthy weight")
  ) %>%
  select(SEQN, any_of("BMI_category"))

# OUTWAS factorization (disease outcomes)
SEQN_outWAS_control <- sample_data_outWAS_vars_subset %>%
  filter(if_all(all_of(setdiff(outWAS_vars, "SEQN")), ~ !is.na(.x) & .x == 0)) %>%
  pull(SEQN)

outdWAS_vars_subset_factorized_for_heatmapheatmap <- sample_data_outWAS_vars_subset %>%
  mutate(
    Asthma = case_when(
      get("ASTHMA", .) == 1 ~ "Asthma",
      SEQN %in% SEQN_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Asthma")),
    
    Diabetes = case_when(
      get("DIABETES", .) == 1 ~ "Diabetes",
      SEQN %in% SEQN_outWAS_control ~ "control",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("control", "Diabetes"))
  ) %>%
  select(SEQN, any_of(c("Asthma", "Diabetes")))

# EXWAS factorization (smoking status)
exWAS_vars_subset_factorized_for_heatmap <- sample_data_exWAS_vars_subset %>%
  mutate(
    Smoking_status = case_when(
      get("SMQ_current_ever_never", .) == 0 ~ "Never smoker",
      get("SMQ_current_ever_never", .) == 2 ~ "Former smoker",
      get("SMQ_current_ever_never", .) == 1 ~ "Current smoker",
      TRUE ~ NA_character_
    ) %>% factor(levels = c("Never smoker", "Former smoker", "Current smoker")) %>% fct_relevel("Never smoker")
  ) %>%
  select(SEQN, any_of("Smoking_status"))

# Combine into datasets list
datasets <- list(
  demoWAS = demoWAS_vars_subset_factorized_for_heatmapheatmap,
  oradWAS = oradWAS_vars_subset_factorized_for_heatmapheatmap,
  exWAS   = exWAS_vars_subset_factorized_for_heatmap,
  pheWAS  = pheWAS_vars_subset_factorized_for_heatmap,
  outWAS  = outdWAS_vars_subset_factorized_for_heatmapheatmap
)

cat("Factorized datasets created:\n")
for (ds_name in names(datasets)) {
  cat("  ", ds_name, ":", ncol(datasets[[ds_name]]) - 1, "categorical vars\n")
}

###############################################################################
## 5. PHYLUM COLORS & HEATMAP PALETTE                                        ##
###############################################################################

phylum_base_col <- c(
  Firmicutes = "#F38400", Bacteroidetes = "#0067A5", Actinobacteria = "#8DB600",
  Proteobacteria = "#E68FAC", Fusobacteria = "#BE0032", Spirochaetes = "#F3C300",
  Cyanobacteria = "#875692", Acidobacteria = "#F6A600", `Candidate division SR1` = "#2B3D26",
  Planctomycetes = "#332288", Saccharibacteria = "#B3446C", Synergistetes = "#A1CAF1",
  Tenericutes = "#654522", Verrucomicrobia = "#C2B280", unclassified = "#DDDDDD",
  unknown = "#999999", TM7 = "#000000"
)

colors <- colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(100)

breaks_fixed <- list(
  clr     = seq(-2.5, 2.5, length.out = 101),
  none    = seq(-1.5, 1.5, length.out = 101),
  lognorm = seq(-1, 1, length.out = 101)
)

###############################################################################
## 6. HELPER FUNCTIONS                                                       ##
###############################################################################

nice <- function(x) gsub("_", " ", x)

tax_tbl <- function(ps){
  df <- as.data.frame(unclass(tax_table(ps)), stringsAsFactors = FALSE)
  df$OTU <- rownames(df)
  df
}

keep_taxa <- function(ps, min_prev = 0.01, trans = "none"){
  m <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) m <- t(m)
  ok1 <- if (trans %in% c("clr", "lognorm"))
    apply(m, 1, var, na.rm = TRUE) > 0
  else
    rowSums(m > 0) / ncol(m) >= min_prev
  ok2 <- !is.na(tax_tbl(ps)$Genus)
  which(ok1 & ok2)
}

p_one <- function(x, g){
  tryCatch({
    ok <- !is.na(x) & !is.na(g)
    if (sum(ok) < 4) return(1)
    x <- x[ok]; g <- droplevels(g[ok])
    if (nlevels(g) < 2 || var(x) == 0) return(1)
    if (nlevels(g) == 2) suppressWarnings(wilcox.test(x ~ g)$p.value)
    else kruskal.test(x, g)$p.value
  }, error = function(e) 1)
}

make_phylum_col <- function(v, alpha = 0.65){
  v[is.na(v) | v == ""] <- "unknown"
  need <- unique(v)
  miss <- setdiff(need, names(phylum_base_col))
  extra <- if (length(miss)){
    qual <- RColorBrewer::brewer.pal(max(3, length(miss)), "Set3")
    setNames(qual[seq_along(miss)], miss)
  } else character(0)
  cols <- c(phylum_base_col, extra)[need]
  grDevices::adjustcolor(cols, alpha.f = alpha)
}

draw_hm <- function(mat, title, key, anno_row, anno_cols){
  h <- max(5, ceiling(nrow(mat) * 5 / 16))
  brks <- breaks_fixed[[key]]
  lgd <- c(min(brks), 0, max(brks))
  
  args <- list(
    mat, cluster_rows = nrow(mat) > 1, clustering_distance_rows = "euclidean",
    cluster_cols = FALSE, treeheight_row = h, color = colors, breaks = brks,
    legend_breaks = lgd, legend_labels = format(lgd, trim = TRUE),
    annotation_row = anno_row, annotation_colors = anno_cols,
    cellwidth = 5.5, cellheight = 5.5, fontsize = 6,
    fontsize_row = 6, fontsize_col = 6, main = title, silent = TRUE
  )
  if (nrow(mat) > 3) args$cutree_rows <- 4
  do.call(pheatmap, args)
}

###############################################################################
## 7. ANALYZE ONE FACTOR                                                     ##
###############################################################################

analyze_factor <- function(otu, g, ps, trans, top_n = 30, alpha = 0.05){
  p_vec <- apply(otu, 1, p_one, g = g)
  q_vec <- p.adjust(p_vec, "BH")
  
  levels_g <- levels(g)
  means_all <- sapply(levels_g, \(lv) rowMeans(otu[, g == lv, drop = FALSE], na.rm = TRUE))
  if (is.null(dim(means_all)))
    means_all <- matrix(means_all, nrow = 1, dimnames = list(rownames(otu), levels_g))
  
  tx <- tax_tbl(ps)[match(rownames(otu), tax_tbl(ps)$OTU), ]
  genus <- ifelse(is.na(tx$Genus)|tx$Genus%in%c("", "unclassified"),
                  paste0("Uncl ", tx$OTU), nice(tx$Genus))
  phylum <- ifelse(is.na(tx$Phylum)|tx$Phylum=="", "unknown", tx$Phylum)
  
  res_tab <- data.frame(
    OTU = tx$OTU, Genus = genus, Phylum = phylum,
    p_value = p_vec, p_adj_BH = q_vec, means_all,
    in_top_panel = FALSE, stringsAsFactors = FALSE
  )
  
  sig_idx <- which(q_vec <= alpha)
  if (!length(sig_idx)) return(list(plots=list(), table=res_tab))
  
  means_sig <- means_all[sig_idx, , drop = FALSE]
  rng <- apply(means_sig, 1, \(z) diff(range(z, na.rm = TRUE)))
  keep <- head(order(rng, decreasing = TRUE), min(top_n, length(rng)))
  res_tab$in_top_panel[sig_idx[keep]] <- TRUE
  means_top <- means_sig[keep, , drop = FALSE]
  
  row_lab <- make.unique(genus[sig_idx][keep])
  rownames(means_top) <- row_lab
  phylum_top <- phylum[sig_idx][keep]
  anno_r <- data.frame(Phylum = phylum_top, row.names = row_lab)
  anno_c <- list(Phylum = make_phylum_col(phylum_top))
  
  if (trans == "clr"){
    panel <- sweep(means_top, 1, means_top[,1], "-") / log(2)
    plt <- draw_hm(panel, paste(nice(attr(g,"varname")), "log₂ fold change"),
                   "clr", anno_r, anno_c)
    return(list(plots = list(log2FC = plt), table = res_tab))
  }
  
  if (trans == "lognorm"){
    panel <- sweep(means_top, 1, means_top[,1], "-")
    plt <- draw_hm(panel, paste(nice(attr(g,"varname")), "log₁₀ diff"),
                   "lognorm", anno_r, anno_c)
    return(list(plots = list(log10Diff = plt), table = res_tab))
  }
  
  panel <- t(scale(t(means_top)))
  plt <- draw_hm(panel, paste(nice(attr(g,"varname")), "z-score"),
                 "none", anno_r, anno_c)
  list(plots = list(Abundance = plt), table = res_tab)
}

###############################################################################
## 8. RUN ONE DATASET                                                        ##
###############################################################################

run_dataset <- function(ds_name, meta, ps, trans, out_dir, 
                        min_prev = 0.01, top_n = 30, alpha = 0.05){
  message("== ", nice(ds_name), " (", trans, ")")
  keep <- keep_taxa(ps, min_prev, trans)
  if (!length(keep)){ message("   no taxa left"); return(list()) }
  ps <- prune_taxa(taxa_names(ps)[keep], ps)
  
  otu <- as(otu_table(ps), "matrix"); if (!taxa_are_rows(ps)) otu <- t(otu)
  common <- intersect(colnames(otu), meta$SEQN)
  if (length(common) < 30){ message("   too few samples"); return(list()) }
  otu <- otu[, common]
  meta <- meta |> filter(SEQN %in% common) |> slice(match(common, SEQN))
  
  fac_vars <- setdiff(names(Filter(is.factor, meta)), "SEQN")
  if (!length(fac_vars)){ message("   no factor vars"); return(list()) }
  
  pl_all <- list()
  for (v in fac_vars){
    g <- droplevels(meta[[v]])
    if (nlevels(g) < 2 || sum(!is.na(g)) < 30){
      message("   ", nice(v), " skipped"); next
    }
    attr(g,"varname") <- v
    ana <- analyze_factor(otu, g, ps, trans, top_n, alpha)
    
    write.csv(ana$table,
              file.path(out_dir, paste0(ds_name, "_", v, "_", trans, "_table.csv")),
              row.names = FALSE)
    
    if (length(ana$plots)){
      nm <- names(ana$plots)[1]
      pl_all[[paste0(v, "_", nm)]] <- ana$plots[[1]]
      message("   ", nice(v), " ✔")
    } else message("   ", nice(v), " – none")
  }
  pl_all
}

###############################################################################
## 9. MASTER EXECUTION LOOP                                                  ##
###############################################################################

grand_total <- 0

for (tr in names(phylo_sets)){
  ps <- phylo_sets[[tr]]
  tr_dir <- file.path(out_root, tr)
  if (!dir.exists(tr_dir)) dir.create(tr_dir, recursive = TRUE)
  booklet <- list()
  
  for (ds in names(datasets)){
    pl <- run_dataset(ds, datasets[[ds]], ps, tr, tr_dir)
    if (length(pl)){
      for (nm in names(pl)){
        pdf(file.path(tr_dir, paste0(ds, "_", nm, "_", tr, ".pdf")),
            10, 7, family = "Helvetica")
        print(pl[[nm]]); dev.off()
      }
      booklet <- c(booklet, pl)
      grand_total <- grand_total + length(pl)
    }
  }
  
  if (length(booklet)){
    flat <- Filter(function(z) inherits(z, "gtable"),
                   unlist(booklet, recursive = TRUE))
    
    pdf(file.path(tr_dir, paste0("ALL_", tr, "_booklet.pdf")),
        12, 18, family = "Helvetica")
    print(marrangeGrob(flat, ncol = 2, nrow = 3))
    dev.off()
  }
}

cat("\n======== SUMMARY =========\n")
cat("Total individual PDFs:", grand_total, "\n")
cat("Outputs in:", out_root, "\n")
cat("==========================\n")

