#!/usr/bin/env Rscript
################################################################################
# nhanes_omp_transformation.R
#
# Step 1: build a transformed microbiome SQLite database with four
# normalized versions of each abundance table (none, hellinger, clr, lognorm).
# No prevalence pre-filtering -- all taxa retained.
#
# Environment: SLURM module `gcc/14.2.0 R/4.4.2` (see run_transformation.sh).
# Locally the project conda env at envs/.conda/envs/nhanes-analysis also works.
#
# Run via the SLURM wrapper:
#   sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh
# or manually:
#   Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_omp_transformation.R \
#     --in_db  data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite \
#     --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(fs)
})

opt_list <- list(
  make_option(c("-i", "--in_db"),  type = "character",
              default = "data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite",
              help    = "input SQLite DB"),
  make_option(c("-o", "--out_db"), type = "character",
              default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite",
              help    = "output SQLite DB")
)
opt <- parse_args(OptionParser(option_list = opt_list))

if (!file.exists(opt$in_db)) stop("Input database not found: ", opt$in_db)
dir_create(dirname(opt$out_db), recurse = TRUE)
if (!file_exists(opt$out_db)) file_copy(opt$in_db, opt$out_db)

EPS <- 1  # integer pseudocount

trans_none <- function(rel_mat) {
  row_sums <- rowSums(rel_mat, na.rm = TRUE)
  if (any(abs(row_sums - 1) > 1e-6, na.rm = TRUE)) {
    n_unclosed <- sum(abs(row_sums - 1) > 1e-6, na.rm = TRUE)
    warning("CLOSURE: ", n_unclosed, " rows not closed to 1 (max deviation: ",
            round(max(abs(row_sums - 1), na.rm = TRUE), 6), ")")
  }
  rel_mat
}

trans_hellinger <- function(rel_mat) sqrt(rel_mat)

trans_lognorm <- function(cnt_mat) {
  cnt_mat  <- matrix(as.numeric(cnt_mat), nrow = nrow(cnt_mat))
  D        <- ncol(cnt_mat)
  n_i      <- rowSums(cnt_mat, na.rm = TRUE)
  mean_lib <- mean(n_i, na.rm = TRUE)
  zero_lib <- n_i == 0
  if (any(zero_lib)) {
    cnt_mat[zero_lib, ] <- NA_real_
    n_i[zero_lib]       <- NA_real_
  }
  lib_size  <- n_i + EPS * D
  scale_mat <- (cnt_mat + EPS) / lib_size
  log10(scale_mat * mean_lib)
}

trans_clr <- function(cnt_mat) {
  cnt_mat <- matrix(as.numeric(cnt_mat), nrow = nrow(cnt_mat))
  n_i <- rowSums(cnt_mat, na.rm = TRUE)
  zero_lib <- n_i == 0
  if (any(zero_lib)) cnt_mat[zero_lib, ] <- NA_real_
  log_cnt_eps <- log(cnt_mat + EPS)
  gm <- exp(rowMeans(log_cnt_eps, na.rm = TRUE))
  log_cnt_eps - log(gm)
}

add_seqn <- function(seqn_vec, mat, original_colnames = NULL) {
  if (is.matrix(mat)) {
    colnames_to_use <- if (!is.null(original_colnames)) original_colnames else colnames(mat)
    df_mat <- as.data.frame(mat, check.names = FALSE)
    if (!is.null(colnames_to_use) && length(colnames_to_use) == ncol(df_mat)) {
      colnames(df_mat) <- colnames_to_use
    }
  } else {
    df_mat <- as.data.frame(mat, check.names = FALSE)
  }
  cbind(SEQN = seqn_vec, df_mat)
}

clone_metadata_rows <- function(con, base_tbl, suffix) {
  new_tbl <- paste0(base_tbl, "_", suffix)
  meta    <- tbl(con, "variable_names_epcf") %>%
               filter(Data.File.Name == base_tbl) %>% collect()
  if (nrow(meta) == 0) return(invisible(FALSE))
  transformation_desc <- case_when(
    suffix == "clr"       ~ " (CLR ln-ratio transformed)",
    suffix == "lognorm"   ~ " (log10-normalized)",
    suffix == "hellinger" ~ " (Hellinger sqrt-transformed)",
    suffix == "none"      ~ " (untransformed copy)",
    TRUE                  ~ paste0(" (", suffix, " transformed)")
  )
  meta_new <- meta %>%
    mutate(
      Data.File.Name        = new_tbl,
      Data.File.Description = paste0(Data.File.Description, "_", suffix, transformation_desc)
    )
  if (tbl(con, "variable_names_epcf") %>%
        filter(Data.File.Name == new_tbl) %>% head(1) %>% collect() %>% nrow() == 0) {
    dbWriteTable(con, "variable_names_epcf", meta_new, append = TRUE)
  }
  TRUE
}

create_view_for_none <- function(con, original_tbl) {
  view_name <- paste0(original_tbl, "_none")
  existing_views <- dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type='view'")$name
  if (!view_name %in% existing_views) {
    dbExecute(con, sprintf("CREATE VIEW IF NOT EXISTS %s AS SELECT * FROM %s",
                           view_name, original_tbl))
    return(TRUE)
  }
  FALSE
}

in_con  <- dbConnect(SQLite(), opt$in_db)
out_con <- dbConnect(SQLite(), opt$out_db)

rel_tables <- c(
  "DADA2RSV_GENUS_RELATIVE_F",  "DADA2RSV_GENUS_RELATIVE_G",
  "DADA2RSV_FAMILY_RELATIVE_F", "DADA2RSV_FAMILY_RELATIVE_G",
  "DADA2RSV_ORDER_RELATIVE_F",  "DADA2RSV_ORDER_RELATIVE_G",
  "DADA2RSV_CLASS_RELATIVE_F",  "DADA2RSV_CLASS_RELATIVE_G",
  "DADA2RSV_PHYLUM_RELATIVE_F", "DADA2RSV_PHYLUM_RELATIVE_G"
)

count_tables <- c(
  "DADA2RSV_GENUS_COUNT_F",  "DADA2RSV_GENUS_COUNT_G",
  "DADA2RSV_FAMILY_COUNT_F", "DADA2RSV_FAMILY_COUNT_G",
  "DADA2RSV_ORDER_COUNT_F",  "DADA2RSV_ORDER_COUNT_G",
  "DADA2RSV_CLASS_COUNT_F",  "DADA2RSV_CLASS_COUNT_G",
  "DADA2RSV_PHYLUM_COUNT_F", "DADA2RSV_PHYLUM_COUNT_G"
)

for (tbl in rel_tables) {
  if (!tbl %in% dbListTables(in_con)) {
    warning("Relative table not found: ", tbl)
    next
  }
  df  <- dbReadTable(in_con, tbl)
  seq <- df$SEQN
  rel <- as.matrix(df[, -1, drop = FALSE])
  original_colnames <- colnames(df)[-1]

  if (create_view_for_none(out_con, tbl)) clone_metadata_rows(out_con, tbl, "none")

  out_hel <- add_seqn(seq, trans_hellinger(rel), original_colnames)
  dbWriteTable(out_con, paste0(tbl, "_hellinger"), out_hel, overwrite = TRUE)
  clone_metadata_rows(out_con, tbl, "hellinger")

  message(tbl, " -> none (view) / hellinger")
}

for (tbl in count_tables) {
  if (!tbl %in% dbListTables(in_con)) {
    warning("Count table not found: ", tbl)
    next
  }
  df  <- dbReadTable(in_con, tbl)
  seq <- df$SEQN
  cnt <- as.matrix(df[, -1, drop = FALSE])
  original_colnames <- colnames(df)[-1]
  relative_colnames <- gsub("_count$", "_relative", original_colnames, ignore.case = TRUE)
  rel_name          <- gsub("_COUNT_", "_RELATIVE_", tbl)

  out_clr <- add_seqn(seq, trans_clr(cnt), relative_colnames)
  dbWriteTable(out_con, paste0(rel_name, "_clr"), out_clr, overwrite = TRUE)
  clone_metadata_rows(out_con, rel_name, "clr")

  out_ln <- add_seqn(seq, trans_lognorm(cnt), relative_colnames)
  dbWriteTable(out_con, paste0(rel_name, "_lognorm"), out_ln, overwrite = TRUE)
  clone_metadata_rows(out_con, rel_name, "lognorm")

  message(tbl, " -> clr / lognorm")
}

dbDisconnect(in_con)
dbDisconnect(out_con)
