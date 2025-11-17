#!/usr/bin/env Rscript
# -----------------------------------------------------------
# nhanes_omp_transformation.R
# Stage 0 ▸ build a transformed microbiome DB with exact schema clones
#   • copies RSV relative‑abundance tables (genus→phylum, cycles F & G)
#   • applies three transforms: none | clr | lognorm
#   • keeps taxa with prevalence ≥ --prev (default 0.001)
#   • **adds matching metadata rows** into `variable_names_epcf`, so that
#     downstream code can resolve the new tables automatically.
#   • writes output DB to   data/00_nhanes_omp_transformed_db/…
# -----------------------------------------------------------
# Usage:
#   Rscript nhanes_omp_transformation.R \
#     --in_db  path/to/nhanes_031725.sqlite \
#     --out_db path/to/nhanes_oral_transformed.sqlite \
#     [--prev  0.001]
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(fs)
})

# ── CLI options ─────────────────────────────────────────────
opt_list <- list(
  make_option(c("-i","--in_db"),  type = "character",
              default = "data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite",
              help    = "input SQLite DB"),
  make_option(c("-o","--out_db"), type = "character",
              default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite",
              help    = "output SQLite DB"),
  make_option(c("-p","--prev"),   type = "double",
              default = 0.001,
              help    = "prevalence threshold")
)
opt <- parse_args(OptionParser(option_list = opt_list))

message("⤷ in_db:  ", opt$in_db)
message("⤷ out_db: ", opt$out_db)
message("⤷ prev:   ", opt$prev)

# ── seed the output DB with an exact copy of the input DB ──
dir_create(dirname(opt$out_db), recurse = TRUE)
if (!file_exists(opt$out_db)) file_copy(opt$in_db, opt$out_db)

# ── transformation helper ──────────────────────────────────
transform_microbiome_data <- function(df, method, prevalence_threshold) {
  stopifnot("SEQN" %in% names(df))
  seqn <- df$SEQN
  otu  <- as.matrix(df[ , -1, drop = FALSE])
  keep <- colMeans(otu > 0) >= prevalence_threshold
  if (!any(keep)) stop("No taxa remain after prevalence filtering")
  otu  <- otu[, keep, drop = FALSE]
  pseudo <- min(otu[otu > 0]) / 1e6
  otu[otu == 0] <- pseudo
  if (method == "clr") {
    if (!requireNamespace("compositions", quietly = TRUE)) {
      install.packages("compositions", repos = "https://cloud.r-project.org")
    }
    otu <- compositions::clr(otu)
  } else if (method == "lognorm") {
    otu <- log(otu)
  }
  cbind(SEQN = seqn, as.data.frame(otu))
}

# ── helper: duplicate metadata rows into variable_names_epcf ─
clone_metadata_rows <- function(con, base_tbl, method) {
  new_tbl <- paste0(base_tbl, "_", method)
  meta <- tbl(con, "variable_names_epcf") %>%
            filter(Data.File.Name == base_tbl) %>%
            collect()
  if (nrow(meta) == 0) return(invisible(FALSE))
  meta_new <- meta %>%
    mutate(Data.File.Name        = new_tbl,
           Data.File.Description = paste0(Data.File.Description, "_", method))
  # only insert if not present
  exists <- tbl(con, "variable_names_epcf") %>%
              filter(Data.File.Name == new_tbl) %>%
              head(1) %>% collect()
  if (nrow(exists) == 0) {
    dbWriteTable(con, "variable_names_epcf", meta_new, append = TRUE)
  }
  TRUE
}

# ── open connections ───────────────────────────────────────
in_con  <- dbConnect(SQLite(), opt$in_db)
out_con <- dbConnect(SQLite(), opt$out_db)

# ── read, transform, and write new tables per taxon level ──
source_tables <- c(
  "DADA2RSV_GENUS_RELATIVE_F", "DADA2RSV_GENUS_RELATIVE_G",
  "DADA2RSV_FAMILY_RELATIVE_F", "DADA2RSV_FAMILY_RELATIVE_G",
  "DADA2RSV_ORDER_RELATIVE_F", "DADA2RSV_ORDER_RELATIVE_G",
  "DADA2RSV_CLASS_RELATIVE_F", "DADA2RSV_CLASS_RELATIVE_G",
  "DADA2RSV_PHYLUM_RELATIVE_F", "DADA2RSV_PHYLUM_RELATIVE_G"
)
methods <- c("none", "clr", "lognorm")

for (tbl in source_tables) {
  if (!tbl %in% dbListTables(in_con)) {
    warning("Table not found in input DB: ", tbl)
    next
  }
  raw <- dbReadTable(in_con, tbl)
  for (m in methods) {
    new_tbl <- paste0(tbl, "_", m)
    message(sprintf("[%s] %s → %s", Sys.time(), tbl, m))
    tr <- transform_microbiome_data(raw, m, opt$prev)
    dbWriteTable(out_con, new_tbl, tr, overwrite = TRUE)
    clone_metadata_rows(out_con, tbl, m)
    message("    • wrote ", nrow(tr), " rows; metadata updated")
  }
}

# ── cleanup ────────────────────────────────────────────────
dbDisconnect(in_con)
dbDisconnect(out_con)