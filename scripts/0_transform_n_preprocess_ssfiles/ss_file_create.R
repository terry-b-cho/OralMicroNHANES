#!/usr/bin/env Rscript
################################################################################
# ss_file_create.R
#
# Step 3 worker: build a single `schema_structure` CSV mapping each
# microbiome variable to each NHANES variable with pairwise sample counts,
# for one (pipeline x transformation) combination. The runner script
# `run_ss_file_create.sh` submits 24 SLURM jobs covering all combinations.
#
# Environment: project conda env at envs/.conda/envs/nhanes-analysis (R >= 4.5).
# Activated by the SLURM wrapper.
#
# Usage:
#   Rscript scripts/0_transform_n_preprocess_ssfiles/ss_file_create.R \
#     --db        data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite \
#     --otu_F     DADA2RSV_GENUS_RELATIVE_F_clr \
#     --otu_G     DADA2RSV_GENUS_RELATIVE_G_clr \
#     --vars_file configs/3_exWAS_vars.txt \
#     --otu_role  dep \
#     --pipeline  3_exWAS \
#     --transform clr \
#     --out_dir   results/0_ss_files
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(readr)
  library(purrr)
  library(fs)
  library(glue)
  library(stringr)
})

opt_parser <- OptionParser(option_list = list(
  make_option("--db",        type = "character",
              default = "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite",
              help    = "path to transformed DB"),
  make_option("--otu_F",     type = "character", help = "cycle-F microbiome table"),
  make_option("--otu_G",     type = "character", help = "cycle-G microbiome table"),
  make_option("--vars_file", type = "character",
              help = "file listing NHANES variable names (one per line)"),
  make_option("--otu_role",  type = "character", default = "indep",
              help = "indep | dep (is microbiome independent?)"),
  make_option("--pipeline",  type = "character", default = "3_exWAS"),
  make_option("--transform", type = "character", default = "none"),
  make_option("--out_dir",   type = "character", default = "results/0_ss_files")
))
opt <- parse_args(opt_parser)

stopifnot(opt$otu_role %in% c("indep", "dep"))
dir_create(opt$out_dir, recurse = TRUE)

con        <- dbConnect(SQLite(), opt$db)
all_tables <- dbListTables(con)

if (!"variable_names_epcf" %in% all_tables) stop("Mapping table 'variable_names_epcf' not found in DB")
if (!opt$otu_F %in% all_tables) stop("OTU table not found: ", opt$otu_F)
if (!opt$otu_G %in% all_tables) stop("OTU table not found: ", opt$otu_G)

get_var_table_for_cycle <- function(var, cycle) {
  lc <- tolower(cycle)
  qry <- dbGetQuery(con, sprintf(
    "SELECT \"Data.File.Name\" AS df
     FROM variable_names_epcf
     WHERE \"Variable.Name\" = '%s'
       AND lower(\"Data.File.Name\") LIKE '%%%s'
     LIMIT 1",
    var, paste0("_", lc)))
  if (nrow(qry) == 0) return(NA_character_)
  qry$df[1]
}

read_otu <- function(tbl) {
  if (!tbl %in% all_tables) stop("Missing OTU table: ", tbl)
  dbReadTable(con, tbl)
}

compute_counts <- function(dep_df, indep_df,
                           dep_prefix, indep_prefix,
                           dep_lookup, indep_lookup) {
  dep_vars   <- setdiff(names(dep_df), "SEQN")
  indep_vars <- setdiff(names(indep_df), "SEQN")
  merged     <- inner_join(dep_df, indep_df, by = "SEQN")
  map_dfr(dep_vars, function(dv) {
    dep_tbl <- dep_lookup[[dv]]
    d_vec   <- merged[[dv]]
    map_dfr(indep_vars, function(iv) {
      indep_tbl <- indep_lookup[[iv]]
      tibble(
        !!paste0(dep_prefix,   "var")   := dv,
        !!paste0(dep_prefix,   "table") := dep_tbl,
        !!paste0(indep_prefix, "var")   := iv,
        !!paste0(indep_prefix, "table") := indep_tbl,
        n = sum(!is.na(d_vec) & !is.na(merged[[iv]]))
      )
    })
  })
}

varlist <- read_lines(opt$vars_file) %>% str_trim() %>% discard(~ .x == "")

otu_F <- read_otu(opt$otu_F)
otu_G <- read_otu(opt$otu_G)

var2table <- list(
  F = setNames(map_chr(varlist, ~ get_var_table_for_cycle(.x, "F")), varlist),
  G = setNames(map_chr(varlist, ~ get_var_table_for_cycle(.x, "G")), varlist)
)

missing_both <- varlist[is.na(var2table$F) & is.na(var2table$G)]
if (length(missing_both) > 0) {
  warning("Variables not found in any cycle: ", paste(missing_both, collapse = ", "))
  varlist     <- setdiff(varlist, missing_both)
  var2table$F <- var2table$F[varlist]
  var2table$G <- var2table$G[varlist]
}

load_vars_cycle <- function(cycle) {
  dfs <- imap(var2table[[cycle]], function(tbl_nm, var) {
    if (is.na(tbl_nm) || !(tbl_nm %in% all_tables)) return(NULL)
    dbReadTable(con, tbl_nm) %>% select(SEQN, !!sym(var))
  })
  reduce(compact(dfs), left_join, by = "SEQN")
}

vars_F <- load_vars_cycle("F")
vars_G <- load_vars_cycle("G")

if (opt$otu_role == "indep") {
  dep_prefix   <- "dep_";   indep_prefix <- "indep_"
  dep_lookup_F <- var2table$F
  dep_lookup_G <- var2table$G
  indep_lookup_F <- setNames(rep(opt$otu_F, ncol(otu_F) - 1), setdiff(names(otu_F), "SEQN"))
  indep_lookup_G <- setNames(rep(opt$otu_G, ncol(otu_G) - 1), setdiff(names(otu_G), "SEQN"))
  dep_df_F   <- vars_F; dep_df_G   <- vars_G
  indep_df_F <- otu_F;  indep_df_G <- otu_G
} else {
  dep_prefix   <- "dep_";   indep_prefix <- "indep_"
  dep_lookup_F <- setNames(rep(opt$otu_F, ncol(otu_F) - 1), setdiff(names(otu_F), "SEQN"))
  dep_lookup_G <- setNames(rep(opt$otu_G, ncol(otu_G) - 1), setdiff(names(otu_G), "SEQN"))
  indep_lookup_F <- var2table$F
  indep_lookup_G <- var2table$G
  dep_df_F   <- otu_F; dep_df_G   <- otu_G
  indep_df_F <- vars_F; indep_df_G <- vars_G
}

ss_F <- compute_counts(dep_df_F, indep_df_F, dep_prefix, indep_prefix,
                      dep_lookup_F, indep_lookup_F) %>% mutate(cycle = "F")
ss_G <- compute_counts(dep_df_G, indep_df_G, dep_prefix, indep_prefix,
                      dep_lookup_G, indep_lookup_G) %>% mutate(cycle = "G")
ss   <- bind_rows(ss_F, ss_G)

outfile <- file.path(opt$out_dir,
                     glue("{opt$pipeline}_{opt$transform}_schema_structure.csv"))
write_csv(ss, outfile)
dbDisconnect(con)
