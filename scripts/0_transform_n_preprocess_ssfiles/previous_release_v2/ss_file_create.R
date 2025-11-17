#!/usr/bin/env Rscript
# -----------------------------------------------------------
#  Schema Structure File Creation
# Build a schema_structure ("ss_file") matrix for 
# transformed microbiome tables with proper naming convention
# -----------------------------------------------------------

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

# ---------- CLI ----------
opt_parser <- OptionParser(option_list = list(
  make_option("--db",        type="character",
              default="data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite",
              help="path to  transformed DB"),
  make_option("--otu_F",     type="character",
              help="cycle-F microbiome table"),
  make_option("--otu_G",     type="character",
              help="cycle-G microbiome table"),
  make_option("--vars_file", type="character",
              help="file listing NHANES variable names (one per line)"),
  make_option("--otu_role",  type="character", default="indep",
              help="indep | dep (is microbiome independent?)"),
  make_option("--pipeline",  type="character", default="3_exWAS"),
  make_option("--transform", type="character", default="none"),
  make_option("--out_dir",   type="character", default="results/0_ss_files")
))
opt <- parse_args(opt_parser)

stopifnot(opt$otu_role %in% c("indep","dep"))
dir_create(opt$out_dir, recurse = TRUE)

cat(" Schema Structure Creation\n")
cat("   Database:", opt$db, "\n")
cat("   OTU F:", opt$otu_F, "\n")
cat("   OTU G:", opt$otu_G, "\n")
cat("   Variables:", opt$vars_file, "\n")
cat("   Pipeline:", opt$pipeline, "\n")
cat("   Transform:", opt$transform, "\n")

# ---------- open DB ----------
con        <- dbConnect(SQLite(), opt$db)
all_tables <- dbListTables(con)

cat("   Available tables:", length(all_tables), "\n")

if (!"variable_names_epcf" %in% all_tables) {
  stop("Mapping table 'variable_names_epcf' not found in DB")
}

# Check if  OTU tables exist
if (!opt$otu_F %in% all_tables) {
  stop(" OTU table not found: ", opt$otu_F)
}
if (!opt$otu_G %in% all_tables) {
  stop(" OTU table not found: ", opt$otu_G)
}

# ---------- helper: map one var + cycle to its Data.File.Name entry ----------
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

# ---------- read OTU tables directly ----------
read_otu <- function(tbl) {
  if (!tbl %in% all_tables) stop("Missing OTU table: ", tbl)
  df <- dbReadTable(con, tbl)
  cat("   Loaded", tbl, ":", nrow(df), "samples ×", ncol(df)-1, "taxa\n")
  df
}

# ---------- compute pairwise schema_structure counts ----------
compute_counts <- function(dep_df, indep_df,
                           dep_prefix, indep_prefix,
                           dep_lookup, indep_lookup) {
  dep_vars   <- setdiff(names(dep_df), "SEQN")
  indep_vars <- setdiff(names(indep_df), "SEQN")
  merged     <- inner_join(dep_df, indep_df, by = "SEQN")
  
  cat("   Computing counts for", length(dep_vars), "×", length(indep_vars), "pairs...\n")
  
  map_dfr(dep_vars, function(dv) {
    dep_tbl <- dep_lookup[[dv]]
    d_vec   <- merged[[dv]]
    map_dfr(indep_vars, function(iv) {
      indep_tbl <- indep_lookup[[iv]]
      tibble(
        !!paste0(dep_prefix,  "var")   := dv,
        !!paste0(dep_prefix,  "table") := dep_tbl,
        !!paste0(indep_prefix, "var")   := iv,
        !!paste0(indep_prefix, "table") := indep_tbl,
        n = sum(!is.na(d_vec) & !is.na(merged[[iv]]))
      )
    })
  })
}

# ---------- load inputs ----------
varlist <- read_lines(opt$vars_file) %>%
           str_trim() %>%
           discard(~ .x == "")

cat("   Loading", length(varlist), "variables from", opt$vars_file, "\n")

otu_F <- read_otu(opt$otu_F)
otu_G <- read_otu(opt$otu_G)

# ---------- build var->table mapping for each cycle ----------
var2table <- list(
  F = setNames(map_chr(varlist, ~ get_var_table_for_cycle(.x, "F")), varlist),
  G = setNames(map_chr(varlist, ~ get_var_table_for_cycle(.x, "G")), varlist)
)

# check that every var appears in at least one cycle
missing_both <- varlist[is.na(var2table$F) & is.na(var2table$G)]
if (length(missing_both) > 0) {
  warning("Variables not found in any cycle: ",
       paste(missing_both, collapse = ", "))
  # Remove missing variables instead of stopping
  varlist <- setdiff(varlist, missing_both)
  var2table$F <- var2table$F[varlist]
  var2table$G <- var2table$G[varlist]
}

# ---------- assemble data frames of NHANES vars per cycle ----------
load_vars_cycle <- function(cycle) {
  dfs <- imap(var2table[[cycle]], function(tbl_nm, var) {
    if (is.na(tbl_nm) || !(tbl_nm %in% all_tables)) return(NULL)
    dbReadTable(con, tbl_nm) %>% select(SEQN, !!sym(var))
  })
  reduce(compact(dfs), left_join, by = "SEQN")
}

vars_F <- load_vars_cycle("F")
vars_G <- load_vars_cycle("G")

cat("   Loaded variables - F:", ncol(vars_F)-1, "vars,", nrow(vars_F), "samples\n")
cat("   Loaded variables - G:", ncol(vars_G)-1, "vars,", nrow(vars_G), "samples\n")

# ---------- decide roles and prefixes ----------
if (opt$otu_role == "indep") {
  dep_prefix   <- "dep_";   indep_prefix <- "indep_"
  dep_lookup_F <- var2table$F
  dep_lookup_G <- var2table$G
  indep_lookup_F <- setNames(rep(opt$otu_F, ncol(otu_F)-1), setdiff(names(otu_F), "SEQN"))
  indep_lookup_G <- setNames(rep(opt$otu_G, ncol(otu_G)-1), setdiff(names(otu_G), "SEQN"))
  dep_df_F     <- vars_F;    dep_df_G     <- vars_G
  indep_df_F   <- otu_F;     indep_df_G   <- otu_G
} else {
  dep_prefix   <- "dep_";   indep_prefix <- "indep_"
  dep_lookup_F <- setNames(rep(opt$otu_F, ncol(otu_F)-1), setdiff(names(otu_F), "SEQN"))
  dep_lookup_G <- setNames(rep(opt$otu_G, ncol(otu_G)-1), setdiff(names(otu_G), "SEQN"))
  indep_lookup_F <- var2table$F
  indep_lookup_G <- var2table$G
  dep_df_F     <- otu_F;     dep_df_G     <- otu_G
  indep_df_F   <- vars_F;    indep_df_G   <- vars_G
}

# ---------- compute schema_structure tables and bind cycles ----------
cat("   Computing schema structure for cycle F...\n")
ss_F <- compute_counts(
          dep_df_F,  indep_df_F,
          dep_prefix, indep_prefix,
          dep_lookup_F, indep_lookup_F
        ) %>% mutate(cycle = "F")

cat("   Computing schema structure for cycle G...\n")
ss_G <- compute_counts(
          dep_df_G,  indep_df_G,
          dep_prefix, indep_prefix,
          dep_lookup_G, indep_lookup_G
        ) %>% mutate(cycle = "G")

ss <- bind_rows(ss_F, ss_G)

# ---------- write out with  prefix ----------
outfile <- file.path(opt$out_dir,
                     glue("{opt$pipeline}_{opt$transform}_schema_structure.csv"))
write_csv(ss, outfile)

cat("✅  Schema Structure Complete\n")
cat("   Wrote", nrow(ss), "rows to", outfile, "\n")
dbDisconnect(con) 