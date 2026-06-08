#!/usr/bin/env Rscript

## -----------------------------------------------------------------------------
## Terry Plot Intermediate Data Preparation
## -----------------------------------------------------------------------------
## This script extracts every intermediate object required to reproduce the
## “Terry plots” from `scripts/4_association_phyloseq_analyses/4_association_phyloseq_analyses.Rmd`.
## **REVIEWER_EDIT:** All derived datasets are written ONLY to
## `results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/intermediate/`
## (isolated from published `4_otu_host_plot_out/`). Do not use the frozen
## `1_otu_host_and_indicies_intermediate_processing.R` for reviewer work.
## Figure typography (compact-layer genus axis/repel sizes) lives in Stage 2:
## `2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`.
## -----------------------------------------------------------------------------

set.seed(42)

## -----------------------------------------------------------------------------
## Helper: resolve repository base path from script location (base R only)
## -----------------------------------------------------------------------------
resolve_base_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
  if (length(script_path) == 0) {
    return(normalizePath("."))
  }
  normalizePath(file.path(dirname(script_path), "..", ".."))
}

base_path <- resolve_base_path()

## -----------------------------------------------------------------------------
## Directory setup — runs before library() so output tree exists even if
## package load fails (otherwise users never see 4_otu_host_plot_out_REVIEWER_EDIT/)
## -----------------------------------------------------------------------------
output_base_path <- file.path(base_path, "results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT")
viz_out_path <- file.path(output_base_path, "figures_out")
intermediate_files_path <- file.path(output_base_path, "intermediate")
supplementary_files_path <- file.path(output_base_path, "supplementary")

dir.create(viz_out_path, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_files_path, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_files_path, recursive = TRUE, showWarnings = FALSE)
for (m in c("compact", "comprehensive")) {
  for (t in c("none", "hellinger", "clr")) {
    dir.create(file.path(viz_out_path, m, t), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(viz_out_path, m, "indices", t), recursive = TRUE, showWarnings = FALSE)
  }
}

message("### Terry Plot Intermediate Processing (REVIEWER_EDIT) ###")
message("Base path: ", base_path)
message("Output root (created if missing): ", output_base_path)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(stringr)
  library(DBI)
  library(RSQLite)
  library(phyloseq)
})

## -----------------------------------------------------------------------------
## Paths to data sources
## -----------------------------------------------------------------------------
aggregated_association_res_path <- file.path(base_path, "results")
phyloseq_obj_files_path <- file.path(base_path, "results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate")
gold_db_ubiome_genus_mapping_path <- file.path(base_path, "results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv")
config_dir_path <- file.path(base_path, "configs")
db_primary <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
db_fallback <- file.path(base_path, "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite")
db_path <- db_primary
if (!file.exists(db_path) || file.info(db_path)$size == 0L) {
  message("Using fallback NHANES SQLite (processed DB missing or empty): ", db_fallback)
  db_path <- db_fallback
}
if (!file.exists(db_path) || file.info(db_path)$size == 0L) {
  stop("SQLite database not found or empty. Tried: ", db_primary, " and ", db_fallback)
}

## -----------------------------------------------------------------------------
## Database connection and metadata lookup
## -----------------------------------------------------------------------------
message("Connecting to SQLite database...")
con <- dbConnect(RSQLite::SQLite(), dbname = db_path)
on.exit({
  try(dbDisconnect(con), silent = TRUE)
}, add = TRUE)

message("Gathering table metadata...")
db_tables <- dbListTables(con)
years_F_or_G <- db_tables[grepl("_F$|_f$|_G$|_g$", db_tables)]
years_F_or_G <- c(
  years_F_or_G,
  "d_outcome_mcq_f", "d_outcome_mcq_g",
  "d_outcome_oh_f", "d_outcome_oh_g",
  "OralDisease_F", "OralDisease_G",
  "phenotype_vars_f",
  "table_names_epcf", "testosterone_vars_f", "testosterone_vars_g", "variable_names_epcf"
)

table_description <- tbl(con, "table_names_epcf") %>%
  collect() %>%
  as.data.table() %>%
  filter(Data.File.Name %in% years_F_or_G)

variable_description <- tbl(con, "variable_names_epcf") %>%
  collect() %>%
  as.data.table() %>%
  filter(Data.File.Name %in% years_F_or_G)

saveRDS(table_description, file.path(intermediate_files_path, "table_description.rds"))
saveRDS(variable_description, file.path(intermediate_files_path, "variable_description.rds"))

## -----------------------------------------------------------------------------
## Load GOLD mapping
## -----------------------------------------------------------------------------
message("Loading GOLD genus mapping...")
ubiome_genus_mapping_complete <- readr::read_csv(
  gold_db_ubiome_genus_mapping_path,
  show_col_types = FALSE
) %>%
  as.data.table()

saveRDS(ubiome_genus_mapping_complete, file.path(intermediate_files_path, "ubiome_genus_mapping_complete.rds"))

## -----------------------------------------------------------------------------
## Load phyloseq objects
## -----------------------------------------------------------------------------
message("Loading phyloseq objects...")
phyloseq_files <- c(
  "ubiome_counts" = "ubiome_counts.rds",
  "ubiome_relative" = "ubiome_relative.rds",
  "ubiome_relative_none" = "ubiome_relative_none.rds",
  "ubiome_relative_clr" = "ubiome_relative_clr.rds",
  "ubiome_relative_hellinger" = "ubiome_relative_hellinger.rds",
  "ubiome_relative_lognorm" = "ubiome_relative_lognorm.rds"
)

phyloseq_objects <- list()
for (name in names(phyloseq_files)) {
  file_path <- file.path(phyloseq_obj_files_path, phyloseq_files[[name]])
  if (!file.exists(file_path)) {
    stop("Missing phyloseq object: ", file_path)
  }
  phyloseq_objects[[name]] <- readRDS(file_path)
  message(sprintf(
    "  Loaded %s (samples: %d, taxa: %d)",
    name,
    nsamples(phyloseq_objects[[name]]),
    ntaxa(phyloseq_objects[[name]])
  ))
}

saveRDS(phyloseq_objects, file.path(intermediate_files_path, "phyloseq_objects.rds"))

## -----------------------------------------------------------------------------
## Load regression result tables (tidied, glanced, rsq)
## -----------------------------------------------------------------------------
message("Loading association result tables...")

analysis_specs <- list(
  demoWAS = list(prefix = "1_demoWAS", transformations = c("none", "clr", "hellinger")),
  oradWAS = list(prefix = "2_oradWAS", transformations = c("none", "clr", "hellinger")),
  exWAS   = list(prefix = "3_exWAS",   transformations = c("none", "clr", "hellinger")),
  pheWAS  = list(prefix = "4_pheWAS",  transformations = c("none", "clr", "hellinger")),
  outWAS  = list(prefix = "5_outWAS",  transformations = c("none", "clr", "hellinger"))
)

load_result <- function(path) {
  if (file.exists(path)) {
    readRDS(path)
  } else {
    tibble()
  }
}

tidied_results <- list()
glanced_results <- list()
rsq_results <- list()

for (analysis in names(analysis_specs)) {
  spec <- analysis_specs[[analysis]]
  for (trans in spec$transformations) {
    prefix <- spec$prefix
    dir_path <- file.path(aggregated_association_res_path, paste0(prefix, "_out"), paste0("result_", trans))
    key <- paste0(analysis, "_", trans)
    tidied_results[[key]] <- load_result(file.path(dir_path, paste0(prefix, "_", trans, "_tidied_complete.rds")))
    glanced_results[[key]] <- load_result(file.path(dir_path, paste0(prefix, "_", trans, "_glanced_complete.rds")))
    rsq_results[[key]] <- load_result(file.path(dir_path, paste0(prefix, "_", trans, "_rsq_complete.rds")))
  }
}

saveRDS(tidied_results, file.path(intermediate_files_path, "tidied_results_raw.rds"))
saveRDS(glanced_results, file.path(intermediate_files_path, "glanced_results_raw.rds"))
saveRDS(rsq_results, file.path(intermediate_files_path, "rsq_results_raw.rds"))

## -----------------------------------------------------------------------------
## Identify OTU column based on analysis type
## -----------------------------------------------------------------------------
message("Deriving OTU columns for each dataset...")

use_dependent <- c("demoWAS", "exWAS")
use_independent <- c("oradWAS", "pheWAS", "outWAS")

processed_tidied <- tidied_results

for (name in names(processed_tidied)) {
  dat <- processed_tidied[[name]]
  if (nrow(dat) == 0) {
    next
  }
  if (any(startsWith(name, use_dependent))) {
    dat <- dat %>%
      mutate(
        otu = dependent_var,
        otu_generic = sub("_relative$", "", dependent_var)
      )
  } else if (any(startsWith(name, use_independent))) {
    dat <- dat %>%
      mutate(
        otu = independent_var,
        otu_generic = sub("_relative$", "", independent_var)
      )
  } else {
    dat <- dat %>% mutate(otu = NA_character_, otu_generic = NA_character_)
  }
  processed_tidied[[name]] <- dat
}

saveRDS(processed_tidied, file.path(intermediate_files_path, "tidied_results_with_otu.rds"))

## -----------------------------------------------------------------------------
## Case-imbalance filtering for binary outcomes
## -----------------------------------------------------------------------------
message("Computing binary case-count filters...")

get_binary_case_counts <- function(var_list, base_name, con, genus_f, genus_g) {
  table_f <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_f" else paste0(base_name, "_F")
  table_g <- if (base_name == "d_outcome_mcq") "d_outcome_mcq_g" else paste0(base_name, "_G")
  db_tables <- dbListTables(con)
  if (!(table_f %in% db_tables && table_g %in% db_tables)) {
    warning(sprintf("Missing tables for %s: %s or %s", base_name, table_f, table_g))
    return(NULL)
  }
  df_f <- tbl(con, table_f) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  df_g <- tbl(con, table_g) %>% collect() %>% mutate(SEQN = as.character(SEQN))
  df_all <- bind_rows(
    df_f %>% filter(SEQN %in% genus_f$SEQN),
    df_g %>% filter(SEQN %in% genus_g$SEQN)
  )
  map_dfr(var_list, function(v) {
    if (v %in% names(df_all)) {
      tibble(var_name = v, cases_count = sum(df_all[[v]] %in% c(1, "1", TRUE), na.rm = TRUE))
    } else {
      tibble(var_name = v, cases_count = NA_integer_)
    }
  })
}

genus_f <- tbl(con, "DADA2RSV_GENUS_RELATIVE_F") %>% select(SEQN) %>% collect() %>% mutate(SEQN = as.character(SEQN))
genus_g <- tbl(con, "DADA2RSV_GENUS_RELATIVE_G") %>% select(SEQN) %>% collect() %>% mutate(SEQN = as.character(SEQN))

orad_vars <- readLines(file.path(config_dir_path, "2_oradWAS_vars.txt"), warn = FALSE)
outwas_vars <- readLines(file.path(config_dir_path, "5_outWAS_vars.txt"), warn = FALSE)

oradWAS_case_counts <- get_binary_case_counts(orad_vars, "OralDisease", con, genus_f, genus_g)
outWAS_case_counts  <- get_binary_case_counts(outwas_vars, "d_outcome_mcq", con, genus_f, genus_g)

binary_case_threshold <- 0.005
case_cutoff <- 9847 * binary_case_threshold

oradWAS_case_counts_pass_list <- oradWAS_case_counts %>% filter(cases_count > case_cutoff) %>% pull(var_name)
outWAS_case_counts_pass_list  <- outWAS_case_counts %>% filter(cases_count > case_cutoff) %>% pull(var_name)

case_filters <- list(
  oradWAS_case_counts = oradWAS_case_counts,
  outWAS_case_counts = outWAS_case_counts,
  oradWAS_case_counts_pass_list = oradWAS_case_counts_pass_list,
  outWAS_case_counts_pass_list = outWAS_case_counts_pass_list,
  binary_case_threshold = binary_case_threshold
)

saveRDS(case_filters, file.path(intermediate_files_path, "binary_case_filters.rds"))

## -----------------------------------------------------------------------------
## OTU prevalence filtering
## -----------------------------------------------------------------------------
message("Computing OTU prevalence filters...")

rsv_genus_relative <- bind_rows(
  tbl(con, "DADA2RSV_GENUS_RELATIVE_F") %>% collect(),
  tbl(con, "DADA2RSV_GENUS_RELATIVE_G") %>% collect()
)

otu_non_zero <- rsv_genus_relative %>%
  summarise(across(where(is.numeric), ~ sum(. > 0), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "otu", values_to = "non_zero_count")

otu_pass_prevalance_list <- otu_non_zero %>%
  filter(non_zero_count > 9847 * 0.01 & otu != "SEQN") %>%
  pull(otu)

prevalence_filters <- list(
  otu_non_zero = otu_non_zero,
  otu_pass_prevalance_list = otu_pass_prevalance_list,
  prevalence_threshold = 0.01
)

saveRDS(prevalence_filters, file.path(intermediate_files_path, "prevalence_filters.rds"))

## -----------------------------------------------------------------------------
## Filter significant association results
## -----------------------------------------------------------------------------
message("Filtering significant association results...")

was_specs <- c(
  "demoWAS_none", "demoWAS_clr", "demoWAS_hellinger",
  "oradWAS_none", "oradWAS_clr", "oradWAS_hellinger",
  "exWAS_none",  "exWAS_clr",  "exWAS_hellinger",
  "pheWAS_none", "pheWAS_clr", "pheWAS_hellinger",
  "outWAS_none", "outWAS_clr", "outWAS_hellinger"
)

significant_results <- list()

for (base in was_specs) {
  dat <- processed_tidied[[base]]
  if (is.null(dat) || nrow(dat) == 0) {
    next
  }
  filtered <- dat %>%
    filter(str_starts(term, independent_var)) %>%
    filter(std.error != 0) %>%
    mutate(dependent_var = map_chr(dependent_var, ~ as.character(.x)[1]))

  if (startsWith(base, "outWAS")) {
    filtered <- filtered %>% filter(dependent_var %in% outWAS_case_counts_pass_list)
  } else if (startsWith(base, "oradWAS")) {
    filtered <- filtered %>% filter(dependent_var %in% oradWAS_case_counts_pass_list)
  }

  # Combined significance: Benjamini–Hochberg FDR (p.value.fdr) and Storey q (q.value), both <= 0.05
  filtered <- filtered %>%
    filter(otu %in% otu_pass_prevalance_list) %>%
    filter(!is.na(p.value.fdr), p.value.fdr <= 0.05) %>%
    filter(!is.na(q.value), q.value <= 0.05) %>%
    arrange(p.value.fdr) %>%
    select(
      otu, term, dependent_var, independent_var, estimate, statistic,
      std.error, p.value, p.value.fdr, q.value, otu_generic,
      available_cycles, effect_scale, interpretation_note, fdr_corrected
    )

  significant_results[[paste0(base, "_sig_res")]] <- filtered
  message(sprintf("  %s: %d significant associations", base, nrow(filtered)))
}

saveRDS(significant_results, file.path(intermediate_files_path, "significant_results.rds"))

## -----------------------------------------------------------------------------
## Dataset list for Terry plots
## -----------------------------------------------------------------------------
datasets <- list(
  demowas_none      = significant_results[["demoWAS_none_sig_res"]],
  demowas_clr       = significant_results[["demoWAS_clr_sig_res"]],
  demowas_hellinger = significant_results[["demoWAS_hellinger_sig_res"]],
  oradwas_none      = significant_results[["oradWAS_none_sig_res"]],
  oradwas_clr       = significant_results[["oradWAS_clr_sig_res"]],
  oradwas_hellinger = significant_results[["oradWAS_hellinger_sig_res"]],
  exwas_none        = significant_results[["exWAS_none_sig_res"]],
  exwas_clr         = significant_results[["exWAS_clr_sig_res"]],
  exwas_hellinger   = significant_results[["exWAS_hellinger_sig_res"]],
  phewas_none       = significant_results[["pheWAS_none_sig_res"]],
  phewas_clr        = significant_results[["pheWAS_clr_sig_res"]],
  phewas_hellinger  = significant_results[["pheWAS_hellinger_sig_res"]],
  outwas_none       = significant_results[["outWAS_none_sig_res"]],
  outwas_clr        = significant_results[["outWAS_clr_sig_res"]],
  outwas_hellinger  = significant_results[["outWAS_hellinger_sig_res"]]
)

saveRDS(datasets, file.path(intermediate_files_path, "datasets_for_plots.rds"))

## -----------------------------------------------------------------------------
## Taxonomy annotations for plotting
## -----------------------------------------------------------------------------
phyloseq_reference <- phyloseq_objects[["ubiome_relative"]]

phylum_info <- data.frame(
  otu = rownames(tax_table(phyloseq_reference)),
  Phylum = as.character(tax_table(phyloseq_reference)[, "Phylum"]),
  stringsAsFactors = FALSE
)

genus_mapping <- data.frame(
  otu = rownames(tax_table(phyloseq_reference)),
  Genus = tax_table(phyloseq_reference)[, "Genus"],
  stringsAsFactors = FALSE
)

saveRDS(list(phylum_info = phylum_info, genus_mapping = genus_mapping),
        file.path(intermediate_files_path, "taxonomy_annotations.rds"))

## -----------------------------------------------------------------------------
## Variable annotation mapping
## -----------------------------------------------------------------------------
message("Creating variable annotation mapping...")

var_files <- file.path(config_dir_path, c(
  "1_demoWAS_vars.txt",
  "2_oradWAS_vars.txt",
  "3_exWAS_vars.txt",
  "4_pheWAS_vars.txt",
  "5_outWAS_vars.txt"
))

var_names <- unique(unlist(lapply(var_files, readLines)))

ubiome_variable_mapping <- data.frame(var_name = var_names, stringsAsFactors = FALSE) %>%
  left_join(
    variable_description %>%
      select(Variable.Name, Variable.Description) %>%
      distinct(Variable.Name, .keep_all = TRUE),
    by = c("var_name" = "Variable.Name")
  ) %>%
  rename(var_description = Variable.Description) %>%
  select(var_name, var_description) %>%
  distinct(var_name, .keep_all = TRUE) %>%
  as.data.table()

annotations <- data.table(
  var_name = c(
    "LBXP1", "LBXP2", "LBDP3", "LBDRFO", "LBXSF6SI", "LBXSCK", "LBXFOLSI",
    "LBXPS4", "LBXTBM", "BMXSUB", "LBDFOT", "LBXMMASI", "BMXTRI", "URX1DC",
    "RIAGENDR_01", "ASTHMA", "BRONCHITIS", "EMPHYSEMA", "ANGINA", "HEART_FAILURE", "HEART_ATTACK",
    "STROKE", "CHD", "CVD", "CANCER_BREAST", "CANCER_COLON", "CANCER_LUNG",
    "CANCER_ESOPHAGEAL", "CANCER_PROSTATE", "CANCER_MOUTH", "DIABETES",
    "AGE_SQUARED", "EDUCATION_LESS9", "EDUCATION_9_11", "EDUCATION_AA",
    "EDUCATION_COLLEGEGRAD", "ETHNICITY_MEXICAN", "ETHNICITY_OTHERHISPANIC",
    "ETHNICITY_OTHER", "ETHNICITY_NONHISPANICBLACK", "BORN_INUSA"
  ),
  var_description = c(
    "Total prostate specific antigen (ng/mL)",
    "Free prostate specific antigen (ng/mL)",
    "Prostate specific antigen ratio (%)",
    "RBC folate (ng/mL)",
    "Mefox oxidation product (nmol/L)",
    "Creatine kinase (U/L)",
    "Serum folate (nmol/L)",
    "Complex prostate specific antigen (ng/mL)",
    "TB Mitogen control result (Positive)",
    "Subscapular Skinfold (mm)",
    "Serum total folate (ng/mL)",
    "Mono-methyl arsenic (ng/mL)",
    "Triceps Skinfold (mm)",
    "NAC-(1,2-dichlorovinyl)-L-cys (ng/mL)",
    "Gender (Female)", "Asthma", "Bronchitis", "Emphysema", "Angina", "Heart failure", "Heart attack",
    "Stroke", "Coronary heart disease", "Cardiovascular disease", "Breast cancer",
    "Colon cancer", "Lung cancer", "Esophageal cancer", "Prostate cancer",
    "Oral cancer", "Diabetes", "Age²", "Edu: <9th", "Edu: 9–11th", "Edu: Associate",
    "Edu: College", "Mexican", "Other Hispanic", "Other ethnicity", "Black", "U.S. Born"
  )
)

missing_vars <- ubiome_variable_mapping[is.na(var_description), var_name]
annotations_needed <- annotations[var_name %in% missing_vars]

ubiome_variable_mapping[
  annotations_needed,
  on = "var_name",
  var_description := i.var_description
]

new_annotations <- annotations[!var_name %in% ubiome_variable_mapping$var_name]
if (nrow(new_annotations) > 0) {
  ubiome_variable_mapping <- rbindlist(
    list(ubiome_variable_mapping, new_annotations),
    use.names = TRUE,
    fill = TRUE
  )
}

ubiome_variable_mapping[, var_description := str_replace_all(
  var_description,
  c(
    " \\(?30 sec\\. pulse \\* 2\\):" = "",
    " concentration \\(" = " \\(",
    " percent \\(" = " \\(",
    " count \\(" = " \\(",
    " number \\(" = " \\(",
    "1000" = "10³",
    "million" = "10⁶"
  )
)]

shortenings <- data.table(
  var_name = c(
    "ENXMEAN", "LBDNENO", "LBXSASSI", "LBXTPO", "LBXTSH1", "URXMC1", "URX4FP",
    "URXMOH", "LBXVME", "URXUMMA", "URXUTM", "DR2TMFAT", "DR2TVARA", "URXCNP",
    "URXMX4", "DR2TPFAT", "URXAMU", "URXMHH", "URXMX5", "URXECP", "URXMX6",
    "URXTCC", "DR2TNUMF", "LBXV2T", "LBXVCM", "LBXVTC", "URXMHP", "URXMX7",
    "DR1TMFAT", "DR1TVARA", "LBXMPAH", "LBXVCT", "SMQ720", "URXBP3", "URXTRS",
    "DR1TPFAT", "DR2TATOA", "LBXGLT", "LBXTT3", "SPXNFEV1", "SPXNFVC", "LBXEPAH",
    "LBXVBM", "DR1TNUMF", "DR2TFDFE", "LBXVDB", "DS1DSCNT", "DR1TATOA", "DR1TSODI",
    "DR1TFDFE", "DS2DSCNT", "SMQ720", "LBXV4C", "DR2TSODI", "LBXV4E", "RIDAGEYR", "INDFMPIR"
  ),
  var_description = c(
    "Mean reproducible FENO (ppb)", "Segmented neutrophils (10³/µL)",
    "Aspartate aminotransferase (U/L)", "Thyroid peroxidase Ab (IU/mL)",
    "Thyroid stimulating hormone (µIU/mL)", "Mono-(3-carboxypropyl) phthalate",
    "4-fluoro-3-phenoxybenzoic acid", "Mono-(2-ethyl-5-oxohexyl) phthalate",
    "Blood methyl t-butyl ether (MTBE)", "Urine monomethylarsonic acid",
    "Urine trimethylarsine oxide", "Total monounsaturated FA (g)",
    "Vitamin A (retinol activity equiv.)", "Mono(carboxynonyl) phthalate",
    "1,3-dimethylxanthine (theophylline)", "Total polyunsaturated FA (g)",
    "5-acetylamino-6-amino-3-methyluracil", "Mono-2-ethyl-5-hydroxyhexyl phthalate",
    "1,7-dimethylxanthine (paraxanthine)", "Mono-2-ethyl-5-carboxypentyl phthalate",
    "3,7-dimethylxanthine (theobromine)", "trans-DCCA",
    "Total foods reported (foods file)", "Blood trans-1,2-Dichloroethene",
    "Blood dibromochloromethane (pg/mL)", "Blood trichloroethene (ng/mL)",
    "Mono-(2-ethyl)-hexyl phthalate", "1,3,7-trimethylxanthine (caffeine)",
    "Total monounsaturated FA (g)", "Vitamin A (retinol activity equiv.)",
    "2-(N-Methyl-PFOS) acetic acid", "Blood carbon tetrachloride (ng/mL)",
    "Cigarettes smoked past 5 days (daily)", "Urinary benzophenone-3",
    "Urinary triclosan", "Total polyunsaturated FA (g)",
    "Added α-tocopherol (Vit E, mg)", "2-hr oral glucose tolerance (mg/dL)",
    "Total triiodothyronine (T3, ng/dL)", "Spirometry baseline FEV₁ (mL)",
    "Spirometry baseline FVC (mL)", "2-(N-Ethyl-PFOS) acetic acid",
    "Blood bromodichloromethane (pg/mL)", "Total foods",
    "Dietary folate equivalents (mcg)", "Blood 1,4-Dichlorobenzene (ng/mL)",
    "Dietary supplements", "Added α-tocopherol (Vit E, mg)",
    "Sodium adjusted for salt (mg)", "Dietary folate equivalents (mcg)",
    "Dietary supplements", "Daily cigarettes past 5 days",
    "Blood tetrachloroethene (ng/mL)", "Sodium adjusted for salt (mg)",
    "Blood 1,1,2,2-Tetrachloroethane (ng/mL)", "Age at screening", "Poverty income ratio"
  )
)

long_vars <- ubiome_variable_mapping[nchar(var_description) > 35, var_name]
shortenings_needed <- shortenings[var_name %in% long_vars]

ubiome_variable_mapping[
  shortenings_needed,
  on = "var_name",
  var_description := i.var_description
]

ubiome_variable_mapping <- unique(ubiome_variable_mapping, by = "var_name")

saveRDS(ubiome_variable_mapping, file.path(intermediate_files_path, "ubiome_variable_mapping.rds"))
fwrite(ubiome_variable_mapping, file.path(intermediate_files_path, "ubiome_variable_mapping.csv"))

## -----------------------------------------------------------------------------
## Bucket definitions (exactly as in original Rmd)
## -----------------------------------------------------------------------------
message("Recording bucket definitions...")

LBD_antigens_antibody <- c(
  "LBDHBG", "LBDHCV", "LBDHD", "LBDHEG", "LBDHEM", "LBDRPCR", "LBDVWCGP"
)

energy_vars <- c("DR1TKCAL", "DR2TKCAL", "DS1TKCAL", "DS2TKCAL", "DSQTKCAL")

alcohol_methylxanthines <- c("DR1TALCO", "DR2TALCO", "DR1TCAFF", "DR2TCAFF", "DR1TTHEO", "DR2TTHEO")

macro_totals <- c(
  "DR1TCARB", "DR2TCARB", "DR1TPROT", "DR2TPROT", "DR1TTFAT", "DR2TTFAT",
  "DR1TSUGR", "DR2TSUGR", "DSQTSUGR", "DR1TFIBE", "DR2TFIBE", "DR1TCHOL",
  "DR2TCHOL", "DR1TMOIS", "DR2TMOIS"
)

fa_class_totals <- c("DR1TSFAT", "DR2TSFAT", "DR1TMFAT", "DR2TMFAT", "DR1TPFAT", "DR2TPFAT")

fa_sfa <- c(
  "DR1TS040", "DR2TS040", "DR1TS060", "DR2TS060", "DR1TS080", "DR2TS080",
  "DR1TS100", "DR2TS100", "DR1TS120", "DR2TS120", "DR1TS140", "DR2TS140",
  "DR1TS160", "DR2TS160", "DR1TS180", "DR2TS180"
)

fa_mufa <- c("DR1TM161", "DR2TM161", "DR1TM181", "DR2TM181", "DR1TM201", "DR2TM201", "DR1TM221", "DR2TM221")

fa_pufa <- c(
  "DR1TP182", "DR2TP182", "DR1TP183", "DR2TP183", "DR1TP184", "DR2TP184",
  "DR1TP204", "DR2TP204", "DR1TP205", "DR2TP205", "DR1TP225", "DR2TP225",
  "DR1TP226", "DR2TP226"
)

carotenoids <- c(
  "DR1TACAR", "DR2TACAR", "DR1TBCAR", "DR2TBCAR", "DR1TCRYP", "DR2TCRYP",
  "DR1TLZ", "DR2TLZ", "DR1TLYCO", "DR2TLYCO"
)

retinoids_tocopherols <- c(
  "DR1TRET", "DR2TRET", "DR1TVARA", "DR2TVARA", "DR1TATOC", "DR2TATOC",
  "DR1TATOA", "DR2TATOA"
)

b_complex <- c(
  "DR1TVB1", "DR2TVB1", "DS1TVB1", "DS2TVB1", "DSQTVB1",
  "DR1TVB2", "DR2TVB2", "DS1TVB2", "DS2TVB2", "DSQTVB2",
  "DR1TNIAC", "DR2TNIAC", "DS1TNIAC", "DS2TNIAC", "DSQTNIAC",
  "DR1TVB6", "DR2TVB6", "DS1TVB6", "DS2TVB6", "DSQTVB6",
  "DR1TVB12", "DR2TVB12", "DS1TVB12", "DS2TVB12", "DSQTVB12",
  "DR1TB12A", "DR2TB12A", "DR2TCHL"
)

folate_group <- unique(c(
  "DR1TFA", "DR2TFA", "DS1TFA", "DS2TFA", "DSQTFA",
  "DR1TFF", "DR2TFF",
  "DR1TFOLA", "DR2TFOLA",
  "DR1TFDFE", "DR2TFDFE", "DS1TFDFE", "DS2TFDFE", "DSQTFDFE"
))

vitamin_c_d_k <- c("DR1TVC", "DR2TVC", "DS1TVC", "DS2TVC", "DSQTVC", "DR1TVD", "DR2TVD", "DS1TVD", "DS2TVD", "DSQTVD", "DR1TVK", "DR2TVK", "DS1TVK", "DS2TVK", "DSQTVK")

minerals <- c(
  "DR1TCALC", "DR2TCALC", "DS1TCALC", "DS2TCALC", "DSQTCALC",
  "DR1TPHOS", "DR2TPHOS", "DSQTPHOS",
  "DR1TMAGN", "DR2TMAGN", "DS1TMAGN", "DS2TMAGN", "DSQTMAGN",
  "DR1TPOTA", "DR2TPOTA", "DS1TPOTA", "DS2TPOTA", "DSQTPOTA",
  "DR1TSODI", "DR2TSODI",
  "DR1TIRON", "DR2TIRON", "DS1TIRON", "DS2TIRON", "DSQTIRON",
  "DR1TZINC", "DR2TZINC", "DS1TZINC", "DS2TZINC", "DSQTZINC",
  "DR1TCOPP", "DR2TCOPP", "DS1TCOPP", "DS2TCOPP", "DSQTCOPP",
  "DR1TSELE", "DR2TSELE", "DS1TSELE", "DS2TSELE", "DSQTSELE",
  "DS2TIODI", "DSQTIODI"
)

intake_counts <- c("DR1TNUMF", "DR2TNUMF", "DS1DSCNT", "DS2DSCNT")

supplement_flags <- c("DS2AN", "DS2ANCNT", "DS1DSCNT", "DS2DS", "DS2DSCNT")

macro_energy_vars       <- c(energy_vars, macro_totals)
fatty_acid_vars         <- c(fa_class_totals, fa_sfa, fa_mufa, fa_pufa)
vitamin_vars            <- c(carotenoids, retinoids_tocopherols, b_complex, folate_group, vitamin_c_d_k)
mineral_vars            <- minerals
alcohol_stims_vars      <- alcohol_methylxanthines
intake_count_vars       <- intake_counts
supplement_indicator_vars <- supplement_flags

chlorinated_phenols_herbicides <- c("URX14D", "URXDCB", "URX1TB", "URX3TB", "URX24D", "URX25T")

phenolic_edcs <- c(
  "URXBP3", "URXBPH", "URX4TO", "URXOPP", "URXTRS", "URXBUP", "URXEPB", "URXMPB", "URXPPB"
)

phthalate_metabolites <- c(
  "URXCNP", "URXCOP", "URXECP", "URXMBP", "URXMC1", "URXMCP", "URXMEP", "URXMHH",
  "URXMHP", "URXMIB", "URXMNM", "URXMNP", "URXMOH", "URXMOP", "URXMZP"
)

pesticide_metabolites <- c(
  "URX4FP", "URXCB3", "URXOPM", "URXTCC", "URXCPM", "URXMAL", "URXPAR", "URXOP1", "URXOP2",
  "URXOP3", "URXOP4", "URXOP5", "URXOP6", "URXDEE", "URXDEA", "URXDHD"
)

pah_metabolites <- c("URXP01", "URXP02", "URXP03", "URXP04", "URXP05", "URXP06", "URXP07", "URXP10", "URXP17")

phytoestrogens_lignans <- c("URXDAZ", "URXEQU", "URXGNS", "URXDMA", "URXETD", "URXETL")

methylxanthines_metabolites <- c(
  "URXAMU", "URXMU1", "URXMU2", "URXMU3", "URXMU4", "URXMU5", "URXMU6", "URXMU7",
  "URXMX1", "URXMX2", "URXMX3", "URXMX4", "URXMX5", "URXMX6", "URXMX7", "URXOXY"
)

tobacco_biomarkers <- c("URXNAL", "URXSCN")

inorganic_anions <- c("URXNO3", "URXUP8")

urinary_metals <- c(
  "URXUAB", "URXUAC", "URXUAS", "URXUAS3", "URXUAS5", "URXUDMA", "URXUMMA", "URXUTM",
  "URXUBA", "URXUBE", "URXUCD", "URXUCO", "URXUCS", "URXUMO", "URXUMN", "URXUPB",
  "URXUPT", "URXUSB", "URXUTL", "URXUTU", "URXUSN", "URXUSR", "URXUHG"
)

voc_solvent_metabolites <- c(
  "URX1DC", "URX2DC", "URX2MH", "URX34M", "URXAAM", "URXAMC", "URXATC", "URXBMA", "URXBPM",
  "URXCEM", "URXCYM", "URXDHB", "URXDPM", "URXGAM", "URXHEM", "URXHP2", "URXHPM", "URXMAD",
  "URXMB1", "URXMB2", "URXMB3", "URXMHNC", "URXPHE", "URXPHG", "URXPMA", "URXPMM", "URXTCV",
  "URXTTC"
)

urinary_chlamydia <- c("URXUCL")

urine_pesticide_herbicide <- c(chlorinated_phenols_herbicides, pesticide_metabolites)
urine_consumer_edcs <- c(phenolic_edcs, phthalate_metabolites)
urine_pah <- pah_metabolites
urine_bioactive_metabolites <- c(phytoestrogens_lignans, methylxanthines_metabolites)
urine_smoke_anion_markers <- c(tobacco_biomarkers, inorganic_anions)
urine_metals <- urinary_metals
urine_infection_marker <- urinary_chlamydia
urine_voc_solvent_metabolites <- voc_solvent_metabolites

hpv_dna <- c("LBX06", "LBX11", "LBX16", "LBX18")

infectious_serology <- c("LBXHA", "LBXHBC", "LBXHBS", "LBXHE1", "LBXHE2", "LBXMEA", "LBXMUM", "LBXRUB")

tb_igra_control <- c("LBXTBN")
autoimmune_serology <- c("LBXTTG")

vitamin_b6_b12 <- c("LBX4PA", "LBXPLP", "LBXB12", "LBDB12SI")
vitamin_d <- c("LBXVD2MS", "LBXVD3MS", "LBXVE3MS", "LBXVIDMS")

fatty_acids <- c(
  "LBXALN", "LBXAR1", "LBXARA", "LBXCAP", "LBXDA1", "LBXDHA", "LBXDP3", "LBXDP6",
  "LBXDTA", "LBXED1", "LBXEN1", "LBXEPA", "LBXET1", "LBXGLA", "LBXHDT", "LBXHGL",
  "LBXLAR", "LBXLG1", "LBXLNA", "LBXML1", "LBXMR1", "LBXMRG", "LBXNR1", "LBXOL1",
  "LBXPEN", "LBXPL1", "LBXPM1", "LBXSD1", "LBXST1", "LBXTSA", "LBXVC1", "LBXOD1",
  "LBXOD9", "LBXOTT"
)

smoke_biomarkers <- c("LBXCOT", "LBX2DF", "LBXNM", "LBXVFN")

toxic_metals <- c("LBXBCD", "LBXBPB", "LBXTHG", "LBXIHG", "LBXBGE", "LBXBGM")

essential_trace_metals <- c("LBXBMN", "LBXBSE")

poly_halogenated_aliphatic <- c(
  "LBXV4E", "LBXVCT", "LBXV4C", "LBXVTC", "LBXVTE", "LBXV2E", "LBXVTP", "LBXVBF",
  "LBXVBM", "LBXVCM", "LBXVCF", "LBXV2P", "LBXVHE",
  "LBXPFBS", "LBXPFDE", "LBXPFDO", "LBXPFHP", "LBXPFHS", "LBXPFNA", "LBXPFOA", "LBXPFOS",
  "LBXPFSA", "LBXPFUA", "LBXEPAH", "LBXMPAH"
)

mono_di_halogenated_aliphatic <- c(
  "LBXV1A", "LBXV1E", "LBXV2A", "LBXV2C", "LBXV2T", "LBXVDE", "LBXVDM", "LBXVDP"
)

halogenated_aromatic <- c("LBXV1D", "LBXV3B", "LBXVDB", "LBXVCB")

aromatic_hydrocarbons <- c("LBXVBZ", "LBXVEB", "LBXVIPB", "LBXVTO", "LBXVOX", "LBXVXY", "LBXVST")

non_halogenated_aliphatics <- c("LBXV06", "LBXVME", "LBXVDX")

nitro_aromatic_voc <- c("LBXVNB")

hematology_indices <- c(
  "LBXBAPCT", "LBXEOPCT", "LBXHCT", "LBXHGB", "LBXLYPCT", "LBXMC", "LBXMCHSI", "LBXMCVSI",
  "LBXMOPCT", "LBXMPSI", "LBXNEPCT", "LBXRBCSI", "LBXRDW", "LBXPLTSI", "LBXWBCSI"
)

rbc_mass_indices <- c("LBXRBCSI", "LBXHGB", "LBXHCT")
rbc_morphology_indices <- c("LBXMC", "LBXMCHSI", "LBXMCVSI", "LBXRDW")
wbc_differential_indices <- c("LBXWBCSI", "LBXBAPCT", "LBXEOPCT", "LBXLYPCT", "LBXMOPCT", "LBXNEPCT")
platelet_count <- c("LBXPLTSI")
platelet_volume <- c("LBXMPSI")

iron_metabolism_markers <- c("LBXFER", "LBXSIR", "LBXTFR")

lipid_markers <- c("LBXAPB", "LBXSCH", "LBXTC", "LBXSTR", "LBXTR")

hepatic_muscle_enzymes <- c(
  "LBXSAL", "LBXSGB", "LBXSTP", "LBXSAPSI", "LBXSASSI", "LBXSATSI",
  "LBXSGTSI", "LBXSLDSI", "LBXSTB", "LBXSCK"
)

renal_electrolyte_markers <- c(
  "LBXSBU", "LBXSCR", "LBXSKSI", "LBXSNASI", "LBXSCLSI", "LBXSC3SI",
  "LBXSCA", "LBXSPH", "LBXSOSSI", "LBXSUA"
)

glucose_insulin_markers <- c("LBXGH", "LBXGLU", "LBXGLT", "LBXSGL", "LBXIN")

thyroid_markers <- c("LBXT3F", "LBXT4F", "LBXTT3", "LBXTT4", "LBXTSH1", "LBXTGN", "LBXATG", "LBXTPO")

folate_onecarbon_markers <- c("LBXFOLSI", "LBXRBFSI", "LBXMMASI", "LBXSF1SI", "LBXSF2SI", "LBXSF3SI", "LBXSF4SI", "LBXSF5SI", "LBXSF6SI")

inflammation_immune_markers <- c("LBXCRP", "LBXTBA", "LBXTBM")

prostate_markers <- c("LBXP1", "LBXP2", "LBXPS4")

reproductive_hormone_markers <- c("LBXTST")

blood_hematology_vars <- c(hematology_indices, iron_metabolism_markers)
blood_metabolic_vars <- c(renal_electrolyte_markers, hepatic_muscle_enzymes, lipid_markers, glucose_insulin_markers)
blood_endocrine_vars <- c(thyroid_markers, reproductive_hormone_markers)
blood_onecarbon_folate_vars <- folate_onecarbon_markers
blood_immune_inflammation_vars <- inflammation_immune_markers
blood_oncology_prostate_vars <- prostate_markers

blood_cell_count_vars <- c("LBDBANO", "LBDEONO", "LBDLYMNO", "LBDMONO", "LBDNENO")

serum_lipids <- c("LBDHDD", "LBDLDL")
prostate_marker <- c("LBDP3")
folate_status <- c("LBDFOT", "LBDRFO")
urine_kidney_markers <- c("URXUCR", "URXUMA")
urine_iodine <- c("URXUIO")
exhaled_nitric_oxide <- c("ENXMEAN")

serum_clinical_vars <- c(serum_lipids, prostate_marker)
folate_status_vars <- folate_status
urine_kidney_iodine_vars <- c(urine_kidney_markers, urine_iodine)

airway_inflammation_var <- c(exhaled_nitric_oxide, "SPXNFEV1", "SPXNFVC")

bucket_definitions <- list(
  LBD_antigens_antibody = LBD_antigens_antibody,
  energy_vars = energy_vars,
  alcohol_methylxanthines = alcohol_methylxanthines,
  macro_totals = macro_totals,
  fa_class_totals = fa_class_totals,
  fa_sfa = fa_sfa,
  fa_mufa = fa_mufa,
  fa_pufa = fa_pufa,
  carotenoids = carotenoids,
  retinoids_tocopherols = retinoids_tocopherols,
  b_complex = b_complex,
  folate_group = folate_group,
  vitamin_c_d_k = vitamin_c_d_k,
  minerals = minerals,
  intake_counts = intake_counts,
  supplement_flags = supplement_flags,
  chlorinated_phenols_herbicides = chlorinated_phenols_herbicides,
  phenolic_edcs = phenolic_edcs,
  phthalate_metabolites = phthalate_metabolites,
  pesticide_metabolites = pesticide_metabolites,
  pah_metabolites = pah_metabolites,
  phytoestrogens_lignans = phytoestrogens_lignans,
  methylxanthines_metabolites = methylxanthines_metabolites,
  tobacco_biomarkers = tobacco_biomarkers,
  inorganic_anions = inorganic_anions,
  urinary_metals = urinary_metals,
  voc_solvent_metabolites = voc_solvent_metabolites,
  urinary_chlamydia = urinary_chlamydia,
  urine_pesticide_herbicide = urine_pesticide_herbicide,
  urine_consumer_edcs = urine_consumer_edcs,
  urine_pah = urine_pah,
  urine_bioactive_metabolites = urine_bioactive_metabolites,
  urine_smoke_anion_markers = urine_smoke_anion_markers,
  urine_metals = urine_metals,
  urine_infection_marker = urine_infection_marker,
  urine_voc_solvent_metabolites = urine_voc_solvent_metabolites,
  hpv_dna = hpv_dna,
  infectious_serology = infectious_serology,
  tb_igra_control = tb_igra_control,
  autoimmune_serology = autoimmune_serology,
  vitamin_b6_b12 = vitamin_b6_b12,
  vitamin_d = vitamin_d,
  fatty_acids = fatty_acids,
  smoke_biomarkers = smoke_biomarkers,
  toxic_metals = toxic_metals,
  essential_trace_metals = essential_trace_metals,
  poly_halogenated_aliphatic = poly_halogenated_aliphatic,
  mono_di_halogenated_aliphatic = mono_di_halogenated_aliphatic,
  halogenated_aromatic = halogenated_aromatic,
  aromatic_hydrocarbons = aromatic_hydrocarbons,
  non_halogenated_aliphatics = non_halogenated_aliphatics,
  nitro_aromatic_voc = nitro_aromatic_voc,
  hematology_indices = hematology_indices,
  rbc_mass_indices = rbc_mass_indices,
  rbc_morphology_indices = rbc_morphology_indices,
  wbc_differential_indices = wbc_differential_indices,
  platelet_count = platelet_count,
  platelet_volume = platelet_volume,
  iron_metabolism_markers = iron_metabolism_markers,
  lipid_markers = lipid_markers,
  hepatic_muscle_enzymes = hepatic_muscle_enzymes,
  renal_electrolyte_markers = renal_electrolyte_markers,
  glucose_insulin_markers = glucose_insulin_markers,
  thyroid_markers = thyroid_markers,
  folate_onecarbon_markers = folate_onecarbon_markers,
  inflammation_immune_markers = inflammation_immune_markers,
  prostate_markers = prostate_markers,
  reproductive_hormone_markers = reproductive_hormone_markers,
  blood_cell_count_vars = blood_cell_count_vars,
  serum_lipids = serum_lipids,
  prostate_marker = prostate_marker,
  folate_status = folate_status,
  urine_kidney_markers = urine_kidney_markers,
  urine_iodine = urine_iodine,
  exhaled_nitric_oxide = exhaled_nitric_oxide,
  airway_inflammation_var = airway_inflammation_var
)

saveRDS(bucket_definitions, file.path(intermediate_files_path, "bucket_definitions.rds"))

## -----------------------------------------------------------------------------
## Extract Terry plot configuration lists directly from the R Markdown
## -----------------------------------------------------------------------------
message("Extracting Terry plot configurations from R Markdown...")

rmd_path <- file.path(base_path, "scripts/4_association_phyloseq_analyses/4_association_phyloseq_analyses.Rmd")
if (!file.exists(rmd_path)) {
  stop("Expected R Markdown file not found at ", rmd_path)
}

rmd_lines <- readLines(rmd_path)

extract_list_block <- function(lines, start_idx) {
  buffer <- character()
  idx <- start_idx
  repeat {
    buffer <- c(buffer, lines[idx])
    if (trimws(lines[idx]) == ")") break
    idx <- idx + 1
    if (idx > length(lines)) {
      stop("Unterminated list encountered when parsing Rmd (start index ", start_idx, ")")
    }
  }
  paste(buffer, collapse = "\n")
}

plot_config_starts <- grep("^plot_configs <- list\\(", rmd_lines)
if (length(plot_config_starts) < 3) {
  stop("Unable to locate all plot configuration lists in the Rmd.")
}

plot_configs_none <- NULL
plot_configs_hellinger <- NULL
plot_configs_clr <- NULL

plot_config_snippets <- lapply(plot_config_starts[1:3], function(idx) extract_list_block(rmd_lines, idx))

for (i in seq_along(plot_config_snippets)) {
  eval(parse(text = plot_config_snippets[[i]]))
  if (i == 1) plot_configs_none <- plot_configs
  if (i == 2) plot_configs_hellinger <- plot_configs
  if (i == 3) plot_configs_clr <- plot_configs
  rm(plot_configs)
}

indices_config_starts <- grep("^indices_plot_configs <- list\\(", rmd_lines)
if (length(indices_config_starts) < 3) {
  stop("Unable to locate all indices plot configuration lists in the Rmd.")
}

indices_plot_configs_none <- NULL
indices_plot_configs_hellinger <- NULL
indices_plot_configs_clr <- NULL

indices_config_snippets <- lapply(indices_config_starts[1:3], function(idx) extract_list_block(rmd_lines, idx))

for (i in seq_along(indices_config_snippets)) {
  eval(parse(text = indices_config_snippets[[i]]))
  if (i == 1) indices_plot_configs_none <- indices_plot_configs
  if (i == 2) indices_plot_configs_hellinger <- indices_plot_configs
  if (i == 3) indices_plot_configs_clr <- indices_plot_configs
  rm(indices_plot_configs)
}

saveRDS(plot_configs_none, file.path(intermediate_files_path, "plot_configs_none.rds"))
saveRDS(plot_configs_hellinger, file.path(intermediate_files_path, "plot_configs_hellinger.rds"))
saveRDS(plot_configs_clr, file.path(intermediate_files_path, "plot_configs_clr.rds"))

saveRDS(indices_plot_configs_none, file.path(intermediate_files_path, "indices_plot_configs_none.rds"))
saveRDS(indices_plot_configs_hellinger, file.path(intermediate_files_path, "indices_plot_configs_hellinger.rds"))
saveRDS(indices_plot_configs_clr, file.path(intermediate_files_path, "indices_plot_configs_clr.rds"))

## -----------------------------------------------------------------------------
## Final message
## -----------------------------------------------------------------------------
message("Intermediate processing complete. Objects saved to ", intermediate_files_path)


