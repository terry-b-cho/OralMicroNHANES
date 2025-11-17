#!/usr/bin/env Rscript
# -----------------------------------------------------------
# nhanes_db_filling_missing_data.R
# Assumes nhanes_oral_transformed.sqlite exists.
# Derives RIAGENDR_01 and fills missing rows in variable_names_epcf
# by mapping each new Variable→Table, using manual descriptions and
# borrowing table‐level fields from existing metadata.
# Usage:
#   Rscript scripts/0_transform_n_preprocess_ssfiles/nhanes_db_filling_missing_data.R \
#     --in_db  data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite \
#     --out_db data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(fs)
})

# ── CLI options ─────────────────────────────────────────────
opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-i","--in_db"), type="character",
    default="data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite",
    help="input transformed DB"),
  make_option(c("-o","--out_db"), type="character",
    default="data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite",
    help="output complete DB")
)))
message("Reading from: ", opt$in_db)
message("Writing to:  ", opt$out_db)

# ── copy input → output ─────────────────────────────────────
dir_create(path_dir(opt$out_db), recurse=TRUE)
file_copy(opt$in_db, opt$out_db, overwrite=TRUE)

con <- dbConnect(SQLite(), opt$out_db)

# ── CRITICAL FIX: Recreate SQL views for "_none" in complete DB ─────
# The SLURM jobs read the complete DB, but views are only in intermediate DB
rel_tables <- c(
  "DADA2RSV_GENUS_RELATIVE_F", "DADA2RSV_GENUS_RELATIVE_G",
  "DADA2RSV_FAMILY_RELATIVE_F","DADA2RSV_FAMILY_RELATIVE_G",
  "DADA2RSV_ORDER_RELATIVE_F", "DADA2RSV_ORDER_RELATIVE_G",
  "DADA2RSV_CLASS_RELATIVE_F", "DADA2RSV_CLASS_RELATIVE_G",
  "DADA2RSV_PHYLUM_RELATIVE_F","DADA2RSV_PHYLUM_RELATIVE_G"
)

message("Creating SQL views for '_none' tables in complete database...")
for (tbl in rel_tables) {
  if (tbl %in% dbListTables(con)) {
    view_name <- paste0(tbl, "_none")
    view_sql <- sprintf("CREATE VIEW IF NOT EXISTS %s AS SELECT * FROM %s", 
                       view_name, tbl)
    dbExecute(con, view_sql)
    message("   Created view: ", view_name, " -> ", tbl)
  }
}

# ── 1) DERIVE RIAGENDR_01 IN‐PLACE ───────────────────────────
for (series in c("F","G")) {
  tbln <- paste0("DEMO_", series)
  if (!tbln %in% dbListTables(con)) next
  cols <- dbListFields(con, tbln)
  if (!"RIAGENDR_01" %in% cols) {
    dbExecute(con,
      sprintf("ALTER TABLE %s ADD COLUMN RIAGENDR_01 INTEGER", tbln))
    dbExecute(con,
      sprintf("UPDATE %s SET RIAGENDR_01 = RIAGENDR - 1", tbln))
  }
}

# ── 2) MANUAL VARIABLE DESCRIPTIONS ─────────────────────────
# Only the vars we know live in tables but aren’t yet in metadata
annotations <- tribble(
  ~var_name,              ~var_description,
  "AGE_SQUARED",          "Age squared (RIDAGEYR²)",
  "BORN_INUSA",           "Born in USA (DMDBORN4=1)",
  "EDUCATION_LESS9",      "Edu: <9th grade (DMDEDUC2=1)",
  "EDUCATION_9_11",       "Edu: 9–11th grade (DMDEDUC2=2)",
  "EDUCATION_AA",         "Edu: AA/some college (DMDEDUC2=4)",
  "EDUCATION_COLLEGEGRAD","Edu: college graduate or above (DMDEDUC2=5)",
  "ETHNICITY_MEXICAN",    "Mexican American (RIDRETH1=1)",
  "ETHNICITY_OTHERHISPANIC","Other Hispanic (RIDRETH1=2)",
  "ETHNICITY_OTHER",      "Other race/multi-racial (RIDRETH1=5)",
  "ETHNICITY_NONHISPANICBLACK","Non-Hispanic Black (RIDRETH1=4)",
  "RIAGENDR_01",          "Gender (0=Male,1=Female)",
  "LBXP1",           "Total prostate specific antigen (ng/mL)",
  "LBXP2",           "Free prostate specific antigen (ng/mL)",
  "LBDP3",           "Prostate specific antigen ratio (%)",
  "LBDRFO",          "RBC folate (ng/mL)",
  "LBXSF6SI",        "Mefox oxidation product (nmol/L)",
  "LBXSCK",          "Creatine kinase (U/L)",
  "LBXFOLSI",        "Serum folate (nmol/L)",
  "LBXPS4",          "Complex prostate specific antigen (ng/mL)",
  "LBXTBM",          "TB Mitogen control result (Positive)",
  "BMXSUB",          "Subscapular Skinfold (mm)",
  "LBDFOT",          "Serum total folate (ng/mL)",
  "LBXMMASI",        "Mono-methyl arsenic (ng/mL)",
  "BMXTRI",          "Triceps Skinfold (mm)",
  "URX1DC",          "NAC-(1,2-dichlorovinyl)-L-cys (ng/mL)",
  "ASTHMA",          "Asthma (self-reported)",
  "BRONCHITIS",      "Bronchitis (self-reported)",
  "EMPHYSEMA",       "Emphysema (self-reported)",
  "ANGINA",          "Angina (self-reported)",
  "HEART_FAILURE",   "Heart failure (self-reported)",
  "HEART_ATTACK",    "Heart attack (self-reported)",
  "STROKE",          "Stroke (self-reported)",
  "CHD",             "Coronary heart disease (self-reported)",
  "CVD",             "Cardiovascular disease (self-reported)",
  "CANCER_BREAST",   "Breast cancer history",
  "CANCER_COLON",    "Colon cancer history",
  "CANCER_LUNG",     "Lung cancer history",
  "CANCER_ESOPHAGEAL","Esophageal cancer history",
  "CANCER_PROSTATE", "Prostate cancer history",
  "CANCER_MOUTH",    "Oral cancer history",
  "DIABETES",        "Diabetes (self-reported)"
)

# ── 3) TABLE‐LEVEL METADATA TEMPLATE ────────────────────────
template <- tbl(con, "variable_names_epcf") %>%
  select(
    Data.File.Name,
    Data.File.Description,
    Begin.Year,
    EndYear,
    Component,
    Use.Constraints
  ) %>%
  distinct() %>%
  collect()

# ── 4) BUILD VARIABLE→TABLE MAPPING ─────────────────────────
all_vars   <- annotations$var_name
all_tables <- dbListTables(con)
fields_map <- setNames(
  lapply(all_tables, function(tbl) dbListFields(con, tbl)),
  all_tables
)

mapping_long <- tibble(Variable = all_vars) %>%
  mutate(Data.File.Name = map(Variable, function(v) {
    names(fields_map)[
      vapply(fields_map, function(flds) v %in% flds, logical(1))
    ]
  })) %>%
  unnest_longer(Data.File.Name) %>%
  filter(!is.na(Data.File.Name)) %>%
  distinct(Variable, Data.File.Name)

# ── 5) FILTER OUT EXISTING METADATA ROWS ────────────────────
existing_meta <- tbl(con, "variable_names_epcf") %>%
  select(Variable.Name, Data.File.Name) %>%
  collect()

to_add <- anti_join(
  mapping_long,
  existing_meta,
  by = c("Variable"="Variable.Name", "Data.File.Name")
)

# ── 6) ASSEMBLE & APPEND NEW METADATA ROWS ──────────────────
if (nrow(to_add)>0) {
  new_meta <- to_add %>%
    left_join(annotations, by=c("Variable"="var_name")) %>%
    rename(Variable.Description = var_description) %>%
    left_join(template, by="Data.File.Name") %>%
    select(
      Variable.Name = Variable,
      Variable.Description,
      Use.Constraints,
      Data.File.Name,
      Data.File.Description,
      Begin.Year,
      EndYear,
      Component
    )

  dbBegin(con)
  tryCatch({
    old <- tbl(con, "variable_names_epcf") %>% collect()
    dbWriteTable(
      con, "variable_names_epcf",
      bind_rows(old, new_meta),
      overwrite=TRUE
    )
    dbCommit(con)
    message("Appended ", nrow(new_meta), " metadata rows.")
  }, error = function(e) {
    dbRollback(con)
    stop("Error appending metadata: ", e$message)
  })
} else {
  message("No new metadata to append.")
}

dbDisconnect(con)
message("✅ Completed: ", opt$out_db)
