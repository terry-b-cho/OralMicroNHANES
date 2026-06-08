################################################################################
# gold_db_process_n_map2genus.R
#
# Aggregates the GOLD microbial phenotype database by genus, parses
# phyloseq genus names, and produces the genus-to-GOLD-feature mapping
# used by downstream NHANES oral microbiome analyses.
#
# Inputs (paths relative to PROJECT_ROOT):
#   - data/00_GOLDdb/goldData_ubiome_selected.csv
#   - data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite
#   - results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_relative.rds
#
# Outputs (under results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/):
#   - gold_db_genus.csv
#   - ubiome_genus_mapping_complete.csv
#   - mapping_summary_stats.csv
#
# Environment: R >= 4.5 with data.table, dplyr, dbplyr, RSQLite, readr, tidyr,
# tibble, purrr, phyloseq. Exact versions in module3_tool_version_list.txt.
#
# Run:
#   Rscript scripts/3_gold_db_microbial_phenotype/gold_db_process_n_map2genus.R
################################################################################

# === USER SETTING =============================================================
# Set PROJECT_ROOT to the absolute path of your local clone of this repository.
PROJECT_ROOT <- "/n/groups/patel/terry/nhanes_oral_mirco_cho"
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(dbplyr)
  library(RSQLite)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(phyloseq)
})

# ---- paths ----
intermediate_files_path <- file.path(PROJECT_ROOT,
  "results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate")
phyloseq_obj_files_path <- file.path(PROJECT_ROOT,
  "results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate")
db_path <- file.path(PROJECT_ROOT,
  "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
gold_db_subset_path <- file.path(PROJECT_ROOT,
  "data/00_GOLDdb/goldData_ubiome_selected.csv")

dir.create(intermediate_files_path, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(db_path),
          file.exists(gold_db_subset_path),
          file.exists(file.path(phyloseq_obj_files_path, "ubiome_relative.rds")))

con <- dbConnect(RSQLite::SQLite(), dbname = db_path)

# ---- Step 1: load and clean GOLD database ----
gold_csv <- fread(gold_db_subset_path,
                  colClasses = "character", encoding = "UTF-8")

gold_csv[] <- lapply(gold_csv, function(x) {
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x[tolower(x) %in% c("na", "unclassified", "undefined", "unknown", "")] <- NA_character_
  x
})

names(gold_csv) <- iconv(names(gold_csv), from = "latin1", to = "UTF-8")
names(gold_csv) <- gsub("<a1>", "°", names(gold_csv))
names(gold_csv) <- gsub("^ORGANISM ", "", names(gold_csv))
names(gold_csv) <- gsub(" ", "_", names(gold_csv))

required_columns <- c("NCBI_KINGDOM", "GOLD_PHYLUM", "GENUS", "SPECIES", "BIOTIC_RELATIONSHIPS",
                      "OXYGEN_REQUIREMENT", "METABOLISM", "ENERGY_SOURCES", "GRAM_STAIN",
                      "SAMPLE_COLLECTION_DATE", "SALINITY", "CELL_SHAPE", "MOTILITY", "SPORULATION",
                      "TEMPERATURE_RANGE")
existing_columns <- intersect(required_columns, names(gold_csv))
gold_csv <- gold_csv[, ..existing_columns]

# ---- Step 2: categorical cleanup + scores ----
gold_csv[, OXYGEN_REQUIREMENT := fcase(
  OXYGEN_REQUIREMENT %in% c("Anaerobe", "Obligate anaerobe"), "Anaerobe",
  OXYGEN_REQUIREMENT %in% c("Aerobe", "Obligate aerobe"), "Aerobe",
  OXYGEN_REQUIREMENT %in% c("Facultative", "Facultative anaerobe", "Facultative aerobe"), "Facultative",
  OXYGEN_REQUIREMENT %in% c("Microaerophilic"), "Microaerophilic",
  default = NA_character_
)]
gold_csv[, OXYGEN_SCORE := fcase(
  OXYGEN_REQUIREMENT == "Anaerobe", 0,
  OXYGEN_REQUIREMENT == "Facultative", 0.5,
  OXYGEN_REQUIREMENT == "Microaerophilic", 0.25,
  OXYGEN_REQUIREMENT == "Aerobe", 1,
  default = NA_real_
)]

gold_csv[, GRAM_STAIN := fcase(
  GRAM_STAIN == "Gram-", "Gram-",
  GRAM_STAIN == "Gram+", "Gram+",
  default = NA_character_
)]
gold_csv[, STAIN_SCORE := fcase(
  GRAM_STAIN == "Gram-", 0,
  GRAM_STAIN == "Gram+", 1,
  default = NA_real_
)]

gold_csv[, MOTILITY := fcase(
  MOTILITY == "Nonmotile", "Nonmotile",
  MOTILITY %in% c("Motile", "Chemotactic"), "Motile",
  default = NA_character_
)]
gold_csv[, MOTILITY_SCORE := fcase(
  MOTILITY == "Nonmotile", 0,
  MOTILITY == "Motile", 1,
  default = NA_real_
)]

gold_csv[, SPORULATION := fcase(
  SPORULATION == "Nonsporulating", "Nonsporulating",
  SPORULATION == "Sporulating", "Sporulating",
  default = NA_character_
)]
gold_csv[, SPORULATION_SCORE := fcase(
  SPORULATION == "Nonsporulating", 0,
  SPORULATION == "Sporulating", 1,
  default = NA_real_
)]

# ---- Step 3: aggregation helpers ----
skip_cols       <- c("OXYGEN_SCORE", "STAIN_SCORE", "MOTILITY_SCORE", "SPORULATION_SCORE")
text_cols       <- c("METABOLISM", "ENERGY_SOURCES")
char_cols       <- names(gold_csv)[sapply(gold_csv, is.character)]
categorical_cols<- setdiff(char_cols, c(text_cols, skip_cols))
features        <- setdiff(names(gold_csv), c(skip_cols, "GENUS"))

feature_types <- sapply(features, function(f) {
  if (f %in% text_cols)        "text"
  else if (f %in% categorical_cols) "categorical"
  else                         "skip"
})

summarize_categorical <- function(values, threshold = 0.7) {
  total_count  <- length(values)
  valid_values <- na.omit(values)
  Non_NA_count <- as.integer(length(valid_values))
  if (Non_NA_count == 0) {
    return(list(aggregation_value = "unknown (0%)", agreement = "none",
                Non_NA_count = 0L, NA_count = total_count))
  }
  freq <- table(valid_values) / Non_NA_count
  aggregation_value <- paste(names(freq), " (", round(freq * 100, 1), "%)",
                             sep = "", collapse = ", ")
  agreement <- if (length(unique(valid_values)) == 1) "high"
               else if (max(freq) >= threshold) "medium"
               else "low"
  list(aggregation_value = aggregation_value, agreement = agreement,
       Non_NA_count = Non_NA_count, NA_count = total_count - Non_NA_count)
}

summarize_text <- function(values, n_top = 3) {
  total_count  <- length(values)
  valid_values <- na.omit(values)
  Non_NA_count <- as.integer(length(valid_values))
  if (Non_NA_count == 0) {
    return(list(aggregation_value = "unknown (0%)",
                Non_NA_count = 0L, NA_count = total_count))
  }
  tokens <- trimws(unlist(strsplit(valid_values, "[|,;]")))
  freq <- table(tokens) / length(tokens)
  top_terms <- head(sort(freq, decreasing = TRUE), n_top)
  aggregation_value <- paste(names(top_terms), " (", round(top_terms * 100, 1), "%)",
                             sep = "", collapse = ", ")
  list(aggregation_value = aggregation_value,
       Non_NA_count = Non_NA_count, NA_count = total_count - Non_NA_count)
}

summarize_index <- function(scores) {
  total_count  <- length(scores)
  scores       <- na.omit(scores)
  Non_NA_count <- as.integer(length(scores))
  if (Non_NA_count == 0) {
    return(list(aggregation_value = "unknown (0%)",
                Non_NA_count = 0L, NA_count = total_count))
  }
  list(aggregation_value = mean(scores),
       Non_NA_count = Non_NA_count, NA_count = total_count - Non_NA_count)
}

summarize_feature <- function(data, feature, type) {
  if (type == "skip") return(NULL)
  values <- data[[feature]]
  if (type == "categorical") {
    s <- summarize_categorical(values)
    return(data.table(feature = feature, aggregation_value = s$aggregation_value,
                      agreement = s$agreement,
                      Non_NA_count = s$Non_NA_count, NA_count = s$NA_count))
  }
  if (type == "text") {
    s <- summarize_text(values)
    return(data.table(feature = feature, aggregation_value = s$aggregation_value,
                      agreement = "medium",
                      Non_NA_count = s$Non_NA_count, NA_count = s$NA_count))
  }
  NULL
}

# ---- Step 4: aggregate by genus ----
stopifnot("GENUS" %in% names(gold_csv))

gold_db_genus <- gold_csv[, {
  feature_summaries <- lapply(features, function(f)
    summarize_feature(.SD, f, feature_types[f]))
  indices <- list(
    data.table(feature = "STAIN_INDEX",
               aggregation_value = summarize_index(.SD$STAIN_SCORE)$aggregation_value,
               agreement = "",
               Non_NA_count = summarize_index(.SD$STAIN_SCORE)$Non_NA_count,
               NA_count     = summarize_index(.SD$STAIN_SCORE)$NA_count),
    data.table(feature = "OXYGEN_INDEX",
               aggregation_value = summarize_index(.SD$OXYGEN_SCORE)$aggregation_value,
               agreement = "",
               Non_NA_count = summarize_index(.SD$OXYGEN_SCORE)$Non_NA_count,
               NA_count     = summarize_index(.SD$OXYGEN_SCORE)$NA_count),
    data.table(feature = "MOTILITY_INDEX",
               aggregation_value = summarize_index(.SD$MOTILITY_SCORE)$aggregation_value,
               agreement = "",
               Non_NA_count = summarize_index(.SD$MOTILITY_SCORE)$Non_NA_count,
               NA_count     = summarize_index(.SD$MOTILITY_SCORE)$NA_count),
    data.table(feature = "SPORULATION_INDEX",
               aggregation_value = summarize_index(.SD$SPORULATION_SCORE)$aggregation_value,
               agreement = "",
               Non_NA_count = summarize_index(.SD$SPORULATION_SCORE)$Non_NA_count,
               NA_count     = summarize_index(.SD$SPORULATION_SCORE)$NA_count)
  )
  rbind(do.call(rbind, feature_summaries),
        do.call(rbind, indices),
        fill = TRUE)
}, by = "GENUS"]

fwrite(gold_db_genus, file.path(intermediate_files_path, "gold_db_genus.csv"))

# ---- Step 5: load phyloseq ----
ubiome_relative <- readRDS(file.path(phyloseq_obj_files_path, "ubiome_relative.rds"))
stopifnot(!is.null(tax_table(ubiome_relative)))

# ---- Step 6: build genus mapping ----
tax_data <- tax_table(ubiome_relative)
stopifnot("Genus" %in% colnames(tax_data))

ubiome_genus_mapping <- data.frame(
  otu   = rownames(tax_data),
  Genus = tax_data[, "Genus"],
  stringsAsFactors = FALSE
)

# ---- Step 7: parse genus names ----
parse_genus <- function(genus_name) {
  if (is.na(genus_name)) return("unclassified")
  genus <- genus_name

  if (genus_name == "Incertae Sedis")     return("unclassified")
  if (genus_name == "Escherichia/Shigella") return("unclassified")

  if (grepl("^\\[.*\\]", genus_name)) {
    genus <- sub("^\\[([^]]+)\\].*$", "\\1", genus_name)
  }
  if (grepl("^Candidatus\\s+", genus)) {
    genus <- sub("^Candidatus\\s+(\\S+).*$", "\\1", genus)
  }
  if (grepl("sensu stricto", genus)) {
    genus <- sub("^(.*?)\\s+sensu stricto\\s+\\d+.*$", "\\1", genus)
  }
  if (grepl("\\s+", genus) && !grepl("sensu stricto", genus)) {
    if (grepl("\\s+.*group$", genus)) {
      if (grepl("^[A-Z]+\\d+[-]?\\d*", genus)) {
        genus <- sub("^(.*?\\d+[-]?\\d*)\\s+.*$", "\\1", genus)
      } else {
        genus <- sub("^(.*?)\\s+.*group$", "\\1", genus)
      }
    } else {
      genus <- sub("^(.*?)\\s+(\\S+)$", "\\1", genus)
    }
  }
  if (grepl("aceae\\s+.*group$", genus)) {
    genus <- sub("^(.*aceae).*$", "\\1", genus)
  }
  if (genus == "" || (genus == genus_name && !grepl("^[A-Za-z0-9-]+$", genus))) {
    genus <- genus_name
  }
  genus
}

ubiome_genus_mapping$parsed_genus <- sapply(ubiome_genus_mapping$Genus, parse_genus)

# ---- Step 8: merge phyloseq genera with GOLD ----
features_to_merge <- c(
  "BIOTIC_RELATIONSHIPS", "OXYGEN_REQUIREMENT", "METABOLISM", "ENERGY_SOURCES",
  "GRAM_STAIN", "SAMPLE_COLLECTION_DATE", "SALINITY", "CELL_SHAPE", "MOTILITY",
  "SPORULATION", "TEMPERATURE_RANGE", "STAIN_INDEX", "OXYGEN_INDEX",
  "MOTILITY_INDEX", "SPORULATION_INDEX"
)

gold_db_genus_wide <- dcast(
  gold_db_genus,
  GENUS ~ feature,
  value.var    = "aggregation_value",
  fun.aggregate = function(x) if (length(x) > 0) x[1] else NA_character_
)

available_features <- intersect(features_to_merge, names(gold_db_genus_wide))
gold_columns <- c("GENUS", available_features)
gold_db_genus_subset <- gold_db_genus_wide[, ..gold_columns]

ubiome_genus_mapping_complete <- merge(
  ubiome_genus_mapping,
  gold_db_genus_subset,
  by.x = "parsed_genus",
  by.y = "GENUS",
  all.x = TRUE,
  all.y = FALSE
)
setDT(ubiome_genus_mapping_complete)

missing_cols <- setdiff(available_features, names(ubiome_genus_mapping_complete))
if (length(missing_cols) > 0) {
  ubiome_genus_mapping_complete[, (missing_cols) := NA_character_]
}

ubiome_genus_mapping_complete <- ubiome_genus_mapping_complete[,
  lapply(.SD, function(x) fifelse(x == "unknown (0%)", NA_character_, x))
]

# ---- Step 9: mapping statistics ----
total_genera           <- nrow(ubiome_genus_mapping_complete)
unique_original_genera <- length(unique(ubiome_genus_mapping_complete$Genus))
unique_parsed_genera   <- length(unique(ubiome_genus_mapping_complete$parsed_genus))

if (length(available_features) > 0) {
  annotated_rows <- ubiome_genus_mapping_complete[,
    any(!is.na(.SD)), .SDcols = available_features, by = .(otu, Genus)]
  annotated_count <- sum(annotated_rows$V1)

  unclassified_count <- sum(ubiome_genus_mapping_complete$parsed_genus == "unclassified")

  na_counts <- ubiome_genus_mapping_complete[,
    lapply(.SD, is.na), .SDcols = available_features]
  unannotated_count <- sum(rowSums(na_counts) == length(available_features))
} else {
  annotated_count    <- 0
  unclassified_count <- sum(ubiome_genus_mapping_complete$parsed_genus == "unclassified")
  unannotated_count  <- total_genera
}

total_excluding_unclassified <- total_genera - unclassified_count
annotation_coverage <- ifelse(total_genera > 0,
                              round((annotated_count / total_genera) * 100, 2), 0)
annotation_coverage_excl_unclassified <- ifelse(
  total_excluding_unclassified > 0,
  round((annotated_count / total_excluding_unclassified) * 100, 2),
  0
)

# ---- Step 10: save outputs ----
fwrite(ubiome_genus_mapping_complete,
       file.path(intermediate_files_path, "ubiome_genus_mapping_complete.csv"))

summary_stats <- data.table(
  metric = c("total_otus", "unique_original_genera", "unique_parsed_genera",
             "annotated_otus", "unclassified_otus", "unannotated_otus",
             "annotation_coverage_percent",
             "annotation_coverage_excl_unclassified_percent"),
  value  = c(total_genera, unique_original_genera, unique_parsed_genera,
             annotated_count, unclassified_count, unannotated_count,
             annotation_coverage, annotation_coverage_excl_unclassified)
)
fwrite(summary_stats, file.path(intermediate_files_path, "mapping_summary_stats.csv"))

dbDisconnect(con)
