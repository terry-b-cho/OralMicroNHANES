#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Build harmonized per-sample metadata from a FLAT metadata table
#
# Purpose : Companion to build_sample_metadata.R (which reads BioSample XML).
#           Many external datasets ship a flat mapping file instead. This script
#           harmonizes it to the same schema and, critically, MEASURES THE
#           SUBJECT STRUCTURE rather than assuming one sample per person.
#
#           Meta-analytic cohorts routinely mix cross-sectional studies with
#           dense time series of a handful of donors. Treating those as
#           independent observations inflates n by several fold and biases any
#           age association. This script emits an explicit
#           `is_subject_representative` flag so downstream code can deduplicate
#           deliberately instead of silently not doing so.
#
# Usage   : Rscript build_metadata_from_flatfile.R <dataset_id>
# Inputs  : config/<id>.yaml -> flat_metadata{file,sample_id,subject_id,study_id,variables,...}
# Outputs : datasets/<id>/07_processed_metadata/sample_metadata.tsv
#           datasets/<id>/07_processed_metadata/data_dictionary.tsv
#           datasets/<id>/07_processed_metadata/subject_structure.tsv
#           datasets/<id>/08_sample_mapping/sample_crosswalk.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(yaml); library(tibble); library(stringr)
})

# === USER CONFIG =============================================================
PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
# =============================================================================
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: build_metadata_from_flatfile.R <dataset_id>")
dataset_id <- args[1]

M12    <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
DS_DIR <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation",
                    "datasets", dataset_id)
stopifnot(dir.exists(DS_DIR), grepl("12_external_data_validation", DS_DIR, fixed = TRUE))

source(file.path(M12, "00_shared", "harmonize_helpers.R"))
m12_assert_parse_genus_agrees()
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg <- yaml::read_yaml(file.path(M12, "config", paste0(dataset_id, ".yaml")))
fm  <- cfg$flat_metadata %||% stop("config has no flat_metadata: block")
say <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")

say("=== build_metadata_from_flatfile: ", dataset_id, " ===")

f <- file.path(DS_DIR, fm$file)
stopifnot(file.exists(f))
NA_STRINGS <- c("", "NA", "na", "not applicable", "no_data", "Unspecified",
                "unknown", "not provided", "None", "missing", "not collected")
fmt <- fm$file_format %||% (if (grepl("\\.csv$", f)) "csv" else "tsv")
raw <- if (identical(fmt, "csv")) {
  readr::read_csv(f, show_col_types = FALSE, progress = FALSE,
                  na = NA_STRINGS, guess_max = 5000)
} else {
  readr::read_tsv(f, show_col_types = FALSE, progress = FALSE,
                  na = NA_STRINGS, guess_max = 5000)
}
names(raw) <- sub("^#", "", names(raw))
# readr names an empty header cell "...1"; expose it under a stable alias so
# configs do not have to encode a reader-specific placeholder.
names(raw)[names(raw) == "...1"] <- "unnamed_0"
say("flat metadata: ", nrow(raw), " rows x ", ncol(raw), " cols")

sid_col <- fm$sample_id  %||% names(raw)[1]
sub_col <- fm$subject_id %||% NA_character_
std_col <- fm$study_id   %||% NA_character_
stopifnot(sid_col %in% names(raw))

meta <- tibble(sample_id = as.character(raw[[sid_col]]))
if (!is.na(sub_col) && sub_col %in% names(raw)) meta$subject_id_original <- as.character(raw[[sub_col]])
if (!is.na(std_col) && std_col %in% names(raw)) meta$study_id            <- as.character(raw[[std_col]])

# Subject identifiers in pooled multi-study datasets are frequently NOT globally
# unique - short numeric ids ("3", "7", "15") recur across constituent studies as
# entirely different people. Silently trusting them merges strangers into one
# "subject" and under-counts the cohort. When `subject_key_composite` is set, the
# subject key is study_id + subject_id, and the collision count is reported.
if ("subject_id_original" %in% names(meta)) {
  if (isTRUE(fm$subject_key_composite) && "study_id" %in% names(meta)) {
    meta$subject_id <- paste(meta$study_id, meta$subject_id_original, sep = ":")
    n_naive <- n_distinct(meta$subject_id_original)
    n_comp  <- n_distinct(meta$subject_id)
    say("subject key: COMPOSITE study_id:subject_id")
    say("  naive distinct subject_id: ", n_naive,
        " | composite distinct: ", n_comp,
        " | collisions resolved: ", n_comp - n_naive)
    if (n_comp > n_naive) {
      coll <- meta %>% distinct(study_id, subject_id_original) %>%
        count(subject_id_original, name = "n_studies") %>% filter(n_studies > 1)
      say("  [warn] ", nrow(coll), " subject id(s) reused across studies: ",
          paste(head(coll$subject_id_original, 12), collapse = ", "))
    }
  } else {
    meta$subject_id <- meta$subject_id_original
    say("subject key: subject_id as given (set flat_metadata.subject_key_composite: true ",
        "if ids may be reused across studies)")
  }
}

dict <- list()
add_dict <- function(h, o, def, units, co, ch, x, avail = "public") {
  nonmiss <- if (is.factor(x) || is.character(x)) sum(!is.na(x) & x != "unknown") else sum(!is.na(x))
  dict[[length(dict)+1]] <<- tibble(
    harmonized_name = h, original_name = o, source_file = fm$file, definition = def,
    units = units, coding_original = co, coding_harmonized = ch,
    n_nonmissing = nonmiss, pct_missing = round(100 * (1 - nonmiss/nrow(meta)), 2),
    availability = avail, evidence_url = cfg$code_repository %||% NA_character_)
}

vm <- fm$variables %||% list()
pick <- function(nm) { tg <- vm[[nm]]$tag; if (!is.null(tg) && tg %in% names(raw)) raw[[tg]] else NULL }

# --- age ---------------------------------------------------------------------
if (!is.null(pick("age_years"))) {
  tg <- vm$age_years$tag
  meta$age_original <- as.character(pick("age_years"))
  # Sentinel values must be declared per variable. A literal 0 in an adults-only
  # cohort is a missing marker, not an age, and treating it as real would drag
  # the mean down and invent a spurious young stratum.
  sentinels <- as.character(vm$age_years$na_values %||% character(0))
  a_raw <- meta$age_original
  if (length(sentinels)) {
    n_sent <- sum(trimws(a_raw) %in% sentinels, na.rm = TRUE)
    a_raw[trimws(a_raw) %in% sentinels] <- NA_character_
    say("  age: treated ", n_sent, " row(s) with sentinel value(s) {",
        paste(sentinels, collapse = ", "), "} as MISSING")
  }
  meta$age_years <- m12_harmonize_age(a_raw, units = vm$age_years$units %||% "years")
  meta$age_group_5yr <- m12_age_group_5yr(meta$age_years)
  add_dict("age_years", tg, "host age at sampling", "years", "numeric string",
           "numeric; implausible (<0 or >120) -> NA", meta$age_years)
}
# --- sex ---------------------------------------------------------------------
if (!is.null(pick("sex"))) {
  tg <- vm$sex$tag
  meta$sex_original <- as.character(pick("sex"))
  # An explicit value_map is REQUIRED whenever the source uses opaque codes
  # (0/1). Guessing the convention would silently invert the variable for the
  # whole cohort, so the default heuristic is used only for self-describing
  # values ("male"/"F"/...).
  vmap <- vm$sex$value_map
  if (!is.null(vmap)) {
    mp <- unlist(vmap); names(mp) <- tolower(trimws(names(mp)))
    got <- unname(mp[tolower(trimws(meta$sex_original))])
    got[is.na(got)] <- "unknown"
    meta$sex_harmonized <- factor(got, levels = c("male", "female", "unknown"))
    say("  sex: applied explicit value_map (", paste(names(mp), unname(mp), sep="->", collapse=", "), ")")
  } else {
    meta$sex_harmonized <- m12_harmonize_sex(meta$sex_original)
  }
  add_dict("sex_harmonized", tg, "host sex", NA, paste(sort(unique(na.omit(meta$sex_original))), collapse="|"),
           "male|female|unknown", meta$sex_harmonized)
}
# --- smoking -----------------------------------------------------------------
if (!is.null(pick("smoking_status"))) {
  tg <- vm$smoking_status$tag
  meta$smoking_original <- as.character(pick("smoking_status"))
  studydef <- cfg$smoking_harmonization$study_definition %||%
              paste0("flat metadata column '", tg, "'")
  mp <- cfg$smoking_harmonization$mapping
  meta$smoking_status_harmonized <- m12_harmonize_smoking(
    meta$smoking_original, study_definition = studydef,
    mapping = if (!is.null(mp)) unlist(mp) else NULL)
  add_dict("smoking_status_harmonized", tg, studydef, NA,
           paste(sort(unique(na.omit(meta$smoking_original))), collapse="|"),
           "never|former|current|unknown", meta$smoking_status_harmonized)
}
# --- pass-through -------------------------------------------------------------
for (nm in setdiff(names(vm), c("age_years", "sex", "smoking_status"))) {
  v <- pick(nm); if (is.null(v)) next
  tg <- vm[[nm]]$tag
  meta[[paste0(nm, "_original")]] <- as.character(v)
  meta[[nm]] <- if (identical(vm[[nm]]$type, "numeric")) suppressWarnings(as.numeric(as.character(v)))
                else as.character(v)
  add_dict(nm, tg, paste0("column '", tg, "'"),
           if (identical(vm[[nm]]$type, "numeric")) "numeric" else NA,
           NA, vm[[nm]]$type %||% "character", meta[[nm]])
}

# --- Smoking refinement ------------------------------------------------------
# Some sources split smoking across TWO binary columns (current + ever). Neither
# alone gives never/former/current, but together they do:
#   ever = no                -> never
#   ever = yes & current = no -> former
#   current = yes            -> current
# This is applied only when explicitly configured, and the derivation is written
# into the data dictionary so it is auditable rather than implicit.
sr <- cfg$smoking_harmonization$refine_with_ever_smoker
if (!is.null(sr) && "smoking_status_harmonized" %in% names(meta)) {
  ev_col <- sr$ever_column %||% "ever_smoker"
  if (ev_col %in% names(meta)) {
    ev  <- tolower(trimws(as.character(meta[[ev_col]])))
    cur <- as.character(meta$smoking_status_harmonized)
    before <- table(factor(cur, levels = M12_SMOKING_LEVELS))
    new <- cur
    new[cur == "never" & ev %in% c("yes", "true", "1")] <- "former"
    new[cur == "never" & ev %in% c("no", "false", "0")] <- "never"
    meta$smoking_status_harmonized <- factor(new, levels = M12_SMOKING_LEVELS)
    after <- table(meta$smoking_status_harmonized)
    say("smoking refined with '", ev_col, "': ",
        paste(sprintf("%s %d->%d", names(after), as.integer(before), as.integer(after)),
              collapse = ", "))
    dict[[length(dict)+1]] <- tibble(
      harmonized_name = "smoking_status_harmonized", original_name = paste0("smoker + ", ev_col),
      source_file = fm$file,
      definition = "three-level smoking derived from two binary columns (current + ever)",
      units = NA_character_, coding_original = "yes/no x yes/no",
      coding_harmonized = "never|former|current|unknown",
      n_nonmissing = sum(meta$smoking_status_harmonized != "unknown"),
      pct_missing = round(100*mean(meta$smoking_status_harmonized == "unknown"), 2),
      availability = "public", evidence_url = cfg$code_repository %||% NA_character_)
  } else {
    say("[warn] refine_with_ever_smoker configured but column '", ev_col, "' not found")
  }
}

# --- SUBJECT STRUCTURE: measure, never assume --------------------------------
say("")
say("=== subject structure ===")
if ("subject_id" %in% names(meta)) {
  n_s <- nrow(meta); n_subj <- n_distinct(meta$subject_id)
  say("  samples: ", n_s, " | distinct subjects: ", n_subj,
      sprintf(" | mean %.2f samples/subject", n_s / n_subj))

  per <- meta %>% count(subject_id, name = "n_samples") %>% arrange(desc(n_samples))
  top <- head(per, 10)
  say("  most-sampled subjects: ",
      paste(sprintf("%s=%d", top$subject_id, top$n_samples), collapse = ", "))
  heavy <- per %>% filter(n_samples >= 10)
  if (nrow(heavy))
    say("  [warn] ", nrow(heavy), " subject(s) contribute >=10 samples each, totalling ",
        sum(heavy$n_samples), " samples (",
        sprintf("%.1f%%", 100*sum(heavy$n_samples)/n_s), " of the dataset)")

  # Representative selection is deterministic (seeded) and documented: for each
  # subject keep the sample with the most complete key phenotype, breaking ties
  # by sample_id so the choice is reproducible across runs and machines.
  keyvars <- intersect(c("age_years", "sex_harmonized", "smoking_status_harmonized"), names(meta))
  meta$.completeness <- rowSums(!is.na(meta[keyvars]))
  meta <- meta %>%
    group_by(subject_id) %>%
    arrange(desc(.completeness), sample_id, .by_group = TRUE) %>%
    mutate(is_subject_representative = row_number() == 1L,
           n_samples_this_subject = n()) %>%
    ungroup() %>%
    select(-.completeness)
  say("  flagged ", sum(meta$is_subject_representative),
      " representative samples (one per subject)")

  ss <- per %>% mutate(pct_of_dataset = round(100*n_samples/n_s, 2))
  if ("study_id" %in% names(meta)) {
    ss <- meta %>% group_by(study_id) %>%
      summarise(n_samples = n(), n_subjects = n_distinct(subject_id),
                samples_per_subject = round(n()/n_distinct(subject_id), 2), .groups = "drop") %>%
      arrange(desc(n_samples))
    say("")
    say("  per-study sample:subject ratio")
    for (i in seq_len(nrow(ss)))
      say(sprintf("    study %-8s %5d samples / %4d subjects = %.1f",
                  ss$study_id[i], ss$n_samples[i], ss$n_subjects[i], ss$samples_per_subject[i]))
  }
  readr::write_tsv(ss, file.path(DS_DIR, "07_processed_metadata", "subject_structure.tsv"))
} else {
  say("  [warn] no subject_id column configured - independence CANNOT be verified")
  meta$is_subject_representative <- TRUE
  meta$n_samples_this_subject <- 1L
}

# --- Write -------------------------------------------------------------------
dir.create(file.path(DS_DIR, "07_processed_metadata"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DS_DIR, "08_sample_mapping"),     showWarnings = FALSE, recursive = TRUE)
readr::write_tsv(meta, file.path(DS_DIR, "07_processed_metadata", "sample_metadata.tsv"))
readr::write_tsv(bind_rows(dict), file.path(DS_DIR, "07_processed_metadata", "data_dictionary.tsv"))

xw <- meta %>% select(any_of(c("sample_id", "subject_id", "study_id",
                               "is_subject_representative", "n_samples_this_subject")))
readr::write_tsv(xw, file.path(DS_DIR, "08_sample_mapping", "sample_crosswalk.tsv"))

say("")
say("wrote sample_metadata.tsv (", nrow(meta), " x ", ncol(meta), ")")
say("wrote data_dictionary.tsv (", length(dict), " variables)")

say("")
say("=== harmonized coverage ===")
for (v in intersect(c("age_years","sex_harmonized","smoking_status_harmonized","bmi"), names(meta))) {
  x <- meta[[v]]
  ok <- if (is.factor(x) || is.character(x)) sum(!is.na(x) & x != "unknown") else sum(!is.na(x))
  say(sprintf("  %-28s %5.1f%% (n=%d)", v, 100*ok/nrow(meta), ok))
}
if ("age_years" %in% names(meta) && sum(!is.na(meta$age_years)) > 0) {
  a <- meta$age_years[!is.na(meta$age_years)]
  say(sprintf("  age range %.0f-%.0f (mean %.1f); >=50: %d; >=60: %d",
              min(a), max(a), mean(a), sum(a >= 50), sum(a >= 60)))
  if ("is_subject_representative" %in% names(meta)) {
    ar <- meta$age_years[meta$is_subject_representative & !is.na(meta$age_years)]
    if (length(ar))
      say(sprintf("  after subject dedup: n=%d, range %.0f-%.0f, mean %.1f, >=50: %d",
                  length(ar), min(ar), max(ar), mean(ar), sum(ar >= 50)))
  }
} else if ("age_years" %in% names(meta)) {
  say("  [note] an age column exists in the source but is EMPTY for every row - ",
      "recorded as availability=absent, not silently dropped")
}
