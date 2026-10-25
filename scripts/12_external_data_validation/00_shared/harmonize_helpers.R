#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Shared harmonization helpers: genus parsing + phenotype recoding
#
# Purpose : Make external datasets joinable to the NHANES analytic tables
#           (956 unique parsed genera across 1,349 OTUs) WITHOUT ever
#           overwriting an original value. Every harmonized field is created
#           alongside its source (`age_original` + `age_years`).
#
# Usage   : source("scripts/12_external_data_validation/00_shared/harmonize_helpers.R")
# Env     : R 4.5.1 (conda env nhanes-analysis)
# =============================================================================

# --- Genus name parser -------------------------------------------------------
# PROVENANCE: lifted VERBATIM from
#   scripts/3_gold_db_microbial_phenotype/gold_db_process_n_genus_mapping.Rmd
#   chunk `parse_genus_names`, function defined at line 558.
# That definition lives inside an Rmd chunk and therefore cannot be source()d,
# so it is duplicated here. `m12_assert_parse_genus_agrees()` below guards
# against the two copies drifting apart. DO NOT EDIT without updating module 3
# and re-running the drift assertion.
parse_genus <- function(genus_name) {
  if (is.na(genus_name)) return("unclassified")

  genus <- genus_name

  if (genus_name == "Incertae Sedis") return("unclassified")
  if (genus_name == "Escherichia/Shigella") return("unclassified")

  # 1. "[<Genus>] <detail> group" or "[<Genus>] <detail>" -> <Genus>
  if (grepl("^\\[.*\\]", genus_name)) {
    genus <- sub("^\\[([^]]+)\\].*$", "\\1", genus_name)
  }
  # 2. "Candidatus <Genus>" -> <Genus>
  if (grepl("^Candidatus\\s+", genus)) {
    genus <- sub("^Candidatus\\s+(\\S+).*$", "\\1", genus)
  }
  # 3. "<Genus> sensu stricto <n>" -> <Genus>
  if (grepl("sensu stricto", genus)) {
    genus <- sub("^(.*?)\\s+sensu stricto\\s+\\d+.*$", "\\1", genus)
  }
  # 4. group patterns
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
  # 5. "aceae" families with "group" suffix
  if (grepl("aceae\\s+.*group$", genus)) {
    genus <- sub("^(.*aceae).*$", "\\1", genus)
  }

  if (genus == "" || (genus == genus_name && !grepl("^[A-Za-z0-9-]+$", genus))) {
    genus <- genus_name
  }
  genus
}

# Drift guard: fixed cases whose expected output is pinned. Run by 07_qc.
M12_PARSE_GENUS_CASES <- list(
  list(input = "Candidatus Saccharimonas",         expect = "Saccharimonas"),
  list(input = "[Eubacterium] brachy group",       expect = "Eubacterium"),
  list(input = "Incertae Sedis",                   expect = "unclassified"),
  list(input = "Escherichia/Shigella",             expect = "unclassified"),
  list(input = "Clostridium sensu stricto 1",      expect = "Clostridium"),
  list(input = "Prevotella",                       expect = "Prevotella"),
  list(input = "Veillonella",                      expect = "Veillonella")
)

m12_assert_parse_genus_agrees <- function(verbose = TRUE) {
  fails <- character(0)
  for (tc in M12_PARSE_GENUS_CASES) {
    got <- parse_genus(tc$input)
    if (!identical(got, tc$expect)) {
      fails <- c(fails, sprintf("  '%s' -> '%s' (expected '%s')", tc$input, got, tc$expect))
    }
  }
  if (length(fails)) {
    stop("ASSERTION FAILED: parse_genus() drifted from module 3:\n",
         paste(fails, collapse = "\n"))
  }
  if (verbose) cat("[OK] parse_genus() matches module 3 on", length(M12_PARSE_GENUS_CASES), "pinned cases\n")
  invisible(TRUE)
}

# --- Genus derivation from lineage strings -----------------------------------
# Many external tables encode taxonomy as a single delimited lineage
# ("Bacteria_Firmicutes_...Streptococcus_oralis" or "k__..;p__..;g__Rothia").
# Returns the genus token, or NA when absent.
m12_genus_from_lineage <- function(x, style = c("auto", "qiime", "underscore")) {
  style <- match.arg(style)
  out <- vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    if (style %in% c("auto", "qiime") && grepl("g__", s, fixed = TRUE)) {
      g <- sub(".*g__([^;|]*).*", "\\1", s)
      g <- trimws(g)
      return(if (nzchar(g)) g else NA_character_)
    }
    if (style %in% c("auto", "underscore") && grepl("_", s, fixed = TRUE)) {
      # Convention: Domain_Phylum_Class_Order_Family_Genus[_species...]
      parts <- strsplit(s, "_", fixed = TRUE)[[1]]
      return(if (length(parts) >= 6) parts[6] else NA_character_)
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
  unname(out)
}

# --- NHANES reference genus universe -----------------------------------------
# The join target: 956 unique parsed genera across 1,349 OTUs.
#
# GOTCHA (verified): the module 3 CSV writes true NA Genus values as an empty
# field. Reading it without `na.strings = c("NA", "")` yields "" which
# parse_genus() passes through instead of mapping to "unclassified", producing
# 384 spurious mismatches. With empty-as-NA, re-parsing reproduces module 3's
# `parsed_genus` for all 1,349 rows exactly (verified 2026-08-30).
M12_NHANES_GENUS_MAP <- file.path(
  "results", "analyses_results", "03_gold_db_microbial_phenotype_out",
  "intermediate", "ubiome_genus_mapping_complete.csv")

m12_nhanes_genus_universe <- function(project_root = Sys.getenv("ORALMICRO_ROOT", getwd())) {
  p <- file.path(project_root, M12_NHANES_GENUS_MAP)
  if (!file.exists(p)) stop("NHANES genus map not found: ", p)
  m <- utils::read.csv(p, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  sort(setdiff(unique(m$parsed_genus), c("unclassified", NA)))
}

# --- Cross-database genus comparison key -------------------------------------
# NHANES uses SILVA nomenclature, which splits several genera into numbered
# clades: Prevotella_2/_6/_7, Selenomonas_3/_4, Treponema_2, Clostridium_sensu
# _stricto_1..10. External datasets built on eHOMD or Greengenes report the
# plain base name. Comparing the two literally makes highly prevalent NHANES
# genera look absent: Prevotella_7 is in 98.3% of participants and would never
# match an eHOMD "Prevotella".
#
# This produces a COMPARISON KEY only. It is a pragmatic cross-database
# harmonization, NOT a taxonomic claim: in SILVA, Prevotella_7 and Prevotella
# are distinct genus-level groups, and collapsing them loses that distinction.
# The original name is always retained alongside, and overlap is reported both
# strictly and normalized so the effect of this choice is visible.
m12_genus_crossdb_key <- function(x) {
  g <- as.character(x)
  g <- sub("_sensu_stricto(_\\d+)?$", "", g)   # Clostridium_sensu_stricto_1 -> Clostridium
  g <- sub("_\\d+$", "", g)                     # Prevotella_7 -> Prevotella
  g <- sub("_(group|subgroup)$", "", g, ignore.case = TRUE)
  g <- gsub("^\\[|\\]$", "", g)                 # [Prevotella] -> Prevotella
  trimws(g)
}

#' Load the per-cycle NHANES genus reference built by
#' 09_reports/build_nhanes_genus_reference.R. Prefer this over the flat
#' taxonomy list: it distinguishes detected from merely-catalogued genera and
#' carries per-cycle prevalence.
m12_nhanes_genus_reference <- function(project_root = Sys.getenv("ORALMICRO_ROOT", getwd())) {
  p <- file.path(project_root, "results/analyses_results/12_external_data_validation",
                 "nhanes_genus_reference.tsv")
  if (!file.exists(p))
    stop("NHANES genus reference not found. Run 09_reports/build_nhanes_genus_reference.R first: ", p)
  utils::read.delim(p, stringsAsFactors = FALSE)
}

#' Genus overlap between an external dataset and NHANES.
#' Reported per dataset in the readiness report — validation is scoped to the
#' abundant oral core, never implied to be whole-microbiome replication.
m12_genus_overlap <- function(external_genera, project_root = Sys.getenv("ORALMICRO_ROOT", getwd())) {
  ext <- sort(unique(stats::na.omit(external_genera)))
  ext <- setdiff(ext, "unclassified")
  nh  <- m12_nhanes_genus_universe(project_root)
  shared <- intersect(ext, nh)
  list(n_external = length(ext), n_nhanes = length(nh), n_shared = length(shared),
       pct_of_external = if (length(ext)) 100 * length(shared) / length(ext) else NA_real_,
       shared = shared, external_only = setdiff(ext, nh))
}

# --- Phenotype harmonization -------------------------------------------------
# NEVER overwrite an original column. These return NEW vectors; callers keep
# both (`smoking_original` and `smoking_status_harmonized`).

M12_SMOKING_LEVELS <- c("never", "former", "current", "unknown")

#' Harmonize smoking to {never, former, current, unknown}.
#' `study_definition` is REQUIRED and recorded in the data dictionary, because
#' study definitions genuinely differ (e.g. some cohorts encode only
#' never/current, leaving former smokers indistinguishable from never).
m12_harmonize_smoking <- function(x, study_definition, mapping = NULL) {
  if (missing(study_definition) || !nzchar(study_definition)) {
    stop("m12_harmonize_smoking(): study_definition is required and must be recorded.")
  }
  v <- tolower(trimws(as.character(x)))
  if (!is.null(mapping)) {
    # Values are lowercased for matching, so the mapping keys must be too.
    # Without this, a config key like `Current:` never matches and its samples
    # silently fall through to "unknown" - a misclassification that is invisible
    # unless you tally the levels. Collisions are reported, never resolved silently.
    names(mapping) <- tolower(trimws(names(mapping)))
    # Redundant case variants (Never/never -> never) are harmless and deduped.
    # A key mapping to two DIFFERENT targets is a real config error and stops.
    conflict <- tapply(unname(mapping), names(mapping),
                       function(z) length(unique(z)) > 1)
    if (any(conflict)) {
      stop("m12_harmonize_smoking(): mapping key(s) assign conflicting values: ",
           paste(names(conflict)[conflict], collapse = ", "))
    }
    mapping <- mapping[!duplicated(names(mapping))]
    out <- unname(mapping[v])
    unmapped <- setdiff(unique(v[!is.na(v) & nzchar(v)]), names(mapping))
    if (length(unmapped)) {
      warning("m12_harmonize_smoking(): ", length(unmapped),
              " unmapped level(s) -> unknown: ", paste(unmapped, collapse = ", "),
              call. = FALSE)
    }
  } else {
    out <- dplyr::case_when(
      v %in% c("never", "never smoker", "nonsmoker", "non-smoker", "no", "0")   ~ "never",
      v %in% c("former", "ex", "ex-smoker", "former smoker", "quit", "past")    ~ "former",
      v %in% c("current", "current smoker", "smoker", "yes", "daily", "1")      ~ "current",
      TRUE                                                                       ~ "unknown"
    )
  }
  out[is.na(out)] <- "unknown"
  bad <- setdiff(unique(out), M12_SMOKING_LEVELS)
  if (length(bad)) stop("m12_harmonize_smoking(): illegal level(s): ", paste(bad, collapse = ", "))
  factor(out, levels = M12_SMOKING_LEVELS)
}

m12_harmonize_sex <- function(x) {
  v <- tolower(trimws(as.character(x)))
  out <- dplyr::case_when(
    v %in% c("m", "male", "1")   ~ "male",
    v %in% c("f", "female", "2") ~ "female",
    TRUE                          ~ "unknown"
  )
  factor(out, levels = c("male", "female", "unknown"))
}

# NHANES-aligned 5-year age bins, identical to scripts/5_age_analyses (lines 78-79)
M12_AGE_BREAKS <- c(14, seq(20, 85, by = 5))
M12_AGE_LABELS <- c("14-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                    "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84")

m12_age_group_5yr <- function(age_years) {
  cut(as.numeric(age_years), breaks = M12_AGE_BREAKS,
      labels = M12_AGE_LABELS, right = FALSE)
}

m12_harmonize_age <- function(x, units = c("years", "months")) {
  units <- match.arg(units)
  a <- suppressWarnings(as.numeric(as.character(x)))
  if (units == "months") a <- a / 12
  a[!is.na(a) & (a < 0 | a > 120)] <- NA_real_   # implausible -> NA, never silently kept
  a
}

# --- Join safety -------------------------------------------------------------
#' Assert a join lost no rows. Module 12 forbids silent row loss.
m12_assert_no_row_loss <- function(before_n, after_n, context = "") {
  if (after_n != before_n) {
    stop(sprintf("ASSERTION FAILED: row count changed in join%s: %d -> %d",
                 if (nzchar(context)) paste0(" [", context, "]") else "",
                 before_n, after_n))
  }
  invisible(TRUE)
}
