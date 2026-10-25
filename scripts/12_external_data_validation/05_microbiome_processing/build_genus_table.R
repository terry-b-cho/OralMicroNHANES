#!/usr/bin/env Rscript
# =============================================================================
# Module 12 — External Data Validation
# Build a genus-level abundance table from a species/ASV table + taxonomy key
#
# Purpose : Collapse an external feature table to GENUS so it can be compared
#           with the NHANES genus universe (956 parsed genera / 1,349 OTUs).
#
#           Collapsing is done by POSITIONAL JOIN against the dataset's taxonomy
#           key, never by splitting lineage strings on a delimiter: eHOMD and
#           SILVA rank names themselves contain underscores and brackets
#           (Lachnospiraceae_[XIV], Saccharibacteria_(TM7)_[G-1]), so a
#           fixed-field split silently corrupts a fraction of features.
#
# Usage   : Rscript build_genus_table.R <dataset_id>
# Inputs  : config/<id>.yaml  -> microbiome.{abundance_file,taxonomy_key_file,genus_column}
#           datasets/<id>/05_processed_microbiome/<abundance>
#           datasets/<id>/05_processed_microbiome/<taxonomy key>
# Outputs : datasets/<id>/05_processed_microbiome/genus_table.tsv        (counts)
#           datasets/<id>/05_processed_microbiome/genus_table_relative.tsv
#           datasets/<id>/05_processed_microbiome/taxonomy_map.tsv
#           datasets/<id>/05_processed_microbiome/genus_overlap_nhanes.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(yaml); library(tibble) })

# === USER CONFIG =============================================================
PROJECT_ROOT <- Sys.getenv("ORALMICRO_ROOT", getwd())
# =============================================================================
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: build_genus_table.R <dataset_id>")
dataset_id <- args[1]

M12    <- file.path(PROJECT_ROOT, "scripts/12_external_data_validation")
DS_DIR <- file.path(PROJECT_ROOT, "results/analyses_results/12_external_data_validation",
                    "datasets", dataset_id)
MIC    <- file.path(DS_DIR, "05_processed_microbiome")
stopifnot(dir.exists(MIC), grepl("12_external_data_validation", DS_DIR, fixed = TRUE))

source(file.path(M12, "00_shared", "harmonize_helpers.R"))
m12_assert_parse_genus_agrees()

cfg <- yaml::read_yaml(file.path(M12, "config", paste0(dataset_id, ".yaml")))
mc  <- cfg$microbiome
say <- function(...) cat("[", format(Sys.time(), tz = "UTC", "%H:%M:%SZ"), "] ", ..., "\n", sep = "")

abund_f <- file.path(MIC, mc$abundance_file    %||% "species.csv")
taxo_f  <- file.path(MIC, mc$taxonomy_key_file %||% "Species_OTU_diversity.csv")
`%||%` <- function(a, b) if (is.null(a)) b else a
stopifnot(file.exists(abund_f), file.exists(taxo_f))

say("=== build_genus_table: ", dataset_id, " ===")

# --- Abundance: samples x features -------------------------------------------
ab <- if (grepl("\\.tsv$|\\.txt$", abund_f)) {
  readr::read_tsv(abund_f, show_col_types = FALSE, progress = FALSE, guess_max = 5000)
} else {
  readr::read_csv(abund_f, show_col_types = FALSE, progress = FALSE, guess_max = 5000)
}
id_col <- names(ab)[1]
# Orientation is NOT guessable from shape alone - a 133x1304 table could be
# either. Declare it in config; guessing here would silently transpose a matrix
# and swap every sample for a taxon.
orientation <- mc$orientation %||% "samples_as_rows"
if (identical(orientation, "features_as_rows")) {
  feature_ids <- as.character(ab[[1]])
  mat <- t(as.matrix(ab[, -1, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  sample_ids <- rownames(mat)
  colnames(mat) <- feature_ids
  say("abundance: TRANSPOSED from features_as_rows (", length(feature_ids),
      " features x ", length(sample_ids), " samples)")
} else {
  sample_ids <- as.character(ab[[1]])
  mat <- as.matrix(ab[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
}
say("abundance: ", nrow(mat), " samples x ", ncol(mat), " features (id column '", id_col, "')")

# --- Taxonomy key ------------------------------------------------------------
tx <- if (grepl("\\.tsv$|\\.txt$", taxo_f)) {
  readr::read_tsv(taxo_f, show_col_types = FALSE, progress = FALSE)
} else {
  readr::read_csv(taxo_f, show_col_types = FALSE, progress = FALSE)
}
say("taxonomy key: ", nrow(tx), " rows x ", ncol(tx), " cols")

# Two join modes, because external datasets use both and confusing them
# silently mislabels every taxon:
#   positional     - taxonomy row i describes abundance column i
#   by_feature_id  - taxonomy has an explicit feature-id column (e.g. an ASV
#                    sequence or OTU id) matched against the column names
join_mode <- mc$taxonomy_join_mode %||% "positional"
say("taxonomy join mode: ", join_mode)

if (identical(join_mode, "by_feature_id")) {
  fid_col <- mc$taxonomy_feature_id_column %||% names(tx)[1]
  lin_col <- mc$taxonomy_lineage_column     %||% names(tx)[2]
  stopifnot(fid_col %in% names(tx), lin_col %in% names(tx))
  say("  feature id column: '", fid_col, "' | lineage column: '", lin_col, "'")

  key <- setNames(as.character(tx[[lin_col]]), as.character(tx[[fid_col]]))
  lineage <- unname(key[colnames(mat)])
  n_hit <- sum(!is.na(lineage))
  say(sprintf("  matched %d/%d features (%.1f%%) to the taxonomy key",
              n_hit, ncol(mat), 100*n_hit/ncol(mat)))
  if (n_hit == 0)
    stop("FEATURE-ID JOIN FAILED: no abundance column matched the taxonomy key. ",
         "Check taxonomy_feature_id_column.")
  if (n_hit / ncol(mat) < 0.90)
    warning("Only ", round(100*n_hit/ncol(mat)), "% of features matched taxonomy; ",
            "unmatched features become 'unclassified' and will be reported.")

  # Some taxonomy keys carry a plain genus name rather than a full lineage
  # string. Parsing "Abiotrophia" as a lineage yields NA and would silently
  # mark every feature unclassified, so this must be declared, not inferred.
  if (isTRUE(mc$taxonomy_lineage_is_genus)) {
    say("  lineage column holds plain genus names (no lineage parsing)")
    genus_raw <- lineage
  } else {
    genus_raw <- m12_genus_from_lineage(lineage, style = mc$lineage_style %||% "qiime")
  }
  genus_raw[is.na(genus_raw) | !nzchar(genus_raw)] <- "unclassified"
  n_unc <- sum(genus_raw == "unclassified")
  say(sprintf("  genus assigned for %d/%d features; %d unclassified (%.1f%%)",
              ncol(mat) - n_unc, ncol(mat), n_unc, 100*n_unc/ncol(mat)))

  tx_aligned <- tibble(feature_id = colnames(mat), lineage = lineage)

} else {
  genus_col <- mc$genus_column %||% {
    cand <- names(tx)[tolower(names(tx)) %in% c("genus", "g", "rank6")]
    if (!length(cand)) stop("No genus column found in taxonomy key; set microbiome.genus_column")
    cand[1]
  }
  say("genus column: '", genus_col, "'")

  # --- The positional-join contract ------------------------------------------
  # Row i of the taxonomy key MUST describe column i of the abundance matrix.
  # Assert it rather than trust it: a silent off-by-one mislabels every taxon.
  if (nrow(tx) != ncol(mat)) {
    stop("POSITIONAL JOIN VIOLATED: taxonomy key has ", nrow(tx),
         " rows but abundance matrix has ", ncol(mat), " feature columns. ",
         "These must correspond 1:1 and in order.")
  }
  say("[OK] positional contract holds: ", nrow(tx), " taxonomy rows == ", ncol(mat), " feature columns")

  sp_col <- names(tx)[tolower(names(tx)) %in% c("species", "s")]
  if (length(sp_col)) {
    sp <- as.character(tx[[sp_col[1]]])
    hits <- mapply(function(cn, s) {
      if (is.na(s) || !nzchar(s)) NA
      else grepl(gsub("[^A-Za-z0-9]", "", s), gsub("[^A-Za-z0-9]", "", cn), fixed = TRUE)
    }, colnames(mat), sp)
    ok <- mean(hits, na.rm = TRUE)
    say(sprintf("alignment corroboration: %.1f%% of feature names contain their key's species epithet",
                100 * ok))
    if (!is.na(ok) && ok < 0.5)
      warning("Positional alignment is weakly corroborated (", round(100*ok), "%). Verify before use.")
  }

  genus_raw <- as.character(tx[[genus_col]])
  genus_raw[is.na(genus_raw) | !nzchar(genus_raw)] <- "unclassified"
  tx_aligned <- tx
}

# Apply the SAME parser NHANES used, so the two vocabularies are comparable.
genus_parsed <- vapply(genus_raw, parse_genus, character(1), USE.NAMES = FALSE)

taxo_map <- tibble(
  feature_index = seq_len(ncol(mat)),
  feature_name  = colnames(mat),
  genus_original = genus_raw,
  genus_parsed   = genus_parsed
)
for (r in names(tx_aligned)) taxo_map[[paste0("key_", r)]] <- as.character(tx_aligned[[r]])

# --- Collapse to genus -------------------------------------------------------
total_before <- sum(mat, na.rm = TRUE)
g <- factor(genus_parsed)
genus_mat <- t(rowsum(t(mat), group = g, na.rm = TRUE))
total_after <- sum(genus_mat, na.rm = TRUE)
say("collapsed ", ncol(mat), " features -> ", ncol(genus_mat), " genera")

# Count preservation is a hard assertion: collapsing must move no mass.
if (!isTRUE(all.equal(total_before, total_after))) {
  stop("MASS NOT PRESERVED by genus collapse: ", total_before, " -> ", total_after)
}
say("[OK] total abundance preserved: ", format(round(total_after, 4), big.mark = ","))

genus_df <- as_tibble(genus_mat) %>% mutate(sample_id = sample_ids, .before = 1)
rel <- genus_mat / pmax(rowSums(genus_mat, na.rm = TRUE), 1)
rel_df <- as_tibble(rel) %>% mutate(sample_id = sample_ids, .before = 1)

readr::write_tsv(genus_df, file.path(MIC, "genus_table.tsv"))
readr::write_tsv(rel_df,   file.path(MIC, "genus_table_relative.tsv"))
readr::write_tsv(taxo_map, file.path(MIC, "taxonomy_map.tsv"))
say("wrote genus_table.tsv (", nrow(genus_df), " x ", ncol(genus_df) - 1, ")")

# --- Overlap with the NHANES genus universe ---------------------------------
ext_genera <- setdiff(colnames(genus_mat), "unclassified")
ov <- m12_genus_overlap(ext_genera, PROJECT_ROOT)
say("")
say("=== genus overlap with NHANES ===")
say("  external genera : ", ov$n_external)
say("  NHANES genera   : ", ov$n_nhanes)
say("  shared          : ", ov$n_shared, sprintf(" (%.1f%% of external)", ov$pct_of_external))
if (length(ov$external_only))
  say("  external-only   : ", paste(head(ov$external_only, 12), collapse = ", "),
      if (length(ov$external_only) > 12) " ..." else "")

ov_tbl <- tibble(
  genus = sort(ext_genera),
  in_nhanes = genus %in% ov$shared,
  total_count = as.numeric(colSums(genus_mat, na.rm = TRUE)[sort(ext_genera)]),
  prevalence = as.numeric(colMeans(genus_mat[, sort(ext_genera), drop = FALSE] > 0, na.rm = TRUE))
) %>% arrange(desc(total_count))
readr::write_tsv(ov_tbl, file.path(MIC, "genus_overlap_nhanes.tsv"))
say("wrote genus_overlap_nhanes.tsv")

say("")
say("NOTE: this feature space is a cross-study INTERSECTION, heavily pruned.")
say("An NHANES genus absent here was filtered out upstream, not absent from saliva.")
