#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 2 — Step 1
# Extract NHANES genus-level AGE and SMOKING effects
#
# Purpose : Turn the module-1 WAS output into a clean per-genus effect table
#           for age and smoking, restricted to the comparison universe
#           (FDR < 0.05 and NHANES prevalence >= 10%).
#
#           Three traps make a naive extraction silently wrong, and each is
#           asserted here rather than commented:
#
#           1. `term` holds the literal model-matrix name. There is NO
#              "indep_var" value. Main effects are the rows flagged by
#              `fdr_corrected`.
#           2. `fdr_corrected` is a "was this row BH-corrected" flag, NOT a
#              significance flag. Covariate and intercept rows carry
#              `p.value.fdr == p.value` (uncorrected), so filtering on FDR
#              before subsetting to main effects returns garbage.
#           3. `SMQ_current_ever_never2` is CURRENT vs never. Two scripts in
#              this repo (modules 6 and 9) label the levels backwards. The
#              SQLite SMQ_F/SMQ_G tables are authoritative and are checked.
#
# Usage   : Rscript 01_extract_nhanes_effects.R
# Outputs : validation/nhanes_effects_{age,smoking,smoking_cotinine}.tsv
#           validation/nhanes_extraction_audit.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(DBI); library(RSQLite)
})
HERE <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/10_validation_analysis"
source(file.path(HERE, "00_validation_config.R"))

vsay("=== Step 1: NHANES genus-level effects ===")

# --- Reference tables --------------------------------------------------------
map <- read.csv(file.path(PROJECT_ROOT,
  "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv"),
  stringsAsFactors = FALSE, na.strings = c("NA", "")) %>%
  select(otu, parsed_genus, Genus)
vsay("otu->genus map: ", nrow(map), " OTUs")

ref <- m12_nhanes_genus_reference(PROJECT_ROOT) %>%
  select(genus, prevalence_F, prevalence_G, prevalence_overall,
         detected_F, detected_G, detected_both, mean_abund_overall, prevalence_tier)

# Phylum from the phyloseq tax_table - the GOLD CSV has no rank columns.
phy_path <- file.path(PROJECT_ROOT,
  "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/nhanes_phyloseq_objects_all.rds")
stopifnot(file.exists(phy_path))
suppressPackageStartupMessages(library(phyloseq))
ps_all <- readRDS(phy_path)
# The saved object is a LIST whose elements are named `counts`, `relative`,
# `relative_clr`, ... - NOT the `ubiome_*` names used inside other modules.
# Resolve by trying the known aliases and fail loudly rather than passing NULL
# into tax_table(), which errors with a misleading "tax_table slot is empty".
ps <- if (inherits(ps_all, "phyloseq")) {
  ps_all
} else {
  cand <- c("relative", "ubiome_relative", "counts", "ubiome_counts")
  hit  <- cand[cand %in% names(ps_all)][1]
  if (is.na(hit))
    stop("no usable phyloseq element in ", basename(phy_path),
         " (have: ", paste(names(ps_all), collapse = ", "), ")")
  vsay("phyloseq element used for taxonomy: '", hit, "'")
  ps_all[[hit]]
}
stopifnot(inherits(ps, "phyloseq"))
tax <- as.data.frame(unclass(phyloseq::tax_table(ps)), stringsAsFactors = FALSE)
tax$otu <- rownames(tax)
phylum_map <- tax %>% transmute(otu, Phylum = as.character(Phylum))
vsay("phylum map: ", nrow(phylum_map), " OTUs, ",
     n_distinct(phylum_map$Phylum), " distinct phyla")
rm(ps_all, ps); invisible(gc())

# --- Trap 3: verify the SMQ level coding against the database ----------------
db <- file.path(PROJECT_ROOT,
  "data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite")
smq_counts <- NULL
if (file.exists(db)) {
  con <- dbConnect(RSQLite::SQLite(), db); on.exit(dbDisconnect(con), add = TRUE)
  got <- lapply(c("SMQ_F", "SMQ_G"), function(t) {
    tryCatch(dbGetQuery(con, sprintf(
      "SELECT SMQ_current_ever_never AS lvl, COUNT(*) AS n FROM %s
        WHERE SMQ_current_ever_never IS NOT NULL GROUP BY lvl", t)) %>% mutate(tbl = t),
      error = function(e) NULL)
  })
  smq_counts <- bind_rows(got)
  if (nrow(smq_counts)) {
    tot <- smq_counts %>% group_by(lvl) %>% summarise(n = sum(n), .groups = "drop") %>%
      mutate(label = SMQ_LEVELS[as.character(lvl)])
    vsay("SMQ_current_ever_never levels (F+G): ",
         paste(sprintf("%s=%s (%d)", tot$lvl, tot$label, tot$n), collapse = ", "))
    # Assert level 2 is the SMALLEST group and level 0 the largest - the
    # never > former > current ordering that identifies 2 as CURRENT.
    n0 <- tot$n[tot$lvl == 0]; n1 <- tot$n[tot$lvl == 1]; n2 <- tot$n[tot$lvl == 2]
    if (!(length(n0) && length(n1) && length(n2) && n0 > n1 && n1 > n2))
      stop("SMQ CODING ASSERTION FAILED: expected never > former > current ",
           "(0 > 1 > 2). Got: ", paste(tot$lvl, tot$n, collapse = " "))
    vsay("[OK] SMQ coding verified: level 2 = CURRENT (smallest group), level 0 = never")
  }
} else {
  vsay("[warn] processed SQLite not found; SMQ coding NOT independently verified")
}

# --- Extraction --------------------------------------------------------------
audit <- list()

extract_one <- function(spec) {
  vsay("--- ", spec$variable, " (", spec$independent_var, " / ", spec$term, ") ---")
  rds <- file.path(PROJECT_ROOT, spec$rds_rel)
  stopifnot(file.exists(rds))
  vsay("  reading ", basename(rds), " (", round(file.size(rds) / 1e6), " MB)")
  d <- readRDS(rds)

  # TRAP 1+2: main effects are the fdr_corrected rows, and only then is
  # p.value.fdr a corrected value.
  if (!"fdr_corrected" %in% names(d))
    stop("expected column `fdr_corrected` not present in ", basename(rds))
  if (any(d$term == "indep_var", na.rm = TRUE))
    stop("unexpected literal term 'indep_var' - the schema has changed")

  main <- d %>% filter(fdr_corrected, independent_var == spec$independent_var,
                       term == spec$term)
  # oradWAS puts the PHENOTYPE on the left and the genus on the right, so the
  # selector above must key off dependent_var instead. Detect and re-filter
  # rather than assume; the row count is asserted below either way.
  if (identical(spec$genus_role, "independent")) {
    main <- d %>% filter(fdr_corrected, dependent_var == spec$independent_var)
  }
  rm(d); invisible(gc())

  n_main <- nrow(main)
  vsay("  main-effect rows: ", n_main)
  if (n_main != spec$expected_n_main)
    warning(sprintf("%s: expected %d main-effect rows, got %d",
                    spec$variable, spec$expected_n_main, n_main), call. = FALSE)

  # Assert the covariate structure so the adjustment asymmetry cannot silently change.
  nc <- unique(main$n_covariates)
  if (!identical(as.integer(nc), as.integer(spec$n_covariates_expected)))
    warning(sprintf("%s: expected n_covariates=%d, got {%s}",
                    spec$variable, spec$n_covariates_expected,
                    paste(nc, collapse = ",")), call. = FALSE)
  vsay("  n_covariates: ", paste(nc, collapse = ","),
       if (identical(as.integer(nc), 0L)) "  (UNADJUSTED)" else "  (adjusted)")

  n_sig_all <- sum(main$p.value.fdr < FDR_THRESHOLD, na.rm = TRUE)
  vsay("  FDR<", FDR_THRESHOLD, " before prevalence filter: ", n_sig_all,
       " of ", n_main)
  if (!is.na(spec$expected_fdr_sig) && n_sig_all != spec$expected_fdr_sig)
    warning(sprintf("%s: expected %d FDR-significant, got %d",
                    spec$variable, spec$expected_fdr_sig, n_sig_all), call. = FALSE)

  # Genus mapping: strip the _relative suffix, drop unclassified
  # Resolve the genus column per this scheme's orientation.
  gcol <- if (identical(spec$genus_role, "independent")) "independent_var" else "dependent_var"
  tb <- main %>%
    mutate(otu = sub("_relative$", "", .data[[gcol]])) %>%
    left_join(map, by = "otu") %>%
    left_join(phylum_map, by = "otu")

  n_unmapped <- sum(is.na(tb$parsed_genus))
  if (n_unmapped) stop(n_unmapped, " taxa did not map to the otu->genus table")
  n_unclass <- sum(tb$parsed_genus == "unclassified")
  tb <- tb %>% filter(parsed_genus != "unclassified")
  vsay("  dropped ", n_unclass, " unclassified OTUs -> ", nrow(tb), " named genera")

  # Split genera: Incertae_Sedis is not a genus; Eubacterium is backed by 7 RSVs.
  tb <- tb %>% filter(parsed_genus != "Incertae_Sedis")
  dup <- tb %>% count(parsed_genus, name = "n_rsv") %>% filter(n_rsv > 1)
  if (nrow(dup)) {
    tb <- tb %>%
      left_join(ref %>% select(genus, prevalence_overall), by = c("parsed_genus" = "genus")) %>%
      group_by(parsed_genus) %>%
      arrange(p.value.fdr, .by_group = TRUE) %>%   # deterministic: smallest FDR wins
      mutate(.n_rsv = n()) %>% slice(1) %>% ungroup()
    for (i in seq_len(nrow(dup)))
      vsay("  [note] ", dup$parsed_genus[i], " backed by ", dup$n_rsv[i],
           " RSVs; kept the lowest-FDR RSV")
  } else {
    tb <- tb %>% left_join(ref %>% select(genus, prevalence_overall),
                           by = c("parsed_genus" = "genus")) %>% mutate(.n_rsv = 1L)
  }

  # Attach the rest of the reference and apply the comparison filters
  out <- tb %>%
    left_join(ref %>% select(-prevalence_overall), by = c("parsed_genus" = "genus")) %>%
    transmute(
      variable = spec$variable, role = spec$role,
      genus = parsed_genus, otu,
      phylum = m12_phylum_std(Phylum), phylum_original = Phylum,
      estimate, std.error, statistic, p.value,
      p.value.fdr, q.value, n_obs,
      direction = ifelse(estimate > 0, "up", ifelse(estimate < 0, "down", NA_character_)),
      sig_fdr = p.value.fdr < FDR_THRESHOLD,
      prevalence_F, prevalence_G, prevalence_overall, detected_F, detected_G,
      detected_both, prevalence_tier, mean_abund_overall,
      n_rsv_backing = .n_rsv,
      independent_var = spec$independent_var, term = spec$term,
      n_covariates = n_covariates,
      adjustment_set = ifelse(n_covariates == 0, "UNADJUSTED (crude)",
                              "age, age^2, sex, PIR, education(4), ethnicity(4), US-born"))

  n_pre <- nrow(out)
  hits <- out %>% filter(sig_fdr, !is.na(prevalence_overall),
                         prevalence_overall >= PREVALENCE_THRESHOLD)
  vsay("  comparison set (FDR<", FDR_THRESHOLD, " & prevalence>=",
       PREVALENCE_THRESHOLD, "): ", nrow(hits), " genera",
       "  [up=", sum(hits$direction == "up"), " down=", sum(hits$direction == "down"), "]")

  audit[[spec$variable]] <<- tibble(
    variable = spec$variable, role = spec$role,
    n_main_effect_rows = n_main, n_covariates = paste(nc, collapse = ","),
    adjustment = ifelse(identical(as.integer(nc), 0L), "UNADJUSTED", "adjusted"),
    n_fdr_sig_all = n_sig_all, n_named_genera = n_pre,
    n_comparison_set = nrow(hits),
    n_up = sum(hits$direction == "up"), n_down = sum(hits$direction == "down"),
    n_obs_model = unique(hits$n_obs)[1])

  # Write the FULL named-genus table (needed later for the reverse direction:
  # external hits that NHANES tested but did not call significant).
  list(all = out, hits = hits)
}

res <- list()
for (i in seq_len(nrow(NHANES_SPEC))) {
  spec <- as.list(NHANES_SPEC[i, ])
  r <- extract_one(spec)
  res[[spec$variable]] <- r
  readr::write_tsv(r$all,  file.path(VAL, paste0("nhanes_effects_", spec$variable, "_all.tsv")))
  readr::write_tsv(r$hits, file.path(VAL, paste0("nhanes_effects_", spec$variable, ".tsv")))
}

aud <- bind_rows(audit)
readr::write_tsv(aud, file.path(VAL, "nhanes_extraction_audit.tsv"))

vsay("")
vsay("=== audit ===")
print(as.data.frame(aud))
vsay("")
vsay("NOTE: age is UNADJUSTED; smoking and gum disease are adjusted for 13 covariates.")
vsay("      Gum disease is modelled phenotype ~ genus (genus is the PREDICTOR), so its")
vsay("      effect scale differs from age/smoking. Directions are comparable; magnitudes are not.")
vsay("      This asymmetry is carried into every downstream table and figure.")
