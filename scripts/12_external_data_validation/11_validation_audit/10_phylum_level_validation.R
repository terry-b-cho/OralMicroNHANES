#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 10
# Phylum-level validation against SLL2's published phylum rows
#
# WHY THIS EXISTS: SLL2's main-text Tables 1 and 2 report directions at THREE
# levels - genus, phylum and alpha diversity. Phases 2 and 3 used only the
# genus rows. The phylum rows are an independent, coarser line of evidence that
# was sitting unused in a file we already held.
#
# NHANES has no phylum-level regression, so the NHANES phylum direction is
# DERIVED: a sign vote across that phylum's primary genus associations. This is
# a summary of the genus-level evidence, NOT a phylum-level model, and it is
# labelled that way everywhere. A phylum is called only when the vote is not a
# tie; ties are reported as indeterminate rather than forced.
#
# Usage   : Rscript 10_phylum_level_validation.R
# Outputs : audit/phylum_level_validation.tsv
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 10: phylum-level validation ===")

DSD <- file.path(M12_OUT, "datasets", "WillisJR_2022_SacaLaLengua2",
                 "02_supplementary", "derived")
prim <- read_tsv(file.path(AUD, "nhanes_primary_associations.tsv"), show_col_types = FALSE)

# --- Parse the published phylum rows -----------------------------------------
# Both tables use a blocked layout: the level/group appears once and the rows
# beneath it inherit it, so column 1 is forward-filled before filtering.
arrow_dir <- function(x) dplyr::case_when(
  grepl("↗", x) ~ "up", grepl("↘", x) ~ "down", TRUE ~ NA_character_)

parse_tbl <- function(file, block_col, name_col, dir_col, want_block) {
  t <- suppressWarnings(read_tsv(file.path(DSD, file), show_col_types = FALSE,
                                 progress = FALSE, col_types = readr::cols(.default = "c")))
  names(t)[seq_len(5)] <- c("c1", "c2", "c3", "c4", "c5")
  t <- t[-1, ]                                    # drop the header-in-data row
  t$block <- t$c1
  for (i in seq_len(nrow(t)))                     # forward-fill the block label
    if (is.na(t$block[i]) || !nzchar(trimws(t$block[i])))
      t$block[i] <- if (i > 1) t$block[i - 1] else NA_character_
  t %>% filter(grepl(want_block, block, ignore.case = TRUE) |
                 grepl("^Phylum:", trimws(c2))) %>%
    mutate(taxon = trimws(sub("^Phylum:\\s*", "", c2)),
           ext_direction = arrow_dir(c3),
           is_phylum = grepl("^Phylum", block, ignore.case = TRUE) |
                       grepl("^Phylum:", trimws(c2))) %>%
    filter(is_phylum, !is.na(ext_direction), nzchar(taxon)) %>%
    select(taxon, ext_direction)
}

age_ph  <- parse_tbl("maintext_Table1_associations.tsv", "c1", "c2", "c3", "^Phylum$")
smk_ph  <- parse_tbl("maintext_Table2_associations.tsv", "c1", "c2", "c3", "^Smokers$")
asay("published phylum rows: age ", nrow(age_ph), ", smoking ", nrow(smk_ph))

# SLL2 uses current phylum nomenclature; NHANES's SILVA release uses the older
# names. Bridge them explicitly - an unmapped name must NOT silently become a
# non-match that looks like "not tested".
PHYLUM_SYNONYM <- c(Patescibacteria = "Saccharibacteria", Bacteroidota = "Bacteroidetes",
                    Bacillota = "Firmicutes", Pseudomonadota = "Proteobacteria",
                    Fusobacteriota = "Fusobacteria", Actinomycetota = "Actinobacteria",
                    Synergistota = "Synergistetes", Spirochaetota = "Spirochaetae")
harmonise_phylum <- function(x) ifelse(x %in% names(PHYLUM_SYNONYM),
                                       unname(PHYLUM_SYNONYM[x]), x)

# --- Derive the NHANES phylum direction by sign vote -------------------------
vote <- function(v) {
  prim %>% filter(nhanes_variable == v) %>%
    mutate(phylum = m12_phylum_std(phylum)) %>%
    count(phylum, direction) %>%
    tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
    mutate(n_up = if ("up" %in% names(.)) up else 0L,
           n_down = if ("down" %in% names(.)) down else 0L,
           n_genera = n_up + n_down,
           nhanes_direction = dplyr::case_when(n_up > n_down ~ "up",
                                               n_down > n_up ~ "down",
                                               TRUE ~ NA_character_),
           vote = sprintf("%d up / %d down", n_up, n_down),
           unanimous = (n_up == 0 | n_down == 0)) %>%
    select(phylum, n_genera, n_up, n_down, nhanes_direction, vote, unanimous)
}

build <- function(pub, v, label, cohort_label) {
  nh <- vote(v)
  pub %>%
    mutate(phylum_published = taxon, phylum = harmonise_phylum(taxon)) %>%
    left_join(nh, by = "phylum") %>%
    transmute(nhanes_variable = v, phenotype_label = label,
              external_dataset = "WillisJR_2022_SacaLaLengua2", cohort = cohort_label,
              phylum_published, phylum_harmonised = phylum,
              nhanes_n_genera = n_genera, nhanes_vote = vote,
              nhanes_direction, nhanes_vote_unanimous = unanimous,
              ext_direction,
              tested = !is.na(nhanes_direction),
              same_direction = ifelse(is.na(nhanes_direction), NA,
                                      nhanes_direction == ext_direction),
              validation_classification = dplyr::case_when(
                is.na(nhanes_direction) & is.na(n_genera) ~ "NOT_TESTABLE_TAXON",
                is.na(nhanes_direction) ~ "MODEL_NOT_COMPARABLE",
                nhanes_direction == ext_direction ~ "DIRECTIONALLY_CONCORDANT_NOT_SIGNIFICANT",
                TRUE ~ "DIRECTIONALLY_DISCORDANT"),
              evidence_level = "phylum-level direction (NHANES side is a sign vote over genus associations, not a phylum regression)")
}

ph <- bind_rows(
  build(age_ph, "RIDAGEYR", "Age (per year)", paste0(COHORT_LABELS[["sll2_age"]], " phylum rows")),
  build(smk_ph, "SMQ_current_ever_never", "Smoking (current vs never)",
        paste0(COHORT_LABELS[["sll2_smoke"]], " phylum rows")))
write_tsv(ph, file.path(AUD, "phylum_level_validation.tsv"))

t <- ph %>% filter(tested)
b <- m12_binom(sum(t$same_direction, na.rm = TRUE), nrow(t))
asay("")
asay("=== phylum-level validation ===")
print(as.data.frame(ph %>% select(phenotype_label, phylum_published, phylum_harmonised,
                                  nhanes_vote, nhanes_direction, ext_direction,
                                  same_direction, validation_classification)))
asay("")
asay(sprintf("testable %d of %d published phylum rows; concordant %d/%d (%s), binomial p = %s",
             nrow(t), nrow(ph), b$k, b$n,
             ifelse(is.na(b$pct), "n/a", sprintf("%.0f%%", b$pct)),
             ifelse(is.na(b$p_value), "n/a", format.pval(b$p_value, digits = 2))))
