#!/usr/bin/env Rscript
# =============================================================================
# Module 12 Phase 3 — Step 1
# Inventory EVERY OralMicroNHANES association across all five WAS schemes
#
# Purpose : Phase 2 validated age and smoking. Those were the first two
#           implemented, not the maximum possible. This step enumerates every
#           phenotype in all five schemes with its FDR-significant,
#           prevalence-qualified genus association count, so no validation
#           opportunity can be missed merely because nobody looked.
#
#           Applies the project's primary reproducibility criteria:
#             NORMALIZATION == "clr"  (all five files are the clr variants)
#             FDR < 0.05
#             pooled NHANES genus prevalence > 10%
#           and preserves Storey q where present.
#
# Usage   : Rscript 01_inventory_nhanes_associations.R
# Outputs : audit/nhanes_variable_inventory.tsv   one row per variable
#           audit/nhanes_primary_associations.tsv one row per variable x genus
# Env     : R 4.5.1 (conda nhanes-analysis)
# =============================================================================
HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "scripts/12_external_data_validation/11_validation_audit"
source(file.path(HERE, "00_audit_config.R"))
asay("=== Step 1: inventory all five WAS schemes ===")

# --- Genus reference: mapping, prevalence, phylum ----------------------------
map <- read.csv(file.path(PROJECT_ROOT,
  "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv"),
  stringsAsFactors = FALSE, na.strings = c("NA", "")) %>% select(otu, parsed_genus)
ref <- m12_nhanes_genus_reference(PROJECT_ROOT) %>%
  select(genus, prevalence_overall, prevalence_tier)

suppressPackageStartupMessages(library(phyloseq))
ps_all <- readRDS(file.path(PROJECT_ROOT,
  "results/analyses_results/02_preprocess_db_n_phyloseq_out/intermediate/nhanes_phyloseq_objects_all.rds"))
nm <- intersect(c("counts", "relative", "ubiome_counts"), names(ps_all))
stopifnot(length(nm) >= 1)
tax <- as.data.frame(phyloseq::tax_table(ps_all[[nm[1]]]))
tax$otu <- rownames(tax)
phy <- tax %>% transmute(otu, phylum = m12_phylum_std(Phylum))

genus_meta <- map %>% left_join(phy, by = "otu") %>%
  left_join(ref, by = c("parsed_genus" = "genus"))
asay("genus reference: ", nrow(map), " OTUs, ", n_distinct(map$parsed_genus), " parsed genera")

# --- Walk each scheme --------------------------------------------------------
# The files are large (exWAS 501 MB, pheWAS 117 MB). Subset to main effects on
# read and gc() immediately; never hold two schemes at once.
inv <- list(); prim <- list()

for (i in seq_len(nrow(WAS_SPEC))) {
  sp <- WAS_SPEC[i, ]
  asay("")
  asay("--- ", sp$analysis_type, " (genus is the ", sp$genus_role, " variable) ---")
  d <- readRDS(file.path(PROJECT_ROOT, sp$rds_rel))

  # Main effects only. `fdr_corrected` is a "was corrected" flag, not a
  # significance flag - covariate rows carry an UNCORRECTED p.value.fdr, so
  # filtering on FDR before this subset returns garbage.
  d <- d %>% filter(fdr_corrected)
  if (!"q.value" %in% names(d)) d$q.value <- NA_real_

  # Resolve which column carries the genus, per this scheme's orientation.
  gcol <- if (sp$genus_role == "dependent") "dependent_var" else "independent_var"
  pcol <- if (sp$genus_role == "dependent") "independent_var" else "dependent_var"
  stopifnot(all(grepl("RSV_genus|_relative", head(d[[gcol]], 50))))

  d <- d %>%
    mutate(otu = sub("_relative$", "", .data[[gcol]]),
           nhanes_variable = .data[[pcol]],
           analysis_type = sp$analysis_type,
           genus_role = sp$genus_role) %>%
    left_join(genus_meta, by = "otu")

  # Drop unclassified genus names: they cannot be matched to any external
  # taxonomy and would inflate every denominator.
  n_before <- nrow(d)
  d <- d %>% filter(!is.na(parsed_genus), !grepl("^unclassified", parsed_genus),
                    parsed_genus != "Incertae_Sedis")
  asay("  main effects: ", format(n_before, big.mark = ","),
       " -> ", format(nrow(d), big.mark = ","), " with a named genus")

  d <- d %>%
    mutate(fdr_sig   = !is.na(p.value.fdr) & p.value.fdr < FDR_THRESHOLD,
           prev_ok   = !is.na(prevalence_overall) & prevalence_overall > PREVALENCE_THRESHOLD,
           primary   = fdr_sig & prev_ok)

  inv[[sp$analysis_type]] <- d %>%
    group_by(analysis_type, genus_role, nhanes_variable) %>%
    summarise(
      n_terms                  = n(),
      n_genus_tested           = n_distinct(parsed_genus),
      n_fdr_significant        = sum(fdr_sig),
      n_fdr_sig_prev_gt10      = n_distinct(parsed_genus[primary]),
      n_fdr_sig_prev_gt10_rsv  = sum(primary),
      n_up   = n_distinct(parsed_genus[primary & estimate > 0]),
      n_down = n_distinct(parsed_genus[primary & estimate < 0]),
      median_abs_beta = median(abs(estimate[primary]), na.rm = TRUE),
      min_fdr = suppressWarnings(min(p.value.fdr, na.rm = TRUE)),
      n_covariates = suppressWarnings(as.integer(names(sort(table(n_covariates), decreasing = TRUE))[1])),
      n_obs_median = median(n_obs, na.rm = TRUE),
      .groups = "drop")

  # The primary set is defined per (variable x GENUS), but several RSVs can
  # parse to the same genus name (Eubacterium is backed by 7). Collapse to one
  # row per genus, keeping the most significant RSV, and record how many RSVs
  # back it so the choice is visible rather than silent.
  prim[[sp$analysis_type]] <- d %>% filter(primary) %>%
    group_by(nhanes_variable, parsed_genus) %>%
    arrange(p.value.fdr, .by_group = TRUE) %>%
    mutate(n_rsv_backing = dplyr::n()) %>%
    slice(1) %>% ungroup() %>%
    transmute(analysis_type, genus_role, nhanes_variable, term,
              genus = parsed_genus, phylum, otu, n_rsv_backing,
              estimate, std.error, p.value, p.value.fdr, q.value,
              n_covariates, n_obs, prevalence_overall, prevalence_tier,
              direction = ifelse(estimate > 0, "up", "down"))
  rm(d); gc(verbose = FALSE)
}

inv_all  <- bind_rows(inv)  %>% arrange(desc(n_fdr_sig_prev_gt10))
prim_all <- bind_rows(prim)

write_tsv(inv_all,  file.path(AUD, "nhanes_variable_inventory.tsv"))
write_tsv(prim_all, file.path(AUD, "nhanes_primary_associations.tsv"))

asay("")
asay("=== scheme totals ===")
print(as.data.frame(inv_all %>% group_by(analysis_type, genus_role) %>%
  summarise(variables = n(),
            with_any_primary = sum(n_fdr_sig_prev_gt10 > 0),
            primary_associations = sum(n_fdr_sig_prev_gt10), .groups = "drop")))
asay("")
asay("PRIMARY VALIDATION SET: ", format(nrow(prim_all), big.mark = ","),
     " variable x genus associations (FDR<0.05 & prevalence>10%) across ",
     n_distinct(prim_all$nhanes_variable), " variables")
asay("")
asay("=== top 25 variables by primary association count ===")
print(as.data.frame(inv_all %>% filter(n_fdr_sig_prev_gt10 > 0) %>%
  select(analysis_type, nhanes_variable, n_genus_tested, n_fdr_significant,
         n_fdr_sig_prev_gt10, n_up, n_down, n_covariates) %>% head(25)))
