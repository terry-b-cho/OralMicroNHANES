#!/usr/bin/env Rscript

# =============================================================================
# MODULE 11 — STEP 04: variable availability + regression-schema documentation
# 04_variable_availability.R      (run AFTER 03_aggregate_per_year_results.sh)
# =============================================================================
# For every association variable (the non-microbiome side of each WAS type),
# classifies availability across cycles and documents the regression schema
# actually fitted, so that no missing-variable situation can ever be mislabeled
# as failed replication downstream.
#
# Availability classes:
#   F_AND_G                    fitted in both cycles with the same coding
#   F_ONLY / G_ONLY            variable exists (and was fitted) in one cycle only
#   UNAVAILABLE                in the pooled schema but produced no per-year
#                              main-effect result in either cycle
#   AVAILABLE_BUT_NOT_COMPARABLE  fitted in both cycles but with materially
#                              different coding (e.g. quartiles in one cycle,
#                              continuous in the other) -> excluded from direct
#                              beta comparison, with the reason recorded
#
# Outputs (all under OUT_ROOT/comparisons/variable_availability/ + tables/):
#   variable_availability.csv     one row per (analysis, variable)
#   regression_schema_table.csv   one row per (analysis, cycle, variable)
#   tables/table1_cycle_characteristics.csv
#   tables/table2_variable_availability.csv   (copy of the availability table)
# =============================================================================

source("scripts/11_per_year_association/00_module11_config.R")  # run from the repo root

suppressPackageStartupMessages({ library(DBI); library(RSQLite) })

cat("MODULE 11 STEP 04: variable availability + regression schema\n")

# The non-microbiome "association variable" per analysis:
#   1_demoWAS / 3_exWAS: microbiome is the DEPENDENT variable -> variable = exposure
#   2_oradWAS / 4_pheWAS / 5_outWAS: microbiome is the INDEPENDENT -> variable = dependent_var
var_col_for <- function(analysis) {
  if (analysis %in% c("1_demoWAS", "3_exWAS")) "exposure" else "dependent_var"
}

# Summarize the fitted coding of one variable in one cycle from its main-effect
# term structure (may vary across the 1,349 OTU pairs when binning falls back).
coding_summary <- function(d, varname) {
  if (nrow(d) == 0) return(list(coding = NA_character_, n_patterns = 0L))
  pat <- d %>%
    group_by(pair_id) %>%
    summarise(terms = paste(sort(unique(term)), collapse = "|"), .groups = "drop") %>%
    count(terms, sort = TRUE)
  lab_of <- function(terms) {
    tl <- strsplit(terms, "\\|")[[1]]
    if (identical(tl, varname)) return("continuous")
    if (all(grepl("Q[2-9]$", tl))) return(sprintf("quartile(%d contrasts)", length(tl)))
    if (all(grepl("T[2-9]$", tl))) return(sprintf("tertile(%d contrasts)", length(tl)))
    sprintf("categorical(%d levels: %s)", length(tl),
            paste(str_remove(tl, paste0("^", varname)), collapse = ","))
  }
  dominant <- lab_of(pat$terms[1])
  list(coding = dominant, n_patterns = nrow(pat))
}

avail_rows  <- list()
schema_rows <- list()

for (an in M11_ANALYSES) {
  vcol <- var_col_for(an)

  # Variables named in the per-cycle schemas (existence in the cycle's data)
  sch_var_col <- if (vcol == "exposure") "indep_var" else "dep_var"
  sch <- lapply(M11_CYCLES, function(cyc) {
    f <- file.path(M11_OUT, "schema_per_year", cyc,
                   paste0(an, "_", M11_NORM, "_schema_structure.csv"))
    unique(read_csv(f, show_col_types = FALSE)[[sch_var_col]])
  })
  names(sch) <- M11_CYCLES
  sch_pooled <- unique(read_csv(
    file.path(M11_BASE, "results/0_ss_files", paste0(an, "_", M11_NORM, "_schema_structure.csv")),
    show_col_types = FALSE)[[sch_var_col]])

  # Fitted main effects per cycle
  fitted <- lapply(M11_CYCLES, function(cyc) {
    d <- m11_load_main_effects(an, cyc)
    if (is.null(d)) stop("Missing aggregated per-year results for ", an, " cycle ", cyc,
                         " - run 03_aggregate_per_year_results.sh first.")
    d$pair_id <- paste(d$dependent_var, d$exposure, sep = "~")
    d$variable <- d[[vcol]]
    d
  })
  names(fitted) <- M11_CYCLES

  all_vars <- sort(unique(c(sch_pooled, sch$F, sch$G)))

  for (v in all_vars) {
    dF <- fitted$F %>% filter(variable == v)
    dG <- fitted$G %>% filter(variable == v)
    in_schema_F <- v %in% sch$F
    in_schema_G <- v %in% sch$G

    cF <- coding_summary(dF, v)
    cG <- coding_summary(dG, v)

    fitted_F <- nrow(dF) > 0
    fitted_G <- nrow(dG) > 0
    identical_coding <- fitted_F && fitted_G &&
      !is.na(cF$coding) && !is.na(cG$coding) && cF$coding == cG$coding

    cls <- if (fitted_F && fitted_G) {
      if (identical_coding) "F_AND_G" else "AVAILABLE_BUT_NOT_COMPARABLE"
    } else if (fitted_F) "F_ONLY"
    else if (fitted_G) "G_ONLY"
    else "UNAVAILABLE"

    reason <- dplyr::case_when(
      cls == "AVAILABLE_BUT_NOT_COMPARABLE" ~
        sprintf("dominant coding differs: F=%s vs G=%s", cF$coding, cG$coding),
      cls == "F_ONLY" & !in_schema_G ~ "variable not measured in cycle G",
      cls == "F_ONLY" &  in_schema_G ~ "in G schema but produced no fit in G",
      cls == "G_ONLY" & !in_schema_F ~ "variable not measured in cycle F",
      cls == "G_ONLY" &  in_schema_F ~ "in F schema but produced no fit in F",
      cls == "UNAVAILABLE" ~ "no per-year main-effect result in either cycle",
      TRUE ~ NA_character_)

    avail_rows[[length(avail_rows) + 1]] <- tibble(
      analysis = an, variable = v,
      availability = cls,
      in_schema_F = in_schema_F, in_schema_G = in_schema_G,
      in_schema_pooled = v %in% sch_pooled,
      fitted_F = fitted_F, fitted_G = fitted_G,
      coding_F = cF$coding, coding_G = cG$coding,
      n_coding_patterns_F = cF$n_patterns, n_coding_patterns_G = cG$n_patterns,
      identical_definition = identical_coding,
      direct_FG_comparison_valid = cls == "F_AND_G",
      reason_not_comparable = reason)

    for (cyc in M11_CYCLES) {
      d <- if (cyc == "F") dF else dG
      if (nrow(d) == 0) next
      dom_formula <- names(sort(table(d$formula_used), decreasing = TRUE))[1]
      schema_rows[[length(schema_rows) + 1]] <- tibble(
        analysis = an, cycle = cyc, variable = v,
        microbiome_role = if (vcol == "exposure") "dependent" else "independent",
        dominant_formula = dom_formula,
        n_formula_variants = length(unique(d$formula_used)),
        n_covariates_range = paste(range(d$n_covariates), collapse = "-"),
        survey_weight = "WTMEC2YR", strata = "SDMVSTRA", psu = "SDMVPSU",
        lonely_psu = "certainty", nest = TRUE,
        participant_n_min = min(d$n_obs), participant_n_max = max(d$n_obs),
        taxa_tested = length(unique(if (vcol == "exposure") d$dependent_var else d$exposure)),
        n_main_effect_tests = nrow(d))
    }
  }
  cat(sprintf("  %-10s: %d variables classified\n", an, length(all_vars)))
}

availability <- bind_rows(avail_rows) %>% arrange(analysis, availability, variable)
schema_tbl   <- bind_rows(schema_rows) %>% arrange(analysis, cycle, variable)

write_csv(availability, file.path(M11_AVAIL, "variable_availability.csv"))
write_csv(availability, file.path(M11_FIGURES_ACCESSORY, "optional_table2_variable_availability.csv"))
write_csv(schema_tbl,   file.path(M11_AVAIL, "regression_schema_table.csv"))

cat("\nAvailability class counts:\n")
print(count(availability, analysis, availability))

# ---- Table 1: cycle-level analytic characteristics --------------------------
prev <- m11_load_prevalence()
tax  <- m11_load_taxonomy()

t1 <- lapply(M11_CYCLES, function(cyc) {
  ncol_s <- prev[[paste0("n_samples_", cyc)]][1]
  pcol   <- prev[[paste0("prevalence_", cyc)]]
  all_n  <- bind_rows(lapply(M11_ANALYSES, function(an) {
    d <- m11_load_main_effects(an, cyc)
    tibble(analysis = an, n_min = min(d$n_obs), n_max = max(d$n_obs),
           n_tests = nrow(d))
  }))
  tibble(
    cycle = cyc,
    nhanes_years = if (cyc == "F") "2009-2010" else "2011-2012",
    microbiome_samples = ncol_s,
    genera_total = nrow(prev),
    genera_prevalence_gt10pct = sum(pcol > M11_PREV_THRESHOLD, na.rm = TRUE),
    regression_participant_n_range = paste(min(all_n$n_min), max(all_n$n_max), sep = "-"),
    total_main_effect_tests = sum(all_n$n_tests),
    survey_design = "svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE); survey.lonely.psu='certainty'",
    weight_note = "unscaled 2-year MEC weights (single-cycle NCHS design)",
    limitations = paste("~half the pooled sample size;",
                        sum(availability$availability == paste0(cyc, "_ONLY")),
                        "variables exclusive to this cycle"))
}) %>% bind_rows()

write_csv(t1, file.path(M11_TABLES, "table1_cycle_characteristics.csv"))
cat("\nTable 1 written. STEP 04 complete.\n")
