# Terry plot standalone pipeline — data, I/O, and engineer contract

Reviewer-response Terry figures plot each retained association as a point at its regression **estimate** on the horizontal axis against categorical **rows** (host-variable descriptions or genus-based OTU groupings) on the vertical axis, with vertical jitter only (zero horizontal jitter), point **fill** mapping host-variable prefix or phylum depending on plot orientation, **size** binned by Benjamini–Hochberg FDR (discrete thresholds from 0.05 down to 1×10^-6), and **color** encoding positive versus negative effect; Stage 1 supplies rows meeting `str_starts(term, independent_var)`, `std.error ≠ 0`, OTU membership in a ≥1% non-zero prevalence list evaluated with a **9847**-participant denominator in this script, non-missing FDR and q both ≤0.05, and oral/outcome binary filters where applicable. **Comprehensive** PDFs under `figures_out/comprehensive/<none|hellinger|clr>/` (and matching **indices** subtrees) show every row passing those gates for each Rmd-derived panel configuration, with optional prefix-based subsetting and up to **four** `geom_text_repel` annotations per row ranked by |estimate| and point size in the default path. **Compact** PDFs under `figures_out/compact/...` apply, within each rendered row stratum, at most **five** β>0 and **five** β<0 points chosen by descending |estimate| with ties broken by `p.value.fdr` then `otu_clean`, where the stratum is **Genus** (or `NA; otu_clean` if genus is missing/unclassified) for OTU-facing panels, **host_variable_column** for standard host panels, or the **(host_prefix, row_label)** pair for prefix-faceted host panels; **indices** compact figures instead take a single global slice of **five** positive and **five** negative estimates. Compact mode uses larger default axis text (**7** pt versus **6** pt on the horizontal axis and on host-variable y-axes), with **genus-bearing** text one point smaller in compact only: OTU Terry y-axis (genus rows), OTU indices y-axis (`parsed_genus`), and host Terry `geom_text_repel` genus labels (mm size reduced by one typographic pt equivalent); stronger vertical jitter; a higher effective annotation cap so essentially all retained points can receive labels; and host-annotation filtering that drops NA cleaned genus strings and incertae/sedis-like names; all outputs are isolated under `4_otu_host_plot_out_REVIEWER_EDIT` separately for each transformation.

This document is the **single source of truth** for paths, inputs, outputs, object shapes, and runtime environment for:

| Role | Path | Editable? |
|------|------|-----------|
| **Original (frozen)** | `1_otu_host_and_indicies_intermediate_processing.R` | **Never edit.** |
| **Original (frozen)** | `2_otu_host_n_indicies_plotting_n_saving.R` | **Never edit.** |
| **Reviewer / engineer** | `1_otu_host_and_indicies_intermediate_processing_REVIEWER_EDIT.R` | Yes — only here for Stage 1 changes. |
| **Reviewer / engineer** | `2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R` | Yes — only here for Stage 2 changes. |

High-level workflow (unchanged from originals):

1. **Stage 1** builds serialised R objects (mostly `.rds`) used by Stage 2.
2. **Stage 2** reads those objects and writes **PDF** Terry plots (host, OTU, and `*_indices`).

Manuscript-oriented narrative also exists in `0_data_description_method_result_explanation.md`; **this** file focuses on **engineering**: exact directories, filenames, column expectations, and **non-overwrite rules** for `*_REVIEWER_EDIT.R`.

---

## Reproducibility: conda environment and YAML

The **same** O2 conda prefix used for NHANES analyses (networks, Terry plots, etc.) is:

- `/n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis`

A **prefix export** is committed at (refresh after `conda install`):

- [`envs/nhanes-analysis_for_reviewers.yml`](../../envs/nhanes-analysis_for_reviewers.yml)

**Typical O2 session:**

```bash
cd /n/groups/patel/terry/nhanes_oral_mirco_cho
module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis
export R_LIBS="$CONDA_PREFIX/lib/R/library"
```

**Reviewer Terry pipeline (Stage 1 then Stage 2, vanilla R):**

```bash
Rscript --vanilla scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing_REVIEWER_EDIT.R
Rscript --vanilla scripts/4_otu_host_plotting_standalone/2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R
```

**Regenerate the YAML** after changing the env:

```bash
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis
conda env export --prefix /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis --no-builds \
  > envs/nhanes-analysis_for_reviewers.yml
```

**PDF font (Stage 2):** `Sys.getenv("NHANES_TERRY_PDF_FONT", "sans")` sets the ggplot base family (`sans` avoids Arial metric issues on some clusters; set to `ArialMT` when Arial is fully registered via `extrafont`).

---

## Reviewer `_REVIEWER_EDIT`: filtering sequence and plot semantics

The subsections below list **every gate in order**, what is **removed**, and what appears **on the PDF**. Terry plots are **not** the same object as network plots: they are **coefficient / effect-size** panels (rows × horizontal effect axis), not pairwise graphs.

### Stage 1 — `1_otu_host_and_indicies_intermediate_processing_REVIEWER_EDIT.R`

**Output root:** `results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/intermediate/`.

**SQLite:** `nhanes_oral_transformed_complete_processed.sqlite` if present and non-empty; otherwise `nhanes_oral_transformed_complete.sqlite` (console message).

**Load upstream WAS results:** for each of `demoWAS`, `oradWAS`, `exWAS`, `pheWAS`, `outWAS` and transforms `none`, `clr`, `hellinger`, read `results/<prefix>_out/result_<trans>/<prefix>_<trans>_tidied_complete.rds` (missing files → empty tibble).

**Derive `otu` on each tidied table (before significance filtering):**

| Analysis groups | `otu` source | `otu_generic` |
|-----------------|--------------|----------------|
| `demoWAS`, `exWAS` | `dependent_var` | `sub("_relative$", "", dependent_var)` |
| `oradWAS`, `pheWAS`, `outWAS` | `independent_var` | same on `independent_var` |

**Binary case-count whitelists (oradWAS / outWAS only):**

1. Variables from `configs/2_oradWAS_vars.txt` / `configs/5_outWAS_vars.txt`.
2. Count cases (`1`, `"1"`, `TRUE`) among participants in `DADA2RSV_GENUS_RELATIVE_F/G` SEQN lists, using `OralDisease_F/G` or `d_outcome_mcq_f/g`.
3. **Pass** if `cases_count > 9847 * 0.005` (literal **9847** in script ≈ 0.5% threshold).  
   **Removed:** binary outcomes with too few cases.

**OTU prevalence whitelist (global):**

4. Stack genus relative abundance F+G tables; count non-zero samples per OTU column.
5. **Pass** if `non_zero_count > 9847 * 0.01` and `otu != "SEQN"` (literal **9847**).  
   **Removed:** taxa with ≤ 1% non-zero prevalence on this fixed denominator.

> **Note (vs network reviewer Stage 1):** network scripts use `nrow(rsv_genus_relative)` for the same 1% and 0.5% **fractions**; Terry reviewer Stage 1 hard-codes **9847** in those inequalities. If cohort sizes ever diverge, prevalence/case lists could differ slightly between pipelines.

**Per `was_specs` row (each `*_tidied` block) — sequence:**

6. Start from `processed_tidied[[base]]`; skip if empty.
7. `filter(str_starts(term, independent_var))` — term matches host-side variable naming.  
   **Removed:** other model terms.
8. `filter(std.error != 0)` — zero standard errors dropped (implicitly allows `NA`? only rows with `std.error != 0`; `NA != 0` is `NA` and dplyr **drops** those rows).
9. `mutate(dependent_var = map_chr(dependent_var, ~ as.character(.x)[1]))` — scalarize list-like cells.
10. **outWAS\***: `filter(dependent_var %in% outWAS_case_counts_pass_list)`. **oradWAS\***: `filter(dependent_var %in% oradWAS_case_counts_pass_list)`.  
    **Removed:** binary outcomes failing case rule.
11. `filter(otu %in% otu_pass_prevalance_list)`.
12. `filter(!is.na(p.value.fdr), p.value.fdr <= 0.05)`.
13. `filter(!is.na(q.value), q.value <= 0.05)`.
14. `arrange(p.value.fdr)` then **select** the columns stored in `significant_results` / `datasets_for_plots`.

**What Stage 1 does not do:** it does **not** apply Terry Stage 2-only rules (prefix subsets, `distinct(otu_clean, estimate)`, compact top-5, annotation suppression). Those are strictly Stage 2.

---

### Stage 2 — `2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`

**Reads:** only `4_otu_host_plot_out_REVIEWER_EDIT/intermediate/*.rds`.  
**Writes:** `figures_out/<compact|comprehensive>/<none|hellinger|clr>/...` and `figures_out/<mode>/indices/<trans>/...`.

**Constants (reviewer):** `COMPACT_MAX_POS = 5`, `COMPACT_MAX_NEG = 5` define the compact **emphasis** set per group (see below); **`ANNOTATION_MAX_POS = 2`**, **`ANNOTATION_MAX_NEG = 2`** (`geom_text_repel` labels). **`JITTER_H_COMPACT = 0.18 × 1.1`** (10% more than the prior 0.18) vs `JITTER_H_COMPREHENSIVE = 0.09`. **`AXIS_TEXT_PT = 6`**, **`POINT_SIZE_SCALE = 1`**, Terry repel **`LABEL_SIZE_MM = 1.8`** (compact repel uses **`LABEL_SIZE_MM_COMPACT`** = one typographic pt smaller in mm). **`ggsave(..., pointsize = 5, family = "ArialMT", colormodel = "rgb")`** matches the frozen Stage 2 script for PDF overlays. **`COMPACT_HIGHLIGHT_ALPHA = 0.65`**, **`COMPACT_FADE_ALPHA = 0.1`** for Terry compact point drawing.

**Compact Terry points (not indices):** Compact **keeps all** association rows (same filtered table as comprehensive). **`mark_assoc_compact_highlight`** flags the top **5** positive and top **5** negative `estimate` rows per group (same grouping as the legacy `subset_assoc_rows_for_compact`); those points use **`COMPACT_HIGHLIGHT_ALPHA`**; all other real associations use **`COMPACT_FADE_ALPHA`**. Dummy padding rows stay at highlight alpha. `*_indices` compact PDFs still use **row subsetting** only (`subset_indices_df_compact`).

**Compact y-axis row order (Terry host/OTU plots):** Y-axis `row_label` order matches the comprehensive figure. Association counts and `dummy_top_bottom` ordering come from the full table (`df2_core` / `df_full_dummy`); factor levels intersect with labels present in the plot data (`terry_factor_row_order`, `terry_host_apply_split_y_order`).

#### Shared pipeline for each Terry config (`plot_and_save`) — order

1. **`prepare_dataset`:** Infer orientation: if `sum(dependent_var == otu) > sum(independent_var == otu)` then microbe is `dependent_var`, host is `independent_var`; else swap. Add `otu_clean = str_remove(otu, "_relative$")`.  
   **Effect:** every row has consistent `host_variable_column` / `otu_column` for plotting.

2. **Optional `prefix_subsets`:** If the plot config lists NHANES variable-name prefixes, **keep** rows where `host_variable_column` starts with any listed prefix (case-insensitive).  
   **Removed:** associations outside the requested domain buckets.

3. **`distinct(otu_clean, estimate, .keep_all = TRUE)`** — one row per OTU×effect pair (first row wins duplicates).  
   **Removed:** duplicate keys after join prep.

4. **Join** `genus_mapping` on `otu_clean` = `otu` (distinct on genus mapping side).

5. **`dot_size = calc_dot_size(p.value.fdr)`** from FDR magnitude bins; **`filter(!is.na(dot_size))`**.  
   **Removed:** rows with FDR outside the binned mapping (should not happen after Stage 1, but enforced here).

6. **Optional `remove_na_unclassified`:** drop `Genus == "unclassified"` or `NA`.

7. **Compact Terry only — emphasis flags:** If `fig_mode == "compact"`, apply **`mark_assoc_compact_highlight`** (uses **`top_pos_neg_slice`** with `COMPACT_MAX_POS` / `COMPACT_MAX_NEG`): **OTU** and **host (no split)** on `df2_core`; **host + `split_by_prefix`** after **`terry_host_join_meta`**. Sets logical **`compact_highlight`**; **no rows are dropped** from the plot table.

8. **Dummy padding rows** (constant label string) are **appended** for layout consistency; they are excluded from annotation tables by filtering `row_label` matching `DUMMY_`.

**What is drawn (Terry host / OTU PDFs):**

- **Y-axis:** one row per `row_label` (genus or `NA; otu` for OTU plots; host variable or facet label for host plots), ordered by association count / configuration, plus dummy spacer rows.
- **X-axis:** model **`estimate`** (effect size); symmetric limits often set from config `xlim`.
- **Points:** one per retained association (compact Terry: **all** associations, with alpha by emphasis); **vertical jitter only** (`width = 0` in position jitter) so x remains the true estimate.
- **Point size:** maps to binned **FDR** (`calc_dot_size`).
- **Color:** **Positive** vs **negative** `estimate` (`association_colors`).
- **Optional text annotations:** `create_annotations` builds repelled `geom_text_repel` labels. Per `row_label`, **`annotation_pos_neg_slice`** keeps up to **2** positive and **2** negative rows (ranked by `\|estimate\|`, **`dot_size`**, **`otu_clean`**). **Host** plots drop cleaned genus `NA` / **incertae** / **sedis**. **OTU** plots label `host_prefix (estimate)`.

**Indices PDFs (`plot_indices` / compact branch):**

- Same upstream `prepare_dataset` / prefix / joins as appropriate for indices configs.
- **Compact:** `subset_indices_df_compact` = global **`top_pos_neg_slice`** on the prepared indices data frame (up to **5** positive + **5** negative by `abs(estimate)`, tie-break `p.value.fdr`, `otu_clean`) — **not** grouped by facet key unless the indices path groups differently (code uses whole-df slice). **PDF height** uses the **same** y-axis row count as comprehensive (`indices_*_compact_dim_row_count` on the full pre-subset `df`) so compact indices pages match comprehensive dimensions while plotting the subset.
- Output filenames: `<file_prefix>_host_indices.pdf` / `_otu_indices.pdf` under `figures_out/<mode>/indices/<trans>/`.

**What differs in compact Terry vs comprehensive:** same rows pass steps 1–6; compact adds **`compact_highlight`** and **`point_alpha`** (no association-row drops). **Indices compact** still drops rows beyond the top 5+5 global slice (`subset_indices_df_compact`).

---

## Critical warnings (read before any edit)

### 1. Do not modify the original scripts

- **`1_otu_host_and_indicies_intermediate_processing.R`** and **`2_otu_host_n_indicies_plotting_n_saving.R`** are the reproducible baseline tied to published / frozen outputs.
- All reviewer-driven or experimental work happens in the `*_REVIEWER_EDIT.R` pair only.

### 2. Never overwrite original pipeline outputs

The originals write under:

```text
<repo_root>/results/analyses_results/4_otu_host_plot_out/
```

**Engineering rule for `*_REVIEWER_EDIT.R`:** every **write** path (intermediates, figures, supplementary artefacts, logs, temp exports) must live under a **separate root** whose path **includes the token `REVIEWER_EDIT`** in a way that cannot collide with `4_otu_host_plot_out/`.

**Implemented canonical layout for `*_REVIEWER_EDIT.R`:**

```text
<repo_root>/results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/
├── intermediate/          # all Stage-1 .rds / .csv (Stage 2 reads ONLY this path)
├── figures_out/
│   ├── compact/             # main-figure layer: per rendered row, keep top 5 positive + top 5 negative associations
│   │   ├── none/
│   │   ├── hellinger/
│   │   ├── clr/
│   │   └── indices/
│   │       ├── none/
│   │       ├── hellinger/
│   │       └── clr/
│   ├── comprehensive/       # supplementary-style: full filtered set (same filters; no compact down-selection)
│   │   ├── none/
│   │   ├── hellinger/
│   │   ├── clr/
│   │   └── indices/
│   │       ├── none/
│   │       ├── hellinger/
│   │       └── clr/
│   └── (no direct PDFs under figures_out root)
└── supplementary/
```

**Reads** (SQLite, upstream WAS `.rds`, phyloseq objects, configs, GOLD CSV, Rmd for plot config extraction) match the originals; only **write** paths use `4_otu_host_plot_out_REVIEWER_EDIT/`.

### 3. Compact vs comprehensive (Stage 2)

- **Comprehensive:** uses the full set of association rows that pass Stage-1 filtering and in-plot prefix/NA rules.
- **Compact:** same inputs, but within each rendered row keeps at most **5** associations with `estimate > 0` and **5** with `estimate < 0`, ranked by `abs(estimate)` (ties broken by `p.value.fdr`, `otu_clean`). Rows are defined by the y-axis grouping: OTU panels group by **Genus** (or `NA; <otu>`), host panels group by **host variable** (and if facetted by `host_prefix`, the rule applies within each facet-row). Slightly larger axis text (`7` pt vs `6` pt).

### 4. Mandatory workflow (every edit to `*_REVIEWER_EDIT.R`)

1. Re-read **this** document before changing code paths.
2. Make a **small** logical edit; do not stack unverified changes.
3. **Grep** reviewer scripts: no runtime `save*` target under `4_otu_host_plot_out/` (without `REVIEWER_EDIT`). Confirm frozen originals are untouched (`git status`).
4. **Run** Stage 1 then Stage 2 top-to-bottom on O2 (or the project `nhanes-analysis` conda env) with repo root as cwd.
5. Append **`0_data_description_method_result_explanation_REVIEWER_EDIT.md`** with what changed and what was generated.
6. If I/O or behaviour changed, update **this** file in the same session.

### 5. Significance rule (Stage 1 reviewer script)

Terry-ready rows require **non-missing** Benjamini–Hochberg **FDR** (`p.value.fdr <= 0.05`) and Storey **q** (`q.value <= 0.05`), in addition to prevalence and case-balance filters (see code).

### 6. Fonts (Stage 2)

Script 2 uses **`extrafont`** and expects Arial TTFs under:

```text
<repo_root>/data/fonts/
```

(pattern `arial.*\.ttf$` per `extrafont::font_import`). If fonts are missing, PDFs may still render but with fallback warnings.

---

## Repository root resolution

Both scripts define:

```r
resolve_base_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", args[grepl(file_arg, args)])
  if (length(script_path) == 0) {
    return(normalizePath("."))
  }
  normalizePath(file.path(dirname(script_path), "..", ".."))
}
```

So when launched as:

```bash
Rscript scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing.R
```

`base_path` is the **repository root** (two levels up from `scripts/4_otu_host_plotting_standalone/`).

If `Rscript` is invoked **without** `--file=...` in `commandArgs`, `base_path` falls back to **`normalizePath(".")`** (current working directory). **Always run from repo root** or use explicit `Rscript path/to/script.R` so `--file=` is present.

---

## Stage 1 — original: inputs (`1_otu_host_and_indicies_intermediate_processing.R`)

All paths below are relative to `base_path` unless noted.

### Upstream data sources

| Kind | Path | Notes |
|------|------|--------|
| **SQLite (NHANES transformed)** | `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` | **Hard fail** if missing (`stop()`). Used for metadata tables, case counts, prevalence, RSV tables. **`1_*_REVIEWER_EDIT.R` only:** if this file is missing **or 0 bytes**, falls back to `nhanes_oral_transformed_complete.sqlite` with a console message. |
| **Phyloseq objects (genus)** | `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/` | Required files (each `readRDS`): `ubiome_counts.rds`, `ubiome_relative.rds`, `ubiome_relative_none.rds`, `ubiome_relative_clr.rds`, `ubiome_relative_hellinger.rds`, `ubiome_relative_lognorm.rds`. |
| **GOLD genus mapping** | `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv` | `readr::read_csv`; columns documented below. |
| **WAS aggregated results** | `results/<N>_<name>WAS_out/result_<trans>/` | See next subsection. |
| **Config variable lists** | `configs/1_demoWAS_vars.txt` … `configs/5_outWAS_vars.txt` | `readLines`; also `2_oradWAS_vars.txt`, `5_outWAS_vars.txt` for binary case filters. |

### WAS result directory pattern

For each analysis `prefix` and transformation `trans`:

```text
results/<prefix>_out/result_<trans>/<prefix>_<trans>_tidied_complete.rds
results/<prefix>_out/result_<trans>/<prefix>_<trans>_glanced_complete.rds
results/<prefix>_out/result_<trans>/<prefix>_<trans>_rsq_complete.rds
```

Where `analysis_specs` defines:

| Analysis key | `prefix` | `trans` values |
|----------------|----------|----------------|
| demoWAS | `1_demoWAS` | none, clr, hellinger |
| oradWAS | `2_oradWAS` | none, clr, hellinger |
| exWAS | `3_exWAS` | none, clr, hellinger |
| pheWAS | `4_pheWAS` | none, clr, hellinger |
| outWAS | `5_outWAS` | none, clr, hellinger |

RDS keys in memory: `"<analysis>_<trans>"` (e.g. `demoWAS_none`). Missing files become **empty tibbles** via `load_result()`.

### SQLite tables referenced (non-exhaustive but script-critical)

- `table_names_epcf`, `variable_names_epcf` — filtered to F/G cycle tables and special outcome tables (see script: `years_F_or_G` construction).
- `DADA2RSV_GENUS_RELATIVE_F`, `DADA2RSV_GENUS_RELATIVE_G` — SEQN + genus columns for prevalence and genus SEQN lists.
- `OralDisease_F` / `OralDisease_G` (or `d_outcome_mcq_f` / `g`) — binary case counts for filters.

Metadata objects saved to RDS retain **whatever columns** exist in those EPC tables after `collect()`; Stage 1 saves:

- `table_description.rds` — from `table_names_epcf` filtered.
- `variable_description.rds` — from `variable_names_epcf` filtered.

Downstream joins use **`Variable.Name`**, **`Variable.Description`** from `variable_description` when building `ubiome_variable_mapping`.

### GOLD CSV: `ubiome_genus_mapping_complete.csv`

Header (single line, comma-separated):

```text
parsed_genus,otu,Genus,BIOTIC_RELATIONSHIPS,OXYGEN_REQUIREMENT,METABOLISM,ENERGY_SOURCES,GRAM_STAIN,SAMPLE_COLLECTION_DATE,SALINITY,CELL_SHAPE,MOTILITY,SPORULATION,TEMPERATURE_RANGE,STAIN_INDEX,OXYGEN_INDEX,MOTILITY_INDEX,SPORULATION_INDEX
```

Stage 1 saves a copy as `ubiome_genus_mapping_complete.rds` under intermediate. Stage 2 indices plots merge on **`otu`** (and use index columns **`STAIN_INDEX`**, **`OXYGEN_INDEX`**, **`MOTILITY_INDEX`**, **`SPORULATION_INDEX`**, plus **`parsed_genus`** for OTU-side indices logic).

---

## Stage 1 — original: outputs

**Root (original only):**

```text
results/analyses_results/4_otu_host_plot_out/
├── intermediate/     # primary Stage-1 artefacts
├── figures_out/      # created early; PDFs written by Stage 2, not Stage 1
└── supplementary/    # directory created; optional future use in baseline
```

### Intermediate `.rds` / `.csv` written by Stage 1

| File | Approximate content |
|------|---------------------|
| `table_description.rds` | `data.table` from EPC table metadata. |
| `variable_description.rds` | `data.table` from EPC variable metadata. |
| `ubiome_genus_mapping_complete.rds` | Full GOLD mapping as loaded from CSV. |
| `phyloseq_objects.rds` | Named list of six phyloseq objects (keys: `ubiome_counts`, `ubiome_relative`, …). |
| `tidied_results_raw.rds` | Named list of tidied tibbles per analysis×transform. |
| `glanced_results_raw.rds` | Named list of glanced tibbles. |
| `rsq_results_raw.rds` | Named list of rsq tibbles. |
| `tidied_results_with_otu.rds` | Same as tidied, plus columns **`otu`**, **`otu_generic`** (orientation rules in script). |
| `binary_case_filters.rds` | List: `oradWAS_case_counts`, `outWAS_case_counts`, pass lists, `binary_case_threshold`. |
| `prevalence_filters.rds` | List: `otu_non_zero` (long: `otu`, `non_zero_count`), `otu_pass_prevalance_list`, `prevalence_threshold`. |
| `significant_results.rds` | Named list `*_sig_res` per WAS×transform; each table columns below. |
| `datasets_for_plots.rds` | List of 18 tibbles (demowas/oradwas/… × transform), keys like `demowas_none`, … — **same rows as** corresponding `*_sig_res`. |
| `taxonomy_annotations.rds` | List: `phylum_info` (`otu`, `Phylum`), `genus_mapping` (`otu`, `Genus`) from phyloseq tax table. |
| `ubiome_variable_mapping.rds` | `data.table`: **`var_name`**, **`var_description`**. |
| `ubiome_variable_mapping.csv` | Same as RDS, `data.table::fwrite`. |
| `bucket_definitions.rds` | Named list of character vectors (host-variable buckets for prefixes / facets). |
| `plot_configs_none.rds` | Parsed from Rmd — list of plot configs for Terry plots. |
| `plot_configs_hellinger.rds` | Same. |
| `plot_configs_clr.rds` | Same. |
| `indices_plot_configs_none.rds` | Parsed from Rmd — list for indices PDFs. |
| `indices_plot_configs_hellinger.rds` | Same. |
| `indices_plot_configs_clr.rds` | Same. |

### `significant_results` / `datasets_for_plots` row columns

After filtering, significant association tables retain these **selected** columns (explicit in Stage 1 script):

- `otu`, `term`, `dependent_var`, `independent_var`, `estimate`, `statistic`, `std.error`, `p.value`, `p.value.fdr`, `q.value`, `otu_generic`, `available_cycles`, `effect_scale`, `interpretation_note`, `fdr_corrected`

**Filters (summary):** term starts with `independent_var`; `std.error != 0`; outcome-specific case lists for outWAS/oradWAS; OTU in prevalence pass list; `p.value.fdr` and `q.value` both non-NA and `< 0.05`.

---

## Stage 2 — original: inputs (`2_otu_host_n_indicies_plotting_n_saving.R`)

**Intermediate directory (original):**

```text
results/analyses_results/4_otu_host_plot_out/intermediate/
```

**Required files** (`ensure_exists()` — hard stop if any missing):

| File |
|------|
| `datasets_for_plots.rds` |
| `taxonomy_annotations.rds` |
| `ubiome_variable_mapping.rds` |
| `bucket_definitions.rds` |
| `ubiome_genus_mapping_complete.rds` |
| `plot_configs_none.rds`, `plot_configs_hellinger.rds`, `plot_configs_clr.rds` |
| `indices_plot_configs_none.rds`, `indices_plot_configs_hellinger.rds`, `indices_plot_configs_clr.rds` |

`bucket_definitions` is injected with `list2env(bucket_definitions, envir = environment())` so bare names (e.g. `energy_vars`) resolve inside plotting functions.

---

## Stage 2 — original: outputs (PDFs only)

**Figures root (original):**

```text
results/analyses_results/4_otu_host_plot_out/figures_out/
```

**Subdirectories created at start:**

- `figures_out/none/`, `figures_out/hellinger/`, `figures_out/clr/`
- `figures_out/indices/none/`, `figures_out/indices/hellinger/`, `figures_out/indices/clr/`

### Terry plot (host / OTU) PDF naming

`save_plot_object()` passes `file.path(transformation, config$file_prefix)` into `plot_and_save()`, which builds:

```text
figures_out/<none|hellinger|clr>/<file_prefix>.pdf
```

`file_prefix` and `plot_type` come **entirely** from the parsed Rmd `plot_configs_*` lists (each list element is typically a list with at least `df`, `plot_type`, `file_prefix`, plus optional knobs like `width_mm`, `prefix_subsets`, etc.).

### Indices PDF naming

Indices loop uses:

```text
file.path("indices", "<none|hellinger|clr>", config$file_prefix)
```

then `plot_indices()` appends:

- **`_host_indices.pdf`** when `plot_type == "host"`
- **`_otu_indices.pdf`** when `plot_type == "otu"`

Full pattern:

```text
figures_out/indices/<none|hellinger|clr>/<file_prefix>_host_indices.pdf
figures_out/indices/<none|hellinger|clr>/<file_prefix>_otu_indices.pdf
```

### REVIEWER_EDIT Stage 2 (`2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`)

- **Reads** intermediates from: `4_otu_host_plot_out_REVIEWER_EDIT/intermediate/` only.
- **Terry (host/OTU) PDFs:**

  ```text
  figures_out/<compact|comprehensive>/<none|hellinger|clr>/<file_prefix>.pdf
  ```

- **Indices PDFs:**

  ```text
  figures_out/<compact|comprehensive>/indices/<none|hellinger|clr>/<file_prefix>_host_indices.pdf
  figures_out/<compact|comprehensive>/indices/<none|hellinger|clr>/<file_prefix>_otu_indices.pdf
  ```

**Device:** PDF via `ggplot2::ggsave(..., device = "pdf", family = "ArialMT", ...)`.

---

## R packages loaded

### Stage 1

`data.table`, `dplyr`, `tidyr`, `tibble`, `purrr`, `readr`, `stringr`, `DBI`, `RSQLite`, `phyloseq`

### Stage 2

Stage 1 list **plus**: `forcats`, `ggplot2`, `egg`, `grid`, `gridExtra`, `extrafont` (the **REVIEWER_EDIT** Stage 2 script does not attach `ggrepel`; association annotations use `geom_text` for PDF stability with **ggplot2 4**).

Ensure the conda environment used on O2 includes **RSQLite**, **phyloseq** (e.g. `mamba install -c conda-forge -c bioconda bioconductor-phyloseq`), **extrafont**, **egg**, **gridExtra**, **ggrepel** (optional if you restore repel layers locally).

---

## O2 / SLURM — recommended runtime sequence

Activation order and reviewer `Rscript` lines are in [Reproducibility: conda environment and YAML](#reproducibility-conda-environment-and-yaml). Repeat here for copy-paste with an interactive allocation:

```bash
# Request resources as needed, e.g.:
# srun --pty -p interactive -t 4:00:00 --mem=32G bash

module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis
cd /n/groups/patel/terry/nhanes_oral_mirco_cho
export R_LIBS="$CONDA_PREFIX/lib/R/library"
```

**Order:** `module purge` → **gcc** → **miniforge conda** → **`eval "$(conda shell.bash hook)"`** → **`conda activate`** → **`export R_LIBS`** → **`cd` repo root** → **`Rscript --vanilla ...`**.

**Memory / time:** Stage 1 is I/O- and SQLite-heavy; Stage 2 is CPU + disk for many PDFs. Interactive `32G` / a few hours is a reasonable starting point; increase for full HPC batch reruns.

### Example SLURM batch skeleton

```bash
#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-08:00
#SBATCH --mem=32G
#SBATCH -p short

set -euo pipefail
module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis

cd /n/groups/patel/terry/nhanes_oral_mirco_cho
export R_LIBS="$CONDA_PREFIX/lib/R/library"
Rscript --vanilla scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing_REVIEWER_EDIT.R
Rscript --vanilla scripts/4_otu_host_plotting_standalone/2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R
```

(Adjust partition/time to your cluster policy.)

---

## Engineer checklist — before merging `*_REVIEWER_EDIT.R` changes

- [ ] **Output root** includes `REVIEWER_EDIT` and is **not** `results/analyses_results/4_otu_host_plot_out/` (unless you use a dedicated subfolder that cannot clash — sibling `4_otu_host_plot_out_REVIEWER_EDIT/` is preferred).
- [ ] **Stage 2** `intermediate_files_path` points to the **reviewer Stage-1** intermediate directory (if you regenerated intermediates); never silently read original intermediates while writing reviewer PDFs to a mixed path.
- [ ] No `saveRDS`, `fwrite`, `ggsave`, or `write*` targets under the original `4_otu_host_plot_out/` tree.
- [ ] Upstream **read-only** inputs unchanged unless the task explicitly documents a new SQLite / WAS path.
- [ ] `git status` shows no accidental edits to **`1_otu_host_and_indicies_intermediate_processing.R`** or **`2_otu_host_n_indicies_plotting_n_saving.R`**.

---

## Quick reference — path strings

| Purpose | Path |
|---------|------|
| **Frozen** Stage 1 & 2 output base | `results/analyses_results/4_otu_host_plot_out` |
| **Reviewer** Stage 1 & 2 output base | `results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT` |
| **Reviewer** intermediates | `.../4_otu_host_plot_out_REVIEWER_EDIT/intermediate` |
| **Reviewer** figures | `.../4_otu_host_plot_out_REVIEWER_EDIT/figures_out/{compact,comprehensive}/...` |
| SQLite DB | `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` |
| Phyloseq dir | `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate` |
| GOLD CSV | `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv` |
| Config dir | `configs` |
| Plot list source Rmd | `scripts/4_association_phyloseq_analyses/4_association_phyloseq_analyses.Rmd` |
| Fonts | `data/fonts` (Arial TTFs) |

---

## Document history

- **Created** as engineer I/O and safety contract for `*_REVIEWER_EDIT.R` work; originals remain the reference implementation for behaviour and science logic.
- **Updated** to reflect implemented `4_otu_host_plot_out_REVIEWER_EDIT` writes, `compact` / `comprehensive` figure trees, mandatory per-edit workflow, and explicit FDR + q significance lines in the reviewer Stage 1 script.
- **2026-04-24:** Added conda YAML path (`envs/nhanes-analysis_for_reviewers.yml`), full **reviewer** filter/plot semantics (Stage 1 sequence, Stage 2 pipeline, compact vs comprehensive, Terry vs network prevalence denominator note), and `NHANES_TERRY_PDF_FONT`.
