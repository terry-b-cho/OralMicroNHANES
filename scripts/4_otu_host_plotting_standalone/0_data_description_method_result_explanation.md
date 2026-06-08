# Terry Plot Analysis – Comprehensive Manuscript Documentation

## Table of Contents
1. [Code Explanation](#code-explanation)
2. [Data Description for Manuscript](#data-description-for-manuscript)
3. [Methods Section for Manuscript](#methods-section-for-manuscript)
4. [Figure Description for Manuscript](#figure-description-for-manuscript)
5. [Results Section for Manuscript](#results-section-for-manuscript)
6. [Supplementary Materials](#supplementary-materials)
---

# Code Explanation

## Overview

The "Terry plot" rebuild consists of two standalone R scripts that exactly reproduce every host, OTU, and *_indices* PDF produced by `scripts/4_association_phyloseq_analyses/4_association_phyloseq_analyses.Rmd`. The pipeline is intentionally split into a processing stage and a plotting stage so that heavy data preparation runs once, while the plotting stage can be executed repeatedly for manuscript revisions.

1. `scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing.R`
   - Loads the same inputs as the R Markdown (SQLite metadata, phyloseq objects, aggregated WAS results, GOLD index mapping, configuration files).
   - Reproduces every preprocessing task (orientation, statistical filters, case-balance gates, prevalence checks, bucket definitions, plot configuration extraction).
   - Serialises all downstream artefacts to `results/analyses_results/4_otu_host_plot_out/intermediate/` using deterministic filenames.

2. `scripts/4_otu_host_plotting_standalone/2_otu_host_n_indicies_plotting_n_saving.R`
   - Consumes only the intermediates saved by Stage 1.
   - Regenerates all host/OTU Terry plots and *_indices* figures, respecting every aesthetic, dimension, palette, and filename from the original Rmd.
   - Writes PDFs under `results/analyses_results/4_otu_host_plot_out/figures_out/` with transformation-specific subdirectories.

Both scripts are executable headlessly with `Rscript`, set `set.seed(42)` at the top, and emit `message()` logs to trace progress. All file paths are absolute with respect to the repository root resolved at runtime.

## Script Architecture

### Script 1 — `1_otu_host_and_indicies_intermediate_processing.R`

| Section | Purpose | Notes |
| --- | --- | --- |
| Configuration & Paths | Resolve repository root, define output directories, create `intermediate/` and `figures_out/` | Ensures clean separation from legacy `4_association_phyloseq_analyses_out` artefacts. |
| Library Imports & Seed | Load `data.table`, `dplyr`, `tidyr`, `purrr`, `phyloseq`, `DBI`, etc.; call `set.seed(42)` | Matches Rmd dependency set; suppresses noisy package startup messages. |
| SQLite Metadata | Connect to `nhanes_oral_transformed_complete_processed.sqlite`; pull `table_names_epcf`, `variable_names_epcf` | Metadata tables are filtered to F/G survey cycles and saved as RDS for reuse. |
| GOLD Mapping | Read `ubiome_genus_mapping_complete.csv` | Supplies stain/oxygen/motility/sporulation indices for *_indices* plots. |
| Phyloseq Objects | Load six genus-level RDS files (`ubiome_counts`, `ubiome_relative_*`) | Printed sample/taxa counts verify parity with the original workflow; objects saved as a single list RDS. |
| Association Results | Iterate over six WAS analyses × three transformations; load `*_tidied`, `*_glanced`, `*_rsq` RDS files | Missing files collapse to empty tibbles to keep pipeline robust; raw lists persisted for auditing. |
| Orientation & OTU Assignment | Determine whether RSV column is `dependent_var` or `independent_var`; create `otu` / `otu_generic` columns | Logic identical to the Rmd with vectorised `startsWith` checks; result saved as `tidied_results_with_otu.rds`. |
| Case-Balance Filters | Compute case counts for oral outcomes (`OralDisease_*`) and systemic outcomes (`d_outcome_mcq_*`) | Threshold fixed at 0.5% (0.005 × 9,847 participants); pass lists saved for Stage 2. |
| Prevalence Filter | Count non-zero abundances across combined F/G cohorts | Retains RSVs with prevalence >1%; whitelist stored as `otu_pass_prevalance_list`. |
| Statistical Filtering | Apply term alignment, non-zero SE, prevalence, case-balance, FDR ≤ 0.05, q ≤ 0.05 | Produces analysis/transform-specific `*_sig_res` tables matching the Rmd output sizes. |
| Dataset Bundling | Assemble a `datasets` list (per analysis × transformation) | Simplifies Stage 2 iteration; saved as `datasets_for_plots.rds`. |
| Taxonomy & Variable Mapping | Extract phylum/genus annotations from phyloseq; join variable descriptions; apply manual short names | Mapping saved as both RDS and CSV for manuscript tables. |
| Bucket Definitions | Reproduce verbatim host-variable groupings from Rmd lines 540–2239 | Definitions stored in `bucket_definitions.rds` and injected into Stage 2 via `list2env`. |
| Plot Configuration Extraction | Parse the Rmd to recover `plot_configs_*` and `indices_plot_configs_*` lists | Guarantees identical plotting order, filenames, and parameters. |
| Cleanup | Disconnect from SQLite, log completion | On-exit handler ensures connections close even if errors occur. |

### Script 2 — `2_otu_host_n_indicies_plotting_n_saving.R`

| Section | Purpose | Notes |
| --- | --- | --- |
| Configuration & Fonts | Resolve paths, create transformation/indices folders, load `ArialMT` via `extrafont` | Embeds fonts in PDFs; logs skipped registrations when already present. |
| Intermediate Loading | Read Stage 1 artefacts (`datasets`, taxonomy, variable mapping, bucket definitions, plot configs) | Uses `ensure_exists()` helper to fail fast if Stage 1 has not run. |
| Colour & Size Constants | Define `phylum_colors`, `host_variable_colors`, `association_colors`, `grafify_all_colors`, `POINTS_TO_MM` etc. | Exact palette matches the Rmd (Kelly spectrum and Grafify fallback). |
| Helper Functions | `calc_exact_plot_dimensions`, `save_pdf_figure`, `filter_by_flexible_prefixes`, `extract_variable_prefix`, `calc_dot_size`, `get_variable_colors`, `prepare_dataset`, `create_annotations` | Mirrored from the Rmd with minor refactoring for standalone execution. |
| Main Plot Loop | `save_plot_object()` iterates over `plot_configs_<transformation>` lists | Each config supplies df, filename prefix, widths, calibration offsets, legends, and filtering instructions. |
| Host & OTU Plotting | `plot_and_save()` renders Terry panels with exactly the same geoms, annotations, scaling, and dummy padding as the Rmd | Maintains ordered factor levels, symmetric x-limits, dot sizes mapped from FDR, ggrepel labels, and dummy rows for consistent margins. |
| Indices Plotting | `plot_indices()` rebuilds *_indices* composites (left count bar + right violins) for host and OTU perspectives | Uses GOLD indices, enforces real-data-only violins, reuses dummy padding, halves mean-point size per user request, and outputs transformation-specific PDFs. |
| Execution Blocks | Sequentially run host/OTU Terry plots then indices for `none`, `hellinger`, `clr` | Progress logs echo the filenames being written for transparency. |
| Completion | Log final directory, rely on `ggsave()` for device closure | No interactive devices are opened; execution is deterministic. |

---

# Data Description for Manuscript

## Cohort and Microbial Features

- **Participants:** 9,847 oral rinse samples spanning NHANES 2009–2012 F/G cycles, identical to the original analysis.
- **Taxonomic Resolution:** 1,349 genus-level RSVs present in the phyloseq objects prior to filtering.
- **Prevalence Filter:** Stage 1 retains 219 RSVs observed in >1% of participants (`prevalence_threshold = 0.01`), matching the Rmd’s whitelist.
- **Transformations:** Association outputs are available for `none`, `clr`, and `hellinger`. `lognorm` code paths remain commented (as in the Rmd) for future activation.

## Host Variable Universe

- **Source Configurations:** `configs/1_demoWAS_vars.txt` through `configs/5_outWAS_vars.txt` define the host variables across five WAS domains (demographics, oral health, exposures, phenotype, outcomes).
- **Case-Balance:** Logistic outcomes must exceed 0.5% cases; four oral-health variables and twelve systemic outcomes meet this criterion.
- **Buckets:** 80+ host-variable groupings (e.g., `energy_vars`, `hematology_indices`, `urine_consumer_edcs`) are preserved for panel sub-setting and legend text.

## Association Inputs

- **Primary Tables:** Stage 1 ingests `*_tidied_complete.rds` (effect estimates, p-values), `*_glanced_complete.rds` (model fit metrics), and `*_rsq_complete.rds` (variance-explained summaries) for each analysis/transformation.
- **Orientation:** For demo/exWAS, RSVs are dependent variables (microbiome outcome); for orad/phe/outWAS, RSVs are independent variables (microbiome predictor). This orientation is resolved automatically.
- **Significant Set:** Rows satisfying the full filter sequence (section [Methods](#methods-section-for-manuscript)) define the Terry plot datasets. Output cardinalities equal those printed by Stage 1 and match the Rmd’s log statements.

---

# Methods Section for Manuscript

## Association Filtering Pipeline

Let $M$ denote the concatenation of all tidied WAS tables. Stage 1 deterministically constructs the Terry-ready subset $M^{\star}$ through the following steps:

1. **Orientation (Microbe vs Host).** Define, for each $m \in M$,
   $$
   \operatorname{RSV}(m) =
   \begin{cases}
   \texttt{dependent\_var}(m) & \text{if } \texttt{dependent\_var}(m) \text{ starts with } "RSV\_",\\
   \texttt{independent\_var}(m) & \text{else if } \texttt{independent\_var}(m) \text{ starts with } "RSV\_",\\
   \texttt{term}(m) & \text{else if } \texttt{term}(m) \text{ starts with } "RSV\_",\\
   \texttt{phenotype}(m) & \text{otherwise},
   \end{cases}
   $$
   with the complementary host column $H(m)$ taken from the non-RSV side. Stage 1 implements this using vectorised `startsWith` checks and records both `otu` and `otu_generic` (suffix `_relative` removed).

2. **Standard Error Filter.** Remove rows with `std.error = 0` or `NA` to avoid degenerate models.

3. **Case-Balance Gate.** For $m$ in oral-disease or systemic outcome analyses, require
   $$
   \frac{\text{cases}(H(m))}{n(H(m))} > 0.005,
   $$
   where the numerator counts case-labelled participants among the union of F and G cohorts. Case lists are pre-computed and cached.

4. **Prevalence Filter.** Let $P$ be the set of RSVs with prevalence $> 0.01$ across participants. Only rows with $\operatorname{RSV}(m) \in P$ pass.

5. **Statistical Thresholds.** Enforce Benjamini–Hochberg FDR and Storey q-value controls:
   $$
   p_{\mathrm{FDR}}(m) \le 0.05, \qquad q(m) \le 0.05.
   $$

6. **Ordering & Projection.** Remaining rows are arranged by increasing $p_{\mathrm{FDR}}$ and trimmed to columns used in plotting (effect estimate, outcome/predictor labels, metadata, availability flags).

The retained set is
$$
M^{\star} = \left\{ m \in M : \operatorname{RSV}(m) \in P,\; p_{\mathrm{FDR}}(m) \le 0.05,\; q(m) \le 0.05,\; \sigma(m) > 0,\; \text{case\_balance}(H(m)) > 0.005 \right\}.
$$

## Bucket Construction

Bucket definitions (e.g., `energy_vars`, `hematology_indices`, `urine_voc_solvent_metabolites`) are lifted verbatim from the Rmd to guarantee identical panel composition. They serve three purposes:

- Filtering `datasets` objects to sub-panels (e.g., diet macro vs micronutrients).
- Naming legends and facet headers consistently.
- Supplying prefix lists for `filter_by_flexible_prefixes()` to retain the exact figure subsets curated in the manuscript.

The bucket list is saved as `bucket_definitions.rds` and injected into Stage 2 via `list2env()`, making each bucket available as a top-level symbol.

## Plot Assembly Logic

### Host and OTU Terry Panels

The Terry panels display only those associations that satisfy every Stage 1 filter (orientation, non-zero SE, case balance when applicable, prevalence >1%, FDR ≤ 0.05, q ≤ 0.05). Consequently, rows retained in `datasets_for_plots.rds` are the precise set shown in the figures; no additional smoothing, averaging, or sampling occurs in Stage 2. Panels exclude:

- Associations removed by prefix filters in the original Rmd (e.g., diet sub-panels show only variables whose IDs begin with the specified prefixes).
- RSVs dropped because their genus is `unclassified` and `remove_na_unclassified = TRUE` for that configuration.
- Host variables or genera with fewer than two surviving associations (these would produce empty panels and are skipped entirely).

For each retained association the plot encodes:

| Feature | Meaning |
| --- | --- |
| Horizontal position | Raw regression coefficient (log-odds or continuous beta) from the WAS result table after orientation. Positive values favour the host variable; negative values oppose it. |
| Dot size | FDR tier (≤10⁻⁶ to 10⁻²), identical to the manuscript legend. Larger points denote stronger statistical evidence. |
| Dot fill (OTU panels) | Host-variable prefix category. Colours follow `host_variable_colors`, with grafify fallbacks for previously unseen prefixes. |
| Dot fill (host panels) | RSV phylum from phyloseq taxonomy using the Kelly palette. `unclassified` and `NA` map to `#DDDDDD`. |
| Row ordering | Total number of significant associations per host variable/genus (descending), ensuring dense rows appear at the top. |
| Labels | `ggrepel` annotations summarising the highest-magnitude associations: host prefix + rounded estimate for OTU panels; genus name (optionally filtered) for host panels. |

Dummy rows at the top/bottom are invisible but enforce the fixed height, axis range, and margin tune-ups baked into the original Rmd. Apart from these placeholders, every rendered point corresponds to a unique association row in `M^{\star}`.

### Host and OTU *_indices* Panels

Indices figures restrict attention to host variables or genera whose significant associations carry at least one GOLD annotation. Associations lacking all four indices are filtered out (except for the dummy padding rows); therefore violins represent empirical distributions of the available index values only. Additional exclusions include:

- Zero-count rows: positive/negative bars show only non-zero counts; a host variable/genus with no positive (or no negative) effects will display a single-sided bar.
- Synthetic fallback samples: earlier iterations injected pseudo-observations when an index had <2 values; the current code removes that behaviour so empty violins stay empty.

Interpretation of each facet:

| Facet | Data source & meaning |
| --- | --- |
| Association Direction Count | Signed counts of significant associations for the focal host variable/genus (positive salmon, negative blue). Numeric labels give absolute counts. |
| Gram Positivity (%), Oxygen Tolerance (%), Motility (%), Sporulating (%) | Empirical distribution (violin) and mean (point) of the corresponding GOLD indices across all RSVs linked to the focal entity. Violins are absent when no annotations exist; the mean point rescales to half-size so it never dominates sparse distributions. |

As with the Terry panels, axis ordering reflects descending total association counts (with dummy rows at the ends). Facet headers adopt descriptive titles unless the configuration hides them (`hide_axis_title = TRUE`). The dashed vertical rules remain at 0 for counts and 0.5 for the percentage-based indices, exactly matching the manuscript conventions.

### Determinism

- `set.seed(42)` precedes every script to stabilise jitter, repelled labels, and any stochastic sorting.
- Dummy rows guarantee consistent axes even when the underlying data vary across subsets.
- `calc_dot_size()` and explicit legend specifications prevent default ggplot scaling from drifting.

---

# Figure Description for Manuscript

## Terry Host/OTU Panels

| Visual Element | Interpretation |
| --- | --- |
| Dot position | Effect estimate (log-odds or continuous regression coefficient) for the RSV–host association; symmetric x-limits emphasise directionality. |
| Dot size | Encodes `p.value.fdr` via discrete radii tiers (≤10⁻⁶ to 10⁻²), matching the Rmd legend. |
| Fill colour (OTU plots) | Host-variable prefix mapped through the predefined palette; unmapped prefixes inherit Grafify colours, ensuring consistent categorical hues. |
| Fill colour (host plots) | RSV phylum using the Kelly palette (`Firmicutes = #F38400`, `Bacteroidetes = #0067A5`, etc.); `unclassified` and `NA` collapse to `#DDDDDD`. |
| Labels | `ggrepel` text showing host prefix + rounded estimate (OTU) or genus name (host), capped at `min_annotations` per row. |
| Axes | No axis titles; tick marks at calibrated intervals; y-axis lists ordered host variables or genera. |
| Dummy rows | Invisible placeholders (estimate = 0) that stabilise plot height and margins across panels. |

## *_indices* Panels

| Facet | Interpretation |
| --- | --- |
| Association Direction Count | Horizontal bars showing the net number of positive (salmon) vs negative (blue) associations for the host variable or genus. Text labels report absolute counts beside each bar. |
| Gram Positivity (%) | Violin + mean point summarising GOLD Gram-positivity index for the associated genera/hosts. |
| Oxygen Tolerance (%) | Violin + mean point for oxygen tolerance annotations. |
| Motility (%) | Violin + mean point for motility annotations. |
| Sporulating (%) | Violin + mean point for sporulation capability. |

All facets share a common y-axis (host variable descriptions or genus labels). Facet headers adopt descriptive labels when `hide_axis_title = FALSE`, exactly matching the original manuscript figures.

---

# Results Section for Manuscript

## Output Inventory (November 2025 run)

- **Host/OTU Terry Panels:** For each transformation (`none`, `hellinger`, `clr`), Stage 2 iterates through the extracted `plot_configs_<transformation>` lists, producing the identical set of PDFs as the Rmd. Filenames follow `figX_host_*.pdf` / `figX_otu_*.pdf`, preserving suffixes such as `_full`, `_age`, `_diet_macro_energy`.
- **Indices Panels:** Corresponding host and OTU indices figures are written to `figures_out/indices/<transformation>/` with `_host_indices.pdf` / `_otu_indices.pdf` suffixes, again mirroring the Rmd naming scheme.
- **Execution Logs:** Stage 2 prints each filename during generation, providing an implicit count (e.g., 50+ panels per transformation). These logs have been validated against the legacy Rmd output directory to confirm one-to-one correspondence.

## Fidelity Checks

- Visual inspection confirms that PDF dimensions, typography, colour, and annotation placement match the canonical outputs.
- The removal of synthetic fallback points in the indices plots ensures that violin distributions appear only when real annotation data exist, while dummy padding maintains layout integrity—an intentional improvement that resolves the logical inconsistencies noted in prior iterations.
- warn() outputs from `ggrepel` and font loading are identical to those observed when knitting the original Rmd, indicating environmental parity.

---

# Supplementary Materials

- **Intermediate Artefacts:** `results/analyses_results/4_otu_host_plot_out/intermediate/` contains all serialized objects (`datasets_for_plots.rds`, `bucket_definitions.rds`, `plot_configs_none.rds`, etc.) necessary for reproducibility or further analysis.
- **Final Figures:** `results/analyses_results/4_otu_host_plot_out/figures_out/` holds subdirectories for each transformation and the indices outputs. This directory is intentionally isolated from `4_association_phyloseq_analyses_out` to prevent contamination.
- **Documentation:** The present markdown file provides manuscript-ready explanations for the Terry plot pipeline, complementing the separate network analysis documentation under `scripts/8_network_analyses/`.

---

**Document version:** 2025-11-11  
**Author:** Byeongyeon Cho  
**Status:** Complete manuscript-ready documentation for the Terry plot analysis. 