# Reviewer response — Terry / OTU / host plots (REVIEWER_EDIT pipeline)

This document is the **reviewer-facing** companion to the engineering I/O spec (`xx_documentation_of_data_files_i_o_for_engineers.md`) and the baseline manuscript narrative (`0_data_description_method_result_explanation.md`).

**Update it after every successful end-to-end run** of the reviewer scripts: date, branch/commit if applicable, code changes, and new or changed outputs under `results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/`.

---

## 2026-04-29 — Overlay parity: typography, PDF device, indices page height

### Code

- **Stage 2** (`2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`): Single **`AXIS_TEXT_PT = 6`** for compact and comprehensive (matches frozen script). **`POINT_SIZE_SCALE = 1`**, **`LABEL_SIZE_MM = 1.8`**. **`ggsave`** adds **`family = "ArialMT"`**, **`colormodel = "rgb"`**, compact uses **`pointsize = 5`** (same as comprehensive). Compact indices PDF **`height_mm`** uses row count from the **full** association table before **`subset_indices_df_compact`** (`indices_host_compact_dim_row_count`, `indices_otu_compact_dim_row_count`). Compact Terry repel uses **`LABEL_SIZE_MM_COMPACT`** (one pt smaller in mm than comprehensive).

### Outputs

Re-run Stage 2 to refresh all PDFs.

---

## 2026-04-29 — Compact Terry: full table + emphasis alpha + jitter

### Code

- **Stage 2** (`2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`): Compact Terry plots **no longer drop** associations outside the top five per sign; **`mark_assoc_compact_highlight`** flags emphasis rows (**`COMPACT_HIGHLIGHT_ALPHA = 0.65`**), others **`COMPACT_FADE_ALPHA = 0.1`**. **`JITTER_H_COMPACT`** increased by **10%** (was 0.18). **`geom_point`** uses **`aes(alpha = point_alpha)`** + **`scale_alpha_identity`**. Indices compact unchanged (`subset_indices_df_compact`). **`subset_assoc_rows_for_compact`** retained for reference / indices-related use but Terry compact paths use marking instead of subsetting.

### Outputs

Re-run Stage 2. Engineer I/O doc updated.

---

## 2026-04-29 — Terry compact y-axis order matches comprehensive

### Code

- **Stage 2** (`2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`): Compact Terry plots reuse **y-axis `row_label` order** derived from the **same comprehensive (full filtered) association table** before compact point subsetting—helpers `terry_otu_*`, `terry_host_*`, `terry_factor_row_order`. Saving comprehensive PDFs first is **not** required; order is computed from `df2_core` on every call.

### Outputs

Re-run Stage 2 to refresh Terry PDFs.

---

## 2026-04-29 — Terry annotations: top two per sign per row

### Code

- **Stage 2** (`2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`): Association **`geom_text_repel`** labels now use **`ANNOTATION_MAX_POS = 2`** and **`ANNOTATION_MAX_NEG = 2`** via **`annotation_pos_neg_slice`** inside **`create_annotations`** (ranked by `\|estimate\|`, then FDR dot size, then `otu_clean`). Applies to **both** comprehensive and compact PDFs. Removed the old **`min_annotations`** path (formerly four strongest overall per row, or “label all” in compact).

### Outputs

Re-run Stage 2 to refresh Terry PDFs. Stage 1 unchanged.

---

## 2026-04-29 — Compact figures: genus text one point smaller

### Code

- **Stage 2** (`2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R`): In `figures_out/compact/...` only, genus on the **y-axis** uses **one fewer pt** than the compact base (OTU Terry plots and OTU indices panels). **Host** Terry plots keep host-variable y-axis size unchanged but shrink **`geom_text_repel`** cleaned-genus labels by **one pt** in mm (`LABEL_SIZE_MM - 25.4/72`). Constants: `COMPACT_GENUS_AXIS_PT_DELTA`, `LABEL_SIZE_MM_COMPACT`, `PT_TO_MM`.
- **Stage 1** header documents that typography is owned by Stage 2. **Engineer I/O** doc updated in the compact-mode paragraph.

### Outputs

Re-run Stage 2 (`Rscript --vanilla ...2_..._REVIEWER_EDIT.R`) to refresh all compact (and shared) PDFs; Stage 1 unchanged unless intermediates are missing.

---

## 2026-04-24 — End-to-end O2 run: conda stack, SQLite fallback, Stage 2 PDF rendering

### Environment and commands

- **Conda env:** `envs/.conda/envs/nhanes-analysis`. **`phyloseq`** was installed with `mamba install -c conda-forge -c bioconda bioconductor-phyloseq` (this also brought **R 4.5.1** and dependencies in line with that solve). **Stage 2** additionally used `mamba install r-egg r-gridextra r-extrafont` (`r-ggrepel` is optional; the current **REVIEWER_EDIT** Stage 2 script does not `library(ggrepel)`).
- **R library path:** from the project root, use a **single** library tree and **ignore** repo `.Renviron` mixing:  
  `export R_LIBS="$CONDA_PREFIX/lib/R/library"` and run **`Rscript --vanilla`** for both stages.
- **Inputs:** if `nhanes_oral_transformed_complete_processed.sqlite` is **missing or 0 bytes**, Stage 1 **REVIEWER_EDIT** automatically uses `nhanes_oral_transformed_complete.sqlite` (message printed to the log).

### Code adjustments this cycle (REVIEWER_EDIT only)

- **Stage 1:** SQLite fallback as above.
- **Stage 2:** association annotations use `ggrepel::geom_text_repel()` again (labels with leader lines; tuned `box.padding`, `point.padding`, etc.). **Default PDF font** is **sans** via `Sys.getenv("NHANES_TERRY_PDF_FONT", "sans")`; set `NHANES_TERRY_PDF_FONT=ArialMT` to try embedded Arial. **Size legend:** simplified `scale_size_continuous()` (no custom six-step FDR labels) for **ggplot2 4** compatibility (`breaks`/`labels` length check). **Guards** for `xlim` when `|estimate|` is zero/empty, and **minimum** PDF width/height in `save_pdf_figure()`.

### Outputs verified

- **Stage 1** log: `4_otu_host_plot_out_REVIEWER_EDIT/stage1_run.log` (~12 min in this run).
- **Stage 2** log: `4_otu_host_plot_out_REVIEWER_EDIT/stage2_run.log` (~10 min in this run; **728** `*.pdf` under `figures_out/`).
- Warnings: many `mbcsToSbcs` / Unicode substitution notes for variable descriptions in PDFs (en-dashes, etc.); no abort after Stage 2 completion.

### How to re-run (recommended)

```bash
module purge
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
eval "$(conda shell.bash hook)"
conda activate /n/groups/patel/terry/nhanes_oral_mirco_cho/envs/.conda/envs/nhanes-analysis
cd /n/groups/patel/terry/nhanes_oral_mirco_cho
export R_LIBS="$CONDA_PREFIX/lib/R/library"
Rscript --vanilla scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing_REVIEWER_EDIT.R 2>&1 | tee results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/stage1_run.log
Rscript --vanilla scripts/4_otu_host_plotting_standalone/2_otu_host_n_indicies_plotting_n_saving_REVIEWER_EDIT.R 2>&1 | tee results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/stage2_run.log
```

---

## 2026-04-23 — Initial implementation of isolated outputs and compact / comprehensive layers

### Scope

- **In scope:** standalone Terry, host, OTU, and `*_indices` PDFs from `*_REVIEWER_EDIT.R` only.
- **Out of scope:** network plots; frozen originals `1_otu_host_and_indicies_intermediate_processing.R` and `2_otu_host_n_indicies_plotting_n_saving.R` (unchanged).

### Reviewer concern addressed (design intent)

- Main association panels were too dense for primary figures. The pipeline now emits two **presentation** layers from the **same** regression-filtered association tables (not a re-analysis):
  - **Compact** — for main-figure use: **per rendered row**, keep at most 5 positive and 5 negative associations by largest `|estimate|` (tie-breaks: `p.value.fdr`, `otu_clean`). For facetted host figures (`split_by_prefix`), the rule applies within each facet-row.
  - **Comprehensive** — full filtered set (same rules as before this layer) for supplement-style use.

### Significance rule (unchanged science, made explicit in code)

Stage 1 (reviewer) keeps association rows with:

- Benjamini–Hochberg FDR: `p.value.fdr` present and `<= 0.05`
- Storey q: `q.value` present and `<= 0.05`

(Plus existing prevalence, case-balance, and term filters as in the original pipeline.)

### Output locations (all under `REVIEWER_EDIT` root)

- **Stage 1 intermediates:** `results/analyses_results/4_otu_host_plot_out_REVIEWER_EDIT/intermediate/`
- **Terry / host / OTU PDFs:** `.../figures_out/{compact,comprehensive}/{none,hellinger,clr}/<name>.pdf`
- **Indices PDFs:** `.../figures_out/{compact,comprehensive}/indices/{none,hellinger,clr}/<name>_{host,otu}_indices.pdf`

### What changed from the original R scripts (behavioural)

- All **writes** use `4_otu_host_plot_out_REVIEWER_EDIT` instead of `4_otu_host_plot_out`.
- Stage 2 runs each plot **twice** (comprehensive then compact). Terry compact uses alpha emphasis and jitter vs comprehensive; **axis text is 6 pt for both** (same as the frozen Stage 2 script). Indices compact subsets plotted rows but **PDF height** matches comprehensive via dimension helpers.
- **Indices** use the same compact row-subset rule on the prepared association `df` before building index panels.

### How to re-run (O2 / project conda)

See **2026-04-24** section above for the current recommended command block (`R_LIBS`, `Rscript --vanilla`, `tee` logs).

---

## Template for future entries

| Field | Value |
|-------|--------|
| Date | YYYY-MM-DD |
| Branch / commit |  |
| Code / doc touched |  |
| Reviewer point |  |
| Method / logic |  |
| New or updated PDF paths |  |

Add a new subsection (or new dated section) for each material change; keep a short entry even when no plots were regenerated (e.g. comment-only edits), noting “no figure regeneration this cycle”.
