-# Microbial Signature Correlation Heatmaps
+
+Pipeline for rebuilding the microbial signature correlation heatmaps outside the
+original R Markdown. Everything lives in `scripts/7_microbial_signature_heatmap`
+and writes to `results/analyses_results/7_microbial_signature_heatmap_out/`.
+
+## Workflow
+
+1. **Build intermediates**
+   ```bash
+   Rscript scripts/7_microbial_signature_heatmap/1_correlation_signature_matrix_intermediate_processing.R
+   ```
+   - Reads the NHANES SQLite database under `data/00_nhanes_omp_transformed_db/`
+     plus the aggregated svyglm outputs.
+   - Applies prevalence, case-count, and FDR/q filters.
+   - Saves per-analysis signature matrices, weighted-correlation matrices,
+     summaries, and `heatmap_significant_results.rds` into
+     `…/intermediate/microbial_signature/`.
+
+2. **Render legacy-sized heatmaps**
+   ```bash
+   Rscript scripts/7_microbial_signature_heatmap/2_microbial_signature_heatmap_sized.R
+   ```
+   - Consumes only the intermediates above.
+   - Recreates the historical “sized” heatmaps (currently CLR only—lognorm
+     hooks remain but require matching intermediates).
+   - Outputs `{analysis}_{transformation}_all_results_conventional_heatmap_sized.pdf`
+     plus an updated `dataset_summary.csv`.
+
+3. **Render optimal-k heatmaps (silhouette-driven splits)**
+   ```bash
+   Rscript scripts/7_microbial_signature_heatmap/2_microbial_signature_heatmap_sized_cutree_optimal_k_sigres.R \
+     otu_threshold=2 fdr_threshold=0.01 prevalence_threshold=0.01
+   ```
+   - Draws only from intermediates.
+   - Enforces configurable OTU/FDR/prevalence filters using
+     `heatmap_significant_results.rds`.
+   - Saves PDFs prefixed `otu<xx>_fdr<yy>_prev<zz>_…` plus a matching summary CSV.
+   - CLI flags: `otu_threshold`, `fdr_threshold`, `prevalence_threshold`,
+     `k_min`, `k_max`.
+
+## Outputs
+
+- `results/…/intermediate/microbial_signature/**` – RDS matrices, CSV summaries,
+  merged `allWAS` resources, verification tables.
+- `results/…/*.pdf` – sized legacy heatmaps.
+- `results/…/otu*_weighted_correlation_legend.pdf` – shared legend.
+- `results/…/otu*_*.pdf` – optimal-k heatmaps (and placeholder PDFs when filters
+  remove all variables).
+- `results/…/otu*_dataset_summary.csv` – per-run stats including selected k.
+
+## Core logic
+
+- Signature rows contain svyglm β estimates + standard errors.
+- Weighted correlations use inverse-variance weights
+  `w = 1/(se₁² + se₂²)` with effective sample size
+  `(Σw)² / Σw²`, Benjamini-Hochberg FDR, and optional duplicate-row removal.
+- Heatmaps size each cell by FDR tier and cluster with `1 - |r|` distances.
+
+## Dependencies
+
+`data.table`, `dplyr`, `tidyr`, `purrr`, `ComplexHeatmap`, `circlize`,
+`ggplot2`, `grid`, `cluster`, `RSQLite`, `DBI`, `extrafont` (Arial). Load fonts
+before first plot (`extrafont::loadfonts(device = "pdf")`).
