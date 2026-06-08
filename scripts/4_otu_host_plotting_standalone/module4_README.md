# Module 4 — Terry plots (host/OTU effect-size panels + indices)

## What this module does
Renders "Terry plot" PDFs from the WAS regression results — each retained association is a point at its regression `estimate` (x-axis) against a categorical row (host variable or genus) on the y-axis, with point size mapped to BH-FDR magnitude, fill mapped to phylum or host-variable prefix, and colour by effect direction. Two parallel pipelines:
- **Baseline** (`*.R`) — writes to `4_otu_host_plot_out/`
- **Additional** (`*_additional.R`) — writes to `4_otu_host_plot_out_additional/`. Adds dual-threshold significance gating (BH-FDR ≤ 0.05 AND Storey q ≤ 0.05), a compact figure variant (top 5 positive + 5 negative per row at full opacity, rest faded), and a `comprehensive` variant.

## Inputs (relative to `PROJECT_ROOT`)
- `data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete_processed.sqlite` *(`...complete.sqlite` fallback in the `_additional` pipeline)*
- `results/analyses_results/2_preprocess_db_n_phyloseq_out/intermediate/ubiome_{counts,relative,relative_none,relative_clr,relative_hellinger,relative_lognorm}.rds`
- `results/analyses_results/3_gold_db_microbial_phenotype_out/intermediate/ubiome_genus_mapping_complete.csv`
- `configs/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_vars.txt`
- `results/{1_demoWAS,2_oradWAS,3_exWAS,4_pheWAS,5_outWAS}_out/result_{none,clr,hellinger}/*_{tidied,glanced,rsq}_complete.rds`
- `data/fonts/arial*.ttf` *(for the ArialMT PDF family; falls back to a system sans family if absent)*

## Outputs

### Baseline pipeline → `results/analyses_results/4_otu_host_plot_out/`
- `intermediate/` — `*_results_raw.rds`, `tidied_results_with_otu.rds`, `significant_results.rds`, `datasets_for_plots.rds`, `phyloseq_objects.rds`, `taxonomy_annotations.rds`, `ubiome_variable_mapping.{rds,csv}`, `bucket_definitions.rds`, `plot_configs_{none,hellinger,clr}.rds`, `indices_plot_configs_{none,hellinger,clr}.rds`, `binary_case_filters.rds`, `prevalence_filters.rds`, `table_description.rds`, `variable_description.rds`, `ubiome_genus_mapping_complete.rds`
- `figures_out/{none,hellinger,clr}/<file_prefix>.pdf` — host & OTU Terry plot PDFs
- `figures_out/indices/{none,hellinger,clr}/<file_prefix>_{host,otu}_indices.pdf`

### Additional pipeline → `results/analyses_results/4_otu_host_plot_out_additional/`
Same intermediates as above plus:
- `figures_out/{compact,comprehensive}/{none,hellinger,clr}/<file_prefix>.pdf`
- `figures_out/{compact,comprehensive}/indices/{none,hellinger,clr}/<file_prefix>_{host,otu}_indices.pdf`

## Scripts
- `1_otu_host_and_indicies_intermediate_processing.R` / `1_otu_host_and_indicies_intermediate_processing_additional.R` — build intermediates (WAS result loading + prevalence/case-balance filters + significance gating + plot-config extraction)
- `2_otu_host_n_indicies_plotting_n_saving.R` / `2_otu_host_n_indicies_plotting_n_saving_additional.R` — render the Terry plot + indices PDFs

## Environment
R >= 4.5 with: `data.table`, `dplyr`, `tidyr`, `tibble`, `purrr`, `readr`, `stringr`, `forcats`, `DBI`, `RSQLite`, `phyloseq`, `ggplot2`, `ggrepel`, `egg`, `grid`, `gridExtra`, `extrafont`.

Conda spec: `envs/nhanes-analysis_for_reviewers.yml`.

## How to run
1. Open each `.R` file and set `PROJECT_ROOT` (top) to your local clone of this repository.
2. From the repo root — choose either pipeline:

```bash
# Baseline pipeline (writes to 4_otu_host_plot_out/)
Rscript scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing.R
Rscript scripts/4_otu_host_plotting_standalone/2_otu_host_n_indicies_plotting_n_saving.R

# Additional pipeline (writes to 4_otu_host_plot_out_additional/)
Rscript scripts/4_otu_host_plotting_standalone/1_otu_host_and_indicies_intermediate_processing_additional.R
Rscript scripts/4_otu_host_plotting_standalone/2_otu_host_n_indicies_plotting_n_saving_additional.R
```

## Run order
Within each pipeline: `1 → 2`. Script 1 builds all intermediates that script 2 reads from `intermediate/`. The baseline and additional pipelines are independent and write to separate output dirs.
