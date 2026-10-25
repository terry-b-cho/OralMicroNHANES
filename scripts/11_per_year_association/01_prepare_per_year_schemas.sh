#!/bin/bash

# =============================================================================
# MODULE 11 — STEP 01: prepare per-year (single-cycle) schema files, CLR only
# 01_prepare_per_year_schemas.sh
# =============================================================================
# Splits each existing CLR schema file in results/0_ss_files (READ-ONLY input,
# verified byte-identical to the production originals) by its `cycle` column
# (column 6: "F" or "G") into two single-cycle schema files.
#
# WHY: universal_was_analysis.R pools cycles iff a (dep_var, indep_var) pair has
# BOTH an F row and a G row in the schema. Feeding it a single-cycle schema makes
# every pair run in single_F / single_G mode automatically:
#   svydesign(ids=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
# i.e., the identical design/weights/covariates/formulas, just not pooled.
#
# Outputs (ONLY inside the module-11 output root — originals never touched):
#   results/analyses_results/11_per_year_association/schema_per_year/F/<analysis>_clr_schema_structure.csv
#   results/analyses_results/11_per_year_association/schema_per_year/G/<analysis>_clr_schema_structure.csv
# Filenames keep the <analysis>_<norm>_schema_structure.csv convention so the
# step-02 submitter can reuse the production naming logic unchanged.
#
# Resume-safe: existing output files are skipped, never overwritten.
# =============================================================================

set -u

# PROJECT_ROOT: path of your local clone; export it or run from the repo root.
BASE_DIR="${PROJECT_ROOT:-$(pwd)}"
SS_DIR="${BASE_DIR}/results/0_ss_files"                                        # READ-ONLY
OUT_ROOT="${BASE_DIR}/results/analyses_results/11_per_year_association"
SCHEMA_OUT="${OUT_ROOT}/schema_per_year"

ANALYSIS_TYPES=("1_demoWAS" "2_oradWAS" "3_exWAS" "4_pheWAS" "5_outWAS")
NORM="clr"

mkdir -p "${SCHEMA_OUT}/F" "${SCHEMA_OUT}/G"

echo "MODULE 11 STEP 01: splitting CLR schema files by cycle"
echo "  Input : $SS_DIR (read-only)"
echo "  Output: $SCHEMA_OUT/{F,G}/"
echo ""

fail=0
for analysis in "${ANALYSIS_TYPES[@]}"; do
    src="${SS_DIR}/${analysis}_${NORM}_schema_structure.csv"
    if [[ ! -f "$src" ]]; then
        echo "ERROR: missing source schema: $src"
        fail=1
        continue
    fi

    for cyc in F G; do
        dst="${SCHEMA_OUT}/${cyc}/${analysis}_${NORM}_schema_structure.csv"
        if [[ -f "$dst" ]]; then
            echo "  SKIP (exists): $dst"
            continue
        fi
        # Keep header; keep rows whose 6th field (cycle) matches.
        # Schema fields contain no embedded commas (verified: dep_var, dep_table,
        # indep_var, indep_table, n, cycle are all comma-free tokens).
        awk -F',' -v c="$cyc" 'NR==1 || $6==c' "$src" > "${dst}.tmp" && mv "${dst}.tmp" "$dst"
    done

    # Integrity check: data rows of F + G must equal data rows of the original
    n_src=$(( $(wc -l < "$src") - 1 ))
    n_f=$(( $(wc -l < "${SCHEMA_OUT}/F/${analysis}_${NORM}_schema_structure.csv") - 1 ))
    n_g=$(( $(wc -l < "${SCHEMA_OUT}/G/${analysis}_${NORM}_schema_structure.csv") - 1 ))
    if [[ $(( n_f + n_g )) -eq $n_src ]]; then
        printf "  OK   %-10s original=%8d  F=%8d  G=%8d  (F+G == original)\n" "$analysis" "$n_src" "$n_f" "$n_g"
    else
        printf "  FAIL %-10s original=%8d  F=%8d  G=%8d  (MISMATCH)\n" "$analysis" "$n_src" "$n_f" "$n_g"
        fail=1
    fi
done

echo ""
if [[ $fail -eq 0 ]]; then
    echo "STEP 01 complete: all per-year schema files ready under $SCHEMA_OUT"
    exit 0
else
    echo "STEP 01 FAILED for at least one analysis — see messages above"
    exit 1
fi
