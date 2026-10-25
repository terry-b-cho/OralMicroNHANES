#!/bin/bash
# =============================================================================
# Module 12 Phase 3 — validation audit, ordered driver
#
# Each step depends on the previous one, so the order is fixed and the script
# stops on the first error.
#
# Usage : bash scripts/12_external_data_validation/11_validation_audit/run_all.sh
# Env   : R 4.5.1 (conda nhanes-analysis); pdftoppm for the report PNGs
# =============================================================================
set -euo pipefail
PROJECT_ROOT="${ORALMICRO_ROOT:-$(pwd)}"
RS="${M12_RSCRIPT:-Rscript}"
HERE="${PROJECT_ROOT}/scripts/12_external_data_validation/11_validation_audit"
LOGS="${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/audit/logs"
mkdir -p "$LOGS"; cd "$PROJECT_ROOT"
command -v "$RS" >/dev/null || { echo "[FATAL] Rscript not found: $RS" >&2; exit 1; }
command -v pdftoppm >/dev/null || echo "[warn] pdftoppm missing - report PNGs will not be generated" >&2
# The workbook is packaged by python3 + openpyxl; the R library has no Excel
# writer and that env cannot build one (its Makeconf defines no compiler).
PY="${M12_PYTHON:-/usr/bin/python3}"
"$PY" -c "import openpyxl" 2>/dev/null || echo "[warn] openpyxl missing for $PY - the supplementary workbook step will fail" >&2

STEPS=(
  "01_inventory_nhanes_associations.R:inventory all five WAS schemes"
  "02_inventory_external_phenotypes.R:external phenotype dictionary + cohort capability"
  "03_build_opportunity_matrix.R:NHANES -> external opportunity matrix"
  "04_run_new_validations.R:fit the newly identified phenotypes"
  "05_master_validation_table.R:canonical association-level evidence base"
  "06_audit_analyses.R:concordance, correlations, enrichment, multi-cohort, strata"
  "08_audit_tables.R:tables 1-6"
  "10_phylum_level_validation.R:phylum-level directions vs SLL2"
  "11_alpha_diversity_validation.R:alpha-diversity directions vs SLL2"
  "13_supplementary_workbook.R:supplementary workbook (3 sheets, one .xlsx)"
)
TS="$(date -u +%Y%m%dT%H%M%SZ)"
for entry in "${STEPS[@]}"; do
  script="${entry%%:*}"; desc="${entry#*:}"
  echo; echo "=== ${script} — ${desc} ==="
  log="${LOGS}/${script%.R}_${TS}.log"
  "$RS" "${HERE}/${script}" > "$log" 2>&1 || true
  grep -E "^\[|FATAL|Error|ASSERTION" "$log" | head -25 || true
  if grep -qE "Execution halted|FATAL|ASSERTION FAILED" "$log"; then
    echo "[FATAL] ${script} failed - see ${log}" >&2; exit 1
  fi
done
echo; echo "=== done ==="
echo "tables:   ${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/audit"
echo "workbook: ${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/S10_summary_of_external_concordance_validation_data_and_results.xlsx"
