#!/bin/bash
# =============================================================================
# Module 12 Phase 2 — ordered driver
#
# Runs the external-validation analysis end to end. Each step depends on the
# previous one, so the order is fixed and the script stops on the first error.
#
# Usage : bash scripts/12_external_data_validation/10_validation_analysis/run_all.sh
# Env   : R 4.5.1 (conda nhanes-analysis); pdftoppm for the report PNGs
# =============================================================================
set -euo pipefail

PROJECT_ROOT="${ORALMICRO_ROOT:-$(pwd)}"
RS="${M12_RSCRIPT:-Rscript}"
HERE="${PROJECT_ROOT}/scripts/12_external_data_validation/10_validation_analysis"
LOGS="${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/validation/logs"
mkdir -p "$LOGS"
cd "$PROJECT_ROOT"

command -v "$RS" >/dev/null || { echo "[FATAL] Rscript not found: $RS" >&2; exit 1; }
command -v pdftoppm >/dev/null || echo "[warn] pdftoppm missing - report PNGs will not be generated" >&2

# The NHANES genus reference is built in Phase 1 and is the denominator for
# every prevalence filter here. Build it if it is absent.
REF="${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/nhanes_genus_reference.tsv"
if [[ ! -f "$REF" ]]; then
  echo "=== prerequisite: building the NHANES genus reference ==="
  "$RS" "${PROJECT_ROOT}/scripts/12_external_data_validation/09_reports/build_nhanes_genus_reference.R"
fi

STEPS=(
  "01_extract_nhanes_effects.R:NHANES age/smoking effect tables"
  "02_reanalyse_external.R:refit-cohort external model fitting"
  "03_assemble_comparisons.R:comparison assembly + concordance"
  "04_figures.R:the twelve panels of Supplementary Figure S15"
  "05_tables.R:tables"
)

TS="$(date -u +%Y%m%dT%H%M%SZ)"
for entry in "${STEPS[@]}"; do
  script="${entry%%:*}"; desc="${entry#*:}"
  echo
  echo "=== ${script} — ${desc} ==="
  log="${LOGS}/${script%.R}_${TS}.log"
  if ! "$RS" "${HERE}/${script}" 2>&1 | tee "$log" | grep -E "^\[|FATAL|Error|ASSERTION"; then
    :   # grep finding nothing is not a failure
  fi
  if grep -qE "Execution halted|FATAL|ASSERTION FAILED" "$log"; then
    echo "[FATAL] ${script} failed - see ${log}" >&2
    exit 1
  fi
done

echo
echo "=== done ==="
echo "figures: ${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/validation/figures/pdf  (Supplementary Figure S15 panels)"
echo "tables:  ${PROJECT_ROOT}/results/analyses_results/12_external_data_validation/validation"
