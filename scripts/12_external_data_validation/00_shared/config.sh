#!/bin/bash
# =============================================================================
# Module 12 — External Data Validation
# Shared configuration: paths, dataset roster, network policy
#
# Purpose : Single source of truth for every Module 12 shell script.
# Usage   : source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
# Inputs  : none (may be overridden by env var ORALMICRO_ROOT)
# Outputs : exported shell variables; creates no files
# Env     : bash; curl/wget on PATH
# =============================================================================

# === USER CONFIG =============================================================
# PROJECT_ROOT resolves from the environment first so this module is portable,
# unlike modules 4-11 which hardcode the old clone path (nhanes_oral_mirco_cho).
PROJECT_ROOT="${ORALMICRO_ROOT:-$(pwd)}"
# =============================================================================

set -euo pipefail

export PROJECT_ROOT
export M12_SCRIPTS="${PROJECT_ROOT}/scripts/12_external_data_validation"
export M12_OUT="${PROJECT_ROOT}/results/analyses_results/12_external_data_validation"
export M12_DATASETS="${M12_OUT}/datasets"
export M12_SHARED="${M12_SCRIPTS}/00_shared"
export M12_CONFIG="${M12_SCRIPTS}/config"

# --- Dataset roster ----------------------------------------------------------
# Processing order is significant: one dataset completes Phases A-I before the
# next begins. Datasets 1-2 carry the age/smoking payload.
#
# This roster is the module's single source of truth for which datasets are
# PRESENT. Datasets assessed and then removed entirely are absent from it; the
# reasons live in module12_README.md and are recorded nowhere else.
#
# TERMINOLOGY - "retained" means CONTRIBUTED EVIDENCE TO THE VALIDATION, and is
# never used for a dataset that was merely acquired. Per-dataset disposition is
# declared in config/<id>.yaml (`validation_disposition`); these arrays mirror
# it so shell callers do not have to parse YAML. 09_reports/build_catalog.R
# asserts the two agree.
export M12_DATASET_IDS=(
  HuangS_2020_AgePrediction
  WillisJR_2022_SacaLaLengua2
  GhensiP_2020_ImplantDisease
  WellsPM_2022_TwinsUKSaliva
)

# Contributed evidence to the external validation.
export M12_RETAINED_IDS=(
  HuangS_2020_AgePrediction
  WillisJR_2022_SacaLaLengua2
  GhensiP_2020_ImplantDisease
  WellsPM_2022_TwinsUKSaliva
)
# Datasets that were assessed and excluded are not present in the module and are
# not named in any config or script. The evidence that disqualified each is in
# module12_README.md and notebooks/02_qc_and_exclusions.html.

# --- Network policy ----------------------------------------------------------
# nature.com and static-content.springer.com return 403 / redirect to an IdP for
# programmatic clients. PMC and the Europe PMC REST API are the working routes
# for full text and supplementary files. See 00_shared/pmc_fetch.sh.
export M12_USER_AGENT="OralMicroNHANES-Module12/1.0 (research use; contact via repository)"
export M12_CURL_RETRY=5
export M12_CURL_RETRY_DELAY=3
export M12_CURL_TIMEOUT=600
export M12_CURL_CONNECT_TIMEOUT=30

export EPMC_REST="https://www.ebi.ac.uk/europepmc/webservices/rest"
export ENA_PORTAL="https://www.ebi.ac.uk/ena/portal/api"
export ENA_BROWSER="https://www.ebi.ac.uk/ena/browser/api"
export NCBI_EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# --- Acquisition tiers -------------------------------------------------------
# Tier 1: published processed tables (supplement / GitHub / Zenodo / cMD)
# Tier 2: public repository metadata (ENA/SRA/BioSample attributes)
# Tier 3: raw reads + reprocessing. NEVER run without an explicit go/no-go,
#         and never before Tiers 1-2 have been shown insufficient.
export M12_ALLOW_TIER3="${M12_ALLOW_TIER3:-0}"

# --- Per-dataset directory helper -------------------------------------------
m12_dataset_dir() {
  local ds="${1:?dataset_id required}"
  echo "${M12_DATASETS}/${ds}"
}

m12_require_dataset() {
  local ds="${1:?dataset_id required}"
  local d; d="$(m12_dataset_dir "$ds")"
  [[ -d "$d" ]] || { echo "[FATAL] Unknown dataset_id '${ds}' (no dir at ${d})" >&2; return 1; }
}
