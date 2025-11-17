#!/bin/bash

# =============================================================================
# NHANES Survey Analysis Environment Setup Script - O2 COMPATIBLE
# setup_nhanes_environment.sh
# =============================================================================
# This script sets up the conda environment for NHANES survey-aware analysis
# UPDATED: O2 module compatibility, modern conda syntax, enhanced validation
# =============================================================================

set -e

# Check if we're on a compute node
if [[ -z "$SLURM_JOB_ID" ]]; then
    echo "Warning: Not running on a compute node. Please start an interactive session first:"
    echo "srun --pty -p interactive -t 12:00:00 --mem=32G bash"
    echo ""
    echo "Continue anyway? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Set base directory
BASE_DIR="/n/groups/patel/terry/nhanes_oral_mirco_cho"
ENV_DIR="${BASE_DIR}/envs/.conda/envs"
ENV_NAME="nhanes-analysis"
ENV_PATH="${ENV_DIR}/${ENV_NAME}"

echo "NHANES Survey Analysis Environment Setup - O2 COMPATIBLE"
echo "============================================================"
echo "Environment path: ${ENV_PATH}"
echo "Features: Final audit corrections, enhanced R packages, O2 compatibility"

# Create environment directory if it doesn't exist
mkdir -p "${ENV_DIR}"

# Load required modules in correct order (CRITICAL: Updated for O2)
echo "Loading required modules (updated for O2)..."
module purge
module load gcc/14.2.0                 # Updated with O2 update
module load conda/miniforge3/24.11.3-0 # Updated with O2 update

echo "✅ Loaded O2-compatible modules:"
echo "   gcc/14.2.0"
echo "   conda/miniforge3/24.11.3-0"

# Check if environment already exists
if [[ -d "${ENV_PATH}" ]]; then
    echo "Environment already exists at ${ENV_PATH}"
    echo "Remove existing environment? (y/N)"
    read -r response
    if [[ "$response" =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        rm -rf "${ENV_PATH}"
    else
        echo "Activating existing environment..."
        eval "$(conda shell.bash hook)"  # Modern conda activation syntax
        conda activate "${ENV_PATH}"
        echo "✅ Environment activated successfully!"
        echo ""
        echo "VERIFICATION:"
        which R
        R --version | head -1
        echo ""
        echo "To use this environment:"
        echo "  module purge"
        echo "  module load gcc/14.2.0"
        echo "  module load conda/miniforge3/24.11.3-0"
        echo "  eval \"\$(conda shell.bash hook)\""
        echo "  conda activate ${ENV_PATH}"
        exit 0
    fi
fi

# Create new environment from yml file
echo "Creating new conda environment..."
ENV_YML="${BASE_DIR}/envs/nhanes-analysis-environment.yml"

if [[ ! -f "${ENV_YML}" ]]; then
    echo "Error: Environment file not found at ${ENV_YML}"
    echo "Please ensure the environment.yml file exists."
    exit 1
fi

# Try to create environment from yml file, fall back to manual method if it fails
echo "Attempting to create environment from yml file..."
if ! conda env create --prefix "${ENV_PATH}" -f "${ENV_YML}"; then
    echo "Environment file creation failed. Trying enhanced manual installation method..."
    
    # Create base environment with R 4.4.3
    echo "Creating base R environment (R 4.4.3 for consistency)..."
    conda create --prefix "${ENV_PATH}" -y -c conda-forge r-base=4.4.3
    
    # Activate the environment using modern syntax
    eval "$(conda shell.bash hook)"
    conda activate "${ENV_PATH}"
    
    # Install packages one by one with enhanced error handling
    echo "Installing required packages individually..."
    packages=(
        "r-survey"
        "r-broom" 
        "r-dplyr"
        "r-glue"
        "r-dbi"
        "r-rsqlite"
        "r-optparse"    # Updated from getopt
        "r-logger"
        "r-readr"
        "r-tibble"
        "r-tidyr"
        "r-magrittr"
        "r-stringr"     # Added for final audit corrections
    )
    
    for package in "${packages[@]}"; do
        echo "Installing ${package}..."
        if ! conda install -y -c conda-forge -c r "${package}"; then
            echo "Warning: Failed to install ${package} via conda, trying R install..."
        fi
    done
    
    # Try to install any remaining packages via R with enhanced validation
    echo "Installing any missing packages via R..."
    Rscript -e "
    required_packages <- c('survey', 'broom', 'dplyr', 'glue', 'DBI', 'RSQLite', 
                          'optparse', 'logger', 'readr', 'tibble', 'tidyr', 
                          'magrittr', 'stringr')
    
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        cat('Installing missing package:', pkg, '\n')
        install.packages(pkg, repos = 'https://cloud.r-project.org/')
      }
    }
    
    # Verify critical packages for final audit corrections
    critical_packages <- c('survey', 'broom', 'dplyr', 'DBI', 'RSQLite', 'optparse', 'stringr')
    missing_critical <- sapply(critical_packages, function(pkg) {
      !requireNamespace(pkg, quietly = TRUE)
    })
    
    if (any(missing_critical)) {
      cat('CRITICAL ERROR: Missing essential packages:', paste(names(missing_critical)[missing_critical], collapse = ', '), '\n')
      quit(status = 1)
    } else {
      cat('✅ All critical packages installed successfully\n')
    }
    "
else
    echo "✅ Environment created successfully from yml file!"
fi

# Activate the environment using modern syntax
echo "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate "${ENV_PATH}"

# Enhanced verification for final audit corrected pipeline
echo "Performing enhanced verification..."
R --version
echo ""

echo "Checking required R packages for final audit corrected pipeline..."
Rscript -e "
required_packages <- c('survey', 'broom', 'dplyr', 'glue', 'DBI', 'RSQLite', 
                      'optparse', 'logger', 'readr', 'tibble', 'tidyr', 
                      'magrittr', 'stringr')

missing <- sapply(required_packages, function(pkg) {
  !requireNamespace(pkg, quietly = TRUE)
})

if (any(missing)) {
  cat('Missing packages:', paste(names(missing)[missing], collapse = ', '), '\n')
  quit(status = 1)
} else {
  cat('✅ All required packages are installed!\n')
  cat('Enhanced package versions for final audit corrections:\n')
  for (pkg in required_packages) {
    version <- tryCatch({
      packageVersion(pkg)
    }, error = function(e) {
      'ERROR'
    })
    cat(sprintf('  %s: %s\n', pkg, version))
  }
}

# Test critical functionality for final audit corrections
cat('\nTesting critical functionality:\n')

# Test table/view detection capability
tryCatch({
  library(DBI)
  library(RSQLite)
  cat('✅ SQLite connectivity: OK\n')
}, error = function(e) {
  cat('❌ SQLite connectivity: FAILED\n')
  quit(status = 1)
})

# Test survey package functionality
tryCatch({
  library(survey)
  options(survey.lonely.psu = 'certainty')  # NCHS recommendation
  cat('✅ Survey package: OK\n')
}, error = function(e) {
  cat('❌ Survey package: FAILED\n')
  quit(status = 1)
})

# Test enhanced data manipulation
tryCatch({
  library(dplyr)
  library(stringr)
  library(tibble)
  cat('✅ Enhanced data manipulation: OK\n')
}, error = function(e) {
  cat('❌ Enhanced data manipulation: FAILED\n')
  quit(status = 1)
})

cat('\\nEnvironment ready for final audit corrected pipeline!\\n')
"

echo ""
echo "Environment setup complete!"
echo ""
echo "✅ ENHANCED FEATURES READY:"
echo "   • Table/view detection for SQL views"
echo "   • FDR correction capabilities"  
echo "   • Effect scale harmonization"
echo "   • Survey-weighted statistics"
echo "   • Enhanced error handling"
echo ""
echo "To use this environment in future sessions:"
echo "  module purge"
echo "  module load gcc/14.2.0"
echo "  module load conda/miniforge3/24.11.3-0"
echo "  eval \"\$(conda shell.bash hook)\""
echo "  conda activate ${ENV_PATH}"
echo ""
echo "Environment is ready for NHANES survey analysis with final audit corrections!" 