#!/bin/bash

# =============================================================================
# NHANES Survey Analysis Environment Setup Script
# setup_nhanes_environment.sh
# =============================================================================
# This script sets up the conda environment for NHANES survey-aware analysis
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

echo "Setting up NHANES survey analysis environment..."
echo "Environment path: ${ENV_PATH}"

# Create environment directory if it doesn't exist
mkdir -p "${ENV_DIR}"

# Load required modules in correct order
echo "Loading required modules..."
module purge
module load gcc/9.2.0
module load miniconda3/23.1.0

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
        source activate "${ENV_PATH}"
        echo "Environment activated successfully!"
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

# Try to create environment from yml file, fall back to simple method if it fails
echo "Attempting to create environment from yml file..."
if ! conda env create --prefix "${ENV_PATH}" -f "${ENV_YML}"; then
    echo "Environment file creation failed. Trying simple installation method..."
    
    # Create base environment
    echo "Creating base R environment..."
    conda create --prefix "${ENV_PATH}" -y -c conda-forge -c r r-base
    
    # Activate the environment
    source activate "${ENV_PATH}"
    
    # Install packages one by one
    echo "Installing required packages individually..."
    packages=(
        "r-survey"
        "r-broom" 
        "r-dplyr"
        "r-glue"
        "r-dbi"
        "r-rsqlite"
        "r-getopt"
        "r-logger"
        "r-readr"
        "r-tibble"
        "r-tidyr"
        "r-magrittr"
    )
    
    for package in "${packages[@]}"; do
        echo "Installing ${package}..."
        conda install -y -c conda-forge -c r "${package}" || {
            echo "Warning: Failed to install ${package} via conda, trying R install..."
        }
    done
    
    # Try to install any remaining packages via R
    echo "Installing any missing packages via R..."
    Rscript -e "
    required_packages <- c('survey', 'broom', 'dplyr', 'glue', 'DBI', 'RSQLite', 'getopt', 'logger', 'readr')
    
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        cat('Installing missing package:', pkg, '\n')
        install.packages(pkg, repos = 'https://cloud.r-project.org/')
      }
    }
    "
else
    echo "Environment created successfully from yml file!"
fi

# Activate the environment
echo "Activating environment..."
source activate "${ENV_PATH}"

# Verify installation
echo "Verifying package installation..."
R --version
echo ""

echo "Checking required R packages..."
Rscript -e "
required_packages <- c('survey', 'broom', 'dplyr', 'glue', 'DBI', 'RSQLite', 'getopt', 'logger', 'readr')
missing <- sapply(required_packages, function(pkg) {
  !requireNamespace(pkg, quietly = TRUE)
})

if (any(missing)) {
  cat('Missing packages:', paste(names(missing)[missing], collapse = ', '), '\n')
  quit(status = 1)
} else {
  cat('All required packages are installed!\n')
  cat('Package versions:\n')
  for (pkg in required_packages) {
    version <- packageVersion(pkg)
    cat(sprintf('  %s: %s\n', pkg, version))
  }
}
"

echo ""
echo "Environment setup complete!"
echo ""
echo "To use this environment in future sessions:"
echo "  module load gcc/9.2.0"
echo "  module load miniconda3/23.1.0"
echo "  source activate ${ENV_PATH}"
echo ""
echo "Environment is ready for NHANES survey analysis!" 