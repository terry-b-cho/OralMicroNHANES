#!/bin/bash
# Master script to run all three modes of integrated diversity analysis with HALF PLOTS
# Using HPC resources: 16 cores, 128GB memory

set -e  # Exit on any error

echo "=== RUNNING ALL THREE MODES OF INTEGRATED DIVERSITY ANALYSIS (HALF PLOTS) ==="
echo "HPC Resources: 16 cores, 128GB memory"
echo "Full dataset: 9,349 samples"
echo "Half violin plots with box+jitter on left side"
echo ""

# Load modules and activate environment
module load conda/miniforge3/24.11.3-0
conda activate o2-nhanes_oral_env

# Set working directory
cd /n/groups/patel/terry/nhanes_oral_mirco_cho/scripts/6_alpha_beta_analyses

# Function to run a specific mode
run_mode() {
    local mode_name=$1
    local no_tick=$2
    local no_tick_asterisk=$3
    
    echo "=========================================="
    echo "RUNNING MODE: $mode_name"
    echo "NO_TICK_AND_BRACKET: $no_tick"
    echo "NO_TICK_AND_BRACKET_ASTERISK: $no_tick_asterisk"
    echo "=========================================="
    
    # Create a temporary script with the correct settings
    cp 2_integrated_diversity_FINAL_HPC_fixed.R temp_mode_script.R
    
    # Update the settings in the temporary script
    sed -i "s/NO_TICK_AND_BRACKET <- .*/NO_TICK_AND_BRACKET <- $no_tick/" temp_mode_script.R
    sed -i "s/NO_TICK_AND_BRACKET_ASTERISK <- .*/NO_TICK_AND_BRACKET_ASTERISK <- $no_tick_asterisk/" temp_mode_script.R
    sed -i "s/HALF_PLOTS <- .*/HALF_PLOTS <- TRUE/" temp_mode_script.R
    
    # Run the analysis with logging
    echo "Starting analysis for $mode_name..."
    Rscript temp_mode_script.R 2>&1 | tee "analysis_log_${mode_name}.txt"
    
    # Check if analysis completed successfully
    if [ $? -eq 0 ]; then
        echo "✓ $mode_name completed successfully"
        
        # Count generated PDF files
        pdf_count=$(find /n/groups/patel/terry/nhanes_oral_mirco_cho/results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/ -name "*.pdf" | wc -l)
        echo "Generated $pdf_count PDF files"
    else
        echo "✗ $mode_name failed"
        exit 1
    fi
    
    # Clean up temporary script
    rm -f temp_mode_script.R
    
    echo ""
}

# Run all three modes
echo "Starting analysis with full dataset (9,349 samples)..."
echo ""

# Mode 1: Normal full mode
run_mode "normal_full_mode" "FALSE" "FALSE"

# Mode 2: no_tick_and_bracket_asterisk mode  
run_mode "no_tick_and_bracket_asterisk" "FALSE" "TRUE"

# Mode 3: no_tick_and_bracket mode
run_mode "no_tick_and_bracket" "TRUE" "FALSE"

echo "=========================================="
echo "ALL THREE MODES COMPLETED SUCCESSFULLY!"
echo "=========================================="

# Final summary
total_pdfs=$(find /n/groups/patel/terry/nhanes_oral_mirco_cho/results/analyses_results/6_alpha_beta_analyses_out/integrated_diversity_plots/ -name "*.pdf" | wc -l)
echo "Total PDF files generated: $total_pdfs"

# List all log files
echo ""
echo "Log files created:"
ls -la analysis_log_*.txt

echo ""
echo "Analysis complete!"
