#!/usr/bin/env Rscript

# NHANES Oral Microbiome GOLD Database Integration Runner
# This script runs the RMarkdown file and saves logs

# Load required libraries
library(rmarkdown)

# Setup paths
base_path <- "/Users/byeongyeoncho/main/github/nhanes_oral_mirco_cho"
script_dir <- file.path(base_path, "scripts/3_gold_db_microbial_phenotype")
logs_dir <- file.path(script_dir, "logs")

# Create logs directory if it doesn't exist
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# Generate timestamp for log files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(logs_dir, paste0("gold_db_integration_", timestamp, ".log"))

# Function to log messages
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste0("[", timestamp, "] ", message)
  cat(log_entry, "\n")
  cat(log_entry, "\n", file = log_file, append = TRUE)
}

# Start logging
log_message("Starting GOLD Database Integration Pipeline")
log_message(paste("Working directory:", getwd()))
log_message(paste("Script directory:", script_dir))
log_message(paste("Log file:", log_file))

# Change to script directory
setwd(script_dir)
log_message(paste("Changed working directory to:", getwd()))

# Check if RMarkdown file exists
rmd_file <- "gold_db_process_n_genus_mapping.Rmd"
if (!file.exists(rmd_file)) {
  log_message(paste("ERROR: RMarkdown file not found:", rmd_file))
  stop("RMarkdown file not found")
}

log_message(paste("Found RMarkdown file:", rmd_file))

# Run the RMarkdown file
tryCatch({
  log_message("Starting RMarkdown rendering...")
  
  # Render the RMarkdown file
  output_file <- render(
    rmd_file,
    output_format = "html_document",
    output_dir = file.path(base_path, "results/analyses_results/03_gold_db_microbial_phenotype_out"),
    quiet = FALSE
  )
  
  log_message(paste("RMarkdown rendering completed successfully"))
  log_message(paste("Output file:", output_file))
  
  # Check if output files were created
  intermediate_dir <- file.path(base_path, "results/analyses_results/03_gold_db_microbial_phenotype_out/intermediate")
  expected_files <- c(
    "gold_db_genus.csv",
    "ubiome_genus_mapping_complete.csv", 
    "mapping_summary_stats.csv"
  )
  
  log_message("Checking for expected output files...")
  for (file in expected_files) {
    file_path <- file.path(intermediate_dir, file)
    if (file.exists(file_path)) {
      file_size <- file.size(file_path)
      log_message(paste("✅", file, "created successfully (", file_size, "bytes)"))
    } else {
      log_message(paste("❌", file, "not found"))
    }
  }
  
  log_message("GOLD Database Integration Pipeline completed successfully!")
  
}, error = function(e) {
  log_message(paste("ERROR during RMarkdown rendering:", e$message))
  log_message("Pipeline failed. Check the log file for details.")
  stop(e)
})

# Final status
log_message("Pipeline execution finished")
log_message(paste("Log saved to:", log_file))

cat("\nGOLD Database Integration Pipeline Complete!\n")
cat("Check the log file for detailed execution information:\n")
cat("   ", log_file, "\n") 