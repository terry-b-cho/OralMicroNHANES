# NHANES Oral Microbiome M-WAS Pipeline - Troubleshooting & Appendix

<div align="center">
  <img src="assets/logo/oral_micro_logo.svg" alt="NHANES Oral Microbiome Analysis Logo" width="400"/>
</div>

---

## **Quick Start (TL;DR)**

```bash
# 1. Setup environment (one-time)
srun --pty -p interactive -t 12:00:00 --mem=32G bash
bash scripts/1_association_pipeline/setup_nhanes_environment.sh

# 2. Run preprocessing (if not done)
sbatch scripts/0_transform_n_preprocess_ssfiles/run_transformation.sh
# Wait for completion, then run missing data filling and schema creation

# 3. Run ALL analyses (6 types × 3 normalizations = 18 total)
./scripts/1_association_pipeline/run_all_was_analyses.sh

# 4. Check results
ls results/*/result_*/
```

**Results will be in `results/<analysis_type>_out/result_<normalization>/` with FDR-corrected p-values ready for downstream analyses!**

---

## 🔧 **Comprehensive Troubleshooting Guide**

### **Environment Setup Issues**

#### **Problem: Conda Environment Creation Fails**
```bash
# Solution 1: Try backup method
bash scripts/1_association_pipeline/setup_nhanes_environment_simple.sh

# Solution 2: Manual package installation
R -e "install.packages(c('survey', 'broom', 'dplyr'), repos='https://cloud.r-project.org/')"

# Solution 3: Check module loading order
module purge
module load gcc/14.2.0          # Updated with O2 update
module load conda/miniforge3/24.11.3-0   # Updated with O2 update
```

#### **Problem: R Package Installation Errors**
```bash
# Check R version compatibility
R --version

# Install packages individually
R -e "install.packages('survey', repos='https://cloud.r-project.org/')"
R -e "install.packages('broom', repos='https://cloud.r-project.org/')"
R -e "install.packages('dplyr', repos='https://cloud.r-project.org/')"

# For persistent issues, try different CRAN mirror
R -e "install.packages('survey', repos='https://cran.rstudio.com/')"
```

### **Job Submission and Execution Issues**

#### **Problem: SLURM Jobs Fail Due to Memory**
```bash
# Increase memory allocation
./scripts/1_association_pipeline/run_complete_was_analysis.sh 4_pheWAS clr "" 120G 24:00:00

# Check memory usage in logs
cat results/4_pheWAS_out/logs/slurm_*.err | grep -i "memory\|oom"

# For very large analyses, use high-memory partition
# Edit SLURM scripts to use: #SBATCH -p highmem
```

#### **Problem: Jobs Timeout**
```bash
# Increase time limit
./scripts/1_association_pipeline/run_complete_was_analysis.sh 2_oradWAS clr "" 64G 48:00:00

# Check job progress
squeue -u $USER --format="%.18i %.9P %.50j %.8u %.8T %.10M %.9l %.6D %R"

# For test runs, use smaller datasets
./scripts/1_association_pipeline/run_complete_was_analysis.sh 1_demoWAS clr test
```

#### **Problem: Permission Denied Errors**
```bash
# Make scripts executable
chmod +x scripts/1_association_pipeline/*.sh
chmod +x scripts/0_transform_n_preprocess_ssfiles/*.sh

# Check file ownership
ls -la scripts/1_association_pipeline/

# Fix permissions if needed
chmod 755 scripts/1_association_pipeline/run_all_was_analyses.sh
```

### **Data and Database Issues**

#### **Problem: Input Database Not Found**
```bash
# Check database exists
ls -la data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite

# Verify database integrity
sqlite3 data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite ".tables"

# Check file permissions
file data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite
```

#### **Problem: Schema Structure Files Missing**
```bash
# Check if preprocessing completed
ls -la results/0_ss_files/

# Re-run schema creation if needed
bash scripts/0_transform_n_preprocess_ssfiles/run_ss_file_create.sh

# Verify schema file content
head -5 results/0_ss_files/1_demoWAS_clr_schema_structure.csv
```

### **Analysis-Specific Issues**

#### **Problem: No Results Generated**
```bash
# Check if aggregation completed
ls -la results/*/result_*/

# Manual aggregation if needed
Rscript scripts/1_association_pipeline/aggregate_was_results.R 4_pheWAS clr results

# Check for empty result files
find results/ -name "*.rds" -size 0
```

#### **Problem: Statistical Errors in R**
```bash
# Check R session info in logs
grep -A 20 "sessionInfo" results/*/logs/slurm_*.out

# Verify survey package version
R -e "packageVersion('survey')"

# Test survey implementation
Rscript scripts/1_association_pipeline/test_survey_implementation.R
```

### **Performance Optimization**

#### **Problem: Jobs Running Too Slowly**
```bash
# Use more CPU cores
# Edit SLURM scripts: #SBATCH -c 16

# Use faster storage if available
# Check if /tmp has more space
df -h /tmp

# Optimize R memory usage
# Add to R scripts: options(mc.cores = parallel::detectCores())
```

#### **Problem: Too Many Jobs in Queue**
```bash
# Check queue status
squeue -u $USER | wc -l

# Cancel all jobs if needed
scancel -u $USER

# Submit jobs with dependencies
# Edit scripts to use: #SBATCH --dependency=afterok:$JOBID
```

---

## 📝 **Publication Methods Templates**

### **Survey Design Methods:**
> "All analyses accounted for the complex survey design of NHANES using appropriate stratification (`SDMVSTRA`), clustering (`SDMVPSU`), and survey weights (`WTINT2YR`/`WTMEC2YR`) scaled across multiple 2-year cycles. Models were fitted using survey-weighted generalized linear models (`svyglm`) from the R `survey` package to ensure proper variance estimation and population-level inference."

### **Multiple Comparisons Methods:**
> "P-values were adjusted for multiple comparisons using the Benjamini-Hochberg false discovery rate (FDR) procedure to control the expected proportion of false discoveries. Associations with FDR-adjusted p-values < 0.05 were considered statistically significant."

### **Microbiome Transformation Methods:**
> "Microbiome OTU relative abundance data were transformed using three methods: (1) centered log-ratio (CLR) transformation to address compositional nature, (2) log-normal transformation with pseudo-count addition, and (3) untransformed relative abundances. CLR transformation was used as the primary method for reporting results."

### **Statistical Analysis Methods:**
> "Linear regression models were used for continuous outcomes and logistic regression for binary outcomes. All models included age, sex, and race/ethnicity as covariates. Survey-weighted variance estimation was performed using Taylor series linearization to account for the complex sampling design of NHANES."

### **Data Processing Methods:**
> "Microbiome data preprocessing included prevalence filtering (retaining taxa present in >0.1% of samples) and three normalization approaches. Quality control procedures removed samples with insufficient sequencing depth and taxa with excessive missing data."

---

## **Result Interpretation Guide**

### **Understanding Output Files**
```r
# Load main results
results <- readRDS("results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds")

# Key columns for interpretation:
# - term: Model coefficient ("indep_var" = main effect of interest)
# - estimate: Effect size (log-odds for logistic, coefficient for linear)
# - std.error: Survey-adjusted standard error
# - p.value: Unadjusted p-value
# - p.value.fdr: FDR-corrected p-value (USE THIS FOR SIGNIFICANCE)
# - p.value.bonferroni: Bonferroni-corrected p-value (very conservative)
# - phenotype: Dependent variable (outcome)
# - exposure: Independent variable (predictor)
# - n_obs: Effective sample size
```

### **Statistical Significance Criteria**
- **Primary**: Use `p.value.fdr < 0.05` for reporting significant associations
- **Conservative**: Use `p.value.bonferroni < 0.05` for very stringent control
- **Effect Size**: Consider biological significance alongside statistical significance
- **Sample Size**: Ensure adequate power (n_obs > 100 recommended)

### **Common Result Patterns**
```r
# Filter significant results
significant_results <- results %>%
  filter(term == "indep_var", p.value.fdr < 0.05)

# Check for multiple significant associations
table(significant_results$phenotype)

# Examine effect sizes
summary(significant_results$estimate)
```

---

## 🏆 **Pipeline Features & Validation**

### **✅ Statistical Rigor Checklist:**
- [ ] Proper NHANES survey design implementation
- [ ] Multiple comparisons correction (FDR + Bonferroni)
- [ ] Survey-weighted variance estimation
- [ ] Complete case analysis with weight scaling
- [ ] Appropriate model selection (linear vs. logistic)

### **✅ Computational Efficiency Features:**
- [ ] Automated SLURM job submission and monitoring
- [ ] Parallel processing of individual dependent variables
- [ ] Resource optimization based on analysis size
- [ ] Comprehensive error handling and logging
- [ ] Automatic result aggregation with summary statistics

### **✅ Reproducibility Standards:**
- [ ] Version-controlled conda environments
- [ ] Comprehensive documentation with examples
- [ ] Standardized file naming and output structure
- [ ] Built-in validation and testing scripts
- [ ] Complete audit trail of analysis parameters

### **✅ Flexibility Options:**
- [ ] Six different analysis types for various research questions
- [ ] Three normalization methods for microbiome data
- [ ] Configurable variable sets via text files
- [ ] Test mode for pipeline validation
- [ ] Custom resource allocation for different analysis sizes

---

## **Advanced Diagnostics**

### **Checking Analysis Completeness**
```bash
# Count expected vs actual result files
expected_files=18  # 6 analyses × 3 normalizations
actual_files=$(find results/ -name "*_tidied_complete.rds" | wc -l)
echo "Expected: $expected_files, Found: $actual_files"

# Check for failed analyses
find results/ -name "slurm_*.err" -exec grep -l "ERROR\|FAILED" {} \;

# Verify aggregation summaries
find results/ -name "*_aggregation_summary.txt" -exec wc -l {} \;
```

### **Performance Monitoring**
```bash
# Check job resource usage
sacct -u $USER --format=JobID,JobName,MaxRSS,Elapsed,State

# Monitor disk usage
du -sh results/*/

# Check database sizes
ls -lh data/00_nhanes_omp_transformed_db/*.sqlite
```

### **Quality Control Checks**
```r
# Load and inspect results
library(dplyr)
results <- readRDS("results/1_demoWAS_out/result_clr/1_demoWAS_clr_tidied_complete.rds")

# Check for unusual patterns
summary(results$p.value)
hist(results$p.value, breaks=50)

# Verify FDR correction
sum(results$p.value.fdr < 0.05, na.rm=TRUE)
sum(results$p.value < 0.05, na.rm=TRUE)

# Check sample sizes
summary(results$n_obs)
table(results$series)  # Should show F and G cycles
```

---

## 📞 **Support & Resources**

### **Getting Help**

1. **📖 Check Documentation First**
   - Main README: [README.md](README.md)
   - Chapter 1: [0_transform_n_preprocess_ssfiles_documentation.md](scripts/0_transform_n_preprocess_ssfiles/0_transform_n_preprocess_ssfiles_documentation.md)
   - Chapter 2: [1_association_pipeline_documentation.md](scripts/1_association_pipeline/1_association_pipeline_documentation.md)

2. **Search Existing Issues**
   - Check GitHub Issues for similar problems
   - Look for closed issues with solutions
   - Review discussion threads

3. **🆕 Create New Issue**
   - Provide detailed error messages
   - Include system information (R version, OS, etc.)
   - Attach relevant log files
   - Describe steps to reproduce the problem

4. **👥 Contact Development Team**
   - For collaboration inquiries
   - For custom analysis requests
   - For integration with other pipelines

### **Useful External Resources**

- **NHANES Survey Design**: [CDC Analytic Guidelines](https://wwwn.cdc.gov/nchs/nhanes/analyticguidelines.aspx)
- **R Survey Package**: [Survey Package Documentation](https://cran.r-project.org/web/packages/survey/index.html)
- **Microbiome Analysis**: [Compositional Data Analysis Methods](https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full)
- **SLURM Documentation**: [SLURM Workload Manager](https://slurm.schedmd.com/documentation.html)
- **O2 Cluster Guide**: [HMS Research Computing](https://harvardmed.atlassian.net/wiki/spaces/O2/overview)

### **Community Resources**

- **Microbiome Analysis**: [QIIME2 Forum](https://forum.qiime2.org/)
- **R Statistical Computing**: [R-help Mailing List](https://stat.ethz.ch/mailman/listinfo/r-help)
- **Survey Statistics**: [Complex Surveys in R](https://r-survey.r-forge.r-project.org/)

---

## **Appendix: File Structure Reference**

### **Input Files Required**
```
data/00_nhanes_omp_abundance_db/nhanes_031725.sqlite    # Raw microbiome data
configs/1_demoWAS_vars.txt                              # Demographics variables
configs/2_oradWAS_vars.txt                              # Oral health variables
configs/3_exWAS_vars.txt                                # Exposure variables
configs/4_pheWAS_vars.txt                               # Phenotype variables
configs/5_outWAS_vars.txt                               # Disease outcome variables
configs/6_zimWAS_vars.txt                               # Lab measurement variables
```

### **Intermediate Files Generated**
```
data/00_nhanes_omp_transformed_db/nhanes_oral_transformed.sqlite         # Transformed data
data/00_nhanes_omp_transformed_db/nhanes_oral_transformed_complete.sqlite # Complete data
results/0_ss_files/*.csv                                                 # Schema files (18 total)
```

### **Final Output Files**
```
results/<analysis>_out/result_<norm>/<analysis>_<norm>_tidied_complete.rds    # Main results
results/<analysis>_out/result_<norm>/<analysis>_<norm>_glanced_complete.rds   # Model stats
results/<analysis>_out/result_<norm>/<analysis>_<norm>_rsq_complete.rds       # R-squared
results/<analysis>_out/result_<norm>/<analysis>_<norm>_aggregation_summary.txt # Summary
results/<analysis>_out/logs/slurm_*.out                                       # Job logs
results/<analysis>_out/logs/slurm_*.err                                       # Error logs
```

### **Log File Locations**
```
results/1_demoWAS_out/logs/     # Demographics analysis logs
results/2_oradWAS_out/logs/     # Oral health analysis logs
results/3_exWAS_out/logs/       # Exposure analysis logs
results/4_pheWAS_out/logs/      # Phenotype analysis logs
results/5_outWAS_out/logs/      # Disease outcome analysis logs
results/6_zimWAS_out/logs/      # Lab measurement analysis logs
```

---

<div align="center">
  <strong>🔧 Need more help? Check the main documentation or create an issue on GitHub!</strong>
</div> 