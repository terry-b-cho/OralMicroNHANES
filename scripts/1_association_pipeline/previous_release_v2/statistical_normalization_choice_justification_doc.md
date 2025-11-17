## **Statistical Justifications for Survey-Weighted Regression**

### **2_oradWAS, 5_outWAS (Microbiome → Binary Outcomes)**
```r
# Survey-weighted logistic regression with CLR normalization:
E[Y|X] = Σᵢ wᵢ × [1 / (1 + exp(-(β₀ + β₁ × CLR(microbiome_i + pseudocount) + Σⱼ βⱼ × covariate_ij)))]

where:
- wᵢ = WTMEC2YR (survey weights)
- Strata: SDMVSTRA  
- PSU: SDMVPSU
- CLR(microbiome) handles compositional constraints
```
- **Binary health outcomes** → survey-weighted logistic regression required
- **CLR handles compositional nature** of microbiome predictors  
- **Survey weights account for NHANES sampling design**
- **Robust to outliers** in abundance distributions

### **1_demoWAS, 3_exWAS (Exposures → Microbiome Abundances)**
```r
# Survey-weighted linear regression with log normalization:
E[log(abundance + pseudocount)|X] = Σᵢ wᵢ × (β₀ + β₁ × exposure_i + Σⱼ βⱼ × covariate_ij)

where:
- wᵢ = WTMEC2YR (survey weights)
- Strata: SDMVSTRA
- PSU: SDMVPSU  
- log(abundance) normalizes right-skewed distributions
```
- **Continuous abundance outcomes** → survey-weighted linear regression appropriate
- **Log transformation normalizes** right/left-skewed abundance distributions
- **Stabilizes variance** across abundance ranges
- **Survey design preserves population representativeness**

### **4_pheWAS, 6_zimWAS (Microbiome → Continuous Outcomes)**
```r
# Survey-weighted linear regression with log normalization:
E[phenotype|X] = Σᵢ wᵢ × (β₀ + β₁ × log(microbiome_i + pseudocount) + Σⱼ βⱼ × covariate_ij)

where:
- wᵢ = WTMEC2YR (survey weights)
- Strata: SDMVSTRA
- PSU: SDMVPSU
- log(microbiome) normalizes right/left-skewed microbiome predictors
```
- **Continuous health outcomes** → survey-weighted linear regression appropriate
- **Log normalization of microbiome predictors** handles right/left-skewed distributions
- **Stabilizes variance** in microbiome abundance ranges
- **Maintains interpretability** as fold-change effects
- **Survey weights ensure population-level inference**

## **Key Survey Design Elements:**

### **All Pipelines Include:**
- **Survey weights**: `WTMEC2YR` (2-year mobile examination center weights)
- **Stratification**: `SDMVSTRA` (masked variance stratum)  
- **Primary sampling units**: `SDMVPSU` (masked variance PSU)
- **Nested design**: Accounts for complex NHANES sampling
- **Population inference**: Results generalizable to US population

### **Normalization Rationale:**
- **CLR for compositional predictors** (microbiome in logistic models)
- **Log for abundance outcomes/predictors** (microbiome in linear models)
- **Survey weighting preserves** population representativeness regardless of normalization

This approach ensures both **statistical validity** and **population-level generalizability** of microbiome associations!