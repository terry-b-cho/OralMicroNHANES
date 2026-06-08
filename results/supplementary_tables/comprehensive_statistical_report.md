# Comprehensive Statistical Report: 20-Schema NHANES Oral Microbiome WAS Analysis

## Executive Summary

This report presents a comprehensive statistical analysis of **20 Wide Association Study (WAS) schemas** investigating associations between NHANES oral microbiome data and various health outcomes. The analysis encompasses **5 analysis types** across **4 transformation methods** with appropriate False Discovery Rate (FDR) correction.

## Methodology

### Analysis Framework
- **5 Analysis Types**:
  - `1_demoWAS`: Demographics → Microbiome (13 demographic predictors)
  - `2_oradWAS`: Microbiome → Oral health outcomes (4 outcomes)  
  - `3_exWAS`: Environmental exposures → Microbiome (473 variables)
  - `4_pheWAS`: Microbiome → Phenotypic measures (133 variables)
  - `5_outWAS`: Microbiome → Disease outcomes (16 variables)

- **4 Transformation Methods**:
  - `clr`: Centered log-ratio transformation
  - `lognorm`: Log-normal transformation  
  - `none`: Raw relative abundances
  - `hellinger`: Square-root transformation

### Statistical Correction
- **Method**: Benjamini-Hochberg False Discovery Rate (FDR) correction
- **Scope**: Scheme-wise correction (each of 20 schemas treated independently)
- **Target**: Main effect tests only (covariates excluded from correction)
- **Note**: Storey's q-value method was not used (use_qvalue=FALSE)

## Results by Schema

### 1. Demographics → Microbiome (1_demoWAS)

#### 1_demoWAS - CLR Transformation
- **Total Tests**: 17,537
- **Nominal Significance**: 
  - p < 0.05: 11,493 (65.54%)
  - p < 0.01: 9,775 (55.74%)
  - p < 0.001: 6,122 (34.91%)
- **FDR-Corrected Significance**:
  - p < 0.05: 11,055 (63.04%)
  - p < 0.01: 8,443 (48.14%)  
  - p < 0.001: 5,370 (30.62%)

#### 1_demoWAS - Log-Normal Transformation
- **Total Tests**: 17,537
- **Nominal Significance**: 
  - p < 0.05: 1,023 (5.83%)
  - p < 0.01: 810 (4.62%)
  - p < 0.001: 619 (3.53%)
- **FDR-Corrected Significance**:
  - p < 0.05: 669 (3.81%)
  - p < 0.01: 553 (3.15%)
  - p < 0.001: 404 (2.30%)

#### 1_demoWAS - None (Raw Abundances)
- **Total Tests**: 16,479
- **Nominal Significance**: 
  - p < 0.05: 1,170 (7.10%)
  - p < 0.01: 705 (4.28%)
  - p < 0.001: 466 (2.83%)
- **FDR-Corrected Significance**:
  - p < 0.05: 505 (3.06%)
  - p < 0.01: 363 (2.20%)
  - p < 0.001: 227 (1.38%)

#### 1_demoWAS - Hellinger Transformation
- **Total Tests**: 16,479
- **Nominal Significance**: 
  - p < 0.05: 1,721 (10.44%)
  - p < 0.01: 1,112 (6.75%)
  - p < 0.001: 764 (4.64%)
- **FDR-Corrected Significance**:
  - p < 0.05: 885 (5.37%)
  - p < 0.01: 662 (4.02%)
  - p < 0.001: 460 (2.79%)

### 2. Microbiome → Oral Health (2_oradWAS)

#### 2_oradWAS - CLR Transformation
- **Total Tests**: 5,396
- **Nominal Significance**: 
  - p < 0.05: 3,893 (72.15%)
  - p < 0.01: 2,910 (53.93%)
  - p < 0.001: 1,909 (35.38%)
- **FDR-Corrected Significance**:
  - p < 0.05: 3,388 (62.79%)
  - p < 0.01: 2,672 (49.52%)
  - p < 0.001: 1,753 (32.49%)

#### 2_oradWAS - Log-Normal Transformation
- **Total Tests**: 5,396
- **Nominal Significance**: 
  - p < 0.05: 1,450 (26.87%)
  - p < 0.01: 206 (3.82%)
  - p < 0.001: 160 (2.97%)
- **FDR-Corrected Significance**:
  - p < 0.05: 166 (3.08%)
  - p < 0.01: 132 (2.45%)
  - p < 0.001: 101 (1.87%)

#### 2_oradWAS - None (Raw Abundances)
- **Total Tests**: 4,642
- **Nominal Significance**: 
  - p < 0.05: 2,553 (55.00%)
  - p < 0.01: 2,261 (48.71%)
  - p < 0.001: 2,065 (44.49%)
- **FDR-Corrected Significance**:
  - p < 0.05: 2,394 (51.57%)
  - p < 0.01: 2,171 (46.77%)
  - p < 0.001: 2,017 (43.45%)

#### 2_oradWAS - Hellinger Transformation
- **Total Tests**: 4,642
- **Nominal Significance**: 
  - p < 0.05: 2,516 (54.20%)
  - p < 0.01: 2,263 (48.75%)
  - p < 0.001: 2,102 (45.28%)
- **FDR-Corrected Significance**:
  - p < 0.05: 2,384 (51.36%)
  - p < 0.01: 2,187 (47.11%)
  - p < 0.001: 2,052 (44.21%)

### 3. Environmental Exposures → Microbiome (3_exWAS)

#### 3_exWAS - CLR Transformation
- **Total Tests**: 638,077
- **Nominal Significance**: 
  - p < 0.05: 95,122 (14.91%)
  - p < 0.01: 42,970 (6.73%)
  - p < 0.001: 18,470 (2.89%)
- **FDR-Corrected Significance**:
  - p < 0.05: 22,108 (3.46%)
  - p < 0.01: 12,425 (1.95%)
  - p < 0.001: 7,709 (1.21%)

#### 3_exWAS - Log-Normal Transformation
- **Total Tests**: 637,604
- **Nominal Significance**: 
  - p < 0.05: 70,931 (11.12%)
  - p < 0.01: 32,604 (5.11%)
  - p < 0.001: 20,153 (3.16%)
- **FDR-Corrected Significance**:
  - p < 0.05: 22,403 (3.51%)
  - p < 0.01: 17,998 (2.82%)
  - p < 0.001: 11,615 (1.82%)

#### 3_exWAS - None (Raw Abundances)
- **Total Tests**: 468,080
- **Nominal Significance**: 
  - p < 0.05: 11,195 (2.39%)
  - p < 0.01: 3,940 (0.84%)
  - p < 0.001: 1,342 (0.29%)
- **FDR-Corrected Significance**:
  - p < 0.05: 512 (0.11%)
  - p < 0.01: 282 (0.06%)
  - p < 0.001: 162 (0.03%)

#### 3_exWAS - Hellinger Transformation
- **Total Tests**: 468,080
- **Nominal Significance**: 
  - p < 0.05: 17,526 (3.74%)
  - p < 0.01: 6,478 (1.38%)
  - p < 0.001: 2,357 (0.50%)
- **FDR-Corrected Significance**:
  - p < 0.05: 1,181 (0.25%)
  - p < 0.01: 698 (0.15%)
  - p < 0.001: 370 (0.08%)

### 4. Microbiome → Phenotypes (4_pheWAS)

#### 4_pheWAS - CLR Transformation
- **Total Tests**: 179,417
- **Nominal Significance**: 
  - p < 0.05: 15,935 (8.88%)
  - p < 0.01: 4,985 (2.78%)
  - p < 0.001: 1,016 (0.57%)
- **FDR-Corrected Significance**:
  - p < 0.05: 172 (0.10%)
  - p < 0.01: 12 (0.01%)
  - p < 0.001: 0 (0.00%)

#### 4_pheWAS - Log-Normal Transformation
- **Total Tests**: 179,417
- **Nominal Significance**: 
  - p < 0.05: 10,555 (5.88%)
  - p < 0.01: 3,992 (2.22%)
  - p < 0.001: 377 (0.21%)
- **FDR-Corrected Significance**:
  - p < 0.05: 32 (0.02%)
  - p < 0.01: 11 (0.01%)
  - p < 0.001: 0 (0.00%)

#### 4_pheWAS - None (Raw Abundances)
- **Total Tests**: 146,593
- **Nominal Significance**: 
  - p < 0.05: 68,898 (47.00%)
  - p < 0.01: 51,295 (34.99%)
  - p < 0.001: 36,019 (24.57%)
- **FDR-Corrected Significance**:
  - p < 0.05: 58,290 (39.76%)
  - p < 0.01: 41,405 (28.24%)
  - p < 0.001: 31,894 (21.76%)

#### 4_pheWAS - Hellinger Transformation
- **Total Tests**: 146,593
- **Nominal Significance**: 
  - p < 0.05: 57,022 (38.90%)
  - p < 0.01: 40,951 (27.94%)
  - p < 0.001: 28,189 (19.23%)
- **FDR-Corrected Significance**:
  - p < 0.05: 44,572 (30.41%)
  - p < 0.01: 31,207 (21.29%)
  - p < 0.001: 24,245 (16.54%)

### 5. Microbiome → Disease Outcomes (5_outWAS)

#### 5_outWAS - CLR Transformation
- **Total Tests**: 21,584
- **Nominal Significance**: 
  - p < 0.05: 5,394 (24.99%)
  - p < 0.01: 1,876 (8.69%)
  - p < 0.001: 642 (2.97%)
- **FDR-Corrected Significance**:
  - p < 0.05: 845 (3.91%)
  - p < 0.01: 285 (1.32%)
  - p < 0.001: 86 (0.40%)

#### 5_outWAS - Log-Normal Transformation
- **Total Tests**: 21,584
- **Nominal Significance**: 
  - p < 0.05: 1,838 (8.52%)
  - p < 0.01: 203 (0.94%)
  - p < 0.001: 51 (0.24%)
- **FDR-Corrected Significance**:
  - p < 0.05: 6 (0.03%)
  - p < 0.01: 5 (0.02%)
  - p < 0.001: 4 (0.02%)

#### 5_outWAS - None (Raw Abundances)
- **Total Tests**: 19,188
- **Nominal Significance**: 
  - p < 0.05: 12,356 (64.39%)
  - p < 0.01: 11,546 (60.17%)
  - p < 0.001: 10,785 (56.21%)
- **FDR-Corrected Significance**:
  - p < 0.05: 12,038 (62.74%)
  - p < 0.01: 11,379 (59.30%)
  - p < 0.001: 10,623 (55.36%)

#### 5_outWAS - Hellinger Transformation
- **Total Tests**: 19,188
- **Nominal Significance**: 
  - p < 0.05: 12,303 (64.12%)
  - p < 0.01: 11,459 (59.72%)
  - p < 0.001: 10,727 (55.90%)
- **FDR-Corrected Significance**:
  - p < 0.05: 11,963 (62.35%)
  - p < 0.01: 11,275 (58.76%)
  - p < 0.001: 10,595 (55.22%)

## Statistical Interpretation

### Q-Values
**Note**: Storey's q-value method was **not** used in this analysis (parameter `use_qvalue=FALSE`). Instead, the Benjamini-Hochberg procedure was applied for FDR correction. Q-values would provide additional information about the proportion of true null hypotheses (π₀), but were not computed in this analysis.

### Key Findings

1. **Transformation Effects**: 
   - CLR and raw abundances generally showed highest significance rates
   - Log-normal transformation was most conservative
   - Hellinger showed intermediate patterns

2. **Analysis Type Patterns**:
   - **Demographics → Microbiome (1_demoWAS)**: High significance with CLR (63% FDR < 0.05)
   - **Microbiome → Oral Health (2_oradWAS)**: Strong associations across all transformations
   - **Environmental → Microbiome (3_exWAS)**: Moderate significance after correction (~3.5% FDR < 0.05)
   - **Microbiome → Phenotypes (4_pheWAS)**: Variable by transformation (39.76% with raw, <0.1% with CLR)
   - **Microbiome → Disease Outcomes (5_outWAS)**: Strong with raw/Hellinger (>60% FDR < 0.05)

3. **FDR Control Effectiveness**:
   - Substantial reduction from nominal to corrected significance in most schemas
   - Maintained biological signal while controlling false discoveries
   - Schema-wise correction preserved meaningful associations

### Recommendations

1. **Priority Analysis**: Focus on high-significance schemas (2_oradWAS, 4_pheWAS with raw/Hellinger, 5_outWAS)
2. **Transformation Selection**: Consider biological interpretability alongside statistical significance
3. **Effect Size Analysis**: Examine magnitude of significant associations for practical relevance
4. **Replication**: Validate findings in independent datasets

## Conclusion

This comprehensive analysis of 20 WAS schemas reveals diverse patterns of microbiome associations across different health domains and statistical transformations. The scheme-wise FDR correction successfully controlled false discoveries while preserving biologically meaningful signals, particularly for oral health, phenotypic measures, and disease outcomes. 
