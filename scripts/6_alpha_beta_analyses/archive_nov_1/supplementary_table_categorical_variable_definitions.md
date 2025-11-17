# Supplementary Table: Categorical Variable Definitions

## Complete Definitions for All Categorical Variables Used in Diversity Analysis

This table provides comprehensive documentation of all categorical variables used in the integrated diversity analysis. Each variable includes:
- **Categorical Variable Name**: Variable name as used in analysis scripts
- **Definition in R Code**: Exact R code used to create the categorical variable
- **Definition Description**: Human-readable description of the variable and its categories
- **Reference Level**: The reference (baseline) category used in statistical comparisons
- **NA Count**: Number of participants with missing data for this variable
- **Count Per Category**: Sample size for each category level
- **Table Source**: Original NHANES table(s) and variable(s) used to derive this categorical variable

---

| Categorical Variable Name | Definition Description | Reference Level | NA Count | Count Per Category | Table Source |
|---------------------------|------------------------|-----------------|----------|---------------------|--------------|
| Gender | Biological sex of participant from NHANES demographics (RIAGENDR). Male=1, Female=2 | Female | 0 | Female = 4809; Male = 4853 | NHANES Demographics (DEMO_F, DEMO_G) |
| Age_group | Age categorized into 10-year groups: 18-29, 30-39, 40-49, 50-59, 60-69, 70-79, 80+ years from RIDAGEYR | 30-39 | 0 | 30-39 = 1632; 14-19 = 1658; 20-29 = 1711; 40-49 = 1648; 50-59 = 1496; 60-69 = 1517 | NHANES Demographics (DEMO_F, DEMO_G) - derived from RIDAGEYR |
| Education_level | Highest education level completed from DMDEDUC2. Categories: <9th Grade, 9-11th Grade, High School, College/AA (Some college or Associate's degree), College Graduate (Bachelor's degree or higher) | College/AA | 188 | College/AA = 2414; < 9th Grade = 1141; 9-11th Grade = 2108; High School = 1897; College Graduate = 1914 | NHANES Demographics (DEMO_F, DEMO_G) - DMDEDUC2 |
| Ethnicity | Race/ethnicity from RIDRETH3. Categories: White (Non-Hispanic White), Black (Non-Hispanic Black), Hispanic (Mexican American, Other Hispanic), Asian (Non-Hispanic Asian), Other (Other Race including Multi-racial) | White | 0 | White = 3472; Black = 2344; Mexican = 1684; Other = 1168; Other Hispanic = 994 | NHANES Demographics (DEMO_F, DEMO_G) - RIDRETH3 |
| US_born | Country of birth from DMDBORN4. US Born vs Non-US Born (1=Born in US, 2=Not born in US) | US Born | 5 | US Born = 6943; non-US Born = 2714 | NHANES Demographics (DEMO_F, DEMO_G) - DMDBORN4 |
| Household_size | Number of people in household from DMDHHSIZ. Categories: 1, 2, 3, 4, 5, 6, 7+ members | 4 | 0 | 4 = 2009; 1 = 882; 2 = 2224; 3 = 1799; 5 = 1342; 6 = 673; 7+ = 733 | NHANES Demographics (DEMO_F, DEMO_G) - DMDHHSIZ |
| Marital_status | Marital status from DMDMARTL. Categories: Married, Widowed, Divorced, Separated, Never married, Living with partner | Married | 1663 | Married = 3959; Divorced = 885; Living with Partner = 748; Never Married = 1830; Separated = 309; Widowed = 268 | NHANES Demographics (DEMO_F, DEMO_G) - DMDMARTL |
| Interview_language | Language of interview from SIALANG. English vs Spanish (1=English, 2=Spanish) | English | 0 | English = 8356; Spanish = 1306 | NHANES Demographics (DEMO_F, DEMO_G) - SIALANG |
| Income_to_poverty_ratio | Family income to poverty threshold ratio from INDFMPIR. Categorized into 7 levels based on percentage of federal poverty level: Below 50%, 50-99%, 100-124%, 125-149%, 150-184%, 185-199%, 200% and Over. Reference level chosen as median category (150-184%) | 150-184% | 805 | 150-184% = 612; Below 50% = 828; 50-99% = 1535; 100-124% = 830; 125-149% = 660; 185-199% = 222; 200% and Over = 4170 | NHANES Demographics (DEMO_F, DEMO_G) - INDFMPIR |
| Denture | Presence of dentures from OHAROCDE. Binary: Denture (OHAROCDE=1) vs control (all oral health indicators=0) | control | 5259 | control = 3917; Denture = 486 | NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCDE |
| Gum_disease | Self-reported gum disease from OHAROCGP. Binary: Gum disease (OHAROCGP=1) vs control (all oral health indicators=0) | control | 3142 | control = 3917; Gum disease = 2603 | NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCGP |
| Oral_hygiene | Poor oral hygiene status from OHAROCOH. Binary: Poor oral hygiene (OHAROCOH=1) vs control (all oral health indicators=0) | control | 2921 | control = 3917; Poor oral hygiene = 2824 | NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCOH |
| Tooth_decay | Presence of tooth decay from OHAROCDT. Binary: Tooth decay (OHAROCDT=1) vs control (all oral health indicators=0) | control | 3344 | control = 3917; Tooth decay = 2401 | NHANES Oral Health (OHQ_F, OHQ_G) - derived variable OHAROCDT |
| Smoking_status | Smoking status derived from SMQ020 (ever smoked 100 cigarettes) and SMQ040 (current smoking). Categories: Never smoker (SMQ020=2), Former smoker (SMQ020=1 & SMQ040 in 1,2), Current smoker (SMQ020=1 & SMQ040=3) | Never smoker | 1662 | Never smoker = 4523; Former smoker = 1890; Current smoker = 1587 | NHANES Smoking (SMQ_F, SMQ_G) - derived from SMQ020, SMQ040 |
| Hepatitis_C_antibody | Hepatitis C antibody test result from LBDHCV. Binary: Positive (LBDHCV=1) vs Negative (LBDHCV=2). Indeterminate results (LBDHCV=5) treated as NA | Negative | 564 | Negative = 8936; Positive = 162 | NHANES Laboratory (HEPC_F, HEPC_G) - LBDHCV |
| HPV_PCR_summary | Human Papillomavirus (HPV) PCR summary result from LBDRPCR. Binary: Positive (LBDRPCR=1) vs Negative (LBDRPCR=2). Inadequate samples (LBDRPCR=3) treated as NA | Negative | 9662 | Negative = 0; Positive = 0 | NHANES Laboratory (HPV_F, HPV_G) - LBDRPCR |
| BMI_category | Body Mass Index category from BMXBMI using CDC criteria. Categories: Underweight (<18.5), Healthy weight (18.5-24.9), Overweight (25.0-29.9), Class 1 Obesity (30.0-34.9), Class 2-3 Obesity (≥35.0 kg/m²) | Healthy weight | 75 | Healthy weight = 3119; Underweight = 283; Overweight = 2913; Class 1 Obesity = 1798; Class 2-3 Obesity = 1474 | NHANES Examination (BMX_F, BMX_G) - BMXBMI |
| Blood_pressure | Blood pressure category using AHA/ACC 2017 guidelines. MSYSTOLIC and MDIASTOLIC are means of three BP readings (BPXSY1-3, BPXDI1-3). Categories: Normal (<120/<80), Elevated (120-129/<80), Stage 1 (130-139/80-89), Stage 2 (≥140/≥90), Crisis (>180/>120 mmHg) | Normal | 320 | Normal = 5126; Elevated = 1460; Hypertension Stage 1 = 1638; Hypertension Stage 2 = 1073; Hypertensive Crisis = 45 | NHANES Examination (BPX_F, BPX_G) - derived from BPXSY1-3, BPXDI1-3 |
| Pulse_category | Pulse rate category from BPXPLS (60-second radial pulse). Categories based on quintile distribution: <60 bpm, 60-70 bpm, 70-75 bpm (reference, median), 75-85 bpm, 85+ bpm | 70-75 bpm | 307 | 70-75 bpm = 1975; <60 bpm = 930; 60-70 bpm = 2635; 75-85 bpm = 2374; 85+ bpm = 1441 | NHANES Examination (BPX_F, BPX_G) - BPXPLS |
| Asthma | Ever diagnosed with asthma from MCQ010. Binary: Asthma (ASTHMA=1) vs control (all outcome indicators=0) | control | 2798 | control = 5333; Asthma = 1531 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable ASTHMA |
| Bronchitis | Ever diagnosed with chronic bronchitis from MCQ160K. Binary: Bronchitis (BRONCHITIS=1) vs control (all outcome indicators=0) | control | 3925 | control = 5333; Bronchitis = 404 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable BRONCHITIS |
| Emphysema | Ever diagnosed with emphysema from MCQ160L. Binary: Emphysema (EMPHYSEMA=1) vs control (all outcome indicators=0) | control | 4215 | control = 5333; Emphysema = 114 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable EMPHYSEMA |
| Angina | Ever diagnosed with angina pectoris from MCQ160C. Binary: Angina (ANGINA=1) vs control (all outcome indicators=0) | control | 4203 | control = 5333; Angina = 126 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable ANGINA |
| Heart_failure | Ever diagnosed with congestive heart failure from MCQ160B. Binary: Heart failure (HEART_FAILURE=1) vs control (all outcome indicators=0) | control | 4149 | control = 5333; Heart failure = 180 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable HEART_FAILURE |
| Heart_attack | Ever diagnosed with heart attack (myocardial infarction) from MCQ160E. Binary: Heart attack (HEART_ATTACK=1) vs control (all outcome indicators=0) | control | 4124 | control = 5333; Heart attack = 205 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable HEART_ATTACK |
| Stroke | Ever diagnosed with stroke from MCQ160F. Binary: Stroke (STROKE=1) vs control (all outcome indicators=0) | control | 4158 | control = 5333; Stroke = 171 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable STROKE |
| Diabetes | Ever diagnosed with diabetes from DIQ010. Binary: Diabetes (DIABETES=1) vs control (all outcome indicators=0) | control | 3352 | control = 5333; Diabetes = 977 | NHANES Questionnaire (DIQ_F, DIQ_G) - derived variable DIABETES |
| CHD | Ever diagnosed with coronary heart disease from MCQ160D. Binary: CHD (CHD=1) vs control (all outcome indicators=0) | control | 4192 | control = 5333; CHD = 137 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CHD |
| CVD | Any cardiovascular disease (composite of heart attack, CHD, angina, stroke, heart failure). Binary: CVD (CVD=1) vs control (all outcome indicators=0) | control | 3818 | control = 5333; CVD = 511 | NHANES Questionnaire (MCQ_F, MCQ_G) - composite derived variable |
| Breast_cancer | Ever diagnosed with breast cancer from MCQ220. Binary: Breast cancer (CANCER_BREAST=1) vs control (all outcome indicators=0). Female participants only | control | 4273 | control = 5333; Breast cancer = 56 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_BREAST |
| Colon_cancer | Ever diagnosed with colon cancer from MCQ220. Binary: Colon cancer (CANCER_COLON=1) vs control (all outcome indicators=0) | control | 4308 | control = 5333; Colon cancer = 21 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_COLON |
| Lung_cancer | Ever diagnosed with lung cancer from MCQ220. Binary: Lung cancer (CANCER_LUNG=1) vs control (all outcome indicators=0) | control | 4324 | control = 5333; Lung cancer = 5 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_LUNG |
| Esophageal_cancer | Ever diagnosed with esophageal cancer from MCQ220. Binary: Esophageal cancer (CANCER_ESOPHAGEAL=1) vs control (all outcome indicators=0) | control | 4327 | control = 5333; Esophageal cancer = 2 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_ESOPHAGEAL |
| Prostate_cancer | Ever diagnosed with prostate cancer from MCQ220. Binary: Prostate cancer (CANCER_PROSTATE=1) vs control (all outcome indicators=0). Male participants only | control | 4276 | control = 5333; Prostate cancer = 53 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_PROSTATE |
| Mouth_cancer | Ever diagnosed with oral/pharyngeal cancer from MCQ220. Binary: Mouth cancer (CANCER_MOUTH=1) vs control (all outcome indicators=0) | control | 4329 | control = 5333; Mouth cancer = 0 | NHANES Questionnaire (MCQ_F, MCQ_G) - derived variable CANCER_MOUTH |

---

## R Code Definitions

For reproducibility, the exact R code used to create each categorical variable is provided below:

### Gender
```r
if("Gender" %in% colnames(.)) fct_relevel(as.factor(Gender), "Female")
```

### Age_group
```r
if("AgeGroup" %in% colnames(.)) fct_relevel(as.factor(AgeGroup), "30-39")
```

### Education_level
```r
factor(EducationLevel, levels = c("< 9th Grade", "9-11th Grade", "High School", "College/AA", "College Graduate")) %>% fct_relevel("College/AA")
```

### Ethnicity
```r
if("Ethnicity" %in% colnames(.)) fct_relevel(as.factor(Ethnicity), "White")
```

### US_born
```r
if("US_Born" %in% colnames(.)) fct_relevel(as.factor(US_Born), "US Born")
```

### Household_size
```r
factor(Household_Size_Factor, levels = c("1","2","3","4","5","6","7+")) %>% fct_relevel("4")
```

### Marital_status
```r
if("Marital_Status" %in% colnames(.)) fct_relevel(as.factor(Marital_Status), "Married")
```

### Interview_language
```r
if("Interview_Language" %in% colnames(.)) fct_relevel(as.factor(Interview_Language), "English")
```

### Income_to_poverty_ratio
```r
case_when(INDFMPIR < 0.50 ~ "Below 50%", INDFMPIR >= 0.50 & INDFMPIR < 1.00 ~ "50-99%", INDFMPIR >= 1.00 & INDFMPIR < 1.25 ~ "100-124%", INDFMPIR >= 1.25 & INDFMPIR < 1.50 ~ "125-149%", INDFMPIR >= 1.50 & INDFMPIR < 1.85 ~ "150-184%", INDFMPIR >= 1.85 & INDFMPIR < 2.00 ~ "185-199%", INDFMPIR >= 2.00 ~ "200% and Over", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("150-184%")
```

### Denture
```r
case_when(DENTURE_OHAROCDE == 1 ~ "Denture", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Denture"))
```

### Gum_disease
```r
case_when(GUM_DISEASE_OHAROCGP == 1 ~ "Gum disease", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Gum disease"))
```

### Oral_hygiene
```r
case_when(ORAL_HYGIENE_OHAROCOH == 1 ~ "Poor oral hygiene", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Poor oral hygiene"))
```

### Tooth_decay
```r
case_when(TOOTH_DECAY_OHAROCDT == 1 ~ "Tooth decay", SEQN %in% SEQN_oradWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Tooth decay"))
```

### Smoking_status
```r
case_when(SMQ_current_ever_never == 0 ~ "Never smoker", SMQ_current_ever_never == 2 ~ "Former smoker", SMQ_current_ever_never == 1 ~ "Current smoker", TRUE ~ NA_character_) %>% factor(levels = c("Never smoker", "Former smoker", "Current smoker")) %>% fct_relevel("Never smoker")
```

### Hepatitis_C_antibody
```r
case_when(LBDHCV == 1 ~ "Positive", LBDHCV == 2 ~ "Negative", TRUE ~ NA_character_) %>% factor(levels = c("Negative", "Positive"))
```

### HPV_PCR_summary
```r
case_when(LBDRPCR == 1 ~ "Positive", LBDRPCR == 2 ~ "Negative", TRUE ~ NA_character_) %>% factor(levels = c("Negative", "Positive"))
```

### BMI_category
```r
case_when(BMXBMI < 18.5 ~ "Underweight", BMXBMI >= 18.5 & BMXBMI < 25 ~ "Healthy weight", BMXBMI >= 25 & BMXBMI < 30 ~ "Overweight", BMXBMI >= 30 & BMXBMI < 35 ~ "Class 1 Obesity", BMXBMI >= 35 ~ "Class 2-3 Obesity", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("Healthy weight")
```

### Blood_pressure
```r
case_when(MSYSTOLIC > 180 | MDIASTOLIC > 120 ~ "Hypertensive Crisis", MSYSTOLIC >= 140 | MDIASTOLIC >= 90 ~ "Hypertension Stage 2", (MSYSTOLIC >= 130 & MSYSTOLIC < 140) | (MDIASTOLIC >= 80 & MDIASTOLIC < 90) ~ "Hypertension Stage 1", (MSYSTOLIC >= 120 & MSYSTOLIC < 130) & (MDIASTOLIC < 80) ~ "Elevated", MSYSTOLIC < 120 & MDIASTOLIC < 80 ~ "Normal", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("Normal")
```

### Pulse_category
```r
case_when(BPXPLS < 60 ~ "<60 bpm", BPXPLS >= 60 & BPXPLS < 70 ~ "60-70 bpm", BPXPLS >= 70 & BPXPLS < 75 ~ "70-75 bpm", BPXPLS >= 75 & BPXPLS < 85 ~ "75-85 bpm", BPXPLS >= 85 ~ "85+ bpm", TRUE ~ NA_character_) %>% factor() %>% fct_relevel("70-75 bpm")
```

### Asthma
```r
case_when(ASTHMA == 1 ~ "Asthma", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Asthma"))
```

### Bronchitis
```r
case_when(BRONCHITIS == 1 ~ "Bronchitis", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Bronchitis"))
```

### Emphysema
```r
case_when(EMPHYSEMA == 1 ~ "Emphysema", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Emphysema"))
```

### Angina
```r
case_when(ANGINA == 1 ~ "Angina", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Angina"))
```

### Heart_failure
```r
case_when(HEART_FAILURE == 1 ~ "Heart failure", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Heart failure"))
```

### Heart_attack
```r
case_when(HEART_ATTACK == 1 ~ "Heart attack", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Heart attack"))
```

### Stroke
```r
case_when(STROKE == 1 ~ "Stroke", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Stroke"))
```

### Diabetes
```r
case_when(DIABETES == 1 ~ "Diabetes", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Diabetes"))
```

### CHD
```r
case_when(CHD == 1 ~ "CHD", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "CHD"))
```

### CVD
```r
case_when(CVD == 1 ~ "CVD", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "CVD"))
```

### Breast_cancer
```r
case_when(CANCER_BREAST == 1 ~ "Breast cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Breast cancer"))
```

### Colon_cancer
```r
case_when(CANCER_COLON == 1 ~ "Colon cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Colon cancer"))
```

### Lung_cancer
```r
case_when(CANCER_LUNG == 1 ~ "Lung cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Lung cancer"))
```

### Esophageal_cancer
```r
case_when(CANCER_ESOPHAGEAL == 1 ~ "Esophageal cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Esophageal cancer"))
```

### Prostate_cancer
```r
case_when(CANCER_PROSTATE == 1 ~ "Prostate cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Prostate cancer"))
```

### Mouth_cancer
```r
case_when(CANCER_MOUTH == 1 ~ "Mouth cancer", SEQN %in% SEQN_outWAS_control ~ "control", TRUE ~ NA_character_) %>% factor(levels = c("control", "Mouth cancer"))
```


## Notes

1. **Control Groups**: For disease outcome variables (Oral-WAS, Outcome-WAS), the 'control' group consists of participants with ALL outcome indicators equal to 0 (no missing values allowed). This ensures a true disease-free comparison group.

2. **Reference Levels**: Reference levels were chosen based on:
   - Most common category (e.g., 'Married' for marital status)
   - Clinical baseline (e.g., 'Healthy weight' for BMI, 'Normal' for blood pressure)
   - Unexposed group (e.g., 'Never smoker' for smoking status, 'control' for disease outcomes)

3. **Missing Data**: Variables with 'Indeterminate' or insufficient sample results are coded as NA and excluded from analysis.

4. **Data Source**: All variables derived from NHANES 2009-2010 (cycles F) and 2011-2012 (cycles G).

---

*Table generated on: 2025-10-14*
*Total variables documented: 35*
*Total participants: 9662*
