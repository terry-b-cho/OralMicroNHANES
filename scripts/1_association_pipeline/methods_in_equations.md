# NHANES Oral Microbiome Analysis Methods

## Methods

### 1. Study Population & Data Sources

We analyzed data from the National Health and Nutrition Examination Survey (NHANES), a continuous program designed by the National Center for Health Statistics (NCHS) that uses a complex, stratified, multistage probability design to obtain representative samples of the civilian, non-institutionalized United States population. Oral‐rinse 16S rRNA gene sequencing was conducted during two survey cycles: 2009–2010 (cycle "F") and 2011–2012 (cycle "G"). All raw demultiplexed sequencing reads, NHANES examination data, questionnaire responses, and laboratory measurements were downloaded from the publicly available NCHS website and linked by the unique participant identifier SEQN.

Our analytic population included participants from both cycles who possessed complete data for oral microbiome analysis. After excluding participants who lacked mobile examination center weights (WTMEC2YR), masked variance strata identifiers (SDMVSTRA), masked primary sampling unit identifiers (SDMVPSU), or any mandatory covariates described below, our final analytic cohort comprised 9,847 participants with oral microbiome data from both cycles combined. The exact participant counts per cycle are automatically determined by the database query system and reported in supplementary materials to ensure accuracy across different analysis runs.

All microbiome count tables, demographic variables, and original NHANES XPT files were integrated into a single SQLite database to ensure reproducible, queryable workflows. This database architecture enabled consistent data access across all analytical pipelines while maintaining data integrity and version control.

### 2. Microbiome Processing

Paired-end 16S rRNA gene sequencing reads were processed using the DADA2 algorithm version 1.24.0 with default parameters and maximum expected errors set to 2. Reads were quality-filtered, denoised, and checked for chimeric sequences before taxonomic assignment using the Greengenes reference database (release 13-8). Amplicon sequence variants (ASVs) were systematically collapsed to five taxonomic ranks: phylum, class, order, family, and genus level, yielding comprehensive taxonomic profiles at each hierarchical level.

Importantly, we retained all detected taxa without implementing prevalence thresholds or abundance filtering. This approach was essential for preserving the compositional structure required for centered log-ratio (CLR) transformation, which demands the complete compositional basis for mathematically valid inference. The retention of rare taxa ensures that our analytical framework can accommodate the full spectrum of microbiome diversity while maintaining the geometric properties of compositional data.

For each taxonomic rank and survey cycle, we generated both raw count matrices and relative abundance matrices, which were subsequently transformed using four distinct normalization approaches as detailed below.

### 3. Normalization & Transformations

Let $i \in \{1, \ldots, N\}$ denote individual participants and $j \in \{1, \ldots, D\}$ denote taxa. We define $C_{ij}$ as the raw count for participant $i$ and taxon $j$, with total library size $n_i = \sum_{j=1}^{D} C_{ij}$ and relative abundance $P_{ij} = C_{ij}/n_i$. A fixed integer pseudo-count $\varepsilon = 1$ was added where logarithmic transformations were applied to avoid numerical instabilities from zero counts.

Zero-library samples were identified using the criterion $n_i = 0$ and systematically replaced with missing values across all taxa before transformation to prevent computational failures. We verified closure properties of relative abundance matrices by confirming that row sums deviated from unity by less than $10^{-6}$ for all non-missing samples.

We implemented four distinct normalization strategies to accommodate different analytical assumptions and interpretative frameworks:

**Proportional (none) transformation:** $T_{ij}^{\mathrm{none}} = P_{ij}$. This transformation preserves data on the probability simplex, maintaining direct interpretability as relative abundances. Effect estimates represent absolute percentage-point changes in microbial proportions, facilitating comparison with traditional epidemiological measures expressed on proportion scales. For computational efficiency, these matrices were implemented as SQL database views rather than materialized tables.

**Hellinger transformation:** $T_{ij}^{\mathrm{hel}} = \sqrt{P_{ij}}$. The square-root transformation provides variance stabilization while preserving Euclidean geometry. This approach down-weights extremely rare taxa while maintaining compositional relationships, with Euclidean distances between transformed samples corresponding to Hellinger divergences between the original compositional vectors.

**Centered log-ratio (CLR) transformation:** We first computed the geometric mean $g_i = \left(\prod_{j=1}^{D}(C_{ij} + \varepsilon)\right)^{1/D}$ for each sample, then applied the transformation $T_{ij}^{\mathrm{clr}} = \ln\left(\frac{C_{ij} + \varepsilon}{g_i}\right)$. CLR transformation removes the unit-sum constraint inherent in compositional data, placing transformed values in $\mathbb{R}^{D-1}$ space. Effect estimates represent natural log-fold changes relative to the geometric mean of all taxa within each sample, providing compositionally coherent inference.

**Log-normal transformation:** $T_{ij}^{\mathrm{ln}} = \log_{10}\left(\frac{C_{ij} + \varepsilon}{n_i + D\varepsilon} \cdot \bar{n}\right)$, where $\bar{n} = N^{-1}\sum_{i=1}^{N} n_i$ represents the global mean library size. This approach first normalizes by individual sequencing depth, then rescales to the population-average library size. Base-10 logarithms yield effect estimates interpretable as log₁₀ fold-changes per mean-depth unit, with the intercept remaining invariant to library size variations.

All transformations were implemented with optimized numerical precision using vectorized operations in R version 4.3.3, with comprehensive error handling for edge cases including infinite values and numerical overflow conditions.

### 4. NHANES Survey Design Specification

Our analytical approach incorporated the complex survey design features of NHANES following NCHS Technical Documentation guidelines (2006). We implemented both single-cycle and pooled-cycle analysis strategies depending on data availability for each variable pair, with automatic detection of optimal analysis modes using a cycle availability algorithm.

**Cycle availability detection:** For each dependent-independent variable pair $(Y, X)$, we computed the set of available cycles $\mathcal{C}_{Y,X} = \{c \in \{\text{F}, \text{G}\} : \text{data available for both } Y \text{ and } X \text{ in cycle } c\}$. Pooled analysis was implemented when $|\mathcal{C}_{Y,X}| = 2$, otherwise single-cycle analysis was performed using $c = \mathcal{C}_{Y,X}$.

**Weight validation procedure:** Before each analysis, we verified the integrity of survey weights using the validation function $V(w) = \mathbf{1}_{w > 0} \cap \mathbf{1}_{w \neq \text{NA}}$ for all weight vectors $w \in \{\text{WTMEC2YR}_F, \text{WTMEC2YR}_G\}$. Analysis proceeded only when $\sum V(w) > 0$ for all required demographic tables.

**Single-cycle analysis:** When $|\mathcal{C}_{Y,X}| = 1$, we created survey design objects using the original 2-year mobile examination center weights (WTMEC2YR), masked variance strata (SDMVSTRA), and masked primary sampling units (SDMVPSU). The survey design was specified as:

$$\texttt{svydesign}(\text{ids} = \sim\text{SDMVPSU}, \text{strata} = \sim\text{SDMVSTRA}, \text{weights} = \sim\text{WTMEC2YR}, \text{nest} = \text{TRUE})$$

**Pooled-cycle analysis:** When $|\mathcal{C}_{Y,X}| = 2$, we implemented NCHS-compliant 4-year pooled analysis following established protocols for combining non-overlapping survey cycles. We computed 4-year weights as $\text{WTMEC4YR} = \text{WTMEC2YR}/2$, representing the average of the 2-year weights across the pooled period. To ensure unique identification across cycles, we created composite design identifiers using the transformations:

$$\text{unique\_psu} = f(\text{cycle}) \times 1000 + \text{SDMVPSU}$$
$$\text{unique\_strata} = f(\text{cycle}) \times 1000 + \text{SDMVSTRA}$$

where $f(\text{cycle})$ represents a numeric encoding function with $f(\text{F}) = 1$ and $f(\text{G}) = 2$.

The pooled survey design was specified as:

$$\texttt{svydesign}(\text{ids} = \sim\text{unique\_psu}, \text{strata} = \sim\text{unique\_strata}, \text{weights} = \sim\text{WTMEC4YR}, \text{nest} = \text{TRUE})$$

Variance estimation employed Taylor series linearization with the "certainty" option for handling strata containing single primary sampling units, following NCHS recommendations for mobile examination center weights. Our pooled-cycle approach was mathematically justified by the non-overlapping nature of the survey cycles and the stability of the sampling frame between 2009 and 2012, with the exception of enhanced Asian-American oversampling introduced in cycle G.

### 5. Regression Pipelines

We implemented six distinct Weighted Association Study (WAS) frameworks, each addressing specific epidemiological questions about relationships between microbiome composition and health outcomes. Let $T_{ij}^{(q)}$ denote the transformed abundance of taxon $j$ for participant $i$ under normalization $q \in \{\text{none}, \text{hel}, \text{clr}, \text{ln}\}$, $X_i$ represent non-microbial exposures, $Y_i$ denote health outcomes, and $\mathbf{C}_i$ represent the covariate vector.

**Progressive covariate selection algorithm:** We implemented a hierarchical covariate selection procedure to ensure model stability while maximizing statistical power. Let $\mathbf{F} = \{\text{RIDAGEYR}, \text{AGE\_SQUARED}, \text{RIAGENDR}, \text{INDFMPIR}\}$ denote the full covariate set and $\mathbf{E} = \{\text{RIDAGEYR}, \text{RIAGENDR}\}$ denote essential covariates. For each covariate $c \in \mathbf{F}$, we computed the missing rate $m_c = N^{-1}\sum_{i=1}^{N} \mathbf{1}_{c_i = \text{NA}}$ and variability criterion $v_c$ defined as the number of unique non-missing values. Covariates were retained if $m_c \leq 0.3$ and $v_c > 1$, with essential covariates always included regardless of these criteria.

**Adaptive binning strategy:** For continuous variables with excessive ties, we implemented a hierarchical binning procedure. Let $x$ be a continuous predictor and $Q_k(x, \alpha)$ denote the quantile-based binning function attempting to create $k$ bins with minimum bin size $\alpha n$. We applied the sequence: quartiles $Q_4(x, 0.05)$, tertiles $Q_3(x, 0.05)$, and finally log-continuous transformation $\log_{10}(x + \min(x[x > 0]))$ if binning failed. This ensures robust handling of tied observations while preserving ordinal relationships.

**Minimum sample size criteria:** We enforced dynamic minimum sample sizes based on model complexity. For linear models, we required $n \geq \max(15, k + 5)$ where $k$ is the number of covariates. For logistic models, we required $n \geq \max(10, k + 5)$. These criteria ensure adequate power for coefficient estimation and convergence stability.

**Model framework implementations:**

**Demographic determinants of microbial abundance (1\_demoWAS):** This pipeline examined how demographic characteristics influence microbiome composition using linear models with identity link:

$$\mathbb{E}[T_{ij}^{(q)}] = \beta_{0j}^{(q)} + \beta_{1j}^{(q)} X_i + \mathbf{C}_i^{\top} \boldsymbol{\gamma}_j^{(q)}$$

where demographic variables $X_i$ included age, sex, race/ethnicity, education, and income measures, with separate models fitted for each taxon as the dependent variable.

**Microbial predictors of oral health outcomes (2\_oradWAS):** This framework investigated associations between individual taxa and binary oral health endpoints using logistic regression:

$$\text{Pr}(Y_i = 1) = \left[1 + \exp\left(-\beta_0^{(q)} - \beta_1^{(q)} T_{ij}^{(q)} - \mathbf{C}_i^{\top} \boldsymbol{\gamma}^{(q)}\right)\right]^{-1}$$

where $Y_i$ represented binary oral health outcomes such as tooth decay, denture use, or gum disease.

**Environmental exposures influencing microbial abundance (3\_exWAS):** This pipeline used the same mathematical framework as 1\_demoWAS but with environmental and dietary exposures $X_i$ replacing demographic predictors:

$$\mathbb{E}[T_{ij}^{(q)}] = \beta_{0j}^{(q)} + \beta_{1j}^{(q)} X_i + \mathbf{C}_i^{\top} \boldsymbol{\gamma}_j^{(q)}$$

**Microbial predictors of continuous phenotypes (4\_pheWAS):** This framework examined associations between taxa and continuous health measures using linear models:

$$\mathbb{E}[Y_i] = \beta_0^{(q)} + \beta_1^{(q)} T_{ij}^{(q)} + \mathbf{C}_i^{\top} \boldsymbol{\gamma}^{(q)}$$

where $Y_i$ represented continuous phenotypes such as body mass index, blood pressure, or laboratory values.

**Microbial predictors of systemic disease outcomes (5\_outWAS):** This pipeline employed the same logistic framework as 2\_oradWAS but focused on systemic disease endpoints:

$$\text{Pr}(Y_i = 1) = \left[1 + \exp\left(-\beta_0^{(q)} - \beta_1^{(q)} T_{ij}^{(q)} - \mathbf{C}_i^{\top} \boldsymbol{\gamma}^{(q)}\right)\right]^{-1}$$

where $Y_i$ represented binary indicators for diseases such as diabetes, cardiovascular disease, or cancer.

**Microbial predictors of laboratory biomarkers (6\_zimWAS):** This framework used linear models identical to 4\_pheWAS but specifically targeted laboratory measurements:

$$\mathbb{E}[Y_i] = \beta_0^{(q)} + \beta_1^{(q)} T_{ij}^{(q)} + \mathbf{C}_i^{\top} \boldsymbol{\gamma}^{(q)}$$

**Model convergence and quality control:** We implemented enhanced convergence monitoring using adaptive iteration limits. Initial models used standard convergence criteria (maximum 100 iterations, tolerance $10^{-8}$). For non-convergent models, we applied relaxed criteria (maximum 200 iterations, tolerance $10^{-6}$) before declaring convergence failure. Infinite values in predictors or outcomes were systematically replaced with missing values before model fitting.

**Effect size calculations:** We computed pseudo-R² for logistic models using $R^2_{\text{pseudo}} = 1 - \text{deviance}/\text{null.deviance}$ and survey-weighted R² for linear models using the formula $R^2_{\text{survey}} = 1 - \text{RSS}_w/\text{TSS}_w$ where weighted sums of squares were computed using survey design weights extracted via the survey package's weight accessor functions.

**Multiple comparisons correction:** We applied the Benjamini-Hochberg false discovery rate (FDR) procedure within each dependent variable to control the expected proportion of false discoveries. For each dependent variable $Y$ and transformation $q$, let $\{p_1, p_2, \ldots, p_m\}$ be the set of p-values for all tested associations. The FDR-adjusted p-values were computed as $\tilde{p}_i = \min\left(1, \min_{j \geq i} \frac{m \cdot p_{(j)}}{j}\right)$ where $p_{(1)} \leq p_{(2)} \leq \ldots \leq p_{(m)}$ are the ordered p-values. We considered associations statistically significant when $\tilde{p}_i < 0.1$.

**Effect scale harmonization:** To facilitate comparison across transformations, we implemented systematic effect scale documentation. Each result was annotated with transformation-specific scale information: proportion (0-1) for none transformation, sqrt-proportion (0-1) for Hellinger, ln-ratio (centered) for CLR, and log₁₀-CPM for log-normal transformation. Interpretation guides were automatically generated to specify the units and mathematical meaning of each effect size estimate.

### 6. Computational Reproducibility

Our complete analytical workflow was implemented in R version 4.3.3 using the survey package version 4.4-2 for complex survey analysis, with all code maintained under version control in a public Git repository. The workflow employed a modular architecture with three main computational stages executed through SLURM job scheduling on high-performance computing infrastructure.

**Data transformation stage:** A dedicated transformation script implements the four normalization strategies using fully vectorized operations for computational efficiency. The script employs a systematic approach: reading raw count and relative abundance matrices from the SQLite database, applying each normalization method with comprehensive numerical precision safeguards, and writing transformed matrices back to the database with systematic naming conventions. For the proportional transformation, SQL views are created rather than duplicate tables, achieving approximately 50% storage reduction while maintaining data access consistency. Zero-library samples are identified and handled gracefully through the transformation pipeline without causing computational failures.

**Schema generation stage:** Twenty-four parallel jobs generate comprehensive variable-pair mapping files, representing each combination of the six analytical pipelines and four normalization methods. These jobs create detailed specifications of all dependent-independent variable combinations to be analyzed, incorporating metadata about data availability across cycles, variable type classifications, and expected analytical frameworks. The schema generation process implements the cycle availability detection algorithm described above, enabling automatic selection between pooled and single-cycle analysis modes for each variable pair.

**Main analytical stage:** The primary analytical pipeline distributes 11,448 individual survey-weighted regression models across the computing cluster infrastructure. Each job processes a single dependent variable across all relevant independent variables and transformations using the progressive covariate selection and adaptive binning algorithms described above. Models are fitted using the survey-weighted generalized linear model framework with comprehensive error handling including convergence monitoring, infinite value replacement, and progressive fallback strategies for problematic model specifications.

**Output standardization:** Each analytical job generates a structured RDS file containing three primary data components following a standardized format. The pe_tidied component contains taxon-level coefficient estimates, standard errors, test statistics, and both raw and FDR-corrected p-values with transformation-specific effect scale annotations. The pe_glanced component provides model-level diagnostic statistics including sample sizes, degrees of freedom, convergence indicators, and goodness-of-fit measures. The rsq component contains survey-weighted or pseudo-R² values computed using the appropriate methods for each model family, along with effect scale harmonization metadata.

**Quality assurance framework:** Comprehensive logging systems capture all computational decisions including sample size checks, convergence diagnostics, covariate reduction sequences, and adaptive binning outcomes. Each analysis records the exact formula used, number of covariates successfully included, transformation method applied, cycle mode implemented, and any fallback strategies invoked during model fitting. A fixed random seed is established for all stochastic processes to ensure reproducibility across computational environments.

**Environmental documentation:** Complete computational environment information is captured using package management tools to document exact software versions, system dependencies, and hardware configurations. The integrated SQLite database architecture ensures that identical data matrices are accessed consistently across all analytical components, eliminating potential discrepancies from file-based data sharing approaches. Combined with version-controlled analysis scripts and comprehensive logging systems, this framework enables complete reproduction of results from raw NHANES XPT files through final publication-ready statistical outputs.

**Computational efficiency optimizations:** The workflow incorporates several performance optimizations including vectorized transformation operations, efficient SQL query patterns for large-scale data access, and intelligent job scheduling to minimize redundant computations. Memory usage is optimized through strategic data loading patterns and garbage collection procedures, enabling analysis of large-scale microbiome datasets within standard high-performance computing resource constraints. 