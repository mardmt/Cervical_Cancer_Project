# CERVICAL CANCER ANALYSIS - METHODOLOGY & INTERPRETATION GUIDE

---

## Statistical Methods Overview

This analysis uses standard epidemiological and biostatistical methods to (1) describe cervical cancer burden and risk factors, and (2) build a predictive model for cancer diagnosis.

---

## PART 1: EXPLORATORY DATA ANALYSIS (EDA)

### 1.1 Missing Data Assessment

**Method:** Count and percentage of missing values per column

**Interpretation:**
- Missing < 5%: Negligible; can ignore or use listwise deletion
- Missing 5–20%: Moderate; consider multiple imputation or variable exclusion
- Missing > 20%: Substantial; use domain knowledge to decide on inclusion
- Missing = "?": Original dataset indicator; treat as missing

**Output:** `01_missing_values.csv`

---

### 1.2 Descriptive Statistics (Continuous Variables)

**Method:**
- **Central Tendency:** Mean, Median
- **Dispersion:** Standard Deviation (SD), Interquartile Range (IQR), Min–Max Range
- **Interpretation:** 
  - Mean ± SD describes the typical patient
  - Large SD relative to mean indicates high variability
  - Median more robust to outliers than mean

**Example Interpretation:**
- Age: Mean = 26.9 years, SD = 8.4 → Most patients aged 18–35
- Sexual partners: Mean = 2.5, SD = 3.2 → Right-skewed; median more informative

**Output:** `02_descriptive_statistics.csv`

---

### 1.3 Categorical Data Distribution

**Method:** Frequency counts and percentages

**Interpretation:**
- Prevalence of sexual history, STD status, contraceptive use
- Baseline cancer diagnosis rate in sample

**Output:** Included in EDA plots and summary files

---

### 1.4 Visual Exploratory Analysis

**Plots:**

1. **Age Distribution (Histogram)**
   - Interpretation: Normal-like vs. skewed distribution
   - Check for outliers or unusual clusters

2. **Global Incidence Rate Distribution (Histogram)**
   - Interpretation: Right-skewed; most countries have low rates, few have very high rates
   - Expected range: 2–25 per 100,000 women

3. **Top 10 Countries by Burden (Bar Chart)**
   - Interpretation: Geographic disparities
   - Pattern: Sub-Saharan Africa and South Asia dominate

4. **Cancer Diagnosis Proportion (Pie Chart)**
   - Interpretation: Class imbalance (cancer is rare in this sample)
   - Expected: ~10% cancer, ~90% no cancer (imbalanced)

**Output:** `07_eda_plots.pdf`

---

## PART 2: RISK STRATIFICATION MODELING

### 2.1 Univariate Logistic Regression

**Model:**
```
log(odds of cancer) = β₀ + β₁ × Risk_Factor
```

**Interpretation of Results:**

| Term | Meaning | Example |
|------|---------|---------|
| **Odds Ratio (OR)** | Multiplicative change in odds per unit increase | OR = 1.5 → 50% increase in cancer odds per unit |
| **OR > 1** | Risk factor increases cancer odds | Smoking increases odds |
| **OR < 1** | Protective factor | Screening programs decrease odds |
| **OR = 1** | No association | p-value will be high |
| **95% CI** | Confidence range; if crosses 1.0, not significant at p<0.05 | CI: [0.9–1.1] crosses 1 → not significant |
| **p-value** | Statistical significance | p < 0.05 → significant association |

**Example Interpretation:**
- Variable: "Age"
- OR = 1.02, 95% CI [1.00–1.04], p = 0.04
- *Each additional year of age increases cancer odds by 2%; significant at p < 0.05*

**Output:** `08_univariate_logistic_regression.csv`

---

### 2.2 Multivariate Logistic Regression

**Model:**
```
log(odds of cancer) = β₀ + β₁×Age + β₂×STD + β₃×Sexual_Partners + ...
```

**Why Multivariate?**
- Controls for confounding (e.g., age confounds the STD–cancer relationship)
- Provides independent effect of each variable
- Builds a predictive model

**Model Selection (Stepwise AIC):**
- Start with all significant univariate variables
- Iteratively remove/add variables to minimize AIC
- AIC balances fit and simplicity (penalizes overfitting)
- Lower AIC = Better model

**Output:** `09_multivariate_logistic_regression.csv`

**Example Interpretation:**
- Model 1 (Initial): AIC = 450
- Model 2 (Stepwise): AIC = 420
- *Model 2 is better (lower AIC); dropped 3 less-important variables*

---

### 2.3 Diagnostic Test Accuracy

**Context:** Four diagnostic tests for cervical cancer are evaluated:
1. Hinselmann (colposcopy-based)
2. Schiller (staining test)
3. Cytology (Pap smear)
4. Biopsy (gold standard pathology)

**Metrics:**

| Metric | Formula | Interpretation |
|--------|---------|-----------------|
| **Sensitivity** | TP/(TP+FN) | "Of patients WITH cancer, what % test positive?" |
| **Specificity** | TN/(TN+FP) | "Of patients WITHOUT cancer, what % test negative?" |
| **PPV (Pos Pred Value)** | TP/(TP+FP) | "If test is positive, probability of cancer?" |
| **NPV (Neg Pred Value)** | TN/(TN+FN) | "If test is negative, probability of NO cancer?" |
| **AUC (ROC)** | Area under curve | Overall discriminatory ability; 0.5=random, 1.0=perfect |

**Clinical Interpretation:**

- **High Sensitivity + High Specificity = Good Test**
  - Few false negatives (diagnoses most cancers)
  - Few false positives (avoids unnecessary treatment)

- **High Sensitivity, Low Specificity = Screening Test**
  - Catches most cancers (few miss) but many false alarms
  - Example: Cytology (Pap smear)

- **Low Sensitivity, High Specificity = Confirmation Test**
  - Few false positives; confirms diagnosis
  - Example: Biopsy

**Example Interpretation:**
- Cytology: Sensitivity = 0.95, Specificity = 0.70, AUC = 0.85
- *Cytology detects 95% of cancers but has many false positives; good for screening but needs confirmation*

- Biopsy: Sensitivity = 0.88, Specificity = 0.98, AUC = 0.95
- *Biopsy is highly specific (few false positives); best for confirmation*

**ROC Curve:**
- X-axis: 1 − Specificity (False Positive Rate)
- Y-axis: Sensitivity (True Positive Rate)
- **AUC = 0.5:** Random guessing
- **AUC = 0.7–0.8:** Acceptable discrimination
- **AUC = 0.8–0.9:** Excellent discrimination
- **AUC > 0.9:** Outstanding discrimination

**Output:** 
- `11_diagnostic_test_accuracy.csv`
- `10_diagnostic_roc_curves.pdf`

---

## PART 3: INTERPRETATION OF KEY FINDINGS

### Geographic Disparities

**Finding:** Cervical cancer incidence ranges from ~2 to ~25 per 100,000 women across countries.

**Interpretation:**
- 12-fold difference between highest and lowest burden countries
- Disparities linked to:
  - Screening program coverage (low in low-income countries)
  - HPV vaccination uptake (high in high-income countries)
  - Healthcare infrastructure

**Data Source:** GLOBOCAN 2022; Our World in Data

---

### Screening Program Coverage

**Finding:** ~46% of countries report organized national screening programs.

**Interpretation:**
- Screening programs concentrated in high-income countries
- Disparity: Low-income countries have 85% of cervical cancer burden but only ~30% of screening programs
- Consequence: Preventable cancers go undetected in under-resourced settings

**Data Source:** WHO Global Health Observatory

---

### Patient-Level Risk Factors

**Expected Key Findings (from literature):**

1. **Sexual History**
   - Early age of first intercourse → increased risk (earlier HPV exposure)
   - Multiple sexual partners → increased risk (more HPV exposure)
   - Interpretation: HPV transmission risk accumulates with sexual exposure

2. **STD History**
   - HPV infection → strong predictor of cancer (necessary but not sufficient cause)
   - Other STDs (Herpes, Chlamydia) → associated but indirect effect
   - Interpretation: HPV is the causal factor; other STDs may be markers of risky behavior

3. **Smoking**
   - Smoking → increased cancer risk (may impair immune response)
   - Interpretation: Immunosuppression enables HPV persistence

4. **Contraceptive Use**
   - Hormonal contraceptives → possible slight increase in risk (debated)
   - IUD use → possible protective effect (inflammatory response)
   - Interpretation: Biological mechanisms not fully understood; findings inconsistent across studies

---

### Diagnostic Test Comparison

**Expected Ranking (from literature):**

1. **Most Sensitive:** Cytology (Pap smear) – detects most cancers
2. **Most Specific:** Biopsy – confirms diagnosis with high certainty
3. **Balanced:** Hinselmann and Schiller – intermediate on both metrics

**Clinical Implication:**
- Screening strategy: Use Cytology to cast wide net, then Biopsy to confirm
- Algorithm: Cytology → HPV test (reflex) → Colposcopy → Biopsy if high-grade lesion

---

## PART 4: STATISTICAL ASSUMPTIONS & LIMITATIONS

### Logistic Regression Assumptions

1. **Binary Outcome:** Cancer diagnosis is Yes/No ✓
2. **Independence:** Each patient is independent observation ✓
3. **No Perfect Separation:** Both outcomes present in data ✓
4. **Linearity of Log-Odds:** Assumes linear relationship between log-odds and predictors
   - Check with Box-Tidwell test (not done here but recommended)
5. **No Multicollinearity:** Predictors not highly correlated
   - Check with VIF (Variance Inflation Factor; not done here but recommended)

### Limitations of This Dataset

1. **Class Imbalance:** ~10% cancer, ~90% non-cancer
   - Consequence: Model biased toward predicting "no cancer"
   - Solution: Use class weights or SMOTE (not done here)

2. **Missing Data:** Multiple missing values, especially in STD variables
   - Consequence: Reduced sample size for those variables
   - Solution: Multiple imputation (not done here; listwise deletion used)

3. **Single Hospital Sample:** Data from Hospital Universitario de Caracas (Venezuela)
   - Consequence: Results may not generalize to other populations
   - Interpretation: Use as proof-of-concept, not for clinical decision-making

4. **Cross-Sectional Data:** No temporal follow-up
   - Consequence: Cannot establish causality, only association
   - Interpretation: Associations are suggestive but not causal

---

## PART 5: NEXT STEPS & EXTENSIONS

### For Portfolio Enhancement

1. **Address Class Imbalance**
   - Use SMOTE (Synthetic Minority Over-Sampling) or undersampling
   - Refit models and compare performance

2. **Model Diagnostics**
   - Residual plots to check linearity assumption
   - VIF to check multicollinearity
   - Influential point detection (Cook's distance)

3. **Predictive Performance**
   - Train/test split; evaluate AUC, accuracy, precision, recall
   - Cross-validation (10-fold CV)

4. **Feature Importance**
   - Plot absolute standardized coefficients
   - Identify top 5 most predictive risk factors

5. **Geographic Link**
   - Correlate screening program status with incidence rates
   - Hypothesis: Countries with screening have lower incidence

### For Clinical Application (Future)

1. **Temporal Trends**
   - Estimate Annual Percentage Change (EAPC) in incidence
   - Forecast incidence to 2030

2. **Spatial Analysis**
   - Identify high-burden clusters using spatial statistics
   - Moran's I test for spatial autocorrelation

3. **Bayesian Inference**
   - Credible intervals for diagnostic test accuracy
   - Posterior predictive distributions

---

## REFERENCES FOR FURTHER READING

### Logistic Regression & Interpretation

- Hosmer Jr, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied logistic regression* (3rd ed.). John Wiley & Sons.

- Sperandei, S. (2014). Understanding logistic regression analysis. *Biochemia Medica*, 24(1), 12–18.

### Diagnostic Test Accuracy

- Greenberg, R. S., Daniels, S. R., Flanders, W. D., et al. (2001). *Medical epidemiology* (3rd ed.). Lange Medical Books/McGraw-Hill.

- Zhou, X. H., Obuchowski, N. A., & McClish, D. K. (2011). *Statistical methods in diagnostic medicine* (2nd ed.). John Wiley & Sons.

### Cervical Cancer Epidemiology

- Bray, F., Laversanne, M., Sung, H., et al. (2024). Global cancer statistics 2022: GLOBOCAN estimates. *CA Cancer J Clin*, 74(3), 229–263.

- Walboomers, J. M., Jacobs, M. V., Manos, M. M., et al. (1999). Human papillomavirus is a necessary cause of invasive cervical cancer worldwide. *J Pathol*, 189(1), 12–19.

---

**Document Version:** 1.0  
**Date:** May 2026  
**Author:** M (Portfolio Project)
