# CERVICAL CANCER ANALYSIS
# Load libraries

library(tidyverse)
library(ggplot2)
library(pROC)

# Set working directory
setwd("F:/GITHUB/Cervical_Cancer_Project")

# Load data
patients <- read.csv("data/cervical_cancer_patients.csv", sep = ";")
epidemiology <- read.csv("data/cervical_cancer_incidence_by_country_2022.csv")
screening <- read.csv("data/countries_with_cervical_cancer_screening_programs.csv")

# Convert "?" to NA
patients[patients == "?"] <- NA

# Check what you loaded
cat("Patients:", nrow(patients), "rows\n")
cat("Epidemiology:", nrow(epidemiology), "rows\n")
cat("Screening:", nrow(screening), "rows\n")
# ====== EXPLORATORY DATA ANALYSIS (EDA) ======

# Missing values summary
cat("\n--- MISSING VALUES ---\n")
missing_summary <- data.frame(
  Column = names(patients),
  Missing_Count = colSums(is.na(patients)),
  Missing_Percent = round(100 * colSums(is.na(patients)) / nrow(patients), 2)
)
print(missing_summary)

# Descriptive statistics
cat("\n--- DESCRIPTIVE STATISTICS ---\n")
numeric_cols <- sapply(patients, is.numeric)
print(summary(patients[, numeric_cols]))

# Check missing values in Age column
cat("Missing values in Age:", sum(is.na(patients$Age)), "\n")

# Check missing values in Sexual Partners column
cat("Missing values in Sexual Partners:", sum(is.na(patients$Number.of.sexual.partners)), "\n")

mean_age <- mean(patients$Age, na.rm = TRUE)
cat("Mean Age:", mean_age, "\n")

sd_age <- sd(patients$Age, na.rm = TRUE)
median_age <- median(patients$Age, na.rm = TRUE)
min_age <- min(patients$Age, na.rm = TRUE)
max_age <- max(patients$Age, na.rm = TRUE)

cat("Age Statistics:\n")
cat("  Mean:", mean_age, "\n")
cat("  SD:", sd_age, "\n")
cat("  Median:", median_age, "\n")
cat("  Min:", min_age, "\n")
cat("  Max:", max_age, "\n")
# Visualize age distribution
hist(patients$Age, 
     main = "Age Distribution", 
     xlab = "Age (years)", 
     ylab = "Frequency",
     col = "steelblue", 
     breaks = 20)

# ====== LOGISTIC REGRESSION (RISK MODELING) ======

cat("\n--- UNIVARIATE LOGISTIC REGRESSION ---\n")

# Build model: Does AGE predict cancer?
model_age <- glm(Dx.Cancer ~ Age, data = patients, family = "binomial")
summary(model_age)
# Extract odds ratio for Age
or_age <- exp(coef(model_age)["Age"])
cat("Odds Ratio for Age:", round(or_age, 3), "\n")

# 95% Confidence Interval
ci_age <- exp(confint(model_age)["Age", ])
cat("95% CI:", round(ci_age[1], 3), "–", round(ci_age[2], 3), "\n")

# Test multiple risk factors
risk_factors <- c("Age", "Number.of.sexual.partners", "Num.of.pregnancies", "Smokes", "STDs")

cat("\n--- Testing Multiple Risk Factors ---\n")

for (factor in risk_factors) {
  formula_str <- paste("Dx.Cancer ~", factor)
  model <- glm(as.formula(formula_str), data = patients, family = "binomial", na.action = na.omit)
  
  or <- exp(coef(model)[2])
  p_val <- summary(model)$coefficients[2, 4]
  
  cat(factor, "| OR:", round(or, 3), "| p-value:", round(p_val, 4), "\n")
}

model_final <- glm(Dx.Cancer ~ Age + STDs, data = patients, family = "binomial", na.action = na.omit)
summary(model_final)

model_final <- glm(Dx.Cancer ~ Age, data = patients, family = "binomial")
summary(model_final)

# ====== DIAGNOSTIC TEST ACCURACY ======

cat("\n--- ROC CURVES: Diagnostic Tests ---\n")

diagnostic_tests <- c("Hinselmann", "Schiller", "Citology", "Biopsy")

for (test in diagnostic_tests) {
  test_results <- patients[[test]]
  valid_idx <- !is.na(test_results) & !is.na(patients$Dx.Cancer)
  
  if (sum(valid_idx) > 0) {
    roc_obj <- roc(patients$Dx.Cancer[valid_idx], test_results[valid_idx], quiet = TRUE)
    cat(test, "| AUC:", round(roc_obj$auc, 3), "\n")
  }
}