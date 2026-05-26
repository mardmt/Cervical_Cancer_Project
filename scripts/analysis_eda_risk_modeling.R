# ==============================================================================
# CERVICAL CANCER EPIDEMIOLOGY & RISK STRATIFICATION ANALYSIS
# EDA + Risk Modeling (Focused Scope)
# ==============================================================================

library(DBI)
library(RSQLite)
library(tidyverse)
library(ggplot2)
library(pROC)

# Create results folder if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

cat("\n========== CERVICAL CANCER ANALYSIS: EDA + RISK MODELING ==========\n")

# ==============================================================================
# PHASE 1: LOAD DATA
# ==============================================================================

cat("\n--- LOADING DATA ---\n")

patients <- read.csv("data/cervical_cancer_patients.csv", sep = ";", stringsAsFactors = FALSE)
epidemiology <- read.csv("data/cervical_cancer_incidence_by_country_2022.csv", stringsAsFactors = FALSE)
screening <- read.csv("data/countries_with_cervical_cancer_screening_programs.csv", stringsAsFactors = FALSE)

cat("Patients loaded:", nrow(patients), "rows,", ncol(patients), "columns\n")
cat("Epidemiology loaded:", nrow(epidemiology), "rows\n")
cat("Screening loaded:", nrow(screening), "rows\n")

# Display column names
cat("\nPatient data columns:\n")
print(head(names(patients), 20))

# ==============================================================================
# PHASE 1: EXPLORATORY DATA ANALYSIS (EDA)
# ==============================================================================

cat("\n========== PHASE 1: EXPLORATORY DATA ANALYSIS ==========\n")

# --- 1.1: DATA QUALITY CHECK ---

cat("\n--- Missing Values Summary ---\n")

missing_summary <- data.frame(
  Column = names(patients),
  Missing_Count = colSums(is.na(patients)),
  Missing_Percent = round(100 * colSums(is.na(patients)) / nrow(patients), 2),
  row.names = NULL
)

print(missing_summary)
write.csv(missing_summary, "results/01_missing_values.csv", row.names = FALSE)

# --- 1.2: DESCRIPTIVE STATISTICS (CONTINUOUS VARIABLES) ---

cat("\n--- Descriptive Statistics (Numeric Variables) ---\n")

numeric_cols <- sapply(patients, is.numeric)
numeric_data <- patients[, numeric_cols]

summary_stats <- data.frame(
  Variable = names(numeric_data),
  N = colSums(!is.na(numeric_data)),
  Mean = round(colMeans(numeric_data, na.rm = TRUE), 2),
  SD = round(apply(numeric_data, 2, sd, na.rm = TRUE), 2),
  Median = round(apply(numeric_data, 2, median, na.rm = TRUE), 2),
  Min = round(apply(numeric_data, 2, min, na.rm = TRUE), 2),
  Max = round(apply(numeric_data, 2, max, na.rm = TRUE), 2),
  row.names = NULL
)

print(summary_stats)
write.csv(summary_stats, "results/02_descriptive_statistics.csv", row.names = FALSE)

# --- 1.3: EPIDEMIOLOGY SUMMARY ---

cat("\n--- Epidemiology Data Summary (2022) ---\n")

# Find incidence rate column
incidence_col <- grep("incidence|rate|per", colnames(epidemiology), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(incidence_col)) {
  epi_numeric <- as.numeric(epidemiology[[incidence_col]])
  
  epi_summary <- data.frame(
    Metric = c("N_Countries", "Mean_Incidence", "SD_Incidence", "Median_Incidence", 
               "Min_Incidence", "Max_Incidence", "Q1_Incidence", "Q3_Incidence"),
    Value = c(
      nrow(epidemiology),
      round(mean(epi_numeric, na.rm = TRUE), 2),
      round(sd(epi_numeric, na.rm = TRUE), 2),
      round(median(epi_numeric, na.rm = TRUE), 2),
      round(min(epi_numeric, na.rm = TRUE), 2),
      round(max(epi_numeric, na.rm = TRUE), 2),
      round(quantile(epi_numeric, 0.25, na.rm = TRUE), 2),
      round(quantile(epi_numeric, 0.75, na.rm = TRUE), 2)
    )
  )
  
  print(epi_summary)
  write.csv(epi_summary, "results/03_epidemiology_summary.csv", row.names = FALSE)
  
  # Top 10 highest burden countries
  top_10 <- epidemiology[order(epi_numeric, decreasing = TRUE), ][1:10, ]
  cat("\nTop 10 Countries by Cervical Cancer Incidence:\n")
  print(top_10[, c("Entity", incidence_col)])
  write.csv(top_10, "results/04_top_10_countries.csv", row.names = FALSE)
}

# --- 1.4: SCREENING PROGRAMS SUMMARY ---

cat("\n--- Screening Programs Distribution ---\n")

screening_col <- grep("screening", colnames(screening), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(screening_col)) {
  screening_summary <- data.frame(
    Status = c("Has Screening Program", "No Screening Program", "Unknown/Missing"),
    Count = c(
      sum(screening[[screening_col]] == "Yes", na.rm = TRUE),
      sum(screening[[screening_col]] == "No", na.rm = TRUE),
      sum(is.na(screening[[screening_col]]))
    )
  )
  
  screening_summary$Percent <- round(100 * screening_summary$Count / sum(screening_summary$Count), 1)
  
  print(screening_summary)
  write.csv(screening_summary, "results/05_screening_distribution.csv", row.names = FALSE)
}

# --- 1.5: CANCER DIAGNOSIS DISTRIBUTION ---

cat("\n--- Cancer Diagnosis Distribution ---\n")

cancer_col <- grep("cancer", colnames(patients), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(cancer_col)) {
  analysis_patients <- patients[!is.na(patients[[cancer_col]]), ]
  
  cancer_dist <- data.frame(
    Diagnosis = c("No Cancer", "Cancer"),
    Count = c(
      sum(analysis_patients[[cancer_col]] == 0 | analysis_patients[[cancer_col]] == "No", na.rm = TRUE),
      sum(analysis_patients[[cancer_col]] == 1 | analysis_patients[[cancer_col]] == "Yes", na.rm = TRUE)
    )
  )
  
  cancer_dist$Percent <- round(100 * cancer_dist$Count / sum(cancer_dist$Count), 1)
  
  print(cancer_dist)
  write.csv(cancer_dist, "results/06_cancer_diagnosis_distribution.csv", row.names = FALSE)
}

# --- 1.6: VISUALIZATIONS (EDA PLOTS) ---

cat("\n--- Generating EDA Plots (PDF) ---\n")

pdf("results/07_eda_plots.pdf", width = 14, height = 10)
par(mfrow = c(2, 2), mar = c(5, 4, 2, 1))

# Plot 1: Age distribution
if ("Age" %in% colnames(patients)) {
  age_data <- as.numeric(patients$Age)
  hist(age_data[!is.na(age_data)], 
       main = "Distribution of Patient Age", 
       xlab = "Age (years)", 
       col = "steelblue", 
       breaks = 20, 
       border = "white")
  grid(axis = "y", lty = 2, col = "gray80")
}

# Plot 2: Incidence rate distribution by country
if (!is.na(incidence_col)) {
  epi_numeric <- as.numeric(epidemiology[[incidence_col]])
  hist(epi_numeric[!is.na(epi_numeric)], 
       main = "Global Cervical Cancer Incidence Rates (2022)", 
       xlab = "Age-standardized rate per 100,000 women", 
       col = "coral", 
       breaks = 20, 
       border = "white")
  grid(axis = "y", lty = 2, col = "gray80")
}

# Plot 3: Top 10 countries by incidence
if (!is.na(incidence_col)) {
  top_10_plot <- epidemiology[order(epi_numeric, decreasing = TRUE), ][1:10, ]
  top_10_rates <- as.numeric(top_10_plot[[incidence_col]])
  top_10_names <- top_10_plot$Entity
  
  barplot(top_10_rates, 
          names.arg = top_10_names,
          main = "Top 10 Countries: Cervical Cancer Incidence (2022)", 
          ylab = "Rate per 100,000 women", 
          col = "lightgreen", 
          las = 2, 
          cex.names = 0.8, 
          border = "white")
  grid(axis = "y", lty = 2, col = "gray80")
}

# Plot 4: Cancer diagnosis proportions
if (!is.na(cancer_col)) {
  analysis_patients <- patients[!is.na(patients[[cancer_col]]), ]
  cancer_counts <- c(
    sum(analysis_patients[[cancer_col]] == 0 | analysis_patients[[cancer_col]] == "No", na.rm = TRUE),
    sum(analysis_patients[[cancer_col]] == 1 | analysis_patients[[cancer_col]] == "Yes", na.rm = TRUE)
  )
  
  pie(cancer_counts, 
      labels = c("No Cancer", "Cancer"),
      main = "Cancer Diagnosis Distribution",
      col = c("lightblue", "salmon"),
      border = "white",
      cex = 1.2)
}

par(mfrow = c(1, 1))
dev.off()

cat("EDA plots saved to: results/07_eda_plots.pdf\n")

# ==============================================================================
# PHASE 2: RISK STRATIFICATION (LOGISTIC REGRESSION)
# ==============================================================================

cat("\n========== PHASE 2: RISK STRATIFICATION ==========\n")

# Prepare data for modeling
cat("\n--- Preparing Data for Modeling ---\n")

# Create clean analysis dataset
model_data <- patients

# Identify cancer diagnosis column and convert to binary
cancer_col <- grep("cancer", colnames(model_data), ignore.case = TRUE, value = TRUE)[1]

if (!is.na(cancer_col)) {
  # Remove rows with missing cancer diagnosis
  model_data <- model_data[!is.na(model_data[[cancer_col]]), ]
  
  # Convert to binary (0/1)
  model_data$cancer_binary <- ifelse(
    model_data[[cancer_col]] == 1 | model_data[[cancer_col]] == "Yes", 
    1, 
    0
  )
  
  cat("Analysis dataset:", nrow(model_data), "patients\n")
  cat("Cancer cases:", sum(model_data$cancer_binary, na.rm = TRUE), "\n")
  cat("Non-cancer:", nrow(model_data) - sum(model_data$cancer_binary, na.rm = TRUE), "\n")
  
  # --- 2.1: UNIVARIATE LOGISTIC REGRESSION ---
  
  cat("\n--- Univariate Logistic Regression ---\n")
  
  # Select numeric columns for univariate analysis (exclude cancer diagnosis and diagnostic tests)
  exclude_cols <- c("Dx.Cancer", "Dx.CIN", "Dx.HPV", "Dx", "Hinselmann", "Schiller", "Citology", "Biopsy", "cancer_binary")
  
  numeric_predictors <- names(model_data)[
    sapply(model_data, is.numeric) & 
    !(names(model_data) %in% exclude_cols)
  ]
  
  univariate_results <- data.frame()
  
  for (pred in numeric_predictors) {
    pred_data <- model_data[[pred]]
    
    # Skip if all missing
    if (all(is.na(pred_data))) next
    
    # Skip if no variation
    if (sd(pred_data, na.rm = TRUE) == 0 | is.na(sd(pred_data, na.rm = TRUE))) next
    
    # Fit univariate model
    try({
      fit <- glm(cancer_binary ~ pred_data, 
                 data = model_data, 
                 family = "binomial", 
                 na.action = na.omit)
      
      summary_fit <- summary(fit)
      coef_table <- summary_fit$coefficients
      
      if (nrow(coef_table) > 1) {
        est <- coef_table[2, 1]
        se <- coef_table[2, 2]
        pval <- coef_table[2, 4]
        or <- exp(est)
        or_ci_lower <- exp(est - 1.96 * se)
        or_ci_upper <- exp(est + 1.96 * se)
        
        univariate_results <- rbind(univariate_results, data.frame(
          Variable = pred,
          OR = round(or, 3),
          CI_Lower = round(or_ci_lower, 3),
          CI_Upper = round(or_ci_upper, 3),
          P_Value = round(pval, 4),
          Significant = ifelse(pval < 0.05, "Yes", "No")
        ))
      }
    }, silent = TRUE)
  }
  
  print(univariate_results)
  write.csv(univariate_results, "results/08_univariate_logistic_regression.csv", row.names = FALSE)
  
  # --- 2.2: MULTIVARIATE LOGISTIC REGRESSION (STEPWISE SELECTION) ---
  
  cat("\n--- Multivariate Logistic Regression (Stepwise AIC) ---\n")
  
  # Select top significant variables from univariate analysis
  if (nrow(univariate_results) > 0) {
    sig_vars <- univariate_results$Variable[univariate_results$Significant == "Yes"]
    
    if (length(sig_vars) > 0) {
      # Build formula for multivariate model
      formula_str <- paste("cancer_binary ~", paste(sig_vars, collapse = " + "))
      
      # Fit initial multivariate model
      try({
        initial_fit <- glm(as.formula(formula_str), 
                          data = model_data, 
                          family = "binomial", 
                          na.action = na.omit)
        
        # Stepwise selection using AIC
        step_fit <- step(initial_fit, direction = "both", trace = 0)
        
        cat("\nFinal Multivariate Model:\n")
        print(summary(step_fit))
        
        # Extract model coefficients
        multivar_coef <- data.frame(
          Variable = names(coef(step_fit))[-1],
          Coefficient = round(coef(step_fit)[-1], 4),
          OR = round(exp(coef(step_fit)[-1]), 3),
          P_Value = round(summary(step_fit)$coefficients[-1, 4], 4)
        )
        
        write.csv(multivar_coef, "results/09_multivariate_logistic_regression.csv", row.names = FALSE)
        
        # Model comparison
        cat("\nModel Comparison (AIC):\n")
        cat("Initial AIC:", AIC(initial_fit), "\n")
        cat("Final AIC:", AIC(step_fit), "\n")
      }, silent = FALSE)
    } else {
      cat("No significant variables found for multivariate model.\n")
    }
  }
  
  # --- 2.3: DIAGNOSTIC TEST ACCURACY (SENSITIVITY, SPECIFICITY, ROC) ---
  
  cat("\n--- Diagnostic Test Accuracy Metrics ---\n")
  
  diagnostic_tests <- c("Hinselmann", "Schiller", "Citology", "Biopsy")
  diagnostic_results <- data.frame()
  
  pdf("results/10_diagnostic_roc_curves.pdf", width = 12, height = 8)
  par(mfrow = c(2, 2), mar = c(5, 4, 2, 1))
  
  for (test in diagnostic_tests) {
    if (test %in% colnames(model_data)) {
      test_col <- model_data[[test]]
      
      # Remove missing values
      valid_idx <- !is.na(test_col) & !is.na(model_data$cancer_binary)
      
      if (sum(valid_idx) > 0) {
        test_results <- as.numeric(test_col[valid_idx])
        cancer_status <- model_data$cancer_binary[valid_idx]
        
        # Calculate sensitivity and specificity
        sensitivity <- sum(test_results == 1 & cancer_status == 1) / sum(cancer_status == 1)
        specificity <- sum(test_results == 0 & cancer_status == 0) / sum(cancer_status == 0)
        ppv <- sum(test_results == 1 & cancer_status == 1) / sum(test_results == 1)
        npv <- sum(test_results == 0 & cancer_status == 0) / sum(test_results == 0)
        
        # Calculate ROC curve and AUC
        roc_obj <- roc(cancer_status, test_results, quiet = TRUE)
        
        diagnostic_results <- rbind(diagnostic_results, data.frame(
          Test = test,
          Sensitivity = round(sensitivity, 3),
          Specificity = round(specificity, 3),
          PPV = round(ppv, 3),
          NPV = round(npv, 3),
          AUC = round(roc_obj$auc, 3),
          N = sum(valid_idx)
        ))
        
        # Plot ROC curve
        plot(roc_obj, 
             main = paste("ROC Curve:", test),
             col = "steelblue", 
             lwd = 2,
             print.auc = TRUE,
             print.auc.x = 0.5,
             print.auc.y = 0.3)
      }
    }
  }
  
  par(mfrow = c(1, 1))
  dev.off()
  
  print(diagnostic_results)
  write.csv(diagnostic_results, "results/11_diagnostic_test_accuracy.csv", row.names = FALSE)
  
  cat("\nDiagnostic ROC curves saved to: results/10_diagnostic_roc_curves.pdf\n")
  
} else {
  cat("Cancer diagnosis column not found. Risk modeling skipped.\n")
}

# ==============================================================================
# ANALYSIS COMPLETE
# ==============================================================================

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("\nOutput files saved to results/ folder:\n")
cat("  01_missing_values.csv\n")
cat("  02_descriptive_statistics.csv\n")
cat("  03_epidemiology_summary.csv\n")
cat("  04_top_10_countries.csv\n")
cat("  05_screening_distribution.csv\n")
cat("  06_cancer_diagnosis_distribution.csv\n")
cat("  07_eda_plots.pdf\n")
cat("  08_univariate_logistic_regression.csv\n")
cat("  09_multivariate_logistic_regression.csv\n")
cat("  10_diagnostic_roc_curves.pdf\n")
cat("  11_diagnostic_test_accuracy.csv\n")
cat("\n✓ Ready for GitHub push\n")
