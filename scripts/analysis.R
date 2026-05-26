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

# ====== PUBLICATION-READY VISUALIZATIONS ======

# Create a multi-panel PDF with 4 plots

pdf("results/cervical_cancer_analysis_figures.pdf", width = 14, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

# PLOT 1: Age distribution
hist(patients$Age, 
     main = "A. Age Distribution of Patients (n=858)",
     xlab = "Age (years)", 
     ylab = "Frequency",
     col = "#2E86AB", 
     breaks = 20,
     border = "white",
     cex.main = 1.3,
     cex.lab = 1.1)
grid(NA, NA, lty = 2, col = "gray80")

# PLOT 2: Cancer diagnosis
cancer_counts <- table(patients$Dx.Cancer)
barplot(cancer_counts,
        main = "B. Cervical Cancer Diagnosis",
        names.arg = c("No Cancer", "Cancer"),
        ylab = "Number of Patients",
        col = c("#A23B72", "#F18F01"),
        border = "white",
        cex.main = 1.3,
        cex.lab = 1.1)


# PLOT 3: Top 10 countries by incidence
incidence_col <- grep("incidence", colnames(epidemiology), ignore.case = TRUE, value = TRUE)[1]
epi_numeric <- as.numeric(epidemiology[[incidence_col]])
top_10 <- epidemiology[order(epi_numeric, decreasing = TRUE), ][1:10, ]
top_10_rates <- as.numeric(top_10[[incidence_col]])
top_10_names <- top_10$Entity

barplot(top_10_rates,
        names.arg = top_10_names,
        main = "C. Top 10 Countries: Cervical Cancer Incidence (2022)",
        ylab = "Age-standardized rate per 100,000 women",
        col = "#06A77D",
        las = 2,
        cex.names = 0.85,
        cex.main = 1.3,
        cex.lab = 1.1,
        border = "white")


# PLOT 4: Global incidence distribution
hist(epi_numeric[!is.na(epi_numeric)],
     main = "D. Global Cervical Cancer Incidence Distribution",
     xlab = "Age-standardized rate per 100,000 women",
     ylab = "Number of Countries",
     col = "#E63946",
     breaks = 15,
     border = "white",
     cex.main = 1.3,
     cex.lab = 1.1)


par(mfrow = c(1, 1))
dev.off()

cat("\n✓ Publication-ready figures saved to: results/cervical_cancer_analysis_figures.pdf\n")

# ====== GEOGRAPHIC DISPARITIES ======

cat("\n--- TOP 10 COUNTRIES: HIGHEST BURDEN ---\n")

incidence_col <- grep("incidence", colnames(epidemiology), ignore.case = TRUE, value = TRUE)[1]
epi_numeric <- as.numeric(epidemiology[[incidence_col]])

top_10 <- epidemiology[order(epi_numeric, decreasing = TRUE), ][1:10, c("Entity", incidence_col)]
print(top_10)
write.csv(top_10, "results/top_10_highest_burden.csv", row.names = FALSE)

cat("\n--- BOTTOM 10 COUNTRIES: LOWEST BURDEN ---\n")

bottom_10 <- epidemiology[order(epi_numeric, decreasing = FALSE), ][1:10, c("Entity", incidence_col)]
print(bottom_10)
write.csv(bottom_10, "results/bottom_10_lowest_burden.csv", row.names = FALSE)

cat("\n--- CORRELATION: Screening vs Incidence ---\n")

# Merge screening data with epidemiology
screening_recent <- screening[screening$Year == max(screening$Year), ]
merged <- merge(epidemiology, screening_recent, by = "Entity", all = FALSE)

# Count screening programs
screening_yes <- sum(merged[, grep("screening", colnames(merged), ignore.case = TRUE)] == "Yes", na.rm = TRUE)
screening_no <- sum(merged[, grep("screening", colnames(merged), ignore.case = TRUE)] == "No", na.rm = TRUE)

cat("Countries WITH screening programs:", screening_yes, "\n")
cat("Countries WITHOUT screening programs:", screening_no, "\n")

# Correlation
incidence_col_epi <- grep("incidence", colnames(merged), ignore.case = TRUE, value = TRUE)[1]
screening_col <- grep("screening", colnames(merged), ignore.case = TRUE, value = TRUE)[1]

# Compare mean incidence
merged_numeric <- as.numeric(merged[[incidence_col_epi]])
screening_status <- merged[[screening_col]]

mean_with <- mean(merged_numeric[screening_status == "Yes"], na.rm = TRUE)
mean_without <- mean(merged_numeric[screening_status == "No"], na.rm = TRUE)

cat("\nMean incidence WITH screening:", round(mean_with, 2), "\n")
cat("Mean incidence WITHOUT screening:", round(mean_without, 2), "\n")
cat("Difference:", round(mean_without - mean_with, 2), "per 100,000\n")

# Visualization: Screening vs Incidence
pdf("results/screening_vs_incidence.pdf", width = 10, height = 6)

# Boxplot comparison
boxplot(merged_numeric ~ screening_status,
        main = "Cervical Cancer Incidence: Impact of National Screening Programs",
        names = c("No Screening", "Screening Program"),
        ylab = "Age-standardized incidence rate per 100,000 women",
        col = c("#E63946", "#06A77D"),
        border = "black",
        cex.main = 1.3,
        cex.lab = 1.1)

# Add mean points
points(c(1, 2), c(mean_without, mean_with), 
       pch = 18, cex = 3, col = "black")

dev.off()

cat("\n✓ Screening impact visualization saved\n")


# ====== SAVE ALL RESULTS TO CSV ======

cat("\n--- EXPORTING SUMMARY TABLES ---\n")

# 1. Missing values summary
write.csv(missing_summary, "results/missing_values_summary.csv", row.names = FALSE)
cat("✓ missing_values_summary.csv\n")

# 2. Descriptive statistics
desc_stats <- data.frame(
  Statistic = c("Mean", "SD", "Median", "Min", "Max"),
  Age = c(mean_age, sd_age, median_age, min_age, max_age)
)
write.csv(desc_stats, "results/age_descriptive_statistics.csv", row.names = FALSE)
cat("✓ age_descriptive_statistics.csv\n")

# 3. Univariate logistic regression results
risk_results <- data.frame(
  Risk_Factor = risk_factors,
  OR = NA,
  P_Value = NA
)
for (i in 1:nrow(risk_results)) {
  factor <- risk_factors[i]
  formula_str <- paste("Dx.Cancer ~", factor)
  model <- glm(as.formula(formula_str), data = patients, family = "binomial", na.action = na.omit)
  risk_results$OR[i] <- round(exp(coef(model)[2]), 3)
  risk_results$P_Value[i] <- round(summary(model)$coefficients[2, 4], 4)
}
write.csv(risk_results, "results/univariate_logistic_regression.csv", row.names = FALSE)
cat("✓ univariate_logistic_regression.csv\n")

# 4. Final model summary
final_model_summary <- data.frame(
  Model = "Age Only (Logistic Regression)",
  Coefficient = round(coef(model_final)["Age"], 4),
  Odds_Ratio = round(exp(coef(model_final)["Age"]), 3),
  P_Value = round(summary(model_final)$coefficients["Age", 4], 4),
  AIC = round(model_final$aic, 2),
  N_Observations = nrow(patients)
)
write.csv(final_model_summary, "results/final_model_summary.csv", row.names = FALSE)
cat("✓ final_model_summary.csv\n")

# 5. Diagnostic test accuracy (ROC AUC)
diagnostic_results <- data.frame(
  Test = diagnostic_tests,
  AUC = NA
)
for (i in 1:nrow(diagnostic_results)) {
  test <- diagnostic_tests[i]
  test_results <- patients[[test]]
  valid_idx <- !is.na(test_results) & !is.na(patients$Dx.Cancer)
  if (sum(valid_idx) > 0) {
    roc_obj <- roc(patients$Dx.Cancer[valid_idx], test_results[valid_idx], quiet = TRUE)
    diagnostic_results$AUC[i] <- round(roc_obj$auc, 3)
  }
}
write.csv(diagnostic_results, "results/diagnostic_test_accuracy.csv", row.names = FALSE)
cat("✓ diagnostic_test_accuracy.csv\n")

# 6. Screening impact summary
screening_summary <- data.frame(
  Group = c("WITH Screening", "WITHOUT Screening"),
  Mean_Incidence = c(round(mean_with, 2), round(mean_without, 2)),
  Difference = c(NA, round(mean_without - mean_with, 2)),
  Percent_Reduction = c(NA, round(100 * (mean_without - mean_with) / mean_without, 1))
)
write.csv(screening_summary, "results/screening_impact_summary.csv", row.names = FALSE)
cat("✓ screening_impact_summary.csv\n")

# 7. Geographic burden rankings
geographic_summary <- data.frame(
  Category = c("Highest Burden (Top)", "Lowest Burden (Bottom)"),
  Country = c(top_10$Entity[1], bottom_10$Entity[1]),
  Rate_Per_100k = c(
    round(as.numeric(top_10[[incidence_col]][1]), 2),
    round(as.numeric(bottom_10[[incidence_col]][1]), 2)
  )
)
write.csv(geographic_summary, "results/geographic_burden_rankings.csv", row.names = FALSE)
cat("✓ geographic_burden_rankings.csv\n")

# ====== TEMPORAL TRENDS ANALYSIS ======

cat("\n--- TEMPORAL TRENDS: 2015-2021 ---\n")

# Check if epidemiology has Year column
if ("Year" %in% colnames(epidemiology)) {
  
  # Filter for years 2015-2021
  epi_trends <- epidemiology[epidemiology$Year >= 2015 & epidemiology$Year <= 2021, ]
  
  # Calculate mean incidence by year (global)
  incidence_col <- grep("incidence", colnames(epi_trends), ignore.case = TRUE, value = TRUE)[1]
  year_trends <- aggregate(
    as.numeric(epi_trends[[incidence_col]]) ~ Year,
    data = epi_trends,
    FUN = mean,
    na.rm = TRUE
  )
  colnames(year_trends) <- c("Year", "Mean_Incidence")
  
  # Calculate EAPC (Estimated Annual Percentage Change)
  if (nrow(year_trends) > 1) {
    years <- year_trends$Year
    incidence <- year_trends$Mean_Incidence
    
    # Linear regression on log scale
    log_incidence <- log(incidence)
    slope <- coef(lm(log_incidence ~ years))[2]
    eapc <- (exp(slope) - 1) * 100
    
    cat("Years available:", min(years), "-", max(years), "\n")
    cat("EAPC:", round(eapc, 2), "%\n")
    
    # Save trends
    write.csv(year_trends, "results/temporal_trends_2015_2021.csv", row.names = FALSE)
    cat("✓ temporal_trends_2015_2021.csv\n")
  }
  
} else {
  cat("⚠ No Year column in epidemiology data.\n")
  cat("Need historical data from Our World in Data (2015-2021).\n")
  cat("Download from: https://ourworldindata.org/grapher/rate-of-new-cervical-cancer-cases-gco\n")
}

cat("\n✓✓✓ ALL RESULTS EXPORTED ✓✓✓\n")