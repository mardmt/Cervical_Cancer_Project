# Cervical Cancer Epidemiology & Risk Stratification Analysis

## Overview

Epidemiological analysis of cervical cancer risk factors and population screening impact. Combines patient-level data (858 individuals), population incidence data (185 countries), and national screening program status to evaluate risk prediction and health disparities.

## Key Results

- Age predicts cancer risk: OR=1.061 per year (p=0.0018, 95% CI: 1.019–1.099)
- Countries with screening programs: 16.6 per 100,000 incidence
- Countries without screening: 21.77 per 100,000 (5.17 per 100k reduction)
- Highest burden: Eswatini (95.89/100k) vs lowest: Yemen (2.14/100k)
- Best diagnostic test: Schiller cytology (AUC=0.655)

## Methods

Logistic regression, ROC analysis, geographic disparities ranking. All analysis reproducible from source data via `scripts/analysis.R`.

## Data

- Patient factors: UCI ML Repository (858 patients, 36 variables)
- Incidence: Our World in Data / GLOBOCAN 2022 (185 countries)
- Screening programs: WHO / Our World in Data (776 countries)

## Code

Language: R (tidyverse, ggplot2, pROC)

## Files

- `scripts/analysis.R` — complete analysis
- `data/` — source datasets
- `results/` — outputs (figures, tables)
- `documentation/METHODOLOGY.md` — technical details

