# Grade II–III Intracranial Meningioma Survival Estimator

This repository contains code used to develop and internally validate a diagnosis-time survival prediction model for adults with WHO grade II–III intracranial meningioma using National Cancer Database (NCDB) data, along with a Shiny application for survival estimation.

## Project Overview

- **Study design:** Retrospective population-based cohort study
- **Data source:** National Cancer Database (NCDB)
- **Population:** Adults with grade II–III intracranial meningioma
- **Outcome:** Overall survival
- **Model:** Multivariable Cox proportional hazards regression
- **Validation:** Temporal holdout testing and bootstrap internal validation with optimism correction
- **Prediction horizons:** 1-, 3-, and 5-year overall survival

The model uses variables available at the time of diagnosis and is intended to balance interpretability, stability, and clinical usability.

## Predictors Included

- Age
- Sex
- Charlson–Deyo comorbidity score
- WHO grade (II vs III)
- Tumor size (continuous, mm)
- Race
- Ethnicity

Age and tumor size are modeled with restricted cubic spline terms to allow for nonlinear effects.

## Repository Structure

```text
R/
  config.R
  utils_recode.R
  utils_survival.R
  01_build_cohort.R
  02_split_data.R
  03_fit_model.R
  04_validate_model.R
  05_make_tables_figures.R

app/
  app.R

data/
  raw/
  processed/

tables/
figures/
