# Survival_QE_Ricardo

Supplementary code, data, and generated figures for the BIOSTATS 748 survival analysis qualifying exam project.

The analysis studies survival after discharge among simulated acute heart failure patients. The workflow moves from data cleaning and exploratory summaries to Kaplan-Meier analysis, Cox modeling, model diagnostics, alternative survival models, and gene-expression extensions.

## Contents

- `data/hf_survival-r.csv`: raw simulated cohort data.
- `output/dat_clean.rds`: cleaned analysis dataset used by the modeling R Markdown files.
- `Survival_QE_Data_Exploration.Rmd`: data checking, recoding, exploratory plots, and creation of `output/dat_clean.rds`.
- `Survival_QE_KM.Rmd`: Kaplan-Meier curves and log-rank comparison.
- `Survival_QE_Cox.Rmd`: Cox models, proportional hazards diagnostics, time-varying Cox models, and Royston-Parmar modeling.
- `Survival_QE_AFT.Rmd`: accelerated failure time model comparison.
- `Survival_QE_Gene.Rmd`: elastic-net Cox, PCA, supervised PCA, and internal validation for gene-expression predictors.
- `scripts/figure_survival_modeling.R`: script for the survival modeling composite figure.
- `scripts/goodness_of_fit_checks.R`: script for functional-form, PH, residual, and model-fit checks.
- `figures/`: generated report figures.
- `tables/`: generated model summaries and validation outputs.

## Suggested Run Order

Run the files from the repository root:

1. `Survival_QE_Data_Exploration.Rmd`
2. `Survival_QE_KM.Rmd`
3. `Survival_QE_Cox.Rmd`
4. `Survival_QE_AFT.Rmd`
5. `Survival_QE_Gene.Rmd`
6. `scripts/figure_survival_modeling.R`
7. `scripts/goodness_of_fit_checks.R`

The checked-in `figures/` and `tables/` folders contain the generated outputs used to support the final report.

## Notes

The gene-expression analysis is exploratory and internally validated only. The selected genes should not be interpreted as externally validated biomarkers.

