suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(splines)
})

cfg <- source("R/config.R")$value
source("R/utils_survival.R")

df <- readRDS(file.path(cfg$processed_dir, "01_analytic_cohort.rds"))

fit <- coxph(
  formula = cox_formula(cfg),
  data = df,
  x = TRUE,
  y = TRUE,
  model = TRUE
)

saveRDS(df,  file.path(cfg$processed_dir, "02_model_data.rds"))
saveRDS(fit, file.path(cfg$processed_dir, "03_cox_model.rds"))

pred <- bind_cols(
  df |> dplyr::select(all_of(cfg$id_col), time_months, event),
  predict_survival_probs(fit, df, cfg$horizons_months)
)

saveRDS(pred, file.path(cfg$processed_dir, "03_full_predictions.rds"))
