suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(splines)
})

cfg <- source("R/config.R")$value
source("R/utils_survival.R")

df <- readRDS(file.path(cfg$processed_dir, "02_analytic_cohort_split.rds"))

train <- df |> filter(split == "train")
test  <- df |> filter(split == "test")

fit <- coxph(
  formula = cox_formula(),
  data = train,
  x = TRUE,
  y = TRUE
)

saveRDS(fit, file.path(cfg$processed_dir, "03_cox_model.rds"))

times <- cfg$horizons_months

train_pred <- bind_cols(
  train |> select(all_of(cfg$id_col), time_months, event, split),
  predict_survival_probs(fit, train, times)
)

test_pred <- bind_cols(
  test |> select(all_of(cfg$id_col), time_months, event, split),
  predict_survival_probs(fit, test, times)
)

saveRDS(train_pred, file.path(cfg$processed_dir, "03_train_predictions.rds"))
saveRDS(test_pred, file.path(cfg$processed_dir, "03_test_predictions.rds"))
