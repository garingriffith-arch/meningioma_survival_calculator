suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(splines)
})

set.seed(1)

cfg <- source("R/config.R")$value
source("R/utils_survival.R")

fit       <- readRDS(file.path(cfg$processed_dir, "03_cox_model.rds"))
test_pred <- readRDS(file.path(cfg$processed_dir, "03_test_predictions.rds"))
split_df  <- readRDS(file.path(cfg$processed_dir, "02_analytic_cohort_split.rds"))
full_df   <- readRDS(file.path(cfg$processed_dir, "01_analytic_cohort.rds"))

times <- cfg$horizons_months

concordance_at_time <- function(df, t, risk_col) {
  tt <- pmin(df$time_months, t)
  status_t <- as.integer(df$event == 1 & df$time_months <= t)
  cc <- survival::concordance(
    Surv(tt, status_t) ~ df[[risk_col]],
    data = df,
    reverse = TRUE
  )
  as.numeric(cc$concordance)
}

# Test-set discrimination
test_df <- test_pred
test_newdata <- split_df |> filter(split == "test")

test_df$lp <- as.numeric(
  predict(fit, newdata = test_newdata, type = "lp")
)

test_harrell_c <- harrell_c(test_df$time_months, test_df$event, test_df$lp)

# Horizon-specific concordance to preserve original manuscript outputs
for (t in times) {
  test_df[[paste0("risk_", t, "m")]] <- 1 - test_df[[paste0("S_", t, "m")]]
}

c_horizon <- sapply(times, function(t) {
  concordance_at_time(test_df, t, paste0("risk_", t, "m"))
})
names(c_horizon) <- paste0("C@_", times, "m")

# Calibration at prespecified horizons
calibration <- lapply(times, function(t) {
  calibration_by_group(test_df, t, paste0("S_", t, "m"), groups = 10)
})
names(calibration) <- paste0(times, "m")

cal_dir <- file.path(cfg$processed_dir, "calibration")
dir.create(cal_dir, recursive = TRUE, showWarnings = FALSE)

for (nm in names(calibration)) {
  write_csv(calibration[[nm]], file.path(cal_dir, paste0("calibration_", nm, ".csv")))
}

# Bootstrap optimism on full analytic cohort
B <- 1000
n <- nrow(full_df)

optim_c <- numeric(B)
optim_s <- numeric(B)

full_fit <- coxph(cox_formula(), data = full_df, x = TRUE, y = TRUE)
full_lp  <- as.numeric(predict(full_fit, newdata = full_df, type = "lp"))

app_c <- harrell_c(full_df$time_months, full_df$event, full_lp)
app_s <- calibration_slope(full_df$time_months, full_df$event, full_lp)

for (b in seq_len(B)) {
  idx <- sample.int(n, size = n, replace = TRUE)
  boot_df <- full_df[idx, , drop = FALSE]
  
  fit_b <- try(coxph(cox_formula(), data = boot_df, x = TRUE, y = TRUE), silent = TRUE)
  
  if (inherits(fit_b, "try-error")) {
    optim_c[b] <- NA_real_
    optim_s[b] <- NA_real_
  } else {
    lp_in  <- as.numeric(predict(fit_b, newdata = boot_df, type = "lp"))
    lp_out <- as.numeric(predict(fit_b, newdata = full_df, type = "lp"))
    
    c_in  <- harrell_c(boot_df$time_months, boot_df$event, lp_in)
    c_out <- harrell_c(full_df$time_months, full_df$event, lp_out)
    
    s_in  <- calibration_slope(boot_df$time_months, boot_df$event, lp_in)
    s_out <- calibration_slope(full_df$time_months, full_df$event, lp_out)
    
    optim_c[b] <- c_in - c_out
    optim_s[b] <- s_in - s_out
  }
  
  if (b %% 50 == 0) {
    cat("Bootstrap", b, "of", B, "\n")
  }
}

results <- list(
  test = list(
    n = nrow(test_df),
    events = sum(test_df$event == 1, na.rm = TRUE),
    harrell_c = test_harrell_c
  ),
  c_horizon = c_horizon,
  calibration = calibration,
  bootstrap = list(
    B = B,
    apparent = list(C = app_c, slope = app_s),
    optimism = list(
      C = optim_c,
      slope = optim_s,
      mean_C = mean(optim_c, na.rm = TRUE),
      mean_slope = mean(optim_s, na.rm = TRUE)
    ),
    corrected = list(
      C = app_c - mean(optim_c, na.rm = TRUE),
      slope = app_s - mean(optim_s, na.rm = TRUE)
    )
  )
)

saveRDS(results, file.path(cfg$processed_dir, "04_validation_results.rds"))

summary_tbl <- data.frame(
  model = "Cox PH",
  test_harrell_c = results$test$harrell_c,
  C_12m = unname(results$c_horizon["C@_12m"]),
  C_36m = unname(results$c_horizon["C@_36m"]),
  C_60m = unname(results$c_horizon["C@_60m"]),
  apparent_C = results$bootstrap$apparent$C,
  optimism_C = results$bootstrap$optimism$mean_C,
  corrected_C = results$bootstrap$corrected$C,
  apparent_slope = results$bootstrap$apparent$slope,
  optimism_slope = results$bootstrap$optimism$mean_slope,
  corrected_slope = results$bootstrap$corrected$slope
)

write_csv(summary_tbl, file.path(cfg$processed_dir, "04_validation_summary.csv"))
