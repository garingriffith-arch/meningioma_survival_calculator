suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(splines)
  library(tibble)
})

set.seed(1)

cfg <- source("R/config.R")$value
source("R/utils_survival.R")

set.seed(cfg$bootstrap_seed)

df  <- readRDS(file.path(cfg$processed_dir, "01_analytic_cohort.rds"))
fit <- readRDS(file.path(cfg$processed_dir, "03_cox_model.rds"))

times <- cfg$horizons_months
B <- cfg$bootstrap_B
n <- nrow(df)

# Apparent horizon-specific IPCW metrics on the full analytic cohort
app_metrics <- horizon_metrics(fit, df, times)

# Apparent global metrics
lp_full <- as.numeric(stats::predict(fit, newdata = df, type = "lp"))
app_c_global <- harrell_c(df$time_months, df$event, lp_full)
app_slope_global <- calibration_slope(df$time_months, df$event, lp_full)

# Calibration-by-group files
cal_dir <- file.path(cfg$processed_dir, "calibration")
dir.create(cal_dir, recursive = TRUE, showWarnings = FALSE)

calibration <- lapply(times, function(t) {
  calibration_by_group(fit, df, t, groups = cfg$calibration_groups)
})
names(calibration) <- paste0(times, "m")

for (nm in names(calibration)) {
  write_csv(calibration[[nm]], file.path(cal_dir, paste0("calibration_", nm, ".csv")))
}

# Bootstrap optimism correction at each horizon
optim_auc   <- matrix(NA_real_, nrow = B, ncol = length(times))
optim_brier <- matrix(NA_real_, nrow = B, ncol = length(times))
colnames(optim_auc) <- paste0(times, "m")
colnames(optim_brier) <- paste0(times, "m")

optim_c_global <- numeric(B)
optim_s_global <- numeric(B)

for (b in seq_len(B)) {
  idx <- sample.int(n, size = n, replace = TRUE)
  boot_df <- df[idx, , drop = FALSE]
  
  fit_b <- try(
    survival::coxph(
      formula = cox_formula(cfg),
      data = boot_df,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    silent = TRUE
  )
  
  if (!inherits(fit_b, "try-error")) {
    met_in  <- horizon_metrics(fit_b, boot_df, times)
    met_out <- horizon_metrics(fit_b, df, times)
    
    optim_auc[b, ]   <- met_in$auc   - met_out$auc
    optim_brier[b, ] <- met_in$brier - met_out$brier
    
    lp_in  <- as.numeric(stats::predict(fit_b, newdata = boot_df, type = "lp"))
    lp_out <- as.numeric(stats::predict(fit_b, newdata = df, type = "lp"))
    
    c_in  <- harrell_c(boot_df$time_months, boot_df$event, lp_in)
    c_out <- harrell_c(df$time_months, df$event, lp_out)
    
    s_in  <- calibration_slope(boot_df$time_months, boot_df$event, lp_in)
    s_out <- calibration_slope(df$time_months, df$event, lp_out)
    
    optim_c_global[b] <- c_in - c_out
    optim_s_global[b] <- s_in - s_out
  } else {
    optim_c_global[b] <- NA_real_
    optim_s_global[b] <- NA_real_
  }
  
  if (b %% 50 == 0) {
    cat("Bootstrap", b, "of", B, "\n")
  }
}

mean_auc_optimism   <- colMeans(optim_auc, na.rm = TRUE)
mean_brier_optimism <- colMeans(optim_brier, na.rm = TRUE)

horizon_tbl <- app_metrics |>
  dplyr::mutate(
    auc_apparent = auc,
    auc_optimism = unname(mean_auc_optimism[paste0(horizon_months, "m")]),
    auc_corrected = auc_apparent - auc_optimism,
    brier_apparent = brier,
    brier_optimism = unname(mean_brier_optimism[paste0(horizon_months, "m")]),
    brier_corrected = brier_apparent - brier_optimism
  ) |>
  dplyr::select(
    horizon_months,
    sample_size,
    deaths_by_horizon,
    auc_apparent,
    auc_optimism,
    auc_corrected,
    brier_apparent,
    brier_optimism,
    brier_corrected
  )

results <- list(
  apparent = list(
    harrell_c = app_c_global,
    calibration_slope = app_slope_global,
    horizon_metrics = app_metrics
  ),
  calibration = calibration,
  bootstrap = list(
    B = B,
    optimism = list(
      auc = optim_auc,
      brier = optim_brier,
      c_global = optim_c_global,
      slope_global = optim_s_global,
      mean_auc = mean_auc_optimism,
      mean_brier = mean_brier_optimism,
      mean_c_global = mean(optim_c_global, na.rm = TRUE),
      mean_slope_global = mean(optim_s_global, na.rm = TRUE)
    ),
    corrected = list(
      horizon_metrics = horizon_tbl,
      harrell_c = app_c_global - mean(optim_c_global, na.rm = TRUE),
      calibration_slope = app_slope_global - mean(optim_s_global, na.rm = TRUE)
    )
  )
)

saveRDS(results, file.path(cfg$processed_dir, "04_validation_results.rds"))
write_csv(horizon_tbl, file.path(cfg$processed_dir, "04_validation_summary.csv"))
