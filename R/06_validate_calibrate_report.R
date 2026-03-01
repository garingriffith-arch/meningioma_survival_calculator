# R/06_validate_calibrate_report.R
# Internal validation, calibration, and diagnostics for the final Cox model (05)
# Outputs: optimism-corrected C-index, calibration plots (1/3/5y), PH checks, tables

library(dplyr)
library(survival)
library(ggplot2)

set.seed(20260227)

# ----------------------------
# load model + analysis data
# ----------------------------
fit <- readRDS("model/cox_simple_fit.rds")
df  <- readRDS("data/processed/analysis_df.rds")

# recreate analysis frame used in 05 (must match exactly)
analysis_df <- df %>%
  mutate(
    time_months = as.numeric(time_months),
    event = as.integer(event),

    AGE = as.numeric(AGE),
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    GRADE = factor(GRADE, levels = c(2, 3), labels = c("II", "III")),
    CDCC_TOTAL_BEST = as.numeric(CDCC_TOTAL_BEST),
    YEAR_OF_DIAGNOSIS = as.numeric(YEAR_OF_DIAGNOSIS),

    # tumor size (collapsed, same logic as 05)
    TUMOR_SIZE_MM = ifelse(TUMOR_SIZE >= 999, NA, TUMOR_SIZE),
    TUMOR_SIZE_CAT = case_when(
      is.na(TUMOR_SIZE_MM) ~ "Unknown",
      TUMOR_SIZE_MM <= 30  ~ "<=3cm",
      TUMOR_SIZE_MM <= 60  ~ "3–6cm",
      TUMOR_SIZE_MM >  60  ~ ">6cm"
    ),
    TUMOR_SIZE_CAT = factor(
      TUMOR_SIZE_CAT,
      levels = c("<=3cm", "3–6cm", ">6cm", "Unknown")
    ),

    # race (collapsed)
    RACE_CAT = case_when(
      RACE == 1 ~ "White",
      RACE == 2 ~ "Black",
      TRUE      ~ "Other/Unknown"
    ),
    RACE_CAT = factor(RACE_CAT)
  ) %>%
  filter(
    !is.na(time_months),
    !is.na(event),
    !is.na(AGE),
    !is.na(SEX),
    !is.na(GRADE),
    !is.na(CDCC_TOTAL_BEST),
    !is.na(YEAR_OF_DIAGNOSIS)
  )

dir.create("model", showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) apparent performance
# ----------------------------
lp_all <- predict(fit, newdata = analysis_df, type = "lp")
c_app <- survConcordance(
  Surv(time_months, event) ~ lp_all,
  data = analysis_df
)$concordance

# ----------------------------
# 2) bootstrap internal validation
# ----------------------------
B <- 200  # chosen for stability vs runtime
n <- nrow(analysis_df)

boot_c_app  <- numeric(B)
boot_c_test <- numeric(B)

for (b in seq_len(B)) {
  idx_b <- sample.int(n, size = n, replace = TRUE)
  d_b   <- analysis_df[idx_b, ]

  fit_b <- coxph(
    Surv(time_months, event) ~
      AGE + SEX + GRADE + CDCC_TOTAL_BEST +
      TUMOR_SIZE_CAT + RACE_CAT + YEAR_OF_DIAGNOSIS,
    data = d_b
  )

  # apparent C in bootstrap sample
  lp_b_in <- predict(fit_b, newdata = d_b, type = "lp")
  boot_c_app[b] <- survConcordance(
    Surv(time_months, event) ~ lp_b_in,
    data = d_b
  )$concordance

  # test C when applying bootstrap-fit to original sample
  lp_b_out <- predict(fit_b, newdata = analysis_df, type = "lp")
  boot_c_test[b] <- survConcordance(
    Surv(time_months, event) ~ lp_b_out,
    data = analysis_df
  )$concordance
}

optimism <- mean(boot_c_app - boot_c_test, na.rm = TRUE)
c_optimism_corrected <- c_app - optimism

message("Apparent C-index: ", round(c_app, 3))
message("Bootstrap optimism: ", round(optimism, 3))
message("Optimism-corrected C-index: ", round(c_optimism_corrected, 3))

write.csv(
  data.frame(
    metric = c("C_index_apparent",
               "C_index_optimism_corrected",
               "bootstrap_optimism",
               "B"),
    value  = c(c_app, c_optimism_corrected, optimism, B)
  ),
  "model/performance_summary.csv",
  row.names = FALSE
)

# ----------------------------
# 3) proportional hazards diagnostics
# ----------------------------
zph <- cox.zph(fit)
capture.output(zph, file = "model/cox_zph.txt")

png("model/cox_zph_plots.png", width = 1200, height = 800)
par(mfrow = c(2, 3))
plot(zph)
dev.off()

# ----------------------------
# 4) calibration at 1 / 3 / 5 years
# ----------------------------
get_calibration_df <- function(data, horizon_months, n_groups = 10) {
  sf <- survfit(fit, newdata = data)

  t_index <- max(which(sf$time <= horizon_months))
  predS <- if (length(t_index) == 0) {
    rep(NA_real_, ncol(sf$surv))
  } else {
    as.numeric(sf$surv[t_index, ])
  }

  tmp <- data %>%
    mutate(pred_surv = predS) %>%
    filter(!is.na(pred_surv)) %>%
    mutate(group = ntile(pred_surv, n_groups))

  out <- tmp %>%
    group_by(group) %>%
    summarise(
      n = n(),
      pred_mean = mean(pred_surv),
      obs = {
        sfg <- survfit(Surv(time_months, event) ~ 1, data = cur_data())
        idx <- max(which(sfg$time <= horizon_months))
        if (length(idx) == 0) NA_real_ else as.numeric(sfg$surv[idx])
      },
      .groups = "drop"
    )

  out$horizon <- horizon_months
  out
}

cal12 <- get_calibration_df(analysis_df, 12)
cal36 <- get_calibration_df(analysis_df, 36)
cal60 <- get_calibration_df(analysis_df, 60)

cal_all <- bind_rows(cal12, cal36, cal60)
write.csv(cal_all, "model/calibration_table.csv", row.names = FALSE)

plot_cal <- function(cal_df, horizon) {
  ggplot(
    cal_df %>% filter(horizon == horizon),
    aes(x = pred_mean, y = obs)
  ) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = paste0("Calibration at ", horizon, " months"),
      x = "Mean predicted survival",
      y = "Observed KM survival"
    ) +
    theme_minimal()
}

ggsave("model/calibration_12m_simple.png",
       plot_cal(cal_all, 12),
       width = 7, height = 5, dpi = 200)

ggsave("model/calibration_36m_simple.png",
       plot_cal(cal_all, 36),
       width = 7, height = 5, dpi = 200)

ggsave("model/calibration_60m_simple.png",
       plot_cal(cal_all, 60),
       width = 7, height = 5, dpi = 200)

# ----------------------------
# 5) coefficient table (for reference / manuscript)
# ----------------------------
hr <- summary(fit)$coef
ci <- summary(fit)$conf.int

coef_table <- data.frame(
  term   = rownames(hr),
  beta  = hr[, "coef"],
  HR    = ci[, "exp(coef)"],
  CI_lo = ci[, "lower .95"],
  CI_hi = ci[, "upper .95"],
  p     = hr[, "Pr(>|z|)"],
  row.names = NULL
)

write.csv(coef_table, "model/cox_hr_table.csv", row.names = FALSE)

message("Saved outputs:")
message(" - model/performance_summary.csv")
message(" - model/cox_zph.txt and model/cox_zph_plots.png")
message(" - model/calibration_table.csv and calibration_*_simple.png")
message(" - model/cox_hr_table.csv")

invisible(TRUE)
