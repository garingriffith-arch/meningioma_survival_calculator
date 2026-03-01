# R/06_validate_calibrate_report.R
# Internal validation + calibration + diagnostics for the SIMPLE Cox model from 05

library(dplyr)
library(survival)
library(ggplot2)

set.seed(20260227)

# ------------------------------------------------------------
# load model + raw data
# ------------------------------------------------------------
fit <- readRDS("model/cox_simple_fit.rds")
df  <- readRDS("data/processed/analysis_df.rds")

dir.create("model", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# recreate analysis_df EXACTLY like 05
# notes-to-self: if I change recoding in 05, I MUST copy it here too
# ------------------------------------------------------------
analysis_df <- df %>%
  mutate(
    # outcome
    time_months = as.numeric(time_months),
    event       = as.integer(event),

    # demographics
    AGE = as.numeric(AGE),
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),

    # grade
    GRADE = factor(GRADE, levels = c(2, 3), labels = c("II", "III")),

    # comorbidity
    CDCC_TOTAL_BEST = as.numeric(CDCC_TOTAL_BEST),

    # secular trend
    YEAR_OF_DIAGNOSIS = as.numeric(YEAR_OF_DIAGNOSIS),

    # tumor size (collapse hard)
    TUMOR_SIZE_MM = if_else(TUMOR_SIZE >= 999, NA_real_, as.numeric(TUMOR_SIZE)),
    TUMOR_SIZE_CAT = case_when(
      is.na(TUMOR_SIZE_MM) ~ "Unknown",
      TUMOR_SIZE_MM <= 30  ~ "<=3cm",
      TUMOR_SIZE_MM <= 60  ~ "3–6cm",
      TRUE                 ~ ">6cm"
    ),
    TUMOR_SIZE_CAT = factor(TUMOR_SIZE_CAT,
                            levels = c("<=3cm", "3–6cm", ">6cm", "Unknown")),

    # race (collapse hard)
    RACE_CAT = case_when(
      RACE == 1 ~ "White",
      RACE == 2 ~ "Black",
      TRUE      ~ "Other/Unknown"
    ),
    RACE_CAT = factor(RACE_CAT)
  ) %>%
  filter(
    !is.na(time_months), !is.na(event),
    !is.na(AGE), !is.na(SEX), !is.na(GRADE),
    !is.na(CDCC_TOTAL_BEST), !is.na(YEAR_OF_DIAGNOSIS)
  ) %>%
  droplevels()

cat("Analysis rows:", nrow(analysis_df), "\n")

# ------------------------------------------------------------
# 1) apparent performance (C-index)
# ------------------------------------------------------------
lp_all <- predict(fit, newdata = analysis_df, type = "lp")
c_app  <- survConcordance(Surv(time_months, event) ~ lp_all, data = analysis_df)$concordance

# ------------------------------------------------------------
# 2) bootstrap internal validation (optimism-corrected C-index)
# ------------------------------------------------------------
B <- 200  # notes-to-self: for final paper bump to 500–1000
n <- nrow(analysis_df)

boot_c_app  <- numeric(B)
boot_c_test <- numeric(B)

form <- Surv(time_months, event) ~
  AGE + SEX + GRADE + CDCC_TOTAL_BEST + TUMOR_SIZE_CAT + RACE_CAT + YEAR_OF_DIAGNOSIS

for (b in seq_len(B)) {
  idx_b <- sample.int(n, size = n, replace = TRUE)
  d_b   <- analysis_df[idx_b, , drop = FALSE] %>% droplevels()

  # refit in bootstrap sample
  fit_b <- coxph(form, data = d_b)

  # apparent (in-bootstrap) C
  lp_in <- predict(fit_b, newdata = d_b, type = "lp")
  boot_c_app[b] <- survConcordance(Surv(time_months, event) ~ lp_in, data = d_b)$concordance

  # test (apply bootstrap-fit to original full sample) C
  lp_out <- predict(fit_b, newdata = analysis_df, type = "lp")
  boot_c_test[b] <- survConcordance(Surv(time_months, event) ~ lp_out, data = analysis_df)$concordance
}

optimism <- mean(boot_c_app - boot_c_test, na.rm = TRUE)
c_corr   <- c_app - optimism

cat("\nApparent C-index:", round(c_app, 3), "\n")
cat("Bootstrap optimism:", round(optimism, 3), "\n")
cat("Optimism-corrected C-index:", round(c_corr, 3), "\n")

write.csv(
  data.frame(
    metric = c("C_index_apparent", "C_index_optimism_corrected", "bootstrap_optimism", "B"),
    value  = c(c_app, c_corr, optimism, B)
  ),
  "model/performance_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 3) PH assumption checks (global + per covariate)
# notes-to-self:
#   - you already have a zph output file in model/cox_zph.txt
#   - in that output, YEAR_OF_DIAGNOSIS and GLOBAL were extremely significant,
#     and TUMOR_SIZE_CAT showed signal too. (double-check interpretation there)
# ------------------------------------------------------------
zph <- cox.zph(fit)
capture.output(zph, file = "model/cox_zph.txt")

# plot layout: adapt to number of terms
png("model/cox_zph_plots.png", width = 1400, height = 900)
k <- nrow(zph$table) - 1  # drop GLOBAL
ncol_plot <- 3
nrow_plot <- ceiling((k + 1) / ncol_plot)  # +1 if GLOBAL included by plot()
par(mfrow = c(nrow_plot, ncol_plot))
plot(zph)
dev.off()

# ------------------------------------------------------------
# helper: predicted survival at a horizon for each subject
# notes-to-self:
#   - safest extraction is summary(survfit(...), times=...)
#   - returns a vector of length (n_subjects * length(times)) in subject blocks
# ------------------------------------------------------------
pred_surv_at <- function(fit_obj, newdata, horizon_months) {
  sf <- survfit(fit_obj, newdata = newdata)

  # summary() handles interpolation + edge cases better than manual indexing
  s <- summary(sf, times = horizon_months)

  # When newdata has N rows and times has length 1:
  # s$surv comes back length N (in order), or length 1 if N=1.
  as.numeric(s$surv)
}

# helper: observed KM survival at horizon for a dataset
obs_km_at <- function(data, horizon_months) {
  sfg <- survfit(Surv(time_months, event) ~ 1, data = data)
  s   <- summary(sfg, times = horizon_months)

  # if horizon precedes first event time, survival is basically 1
  if (length(s$surv) == 0) return(1)
  as.numeric(s$surv)
}

# ------------------------------------------------------------
# 4) calibration @ 12/36/60 months
# approach:
#   - compute predicted S_i(t) for each subject
#   - bin into deciles by predicted survival (low S = higher risk)
#   - within each bin: mean predicted vs KM observed at same horizon
# ------------------------------------------------------------
get_calibration_df <- function(data, horizon_months, n_groups = 10) {
  predS <- pred_surv_at(fit, data, horizon_months)

  tmp <- data %>%
    mutate(pred_surv = predS) %>%
    filter(!is.na(pred_surv)) %>%
    # notes-to-self: group 1 = lowest predicted survival (highest risk)
    mutate(group = ntile(pred_surv, n_groups))

  out <- tmp %>%
    group_by(group) %>%
    summarise(
      n = n(),
      pred_mean = mean(pred_surv),
      obs = obs_km_at(cur_data(), horizon_months),
      .groups = "drop"
    ) %>%
    mutate(horizon = horizon_months)

  out
}

cal_all <- bind_rows(
  get_calibration_df(analysis_df, 12),
  get_calibration_df(analysis_df, 36),
  get_calibration_df(analysis_df, 60)
)

write.csv(cal_all, "model/calibration_table.csv", row.names = FALSE)

plot_cal <- function(cal_df, horizon_val) {
  ggplot(cal_df %>% filter(.data$horizon == horizon_val),
         aes(x = pred_mean, y = obs)) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = paste0("Calibration at ", horizon_val, " months"),
      x = "Mean predicted survival",
      y = "Observed KM survival"
    ) +
    theme_minimal()
}

ggsave("model/calibration_12m_simple.png", plot_cal(cal_all, 12), width = 7, height = 5, dpi = 200)
ggsave("model/calibration_36m_simple.png", plot_cal(cal_all, 36), width = 7, height = 5, dpi = 200)
ggsave("model/calibration_60m_simple.png", plot_cal(cal_all, 60), width = 7, height = 5, dpi = 200)

# ------------------------------------------------------------
# 5) coefficient table for manuscript (HR + 95% CI)
# ------------------------------------------------------------
s <- summary(fit)

hr <- s$coef
ci <- s$conf.int

coef_table <- data.frame(
  term    = rownames(hr),
  beta    = hr[, "coef"],
  HR      = ci[, "exp(coef)"],
  CI_low  = ci[, "lower .95"],
  CI_high = ci[, "upper .95"],
  p       = hr[, "Pr(>|z|)"],
  row.names = NULL
)

write.csv(coef_table, "model/cox_hr_table.csv", row.names = FALSE)

cat("\nSaved:\n",
    "- model/performance_summary.csv\n",
    "- model/cox_zph.txt and model/cox_zph_plots.png\n",
    "- model/calibration_table.csv and calibration_*_simple.png\n",
    "- model/cox_hr_table.csv\n")
