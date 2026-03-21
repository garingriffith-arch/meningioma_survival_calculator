# R/04_make_tables_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(survival)
  library(tibble)
})

cfg <- source("R/config.R")$value
source("R/utils_survival.R")

dir.create(cfg$tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$tables_dir, "manuscript"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$figures_dir, "manuscript"), recursive = TRUE, showWarnings = FALSE)

df      <- readRDS(file.path(cfg$processed_dir, "01_analytic_cohort.rds"))
fit     <- readRDS(file.path(cfg$processed_dir, "03_cox_model.rds"))
val_res <- readRDS(file.path(cfg$processed_dir, "04_validation_results.rds"))
flow    <- read_csv(file.path(cfg$processed_dir, "01_cohort_flow.csv"), show_col_types = FALSE)
val_sum <- read_csv(file.path(cfg$processed_dir, "04_validation_summary.csv"), show_col_types = FALSE)

fmt_n_pct <- function(x) {
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
}

fmt_mean_sd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
}

fmt_med_iqr <- function(x) {
  q <- stats::quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE)
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

fmt_p <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

# -----------------------------
# Table 1: Cohort characteristics
# -----------------------------
tab1 <- tibble(
  group = "Overall",
  N = nrow(df),
  age_mean_sd = fmt_mean_sd(df$age),
  age_median_iqr = fmt_med_iqr(df$age),
  male = fmt_n_pct(df$SEX == "Male"),
  female = fmt_n_pct(df$SEX == "Female"),
  cdcc_0 = fmt_n_pct(df$cdcc == 0),
  cdcc_1 = fmt_n_pct(df$cdcc == 1),
  cdcc_2plus = fmt_n_pct(df$cdcc >= 2),
  grade_II = fmt_n_pct(df$grade23 == "2"),
  grade_III = fmt_n_pct(df$grade23 == "3"),
  race_white = fmt_n_pct(df$race == "White"),
  race_black = fmt_n_pct(df$race == "Black"),
  race_other = fmt_n_pct(df$race == "Other"),
  ethnicity_non_hispanic = fmt_n_pct(df$ethnicity == "Non-Hispanic"),
  ethnicity_hispanic = fmt_n_pct(df$ethnicity == "Hispanic"),
  tumor_size_mm_median_iqr = fmt_med_iqr(df$tumor_size_mm),
  followup_months_median_iqr = fmt_med_iqr(df$time_months),
  deaths = sum(df$event == 1, na.rm = TRUE)
)

write_csv(
  tab1,
  file.path(cfg$tables_dir, "manuscript", "Table1_cohort_characteristics.csv")
)

# -----------------------------
# Table 2: Cohort flow
# -----------------------------
write_csv(
  flow,
  file.path(cfg$tables_dir, "manuscript", "Table2_cohort_flow.csv")
)

# -----------------------------
# Table 3: Multivariable model results
# -----------------------------
ref <- make_reference_profile(df)

ref_age50 <- ref
ref_age50$age <- 50
new_age70 <- ref
new_age70$age <- 70

ref_size30 <- ref
ref_size30$tumor_size_mm <- 30
new_size60 <- ref
new_size60$tumor_size_mm <- 60

age_contrast <- contrast_hr_ci(fit, ref_age50, new_age70) |>
  mutate(term = "Age: 70 vs 50 years")

size_contrast <- contrast_hr_ci(fit, ref_size30, new_size60) |>
  mutate(term = "Tumor size: 60 vs 30 mm")

contrast_tbl <- bind_rows(
  age_contrast |>
    transmute(
      term,
      HR = round(HR, 3),
      lo = round(lo, 3),
      hi = round(hi, 3),
      p = NA_character_
    ),
  size_contrast |>
    transmute(
      term,
      HR = round(HR, 3),
      lo = round(lo, 3),
      hi = round(hi, 3),
      p = NA_character_
    )
)

sfit <- summary(fit)
coef_mat <- as.data.frame(sfit$coefficients)
ci_mat   <- as.data.frame(sfit$conf.int)

cat_tbl <- tibble(
  term = rownames(coef_mat),
  HR = ci_mat$`exp(coef)`,
  lo = ci_mat$`lower .95`,
  hi = ci_mat$`upper .95`,
  p  = coef_mat$`Pr(>|z|)`
) |>
  filter(!grepl("^splines::ns\\(", term)) |>
  mutate(
    term = recode(
      term,
      "SEXFemale" = "Female sex (vs male)",
      "cdcc" = "Charlson-Deyo score (per point)",
      "grade233" = "WHO grade III (vs II)",
      "raceBlack" = "Black race (vs White)",
      "raceOther" = "Other race (vs White)",
      "ethnicityHispanic" = "Hispanic ethnicity (vs non-Hispanic)",
      .default = term
    ),
    HR = round(HR, 3),
    lo = round(lo, 3),
    hi = round(hi, 3),
    p = fmt_p(p)
  )

# robust extraction of overall and nonlinearity p-values
drop_tab <- stats::drop1(fit, test = "Chisq")
drop_p_col <- grep("^Pr\\(", colnames(drop_tab), value = TRUE)[1]

overall_age_p <- drop_tab[
  grep("^splines::ns\\(age", rownames(drop_tab)),
  drop_p_col
][1]

overall_size_p <- drop_tab[
  grep("^splines::ns\\(tumor_size_mm", rownames(drop_tab)),
  drop_p_col
][1]

fit_age_linear <- stats::update(
  fit,
  . ~ . - splines::ns(age, df = cfg$age_spline_df) + age
)

fit_size_linear <- stats::update(
  fit,
  . ~ . - splines::ns(tumor_size_mm, df = cfg$size_spline_df) + tumor_size_mm
)

an_age <- stats::anova(fit_age_linear, fit, test = "Chisq")
an_size <- stats::anova(fit_size_linear, fit, test = "Chisq")

anova_p_col_age <- grep("^P\\(", colnames(an_age), value = TRUE)[1]
anova_p_col_size <- grep("^P\\(", colnames(an_size), value = TRUE)[1]

nonlinear_age_p <- an_age[[anova_p_col_age]][2]
nonlinear_size_p <- an_size[[anova_p_col_size]][2]

p_tbl <- tibble(
  term = c("Age spline term", "Tumor size spline term"),
  HR = NA_real_,
  lo = NA_real_,
  hi = NA_real_,
  p = c(
    paste0(
      "overall=", fmt_p(overall_age_p),
      "; nonlinear=", fmt_p(nonlinear_age_p)
    ),
    paste0(
      "overall=", fmt_p(overall_size_p),
      "; nonlinear=", fmt_p(nonlinear_size_p)
    )
  )
)

tab3 <- bind_rows(contrast_tbl, cat_tbl, p_tbl)

write_csv(
  tab3,
  file.path(cfg$tables_dir, "manuscript", "Table3_cox_model_effects.csv")
)

# -----------------------------
# Table 4: Horizon-specific validation summary
# -----------------------------
tab4 <- val_sum |>
  mutate(
    across(
      c(
        auc_apparent, auc_optimism, auc_corrected,
        brier_apparent, brier_optimism, brier_corrected
      ),
      ~ round(.x, 3)
    )
  )

write_csv(
  tab4,
  file.path(cfg$tables_dir, "manuscript", "Table4_validation_12_36_60m.csv")
)

# -----------------------------
# Figure 1: Calibration plots
# -----------------------------
cal_dir <- file.path(cfg$processed_dir, "calibration")
cal_files <- list.files(cal_dir, pattern = "^calibration_\\d+m\\.csv$", full.names = TRUE)

cal_all <- bind_rows(lapply(cal_files, function(f) {
  horizon <- as.numeric(gsub("calibration_(\\d+)m\\.csv", "\\1", basename(f)))
  read_csv(f, show_col_types = FALSE) |>
    mutate(horizon = horizon)
}))

make_cal_plot <- function(h) {
  dat <- filter(cal_all, horizon == h)
  
  ggplot(dat, aes(x = pred_surv_mean, y = km_surv)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 2) +
    geom_line() +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(labels = percent) +
    scale_y_continuous(labels = percent) +
    labs(
      title = paste("Calibration at", h, "months"),
      x = "Predicted overall survival",
      y = "Observed overall survival"
    ) +
    theme_minimal()
}

for (h in cfg$horizons_months) {
  ggsave(
    filename = file.path(cfg$figures_dir, "manuscript", paste0("Figure1_calibration_", h, "m.png")),
    plot = make_cal_plot(h),
    width = 6,
    height = 5,
    dpi = 300
  )
}

# -----------------------------
# Figure 2: Predicted survival by risk group
# -----------------------------
lp <- as.numeric(predict(fit, newdata = df, type = "lp"))
q <- quantile(lp, probs = c(0.2, 0.8), na.rm = TRUE)

df_risk <- df |>
  mutate(
    risk_group = case_when(
      lp <= q[1] ~ "Low risk",
      lp >= q[2] ~ "High risk",
      TRUE ~ "Intermediate risk"
    )
  )

profiles <- df_risk |>
  group_by(risk_group) |>
  summarise(
    age = median(age, na.rm = TRUE),
    tumor_size_mm = median(tumor_size_mm, na.rm = TRUE),
    cdcc = median(cdcc, na.rm = TRUE),
    grade23 = safe_mode(grade23),
    SEX = safe_mode(SEX),
    race = safe_mode(race),
    ethnicity = safe_mode(ethnicity),
    .groups = "drop"
  ) |>
  mutate(
    grade23 = factor(grade23, levels = levels(df$grade23)),
    SEX = factor(SEX, levels = levels(df$SEX)),
    race = factor(race, levels = levels(df$race)),
    ethnicity = factor(ethnicity, levels = levels(df$ethnicity))
  )

sf_df <- bind_rows(lapply(seq_len(nrow(profiles)), function(i) {
  sf <- survfit(fit, newdata = profiles[i, ])
  data.frame(
    time = sf$time,
    surv = sf$surv,
    lower = sf$lower,
    upper = sf$upper,
    risk_group = profiles$risk_group[i]
  )
}))

p_surv <- ggplot(sf_df, aes(x = time, y = surv, color = risk_group, fill = risk_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
  geom_step(linewidth = 1.1) +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Time (months)",
    y = "Predicted overall survival",
    color = "Risk group",
    fill = "Risk group"
  ) +
  coord_cartesian(xlim = c(0, max(cfg$horizons_months))) +
  theme_minimal()

ggsave(
  filename = file.path(cfg$figures_dir, "manuscript", "Figure2_predicted_survival_risk_groups.png"),
  plot = p_surv,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------
# Figure 3: Adjusted spline effect for age
# -----------------------------
age_df <- make_spline_effect_df(
  fit = fit,
  df = df,
  variable = "age",
  ref_value = cfg$age_ref_years
)

p_age <- ggplot(age_df, aes(x = x, y = HR)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    x = "Age (years)",
    y = paste0("Hazard ratio (reference = ", cfg$age_ref_years, " years)"),
    title = "Adjusted nonlinear association of age with overall survival"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(cfg$figures_dir, "manuscript", "Figure3_age_spline_effect.png"),
  plot = p_age,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# Figure 4: Adjusted spline effect for tumor size
# -----------------------------
size_df <- make_spline_effect_df(
  fit = fit,
  df = df,
  variable = "tumor_size_mm",
  ref_value = cfg$size_ref_mm
)

p_size <- ggplot(size_df, aes(x = x, y = HR)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    x = "Tumor size (mm)",
    y = paste0("Hazard ratio (reference = ", cfg$size_ref_mm, " mm)"),
    title = "Adjusted nonlinear association of tumor size with overall survival"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(cfg$figures_dir, "manuscript", "Figure4_tumor_size_spline_effect.png"),
  plot = p_size,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# Table 5: Global validation summary
# -----------------------------
tab5 <- tibble(
  metric = c("Harrell's C", "Calibration slope"),
  apparent = c(
    val_res$apparent$harrell_c,
    val_res$apparent$calibration_slope
  ),
  optimism = c(
    val_res$bootstrap$optimism$mean_c_global,
    val_res$bootstrap$optimism$mean_slope_global
  ),
  corrected = c(
    val_res$bootstrap$corrected$harrell_c,
    val_res$bootstrap$corrected$calibration_slope
  )
) |>
  mutate(
    apparent = round(apparent, 3),
    optimism = round(optimism, 3),
    corrected = round(corrected, 3)
  )

write_csv(
  tab5,
  file.path(cfg$tables_dir, "manuscript", "Table5_global_validation_summary.csv")
)
