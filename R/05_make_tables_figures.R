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

df       <- readRDS(file.path(cfg$processed_dir, "01_analytic_cohort.rds"))
df_split <- readRDS(file.path(cfg$processed_dir, "02_analytic_cohort_split.rds"))
fit      <- readRDS(file.path(cfg$processed_dir, "03_cox_model.rds"))
val_res  <- readRDS(file.path(cfg$processed_dir, "04_validation_results.rds"))
flow     <- read_csv(file.path(cfg$processed_dir, "01_cohort_flow.csv"), show_col_types = FALSE)
val_sum  <- read_csv(file.path(cfg$processed_dir, "04_validation_summary.csv"), show_col_types = FALSE)

fmt_n_pct <- function(x) {
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
}

fmt_mean_sd <- function(x) {
  sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_med_iqr <- function(x) {
  q <- quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE)
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}

make_table1 <- function(dat, label) {
  tibble(
    group = label,
    N = nrow(dat),
    age_mean_sd = fmt_mean_sd(dat$age),
    age_median_iqr = fmt_med_iqr(dat$age),
    male = fmt_n_pct(dat$SEX == "Male"),
    female = fmt_n_pct(dat$SEX == "Female"),
    cdcc_0 = fmt_n_pct(dat$cdcc == 0),
    cdcc_1 = fmt_n_pct(dat$cdcc == 1),
    cdcc_2plus = fmt_n_pct(dat$cdcc >= 2),
    grade_II = fmt_n_pct(dat$grade23 == "2"),
    grade_III = fmt_n_pct(dat$grade23 == "3"),
    race_white = fmt_n_pct(dat$race == "White"),
    race_black = fmt_n_pct(dat$race == "Black"),
    race_other = fmt_n_pct(dat$race == "Other"),
    ethnicity_non_hispanic = fmt_n_pct(dat$ethnicity == "Non-Hispanic"),
    ethnicity_hispanic = fmt_n_pct(dat$ethnicity == "Hispanic"),
    tumor_size_mm_median_iqr = fmt_med_iqr(dat$tumor_size_mm),
    followup_months_median_iqr = fmt_med_iqr(dat$time_months),
    deaths = sum(dat$event == 1, na.rm = TRUE)
  )
}

tab1 <- bind_rows(
  make_table1(df, "Overall"),
  make_table1(filter(df_split, split == "train"), "Train"),
  make_table1(filter(df_split, split == "test"), "Test")
)

write_csv(
  tab1,
  file.path(cfg$tables_dir, "manuscript", "Table1_cohort_characteristics.csv")
)

write_csv(
  flow,
  file.path(cfg$tables_dir, "manuscript", "Table2_cohort_flow.csv")
)

# Table 3: Cox model effects
ref <- df %>%
  summarise(
    age = median(age, na.rm = TRUE),
    tumor_size_mm = median(tumor_size_mm, na.rm = TRUE),
    cdcc = median(cdcc, na.rm = TRUE),
    grade23 = safe_mode(grade23),
    SEX = safe_mode(SEX),
    race = safe_mode(race),
    ethnicity = safe_mode(ethnicity)
  ) %>%
  mutate(
    grade23 = factor(grade23, levels = levels(df$grade23)),
    SEX = factor(SEX, levels = levels(df$SEX)),
    race = factor(race, levels = levels(df$race)),
    ethnicity = factor(ethnicity, levels = levels(df$ethnicity))
  )

hr_contrast <- function(new_age = ref$age, new_tsize = ref$tumor_size_mm) {
  new <- ref
  new$age <- new_age
  new$tumor_size_mm <- new_tsize
  
  lp_ref <- predict(fit, newdata = ref, type = "lp")
  lp_new <- predict(fit, newdata = new, type = "lp")
  
  as.numeric(exp(lp_new - lp_ref))
}

contrast_tbl <- tibble(
  term = c(
    "Age: 70 vs 50 years",
    "Tumor size: 60 vs 30 mm"
  ),
  HR = c(
    hr_contrast(new_age = 70) / hr_contrast(new_age = 50),
    hr_contrast(new_tsize = 60) / hr_contrast(new_tsize = 30)
  ),
  lo = NA_real_,
  hi = NA_real_,
  p = NA_character_
) %>%
  mutate(
    HR = round(HR, 3)
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
) %>%
  filter(!grepl("^ns\\(", term)) %>%
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
    p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

tab3 <- bind_rows(contrast_tbl, cat_tbl)

write_csv(
  tab3,
  file.path(cfg$tables_dir, "manuscript", "Table3_cox_model_effects.csv")
)

# Table 4: Validation summary
write_csv(
  val_sum,
  file.path(cfg$tables_dir, "manuscript", "Table4_validation_summary.csv")
)

# Table 5: Bootstrap optimism summary
boot_sum <- tibble(
  metric = c("Harrell's C", "Calibration slope"),
  apparent = c(
    val_res$bootstrap$apparent$C,
    val_res$bootstrap$apparent$slope
  ),
  optimism = c(
    val_res$bootstrap$optimism$mean_C,
    val_res$bootstrap$optimism$mean_slope
  ),
  corrected = c(
    val_res$bootstrap$corrected$C,
    val_res$bootstrap$corrected$slope
  )
)

write_csv(
  boot_sum,
  file.path(cfg$tables_dir, "manuscript", "Table5_bootstrap_optimism_summary.csv")
)

# Calibration plots
cal_dir <- file.path(cfg$processed_dir, "calibration")
cal_files <- list.files(cal_dir, pattern = "^calibration_\\d+m\\.csv$", full.names = TRUE)

cal_all <- bind_rows(lapply(cal_files, function(f) {
  horizon <- as.numeric(gsub("calibration_(\\d+)m\\.csv", "\\1", basename(f)))
  read_csv(f, show_col_types = FALSE) |> mutate(horizon = horizon)
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

# Figure 2: Predicted survival curves for representative risk groups
lp <- as.numeric(predict(fit, newdata = df, type = "lp"))
q <- quantile(lp, probs = c(0.2, 0.8), na.rm = TRUE)

df_risk <- df %>%
  mutate(
    risk_group = case_when(
      lp <= q[1] ~ "Low risk",
      lp >= q[2] ~ "High risk",
      TRUE ~ "Intermediate risk"
    )
  )

profiles <- df_risk %>%
  group_by(risk_group) %>%
  summarise(
    age = median(age, na.rm = TRUE),
    tumor_size_mm = median(tumor_size_mm, na.rm = TRUE),
    cdcc = median(cdcc, na.rm = TRUE),
    grade23 = safe_mode(grade23),
    SEX = safe_mode(SEX),
    race = safe_mode(race),
    ethnicity = safe_mode(ethnicity),
    .groups = "drop"
  ) %>%
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

# Figure 3: Discrimination across prediction horizons
cvec <- val_res$c_horizon

dfc <- tibble(
  horizon = as.numeric(gsub("C@_|m", "", names(cvec))),
  c_index = as.numeric(cvec)
)

p_disc <- ggplot(dfc, aes(x = factor(horizon), y = c_index)) +
  geom_col(fill = "#198754") +
  geom_text(
    aes(label = sprintf("%.3f", c_index)),
    vjust = -0.5
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    name = "Time-dependent C-index"
  ) +
  scale_x_discrete(
    name = "Prediction Horizon (months)"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(cfg$figures_dir, "manuscript", "Figure3_discrimination.png"),
  plot = p_disc,
  width = 7,
  height = 5,
  dpi = 300
)

# Figure 4A: Bootstrap optimism in Harrell's C
df_opt_c <- tibble(
  optimism_C = as.numeric(val_res$bootstrap$optimism$C)
) %>%
  filter(is.finite(optimism_C))

p4a <- ggplot(df_opt_c, aes(x = optimism_C)) +
  geom_histogram(bins = 30, fill = "#fd7e14", color = "white") +
  labs(
    title = "Bootstrap optimism in Harrell's C",
    x = "Optimism (C-index)",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(cfg$figures_dir, "manuscript", "Figure4A_bootstrap_optimism_C.png"),
  plot = p4a,
  width = 7,
  height = 5,
  dpi = 300
)

# Figure 4B: Bootstrap optimism in calibration slope
df_opt_s <- tibble(
  optimism_slope = as.numeric(val_res$bootstrap$optimism$slope)
) %>%
  filter(is.finite(optimism_slope))

p4b <- ggplot(df_opt_s, aes(x = optimism_slope)) +
  geom_histogram(bins = 30, fill = "#6f42c1", color = "white") +
  labs(
    title = "Bootstrap optimism in calibration slope",
    x = "Optimism (slope)",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(cfg$figures_dir, "manuscript", "Figure4B_bootstrap_optimism_slope.png"),
  plot = p4b,
  width = 7,
  height = 5,
  dpi = 300
)
