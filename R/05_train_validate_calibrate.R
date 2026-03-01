# R/05_train_validate_calibrate.R
# SIMPLE, STABLE Cox model for Grade II–III intracranial meningioma

library(dplyr)
library(survival)

set.seed(20260227)

# ------------------------------------------------------------
# load analytic dataset
# ------------------------------------------------------------
df <- readRDS("data/processed/analysis_df.rds")

# ------------------------------------------------------------
# minimal, defensible recoding
# ------------------------------------------------------------
analysis_df <- df %>%
  mutate(
    # outcome
    time_months = as.numeric(time_months),
    event       = as.integer(event),

    # demographics
    AGE = as.numeric(AGE),
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),

    # grade (this analysis is only II–III anyway)
    GRADE = factor(GRADE, levels = c(2, 3), labels = c("II", "III")),

    # comorbidity
    CDCC_TOTAL_BEST = as.numeric(CDCC_TOTAL_BEST),

    # year (just to catch secular drift)
    YEAR_OF_DIAGNOSIS = as.numeric(YEAR_OF_DIAGNOSIS),

    # --------------------------------------------------------
    # tumor size: collapse hard (NCDB-ish 999 codes etc.)
    # notes-to-self: treat >=999 as unknown/missing
    # --------------------------------------------------------
    TUMOR_SIZE_MM = if_else(TUMOR_SIZE >= 999, NA_real_, as.numeric(TUMOR_SIZE)),
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

    # --------------------------------------------------------
    # race: collapse hard (keep it simple + stable)
    # --------------------------------------------------------
    RACE_CAT = case_when(
      RACE == 1 ~ "White",
      RACE == 2 ~ "Black",
      TRUE      ~ "Other/Unknown"
    ),
    RACE_CAT = factor(RACE_CAT)
  ) %>%
  # notes-to-self: drop rows missing core predictors/outcome
  filter(
    !is.na(time_months),
    !is.na(event),
    !is.na(AGE),
    !is.na(SEX),
    !is.na(GRADE),
    !is.na(CDCC_TOTAL_BEST),
    !is.na(YEAR_OF_DIAGNOSIS)
  )

cat("Final analytic rows:", nrow(analysis_df), "\n")

# ------------------------------------------------------------
# train / test split (70/30)
# ------------------------------------------------------------
n <- nrow(analysis_df)
idx <- sample.int(n)

n_train <- floor(0.7 * n)
train <- analysis_df[idx[seq_len(n_train)], ]
test  <- analysis_df[idx[(n_train + 1):n], ]

cat("Train rows:", nrow(train), " | Test rows:", nrow(test), "\n")

# ------------------------------------------------------------
# Cox proportional hazards model
# ------------------------------------------------------------
cox_fit <- coxph(
  Surv(time_months, event) ~
    AGE + SEX + GRADE + CDCC_TOTAL_BEST +
    TUMOR_SIZE_CAT + RACE_CAT + YEAR_OF_DIAGNOSIS,
  data = train,
  x = TRUE,
  y = TRUE
)

print(summary(cox_fit))

# ------------------------------------------------------------
# test-set discrimination (Harrell's C)
# notes-to-self: basic sanity check; not the full story
# ------------------------------------------------------------
lp_test <- predict(cox_fit, newdata = test, type = "lp")

c_index <- survConcordance(
  Surv(time_months, event) ~ lp_test,
  data = test
)$concordance

cat("Test-set Harrell C-index:", round(c_index, 3), "\n")

# ------------------------------------------------------------
# save model (so downstream scripts don't refit)
# ------------------------------------------------------------
dir.create("model", showWarnings = FALSE, recursive = TRUE)
saveRDS(cox_fit, file = "model/cox_simple_fit.rds")

# ------------------------------------------------------------
# survival probabilities at 1 / 3 / 5 years
# notes-to-self:
#   - survfit(cox_fit, newdata=test) gives one curve per row
#   - summary(sf, times=...) is the safest way to extract S(t)
# ------------------------------------------------------------
times <- c(12, 36, 60)  # months = 1/3/5 years
sf <- survfit(cox_fit, newdata = test)

# returns a list; $surv is a vector of length (nrow(test)*length(times))
sf_sum <- summary(sf, times = times)

# reshape into nrow(test) x length(times)
# summary() repeats times for each subject, in order
s_mat <- matrix(sf_sum$surv, nrow = nrow(test), ncol = length(times), byrow = TRUE)
colnames(s_mat) <- paste0("S_", times)

pred_df <- bind_cols(
  test,
  as.data.frame(s_mat)
)

write.csv(pred_df, "model/test_predictions_simple.csv", row.names = FALSE)

cat("Saved model + test predictions to model/.\n")

