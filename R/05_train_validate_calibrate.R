# R/05_train_validate_calibrate.R
# SIMPLE, STABLE Cox survival model for Grade II–III intracranial meningioma
# No rms, no splines, no fragile factor handling

library(dplyr)
library(survival)

set.seed(20260227)

# ----------------------------
# Load analytic dataset
# ----------------------------
df <- readRDS("data/processed/analysis_df.rds")

# ----------------------------
# Minimal, defensible recoding
# ----------------------------
analysis_df <- df %>%
  mutate(
    # outcome
    time_months = as.numeric(time_months),
    event = as.integer(event),
    
    # demographics
    AGE = as.numeric(AGE),
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    
    # grade
    GRADE = factor(GRADE, levels = c(2, 3), labels = c("II", "III")),
    
    # comorbidity
    CDCC_TOTAL_BEST = as.numeric(CDCC_TOTAL_BEST),
    
    # year (secular trend)
    YEAR_OF_DIAGNOSIS = as.numeric(YEAR_OF_DIAGNOSIS),
    
    # ----------------------------
    # Tumor size: collapse HARD
    # ----------------------------
    TUMOR_SIZE_MM = ifelse(TUMOR_SIZE >= 999, NA, TUMOR_SIZE),
    TUMOR_SIZE_CAT = case_when(
      is.na(TUMOR_SIZE_MM)        ~ "Unknown",
      TUMOR_SIZE_MM <= 30         ~ "<=3cm",
      TUMOR_SIZE_MM <= 60         ~ "3–6cm",
      TUMOR_SIZE_MM > 60          ~ ">6cm"
    ),
    TUMOR_SIZE_CAT = factor(
      TUMOR_SIZE_CAT,
      levels = c("<=3cm", "3–6cm", ">6cm", "Unknown")
    ),
    
    # ----------------------------
    # Race: collapse HARD
    # ----------------------------
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

cat("Final analytic rows:", nrow(analysis_df), "\n")

# ----------------------------
# Train / test split (70/30)
# ----------------------------
n <- nrow(analysis_df)
idx <- sample(seq_len(n))

train <- analysis_df[idx[1:floor(0.7 * n)], ]
test  <- analysis_df[idx[(floor(0.7 * n) + 1):n], ]

# ----------------------------
# Cox proportional hazards model
# ----------------------------
cox_fit <- coxph(
  Surv(time_months, event) ~
    AGE + SEX + GRADE + CDCC_TOTAL_BEST +
    TUMOR_SIZE_CAT + RACE_CAT + YEAR_OF_DIAGNOSIS,
  data = train,
  x = TRUE,
  y = TRUE
)


summary(cox_fit)

# ----------------------------
# Test-set discrimination (C-index)
# ----------------------------
lp_test <- predict(cox_fit, newdata = test, type = "lp")

c_index <- survConcordance(
  Surv(test$time_months, test$event) ~ lp_test
)$concordance

cat("Test-set Harrell C-index:", round(c_index, 3), "\n")

# ----------------------------
# Save model for later use
# ----------------------------
dir.create("model", showWarnings = FALSE, recursive = TRUE)
saveRDS(cox_fit, "model/cox_simple_fit.rds")

# ----------------------------
# Survival probabilities at 1 / 3 / 5 years
# ----------------------------
times <- c(12, 36, 60)
sf <- survfit(cox_fit, newdata = test)

get_surv <- function(sf_obj, t) {
  idx <- max(which(sf_obj$time <= t))
  if (length(idx) == 0) return(rep(NA, ncol(sf_obj$surv)))
  sf_obj$surv[idx, ]
}

pred_df <- test %>%
  mutate(
    S_12 = get_surv(sf, 12),
    S_36 = get_surv(sf, 36),
    S_60 = get_surv(sf, 60)
  )

write.csv(pred_df, "model/test_predictions_simple.csv", row.names = FALSE)

cat("Saved model and predictions.\n")
