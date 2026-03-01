# R/04_data_dictionary.R
# Step 4: Empirical data dictionary for key variables

library(dplyr)

df <- readRDS("data/processed/analysis_df.rds")

cat("ROWS:", nrow(df), "\n")
cat("COLS:", ncol(df), "\n\n")

# -----------------------------
# Helper functions
# -----------------------------
summarize_numeric <- function(x) {
  tibble(
    type = "numeric",
    n_nonmissing = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    min = suppressWarnings(min(x, na.rm = TRUE)),
    median = suppressWarnings(median(x, na.rm = TRUE)),
    mean = suppressWarnings(mean(x, na.rm = TRUE)),
    max = suppressWarnings(max(x, na.rm = TRUE))
  )
}

summarize_categorical <- function(x) {
  top_vals <- sort(table(x), decreasing = TRUE)
  tibble(
    type = "categorical",
    n_nonmissing = sum(!is.na(x)),
    n_missing = sum(is.na(x)),
    n_levels = length(unique(x)),
    top_values = paste(
      names(top_vals)[1:min(5, length(top_vals))],
      collapse = ", "
    )
  )
}

# -----------------------------
# Choose important variables
# (baseline + outcome only)
# -----------------------------
important_vars <- c(
  # outcomes
  "time_months", "event",
  
  # demographics
  "AGE", "SEX", "RACE", "SPANISH_HISPANIC_ORIGIN", "INSURANCE_STATUS",
  
  # tumor
  "PRIMARY_SITE", "LATERALITY", "GRADE", "BEHAVIOR",
  "TUMOR_SIZE", "TUMOR_SIZE_SUMMARY_16",
  
  # comorbidity
  "CDCC_TOTAL_BEST",
  
  # staging
  "TNM_CLIN_T", "TNM_CLIN_N", "TNM_CLIN_M", "TNM_CLIN_STAGE_GROUP",
  "ANALYTIC_STAGE_GROUP",
  
  # facility / SES
  "FACILITY_TYPE_CD", "FACILITY_LOCATION_CD",
  "MED_INC_QUAR_16", "NO_HSD_QUAR_16",
  
  # time
  "YEAR_OF_DIAGNOSIS"
)

important_vars <- intersect(important_vars, names(df))

# -----------------------------
# Build dictionary
# -----------------------------
dictionary <- lapply(important_vars, function(v) {
  x <- df[[v]]
  
  if (is.numeric(x)) {
    cbind(variable = v, summarize_numeric(x))
  } else {
    cbind(variable = v, summarize_categorical(x))
  }
}) %>%
  bind_rows()

print(dictionary)

write.csv(dictionary,
          "data/processed/data_dictionary.csv",
          row.names = FALSE)

cat("\nData dictionary saved to data/processed/data_dictionary.csv\n")