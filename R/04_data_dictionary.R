# R/04_data_dictionary.R
# Step 4: quick-and-dirty empirical data dictionary (baseline + outcome vars)
#
# Goal: for a curated set of variables, spit out:
#   - missingness
#   - basic numeric summaries OR quick categorical overview
# Then save a CSV so I can eyeball what I'm working with.

library(dplyr)

df <- readRDS("data/processed/analysis_df.rds")

cat("ROWS:", nrow(df), "\n")
cat("COLS:", ncol(df), "\n\n")

# -------------------------------------------------------------------
# notes-to-self: helpers
# -------------------------------------------------------------------
# For numeric columns: basic distribution + missingness.
# Guard rails:
#   - if everything is NA, don't let min/max turn into Inf/-Inf
summarize_numeric <- function(x) {
  n_nonmissing <- sum(!is.na(x))
  n_missing    <- sum(is.na(x))

  if (n_nonmissing == 0) {
    return(tibble(
      type        = "numeric",
      n_nonmissing = n_nonmissing,
      n_missing    = n_missing,
      min          = NA_real_,
      median       = NA_real_,
      mean         = NA_real_,
      max          = NA_real_
    ))
  }

  tibble(
    type         = "numeric",
    n_nonmissing = n_nonmissing,
    n_missing    = n_missing,
    min          = suppressWarnings(min(x, na.rm = TRUE)),
    median       = suppressWarnings(median(x, na.rm = TRUE)),
    mean         = suppressWarnings(mean(x, na.rm = TRUE)),
    max          = suppressWarnings(max(x, na.rm = TRUE))
  )
}

# For categorical-ish columns: missingness + number of distinct values
# + top few most common values.
#
# notes-to-self:
#   - treat NA as missing (don’t include in "top values")
#   - coerce to character so factors don’t surprise me later
summarize_categorical <- function(x, top_n = 5) {
  n_nonmissing <- sum(!is.na(x))
  n_missing    <- sum(is.na(x))

  x_chr <- as.character(x)

  # distinct levels, excluding NA
  n_levels <- length(unique(x_chr[!is.na(x_chr)]))

  # top values, excluding NA
  top_tbl <- sort(table(x_chr, useNA = "no"), decreasing = TRUE)

  top_values <- if (length(top_tbl) == 0) {
    NA_character_
  } else {
    paste(names(top_tbl)[seq_len(min(top_n, length(top_tbl)))], collapse = ", ")
  }

  tibble(
    type         = "categorical",
    n_nonmissing = n_nonmissing,
    n_missing    = n_missing,
    n_levels     = n_levels,
    top_values   = top_values
  )
}

# -------------------------------------------------------------------
# notes-to-self: the “important vars” list
# -------------------------------------------------------------------
# Baseline + outcome only (keep this small-ish so the CSV is readable)
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

# Only keep the ones that actually exist in df (prevents dumb crashes)
important_vars <- intersect(important_vars, names(df))

# -------------------------------------------------------------------
# build the dictionary table
# -------------------------------------------------------------------
# notes-to-self:
#   - treat integers as numeric too
#   - everything else gets the categorical summary
dictionary <- lapply(important_vars, function(v) {
  x <- df[[v]]

  out <- if (is.numeric(x) || is.integer(x)) {
    summarize_numeric(x)
  } else {
    summarize_categorical(x)
  }

  # stick variable name in front
  mutate(out, variable = v, .before = 1)
}) %>%
  bind_rows() %>%
  # nice-to-have: keep output stable and readable
  arrange(variable)

print(dictionary)

# -------------------------------------------------------------------
# save for later eyeballing
# -------------------------------------------------------------------
out_path <- "data/processed/data_dictionary.csv"
write.csv(dictionary, out_path, row.names = FALSE)

cat("\nData dictionary saved to:", out_path, "\n")
