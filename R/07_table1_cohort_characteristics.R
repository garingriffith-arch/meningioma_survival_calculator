# R/07_table1_cohort_characteristics.R
# Table 1: Baseline Cohort Characteristics (Grade II–III intracranial meningioma)

pkgs <- c("dplyr", "tidyr", "stringr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)

library(dplyr)
library(tidyr)
library(stringr)

# ---- Load cohort ----
df <- readRDS("data/processed/analysis_df.rds")

# ---- Recreate the exact modeling-derived categories (to stay consistent) ----
cohort <- df %>%
  mutate(
    AGE = as.numeric(AGE),
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    GRADE = factor(GRADE, levels = c(2, 3), labels = c("II", "III")),
    CDCC_TOTAL_BEST = as.numeric(CDCC_TOTAL_BEST),
    time_months = as.numeric(time_months),
    
    # Tumor size category (mm -> cm bins); 999 treated as unknown
    TUMOR_SIZE_MM = ifelse(TUMOR_SIZE >= 999, NA, as.numeric(TUMOR_SIZE)),
    TUMOR_SIZE_CAT = case_when(
      is.na(TUMOR_SIZE_MM) ~ "Unknown",
      TUMOR_SIZE_MM <= 30  ~ "<=3cm",
      TUMOR_SIZE_MM <= 60  ~ "3–6cm",
      TRUE                 ~ ">6cm"
    ),
    TUMOR_SIZE_CAT = factor(TUMOR_SIZE_CAT, levels = c("<=3cm", "3–6cm", ">6cm", "Unknown")),
    
    # Race collapsed
    RACE_CAT = case_when(
      RACE == 1 ~ "White",
      RACE == 2 ~ "Black",
      TRUE      ~ "Other/Unknown"
    ),
    RACE_CAT = factor(RACE_CAT, levels = c("White", "Black", "Other/Unknown"))
  ) %>%
  filter(
    !is.na(time_months),
    !is.na(AGE),
    !is.na(SEX),
    !is.na(GRADE),
    !is.na(CDCC_TOTAL_BEST)
  )

N <- nrow(cohort)

# ---- Helper formatters ----
fmt_n_pct <- function(x, denom) {
  n <- sum(x, na.rm = TRUE)
  pct <- 100 * n / denom
  sprintf("%d (%.1f%%)", n, pct)
}

fmt_median_iqr <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.1f (%.1f–%.1f)", q[[2]], q[[1]], q[[3]])
}

# ---- Build Table 1 rows ----
table1 <- tibble::tibble(
  Characteristic = c(
    "Patients, n",
    "Age, years, median (IQR)",
    "Sex, n (%)",
    "  Male",
    "  Female",
    "WHO Grade, n (%)",
    "  II",
    "  III",
    "Charlson–Deyo score, median (IQR)",
    "  0",
    "  1",
    "  2",
    "  3",
    "Tumor size category, n (%)",
    "  <=3 cm",
    "  3–6 cm",
    "  >6 cm",
    "  Unknown",
    "Race, n (%)",
    "  White",
    "  Black",
    "  Other/Unknown",
    "Follow-up time, months, median (IQR)"
  ),
  Value = c(
    sprintf("%d", N),
    fmt_median_iqr(cohort$AGE),
    
    "",  # section header
    fmt_n_pct(cohort$SEX == "Male", N),
    fmt_n_pct(cohort$SEX == "Female", N),
    
    "",  # section header
    fmt_n_pct(cohort$GRADE == "II", N),
    fmt_n_pct(cohort$GRADE == "III", N),
    
    fmt_median_iqr(cohort$CDCC_TOTAL_BEST),
    fmt_n_pct(cohort$CDCC_TOTAL_BEST == 0, N),
    fmt_n_pct(cohort$CDCC_TOTAL_BEST == 1, N),
    fmt_n_pct(cohort$CDCC_TOTAL_BEST == 2, N),
    fmt_n_pct(cohort$CDCC_TOTAL_BEST == 3, N),
    
    "",  # section header
    fmt_n_pct(cohort$TUMOR_SIZE_CAT == "<=3cm", N),
    fmt_n_pct(cohort$TUMOR_SIZE_CAT == "3–6cm", N),
    fmt_n_pct(cohort$TUMOR_SIZE_CAT == ">6cm", N),
    fmt_n_pct(cohort$TUMOR_SIZE_CAT == "Unknown", N),
    
    "",  # section header
    fmt_n_pct(cohort$RACE_CAT == "White", N),
    fmt_n_pct(cohort$RACE_CAT == "Black", N),
    fmt_n_pct(cohort$RACE_CAT == "Other/Unknown", N),
    
    fmt_median_iqr(cohort$time_months)
  )
)

# ---- Save outputs ----
dir.create("model", showWarnings = FALSE, recursive = TRUE)
write.csv(table1, "model/table1_cohort_characteristics.csv", row.names = FALSE)

# Optional: also write a nice Word table (if you want)
if (!"flextable" %in% rownames(installed.packages())) install.packages("flextable")
if (!"officer" %in% rownames(installed.packages())) install.packages("officer")
library(flextable)
library(officer)

ft <- flextable(table1) |>
  autofit()

doc <- read_docx() |>
  body_add_par("Table 1. Baseline Cohort Characteristics", style = "heading 1") |>
  body_add_flextable(ft)

print(doc, target = "model/table1_cohort_characteristics.docx")

cat("Saved:\n- model/table1_cohort_characteristics.csv\n- model/table1_cohort_characteristics.docx\n")