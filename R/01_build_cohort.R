suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

cfg <- source("R/config.R")$value
source("R/utils_recode.R")

dir.create(cfg$processed_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(cfg$raw_csv))

raw_dt <- fread(
  cfg$raw_csv,
  na.strings = c("", "NA", "NaN"),
  showProgress = TRUE
)

raw_df <- as_tibble(raw_dt)

flow <- tibble(step = character(), n_remaining = integer())
add_flow <- function(df, label) {
  flow <<- bind_rows(flow, tibble(step = label, n_remaining = nrow(df)))
  df
}

saveRDS(raw_df, file.path(cfg$processed_dir, "00_raw.rds"))

qc <- list(
  n_rows = nrow(raw_df),
  n_cols = ncol(raw_df),
  n_dup_id = sum(duplicated(raw_df[[cfg$id_col]])),
  n_time_missing = sum(is.na(raw_df[[cfg$time_col]])),
  n_time_zero = sum(raw_df[[cfg$time_col]] == 0, na.rm = TRUE),
  n_vital_missing = sum(is.na(raw_df[[cfg$event_col]])),
  n_site_missing = sum(is.na(raw_df[[cfg$site_col]])),
  n_hist_missing = sum(is.na(raw_df[[cfg$hist_col]])),
  n_grade_missing = sum(is.na(raw_df[[cfg$grade_col]]))
)

saveRDS(qc, file.path(cfg$processed_dir, "00_qc.rds"))

df <- raw_df |>
  add_flow("Raw NCDB extract") |>
  dplyr::filter(
    !is.na(.data[[cfg$time_col]]),
    !is.na(.data[[cfg$event_col]])
  ) |>
  add_flow("Complete survival fields") |>
  dplyr::mutate(
    primary_site_std = normalize_site(.data[[cfg$site_col]]),
    histology_num = as.integer(substr(as.character(.data[[cfg$hist_col]]), 1, 4)),
    grade_chr = trimws(as.character(.data[[cfg$grade_col]])),
    time_months = pmax(as.numeric(.data[[cfg$time_col]]), cfg$time0_epsilon_months),
    event = as.integer(.data[[cfg$event_col]] == 0),
    age = as.numeric(AGE),
    year_dx = as.integer(YEAR_OF_DIAGNOSIS),
    tumor_size_mm = as.numeric(TUMOR_SIZE_SUMMARY_16),
    SEX = factor(recode_sex(SEX), levels = c("Male", "Female")),
    cdcc = as.integer(CDCC_TOTAL_BEST),
    grade23 = factor(grade_chr, levels = c("2", "3")),
    race = factor(recode_race_wbo(RACE), levels = c("White", "Black", "Other")),
    ethnicity = factor(
      recode_ethnicity_binary(SPANISH_HISPANIC_ORIGIN),
      levels = c("Non-Hispanic", "Hispanic")
    )
  ) |>
  dplyr::filter(primary_site_std %in% cfg$strict_sites) |>
  add_flow("Primary site: intracranial meninges") |>
  dplyr::filter(histology_num %in% cfg$strict_histology) |>
  add_flow("Histology: meningioma 9530-9539") |>
  dplyr::filter(grade_chr %in% cfg$grade_keep) |>
  add_flow("Grade II-III") |>
  dplyr::filter(age >= 18) |>
  add_flow("Age >=18 years") |>
  dplyr::filter(
    !is.na(tumor_size_mm),
    tumor_size_mm > 0,
    tumor_size_mm <= cfg$max_tumor_size_mm
  ) |>
  add_flow("Known tumor size within plausible range") |>
  tidyr::drop_na(
    time_months, event, age, year_dx, SEX, cdcc, grade23, tumor_size_mm, race, ethnicity
  ) |>
  add_flow("Complete-case analytic cohort")

saveRDS(df, file.path(cfg$processed_dir, "01_analytic_cohort.rds"))
write_csv(flow, file.path(cfg$processed_dir, "01_cohort_flow.csv"))
