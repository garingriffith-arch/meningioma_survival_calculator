suppressPackageStartupMessages({
  library(dplyr)
})

cfg <- source("R/config.R")$value
df <- readRDS(file.path(cfg$processed_dir, "01_analytic_cohort.rds"))

years <- sort(unique(df$year_dx))
stopifnot(length(years) >= 1)

if (length(years) < 3) {
  set.seed(1)
  df <- df |>
    mutate(split = ifelse(runif(n()) < 0.8, "train", "test"))
  message("Fewer than 3 diagnosis years available; using random 80/20 split.")
} else {
  n_total <- nrow(df)
  target_test_n <- max(100L, floor(0.20 * n_total))
  
  n_by_year <- df |>
    count(year_dx, name = "n") |>
    arrange(year_dx)
  
  k <- 1L
  while (k < length(years)) {
    test_years <- tail(years, k)
    n_test <- sum(n_by_year$n[n_by_year$year_dx %in% test_years])
    if (n_test >= target_test_n && k < length(years)) break
    k <- k + 1L
  }
  
  test_years <- tail(years, k)
  cut_year <- min(test_years)
  
  df <- df |>
    mutate(split = ifelse(year_dx >= cut_year, "test", "train"))
  
  message(
    "Temporal split used. Test years: ",
    paste(test_years, collapse = ", "),
    " | n_test = ",
    sum(df$split == "test")
  )
}

tab <- table(df$split)
print(tab)
stopifnot(all(c("train", "test") %in% names(tab)))

saveRDS(df, file.path(cfg$processed_dir, "02_analytic_cohort_split.rds"))
