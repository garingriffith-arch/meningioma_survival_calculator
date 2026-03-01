# R/03_build_analysis_df.R
# Build final analysis dataset with survival outcome
# This is the dataset used for all downstream modeling.

library(arrow)
library(dplyr)
library(stringr)

# ---- load processed data ----
df <- read_parquet("data/processed/NCBD_Brain.parquet")

# ---- construct analysis dataset ----
analysis_df <- df %>%
  mutate(
    # primary site as character for regex filtering
    PRIMARY_SITE_chr = as.character(PRIMARY_SITE),

    # survival time (months)
    time_months = as.numeric(DX_LASTCONTACT_DEATH_MONTHS),

    # event indicator
    event = case_when(
      PUF_VITAL_STATUS %in% c(1, "1", "Dead", "DEAD")   ~ 1,
      PUF_VITAL_STATUS %in% c(0, "0", "Alive", "ALIVE") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    # intracranial sites
    str_detect(PRIMARY_SITE_chr, "^C7(0|1)"),

    # grade II–III only
    GRADE %in% c(2, 3),

    # valid survival data
    !is.na(time_months),
    !is.na(event),
    time_months >= 0
  )

# ---- save analysis dataset ----
saveRDS(analysis_df, "data/processed/analysis_df.rds")

message("Final analysis dataset saved:")
message(" - data/processed/analysis_df.rds")
message("Rows: ", nrow(analysis_df))

invisible(TRUE)
