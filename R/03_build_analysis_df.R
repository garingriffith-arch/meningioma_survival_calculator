# R/03_build_analysis_df.R
# Step 3: Final cohort + survival outcome

library(arrow)
library(dplyr)
library(stringr)

df <- read_parquet("data/processed/NCBD_Brain.parquet")

analysis_df <- df %>%
  mutate(
    PRIMARY_SITE_chr = as.character(PRIMARY_SITE),
    time_months = as.numeric(DX_LASTCONTACT_DEATH_MONTHS),
    event = case_when(
      PUF_VITAL_STATUS %in% c(1, "1", "Dead", "DEAD") ~ 1,
      PUF_VITAL_STATUS %in% c(0, "0", "Alive", "ALIVE") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    str_detect(PRIMARY_SITE_chr, "^C7(0|1)"),
    GRADE %in% c(2, 3),
    !is.na(time_months),
    !is.na(event),
    time_months >= 0
  )

saveRDS(analysis_df, "data/processed/analysis_df.rds")

cat(
  "Final analysis dataset saved\n",
  "Rows:", nrow(analysis_df), "\n"
)