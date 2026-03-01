# R/01_prepare_data.R
# Step 1: load Excel -> save Parquet -> basic cohort filter stub

# --- Packages ---
pkgs <- c("readxl", "dplyr", "stringr", "arrow")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)

library(readxl)
library(dplyr)
library(stringr)
library(arrow)

# --- Paths (run from the PROJECT ROOT in RStudio) ---
infile <- "data/raw/NCBD Brain.xlsx"

# --- 0) Confirm file exists ---
stopifnot(file.exists(infile))

# --- 1) Identify sheets ---
sheets <- excel_sheets(infile)
print(sheets)

# --- 2) Read FIRST sheet (change sheets[1] to another index if needed) ---
# .name_repair ensures column names are unique even if Excel has duplicates
df <- read_excel(infile, sheet = sheets[1], .name_repair = "unique")

# --- 3) Quick sanity checks ---
cat("\nRows:", nrow(df), " | Cols:", ncol(df), "\n\n")
print(names(df)[1:min(30, length(names(df)))])  # show first 30 column names

# --- 4) Save a fast, analysis-friendly copy (Parquet) ---
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write_parquet(df, "data/processed/NCBD_Brain.parquet")

# Optional: Save a small preview CSV (first 200 rows) so you can inspect quickly
write.csv(head(df, 200), "data/processed/peek_200rows.csv", row.names = FALSE)

# --- 5) First-pass cohort filter (stub) ---
# NOTE: This is a generic filter. We'll tighten coding once we inspect values.
cohort_stub <- df %>%
  mutate(
    PRIMARY_SITE_chr = as.character(.data$PRIMARY_SITE),
    GRADE_chr        = as.character(.data$GRADE)
  ) %>%
  filter(
    str_detect(PRIMARY_SITE_chr, "^C7(0|1)"),
    GRADE_chr %in% c("2", "3", "II", "III")
  )

saveRDS(cohort_stub, "data/processed/cohort_stub.rds")

cat("\nSaved:\n",
    "- data/processed/NCBD_Brain.parquet\n",
    "- data/processed/peek_200rows.csv\n",
    "- data/processed/cohort_stub.rds\n",
    "\nCohort stub rows:", nrow(cohort_stub), "\n")

print("finished with 01_prepare_data.R")