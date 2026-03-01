# R/01_prepare_data.R
# Prep: read Excel -> write Parquet -> quick first-pass cohort stub (intracranial + grade II/III)

# ---- packages ----
pkgs <- c("readxl", "dplyr", "stringr", "arrow")
need <- setdiff(pkgs, rownames(installed.packages()))
if (length(need)) install.packages(need)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(arrow)
})

# ---- paths ----
# run this from project root (meningioma-survival-calculator/)
infile <- file.path("data", "raw", "NCBD Brain.xlsx")

# ---- sanity check ----
stopifnot(file.exists(infile))

# ---- sheet discovery ----
sheets <- excel_sheets(infile)
message("Sheets found: ", paste(sheets, collapse = ", "))

# NOTE to self: if the wrong sheet loads, change sheets[1] to the right one
df <- read_excel(infile, sheet = sheets[1], .name_repair = "unique")

# ---- quick peek ----
message("Rows: ", nrow(df), " | Cols: ", ncol(df))
print(names(df)[1:min(30, length(names(df)))])

# ---- write processed outputs ----
dir.create(file.path("data", "processed"), showWarnings = FALSE, recursive = TRUE)

write_parquet(df, file.path("data", "processed", "NCBD_Brain.parquet"))
write.csv(head(df, 200),
          file.path("data", "processed", "peek_200rows.csv"),
          row.names = FALSE)

# ---- cohort stub (just to confirm we have expected rows) ----
# NOTE to self: this is intentionally loose; tighten once we verify coding in 02
cohort_stub <- df %>%
  mutate(
    PRIMARY_SITE_chr = as.character(PRIMARY_SITE),
    GRADE_chr = as.character(GRADE)
  ) %>%
  filter(
    str_detect(PRIMARY_SITE_chr, "^C7(0|1)"),
    GRADE_chr %in% c("2", "3", "II", "III")
  )

saveRDS(cohort_stub, file.path("data", "processed", "cohort_stub.rds"))

message("Saved outputs:")
message(" - data/processed/NCBD_Brain.parquet")
message(" - data/processed/peek_200rows.csv")
message(" - data/processed/cohort_stub.rds")
message("Cohort stub rows: ", nrow(cohort_stub))

invisible(TRUE)
