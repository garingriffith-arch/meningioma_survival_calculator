# R/02_inspect_coding.R
# Inspect coding for PRIMARY_SITE, GRADE, and HISTOLOGY
# Goal here is just to understand how things are coded before filtering hard.

library(arrow)
library(dplyr)
library(stringr)

# ---- load processed data ----
df <- read_parquet("data/processed/NCBD_Brain.parquet")

cat("\nTOTAL ROWS:", nrow(df), "\n")

# ===============================
# 1) PRIMARY_SITE coding
# ===============================
cat("\n=== PRIMARY_SITE: most common values ===\n")
df %>%
  count(PRIMARY_SITE, sort = TRUE) %>%
  head(20) %>%
  print()

# ===============================
# 2) GRADE coding
# ===============================
cat("\n=== GRADE: all observed values ===\n")
df %>%
  count(GRADE, sort = TRUE) %>%
  print()

# ===============================
# 3) HISTOLOGY coding (meningioma check)
# ICD-O-3 meningioma histologies typically start with 953*
# ===============================
cat("\n=== HISTOLOGY codes starting with 953 (meningioma family) ===\n")
df %>%
  mutate(HISTOLOGY_chr = as.character(HISTOLOGY)) %>%
  filter(str_detect(HISTOLOGY_chr, "^953")) %>%
  count(HISTOLOGY_chr, sort = TRUE) %>%
  print()

# ===============================
# 4) Cross-check: PRIMARY_SITE x HISTOLOGY
# Make sure meningioma histologies are showing up where we expect
# ===============================
cat("\n=== PRIMARY_SITE distribution among 953* histologies ===\n")
df %>%
  mutate(
    PRIMARY_SITE_chr = as.character(PRIMARY_SITE),
    HISTOLOGY_chr    = as.character(HISTOLOGY)
  ) %>%
  filter(str_detect(HISTOLOGY_chr, "^953")) %>%
  count(PRIMARY_SITE_chr, sort = TRUE) %>%
  head(20) %>%
  print()

cat("\nStep 2 inspection complete.\n")
invisible(TRUE)
