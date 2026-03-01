# R/02_inspect_coding.R
# Step 2: Inspect coding for SITE, GRADE, HISTOLOGY

library(arrow)
library(dplyr)
library(stringr)

# --- Load fast Parquet copy ---
df <- read_parquet("data/processed/NCBD_Brain.parquet")

cat("\nTOTAL ROWS:", nrow(df), "\n")

# ===============================
# 1) PRIMARY SITE coding
# ===============================
cat("\n=== PRIMARY_SITE (top values) ===\n")
print(
  df %>%
    count(PRIMARY_SITE, sort = TRUE) %>%
    head(20)
)

# ===============================
# 2) GRADE coding
# ===============================
cat("\n=== GRADE (all observed values) ===\n")
print(
  df %>%
    count(GRADE, sort = TRUE)
)

# ===============================
# 3) HISTOLOGY coding (meningioma check)
# ICD-O-3 meningioma codes typically start with 953*
# ===============================
cat("\n=== HISTOLOGY values starting with 953 (meningioma family) ===\n")
print(
  df %>%
    mutate(HISTOLOGY_chr = as.character(HISTOLOGY)) %>%
    filter(str_detect(HISTOLOGY_chr, "^953")) %>%
    count(HISTOLOGY_chr, sort = TRUE)
)

# ===============================
# 4) Cross-check: SITE x HISTOLOGY
# ===============================
cat("\n=== Cross-tab: PRIMARY_SITE x meningioma histology ===\n")
print(
  df %>%
    mutate(
      PRIMARY_SITE_chr = as.character(PRIMARY_SITE),
      HISTOLOGY_chr    = as.character(HISTOLOGY)
    ) %>%
    filter(str_detect(HISTOLOGY_chr, "^953")) %>%
    count(PRIMARY_SITE_chr, sort = TRUE) %>%
    head(20)
)

cat("\nStep 2 inspection complete.\n")