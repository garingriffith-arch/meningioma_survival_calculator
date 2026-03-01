# R/08_table2_multivariable_cox.R
library(survival)
library(dplyr)
library(tibble)

fit <- readRDS("model/cox_simple_fit.rds")
s <- summary(fit)

table2 <- tibble(
  Variable = rownames(s$coef),
  HR = s$conf.int[, "exp(coef)"],
  CI_lower = s$conf.int[, "lower .95"],
  CI_upper = s$conf.int[, "upper .95"],
  p_value = s$coef[, "Pr(>|z|)"]
)

# ---- Rescale YEAR_OF_DIAGNOSIS to per 5-year increase (presentation only) ----
yr_row <- table2$Variable == "YEAR_OF_DIAGNOSIS"
table2$HR[yr_row]       <- table2$HR[yr_row]^5
table2$CI_lower[yr_row] <- table2$CI_lower[yr_row]^5
table2$CI_upper[yr_row] <- table2$CI_upper[yr_row]^5
table2$Variable[yr_row] <- "Year of diagnosis (per 5-year increase)"

# ---- Drop reference-level race row (presentation fix) ----
table2 <- table2 %>%
  filter(Variable != "RACE_CATWhite")

# Drop year of diagnosis from Table 2 (adjustment variable only)
table2 <- table2 %>%
  filter(!grepl("^Year of diagnosis", Variable))

# ---- Clean labels ----
clean_label <- function(x) {
  x <- gsub("^AGE$", "Age (per year)", x)
  x <- gsub("^SEXFemale$", "Sex: Female vs Male", x)
  x <- gsub("^GRADEIII$", "WHO Grade III vs II", x)
  x <- gsub("^CDCC_TOTAL_BEST$", "Charlson–Deyo score (per point)", x)
  x <- gsub("^TUMOR_SIZE_CAT3–6cm$", "Tumor size 3–6 cm vs ≤3 cm", x)
  x <- gsub("^TUMOR_SIZE_CAT>6cm$", "Tumor size >6 cm vs ≤3 cm", x)
  x <- gsub("^TUMOR_SIZE_CATUnknown$", "Tumor size unknown vs ≤3 cm", x)
  x <- gsub("^RACE_CATOther/Unknown$", "Race: Other/Unknown vs White", x)
  x
}

table2$Variable <- clean_label(table2$Variable)

# ---- Format numbers for publication ----
table2 <- table2 %>%
  mutate(
    HR = round(HR, 3),
    CI_lower = round(CI_lower, 3),
    CI_upper = round(CI_upper, 3),
    p_value = ifelse(p_value < 0.001, "<0.001", signif(p_value, 3))
  )

write.csv(table2, "model/table2_multivariable_cox.csv", row.names = FALSE)
print(table2)