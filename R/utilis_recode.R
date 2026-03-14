normalize_site <- function(x) {
  gsub("\\.", "", toupper(trimws(as.character(x))))
}

recode_sex <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x == "1" ~ "Male",
    x == "2" ~ "Female",
    TRUE ~ NA_character_
  )
}

recode_race_wbo <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x == "1" ~ "White",
    x == "2" ~ "Black",
    x %in% c("98", "99", "", " ", NA) ~ NA_character_,
    TRUE ~ "Other"
  )
}

recode_ethnicity_binary <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x == "0" ~ "Non-Hispanic",
    x %in% c("1", "2", "3", "4", "5", "6", "7", "8") ~ "Hispanic",
    TRUE ~ NA_character_
  )
}
