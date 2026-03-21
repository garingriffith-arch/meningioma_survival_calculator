cfg <- list(
  raw_csv = file.path("data", "raw", "NCDB CNS.csv"),
  
  id_col    = "PUF_CASE_ID",
  time_col  = "DX_LASTCONTACT_DEATH_MONTHS",
  event_col = "PUF_VITAL_STATUS",
  
  site_col  = "PRIMARY_SITE",
  hist_col  = "HISTOLOGY",
  grade_col = "GRADE",
  
  strict_sites     = c("C700", "C709"),
  strict_histology = 9530:9539,
  grade_keep       = c("2", "3"),
  
  time0_epsilon_months = 0.5,
  max_tumor_size_mm    = 150,
  
  horizons_months = c(12, 36, 60),
  
  age_spline_df   = 4,
  size_spline_df  = 3,
  
  age_ref_years   = 60,
  size_ref_mm     = 30,
  
  bootstrap_B = 300,
  bootstrap_seed = 1,
  
  calibration_groups = 10,
  
  processed_dir = file.path("data", "processed"),
  tables_dir    = "tables",
  figures_dir   = "figures"
)

cfg
