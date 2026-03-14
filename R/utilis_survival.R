cox_formula <- function() {
  survival::Surv(time_months, event) ~
    splines::ns(age, df = 4) +
    splines::ns(tumor_size_mm, df = 3) +
    SEX + cdcc + grade23 + race + ethnicity
}

predict_survival_probs <- function(fit, newdata, times) {
  sf <- survival::survfit(fit, newdata = newdata)
  
  surv_list <- lapply(seq_len(nrow(newdata)), function(i) {
    s <- summary(sf[i], times = times, extend = TRUE)$surv
    as.numeric(s)
  })
  
  out <- do.call(rbind, surv_list)
  out <- as.data.frame(out)
  names(out) <- paste0("S_", times, "m")
  out
}

harrell_c <- function(time, event, score) {
  cc <- survival::concordance(survival::Surv(time, event) ~ score, reverse = TRUE)
  as.numeric(cc$concordance)
}

calibration_slope <- function(time, event, lp) {
  fit <- survival::coxph(survival::Surv(time, event) ~ lp)
  as.numeric(stats::coef(fit))
}

calibration_by_group <- function(df, horizon, surv_col, groups = 10) {
  df <- df |>
    dplyr::mutate(
      pred_surv = .data[[surv_col]],
      pred_risk = 1 - pred_surv
    ) |>
    dplyr::filter(is.finite(pred_surv), is.finite(pred_risk)) |>
    dplyr::mutate(group = dplyr::ntile(pred_risk, groups))
  
  df |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      n = dplyr::n(),
      pred_surv_mean = mean(pred_surv, na.rm = TRUE),
      pred_risk_mean = mean(pred_risk, na.rm = TRUE),
      km_surv = {
        sf <- survival::survfit(survival::Surv(time_months, event) ~ 1, data = dplyr::cur_data())
        s <- summary(sf, times = horizon, extend = TRUE)$surv
        as.numeric(s[1])
      },
      km_risk = 1 - km_surv,
      .groups = "drop"
    )
}

safe_mode <- function(x) {
  ux <- stats::na.omit(x)
  names(sort(table(ux), decreasing = TRUE))[1]
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}
