cox_formula <- function(cfg) {
  survival::Surv(time_months, event) ~
    splines::ns(age, df = cfg$age_spline_df) +
    splines::ns(tumor_size_mm, df = cfg$size_spline_df) +
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

safe_mode <- function(x) {
  ux <- stats::na.omit(x)
  names(sort(table(ux), decreasing = TRUE))[1]
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}

harrell_c <- function(time, event, score) {
  cc <- survival::concordance(survival::Surv(time, event) ~ score, reverse = TRUE)
  as.numeric(cc$concordance)
}

calibration_slope <- function(time, event, lp) {
  fit <- survival::coxph(survival::Surv(time, event) ~ lp)
  as.numeric(stats::coef(fit))
}

make_reference_profile <- function(df) {
  df |>
    dplyr::summarise(
      age = stats::median(age, na.rm = TRUE),
      tumor_size_mm = stats::median(tumor_size_mm, na.rm = TRUE),
      cdcc = stats::median(cdcc, na.rm = TRUE),
      grade23 = safe_mode(grade23),
      SEX = safe_mode(SEX),
      race = safe_mode(race),
      ethnicity = safe_mode(ethnicity)
    ) |>
    dplyr::mutate(
      grade23 = factor(grade23, levels = levels(df$grade23)),
      SEX = factor(SEX, levels = levels(df$SEX)),
      race = factor(race, levels = levels(df$race)),
      ethnicity = factor(ethnicity, levels = levels(df$ethnicity))
    )
}

model_matrix_from_fit <- function(fit, newdata) {
  tt <- delete.response(terms(fit))
  mf <- model.frame(tt, data = newdata, xlev = fit$xlevels, na.action = stats::na.pass)
  X <- model.matrix(tt, mf, contrasts.arg = fit$contrasts)
  
  if ("(Intercept)" %in% colnames(X)) {
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  }
  
  beta_names <- names(stats::coef(fit))
  X <- X[, beta_names, drop = FALSE]
  
  X
}

contrast_hr_ci <- function(fit, refdata, newdata) {
  X0 <- model_matrix_from_fit(fit, refdata)
  X1 <- model_matrix_from_fit(fit, newdata)
  
  dX <- X1 - X0
  beta <- stats::coef(fit)
  V <- stats::vcov(fit)
  
  lp_diff <- as.numeric(dX %*% beta)
  se_lp   <- sqrt(as.numeric(dX %*% V %*% t(dX)))
  
  tibble::tibble(
    log_hr = lp_diff,
    se_log_hr = se_lp,
    HR = exp(lp_diff),
    lo = exp(lp_diff - 1.96 * se_lp),
    hi = exp(lp_diff + 1.96 * se_lp)
  )
}

overall_and_nonlinear_p <- function(fit, cfg) {
  drop_tab <- stats::drop1(fit, test = "Chisq")
  
  overall_age_p <- drop_tab[grep("^splines::ns\\(age", rownames(drop_tab)), "Pr(>Chi)"][1]
  overall_size_p <- drop_tab[grep("^splines::ns\\(tumor_size_mm", rownames(drop_tab)), "Pr(>Chi)"][1]
  
  fit_age_linear <- stats::update(
    fit,
    . ~ . - splines::ns(age, df = cfg$age_spline_df) + age
  )
  fit_size_linear <- stats::update(
    fit,
    . ~ . - splines::ns(tumor_size_mm, df = cfg$size_spline_df) + tumor_size_mm
  )
  
  an_age  <- stats::anova(fit_age_linear, fit, test = "Chisq")
  an_size <- stats::anova(fit_size_linear, fit, test = "Chisq")
  
  nonlinear_age_p  <- an_age$`P(>|Chi|)`[2]
  nonlinear_size_p <- an_size$`P(>|Chi|)`[2]
  
  tibble::tibble(
    predictor = c("Age", "Tumor size"),
    overall_p = c(overall_age_p, overall_size_p),
    nonlinear_p = c(nonlinear_age_p, nonlinear_size_p)
  )
}

km_censor_surv <- function(time, event) {
  censor_fit <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
  times <- c(0, censor_fit$time)
  survs <- c(1, censor_fit$surv)
  
  list(
    G_t = stats::stepfun(times, c(survs, tail(survs, 1))),
    G_left = function(x) {
      idx <- findInterval(x - 1e-10, times)
      idx[idx < 1] <- 1
      survs_ext <- c(1, censor_fit$surv)
      survs_ext[idx]
    }
  )
}

ipcw_brier <- function(time, event, pred_surv, horizon) {
  kmc <- km_censor_surv(time, event)
  Gt <- kmc$G_t(horizon)
  Gt <- ifelse(Gt <= 0, NA_real_, Gt)
  
  term1 <- ifelse(
    event == 1 & time <= horizon,
    (0 - pred_surv)^2 / pmax(kmc$G_left(time), 1e-8),
    0
  )
  
  term2 <- ifelse(
    time > horizon,
    (1 - pred_surv)^2 / pmax(Gt, 1e-8),
    0
  )
  
  mean(term1 + term2, na.rm = TRUE)
}

ipcw_auc <- function(time, event, risk, horizon) {
  kmc <- km_censor_surv(time, event)
  
  case_idx <- which(event == 1 & time <= horizon)
  ctrl_idx <- which(time > horizon)
  
  if (length(case_idx) == 0 || length(ctrl_idx) == 0) {
    return(NA_real_)
  }
  
  w_case <- 1 / pmax(kmc$G_left(time[case_idx]), 1e-8)
  
  auc_case <- vapply(seq_along(case_idx), function(i) {
    r_case <- risk[case_idx[i]]
    r_ctrl <- risk[ctrl_idx]
    mean((r_case > r_ctrl) + 0.5 * (r_case == r_ctrl))
  }, numeric(1))
  
  weighted.mean(auc_case, w = w_case, na.rm = TRUE)
}

horizon_metrics <- function(fit, data, horizons) {
  surv_mat <- predict_survival_probs(fit, data, horizons)
  
  out <- lapply(seq_along(horizons), function(i) {
    h <- horizons[i]
    s_col <- surv_mat[[paste0("S_", h, "m")]]
    r_col <- 1 - s_col
    
    tibble::tibble(
      horizon_months = h,
      sample_size = nrow(data),
      deaths_by_horizon = sum(data$event == 1 & data$time_months <= h, na.rm = TRUE),
      auc = ipcw_auc(
        time = data$time_months,
        event = data$event,
        risk = r_col,
        horizon = h
      ),
      brier = ipcw_brier(
        time = data$time_months,
        event = data$event,
        pred_surv = s_col,
        horizon = h
      )
    )
  })
  
  dplyr::bind_rows(out)
}

calibration_by_group <- function(fit, data, horizon, groups = 10) {
  pred <- predict_survival_probs(fit, data, horizon)
  
  df <- data |>
    dplyr::mutate(
      pred_surv = pred[[paste0("S_", horizon, "m")]],
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
    ) |>
    dplyr::mutate(horizon = horizon)
}

make_spline_effect_df <- function(fit, df, variable, ref_value, grid = 200) {
  ref <- make_reference_profile(df)
  
  x_rng <- stats::quantile(df[[variable]], probs = c(0.01, 0.99), na.rm = TRUE)
  x_seq <- seq(x_rng[1], x_rng[2], length.out = grid)
  
  ref_profile <- ref[rep(1, length(x_seq)), , drop = FALSE]
  new_profile <- ref_profile
  
  ref_profile[[variable]] <- ref_value
  new_profile[[variable]] <- x_seq
  
  X0 <- model_matrix_from_fit(fit, ref_profile)
  X1 <- model_matrix_from_fit(fit, new_profile)
  
  dX <- X1 - X0
  beta <- stats::coef(fit)
  V <- stats::vcov(fit)
  
  lp <- as.numeric(dX %*% beta)
  se <- sqrt(diag(dX %*% V %*% t(dX)))
  
  tibble::tibble(
    x = x_seq,
    HR = exp(lp),
    lo = exp(lp - 1.96 * se),
    hi = exp(lp + 1.96 * se)
  )
}
