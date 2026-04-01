suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(splines)
  library(bslib)
  library(ggplot2)
})

fit <- readRDS(file.path("..", "data", "processed", "03_cox_model.rds"))
xlev <- fit$xlevels
model_n <- fit$n

safe_factor <- function(val, levels) {
  if (is.null(levels) || length(levels) == 0) return(factor(val))
  if (!val %in% levels) val <- levels[1]
  factor(val, levels = levels)
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}

median_stats <- function(sf) {
  out <- list(med = NA_real_, lower = NA_real_, upper = NA_real_)
  q <- tryCatch(quantile(sf, probs = 0.5, conf.int = TRUE), error = function(e) NULL)
  if (!is.null(q)) {
    out$med <- suppressWarnings(as.numeric(q$quantile[1]))
    if (!is.null(q$lower)) out$lower <- suppressWarnings(as.numeric(q$lower[1]))
    if (!is.null(q$upper)) out$upper <- suppressWarnings(as.numeric(q$upper[1]))
    return(out)
  }
  if (!is.null(sf$time) && !is.null(sf$surv) && length(sf$time) > 0) {
    idx <- which(sf$surv <= 0.5)[1]
    if (!is.na(idx)) out$med <- as.numeric(sf$time[idx])
  }
  out
}

pretty_ins <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " / ", x, fixed = TRUE)
  x <- gsub("/", " / ", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

pretty_income <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

ins_choices <- setNames(xlev$insurance_bin4, pretty_ins(xlev$insurance_bin4))
inc_choices <- setNames(xlev$income_lowhigh, pretty_income(xlev$income_lowhigh))

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    primary = "#1f4e79",
    bg = "#f4f7fb",
    fg = "#243447"
  ),

  tags$head(
    tags$style(HTML("
      :root {
        --page-max: 1320px;
        --card-radius: 24px;
        --shadow-soft: 0 8px 28px rgba(31, 52, 73, 0.07);
        --border-soft: #e7edf5;
        --text-main: #243447;
        --text-muted: #5b6b7f;
        --bg-soft: #f4f7fb;
        --accent: #1f4e79;
      }

      body {
        background: var(--bg-soft);
      }

      .app-container {
        max-width: var(--page-max);
        margin: 0 auto;
        padding: 24px 22px 36px 22px;
      }

      .app-header {
        background: #ffffff;
        border-radius: 28px;
        padding: clamp(18px, 2.2vw, 30px);
        margin-bottom: 24px;
        box-shadow: var(--shadow-soft);
        border: 1px solid var(--border-soft);
      }

      .header-grid {
        display: grid;
        grid-template-columns: minmax(70px, 96px) 1fr;
        gap: 20px;
        align-items: center;
      }

      .logo-wrap {
        display: flex;
        align-items: center;
        justify-content: center;
      }

      .ohsu-logo {
        width: clamp(58px, 6vw, 92px);
        height: auto;
        display: block;
      }

      .header-title {
        margin: 0 0 8px 0;
        font-weight: 800;
        line-height: 1.04;
        font-size: clamp(2rem, 3.7vw, 3.4rem);
        color: var(--text-main);
        max-width: 900px;
      }

      .ohsu-subtitle {
        color: var(--text-muted);
        margin: 0 0 3px 0;
        font-size: 1.05rem;
      }

      .ohsu-dept {
        color: #738396;
        margin: 0;
        font-size: 0.98rem;
      }

      .input-card, .metric-card, .plot-card, .detail-card {
        background: #ffffff;
        border: 1px solid var(--border-soft) !important;
        border-radius: var(--card-radius) !important;
        box-shadow: var(--shadow-soft);
      }

      .metric-card .card-body,
      .plot-card .card-body,
      .detail-card .card-body {
        padding: 22px;
      }

      .input-card .card-body {
        padding: 18px 18px 16px 18px;
      }

      .sticky-panel {
        position: sticky;
        top: 24px;
      }

      .section-title {
        font-weight: 800;
        color: var(--text-main);
        margin-bottom: 14px;
        line-height: 1.06;
        font-size: clamp(1.55rem, 2vw, 2rem);
      }

      .plot-title {
        font-weight: 800;
        color: var(--text-main);
        margin-bottom: 10px;
        font-size: 1.15rem;
      }

      .form-label {
        font-weight: 650;
        color: #2f4257;
        margin-bottom: 5px;
        font-size: 0.97rem;
      }

      .shiny-input-container {
        margin-bottom: 10px;
      }

      .form-control, .form-select {
        border-radius: 14px !important;
        border: 1px solid #d4dde8 !important;
        min-height: 44px;
        box-shadow: none !important;
      }

      .form-control:focus, .form-select:focus {
        border-color: #8db1d5 !important;
        box-shadow: 0 0 0 0.16rem rgba(31, 78, 121, 0.10) !important;
      }

      .btn-primary {
        background-color: #245789 !important;
        border-color: #245789 !important;
        border-radius: 14px !important;
        font-weight: 750;
        min-height: 46px;
        margin-top: 6px;
      }

      .btn-primary:hover {
        background-color: #1d476f !important;
        border-color: #1d476f !important;
      }

      .metric-grid {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: 16px;
        margin-bottom: 18px;
      }

      .metric-card {
        min-height: 126px;
      }

      .metric-value {
        font-size: clamp(1.8rem, 2.4vw, 2.45rem);
        line-height: 1;
        font-weight: 800;
        color: var(--accent);
        margin-bottom: 10px;
      }

      .metric-label {
        font-size: 0.96rem;
        color: var(--text-muted);
        line-height: 1.35;
      }

      .detail-card h3 {
        font-size: 1.08rem;
        font-weight: 750;
        color: var(--text-main);
        margin-top: 0;
        margin-bottom: 0.8rem;
      }

      .detail-card ul {
        margin-bottom: 0;
        padding-left: 1.15rem;
      }

      .detail-card li {
        color: #425466;
        margin-bottom: 0.48rem;
        line-height: 1.5;
      }

      .block-gap {
        height: 18px;
      }

      .checkbox {
        margin-top: 6px;
        margin-bottom: 0;
      }

      .plot-card .shiny-plot-output {
        margin-top: 2px;
      }

      .detail-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 26px 38px;
      }

      .detail-section {
        min-width: 0;
      }

      @media (max-width: 1199px) {
        .metric-grid {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
        .sticky-panel {
          position: static;
        }
      }

      @media (max-width: 767px) {
        .app-container {
          padding: 18px 14px 28px 14px;
        }
        .header-grid {
          grid-template-columns: 1fr;
          gap: 14px;
          text-align: center;
        }
        .logo-wrap {
          justify-content: center;
        }
        .metric-grid,
        .detail-grid {
          grid-template-columns: 1fr;
        }
        .input-card .card-body,
        .metric-card .card-body,
        .plot-card .card-body,
        .detail-card .card-body {
          padding: 18px;
        }
      }
    "))
  ),

  div(
    class = "app-container",

    div(
      class = "app-header",
      div(
        class = "header-grid",
        div(
          class = "logo-wrap",
          img(src = "ohsu_logo.png", class = "ohsu-logo")
        ),
        div(
          h1("Grade II–III Intracranial Meningioma Survival Estimator", class = "header-title"),
          p("Oregon Health & Science University", class = "ohsu-subtitle"),
          p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),

    layout_columns(
      col_widths = c(4, 8),

      div(
        class = "sticky-panel",
        card(
          class = "input-card",
          card_body(
            h2("Patient characteristics", class = "section-title"),

            numericInput("age", "Age (years)", value = 55, min = 18, max = 90, step = 1),

            selectInput(
              "sex", "Sex",
              choices = xlev$SEX,
              selected = xlev$SEX[1],
              selectize = FALSE
            ),

            selectInput(
              "race", "Race",
              choices = xlev$race,
              selected = xlev$race[1],
              selectize = FALSE
            ),

            selectInput(
              "ethnicity", "Ethnicity",
              choices = xlev$ethnicity,
              selected = xlev$ethnicity[1],
              selectize = FALSE
            ),

            numericInput("tsize_mm", "Tumor size (mm)", value = 45, min = 1, max = 150, step = 1),

            selectInput(
              "grade", "WHO grade",
              choices = c("Grade II" = "2", "Grade III" = "3"),
              selected = "2",
              selectize = FALSE
            ),

            selectInput(
              "cdcc", "Charlson-Deyo score",
              choices = 0:3,
              selected = 0,
              selectize = FALSE
            ),

            selectInput(
              "insurance", "Insurance",
              choices = ins_choices,
              selected = xlev$insurance_bin4[1],
              selectize = FALSE
            ),

            selectInput(
              "income", "Neighborhood income group",
              choices = inc_choices,
              selected = xlev$income_lowhigh[1],
              selectize = FALSE
            ),

            actionButton("calc", "Estimate survival", class = "btn-primary w-100"),
            checkboxInput("show_ci", "Show confidence bands", value = FALSE)
          )
        )
      ),

      div(
        div(
          class = "metric-grid",

          card(
            class = "metric-card",
            card_body(
              div(textOutput("s1"), class = "metric-value"),
              div("1-year overall survival", class = "metric-label")
            )
          ),

          card(
            class = "metric-card",
            card_body(
              div(textOutput("s3"), class = "metric-value"),
              div("3-year overall survival", class = "metric-label")
            )
          ),

          card(
            class = "metric-card",
            card_body(
              div(textOutput("s5"), class = "metric-value"),
              div("5-year overall survival", class = "metric-label")
            )
          ),

          card(
            class = "metric-card",
            card_body(
              div(textOutput("med"), class = "metric-value"),
              div(
                tags$span(
                  "Median predicted survival",
                  class = "metric-label",
                  title = "‘Not reached’ means the predicted survival curve does not fall below 50% within available follow-up, so the median time cannot be estimated."
                )
              )
            )
          )
        ),

        div(class = "block-gap"),

        card(
          class = "plot-card",
          card_body(
            h2("Estimated overall survival", class = "plot-title"),
            plotOutput("survplot", height = "560px")
          )
        )
      )
    ),

    div(class = "block-gap"),

    card(
      class = "detail-card",
      card_body(
        h2("Model details, analysis summary, and intended use", class = "section-title"),

        div(
          class = "detail-grid",

          div(
            class = "detail-section",
            h3("Model cohort and intended use"),
            tags$ul(
              tags$li(paste0("Model cohort: n = ", format(model_n, big.mark = ","), " patients.")),
              tags$li("Intended use: this tool provides population-level survival estimates derived from the National Cancer Database."),
              tags$li("It is intended to support clinician-patient discussion and does not replace individualized clinical judgment.")
            )
          ),

          div(
            class = "detail-section",
            h3("Cohort and variables"),
            tags$ul(
              tags$li("Data source: National Cancer Database (NCDB)."),
              tags$li("Study population: adults with intracranial meningioma limited to primary site codes C70.0 and C70.9, histology codes 9530–9539, and WHO grade II or III disease."),
              tags$li("Outcome: overall survival, measured in months from diagnosis."),
              tags$li("Predictors included in the model: age, tumor size, sex, race, ethnicity, Charlson-Deyo comorbidity score, WHO grade, insurance category, and neighborhood income group.")
            )
          ),

          div(
            class = "detail-section",
            h3("Statistical analysis"),
            tags$ul(
              tags$li("Model type: multivariable Cox proportional hazards regression."),
              tags$li("Continuous predictors were modeled flexibly using spline terms where specified in the final model."),
              tags$li("Predicted overall survival probabilities are displayed at 1, 3, and 5 years, along with the estimated survival curve."),
              tags$li("Internal validation was performed with bootstrap optimism correction and horizon-specific performance assessment.")
            )
          ),

          div(
            class = "detail-section",
            h3("Performance and interpretation"),
            tags$ul(
              tags$li("This calculator is intended to provide population-level risk estimates rather than a deterministic patient-specific prognosis."),
              tags$li("Estimates depend on the quality and coding of registry data and should be interpreted in the context of pathology, imaging, treatment details, and clinical judgment."),
              tags$li("The median predicted survival is reported only when estimable from the predicted survival function.")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  observe({
    current_age <- suppressWarnings(as.numeric(input$age))
    if (!is.na(current_age) && current_age > 90) updateNumericInput(session, "age", value = 90)
    if (!is.na(current_age) && current_age < 18) updateNumericInput(session, "age", value = 18)
  })

  observe({
    current_val <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.na(current_val) && current_val > 150) updateNumericInput(session, "tsize_mm", value = 150)
    if (!is.na(current_val) && current_val < 1) updateNumericInput(session, "tsize_mm", value = 1)
  })

  newdata <- eventReactive(input$calc, {
    age_val <- suppressWarnings(as.numeric(input$age))
    if (!is.finite(age_val)) age_val <- 55
    age_val <- min(max(age_val, 18), 90)

    tum <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.finite(tum)) tum <- 45
    tum <- min(max(tum, 1), 150)

    data.frame(
      age = age_val,
      SEX = safe_factor(input$sex, xlev$SEX),
      cdcc = as.numeric(input$cdcc),
      grade23 = safe_factor(input$grade, xlev$grade23),
      tumor_size_mm = tum,
      race = safe_factor(input$race, xlev$race),
      ethnicity = safe_factor(input$ethnicity, xlev$ethnicity),
      insurance_bin4 = safe_factor(input$insurance, xlev$insurance_bin4),
      income_lowhigh = safe_factor(input$income, xlev$income_lowhigh),
      check.names = FALSE
    )
  })

  surv_obj <- eventReactive(input$calc, {
    req(newdata())
    survfit(fit, newdata = newdata())
  })

  output$s1 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 12))
  })

  output$s3 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 36))
  })

  output$s5 <- renderText({
    req(surv_obj())
    sprintf("%.1f%%", 100 * surv_at(surv_obj(), 60))
  })

  output$med <- renderText({
    req(surv_obj())
    ms <- median_stats(surv_obj())
    if (!is.finite(ms$med)) return("Not reached")
    if (is.finite(ms$lower) && is.finite(ms$upper)) {
      pm <- as.integer(round((ms$upper - ms$lower) / 2))
      med_round <- as.integer(round(ms$med))
      if (is.finite(pm) && pm > 0) return(sprintf("%d \u00B1 %d months", med_round, pm))
    }
    sprintf("%d months", as.integer(round(ms$med)))
  })

  output$survplot <- renderPlot({
    req(surv_obj())
    sf <- surv_obj()

    df <- data.frame(time = sf$time, surv = sf$surv)
    if (!is.null(sf$lower) && !is.null(sf$upper)) {
      df$lower <- sf$lower
      df$upper <- sf$upper
    } else {
      df$lower <- NA_real_
      df$upper <- NA_real_
    }

    df <- df[df$time <= 120, , drop = FALSE]
    if (nrow(df) == 0) {
      df <- data.frame(time = c(0, 60), surv = c(1, 1), lower = c(1, 1), upper = c(1, 1))
    } else if (min(df$time) > 0) {
      df <- rbind(data.frame(time = 0, surv = 1, lower = 1, upper = 1), df)
    }

    pts <- c(12, 36, 60, 120)
    s_main <- summary(sf, times = pts, extend = TRUE)
    pts_df <- data.frame(time = s_main$time, surv = s_main$surv)
    guide_df <- data.frame(time = c(12, 36, 60, 120))

    ggplot(df, aes(x = time, y = surv)) +
      {
        if (isTRUE(input$show_ci)) {
          geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#8fb3d9", alpha = 0.22)
        }
      } +
      geom_vline(
        data = guide_df,
        aes(xintercept = time),
        linetype = "dashed",
        linewidth = 0.55,
        color = "#cbd6e2"
      ) +
      geom_step(color = "#1f6feb", linewidth = 1.5, direction = "hv") +
      geom_point(
        data = pts_df,
        aes(x = time, y = surv),
        inherit.aes = FALSE,
        color = "#1f6feb",
        size = 3.2
      ) +
      scale_x_continuous(
        limits = c(0, 120.8),
        breaks = seq(0, 120, by = 12),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      scale_y_continuous(
        limits = c(0, 1.04),
        breaks = seq(0, 1, by = 0.2),
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      labs(x = "Months", y = "Overall survival") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e7edf5", linewidth = 0.7),
        axis.title = element_text(color = "#2f4257", face = "bold"),
        axis.text = element_text(color = "#425466"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        panel.background = element_rect(fill = "#ffffff", color = NA),
        plot.margin = margin(10, 10, 8, 8)
      )
  }, res = 120)
}

shinyApp(ui, server)
