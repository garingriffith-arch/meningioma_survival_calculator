suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(splines)
  library(bslib)
  library(ggplot2)
})

fit <- readRDS(file.path("..", "data", "processed", "03_cox_model.rds"))
xlev <- fit$xlevels

safe_factor <- function(val, levels) {
  if (is.null(levels) || length(levels) == 0) {
    return(factor(val))
  }
  if (!val %in% levels) {
    val <- levels[1]
  }
  factor(val, levels = levels)
}

surv_at <- function(sf, t) {
  s <- summary(sf, times = t, extend = TRUE)$surv
  as.numeric(s[1])
}

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
      .app-container {
        max-width: 1400px;
        margin: 0 auto;
        padding: 18px 18px 28px 18px;
      }

      .app-header {
        background: #ffffff;
        border-radius: 22px;
        padding: clamp(16px, 2vw, 26px);
        margin-bottom: 20px;
        box-shadow: 0 6px 24px rgba(31, 52, 73, 0.08);
      }

      .header-grid {
        display: grid;
        grid-template-columns: minmax(70px, 10vw) 1fr;
        gap: clamp(14px, 2vw, 26px);
        align-items: center;
      }

      .logo-wrap {
        display: flex;
        align-items: center;
        justify-content: center;
      }

      .ohsu-logo {
        width: clamp(58px, 7vw, 110px);
        height: auto;
        display: block;
      }

      .header-title {
        margin: 0 0 6px 0;
        font-weight: 800;
        line-height: 1.08;
        font-size: clamp(1.9rem, 3.8vw, 3.7rem);
        color: #243447;
      }

      .ohsu-subtitle {
        color: #5b6b7f;
        margin: 0;
        font-size: clamp(0.98rem, 1.4vw, 1.15rem);
      }

      .ohsu-dept {
        color: #738396;
        margin: 0;
        font-size: clamp(0.92rem, 1.2vw, 1.02rem);
      }

      .input-card, .metric-card, .plot-card, .info-card {
        background: #ffffff;
        border: 1px solid #e7edf5 !important;
        border-radius: 22px !important;
        box-shadow: 0 6px 24px rgba(31, 52, 73, 0.06);
      }

      .section-title {
        font-weight: 700;
        color: #243447;
        margin-bottom: 14px;
        line-height: 1.05;
      }

      .metric-card {
        min-height: 138px;
      }

      .metric-value {
        font-size: clamp(2rem, 2.8vw, 2.6rem);
        line-height: 1;
        font-weight: 800;
        color: #1f4e79;
        margin-bottom: 10px;
      }

      .metric-label {
        font-size: 0.98rem;
        color: #5b6b7f;
      }

      .plot-title {
        font-weight: 700;
        color: #243447;
        margin-bottom: 12px;
      }

      .form-label {
        font-weight: 600;
        color: #2f4257;
        margin-bottom: 6px;
      }

      .shiny-input-container {
        margin-bottom: 14px;
      }

      .form-control, .form-select {
        border-radius: 12px !important;
        border: 1px solid #d4dde8 !important;
        min-height: 46px;
      }

      .btn-primary {
        background-color: #245789 !important;
        border-color: #245789 !important;
        border-radius: 14px !important;
        font-weight: 700;
        min-height: 46px;
      }

      .btn-primary:hover {
        background-color: #1d476f !important;
        border-color: #1d476f !important;
      }

      .info-text {
        font-size: 0.98rem;
        color: #425466;
        margin-bottom: 0;
      }

      .info-label {
        font-weight: 700;
        color: #243447;
      }

      .checkbox {
        margin-top: 4px;
      }

      @media (max-width: 768px) {
        .section-title {
          font-size: 1.8rem;
        }

        .metric-card {
          min-height: 118px;
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
          img(
            src = "ohsu_logo.png",
            class = "ohsu-logo"
          )
        ),
        div(
          h1(
            "Grade II–III Intracranial Meningioma Survival Estimator",
            class = "header-title"
          ),
          p("Oregon Health & Science University", class = "ohsu-subtitle"),
          p("Department of Neurological Surgery", class = "ohsu-dept")
        )
      )
    ),
    
    layout_columns(
      col_widths = c(4, 8),
      
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
          
          numericInput(
            "tsize_mm",
            "Tumor size (mm)",
            value = 45,
            min = 1,
            max = 150,
            step = 1
          ),
          
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
          
          div(style = "margin-top: 10px;"),
          
          actionButton(
            "calc",
            "Estimate survival",
            class = "btn-primary w-100"
          ),
          
          checkboxInput("show_ci", "Show confidence bands", value = FALSE)
        )
      ),
      
      div(
        layout_columns(
          col_widths = c(4, 4, 4),
          
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
          )
        ),
        
        div(style = "height: 16px;"),
        
        card(
          class = "plot-card",
          card_body(
            h2("Estimated overall survival", class = "plot-title"),
            plotOutput("survplot", height = "450px")
          )
        ),
        
        div(style = "height: 16px;"),
        
        card(
          class = "info-card",
          card_body(
            p(
              class = "info-text",
              span("Intended use: ", class = "info-label"),
              "This tool provides population-level survival estimates derived from the National Cancer Database. It is intended to support clinician-patient discussion and does not replace individualized clinical judgment."
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
    if (!is.na(current_age) && current_age > 90) {
      updateNumericInput(session, "age", value = 90)
    }
    if (!is.na(current_age) && current_age < 18) {
      updateNumericInput(session, "age", value = 18)
    }
  })
  
  observe({
    current_val <- suppressWarnings(as.numeric(input$tsize_mm))
    if (!is.na(current_val) && current_val > 150) {
      updateNumericInput(session, "tsize_mm", value = 150)
    }
    if (!is.na(current_val) && current_val < 1) {
      updateNumericInput(session, "tsize_mm", value = 1)
    }
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
  
  output$survplot <- renderPlot({
    req(surv_obj())
    sf <- surv_obj()
    
    df <- data.frame(
      time = sf$time,
      surv = sf$surv
    )
    
    if (!is.null(sf$lower) && !is.null(sf$upper)) {
      df$lower <- sf$lower
      df$upper <- sf$upper
    } else {
      df$lower <- NA_real_
      df$upper <- NA_real_
    }
    
    df <- df[df$time <= 60, , drop = FALSE]
    
    if (nrow(df) == 0) {
      df <- data.frame(
        time = c(0, 60),
        surv = c(1, 1),
        lower = c(1, 1),
        upper = c(1, 1)
      )
    } else if (min(df$time) > 0) {
      df <- rbind(
        data.frame(time = 0, surv = 1, lower = 1, upper = 1),
        df
      )
    }
    
    pts <- c(12, 36, 60)
    s_main <- summary(sf, times = pts, extend = TRUE)
    
    pts_df <- data.frame(
      time = s_main$time,
      surv = s_main$surv
    )
    
    ann_df <- data.frame(
      time = c(12, 36, 60),
      label = c("1-year", "3-year", "5-year")
    )
    
    p <- ggplot(df, aes(x = time, y = surv)) +
      {
        if (isTRUE(input$show_ci)) {
          geom_ribbon(
            aes(ymin = lower, ymax = upper),
            fill = "#8fb3d9",
            alpha = 0.22
          )
        }
      } +
      geom_vline(
        data = ann_df,
        aes(xintercept = time),
        linetype = "dashed",
        linewidth = 0.55,
        color = "#cbd6e2"
      ) +
      geom_step(
        color = "#1f6feb",
        linewidth = 1.5,
        direction = "hv"
      ) +
      geom_point(
        data = pts_df,
        aes(x = time, y = surv),
        inherit.aes = FALSE,
        color = "#1f6feb",
        size = 3.2
      ) +
      annotate(
        "text",
        x = c(12, 36, 58.5),
        y = c(1.015, 1.015, 1.015),
        label = c("1-year", "3-year", "5-year"),
        size = 3.8,
        fontface = "bold",
        color = "#5b6b7f",
        vjust = 0
      ) +
      scale_x_continuous(
        limits = c(0, 60.8),
        breaks = seq(0, 60, by = 12),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      scale_y_continuous(
        limits = c(0, 1.04),
        breaks = seq(0, 1, by = 0.2),
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      labs(
        x = "Months",
        y = "Overall survival"
      ) +
      coord_cartesian(clip = "off") +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e7edf5", linewidth = 0.7),
        axis.title = element_text(color = "#2f4257", face = "bold"),
        axis.text = element_text(color = "#425466"),
        plot.background = element_rect(fill = "#ffffff", color = NA),
        panel.background = element_rect(fill = "#ffffff", color = NA),
        plot.margin = margin(18, 30, 10, 10)
      )
    
    p
  }, res = 120)
}

shinyApp(ui, server)
