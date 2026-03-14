suppressPackageStartupMessages({
  library(shiny)
  library(survival)
  library(splines)
  library(bslib)
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

ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Roboto"),
    heading_font = font_google("Roboto")
  ),
  
  fluidRow(
    style = "
      background-color: #f8f9fa;
      padding: 20px 15px;
      border-bottom: 1px solid #ddd;
      align-items: center;
      display: flex;
    ",
    column(
      width = 3,
      img(
        src = "ohsu_logo.png",
        style = "height: 100px; width: auto; display: block;"
      )
    ),
    column(
      width = 9,
      div(
        style = "display: flex; flex-direction: column; justify-content: center;",
        h3("Grade II-III Intracranial Meningioma Survival Estimator", style = "margin-bottom: 5px;"),
        h5("Oregon Health & Science University", style = "color: #555; margin-top: 0;"),
        h6("Department of Neurological Surgery", style = "color: #666; margin-top: 2px;")
      )
    )
  ),
  
  div(style = "height: 12px;"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient Characteristics"),
      
      numericInput("age", "Age (years)", value = 55, min = 18, max = 90),
      
      selectInput("sex", "Sex", choices = xlev$SEX, selected = xlev$SEX[1]),
      
      selectInput(
        "race", "Race",
        choices = xlev$race,
        selected = xlev$race[1]
      ),
      
      selectInput(
        "ethnicity", "Ethnicity",
        choices = xlev$ethnicity,
        selected = xlev$ethnicity[1]
      ),
      
      numericInput("tsize_mm", "Tumor size (mm)", value = 45, min = 1, max = 150),
      
      selectInput(
        "grade", "WHO Grade",
        choices = c("Grade II" = "2", "Grade III" = "3"),
        selected = "2"
      ),
      
      numericInput("cdcc", "Charlson-Deyo score", value = 0, min = 0, max = 3),
      
      hr(),
      
      actionButton("calc", "Estimate Survival", class = "btn-primary", style = "width: 100%;"),
      checkboxInput("show_ci", "Show confidence bands", value = FALSE)
    ),
    
    mainPanel(
      width = 9,
      h4("Estimated Overall Survival"),
      
      fluidRow(
        column(4, h5(textOutput("s1"))),
        column(4, h5(textOutput("s3"))),
        column(4, h5(textOutput("s5")))
      ),
      
      plotOutput("survplot", height = "420px"),
      
      br(),
      
      wellPanel(
        p(
          strong("Intended use: "),
          "This tool provides population-level survival estimates derived from the National Cancer Database. ",
          "It is intended to support clinician-patient discussion and does not replace individualized clinical judgment."
        )
      )
    )
  )
)

server <- function(input, output, session) {
  newdata <- eventReactive(input$calc, {
    tum <- as.numeric(input$tsize_mm)
    if (!is.finite(tum) || tum <= 0 || tum > 150) {
      tum <- 45
    }
    
    data.frame(
      age = as.numeric(input$age),
      SEX = safe_factor(input$sex, xlev$SEX),
      cdcc = as.numeric(input$cdcc),
      grade23 = safe_factor(input$grade, xlev$grade23),
      tumor_size_mm = tum,
      race = safe_factor(input$race, xlev$race),
      ethnicity = safe_factor(input$ethnicity, xlev$ethnicity),
      check.names = FALSE
    )
  })
  
  output$s1 <- renderText({
    req(newdata())
    sf <- survfit(fit, newdata = newdata())
    sprintf("1-year OS: %.1f%%", 100 * surv_at(sf, 12))
  })
  
  output$s3 <- renderText({
    req(newdata())
    sf <- survfit(fit, newdata = newdata())
    sprintf("3-year OS: %.1f%%", 100 * surv_at(sf, 36))
  })
  
  output$s5 <- renderText({
    req(newdata())
    sf <- survfit(fit, newdata = newdata())
    sprintf("5-year OS: %.1f%%", 100 * surv_at(sf, 60))
  })
  
  output$survplot <- renderPlot({
    req(newdata())
    sf <- survfit(fit, newdata = newdata())
    
    plot(
      sf,
      xlab = "Months",
      ylab = "Overall survival",
      xlim = c(0, 60),
      ylim = c(0, 1),
      lwd = 3,
      col = "#0d6efd",
      conf.int = FALSE,
      mark.time = FALSE
    )
    grid()
    
    pts <- c(12, 36, 60)
    s_main <- summary(sf, times = pts, extend = TRUE)
    
    points(s_main$time, s_main$surv, pch = 16, cex = 0.9, col = "#0d6efd")
    
    if (isTRUE(input$show_ci)) {
      points(s_main$time, s_main$lower, pch = 15, cex = 0.85, col = "#6c757d")
      points(s_main$time, s_main$upper, pch = 15, cex = 0.85, col = "#6c757d")
      segments(s_main$time, s_main$lower, s_main$time, s_main$upper, col = "#6c757d", lwd = 1)
    }
  })
}

shinyApp(ui, server)
