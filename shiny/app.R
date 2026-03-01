# shiny/app.R
# Grade II–III Intracranial Meningioma Survival Estimator
# Oregon Health & Science University

library(shiny)
library(survival)
library(dplyr)
library(bslib)

# ------------------------------------------------------------
# load frozen model (built in 05)
# ------------------------------------------------------------
fit <- readRDS("../model/cox_simple_fit.rds")

# ------------------------------------------------------------
# UI
# ------------------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Roboto"),
    heading_font = font_google("Roboto")
  ),

  # header strip (logo + title)
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
        h3(
          "Grade II–III Intracranial Meningioma Survival Estimator",
          style = "margin-bottom: 5px;"
        ),
        h5(
          "Oregon Health & Science University",
          style = "color: #555; margin-top: 0;"
        )
      )
    )
  ),

  div(style = "height: 10px;"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient Characteristics"),

      numericInput("age", "Age (years)", value = 55, min = 18, max = 90),

      selectInput("sex", "Sex", choices = c("Male", "Female")),

      selectInput("grade", "WHO Grade", choices = c("II", "III")),

      numericInput(
        "cdcc",
        "Charlson–Deyo comorbidity score",
        value = 0, min = 0, max = 3
      ),

      selectInput(
        "size",
        "Tumor size category",
        choices = c("<=3cm", "3–6cm", ">6cm", "Unknown")
      ),

      selectInput(
        "race",
        "Race",
        choices = c("White", "Black", "Other/Unknown")
      ),

      actionButton(
        "calc",
        "Estimate Survival",
        class = "btn-primary",
        style = "width: 100%; margin-top: 15px;"
      )
    ),

    mainPanel(
      width = 9,

      h4("Estimated Overall Survival"),

      fluidRow(
        column(4, div(style = "font-size: 22px; font-weight: bold;", textOutput("s1"))),
        column(4, div(style = "font-size: 22px; font-weight: bold;", textOutput("s3"))),
        column(4, div(style = "font-size: 22px; font-weight: bold;", textOutput("s5")))
      ),

      br(),
      plotOutput("survplot", height = "420px"),

      br(),
      wellPanel(
        p(
          strong("Intended use: "),
          "This tool estimates overall survival based on population-level registry data ",
          "and is intended to support clinical discussion. Predictions assume contemporary ",
          "management approximated by the most recent era available in the dataset and are ",
          "not a substitute for individualized clinical judgment."
        )
      )
    )
  )
)

# ------------------------------------------------------------
# Server
# ------------------------------------------------------------
server <- function(input, output, session) {

  # fixed model-era assumption (remove from UI; keep consistent)
  YEAR_FIXED <- 2017

  # notes-to-self: stable single-time survival extraction
  surv_at <- function(sf, t_months) {
    s <- summary(sf, times = t_months)
    if (length(s$surv) == 0) return(NA_real_)
    as.numeric(s$surv)
  }

  # one click => build newdata + run survfit once + cache outputs
  res <- eventReactive(input$calc, {

    nd <- data.frame(
      AGE = as.numeric(input$age),
      SEX = factor(input$sex, levels = c("Male", "Female")),
      GRADE = factor(input$grade, levels = c("II", "III")),
      CDCC_TOTAL_BEST = as.numeric(input$cdcc),
      TUMOR_SIZE_CAT = factor(input$size, levels = c("<=3cm", "3–6cm", ">6cm", "Unknown")),
      RACE_CAT = factor(input$race, levels = c("White", "Black", "Other/Unknown")),
      YEAR_OF_DIAGNOSIS = YEAR_FIXED
    )

    sf <- survfit(fit, newdata = nd)

    list(
      nd = nd,
      sf = sf,
      s12 = surv_at(sf, 12),
      s36 = surv_at(sf, 36),
      s60 = surv_at(sf, 60)
    )
  }, ignoreInit = TRUE)

  output$s1 <- renderText({
    r <- res(); req(r)
    sprintf("1-year OS: %.1f%%", 100 * r$s12)
  })

  output$s3 <- renderText({
    r <- res(); req(r)
    sprintf("3-year OS: %.1f%%", 100 * r$s36)
  })

  output$s5 <- renderText({
    r <- res(); req(r)
    sprintf("5-year OS: %.1f%%", 100 * r$s60)
  })

  output$survplot <- renderPlot({
    r <- res(); req(r)

    plot(
      r$sf,
      xlab = "Months since diagnosis",
      ylab = "Overall survival probability",
      lwd = 3,
      col = "#003B5C",
      conf.int = FALSE
    )
    grid()
  })
}

shinyApp(ui, server)
