# shiny/app.R
# Grade II–III Intracranial Meningioma Survival Estimator
# Oregon Health & Science University

library(shiny)
library(survival)
library(dplyr)
library(bslib)

# Load frozen model (Shiny runs from /shiny)
fit <- readRDS("../model/cox_simple_fit.rds")

# ----------------------------
# UI
# ----------------------------
ui <- fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Roboto"),
    heading_font = font_google("Roboto")
  ),
  
  # ---- Header ----
  # ---- Header ----
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
        column(
          width = 4,
          div(style = "font-size: 22px; font-weight: bold;",
              textOutput("s1"))
        ),
        column(
          width = 4,
          div(style = "font-size: 22px; font-weight: bold;",
              textOutput("s3"))
        ),
        column(
          width = 4,
          div(style = "font-size: 22px; font-weight: bold;",
              textOutput("s5"))
        )
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

# ----------------------------
# Server
# ----------------------------
server <- function(input, output, session) {
  
  newdata <- eventReactive(input$calc, {
    data.frame(
      AGE = as.numeric(input$age),
      SEX = factor(input$sex, levels = c("Male", "Female")),
      GRADE = factor(input$grade, levels = c("II", "III")),
      CDCC_TOTAL_BEST = as.numeric(input$cdcc),
      TUMOR_SIZE_CAT = factor(
        input$size,
        levels = c("<=3cm", "3–6cm", ">6cm", "Unknown")
      ),
      RACE_CAT = factor(
        input$race,
        levels = c("White", "Black", "Other/Unknown")
      ),
      YEAR_OF_DIAGNOSIS = 2017
    )
  })
  
  get_surv <- function(sf, t) {
    idx <- max(which(sf$time <= t))
    if (length(idx) == 0) NA_real_ else sf$surv[idx]
  }
  
  output$s1 <- renderText({
    nd <- newdata()
    sf <- survfit(fit, newdata = nd)
    sprintf("1-year OS: %.1f%%", 100 * get_surv(sf, 12))
  })
  
  output$s3 <- renderText({
    nd <- newdata()
    sf <- survfit(fit, newdata = nd)
    sprintf("3-year OS: %.1f%%", 100 * get_surv(sf, 36))
  })
  
  output$s5 <- renderText({
    nd <- newdata()
    sf <- survfit(fit, newdata = nd)
    sprintf("5-year OS: %.1f%%", 100 * get_surv(sf, 60))
  })
  
  output$survplot <- renderPlot({
    nd <- newdata()
    sf <- survfit(fit, newdata = nd)
    
    plot(
      sf,
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