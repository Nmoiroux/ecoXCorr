library(shiny)
library(ecoXCorr)
library(readr)

ui <- fluidPage(

  titlePanel("ecoXCorr: Lagged Cross-Correlation Analysis"),

  sidebarLayout(
    sidebarPanel(

      h4("Data source"),

      radioButtons(
        "data_source",
        "Choose data source:",
        choices = c(
          "Example datasets (package)" = "example",
          "Upload my own data" = "upload"
        )
      ),

      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("meteo_file", "Meteorological data (CSV)"),
        fileInput("resp_file", "Response data (CSV)")
      ),

      hr(),

      h4("Variables"),

      textInput("date_meteo", "Meteo date column", value = "date"),
      textInput("date_resp",  "Response date column", value = "date"),

      textInput("value_cols", "Meteorological variable", value = "rain_sum"),
      textInput("response",   "Response variable", value = "individualCount"),

      selectInput(
        "agg_fun",
        "Aggregation function",
        choices = c("mean", "sum", "min", "max"),
        selected = "sum"
      ),

      hr(),

      h4("Lag structure"),

      numericInput("lag_unit", "Lag unit (days)", value = 7, min = 1),
      numericInput("max_lag",  "Maximum lag", value = 8, min = 1),

      hr(),

      h4("Model"),

      selectInput(
        "family",
        "Model family",
        choices = c(
          "gaussian",
          "binomial",
          "poisson",
          "quasipoisson",
          "nbinom2",
          "truncated_nbinom2"
        ),
        selected = "gaussian"
      ),

      textInput(
        "random",
        "Random effects (optional)",
        value = ""
      ),

      hr(),

      h4("Plot options"),

      selectInput(
        "model_outcome",
        "Model outcome to display",
        choices = c("d_aic","R2sign","R2","betas","weight"),
        selected = "estimate"
      ),

      numericInput(
        "threshold_p",
        "Significance threshold (p-value)",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      ),


      actionButton("run", "Run ecoXCorr", class = "btn-primary")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Cross-correlation map",
                 plotOutput("ccm_plot", height = "600px")),
        tabPanel("Messages",
                 verbatimTextOutput("messages")),
        tabPanel("Results table",
                 tableOutput("results_head"))
      )
    )
  )
)



server <- function(input, output, session) {

  # Reactive datasets --------------------------------------------------

  meteo_data <- reactive({
    if (input$data_source == "example") {
      meteoMPL2023
    } else {
      req(input$meteo_file)
      read_csv(input$meteo_file$datapath)
    }
  })

  response_data <- reactive({
    if (input$data_source == "example") {
      albopictusMPL2023
    } else {
      req(input$resp_file)
      read_csv(input$resp_file$datapath)
    }
  })

  # Run ecoXCorr -------------------------------------------------------

  results <- eventReactive(input$run, {

    ecoXCorr(
      meteo_data    = meteo_data(),
      response_data = response_data(),
      date_col_meteo = input$date_meteo,
      date_col_resp  = input$date_resp,
      value_cols     = input$value_cols,
      agg_fun        = input$agg_fun,
      response       = input$response,
      lag_unit       = input$lag_unit,
      max_lag        = input$max_lag,
      random         = input$random,
      family         = input$family
    )
  })



  # Plot + message capture --------------------------------------------

  output$ccm_plot <- renderPlot({

    req(results())

    msg <- NULL
    plot <- withCallingHandlers(
      plotCCM(
        results(),
        model_outcome = input$model_outcome,
        threshold_p   = input$threshold_p
      ),
      message = function(m) {
        msg <<- conditionMessage(m)
        invokeRestart("muffleMessage")
      }
    )

    output$messages <- renderText(msg)
    plot
  })

  # Table preview ------------------------------------------------------

  output$results_head <- renderTable({
    req(results(), input$model_outcome)

    res <- results()

    if (input$model_outcome == "R2"){
      ord <- names(res)[which(grepl('^R2(?!s)',names(res), perl=TRUE) ==TRUE)]
    } else {
      ord <- input$model_outcome
    }


    head(res[order(res[[ord]], decreasing = TRUE), ])
  })
}


shinyApp(ui, server)
