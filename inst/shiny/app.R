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

      numericInput("interval", "Lag unit (days)", value = 7, min = 1),
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
        tabPanel("Cross-Correlation Map",
                 plotOutput("ccm_plot", height = "600px")),
        tabPanel("Figure description",
                 textOutput("messages")),
        tabPanel("Results table (preview)",
                 tableOutput("results_head")),
        tabPanel("R Code", verbatimTextOutput("rcode")),
        tabPanel("About", value = "panel4", # Tab About----
                 h3("About"),
                 tags$p("This shiny app is a friendly GUI for the ecoXCross package : Moiroux Nicolas, Bartholomée Colombine, and Taconet Paul. 2006. ecoXCorr: An R Package to Explore Lagged Associations between Environmental Time Series and Ecological Responses. doi:10.5281/zenodo.18600567"),
                  tags$a(href="https://github.com/Nmoiroux/ecoXCorr", "https://github.com/Nmoiroux/ecoXCorr"),

                 h4("References"),
                 tags$p("Curriero, Frank C., Scott M. Shone, and Gregory E. Glass. ‘Cross Correlation Maps: A Tool for Visualizing and Modeling Time Lagged Associations’. Vector Borne and Zoonotic Diseases (Larchmont, N.Y.) 5, no. 3 (2005): 267–75. https://doi.org/10.1089/vbz.2005.5.267."),
                 tags$p("Pol, Martijn van de, Liam D. Bailey, Nina McLean, Laurie Rijsdijk, Callum R. Lawson, and Lyanne Brouwer. ‘Identifying the Best Climatic Predictors in Ecology and Evolution’. Methods in Ecology and Evolution 7, no. 10 (2016): 1246–57. https://doi.org/10.1111/2041-210X.12590."),

                 )
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
      interval       = input$interval,
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


    head(res[order(abs(res[[ord]]), decreasing = TRUE), ])
  })

  # R code ----------------------------------------------------------
  output$rcode <- renderText({
    paste(
      "# Code used\n",
      "# install.packages(\"devtools\")",
      "# devtools::install_github(\"Nmoiroux/ecoXCorr\")\n",
      "library(ecoXCorr)\n",
      "results <- ecoXCorr(",
      "  meteo_data     = meteo_data,",
      "  response_data  = response_data,",
      sprintf("  date_col_meteo = %s,", input$date_meteo),
      sprintf("  date_col_resp  = %s,", input$date_resp),
      sprintf("  value_cols     = %s,", input$value_cols),
      sprintf("  agg_fun        = %s,", input$agg_fun),
      sprintf("  response       = %s,", input$response),
      sprintf("  interval       = %d,", input$interval),
      sprintf("  max_lag        = %d,", input$max_lag),
      sprintf("  random         = %s,", input$random),
      sprintf("  family         = %s", input$family),")\n",
      "plotCCM(results,",
      sprintf("  model_outcome = %s,", input$model_outcome),
      sprintf("  threshold_p   = %.2f,", input$threshold_p),")",
      sep = "\n"
    )
  })
}


shinyApp(ui, server)
