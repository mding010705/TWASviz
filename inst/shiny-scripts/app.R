library(shiny)
library(data.table)
library(TWASviz)

ui <- fluidPage(
  titlePanel("Sparse-group lasso Transcriptome Wide Association Analysis (TWAS)"),
  sidebarLayout(
    sidebarPanel(
      tags$p("This is a Shiny App that is part of TWASviz in R."),
      br(),

      tags$b("Description: ."),

      # br() element to introduce extra vertical spacing ----
      br(),
      br(),

      # input
      tags$p("Instructions: Upload one or several PrediXcan predicted gene
      expression files to run a TWAS on each file.
      Then press 'Run Analysis'.
      Navigate through the different tabs on the right to visualize your results."),

      # br() element to introduce extra vertical spacing ----
      br(),

      fileInput(
        inputId = "Predicted gene expression files",
        label = "Choose your Predicted gene expression files",
        multiple = TRUE,
        accept = c(".txt")
      ),
      tags$p("Note: if you do not upload any files and then press Run Analysis,
             we will use the files found in ./inst/extdata/ .
             The files small_cell1.txt, small_cell2.txt, ..., small_cell6.txt
             will be used as predicted gene expression input, and phenotype_covars.txt
             will be used for phenotype and covariate input."),

      actionButton("predGEBtn", "Example predicted gene expression file format"),  # Changed to actionButton for simplicity
      br(), br(),
      fileInput("phenotypeFile", "Upload a file with your phenotype data.", accept = c(".txt")),
      actionButton("alphaMissenseSampleBtn", "Example AlphaMissense file format"),  # Changed to actionButton for simplicity
      br(),
      # Side note for downloading AlphaMissense Predictions
      uiOutput("downloadPredictedGeneExpressionAndPhenoCovar"),

      actionButton("runTWAS", "Run TWAS")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Class Summary",
                           DT::dataTableOutput("classSummaryTable")),
                  tabPanel("Pathway Visualization",
                           tags$p("Please wait for the result to complete, it could take a few minutes."),
                           br(),
                           selectInput("pathwayType", "Select Pathway Type:", choices = c("all significant genes"= "sig", "up-regulated genes" = "up", "down-regulated genes" = "down")),
                           plotOutput("pathwayPlot")))
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$vcfSampleBtn, {
    showModal(modalDialog(
      title = "Example VCF Sample Format",
      DT::dataTableOutput("vcfSampleTable")
    ))
  })

  output$vcfSampleTable <- DT::renderDataTable({
    DT::datatable(head(MissensePathoR::vcfSample))
  })

  observeEvent(input$alphaMissenseSampleBtn, {
    showModal(modalDialog(
      title = "Example AlphaMissense Sample Format",
      DT::dataTableOutput("alphaMissenseSampleTable")
    ))
  })

  output$alphaMissenseSampleTable <- DT::renderDataTable({
    DT::datatable(head(MissensePathoR::AlphaMissenseSample))
  })


  # For classSummary
  runAnalysis <- eventReactive(input$runAnalysis, {
    if (is.null(input$vcfFile) || is.null(input$alphaMissenseFile)) {
      return(NULL)
    }
    # Read the files
    vcfData <- fread(input$vcfFile$datapath)
    alphaMissenseData <- fread(input$alphaMissenseFile$datapath)

    # Run predictPathoScore (assuming this function is defined in your package)
    predScoreSample <- predictPathoScore(vcfData, alphaMissenseData)

    # Run classSummary
    classSummary(predScoreSample)
  })

  # Output for classSummary
  output$classSummaryTable <- DT::renderDataTable({
    runAnalysis()
  })

  pathwayData <- eventReactive(input$runAnalysis, {
    req(input$vcfFile)
    req(input$alphaMissenseFile)

    # Read the files
    vcfData <- fread(input$vcfFile$datapath)
    alphaMissenseData <- fread(input$alphaMissenseFile$datapath)
    # Run predictPathoScore
    predScoreSample <- predictPathoScore(vcfData, alphaMissenseData)

    # Continuing the analysis pipeline
    result <- mapGene(predScoreSample)
    result$sample_name <- paste0(result$sample_name, "_", result$group)
    set.seed(1)
    diffOut <- diffSNPs(result, "0h")
    enrichOut <- enrichSNP(diffOut)
    enrichOut
  })

  # Plot output for pathwayViz
  output$pathwayPlot <- renderPlot({
    req(pathwayData())
    enrichOut <- pathwayData()
    pathwayViz(enrichOut, input$pathwayType)
  })

  # example data download
  # URL for downloading VCF sample data
  output$downloadLinkVcfSample <- renderUI({
    a("Download VCF Sample Dataset", href = "https://github.com/Lola-W/MissensePathoR/raw/main/inst/extdata/vcfSample.csv", target = "_blank")
  })

  # Modal for VCF sample data description
  observeEvent(input$downloadLinkVcfSample, {
    shinyalert(title = "VCF Sample Dataset",
               text = "This dataset contains combined VCF data for Hela cell replicates across four time points (0, 1, 4, and 8 hours) after introducing H2O2, processed with the `readVCF` function from the MissensePathoR package. Citation: Rendleman J, Cheng Z, Maity S, et al. New insights into the cellular temporal response to proteostatic stress. Elife. 2018;7:e39054. doi: 10.7554/eLife.39054.",
               type = "info")
  })

  # URL for downloading AlphaMissense sample data
  output$downloadLinkAlphaMissenseSample <- renderUI({
    a("Download AlphaMissense Sample Dataset", href = "https://github.com/Lola-W/MissensePathoR/raw/main/inst/extdata/AlphaMissenseSample.tsv", target = "_blank")
  })

  # Modal for AlphaMissense sample data description
  observeEvent(input$downloadLinkAlphaMissenseSample, {
    shinyalert(title = "AlphaMissense Sample Dataset",
               text = "A sample dataset containing missense variant predictions from the AlphaMissense community dataset resource. Citation: Jun Cheng et al., 'Accurate proteome-wide missense variant effect prediction with AlphaMissense.' Science 381, eadg7492 (2023). DOI: 10.1126/science.adg7492.",
               type = "info")
  })

  # Example data display (optional, if you want to show example data by default)
  output$classSummaryTable <- DT::renderDataTable({
    if (is.null(runAnalysis())) {
      classSummary(predictPathoScore(MissensePathoR::vcfSample, MissensePathoR::AlphaMissenseSample))
    } else {
      runAnalysis()
    }
  })

}

shinyApp(ui = ui, server = server)
