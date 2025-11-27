library(shiny)

ui <- fluidPage(
  titlePanel("Transcriptome Wide Association Analysis (TWAS) Visualizations"),
  sidebarLayout(
    sidebarPanel(
      tags$p("This is a Shiny App that is part of TWASviz in R."),
      br(),

      tags$b("Description: ."),

      # br() element to introduce extra vertical spacing
      br(),
      br(),

      # input
      tags$p("Instructions: Upload one or several PrediXcan association files
      or file output from TWASviz sparse group lasso.
      Then press 'Make Plots'.
      Navigate through the different tabs on the right to visualize your results."),
      br(),

      fileInput(
        inputId = "assoc_files",
        label = "Choose your TWAS result files",
        multiple = TRUE,
        accept = c(".txt")
      ),
      tags$p("Note: if you do not upload any TWAS result files and then
      press Make Plots,
             we will use the files found in ./inst/extdata/ .
             The files predixcan_twas1.txt, predixcan_twas2.txt, ...,
             predixcan_twas6.txt will be used as TWAS result files, and
             all_pathways.rds will be used as pathway input."),

      br(), br(),
      fileInput("pathway_file", "Upload a file with your pathway list (Optional)", accept = c(".rds", ".RDS", ".Rds")),
      br(),
      fileInput("bkgrnd_file",
                "Upload files in order with your background gene set data for each gene expression input (Optional)",
                multiple = TRUE,
                accept = c(".txt")),
      br(),
      textInput(
        inputId = "gene_colname",
        label = "Type the gene column name",
        value = "gene"),
      br(),
      textInput(
        inputId = "eff_size_colname",
        label = "Type the effect size column name",
        value = "zscore"),
      br(),
      textInput(
        inputId = "pval_colname",
        label = "Type the p-value column name",
        value = "pvalue"),
      br(),
      numericInput(
        inputId = "pval_inc",
        label = "Type the p value at which your given gene associations are
        significant",
        value = 0.05,
        min = 0,
        max = 1),
      br(),
      selectInput(
        inputId = "gene_nom",
        label = "Select the gene nomenclature in your TWAS files",
        choices = c("ENSEMBL", "SYMBOL", "ENTREZID"),
        selected = NULL,
        multiple = FALSE),
      br(),
      selectInput(
        inputId = "go_sub",
        label = "Select the gene ontology subontology in which to find
        enrichments",
        choices = c("BP", "MF", "CC", "ALL"),
        selected = "BP",
        multiple = FALSE),
      br(),
      selectInput(
        inputId = "organism",
        label = "Select the organism",
        choices = c("Anopheles", "Bovine", "Canine", "Chicken", "Chimp",
                    "E coli strain K12, E coli strain Sakai",
                    "Fly", "Human", "Mouse", "Pig", "Rat", "Rhesus",
                    "Streptomyces coelicolor", 'Worm', 'Xenopus',
                    "Yeast", "Zebrafish"),
        selected = "Human",
        multiple = FALSE),
      br(),
      numericInput(
        inputId = "pval_go",
        label = "Type the p value at which gene ontology enrichments are
        significant",
        value = 0.05,
        min = 0,
        max = 1),
      # Side note for downloading sample files
      uiOutput("dl_TWASs_pathways"),

      actionButton("run_plots", "Make Plots")
    )),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Instructions",
                    HTML('
    <div style="max-width: 800px; line-height: 1.6;">

      <h3>This is the TWASviz app</h3>

      <h4>About the TWASviz app</h4>

      <p>
        The TWASviz app allows users to take Transcriptome Wide Association Study
        (TWAS) results from either the sparse group lasso association test from
        TWASviz itself in R or PrediXcan in python
        and to effectively visualize these outputs. This app provides
        visualization of the relationship between the p-value and effect sizes,
        the overlap of sets of important gene transcripts
        from multiple cell types or tissues, and the Gene Ontology term
        enrichments from the TWAS.
      </p>

      <p>
        Sample files are available in <code>inst/extdata</code> of
        the TWASviz GitHub repository.
      </p>

      <ul>
        <li><code>A_sgl_TWASgene.txt, E_sgl_TWASgene.txt,
             F_sgl_TWASgene.txt</code> - different sparse group lasso derived
             TWAS result files</li>
             code>predixcan_twas1.txt, predixcan_twas2.txt, ...,
             predixcan_twas6.txt</code> - different PrediXcan derived
             TWAS result files</li>
        <li><code>all_pathways.rds</code> - file with 1077 different pathways
        from KEGG, Reactome, and Biocarta</li>
        <li><code>example_background_gene_set.txt</code> - example file for
        format of background gene set input</li>
      </ul>

      <h4>How to Use This App</h4>

      <ol>
        <li>
          <strong>Run through Application 2 in the Introduction_TWASviz
          vignette to understand how to perform a sparse group lasso based
          TWAS in R with TWASviz, or use the association output from PrediXcan.
          Be careful with your TWAS! It is easy to feel like these high dimensional
          analyses can be automated but it is important that you understand each
          component and parameter specification going each TWAS you perform.
        </li>

        <li>
          <strong>Upload one or more TWAS result files.
        </li>
        <li>
          <strong>Type the gene, effect size, and p-value column name. Leave
          the p-value column name text field empty if your file does not contain
          p-values. The gene column name does not have to contain genes, it could
          contain pathways.
        </li>
        <li>
          <strong>Type the p-value cutoff for significance for the volcano plot,
          and for use in the gene overlap heatmap and Gene Ontology enrichment
          analysis. This will be ignored if you do not specify the p-value
          column name above.
        </li>
        <li>
          <strong>Type the
          enrichment p-value cutoff, and choose the organism type,
          ontology type ("BP", "MF", "CC", "ALL") and the gene nomenclature
          ("ENSEMBL", "SYMBOL", "ENTREZID"). This is for the Gene Ontology heatmap.
        </li>
        <li>
          <strong> (Optional) Upload RDS files containing a pathway list and/or
          a list of vectors containing the background gene set for each of the
          uploaded TWAS files.

        </li>
        </li>
          Find your plots in the tab panels:
          <ul>
            <li><strong>Volcano Plots</strong> - Plot of effect
            size vs -log10(p-value). It will be empty if you leave the p-value
            column name text field empty.
            </li>
            <li><strong>Gene Overlap Heatmap</strong> - Showing the pairwise
            complete correlation and the number of overlapping genes between
            your TWAS files</li>
            <li><strong>Gene Ontology Enrichment Heatmap</strong> - Showing the
            Gene Ontology enrichments for each of your TWAS files</li>

          </ul>
        </li>

    </div>
  ')),
                  tabPanel("Volcano Plots",
                           tags$p("Please wait for the files to read, it could take a few minutes."),
                           br(),
                           plotOutput("vol_plot")),
                  tabPanel("Gene Overlap Heatmap",
                           tags$p("Please wait for the files to read, it could take a few minutes."),
                           br(),
                           plotOutput("overlap_plot")),
                  tabPanel("Gene Ontology Enrichment Heatmap",
                           tags$p("Please wait for the GO enrichment to complete, it could take a few minutes."),
                           br(),
                           plotOutput("go_plot")),
                  tabPanel("References",HTML('
    <div style="max-width: 800px; line-height: 1.6;">

      <h3>References</h3>

      <ul>
        <li>
          BioRender.com. <em>BioRender</em> [Online]. Available at:
          <a href="https://www.biorender.com"
          target="_blank">https://www.biorender.com</a>
          (accessed 26 October 2025).
        </li>

        <li>
          Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. (2025).
          <em>Shiny: Web Application Framework for R — Tabsets example.</em>
          Shiny Gallery, RStudio.
          <a href="https://shiny.posit.co/r/gallery/application-layout/tabsets"
          target="_blank">
            Link
          </a>
        </li>

        <li>
          Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. (2025).
          <em>Shiny: Web Application Framework for R — File upload example.</em>
          Shiny Gallery, RStudio.
          <a href="https://shiny.posit.co/r/gallery/widgets/file-upload/"
          target="_blank">
            Link
          </a>
        </li>

        <li>
          OpenAI. (2025). ChatGPT (GPT-5) large language model.
          <a href="https://chat.openai.com/" target="_blank">
          https://chat.openai.com/</a>
        </li>

        <li>
          Pagès H, Aboyoun P, Gentleman R & DebRoy S. (2025).
          <em>Biostrings: Efficient manipulation of biological strings</em>
          (R package version 2.77.2).
          <a href="https://bioconductor.org/packages/Biostrings"
          target="_blank">
            Bioconductor
          </a>, doi:10.18129/B9.bioc.Biostrings
        </li>

        <li>
          Park J, Kaufman E, Valdmanis PN & Bafna V. (2023).
          <em>TRviz: A Python Library for decomposing and Visualizing
          Tandem Repeat Sequences.</em>
          Bioinformatics Advances 3.
        </li>

        <li>
          R Core Team. (2025). <em>R: A Language and Environment for Statistical
          Computing.</em>
          Vienna, Austria: R Foundation for Statistical Computing.
          <a href="https://www.R-project.org/" target="_blank">
          https://www.R-project.org/</a>
        </li>

        <li>
          Wickham H. (2016). <em>ggplot2: Elegant Graphics for Data Analysis.
          </em>
          Springer-Verlag, New York.
        </li>

        <li>
          Wickham H. (2007). Reshaping data with the reshape package.
          <em>J. Stat. Softw.</em> 21, 1–20.
        </li>

      </ul>

    </div>
  '))
    )
  )
)
server <- function(input, output, session) {
  # For classSummary
  make_plots <- eventReactive(input$run_plots, {
    # Read the files
    if (is.null(input$assoc_files)) {
      upload_assoc <- list()
      assoc_fns <- paste0("../extdata/predixcan_twas", 1:6, ".txt")
      for(i in seq_along(assoc_fns)){
        upload_assoc[[i]] <- as.data.frame(data.table::fread(
          file = default_files[i], header = TRUE))
      }
      tissue_names <- sub("(.*)\\..*$", "\\1", basename(default_files))
    } else {
      upload_assoc <- list()
      assoc_fns <- seq_along(input$assoc_files[, 1])
      tissue_names <- seq_along(input$assoc_files[, 1])
      for(i in seq_along(input$assoc_files[, 1])){
        upload_assoc[[i]] <- as.data.frame(data.table::fread(
          file = input$assoc_files[[i, 'datapath']], header = TRUE))
        assoc_fns[i] <- input$assoc_files[[i, 'datapath']]
        tissue_names[i] <- sub("(.*)\\..*$", "\\1",
                               basename(input$assoc_files[[i, 'datapath']]))
      }

    }

    if (!is.null(input$pathway_file)) {
      all_pathways <- readRDS(input$pathway_file$datapath)
    } else {
      all_pathways <- readRDS("../extdata/all_pathways.rds")
    }

    if (!is.null(input$bkgrnd_files) &&
        (length(input$bkgrnd_files) == length(upload_assoc))) {
      upload_bkgrnd <- list()
      for(i in seq_along(input$bkgrnd_files[, 1])){
        upload_bkgrnd[[i]] <- as.data.frame(data.table::fread(
          file = input$bkgrnd_files[[i, 'datapath']], header = FALSE))
      }
    } else {
      upload_bkgrnd <- NULL
    }

    if (!(input$gene_colname %in% colnames(upload_assoc[[1]]))){
      stop("Gene column name not in your TWAS input files.")
    }
    if (!(input$eff_size_colname %in% colnames(upload_assoc[[1]]))){
      stop("Effect size column name not in your TWAS input files.")
    }
    if (!is.null(input$pval_colname) &&
        !(input$pval_colname %in% colnames(upload_assoc[[1]]))){
      stop("p-value column name not in your TWAS input files.")
    }

    return(list(assoc_fns = assoc_fns,
                upload_assoc = upload_assoc,
                tissue_names = tissue_names,
                upload_bkgrnd = upload_bkgrnd))
  })

  proc_input <- read_input()

  # Volcano plot tab
  output$vol_plot <- renderPlot({
    req(proc_input)
    if (!is.null(input$pval_colname)){
      vol_plot_list <- lapply(proc_input$upload_assoc, function(x){
        TWASviz::volcano_plot(betas_p_vals = x,
                              gene_colname = input$gene_colname,
                              effect_size_colname = input$eff_size_colname,
                              pvalue_colname = input$pval_colname,
                              p_thresh = input$pval_inc)
      })
    } else {
      vol_plot_list <- NULL
    }
    cowplot::plot_grid(vol_plot_list, labels = proc_input$tissue_names)
  })

  # Overlap plot tab
  output$overlap_plot <- renderPlot({
    req(proc_input)
    adj_df <- TWASviz::predixcan2adj_df(predixcan_assoc_filenames = proc_input$assoc_fns,
                                        gene_colname = input$gene_colname,
                                        effect_size_colname = input$eff_size_colname,
                                        pvalue_colname = input$pval_colname,
                                        pvalue_thresh = input$pval_inc,
                                        tissue_names = proc_input$tissue_names)
    TWASviz::correlation_overlap_heatmap(adj_df,
              tissue_names = proc_input$tissue_names)
  })

  # GO plot tab
  output$go_plot <- renderPlot({
    req(proc_input)
    org2db <- c(Anopheles = "org.Ag.eg.db", Bovine = "org.Bt.eg.db",
                Canine = "org.Cf.eg.db", Chicken = "org.Gg.eg.db",
                Chimp = "org.Pt.eg.db", `E coli strain K12` = "org.EcK12.eg.db",
                `E coli strain Sakai` = "org.EcSakai.eg.db", Fly = "org.Dm.eg.db",
                Human = "org.Hs.eg.db", Mouse = "org.Mm.eg.db",
                Pig = "org.Ss.eg.db", Rat = "org.Rn.eg.db",
                Rhesus = "org.Mmu.eg.db",
                `Streptomyces coelicolor` = "org.Sco.eg.db", Worm = "org.Ce.eg.db",
                Xenopus = "org.Xl.eg.db", Yeast = "org.Sc.sgd.db",
                Zebrafish = "org.Dr.eg.db")
    if (!is.null(input$pval_colname)){
      genes <- lapply(proc_input$upload_assoc, function(x){
        x[x[,input$pval_colname] <= input$pval_inc, ][,input$gene_colname]
      })
    } else {
      genes <- lapply(proc_input$upload_assoc, function(x){
        x[,input$gene_colname]
      })
    }

    names(genes) <- proc_input$tissue_names
    enrich_res <- TWASviz::gene_enrichment(genes,
                                           organism = org2db[input$organism],
                                           background = proc_input$upload_bkgrnd,
                                           p_cutoff = input$pval_go,
                                           ont_type = input$go_sub,
                                           gene_nom = input$gene_nom)

    TWASviz::goenrich_heatmap(enrich_res = enrich_res, x_label = "")
  })

  # example data download
  output$dl_TWASs_pathways <- renderUI({
    a("Download Sample Data", href = "https://github.com/Lola-W/MissensePathoR/raw/main/inst/extdata/vcfSample.csv", target = "_blank")
  })

  # Data description
  observeEvent(input$dl_TWASs_pathways, {
    shinyalert(title = "VCF Sample Dataset",
               text = "This dataset contains combined VCF data for Hela cell replicates across four time points (0, 1, 4, and 8 hours) after introducing H2O2, processed with the `readVCF` function from the MissensePathoR package. Citation: Rendleman J, Cheng Z, Maity S, et al. New insights into the cellular temporal response to proteostatic stress. Elife. 2018;7:e39054. doi: 10.7554/eLife.39054.",
               type = "info")
  })

}

shinyApp(ui = ui, server = server)
