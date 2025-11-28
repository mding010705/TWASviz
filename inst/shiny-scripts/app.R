# Tabsets. Shiny. (2014, July 30). <https://shiny.posit.co/r/gallery/application-layout/tabsets>.
#
# Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J,
# Dipert A, Borges B (2025). _shiny: Web Application Framework for R_. R package
# version 1.11.1, <https://CRAN.R-project.org/package=shiny>.


library(shiny)

ui <- fluidPage(
  titlePanel("Transcriptome Wide Association Study (TWAS) Visualizations"),
  sidebarLayout(
    sidebarPanel(

      tags$p("This is a Shiny App that is part of TWASviz in R."),
      br(),

      # input
      tags$p("Instructions: Upload one or several PrediXcan association files
      or file outputs from TWASviz sparse group lasso.
      Then press 'Make Plots'.
      Navigate through the different tabs on the right to visualize your results."),
      br(),

      fileInput(
        inputId = "assoc_files",
        label = "Choose your TWAS result files",
        multiple = TRUE,
        accept = c(".txt")
      ),
      tags$p('Note: if you do not upload any TWAS result files and then
      press Make Plots,
             we will use the files found in ./inst/extdata/ .
             The files predixcan_twas1.txt, predixcan_twas2.txt, ...,
             predixcan_twas6.txt will be used as TWAS result files. These
             are the same TWAS result files in default_shiny_input_TWASviz.zip,
             the file that you are taken to by the "Download Sample Data"
             button below.'),

      br(), br(),
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
        label = "Type the p-value at which your given gene associations are
        significant",
        value = 0.2,
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
      tags$p("BP = Biological Process, MF = Molecular Function,
             CC = Cellular Component, ALL = BP, MF, and CC"),
      br(),
      selectInput(
        inputId = "organism",
        label = "Select the organism",
        choices = c("Anopheles", "Bovine", "Canine", "Chicken", "Chimp",
                    "E coli strain K12", "E coli strain Sakai",
                    "Fly", "Human", "Mouse", "Pig", "Rat", "Rhesus",
                    'Worm', 'Xenopus',
                    "Yeast", "Zebrafish"),
        selected = "Human",
        multiple = FALSE),
      br(),
      numericInput(
        inputId = "pval_go",
        label = "Type the p-value at which gene ontology enrichments are
        significant",
        value = 0.5,
        min = 0,
        max = 1),
      # Side note for downloading sample files
      uiOutput("dl_TWASs_pathways"),

      actionButton("run_plots", "Make Plots")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel(
                    "Instructions",
                    HTML('
      <h3>About the TWASviz app</h3>

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
        <li><code>predixcan_twas1.txt, predixcan_twas2.txt, ...,
             predixcan_twas6.txt</code> - different PrediXcan derived
             TWAS result files</li>
        <li><code>example_background_gene_set.txt</code> - example file for
        format of background gene set input</li>
      </ul>

      <h4>How to Use this App</h4>

      <ol>
        <li>
          Run through Application 2 in the Introduction_TWASviz
          vignette to understand how to perform a sparse group lasso based
          TWAS in R with TWASviz, or use the association output from PrediXcan.
          Be careful with your TWAS! It is easy to feel like these high dimensional
          analyses can be automated but it is important that you understand each
          component and parameter specification going each TWAS you perform.
        </li>

        <li>
          Upload one or more TWAS result files.
        </li>
        <li>
          Type the gene, effect size, and p-value column name. Leave
          the p-value column name text field empty if your file does not contain
          p-values. The gene column name does not have to contain genes, it could
          contain pathways.
        </li>
        <li>
          Type the p-value cutoff for significance for the volcano plot,
          and for use in the gene overlap heatmap and Gene Ontology enrichment
          analysis. This will be ignored if you do not specify the p-value
          column name above.
        </li>
        <li>
          Type the
          enrichment p-value cutoff, and choose the organism type,
          ontology type ("BP", "MF", "CC", "ALL") and the gene nomenclature
          ("ENSEMBL", "SYMBOL", "ENTREZID").
          This is for the Gene Ontology heatmap.
        </li>
        <li>
           (Optional) Upload an RDS file containing
          a list of vectors containing the background gene set for each of the
          uploaded TWAS files.

        </li>
        </li>
          <strong> Find your plots in the tab panels:</strong>
          <ul>
            <li><strong>Volcano Plots</strong> - Plot of effect
            size vs -log10(p-value). It will be empty if you leave the p-value
            column name text field empty.
            </li>
            <li><strong>Gene Overlap Heatmap</strong> - Showing the pairwise
            complete correlation and the number of overlapping genes between
            your TWAS files</li>
            <li><strong>Gene Ontology Enrichment Heatmap</strong> - Showing the
            Gene Ontology enrichments for each of your TWAS files. There will be
            an "error: No enrichGO results with nonzero rows found" message
            if none of your uploaded TWAS files result in significant ontology
            term enrichments at the specified p-value thresholds.</li>

          </ul>
        </li>
        </ol>')),
                  tabPanel("Volcano Plots",
                           tags$p("Please wait for the files to read, it could take a
                                  few minutes. Do not press *Make Plots* again while
                                  waiting for output."),
                           br(),
                           selectInput("fileView", "Select which file to view:",
                                       choices = c()),
                           plotOutput("vol_plot")),
                  tabPanel("Gene Overlap Heatmap",
                           tags$p("Please wait for the files to read, it could take a
                                  few minutes. Do not press *Make Plots* again while
                                  waiting for output."),
                           br(),
                           plotOutput("overlap_plot")),
                  tabPanel("Gene Ontology Enrichment Heatmap",
                           tags$p("Please wait for the GO enrichment to complete, it could take
                                  a few minutes. Do not press *Make Plots* again while
                                  waiting for output."),
                           br(),
                           plotOutput("go_plot")),
                  tabPanel("References",HTML('
      <h3>References</h3>

      <ul>
        <li>
          Aleksander, S. A., Balhoff, J., Carbon, S., Cherry, J. M., Drabkin, H. J.,
          Ebert, D., Feuermann, M., Gaudet, P., Harris, N. L., Hill, D. P., Lee, R.,
          Mi, H., Moxon, S., Mungall, C. J., Muruganugan, A., Mushayahama, T.,
          Sternberg, P. W., Thomas, P. D., Van Auken, K., … Westerfield, M. (2023).
          The gene ontology knowledgebase in 2023. GENETICS, 224(1).
          https://doi.org/10.1093/genetics/iyad031
        </li>
        <li>
          Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H.,
          Cherry, J. M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T.,
          Harris, M. A., Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S.,
          Matese, J. C., Richardson, J. E., Ringwald, M., Rubin, G. M. & Sherlock, G.
          (2000). Gene ontology: Tool for the unification of biology. Nature Genetics,
          25(1), 25–29. https://doi.org/10.1038/75556
        </li>
        <li>
          Bache, S., Wickham, H. (2022). magrittr: A Forward-Pipe Operator for R.
          R package version 2.0.3. https://CRAN.R-project.org/package=magrittr.
        </li>
        <li>
          Barrett, T., Dowle, M., Srinivasan, A., Gorecki, J., Chirico, M.,
          Hocking, T., Schwendinger, B., Krylov, I. (2025). data.table: Extension
          of `data.frame`. R package version 1.17.8,
          https://CRAN.R-project.org/package=data.table.
        </li>
        <li>
          Bates, D., Maechler, M., Jagan, M. (2025). Matrix: Sparse and Dense
          Matrix Classes and Methods. R package version 1.7-4.
          https://CRAN.R-project.org/package=Matrix.
        </li>
        <li>
          BioRender.com. BioRender. https://www.biorender.com (accessed 3 November 2025)
        </li>
        <li>
          Campitelli, E. (2025). ggnewscale: Multiple Fill and Colour Scales in
          \'ggplot2\'. R package version 0.5.2.
          https://CRAN.R-project.org/package=ggnewscale.
        </li>
        <li>
          Carlson M (2024). org.Hs.eg.db: Genome wide annotation for Human.
          R package version 3.20.0.
        </li>
        <li>
          Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Xie, Y.,
          Allen, J., McPherson, J., Dipert, A., Borges, B. (2025). shiny: Web
          Application Framework for R. R package version 1.11.1.
          https://CRAN.R-project.org/package=shiny.
        </li>
        <li>
          Devailly, G.  (2021, October 14). Plotting heatmaps from big matrices in R.
          https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/
        </li>
        <li>
          Diamant, I., Clarke, D. J. B., Evangelista, J. E., Lingam, N., & Ma’ayan, A.
          (2024). Harmonizome 3.0: Integrated knowledge about genes and proteins from
          diverse multi-omics resources. Nucleic Acids Research, 53(D1).
          https://doi.org/10.1093/nar/gkae1080
        </li>
        <li>
          Gamazon, E. R., Wheeler, H. E., Shah, K. P., Mozaffari, S. V., Aquino-Michaels,
          K., Carroll, R. J., Eyler, A. E., Denny, J. C., Nicolae, D. L.,
          Cox, N. J., & Im, H. K. (2015). A gene-based association method for
          mapping traits using reference transcriptome data. Nature Genetics, 47(9),
          1091–1098. https://doi.org/10.1038/ng.3367
        </li>
        <li>
          Kanehisa, M. (2000). Kegg: Kyoto encyclopedia of genes and genomes.
          Nucleic Acids Research, 28(1), 27–30. https://doi.org/10.1093/nar/28.1.27
        </li>
        <li>
          Liang, X., Cohen, A., Solón Heinsfeld, A., Pestilli, F., McDonald, D.J.
          (2024). “sparsegl: An R Package for Estimating Sparse Group Lasso.” Journal
          of Statistical Software, 110(6), 1-23. doi:10.18637/jss.v110.i06.
          https://doi.org/10.18637/jss.v110.i06.
        </li>
        <li>
          Milacic, M., Beavers, D., Conley, P., Gong, C., Gillespie, M., Griss, J.,
          Haw, R., Jassal, B., Matthews, L., May, B., Petryszak, R., Ragueneau, E.,
          Rothfels, K., Sevilla, C., Shamovsky, V., Stephan, R., Tiwari, K.,
          Varusai, T., Weiser, J., … D’Eustachio, P. (2023). The reactome pathway
          knowledgebase 2024. Nucleic Acids Research, 52(D1).
          https://doi.org/10.1093/nar/gkad1025
        </li>
        <li>
          OpenAI. (2025). ChatGPT (GPT-5) large language model. https://chat.openai.com/
          (accessed 28 November 2025)
        </li>
        <li>
          R Core Team (2024). R: A Language and Environment for Statistical Computing.
          R Foundation for Statistical Computing, Vienna, Austria.
          https://www.R-project.org/.
        </li>
        <li>
          Revelle, W. (2025). psych: Procedures for Psychological, Psychometric,
          and Personality Research. Northwestern University, Evanston, Illinois.
          R package version 2.5.6, https://CRAN.R-project.org/package=psych.
        </li>
        <li>
          Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013).
          A sparse-group lasso. Journal of Computational and Graphical Statistics,
          22(2), 231–245. https://doi.org/10.1080/10618600.2012.681250
        </li>
        <li>
          Slowikowski K (2024). ggrepel: Automatically Position Non-Overlapping
          Text Labels with \'ggplot2\'. R package version 0.9.6.
          https://CRAN.R-project.org/package=ggrepel.
        </li>
        <li>
          Storey, J. D., & Tibshirani, R. (2003). Statistical significance for
          genomewide studies. Proceedings of the National Academy of Sciences,
          100(16), 9440–9445. https://doi.org/10.1073/pnas.1530509100
        </li>
        <li>
          Tabsets. Shiny. (2014, July 30).
          https://shiny.posit.co/r/gallery/application-layout/tabsets
        </li>
        <li>
          Wainberg, M., Sinnott-Armstrong, N., Mancuso, N., Barbeira, A. N.,
          Knowles, D. A., Golan, D., Ermel, R., Ruusalepp, A., Quertermous, T., Hao,
          K., Björkegren, J. L., Im, H. K., Pasaniuc, B., Rivas, M. A., & Kundaje, A.
          (2019). Opportunities and challenges for transcriptome-wide association
          studies. Nature Genetics, 51(4), 592–599.
          https://doi.org/10.1038/s41588-019-0385-z
        </li>
        <li>
          Wang, N., Ye, Z., & Ma, T. (2024). Tips: A novel pathway-guided joint model
          for transcriptome-wide association studies. Briefings in Bioinformatics, 25(6).
          https://doi.org/10.1093/bib/bbae587
        </li>
        <li>
          Wickham, H., François, R., Henry, L., Müller, K., Vaughan, D. (2023).
          dplyr: A Grammar of Data Manipulation. R package version 1.1.4.
          https://CRAN.R-project.org/package=dplyr.
        </li>
        <li>
          Wickham, H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
          New York, 2016.
        </li>
        <li>
          Wickham, H. (2025). stringr: Simple, Consistent Wrappers for Common String
          Operations. R package version 1.6.0. https://CRAN.R-project.org/package=stringr.
        </li>
        <li>
          Wickham, H., Vaughan, D., Girlich, M. (2024). tidyr: Tidy Messy Data.
          R package version 1.3.1. https://CRAN.R-project.org/package=tidyr.
        </li>
        <li>
          Xu, S., Hu, E., Cai, Y., Xie, Z., Luo, X., Zhan, L., Tang, W., Wang, Q.,
          Liu, B., Wang, R., Xie, W., Wu, T., Xie, L., Yu, G. Using clusterProfiler
          to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320.
        </li>

      </ul>'))
    )
  )
))

server <- function(input, output, session) {
  # When run_plots pressed
  save_inputs <- eventReactive(input$run_plots, {
    list(gene_colname = input$gene_colname,
         eff_size_colname = input$eff_size_colname,
         pval_colname = input$pval_colname, pval_inc = input$pval_inc,
         organism = input$organism,
         pval_go = input$pval_go, go_sub = input$go_sub,
         gene_nom = input$gene_nom)})

  proc_input <- eventReactive(input$run_plots, {
    # Read the files
    if (is.null(input$assoc_files)) {
      upload_assoc <- list()
      assoc_fns <- paste0("../extdata/predixcan_twas", 1:6, ".txt")
      for(i in seq_along(assoc_fns)){
        upload_assoc[[i]] <- as.data.frame(data.table::fread(
          file = assoc_fns[i], header = TRUE))
      }
      tissue_names <- sub("(.*)\\..*$", "\\1", basename(assoc_fns))
    } else {
      upload_assoc <- list()
      assoc_fns <- seq_along(input$assoc_files[, 1])
      tissue_names <- seq_along(input$assoc_files[, 1])
      for(i in seq_along(input$assoc_files[, 1])){
        upload_assoc[[i]] <- as.data.frame(data.table::fread(
          file = input$assoc_files[[i, 'datapath']], header = TRUE))
        assoc_fns[i] <- input$assoc_files[[i, 'datapath']]
        tissue_names[i] <- sub("(.*)\\..*$", "\\1",
                               basename(input$assoc_files[[i, 'name']]))
      }
    # colname input validation
      if (!(input$gene_colname %in% colnames(upload_assoc[[i]]))){
        stop("Gene column name not in your TWAS input files.")
      }
      if (!(input$eff_size_colname %in% colnames(upload_assoc[[i]]))){
        stop("Effect size column name not in your TWAS input files.")
      }
      if (!(input$pval_colname == "") &&
          !(input$pval_colname %in% colnames(upload_assoc[[i]]))){
        stop("p-value column name not in your TWAS input files.")
      }

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


    if (!(input$pval_colname == "")){ # update drop down selection to include file names
      vol_plot_choices <- seq_along(assoc_fns)
      names(vol_plot_choices) <- basename(assoc_fns)
      updateSelectInput(session, "fileView",
                        label = "Select which file to view:",
                        choices = vol_plot_choices,
                        selected = tail(vol_plot_choices, 1))
  } else { # no p-values -> no volcano plots
    updateSelectInput(session, "fileView",
                      label = "No p-values, no volcano plot",
                      choices = c())
   }
    return(list(assoc_fns = assoc_fns,
                upload_assoc = upload_assoc,
                tissue_names = tissue_names,
                upload_bkgrnd = upload_bkgrnd))
  })

  # Volcano plot tab
  output$vol_plot <- renderPlot({
    req(proc_input)
    if (!(save_inputs()$pval_colname == "")){
        TWASviz::volcano_plot(betas_p_vals =
                         proc_input()$upload_assoc[[as.numeric(input$fileView)]],
                              gene_colname = save_inputs()$gene_colname,
                              effect_size_colname = save_inputs()$eff_size_colname,
                              pvalue_colname = save_inputs()$pval_colname,
                              p_thresh = save_inputs()$pval_inc)
    } else {
      vol_plot_list <- NULL
      return("p-values not available")
    }
  })

  # Overlap plot tab
  output$overlap_plot <- renderPlot({
    req(proc_input)
    adj_df <- TWASviz::predixcan2adj_df(predixcan_assoc_filenames = proc_input()$assoc_fns,
                                        gene_colname = save_inputs()$gene_colname,
                                        effect_size_colname = save_inputs()$eff_size_colname,
                                        pvalue_colname = save_inputs()$pval_colname,
                                        pvalue_thresh = save_inputs()$pval_inc,
                                        tissue_names = proc_input()$tissue_names)
    TWASviz::correlation_overlap_heatmap(adj_df,
              tissue_names = proc_input()$tissue_names)
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
    if (!(save_inputs()$pval_colname == "")){
      genes <- lapply(proc_input()$upload_assoc, function(x){
        x[x[,save_inputs()$pval_colname] <= save_inputs()$pval_inc, ][,save_inputs()$gene_colname]
      })
    } else {
      genes <- lapply(proc_input()$upload_assoc, function(x){
        x[,save_inputs()$gene_colname]
      })
    }

    enrich_res <- TWASviz::gene_enrichment(genes,
                                           tissue_types = proc_input()$tissue_names,
                                           organism = org2db[save_inputs()$organism],
                                           background = proc_input()$upload_bkgrnd,
                                           p_cutoff = save_inputs()$pval_go,
                                           ont_type = save_inputs()$go_sub,
                                           gene_nom = save_inputs()$gene_nom)

    TWASviz::goenrich_heatmap(enrich_res = enrich_res, x_label = "")
  })

  # Example data download
  output$dl_TWASs_pathways <- renderUI({
    a("Download Sample Data", href = "https://github.com/mding010705/TWASviz/blob/main/inst/extdata/default_shiny_input_TWASviz.zip", target = "_blank")
  })

}

shinyApp(ui = ui, server = server)

# [END]
