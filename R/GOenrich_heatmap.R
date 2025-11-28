#' Perform Gene Ontology Enrichments for Multiple Tissue/Cell Types
#'
#' Runs Gene Ontology enrichment analysis across multiple tissues/cell types
#' using **clusterProfiler::enrichGO()**. Each entry in the input list is treated
#' as a gene set from one tissue/condition, and enrichment results are returned
#' for each.
#'
#' This function optionally accepts a background/universe gene list for each
#' tissue type. If not supplied, enrichGO() uses the entire gene universe
#' defined by the OrgDb annotation.
#'
#' **References**
#' - Yu, G. (2024). *Thirteen years of clusterProfiler.*
#'   The Innovation, 5(6):100722.
#' - Yu, G., Wang, L.-G., Han, Y., & He, Q.-Y. (2012).
#'   *clusterProfiler: an R package for comparing biological themes among
#'   gene clusters.* OMICS: A Journal of Integrative Biology, 16(5), 284–287.
#'
#' @param gene_set List of character vectors. Each element contains significant
#'   gene symbols/IDs for a tissue or cell type.
#' @param organism Character string giving the OrgDb annotation database
#'   (e.g., `"org.Hs.eg.db"`). You must have the annotation package installed.
#' @param tissue_types Optional vector giving names for each entry in gene_set.
#'   Defaults to the names of gene_set.
#' @param background Optional list of background (universe) genes per tissue.
#'   If provided, each entry must match the corresponding entry in gene_set.
#' @param p_adj_method Adjustment method for p-values. One of:
#'   `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`,
#'   `"fdr"`, `"none"`.
#' @param min_gene_set_size Minimum GO term size considered.
#' @param max_gene_set_size Maximum GO term size considered.
#' @param p_cutoff Maximum p-value threshold for enrichment.
#' @param ont_type GO subontology: `"BP"` (Biological Process),
#'   `"MF"` (Molecular Function), `"CC"` (Cellular Component), or `"ALL"`.
#' @param gene_nom Gene ID type used (e.g. `"SYMBOL"`, `"ENSEMBL"`, `"ENTREZID"`).
#'
#' @return A named list of enrichGO objects (one per tissue/cell type).
#'
#' @examples
#' \dontrun{
#'   gene_set <- list(
#'     Liver = c("CHST2", "B3GNT1", "B3GNT2"),
#'     Brain = c("ST3GAL1", "CHST4", "FUT8")
#'   )
#'
#'   enrich_res <- gene_enrichment(
#'     gene_set = gene_set,
#'     organism = "org.Hs.eg.db",
#'     gene_nom = "SYMBOL"
#'   )
#'
#'   enrich_res$Liver
#' }
#'
#' @export
#' @importFrom clusterProfiler enrichGO


gene_enrichment <- function(gene_set,
                            organism = "org.Hs.eg.db",
                            tissue_types = names(gene_set),
                            background = NULL,
                            p_adj_method = "fdr",
                            min_gene_set_size = 3,
                            max_gene_set_size = 800,
                            p_cutoff = 0.05,
                            ont_type = "BP",
                            gene_nom = "ENSEMBL") {

  # Input Validation
  if (!is.list(gene_set))
    stop("gene_set must be a list of character vectors.")

  if (is.null(tissue_types))
    stop("tissue_types must be provided or gene_set must have names.")

  if (!is.null(background) && (length(background) != length(gene_set)))
    stop("background must be NULL or a list of same length as gene_set.")

  if (!gene_nom %in% c("SYMBOL", "ENSEMBL", "ENTREZID"))
    warning("gene_nom is not one of SYMBOL, ENSEMBL, ENTREZID. Ensure it matches OrgDb.")

  if (!ont_type %in% c("BP", "MF", "CC", "ALL"))
    stop("ont_type must be one of 'BP', 'MF', 'CC', or 'ALL'.")

  # Try loading the OrgDb annotation package
  if (!requireNamespace(organism, quietly = TRUE)) {
    stop(paste0(
      "The organism annotation package '", organism,
      "' is not installed. Install using BiocManager::install()."
    ))
  }

  enrich_results <- list()

  # Perform enrichment for each tissue
  for (i in seq_along(gene_set)) {

    genes <- gene_set[[i]]

    if (!is.character(genes))
      stop("Each entry of gene_set must be a character vector of gene IDs.")

    # No background universe supplied
    if (is.null(background) || is.null(background[[i]])) {
      enrich_results[[tissue_types[i]]] <-
        clusterProfiler::enrichGO(
          gene          = genes,
          OrgDb         = organism,
          keyType       = gene_nom,
          ont           = ont_type,
          pAdjustMethod = p_adj_method,
          pvalueCutoff  = p_cutoff,
          minGSSize     = min_gene_set_size,
          maxGSSize     = max_gene_set_size,
          readable      = TRUE
        )

    # Background provided
    } else {
      enrich_results[[tissue_types[i]]] <-
        clusterProfiler::enrichGO(
          gene          = genes,
          universe      = as.character(background[[i]]),
          OrgDb         = organism,
          keyType       = gene_nom,
          ont           = ont_type,
          pAdjustMethod = p_adj_method,
          pvalueCutoff  = p_cutoff,
          minGSSize     = min_gene_set_size,
          maxGSSize     = max_gene_set_size,
          readable      = TRUE
        )
    }
  }

  return(enrich_results)
}


#' Plot Gene Ontology Enrichment Heatmap Across Tissues/Cell Types
#'
#' Visualizes enriched GO terms by generating a heatmap of
#' \eqn{-\log_{10}(adjusted\ p\text{-}values)} across tissues/cell types.
#'
#' Each tissue contributes its *top N* most significant GO terms.
#' Terms are clustered and displayed along the y-axis.
#'
#' @param enrich_res List of enrichGO objects (output from gene_enrichment).
#' @param top_n Integer specifying number of GO terms to plot per tissue.
#' @param x_label Label for x-axis (default: `"Cell Type"`).
#'
#' @return A ggplot2 heatmap object.
#'
#' @examples
#' \dontrun{
#'   enrich_res <- gene_enrichment(
#'     gene_set = list(
#'       TissueA = c("CHST2", "B3GNT1", "B3GNT2"),
#'       TissueB = c("ST3GAL1", "CHST4", "FUT8")
#'     ),
#'     organism = "org.Hs.eg.db",
#'     gene_nom = "SYMBOL"
#'   )
#'
#'   p <- goenrich_heatmap(enrich_res, top_n = 5)
#'   print(p)
#' }
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2


goenrich_heatmap <- function(enrich_res, top_n = 5, x_label = "Cell Type") {

  # Input validation
  if (!is.list(enrich_res))
    stop("enrich_res must be a list of enrichGO result objects.")

  if (!all(sapply(enrich_res, function(x) inherits(x, "enrichResult") || is.null(x))))
    stop("All non-null entries of enrich_res must be 'enrichResult' objects.")

  if (top_n <= 0)
    stop("top_n must be a positive integer.")

  # Remove null or empty results
  filtered_results <- enrich_res[unlist(sapply(enrich_res,
                                        function(x) (!is.null(x)) && (nrow(x) > 0)))]
  if (length(filtered_results) == 0)
    stop("No enrichGO results with nonzero rows found.")


  # Extract top N GO terms per tissue
  plot_data <- lapply(names(filtered_results), function(mod) {
    df <- filtered_results[[mod]]@result
    df$Module <- mod

    df %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n)
  }) %>% dplyr::bind_rows()

  # Shorten GO term names for readability
  plot_data$Term <- stringr::str_trunc(plot_data$Description, 40)

  # Ensure uniqueness of terms per module
  plot_data <- plot_data[!duplicated(plot_data[c("Module", "Term")]), ]

  # Prepare GO term × tissue/cell type matrix
  mat_df <- plot_data %>%
    dplyr::mutate(logp = -log10(p.adjust)) %>%
    dplyr::select(Module, Term, logp) %>%
    tidyr::pivot_wider(names_from = Term,
                       values_from = logp,
                       values_fill = 0)

  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$Module

  # Hierarchical clustering of tissues/cell types
  row_order <- tryCatch(
    hclust(dist(mat))$order,
    error = function(e) seq_len(nrow(mat))
  )

  ordered_modules <- rownames(mat)[row_order]
  plot_data$Module <- factor(plot_data$Module, levels = ordered_modules)

  # Produce heatmap
  return(
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = Module, y = Term)
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = -log10(p.adjust)),
        color = "white"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = Count), size = 3) +
      ggplot2::scale_fill_gradient(
        low = "white",
        high = "red",
        name = "-log10 adj p"
      ) +
      ggplot2::labs(
        x = x_label,
        y = "GO Term",
        title = "GO Enrichment Heatmap"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      ))
}

# [END]
