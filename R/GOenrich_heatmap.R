#' Perform Gene Ontology Enrichments for Multiple Tissue/Cell Types
#'
#' Uses clusterProfiler function enrichGO() to find Gene Ontology (GO)
#' enrichments for multiple tissue types.
#'
#' @param gene_set List where each entry is a set of significant genes for
#' a tissue.
#' @param organism Genome wide annotation for your organsim of choice.
#' You must have your chosen annotation package installed.
#' @param tissue_types Vector of names that correspond to entries of gene_set.
#' @param background List where each entry is a vector of gene names that
#' describe the transcriptomic universe of the corresponding
#' tissue type in gene_set.
#' @param p_adj_method p-value adjustment method. One of "holm", "hochberg",
#' "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min_gene_set_size Integer describing the smallest GO term size to be
#' considered by enrichGO.
#' @param max_gene_set_size Integer describing the largest GO term size to be
#' considered by enrichGO.
#' @param p_cutoff Max p-value that will be considered as significant.
#' @param ont_type Which GO subontology to query.
#' One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param gene_nom Nomenclature/keytype of genes in gene_set inputs.
#'
#' @return Returns a list of enrichGO() results.
#'
#' @examples
#' gene_enrichment(gene_set = list(a=c("CHST2", "B3GNT1", "B3GNT2", "CHST1",
#' "B4GALT1","B3GNT7", "ST3GAL2", "B4GALT3", "ST3GAL1", "B4GALT2", "B4GALT4",
#' "CHST4", "ST3GAL3", "FUT8", "CHST6")),
#' organism = "org.Hs.eg.db", gene_nom = "SYMBOL")
#'
#' @references
#' G Yu. Thirteen years of clusterProfiler. The Innovation. 2024, 5(6):100722
#'
#' @export
#' @import clusterProfiler

gene_enrichment <- function(gene_set, organism = "org.Hs.eg.db",
                            tissue_types = names(gene_set),
                            background = NULL,
                                    p_adj_method = "fdr",
                                    min_gene_set_size = 3,
                                    max_gene_set_size = 800,
                           p_cutoff = 0.05, ont_type = "BP",
                           gene_nom = "ENSEMBL"){

  enrich_results <- list()
    for (i in seq_along(gene_set)){
      if(is.null(background[[i]])){
        enrich_results[[tissue_types[i]]] <-
          clusterProfiler::enrichGO(
            gene = gene_set[[i]],
            OrgDb = organism,
            keyType = gene_nom,
            ont = ont_type,
            pAdjustMethod = p_adj_method,
            pvalueCutoff = p_cutoff,
            readable = TRUE)
      } else {
        enrich_results[[tissue_types[i]]] <-
          clusterProfiler::enrichGO(
            gene = gene_set[[i]],
            universe = as.character(sapply(as.vector(background[[i]]),
                                          as.character)),
            OrgDb = organism,
            keyType = gene_nom,
            ont = ont_type,
            pAdjustMethod = p_adj_method,
            pvalueCutoff = p_cutoff,
            readable = TRUE)
      }

    }
    return(enrich_results)
}


#' Plot Gene Ontology Enrichments for Multiple Tissue/Cell Types
#'
#' Plots the -log10(adjusted p-value) of the Gene Ontology (GO) enrichments
#' for multiple tissue types.
#'
#' @param enrich_res List where each entry is enrichGO result for a tissue.
#' @param top_n Number of GO terms to use per tissue type.
#'
#' @return Returns a heatmap plot of p-value strength of each tissue's
#' enrichment for top GO terms.
#'
#' @examples
#' goenrich_heatmap(enrich_res = gene_enrichment(gene_set = list(a=c("CHST2",
#' "B3GNT1", "B3GNT2", "CHST1",
#' "B4GALT1","B3GNT7", "ST3GAL2", "B4GALT3", "ST3GAL1", "B4GALT2", "B4GALT4",
#' "CHST4", "ST3GAL3", "FUT8", "CHST6")),
#' organism = "org.Hs.eg.db", gene_nom = "SYMBOL"))
#'
#' @export
#' @import dplyr
#' @import stringr
#' @import ggplot2

goenrich_heatmap <- function(enrich_res, top_n = 5, x_label = "Cell Type"){
  # Select top N GO terms per module
  filtered_results <- enrich_res[sapply(enrich_res,
                                        function(x) !is.null(x) && nrow(x) > 0)]

  plot_data <- lapply(names(filtered_results), function(mod) {
    df <- filtered_results[[mod]]@result
    df$Module <- mod
    dplyr::arrange(df, p.adjust) %>%
      dplyr::slice_head(n = top_n)
  }) %>% dplyr::bind_rows()

  # Truncate GO term names
  plot_data$Term <- stringr::str_trunc(plot_data$Description, 40)
  plot_data <-  plot_data[!duplicated(plot_data[c("Module","Term")]),]

  # Create module x GO term matrix for clustering modules
  mat_df <- plot_data %>%
    dplyr::mutate(logp = -log10(p.adjust)) %>%
    dplyr::select(Module, Term, logp) %>%
    tidyr::pivot_wider(names_from = Term, values_from = logp, values_fill = 0)

  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- mat_df$Module

  # Cluster modules
  row_order <- tryCatch(hclust(dist(mat))$order, error=function(e){seq(along = mat)})
  ordered_modules <- rownames(mat)[row_order]
  plot_data$Module <- factor(plot_data$Module, levels = ordered_modules)

  # Plot heatmap with GO terms on y-axis
  return(ggplot2::ggplot(plot_data, ggplot2::aes(x = Module, y = Term)) +
    ggplot2::geom_tile(ggplot2::aes(fill = -log10(p.adjust)), color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = Count), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "red", name = "-log10 adj p") +
    ggplot2::labs(x = x_label, y = "GO Term", title = "GO Enrichment Heatmap") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   axis.text.y = element_text(size = 10),
                   panel.grid = element_blank()))
}

# [END]
