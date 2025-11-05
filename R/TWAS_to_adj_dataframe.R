
#' Format Coefficient Output from Multiple Sparse Group Lassos into Adjacency Matrix
#'
#' Produces an adjacency-matrix-like dataframe where each row is a gene and
#' each column is the summed estimated effect size of the gene from a sparse
#' group lasso TWAS. Also produces another adjacency-matrix-like dataframe
#' but with pathways.
#'
#' @param gene_coefs List of coefficient output from cv.sparsegl object.
#' @param tissue_names Names for the tissues/cell/condition types for the
#' imputed gene expression data used to train each cv.sparsegl object that
#' produced the above coefficients.
#'
#' @return Returns a list(
#' gene = gene adjacency-matrix-like data frame,
#' pathway = pathway adjacency-matrix-like data frame
#' ).
#'
#' @export
#' @import stringr

sgl2adj_df <- function(gene_coefs,
                       tissue_names = seq_along(gene_coefs)){
  gene_adj_df <- data.frame(gene = c(NA))
  pathway_adj_df <- data.frame(pathway = c(NA))

  for(coef in gene_coefs){
    coef <- coef[coef[,1]!= 0,]
    gene_pathway <- stringr::str_split_fixed(rownames(coef), "_", n = 2)
    gene_pathway[gene_pathway[,2]=="", 2] <-
      gene_pathway[gene_pathway[,2]=="", 1]

    gene_df <- data.frame(gene = gene_pathway[,1],
                          betas = coef)
    gene_df <- aggregate(. ~ gene, data = gene_df, FUN = sum)
    pathway_df <- data.frame(pathway = gene_pathway[,2],
                          betas = coef)
    pathway_df <- aggregate(. ~ pathway, data = pathway_df, FUN = sum)

    gene_adj_df <- merge(gene_adj_df, gene_df,
                         by = c("gene"), all = TRUE)
    pathway_adj_df <- merge(pathway_adj_df, gene_df,
                            by = c("pathway"), all = TRUE)

  }
  colnames(gene_adj_df) <- c("gene", tissue_names)
  colnames(pathway_adj_df) <- c("pathway", tissue_names)
  rownames(gene_adj_df) <- gene_adj_df$gene
  rownames(pathway_adj_df) <- pathway_adj_df$gene
  gene_adj_df <- gene_adj_df[, -1]
  pathway_adj_df <- pathway_adj_df[, -1]

  return(list(gene = gene_adj_df, pathway = pathway_adj_df))
}



#' Format Coefficient Output from PrediXcan Association Output.
#'
#' Produces an adjacency-matrix-like dataframe where each row is a gene and
#' each column is the estimated effect size of the gene from the PrediXcan
#' association software.
#'
#' @param predixcan_assoc_filenames List of paths to file output from
#' PrediXcan association software.
#' @param tissue_names Names for the tissues/cell/condition types for the
#' imputed gene expression data used to produce each of the above PrediXcan
#' associations.
#'
#' @return Returns an adjacency-matrix-like dataframe of genes and their
#' effect sizes for multiple tissues/cells/conditions.
#'
#' @export
#' @import data.table
#' @import qvalue

predixcan2adj_df <- function(predixcan_assoc_filenames,
                             tissue_names = seq_along(predixcan_assoc_filenames),
                             use_fdr = FALSE,
                             pvalue_thresh = 0.05){
  gene_adj_df <- data.frame(gene = c(NA))
  for (f in predixcan_assoc_filenames){
    assoc <- as.data.frame(data.table::fread(f, header = TRUE))
    if(use_fdr == TRUE){
      assoc$qvalue <- qvalue::qvalue(assoc$pvalue)
      gene_adj_df <- merge(gene_adj_df, assoc[which(assoc$qvalue < pvalue_thresh),
                                              c("gene", "se")],
                           by = c("gene"), sort = FALSE, all = TRUE)
    } else{
      gene_adj_df <- merge(gene_adj_df, assoc[which(assoc$pvalue < pvalue_thresh),
                                              c("gene", "se")],
                           by = c("gene"), sort = FALSE, all = TRUE)
    }
  }
  gene_adj_df <- gene_adj_df[!is.na(gene_adj_df$gene), ]
  colnames(gene_adj_df) <- c("gene", tissue_names)
  rownames(gene_adj_df) <- gene_adj_df$gene
  gene_adj_df <- gene_adj_df[, -1]
  return(gene_adj_df)
}

# [END]
