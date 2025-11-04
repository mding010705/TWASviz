#' Formatting Gene Expression and Pathway Data to Use with Sparse Group Lasso
#'
#' Duplicates genes in your imputed gene expression matrix that appear multiple
#' times in your specified pathways, renames genes to gene_pathway, subsets and
#' reorders your gene expression and phenotype files to match the participant
#' ID order, and regresses out the covariates from your phenotype if you
#' specify family_func = gaussian(link = "identity") or provides offsets from
#' predicted values from glm(phenotype ~ ., data = covars, family = family_func).
#'
#' @param gene_expression_filename Path to text file with a column called FID
#' containing unique participant IDs and columns named with gene names,
#' containing imputed gene expressions.
#' @param pathways_filenames Vector of paths to RDS files containing a list
#' named by pathway, where each index contains a vector of genes in the pathway.
#' @param phenotype_filename Path to text file with a column called FID
#' containing unique participant IDs and a column with the phenotype
#' of interest.
#' @param phenotype_colname Name of the column containing the phenotype of
#' interest in the phenotype file specified above.
#' @param covariates_filename Path to text file with a column called FID
#' containing unique participant IDs and columns with covariates.
#' @param covariates_colnames Vector of names of the columns containing the
#' covariates in the covariate file specified above.
#' @param family_func Family object specifying the type of model you
#' want to use in the sparse group lasso function. See ?stats::family for more
#' details.
#'
#' @return Returns a list(
#' X = gene expression matrix,
#' y = residuals on phenotype ~ covariates if covariates_filename != NULL
#' and family_func = gaussian(link = "identity"), or just vector of phenotype
#' data otherwise,
#' groups = vector of consecutive integers describing the pathway grouping of
#' the genes,
#' gene2pathway = list with genes as names and the
#' length of pathways that contain the gene,
#' FID = the IDs for the participants in order ,
#' offsets = NULL if the no covariates are specified or
#' family_func = gaussian(link = "identity"), vector of predicted values
#' from generalized linear regression using specified covariates otherwise,
#' processed_pathways = list with pathways as names and the genes in each
#' pathway formatted as gene_pathway).
#'
#' @export
#' @import data.table



preprocess_expressions_pathways <- function(gene_expression_filename,
                                             pathways_filenames,
                                             phenotype_filename = NULL,
                                            phenotype_colname = NULL,
                                            covariates_filename = NULL,
                                            covariates_colnames = NULL,
                                            family_func =
                                              gaussian(link = "identity")){
  gene_expr <- data.table::fread(gene_expression_filename, header = TRUE)
  gene_expr <- as.data.frame(gene_expr)
  gene_expr <- gene_expr[order(gene_expr$FID),]
  rownames(gene_expr) <- gene_expr$FID
  gene_expr <- gene_expr[,-(1:2)]

  if (is.null(phenotype_filename)){
    phenotype <- NULL
  } else {
    phenotype <- data.table::fread(phenotype_filename, header = TRUE)
    phenotype <- phenotype[!(duplicated(phenotype$FID)), ]
    phenotype <- as.data.frame(phenotype)
    phenotype <- phenotype[(match(rownames(gene_expr), phenotype$FID)), ]
    phenotype <- phenotype[!is.na(phenotype$FID), ][, phenotype_colname]

  }

  if (is.null(covariates_filename)){
    offsets <- NULL
  } else {
    covars <- data.table::fread(covariates_filename, header = TRUE)
    covars <- covars[!(duplicated(covars$FID)), ]
    covars <- as.data.frame(covars)
    covars <- covars[(match(phenotype$FID, covars$FID)), ]
    covars <- covars[!is.na(covars$FID), ][, covariates_colnames]
    phenotype <- phenotype[(match(rownames(covars$FID), phenotype$FID)), ]
    fit <- glm(phenotype ~ ., data = covars, family = family_func)
    offsets <- predict(fit, covars)
    if(family_func == gaussian(link = "identity")){
      phenotype <- phenotype - offsets
      offsets <- NULL
    }
  }
  gene_expr <- gene_expr[match(phenotype$FID, rownames(covars$FID)), ]
  # pathway processing
  all_genes <- colnames(gene_expr)

  gene_pathways <- list() # key = gene, value = lengths of pathways gene is in
  proc_pathways <- list()
  all_pathways <- list()
  # pathways
  for (fn in pathways_filenames){
    all_pathways <- c(all_pathways, readRDS(fn))
  }

  for (p in seq_along(all_pathways)){
    genes <- intersect(all_pathways[[p]], all_genes)
    if (length(genes) == 0){
      next
    }
    genes_dup <- genes
    for (g in seq_along(genes)){
      if (genes[g] %in% names(gene_pathways)){
        gene_pathways[[genes[g]]] <- c(gene_pathways[[genes[g]]], length(genes))
        names(gene_pathways[[genes[g]]])[length(gene_pathways[[genes[g]]])] <- (names(all_pathways)[p])
      } else {
        gene_pathways[[genes[g]]] <- c(length(genes))
        names(gene_pathways[[genes[g]]]) <- c(names(all_pathways)[p])
      }
      genes_dup[g] <- paste0(genes[g], "_", names(all_pathways)[p])

    }
    proc_pathways[[names(all_pathways)[p]]] <- genes_dup
  }

  # now add singleton genes that don't belong to a pathway
  for (g in all_genes){
    if (!(g %in% names(gene_pathways))){
      gene_pathways[[g]] <- c(1)
      names(gene_pathways[[g]]) <- g
      proc_pathways[[g]] <- paste0(g, "_", g)
    }
  }
  gene_expr_proc <- gene_expr[, match(sub("_.*", "",
                                         unname(unlist(proc_pathways))),
                                colnames(gene_expr))]
  colnames(gene_expr_proc) <- unname(unlist(proc_pathways))
  grouping <- rep(seq_along(proc_pathways), unname(unlist(sapply(proc_pathways, length))))

  return(list(X = as.matrix(gene_expr_proc), y = phenotype, groups = grouping,
              gene2pathway = gene_pathways, FID = rownames(gene_expr_proc),
              offsets = offsets, processed_pathways = proc_pathways))
}


# [END]
