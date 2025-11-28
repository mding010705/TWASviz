test_that("No error on good input.", {
  expect_no_error(TWASviz::preprocess_expressions_pathways(
    gene_expr = data.frame(FID = c(1,2,3), IID = c(10,11,12), gene1 = 1:3,
                                 gene2 = 1:3),
    all_pathways = list(p1 = c("gene1", "gene2"), p2 = c("gene4"))))
})


test_that("Error on empty input.", {
  expect_error(TWASviz::preprocess_expressions_pathways(
    gene_expr = data.frame(FID = c(), IID = c()),
    all_pathways = list(p1 = c(), p2 = c())))
})


# Write temporary phenotype / covariate files for testing

write_temp_file <- function(df) {
  tf <- tempfile(fileext = ".txt")
  data.table::fwrite(df, tf, sep = "\t")
  return(tf)
}


# Basic gene expression + pathways

gene_expr_base <- data.frame(
  FID = 1:4, IID = c(4,1,3,2),
  G1 = 1:4, G2 = 5:8, G3 = 9:12
)

pathways_base <- list(
  PW1 = c("G1", "G2"),
  PW2 = c("G2", "G3")
)

test_that("Input validation: missing gene_expr fails.", {
  expect_error(preprocess_expressions_pathways(all_pathways = pathways_base))
})

test_that("Input validation: gene_expr must have FID/IID columns.", {
  bad <- gene_expr_base
  colnames(bad)[1:2] <- c("X", "Y")
  expect_error(
    preprocess_expressions_pathways(
      gene_expr = bad,
      all_pathways = pathways_base
    )
  )
})

test_that("Input validation: all_pathways must be a named list.", {
  expect_error(
    preprocess_expressions_pathways(
      gene_expr = gene_expr_base,
      all_pathways = list(c("G1", "G2"))
    )
  )
})

test_that("Input validation: phenotype_colname required when phenotype
          file provided.", {
  pheno_file <- write_temp_file(data.frame(IID=1:4, P=1:4))
  expect_error(
    preprocess_expressions_pathways(
      gene_expr = gene_expr_base,
      all_pathways = pathways_base,
      phenotype_filename = pheno_file
    )
  )
})

test_that("Input validation: unsupported family raises error.", {
  expect_error(
    preprocess_expressions_pathways(
      gene_expr = gene_expr_base,
      all_pathways = pathways_base,
      family_func = "poisson"
    )
  )
})


# Core functionality: no phenotype, no covariates

test_that("Expression matrix is sorted by IID and FID/IID removed.", {
  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base
  )

  expect_equal(rownames(out$X), as.character(sort(gene_expr_base$IID)))
  expect_equal(ncol(out$X), length(unlist(out$processed_pathways)))
})

test_that("Pathways produce duplicated genes with suffix gene_pathway.", {
  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base
  )

  # PW1 should produce G1_PW1, G2_PW1
  expect_true(all(c("G1_PW1", "G2_PW1") %in% colnames(out$X)))

  # Pathways with overlapping G2 = G2_PW1 + G2_PW2
  expect_true(all(c("G2_PW1", "G2_PW2") %in% colnames(out$X)))
})

test_that("Group vector length equals number of duplicated genes.", {
  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base
  )

  expect_equal(length(out$groups), ncol(out$X))
})

test_that("Gene to pathway mapping created correctly.", {
  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base
  )

  # G2 appears in PW1 and PW2
  expect_true("G2" %in% names(out$gene2pathway))
  expect_true(all(c("PW1", "PW2") %in% names(out$gene2pathway$G2)))
})

# Phenotype alignment (Gaussian family)

test_that("Phenotype aligned correctly and residualized for Gaussian.", {

  pheno <- data.frame(IID=c(1,2,3,4), Y=c(10,20,30,40))
  pheno_file <- write_temp_file(pheno)

  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base,
    phenotype_filename = pheno_file,
    phenotype_colname = "Y",
    family_func = "gaussian"
  )

  # Y should be reordered according to sorted IID: 1,2,3,4
  expect_equal(out$IID, c("1","2","3","4"))
  expect_equal(length(out$y), 4)

  # Gaussian: offsets == NULL
  expect_null(out$offsets)
})


# Covariate alignment + GLM fitting (Gaussian)

test_that("Gaussian family residualizes phenotype after GLM covariate
          regression.", {

  pheno <- data.frame(IID=1:4, Y=c(10,20,30,40))
  covars <- data.frame(IID=1:4, age=c(20,30,40,50))
  pheno_file <- write_temp_file(pheno)
  covar_file <- write_temp_file(covars)

  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base,
    phenotype_filename = pheno_file,
    phenotype_colname = "Y",
    covariates_filename = covar_file,
    covariates_colnames = "age",
    family_func = "gaussian"
  )

  # GLM residuals should have mean close to 0
  expect_true(abs(mean(out$y)) < 1e-6)

  # offsets must be NULL for Gaussian
  expect_null(out$offsets)
})


# Binomial family: offsets used and phenotype converted to factor

test_that("Binomial family produces offsets and y is factor.", {

  pheno <- data.frame(IID=1:4, Y=c(0,1,0,1))
  covars <- data.frame(IID=1:4, age=c(1,2,3,4))
  pheno_file <- write_temp_file(pheno)
  covar_file <- write_temp_file(covars)

  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_base,
    phenotype_filename = pheno_file,
    phenotype_colname = "Y",
    covariates_filename = covar_file,
    covariates_colnames = "age",
    family_func = "binomial"
  )

  expect_true(is.factor(out$y))
  expect_true(length(out$offsets) == 4)
})


# Genes not in pathways become single-gene pseudo-pathways

test_that("Genes not in pathways get singleton pathway entries.", {
  pathways_partial <- list(PW1 = c("G1"))

  out <- preprocess_expressions_pathways(
    gene_expr_base,
    pathways_partial
  )

  # G2 and G3 should get own pathway G2/G3
  expect_true("G2" %in% names(out$processed_pathways))
  expect_true("G3" %in% names(out$processed_pathways))

  expect_true("G2_G2" %in% colnames(out$X))
  expect_true("G3_G3" %in% colnames(out$X))
})

