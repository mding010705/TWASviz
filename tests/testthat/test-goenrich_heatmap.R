if (!requireNamespace("org.Hs.eg.db")){
  install.packages("org.Hs.eg.db")
}
library(ggplot2)

test_that("intergration test: no error on good input", {
  expect_no_error(TWASviz::goenrich_heatmap(enrich_res =
                gene_enrichment(gene_set = list(a=c("CHST2", "B3GNT1", "B3GNT2",
                                                    "CHST1","B4GALT1","B3GNT7",
                                                    "ST3GAL2"),
                                                b=c("B4GALT3",
                                                    "ST3GAL1", "B4GALT2",
                                                    "B4GALT4", "CHST4", "ST3GAL3",
                                                    "FUT8", "CHST6"),
                                                d=c("CHST2", "B3GNT1", "B3GNT2",
                                                    "B4GALT3",
                                                    "ST3GAL1", "B4GALT2",
                                                    "B4GALT4", "CHST4", "ST3GAL3",
                                                    "FUT8", "CHST6")),
                                                    organism = "org.Hs.eg.db",
                                                    gene_nom = "SYMBOL")))
})


test_that("error on empty input", {
  expect_error(TWASviz::goenrich_heatmap(enrich_res = list()))
})


# Create a mock enrichResult object
mock_enrich_res <- list(
  Tissue1 = structure(
    list(
      result = data.frame(
        ID = c("GO:1", "GO:2"),
        Description = c("TermA", "TermB"),
        p.adjust = c(0.001, 0.05),
        Count = c(10, 5),
        stringsAsFactors = FALSE
      )
    ),
    class = "enrichResult"
  ),

  Tissue2 = structure(
    list(
      result = data.frame(
        ID = c("GO:3", "GO:4"),
        Description = c("TermC", "TermD"),
        p.adjust = c(0.02, 0.03),
        Count = c(3, 7),
        stringsAsFactors = FALSE
      )
    ),
    class = "enrichResult"
  )
)

test_that("goenrich_heatmap validates inputs", {
  expect_error(goenrich_heatmap("not a list"))
  expect_error(goenrich_heatmap(list(A = "not enrichResult")))
  expect_error(goenrich_heatmap(mock_enrich_res, top_n = 0))
})

test_that("goenrich_heatmap fails gracefully on empty enrichment results", {
  empty_res <- list(
    T1 = structure(list(result = data.frame()), class = "enrichResult")
  )
  expect_error(goenrich_heatmap(empty_res), regexp = "No enrichGO")
})

test_that("goenrich_heatmap returns a ggplot object.", {
  p <- goenrich_heatmap(mock_enrich_res, top_n = 2)

  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% sapply(p$layers, function(x) class(x$geom)[1]))
  expect_true("GeomText" %in% sapply(p$layers, function(x) class(x$geom)[1]))
})

test_that("goenrich_heatmap handles duplicated terms", {
  dupe_res <- mock_enrich_res
  dupe_res$Tissue1@result$Description[2] <- "TermA"  # duplicate name

  p <- goenrich_heatmap(dupe_res, top_n = 2)
  expect_s3_class(p, "ggplot")
})
