if (!requireNamespace("org.Hs.eg.db")){
  install.packages("org.Hs.eg.db")
}
library(ggplot2)

test_that("Intergration test: no error on good input.", {
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


test_that("Error on empty input.", {
  expect_error(TWASviz::goenrich_heatmap(enrich_res = list()))
})


test_that("goenrich_heatmap requires a list.", {
  expect_error(goenrich_heatmap("not a list"),
               "enrich_res must be a list")
})

test_that("goenrich_heatmap rejects non-enrichResult entries.", {
  bad <- list(A = data.frame(x = 1))
  expect_error(goenrich_heatmap(bad),
               "must be 'enrichResult'")
})

test_that("top_n must be > 0.", {
  expect_error(goenrich_heatmap(list(), top_n = 0),
               "top_n must be a positive integer")
})

test_that("goenrich_heatmap validates inputs", {
  expect_error(goenrich_heatmap("not a list"))
  expect_error(goenrich_heatmap(list(A = "not enrichResult")))
  expect_error(goenrich_heatmap(mock_enrich_res, top_n = 0))
})




