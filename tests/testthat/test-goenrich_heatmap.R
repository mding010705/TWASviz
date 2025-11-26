if (!requireNamespace("org.Hs.eg.db")){
  install.packages("org.Hs.eg.db")
}

test_that("no error on good input", {
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
