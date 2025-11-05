test_that("no error on good input", {
  expect_no_error(TWASviz::gene_enrichment(gene_set = list(a=c("ENSG00000106299", "ENSG00000143632",
                                             "ENSG00000241553", "ENSG00000101182",
                                             "ENSG00000162704", "ENSG00000138071",
                                             "ENSG00000061676", "ENSG00000163466",
                                             "ENSG00000136238", "ENSG00000198400",
                                             "ENSG00000158092")),
                                organism = "org.Hs.eg.db"))
})


test_that("no enrichment", {
  expect_equal(length(TWASviz::gene_enrichment(gene_set =
                                                 list(a=c("A")),
                                           organism = "org.Hs.eg.db")), 0)
})

test_that("correct enrichment", {
  expect_equal(TWASviz::gene_enrichment(gene_set =
                                                 list(a=c("ENSG00000131067", "ENSG00000139631",
                                                "ENSG00000099998", "ENSG00000100031",
                                                "ENSG00000128683", "ENSG00000276559",
                                                "ENSG00000136881", "ENSG00000136750",
                                                "ENSG00000167741", "ENSG00000129596",
                                                "ENSG00000181915", "ENSG00000106299",
                                                "ENSG00000143632", "ENSG00000167468",
                                                "ENSG00000241553", "ENSG00000101182",
                                                "ENSG00000162704", "ENSG00000138071",
                                                "ENSG00000061676", "ENSG00000163466",
                                                "ENSG00000136238", "ENSG00000198400",
                                                "ENSG00000158092", "ENSG00000224586",
                                                "ENSG00000117862", "ENSG00000233276",
                                                "ENSG00000176153", "ENSG00000211445",
                                                "ENSG00000182054", "ENSG00000138413",
                                                "ENSG00000178814", "ENSG00000023909",
                                                "ENSG00000167741", "ENSG00000166825",
                                                "ENSG00000131067", "ENSG00000065621",
                                                "ENSG00000099998", "ENSG00000160211",
                                                "ENSG00000115758")),
                                               organism = "org.Hs.eg.db")[["a"]]@result[["ID"]][1],
                                                "GO:0006749")
})
