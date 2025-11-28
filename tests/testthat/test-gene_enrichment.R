# Mock enrichGO to avoid requiring real OrgDb database
fake_enrichGO <- function(...){
  # Return an object with class "enrichResult"
  structure(
    list(
      result = data.frame(
        ID = "GO:0000001",
        Description = "Mock GO Term",
        p.adjust = 0.01,
        Count = 5,
        stringsAsFactors = FALSE
      )
    ),
    class = "enrichResult"
  )
}

test_that("gene_enrichment validates input types.", {
  expect_error(gene_enrichment("not a list"))
  expect_error(gene_enrichment(list(A = c("G1")), tissue_types = NULL))

  # incorrect background length
  expect_error(gene_enrichment(
    gene_set = list(A = c("G1")),
    background = list(1:10, 1:10)
  ))
})

test_that("gene_enrichment stops if OrgDb cannot be loaded.", {
  expect_error(
    gene_enrichment(
      gene_set = list(A = c("G1")),
      organism = "nonexistent.db"
    ),
    regexp = "is not installed"
  )
})

test_that("gene_enrichment runs enrichGO for each tissue (mocked).", {
  # Patch requireNamespace and enrichGO
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) TRUE,
    `clusterProfiler::enrichGO` = fake_enrichGO,

    {
      res <- gene_enrichment(
        gene_set = list(
          T1 = c("GENE1", "GENE2"),
          T2 = c("GENE3", "GENE4")
        ),
        organism = "org.Hs.eg.db",
        gene_nom = "SYMBOL"
      )

      expect_equal(length(res), 2)
      expect_true(all(sapply(res, inherits, what = "enrichResult")))
      expect_equal(names(res), c("T1", "T2"))
    }
  )
})

test_that("gene_enrichment supports background list (mocked).", {
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) TRUE,
    `clusterProfiler::enrichGO` = fake_enrichGO,

    {
      gene_set <- list(
        T1 = c("G1", "G2")
      )
      background <- list(
        T1 = c("G1", "G2", "G3", "G4")
      )

      res <- gene_enrichment(
        gene_set = gene_set,
        organism = "org.Hs.eg.db",
        background = background
      )

      expect_true(inherits(res$T1, "enrichResult"))
    }
  )
})

