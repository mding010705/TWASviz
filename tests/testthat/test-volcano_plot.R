library(ggplot2)

test_that("No error on good input.", {
  expect_no_error(TWASviz::volcano_plot(data.frame(beta = c(1, 2, 3),
                                                   p_value = c(0.5, 0.001, 1)),
                                        add_sig_gene_labels = FALSE,
                                        effect_size_colname = "beta",
                                        pvalue_colname = "p_value"))
})


test_that("Error on empty input.", {
  expect_error(TWASviz::volcano_plot())
})



# Input validation

test_that("volcano_plot returns empty ggplot if betas_p_vals is NULL.", {
  p <- volcano_plot(NULL)
  expect_s3_class(p, "ggplot")
  expect_equal(length(p$layers), 0)  # empty plot
})

test_that("volcano_plot errors on non-data.frame input.", {
  expect_error(volcano_plot(1), "data.frame")
  expect_error(volcano_plot(list(a=1)), "data.frame")
})

test_that("volcano_plot checks for pvalue column.", {
  df <- data.frame(zscore = 1:3)
  expect_error(
    volcano_plot(df, pvalue_colname = "pvalue"),
    "must contain your specified pvalue column"
  )
})

test_that("volcano_plot checks for effect size column.", {
  df <- data.frame(pvalue = c(0.1, 0.2))
  expect_error(
    volcano_plot(df, effect_size_colname = "zscore"),
    "must contain your specified effect size column"
  )
})

# No significant genes

test_that("volcano_plot returns message string if no significant genes.", {
  df <- data.frame(
    gene   = c("A","B","C"),
    zscore = c(1,2,3),
    pvalue = c(0.2, 0.3, 0.5)  # none < .05
  )
  out <- volcano_plot(df)
  expect_type(out, "character")
  expect_match(out, "No significant gene associations")
})

# Gene name handling

test_that("volcano_plot uses rownames as gene names if missing gene column.", {
  df <- data.frame(
    zscore = c(2, -1, 3),
    pvalue = c(0.01, 0.2, 0.001)
  )
  rownames(df) <- c("GeneA","GeneB","GeneC")

  p <- volcano_plot(df, add_sig_gene_labels = FALSE)

  # look at plot data
  d <- layer_data(p)

  expect_true("label" %in% names(d))
  expect_setequal(d$label, rownames(df))
})

# Check that output is a ggplot object w/ correct layers

test_that("volcano_plot returns a ggplot object with expected layers.", {
  df <- data.frame(
    gene   = c("A","B","C"),
    zscore = c(2, 0.5, -3),
    pvalue = c(0.01, 0.9, 0.02)
  )

  p <- volcano_plot(df, add_sig_gene_labels = TRUE)

  expect_s3_class(p, "ggplot")

  # Layers: points, vertical/horizontal lines, repelled labels
  layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))

  expect_true("GeomPoint"     %in% layer_geoms)
  expect_true("GeomVline"     %in% layer_geoms)
  expect_true("GeomHline"     %in% layer_geoms)
  expect_true("GeomTextRepel" %in% layer_geoms)
})



# Test significance colouring logic

test_that("volcano_plot colors significant vs non-significant genes correctly.", {
  df <- data.frame(
    gene   = c("A","B","C"),
    zscore = c(1,2,3),
    pvalue = c(0.001, 0.2, 0.04)  # A & C significant
  )

  p <- volcano_plot(df, add_sig_gene_labels = FALSE)

  d <- layer_data(p)  # extract plot-mapped data

  expect_setequal(unique(d$colour), c("blue","black"))
})

# Labeling behavior

test_that("volcano_plot adds labels only for significant genes.", {
  df <- data.frame(
    gene   = c("A","B","C"),
    zscore = c(1,3,2),
    pvalue = c(0.001, 0.5, 0.02)   # A and C significant
  )

  p <- volcano_plot(df, add_sig_gene_labels = TRUE)

  layer_geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))

  expect_true("GeomTextRepel" %in% layer_geoms)

  # Extract label layer data
  repel_idx <- which(layer_geoms == "GeomTextRepel")
  repel_data <- p$layers[[repel_idx]]$data

  expect_setequal(repel_data$gene, c("A","C"))
})

