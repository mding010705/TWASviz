#' Launch Shiny App for TWASviz
#'
#' Launches the Shiny app for TWASviz, providing
#' a user interface to interactively produce volcano plots, overlap heatmaps,
#' and Gene Ontology term enrichment heatmaps this package.
#' The app's code is located in `"./inst/shiny-scripts"`.
#'
#' @return No return value.
#'
#' @examples
#' \dontrun{
#'   run_TWAS_analysis()
#' }
#'
#' @export
#' @importFrom shiny runApp
#'
#' @references Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Xie, Y., Allen, J.,
#' McPherson, J., Dipert, A., Borges, B. (2025). shiny: Web Application Framework for R.
#' R package version 1.11.1. https://CRAN.R-project.org/package=shiny.


run_TWAS_analysis <- function() {
  appDir <- system.file("shiny-scripts", package = "TWASviz")
  shiny::runApp(appDir, display.mode = "normal")
}
