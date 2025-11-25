#' Launch Shiny App for TWASviz
#'
#' Launches the Shiny app for TWASviz, providing
#' a user interface to interactively work with the this package.
#' The app's code is located in \code{./inst/shiny-scripts}.
#'
#' @return No return value.
#'
#' @examples
#' \dontrun{
#'   run_TWAS_analysis()
#' }
#'
#' @references
#' Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J, McPherson J, Dipert A, Borges B (2025). _shiny: Web Application Framework for R_. R package version 1.11.1, <https://CRAN.R-project.org/package=shiny>.
#'
#' @export
#' @importFrom shiny runApp

run_TWAS_analysis <- function() {
  appDir <- system.file("shiny-scripts", package = "TWASviz")
  shiny::runApp(appDir, display.mode = "normal")
}
