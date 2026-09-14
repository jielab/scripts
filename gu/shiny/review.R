#!/usr/bin/env Rscript
# Standalone Shiny host for the prepared, self-contained evidence review page.
# Uses the existing shiny package only. Does not modify any environment or database.
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  ix <- match(flag, args)
  if (is.na(ix)) return(default)
  if (ix == length(args)) stop("Missing value for ", flag)
  args[[ix + 1L]]
}
data_dir <- value("--data")
if (is.null(data_dir)) stop("Usage: Rscript shiny/app.R --review --data /mnt/d/analysis/gu/final/review [--port 3839] [--host 127.0.0.1]")
port <- suppressWarnings(as.integer(value("--port", "3839")))
if (is.na(port) || port < 1L || port > 65535L) stop("Invalid port")
host <- value("--host", "127.0.0.1")
if (!requireNamespace("shiny", quietly = TRUE)) stop("The existing R session does not have shiny. Activate your working GU R environment; this script will not install packages.")
data_dir <- normalizePath(data_dir, mustWork = TRUE)
if (!file.exists(file.path(data_dir, "review.html"))) stop("review.html missing. Run prepare_review.py first.")
shiny::addResourcePath("gu_review", data_dir)
ui <- shiny::fluidPage(
  shiny::tags$head(shiny::tags$style("html,body{margin:0;padding:0}.container-fluid{padding:0}")),
  shiny::tags$iframe(src = "gu_review/review.html", title = "GU evidence review",
                    style = "display:block;width:100%;height:100vh;border:0")
)
app <- shiny::shinyApp(ui, server = function(input, output, session) {})
shiny::runApp(app, host = host, port = port, launch.browser = FALSE)

