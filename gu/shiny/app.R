#!/usr/bin/env Rscript
# GU Shiny launcher. gu.sh keeps the command, paths and environment settings.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
app_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
if ("--review" %in% commandArgs(trailingOnly = TRUE)) {
  source(file.path(app_dir,"review.R"),local=TRUE)
  quit(status = 0)
}
source(file.path(app_dir,"data.R"),local=TRUE)
source(file.path(app_dir,"ui.R"),local=TRUE)
source(file.path(app_dir,"server.R"),local=TRUE)

app <- shiny::shinyApp(ui, server)
.gu_ui_handler <- app$httpHandler
app$httpHandler <- function(req) {
  response<-gu_igv_reference_response(req,.gu_reference_resources)
  if(is.null(response)).gu_ui_handler(req)else response
}
shiny::runApp(app,
  host = Sys.getenv("GU_SHINY_HOST", "127.0.0.1"),
  port = as.integer(Sys.getenv("GU_SHINY_PORT", "3838")), launch.browser = interactive())
