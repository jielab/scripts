args <- commandArgs(trailingOnly=TRUE)
app <- if(length(args)) args[1] else file.path(getwd(), "shiny")
shiny::runApp(app, host=Sys.getenv("GU_SHINY_HOST","127.0.0.1"), port=as.integer(Sys.getenv("GU_SHINY_PORT","3838")), launch.browser=interactive())
