# Exercise the production move stage on synthetic coordinates, without network.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
root <- normalizePath(file.path(dirname(script), ".."))
exprs <- parse(file.path(root, "f/phe.R"))
stage <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("run_if")) &&
  identical(x[[2]], "move"), as.list(exprs))[[1]][[3]]

run_case <- function(fail_offline = FALSE, bad_coordinates = FALSE) {
  temp <- tempfile("move-map-test-")
  dir.create(temp); dir.create(file.path(temp, "Rdata"))
  old <- setwd(temp)
  on.exit({setwd(old); unlink(temp, recursive = TRUE)}, add = TRUE)
  env <- new.env(parent = globalenv())
  env$indir <- temp
  env$map_tile_cache <- file.path(temp, "tiles")
  env$center.lst <- c(c1 = "Croydon")
  env$phe <- data.frame(eid = c("a", "b"), tdi = c(2, 4), center = "c1",
    birth_east = c(530000, 531000), birth_north = 180000,
    home_east = c(530000, 531000), home_north = 180000)
  if (bad_coordinates) env$phe$birth_east <- c("invalid", "invalid")
  calls <- logical()
  env$plot_map <- function(..., basemap) {
    stopifnot(file.exists(file.path(temp, "Rdata/move.rds")))
    calls <<- c(calls, basemap)
    list(basemap = basemap)
  }
  env$ggsave <- function(filename, plot, ...) {
    writeLines("partial render", filename)
    if (plot$basemap || fail_offline) stop("Zero tiles were loaded for type osm")
    writeLines("completed synthetic map", filename)
  }
  existing <- file.path(temp, "Fig1_UK_all_centers_tdi.png")
  writeLines("old complete map", existing)
  error <- tryCatch({eval(stage, env); NULL}, error = identity)
  if (bad_coordinates) {
    stopifnot(inherits(error, "error"), !file.exists(file.path(temp, "Rdata/move.rds")))
    return(invisible(NULL))
  }
  stopifnot(is.null(error))
  dat <- readRDS(file.path(temp, "Rdata/move.rds"))
  stopifnot(identical(dat$eid, c("a", "b")), all(dat$move == 0),
    identical(dat$tdi.birth, c(2, 4)), all(dat$tdi.birth.nn == 0),
    identical(calls, c(TRUE, FALSE, FALSE, FALSE)),
    length(list.files(temp, pattern = "^\\.move-map-", all.files = TRUE)) == 0)
  stopifnot(readLines(existing) == if (fail_offline) "old complete map" else "completed synthetic map")
}
run_case()
run_case(fail_offline = TRUE)
run_case(bad_coordinates = TRUE)
cat("PASS: basemap fallback, optional-map failure, scientific-error propagation\n")
