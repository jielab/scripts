#!/usr/bin/env Rscript
# bplot launcher and session state. Data, plots and FASTA serving live beside it.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L || any(!grepl("^--", args[seq(1L, length(args), 2L)]))) stop("Expected --option value pairs")
opt <- as.list(args[seq(2L, length(args), 2L)])
names(opt) <- gsub("-", "_", sub("^--", "", args[seq(1L, length(args), 2L)]))
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
app_dir <- dirname(normalizePath(script))
for (package in c("data.table", "jsonlite")) if (!requireNamespace(package, quietly = TRUE)) stop("Install R package: ", package)
library(data.table)
data.table::setDTthreads(2L)
source(file.path(app_dir, "data.R"))
if (!is.null(opt$prepare_config)) {
  bplot_prepare(readRDS(opt$prepare_config), opt$build, opt$progress)
  quit(save = "no", status = 0L)
}
for (package in c("shiny", "plotly", "htmlwidgets", "processx", "bslib"))
  if (!requireNamespace(package, quietly = TRUE)) stop("Install R package: ", package)
library(shiny)
source(file.path(app_dir, "plot.R"))
source(file.path(app_dir, "ref.R"))
source(file.path(app_dir, "ld.R"))
manifest <- jsonlite::read_json(opt$manifest, simplifyVector = TRUE)
opt$source_build <- as.character(manifest$source_build)
opt$tracks <- manifest$tracks
for (path in opt$tracks$file[grepl("[.]thin[.]gz$", opt$tracks$file)]) {
  source_path <- sub("[.]thin[.]gz$", ".gz", path)
  if (!file.exists(path) || !file.exists(paste0(path, ".tbi")))
    stop("Thin GWAS or index missing; regenerate with gwas_format.sh thin: ", path)
  if (file.exists(source_path) && file.info(path)$mtime < file.info(source_path)$mtime)
    stop("Thin GWAS is stale; regenerate with gwas_format.sh thin: ", path)
}
opt$app_dir <- app_dir
opt$max_points <- as.integer(opt$max_points)
opt$ld_max_snps <- as.integer(opt$ld_max_snps)
opt$port <- as.integer(opt$port)
opt$p_threshold <- as.numeric(opt$p_threshold)
opt$tracks$label <- opt$tracks$trait
opt$tracks$race <- toupper(sub("^.*\\.", "", opt$tracks$trait))
for (field in c("labels", "races")) if (!is.null(opt[[field]])) {
  values <- trimws(strsplit(opt[[field]], ",", fixed = TRUE)[[1]])
  if (length(values) != nrow(opt$tracks) || any(!nzchar(values))) stop("--", field, " must match input count")
  opt$tracks[[if (field == "labels") "label" else "race"]] <- values
}
opt$tracks$race <- toupper(opt$tracks$race)
if (any(!grepl("^[A-Z0-9_-]+$", opt$tracks$race))) stop("Invalid ancestry suffix; use --races CSV")
if (!opt$grch %in% c("37", "38")) stop("--grch must be 37 or 38")
dir.create(file.path(opt$output_dir, "cache"), recursive = TRUE, showWarnings = FALSE)
run_dir <- file.path(opt$output_dir, "logs", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid()))
dir.create(run_dir, recursive = TRUE)
config <- file.path(run_dir, "config.rds")
bplot_atomic_rds(opt, config)

assets <- file.path(app_dir, "www")
addResourcePath("bplot_assets", assets)
refs <- bplot_local_igv_references(opt$reference_dir)
# Small annotation files are cached; large indexed references are served in place.
genes <- file.path(opt$output_dir, "cache", "genes")
dir.create(genes, showWarnings = FALSE)
for (build in c("37", "38")) {
  src <- file.path(opt$gene_dir, paste0("glist.", build, ".bed"))
  if (!file.exists(src)) next
  d <- fread(src, header = FALSE, showProgress = FALSE)
  if (ncol(d) < 4L) next
  d <- d[, 1:4]
  d <- d[bplot_chr(d[[1]]) %in% names(bplot_lengths(build))]
  d[[1]] <- paste0("chr", bplot_chr(d[[1]]))
  fwrite(d, file.path(genes, paste0(build, ".bed")), sep = "\t", col.names = FALSE, quote = FALSE)
}
addResourcePath("bplot_genes", genes)
source(file.path(app_dir, "ui.R"))
source(file.path(app_dir, "server.R"))

app <- shinyApp(ui, server)
handler <- app$httpHandler
app$httpHandler <- function(req) {
  response <- bplot_igv_reference_response(req, refs)
  if (is.null(response)) handler(req) else response
}
cat(sprintf("bplot: http://%s:%d · output: %s\n", opt$host, opt$port, opt$output_dir))
runApp(app, host = opt$host, port = opt$port, launch.browser = identical(opt$launch_browser, "TRUE"))
