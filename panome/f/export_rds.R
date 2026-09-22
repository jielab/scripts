# Read-only bridge. No sourcing of phenotype pipelines or installing R packages.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
x <- readRDS(args[1])
if (!is.data.frame(x)) stop("Expected one data.frame in RDS")
cols <- readLines(args[3], warn = FALSE)
if (length(cols)) {
  absent <- setdiff(cols, names(x))
  if (length(absent)) stop("Missing columns: ", paste(absent, collapse = ", "))
  x <- x[, cols, drop = FALSE]
}
for (nm in names(x)) {
  if (inherits(x[[nm]], "Date") || inherits(x[[nm]], "POSIXt"))
    x[[nm]] <- format(x[[nm]], "%Y-%m-%d")
  if (inherits(x[[nm]], "integer64")) x[[nm]] <- as.character(x[[nm]])
}
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::fwrite(x, args[2], na = "")
} else {
  write.csv(x, args[2], row.names = FALSE, na = "")
}
