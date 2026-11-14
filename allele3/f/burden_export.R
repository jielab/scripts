#!/usr/bin/env Rscript
# Optional read-only adapter for R classes that pyreadr cannot safely interpret.
# No LE8/source pipeline execution and no package installation. The Python
# assoc_phenotypes.py stage then applies the exact cohort/outcome contract.
# Usage: Rscript assoc_phenotypes.R --input /UKB/Rdata/all.rds --out cohort.tsv.gz
#        [--columns eid,date_attend,...]
args <- commandArgs(trailingOnly = TRUE)
opt <- list()
if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat("Read-only RDS export: --input FILE --out FILE.tsv[.gz] [--columns CSV]\n")
    quit(save = "no", status = 0)
}
if (length(args) %% 2L != 0L) stop("Each option requires a value")
for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% c("input", "out", "columns")) stop("Unknown option: ", args[[i]])
    opt[[key]] <- args[[i + 1L]]
}
if (is.null(opt$input) || is.null(opt$out)) stop("--input and --out are required")
x <- readRDS(opt$input)
if (!is.data.frame(x)) stop("Expected one R data.frame")
if (anyDuplicated(names(x))) stop("Duplicate source columns")
if (!is.null(opt$columns)) {
    wanted <- trimws(strsplit(opt$columns, ",", fixed = TRUE)[[1L]])
    if (length(setdiff(wanted, names(x)))) stop("Requested columns absent")
    x <- x[, wanted, drop = FALSE]
}
for (name in names(x)) {
    if (inherits(x[[name]], "Date") || inherits(x[[name]], "POSIXt"))
        x[[name]] <- format(x[[name]], "%Y-%m-%d")
    if (inherits(x[[name]], "integer64")) {
        if (!requireNamespace("bit64", quietly = TRUE))
            stop("integer64 source needs installed bit64 for lossless text export")
        x[[name]] <- as.character(x[[name]])
    }
    if (is.factor(x[[name]])) x[[name]] <- as.character(x[[name]])
    if (is.list(x[[name]])) stop("List-valued column cannot be exported: ", name)
}
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
tmp <- paste0(opt$out, ".tmp")
con <- if (grepl("[.]gz$", opt$out)) gzfile(tmp, "wt") else file(tmp, "wt")
tryCatch(write.table(x, con, sep = "\t", row.names = FALSE, quote = TRUE, na = "NA"),
         finally = close(con))
if (!file.rename(tmp, opt$out)) stop("Failed to publish phenotype export")
cat("[ASSOC] Exported existing columns; no phenotype definitions or missing values changed.\n")
