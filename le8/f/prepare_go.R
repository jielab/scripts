#!/usr/bin/env Rscript
# Export installed annotation snapshots without network access, on local or HPC.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('Usage: Rscript f/prepare_go.R OUTPUT_GO_DIRECTORY')
need <- c('AnnotationDbi', 'org.Hs.eg.db', 'GO.db', 'data.table')
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install annotation packages before export: ', paste(missing, collapse = ', '))
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
out <- file.path(args[1], 'gene_go.tsv')
if (file.exists(out)) stop('Database already exists; use a new snapshot directory: ', out)
db <- org.Hs.eg.db::org.Hs.eg.db
x <- suppressMessages(AnnotationDbi::select(db, keys = AnnotationDbi::keys(db, keytype = 'SYMBOL'),
    keytype = 'SYMBOL', columns = c('GO', 'ONTOLOGY')))
x <- unique(x[!is.na(x$SYMBOL) & !is.na(x$GO) & x$ONTOLOGY %in% c('BP','MF','CC'), ])
y <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = unique(x$GO), keytype = 'GOID', columns = 'TERM'))
x$TERM <- y$TERM[match(x$GO, y$GOID)]
x <- x[order(x$SYMBOL, x$GO, x$ONTOLOGY), ]
tmp <- tempfile(tmpdir = args[1])
data.table::fwrite(x, tmp, sep = '\t', quote = FALSE, na = '')
if (!file.rename(tmp, out)) { unlink(tmp); stop('Cannot publish GO snapshot') }
writeLines(c(paste('R:', R.version.string),
    vapply(need, function(p) paste(p, as.character(utils::packageVersion(p))), character(1)),
    paste('gene_go.tsv md5:', unname(tools::md5sum(out)))), file.path(args[1], 'gene_go.provenance.txt'))
cat('Offline GO snapshot:', out, 'rows:', nrow(x), '\n')
