#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value after ", key, call. = FALSE)
  args[[i + 1L]]
}

infile <- get_arg("--in")
outdir <- get_arg("--out-dir")
pop <- toupper(get_arg("--pop", ""))
if (is.null(infile) || is.null(outdir) || !nzchar(pop)) {
  stop("Required: --in FILE --out-dir DIR --pop POP", call. = FALSE)
}
if (!file.exists(infile)) stop("Missing input: ", infile, call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Read only the header. The full GWAS is streamed once through awk, rather than
# being loaded 22 times (once per chromosome).
hdr <- names(fread(infile, nrows = 0L, showProgress = FALSE))
if (!length(hdr)) stop("Cannot read the header from ", infile, call. = FALSE)
canon <- toupper(gsub("[^A-Z0-9]", "", hdr))
pick_index <- function(label, aliases) {
  idx <- match(toupper(gsub("[^A-Z0-9]", "", aliases)), canon, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (!length(idx)) stop("Cannot identify ", label, " column in ", infile,
                         ". Header: ", paste(hdr, collapse = ", "), call. = FALSE)
  idx[[1L]]
}
idx <- c(
  SNP = pick_index("SNP", c("SNP", "RSID", "MARKERNAME", "VARIANTID", "ID")),
  CHR = pick_index("CHR", c("CHR", "CHROM", "CHROMOSOME")),
  A1 = pick_index("effect allele", c("A1", "EA", "EFFECTALLELE", "ALT", "ALLELE1")),
  A2 = pick_index("other allele", c("A2", "NEA", "OTHERALLELE", "REF", "ALLELE2")),
  BETA = pick_index("BETA", c("BETA", "EFFECT", "LOGOR", "B")),
  SE = pick_index("SE", c("SE", "STDERR", "STANDARDERROR", "BETASE"))
)

awk_file <- tempfile(fileext = ".awk")
on.exit(unlink(awk_file), add = TRUE)
# The summary-statistic files used here are whitespace/tab delimited. AWK's
# default FS therefore handles either form. Each chromosome receives a header.
awk_code <- sprintf('
BEGIN {
  OFS="\\t"
  for (c=1; c<=22; c++) {
    f=OUTDIR "/" POP ".chr" c ".dat"
    print "SNP", "A1", "A2", "BETA", "SE" > f
  }
}
NR==1 {next}
{
  chr=toupper($%d); sub(/^CHR/, "", chr)
  if (chr !~ /^[0-9]+$/) next
  chr += 0
  snp=$%d; a1=toupper($%d); a2=toupper($%d); beta=$%d; se=$%d
  num="^[-+]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eE][-+]?[0-9]+)?$"
  if (chr>=1 && chr<=22 &&
      snp!="" && snp!="NA" &&
      a1 ~ /^[ACGT]$/ && a2 ~ /^[ACGT]$/ && a1!=a2 &&
      beta ~ num && se ~ num && (se+0)>0) {
    key=chr SUBSEP snp
    if (!(key in seen)) {
      seen[key]=1
      f=OUTDIR "/" POP ".chr" chr ".dat"
      print snp, a1, a2, beta, se >> f
    }
  }
}
', idx[["CHR"]], idx[["SNP"]], idx[["A1"]], idx[["A2"]], idx[["BETA"]], idx[["SE"]])
writeLines(awk_code, awk_file)

reader <- if (grepl("\\.(gz|bgz|bgzip)$", infile, ignore.case = TRUE)) "gzip -cd" else "cat"
cmd <- paste(
  "set -o pipefail;",
  reader, shQuote(normalizePath(infile)), "|",
  "awk", "-v", paste0("OUTDIR=", shQuote(normalizePath(outdir))),
  "-v", paste0("POP=", shQuote(pop)),
  "-f", shQuote(awk_file)
)
status <- system2("bash", c("-c", shQuote(cmd)))
if (!identical(status, 0L)) stop("Streaming summary-statistic preparation failed for ", infile, call. = FALSE)

qc <- rbindlist(lapply(1:22, function(chr) {
  f <- file.path(outdir, paste0(pop, ".chr", chr, ".dat"))
  if (!file.exists(f)) stop("Missing output: ", f, call. = FALSE)
  wc <- system2("wc", c("-l", f), stdout = TRUE, stderr = TRUE)
  n_line <- suppressWarnings(as.integer(strsplit(trimws(wc[[1L]]), "[[:space:]]+")[[1L]][1L]))
  if (!is.finite(n_line)) stop("Cannot count lines in ", f, call. = FALSE)
  n <- max(0L, n_line - 1L)
  if (n < 100L) stop("Only ", n, " valid variants for chromosome ", chr, " in ", infile, call. = FALSE)
  data.table(pop = pop, chr = chr, variants = n, file = normalizePath(f))
}))
fwrite(qc, file.path(outdir, paste0(pop, ".sumstats_qc.tsv")), sep = "\t", quote = FALSE)
cat("Prepared", sum(qc$variants), "variants for", pop, "across chromosomes 1-22\n")
