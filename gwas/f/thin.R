#!/usr/bin/env Rscript
# Independent display sidecar; the standardized source is never rewritten.
thin_fun <- function(path) {
  expr <- parse(path)
  hit <- vapply(expr, function(e) is.call(e) && is.symbol(e[[1]]) && as.character(e[[1]]) %in% c("<-", "=") &&
    is.symbol(e[[2]]) && as.character(e[[2]]) == "thinP0", logical(1))
  if (sum(hit) != 1L) stop("Expected one thinP0 definition in ", path)
  env <- new.env(parent = globalenv())
  env$is.data.table <- data.table::is.data.table
  eval(expr[[which(hit)]], envir = env)
  env$thinP0
}

thin_rows <- function(data, limit, fun) {
  data <- data.table::copy(data)
  data[, .chr := toupper(sub("^chr", "", as.character(CHR), ignore.case = TRUE))]
  data[.chr == "M", .chr := "MT"]
  data[, .chr := match(.chr, c(as.character(1:22), "X", "Y", "MT"))]
  # Standardized GWAS also use numeric 23/24/25.
  data[is.na(.chr), .chr := suppressWarnings(as.integer(sub("^chr", "", CHR, ignore.case = TRUE)))]
  data[, .pv := suppressWarnings(as.numeric(P))]
  data[, .score := -log10(pmax(.pv, 1e-300))]
  if ("LOG10P" %in% names(data)) {
    lp <- suppressWarnings(as.numeric(data$LOG10P))
    keep <- is.finite(lp) & as.numeric(data$P) <= 1e-300 & lp >= 300
    data$.score[keep] <- lp[keep]
  }
  data.table::setorderv(data, c(".chr", "POS", "SNP"))
  audit <- list(); out <- list()
  for (chr in unique(data$.chr)) {
    part <- data[.chr == chr]
    selected <- data.table::as.data.table(fun(part, P = 1e-3, p_col = "P", seed = 1234L + chr))
    signal <- selected[.pv <= 1e-3]
    bg <- selected[.pv > 1e-3]
    data.table::setorderv(signal, c(".score", "POS", "SNP"), c(-1L, 1L, 1L))
    lost <- max(0L, nrow(signal) - limit)
    signal <- head(signal, limit)
    room <- limit - nrow(signal)
    if (nrow(bg) > room) bg <- if (room) bg[unique(round(seq(1, .N, length.out = room)))] else bg[0]
    kept <- data.table::rbindlist(list(signal, bg))
    stopifnot(nrow(kept) <= limit)
    audit[[length(audit) + 1L]] <- data.table::data.table(CHR = chr, candidates = nrow(part),
      before_cap = nrow(selected), retained = nrow(kept), background = nrow(bg), signal_omitted_by_cap = lost)
    out[[length(out) + 1L]] <- kept
  }
  result <- data.table::rbindlist(out)
  if (nrow(result)) data.table::setorderv(result, c(".chr", "POS", "SNP"))
  result[, c(".chr", ".score", ".pv") := NULL]
  list(data = result, audit = data.table::rbindlist(audit))
}

thin_main <- function(args) {
  if (length(args) %% 2L) stop("Expected --option value pairs")
  opt <- as.list(args[seq(2L, length(args), 2L)])
  names(opt) <- gsub("-", "_", sub("^--", "", args[seq(1L, length(args), 2L)]))
  if (is.null(opt$hm3_file) || !nzchar(opt$hm3_file)) stop("--hm3-file must specify the HM3 reference list")
  library(data.table)
  setDTthreads(2L)
  limit <- as.integer(opt$chr_max)
  if (!is.finite(limit) || limit < 1L) stop("--chr-max must be positive")
  if (normalizePath(opt$input) == normalizePath(opt$output, mustWork = FALSE)) stop("Thin output must differ from input")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
  awk <- file.path(dirname(normalizePath(script)), "thin.awk")
  identity <- function(path) {
    if (is.null(path) || !nzchar(path)) return(NULL)
    info <- file.info(path); if (is.na(info$size)) stop("Missing file: ", path)
    list(path = normalizePath(path), size = info$size, mtime = as.numeric(info$mtime))
  }
  fun <- thin_fun(opt$phe_r)
  # Without persistent metadata, regenerate to honor all current settings.
  source_identity <- identity(opt$input)
  dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
  work <- tempfile("thin-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  cmd <- paste("gzip -cd --", shQuote(opt$input), "| awk -v", shQuote(paste0("hm3=", opt$hm3_file)),
    "-v", shQuote(paste0("hm3pos=", if (is.null(opt$hm3_pos)) "" else opt$hm3_pos)), "-f", shQuote(awk))
  candidates <- file.path(work, "candidates.tsv")
  rc <- system2("bash", c("-o", "pipefail", "-c", shQuote(cmd)), stdout = candidates)
  if (rc != 0L) stop("Cannot read thin candidates from ", opt$input)
  data <- fread(candidates, showProgress = FALSE)
  if (!nrow(data)) stop("No valid thin candidates in ", opt$input)
  data[, POS := as.numeric(POS)]
  result <- thin_rows(data, limit, fun)
  # tabix requires decimal integer coordinates. fwrite may serialize numeric
  # positions such as 105000000 as 1.05e+08, which tabix reads as position 1.
  result$data[, POS := sprintf("%.0f", POS)]
  tsv <- file.path(work, "thin.tsv")
  fwrite(result$data, tsv, sep = "\t", quote = FALSE, na = "NA")
  stage <- paste0(opt$output, ".", Sys.getpid(), ".tmp")
  on.exit(unlink(c(stage, paste0(stage, ".tbi"))), add = TRUE)
  if (system2("bgzip", c("-c", shQuote(tsv)), stdout = stage) != 0L) stop("bgzip failed")
  cols <- match(c("CHR", "POS"), names(result$data))
  if (system2("tabix", c("-f", "-S", "1", "-s", cols[1], "-b", cols[2], "-e", cols[2], shQuote(stage))) != 0L)
    stop("Indexing thin GWAS failed")
  if (!identical(source_identity, identity(opt$input))) stop("Source changed while thinning; retry: ", opt$input)
  if (!file.rename(stage, opt$output) || !file.rename(paste0(stage, ".tbi"), paste0(opt$output, ".tbi"))) stop("Cannot publish thin GWAS")
  # Remove sidecars left by earlier versions after successful publication.
  unlink(paste0(opt$output, c(".grch", ".n.tsv", ".meta.rds")))
  cat(sprintf("Thin GWAS: %s candidates -> %s SNPs; at most %d per chromosome; %s\n",
    nrow(data), nrow(result$data), limit, opt$output))
}
if (sys.nframe() == 0L) thin_main(commandArgs(trailingOnly = TRUE))
