#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop("Usage: gwas_post_mplot.R GWAS OUTPUT TRAIT METHOD PANEL MAGMA_GENES PLOT_F CIS_BED MH_PLOT_BED")
}

input <- args[[1L]]
output <- args[[2L]]
trait <- args[[3L]]
method <- args[[4L]]
panel <- args[[5L]]
magma_file <- args[[6L]]
plot_f <- args[[7L]]
cis_bed <- if (nzchar(args[[8L]])) args[[8L]] else NULL
mh_plot_bed <- if (nzchar(args[[9L]])) args[[9L]] else NULL

if (!file.exists(input) || file.info(input)$size <= 0) stop("GWAS plotting input is missing: ", input)
if (!method %in% c("self", "CMplot", "qqman")) stop("Unsupported plot method: ", method)
if (!panel %in% c("none", "magma")) stop("Unsupported additional panel: ", panel)
if (panel == "magma" && (!file.exists(magma_file) || file.info(magma_file)$size <= 0)) {
  stop("MAGMA panel requested but genes.out is missing: ", magma_file)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
top_png <- tempfile(pattern = paste0(".", basename(output), "."), tmpdir = dirname(output), fileext = ".top.png")
bottom_png <- tempfile(pattern = paste0(".", basename(output), "."), tmpdir = dirname(output), fileext = ".magma.png")
combined_png <- tempfile(pattern = paste0(".", basename(output), "."), tmpdir = dirname(output), fileext = ".combined.png")
on.exit(unlink(c(top_png, bottom_png, combined_png), force = TRUE), add = TRUE)

open_png <- function(path, width, height, res = 300L) {
  a <- list(filename = path, width = width, height = height, units = "px", res = res, bg = "white")
  if (capabilities("cairo")) a$type <- "cairo"
  do.call(grDevices::png, a)
}

read_gwas <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    x <- data.table::fread(path, showProgress = FALSE, data.table = FALSE, check.names = FALSE)
  } else {
    con <- if (grepl("\\.(gz|bgz)$", path, ignore.case = TRUE)) gzfile(path) else path
    x <- utils::read.table(con, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
  }
  names(x) <- toupper(sub("^#", "", names(x)))
  need <- c("CHR", "POS", "P")
  if (!all(need %in% names(x))) stop("GWAS plotting input must contain ", paste(need, collapse = "/"))
  x$CHR <- toupper(sub("^CHR", "", as.character(x$CHR)))
  x$CHR[x$CHR == "X"] <- "23"
  x$CHR[x$CHR == "Y"] <- "24"
  x$CHR[x$CHR %in% c("M", "MT")] <- "25"
  x$CHR <- suppressWarnings(as.integer(x$CHR))
  x$POS <- suppressWarnings(as.numeric(x$POS))
  x$P <- suppressWarnings(as.numeric(x$P))
  positive <- x$P[is.finite(x$P) & x$P > 0]
  p_floor <- if (length(positive)) max(.Machine$double.xmin, min(positive) * 0.1) else .Machine$double.xmin
  x$P[is.finite(x$P) & x$P <= 0] <- p_floor
  x <- x[is.finite(x$CHR) & x$CHR >= 1 & x$CHR <= 24 & is.finite(x$POS) & x$POS > 0 &
           is.finite(x$P) & x$P > 0 & x$P <= 1, , drop = FALSE]
  if (!nrow(x)) stop("No valid CHR/POS/P rows in ", path)
  if (!"SNP" %in% names(x)) x$SNP <- paste0("chr", x$CHR, ":", format(x$POS, scientific = FALSE, trim = TRUE))
  x$SNP <- make.unique(as.character(x$SNP), sep = "_")
  x
}

render_self <- function(path) {
  if (!file.exists(plot_f)) stop("self plot method requires --plot-f: ", plot_f)
  source(plot_f)
  if (!exists("mh_plot", mode = "function")) stop("mh_plot() was not found in ", plot_f)
  cis_gene <- NULL
  if (!is.null(cis_bed) && file.exists(cis_bed)) {
    z <- tryCatch(utils::read.table(cis_bed, header = FALSE, comment.char = "#", stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (!is.null(z) && ncol(z) >= 4L && any(as.character(z[[4L]]) == trait)) cis_gene <- trait
  }
  open_png(path, 4000L, 1900L)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(mh_plot(input, col = c("#4B5563", "#A7B0BD"), cis_gene = cis_gene,
                cis_bed = cis_bed, mh_plot_bed = mh_plot_bed,
                cis_color = "#D62728", other_top_color = "#159947",
                other_top_max_per_chr = 2, locus_size = 1e6, main = trait))
}

render_cmplot <- function(path) {
  if (!requireNamespace("CMplot", quietly = TRUE)) stop("--plot-method CMplot requires the CMplot R package")
  x <- read_gwas(input)
  pmap <- data.frame(SNP = x$SNP, Chromosome = x$CHR, Position = floor(x$POS), P = x$P,
                     stringsAsFactors = FALSE)
  open_png(path, 4000L, 1900L)
  on.exit(grDevices::dev.off(), add = TRUE)
  CMplot::CMplot(pmap, plot.type = "m", LOG10 = TRUE, file.output = FALSE,
                 col = c("#2F6690", "#9FB4C7"), pch = 19, cex = c(0.35, 0.8, 1),
                 points.alpha = 75L, threshold = 5e-8, threshold.col = "#D62728",
                 threshold.lty = 2, threshold.lwd = 1.5, amplify = FALSE,
                 axis.cex = 0.9, lab.cex = 1.1, mar = c(4.5, 6, 3, 2),
                 main = trait, main.cex = 1.25, verbose = FALSE)
}

render_qqman <- function(path) {
  if (!requireNamespace("qqman", quietly = TRUE)) stop("--plot-method qqman requires the qqman R package")
  x <- read_gwas(input)
  q <- data.frame(SNP = x$SNP, CHR = x$CHR, BP = floor(x$POS), P = x$P,
                  stringsAsFactors = FALSE)
  present <- sort(unique(q$CHR))
  chr_names <- c(as.character(1:22), "X", "Y")[present]
  open_png(path, 4000L, 1900L)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 6, 3, 2) + 0.1, mgp = c(3.2, 1, 0), las = 1)
  qqman::manhattan(q, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
                   col = c("#2F6690", "#9FB4C7"), chrlabs = chr_names,
                   suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
                   cex = 0.42, cex.axis = 0.85, main = trait)
}

read_magma <- function(path) {
  x <- utils::read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- toupper(sub("^#", "", names(x)))
  if (!all(c("CHR", "P") %in% names(x))) stop("MAGMA genes.out must contain CHR and P: ", path)
  start_col <- if ("START" %in% names(x)) "START" else if ("POS" %in% names(x)) "POS" else NA_character_
  stop_col <- if ("STOP" %in% names(x)) "STOP" else start_col
  if (is.na(start_col)) stop("MAGMA genes.out must contain START/STOP or POS: ", path)
  gene_col <- intersect(c("GENE", "SYMBOL", "ID"), names(x))
  gene <- if (length(gene_col)) as.character(x[[gene_col[[1L]]]]) else paste0("gene_", seq_len(nrow(x)))
  chr <- toupper(sub("^CHR", "", as.character(x$CHR)))
  chr[chr == "X"] <- "23"; chr[chr == "Y"] <- "24"; chr[chr %in% c("M", "MT")] <- "25"
  out <- data.frame(GENE = gene, CHR = suppressWarnings(as.integer(chr)),
                    START = suppressWarnings(as.numeric(x[[start_col]])),
                    STOP = suppressWarnings(as.numeric(x[[stop_col]])),
                    P = suppressWarnings(as.numeric(x$P)), stringsAsFactors = FALSE)
  positive <- out$P[is.finite(out$P) & out$P > 0]
  p_floor <- if (length(positive)) max(.Machine$double.xmin, min(positive) * 0.1) else .Machine$double.xmin
  out$P[is.finite(out$P) & out$P <= 0] <- p_floor
  out <- out[is.finite(out$CHR) & out$CHR >= 1 & out$CHR <= 24 & is.finite(out$START) &
               is.finite(out$STOP) & is.finite(out$P) & out$P > 0 & out$P <= 1, , drop = FALSE]
  if (!nrow(out)) stop("No valid gene regions in ", path)
  swap <- out$START > out$STOP
  if (any(swap)) {
    z <- out$START[swap]; out$START[swap] <- out$STOP[swap]; out$STOP[swap] <- z
  }
  out
}

render_magma <- function(path) {
  x <- read_magma(magma_file)
  chr_order <- sort(unique(x$CHR))
  chr_len <- vapply(chr_order, function(ch) max(x$STOP[x$CHR == ch], na.rm = TRUE), numeric(1L))
  offsets <- c(0, head(cumsum(chr_len), -1L))
  names(offsets) <- as.character(chr_order)
  x$OFFSET <- offsets[as.character(x$CHR)]
  x$X1 <- x$START + x$OFFSET
  x$X2 <- x$STOP + x$OFFSET
  x$XMID <- (x$X1 + x$X2) / 2
  x$Y <- -log10(x$P)
  total <- sum(chr_len)
  # Preserve the actual interval centers while guaranteeing that very short genes
  # remain visible as a one-pixel horizontal segment at the target resolution.
  visible_half <- total / (2 * 4000)
  half <- pmax((x$X2 - x$X1) / 2, visible_half)
  x$DRAW1 <- x$XMID - half
  x$DRAW2 <- x$XMID + half
  sig <- 0.05 / nrow(x)
  sig_y <- -log10(sig)
  ymax <- max(x$Y, sig_y, na.rm = TRUE)
  ymax <- max(2, ceiling(ymax * 1.10))
  centers <- offsets + chr_len / 2
  labels <- c(as.character(1:22), "X", "Y")[chr_order]
  cols <- rep(c("#3A6EA5", "#A7B7C9"), length.out = length(chr_order))
  names(cols) <- as.character(chr_order)
  seg_col <- unname(cols[as.character(x$CHR)])
  significant <- x$P <= sig

  open_png(path, 4000L, 1200L)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(5, 6, 3.2, 2) + 0.1, mgp = c(3.2, 1, 0), las = 1)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, total), ylim = c(0, ymax), xaxs = "i", yaxs = "i")
  graphics::abline(h = graphics::axTicks(2), col = "#ECEFF3", lwd = 0.8)
  graphics::segments(x$DRAW1[!significant], x$Y[!significant], x$DRAW2[!significant], x$Y[!significant],
                     col = grDevices::adjustcolor(seg_col[!significant], alpha.f = 0.78), lwd = 1.35, lend = 1)
  if (any(significant)) {
    graphics::segments(x$DRAW1[significant], x$Y[significant], x$DRAW2[significant], x$Y[significant],
                       col = "#D62728", lwd = 2.2, lend = 1)
  }
  graphics::abline(h = sig_y, col = "#D62728", lty = 2, lwd = 1.5)
  graphics::axis(1, at = centers, labels = labels, tick = FALSE, cex.axis = 0.85)
  graphics::axis(2, las = 1, cex.axis = 0.85)
  graphics::box(bty = "l", lwd = 1.1)
  graphics::title(main = paste0(trait, " — MAGMA gene regions"), xlab = "Chromosome",
                  ylab = expression(-log[10](italic(P))), cex.main = 1.25, font.main = 2)
  graphics::mtext(sprintf("Bonferroni: %.3g (%s genes)", sig, format(nrow(x), big.mark = ",")),
                  side = 3, line = 0.25, adj = 1, cex = 0.72, col = "#8B1A1A")
  label_n <- min(10L, sum(significant))
  if (label_n > 0L) {
    ii <- head(order(x$P), label_n)
    graphics::text(x$XMID[ii], pmin(ymax * 0.97, x$Y[ii] + ymax * 0.018), labels = x$GENE[ii],
                   cex = 0.58, srt = 30, adj = c(0, 0), col = "#7F1D1D", xpd = NA)
  }
}

switch(method,
       self = render_self(top_png),
       CMplot = render_cmplot(top_png),
       qqman = render_qqman(top_png))
if (!file.exists(top_png) || file.info(top_png)$size <= 0) stop("Top Manhattan panel was not created")

if (panel == "none") {
  if (file.exists(output)) unlink(output)
  if (!file.rename(top_png, output)) stop("Could not install output: ", output)
} else {
  render_magma(bottom_png)
  if (!file.exists(bottom_png) || file.info(bottom_png)$size <= 0) stop("MAGMA panel was not created")
  if (requireNamespace("magick", quietly = TRUE)) {
    images <- c(magick::image_read(top_png), magick::image_read(bottom_png))
    magick::image_write(magick::image_append(images, stack = TRUE), path = combined_png, format = "png")
  } else if (nzchar(Sys.which("convert"))) {
    status <- system2(Sys.which("convert"), c(shQuote(top_png), shQuote(bottom_png), "-append", shQuote(combined_png)))
    if (!identical(status, 0L)) stop("ImageMagick convert failed with exit status ", status)
  } else {
    stop("Combining --add-panel magma requires the R magick package or ImageMagick convert")
  }
  if (!file.exists(combined_png) || file.info(combined_png)$size <= 0) stop("Combined Manhattan plot was not created")
  if (file.exists(output)) unlink(output)
  if (!file.rename(combined_png, output)) stop("Could not install output: ", output)
}

message("mplot: ", output, " [method=", method, ", panel=", panel, "]")
