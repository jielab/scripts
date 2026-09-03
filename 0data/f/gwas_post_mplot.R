#!/usr/bin/env Rscript

# gwas_post.sh-specific CLI adapter for the reusable mh_plot() implementation
# supplied by PLOT_F (normally /mnt/d/scripts/0f/mplot.f.R).

main <- function() {
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(11L, 16L, 19L)) {
  stop(paste("Usage: gwas_post_mplot.R GWAS OUTPUT TRAIT SELF PANEL MAGMA_GENES",
             "PLOT_F CIS_BED MH_PLOT_BED GRCH FLAG_OUTPUT [ADD_SIGNAL MATCH_COL",
             "MATCH_VALUE LOCUS_POS DISPLAY_COL [WRITE_SIG SIG_OUTPUT COJO_FILE]]"))
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
grch <- args[[10L]]
flag_output <- args[[11L]]
add_signal <- if (length(args) >= 12L && nzchar(args[[12L]])) args[[12L]] else NULL
match_col <- if (length(args) >= 13L && nzchar(args[[13L]])) args[[13L]] else "Protein"
match_value <- if (length(args) >= 14L && nzchar(args[[14L]])) args[[14L]] else "[trait]"
locus_pos <- if (length(args) >= 15L && nzchar(args[[15L]])) args[[15L]] else "Gene_chr,Gene_Start,Gene_end"
display_col <- if (length(args) >= 16L && nzchar(args[[16L]])) args[[16L]] else "Gene,beta,p-value"
write_sig <- if (length(args) >= 17L) toupper(args[[17L]]) else "FALSE"
sig_output <- if (length(args) >= 18L && nzchar(args[[18L]])) args[[18L]] else sub("\\.png$", ".sig.txt", output, ignore.case = TRUE)
cojo_file <- if (length(args) >= 19L && nzchar(args[[19L]])) args[[19L]] else ""

if (!file.exists(input) || file.info(input)$size <= 0) stop("GWAS plotting input is missing: ", input)
if (!identical(method, "self")) stop("The only supported plot method is self")
if (!panel %in% c("none", "magma")) stop("Unsupported additional panel: ", panel)
if (!write_sig %in% c("TRUE", "FALSE")) stop("WRITE_SIG must be TRUE or FALSE")
write_sig <- identical(write_sig, "TRUE")
if (write_sig && !nzchar(sig_output)) stop("SIG_OUTPUT is required when WRITE_SIG=TRUE")
if (panel == "magma" && (!file.exists(magma_file) || file.info(magma_file)$size <= 0)) {
  stop("MAGMA panel requested but genes.out is missing: ", magma_file)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
top_png <- tempfile(pattern = paste0(".", basename(output), "."), tmpdir = dirname(output), fileext = ".top.png")
sig_tmp <- if (write_sig) tempfile(pattern = paste0(".", basename(sig_output), "."),
                                   tmpdir = dirname(sig_output), fileext = ".tmp") else NULL
display_input <- NULL
on.exit({
  unlink(top_png, force = TRUE)
  if (!is.null(sig_tmp)) unlink(sig_tmp, force = TRUE)
  if (!is.null(display_input)) unlink(display_input, force = TRUE)
}, add = TRUE)

open_png <- function(path, width, height) {
  a <- list(filename = path, width = width, height = height, units = "px", res = 300L, bg = "white")
  if (capabilities("cairo")) a$type <- "cairo"
  do.call(grDevices::png, a)
}

render_self <- function(path) {
  if (!file.exists(plot_f)) stop("Plot function file is missing: ", plot_f)
  source(plot_f)
  if (!exists("mh_plot", mode = "function")) stop("mh_plot() was not found in ", plot_f)

  # Use the shared phenotype helper's display-only P-value compression.  Keep a
  # local fallback so a user-supplied PLOT_F remains portable when 0phe.f.R is
  # not installed next to it.
  phe_f <- file.path(dirname(plot_f), "0phe.f.R")
  if (file.exists(phe_f)) source(phe_f, local = TRUE)
  # thinP0() in 0phe.f.R accepts both data.frames and data.tables but calls
  # is.data.table() without a namespace.  Keep it usable in this lean Rscript
  # process without attaching the whole data.table package.
  if (!exists("is.data.table", mode = "function")) {
    is.data.table <- function(x) inherits(x, "data.table")
  }
  if (!exists("thinP1", mode = "function")) {
    thinP1 <- function(p, min_p = 1e-300, max_log10 = 50, compress_threshold = 10) {
      p <- as.numeric(p)
      p[!is.na(p) & (p < min_p | p == 0)] <- min_p
      log10p <- -log10(p)
      current_max <- max(log10p[is.finite(log10p)], na.rm = TRUE)
      if (current_max <= max_log10) return(p)
      scaling_factor <- (max_log10 - compress_threshold) / (current_max - compress_threshold)
      log10p_new <- ifelse(log10p <= compress_threshold, log10p,
                           compress_threshold + (log10p - compress_threshold) * scaling_factor)
      10^(-log10p_new)
    }
  }
  if (!exists("thinP0", mode = "function")) {
    thinP0 <- function(dat1, P = 1e-3, p_col = "P", base_keep = 0.10,
                       decade_factor = 0.5, min_keep = 0.01, seed = 1234) {
      set.seed(seed)
      p <- suppressWarnings(as.numeric(dat1[[p_col]]))
      ok <- !is.na(p) & p > 0 & p <= 1
      keep <- rep(TRUE, length(p))
      idx <- which(ok & p > P)
      if (length(idx)) {
        decade <- floor(log10(p[idx] / P))
        prob <- pmax(base_keep * (decade_factor ^ decade), min_keep)
        keep[idx] <- runif(length(idx)) < prob
      }
      dat1[keep, , drop = FALSE]
    }
  }

  read_gwas <- function(x) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fread(x, showProgress = FALSE, data.table = FALSE)
    } else {
      con <- if (grepl("\\.gz$", x, ignore.case = TRUE)) gzfile(x) else x
      on.exit(close(con), add = TRUE)
      utils::read.table(con, header = TRUE, sep = "\t", quote = "", comment.char = "",
                        check.names = FALSE)
    }
  }
  write_gwas <- function(x, path) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
    } else {
      con <- gzfile(path, open = "wt")
      on.exit(close(con), add = TRUE)
      utils::write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }
  }

  gwas <- read_gwas(input)
  if (!all(c("CHR", "POS", "P") %in% names(gwas))) {
    stop("GWAS plotting input must contain CHR, POS, and P: ", input)
  }
  raw_p <- suppressWarnings(as.numeric(gwas$P))
  plot_gwas <- thinP0(gwas, P = 1e-3, p_col = "P")
  thinned_n <- nrow(gwas) - nrow(plot_gwas)
  plot_raw_p <- suppressWarnings(as.numeric(plot_gwas$P))
  plot_p <- thinP1(plot_raw_p, min_p = 1e-300, max_log10 = 50, compress_threshold = 10)
  raw_logp <- -log10(pmax(raw_p, 1e-300))
  raw_logp[!is.finite(raw_p) | raw_p > 1] <- NA_real_
  raw_max <- if (any(is.finite(raw_logp))) max(raw_logp, na.rm = TRUE) else NA_real_
  extreme_n <- sum(is.finite(raw_logp) & raw_logp > 50)
  changed <- (is.na(plot_raw_p) != is.na(plot_p)) |
             (!is.na(plot_raw_p) & !is.na(plot_p) & plot_raw_p != plot_p)
  plot_input <- input
  if (thinned_n > 0L || any(changed)) {
    plot_gwas$P <- plot_p
    display_input <<- tempfile(pattern = paste0(".", trait, ".thinP1."),
                               tmpdir = dirname(output), fileext = ".gz")
    write_gwas(plot_gwas, display_input)
    plot_input <- display_input
    message(sprintf("mplot: thinP0 removed %s of %s display rows above P=1e-3; thinP1 compressed %s P values (%s above -log10(P)=50; raw max %.2f -> display max %.2f)",
                    format(thinned_n, big.mark = ","), format(nrow(gwas), big.mark = ","),
                    format(sum(changed), big.mark = ","), format(extreme_n, big.mark = ","), raw_max,
                    max(-log10(plot_p[is.finite(plot_p) & plot_p > 0]), na.rm = TRUE)))
  }
  normalize_chr <- function(x) {
    x <- toupper(sub("^CHR", "", as.character(x)))
    x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x %in% c("M", "MT")] <- "25"
    suppressWarnings(as.integer(x))
  }
  display_chr <- function(x) {
    x <- as.character(x)
    x[x == "23"] <- "X"; x[x == "24"] <- "Y"; x[x == "25"] <- "MT"
    paste0("chr", x)
  }

  split_columns <- function(x, arg_name) {
    out <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
    if (!length(out) || any(!nzchar(out))) stop(arg_name, " must be a comma-separated column list")
    out
  }
  resolve_columns <- function(requested, available, arg_name) {
    canonical <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
    vapply(requested, function(wanted) {
      exact <- which(available == wanted)
      if (length(exact) == 1L) return(available[[exact]])
      relaxed <- which(canonical(available) == canonical(wanted))
      if (length(relaxed) != 1L) {
        stop(arg_name, " column '", wanted, "' was not found uniquely in ", add_signal)
      }
      available[[relaxed]]
    }, character(1L), USE.NAMES = FALSE)
  }
  format_signal_value <- function(x) {
    if (length(x) == 0L || is.na(x)) return("NA")
    if (is.numeric(x)) {
      if (!is.finite(x)) return("NA")
      if (x != 0 && (abs(x) < 1e-3 || abs(x) >= 1e4)) {
        return(format(x, scientific = TRUE, digits = 4L, trim = TRUE))
      }
      return(format(signif(x, 5L), scientific = FALSE, trim = TRUE))
    }
    as.character(x)
  }
  format_report_value <- function(x) {
    if (length(x) == 0L || is.na(x)) return("NA")
    if (is.numeric(x)) {
      if (!is.finite(x)) return("NA")
      return(sprintf("%.15g", x))
    }
    as.character(x)
  }

  signal_rows <- NULL
  if (!is.null(add_signal)) {
    if (!file.exists(add_signal) || file.info(add_signal)$size <= 0) {
      stop("Signal table is missing: ", add_signal)
    }
    signal <- read_gwas(add_signal)
    match_name <- resolve_columns(match_col, names(signal), "--match-col")
    locus_requested <- split_columns(locus_pos, "--locus-pos")
    if (length(locus_requested) != 3L) {
      stop("--locus-pos must specify chromosome,start,end")
    }
    locus_names <- resolve_columns(locus_requested, names(signal), "--locus-pos")
    display_requested <- split_columns(display_col, "--display-col")
    display_names <- resolve_columns(display_requested, names(signal), "--display-col")
    target <- if (match_value %in% c("[trait]", "{trait}")) trait else match_value
    selected <- signal[as.character(signal[[match_name]]) == target, , drop = FALSE]
    if (nrow(selected)) {
      signal_rows <- data.frame(
        CHR = normalize_chr(selected[[locus_names[[1L]]]]),
        START = suppressWarnings(as.numeric(selected[[locus_names[[2L]]]])),
        END = suppressWarnings(as.numeric(selected[[locus_names[[3L]]]])),
        stringsAsFactors = FALSE
      )
      labels <- vapply(seq_len(nrow(selected)), function(i) {
        values <- vapply(display_names, function(nm) format_signal_value(selected[[nm]][i]), character(1L))
        paste(values, collapse = " | ")
      }, character(1L))
      signal_rows$LABEL <- labels
      valid <- is.finite(signal_rows$CHR) & signal_rows$CHR >= 1L & signal_rows$CHR <= 24L &
               is.finite(signal_rows$START) & is.finite(signal_rows$END)
      if (any(!valid)) {
        message("mplot signals: skipping ", sum(!valid),
                " matched row(s) with missing or invalid locus coordinates")
        signal_rows <- signal_rows[valid, , drop = FALSE]
      }
      if (nrow(signal_rows)) {
        swap <- signal_rows$START > signal_rows$END
        if (any(swap)) {
          z <- signal_rows$START[swap]
          signal_rows$START[swap] <- signal_rows$END[swap]
          signal_rows$END[swap] <- z
        }
      } else {
        signal_rows <- NULL
      }
    }
    message("mplot signals: trait=", trait, " match ", match_name, "=", target,
            " rows=", nrow(selected), " plottable=", if (is.null(signal_rows)) 0L else nrow(signal_rows))
  }

  read_bed <- function(path) {
    if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NULL)
    z <- tryCatch(utils::read.table(path, header = FALSE, comment.char = "#",
                                    stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(z) || ncol(z) < 4L) return(NULL)
    z <- z[, 1:4, drop = FALSE]
    names(z) <- c("CHR", "START", "END", "GENE")
    z$CHR <- normalize_chr(z$CHR)
    z$START <- suppressWarnings(as.numeric(z$START))
    z$END <- suppressWarnings(as.numeric(z$END))
    z$GENE <- as.character(z$GENE)
    z[is.finite(z$CHR) & is.finite(z$START) & is.finite(z$END) & nzchar(z$GENE), , drop = FALSE]
  }
  annotation_bed <- read_bed(mh_plot_bed)
  all_cis_bed <- read_bed(cis_bed)
  cis_rows <- if (!is.null(all_cis_bed)) {
    all_cis_bed[all_cis_bed$GENE == trait, , drop = FALSE]
  } else NULL
  cis_gene <- if (!is.null(cis_rows) && nrow(cis_rows)) trait else NULL

  # Every mplot run emits one atomic per-GWAS flag fragment.  The shell
  # coordinator merges these fragments into mplot/0flag.tsv under a lock.
  flag <- data.frame(GWAS = character(), FLAG = character(), stringsAsFactors = FALSE)
  if (!is.null(cis_rows) && nrow(cis_rows)) {
    gwas_chr <- normalize_chr(gwas$CHR)
    gwas_pos <- suppressWarnings(as.numeric(gwas$POS))
    significant <- is.finite(raw_p) & raw_p >= 0 & raw_p <= 5e-8
    cis_chr <- sort(unique(cis_rows$CHR))
    cis_significant <- FALSE
    for (chr_i in cis_chr) {
      r <- cis_rows[cis_rows$CHR == chr_i, , drop = FALSE]
      cis_center <- (min(r$START) + max(r$END)) / 2
      cis_significant <- cis_significant || any(significant & gwas_chr == chr_i &
                                                  gwas_pos >= cis_center - 5e5 &
                                                  gwas_pos <= cis_center + 5e5,
                                                na.rm = TRUE)
    }
    trans_chr <- sort(unique(gwas_chr[significant & !gwas_chr %in% cis_chr &
                                            is.finite(gwas_chr)]))
    if (!cis_significant && length(trans_chr)) {
      msg <- paste0("cis gene on ", paste(display_chr(cis_chr), collapse = ", "),
                    ", not Significant, but Significant on ",
                    paste(display_chr(trans_chr), collapse = ", "))
      flag <- data.frame(GWAS = trait, FLAG = msg, stringsAsFactors = FALSE)
      message("mplot flag: ", trait, ": ", msg)
    }
  }
  dir.create(dirname(flag_output), recursive = TRUE, showWarnings = FALSE)
  flag_tmp <- paste0(flag_output, ".tmp.", Sys.getpid())
  utils::write.table(flag, flag_tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (file.exists(flag_output)) unlink(flag_output, force = TRUE)
  if (!file.rename(flag_tmp, flag_output)) stop("Could not install mplot flag fragment: ", flag_output)

  mirror <- if (panel == "magma") read_magma(magma_file) else NULL
  if (!is.null(mirror)) {
    magma_grch <- attr(mirror, "grch")
    if (nzchar(magma_grch) && !identical(magma_grch, grch)) {
      stop("GWAS/MAGMA build mismatch: plotting GRCh", grch,
           " but ", magma_file, " was generated for GRCh", magma_grch,
           ". Use --grch auto or the matching build.")
    }
  }

  if (write_sig) {
    genomewide_p <- 5e-8
    magma_p <- 2.5e-6
    locus_size <- 1e6
    gwas_chr <- normalize_chr(gwas$CHR)
    gwas_pos <- suppressWarnings(as.numeric(gwas$POS))
    gwas_snp <- if ("SNP" %in% names(gwas)) as.character(gwas$SNP) else paste0(gwas_chr, ":", gwas_pos)
    significant <- is.finite(raw_p) & raw_p >= 0 & raw_p <= genomewide_p &
                   is.finite(gwas_chr) & gwas_chr >= 1L & gwas_chr <= 24L &
                   is.finite(gwas_pos) & gwas_pos > 0

    cis_start <- cis_end <- NA_real_
    cis_chr <- NA_integer_
    cis_mask <- rep(FALSE, nrow(gwas))
    if (!is.null(cis_rows) && nrow(cis_rows)) {
      cis_chr <- cis_rows$CHR[[1L]]
      cis_center <- (min(cis_rows$START) + max(cis_rows$END)) / 2
      cis_start <- cis_center - locus_size / 2
      cis_end <- cis_center + locus_size / 2
      cis_mask <- gwas_chr == cis_chr & gwas_pos >= cis_start & gwas_pos <= cis_end
      cis_mask[is.na(cis_mask)] <- FALSE
    }

    choose_lead <- function(ii) {
      ii[order(raw_p[ii], gwas_pos[ii], gwas_snp[ii], na.last = TRUE)[[1L]]]
    }
    top_idx <- integer()
    top_class <- character()
    cis_eligible <- which(significant & cis_mask)
    if (length(cis_eligible)) {
      top_idx <- c(top_idx, choose_lead(cis_eligible))
      top_class <- c(top_class, "cis")
    }

    # Match mh_plot(): significant background markers are split into loci at
    # >1 Mb gaps, then only the strongest locus per chromosome is colored green.
    trans_eligible <- which(significant & !cis_mask)
    trans_leads <- integer()
    if (length(trans_eligible)) {
      by_chr <- split(trans_eligible, gwas_chr[trans_eligible])
      all_loci <- list()
      k <- 0L
      for (ii in by_chr) {
        ii <- ii[order(gwas_pos[ii])]
        locus_id <- cumsum(c(TRUE, diff(gwas_pos[ii]) > locus_size))
        for (jj in split(ii, locus_id)) {
          lead <- choose_lead(jj)
          k <- k + 1L
          all_loci[[k]] <- data.frame(INDEX = lead, CHR = gwas_chr[[lead]],
                                      POS = gwas_pos[[lead]], P = raw_p[[lead]])
        }
      }
      loci <- do.call(rbind, all_loci)
      loci <- loci[order(loci$CHR, loci$P, loci$POS), , drop = FALSE]
      keep <- unlist(lapply(split(seq_len(nrow(loci)), loci$CHR), head, 1L), use.names = FALSE)
      trans_leads <- as.integer(loci$INDEX[keep])
      top_idx <- c(top_idx, trans_leads)
      top_class <- c(top_class, rep("trans", length(trans_leads)))
    }

    near_genes <- function(b, chr, pos, radius) {
      if (is.null(b) || !nrow(b)) return(NULL)
      x <- b[b$CHR == chr, , drop = FALSE]
      if (!nrow(x)) return(NULL)
      x$DIST <- pmax(x$START - pos, pos - x$END, 0)
      x <- x[x$DIST <= radius, , drop = FALSE]
      x[order(x$DIST, x$START, x$GENE), , drop = FALSE]
    }
    gene_label <- function(chr, pos) {
      x <- near_genes(annotation_bed, chr, pos, 0)
      if (is.null(x) || !nrow(x)) x <- near_genes(all_cis_bed, chr, pos, 0)
      if (!is.null(x) && nrow(x)) return(paste(head(unique(x$GENE), 3L), collapse = ","))
      x <- near_genes(annotation_bed, chr, pos, 1e6)
      if (is.null(x) || !nrow(x)) x <- near_genes(all_cis_bed, chr, pos, 1e6)
      if (!is.null(x) && nrow(x)) return(paste(head(unique(x$GENE), 3L), collapse = ","))
      b <- if (!is.null(annotation_bed) && any(annotation_bed$CHR == chr)) {
        annotation_bed[annotation_bed$CHR == chr, , drop = FALSE]
      } else if (!is.null(all_cis_bed)) {
        all_cis_bed[all_cis_bed$CHR == chr, , drop = FALSE]
      } else NULL
      if (is.null(b) || !nrow(b)) return(paste0("chr", chr, ":", format(pos, scientific = FALSE)))
      up <- b[b$END < pos, , drop = FALSE]
      down <- b[b$START > pos, , drop = FALSE]
      labels <- character()
      if (nrow(up)) {
        u <- up[which.max(up$END), ]
        labels <- c(labels, sprintf("%s (%+.1fMb)", u$GENE, -(pos - u$END) / 1e6))
      }
      if (nrow(down)) {
        d <- down[which.min(down$START), ]
        labels <- c(labels, sprintf("%s (%+.1fMb)", d$GENE, (d$START - pos) / 1e6))
      }
      paste(labels, collapse = ",")
    }

    cojo <- NULL
    if (nzchar(cojo_file) && file.exists(cojo_file) && file.info(cojo_file)$size > 0) {
      cojo <- tryCatch(read_gwas(cojo_file), error = function(e) {
        message("mplot sig: could not read COJO file; using GWAS values: ", conditionMessage(e))
        NULL
      })
    }
    find_column <- function(x, candidates) {
      if (is.null(x)) return(NA_character_)
      available <- toupper(sub("^#", "", names(x)))
      hit <- match(toupper(candidates), available, nomatch = 0L)
      hit <- hit[hit > 0L]
      if (length(hit)) names(x)[hit[[1L]]] else NA_character_
    }
    cojo_snp <- find_column(cojo, "SNP")
    cojo_chr <- find_column(cojo, "CHR")
    cojo_pos <- find_column(cojo, c("BP", "POS"))
    cojo_beta <- find_column(cojo, c("B", "BETA"))
    cojo_se <- find_column(cojo, "SE")
    cojo_p <- find_column(cojo, "P")
    gwas_beta <- find_column(gwas, "BETA")
    gwas_se <- find_column(gwas, "SE")
    value_at <- function(x, column, i) {
      if (is.null(x) || is.na(column) || !length(i)) return("NA")
      format_report_value(x[[column]][i[[1L]]])
    }

    common <- data.frame(SNP = character(), CHR = character(), POS = character(),
                         BETA = character(), SE = character(), P = character(),
                         GENE = character(), CLASS = character(), stringsAsFactors = FALSE)
    if (length(top_idx)) {
      ord <- order(gwas_chr[top_idx], gwas_pos[top_idx])
      top_idx <- top_idx[ord]
      top_class <- top_class[ord]
      common <- do.call(rbind, lapply(seq_along(top_idx), function(j) {
        i <- top_idx[[j]]
        ci <- integer()
        if (!is.null(cojo) && !is.na(cojo_snp)) {
          ci <- which(as.character(cojo[[cojo_snp]]) == gwas_snp[[i]])
          if (length(ci) > 1L && !is.na(cojo_chr) && !is.na(cojo_pos)) {
            exact_match <- normalize_chr(cojo[[cojo_chr]][ci]) == gwas_chr[[i]] &
                           suppressWarnings(as.numeric(cojo[[cojo_pos]][ci])) == gwas_pos[[i]]
            exact <- ci[which(exact_match)]
            if (length(exact)) ci <- exact
          }
        }
        gene <- if (top_class[[j]] == "cis") trait else gene_label(gwas_chr[[i]], gwas_pos[[i]])
        data.frame(
          SNP = gwas_snp[[i]], CHR = as.character(gwas_chr[[i]]),
          POS = format(gwas_pos[[i]], scientific = FALSE, trim = TRUE),
          BETA = if (length(ci)) value_at(cojo, cojo_beta, ci) else value_at(gwas, gwas_beta, i),
          SE = if (length(ci)) value_at(cojo, cojo_se, ci) else value_at(gwas, gwas_se, i),
          P = if (length(ci)) value_at(cojo, cojo_p, ci) else format_report_value(raw_p[[i]]),
          GENE = gene, CLASS = top_class[[j]], stringsAsFactors = FALSE
        )
      }))
    }

    magma_sig <- data.frame(GENE = character(), CHR = character(), REGION = character(),
                            ZSTAT = character(), P = character(), TRAIT = character(),
                            CLASS = character(), stringsAsFactors = FALSE)
    if (!is.null(mirror)) {
      mm <- mirror[is.finite(mirror$P) & mirror$P <= magma_p, , drop = FALSE]
      if (nrow(mm)) {
        mm <- mm[order(mm$CHR, mm$START, mm$P), , drop = FALSE]
        mm_class <- rep("trans", nrow(mm))
        if (is.finite(cis_chr)) {
          overlap <- mm$CHR == cis_chr & mm$STOP >= cis_start & mm$START <= cis_end
          mm_class[overlap] <- "cis"
        }
        magma_sig <- data.frame(
          GENE = as.character(mm$GENE), CHR = as.character(mm$CHR),
          REGION = paste0(format(mm$START, scientific = FALSE, trim = TRUE), "-",
                          format(mm$STOP, scientific = FALSE, trim = TRUE)),
          ZSTAT = vapply(mm$ZSTAT, format_report_value, character(1L)),
          P = vapply(mm$P, format_report_value, character(1L)),
          TRAIT = trait, CLASS = mm_class, stringsAsFactors = FALSE
        )
      }
    }

    dir.create(dirname(sig_tmp), recursive = TRUE, showWarnings = FALSE)
    con <- file(sig_tmp, open = "wt")
    tryCatch({
      writeLines("common:", con)
      utils::write.table(common, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
      writeLines(c("", "magma:"), con)
      utils::write.table(magma_sig, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }, finally = close(con))
    message("mplot sig: ", sig_output, " [common=", nrow(common), ", magma=", nrow(magma_sig),
            ", cojo=", if (is.null(cojo)) "fallback-to-GWAS" else cojo_file, "]")
  }

  # The two data sets retain their own -log10(P) coordinates, but their
  # display extents are fixed at 60% GWAS and 40% MAGMA instead of allowing an
  # extreme SNP P value to squeeze the lower panel.  Expand both spans together
  # when necessary so neither data set nor its labels are clipped.
  gwas_y <- -log10(plot_p[is.finite(plot_p) & plot_p > 0 & plot_p <= 1])
  gwas_threshold_y <- -log10(5e-8)
  top_need <- max(ceiling(max(gwas_y, gwas_threshold_y, na.rm = TRUE) + 3),
                  ceiling(gwas_threshold_y + 16))
  signal_label_room <- if (is.null(signal_rows)) 0 else 14
  if (is.null(mirror)) {
    plot_ylim <- c(-signal_label_room, top_need)
  } else {
    magma_y <- -log10(mirror$P[is.finite(mirror$P) & mirror$P > 0 & mirror$P <= 1])
    magma_threshold_y <- -log10(2.5e-6)
    bottom_need <- ceiling(max(magma_y, magma_threshold_y, na.rm = TRUE)) +
                   if (is.null(signal_rows)) 4 else signal_label_room
    top_span <- max(top_need, bottom_need * 0.60 / 0.40)
    bottom_span <- top_span * 0.40 / 0.60
    plot_ylim <- c(-bottom_span, top_span)
  }
  open_png(path, 4000L, if (is.null(mirror)) 1900L else 2600L)
  on.exit(grDevices::dev.off(), add = TRUE)
  signal_panel <- NULL
  if (!is.null(signal_rows)) {
    signal_panel <- function(x, y, getgenpos, ...) {
      midpoint_mb <- (signal_rows$START + signal_rows$END) / 2e6
      signal_x <- as.numeric(getgenpos(signal_rows$CHR, midpoint_mb))
      valid_x <- is.finite(signal_x)
      if (!any(valid_x)) return(invisible(NULL))
      signal_x <- signal_x[valid_x]
      signal_labels <- signal_rows$LABEL[valid_x]
      limits <- lattice::current.panel.limits()
      y_span <- diff(limits$ylim)
      label_y <- limits$ylim[[1L]] + y_span * (0.015 + 0.018 * ((seq_along(signal_x) - 1L) %% 4L))
      lattice::panel.abline(v = signal_x, lty = 2L, lwd = 1.25, col = "#2563EB")
      x_fraction <- (signal_x - limits$xlim[[1L]]) / diff(limits$xlim)
      label_adj <- ifelse(x_fraction < 0.15, 0, ifelse(x_fraction > 0.85, 1, 0.5))
      for (i in seq_along(signal_x)) {
        lattice::panel.text(signal_x[[i]], label_y[[i]], labels = signal_labels[[i]],
                            adj = c(label_adj[[i]], 0), col = "#2563EB", cex = 0.58, font = 2L)
      }
    }
  }

  print(mh_plot(plot_input, col = c("#4B5563", "#A7B0BD"), cis_gene = cis_gene,
                cis_bed = cis_bed, mh_plot_bed = mh_plot_bed,
                cis_color = "#D62728", other_top_color = "#159947",
                other_top_max_per_chr = 1, locus_size = 1e6, main = trait,
                mirror = mirror, mirror.sig.level = 2.5e-6,
                mirror.label.n = 10L, ylim = plot_ylim,
                threshold.col = "#E67E22", mirror.threshold.col = "#E67E22",
                chr.label.col = "#E67E22", panel.extra = signal_panel,
                xlab = if (is.null(mirror)) "Chromosome" else NULL))
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
  meta_file <- file.path(dirname(path), "magma.meta.tsv")
  magma_grch <- ""
  if (file.exists(meta_file)) {
    meta <- tryCatch(utils::read.table(meta_file, header = TRUE, sep = "\t", quote = "",
                                      comment.char = "", stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (!is.null(meta) && all(c("key", "value") %in% names(meta))) {
      z <- as.character(meta$value[meta$key == "grch"])
      if (length(z)) magma_grch <- z[[1L]]
    }
    gene_loc <- if (!is.null(meta) && all(c("key", "value") %in% names(meta))) {
      z <- as.character(meta$value[meta$key == "gene_loc"])
      if (length(z)) z[[1L]] else ""
    } else ""
    if (nzchar(gene_loc) && !file.exists(gene_loc)) {
      candidates <- unique(c(
        sub("NCBI(37|38)\\.gene\\.loc$", "NCBI.\\1.gene.loc", gene_loc),
        sub("NCBI\\.(37|38)\\.gene\\.loc$", "NCBI\\1.gene.loc", gene_loc)
      ))
      existing <- candidates[file.exists(candidates)]
      if (length(existing)) gene_loc <- existing[[1L]]
    }
    if (nzchar(gene_loc) && file.exists(gene_loc)) {
      loc <- tryCatch(utils::read.table(gene_loc, header = FALSE, stringsAsFactors = FALSE,
                                       comment.char = "", quote = ""),
                      error = function(e) NULL)
      if (!is.null(loc) && ncol(loc) >= 6L) {
        symbol_map <- setNames(as.character(loc[[6L]]), as.character(loc[[1L]]))
        symbol <- unname(symbol_map[gene])
        use_symbol <- !is.na(symbol) & nzchar(symbol)
        gene[use_symbol] <- symbol[use_symbol]
      }
    }
  }
  chr <- toupper(sub("^CHR", "", as.character(x$CHR)))
  chr[chr == "X"] <- "23"; chr[chr == "Y"] <- "24"; chr[chr %in% c("M", "MT")] <- "25"
  out <- data.frame(GENE = gene, CHR = suppressWarnings(as.integer(chr)),
                    START = suppressWarnings(as.numeric(x[[start_col]])),
                    STOP = suppressWarnings(as.numeric(x[[stop_col]])),
                    ZSTAT = if ("ZSTAT" %in% names(x)) suppressWarnings(as.numeric(x$ZSTAT)) else NA_real_,
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
  attr(out, "grch") <- magma_grch
  out
}

render_self(top_png)
if (!file.exists(top_png) || file.info(top_png)$size <= 0) stop("Top Manhattan panel was not created")
if (write_sig && (!file.exists(sig_tmp) || file.info(sig_tmp)$size <= 0)) {
  stop("Significant-hit summary was not created")
}

if (file.exists(output)) unlink(output)
if (!file.rename(top_png, output)) stop("Could not install output: ", output)
if (write_sig) {
  dir.create(dirname(sig_output), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(sig_output)) unlink(sig_output)
  if (!file.rename(sig_tmp, sig_output)) stop("Could not install significant-hit summary: ", sig_output)
}

unlink(top_png, force = TRUE)
message("mplot: ", output, " [panel=", panel, ", res=300",
        ", pixels=4000x", if (panel == "magma") 2600L else 1900L, "]")
}

main()
