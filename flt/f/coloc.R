# Single-signal Wakefield ABF colocalization in native GWAS coordinates.
# Reuse allele/coordinate QC but retain all variants, including the large Z
# values excluded by the composite-null PLACO test.
source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
options(flt.placo.functions_only = TRUE)
source(file.path(Sys.getenv("FLT_CODE"), "s3_2_placo.R"))
library(data.table)
root <- file.path(flt_root, "03placo")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
prefix <- flt_placo_file("")
if (!file.exists(paste0(prefix, "_hits.tsv"))) stop("Run PLACO for G732 and mdd_2025 first")
hits <- fread(paste0(prefix, "_hits.tsv"))
qc <- function(...) invisible(NULL)
width <- as.integer(Sys.getenv("FLT_COLOC_HALF_WINDOW", "250000"))
min_snps <- 100L
prior_sd <- c(as.numeric(Sys.getenv("FLT_COLOC_FLT_PRIOR_SD", "0.15")), 0.20)
stopifnot(is.finite(width), width > 0L, all(is.finite(prior_sd) & prior_sd > 0))
# Greedily choose the strongest candidate within each non-overlapping native
# window. These are physical windows and may be correlated in long-range LD.
setorder(hits, P_PLACO, SNP)
hits <- hits[is.finite(POS1) & !is.na(CHR1)]
loci <- integer()
for (i in seq_len(nrow(hits))) {
    if (!any(hits$CHR1[loci] == hits$CHR1[i] & abs(hits$POS1[loci] - hits$POS1[i]) <= 2 * width))
        loci <- c(loci, i)
}
if (!length(loci)) {
    fwrite(data.table(lead = character(), status = character()), file.path(root, "coloc-results.csv"))
    quit(status = 0)
}
plan <- hits[loci, .(lead = SNP, chr = as.character(CHR1), center = POS1, PLACO_P = P_PLACO)]
plan[, `:=`(start = pmax(1, center - width), end = center + width)]
read_window_inputs <- function(file, label) {
    d <- read_placo(file, label, qc)
    ranges <- as.data.table(plan)[, .(CHR = chr, begin = start, end)]
    d[, BP_END := POS]
    setkey(ranges, CHR, begin, end)
    d <- foverlaps(d[is.finite(POS)], ranges, by.x = c("CHR", "POS", "BP_END"), nomatch = 0L)
    d[, .(SNP, CHR, POS, EA, NEA, BETA, SE, P, Z)]
}
a <- read_window_inputs(files[1], trait1)
b <- read_window_inputs(files[2], trait2)
d <- harmonise_placo(a, b, qc)
rm(a, b)

lbf <- function(beta, se, sd) {
    r <- sd^2 / (sd^2 + se^2)
    .5 * (log1p(-r) + r * (beta / se)^2)
}
lsum <- function(a) {
    m <- max(a)
    m + log(sum(exp(a - m)))
}
coloc_fit <- function(z, p12) {
    x <- lbf(z$BETA1, z$SE1, prior_sd[1])
    y <- lbf(z$BETA2, z$SE2, prior_sd[2])
    xy <- x + y
    # Sum different-SNP configurations stably without subtracting large BFs.
    best <- which.max(y)
    y_other <- rep(lsum(y), length(y)) + log1p(-exp(y - lsum(y)))
    y_other[best] <- lsum(y[-best])
    p1 <- p2 <- 1e-4
    logp <- c(0, log(p1) + lsum(x), log(p2) + lsum(y),
              log(p1 * p2) + lsum(x + y_other), log(p12) + lsum(xy))
    list(pp = exp(logp - lsum(logp)), lead = z$SNP[which.max(xy)],
         pp_snp = exp(max(xy) - lsum(xy)))
}
result <- lapply(seq_len(nrow(plan)), function(i) {
    locus <- plan[i]
    z <- d[CHR1 == locus$chr & between(POS1, locus$start, locus$end)]
    rbindlist(lapply(c(1e-6, 5e-6, 1e-5), function(p12) {
        fit <- if (nrow(z) >= min_snps && locus$lead %in% z$SNP) coloc_fit(z, p12) else NULL
        cbind(locus, data.table(n_shared = nrow(z), p1 = 1e-4, p2 = 1e-4, p12,
          sd_flt = prior_sd[1], sd_mdd = prior_sd[2],
          status = if (is.null(fit)) "insufficient_coverage" else "ok",
          PP_H0 = fit$pp[1] %||% NA_real_, PP_H1 = fit$pp[2] %||% NA_real_,
          PP_H2 = fit$pp[3] %||% NA_real_, PP_H3 = fit$pp[4] %||% NA_real_,
          PP_H4 = fit$pp[5] %||% NA_real_, shared_lead = fit$lead %||% NA_character_,
          PP_SNP_given_H4 = fit$pp_snp %||% NA_real_))
    }))
})
fwrite(rbindlist(result), file.path(root, "coloc-results.csv"))
