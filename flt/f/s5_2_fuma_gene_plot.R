source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
# PLACO / FUMA gene mapping plots following Fig. 6a–b of Jiang et al. (2025).
# Add --fig6a-only to save only the circular plot.
# install.packages(c("data.table", "openxlsx"))

required <- c("data.table", "openxlsx")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(data.table))

# ---- Settings --------------------------------------------------------------
root <- Sys.getenv("FLT_OUT")
fuma_dir <- Sys.getenv("FLT_FUMA_DIR", file.path(root, "03placo/FUMA_job787285"))
smr_file <- file.path(root, "04smr", "matched.tsv")
out_dir <- file.path(Sys.getenv("FLT_PAPER"), "03fuma/gene_mapping")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- "G732-mdd_2025"
food_label <- "Chicken liking"
disease_label <- "MDD"

min_evidence <- 4L     # Fig. 6b: use 2L for the paper's display rule
genes_per_page <- 28L  # Fig. 6b: paginate long gene lists
magma_fdr <- 0.05
dpi <- 300
stopifnot(min_evidence %in% 1:6, genes_per_page >= 1)

# ---- Input / validation -----------------------------------------------------
paths <- c(genes = file.path(fuma_dir, "genes.txt"),
           loci = file.path(fuma_dir, "GenomicRiskLoci.txt"),
           magma = file.path(fuma_dir, "magma.genes.out"),
           smr = smr_file)
if (any(!file.exists(paths[names(paths) != "smr"]))) {stop("Missing files: ", paste(paths[names(paths) != "smr" & !file.exists(paths)], collapse = ", "))}
check_columns <- function(x, expected, path) {
  absent <- setdiff(expected, names(x))
  if (length(absent)) stop(basename(path), " missing columns: ", paste(absent, collapse = ", "))
}

genes <- fread(paths[["genes"]])
loci <- fread(paths[["loci"]])
magma <- fread(paths[["magma"]])
# Matched table is published after both traits pass SMR/HEIDI in the same tissue.
smr_used <- Sys.getenv("FLT_SMR_CURRENT") == "TRUE" && file.exists(paths[["smr"]])
smr <- if (smr_used) fread(paths[["smr"]]) else
  data.table(Tissue = character(), Gene = character(), Gene_id = character())
tissues_file <- file.path(root, "04smr/tissues.txt")
tissues <- if (smr_used && file.exists(tissues_file)) readLines(tissues_file) else character()
check_columns(genes, c("ensg", "symbol", "GenomicLocus", "posMapSNPs", "eqtlMapSNPs", "eqtlMapts"), paths[["genes"]])
check_columns(loci, c("GenomicLocus", "rsID", "chr", "start", "end", "p"), paths[["loci"]])
check_columns(magma, c("GENE", "SYMBOL", "P"), paths[["magma"]])
check_columns(smr, c("Tissue", "Gene", "Gene_id"), paths[["smr"]])

ensembl_id <- function(x) sub("\\.[0-9]+$", "", trimws(as.character(x)))
genes[, ensg := ensembl_id(ensg)]
magma[, GENE := ensembl_id(GENE)]
smr[, Gene_id := ensembl_id(Gene_id)]
loci[, GenomicLocus := as.character(GenomicLocus)]
if (anyDuplicated(loci$GenomicLocus)) stop("Duplicate GenomicLocus in GenomicRiskLoci.txt.")

# ---- Fig. 6a: disease -> food liking -> PLACO locus -> mapped gene ----------
mapped <- genes[!is.na(ensg) & nzchar(ensg) & (posMapSNPs > 0 | eqtlMapSNPs > 0)]
if (!nrow(mapped)) stop("No genes with positional or eQTL mapping.")
mapped[, label := fifelse(is.na(symbol) | !nzchar(symbol), ensg, symbol)]
links <- mapped[, .(GenomicLocus = trimws(strsplit(region_ids, ":", fixed = TRUE)[[1L]])), by = .(Ensembl = ensg, Gene = label, region_ids = as.character(GenomicLocus))]
links[, region_ids := NULL]
links <- unique(links, by = c("Ensembl", "GenomicLocus"))
links <- merge(links, loci[, .(GenomicLocus, chr, start, end, TopLeadSNP = rsID, P_PLACO = p)], by = "GenomicLocus", all.x = TRUE)
if (anyNA(links$TopLeadSNP)) stop("Some mapped genes have no matching FUMA locus.")
setorder(links, P_PLACO, GenomicLocus, Gene)
chosen <- unique(links$GenomicLocus)
circle <- copy(links)

# ---- Fig. 6b: six evidence criteria, joined by Ensembl gene ID --------------
magma <- magma[is.finite(P) & P >= 0 & P <= 1]
magma[, MAGMA_q := p.adjust(P, method = "BH")]
# Use the user's preselected matched results directly, without P-value filtering.
smr_pass <- smr[!is.na(Gene_id) & nzchar(Gene_id)]

gene_sets <- list(
  Position = genes[posMapSNPs > 0, ensg],
  MAGMA = magma[MAGMA_q < magma_fdr, GENE],
  eQTLGen = genes[grepl("eQTLGen_cis_eQTLs", eqtlMapts, fixed = TRUE), ensg],
  GTEx_Whole_Blood = genes[grepl("GTEx_v8/Whole_Blood", eqtlMapts, fixed = TRUE), ensg],
  SMR_Whole_Blood = smr_pass[Tissue == "Whole_Blood", Gene_id],
  SMR_eQTLGen = smr_pass[Tissue == "eQTLGen", Gene_id]
)
ids <- unique(unlist(gene_sets, use.names = FALSE))
evidence <- data.table(Ensembl = ids[!is.na(ids) & nzchar(ids)])
symbols <- unique(rbindlist(list(genes[, .(Ensembl = ensg, Gene = symbol)],
                                 magma[, .(Ensembl = GENE, Gene = SYMBOL)],
                                 smr[, .(Ensembl = Gene_id, Gene)]))[!is.na(Gene) & nzchar(Gene)],
                  by = "Ensembl")
evidence[symbols, on = "Ensembl", Gene := i.Gene]
evidence[is.na(Gene), Gene := Ensembl]
for (method in names(gene_sets)) evidence[, (method) := Ensembl %chin% gene_sets[[method]]]
if (!"Whole_Blood" %in% tissues) evidence[, SMR_Whole_Blood := NA]
if (!"eQTLGen" %in% tissues) evidence[, SMR_eQTLGen := NA]
evidence[, Evidence_count := rowSums(.SD, na.rm = TRUE), .SDcols = names(gene_sets)]
magma_per_gene <- unique(magma[order(P), .(Ensembl = GENE, MAGMA_P = P, MAGMA_q)],
                         by = "Ensembl")
evidence[magma_per_gene, on = "Ensembl", `:=`(MAGMA_P = i.MAGMA_P, MAGMA_q = i.MAGMA_q)]
setorder(evidence, -Evidence_count, MAGMA_P, Gene, na.last = TRUE)
shown <- evidence[Evidence_count >= min_evidence]


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output <- function(suffix) file.path(out_dir, paste0(prefix, "_", suffix))
fwrite(links, output("Fig6a_locus_genes.tsv"), sep = "\t", na = "NA")
fwrite(circle, output("Fig6a_plotted_locus_genes.tsv"), sep = "\t", na = "NA")
fwrite(evidence, output("Fig6b_gene_evidence.tsv"), sep = "\t", na = "NA")
print(evidence[, lapply(.SD, sum, na.rm = TRUE), .SDcols = names(gene_sets)])
message(sprintf("Mapped genes: %d; locus–gene links: %d; genes with >=2 criteria: %d; plotted genes: %d",
                uniqueN(links$Ensembl), nrow(links), sum(evidence$Evidence_count >= 2), nrow(shown)))

# ---- Circular layout (base graphics, as in s8_1_gen_mapping_plot.R) ---------
xy <- function(r, a) c(r * cos(a), r * sin(a))
branch <- function(p, q, col = "#B5B5B5") {
  t <- seq(0, 1, length.out = 50)
  mid <- (sqrt(sum(p*p)) + sqrt(sum(q*q))) / 2
  ctrl <- function(v) if (sum(v*v) < 1e-12) v else v / sqrt(sum(v*v)) * mid
  c1 <- ctrl(p); c2 <- ctrl(q)
  b <- outer((1-t)^3, p) + outer(3*(1-t)^2*t, c1) +
    outer(3*(1-t)*t^2, c2) + outer(t^3, q)
  lines(b, col = col, lwd = 1.1)
}
radial_text <- function(r, a, label, cex, font = 1) {
  deg <- (a * 180/pi) %% 360
  left <- deg > 90 & deg < 270
  p <- xy(r, a)
  text(p[1], p[2], label, srt = if (left) deg + 180 else deg,
       adj = c(if (left) 1 else 0, .5), cex = cex, font = font, xpd = NA)
}
draw_circle <- function() {
  par(mar = c(3, 1, 3, 1), family = "sans")
  color <- "#934953"
  x <- copy(circle)
  x[, locus_order := match(GenomicLocus, chosen)]
  setorder(x, locus_order, Gene, Ensembl)
  # Keep each gene-locus assignment as a separate leaf, as in the reference.
  n <- nrow(x)
  gap <- c(0, if (n > 1L) ifelse(x$GenomicLocus[-1L] !=
                                x$GenomicLocus[-n], 1.2, 0) else numeric())
  pos <- seq_len(n) + cumsum(gap)
  angles <- pi / 2 - 2 * pi * (pos - .5) / (max(pos) + 1.5)
  plot.new()
  plot.window(c(-1.44, 1.44), c(-1.44, 1.44), asp = 1)

  food_angle <- mean(range(angles))
  food_point <- xy(.29, food_angle)
  branch(c(0, 0), food_point)
  nodes <- list(c(0, 0), food_point)
  locus_labels <- list()
  for (loc in unique(x$GenomicLocus)) {
    ii <- which(x$GenomicLocus == loc)
    locus_angle <- mean(range(angles[ii]))
    locus_point <- xy(.63, locus_angle)
    branch(food_point, locus_point)
    nodes[[length(nodes) + 1L]] <- locus_point
    locus_labels[[length(locus_labels) + 1L]] <- list(
      a = locus_angle, label = x$TopLeadSNP[ii[1L]])
    for (j in ii) {
      gene_point <- xy(1, angles[j])
      branch(locus_point, gene_point)
      nodes[[length(nodes) + 1L]] <- gene_point
    }
  }
  points(do.call(rbind, nodes), pch = 21, bg = "white", col = color,
         lwd = 1.5, cex = .75)
  gene_cex <- if (n > 150) .55 else .72
  for (j in seq_len(n)) radial_text(1.025, angles[j], x$Gene[j], gene_cex, 2)
  for (v in locus_labels) radial_text(.65, v$a, v$label, .58, 2)
  radial_text(.32, food_angle, food_label, .95)
  text(0, .08, disease_label, col = color, cex = 2.5, font = 2)
  title(main = sprintf("a   All mapped genes across %d PLACO loci", length(chosen)),
        col.main = "#192637", line = 1)
  mtext("Disease / food liking / genomic locus (lead SNP) / mapped gene",
        side = 1, line = .5, cex = .7)

}

# ---- Evidence dot plot -----------------------------------------------------
method_labels <- c("Positional mapping", "MAGMA gene", "eQTLGen cis-eQTL",
                   "GTEx v8 whole blood", "SMR - GTEx whole blood", "SMR - eQTLGen")
method_colors <- c("#298273", "#456A9B", "#6D76B2", "#9465A2", "#B9783F", "#A84E5A")
dot_pages <- split(seq_len(nrow(shown)), ceiling(seq_len(nrow(shown)) / genes_per_page))

# Fig6b: full portrait-A4 paper width, compact 70 mm height; all text 8 pt.
dot_width <- 210/25.4
dot_height <- 70/25.4
dot_font_pt <- 8
draw_dots <- function(index, page) {
  x <- shown[index]
  par(mai=c(.82,1.48,.48,.12),family="sans",ps=dot_font_pt,
      cex=1,cex.axis=1,cex.lab=1,cex.main=1,cex.sub=1,mgp=c(2,.35,0))
  plot.new()
  plot.window(c(.45,nrow(x)+.55),c(.45,6.55),xaxs="i",yaxs="i")
  for(i in seq_along(gene_sets)) {
    y <- length(gene_sets)+1L-i
    if(i%%2L==1L) rect(.45,y-.49,nrow(x)+.55,y+.49,col="#F1F4F8",border=NA)
    flag <- x[[names(gene_sets)[i]]]
    points(seq_len(nrow(x)),rep(y,nrow(x)),pch=16,
           cex=ifelse(flag,1.1,.38),col=ifelse(flag,method_colors[i],"#BEC9D6"))
  }
  axis(2,at=6:1,labels=method_labels,las=1,tick=FALSE,cex.axis=1)
  axis(1,at=seq_len(nrow(x)),labels=x$Gene,las=2,tick=FALSE,cex.axis=1)
  page_label <- if(length(dot_pages)>1L) paste0(" | ",page,"/",length(dot_pages)) else ""
  title(main=paste0("b   Genes meeting at least ",min_evidence," of six criteria",page_label),
        line=2,col.main="#192637",cex.main=1)
  mtext("Large dot: meets criterion    Small grey dot: does not meet criterion",
        side=3,line=.65,cex=1,col="#526073")
  note <- paste0("MAGMA BH-FDR < ",magma_fdr,
                 "; SMR uses the nominal paired screen. Missing tissues have unavailable evidence.")
  mtext(note,side=1,line=4.8,cex=1,col="#526073",adj=1)
}

# ---- Save PDF and PNG ------------------------------------------------------
save_plot <- function(filename, format, width, height, draw, pointsize=12) {
  if (format == "pdf") {
    pdf(filename, width = width, height = height, pointsize=pointsize, useDingbats = FALSE)
  } else {
    png(filename, width = width, height = height, units = "in", res = dpi, bg = "white", pointsize=pointsize)
  }
  on.exit(dev.off())
  draw()
}

# Increase the canvas for the full gene set while preserving readable labels.
circle_size <- max(14, ceiling(nrow(circle) / 20))
save_plot(output("Fig6a_reference_style.pdf"), "pdf", circle_size, circle_size, draw_circle)
# save_plot(output("Fig6a_reference_style.png"), "png", circle_size, circle_size, draw_circle)
if (nrow(shown) && !("--fig6a-only" %in% commandArgs(trailingOnly = TRUE))) {
for (page in seq_along(dot_pages)) {
  suffix <- if (length(dot_pages) == 1L) "Fig6b.png" else sprintf("Fig6b_%02d.png", page)
  save_plot(output(suffix), "png", dot_width, dot_height, function() draw_dots(dot_pages[[page]], page), pointsize=dot_font_pt)
}
}
message("Outputs saved in: ", out_dir)

