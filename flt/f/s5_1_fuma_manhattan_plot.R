source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
# PLACO Manhattan plot with FUMA lead / independent significant SNPs.
# Run: Rscript scripts/s3_2_placo_fuma_plot.R
# Packages: install.packages(c("data.table", "R.utils", "ggplot2", "ggrepel"))
# Coordinates (GRCh38) and P_PLACO come exclusively from local PLACO_all.
# FUMA supplies only rsID highlight lists; no liftOver or FUMA rerun is needed.

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
  library(ggplot2)
  library(ggrepel)
})

# ---- Settings --------------------------------------------------------------
root <- Sys.getenv("FLT_OUT")
input_dir <- dirname(flt_placo_file(""))
prefix <- "G732-mdd_2025_PLACO"
fuma_dir <- Sys.getenv("FLT_FUMA_DIR")
out_dir <- Sys.getenv("FLT_PAPER")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_title <- "G732 × MDD (PLACO)"
label_mode <- "lead" # "lead", "independent", or "none"
label_top_n <- Inf   # e.g. 40 for a less crowded plot; ranked by PLACO P
threshold <- 1e-5    # Matches params.config: leadP = 1e-5
show_genome_wide_line <- TRUE # Also mark conventional P = 5e-8
width <- 18
height <- 8
dpi <- 300
y_axis_max <- 12
# Thin only background points in display-sized bins. All P <= threshold and
# all FUMA highlighted SNPs are preserved. Set FALSE to plot every SNP.
thin_background <- TRUE
chr_colors <- rep(c("#4594BD", "#E5B63C", "#51468D", "#626E33",
                    "#BD287A", "#D86B32", "#C581B5", "#79BCA1",
                    "#782B44", "#307B68"), length.out = 22)
stopifnot(label_mode %in% c("lead", "independent", "none"))

# ---- Input / validation ----------------------------------------------------
paths <- c(full = flt_placo_file("_all.tsv.gz"),
           lead = file.path(fuma_dir, "leadSNPs.txt"),
           ind = file.path(fuma_dir, "IndSigSNPs.txt"))
if (!file.exists(paths[["full"]])) stop("Missing PLACO summary statistics: ", paths[["full"]])
annotate_fuma <- Sys.getenv("FLT_FUMA_MATCHED", "FALSE") == "TRUE"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Reading local PLACO GRCh38 coordinates and P_PLACO...")
d <- fread(paths[["full"]], select = c("SNP", "CHR1", "POS1", "P_PLACO"))
setnames(d, c("CHR1", "POS1", "P_PLACO"), c("CHR", "BP", "P"))
if (anyDuplicated(d$SNP)) stop("Duplicate SNP IDs in PLACO_all; resolve before matching.")
n_input <- nrow(d)
d[, CHR := suppressWarnings(as.integer(sub("^chr", "", as.character(CHR))))]
bad <- d[is.na(CHR) | !CHR %in% 1:22 | !is.finite(BP) | BP <= 0 | !is.finite(P) | P < 0 | P > 1]
fwrite(bad, file.path(out_dir, paste0(prefix, "_excluded_plot_rows.tsv")), sep = "\t")
d <- d[!is.na(CHR) & CHR %in% 1:22 & is.finite(BP) & BP > 0 & is.finite(P) & P >= 0 & P <= 1]
if (!nrow(d)) stop("No valid autosomal SNPs to plot.")
if (any(d$P == 0)) warning("Zero P values capped at .Machine$double.xmin for plotting.")
d[, logp := -log10(pmax(P, .Machine$double.xmin))]
lead <- if (annotate_fuma) fread(paths[["lead"]], select = "rsID") else data.table(rsID = character())
ind <- if (annotate_fuma) fread(paths[["ind"]], select = "rsID") else data.table(rsID = character())
y_axis_max <- max(y_axis_max, ceiling(max(d$logp)))
d[, type := fifelse(SNP %chin% lead$rsID, "Lead SNPs", fifelse(SNP %chin% ind$rsID, "IndSig SNPs", "Background"))]
unmatched <- rbindlist(list(lead[, .(SNP = rsID, source = "leadSNPs")],ind[, .(SNP = rsID, source = "IndSigSNPs")]))[!SNP %chin% d$SNP]
# fwrite(unmatched, file.path(out_dir, paste0(prefix, "_unmatched_FUMA_SNPs.tsv")), sep = "\t")
if (nrow(unmatched)) warning(nrow(unmatched), " FUMA entries could not be plotted; see unmatched table.")

# ---- Cumulative chromosome positions --------------------------------------
axis <- d[, .(length = max(BP)), by = CHR][order(CHR)]
gap <- sum(axis$length) * 0.003
axis[, offset := c(0, head(cumsum(length + gap), -1))]
axis[, center := offset + length / 2]
d[axis, on = "CHR", x := BP + i.offset]
hi <- d[type != "Background"]
labels <- switch(label_mode, lead = hi[type == "Lead SNPs"],
                 independent = hi[SNP %chin% ind$rsID], none = hi[0])
setorder(labels, P, SNP)
if (is.finite(label_top_n)) labels <- head(labels, label_top_n)
# fwrite(hi[, .(SNP, CHR, BP, P, logp, type)], file.path(out_dir, paste0(prefix, "_highlighted_SNPs.tsv")), sep = "\t")

bg <- d[type == "Background"]
if (thin_background) {
  # Keep strongest point in each x/y bin, preserving every significant point.
  xbin <- max(d$x) / (width * dpi)
  bg[, `:=`(xb = floor(x / xbin), yb = floor(logp / 0.025))]
  bg <- rbindlist(list(bg[P <= threshold], bg[P > threshold][order(P), .SD[1L], by = .(CHR, xb, yb)]), use.names = TRUE, fill = TRUE)
}
message(sprintf("Valid SNPs: %s; drawn background: %s; lead: %s; independent-only: %s; labels: %s", nrow(d), nrow(bg), sum(hi$type == "Lead SNPs"), sum(hi$type == "IndSig SNPs"), nrow(labels)))
label_floor <- -log10(threshold) + 1

# ---- Plot -----------------------------------------------------------------
p <- ggplot() +
  geom_point(data = bg, aes(x, logp, color = factor(CHR)), size = 0.22, alpha = 0.65) +
  scale_color_manual(values = setNames(chr_colors, as.character(1:22)), guide = "none") +
  geom_hline(yintercept = -log10(threshold), linetype = "dotted", linewidth = 0.6, color = "#333333") +
  geom_point(data = hi, aes(x, logp, fill = type), shape = 21, size = 2, stroke = 0.25, color = "white") +
  scale_fill_manual(values = c("Lead SNPs" = "#781B31", "IndSig SNPs" = "#397660"), breaks = c("Lead SNPs", "IndSig SNPs"), name = NULL) +
  scale_x_continuous(breaks = axis$center, labels = axis$CHR, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = seq(0, y_axis_max, by = 5),  expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(0, y_axis_max), clip = "on") +
  labs(title = plot_title, x = "Chromosome", y = expression(-log[10](italic(P))),
       caption = paste0("Dotted: PLACO candidate threshold P = 1e-5",
                        if (show_genome_wide_line) "; dashed: P = 5e-8" else "",
                        if (annotate_fuma) ". FUMA highlights matched by rsID." else ". Native GRCh38 coordinates.")) +
  theme_classic(base_size = 17) +
  theme(axis.text = element_text(size = 15, color = "#222222", face = "bold"),
        axis.title = element_text(size = 19, face = "bold"),
        plot.title = element_text(size = 21, face = "bold"),
        legend.position = "right", legend.text = element_text(size = 16),
        plot.caption = element_text(size = 11, hjust = 0),
        plot.margin = margin(15, 25, 10, 10))
if (show_genome_wide_line) p <- p + geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "#999999", linewidth = 0.4)
if (nrow(labels)) p <- p + ggrepel::geom_text_repel(
  data = labels, aes(x, logp, label = SNP), seed = 732,
  size = 4, fontface = "bold", color = "#444444", min.segment.length = 0,
  segment.color = "#444444", segment.linetype = "dashed", segment.size = 0.3,
  box.padding = 0.35, point.padding = 0.15, max.overlaps = Inf,
  max.iter = 100000, max.time = 20, force = 3,
  ylim = c(label_floor, y_axis_max - 0.5))

base <- file.path(out_dir, "Fig5A")
ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
writeLines(c(paste("Input rows:", n_input), paste("Valid rows:", nrow(d)),
             paste("Excluded rows:", nrow(bad)),
             "Coordinates: PLACO_all CHR1/POS1 (GRCh38); P values: P_PLACO.",
             "FUMA supplies only rsID highlight lists; no liftOver or FUMA rerun required.",
             capture.output(sessionInfo())), paste0(base, "_run_info.txt"))
message("Saved: ", base, ".png")
