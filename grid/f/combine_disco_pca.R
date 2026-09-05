#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value for ", name, call. = FALSE)
  args[i + 1L]
}

dir <- get_arg("--dir")
pca_file <- get_arg("--pca")
out <- get_arg("--out", pca_file)
med_file <- get_arg("--med")
outdir <- get_arg("--outdir")
n_pc <- as.integer(get_arg("--n-pc", "20"))
distance_pcs <- as.integer(get_arg("--distance-pcs", "10"))
if (is.null(outdir)) stop("Missing --outdir", call. = FALSE)
if (!is.finite(n_pc) || n_pc < 5L) stop("--n-pc must be >=5", call. = FALSE)
if (!is.finite(distance_pcs) || distance_pcs < 5L || distance_pcs > n_pc) stop("Invalid --distance-pcs", call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!is.null(dir)) {
  if (!dir.exists(dir)) stop("Missing --dir directory", call. = FALSE)
  if (is.null(out)) out <- file.path(dir, "ukb.discodivas.pca.tsv.gz")
  files <- file.path(dir, paste0("chr", 1:22, ".sscore"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing PCA score files: ", paste(missing, collapse = ", "), call. = FALSE)

  scores <- rbindlist(lapply(seq_along(files), function(chr) {
    x <- fread(files[chr], showProgress = FALSE)
    id_col <- intersect(c("IID", "#IID"), names(x))[1L]
    score_cols <- grep("_SUM$", names(x), value = TRUE)
    if (length(id_col) != 1L || is.na(id_col) || length(score_cols) != n_pc) {
      stop("Unexpected PLINK2 PCA columns in ", files[chr], "; expected ", n_pc,
           " score columns but found ", length(score_cols), call. = FALSE)
    }
    y <- x[, c(id_col, score_cols), with = FALSE]
    setnames(y, id_col, "IID")
    y[, IID := as.character(IID)]
    if (anyDuplicated(y$IID)) stop("Duplicate IID in ", files[chr], call. = FALSE)
    y[, chr := chr]
    y
  }), fill = TRUE)

  coverage <- scores[, .(n_chr = uniqueN(chr)), by = IID]
  incomplete <- coverage[n_chr != 22L]
  if (nrow(incomplete)) {
    warning(nrow(incomplete), " IDs are absent from one or more chromosome PCA files and will be excluded", call. = FALSE)
    scores <- scores[!IID %chin% incomplete$IID]
  }
  score_cols <- grep("_SUM$", names(scores), value = TRUE)
  pca <- scores[, lapply(.SD, sum), by = IID, .SDcols = score_cols]
  setnames(pca, score_cols, paste0("PC", seq_len(n_pc)))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  fwrite(pca, out, sep = "\t", quote = FALSE)
  pca_file <- out
} else {
  if (is.null(pca_file) || !file.exists(pca_file)) stop("Provide --dir or an existing --pca file", call. = FALSE)
  pca <- fread(pca_file, showProgress = FALSE)
}

if (toupper(get_arg("--projection-only", "FALSE")) == "TRUE") quit(status = 0L)
if (is.null(med_file) || !file.exists(med_file)) stop("Missing --med file", call. = FALSE)
med <- fread(med_file, showProgress = FALSE)
id_col <- intersect(c("IID", "#IID", "eid"), names(pca))[1L]
if (length(id_col) != 1L || is.na(id_col)) stop("PCA file has no IID column", call. = FALSE)
setnames(pca, id_col, "IID")
pca[, IID := as.character(IID)]
all_pc_cols <- paste0("PC", seq_len(n_pc))
if (!all(all_pc_cols %in% names(pca))) stop("PCA file lacks ", paste(setdiff(all_pc_cols, names(pca)), collapse = ","), call. = FALSE)
if (anyNA(pca[, ..all_pc_cols])) stop("Target PCA contains missing values", call. = FALSE)

pop_col <- names(med)[1L]
med_available <- intersect(all_pc_cols, names(med))
actual_distance_pcs <- min(distance_pcs, length(med_available))
if (actual_distance_pcs < 5L) stop("Reference-center file contains fewer than five PCs", call. = FALSE)
dist_pc_cols <- paste0("PC", seq_len(actual_distance_pcs))
if (!all(dist_pc_cols %in% names(med))) stop("Reference centers must contain consecutive PC1-PC", actual_distance_pcs, call. = FALSE)
if (anyNA(med[, ..dist_pc_cols])) stop("Reference centers contain missing values", call. = FALSE)

x <- as.matrix(pca[, ..dist_pc_cols]); storage.mode(x) <- "double"
centers <- as.matrix(med[, ..dist_pc_cols]); storage.mode(centers) <- "double"
center_names <- toupper(as.character(med[[pop_col]]))
rownames(centers) <- center_names

dist_mat <- vapply(seq_len(nrow(centers)), function(i) {
  sqrt(rowSums(sweep(x, 2L, centers[i, ], FUN = "-")^2))
}, numeric(nrow(pca)))
if (is.null(dim(dist_mat))) dist_mat <- matrix(dist_mat, ncol = nrow(centers))
colnames(dist_mat) <- center_names
nearest_i <- max.col(-dist_mat, ties.method = "first")
nearest_dist <- dist_mat[cbind(seq_len(nrow(dist_mat)), nearest_i)]
sorted_dist <- t(apply(dist_mat, 1L, sort, partial = 1:2))
distance_margin <- sorted_dist[, 2L] - sorted_dist[, 1L]
relative_margin <- distance_margin / pmax(sorted_dist[, 2L], sqrt(.Machine$double.eps))

out_distance <- data.table(
  IID = pca$IID,
  closest_reference = center_names[nearest_i],
  nearest_euclidean_distance = nearest_dist,
  second_minus_first_distance = distance_margin,
  relative_distance_margin = relative_margin
)
out_distance <- cbind(out_distance, as.data.table(dist_mat))
fwrite(out_distance, file.path(outdir, "ukb_reference_distances.tsv.gz"), sep = "\t", quote = FALSE)

pc_summary <- rbindlist(lapply(all_pc_cols, function(pc) {
  z <- pca[[pc]]
  q <- quantile(z, c(0, .01, .25, .5, .75, .99, 1), names = FALSE)
  data.table(PC = pc, n = length(z), mean = mean(z), sd = sd(z), min = q[1L], q01 = q[2L],
             q25 = q[3L], median = q[4L], q75 = q[5L], q99 = q[6L], max = q[7L])
}))
nearest_counts <- out_distance[, .(N = .N, percent = 100 * .N / nrow(out_distance)), by = closest_reference][order(-N)]
center_dist <- as.data.table(as.matrix(dist(centers, upper = TRUE, diag = TRUE)), keep.rownames = "from")
setnames(center_dist, c("from", center_names))

match_files <- if (!is.null(dir)) file.path(dir, paste0("chr", 1:22, ".sscore.vars")) else character()
match_qc <- if (length(match_files) && all(file.exists(match_files))) {
  rbindlist(lapply(1:22, function(chr) {
    data.table(chr = chr, matched_pca_variants = nrow(fread(match_files[chr], header = FALSE, showProgress = FALSE)))
  }))
} else data.table(chr = integer(), matched_pca_variants = integer())

run_info <- data.table(
  item = c("target_n", "projected_pcs", "distance_pcs", "reference_groups", "pca_file", "med_file"),
  value = c(nrow(pca), n_pc, actual_distance_pcs, paste(center_names, collapse = ","),
            normalizePath(pca_file), normalizePath(med_file))
)

wb <- createWorkbook()
header_style <- createStyle(fontColour = "white", fgFill = "#1F4E78", textDecoration = "bold", halign = "center")
write_sheet <- function(name, z) {
  name <- substr(name, 1L, 31L)
  addWorksheet(wb, name)
  if (!nrow(z)) { writeData(wb, name, "No rows"); return(invisible(NULL)) }
  writeDataTable(wb, name, z, tableStyle = "TableStyleMedium2")
  addStyle(wb, name, header_style, rows = 1, cols = seq_len(ncol(z)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(z)), widths = "auto")
}
write_sheet("Run_info", run_info)
write_sheet("Target_PC_summary", pc_summary)
write_sheet("Reference_centers", med[, c(pop_col, dist_pc_cols), with = FALSE])
write_sheet("Closest_counts", nearest_counts)
write_sheet("Center_distances", center_dist)
write_sheet("Variant_match", match_qc)
saveWorkbook(wb, file.path(outdir, "pca_qc.xlsx"), overwrite = TRUE)

set.seed(107)
idx <- if (nrow(pca) > 50000L) sample.int(nrow(pca), 50000L) else seq_len(nrow(pca))
plot_dt <- data.table(PC1 = pca$PC1[idx], PC2 = pca$PC2[idx])
center_dt <- data.table(reference = center_names, PC1 = centers[, "PC1"], PC2 = centers[, "PC2"])
theme_qc <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

p1 <- ggplot(plot_dt, aes(PC1, PC2)) + geom_point(alpha = .10, size = .30) +
  geom_point(data = center_dt, aes(PC1, PC2), shape = 21, fill = "white", size = 3, stroke = 1) +
  geom_text(data = center_dt, aes(PC1, PC2, label = reference), vjust = -1, fontface = "bold") +
  labs(title = "A. Target projection and 1KG reference centers",
       subtitle = "Random sample of up to 50,000 target participants", x = "PC1", y = "PC2") + theme_qc
p2 <- ggplot(out_distance, aes(nearest_euclidean_distance)) + geom_histogram(bins = 70) +
  labs(title = "B. Distance to closest reference center", x = paste0("Euclidean distance using PC1-PC", actual_distance_pcs), y = "Participants") + theme_qc
p3 <- ggplot(nearest_counts, aes(reorder(closest_reference, -N), N)) + geom_col() +
  geom_text(aes(label = paste0(format(N, big.mark = ","), " (", sprintf("%.1f", percent), "%)")), vjust = -.3, size = 3.2) +
  labs(title = "C. Closest reference center", subtitle = "QC fallback only; not a formal ancestry assignment", x = NULL, y = "Participants") +
  scale_y_continuous(expand = expansion(mult = c(0, .16))) + theme_qc
heat <- melt(center_dist, id.vars = "from", variable.name = "to", value.name = "distance")
p4 <- ggplot(heat, aes(from, to, fill = distance)) + geom_tile() + geom_text(aes(label = sprintf("%.1f", distance)), size = 3) +
  labs(title = "D. Distances among reference centers", x = NULL, y = NULL, fill = "Distance") + theme_qc
p5 <- if (nrow(match_qc)) {
  ggplot(match_qc, aes(chr, matched_pca_variants)) + geom_col() +
    labs(title = "E. PCA variants matched by chromosome", x = "Chromosome", y = "Matched variants") + theme_qc
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "Variant-match files unavailable") +
    labs(title = "E. PCA variants matched by chromosome") + theme_void()
}
p6 <- ggplot(pc_summary, aes(as.integer(sub("PC", "", PC)), sd)) + geom_line() + geom_point() +
  labs(title = "F. Scale of projected PCs", x = "PC", y = "Standard deviation") + theme_qc
fig <- (p1 | p2) / (p3 | p4) / (p5 | p6) + plot_annotation(title = "DiscoDivas PCA projection quality control")
ggsave(file.path(outdir, "Fig1.PCA_QC.png"), fig, width = 14, height = 13, dpi = 220, bg = "white")

cat("Target PCA samples:", nrow(pca), "\n")
cat("Projected PCs:", n_pc, "\n")
cat("Distance PCs:", actual_distance_pcs, "\n")
print(nearest_counts)
