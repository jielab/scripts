source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
# MAGMA enrichment plots for chicken liking x MDD. No additional enrichment runs.
# Rscript scripts/s5_3_enrichment_plot.R
# Optional environment variables: FLT_ROOT, FLT_FUMA_DIR, FLT_ENRICHMENT_OUT
required <- c("ggplot2", "patchwork", "deTS")
for (pkg in required) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
library(ggplot2)
root <- Sys.getenv("FLT_OUT")
input <- Sys.getenv("FLT_FUMA_DIR", file.path(root, "03placo/FUMA_job787285"))
out <- Sys.getenv("FLT_ENRICHMENT_OUT", file.path(root, "03fuma"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05
read_magma <- function(filename, type) {
  path <- file.path(input, filename)
  d <- read.table(path, header = TRUE, comment.char = "#", stringsAsFactors = FALSE,
                  check.names = FALSE, quote = "")
  stopifnot(all(c("VARIABLE", "TYPE", "NGENES", "BETA", "P") %in% names(d)))
  d <- d[d$TYPE == type, , drop = FALSE]
  valid <- is.finite(d$P) & d$P >= 0 & d$P <= 1
  if (any(!valid)) warning(filename, ": excluded ", sum(!valid), " invalid P values")
  d <- d[valid, , drop = FALSE]
  if (!nrow(d)) stop("No valid tests: ", filename)
  d$name <- if ("FULL_NAME" %in% names(d)) d$FULL_NAME else d$VARIABLE
  if (anyDuplicated(d$name)) stop("Duplicate tests: ", filename)
  d$P_BH <- p.adjust(d$P, "BH")
  d$P_Bonferroni <- p.adjust(d$P, "bonferroni")
  d$n_tests <- nrow(d)
  d$source_file <- filename
  d[order(d$P), ]
}
sets_file <- if (file.exists(file.path(input, "magma.gsa.out"))) "magma.gsa.out" else "magma.sets.out"
sets <- read_magma(sets_file, "SET")
tissue <- read_magma("magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out", "COVAR")
general <- read_magma("magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out", "COVAR")
if (nrow(tissue) != 54L) stop("Expected 54 GTEx specific tissues; found ", nrow(tissue))
pretty <- function(x) {
  x <- sub("^(GOCC_|GOBP_|GOMF_|GO_)", "", x)
  gsub("_", " ", x)
}
palette <- c("Gene set" = "#237C8C", Brain = "#3577B5", Pituitary = "#DF8844", Other = "#8A929B")
make_plot <- function(d, title, subtitle, correction = "P_Bonferroni", tissue_plot = FALSE,
                      significant_only = TRUE) {
  if (significant_only) d <- d[d[[correction]] < alpha, , drop = FALSE]
  if (!nrow(d)) return(ggplot() + annotate("text", 0, 0, label = "No significant results") + theme_void() + labs(title = title))
  d <- d[order(d[[correction]], d$P), ]
  d$label <- vapply(pretty(d$name), function(s) paste(strwrap(s, width = 47), collapse = "\n"), character(1))
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  d$group <- if (tissue_plot) ifelse(grepl("^Brain", d$name), "Brain", ifelse(d$name == "Pituitary", "Pituitary", "Other")) else "Gene set"
  d$score <- -log10(pmax(d[[correction]], .Machine$double.xmin))
  p <- ggplot(d, aes(score, label)) +
    geom_vline(xintercept = -log10(alpha), linetype = 2, colour = "#A4ABB2", linewidth = .45) +
    geom_segment(aes(x = 0, xend = score, yend = label), colour = "#DCE3E8", linewidth = .65)
  if (tissue_plot) p <- p + geom_point(aes(colour = group), size = 3.5) else
    p <- p + geom_point(aes(size = NGENES, colour = group)) + scale_size_area(max_size = 9, name = "Genes in set")
  p + scale_colour_manual(values = palette, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, .12))) +
    labs(title = title, subtitle = subtitle, x = if (correction == "P_BH") expression(-log[10](P[BH])) else expression(-log[10](P[Bonferroni])), y = NULL) +
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(colour = "#263441", size = 10), legend.position = "bottom",
          plot.margin = margin(12, 16, 12, 12)) +
    guides(colour = if (tissue_plot) guide_legend(override.aes = list(size = 3)) else "none")
}
save_plot <- function(p, name, width, height) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out, paste0(name, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
}
a <- make_plot(sets, "A  MAGMA gene sets", sprintf("Bonferroni P < 0.05 | %s tested sets", format(nrow(sets), big.mark = ",")))
b <- make_plot(tissue, "B  GTEx v8 tissues", "MAGMA gene-property analysis | 54 tissues", tissue_plot = TRUE)
save_plot(a, "MAGMA_pathways_Bonferroni", 9, 5)
save_plot(b, "MAGMA_GTEx54_Bonferroni", 9, 7)
all_tissues <- make_plot(tissue, "GTEx v8: all 54 tissues",
                        "Chicken liking x MDD | MAGMA gene-property analysis",
                        tissue_plot = TRUE, significant_only = FALSE) +
  labs(caption = "Dashed line: Bonferroni P = 0.05 (54 tests).\nAdjusted P = 1 is plotted at zero; tissues are ordered by raw P.")
save_plot(all_tissues, "MAGMA_GTEx54_all_tissues", 11, 17)
# All tested sets: raw P preserves the distribution that would otherwise collapse
# to zero after Bonferroni adjustment. Colour uses full-family adjusted P values.
ranked <- sets[order(sets$P), ]
ranked$rank <- seq_len(nrow(ranked))
ranked$score <- -log10(pmax(ranked$P, .Machine$double.xmin))
ranked$status <- factor(ifelse(ranked$P_Bonferroni < alpha, "Bonferroni < 0.05",
                        ifelse(ranked$P_BH < alpha, "BH < 0.05 only", "Not significant")),
                        levels = c("Bonferroni < 0.05", "BH < 0.05 only", "Not significant"))
top <- ranked[ranked$P_Bonferroni < alpha, ]
key <- paste(paste0(top$rank, ". ", pretty(top$name)), collapse = "\n")
all_sets <- ggplot(ranked, aes(rank, score)) +
  geom_hline(yintercept = -log10(alpha / nrow(ranked)), linetype = 2, colour = "#DF8844") +
  geom_point(aes(colour = status), size = 1.8, alpha = .7) +
  geom_text(data = top, aes(label = rank), vjust = -1, size = 3.5, show.legend = FALSE) +
  scale_x_log10(breaks = c(1, 5, 22, 100, 1000, nrow(ranked)),
                labels = c("1", "5", "22", "100", "1,000", format(nrow(ranked), big.mark = ","))) +
  scale_y_continuous(expand = expansion(mult = c(.02, .12))) +
  scale_colour_manual(values = c("Bonferroni < 0.05" = "#237C8C", "BH < 0.05 only" = "#DF8844", "Not significant" = "#A4ABB2"), drop = FALSE) +
  labs(title = "MAGMA: all tested gene sets", subtitle = "Chicken liking x MDD | 17,008 tests; ranked by raw P",
       x = "Gene-set rank (log scale)", y = expression(-log[10](P[raw])), colour = NULL,
       caption = paste("Every tested gene set is included. Dashed line: Bonferroni threshold (0.05 / number of tests).",
                       "Numbered points: Bonferroni-significant gene sets.", key, sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 15),
        plot.caption = element_text(hjust = 0, lineheight = 1.15), plot.margin = margin(15, 20, 15, 15))
save_plot(all_sets, "MAGMA_pathways_all", 11, 8)
combined <- patchwork::wrap_plots(a, b, ncol = 2, widths = c(1, 1.05)) +
  patchwork::plot_annotation(title = "Chicken liking x MDD", subtitle = "PLACO associations: gene-set and tissue-expression enrichment",
                            theme = theme(plot.title = element_text(face = "bold", size = 18)))
save_plot(combined, "MAGMA_enrichment_AB", 16, 7.5)
save_plot(make_plot(sets, "MAGMA gene sets", "BH FDR < 0.05", "P_BH"), "MAGMA_pathways_BH", 10,
          max(5, 2 + sum(sets$P_BH < alpha) * .36))
save_plot(make_plot(general, "GTEx v8 general tissues", paste(nrow(general), "tissues; separately corrected"), tissue_plot = TRUE),
          "MAGMA_GTEx_general_Bonferroni", 9, 5)
results <- list(pathways = sets, GTEx54 = tissue, GTEx_general = general)
summary <- do.call(rbind, lapply(names(results), function(n) {
  d <- results[[n]]
  write.csv(d, file.path(out, paste0("MAGMA_", n, "_all.csv")), row.names = FALSE)
  data.frame(analysis = n, tests = nrow(d), Bonferroni_significant = sum(d$P_Bonferroni < alpha), BH_significant = sum(d$P_BH < alpha))
}))
write.csv(summary, file.path(out, "MAGMA_test_summary.csv"), row.names = FALSE)
writeLines(c("Chicken liking x MDD: existing FUMA MAGMA results.",
 "A: competitive gene-set analysis. Dot area represents NGENES; axis is -log10(Bonferroni P).",
 "B: MAGMA gene-property analysis using GTEx v8. Brain and pituitary are coloured separately.",
 "Bonferroni and BH corrections use all valid tests within each complete input table before selection.",
 "Significant results are displayed except in the all-pathways and all-54-tissues figures; CSV tables contain all tests.",
 "All-pathways plot: raw -log10(P) versus rank (log scale); dashed line is the raw Bonferroni threshold 0.05/n_tests. Colours use adjusted P.",
 "Other plots: dashed line is adjusted P = 0.05.",
 "The all-54-tissues figure includes nonsignificant tissues; adjusted P = 1 is plotted at zero, with ordering by raw P.",
 "Specific (54) and general tissues are corrected separately. No rich factor is inferred.",
 "MSigDB version is not established by these output tables; consult the original FUMA job metadata.",
 "Metascape import status is recorded in metascape/status.txt. deTS methods and results are saved under deTS/. Whole-blood eQTL mapping is separate from expression enrichment.",
 paste("Input:", normalizePath(input, winslash = "/"))), file.path(out, "figure_notes.txt"))
# ---- Metascape: candidate export and native result import -------------------
# Submit candidates_ensembl.txt at https://metascape.org/ (human -> human).
# Extract the downloaded ZIP under metascape/results, retaining Enrichment_GO/
# _FINAL_GO.csv. Alternatively set FLT_METASCAPE_RESULTS to that exact CSV.
# Do not recompute BH using the exported, prefiltered enrichment table.
meta_dir <- file.path(out, "metascape")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
mapped <- read.delim(file.path(input, "genes.txt"), check.names = FALSE)
gene_stats <- read.delim(file.path(input, "magma.genes.out"), check.names = FALSE)
stopifnot(all(c("ensg", "symbol") %in% names(mapped)),
          all(c("GENE", "SYMBOL", "P") %in% names(gene_stats)))
clean_id <- function(x) sub("\\.[0-9]+$", "", trimws(as.character(x)))
mapped$ensg <- clean_id(mapped$ensg)
gene_stats$GENE <- clean_id(gene_stats$GENE)
valid <- is.finite(gene_stats$P) & gene_stats$P >= 0 & gene_stats$P <= 1
gene_stats <- gene_stats[valid, ]
stopifnot(!anyDuplicated(gene_stats$GENE))
gene_stats$P_Bonferroni <- p.adjust(gene_stats$P, "bonferroni")
gene_stats$P_BH <- p.adjust(gene_stats$P, "BH")
gene_method <- Sys.getenv("FLT_METASCAPE_GENE_CORRECTION", "bonferroni")
if (!gene_method %in% c("bonferroni", "BH")) stop("Gene correction must be bonferroni or BH")
gene_p <- if (gene_method == "BH") gene_stats$P_BH else gene_stats$P_Bonferroni
m_ids <- unique(mapped$ensg[!is.na(mapped$ensg) & nzchar(mapped$ensg)])
g_ids <- unique(gene_stats$GENE[gene_p < alpha])
ids <- sort(union(m_ids, g_ids))
symbols <- mapped$symbol[match(ids, mapped$ensg)]
missing_symbol <- is.na(symbols) | !nzchar(symbols)
symbols[missing_symbol] <- gene_stats$SYMBOL[match(ids[missing_symbol], gene_stats$GENE)]
candidates <- data.frame(ensembl = ids, symbol = symbols,
                         FUMA_mapped = ids %in% m_ids, MAGMA_significant = ids %in% g_ids)
write.csv(candidates, file.path(meta_dir, "candidate_gene_evidence.csv"), row.names = FALSE)
writeLines(ids, file.path(meta_dir, "candidates_ensembl.txt"))
writeLines(sort(unique(symbols[!is.na(symbols) & nzchar(symbols)])), file.path(meta_dir, "candidates_symbols.txt"))
write.csv(gene_stats, file.path(meta_dir, "MAGMA_gene_statistics.csv"), row.names = FALSE)
counts <- data.frame(mapped = length(m_ids), magma_significant = length(g_ids),
                     intersection = length(intersect(m_ids, g_ids)), union = length(ids),
                     gene_tests = nrow(gene_stats), gene_correction = gene_method, alpha = alpha)
write.csv(counts, file.path(meta_dir, "candidate_counts.csv"), row.names = FALSE)
print(counts)
meta_default <- file.path(meta_dir, "result/Enrichment_GO/_FINAL_GO.csv")
if (!file.exists(meta_default)) meta_default <- file.path(meta_dir, "results/Enrichment_GO/_FINAL_GO.csv")
meta_file <- Sys.getenv("FLT_METASCAPE_RESULTS", meta_default)
writeLines(c("Metascape workflow (not Enrichr)",
 "1. Submit candidates_ensembl.txt to https://metascape.org/; choose Homo sapiens as input and analysis species.",
 paste("Candidate definition: all FUMA mapped genes union MAGMA", gene_method, "P < 0.05 genes."),
 "2. Record run date, database versions, background, submitted and successfully recognized gene counts in the job records.",
 "Use the documented Metascape default genome background unless specifying a justified custom assay universe; record your choice.",
 "3. Download the complete result ZIP and extract it under metascape/results.",
 "4. Rerun this script. FLT_METASCAPE_RESULTS can point to Enrichment_GO/_FINAL_GO.csv.",
 "The import uses native Log(q-value), Enrichment and #GeneInGOAndHitList values.",
 "It selects the lowest-q significant term per GroupID (ties broken by raw P), then plots up to 20 groups.",
 "The script does not submit a web job or substitute Enrichr for Metascape.",
 "Do not use Excel group-summary rows as term-level enrichment results.",
 "References: https://metascape.org/blog/?p=122 ; https://metascape.org/blog/?p=252"),
 file.path(meta_dir, "README.txt"))
if (file.exists(meta_file)) {
  meta <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"GroupID" %in% names(meta) && "GROUP_ID" %in% names(meta)) meta$GroupID <- meta$GROUP_ID
  needed <- c("Description", "GroupID", "Log(q-value)", "LogP", "Enrichment", "#GeneInGOAndHitList")
  absent <- setdiff(needed, names(meta))
  if (length(absent)) stop("Metascape CSV missing columns: ", paste(absent, collapse = ", "),
                           "; use the native full Enrichment_GO/_FINAL_GO.csv")
  numeric_cols <- needed[3:6]
  for (nm in numeric_cols) meta[[nm]] <- suppressWarnings(as.numeric(meta[[nm]]))
  if (any(!is.finite(meta[["Log(q-value)"]])) || any(meta[["Log(q-value)"]] > 0))
    stop("Invalid native log10 q values in Metascape results")
  if (any(!is.finite(meta$Enrichment) | meta$Enrichment <= 0) ||
      any(!is.finite(meta[["#GeneInGOAndHitList"]]) | meta[["#GeneInGOAndHitList"]] <= 0))
    stop("Invalid enrichment factors or hit counts")
  if (anyNA(meta$GroupID) || any(!nzchar(as.character(meta$GroupID)))) stop("Missing Metascape group IDs")
  meta$neg_log10_q <- -meta[["Log(q-value)"]]
  meta$BH_q <- 10^meta[["Log(q-value)"]]
  # These counts are native Metascape analysis counts, not inferred from the
  # number of genes appearing in significant terms.
  native_counts <- unique(meta[c("#TotalGeneInLibrary", "#GeneInHitList")])
  write.csv(native_counts, file.path(meta_dir, "Metascape_native_gene_counts.csv"), row.names = FALSE)
  write.csv(meta, file.path(meta_dir, "Metascape_all_imported.csv"), row.names = FALSE)
  selected <- meta[meta[["Log(q-value)"]] < log10(alpha), ]
  selected <- selected[order(selected[["Log(q-value)"]], selected$LogP), ]
  selected <- selected[!duplicated(selected$GroupID), ]
  import_summary <- data.frame(imported_terms = nrow(meta),
                               BH_significant_terms = sum(meta$BH_q < alpha),
                               BH_significant_clusters = nrow(selected),
                               displayed_clusters = min(20, nrow(selected)))
  write.csv(import_summary, file.path(meta_dir, "Metascape_result_summary.csv"), row.names = FALSE)
  print(import_summary)
  print(head(selected[c("Description", "Enrichment", "#GeneInGOAndHitList", "BH_q")], 10))
  write.csv(selected, file.path(meta_dir, "Metascape_significant_cluster_representatives.csv"), row.names = FALSE)
  show <- head(selected, 20)
  if (nrow(show)) {
    show$label <- paste0(show$Description, " [", show$GroupID, "]")
    show$label <- vapply(show$label, function(x) paste(strwrap(x, 55), collapse = "\n"), character(1))
    show$label <- factor(show$label, levels = rev(show$label))
    p_meta <- ggplot(show, aes(Enrichment, label, size = `#GeneInGOAndHitList`, colour = neg_log10_q)) +
      geom_point(alpha = .9) + scale_size_area(max_size = 9) +
      scale_colour_gradient(low = "#80B1D3", high = "#D95F5F") +
      labs(title = "Metascape pathway enrichment", subtitle = "Chicken liking x MDD | FUMA mapped + MAGMA significant genes",
           x = "Enrichment factor", y = NULL, size = "Hit genes", colour = expression(-log[10](q)),
           caption = "Native Metascape BH q < 0.05; one representative per cluster; up to 20 clusters shown.") +
      theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold"))
    save_plot(p_meta, "Metascape_pathways", 12, max(5, 2 + nrow(show) * .42))
  }
  writeLines(paste("Imported:", normalizePath(meta_file), "Significant clusters:", nrow(selected)), file.path(meta_dir, "status.txt"))
} else {
  writeLines("PENDING: candidate lists prepared; native Metascape results not found. No Metascape enrichment has been run by this script.", file.path(meta_dir, "status.txt"))
  message("Metascape candidates exported; awaiting native result file: ", meta_file)
}
# ---- deTS candidate-list tissue enrichment ---------------------------------
# Native references bundled in deTS; these are NOT the FUMA GTEx v8 panel.
# Match exact supplied symbols, report unmatched symbols, and avoid ambiguous
# automatic alias replacement. Both references use the same candidate list.
dets_dir <- file.path(out, "deTS")
dir.create(dets_dir, recursive = TRUE, showWarnings = FALSE)
dets_ratio <- 0.05
dets_query <- sort(unique(trimws(symbols[!is.na(symbols) & nzchar(trimws(symbols))])))
if (length(dets_query) < 20) stop("deTS requires at least 20 candidate symbols")
writeLines(dets_query, file.path(dets_dir, "input_symbols.txt"))
ref_env <- new.env()
data("GTEx_t_score", package = "deTS", envir = ref_env)
data("ENCODE_z_score", package = "deTS", envir = ref_env)
refs <- list(GTEx = ref_env$GTEx_t_score, ENCODE = ref_env$ENCODE_z_score)
dets_plots <- list()
dets_counts <- list()
for (panel in names(refs)) {
  score <- refs[[panel]]
  stopifnot(!anyDuplicated(rownames(score)), !anyDuplicated(colnames(score)))
  matched <- intersect(dets_query, rownames(score))
  if (length(matched) < 20) stop(panel, ": fewer than 20 matched candidate genes")
  writeLines(matched, file.path(dets_dir, paste0(panel, "_matched_symbols.txt")))
  writeLines(setdiff(dets_query, matched), file.path(dets_dir, paste0(panel, "_unmatched_symbols.txt")))
  matching <- candidates
  matching$matched_reference <- matching$symbol %in% rownames(score)
  write.csv(matching, file.path(dets_dir, paste0(panel, "_candidate_matching.csv")), row.names = FALSE)
  raw <- deTS::tsea.analysis(dets_query, score, ratio = dets_ratio, p.adjust.method = "none")
  # Independently recover contingency counts for interpretation and verification.
  result <- do.call(rbind, lapply(seq_len(ncol(score)), function(j) {
    values <- as.numeric(score[, j])
    cutoff <- quantile(values, 1 - dets_ratio, na.rm = TRUE)
    tissue_genes <- rownames(score)[which(values > cutoff)]
    hits <- intersect(matched, tissue_genes)
    a <- length(hits); b <- length(matched) - a
    c <- length(tissue_genes) - a; d <- nrow(score) - a - b - c
    check_p <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")$p.value
    stopifnot(isTRUE(all.equal(as.numeric(raw[j, 1]), check_p, tolerance = 1e-12)))
    data.frame(reference = panel, tissue = colnames(score)[j], P = raw[j, 1],
               hit_genes = a, matched_candidates = length(matched),
               tissue_specific_genes = length(tissue_genes), background_genes = nrow(score),
               enrichment_factor = (a / length(matched)) / (length(tissue_genes) / nrow(score)),
               hit_symbols = paste(hits, collapse = ";"))
  }))
  result$P_BH <- p.adjust(result$P, "BH")
  result$P_Bonferroni <- p.adjust(result$P, "bonferroni")
  result <- result[order(result$P), ]
  write.csv(result, file.path(dets_dir, paste0("deTS_", panel, "_all.csv")), row.names = FALSE)
  write.csv(result[result$P_BH < alpha, ], file.path(dets_dir, paste0("deTS_", panel, "_BH_significant.csv")), row.names = FALSE)
  dets_counts[[panel]] <- data.frame(reference = panel, candidate_ensembl = length(ids),
      input_symbols = length(dets_query), matched_symbols = length(matched),
      unmatched_symbols = length(setdiff(dets_query, matched)), reference_genes = nrow(score),
      tissue_tests = ncol(score), BH_significant = sum(result$P_BH < alpha),
      Bonferroni_significant = sum(result$P_Bonferroni < alpha))
  plot_data <- result
  plot_data$tissue <- factor(plot_data$tissue, levels = rev(plot_data$tissue))
  plot_data$significance <- ifelse(plot_data$P_BH < alpha, "BH < 0.05", "Not significant")
  p_dets <- ggplot(plot_data, aes(-log10(pmax(P_BH, .Machine$double.xmin)), tissue)) +
    geom_vline(xintercept = -log10(alpha), linetype = 2, colour = "#A4ABB2") +
    geom_point(aes(size = hit_genes, colour = significance), alpha = .9) +
    scale_size_area(max_size = 7) +
    scale_colour_manual(values = c("BH < 0.05" = "#237C8C", "Not significant" = "#A4ABB2")) +
    labs(title = paste("deTS:", panel, "tissue enrichment"),
         subtitle = sprintf("%d/%d symbols matched | %d tissues | top 5%% specificity scores", length(matched), length(dets_query), ncol(score)),
         x = expression(-log[10](P[BH])), y = NULL, size = "Hit genes", colour = NULL,
         caption = "One-sided Fisher exact test; BH correction within reference. Native deTS reference, not FUMA GTEx v8.") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"), legend.position = "bottom", axis.text.y = element_text(size = 9))
  save_plot(p_dets, paste0("deTS_", panel, "_all_tissues"), 11, max(6, 2 + nrow(result) * .24))
  dets_plots[[panel]] <- p_dets
  print(head(result[c("tissue", "hit_genes", "enrichment_factor", "P_BH")], 5))
}
dets_summary <- do.call(rbind, dets_counts)
write.csv(dets_summary, file.path(dets_dir, "deTS_matching_and_test_summary.csv"), row.names = FALSE)
print(dets_summary)
save_plot(patchwork::wrap_plots(dets_plots, ncol = 2), "deTS_GTEx_ENCODE", 19,
          max(8, 2 + max(vapply(refs, ncol, integer(1))) * .24))
writeLines(c(paste("deTS version:", as.character(packageVersion("deTS"))),
 "Input: same FUMA-mapped union MAGMA-significant candidate genes as Metascape; deduplicated FUMA/MAGMA symbols.",
 "Exact symbol matching; unmatched symbols are exported separately for each panel. No automatic alias substitution.",
 "References: package-bundled GTEx_t_score and ENCODE_z_score; not the FUMA GTEx v8 reference.",
 "Tissue-specific genes: scores strictly above the 95th percentile, following deTS implementation (ratio=0.05).",
 "Background: all genes in each reference matrix; candidate list is intersected with this background.",
 "Test: one-sided Fisher exact test (greater), NOT permutation testing.",
 "BH and Bonferroni adjustments are computed across all tissues separately within each reference.",
 "All tissue results are plotted, including nonsignificant results. Dot area indicates hit gene count.",
 "Raw deTS P values verified against independently constructed Fisher contingency tables.",
 "Reference: https://search.r-project.org/CRAN/refmans/deTS/html/tsea.analysis.html"),
 file.path(dets_dir, "methods.txt"))
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
print(summary)
message("Saved plots and tables to: ", out)

