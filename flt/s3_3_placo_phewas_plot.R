# PheWAS 绘图：先运行 s3_3_placo_phewas.R，再直接运行本脚本。
# 只读取分类结果，不联网、不重新分类；修改字典后先重跑分析脚本。
# 必需 data.table；wordcloud2、htmlwidgets 仅用于 HTML 词云。
# install.packages(c("data.table", "wordcloud2", "htmlwidgets"))
library(data.table)

dir0 <- "D:/"
project <- file.path(dir0, "projects/flt")
data_dir <- file.path(project, "analysis/03phewas")
infile <- file.path(data_dir, "lead_phewas_results.tsv")
outdir <- file.path(data_dir, "plots")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 一、绘图参数
level <- "both"       # both / trait / class：两套图 / 标准表型 / 类别
metric <- "n_loci"    # 词云指标：n_loci / n_leads / n_records
top_n <- Inf          # 不截断并列词条；在下面用关联记录数阈值筛选
trait_min_records <- 20L # 标准表型至少20条记录
draw_cloud <- TRUE    # TRUE 生成彩色词云，FALSE 只生成柱状图
# 二、读取分析结果
if (!file.exists(infile)) stop("请先运行 s3_3_placo_phewas.R：", infile)
plot_audit <- list()
write_tsv <- function(x, name) plot_audit[[name]] <<- copy(x)

d <- as.data.table(read.delim(infile, fileEncoding="UTF-8-BOM", check.names=FALSE, stringsAsFactors=FALSE, na.strings=c("", "NA")))
stopifnot(all(c("trait", "pair", "lead_snp", "locus", "source", "p") %in% names(d)))
d <- d[is.finite(p) & p >= 0 & p < 1e-5]
if (!nrow(d)) stop("No associations with P < 1e-5.")
if (anyNA(d$locus) || anyNA(d$pair)) stop("FUMA locus and pair are required for counting.")
required <- c("canonical_trait", "trait_class", "analysis_scope", "input_trait")
if (!all(required %in% names(d))) stop("缺少分类字段，请先运行 s3_3_placo_phewas.R 。")

# Keep non-phenotype and index-trait records in the full output, with separate scopes.
main <- d[analysis_scope == "phenotype" & !input_trait]
if (!nrow(main)) stop("No classified phenotypes for the main plots; see classification_needs_review.tsv.")
coverage_file <- file.path(data_dir, "source_coverage.tsv")
coverage_note <- "Query coverage unavailable"
if (file.exists(coverage_file)) {
  coverage <- fread(coverage_file)
  write_tsv(coverage[, .(n_snps = .N), by = .(source, status)], "source_coverage_summary.tsv")
  coverage_note <- "Database coverage incomplete"
  if (nrow(coverage) && setequal(unique(coverage$source), c("GWASCatalog", "OpenGWAS", "GWASAtlas")) &&
      all(coverage$status %in% "complete") &&
      all(coverage[, .(n_sources = uniqueN(source)), by = query_snp]$n_sources == 3L)) {coverage_note <- "All three sources queried for each lead SNP"}
}

plot_levels <- if (level == "both") c("class", "trait") else level
metric_labels <- c(n_loci = "Number of distinct FUMA genomic loci", n_leads = "Number of distinct lead SNPs", n_records = "Number of database records")
have_cloud <- draw_cloud && requireNamespace("wordcloud2", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)
if (draw_cloud && !have_cloud) message("For optional HTML word clouds, install wordcloud2 and htmlwidgets.")

# 三、绘图并保存：柱状图 PDF/PNG，词云 HTML
draw_bar <- function(z, pair, by_level, denominator) {
  labels <- vapply(rev(z$word), function(x) paste(strwrap(x, width = 48), collapse = "\n"), character(1))
  values <- rev(z$n_loci)
  par(mar = c(5, min(25, max(12, max(nchar(sub("\n.*", "", labels))) * .47)), 5, 2), las = 1)
  at <- barplot(values, horiz = TRUE, names.arg = labels, cex.names = .83, col = "#3C7180", border = NA, xlim = c(0, max(values) * 1.17), xlab = metric_labels[["n_loci"]], las = 1)
  text(values + max(values) * .018, at, labels = values, pos = 4, cex = .83, xpd = NA)
  title(main = paste(pair, if (by_level == "class") "| Phenotype categories" else "| Canonical phenotypes"))
  mtext(sprintf("%d loci with any database hit at P < 1e-5 | %s", denominator, paste(coverage_note, selection_note, sep=" | ")), side = 3, line = .25, cex = .78)
  mtext("Main plot excludes input phenotypes, molecular traits, composite tests and administrative terms.", side = 1, line = 3.7, cex = .70, col = "#526073")
}

save_bar <- function(path, format, height, draw) {
  if (format == "pdf") pdf(path, width = 12, height = height, useDingbats = FALSE) else
    png(path, width = 12, height = height, units = "in", res = 260, bg = "white")
  on.exit(dev.off())
  draw()
}

for (by_level in plot_levels) {
  column <- if (by_level == "class") "trait_class" else "canonical_trait"
  counts <- main[, .(n_records = .N, n_leads = uniqueN(lead_snp), n_loci = uniqueN(locus), n_sources = uniqueN(source)), by = .(pair, word = get(column))]
  setorder(counts, pair, -n_loci, word)
  write_tsv(counts, paste0(by_level, "_main_counts.tsv"))
  for (pair_i in unique(counts$pair)) {
    z <- copy(counts[pair == pair_i])
    if (by_level == "trait") z <- z[n_records >= trait_min_records]
    if (!nrow(z)) {
      warning("No terms pass the display threshold for ", pair_i, " / ", by_level)
      next
    }
    selection_note <- if (by_level == "trait")
      paste0("n_records >= ", trait_min_records, "; ", nrow(z), " canonical phenotypes") else
      paste0("All ", nrow(z), " phenotype categories")
    write_tsv(z, paste0(by_level, "_display_counts.tsv"))
    message(pair_i, " / ", by_level, ": ", selection_note)
    stem <- paste(gsub("[^A-Za-z0-9_-]+", "_", pair_i), by_level, sep = "_")
    b <- if (is.infinite(top_n)) copy(z) else head(z, top_n)
    height <- max(5.5, nrow(b) * .37 + 2)
    n_loci_any <- uniqueN(d[pair == pair_i, locus])
    draw <- function() draw_bar(b, pair_i, by_level, n_loci_any)
    for (extension in c("png", "pdf")) {save_bar(file.path(outdir, paste0(stem, "_n_loci_barplot.", extension)), extension, height, draw)}
    if (have_cloud) {
      z[, freq := get(metric)]
      setorder(z, -freq, word)
      z <- if (is.infinite(top_n)) copy(z) else head(z, top_n)
      # 每张词云在实际展示的词条中计算分位数；相同频数同色。
      # 颜色与字号均编码当前 metric（默认 n_loci），不是 P 值或表型类别。
      q <- quantile(z$freq, probs = seq(0, 1, 0.1), names = FALSE)
      z[, color := fcase(freq > q[10], "#E64B35", freq > q[9], "#4DBBD5", freq > q[6], "#00A087", freq > q[3], "#3C5488", default = "#F39B7F")]
      write_tsv(z, paste0(stem, "_", metric, "_wordcloud_counts.tsv"))
      cloud <- wordcloud2::wordcloud2(as.data.frame(z[, .(word, freq)]), size = .65, minSize = 2, color = z$color, backgroundColor = "white", shape = "pentagon", rotateRatio = 0, shuffle = FALSE, widgetsize = c(1800, 1200))
      htmlwidgets::saveWidget(cloud, file.path(outdir, paste0(stem, "_", metric, "_wordcloud.html")), selfcontained = FALSE, title = paste(pair_i, by_level, metric, selection_note, coverage_note))
    }
  }
}
dir.create(file.path(data_dir, "cache"), recursive=TRUE, showWarnings=FALSE)
saveRDS(plot_audit, file.path(data_dir, "cache", "plot_audit.rds"))
message("Records: ", nrow(d), "; raw terms: ", uniqueN(d$trait), "; canonical terms: ", uniqueN(d$canonical_trait), "; unclassified terms: ", uniqueN(d[analysis_scope == "unclassified", trait]))
message("Outputs: ", outdir)
