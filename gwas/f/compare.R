#!/usr/bin/env Rscript


# 🚩 Select workflow
.compare_args <- commandArgs(trailingOnly = TRUE)
if (!length(.compare_args) || .compare_args[1L] %in% c("-h", "--help")) {
	cat("Modes: compare, ldsc\n")
	quit(status = 0L)
}
.compare_mode <- .compare_args[1L]
.compare_args <- .compare_args[ - 1L]
if (!.compare_mode %in% c('compare', 'ldsc')) stop("Unknown mode: ", .compare_mode)


# 🚩 compare
if (.compare_mode == "compare") {
# Multi-GWAS QC comparison of standardized summary statistics.
suppressPackageStartupMessages(library(data.table))
args <- .compare_args
if (!length(args) || any(args %in% c('-h', '--help'))) {
	cat('Usage: compare.sh compare --gwas-files A.gz,B.gz[,C.gz] [options]\n',
			'--output-dir DIR (default ./gwas_compare) --labels CSV --grch 37|38\n',
			'--compare-beta TRUE --compare-EAF TRUE\n',
			'--p-threshold 5e-8 --significant first|either|both (default first)\n',
			'Input: SNP CHR POS EA NEA P; BETA/EAF required when compared. Same build required.\n',
			'First GWAS is compared with each follower. Palindromic SNPs and duplicate\n',
			'allele/position keys are excluded from scatter plots and counted in QC.\n',
			'\nExamples:\n',
			'  compare.sh compare --gwas-files A.gz,B.gz,C.gz --labels A,B,C --output-dir qc\n',
			'  compare.sh compare --gwas-files A.gz,B.gz --significant either --p-threshold 5e-8\n',
			'  compare.sh compare --gwas-files A.gz,B.gz --compare-beta TRUE --compare-EAF TRUE\n',
			'Use compare.sh -h for full paths and output descriptions.\n')
	quit(status = 0)
}
opt <- list('output-dir' = 'gwas_compare', 'compare-beta' = 'TRUE',
				'compare-eaf' = 'TRUE', 'p-threshold' = '5e-8', significant = 'first')
allowed <- c(names(opt), 'gwas-files', 'labels', 'grch')
if (length(args) %% 2L) stop('Options require values')
for (i in seq(1, length(args), 2)) {
	k <- tolower(sub('^--', '', args[i]))
	if (!k %in% allowed) stop('Unknown option: ', args[i])
	opt[[k]] <- args[i + 1]
}
flag <- function(x) { if (!toupper(x) %in% c('TRUE', 'FALSE')) stop('Expected TRUE/FALSE'); toupper(x) == 'TRUE' }
cb <- flag(opt[['compare-beta']]); ce <- flag(opt[['compare-eaf']])
if (is.null(opt[['gwas-files']])) stop('--gwas-files required')
files <- trimws(strsplit(opt[['gwas-files']], ',', fixed = TRUE)[[1]])
if (length(files) < 2 || any(!file.exists(files))) stop('Provide at least two existing GWAS files')
labels <- if (is.null(opt$labels)) basename(files) else strsplit(opt$labels, ',', fixed = TRUE)[[1]]
if (length(labels) != length(files) || anyDuplicated(labels)) stop('Labels must be unique and match files')
threshold <- as.numeric(opt[['p-threshold']])
if (!is.finite(threshold) || threshold <= 0 || threshold > 1) stop('Invalid P threshold')
if (!opt$significant %in% c('first', 'either', 'both')) stop('Invalid significant selector')
if (!is.null(opt$grch) && !opt$grch %in% c('37', '38')) stop('Invalid GRCh')
builds <- unique(unlist(lapply(files, function(f) if(file.exists(paste0(f, '.grch'))) trimws(readLines(paste0(f, '.grch'))) else NULL)))
if(length(unique(c(builds, opt$grch))) > 1L) stop('Mixed GRCh builds: liftover before comparison')
out <- opt[['output-dir']]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
comp <- function(x) chartr('ACGT', 'TGCA', x)
read_gwas <- function(f) {
	# fread reads one genome-wide file at a time. gzip shell command is quoted.
	d <- if (grepl('\\.(gz|bgz)$', f)) fread(cmd = paste('gzip -cd --', shQuote(f)), showProgress = FALSE) else fread(f, showProgress = FALSE)
	setnames(d, toupper(sub('^#', '', names(d))))
	required <- c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'P', if(cb) 'BETA', if(ce) 'EAF')
	if (length(setdiff(required, names(d)))) stop(f, ': missing ', paste(setdiff(required, names(d)), collapse = ','))
	d <- d[, ..required]
	for (k in intersect(c('POS', 'P', 'EAF', 'BETA'), names(d))) set(d, j = k, value = suppressWarnings(as.numeric(d[[k]])))
	d[, CHR := toupper(sub('^CHR', '', toupper(as.character(CHR))))]
	d[CHR == 'X', CHR := '23']; d[CHR == 'Y', CHR := '24']; d[CHR %in% c('MT', 'M'), CHR := '25']
	d[, CHR := suppressWarnings(as.integer(CHR))]
	d[, `:=`(EA = toupper(EA), NEA = toupper(NEA))]
	d[, valid := is.finite(P) & P >= 0 & P <= 1 & is.finite(POS) & POS > 0 & POS == floor(POS) & CHR %in% 1:25 & !is.na(EA) & !is.na(NEA) & EA != NEA]
	invalid <- d[valid == FALSE, .N]; d <- d[valid == TRUE]; d[, valid := NULL]
	if(!nrow(d)) stop(f, ': no valid variants')
	# Canonicalize SNP complements; indels require literal matching (no strand inference).
	d[, pair := paste(pmin(EA, NEA), pmax(EA, NEA), sep = '/')]
	d[, pal := pair %in% c('A/T', 'C/G')]
	d[, canonical := pair]
	d[nchar(EA) == 1 & nchar(NEA) == 1 & grepl('^[ACGT]$', EA) & grepl('^[ACGT]$', NEA),
		canonical := pmin(pair, paste(pmin(comp(EA), comp(NEA)), pmax(comp(EA), comp(NEA)), sep = '/'))]
	d[, key := paste(CHR, POS, canonical, sep = ':')]
	d[, duplicate := duplicated(key) | duplicated(key, fromLast = TRUE)]
	attr(d, 'invalid') <- invalid
	d
}
# Keep significant identities first; follower-only significant sites are included on demand.
first <- read_gwas(files[1])
base <- first[!duplicate & !pal]
if(opt$significant %in% c('first', 'both')) base <- base[P <= threshold]
audit <- list(); summary <- list()
for (i in seq_along(files)) {
	d <- if (i == 1L) first else read_gwas(files[i])
	audit[[i]] <- data.table(file = files[i], label = labels[i], invalid_rows = attr(d, 'invalid'), valid_rows = nrow(d), duplicate_rows = sum(d$duplicate), palindromic_rows = sum(d$pal))
	if (i == 1L) { rm(first); gc(verbose = FALSE); next }
	if (!(cb || ce)) next
	if(opt$significant %in% c('first', 'both')) d <- d[key %in% base$key]
	x <- merge(base, d[!duplicate & !pal], by = 'key', suffixes = c('.first', '.other'))
	x <- x[switch(opt$significant, first = P.first <= threshold, either = P.first <= threshold | P.other <= threshold, both = P.first <= threshold & P.other <= threshold)]
	x[, flip := EA.first == NEA.other & NEA.first == EA.other |
			(nchar(EA.first) == 1 & nchar(NEA.first) == 1 & EA.first == comp(NEA.other) & NEA.first == comp(EA.other))]
	if (cb) x[flip == TRUE, BETA.other :=  - BETA.other]
	if (ce) x[flip == TRUE, EAF.other := 1 - EAF.other]
	tag <- sprintf('01_vs_%02d', i)
	matched_file <- file.path(out, paste0(tag, '.harmonized.tsv.gz'))
	if (nrow(x)) {
		fwrite(x, matched_file, sep = '\t')
	} else {
		# Some data.table versions omit the gzip trailer for zero-row tables.
		handle <- gzfile(matched_file, 'wt')
		writeLines(paste(names(x), collapse = '\t'), handle)
		close(handle)
	}
	for (v in c(if(cb) 'BETA', if(ce) 'EAF')) {
		a <- x[[paste0(v, '.first')]]; b <- x[[paste0(v, '.other')]]
		ok <- is.finite(a) & is.finite(b)
		if (v == 'EAF') ok <- ok & a >= 0 & a <= 1 & b >= 0 & b <= 1
		a <- a[ok]; b <- b[ok]
		r <- if(length(a) > 1 && sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_
		summary[[length(summary) + 1]] <- data.table(first = labels[1], other = labels[i], metric = v, matched_significant = nrow(x), n = length(a), flipped = sum(x$flip), r = r, mean_difference = if(length(a)) mean(b - a) else NA_real_)
		png(file.path(out, paste0(tag, '.', v, '.png')), width = 1400, height = 1400, res = 180)
		if (length(a)) {
			lim <- range(c(a, b)); if(diff(lim) == 0) lim <- lim + c( - .01, .01)
			plot(a, b, pch = 16, cex = .45, col = adjustcolor('#286b9e', alpha.f = .35), xlim = lim, ylim = lim,
					 xlab = paste(labels[1], v), ylab = paste(labels[i], v), main = sprintf('n=%d; r=%.3f', length(a), r)); abline(0, 1, col = 'firebrick')
		} else { plot.new(); title(main = paste(v, ': no eligible matched significant variants')) }
		dev.off()
	}
}
fwrite(rbindlist(audit), file.path(out, 'input_qc.tsv'), sep = '\t')
if(length(summary)) fwrite(rbindlist(summary), file.path(out, 'comparison_qc.tsv'), sep = '\t')
writeLines(c(paste('GRCh:', if(is.null(opt$grch)) 'unspecified; caller must ensure same build' else opt$grch),
					paste('Significance:', opt$significant, 'P <=', threshold), capture.output(sessionInfo())), file.path(out, 'compare.log'))
cat('Comparison written to ', normalizePath(out), '\n', sep = '')

}


# 🚩 ldsc
if (.compare_mode == "ldsc") {
# Same pairwise heatmap idea as 0f/plotting.R::plot_rg, using installed ggplot2.
# Keep raw estimates (including rg > 1) in labels/tables; saturate only the colors.
suppressPackageStartupMessages(library(ggplot2))
out <- .compare_args[1]
rg <- read.delim(file.path(out, 'rg.tsv'), check.names = FALSE, na.strings = c('NA', 'nan'))
h2 <- read.delim(file.path(out, 'h2.tsv'), check.names = FALSE)
traits <- h2$trait
m <- expand.grid(p1 = traits, p2 = traits, stringsAsFactors = FALSE)
m$rg <- NA_real_; m$p <- NA_real_
for (i in seq_len(nrow(rg))) {
	hit <- (m$p1 == rg$p1[i] & m$p2 == rg$p2[i]) | (m$p1 == rg$p2[i] & m$p2 == rg$p1[i])
	if (rg$status[i] == 'ESTIMATED') {
		m$rg[hit] <- rg$rg[i]; m$p[hit] <- rg$p[i]
	}
}
available <- unique(c(rg$p1[rg$status == 'ESTIMATED'], rg$p2[rg$status == 'ESTIMATED'],
					h2$trait[h2$status == 'ESTIMATED']))
diag <- m$p1 == m$p2 & m$p1 %in% available
m$rg[diag] <- 1; m$p[diag] <- 0
m$label <- ifelse(is.finite(m$rg), sprintf('%.2f', m$rg), 'NA')
m$label[is.finite(m$p) & m$p < 0.05 & m$p1 != m$p2] <- paste0(m$label[is.finite(m$p) & m$p < 0.05 & m$p1 != m$p2], '*')
m$p1 <- factor(m$p1, levels = traits)
m$p2 <- factor(m$p2, levels = rev(traits))
p <- ggplot(m, aes(p1, p2, fill = rg)) +
	geom_tile(color = 'white', linewidth = 0.4) + geom_text(aes(label = label), size = 3.3) +
	scale_fill_gradient2(low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 0,
				limits = c( - 1, 1), oob = scales::squish, na.value = 'grey90', name = 'rg') +
	coord_fixed() + labs(x = NULL, y = NULL, title = 'LDSC genetic correlation',
		subtitle = 'Chromosomes 1–22 · * nominal P < 0.05 · NA: unavailable',
		caption = 'Raw estimates are shown; color scale is capped at ±1. Diagonal = 1 by definition for available traits.') +
	theme_minimal(base_size = 12) + theme(panel.grid = element_blank(),
		axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = 'bold'))
tmp <- file.path(out, 'rg.tmp.png')
ggsave(tmp, p, width = max(8, 3 + 0.65 * length(traits)), height = max(7, 2 + 0.65 * length(traits)), dpi = 180, bg = 'white')
if (!file.rename(tmp, file.path(out, 'rg.png'))) stop('Could not publish rg.png')

}
