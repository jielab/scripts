#!/usr/bin/env Rscript


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 IBDmix: genome-wide burden + selected-locus summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic genome-wide window burden:
# Rscript ibdmix.R --in len1000.lod5.segments.tsv.gz --outdir report
# With sample groups:
# Rscript ibdmix.R --in len1000.lod5.segments.tsv.gz --outdir report \
# --sample_file samples_v3.ALL.panel --stats_by_group super_pop,gender
# With focal loci:
# Rscript ibdmix.R --in len1000.lod5.segments.tsv.gz --outdir report \
# --sample_file samples_v3.ALL.panel --stats_by_group super_pop,gender \
# --selected_loci gu.loci.bed
if (!requireNamespace("data.table", quietly = TRUE)) {
	stop("Required R package not installed: data.table. Install it with install.packages('data.table').")
}
suppressPackageStartupMessages(library(data.table))

if (!exists("write_xlsx", mode = "function")) {
	write_xlsx <- function(x, path) {
		if (requireNamespace("writexl", quietly = TRUE)) {
			writexl::write_xlsx(x, path)
		} else {
			d <- paste0(path, ".csv_sheets")
			dir.create(d, recursive = TRUE, showWarnings = FALSE)
			for (nm in names(x)) fwrite(x[[nm]], file.path(d, paste0(nm, ".csv")))
			message("writexl is not installed; wrote CSV sheets under: ", d)
		}
	}
}

phe_fun_file <- c("D:/scripts/0f/0phe.f.R", "/mnt/d/scripts/0f/0phe.f.R", "/work/sph-huangj/scripts/0f/0phe.f.R")
phe_fun_file <- phe_fun_file[file.exists(phe_fun_file)][1]
if (!is.na(phe_fun_file)) source(phe_fun_file)
# Inlined helper functions for ibdmix.R

get_arg <- function(key, default = NA_character_) {
	i <- match(key, args)
	if (is.na(i) || i == length(args)) return(default)
	args[i + 1]
}

find_input <- function(f, outdir) {
	bn <- basename(f)
	bn2 <- sub("\\.tsv\\.gz$", ".segments.tsv.gz", bn)
	cand <- unique(c(
		f,
		file.path(outdir, bn),
		sub("\\.tsv\\.gz$", ".segments.tsv.gz", f),
		file.path(outdir, bn2),
		file.path(dirname(f), "report", bn),
		file.path(dirname(f), "report", bn2)
	))
	cand[file.exists(cand)][1]
}

out_file <- function(suffix) file.path(outdir, suffix)

out_dir <- function(dirname) file.path(outdir, dirname)

clean_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

norm_chrom <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

chr_label <- function(x) paste0("chr", norm_chrom(x))

# IBDmix summary.sh writes BED-like intervals: length = end - start.
# Use this convention consistently when summing segment burden.
segment_width <- function(start, end) pmax(0L, as.integer(end) - as.integer(start))

sum_num <- function(x) sum(as.numeric(x), na.rm = TRUE)

reduce_segments <- function(d) {
	if (nrow(d) == 0) return(d)
	key_cols <- intersect(c("ID", "anc", "chrom", "sample_set", "super_pop"), names(d))
	d <- copy(d)
	setorderv(d, c(key_cols, "start", "end"))
	d[, prev_max_end := shift(cummax(end), fill = NA_integer_), by = key_cols]
	d[, new_interval := is.na(prev_max_end) | start > prev_max_end]
	d[, interval_id := cumsum(new_interval), by = key_cols]
	out <- d[, .(
		start = min(start, na.rm = TRUE),
		end = max(end, na.rm = TRUE),
		n_raw_segments = .N,
		raw_bp_before_reduce = sum_num(seg_len),
		slod = max(slod, na.rm = TRUE),
		sites = sum(sites, na.rm = TRUE),
		positive_lods = sum(positive_lods, na.rm = TRUE),
		negative_lods = sum(negative_lods, na.rm = TRUE)
	), by = c(key_cols, "interval_id")]
	out[, `:=`(
		length = segment_width(start, end),
		seg_len = segment_width(start, end),
		bp_removed_by_reduce = raw_bp_before_reduce - segment_width(start, end)
	)]
	out[, interval_id := NULL]
	setcolorder(out, intersect(c("ID", "chrom", "start", "end", "length", "slod", "sites", "positive_lods", "negative_lods", "sample_set", "super_pop", "anc", "seg_len", "n_raw_segments", "raw_bp_before_reduce", "bp_removed_by_reduce"), names(out)))
	setorderv(out, intersect(c("ID", "anc", "chrom", "start", "end"), names(out)))
	out[]
}

count_do <- function(f) {
	cmd <- paste("gzip -dc", shQuote(f), "| awk 'BEGIN{n=0} $1==\"DO\"{n++} END{print n}'")
	as.integer(system(cmd, intern = TRUE))
}

chrom_rank <- function(x) {
	r <- match(x, chrom_order)
	other <- match(x, sort(unique(x)))
	ifelse(is.na(r), length(chrom_order) + other, r)
}

make_windows <- function(d, w) {
	chr.max <- d[, .(chr_end = max(end, na.rm = TRUE)), by = chrom]
	chr.max[, nwin := ceiling(chr_end / w)]
	wins <- chr.max[, .(win_id = seq_len(nwin)), by = .(chrom, chr_end)]
	wins[, start := (win_id - 1L) * w + 1L]
	wins[, end := pmin(win_id * w, chr_end)]
	wins[, window := paste0(chrom, ":", start, "-", end)]
	wins[, chr_rank := chrom_rank(chrom)]
	setorder(wins, chr_rank, start)
	wins[, genome_pos := cumsum(as.numeric(c(1L, head(end - start + 1L, -1L))))]
	wins[]
}

draw_linear_gw_plot <- function(file) {
	png(file, width = 2600, height = 1800, units = "px", res = 220)
	op <- par(mfrow = c(3, 2), mar = c(4.2, 4.4, 3.4, 1.2), oma = c(0, 0, 2, 0))
	for (a in plot.names) {
		if (is.na(a)) {
			plot.new()
			next
		}
		d <- win.plot[anc == a]
		if (nrow(d) == 0) {
			plot.new()
			next
		}
		main <- ifelse(a == "ANY_REF", "Any archaic reference", a)
		plot(d$genome_pos, d$burden_percent, pch = 20, cex = 0.35, xaxt = "n",
			xlab = "Chromosome", ylab = "Samples with segment (%)", main = main)
		axis(1, at = chr.mid$mid, labels = chr.mid$chrom, cex.axis = 0.55)
		q99 <- as.numeric(quantile(d$burden_percent, 0.99, na.rm = TRUE))
		abline(h = q99, lty = 2, lwd = 1)
		mtext(paste0("99th pct = ", round(q99, 2), "%"), side = 3, line = 0.15, cex = 0.75)
	}
	mtext(paste0("Genome-wide IBDmix window burden, window = ", format(window_size, big.mark = ","), " bp"),
		outer = TRUE, cex = 1.1, font = 2)
	par(op)
	dev.off()
}

draw_circular_panel <- function(d, main) {
	r0 <- 0.22
	r1 <- 0.92
	plot.new()
	plot.window(xlim = c(-1.15, 1.15), ylim = c(-1.15, 1.15), asp = 1)
	q99 <- as.numeric(quantile(d$burden_percent, 0.99, na.rm = TRUE))
	if (!is.finite(q99) || q99 <= 0) q99 <- max(d$burden_percent, na.rm = TRUE)
	if (!is.finite(q99) || q99 <= 0) q99 <- 1
	th <- pi / 2 - 2 * pi * (d$genome_start - 1) / genome_total_bp
	y <- pmin(d$burden_percent, q99)
	r <- r0 + (r1 - r0) * y / q99
	cols <- ifelse(as.integer(d$chr_rank) %% 2L == 0L, "grey35", "grey60")
	segments(r0 * cos(th), r0 * sin(th), r * cos(th), r * sin(th), col = cols, lwd = 0.55)
	points(r * cos(th), r * sin(th), pch = 20, cex = 0.18, col = "black")
# inner baseline and outer 99th-percentile ring
	theta_grid <- seq(0, 2 * pi, length.out = 720)
	lines(r0 * cos(theta_grid), r0 * sin(theta_grid), col = "grey75", lwd = 0.6)
	lines(r1 * cos(theta_grid), r1 * sin(theta_grid), col = "grey25", lty = 2, lwd = 0.7)
# chromosome boundaries and labels
	for (i in seq_len(nrow(chr.bounds))) {
		b <- chr.bounds$chr_start[i]
		tb <- pi / 2 - 2 * pi * (b - 1) / genome_total_bp
		segments((r0 - 0.02) * cos(tb), (r0 - 0.02) * sin(tb), (r1 + 0.03) * cos(tb), (r1 + 0.03) * sin(tb),
			col = "grey70", lwd = 0.4)
		if (i %% 2L == 1L || nrow(chr.bounds) <= 16L) {
			mt <- pi / 2 - 2 * pi * (chr.bounds$chr_mid[i] - 1) / genome_total_bp
			text(1.05 * cos(mt), 1.05 * sin(mt), labels = chr.bounds$chrom[i], cex = 0.42)
		}
	}
	text(0, 0.05, main, font = 2, cex = 0.9)
	text(0, -0.08, paste0("99th pct\n", round(q99, 2), "%"), cex = 0.62)
}

draw_circular_gw_plot <- function(file) {
	png(file, width = 2500, height = 2500, units = "px", res = 260)
	op <- par(mfrow = c(3, 2), mar = c(1.2, 1.2, 2.2, 1.2), oma = c(0, 0, 3, 0))
	for (a in plot.names) {
		if (is.na(a)) {
			plot.new()
			next
		}
		d <- win.plot[anc == a]
		if (nrow(d) == 0) {
			plot.new()
			next
		}
		main <- ifelse(a == "ANY_REF", "Any archaic reference", a)
		draw_circular_panel(d, main)
	}
	mtext(paste0("Genome-wide IBDmix window burden, circular density, window = ",
		format(window_size, big.mark = ","), " bp"), outer = TRUE, cex = 1.05, font = 2)
	par(op)
	dev.off()
}

parse_loci <- function(x) {
	if (is.na(x) || !nzchar(x)) return(data.table())
	if (file.exists(x)) {
		lines <- readLines(x, warn = FALSE)
		lines <- trimws(lines)
		lines <- lines[nzchar(lines) & !grepl("^#", lines)]
		fields <- strsplit(lines, "\\s+")
		if (!length(fields) || any(lengths(fields) < 3)) {
			stop("Bad selected_loci file: expected at least 3 whitespace-separated columns: ", x)
		}
		max_cols <- max(lengths(fields))
		bed <- rbindlist(lapply(fields, function(z) {
			length(z) <- max_cols
			as.list(setNames(z, paste0("V", seq_len(max_cols))))
		}), fill = TRUE)
		bed[, `:=`(
			chrom = sub("^chr", "", V1, ignore.case = TRUE),
			locus_start = suppressWarnings(as.integer(gsub(",", "", V2))),
			locus_end = suppressWarnings(as.integer(gsub(",", "", V3))),
			label = if ("V4" %in% names(bed)) trimws(V4) else NA_character_
		)]
		bed <- bed[!is.na(locus_start) & !is.na(locus_end) & nzchar(chrom)]
		if (nrow(bed) == 0) stop("No valid selected loci in file: ", x)
		bed[locus_end < locus_start, c("locus_start", "locus_end") := .(locus_end, locus_start)]
		bed[is.na(label) | label == "", label := paste0("selected", .I)]
		out <- bed[, .(
			locus_id = make.unique(clean_name(label), sep = "_"),
			locus = paste0(chrom, ":", locus_start, "-", locus_end),
			chrom,
			locus_start,
			locus_end
		)]
		out[, locus_len := pmax(0L, locus_end - locus_start)]
		return(out)
	}
	x <- gsub("[\\x{2010}-\\x{2015}\\x{2212}]", "-", x, perl = TRUE)
	parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
	parts <- parts[nzchar(parts)]
	out <- rbindlist(lapply(seq_along(parts), function(i) {
		p <- parts[i]
		m <- regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", p)
		r <- regmatches(p, m)[[1]]
		if (length(r) != 4) stop("Bad locus format: ", p, ". Expected chr:start-end;chr:start-end")
		data.table(
			locus_id = paste0("locus", i),
			locus = p,
			chrom = sub("^chr", "", r[2], ignore.case = TRUE),
			locus_start = as.integer(gsub(",", "", r[3])),
			locus_end = as.integer(gsub(",", "", r[4]))
		)
	}), fill = TRUE)
	out[locus_end < locus_start, c("locus_start", "locus_end") := .(locus_end, locus_start)]
	out[, locus_len := locus_end - locus_start + 1L]
	out
}

plot_selected_panel <- function(d, main_title, show_legend = FALSE) {
	tab <- dcast(d, anc ~ locus, value.var = "locus_burden_percent", fill = 0)
	mat <- as.matrix(tab[, setdiff(names(tab), "anc"), with = FALSE])
	rownames(mat) <- tab$anc
	ymax <- max(mat, na.rm = TRUE) * 1.2
	if (!is.finite(ymax) || ymax <= 0) ymax <- 1
	barplot(mat, beside = TRUE, las = 1, cex.names = 0.82,
		ylab = "Selected-locus burden (% of copy-number-adjusted locus)",
		xlab = "Selected locus",
		main = main_title,
		ylim = c(0, ymax),
		legend.text = if (show_legend) rownames(mat) else NULL,
		args.legend = list(x = "topright", cex = 0.68, title = "Archaic reference"))
}

make_group_percent_table <- function(seg, group_col) {
	if (!(group_col %in% names(seg))) return(data.table())
	d <- seg[!is.na(get(group_col)) & get(group_col) != ""]
	if (nrow(d) == 0) return(data.table())
	out <- d[, .(
		n_samples_total = uniqueN(ID),
		n_samples_with_segment = uniqueN(ID[n_segments > 0]),
		mean_total_bp_per_sample = mean(total_bp, na.rm = TRUE),
		mean_segments_per_sample = as.numeric(mean(n_segments, na.rm = TRUE)),
		median_segments_per_sample = as.numeric(median(n_segments, na.rm = TRUE)),
		mean_total_mb_per_sample = mean(total_mb, na.rm = TRUE),
		median_total_mb_per_sample = median(total_mb, na.rm = TRUE),
		mean_denominator_gb_haploid = mean(denominator_gb_haploid, na.rm = TRUE),
		mean_denominator_gb_diploid = mean(denominator_gb_diploid, na.rm = TRUE),
		segment_percent_haploid = mean(segment_percent_haploid, na.rm = TRUE),
		segment_percent_diploid = mean(segment_percent_diploid, na.rm = TRUE),
		segment_percent = mean(segment_percent, na.rm = TRUE)
	), by = c(group_col, "anc")]
	out[, sample_prevalence_percent := 100 * n_samples_with_segment / n_samples_total]
	setorderv(out, c(group_col, "anc"))
	out[]
}

draw_group_percent_plot <- function(gp, group_col, file, title, subtitle = "", percent_col = "segment_percent_diploid", ylab = "% of diploid genome") {
	if (nrow(gp) == 0 || !(group_col %in% names(gp))) return(invisible(FALSE))
	d <- copy(gp)
	d[, group_value := as.character(get(group_col))]
	tab <- dcast(d, group_value ~ anc, value.var = percent_col, fill = 0)
	anc_cols <- setdiff(names(tab), "group_value")
	anc_cols <- c(setdiff(anc_cols, "ANY_REF"), intersect("ANY_REF", anc_cols))
	if (!length(anc_cols)) return(invisible(FALSE))
	mat <- as.matrix(tab[, anc_cols, with = FALSE])
	rownames(mat) <- tab$group_value
	plot_ylim <- c(0, max(mat, na.rm = TRUE) * 1.20)
	if (!all(is.finite(plot_ylim)) || plot_ylim[2] <= 0) plot_ylim <- c(0, 1)
	png(file, width = if (nrow(mat) > 12) 2200 else 1700, height = 1350, units = "px", res = 220)
	op <- par(mar = c(if (nrow(mat) > 12) 8.5 else 5.2, 5, 1.0, 1.2), oma = c(0, 0, if (nzchar(subtitle)) 3.6 else 2.7, 0))
	if (nrow(mat) > 12) {
		barplot(t(mat), beside = TRUE, las = 2, cex.names = 0.55, ylab = ylab, xlab = group_col, ylim = plot_ylim,
			legend.text = anc_cols, args.legend = list(x = "topright", bty = "n", cex = 0.62, title = "Reference"))
	} else {
		barplot(mat, beside = TRUE, las = 1, ylab = ylab, xlab = "Archaic reference", ylim = plot_ylim,
			legend.text = rownames(mat), args.legend = list(x = "top", horiz = TRUE, bty = "n", cex = 0.78, title = group_col))
	}
	mtext(title, outer = TRUE, line = if (nzchar(subtitle)) 2.4 else 1.4, font = 2)
	if (nzchar(subtitle)) mtext(subtitle, outer = TRUE, line = 1.1, cex = 0.72)
	par(op)
	dev.off()
	invisible(TRUE)
}

draw_denom_sensitivity_plot <- function(gp, group_col, file, title) {
	if (nrow(gp) == 0 || !(group_col %in% names(gp))) return(invisible(FALSE))
	d <- copy(gp)
	d[, group_value := as.character(get(group_col))]
	if (uniqueN(d$group_value) > 12) d <- d[anc == "ANY_REF"]
	png(file, width = 2400, height = 1250, units = "px", res = 220)
	op <- par(mfrow = c(1, 2), mar = c(5.2, 5, 1.0, 1.2), oma = c(0, 0, 3, 0))
	for (pc in c("segment_percent_diploid", "segment_percent_haploid")) {
		tab <- dcast(d, group_value ~ anc, value.var = pc, fill = 0)
		anc_cols <- setdiff(names(tab), "group_value")
		anc_cols <- c(setdiff(anc_cols, "ANY_REF"), intersect("ANY_REF", anc_cols))
		mat <- as.matrix(tab[, anc_cols, with = FALSE])
		rownames(mat) <- tab$group_value
		ymax <- max(mat, na.rm = TRUE) * 1.20
		if (!is.finite(ymax) || ymax <= 0) ymax <- 1
		barplot(mat, beside = TRUE, las = 1, ylab = ifelse(pc == "segment_percent_diploid", "% of diploid genome", "% of haploid genome"),
			xlab = "Archaic reference", ylim = c(0, ymax),
			legend.text = rownames(mat), args.legend = list(x = "top", horiz = TRUE, bty = "n", cex = 0.68))
		mtext(ifelse(pc == "segment_percent_diploid", "Diploid denominator", "Haploid denominator"), side = 3, line = 0.1, font = 2)
	}
	mtext(title, outer = TRUE, line = 1.6, font = 2)
	par(op)
	dev.off()
	invisible(TRUE)
}



ref_group_label <- function(x) {
	fifelse(x %in% c("Altai", "Chagyr", "Vindija"), "Neanderthal_refs",
		fifelse(x %in% c("Denisova", "Denisova25"), "Denisova_refs",
			fifelse(x == "ANY_REF", "ANY_REF", "Other")))
}

selected_overlap_diagnostics <- function(seg, loci, ph, stats_by_group, cover_thresholds = c(0, 0.01, 0.10, 0.50, 0.80)) {
	if (nrow(loci) == 0) return(invisible(NULL))
	if (nrow(seg) == 0) return(invisible(NULL))
	seg0 <- copy(seg)
	seg0[, chrom := as.character(chrom)]
	loc0 <- copy(loci)
	loc0[, chrom := as.character(chrom)]
	x <- rbindlist(lapply(seq_len(nrow(loc0)), function(i) {
		l <- loc0[i]
		z <- seg0[chrom == l$chrom & end > l$locus_start & start < l$locus_end]
		if (nrow(z) == 0) return(data.table())
		z[, `:=`(
			locus_id = l$locus_id,
			locus = l$locus,
			locus_start = l$locus_start,
			locus_end = l$locus_end,
			locus_len = l$locus_len
		)]
		z
	}), fill = TRUE)
	if (nrow(x) > 0) {
		x[, overlap_bp := as.numeric(pmax(0L, pmin(end, locus_end) - pmax(start, locus_start)))]
		x <- x[overlap_bp > 0]
		x[, cover_prop := overlap_bp / locus_len]
		x[, ref_group := ref_group_label(anc)]
	}
	if (nrow(x) == 0) {
		x <- data.table(
			ID = character(), chrom = character(), start = integer(), end = integer(), length = integer(),
			slod = numeric(), sites = numeric(), positive_lods = numeric(), negative_lods = numeric(),
			sample_set = character(), super_pop = character(), anc = character(), locus_id = character(),
			locus = character(), locus_start = integer(), locus_end = integer(), locus_len = integer(),
			overlap_bp = numeric(), cover_prop = numeric(), ref_group = character()
		)
	}
	fwrite(x, out_file("selected_all_overlaps.no_cover_filter.tsv"), sep = "\t")

	threshold_dt <- data.table(cover_threshold = cover_thresholds)
	threshold_dt[, threshold_label := paste0("cover", sprintf("%02d", round(100 * cover_threshold)))]
	if (nrow(x) > 0) {
		xthr <- rbindlist(lapply(seq_len(nrow(threshold_dt)), function(i) {
			th <- threshold_dt[i]
			z <- x[cover_prop >= th$cover_threshold]
			if (nrow(z) == 0) return(data.table())
			z[, `:=`(cover_threshold = th$cover_threshold, threshold_label = th$threshold_label)]
			z
		}), fill = TRUE)
	} else {
		xthr <- data.table()
	}
	n_total <- uniqueN(ph$ID)
	locus_ref_grid <- merge(
		merge(copy(loc0)[, tmp := 1], unique(threshold_dt[, .(cover_threshold, threshold_label)])[, tmp := 1], by = "tmp", allow.cartesian = TRUE),
		data.table(ref_group = c("ANY_REF", "Neanderthal_refs", "Denisova_refs", "Other"), tmp = 1),
		by = "tmp", allow.cartesian = TRUE
	)[, tmp := NULL]
	if (nrow(xthr) > 0) {
		sum_ref <- xthr[, .(
			n_segments = as.integer(.N),
			n_samples = as.integer(uniqueN(ID)),
			median_overlap_bp = as.numeric(median(overlap_bp, na.rm = TRUE)),
			max_overlap_bp = as.numeric(max(overlap_bp, na.rm = TRUE)),
			median_cover_prop = as.numeric(median(cover_prop, na.rm = TRUE)),
			max_cover_prop = as.numeric(max(cover_prop, na.rm = TRUE)),
			median_slod = as.numeric(median(slod, na.rm = TRUE)),
			median_sites = as.numeric(median(sites, na.rm = TRUE))
		), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, cover_threshold, threshold_label, ref_group)]
	} else {
		sum_ref <- locus_ref_grid[0, .(locus_id, locus, chrom, locus_start, locus_end, locus_len, cover_threshold, threshold_label, ref_group,
			n_segments = integer(), n_samples = integer(), median_overlap_bp = numeric(), max_overlap_bp = numeric(),
			median_cover_prop = numeric(), max_cover_prop = numeric(), median_slod = numeric(), median_sites = numeric())]
	}
	sum_ref <- merge(locus_ref_grid, sum_ref,
		by = c("locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "cover_threshold", "threshold_label", "ref_group"), all.x = TRUE)
	for (cc in c("n_segments", "n_samples")) sum_ref[is.na(get(cc)), (cc) := 0L]
	sum_ref[, n_samples_total := n_total]
	sum_ref[, carrier_percent := 100 * n_samples / pmax(1, n_samples_total)]
	setorder(sum_ref, locus_id, cover_threshold, ref_group)
	fwrite(sum_ref, out_file("selected_overlap_by_locus_ref_group.csv"))

	by_anc <- data.table()
	if (nrow(xthr) > 0) {
		by_anc <- xthr[, .(
			n_segments = as.integer(.N),
			n_samples = as.integer(uniqueN(ID)),
			median_cover_prop = as.numeric(median(cover_prop, na.rm = TRUE)),
			max_cover_prop = as.numeric(max(cover_prop, na.rm = TRUE)),
			median_slod = as.numeric(median(slod, na.rm = TRUE))
		), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, cover_threshold, threshold_label, anc, ref_group)]
		by_anc[, n_samples_total := n_total]
		by_anc[, carrier_percent := 100 * n_samples / pmax(1, n_samples_total)]
		setorder(by_anc, locus_id, cover_threshold, ref_group, anc)
	}
	if (ncol(by_anc) == 0) {
		by_anc <- data.table(locus_id = character(), locus = character(), chrom = character(), locus_start = integer(),
			locus_end = integer(), locus_len = integer(), cover_threshold = numeric(), threshold_label = character(),
			anc = character(), ref_group = character(), n_segments = integer(), n_samples = integer(),
			median_cover_prop = numeric(), max_cover_prop = numeric(), median_slod = numeric(),
			n_samples_total = integer(), carrier_percent = numeric())
	}
	fwrite(by_anc, out_file("selected_overlap_by_locus_anc.csv"))

	by_group <- data.table()
	group_cols <- intersect(stats_by_group, names(ph))
	if (nrow(xthr) > 0 && length(group_cols) > 0) {
		group_ph <- unique(ph[, c("ID", group_cols), with = FALSE])
		xthr_group <- copy(xthr)
		dup_group_cols <- intersect(group_cols, setdiff(names(xthr_group), "ID"))
		if (length(dup_group_cols) > 0) xthr_group[, (dup_group_cols) := NULL]
		mg <- merge(xthr_group, group_ph, by = "ID", all.x = FALSE, sort = FALSE)
		if (nrow(mg) > 0) {
			den.g <- group_ph[, .(n_samples_total = uniqueN(ID)), by = group_cols]
			by_group <- mg[, .(
				n_segments = as.integer(.N),
				n_samples = as.integer(uniqueN(ID)),
				median_cover_prop = as.numeric(median(cover_prop, na.rm = TRUE)),
				max_cover_prop = as.numeric(max(cover_prop, na.rm = TRUE)),
				median_slod = as.numeric(median(slod, na.rm = TRUE))
			), by = c(group_cols, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "cover_threshold", "threshold_label", "ref_group")]
			by_group <- merge(by_group, den.g, by = group_cols, all.x = TRUE)
			by_group[, carrier_percent := 100 * n_samples / pmax(1, n_samples_total)]
			setorderv(by_group, c(group_cols, "locus_id", "cover_threshold", "ref_group"))
		}
	}
	if (ncol(by_group) == 0) {
		by_group <- data.table(locus_id = character(), locus = character(), chrom = character(), locus_start = integer(),
			locus_end = integer(), locus_len = integer(), cover_threshold = numeric(), threshold_label = character(),
			ref_group = character(), n_segments = integer(), n_samples = integer(), median_cover_prop = numeric(),
			max_cover_prop = numeric(), median_slod = numeric(), n_samples_total = integer(), carrier_percent = numeric())
	}
	fwrite(by_group, out_file("selected_overlap_by_group_ref_group.csv"))

	if (nrow(sum_ref) > 0) {
		n_loci <- uniqueN(sum_ref$locus_id)
		png(out_file("selected_overlap_thresholds.png"), width = 1900, height = max(1100, 520 * n_loci), units = "px", res = 220)
		op <- par(mfrow = c(n_loci, 1), mar = c(5, 5, 3.2, 1.2), oma = c(0, 0, 3.5, 0))
		for (lid in unique(sum_ref$locus_id)) {
			d <- sum_ref[locus_id == lid]
			d <- d[ref_group %in% c("ANY_REF", "Neanderthal_refs", "Denisova_refs")]
			tab <- dcast(d, threshold_label ~ ref_group, value.var = "carrier_percent", fill = 0)
			ref_cols <- intersect(c("ANY_REF", "Neanderthal_refs", "Denisova_refs"), names(tab))
			mat <- t(as.matrix(tab[, ref_cols, with = FALSE]))
			colnames(mat) <- tab$threshold_label
			ymax <- max(mat, na.rm = TRUE) * 1.20
			if (!is.finite(ymax) || ymax <= 0) ymax <- 1
			barplot(mat, beside = TRUE, las = 1, ylim = c(0, ymax), ylab = "Carrier samples (%)", xlab = "Minimum locus coverage",
				legend.text = rownames(mat), args.legend = list(x = "topright", bty = "n", cex = 0.75))
			mtext(unique(d$locus)[1], side = 3, line = 0.5, font = 2, cex = 0.85)
		}
		mtext("Selected-locus overlap diagnostic before locus_min_cover filtering", outer = TRUE, line = 2.1, font = 2)
		par(op)
		dev.off()
	}

	diag_sheets <- list(
		all_overlaps_no_cover_filter = x,
		by_locus_ref_group = sum_ref
	)
	if (nrow(by_anc) > 0) diag_sheets$by_locus_archaic_reference <- by_anc
	if (nrow(by_group) > 0) diag_sheets$by_group_ref_group <- by_group
	write_xlsx(diag_sheets, out_file("selected_overlap_diagnostics.xlsx"))
	invisible(list(all_overlaps = x, by_locus_ref_group = sum_ref, by_locus_anc = by_anc, by_group = by_group))
}

args <- commandArgs(trailingOnly = TRUE)

deprecated_args <- intersect(args, c("--summary_by", "--stats_by"))
if (length(deprecated_args) > 0) {
	stop("Deprecated option(s) not supported: ", paste(deprecated_args, collapse = ", "),
		". Please use only --stats_by_group and --selected_loci.")
}
valid_args <- c("--in", "--outdir", "--sample_file", "--sample_keep", "--sample_id_col", "--stats_by_group",
    "--selected_loci", "--window_size", "--locus_min_cover", "--top_n_windows", "--plot_style", "--admixed_african_pops", "--genome_build")
unknown_args <- args[grepl("^--", args) & !(args %in% valid_args)]
if (length(unknown_args) > 0) stop("Unknown option(s): ", paste(unknown_args, collapse = ", "))

infile <- get_arg("--in")
outdir <- get_arg("--outdir", ".")
sample_file <- get_arg("--sample_file", NA_character_)
sample_keep_file <- get_arg("--sample_keep", NA_character_)
sample_id_col <- get_arg("--sample_id_col", "sample")
stats_by_group_arg <- get_arg("--stats_by_group", NA_character_)
selected_loci_arg <- get_arg("--selected_loci", NA_character_)
window_size <- as.integer(get_arg("--window_size", "1000000"))
locus_min_cover <- as.numeric(get_arg("--locus_min_cover", "0.5"))
top_n_windows <- as.integer(get_arg("--top_n_windows", "100"))
plot_style <- tolower(get_arg("--plot_style", "both"))
admixed_african_pops_arg <- get_arg("--admixed_african_pops", "ACB,ASW")
genome_build <- tolower(get_arg("--genome_build", "b37"))
if (!(genome_build %in% c("b37", "b38"))) stop("--genome_build must be b37 or b38")

# Internal segment/locus coordinates are treated as BED-like half-open intervals.
# Convert the 1-based inclusive PAR definitions used by bcftools into half-open
# coordinates so chrX locus denominators can account for male copy number.
x_par_half_open <- if (genome_build == "b38") {
    matrix(c(10000, 2781479, 155701382, 156030895), ncol = 2, byrow = TRUE)
} else {
    matrix(c(60000, 2699520, 154931043, 155260560), ncol = 2, byrow = TRUE)
}
interval_overlap_bp <- function(start, end, left, right) {
    pmax(0, pmin(as.numeric(end), as.numeric(right)) - pmax(as.numeric(start), as.numeric(left)))
}
locus_denominator_bp <- function(chrom, start, end, n_female, n_male, n_unknown = 0) {
    chrom <- norm_chrom(chrom)
    len <- pmax(0, as.numeric(end) - as.numeric(start))
    n_female <- as.numeric(n_female); n_male <- as.numeric(n_male); n_unknown <- as.numeric(n_unknown)
    if (!(chrom %in% c("X", "23"))) return((n_female + n_male + n_unknown) * 2 * len)
    if (n_unknown > 0) stop("chrX selected-locus denominator requires known sex for every sample")
    par_bp <- sum(vapply(seq_len(nrow(x_par_half_open)), function(i) {
        interval_overlap_bp(start, end, x_par_half_open[i, 1], x_par_half_open[i, 2])
    }, numeric(1)))
    par_bp <- pmin(len, par_bp)
    nonpar_bp <- pmax(0, len - par_bp)
    n_female * 2 * len + n_male * (2 * par_bp + nonpar_bp)
}
admixed_african_pops <- trimws(unlist(strsplit(admixed_african_pops_arg, ",", fixed = TRUE)))
admixed_african_pops <- admixed_african_pops[nzchar(admixed_african_pops)]

if (is.na(infile)) stop("Please provide --in input.segments.tsv.gz")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (is.na(window_size) || window_size <= 0) stop("Bad --window_size: ", window_size)
if (is.na(locus_min_cover) || locus_min_cover <= 0 || locus_min_cover > 1) stop("Bad --locus_min_cover: ", locus_min_cover)
if (!(plot_style %in% c("both", "circular", "linear"))) stop("Bad --plot_style: ", plot_style, "; use both, circular, or linear")

stats_by_group <- character(0)
if (!is.na(stats_by_group_arg) && nzchar(stats_by_group_arg)) {
	stats_by_group <- trimws(unlist(strsplit(stats_by_group_arg, ",", fixed = TRUE)))
	stats_by_group <- stats_by_group[nzchar(stats_by_group)]
}

in0 <- find_input(infile, outdir)
if (is.na(in0)) stop("Input file not found: ", infile)
infile <- in0

input_prefix <- basename(infile)
input_prefix <- sub("\\.segments\\.tsv\\.gz$", "", input_prefix)
input_prefix <- sub("\\.tsv\\.gz$", "", input_prefix)
input_prefix <- sub("\\.gz$", "", input_prefix)
unlink(file.path(outdir, c(
	"check.tsv",
	"check.csv",
	"burden.tsv.gz",
	"burden.csv.gz",
	"burden.csv",
	"burden_by_group.tsv.gz",
	"burden_by_group.csv.gz",
	"burden_by_group.csv",
	"top_windows.tsv",
	"top_windows.csv",
	"top.csv",
	"burden.png",
	"all_gw.csv",
	"all_gw_by_group.csv",
	"all_gw.png",
	"all_gw_linear.png",
	"all_gw_circular.png",
	"all_gw_peaks.png",
	"all.xlsx",
	"group_burden.tsv",
	"group_burden.csv",
	"group_burden.png",
	"all_percent.csv",
	"all_percent.per_reference.csv",
	"all_percent_by_pop.csv",
	"all_percent_by_pop.png",
	"all_percent_afr_sensitivity.csv",
	"all_percent_afr_sensitivity.png",
	"all_percent_denom_sensitivity.png",
	"all_percent.png",
	"overlap_reduction_diagnostics.csv",
	"overlap_reduction_by_super_pop.csv",
	"overlap_reduction.png",
	"cross_reference_overlap_diagnostics.csv",
	"cross_reference_overlap_by_super_pop.csv",
	"genomewide_segments.any_ref.reduced.tsv.gz",
	"ibdmix_diagnostics.tsv",
	"sample_burden.tsv",
	"sample_burden.csv",
	"sample_boxplot.png",
	"all_length.csv",
	"all_length.png",
	"genomewide_segments.raw.tsv.gz",
	"genomewide_segments.reduced.tsv.gz",
	"genomewide_segments.tsv.gz",
	"locus_matches.tsv",
	"locus_matches.csv",
	"locus_summary.tsv",
	"locus_summary.csv",
	"locus_summary.png",
	"selected_overall.csv",
	"selected_overall.png",
	"selected.png",
	"selected.xlsx",
	"locus_group.tsv",
	"locus_group.csv",
	"locus_group.png",
	"selected_percent.csv",
	"selected_by_group.csv",
	"selected_percent.png",
	"input.check.tsv",
	"group.by_anc.segment_percent.tsv",
	"sample.by_anc.segment_burden.tsv",
	"group.by_anc.segment_percent.png",
	"sample.by_anc.segment_mb.boxplot.png",
	"locus.sample_matches.tsv",
	"locus.by_anc.summary.tsv",
	"locus.by_anc.summary.png",
	"locus.by_super_pop_gender.summary.tsv",
	"locus.by_group_anc.summary.png",
	"locus.by_group_anc.summary.tsv",
	"genomewide.window1000000.burden.tsv.gz",
	"genomewide.window1000000.by_super_pop_gender.burden.tsv.gz",
	"genomewide.window1000000.top_windows.tsv",
	"genomewide.window1000000.burden.3x2.png"
)), force = TRUE)
unlink(c(out_dir("sample.by_anc.segment_mb.boxplot"), out_dir("locus.by_group_anc.summary"), out_dir(paste0("bedGraph.window", window_size))), recursive = TRUE, force = TRUE)

message("Input:  ", infile)
message("Outdir: ", outdir)
message("Input prefix ignored for output names: ", input_prefix)
message("Window size: ", window_size)
if (!is.na(sample_file)) message("Sample file: ", sample_file)
if (!is.na(sample_keep_file)) message("Sample keep: ", sample_keep_file)
if (length(stats_by_group) > 0) message("Stats by group: ", paste(stats_by_group, collapse = ","))
message("Admixed African populations handled separately: ", paste(admixed_african_pops, collapse = ","))
if (!is.na(selected_loci_arg) && nzchar(selected_loci_arg)) message("Selected loci: ", selected_loci_arg)
message("Genome-wide PNG plots: disabled; use all_gw.csv or IGV tracks")

n_do <- tryCatch(count_do(infile), error = function(e) NA_integer_)

cmd <- paste("gzip -dc", shQuote(infile), "| awk 'BEGIN{FS=OFS=\"\\t\"} $1!=\"DO\"'")
dat <- fread(cmd = cmd, sep = "\t", header = TRUE, fill = TRUE, showProgress = TRUE)

need <- c("ID", "chrom", "start", "end", "length", "slod", "sites", "positive_lods", "negative_lods", "anc")
miss <- setdiff(need, names(dat))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

if (!("sample_set" %in% names(dat))) dat[, sample_set := "sample_keep"]

dat <- dat[!is.na(ID) & ID != "" & ID != "ID" & ID != "DO"]
dat[, ID := as.character(ID)]
dat[, chrom := as.character(chrom)]
dat[, chrom := sub("^chr", "", chrom, ignore.case = TRUE)]
dat[, start := as.integer(start)]
dat[, end := as.integer(end)]
dat[, seg_len := suppressWarnings(as.integer(length))]
dat[, slod := as.numeric(slod)]
dat[, sites := as.numeric(sites)]
dat[, positive_lods := as.numeric(positive_lods)]
dat[, negative_lods := as.numeric(negative_lods)]
dat[, anc := as.character(anc)]
dat[, sample_set := as.character(sample_set)]
if (!("super_pop" %in% names(dat))) dat[, super_pop := NA_character_]
dat[, super_pop := as.character(super_pop)]
dat <- dat[!is.na(start) & !is.na(end) & !is.na(anc) & anc != ""]
dat <- dat[end > start]
dat[, coord_len := segment_width(start, end)]
dat[is.na(seg_len) | seg_len <= 0 | abs(seg_len - coord_len) > 1L, seg_len := coord_len]
dat <- dat[seg_len > 0]
dat[, coord_len := NULL]
if (nrow(dat) == 0) stop("No usable segment rows after cleanup")

n0 <- nrow(dat)
dat <- unique(dat, by = c("ID", "anc", "chrom", "start", "end"))
n_dup_removed <- n0 - nrow(dat)
segment.detail <- copy(dat)

seg_ids <- sort(unique(dat$ID))
ids <- seg_ids
if (!is.na(sample_keep_file) && file.exists(sample_keep_file)) {
	ids0 <- fread(sample_keep_file, header = FALSE, showProgress = FALSE)[[1]]
	ids0 <- sort(unique(as.character(ids0[!is.na(ids0) & nzchar(ids0)])))
	if (length(ids0) == 0) stop("sample_keep has no sample IDs: ", sample_keep_file)
	dat <- dat[ID %in% ids0]
	ids <- ids0
	seg_ids <- sort(unique(dat$ID))
} else if (!is.na(sample_keep_file)) {
	stop("sample_keep file not found: ", sample_keep_file)
}
if (nrow(dat) == 0) stop("No segment rows remain after applying sample_keep: ", sample_keep_file)
segment.detail.raw <- copy(dat)
raw_rows_after_sample_keep <- nrow(segment.detail.raw)
raw_bp_after_sample_keep <- segment.detail.raw[, sum(as.numeric(seg_len), na.rm = TRUE)]

# Merge/reduce overlapping intervals within each ID + archaic reference + chromosome.
# The per-reference reduced table is written as genomewide_segments.reduced.tsv.gz.
# A second ANY_REF table sets anc=ANY_REF before reduction, which prevents
# interpreting overlapping Altai/Chagyr/Vindija calls as additive genome burden.
dat_ref <- reduce_segments(segment.detail.raw)
segment.detail <- copy(dat_ref)
reduced_rows <- nrow(segment.detail)
reduced_bp <- segment.detail[, sum(as.numeric(seg_len), na.rm = TRUE)]
reduced_bp_removed <- raw_bp_after_sample_keep - reduced_bp

dat_any_ref_raw <- copy(segment.detail.raw)
dat_any_ref_raw[, anc := "ANY_REF"]
dat_any_ref <- reduce_segments(dat_any_ref_raw)
any_ref_rows <- nrow(dat_any_ref)
any_ref_bp <- dat_any_ref[, sum(as.numeric(seg_len), na.rm = TRUE)]
cross_ref_bp_removed <- reduced_bp - any_ref_bp

per_ref_ancs <- sort(unique(segment.detail$anc))
ancs <- c(per_ref_ancs, "ANY_REF")
dat <- rbindlist(list(segment.detail, dat_any_ref), fill = TRUE, use.names = TRUE)
chrom_order <- c(as.character(1:22), "X", "Y", "MT", "M")

message("Rows after cleanup: ", nrow(segment.detail.raw))
message("Rows after within-reference interval reduction: ", nrow(segment.detail))
message("Rows after cross-reference ANY_REF reduction: ", nrow(dat_any_ref))
message("Raw bp removed by within-reference interval reduction: ", format(round(reduced_bp_removed), big.mark = ","))
message("Per-reference bp removed by cross-reference ANY_REF reduction: ", format(round(cross_ref_bp_removed), big.mark = ","))
message("Samples in segment file: ", length(seg_ids))
message("Samples in denominator: ", length(ids))
message("Archaic references: ", paste(per_ref_ancs, collapse = ", "), "; plus ANY_REF for cross-reference-deduplicated burden")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step1: read sample metadata for group summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ph <- NULL
if (!is.na(sample_file) && file.exists(sample_file)) {
	info <- fread(sample_file, header = TRUE, fill = TRUE, check.names = FALSE, showProgress = FALSE)
	setnames(info, trimws(sub("^\ufeff", "", names(info))))
	if (!(sample_id_col %in% names(info))) {
		stop("sample_id_col not found in sample_file: ", sample_id_col,
			"; columns read: ", paste(names(info), collapse = ", "))
	}
	if (length(stats_by_group) > 0) {
		generated_group_cols <- c("afr_sensitivity_group", "afr_core_status")
		miss <- setdiff(stats_by_group, c(names(info), generated_group_cols))
		if (length(miss) > 0) stop("stats_by_group column(s) not found in sample_file or generated diagnostics: ", paste(miss, collapse = ", "))
	}
		sex_col <- intersect(c("sex", "gender", "Sex", "Gender", "SEX", "GENDER"), names(info))[1]
		keep_cols <- unique(c(sample_id_col, stats_by_group, "pop", "super_pop", sex_col))
		keep_cols <- keep_cols[!is.na(keep_cols) & keep_cols %in% names(info)]
		ph <- unique(info[, keep_cols, with = FALSE])
		setnames(ph, sample_id_col, "ID")
		if (!is.na(sex_col)) ph[, .sex_for_denom := as.character(get(sex_col))]
		ph[, ID := as.character(ID)]
		ph <- ph[ID %in% ids]
		if ("pop" %in% names(ph)) ph[, pop := as.character(pop)]
		if ("super_pop" %in% names(ph)) ph[, super_pop := as.character(super_pop)]
		if (all(c("pop", "super_pop") %in% names(ph))) {
			ph[, afr_sensitivity_group := fifelse(super_pop == "AFR" & pop %in% admixed_african_pops,
				paste0("AFR_admixed_", paste(admixed_african_pops, collapse = "_")),
				fifelse(super_pop == "AFR", paste0("AFR_core_no_", paste(admixed_african_pops, collapse = "_")), super_pop))]
			ph[, afr_core_status := fifelse(super_pop == "AFR" & pop %in% admixed_african_pops,
				paste0("AFR_admixed_", paste(admixed_african_pops, collapse = "_")),
				fifelse(super_pop == "AFR", paste0("AFR_core_no_", paste(admixed_african_pops, collapse = "_")), "non_AFR"))]
		}
		generated_missing <- setdiff(intersect(stats_by_group, generated_group_cols), names(ph))
		if (length(generated_missing) > 0) {
			stop("stats_by_group generated column(s) require pop and super_pop in sample_file: ",
				paste(generated_missing, collapse = ", "))
		}
		for (x in stats_by_group) ph[, (x) := as.character(get(x))]
		if (!(".sex_for_denom" %in% names(ph))) ph[, .sex_for_denom := NA_character_]
		ph[, .sex_for_denom := tolower(trimws(as.character(.sex_for_denom)))]
		for (x in stats_by_group) ph <- ph[!is.na(get(x)) & get(x) != ""]
	if (nrow(ph) == 0) stop("No matching sample IDs between --in and --sample_file")
	message("Matched sample_file IDs: ", uniqueN(ph$ID))
} else if (!is.na(sample_file)) {
	message("Sample file not found, skip group-specific outputs: ", sample_file)
	stats_by_group <- character(0)
}

if (is.null(ph)) {
	ph <- data.table(ID = ids)
	ph[, .sex_for_denom := NA_character_]
}
if (length(stats_by_group) > 0) {
	missing_group_cols <- setdiff(stats_by_group, names(ph))
	if (length(missing_group_cols) > 0) {
		message("Stats by group column(s) unavailable after sample metadata load; skipping: ",
			paste(missing_group_cols, collapse = ", "))
		stats_by_group <- intersect(stats_by_group, names(ph))
	}
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step3: compute genome-wide window inherited-frequency table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
windows <- make_windows(dat, window_size)
genome_bp <- windows[, sum(as.numeric(end - start + 1L), na.rm = TRUE)]
chr_lengths_b37 <- c(
    "1"=249250621,"2"=243199373,"3"=198022430,"4"=191154276,"5"=180915260,"6"=171115067,
    "7"=159138663,"8"=146364022,"9"=141213431,"10"=135534747,"11"=135006516,"12"=133851895,
    "13"=115169878,"14"=107349540,"15"=102531392,"16"=90354753,"17"=81195210,"18"=78077248,
    "19"=59128983,"20"=63025520,"21"=48129895,"22"=51304566,"X"=155270560,"23"=155270560,
    "Y"=59373566,"MT"=16569,"M"=16569)
chr_lengths_b38 <- c(
    "1"=248956422,"2"=242193529,"3"=198295559,"4"=190214555,"5"=181538259,"6"=170805979,
    "7"=159345973,"8"=145138636,"9"=138394717,"10"=133797422,"11"=135086622,"12"=133275309,
    "13"=114364328,"14"=107043718,"15"=101991189,"16"=90338345,"17"=83257441,"18"=80373285,
    "19"=58617616,"20"=64444167,"21"=46709983,"22"=50818468,"X"=156040895,"23"=156040895,
    "Y"=57227415,"MT"=16569,"M"=16569)
chr_lengths <- if (genome_build == "b38") chr_lengths_b38 else chr_lengths_b37
denom_chroms <- sort(unique(windows$chrom))
known_denom_chroms <- denom_chroms[denom_chroms %in% names(chr_lengths)]
unknown_denom_chroms <- setdiff(denom_chroms, known_denom_chroms)
haploid_denom_bp <- sum(as.numeric(chr_lengths[known_denom_chroms]), na.rm = TRUE)
if (length(unknown_denom_chroms) > 0) {
	haploid_denom_bp <- haploid_denom_bp + windows[chrom %in% unknown_denom_chroms, sum(as.numeric(end - start + 1L), na.rm = TRUE)]
}
if (haploid_denom_bp <= 0) haploid_denom_bp <- genome_bp
# Segment totals are summed across both haplotypes per individual, so the
# default percentage uses a diploid genome denominator.  The haploid percentage
# column is retained for comparison with older haploid-scaled summaries.
diploid_genome_bp <- 2 * haploid_denom_bp
sample.denom <- unique(ph[, .(ID, .sex_for_denom)])
sample.denom <- merge(data.table(ID = ids), sample.denom, by = "ID", all.x = TRUE)
sample.denom[is.na(.sex_for_denom), .sex_for_denom := ""]
sample.denom[, .sex_for_denom := fifelse(.sex_for_denom %in% c("1","m","male","xy"), "male",
    fifelse(.sex_for_denom %in% c("2","f","female","xx"), "female", .sex_for_denom))]
x_included <- any(denom_chroms %in% c("X", "23"))
x_bp <- if (x_included) as.numeric(chr_lengths["X"]) else 0
sample.denom[, denominator_bp_haploid := haploid_denom_bp]
sample.denom[, denominator_bp_diploid := diploid_genome_bp]
if (x_included) sample.denom[.sex_for_denom == "male", denominator_bp_diploid := diploid_genome_bp - x_bp]
sample.denom[, denominator_gb_haploid := denominator_bp_haploid / 1e9]
sample.denom[, denominator_gb_diploid := denominator_bp_diploid / 1e9]
wins.map <- windows[, .(chrom, win_id, window, win_start = start, win_end = end, chr_rank, genome_pos)]
setkey(wins.map, chrom, win_id)

# Expand each inherited segment to the window(s) it overlaps. With the default
# 1 Mb window, most segments contribute to only one row here.
dseg <- dat[, .(
	seg_id = .I, ID, anc, chrom,
	start_win = ((start - 1L) %/% window_size) + 1L,
	end_win = ((end - 1L) %/% window_size) + 1L
)]
ov <- dseg[, .(win_id = seq.int(start_win[1], end_win[1])), by = .(seg_id, ID, anc, chrom)]
ov[, seg_id := NULL]
ov <- unique(ov)
ov <- wins.map[ov, on = .(chrom, win_id), nomatch = 0L]
ov <- ov[, .(ID, anc, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)]

den_all <- length(ids)
win.burden <- ov[, .(n_samples_with_segment = uniqueN(ID)), by = .(anc, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)]
win.burden[, n_samples_total := den_all]
win.burden[, burden_fraction := n_samples_with_segment / n_samples_total]
win.burden[, burden_percent := burden_fraction * 100]
setcolorder(win.burden, c("anc", "chrom", "win_id", "window", "win_start", "win_end", "n_samples_with_segment", "n_samples_total", "burden_fraction", "burden_percent", "chr_rank", "genome_pos"))
setorder(win.burden, anc, chr_rank, win_start)

win.file <- out_file("all_gw.csv")
fwrite(win.burden, win.file)
fwrite(segment.detail.raw, out_file("genomewide_segments.raw.tsv.gz"), sep = "\t")
fwrite(segment.detail, out_file("genomewide_segments.reduced.tsv.gz"), sep = "\t")
fwrite(dat_any_ref, out_file("genomewide_segments.any_ref.reduced.tsv.gz"), sep = "\t")

bedgraph.dir <- out_dir("bedgraph_tmp")
dir.create(bedgraph.dir, recursive = TRUE, showWarnings = FALSE)
for (a in sort(unique(win.burden$anc))) {
	bg <- win.burden[anc == a, .(
		chrom = ifelse(chrom %in% c("X", "Y", "M", "MT"), paste0("chr", chrom), paste0("chr", chrom)),
		start0 = win_start - 1L,
		end = win_end,
		burden_percent,
		chr_rank
	)]
	setorder(bg, chr_rank, start0)
	bg[, chr_rank := NULL]
	fwrite(bg, file.path(bedgraph.dir, paste0(a, ".window", window_size, ".burden_percent.bedGraph")), sep = "\t", col.names = FALSE)
}

top.file <- out_file("top.csv")
top.windows <- win.burden[, head(.SD[order(-burden_percent, -n_samples_with_segment)], top_n_windows), by = anc]
fwrite(top.windows, top.file)

# Optional group-specific window burden table, when --sample_file and --stats_by_group are provided.
if (length(stats_by_group) > 0) {
	ov.g <- merge(ov[, .(ID, anc, chrom, win_id, window, win_start, win_end, chr_rank, genome_pos)], ph, by = "ID", all.x = FALSE)
	if (nrow(ov.g) > 0) {
		den.group <- ph[, .(n_samples_total = uniqueN(ID)), by = stats_by_group]
		win.g <- ov.g[, .(n_samples_with_segment = uniqueN(ID)), by = c(stats_by_group, "anc", "chrom", "win_id", "window", "win_start", "win_end", "chr_rank", "genome_pos")]
		win.g <- merge(win.g, den.group, by = stats_by_group, all.x = TRUE)
		win.g[, burden_fraction := n_samples_with_segment / n_samples_total]
		win.g[, burden_percent := burden_fraction * 100]
		setorderv(win.g, c(stats_by_group, "anc", "chr_rank", "win_start"))
		group_tag <- paste(stats_by_group, collapse = "_")
		fwrite(win.g, out_file("all_gw_by_group.csv"))
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step4: compute sample and group inherited-frequency summaries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seg.sample <- dat[, .(
	n_segments = .N,
	n_raw_segments = sum_num(n_raw_segments),
	total_bp = sum_num(seg_len),
	raw_bp_before_reduce = sum_num(raw_bp_before_reduce),
	bp_removed_by_reduce = sum_num(bp_removed_by_reduce),
	total_mb = sum_num(seg_len) / 1e6,
	mean_length = mean(seg_len, na.rm = TRUE),
	mean_slod = mean(slod, na.rm = TRUE),
	median_slod = median(slod, na.rm = TRUE)
), by = .(ID, anc)]

all.sample.anc <- CJ(ID = ids, anc = ancs, unique = TRUE)
seg.sample <- merge(all.sample.anc, seg.sample, by = c("ID", "anc"), all.x = TRUE)
for (cc in c("n_segments", "n_raw_segments", "total_bp", "raw_bp_before_reduce", "bp_removed_by_reduce", "total_mb")) seg.sample[is.na(get(cc)), (cc) := 0]
seg.sample <- merge(seg.sample, ph, by = "ID", all.x = TRUE)
seg.sample <- merge(seg.sample, sample.denom[, .(ID, denominator_bp_haploid, denominator_bp_diploid, denominator_gb_haploid, denominator_gb_diploid)], by = "ID", all.x = TRUE)
seg.sample[, segment_percent_haploid := 100 * total_bp / denominator_bp_haploid]
seg.sample[, segment_percent_diploid := 100 * total_bp / denominator_bp_diploid]
seg.sample[, segment_percent := segment_percent_diploid]
fwrite(seg.sample, out_file("all_length.csv"))

if (length(stats_by_group) > 0) {
	group.cols <- stats_by_group
	seg.sample.group <- copy(seg.sample)
	for (g in group.cols) seg.sample.group <- seg.sample.group[!is.na(get(g)) & get(g) != ""]
	if (nrow(seg.sample.group) == 0) stop("No samples with complete stats_by_group columns: ", paste(group.cols, collapse = ","))
	group.percent <- seg.sample.group[, .(
		n_samples_total = uniqueN(ID),
		n_samples_with_segment = uniqueN(ID[n_segments > 0]),
		mean_total_bp_per_sample = mean(total_bp, na.rm = TRUE),
		mean_segments_per_sample = as.numeric(mean(n_segments, na.rm = TRUE)),
		median_segments_per_sample = as.numeric(median(n_segments, na.rm = TRUE)),
		mean_total_mb_per_sample = mean(total_mb, na.rm = TRUE),
		median_total_mb_per_sample = median(total_mb, na.rm = TRUE),
		mean_denominator_gb_haploid = mean(denominator_gb_haploid, na.rm = TRUE),
		mean_denominator_gb_diploid = mean(denominator_gb_diploid, na.rm = TRUE),
		segment_percent_haploid = mean(segment_percent_haploid, na.rm = TRUE),
		segment_percent_diploid = mean(segment_percent_diploid, na.rm = TRUE),
		segment_percent = mean(segment_percent, na.rm = TRUE)
	), by = c(group.cols, "anc")]
	group.percent[, sample_prevalence_percent := 100 * n_samples_with_segment / n_samples_total]
	setorderv(group.percent, c(group.cols, "anc"))
	fwrite(group.percent, out_file("all_percent.csv"))
	fwrite(group.percent[anc != "ANY_REF"], out_file("all_percent.per_reference.csv"))

	png(out_file("all_percent.png"), width = if (length(group.cols) == 1) 1700 else 1700 * length(group.cols), height = 1350, units = "px", res = 220)
	op <- par(mfrow = c(1, length(group.cols)), mar = c(5, 5, 0.8, 1.2), oma = c(0, 0, 3.7, 0))
	for (g in group.cols) {
		dg <- make_group_percent_table(seg.sample, g)
		dg[, group_value := as.character(get(g))]
		tab <- dcast(dg, group_value ~ anc, value.var = "segment_percent_diploid", fill = 0)
		anc_cols <- setdiff(names(tab), "group_value")
		anc_cols <- c(setdiff(anc_cols, "ANY_REF"), intersect("ANY_REF", anc_cols))
		mat <- as.matrix(tab[, anc_cols, with = FALSE])
		rownames(mat) <- tab$group_value
		group_colors <- rainbow(nrow(mat))
		plot_ylim <- c(0, max(mat, na.rm = TRUE) * 1.18)
		if (!all(is.finite(plot_ylim)) || plot_ylim[2] <= 0) plot_ylim <- c(0, 1)
		barplot(mat, beside = TRUE, las = 1, ylab = "% of diploid genome", col = group_colors,
			xlab = "Archaic reference", ylim = plot_ylim,
			legend.text = rownames(mat), args.legend = list(x = "top", horiz = TRUE, bty = "n", cex = 0.78, title = g))
	}
	mtext("Genome-wide archaic segment burden", outer = TRUE, line = 2.45, font = 2)
	mtext(paste0("Diploid denominator = ", round(diploid_genome_bp / 1e9, 2),
		" Gb; ANY_REF is cross-reference-deduplicated and should be used for total burden"), outer = TRUE, line = 1.15, cex = 0.72)
	par(op)
	dev.off()

	first_group <- group.cols[1]
	draw_denom_sensitivity_plot(make_group_percent_table(seg.sample, first_group), first_group, out_file("all_percent_denom_sensitivity.png"),
		"Denominator sensitivity: diploid versus haploid scaling")
}

if ("pop" %in% names(seg.sample)) {
	percent.by_pop <- make_group_percent_table(seg.sample, "pop")
	if (nrow(percent.by_pop) > 0) {
		fwrite(percent.by_pop, out_file("all_percent_by_pop.csv"))
		draw_group_percent_plot(percent.by_pop, "pop", out_file("all_percent_by_pop.png"),
			"Genome-wide archaic segment burden by 1000G population",
			"Use this plot to check whether AFR signal is driven by ACB/ASW or other populations")
	}
}

if ("afr_sensitivity_group" %in% names(seg.sample)) {
	percent.afr <- make_group_percent_table(seg.sample, "afr_sensitivity_group")
	if (nrow(percent.afr) > 0) {
		fwrite(percent.afr, out_file("all_percent_afr_sensitivity.csv"))
		draw_group_percent_plot(percent.afr, "afr_sensitivity_group", out_file("all_percent_afr_sensitivity.png"),
			"AFR sensitivity: core AFR versus admixed ACB/ASW",
			paste0("Admixed African populations separated: ", paste(admixed_african_pops, collapse = ", ")))
	}
}

# Diagnostics for within-reference overlap reduction and cross-reference double counting.
overlap.raw <- segment.detail.raw[, .(raw_n_segments = .N, raw_total_bp = sum_num(seg_len)), by = .(ID, anc)]
overlap.red <- segment.detail[, .(reduced_n_segments = .N, reduced_total_bp = sum_num(seg_len), bp_removed_by_reduce = sum_num(bp_removed_by_reduce)), by = .(ID, anc)]
overlap.sample <- merge(CJ(ID = ids, anc = per_ref_ancs, unique = TRUE), overlap.raw, by = c("ID", "anc"), all.x = TRUE)
overlap.sample <- merge(overlap.sample, overlap.red, by = c("ID", "anc"), all.x = TRUE)
for (cc in c("raw_n_segments", "raw_total_bp", "reduced_n_segments", "reduced_total_bp", "bp_removed_by_reduce")) overlap.sample[is.na(get(cc)), (cc) := 0]
overlap.sample[, bp_removed_fraction := fifelse(raw_total_bp > 0, bp_removed_by_reduce / raw_total_bp, 0)]
overlap.summary <- overlap.sample[, .(
	n_samples = uniqueN(ID), raw_n_segments = sum_num(raw_n_segments), reduced_n_segments = sum_num(reduced_n_segments),
	raw_total_bp = sum_num(raw_total_bp), reduced_total_bp = sum_num(reduced_total_bp), bp_removed_by_reduce = sum_num(bp_removed_by_reduce),
	bp_removed_percent = 100 * sum_num(bp_removed_by_reduce) / pmax(1, sum_num(raw_total_bp))
), by = anc]
setorder(overlap.summary, anc)
fwrite(overlap.summary, out_file("overlap_reduction_diagnostics.csv"))

if ("super_pop" %in% names(ph)) {
	overlap.by_super <- merge(overlap.sample, unique(ph[, .(ID, super_pop)]), by = "ID", all.x = TRUE)
	overlap.by_super <- overlap.by_super[!is.na(super_pop) & super_pop != "", .(
		n_samples = uniqueN(ID), raw_n_segments = sum_num(raw_n_segments), reduced_n_segments = sum_num(reduced_n_segments),
		raw_total_bp = sum_num(raw_total_bp), reduced_total_bp = sum_num(reduced_total_bp), bp_removed_by_reduce = sum_num(bp_removed_by_reduce),
		bp_removed_percent = 100 * sum_num(bp_removed_by_reduce) / pmax(1, sum_num(raw_total_bp))
	), by = .(super_pop, anc)]
	setorder(overlap.by_super, super_pop, anc)
	fwrite(overlap.by_super, out_file("overlap_reduction_by_super_pop.csv"))
}

png(out_file("overlap_reduction.png"), width = 1700, height = 1200, units = "px", res = 220)
op <- par(mar = c(5, 5, 1, 1), oma = c(0, 0, 2.5, 0))
barplot(setNames(overlap.summary$bp_removed_percent, overlap.summary$anc), las = 1, ylab = "% raw bp removed", xlab = "Archaic reference")
mtext("Within-reference overlap reduction diagnostic", outer = TRUE, line = 1.4, font = 2)
par(op)
dev.off()

perref.sample.total <- segment.detail[, .(per_ref_sum_bp = sum_num(seg_len), per_ref_n_segments = .N), by = ID]
anyref.sample.total <- dat_any_ref[, .(any_ref_bp = sum_num(seg_len), any_ref_n_segments = .N), by = ID]
crossref.sample <- merge(data.table(ID = ids), perref.sample.total, by = "ID", all.x = TRUE)
crossref.sample <- merge(crossref.sample, anyref.sample.total, by = "ID", all.x = TRUE)
for (cc in c("per_ref_sum_bp", "per_ref_n_segments", "any_ref_bp", "any_ref_n_segments")) crossref.sample[is.na(get(cc)), (cc) := 0]
crossref.sample[, cross_ref_bp_removed := per_ref_sum_bp - any_ref_bp]
crossref.sample[, cross_ref_removed_percent := 100 * cross_ref_bp_removed / pmax(1, per_ref_sum_bp)]
fwrite(crossref.sample, out_file("cross_reference_overlap_diagnostics.csv"))
if ("super_pop" %in% names(ph)) {
	crossref.by_super <- merge(crossref.sample, unique(ph[, .(ID, super_pop)]), by = "ID", all.x = TRUE)
	crossref.by_super <- crossref.by_super[!is.na(super_pop) & super_pop != "", .(
		n_samples = uniqueN(ID), per_ref_sum_bp = sum_num(per_ref_sum_bp), any_ref_bp = sum_num(any_ref_bp),
		cross_ref_bp_removed = sum_num(cross_ref_bp_removed),
		cross_ref_removed_percent = 100 * sum_num(cross_ref_bp_removed) / pmax(1, sum_num(per_ref_sum_bp))
	), by = super_pop]
	setorder(crossref.by_super, super_pop)
	fwrite(crossref.by_super, out_file("cross_reference_overlap_by_super_pop.csv"))
}

diag <- data.table(
	item = c("raw_rows_after_sample_keep", "raw_bp_after_sample_keep", "within_ref_reduced_rows", "within_ref_reduced_bp",
		"within_ref_bp_removed", "any_ref_reduced_rows", "any_ref_reduced_bp", "cross_ref_bp_removed",
		"n_samples_denominator", "haploid_denominator_gb", "diploid_denominator_gb", "admixed_african_pops",
		"african_admixed_sample_n", "african_core_sample_n"),
	value = as.character(c(raw_rows_after_sample_keep, round(raw_bp_after_sample_keep), reduced_rows, round(reduced_bp),
		round(reduced_bp_removed), any_ref_rows, round(any_ref_bp), round(cross_ref_bp_removed),
		length(ids), round(haploid_denom_bp / 1e9, 4), round(diploid_genome_bp / 1e9, 4), paste(admixed_african_pops, collapse = ","),
		if (all(c("pop", "super_pop") %in% names(ph))) uniqueN(ph[super_pop == "AFR" & pop %in% admixed_african_pops, ID]) else NA_integer_,
		if (all(c("pop", "super_pop") %in% names(ph))) uniqueN(ph[super_pop == "AFR" & !(pop %in% admixed_african_pops), ID]) else NA_integer_))
)
fwrite(diag, out_file("ibdmix_diagnostics.tsv"), sep = "	")
for (i in seq_len(nrow(diag))) message("DIAG ", diag$item[i], "=", diag$value[i])

all_sheets <- list(all_gw = win.burden, all_length = seg.sample)
if (exists("win.g")) all_sheets$all_gw_by_group <- win.g
if (exists("group.percent")) all_sheets$all_percent <- group.percent
if (exists("percent.by_pop")) all_sheets$all_percent_by_pop <- percent.by_pop
if (exists("percent.afr")) all_sheets$all_percent_afr_sensitivity <- percent.afr
if (exists("overlap.summary")) all_sheets$overlap_reduction <- overlap.summary
if (exists("crossref.by_super")) all_sheets$cross_reference_overlap_by_super_pop <- crossref.by_super
write_xlsx(all_sheets, out_file("all.xlsx"))
# Keep the CSV tables as lightweight debug/replot inputs; they are also copied into all.xlsx


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step5: skip genome-wide inherited-frequency PNG plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unlink(out_file(c("all_gw.png", "all_gw_linear.png", "all_gw_circular.png", "all_gw_peaks.png")), force = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Step6: optionally summarize selected/focal loci
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
loci <- parse_loci(selected_loci_arg)
if (nrow(loci) > 0) {
	loci[, chrom := as.character(chrom)]
	selected_diag <- selected_overlap_diagnostics(dat, loci, ph, stats_by_group)
	lo <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
		l <- loci[i]
		z <- dat[chrom == l$chrom & end > l$locus_start & start < l$locus_end,
			.(ID, anc, chrom, seg_start = start, seg_end = end, seg_len, slod, sites)]
		if (nrow(z) == 0) return(data.table())
		z[, `:=`(
			locus_id = l$locus_id,
			locus = l$locus,
			locus_start = l$locus_start,
			locus_end = l$locus_end,
			locus_len = l$locus_len
		)]
		z
	}), fill = TRUE)
	if (nrow(lo) > 0) {
		lo[, overlap_start := pmax(seg_start, locus_start)]
		lo[, overlap_end := pmin(seg_end, locus_end)]
		lo[, overlap_bp := pmax(0L, overlap_end - overlap_start)]
		lo[, locus_cover := overlap_bp / locus_len]
		lo <- lo[locus_cover >= locus_min_cover]
	}

	if (nrow(lo) == 0) {
		message("No locus matches at locus_min_cover >= ", locus_min_cover)
		matches2 <- data.table(
			locus_id = character(), locus = character(), chrom = character(),
			locus_start = integer(), locus_end = integer(), locus_len = integer(),
			ID = character(), anc = character(), n_matching_segments = integer(),
			max_locus_cover = numeric(), total_overlap_bp = integer(),
			best_segment_start = integer(), best_segment_end = integer(),
			best_segment_length = integer(), best_slod = numeric(), best_sites = numeric()
		)
	} else {
		matches <- lo[, .(
			n_matching_segments = .N,
			max_locus_cover = max(locus_cover, na.rm = TRUE),
			total_overlap_bp = sum(overlap_bp, na.rm = TRUE),
			best_segment_start = seg_start[which.max(locus_cover)],
			best_segment_end = seg_end[which.max(locus_cover)],
			best_segment_length = seg_len[which.max(locus_cover)],
			best_slod = slod[which.max(locus_cover)],
			best_sites = sites[which.max(locus_cover)]
		), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, ID, anc)]
		setorder(matches, locus_id, anc, ID)
		matches2 <- matches

		den.locus <- uniqueN(ph$ID)
		locus.by_anc <- matches[, .(n_samples_match = uniqueN(ID)), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc)]
		locus.by_anc[, n_samples_total := den.locus]
		locus.by_anc[, match_fraction := n_samples_match / n_samples_total]
		locus.by_anc[, match_percent := match_fraction * 100]
		setorder(locus.by_anc, locus_id, anc)
		fwrite(locus.by_anc, out_file("selected_overall.csv"))

		if (length(stats_by_group) > 0) {
			mg <- merge(matches, ph, by = "ID", all.x = FALSE)
			den.g <- ph[, .(n_samples_total = uniqueN(ID)), by = stats_by_group]
			locus.g <- mg[, .(n_samples_match = uniqueN(ID)), by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
			locus.g <- merge(locus.g, den.g, by = stats_by_group, all.x = TRUE)
			locus.g[, match_fraction := n_samples_match / n_samples_total]
			locus.g[, match_percent := match_fraction * 100]
			group_tag <- paste(stats_by_group, collapse = "_")
			setorderv(locus.g, c(stats_by_group, "locus_id", "anc"))
			fwrite(locus.g, out_file("selected_by_group.csv"))
		}
	}

	if (nrow(matches2) > 0) {
		matches2[, chrom := as.character(chrom)]
		matches2[, covered_bp := pmin(total_overlap_bp, locus_len)]
	}
	locus.grid <- merge(copy(loci)[, tmp := 1], data.table(anc = ancs, tmp = 1), by = "tmp", allow.cartesian = TRUE)[, tmp := NULL]
	if (nrow(matches2) > 0) {
		locus.by_anc <- matches2[, .(
			n_samples_carrier = uniqueN(ID),
			n_matching_segments = as.integer(sum(n_matching_segments, na.rm = TRUE)),
			total_covered_bp = as.numeric(sum(covered_bp, na.rm = TRUE)),
			median_max_locus_cover = as.numeric(median(max_locus_cover, na.rm = TRUE)),
			median_best_slod = as.numeric(median(best_slod, na.rm = TRUE)),
			median_best_segment_length = as.numeric(median(best_segment_length, na.rm = TRUE))
		), by = .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc)]
	} else {
		locus.by_anc <- locus.grid[0, .(locus_id, locus, chrom, locus_start, locus_end, locus_len, anc,
			n_samples_carrier = integer(), n_matching_segments = integer(), total_covered_bp = numeric(),
			median_max_locus_cover = numeric(), median_best_slod = numeric(),
			median_best_segment_length = numeric())]
	}
	locus.by_anc <- merge(locus.grid, locus.by_anc,
		by = c("locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc"), all.x = TRUE)
	for (cc in c("n_samples_carrier", "n_matching_segments", "total_covered_bp")) {
		locus.by_anc[is.na(get(cc)), (cc) := 0]
	}
	sex.den.all <- ph[, .(
		n_samples_total = uniqueN(ID),
		n_female = uniqueN(ID[.sex_for_denom == "female"]),
		n_male = uniqueN(ID[.sex_for_denom == "male"]),
		n_unknown_sex = uniqueN(ID[!(.sex_for_denom %in% c("female", "male"))])
	)]
	locus.by_anc[, `:=`(
		n_samples_total = sex.den.all$n_samples_total,
		n_female = sex.den.all$n_female,
		n_male = sex.den.all$n_male,
		n_unknown_sex = sex.den.all$n_unknown_sex
	)]
	locus.by_anc[, locus_denominator_bp := mapply(locus_denominator_bp, chrom, locus_start, locus_end,
		n_female, n_male, n_unknown_sex)]
	# Retain the legacy column name for downstream compatibility; on chrX it is
	# now sex/PAR-aware rather than assuming two copies for every individual.
	locus.by_anc[, diploid_locus_bp := locus_denominator_bp]
	locus.by_anc[, locus_burden_fraction := total_covered_bp / locus_denominator_bp]
	locus.by_anc[, locus_burden_percent := 100 * locus_burden_fraction]
	locus.by_anc[, carrier_fraction := n_samples_carrier / n_samples_total]
	locus.by_anc[, carrier_percent := 100 * carrier_fraction]
	setorder(locus.by_anc, locus_id, anc)
	fwrite(locus.by_anc, out_file("selected_overall.csv"))

	if (length(stats_by_group) > 0) {
		group_tag <- paste(stats_by_group, collapse = "_")
		den.g <- ph[, .(
			n_samples_total = uniqueN(ID),
			n_female = uniqueN(ID[.sex_for_denom == "female"]),
			n_male = uniqueN(ID[.sex_for_denom == "male"]),
			n_unknown_sex = uniqueN(ID[!(.sex_for_denom %in% c("female", "male"))])
		), by = stats_by_group]
		group.levels <- unique(ph[, stats_by_group, with = FALSE])
		group.levels[, tmp := 1]
		locus.grid.g <- merge(
			merge(copy(loci)[, tmp := 1], data.table(anc = ancs, tmp = 1), by = "tmp", allow.cartesian = TRUE),
			group.levels, by = "tmp", allow.cartesian = TRUE
		)[, tmp := NULL]

		if (nrow(matches2) > 0) {
			mg <- merge(matches2, ph, by = "ID", all.x = FALSE)
			if (nrow(mg) > 0) {
				locus.g <- mg[, .(
					n_samples_carrier = uniqueN(ID),
					n_matching_segments = as.integer(sum(n_matching_segments, na.rm = TRUE)),
					total_covered_bp = as.numeric(sum(covered_bp, na.rm = TRUE))
				), by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
			} else {
				locus.g <- locus.grid.g[0, .SD, .SDcols = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
				locus.g[, `:=`(n_samples_carrier = integer(), n_matching_segments = integer(), total_covered_bp = numeric())]
			}
		} else {
			locus.g <- locus.grid.g[0, .SD, .SDcols = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc")]
			locus.g[, `:=`(n_samples_carrier = integer(), n_matching_segments = integer(), total_covered_bp = numeric())]
		}

		locus.g <- merge(locus.grid.g, locus.g,
			by = c(stats_by_group, "locus_id", "locus", "chrom", "locus_start", "locus_end", "locus_len", "anc"), all.x = TRUE)
		for (cc in c("n_samples_carrier", "n_matching_segments", "total_covered_bp")) {
			locus.g[is.na(get(cc)), (cc) := 0]
		}
		locus.g <- merge(locus.g, den.g, by = stats_by_group, all.x = TRUE)
		locus.g[, locus_denominator_bp := mapply(locus_denominator_bp, chrom, locus_start, locus_end,
			n_female, n_male, n_unknown_sex)]
		locus.g[, diploid_locus_bp := locus_denominator_bp]
		locus.g[, locus_burden_fraction := total_covered_bp / locus_denominator_bp]
		locus.g[, locus_burden_percent := 100 * locus_burden_fraction]
		locus.g[, carrier_fraction := n_samples_carrier / n_samples_total]
		locus.g[, carrier_percent := 100 * carrier_fraction]
		setorderv(locus.g, c(stats_by_group, "locus_id", "anc"))
		fwrite(locus.g, out_file("selected_by_group.csv"))
	}

	selected_match <- copy(matches2)
	if (nrow(selected_match) > 0 && exists("ph")) {
		add_cols <- intersect(c("ID", "pop", "super_pop"), names(ph))
		selected_match <- merge(selected_match, unique(ph[, add_cols, with = FALSE]),
			by = "ID", all.x = TRUE, sort = FALSE)
		front <- intersect(c("locus_id", "locus", "chrom", "locus_start", "locus_end",
			"locus_len", "ID", "pop", "super_pop", "anc"), names(selected_match))
		selected_match <- selected_match[, c(front, setdiff(names(selected_match), front)), with = FALSE]
	}
	selected_sheets <- list(selected_overall = locus.by_anc)
	if (exists("locus.g")) selected_sheets$selected_by_group <- locus.g
	selected_sheets$selected_match <- selected_match
	if (exists("selected_diag") && is.list(selected_diag)) {
		selected_sheets$all_overlaps_no_cover_filter <- selected_diag$all_overlaps
		selected_sheets$overlap_by_locus_ref_group <- selected_diag$by_locus_ref_group
		if (nrow(selected_diag$by_locus_anc) > 0) selected_sheets$overlap_by_locus_anc <- selected_diag$by_locus_anc
		if (nrow(selected_diag$by_group) > 0) selected_sheets$overlap_by_group_ref_group <- selected_diag$by_group
	}
	write_xlsx(selected_sheets, out_file("selected.xlsx"))

	panel_data <- list(list(title = "Overall", data = locus.by_anc))
	if (exists("locus.g")) {
		locus.g[, group_label := do.call(paste, c(.SD, sep = " / ")), .SDcols = stats_by_group]
		group_title <- paste(stats_by_group, collapse = " / ")
		for (gv in sort(unique(locus.g$group_label))) {
			panel_data[[length(panel_data) + 1L]] <- list(
				title = paste0(group_title, " = ", gv),
				data = locus.g[group_label == gv]
			)
		}
	}
	png(out_file("selected.png"), width = 1850, height = 780 * length(panel_data), units = "px", res = 220)
	op <- par(mfrow = c(length(panel_data), 1), mar = c(5, 7, 3.5, 1.2), oma = c(0, 0, 3.8, 0))
	for (i in seq_along(panel_data)) {
		plot_selected_panel(panel_data[[i]]$data, panel_data[[i]]$title, show_legend = i == 1L)
	}
	mtext("Selected loci: copy-number-adjusted locus burden", outer = TRUE, line = 2.1, font = 2)
	mtext(paste0("Burden = selected-locus bp covered by IBDmix segments / sex- and PAR-aware locus copy bp; ",
		"segments must cover at least ", round(100 * locus_min_cover), "% of the locus."),
		outer = TRUE, line = 0.7, cex = 0.75)
	par(op)
	dev.off()
	unlink(out_file(c("selected_overall.csv", "selected_percent.csv", "selected_by_group.csv")), force = TRUE)
	}

message("Done:")
unlink(bedgraph.dir, recursive = TRUE, force = TRUE)
message("  ", out_file("all.xlsx"))
message("  ", top.file)
if (nrow(loci) > 0) {
	message("  ", out_file("selected.xlsx"))
	message("  ", out_file("selected.png"))
}
