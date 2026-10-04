source(Sys.getenv("GU_RESULTS_R", "/mnt/d/scripts/0f/results.R"))


# 🚩 phyml Circular phylogeny adapted from Bald/archaic introgression/Fig2.R, panel B.  Shared by the batch exporter
# and Shiny. Tree inference is never done here.
gu_b_read <- function(path) {
	if (!file.exists(path))
		return(data.frame())
	as.data.frame(data.table::fread(path, na.strings = c("", "NA"), showProgress = FALSE))
}
gu_b_value <- function(row, name, default = "") {
	x <- row[[name]]
	if (is.null(x) || !length(x) || is.na(x[[1]]))
		default else x[[1]]
}
gu_b_desc <- function(tr, node) {
	kids <- tr$edge[tr$edge[, 1] == node, 2]
	unlist(lapply(kids, function(k) if (k <= ape::Ntip(tr))
		k else gu_b_desc(tr, k)), use.names = FALSE)
}
gu_b_span <- function(a) {
	a <- sort((a + 2 * pi)%%(2 * pi))
	g <- c(diff(a), a[1] + 2 * pi - tail(a, 1))
	j <- which.max(g)
	s <- a[(j%%length(a)) + 1]
	e <- a[j]
	if (e < s)
		e <- e + 2 * pi
	c(s, e)
}
gu_b_wedge <- function(a0, a1, r0, r1, col) {
	th <- seq(a0, a1, length.out = 40)
	polygon(c(r0 * cos(th), rev(r1 * cos(th))), c(r0 * sin(th), rev(r1 * sin(th))), col = col, border = NA)
}
gu_b_ring <- function(final, locus_id) {
	h <- gu_b_read(file.path(final, "gwas_haplotypes.tsv"))
	if (nrow(h))
		h$n <- h$n_copies else h <- gu_b_read(file.path(final, "evidence_haplotypes.tsv"))
	if (!nrow(h))
		return(data.frame())
	h <- h[h$locus_id == locus_id & !duplicated(paste(h$locus_id, h$hap_id)), , drop = FALSE]
	groups <- c("AFR", "AMR", "EAS", "EUR", "SAS", "UNKNOWN")
	z <- data.frame(label = h$hap_id, n = as.numeric(h$n), stringsAsFactors = FALSE)
	for (g in groups) z[[g]] <- 0
	for (i in seq_len(nrow(h))) {
		tokens <- strsplit(as.character(h$superpopulation_copy_counts[i]), ",", fixed = TRUE)[[1]]
		for (token in tokens) {
			pair <- strsplit(token, ":", fixed = TRUE)[[1]]
			if (length(pair) != 2)
				next
			g <- if (pair[1] %in% groups)
				pair[1] else "UNKNOWN"
			count <- suppressWarnings(as.numeric(pair[2]))
			if (is.finite(count))
				z[i, g] <- z[i, g] + count
		}
		# Count chromosome copies, not distinct individuals (unlike old ring_dat).
		missing <- z$n[i] - sum(as.numeric(z[i, groups]))
		if (is.finite(missing) && missing > 0)
			z$UNKNOWN[i] <- z$UNKNOWN[i] + missing
	}
	z
}
gu_b_bundle <- function(row, final, locus_id) {
	nw <- gu_b_value(row, "tree_newick")
	if (!nzchar(nw))
		stop("No completed Newick tree for this locus")
	tr <- ape::read.tree(text = nw)
	lineage <- gu_b_value(row, "expected_lineage", gu_b_value(row, "candidate_lineage", gu_b_value(row, "lineage_filter",
		"ALL")))
	pass <- as.character(gu_b_value(row, "candidate_clade_pass", 0)) %in% c("1", "TRUE")
	mixed_bs <- suppressWarnings(as.numeric(gu_b_value(row, "mixed_lineage_bootstrap", NA_real_)))
	mixed <- !pass && is.finite(mixed_bs)
	fields <- if (pass)
		c("candidate_tips_in_clade", "control_tips_in_clade", "candidate_context_tips_in_clade", "expected_archaic_tips_in_clade") else if (mixed)
		"mixed_lineage_edge_tips" else character()
	labels <- unique(unlist(strsplit(paste(vapply(fields, function(f) as.character(gu_b_value(row, f)), character(1)), collapse = ","),
		",", fixed = TRUE)))
	if ((pass || mixed) && !length(intersect(labels, tr$tip.label)))
		labels <- strsplit(gu_b_value(row, "candidate_edge_tips", gu_b_value(row, "candidate_clade_tips")), ",", fixed = TRUE)[[1]]
	regional <- grepl("region_", gu_b_value(row, "tree_scope")) || grepl("regional_tree", gu_b_value(row, "candidate_clade_rule"))
	if (regional) {
		labels <- character()
		pass <- FALSE
		mixed <- FALSE
	}
	risk <- gu_b_read(file.path(final, "gwas_haplotypes.tsv"))
	risk_labels <- if (nrow(risk))
		risk$hap_id[risk$locus_id == locus_id & risk$role == "risk"] else character()
	list(risk_labels = risk_labels, tree = tr, ring = gu_b_ring(final, locus_id), locus = locus_id, lineage = lineage, side = intersect(labels,
		tr$tip.label), bootstrap = if (pass) suppressWarnings(as.numeric(gu_b_value(row, "candidate_clade_bootstrap", NA_real_))) else if (mixed) mixed_bs else NA_real_,
		state = if (pass) "Lineage-specific tree support" else if (mixed) "Mixed archaic edge; lineage unresolved" else if (regional) "Regional topology; affinity only" else "Predefined GWAS risk haplotypes; tree unconfirmed",
		mixed = mixed)
}
# Zero-length terminal branches still need distinct label angles. Their inferred branch lengths and topology are left
# untouched.
gu_b_tip_angles <- function(p, nt) {
	x <- p$xx[seq_len(nt)]
	y <- p$yy[seq_len(nt)]
	a <- atan2(y, x)
	zero <- which(sqrt(x * x + y * y) < 1e-12)
	used <- a[setdiff(seq_len(nt), zero)]%%(2 * pi)
	for (i in zero) {
		if (!length(used))
			v <- 0 else {
			z <- sort(used)
			gaps <- diff(c(z, z[1] + 2 * pi))
			j <- which.max(gaps)
			v <- (z[j] + gaps[j]/2)%%(2 * pi)
		}
		a[i] <- v
		used <- c(used, v)
	}
	a
}


# 🚩 Saved tree data used by both PNG and XLSX
gu_b_prepare <- function(bundle, min_copies = 11L) {
	tr <- bundle$tree
	ring <- bundle$ring
	original_n <- ape::Ntip(tr)
	arch <- grepl("Altai|Chagyr|Vindija|Neander|Denis|Archaic", tr$tip.label, ignore.case = TRUE)
	counts <- ring$n[match(tr$tip.label, ring$label)]
	keep <- arch | tr$tip.label == "Ancestral" | (!is.na(counts) & counts >= min_copies)
	modern0 <- !arch & tr$tip.label != "Ancestral"
	risk0 <- tr$tip.label %in% bundle$risk_labels
	if (any(risk0) && (!any(keep & risk0) || !any(keep & modern0 & !risk0))) {
		min_copies <- 2L
		keep <- arch | tr$tip.label == "Ancestral" | (!is.na(counts) & counts >= min_copies)
	}
	if (sum(keep) >= 3 && any(!keep))
		tr <- ape::drop.tip(tr, tr$tip.label[!keep])
	if ("Ancestral" %in% tr$tip.label)
		tr <- tryCatch(ape::root(ape::unroot(tr), "Ancestral", resolve.root = TRUE), error = function(e) tr)
	# Preserve ancestral and all other branch lengths; do not shorten to fit.
	tr <- ape::ladderize(tr)
	nt <- ape::Ntip(tr)
	labs <- tr$tip.label
	list(tree = tr, ring = ring, original_n = original_n, min_copies = min_copies)
}

gu_b_write_workbook <- function(tables, path) {
	if (requireNamespace("openxlsx", quietly = TRUE) && requireNamespace("zip", quietly = TRUE)) {
		result_write_workbook(tables, path)
		return(invisible(path))
	}
	# The analysis environment may omit spreadsheet packages. Keep its plotting runtime and use the existing system R
	# spreadsheet runtime for this export.
	temporary <- tempfile("gu-figure-data-", tmpdir = "/tmp", fileext = ".rds")
	on.exit(unlink(temporary), add = TRUE)
	saveRDS(tables, temporary)
	code <- paste0("source(", deparse(Sys.getenv("GU_RESULTS_R", "/mnt/d/scripts/0f/results.R")), "); result_write_workbook(readRDS(",
		deparse(temporary), "), ", deparse(path), ")")
	status <- system2("env", c("-u", "R_HOME", "-u", "R_LIBS", "-u", "R_LIBS_USER", "-u", "R_LIBS_SITE", "-u", "R_ENVIRON_USER",
		"/usr/bin/Rscript", "-e", shQuote(code)))
	if (status != 0L)
		stop("Figure workbook export failed: ", path)
	invisible(path)
}

gu_b_tables <- function(bundle, min_copies = 11L) {
	prepared <- gu_b_prepare(bundle, min_copies)
	tr <- prepared$tree
	labels <- tr$tip.label
	tips <- data.frame(label = labels, risk_haplotype = labels %in% bundle$risk_labels, candidate_side = labels %in% bundle$side,
		stringsAsFactors = FALSE)
	ring <- prepared$ring
	if (nrow(ring))
		tips <- merge(tips, ring, by = "label", all.x = TRUE, sort = FALSE)
	edges <- data.frame(parent_node = tr$edge[, 1], child_node = tr$edge[, 2], branch_length = tr$edge.length)
	edges$child_label <- ifelse(edges$child_node <= length(labels), labels[edges$child_node], "")
	nodes <- data.frame(node = length(labels) + seq_len(tr$Nnode), bootstrap = if (is.null(tr$node.label))
		rep(NA_character_, tr$Nnode) else tr$node.label)
	evidence <- data.frame(locus = bundle$locus, lineage = bundle$lineage, state = bundle$state, candidate_bootstrap = bundle$bootstrap,
		display_min_copies = prepared$min_copies, displayed_tips = length(labels), inferred_tips = prepared$original_n, bar_scale_max_copies = if (nrow(ring))
			max(ring$n, na.rm = TRUE) else 1)
	list(haplotypes = tips, branches = edges, node_support = nodes, evidence = evidence)
}

gu_draw_panel_b <- function(bundle, min_copies = 11L) {
	sp_col <- c(AFR = "#9D76B1", AMR = "#73B28E", EAS = "#BA4E4F", EUR = "#55728E", SAS = "#CA997A", UNKNOWN = "#BDBDBD")
	lin_col <- c(Neanderthal = "#31859C", Denisovan = "#984807")
	anc_col <- "#2F6B3F"
	prepared <- gu_b_prepare(bundle, min_copies)
	tr <- prepared$tree
	ring <- prepared$ring
	original_n <- prepared$original_n
	min_copies <- prepared$min_copies
	nt <- ape::Ntip(tr)
	labs <- tr$tip.label
	ancestral <- labs == "Ancestral"
	arch <- grepl("Altai|Chagyr|Vindija|Neander|Denis|Archaic", labs, ignore.case = TRUE)
	modern <- !ancestral & !arch
	dep <- ape::node.depth.edgelength(tr)[seq_len(nt)]
	r <- max(dep, na.rm = TRUE)
	if (!is.finite(r) || r <= 0) {
		tr$edge.length <- rep(1, nrow(tr$edge))
		r <- max(ape::node.depth.edgelength(tr))
	}
	rr <- r * c(lab = 1.04, text = 1.11, pop0 = 1.47, pop1 = 1.6, bar0 = 1.72, bar1 = 1.9, lim = 2.08)
	oldpar <- par(mar = c(5.2, 0.3, 3.8, 0.3), xpd = NA, pty = "s")
	on.exit(par(oldpar))
	draw <- function(rotation = 0, visible = TRUE) ape::plot.phylo(tr, "fan", use.edge.length = TRUE, show.tip.label = FALSE,
		no.margin = FALSE, edge.width = 0.65, rotate.tree = rotation, plot = visible, x.lim = c(-rr["lim"], rr["lim"]), y.lim = c(-rr["lim"],
			rr["lim"]))
	draw(visible = FALSE)
	coordinates <- get("last_plot.phylo", envir = getFromNamespace(".PlotPhyloEnv", "ape"))
	ia <- which(ancestral)
	rotation <- if (length(ia))
		90 - gu_b_tip_angles(coordinates, nt)[ia] * 180/pi else 0
	draw(rotation)
	p <- get("last_plot.phylo", envir = getFromNamespace(".PlotPhyloEnv", "ape"))
	a <- gu_b_tip_angles(p, nt)
	side <- intersect(bundle$side, labs)
	hi <- if (length(side) >= 2)
		ape::getMRCA(tr, side) else NA_integer_
	col_hi <- if (bundle$mixed)
		"#9C7A9C" else if (bundle$lineage %in% names(lin_col))
		lin_col[bundle$lineage] else "#909090"
	if (is.finite(hi)) {
		ix <- gu_b_desc(tr, hi)
		# A collapsed/pruned MRCA must never shade unrelated retained tips.
		if (!all(labs[ix] %in% side) || !any(modern[ix]) || any(ancestral[ix]))
			hi <- NA_integer_ else {
			s <- gu_b_span(a[ix])
			r0 <- sqrt(p$xx[hi]^2 + p$yy[hi]^2)
			gu_b_wedge(s[1], s[2], r0, rr["lab"], adjustcolor(col_hi, alpha.f = 0.2))
			if (is.finite(bundle$bootstrap))
				ape::nodelabels(bundle$bootstrap, node = hi, frame = "none", cex = 0.85)
		}
	}
	segments(p$xx[seq_len(nt)], p$yy[seq_len(nt)], rr["lab"] * cos(a), rr["lab"] * sin(a), lty = 3, col = "#909090", lwd = 0.35)
	cex <- max(0.25, min(0.85, 6.2/sqrt(nt)))
	deg <- a * 180/pi
	flip <- deg < -90 | deg > 90
	cols <- ifelse(ancestral, anc_col, ifelse(arch, ifelse(grepl("Denis", labs, TRUE), lin_col["Denisovan"], lin_col["Neanderthal"]),
		"#262626"))
	cols[modern & labs %in% bundle$risk_labels] <- "#B2182B"
	for (i in seq_len(nt)) text(rr["text"] * cos(a[i]), rr["text"] * sin(a[i]), labs[i], srt = if (flip[i])
		deg[i] + 180 else deg[i], adj = if (flip[i])
		c(1, 0.5) else c(0, 0.5), cex = cex, col = cols[i])
	th <- seq(0, 2 * pi, length.out = 720)
	for (rad in rr[c("pop0", "pop1", "bar0", "bar1")]) lines(rad * cos(th), rad * sin(th), col = "#D2D2D2", lwd = 0.4)
	max_n <- if (nrow(ring))
		max(ring$n, na.rm = TRUE) else 1
	da <- 2 * pi/nt * 0.36
	for (i in seq_len(nt)) {
		j <- match(labs[i], ring$label)
		if (modern[i] && !is.na(j)) {
			z <- ring[j, , drop = FALSE]
			total <- sum(as.numeric(z[names(sp_col)]))
			r0 <- rr["pop0"]
			for (sp in names(sp_col)) {
				prop <- if (total > 0)
						as.numeric(z[[sp]])/total else 0
				if (prop > 0)
						gu_b_wedge(a[i] - da, a[i] + da, r0, r0 + prop * (rr["pop1"] - rr["pop0"]), sp_col[sp])
				r0 <- r0 + prop * (rr["pop1"] - rr["pop0"])
			}
			h1 <- rr["bar0"] + z$n/max_n * (rr["bar1"] - rr["bar0"])
			gu_b_wedge(a[i] - da, a[i] + da, rr["bar0"], h1, "#FCB4A4")
			if (labs[i] %in% side || nt <= 35)
				text((h1 + 0.025 * r) * cos(a[i]), (h1 + 0.025 * r) * sin(a[i]), z$n, srt = if (flip[i])
						deg[i] + 180 else deg[i], adj = if (flip[i])
						c(1, 0.5) else c(0, 0.5), cex = cex * 0.85)
		} else if (arch[i])
			points(mean(rr[c("pop0", "pop1")]) * cos(a[i]), mean(rr[c("pop0", "pop1")]) * sin(a[i]), pch = 16, col = cols[i],
				cex = 0.8)
	}
	title(main = paste(bundle$locus, bundle$lineage, sep = " | "), cex.main = 1.15, line = 2)
	mtext(paste0(bundle$state, "; red modern labels = GWAS risk allele"), side = 3, line = 0.6, cex = 0.65)
	try(ape::add.scale.bar(x = -rr["lim"] * 0.91, y = -rr["lim"] * 0.91, length = signif(r/5, 2), lwd = 0.8, cex = 0.65),
		silent = TRUE)
	# A separate footer viewport keeps the legend and count note off the rings.
	par(fig = c(0, 1, 0, 0.105), new = TRUE, mar = c(0, 0, 0, 0), pty = "m")
	plot.new()
	plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
	pops <- if (nrow(ring) && any(ring$UNKNOWN > 0))
		names(sp_col) else names(sp_col)[1:5]
	legend(0.5, 0.68, xjust = 0.5, yjust = 0.5, legend = c(pops, "Haplotype copies", "Neanderthal", "Denisovan", "Ancestral"),
		fill = c(sp_col[pops], "#FCB4A4", lin_col, anc_col), border = NA, bty = "n", ncol = 5, cex = 0.7)
	text(0.5, 0.16, sprintf("Display: modern n >= %d; %d / %d tips. Inference uses the full saved tree. Bars: copies (max %s).",
		min_copies, nt, original_n, format(max_n)), cex = 0.64)
	invisible(list(displayed_tips = nt, inferred_tips = original_n, highlight = hi))
}
gu_export_panel_b <- function(out) {
	final <- file.path(out, "final")
	tables <- lapply(c("evidence_trees.tsv", "region_common_trees.tsv"), function(f) gu_b_read(file.path(final, f)))
	exports <- list()
	for (d in tables) if (nrow(d))
		for (i in seq_len(nrow(d))) {
			row <- d[i, , drop = FALSE]
			if (!nzchar(gu_b_value(row, "tree_newick")))
				next
			phy <- gu_b_value(row, "phy_file")
			if (!nzchar(phy))
				next
			# Export inside the requested run even if a restored summary still names its original absolute path or the alignment
			# itself was not archived.
			local_dir <- file.path(out, "loci")
			if (!file.exists(file.path(local_dir, "sites.tsv")))
				local_dir <- file.path(local_dir, row$locus_id[[1]])
			if (!dir.exists(local_dir))
				stop("Missing locus output directory: ", local_dir)
			phy <- file.path(local_dir, basename(phy))
			bundle <- gu_b_bundle(row, final, row$locus_id[[1]])
			stem <- paste0(phy, "_phyml_tree.", gsub("[^A-Za-z0-9_-]", "_", bundle$lineage), ".panelB")
			for (min_copies in c(11L, 2L)) {
				variant <- if (min_copies == 11L)
						"" else ".full"
				dest <- paste0(stem, variant, ".png")
				tmp <- tempfile("gu-tree-", tmpdir = "/tmp", fileext = ".png")
				grDevices::png(tmp, width = 2400, height = 2500, res = 240, type = "cairo")
				tryCatch(gu_draw_panel_b(bundle, min_copies), finally = grDevices::dev.off())
				if (!file.copy(tmp, dest, overwrite = TRUE))
						stop("Could not install ", dest)
				unlink(tmp)
				gu_b_write_workbook(gu_b_tables(bundle, min_copies), paste0(stem, variant, ".xlsx"))
			}
			exports[[length(exports) + 1L]] <- data.frame(locus_id = bundle$locus, lineage = bundle$lineage, tree_file = gu_b_value(row,
				"tree_file"), panel_b_png = paste0(stem, ".png"), panel_b_xlsx = paste0(stem, ".xlsx"), panel_b_full_png = paste0(stem,
				".full.png"), state = bundle$state)
		}
	data.table::fwrite(data.table::rbindlist(exports, fill = TRUE), file.path(final, "phylogeny_panel_b.tsv"), sep = "\t",
		na = "")
	message("PHYML PANEL B: exported ", length(exports), " trees")
}
if (sys.nframe() == 0L) {
	args <- commandArgs(trailingOnly = TRUE)
	i <- match("--out", args)
	if (is.na(i) || i == length(args))
		stop("Usage: Rscript --vanilla phyml_panel_b.R --out PHYML_RUN_DIR")
	gu_export_panel_b(args[[i + 1L]])
}
