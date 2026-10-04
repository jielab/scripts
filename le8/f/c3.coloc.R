# 🚩 c3.coloc
# C3: Colocalization — locus-specific QTL/outcome colocalization and fine-mapping evidence.
suppressPackageStartupMessages({
	fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(fdir, "0.common.R"))
	source(file.path(fdir, "c1.correlate.R"))

	# Cell-specific analyses are executed only by c5_cellulation.
})
LE8_JOB <- "c3_coloc"
C3_CODE_VERSION <- "2026-09-19.mr-fdr-v1"

# c3 mr selection C3 candidates are significant C2 analyses and their actual retained instruments.  This
# avoids re-running COJO coordinate matching (including PGS-only preparation).
C3_SELECTION_VERSION <- "2026-09-19.mr-fdr-v1"
c3_file_stamp <- function(paths) {
	paths <- unique(paths[!is.na(paths) & nzchar(paths)])
	info <- file.info(paths)
	data.frame(
		path = normalizePath(paths, mustWork = FALSE), size = info$size, mtime = as.numeric(info$mtime),
		ctime = as.numeric(info$ctime)
	)
}
c3_select_mr <- function(mr, layer, fdr = 0.05, max_features = 200L) {
	required <- c("exposure", "analysis", "pval", "FDR_all", "instrument_snps", "instrument_positions")
	if (!all(required %in% names(mr)))
		stop("C3 requires C2 MR results with FDR_all and retained instrument IDs/positions; rerun C2.", call. = FALSE)
	classes <- if (layer == "protein")
		c("cis", "trans") else c("local", "distal")
	z <- mr |>
		filter(is.finite(FDR_all), FDR_all < fdr, is.finite(pval), !is.na(exposure), nzchar(exposure), analysis %in%
			classes) |>
		arrange(FDR_all, pval, exposure, analysis)
	if (any(is.na(z$instrument_snps) | !nzchar(z$instrument_snps)))
		stop("Significant C2 MR rows lack retained instrument IDs; rerun C2.", call. = FALSE)
	features <- unique(z$exposure)
	if (length(features) > max_features)
		message("C3: ", length(features), " MR-significant features; C3_MAX_FEATURES retains ", max_features)
	z |>
		filter(exposure %in% head(features, max_features))
}
c3_mr_instruments <- function(mr) {
	bind_rows(lapply(seq_len(nrow(mr)), function(i) {
		ids <- trimws(strsplit(mr$instrument_snps[i], ";", fixed = TRUE)[[1]])
		pos <- strsplit(ifelse(is.na(mr$instrument_positions[i]), "", mr$instrument_positions[i]), ";", fixed = TRUE)[[1]]
		if (length(pos) != length(ids))
			pos <- rep(NA_character_, length(ids))
		tibble(SNP = ids, analysis = mr$analysis[i], mr_position = pos)
	})) |>
		distinct()
}
c3_read_mr_qtl <- function(mr, qfile, cache_root) {
	wanted <- c3_mr_instruments(mr)
	key <- le8_hash_object(list(version = C3_SELECTION_VERSION, file = c3_file_stamp(qfile), wanted = wanted, metadata = c3_file_stamp(Sys.getenv(
		"LE8_GWAS_MANIFEST",
		""
	)), N = Sys.getenv("C2_SUMSTAT_N", "100000")))
	cache <- file.path(cache_root, paste0(key, ".rds"))
	iv <- read_stage_cache(cache)
	if (is.data.frame(iv)) {
		message("  QTL instrument cache hit: ", nrow(iv), " records")
		return(iv)
	}
	# Positions in C2 are from marginal QTL records, already on the QTL build.  Query indexed intervals first;
	# missing/legacy coordinates get one ID scan.
	positions <- unique(wanted$mr_position)
	positions <- positions[!is.na(positions) & grepl("^[^:]+:[0-9]+$", positions)]
	tabix_ok <- nzchar(Sys.which("tabix")) && any(vapply(paste0(qfile, c(".tbi", ".csi")), function(f) file.exists(f) &&
		file.info(f)$mtime >= file.info(qfile)$mtime, logical(1)))
	q <- tibble(SNP = character())
	if (tabix_ok && length(positions)) {
		message("  QTL indexed lookup: ", length(positions), " MR instrument positions")
		q <- bind_rows(lapply(positions, function(p) {
			s <- strsplit(p, ":", fixed = TRUE)[[1]]
			read_sumstat_region(qfile, s[1], as.numeric(s[2]), as.numeric(s[2]))
		}))
	}
	missing <- setdiff(wanted$SNP, q$SNP)
	if (length(missing)) {
		message("  QTL ID scan: ", length(missing), " instruments absent from indexed lookup")
		q <- bind_rows(q, read_sumstat_snps(qfile, missing))
	}
	if (!nrow(q))
		q <- tibble(SNP = character(), CHR = character(), POS = numeric(), P = numeric())
	q <- q |>
		distinct(SNP, .keep_all = TRUE)
	iv <- wanted |>
		inner_join(q, by = "SNP")
	lost <- setdiff(wanted$SNP, iv$SNP)
	if (length(lost))
		stop("C3 cannot recover retained MR instruments from ", qfile, ": ", paste(lost, collapse = ", "), "; check C2/QTL provenance.",
			call. = FALSE
		)
	if (any(!is.finite(iv$POS) | is.na(iv$CHR) | !is.finite(iv$P)))
		stop("C3 retained MR instruments have invalid QTL coordinates/P values: ", qfile, call. = FALSE)
	write_stage_cache(iv, cache)
	iv
}
c3_locus_settings <- function() list(window = WINDOW_BP, min_snps = MIN_SNPS, p12 = C3_P12, susie = Sys.getenv(
	"C3_RUN_SUSIE",
	"FALSE"
), susie_max = Sys.getenv("C3_SUSIE_MAX_SNPS", "5000"), ld = c3_file_stamp(Sys.getenv(
	"LE8_LD_MANIFEST",
	""
)), metadata = c3_file_stamp(Sys.getenv("LE8_GWAS_MANIFEST", "")), N = Sys.getenv("C2_SUMSTAT_N", "100000"))
c3_cached_locus <- function(cache, legacy_files, feature, layer, chr, pos, cls, qfile, yfile, outtype, sfrac) {
	z <- read_stage_cache(cache)
	if (!is.null(z))
		return(z)
	# Legacy caches have no input signature. Migrate only default ABF settings, unchanged input files and
	# exactly matching geometry/class; never a SuSiE run.
	legacy_ok <- !LE8_REPLACE && WINDOW_BP == 5e+05 && MIN_SNPS == 50 && identical(unname(C3_P12), c(
		1e-06, 1e-05,
		1e-04
	)) && !truthy(Sys.getenv("C3_RUN_SUSIE", "FALSE")) && Sys.getenv("C2_SUMSTAT_N", "100000") == "100000" &&
		!nzchar(Sys.getenv("LE8_GWAS_MANIFEST", "")) && !is.finite(sfrac)
	if (!legacy_ok)
		return(NULL)
	suffix <- paste0("_", gsub("[^A-Za-z0-9._-]", "_", feature), "_chr", chr, "_", format(pos,
		scientific = FALSE,
		trim = TRUE
	), ".rds")
	paths <- legacy_files[endsWith(basename(legacy_files), suffix)]
	for (f in paths) {
		if (any(file.info(c(qfile, yfile))$mtime > file.info(f)$mtime))
			next
		old <- read_stage_cache(f)
		s <- old$summary
		if (is.data.frame(s) && nrow(s) == 1L && all(c("locus_class", "case_fraction_source") %in% names(s)) &&
			identical(as.character(s$feature), feature) && identical(as.character(s$layer), layer) && identical(
			as.character(s$locus_class),
			cls
		) && as.character(s$chr) == as.character(chr) && s$lead_pos == pos && s$start == max(1, pos - WINDOW_BP) &&
			s$end == pos + WINDOW_BP && s$case_fraction_source == "not required for beta/varbeta cc ABF" && all(c(
			"variants",
			"regional"
		) %in% names(old))) {
			write_stage_cache(old, cache)
			return(old)
		}
	}
	NULL
}


C3_MR_FDR <- as.numeric(Sys.getenv("C3_MR_FDR", unset = "0.05"))
if (!LE8_REUSE_RESULTS) suppressPackageStartupMessages(pacman::p_load(coloc))
MAX_FEATURES <- as.integer(Sys.getenv("C3_MAX_FEATURES", unset = "200"))
MAX_LOCI_PER_FEATURE <- as.integer(Sys.getenv("C3_MAX_LOCI_PER_FEATURE", unset = "3"))
WINDOW_BP <- as.numeric(Sys.getenv("COLOC_WINDOW_BP", unset = "500000"))
MIN_SNPS <- as.integer(Sys.getenv("C3_MIN_NSNP", unset = "50"))
H4_STRONG <- as.numeric(Sys.getenv("C3_H4", unset = "0.70"))
C3_P12 <- c(conservative = 1e-06, default = 1e-05, liberal = 1e-04)

run_gpu_coloc_step <- function(layer, rawdir, manifest, ygfile, outtype, gpudir = le8_cache_dir("gpu_coloc", basename(dirname(rawdir)))) {
	dir.create(gpudir, recursive = TRUE, showWarnings = FALSE)
	result_file <- file.path(gpudir, "gpu_coloc.results.tsv")
	if (!truthy(Sys.getenv("RUN_GPU_COLOC", unset = "TRUE")))
		return(invisible(0L))
	if (cache_valid(result_file)) {
		cache_message(paste0("GPU-coloc/", layer), result_file)
		return(invisible(0L))
	}
	if (!nrow(manifest))
		return(invisible(0L))
	sh <- file.path(Sys.getenv("LE8_FDIR"), "c3.coloc_GPU.sh")
	manifest_file <- file.path(rawdir, "qtl_cad_manifest.tsv")
	write_raw_tsv(manifest, basename(manifest_file), rawdir)
	status <- system2("bash", c(
		sh, "--qtl-manifest", manifest_file, "--cad-gwas", ygfile, "--outdir", gpudir,
		"--outcome-type", outtype, "--H4", as.character(H4_STRONG)
	))
	if (status != 0)
		warning("GPU-coloc runner returned status ", status)
	invisible(status)
}

read_gpu_coloc_results <- function(rawdir, gpudir = le8_cache_dir("gpu_coloc", basename(dirname(rawdir)))) {
	rf <- file.path(gpudir, "gpu_coloc.results.tsv")
	sf <- file.path(gpudir, "signal_preparation_status.tsv")
	ans <- list(results = tibble(), status = tibble())
	if (!truthy(Sys.getenv("RUN_GPU_COLOC", unset = "TRUE")))
		return(ans)
	if (file.exists(rf) && file.size(rf) > 0)
		ans$results <- tryCatch(as_tibble(data.table::fread(rf, showProgress = FALSE, check.names = FALSE)), error = function(e) tibble())
	if (file.exists(sf) && file.size(sf) > 0)
		ans$status <- tryCatch(as_tibble(data.table::fread(sf, showProgress = FALSE, check.names = FALSE)), error = function(e) tibble())
	d <- ans$results
	if (nrow(d)) {
		pp <- names(d)[toupper(names(d)) == "PP.H4"][1]
		if (is.na(pp))
			pp <- names(d)[str_detect(toupper(names(d)), "PP.*H4")][1]
		sigcols <- names(d)[vapply(d, function(x) any(str_detect(as.character(x), "QTL__"), na.rm = TRUE), logical(1))]
		if ("trait" %in% names(d))
			d$feature <- as.character(d$trait) else if (length(sigcols)) {
			sig <- as.character(d[[sigcols[1]]])
			trait <- str_match(sig, "^QTL__[^_]+__(.*?)__chr")
			d$feature <- if (ncol(trait) >= 2)
				trait[, 2] else NA_character_
		} else d$feature <- NA_character_
		d$GPU_PP.H4 <- if (length(pp) && !is.na(pp))
			suppressWarnings(as.numeric(d[[pp]])) else NA_real_
		ans$results <- d
	}
	ans
}

plot_gpu_coloc_validation <- function(gpu, abf, outdir) {
	d <- gpu$results
	st <- gpu$status
	if (nrow(d) && "GPU_PP.H4" %in% names(d)) {
		if ("region" %in% names(d) && "locus" %in% names(abf)) {
			z <- d |>
				filter(is.finite(GPU_PP.H4)) |>
				group_by(feature, region) |>
				slice_max(GPU_PP.H4, n = 1, with_ties = FALSE) |>
				ungroup() |>
				left_join(abf |>
					filter(status == "ok") |>
					transmute(feature, region = locus, ABF_PP.H4 = PP.H4), by = c("feature", "region"))
		} else {
			z <- d |>
				filter(is.finite(GPU_PP.H4)) |>
				group_by(feature) |>
				slice_max(GPU_PP.H4, n = 1, with_ties = FALSE) |>
				ungroup() |>
				left_join(abf |>
					filter(status == "ok") |>
					group_by(feature) |>
					summarise(ABF_PP.H4 = max(PP.H4, na.rm = TRUE), .groups = "drop"), by = "feature") |>
				mutate(region = NA_character_)
		}
		z <- z |>
			arrange(desc(GPU_PP.H4)) |>
			slice_head(n = 40) |>
			mutate(
				label = ifelse(is.na(region) | !nzchar(region), feature, paste(feature, region, sep = " | ")),
				label = factor(label, levels = rev(label))
			)
		pa <- if (!nrow(z))
			blank_plot("a. GPU-coloc regional results") else ggplot(z, aes(GPU_PP.H4, label)) +
			geom_vline(xintercept = H4_STRONG, linetype = 2, color = "grey55") +
			geom_segment(aes(x = 0, xend = GPU_PP.H4, yend = label), color = "grey75") +
			geom_point(aes(color = GPU_PP.H4 >=
				H4_STRONG), size = 2.2) +
			scale_color_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "grey65"), guide = "none") +
			scale_x_continuous(limits = c(0, 1), labels = label_percent()) +
			labs(
				title = "a. GPU-coloc regional marginal-signal evidence",
				subtitle = "Exact manifest QTL-outcome pairs only", x = "GPU-coloc PP(H4)", y = NULL
			) +
			theme_5c(8)
		pb <- if (!nrow(z) || !any(is.finite(z$ABF_PP.H4)))
			blank_plot("b. ABF versus GPU-coloc", "No feature could be aligned across methods") else ggplot(z, aes(ABF_PP.H4, GPU_PP.H4)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
			geom_hline(yintercept = H4_STRONG, linetype = 3, color = "grey70") +
			geom_vline(
				xintercept = H4_STRONG,
				linetype = 3, color = "grey70"
			) +
			geom_point(color = "#6A3D9A", size = 2) +
			ggrepel::geom_text_repel(aes(label = label),
				size = 2.5, max.overlaps = 15, seed = 63
			) +
			coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
			labs(
				title = "b. coloc.abf versus GPU-coloc",
				subtitle = "Same marginal locus and p12; disagreement flags input or prior-scale differences", x = "coloc.abf PP(H4)",
				y = "GPU-coloc PP(H4)"
			) +
			theme_5c(9)
	} else {
		pa <- blank_plot("a. GPU-coloc regional results", "GPU-coloc was disabled, unavailable, or produced no PP(H4) column")
		pb <- blank_plot("b. ABF versus GPU-coloc")
	}
	pc <- if (!nrow(st) || !"status" %in% names(st))
		blank_plot("c. GPU-coloc preparation audit", "No status file") else st |>
		count(status) |>
		ggplot(aes(n, fct_reorder(status, n), fill = status)) +
		geom_col() +
		geom_text(aes(label = n),
			hjust =  - 0.1,
			fontface = "bold"
		) +
		scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
		labs(
			title = "c. Signal-preparation audit",
			x = "Loci", y = NULL, fill = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "none")
	save_plot((pa | pb) / pc + plot_layout(heights = c(1, 0.55)), "c3.Fig6.gpu_coloc_validation.png", 15.5, 11, outdir = outdir)
}

infer_outcome_type <- function(y) {
	explicit <- tolower(Sys.getenv("C3_OUTCOME_TYPE", unset = "auto"))
	if (explicit %in% c("cc", "quant"))
		return(explicit)
	if (str_detect(tolower(y), "bmi|height|weight|bp$|ldl|hdl|hba1c|glucose|age"))
		"quant" else "cc"
}


# Retained for reproducibility of the original four-panel output.  The publication function below is the
# active implementation.
plot_coloc_results_legacy <- function(res, regional, variants, layer, outdir) {
	ok <- res |>
		filter(status == "ok") |>
		arrange(desc(PP.H4))
	if (!nrow(ok)) {
		blank_files <- c(
			"c3.Fig1.posterior_lollipop.png", "c3.Fig2.posterior_hypotheses.png", "c3.Fig3.regional_top_loci.png",
			"c3.Fig4.coloc_network.png"
		)
		for (i in seq_along(blank_files)) save_plot(blank_plot(paste0("C3 Figure ", i), "No locus passed the minimum aligned-SNP requirement"),
			blank_files[[i]], 9, 5,
			outdir = outdir
		)
		return(invisible(NULL))
	}
	top <- ok |>
		slice_head(n = 50) |>
		mutate(label = paste(feature, locus, sep = " | "), label = factor(label, levels = rev(label)), tier = case_when(PP.H4 >=
			0.8 ~ "Strong", PP.H4 >= H4_STRONG ~ "Probable", PP.H4 >= 0.5 ~ "Suggestive", TRUE ~ "Weak"))
	p1 <- ggplot(top, aes(PP.H4, label)) +
		geom_segment(aes(x = 0, xend = PP.H4, y = label, yend = label, color = tier),
			linewidth = 0.75
		) +
		geom_point(aes(color = tier, size = n_snps), alpha = 0.9) +
		geom_vline(
			xintercept = H4_STRONG,
			linetype = 2, color = "grey45"
		) +
		scale_color_manual(values = c(
			Strong = "#1B9E77", Probable = "#66A61E",
			Suggestive = "#E6AB02", Weak = "grey70"
		)) +
		scale_x_continuous(limits = c(0, 1), labels = label_percent()) +
		labs(
			title = "Colocalization posterior probability", subtitle = "Each row is a QTL locus tested against the disease GWAS",
			x = "PP(H4): shared causal signal", y = NULL, color = NULL, size = "Aligned SNPs"
		) +
		theme_5c(9) +
		theme(legend.position = "bottom")
	save_plot(p1, "c3.Fig1.posterior_lollipop.png", 10, 10, outdir = outdir)

	post <- top |>
		select(label, starts_with("PP.H")) |>
		pivot_longer(starts_with("PP.H"), names_to = "hypothesis", values_to = "posterior") |>
		mutate(hypothesis = factor(hypothesis, levels = paste0("PP.H", 0 : 4)))
	p2 <- ggplot(post, aes(label, posterior, fill = hypothesis)) +
		geom_col(width = 0.78) +
		coord_flip() +
		scale_fill_brewer(palette = "Set2") +
		scale_y_continuous(labels = label_percent()) +
		labs(
			title = "Competing colocalization hypotheses", x = NULL,
			y = "Posterior probability", fill = "Hypothesis"
		) +
		theme_5c(9) +
		theme(legend.position = "bottom")
	save_plot(p2, "c3.Fig2.posterior_hypotheses.png", 10, 10, outdir = outdir)

	loci <- ok |>
		slice_head(n = 6) |>
		pull(locus)
	rr <- regional |>
		filter(locus %in% loci, is.finite(p), p > 0) |>
		mutate(y =  - log10(pmax(p, 1e-300)), position_mb = position / 1e+06)
	p3 <- ggplot(rr, aes(position_mb, y, color = dataset)) +
		geom_point(size = 0.8, alpha = 0.72) +
		geom_vline(data = ok |>
			filter(locus %in% loci), aes(xintercept = lead_pos / 1e+06), inherit.aes = FALSE, linetype = 2, color = "grey55") +
		facet_grid(dataset ~ feature + locus, scales = "free_x", space = "free_x") +
		scale_color_manual(values = c(
			QTL = "#4C78A8",
			setNames("#E45756", Y)
		)) +
		labs(
			title = "Regional QTL–disease comparison at top loci", x = "Position (Mb)",
			y = expression( - log[10](P)), color = NULL
		) +
		theme_5c(8) +
		theme(legend.position = "bottom", axis.text.x = element_text(
			angle = 45,
			hjust = 1
		))
	save_plot(p3, "c3.Fig3.regional_top_loci.png", 17, 7.8, outdir = outdir)

	# Evidence network: feature nodes on the outside, genomic loci inside; edge width is PP4.
	net <- ok |>
		filter(PP.H4 >= 0.5) |>
		distinct(feature, locus, PP.H4, chr, lead_pos)
	if (!nrow(net))
		p4 <- blank_plot("Colocalization evidence network", "No locus had PP(H4) >= 0.5") else {
		fnode <- tibble(name = unique(net$feature), type = "feature", angle = seq(0, 2 * pi, length.out = n_distinct(net$feature) +
			1)[ - (n_distinct(net$feature) + 1)], r = 1)
		lnode <- tibble(name = unique(net$locus), type = "locus", angle = seq(0, 2 * pi, length.out = n_distinct(net$locus) +
			1)[ - (n_distinct(net$locus) + 1)] + pi / 10, r = 0.48)
		nodes <- bind_rows(fnode, lnode) |>
			mutate(x = r * cos(angle), y = r * sin(angle))
		edges <- net |>
			left_join(nodes |>
				select(feature = name, x1 = x, y1 = y), by = "feature") |>
			left_join(nodes |>
				select(locus = name, x2 = x, y2 = y), by = "locus")
		p4 <- ggplot() +
			geom_curve(data = edges, aes(
				x = x1, y = y1, xend = x2, yend = y2, linewidth = PP.H4,
				color = PP.H4
			), curvature = 0.15, alpha = 0.65) +
			geom_point(
				data = nodes, aes(x, y, shape = type),
				size = 3, fill = "white"
			) +
			ggrepel::geom_text_repel(
				data = nodes, aes(x, y, label = name), size = 2.5,
				max.overlaps = 40, seed = 4
			) +
			scale_linewidth(range = c(0.4, 2.5), guide = "none") +
			scale_color_gradient(
				low = "#FEE8C8",
				high = "#B30000", limits = c(0.5, 1), name = "PP(H4)"
			) +
			coord_equal() +
			theme_void() +
			labs(title = "Shared-locus evidence network") +
			theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
	}
	save_plot(p4, "c3.Fig4.coloc_network.png", 11, 9, outdir = outdir)
}

# Publication figures cover evidence, posterior, regional, and prior diagnostics.
plot_coloc_results <- function(res, regional, variants, layer, outdir) {
	ok <- res |>
		filter(status == "ok") |>
		arrange(desc(PP.H4))
	files <- c(
		"c3.Fig1.evidence_triage.png", "c3.Fig2.posterior_diagnostics.png", "c3.Fig3.regional_top_loci.png",
		"c3.Fig4.credible_sets.png", "c3.Fig5.prior_sensitivity.png"
	)
	if (!nrow(ok)) {
		for (i in seq_along(files)) save_plot(blank_plot(paste0("C3 Figure ", i), "No locus passed the aligned-SNP requirement"),
			files[[i]], 9, 5,
			outdir = outdir
		)
		return(invisible(NULL))
	}
	for (nm in c("PP.H4_p12_conservative", "PP.H4_p12_default", "PP.H4_p12_liberal", "PP.H4_robust_min", "lead_shared_pp")) if (!nm %in%
		names(ok))
		ok[[nm]] <- if (nm == "PP.H4_p12_default")
			ok$PP.H4 else NA_real_
	ok <- ok |>
		arrange(desc(coalesce(PP.H4_robust_min, PP.H4)), desc(PP.H4))
	top <- ok |>
		slice_head(n = 30) |>
		mutate(
			label = paste(feature, str_replace(locus, "^chr", "chr"), sep = " | "), label = factor(label, levels = rev(label)),
			robust = coalesce(PP.H4_robust_min, PP.H4) >= H4_STRONG, lo = pmin(PP.H4_p12_conservative, PP.H4_p12_default,
				PP.H4_p12_liberal,
				na.rm = TRUE
			), hi = pmax(PP.H4_p12_conservative, PP.H4_p12_default, PP.H4_p12_liberal,
				na.rm = TRUE
			)
		)
	p1 <- ggplot(top, aes(PP.H4, label)) +
		geom_vline(xintercept = H4_STRONG, linetype = 2, color = "grey50") +
		geom_errorbarh(aes(xmin = lo, xmax = hi, color = robust), height = 0.1, linewidth = 0.8) +
		geom_point(aes(
			size = lead_shared_pp,
			color = robust
		), alpha = 0.9) +
		scale_color_manual(
			values = c(`TRUE` = "#1B9E77", `FALSE` = "#D95F02"),
			labels = c(`TRUE` = "Robust to p12", `FALSE` = "Prior-sensitive")
		) +
		scale_x_continuous(
			limits = c(0, 1),
			labels = label_percent()
		) +
		scale_size_continuous(range = c(1.8, 5), name = "Lead SNP PP") +
		labs(
			title = "Colocalization evidence triage",
			subtitle = "Point: default PP(H4); interval: p12 = 1e-6 to 1e-4", x = "Posterior probability of a shared signal",
			y = NULL, color = NULL
		) +
		theme_5c(8) +
		theme(legend.position = "bottom")
	save_plot(p1, files[[1]], 11.5, 10, outdir = outdir)

	t25 <- top |>
		slice_head(n = 25)
	post <- t25 |>
		select(label, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4) |>
		pivot_longer(starts_with("PP.H"), names_to = "hypothesis", values_to = "posterior") |>
		mutate(hypothesis = factor(hypothesis, levels = paste0("PP.H", 0 : 4)))
	pa <- ggplot(post, aes(hypothesis, label, fill = posterior)) +
		geom_tile(color = "white", linewidth = 0.25) +
		scale_fill_viridis_c(limits = c(0, 1), labels = label_percent()) +
		labs(
			title = "a. Competing hypotheses",
			x = NULL, y = NULL, fill = "Posterior"
		) +
		theme_5c(8)
	pb <- ggplot(t25, aes(credible_set_n, n_snps, color = PP.H4, size = lead_shared_pp)) +
		geom_point(alpha = 0.82) +
		ggrepel::geom_text_repel(aes(label = as.character(label)), size = 2.3, seed = 52, max.overlaps = 12) +
		scale_x_log10() +
		scale_y_log10() +
		scale_color_viridis_c(limits = c(0, 1), name = "PP(H4)") +
		labs(
			title = "b. Shared-signal posterior resolution",
			subtitle = "Small sets and high lead-SNP posterior are preferred, conditional on the coloc single-signal model",
			x = "95% shared-signal posterior-set size", y = "Aligned SNPs", size = "Lead SNP PP"
		) +
		theme_5c(8)
	save_plot(pa | pb, files[[2]], 16, 9.5, outdir = outdir)

	loc <- ok |>
		slice_head(n = 4) |>
		select(feature, locus, lead_pos)
	rr <- regional |>
		semi_join(loc, by = c("feature", "locus")) |>
		filter(is.finite(p), p > 0) |>
		mutate(y =  - log10(pmax(p, 1e-300)), position_mb = position / 1e+06, panel = paste(feature, locus, sep = " | "))
	p3 <- if (!nrow(rr))
		blank_plot("Regional QTL–disease comparisons") else ggplot(rr, aes(position_mb, y, color = dataset)) +
		geom_line(linewidth = 0.35, alpha = 0.5) +
		geom_point(
			size = 0.85,
			alpha = 0.72
		) +
		geom_vline(
			data = loc |>
				mutate(panel = paste(feature, locus, sep = " | ")), aes(xintercept = lead_pos / 1e+06), inherit.aes = FALSE,
			linetype = 2, color = "grey45"
		) +
		facet_wrap( ~ panel, ncol = 2, scales = "free") +
		scale_color_manual(values = c(
			QTL = "#4C78A8",
			setNames("#E45756", Y)
		)) +
		labs(
			title = "Regional QTL–disease comparisons", subtitle = "Four strongest loci; dashed line is the selected QTL lead",
			x = "Position (Mb)", y = expression( - log[10](P)), color = NULL
		) +
		theme_5c(8) +
		theme(legend.position = "bottom")
	save_plot(p3, files[[3]], 14.5, 10, outdir = outdir)

	vv <- variants |>
		semi_join(loc, by = c("feature", "locus"))
	if (nrow(vv) && all(c("SNP.PP.H4", "snp") %in% names(vv))) {
		vv <- vv |>
			group_by(feature, locus) |>
			arrange(desc(SNP.PP.H4), .by_group = TRUE) |>
			slice_head(n = 40) |>
			ungroup() |>
			mutate(rank = ave( - SNP.PP.H4, interaction(feature, locus), FUN = rank), panel = paste(feature, locus,
				sep = " | "
			), credible95 = coalesce(credible95, FALSE), posterior_plot = pmax(coalesce(
				SNP.PP.H4,
				0
			), 1e-10))
		p4a <- ggplot(vv, aes(rank, posterior_plot, color = credible95)) +
			geom_segment(aes(
				xend = rank, y = 1e-10,
				yend = posterior_plot
			), linewidth = 0.4) +
			geom_point(size = 1.55) +
			facet_wrap( ~ panel, ncol = 2, scales = "free_x") +
			scale_y_log10(limits = c(1e-10, 1), breaks = 10 ^ c( - 10, - 8, - 6, - 4, - 2, 0)) +
			scale_color_manual(values = c(
				`TRUE` = "#1B9E77",
				`FALSE` = "grey70"
			)) +
			labs(
				title = "a. Variant-level shared-signal posterior (log scale)", subtitle = "A posterior near 1 with all alternatives near 0 is visible rather than appearing as a blank panel",
				x = "Variant rank within locus", y = "SNP posterior under H4", color = "H4-conditional 95% set"
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")
	} else p4a <- blank_plot("a. Variant-level shared-signal posterior", "No variant posterior was available")
	cs <- ok |>
		filter(is.finite(credible_set_n), credible_set_n > 0, is.finite(PP.H4)) |>
		mutate(robust = coalesce(PP.H4_robust_min, PP.H4) >= H4_STRONG, degenerate = credible_set_n == 1 & coalesce(
			lead_shared_pp,
			0
		) > 0.999, label = ifelse(feature %in% loc$feature | degenerate & PP.H4 >= H4_STRONG, paste(feature,
			locus,
			sep = " | "
		), NA_character_))
	p4b <- if (!nrow(cs))
		blank_plot("b. Credible-set resolution", "No credible-set summary was available") else ggplot(cs, aes(credible_set_n, coalesce(PP.H4_robust_min, PP.H4), color = degenerate, size = n_snps)) +
		geom_hline(yintercept = H4_STRONG, linetype = 2, color = "grey55") +
		geom_point(alpha = 0.62) +
		ggrepel::geom_text_repel(aes(label = label),
			size = 2.2, max.overlaps = 15, seed = 72
		) +
		scale_x_log10() +
		scale_color_manual(values = c(
			`TRUE` = "#D7301F",
			`FALSE` = "#3F78A8"
		), labels = c(`TRUE` = "1-SNP posterior concentration", `FALSE` = "Other")) +
		scale_y_continuous(limits = c(
			0,
			1
		), labels = label_percent()) +
		labs(
			title = "b. Resolution across every tested locus", x = "H4-conditional 95% set size (log scale)",
			y = "Robust minimum PP(H4)", color = NULL, size = "Aligned SNPs"
		) +
		theme_5c(8) +
		theme(legend.position = "bottom")
	p4c <- if (!nrow(cs))
		blank_plot("c. Credible-set size distribution") else ggplot(cs, aes(credible_set_n, color = robust)) +
		stat_ecdf(linewidth = 0.9) +
		scale_x_log10() +
		scale_color_manual(values = c(
			`TRUE` = "#1B9E77",
			`FALSE` = "grey55"
		), labels = c(`TRUE` = "Robust shared signal", `FALSE` = "Not robust")) +
		labs(
			title = "c. Shared-signal resolution audit",
			subtitle = paste0(sum(cs$degenerate), " / ", nrow(cs), " loci have a one-SNP set with lead PP > 0.999; verify LD and harmonization"),
			x = "H4-conditional 95% set size (log scale)", y = "Cumulative fraction of loci", color = NULL
		) +
		theme_5c(8) +
		theme(legend.position = "bottom")
	p4 <- p4a / plot_spacer() / (p4b | p4c) + plot_layout(heights = c(1.25, 0.07, 1))
	save_plot(p4, files[[4]], 17, 10, outdir = outdir)

	sens <- top |>
		mutate(highlight = row_number() <= 10) |>
		select(label, highlight, conservative = PP.H4_p12_conservative, default = PP.H4_p12_default, liberal = PP.H4_p12_liberal) |>
		pivot_longer(c(conservative, default, liberal), names_to = "prior", values_to = "PP4") |>
		mutate(prior = factor(prior, levels = c("conservative", "default", "liberal"), labels = c(
			"p12=1e-6", "p12=1e-5",
			"p12=1e-4"
		)))
	p5 <- ggplot(sens, aes(prior, PP4, group = label, color = highlight)) +
		geom_hline(
			yintercept = H4_STRONG,
			linetype = 2, color = "grey55"
		) +
		geom_line(alpha = 0.48) +
		geom_point(size = 1.7) +
		scale_color_manual(values = c(
			`TRUE` = "#D95F02",
			`FALSE` = "grey72"
		), guide = "none") +
		scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
		labs(
			title = "Shared-signal prior sensitivity", subtitle = "Highlighted lines are the ten strongest default-prior loci",
			x = "Prior probability that one variant affects both traits", y = "PP(H4)"
		) +
		theme_5c(10)
	save_plot(p5, files[[5]], 11.5, 7.5, outdir = outdir)
}

credible_set_audit <- function(res, variants = tibble()) {
	ok <- res |>
		filter(status == "ok")
	locus <- if (!nrow(variants) || !all(c("feature", "locus", "SNP.PP.H4") %in% names(variants)))
		tibble() else variants |>
		group_by(feature, locus) |>
		arrange(desc(SNP.PP.H4), .by_group = TRUE) |>
		summarise(
			variant_rows = n(), posterior_sum = sum(SNP.PP.H4, na.rm = TRUE), lead_pp = max(SNP.PP.H4, na.rm = TRUE),
			second_pp = ifelse(n() > 1, nth(SNP.PP.H4, 2), NA_real_), zero_pp_fraction = mean(coalesce(
				SNP.PP.H4,
				0
			) <= 1e-12), .groups = "drop"
		)
	by_locus <- ok |>
		select(feature, locus, n_snps, PP.H4, PP.H4_robust_min, lead_shared_pp, credible_set_n) |>
		left_join(locus, by = c("feature", "locus")) |>
		mutate(robust = coalesce(PP.H4_robust_min, PP.H4) >= H4_STRONG, one_snp_concentration = credible_set_n ==
			1 & coalesce(lead_shared_pp, lead_pp, 0) > 0.999, audit_flag = case_when(
			one_snp_concentration ~ "verify LD/harmonization: one-SNP posterior concentration",
			!is.finite(credible_set_n) | credible_set_n < 1 ~ "missing credible set", TRUE ~ "ok"
		))
	overall <- tibble(metric = c(
		"tested loci", "robust shared-signal loci", "median credible-set size", "one-SNP posterior concentration",
		"missing credible set"
	), value = c(nrow(by_locus), sum(by_locus$robust, na.rm = TRUE), median(by_locus$credible_set_n,
		na.rm = TRUE
	), sum(by_locus$one_snp_concentration, na.rm = TRUE), sum(!is.finite(by_locus$credible_set_n) |
		by_locus$credible_set_n < 1)))
	list(overall = overall, by_locus = by_locus)
}

# Locus evidence and optional SuSiE analysis

select_nonoverlapping_leads <- function(iv, max_loci = MAX_LOCI_PER_FEATURE, window_bp = WINDOW_BP) {
	iv <- iv |>
		filter(is.finite(POS), !is.na(CHR), is.finite(P))
	if (!nrow(iv))
		return(iv)
	# Reserve the first slot for a genuine cis/local signal, never a trans proxy.
	iv <- iv |>
		mutate(.priority = ifelse(analysis %in% c("cis", "local"), 0L, 1L)) |>
		arrange(.priority, P, SNP)
	keep <- integer()
	for (i in seq_len(nrow(iv))) {
		overlap <- length(keep) && any(iv$CHR[keep] == iv$CHR[i] & abs(iv$POS[keep] - iv$POS[i]) <= 2 * window_bp)
		if (!overlap)
			keep <- c(keep, i)
		if (length(keep) >= max_loci)
			break
	}
	iv[keep, , drop = FALSE] |>
		select( - .priority)
}

le8_locus_class <- function(feature, layer, chr, pos) {
	if (layer == "metabolite") {
		fs <- find_qtl_files(feature, dir.met.gwas, layer)
		if (is.na(fs$joint))
			return("unknown")
		iv <- read_qtl_instruments(feature, dir.met.gwas, layer)$instruments
		if (!nrow(iv))
			return("unknown")
		lead <- iv[which.min(iv$P), , drop = FALSE]
		return(if (chr == lead$CHR && abs(pos - lead$POS) <= le8_num_env("C2_LOCAL_WINDOW_BP", 1e+06)) "local" else "distal")
	}
	a <- layer_annotation(layer, feature)
	if (!nrow(a) || is.na(a$chr[1]) || !is.finite(a$start[1]) || !is.finite(a$end[1]))
		return("unknown")
	pad <- le8_num_env("C2_CIS_WINDOW_BP", 1e+06)
	if (as.character(chr) == as.character(a$chr[1]) && pos >= a$start[1] - pad && pos <= a$end[1] + pad)
		"cis" else "trans"
}

le8_coloc_dataset <- function(d, which, type, file, case_frac = NA_real_) {
	suffix <- if (which == 1L)
		"x" else "y"
	m <- le8_gwas_metadata(file)
	beta <- d[[paste0("BETA_", suffix)]]
	se <- d[[paste0("SE_", suffix)]]
	z <- list(beta = beta, varbeta = se ^ 2, snp = d$SNP, position = d[[paste0("POS_", suffix)]], type = type)
	n <- d[[paste0("N_", suffix)]]
	n <- n[is.finite(n) & n > 0]
	N <- if (length(n))
		median(n) else m$N
	if (is.finite(N) && N > 0)
		z$N <- N
	maf <- pmin(d[[paste0("EAF_", suffix)]], 1 - d[[paste0("EAF_", suffix)]])
	if (all(is.finite(maf) & maf > 0 & maf <= 0.5))
		z$MAF <- maf
	if (type == "cc") {
		# beta/varbeta case-control ABF does not require s. Never substitute UKB incident prevalence for the
		# discovery GWAS case fraction.
		if (is.finite(case_frac) && case_frac > 0 && case_frac < 1)
			z$s <- case_frac
	} else {
		if (is.finite(m$sdY) && m$sdY > 0)
			z$sdY <- m$sdY
		if (is.null(z$sdY) && (is.null(z$N) || is.null(z$MAF)))
			stop("Quantitative coloc needs documented sdY or real N plus MAF: ", file)
	}
	z
}

le8_dense_ld <- function(trait, chr, start, end, snps, EA, NEA) {
	mf <- Sys.getenv("LE8_LD_MANIFEST", unset = "")
	if (!nzchar(mf) || !file.exists(mf))
		return(list(status = "dense LD manifest unavailable"))
	m <- data.table::fread(mf, showProgress = FALSE)
	req <- c("trait", "chr", "start", "end", "file")
	if (!all(req %in% names(m)))
		return(list(status = "LD manifest requires trait,chr,start,end,file"))
	ii <- which(m$trait == trait & as.character(m$chr) == as.character(chr) & m$start <= start & m$end >= end)
	if (!length(ii))
		return(list(status = paste("no covering dense LD for", trait)))
	j <- ii[which.min(m$end[ii] - m$start[ii])]
	f <- as.character(m$file[j])
	z <- tryCatch(readRDS(f), error = function(e) NULL)
	if (!is.list(z) || is.null(z$R) || is.null(z$alleles))
		return(list(status = "LD RDS must contain R and alleles"))
	R <- z$R
	a <- as.data.frame(z$alleles)
	if (!all(c("SNP", "EA", "NEA") %in% names(a)) || anyDuplicated(a$SNP) || !is.matrix(R) || nrow(R) != ncol(R) ||
		is.null(rownames(R)) || !identical(rownames(R), colnames(R)))
		return(list(status = "invalid LD names/alleles"))
	if (!all(snps %in% rownames(R)) || !all(snps %in% a$SNP))
		return(list(status = "incomplete dense LD coverage; no silent SNP thinning"))
	a <- a[match(snps, a$SNP), , drop = FALSE]
	R <- R[snps, snps, drop = FALSE]
	same <- toupper(a$EA) == toupper(EA) & toupper(a$NEA) == toupper(NEA)
	flip <- toupper(a$EA) == toupper(NEA) & toupper(a$NEA) == toupper(EA)
	if (anyNA(same | flip) || any(!(same | flip)))
		return(list(status = "LD/summary allele mismatch"))
	sign <- ifelse(same, 1, - 1)
	R <- R * outer(sign, sign)
	if (any(!is.finite(R)) || max(abs(R - t(R))) > 1e-06 || max(abs(diag(R) - 1)) > 1e-04)
		return(list(status = "invalid dense LD correlations"))
	# chol permits positive semidefinite external LD only after a tiny numeric ridge. No shrinkage is silently
	# used to rescue a materially indefinite R.
	ev <- tryCatch(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values), error = function(e) NA_real_)
	if (!is.finite(ev) || ev <  - 1e-05)
		return(list(status = "dense LD is not positive semidefinite"))
	list(
		status = "ok", R = R, file = f, n_ref = z$n_ref %||% NA_integer_, build = z$build %||% NA_character_,
		ancestry = z$ancestry %||% NA_character_
	)
}

le8_run_susie <- function(dx, dy, d, feature, chr, start, end, locus) {
	no <- list(status = "disabled", summary = tibble(), pip = tibble())
	if (!truthy(Sys.getenv("C3_RUN_SUSIE", unset = "FALSE")))
		return(no)
	if (!requireNamespace("susieR", quietly = TRUE))
		return(modifyList(no, list(status = "susieR missing")))
	if (is.null(dx$N) || is.null(dy$N))
		return(modifyList(no, list(status = "actual discovery N required for each SuSiE trait")))
	if (length(dx$snp) > le8_num_env("C3_SUSIE_MAX_SNPS", 5000))
		return(modifyList(no, list(status = "locus exceeds configured dense-LD memory cap; not thinned")))
	lx <- le8_dense_ld(feature, chr, start, end, d$SNP, d$EA_x, d$NEA_x)
	ly <- le8_dense_ld(Y, chr, start, end, d$SNP, d$EA_x, d$NEA_x) # harmonization orients both betas to EA_x
	if (lx$status != "ok" || ly$status != "ok")
		return(modifyList(no, list(status = paste(lx$status, ly$status, sep = "; "))))
	dx$LD <- lx$R
	dy$LD <- ly$R
	ans <- tryCatch(
		{
			fx <- coloc::runsusie(dx, maxit = 1000, estimate_residual_variance = FALSE)
			fy <- coloc::runsusie(dy, maxit = 1000, estimate_residual_variance = FALSE)
			if (!isTRUE(fx$converged) || !isTRUE(fy$converged))
				stop("SuSiE did not converge")
			fits <- lapply(C3_P12, function(p12) coloc::coloc.susie(fx, fy, p1 = 1e-04, p2 = 1e-04, p12 = p12))
			sm <- imap_dfr(fits, function(ff, nm) as_tibble(ff$summary) |>
				mutate(prior = nm, feature = feature, locus = locus))
			pp <- tibble(feature, locus, SNP = d$SNP, QTL_PIP = fx$pip, outcome_PIP = fy$pip)
			rd <- .le8_analysis_state$rawdir
			dir.create(file.path(rd, "susie"), showWarnings = FALSE)
			stem <- paste0(gsub("[^A-Za-z0-9_.-]", "_", paste(feature, locus)), ".rds")
			saveRDS(
				list(QTL = fx, outcome = fy, coloc = fits, summary = sm, pip = pp, LD_QTL = lx$file, LD_outcome = ly$file),
				file.path(rd, "susie", stem)
			)
			list(status = "ok; signal-pair posterior available", summary = sm, pip = pp)
		},
		error = function(e) modifyList(no, list(status = conditionMessage(e)))
	)
	ans
}

coloc_one_locus <- function(
	feature, layer, qtl_file, lead_chr, lead_pos, ygwas_file, outcome_type, case_frac,
	locus_class = NULL
) {
	start <- max(1, lead_pos - WINDOW_BP)
	end <- lead_pos + WINDOW_BP
	locus <- paste0("chr", lead_chr, ":", floor(start), "-", ceiling(end))
	cls <- if (is.null(locus_class))
		le8_locus_class(feature, layer, lead_chr, lead_pos) else locus_class
	s <- tibble(layer, feature, locus,
		chr = as.character(lead_chr), start, end, lead_pos, locus_class = cls, status = "failed",
		n_snps = 0L, PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_, PP.H3 = NA_real_, PP.H4 = NA_real_,
		PP.H4_p12_conservative = NA_real_, PP.H4_p12_default = NA_real_, PP.H4_p12_liberal = NA_real_, PP.H4_robust_min = NA_real_,
		lead_shared = NA_character_, lead_shared_pp = NA_real_, credible_set_n = NA_integer_, message = "", credible_set_interpretation = "SNP posterior CONDITIONAL ON H4; not a trait fine-mapping credible set",
		susie_status = "not attempted", case_fraction_source = if (is.finite(case_frac))
			"configured discovery GWAS" else "not required for beta/varbeta cc ABF"
	)
	fail <- function(msg) {
		s$message <- msg
		list(summary = s, variants = tibble(), regional = tibble())
	}
	q <- read_sumstat_region(qtl_file, lead_chr, start, end)
	y <- read_sumstat_region(ygwas_file, lead_chr, start, end)
	d <- harmonize_sumstats(q, y)
	if (!nrow(d))
		return(fail("no allele-aligned overlapping SNPs"))
	d <- d |>
		filter(is.finite(BETA_x), is.finite(SE_x), SE_x > 0, is.finite(BETA_y), is.finite(SE_y), SE_y > 0)
	s$n_snps <- nrow(d)
	if (nrow(d) < MIN_SNPS)
		return(fail("insufficient dense overlap"))
	mq <- le8_gwas_metadata(qtl_file)
	my <- le8_gwas_metadata(ygwas_file)
	if (!is.na(mq$build) && !is.na(my$build) && as.character(mq$build) != as.character(my$build))
		return(fail("declared discovery builds differ; harmonize coordinates before C3"))
	datasets <- tryCatch(list(x = le8_coloc_dataset(d, 1L, "quant", qtl_file), y = le8_coloc_dataset(
		d, 2L, outcome_type,
		ygwas_file, case_frac
	)), error = function(e) e)
	if (inherits(datasets, "condition"))
		return(fail(conditionMessage(datasets)))
	fits <- lapply(C3_P12, function(p12) tryCatch(
		{
			f <- NULL
			invisible(capture.output(f <- coloc::coloc.abf(datasets$x, datasets$y, p1 = 1e-04, p2 = 1e-04, p12 = p12)))
			f
		},
		error = function(e) e
	))
	f <- fits[["default"]]
	if (inherits(f, "condition"))
		return(fail(conditionMessage(f)))
	for (i in 0 : 4) s[[paste0("PP.H", i)]] <- as.numeric(f$summary[[paste0("PP.H", i, ".abf")]])
	pp <- vapply(fits, function(x) if (inherits(x, "condition"))
		NA_real_ else as.numeric(x$summary[["PP.H4.abf"]]), numeric(1))
	s$PP.H4_p12_conservative <- pp["conservative"]
	s$PP.H4_p12_default <- pp["default"]
	s$PP.H4_p12_liberal <- pp["liberal"]
	s$PP.H4_robust_min <- if (all(is.finite(pp)))
		min(pp) else NA_real_ # missing prior is not a pass
	v <- as_tibble(f$results) |>
		arrange(desc(SNP.PP.H4)) |>
		mutate(
			cum_pp = cumsum(SNP.PP.H4), credible95 = lag(cum_pp, default = 0) < 0.95, feature = feature, layer = layer,
			locus = locus, interpretation = "Conditional shared-variant posterior under H4"
		)
	s$lead_shared <- as.character(v$snp[1])
	s$lead_shared_pp <- v$SNP.PP.H4[1]
	s$credible_set_n <- sum(v$credible95)
	s$status <- "ok"
	su <- le8_run_susie(datasets$x, datasets$y, d, feature, lead_chr, start, end, locus)
	s$susie_status <- su$status
	r <- bind_rows(
		d |>
			transmute(layer, feature, locus, dataset = "QTL", SNP, position = POS_x, p = P_x, beta = BETA_x, se = SE_x),
		d |>
			transmute(layer, feature, locus, dataset = Y, SNP, position = POS_y, p = P_y, beta = BETA_y, se = SE_y)
	)
	list(summary = s, variants = v, regional = r, susie = su)
}

le8_c3_sets <- function(res, mr, layer) {
	e <- le8_same_locus_evidence(mr, res, layer)
	g <- unique(e$feature[e$eligible %in% TRUE])
	# Compatibility keys retained; names are NOT a claim of identified causation.
	list(Causal_Tier1 = character(), Causal_Tier2plus = g, Causal_any = g)
}

le8_c3_additions <- function(res, mr, layer, outdir) {
	e <- le8_same_locus_evidence(mr, res, layer)
	rd <- le8_job_dir(outdir, "c3_coloc")
	write_raw_csv(e, "c3.same_locus_evidence.csv", rd)
	ss <- list.files(file.path(rd, "susie"), pattern = "\\.rds$", full.names = TRUE)
	su <- bind_rows(lapply(ss, function(f) {
		z <- readRDS(f)
		as_tibble(z$summary)
	}))
	pi <- bind_rows(lapply(ss, function(f) {
		z <- readRDS(f)
		as_tibble(z$pip)
	}))
	write_raw_csv(su, "c3.susie_signal_pairs.csv", rd)
	write_raw_csv(pi, "c3.susie_trait_PIP.csv", rd)
	status <- if ("susie_status" %in% names(res))
		res |>
			count(susie_status) else tibble(status = "no loci")
	list(same_locus = e, susie_pairs = su, susie_PIP = pi, susie_status = status)
}

run_c3_layer <- function(layer = c("protein", "metabolite")) {
	layer <- match.arg(layer)

	# Cell-specific analyses are executed only by c5_cellulation.
	if (LE8_REUSE_RESULTS) {
		restored <- le8_restore_outputs(layer, "c3_coloc")
		od <- if (layer == "protein")
			out.prot else out.met
		cf <- file.path(od, "c3_coloc", "c3.coloc_summary.csv")
		co <- if (file.exists(cf))
			as.data.frame(data.table::fread(cf)) else data.frame()
		tri <- read_c3_pgs_integration(layer, od, co)
		plot_c3_pgs_integration(tri, od)
		return(restored)
	}
	# Check reusable results and initialize the analysis output directory.
	layer <- match.arg(layer)
	le8_begin_analysis(layer, "c3_coloc")
	.le8_analysis_env <- environment()
	on.exit(le8_finish_analysis(layer, "c3_coloc", .le8_analysis_env), add = TRUE)

	layer <- match.arg(layer)
	outdir <- if (layer == "protein")
		out.prot else out.met
	setwd2(outdir)
	rawdir <- le8_job_dir(outdir, LE8_JOB)
	dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	cache <- file.path(rawdir, "c3.res.rds")
	if (!is.finite(C3_MR_FDR) || C3_MR_FDR <= 0 || C3_MR_FDR >= 1 || !is.finite(MAX_FEATURES) || MAX_FEATURES <
		1 || !is.finite(MAX_LOCI_PER_FEATURE) || MAX_LOCI_PER_FEATURE < 1)
		stop("Invalid C3 MR FDR / feature / locus limits", call. = FALSE)
	c2f <- file.path(le8_job_dir(outdir, "c2_cause"), "c2.res.rds")
	if (!file.exists(c2f))
		stop("C3 requires completed C2 MR results; run C2 first.", call. = FALSE)
	c2 <- readRDS(c2f)
	mr <- c2$MR %||% tibble()
	selected_mr <- c3_select_mr(mr, layer, C3_MR_FDR, MAX_FEATURES)
	candidates <- unique(selected_mr$exposure)
	message(
		"C3/", layer, ": FDR_all < ", C3_MR_FDR, "; ", nrow(selected_mr), " significant MR analyses, ", length(candidates),
		" features; up to ", MAX_LOCI_PER_FEATURE, " loci per feature"
	)
	base <- if (layer == "protein")
		dir.X else dir.met.gwas
	ygfile <- get_y_gwas_file(Y, TRUE)
	qfiles <- setNames(lapply(candidates, function(feature) find_qtl_files(feature, base, layer)), candidates)
	input_files <- c(ygfile, unlist(lapply(qfiles, function(fs) fs$full)))
	selection_signature <- le8_hash_object(list(
		version = C3_SELECTION_VERSION, mr = as.data.frame(selected_mr[
,
			c("exposure", "analysis", "pval", "FDR_all", "instrument_snps", "instrument_positions")
		]), settings = c3_locus_settings(),
		files = c3_file_stamp(input_files), max_loci = MAX_LOCI_PER_FEATURE, fdr = C3_MR_FDR, H4 = H4_STRONG, case_frac = Sys.getenv(
			"COLOC_CASE_FRAC",
			"NA"
		), gpu = Sys.getenv("RUN_GPU_COLOC", "TRUE"), gpu_p12 = Sys.getenv("GPU_COLOC_P12", "1e-5")
	))
	gpudir <- file.path(le8_cache_dir("gpu_coloc", layer), "runs", selection_signature)
	qtl_cache_root <- Sys.getenv("C3_QTL_CACHE_DIR", unset = file.path(le8_cache_root(), "c3_qtl"))
	write_raw_csv(selected_mr, "c3.selected_mr.csv", rawdir)
	if (cache_valid(cache)) {
		old <- tryCatch(readRDS(cache), error = function(e) NULL)
		if (!is.null(old) && identical(old$meta$selection_signature, selection_signature) && all(c(
			"summary", "regional",
			"variants"
		) %in% names(old))) {
			message("C3/", layer, ": reuse locus results and regenerate ", C3_CODE_VERSION, " figures")
			plot_coloc_results(old$summary, old$regional, old$variants, layer, outdir)
			gpu <- old$GPU_coloc %||% read_gpu_coloc_results(rawdir, gpudir)
			plot_gpu_coloc_validation(gpu, old$summary, outdir)
			aud <- credible_set_audit(old$summary, old$variants)
			tri <- read_c3_pgs_integration(layer, outdir, old$summary)
			plot_c3_pgs_integration(tri, outdir)
			write_raw_csv(aud$overall, "c3.credible_set_audit.csv", rawdir)
			write_raw_csv(aud$by_locus, "c3.credible_set_by_locus.csv", rawdir)
			write_raw_csv(tri, "c3.pgs_observed_coloc_triangulation.csv", rawdir)
			old$credible_set_audit <- aud
			old$pgs_triangulation <- tri
			old$meta$code_version <- C3_CODE_VERSION
			lists <- old$causal_lists %||% list()
			gpu <- old$GPU_coloc %||% list(results = tibble(), status = tibble())
			write_xlsx2(list(
				coloc_summary = old$summary, credible_set_audit = aud$overall, credible_set_by_locus = aud$by_locus,
				GPU_results = gpu$results %||% tibble(),
				GPU_status = gpu$status %||% tibble(), GPU_manifest = old$manifest %||% tibble(), causal_sets = if (length(lists)) stack(lists) else tibble(),
				pgs_observed_coloc = tri
			), "c3.out.xlsx")
			saveRDS(old, cache, compress = "xz")
			finalize_outputs(LE8_JOB, outdir)
			return(old)
		}
	}
	outtype <- infer_outcome_type(Y)
	sfrac <- if (outtype == "cc")
		get_case_fraction(Y) else NA_real_
	if (outtype == "cc" && (!is.finite(sfrac) || sfrac <= 0 || sfrac >= 1))
		message("C3: discovery case fraction not supplied; beta/varbeta cc ABF omits s, rather than substituting cohort prevalence")
	coloc_started <- le8_stage_start(paste0("C3/", layer, " coloc"))
	rows <- list()
	vrows <- list()
	rrows <- list()
	manifest <- list()
	k <- 0L
	reused_loci <- 0L
	legacy_files <- list.files(le8_cache_dir("coloc_loci", layer), pattern = "[.]rds$", full.names = TRUE)
	locus_cache_dir <- file.path(le8_cache_dir("coloc_loci", layer), "mr_selected_v1")
	dir.create(locus_cache_dir, recursive = TRUE, showWarnings = FALSE)
	selection_audit <- list()
	write_raw_csv(tibble(
		SNP = character(), analysis = character(), mr_position = character(), feature = character(),
		selected_lead = logical()
	), "c3.mr_locus_selection.csv", rawdir)
	for (feature in candidates) {
		feature_started <- Sys.time()
		message(
			"C3/", layer, ": feature ", match(feature, candidates), "/", length(candidates), " ", feature,
			"; completed loci=", k, ", reused=", reused_loci
		)
		qf <- qfiles[[feature]]$full
		if (length(qf) != 1L || is.na(qf) || !file.exists(qf))
			stop("Missing full QTL for MR-selected feature: ", feature, call. = FALSE)
		fm <- selected_mr |>
			filter(exposure == feature)
		iv <- c3_read_mr_qtl(fm, qf, qtl_cache_root)
		leads <- select_nonoverlapping_leads(iv, MAX_LOCI_PER_FEATURE, WINDOW_BP)
		selection_audit[[feature]] <- iv |>
			mutate(feature = feature, selected_lead = SNP %in% leads$SNP)
		write_raw_csv(bind_rows(selection_audit), "c3.mr_locus_selection.csv", rawdir)
		for (i in seq_len(nrow(leads))) {
			k <- k + 1
			locus_key <- le8_hash_object(list(
				version = C3_SELECTION_VERSION, feature = feature, layer = layer,
				chr = leads$CHR[i], pos = leads$POS[i], class = leads$analysis[i], files = c3_file_stamp(c(
					qf,
					ygfile
				)), settings = c3_locus_settings(), outtype = outtype, sfrac = sfrac
			))
			locus_cache <- file.path(locus_cache_dir, paste0(locus_key, ".rds"))
			z <- c3_cached_locus(
				locus_cache, legacy_files, feature, layer, leads$CHR[i], leads$POS[i], leads$analysis[i],
				qf, ygfile, outtype, sfrac
			)
			reused <- is.list(z) && all(c("summary", "variants", "regional") %in% names(z))
			message(
				"  locus ", i, "/", nrow(leads), " chr", leads$CHR[i], ":", leads$POS[i], " [", leads$analysis[i],
				"] ", if (reused)
					"cache hit" else "computing"
			)
			if (!reused) {
				z <- coloc_one_locus(feature, layer, qf, leads$CHR[i], leads$POS[i], ygfile, outtype, sfrac, locus_class = leads$analysis[i])
				write_stage_cache(z, locus_cache)
			} else reused_loci <- reused_loci + 1L
			rows[[k]] <- z$summary
			vrows[[k]] <- z$variants
			rrows[[k]] <- z$regional
			manifest[[k]] <- tibble(omics = layer, trait = feature, file = qf, type = "quant", region = z$summary$locus)
		}
		message("C3/", layer, ": ", feature, " done in ", round(as.numeric(difftime(Sys.time(), feature_started,
			units = "secs"
		)), 1), " s; loci=", k, ", reused=", reused_loci)
	}
	res <- bind_rows(rows)
	variants <- bind_rows(vrows)
	regional <- bind_rows(rrows)
	mani <- bind_rows(manifest)
	le8_stage_done(paste0("C3/", layer, " coloc"), coloc_started, paste0("loci=", nrow(res)))
	if (!nrow(res)) {
		message("C3/", layer, ": no QTL locus could be constructed; write auditable empty outputs and blank panels")
		res <- tibble(
			layer = character(), feature = character(), locus = character(), chr = character(), start = numeric(),
			end = numeric(), lead_pos = numeric(), status = character(), n_snps = integer(), PP.H0 = numeric(),
			PP.H1 = numeric(), PP.H2 = numeric(), PP.H3 = numeric(), PP.H4 = numeric(), PP.H4_p12_conservative = numeric(),
			PP.H4_p12_default = numeric(), PP.H4_p12_liberal = numeric(), PP.H4_robust_min = numeric(), lead_shared = character(),
			lead_shared_pp = numeric(), credible_set_n = integer(), message = character(), PP4_rank_fraction = numeric(),
			tier = character()
		)
		aud <- credible_set_audit(res, variants)
		write_raw_csv(res, "c3.coloc_summary.csv", rawdir)
		write_raw_csv(variants, "c3.variant_posteriors.csv", rawdir)
		write_raw_csv(regional, "c3.regional_rows.csv", rawdir)
		write_raw_tsv(mani, "qtl_cad_manifest.tsv", rawdir)
		write_raw_csv(aud$overall, "c3.credible_set_audit.csv", rawdir)
		write_raw_csv(aud$by_locus, "c3.credible_set_by_locus.csv", rawdir)
		plot_coloc_results(res, regional, variants, layer, outdir)
		gpu <- list(results = tibble(), status = tibble())
		plot_gpu_coloc_validation(gpu, res, outdir)
		tri <- read_c3_pgs_integration(layer, outdir, res)
		plot_c3_pgs_integration(tri, outdir)
		write_raw_csv(tri, "c3.pgs_observed_coloc_triangulation.csv", rawdir)
		lists <- list(Causal_Tier1 = character(), Causal_Tier2plus = character(), Causal_any = character())
		saveRDS(lists, file.path(rawdir, paste0("c3.causal_", layer, "_lists.rds")), compress = "xz")
		out <- list(
			meta = module_meta(layer, extra = list(
				outcome_type = outtype, case_fraction = sfrac, code_version = C3_CODE_VERSION,
				selection_signature = selection_signature, status = "no QTL locus constructed"
			)), summary = res, variants = variants,
			regional = regional, manifest = mani, GPU_coloc = gpu, causal_lists = lists, credible_set_audit = aud,
			pgs_triangulation = tri
		)
		saveRDS(out, cache, compress = "xz")
		write_xlsx2(list(
			coloc_summary = res, credible_set_audit = aud$overall, credible_set_by_locus = aud$by_locus,
			GPU_results = gpu$results, GPU_status = gpu$status,
			GPU_manifest = mani, causal_sets = stack(lists), pgs_observed_coloc = tri
		), "c3.out.xlsx")
		finalize_outputs(LE8_JOB, outdir)
		return(out)
	}
	res <- res |>
		mutate(PP4_rank_fraction = ifelse(status == "ok" & sum(status == "ok") > 0, rank( - PP.H4,
			ties.method = "min",
			na.last = "keep"
		) / sum(status == "ok"), NA_real_), tier = case_when(status != "ok" ~ "Not tested", PP.H4 >=
			0.8 ~ "Tier 1", PP.H4 >= H4_STRONG ~ "Tier 2", PP.H4 >= 0.5 ~ "Tier 3", TRUE ~ "Tier 4"))
	aud <- credible_set_audit(res, variants)
	write_raw_csv(res, "c3.coloc_summary.csv", rawdir)
	write_raw_csv(variants, "c3.variant_posteriors.csv", rawdir)
	write_raw_csv(regional, "c3.regional_rows.csv", rawdir)
	write_raw_tsv(mani, "qtl_cad_manifest.tsv", rawdir)
	write_raw_csv(aud$overall, "c3.credible_set_audit.csv", rawdir)
	write_raw_csv(aud$by_locus, "c3.credible_set_by_locus.csv", rawdir)
	le8_stage(
		paste0("C3/", layer, " GPU-coloc"), run_gpu_coloc_step(layer, rawdir, mani, ygfile, outtype, gpudir),
		paste0("loci=", nrow(mani))
	)
	plot_coloc_results(res, regional, variants, layer, outdir)
	gpu <- read_gpu_coloc_results(rawdir, gpudir)
	plot_gpu_coloc_validation(gpu, res, outdir)
	tri <- read_c3_pgs_integration(layer, outdir, res)
	plot_c3_pgs_integration(tri, outdir)
	write_raw_csv(tri, "c3.pgs_observed_coloc_triangulation.csv", rawdir)
	lists <- le8_c3_sets(res, mr, layer)
	saveRDS(lists, file.path(rawdir, paste0("c3.causal_", layer, "_lists.rds")), compress = "xz")
	if (layer == "protein")
		saveRDS(lists, file.path(rawdir, "c3.causal_protein_lists.rds"), compress = "xz")
	out <- list(
		meta = module_meta(layer, extra = list(
			outcome_type = outtype, case_fraction = sfrac, code_version = C3_CODE_VERSION,
			selection_signature = selection_signature
		)), summary = res, variants = variants, regional = regional, manifest = mani,
		GPU_coloc = gpu, causal_lists = lists, credible_set_audit = aud, pgs_triangulation = tri
	)
	saveRDS(out, cache, compress = "xz")
	write_xlsx2(list(
		coloc_summary = res, credible_set_audit = aud$overall, credible_set_by_locus = aud$by_locus,
		GPU_results = gpu$results, GPU_status = gpu$status,
		GPU_manifest = mani, causal_sets = stack(lists), pgs_observed_coloc = tri
	), "c3.out.xlsx")
	finalize_outputs(LE8_JOB, outdir)
	out
}

if (prot_DO) invisible(le8_stage("C3/protein", run_c3_layer("protein")))
if (met_DO) invisible(le8_stage("C3/metabolite", run_c3_layer("metabolite")))
