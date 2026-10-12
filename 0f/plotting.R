# 🚩 Shared plotting themes
leg_theme <- ggplot2::theme(legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.65), colour = "grey60"), legend.text = ggplot2::element_text(size = 8), legend.key.size = grid::unit(0.45, "lines"))

fig_theme <- function(base_size = 16) {
	ggplot2::theme_classic(base_size = base_size) +
		ggplot2::theme(
			plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.18, hjust = .5),
			plot.subtitle = ggplot2::element_text(size = base_size * .88, hjust = .5, margin = ggplot2::margin(b = 8)),
			axis.title = ggplot2::element_text(face = "bold"),
			axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
			axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
			axis.text = ggplot2::element_text(face = "bold", color = "black"),
			legend.title = ggplot2::element_text(face = "bold"),
			legend.text = ggplot2::element_text(face = "bold"),
			panel.grid.major = ggplot2::element_line(color = "grey88", linewidth = .45),
			panel.grid.minor = ggplot2::element_blank(),
			plot.margin = ggplot2::margin(12, 18, 12, 18)
		)
}


# 🚩 Forest plots
plot_forest <- function(dat1, lab.X = NA, x = "estimate", xmin = "conf.low", xmax = "conf.high", lab.Y = "lab.Y", col = "color", n_col = "N_event", ref = 1, point_size = 3, line_width = 0.9, err_height = 0.25, text_nudge = 0.20, text_size = 3.0, base_size = 13) {
	stopifnot(all(c(x, xmin, xmax, lab.Y, col, n_col) %in% names(dat1)))
	if (length(lab.X) > 1) {
		plots <- lapply(lab.X, function(exp_name) {
			sub_dat <- dat1[dat1$Exposure == exp_name, ]
			if (nrow(sub_dat) == 0) return(NULL)
			plot_forest(sub_dat, lab.X = exp_name, x = x, xmin = xmin, xmax = xmax, lab.Y = lab.Y, col = col, n_col = n_col, ref = ref, point_size = point_size, line_width = line_width, err_height = err_height, text_nudge = text_nudge, text_size = text_size, base_size = base_size)
		})
		names(plots) <- lab.X
		return(plots)
	}
	if (!is.factor(dat1[[lab.Y]])) dat1[[lab.Y]] <- factor(dat1[[lab.Y]])
	p <- ggplot(dat1, aes(x = .data[[x]], y = as.numeric(.data[[lab.Y]]))) + geom_vline(xintercept = ref, linetype = "dashed", colour = "grey55") +
		geom_errorbar(aes(xmin = .data[[xmin]], xmax = .data[[xmax]]), width = err_height, colour = dat1[[col]], linewidth = line_width, orientation = "y") +
		geom_point(colour = dat1[[col]], size = point_size) +
		 geom_text(aes(label = sprintf("(N=%s, HR=%.2f)", format(.data[[n_col]], big.mark = ",", trim = TRUE), .data[[x]])), y = as.numeric(dat1[[lab.Y]]) - text_nudge, colour = dat1[[col]], size = text_size, vjust = 1) +
		scale_y_continuous(breaks = seq_along(levels(dat1[[lab.Y]])), labels = levels(dat1[[lab.Y]])) +
		labs(x = if (any(is.na(lab.X))) "Hazard ratio" else paste0("Hazard ratio (", lab.X, ")"), y = NULL) +
		theme_classic(base_size = base_size) + theme(axis.text.y = element_text(size = 11), plot.margin = margin(6, 10, 14, 10))
	return(p)
}


# 🚩 Phenome-wide association plots
plot_phewas <- function(dat1, phecode, Xs, varX, output_dir = NULL, resume = TRUE, make_plots = TRUE,
		cores = suppressWarnings(as.integer(Sys.getenv("PHEWAS_CORES", unset = "1")))) {
	if (length(cores) != 1L || is.na(cores) || cores < 1L) {
		warning("Invalid PHEWAS_CORES; using 1.", call. = FALSE)
		cores <- 1L
	}
	message("Running PheWAS with ", cores, " core(s), one genetic variable at a time")
	phenos <- if(length(phecode) == 1 && is.na(phecode)) readRDS(paste0(indir, "/Rdata/PheWAS.rds")) else phecode
	cov <- dplyr::select(dat1, id = eid, all_of(varX))
	checkpoint_dir <- if (is.null(output_dir)) NULL else file.path(output_dir, "phewas_checkpoints")
	if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
	results <- vector("list", length(Xs)); names(results) <- Xs
	plots <- vector("list", length(Xs)); names(plots) <- Xs
	for (i in seq_along(Xs)) {
		x <- Xs[[i]]
		checkpoint <- if (is.null(checkpoint_dir)) NULL else file.path(checkpoint_dir, paste0(x, ".phewas.rds"))
		png_file <- if (is.null(output_dir)) NULL else file.path(output_dir, paste0(x, ".phewas.png"))
		if (isTRUE(resume) && !is.null(checkpoint) && file.exists(checkpoint)) {
			message(sprintf("[%d/%d] Resume %s from checkpoint", i, length(Xs), x))
			res_x <- readRDS(checkpoint)
		} else {
			message(sprintf("[%d/%d] START %s", i, length(Xs), x))
			res_x <- phewas(
				phenotypes = phenos,
				genotypes = dplyr::select(dat1, id = eid, all_of(x)),
				covariates = cov, cores = cores, significance.threshold = "bonferroni"
			) %>% addPhecodeInfo() %>% rename(category = group) %>% mutate(p = thinP1(p)) %>% as.data.frame()
			if (!is.null(checkpoint)) {
				tmp_checkpoint <- paste0(checkpoint, ".tmp")
				saveRDS(res_x, tmp_checkpoint)
				if (!file.rename(tmp_checkpoint, checkpoint)) stop("Could not write checkpoint: ", checkpoint)
			}
		}
		dat2 <- res_x[res_x$snp == x & is.finite(res_x$p), ]
		p <- if (!isTRUE(make_plots) || !nrow(dat2)) NULL else suppressMessages(phewasManhattan(dat2,
			annotate.phenotype.description = dat2[, c("phenotype", "description")], title = x, OR.direction = FALSE, size.x.labels = 14, size.y.labels = 14)
		) + theme(text = element_text(size = 12))
		if (!is.null(png_file) && !is.null(p) && (!isTRUE(resume) || !file.exists(png_file))) {
			ggsave(png_file, p, width = 14, height = 8, dpi = 300, bg = "white", limitsize = FALSE)
		}
		message(sprintf("[%d/%d] DONE %s%s", i, length(Xs), x,
			if (!isTRUE(make_plots) || is.null(png_file) || !file.exists(png_file)) "" else paste0(" -> ", png_file)))
		results[[i]] <- res_x
		plots[i] <- list(if (is.null(output_dir)) p else NULL)
		rm(res_x, dat2, p); invisible(gc())
	}
	list(plots = plots, res = dplyr::bind_rows(results))
}

# Align panels within each row while keeping row-specific y-axis widths.
align_panel_rows <- function(rows, rel_widths = NULL, rel_heights = NULL,
		align = "h", axis = "tb", row_gap = 0.10,
		axis_text_size = 10.5, axis_title_size = 11.5, side_pad = 0.025) {
	if (!requireNamespace("cowplot", quietly = TRUE)) {
		stop("align_panel_rows() requires the cowplot package.", call. = FALSE)
	}
	if (!is.list(rows) || !length(rows) || !all(vapply(rows, is.list, logical(1)))) {
		stop("rows must be a non-empty list of plot lists.", call. = FALSE)
	}
	widths_by_row <- if (is.null(rel_widths)) {
		rep(list(NULL), length(rows))
	} else if (is.list(rel_widths)) {
		if (length(rel_widths) != length(rows)) stop("rel_widths list must match rows.", call. = FALSE)
		rel_widths
	} else {
		rep(list(rel_widths), length(rows))
	}
	# Apply one typography contract before grob construction.  Keeping this here
	# prevents individual figures from accumulating contradictory theme margins.
	rows <- lapply(rows, function(x) lapply(x, function(p) {
		if (inherits(p, "ggplot")) {
			p + ggplot2::theme(
				axis.title = ggplot2::element_text(face = "bold", size = axis_title_size),
				axis.text = ggplot2::element_text(face = "bold", size = axis_text_size, color = "black")
			)
		} else p
	}))
	row_grobs <- Map(function(x, w) {
		args <- list(plotlist = x, nrow = 1, align = align, axis = axis)
		if (!is.null(w)) args$rel_widths <- w
		do.call(cowplot::plot_grid, args)
	}, rows, widths_by_row)
	if (is.finite(side_pad) && side_pad > 0 && side_pad < .25) {
		row_grobs <- lapply(row_grobs, function(g) cowplot::plot_grid(
			NULL, g, NULL, nrow = 1, rel_widths = c(side_pad, 1 - 2 * side_pad, side_pad)
		))
	}
	# Insert actual grid rows between plot rows.  Unlike plot.margin, this spacing
	# is independent of axis-label length and therefore remains stable.
	if (length(row_grobs) > 1L && is.finite(row_gap) && row_gap > 0) {
		with_gaps <- vector("list", 2L * length(row_grobs) - 1L)
		with_gaps[seq(1L, length(with_gaps), 2L)] <- row_grobs
		with_gaps[seq(2L, length(with_gaps), 2L)] <- rep(list(NULL), length(row_grobs) - 1L)
		base_h <- if (is.null(rel_heights)) rep(1, length(row_grobs)) else rel_heights
		gap_h <- row_gap * mean(base_h)
		out_h <- as.numeric(rbind(base_h[ - length(base_h)], rep(gap_h, length(base_h) - 1L)))
		out_h <- c(out_h, base_h[length(base_h)])
		row_grobs <- with_gaps
		rel_heights <- out_h
	}
	args <- list(plotlist = row_grobs, ncol = 1, align = "none")
	if (!is.null(rel_heights)) args$rel_heights <- rel_heights
	do.call(cowplot::plot_grid, args)
}

# Volcano plot.


# 🚩 Volcano and protein association plots
plot_volcano <- function(label, dat, X, BETA, P, sig, label_x, label_y, topN = NA, topP = NA) {
	dat$X <- dat[[X]]; dat$BETA <- dat[[BETA]]; dat$P <- dat[[P]]
	dat <- dat %>% mutate(color = ifelse(P < sig & BETA > 0, "positive", ifelse(P < sig & BETA < 0, "negative", "NS")))
	if (!is.na(topN)) top <- dat %>% arrange(P) %>% head(topN)
	if (!is.na(topP)) top <- dat %>% filter(P <= topP)
	top$X <- top$X
	top$nudge_x <-  - 0.6 * top$BETA # Pull labels toward the center at x = 0.
	ggplot(dat, aes(x = BETA, y =  - log10(P), color = color)) + geom_point(size = 2) +
		scale_color_manual(values = c(positive = "purple", negative = "green", NS = "gray")) +
		ggrepel::geom_text_repel(data = top, aes(label = X),
			nudge_x = top$nudge_x, direction = "y", # Push labels toward zero first, then repel mainly along the y direction.
			color = "darkblue", fontface = "bold", size = 3, box.padding = 0.35, point.padding = 0.25, min.segment.length = 0, segment.color = "black", segment.size = 0.3,
			arrow = grid::arrow(length = grid::unit(0.015, "npc"), type = "open"), max.overlaps = Inf, seed = 1) +
		geom_hline(yintercept =  - log10(sig), linetype = "dashed", color = "red", linewidth = 1.2) +
		geom_vline(xintercept = 0, linetype = "dotted", linewidth = 1.2) +
		labs(x = label_x, y = label_y, title = label) + theme_minimal() +
		theme(axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 14, face = "bold"), axis.line = element_line(linewidth = 1.2),
			  legend.position = "none", plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
}


# Genetic-correlation heatmap.


# 🚩 Genetic correlations
plot_rg <- function(dat1, p1_col, p2_col, rg_col, p_col, alpha = 0.05, cluster_abs = TRUE, diag_value = 1, show_diag_text = FALSE, ns_label = "ns"){
	pal_rg <- c("#d73027", "#fc8d59", "#91cf60", "#1a9850") # Red -> orange -> green, bright palette.
	col_ns <- "grey80"; col_diag <- "skyblue2"; col_na <- "white"
	dat <- dat1 %>% dplyr::mutate(p1 = as.character(.data[[p1_col]]), p2 = as.character(.data[[p2_col]]),
		rg = pmax( - 1, pmin(1, as.numeric(.data[[rg_col]]))), p = as.numeric(.data[[p_col]]),
		a = pmin(p1, p2), b = pmax(p1, p2)) %>% dplyr::group_by(a, b) %>%
		summarise(rg = {x <- rg[!is.na(rg)]; if(!length(x)) NA_real_ else mean(x)}, p = {x <- p[!is.na(p)]; if(!length(x)) NA_real_ else min(x)}, .groups = "drop")
	traits <- sort(unique(c(dat$a, dat$b))); n <- length(traits)
	mat_rg <- matrix(NA_real_, n, n, dimnames = list(traits, traits)); mat_p <- mat_rg
	ij <- cbind(match(dat$a, traits), match(dat$b, traits)); mat_rg[ij] <- dat$rg; mat_p[ij] <- dat$p
	mat_rg <- pmax(mat_rg, t(mat_rg), na.rm = TRUE); mat_p <- pmin(mat_p, t(mat_p), na.rm = TRUE); mat_p[is.infinite(mat_p)] <- NA_real_
	diag(mat_rg) <- diag_value; diag(mat_p) <- 0
	mat_clust <- mat_rg; mat_clust[is.na(mat_clust)] <- 0; dist_mat <- if(cluster_abs) stats::as.dist(1 - abs(mat_clust)) else stats::as.dist(1 - mat_clust)
	sig <- !is.na(mat_rg) & !is.na(mat_p) & mat_p <= alpha; ns <- !is.na(mat_rg) & !is.na(mat_p) & mat_p > alpha; diagm <- row(mat_rg) == col(mat_rg)
	na_cell <- (is.na(mat_rg) | is.na(mat_p)) & !diagm
	lab <- matrix("", n, n, dimnames = dimnames(mat_rg)); lab[sig] <- sprintf("%.2f", mat_rg[sig]); lab[ns] <- ns_label; lab[na_cell] <- "NA"
	diag(lab) <- if(show_diag_text) sprintf("%.2f", diag_value) else ""
	rg_has_neg <- any(mat_rg[sig] < 0, na.rm = TRUE); rg_min <- if(rg_has_neg) - 1 else 0; rg_max <- 1
	col_fun <- circlize::colorRamp2(seq(rg_min, rg_max, length.out = length(pal_rg)), pal_rg)
	ht <- ComplexHeatmap::Heatmap(mat_rg, name = "rg", col = col_fun, na_col = col_na,
		cluster_rows = stats::hclust(dist_mat, method = "average"), cluster_columns = stats::hclust(dist_mat, method = "average"),
		show_row_names = TRUE, show_column_names = TRUE,
		row_names_gp = grid::gpar(fontsize = 14), column_names_gp = grid::gpar(fontsize = 14),
		cell_fun = function(i, j, x, y, w, h, fill){
			if(diagm[j, i]) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = col_diag, col = NA))
			else if(na_cell[j, i]) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = col_na, col = NA))
			else if(ns[j, i]) grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = col_ns, col = NA))
			if(lab[j, i] != "") grid::grid.text(lab[j, i], x, y, gp = grid::gpar(fontsize = 12))
	})
	ComplexHeatmap::draw(ht)
	invisible(list(mat_rg = mat_rg, mat_p = mat_p, labels = lab, traits = traits, dat_sym = dat))
}


# 🚩 Trends and maps


plot_map <- function(dat1, centers = NULL, separate = TRUE, width = 50, nrow = 1,
			dist_outlier = 0.98, transparency = "tdi", center.lst = center.lst,
			cell_size = 1000, alpha_range = c(0.18, 0.95),
			basemap = TRUE, basemap_type = "osm", tile_zoom = NULL,
			tile_cachedir = NULL) {

	req <- c("center", "home_east", "home_north", transparency)
	if (length(setdiff(req, names(dat1))) > 0) stop("Missing columns: ", paste(setdiff(req, names(dat1)), collapse = ", "))
	if (is.null(tile_zoom)) tile_zoom <- ifelse(separate, 11, 6)
	if (!is.numeric(nrow) || length(nrow) != 1 || is.na(nrow) || nrow < 1) stop("nrow must be a positive number.")
	nrow <- as.integer(nrow)
	norm <- function(x) tolower(trimws(as.character(x)))

	key <- tibble(center.id = names(center.lst), center.name = unname(center.lst), center.label = paste0(center.id, " (", center.name, ")"))
	parse_centers <- function(x) {
		if (is.null(x) || length(x) == 0 || all(is.na(x) | x == "")) return(key %>% transmute(center.id, panel = center.label, panel.order = row_number()))
		bind_rows(lapply(seq_along(x), function(i) {
			z <- x[i]; hit <- key %>% filter(norm(center.id) == norm(z) | norm(center.name) == norm(z) | norm(center.label) == norm(z))
			if (nrow(hit) == 0) stop("No center matched: ", z)
			hit %>% transmute(center.id, panel = ifelse(any(norm(center.name) == norm(z)), center.name, center.label), panel.order = i)
		})) %>% distinct(center.id, .keep_all = TRUE)
	}

	d <- dat1 %>%
		mutate(
			center0 = sub("\\.0$", "", trimws(as.character(center))),
			center.id = case_when(
				center0 %in% names(center.lst) ~ center0,
				paste0("c", center0) %in% names(center.lst) ~ paste0("c", center0),
				tolower(center0) %in% tolower(unname(center.lst)) ~ names(center.lst)[match(tolower(center0), tolower(unname(center.lst)))],
				TRUE ~ NA_character_
			)
		) %>%
		transmute(center.id, home_east, home_north, alpha.raw = .data[[transparency]]) %>%
		drop_na(center.id, home_east, home_north, alpha.raw) %>%
		left_join(key, by = "center.id") %>%
		inner_join(parse_centers(centers), by = "center.id")
	if (nrow(d) == 0) stop("No data after filtering centers and missing values.")

	rescale_alpha <- function(x) {
		q <- quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE)
		if (!all(is.finite(q)) || q[1] == q[2]) return(rep(mean(alpha_range), length(x)))
		scales::rescale(pmin(pmax(x, q[1]), q[2]), to = alpha_range, from = q)
	}

	cell <- d %>%
		mutate(x = round(home_east / cell_size) * cell_size, y = round(home_north / cell_size) * cell_size) %>%
		group_by(center.id, center.label, panel, panel.order, x, y) %>%
		summarise(n = n(), alpha.raw = median(alpha.raw), .groups = "drop") %>%
		mutate(alpha01 = rescale_alpha(alpha.raw))

	bbox_q <- function(z) {
		p <- (1 - dist_outlier) / 2
		x <- unname(quantile(z$home_east, c(p, 1 - p), na.rm = TRUE))
		y <- unname(quantile(z$home_north, c(p, 1 - p), na.rm = TRUE))
		c(xmin = x[1] - 500, xmax = x[2] + 500, ymin = y[1] - 500, ymax = y[2] + 500)
	}

	bbox_dense <- function(z) {
		if (is.null(width) || !is.numeric(width) || length(width) != 1 || is.na(width) || width <= 0) stop("width must be a positive number in km.")
		w <- width * 1000
		xs <- seq(min(z$home_east), max(z$home_east), by = cell_size)
		ys <- seq(min(z$home_north), max(z$home_north), by = cell_size)
		best <- c(n =  - 1, xmin = NA, xmax = NA, ymin = NA, ymax = NA)
		for (x in xs) for (y in ys) {
			n <- sum(z$home_east >= x - w / 2 & z$home_east <= x + w / 2 & z$home_north >= y - w / 2 & z$home_north <= y + w / 2)
			if (n > best["n"]) best <- c(n = n, xmin = x - w / 2, xmax = x + w / 2, ymin = y - w / 2, ymax = y + w / 2)
		}
		best[ - 1]
	}

	make_sf <- function(z) {
		st_sf(z, geometry = st_sfc(lapply(seq_len(nrow(z)), function(i) {
			h <- cell_size / 2; x <- z$x[i]; y <- z$y[i]
			st_polygon(list(matrix(c(x - h, y - h, x + h, y - h, x + h, y + h, x - h, y + h, x - h, y - h), ncol = 2, byrow = TRUE)))
		}), crs = 27700))
	}

	bb2ll <- function(bb) st_bbox(st_transform(st_as_sfc(st_bbox(c(
		xmin = unname(bb["xmin"]), ymin = unname(bb["ymin"]),
		xmax = unname(bb["xmax"]), ymax = unname(bb["ymax"])
	), crs = 27700)), 4326))
	scale_bar <- function(bb, km = 10) {
		w <- unname(bb["xmax"] - bb["xmin"])
		h <- unname(bb["ymax"] - bb["ymin"])
		len <- km * 1000
		if (!is.finite(w) || !is.finite(h) || w <= len * 1.15 || h <= 0) return(NULL)
		x0 <- unname(bb["xmin"]) + 0.06 * w
		y0 <- unname(bb["ymin"]) + 0.08 * h
		tick <- 0.025 * h
		line <- st_sf(
			geometry = st_sfc(
				st_linestring(matrix(c(x0, y0, x0 + len, y0), ncol = 2, byrow = TRUE)),
				st_linestring(matrix(c(x0, y0 - tick, x0, y0 + tick), ncol = 2, byrow = TRUE)),
				st_linestring(matrix(c(x0 + len, y0 - tick, x0 + len, y0 + tick), ncol = 2, byrow = TRUE)),
				crs = 27700
			)
		)
		label <- st_sf(
			label = paste0(km, " km"),
			geometry = st_sfc(st_point(c(x0 + len / 2, y0 + tick * 2.3)), crs = 27700)
		)
		list(line = line, label = label)
	}
	uk <- rnaturalearth::ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf") %>% st_transform(4326)
	adm <- tryCatch(rnaturalearth::ne_states(country = "United Kingdom", returnclass = "sf") %>% st_transform(4326), error = function(e) NULL)

	draw <- function(z, bb, ttl, fill_var, cols, legend = TRUE) {
		sb <- if (separate) scale_bar(bb, 10) else NULL
		p <- ggplot()
		if (basemap) p <- p + ggspatial::annotation_map_tile(
			type = basemap_type, zoom = tile_zoom, cachedir = tile_cachedir, progress = "none")
		p <- p + geom_sf(data = uk, fill = alpha("#F4F1E8", 0.35), color = alpha("grey40", 0.55), linewidth = 0.25)
		if (!is.null(adm)) p <- p + geom_sf(data = adm, fill = NA, color = alpha("grey35", 0.40), linewidth = 0.18)
		p <- p + geom_sf(data = make_sf(z), aes(fill = .data[[fill_var]], alpha = alpha01), color = NA)
		if (!is.null(sb)) {
			p <- p +
				geom_sf(data = sb$line, inherit.aes = FALSE, color = "white", linewidth = 1.1) +
				geom_sf(data = sb$line, inherit.aes = FALSE, color = "black", linewidth = 0.45) +
				geom_sf_text(data = sb$label, aes(label = label), inherit.aes = FALSE, size = 3.0, fontface = "bold", color = "black")
		}
		p +
				scale_fill_manual(values = cols, name = "Center") + scale_alpha_identity() +
				coord_sf(
					xlim = c(bb["xmin"], bb["xmax"]), ylim = c(bb["ymin"], bb["ymax"]),
					crs = st_crs(27700), default_crs = st_crs(27700),
					datum = NA, expand = FALSE
				) +
				labs(title = ttl) +
				theme_void(base_size = 11) +
				theme(
					plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
					legend.position = ifelse(legend, "right", "none"),
					legend.title = element_text(size = 11),
					legend.text = element_text(size = 9),
					legend.key.height = unit(0.42, "cm"),
					panel.border = element_blank(),
					plot.background = element_rect(fill = "white", color = NA),
					plot.margin = if (separate) margin(6, 6, 6, 6) else margin(4, 5, 4, 5)
				) +
				guides(fill = guide_legend(ncol = 1, byrow = TRUE))
	}

	if (!separate) {
		z <- cell %>% group_by(x, y) %>% arrange(desc(n), center.label) %>% slice(1) %>% ungroup()
		lv <- key %>% filter(center.label %in% unique(z$center.label)) %>% pull(center.label)
		return(draw(z, st_bbox(st_transform(uk, 27700)), "UK Biobank assessment centers", "center.label", setNames(rainbow(length(lv)), lv), TRUE))
	}

	pan <- d %>% count(panel, panel.order, name = "N") %>% arrange(panel.order)
	cols <- setNames(rainbow(nrow(pan)), pan$panel)
		plist <- lapply(pan$panel, function(pp) {
			di <- d %>% filter(panel == pp)
			bb0 <- if (!is.null(width)) bbox_dense(di) else bbox_q(di)
			ci <- cell %>% filter(panel == pp, x >= bb0["xmin"], x <= bb0["xmax"], y >= bb0["ymin"], y <= bb0["ymax"])
			draw(ci, bb0, paste0(pp, " (N=", scales::comma(pan$N[pan$panel == pp]), ")"), "panel", cols, FALSE)
		})
		if (length(plist) > 1) {
			wrap_plots(plist, nrow = nrow) & theme(plot.margin = margin(6, 6, 6, 6))
		} else {
			wrap_plots(plist, nrow = 1)
		}
	}


# 


# 🚩 Time-series plots


library(lattice)
.mh_plot_vectors <- function(chr, pos, pvalue,
	sig.level = NA, annotate = NULL, ann.default = list(),
	should.thin = T, thin.pos.places = 2, thin.logp.places = 2,
	xlab = "Chromosome", ylab = expression( - log[10](p - value)),
	col = c("gray", "darkgray"), panel.extra = NULL, pch = 20, cex = 0.8, ...) {

	if (length(chr) == 0) stop("chromosome vector is empty")
	if (length(pos) == 0) stop("position vector is empty")
	if (length(pvalue) == 0) stop("pvalue vector is empty")

	#make sure we have an ordered factor
	if(!is.ordered(chr)) {
		chr <- ordered(chr)
	} else {
		chr <- chr[, drop = T]
	}

	#make sure positions are in kbp
	if (any(pos > 1e6)) pos <- pos / 1e6;

	#calculate absolute genomic position
	#from relative chromosomal positions
	posmin <- tapply(pos, chr, min);
	posmax <- tapply(pos, chr, max);
	posshift <- head(c(0, cumsum(posmax)), - 1);
	names(posshift) <- levels(chr)
	genpos <- pos + posshift[chr];
	getGenPos <- function(cchr, cpos) {
		p <- posshift[as.character(cchr)] + cpos
		return(p)
	}

	#parse annotations
	grp <- NULL
	ann.settings <- list()
	label.default <- list(x = "peak", y = "peak", adj = NULL, pos = 3, offset = 0.5,
		col = NULL, fontface = NULL, fontsize = NULL, show = F)
	parse.label <- function(rawval, groupname) {
		r <- list(text = groupname)
		if(is.logical(rawval)) {
			if(!rawval) {r$show <- F}
		} else if (is.character(rawval) || is.expression(rawval)) {
			if(nchar(rawval) >= 1) {
				r$text <- rawval
			}
		} else if (is.list(rawval)) {
			r <- modifyList(r, rawval)
		}
		return(r)
	}

	if(!is.null(annotate)) {
		if (is.list(annotate)) {
			grp <- annotate[[1]]
		} else {
			grp <- annotate
		}
		if (!is.factor(grp)) {
			grp <- factor(grp)
		}
	} else {
		grp <- factor(rep(1, times = length(pvalue)))
	}

	ann.settings <- vector("list", length(levels(grp)))
	ann.settings[[1]] <- list(pch = pch, col = col, cex = cex, fill = col, label = label.default)

	if (length(ann.settings) > 1) {
		lcols <- trellis.par.get("superpose.symbol")$col
		lfills <- trellis.par.get("superpose.symbol")$fill
		for(i in 2:length(levels(grp))) {
			ann.settings[[i]] <- list(pch = pch,
				col = lcols[(i - 2) %% length(lcols) + 1 ],
				fill = lfills[(i - 2) %% length(lfills) + 1 ],
				cex = cex, label = label.default);
			ann.settings[[i]]$label$show <- T
		}
	}
	names(ann.settings) <- levels(grp)
	for(i in 1:length(ann.settings)) {
		if (i > 1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
		ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label,
			parse.label(ann.settings[[i]]$label, levels(grp)[i]))
	}
	if(is.list(annotate) && length(annotate) > 1) {
		user.cols <- 2:length(annotate)
		ann.cols <- c()
		if(!is.null(names(annotate[ - 1])) && all(names(annotate[ - 1]) != "")) {
			ann.cols <- match(names(annotate)[ - 1], names(ann.settings))
		} else {
			ann.cols <- user.cols - 1
		}
		for(i in seq_along(user.cols)) {
			if(!is.null(annotate[[user.cols[i]]]$label)) {
				annotate[[user.cols[i]]]$label <- parse.label(annotate[[user.cols[i]]]$label,
					levels(grp)[ann.cols[i]])
			}
			ann.settings[[ann.cols[i]]] <- modifyList(ann.settings[[ann.cols[i]]],
				annotate[[user.cols[i]]])
		}
	}
	 rm(annotate)

	#reduce number of points plotted
	if(should.thin) {
		thinned <- unique(data.frame(
			logp = round( - log10(pvalue), thin.logp.places),
			pos = round(genpos, thin.pos.places),
			chr = chr,
			grp = grp)
		)
		logp <- thinned$logp
		genpos <- thinned$pos
		chr <- thinned$chr
		grp <- thinned$grp
		rm(thinned)
	} else {
		logp <-  - log10(pvalue)
	}
	rm(pos, pvalue)
	gc()

	#custom axis to print chromosome names
	axis.chr <- function(side, ...) {
		if(side == "bottom") {
			panel.axis(side = side, outside = T,
				at = ((posmax + posmin) / 2 + posshift),
				labels = levels(chr),
				ticks = F, rot = 0,
				check.overlap = F
			)
		} else if (side == "top" || side == "right") {
			panel.axis(side = side, draw.labels = F, ticks = F);
		}
		else {
			axis.default(side = side, ...);
		}
	 }

	#make sure the y-lim covers the range (plus a bit more to look nice)
	prepanel.chr <- function(x, y, ...) {
		A <- list();
		maxy <- ceiling(max(y, ifelse(!is.na(sig.level), - log10(sig.level), 0))) + .5;
		A$ylim = c(0, maxy);
		A;
	}

	xyplot(logp ~ genpos, chr = chr, groups = grp,
		axis = axis.chr, ann.settings = ann.settings,
		prepanel = prepanel.chr, scales = list(axs = "i"),
		panel = function(x, y, ..., getgenpos) {
			if(!is.na(sig.level)) {
				#add significance line (if requested)
				panel.abline(h =  - log10(sig.level), lty = 2);
			}
			panel.superpose(x, y, ..., getgenpos = getgenpos);
			if(!is.null(panel.extra)) {
				panel.extra(x, y, getgenpos, ...)
			}
		},
		panel.groups = function(x, y, ..., subscripts, group.number) {
			A <- list(...)
			#allow for different annotation settings
			gs <- ann.settings[[group.number]]
			A$col.symbol <- gs$col[(as.numeric(chr[subscripts]) - 1) %% length(gs$col) + 1]
			A$cex <- gs$cex[(as.numeric(chr[subscripts]) - 1) %% length(gs$cex) + 1]
			A$pch <- gs$pch[(as.numeric(chr[subscripts]) - 1) %% length(gs$pch) + 1]
			A$fill <- gs$fill[(as.numeric(chr[subscripts]) - 1) %% length(gs$fill) + 1]
			A$x <- x
			A$y <- y
			do.call("panel.xyplot", A)
			#draw labels (if requested)
			if(gs$label$show) {
				gt <- gs$label
				names(gt)[which(names(gt) == "text")] <- "labels"
				gt$show <- NULL
				if(is.character(gt$x) | is.character(gt$y)) {
					peak = which.max(y)
					center = mean(range(x))
					if (is.character(gt$x)) {
						if(gt$x == "peak") {gt$x <- x[peak]}
						if(gt$x == "center") {gt$x <- center}
					}
					if (is.character(gt$y)) {
						if(gt$y == "peak") {gt$y <- y[peak]}
					}
				}
				if(is.list(gt$x)) {
					gt$x <- A$getgenpos(gt$x[[1]], gt$x[[2]])
				}
				do.call("panel.text", gt)
			}
		},
		xlab = xlab, ylab = ylab,
		panel.extra = panel.extra, getgenpos = getGenPos, ...
	);
}

# Draw raw and adj2-adjusted Yin-Yang trajectories on the same b2e time scale.
# Both panels carry the same two-row statistical legend inside the plot.


# Draw a Manhattan plot from a standardized small GWAS file. cis_bed supplies
# CHR, START, END, and GENE; cis_gene identifies the red target gene.


# 🚩 Manhattan plots
mh_plot <- function(smalled_gwas_file, col = c("gray", "darkgray"), cis_gene = NULL, cis_bed = NULL, mh_plot_bed = NULL,
		cis_color = "red", other_top_color = "green", other_top_max_per_chr = 2,
		locus_size = 1e6, genomewide_p = 5e-8, main = NULL, ...) {
	if (!file.exists(smalled_gwas_file)) stop("GWAS file does not exist: ", smalled_gwas_file)
	if (requireNamespace("data.table", quietly = TRUE)) {
		dat <- data.table::fread(smalled_gwas_file, showProgress = FALSE, data.table = FALSE)
	} else {
		con <- if (grepl("\\.gz$", smalled_gwas_file, ignore.case = TRUE)) gzfile(smalled_gwas_file) else smalled_gwas_file
		dat <- read.table(con, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE)
	}
	need <- c("CHR", "POS", "P")
	if (!all(need %in% names(dat))) stop("GWAS file must contain columns: ", paste(need, collapse = ", "))
	dat$CHR <- sub("^chr", "", as.character(dat$CHR), ignore.case = TRUE)
	dat$CHR[dat$CHR == "X"] <- "23"; dat$CHR[dat$CHR == "Y"] <- "24"
	dat$CHR <- suppressWarnings(as.integer(dat$CHR)); dat$POS <- suppressWarnings(as.numeric(dat$POS)); dat$P <- suppressWarnings(as.numeric(dat$P))
	positive <- dat$P[is.finite(dat$P) & dat$P > 0]
	p_floor <- if (length(positive)) max(.Machine$double.xmin, min(positive) * 0.1) else .Machine$double.xmin
	dat$P[is.finite(dat$P) & dat$P <= 0] <- p_floor
	dat <- dat[is.finite(dat$CHR) & dat$CHR >= 1 & dat$CHR <= 24 & is.finite(dat$POS) & dat$POS > 0 & is.finite(dat$P) & dat$P > 0 & dat$P <= 1, , drop = FALSE]
	if (!nrow(dat)) stop("No valid CHR/POS/P rows in ", smalled_gwas_file)
	group <- rep("Background", nrow(dat))
	locus_size <- as.numeric(locus_size)
	if (!is.finite(locus_size) || locus_size <= 0) stop("locus_size must be positive")
	read_gene_bed <- function(x, arg_name) {
		if (is.null(x)) return(NULL)
		if (is.character(x) && length(x) == 1) {
			if (!file.exists(x)) stop(arg_name, " does not exist: ", x)
			b <- read.table(x, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
		} else b <- as.data.frame(x, stringsAsFactors = FALSE)
		if (ncol(b) < 4) stop(arg_name, " must contain CHR, START, END, and GENE")
		b <- b[, 1:4, drop = FALSE]; names(b) <- c("CHR", "START", "END", "GENE")
		b$CHR <- sub("^chr", "", as.character(b$CHR), ignore.case = TRUE)
		b$CHR[b$CHR == "X"] <- "23"; b$CHR[b$CHR == "Y"] <- "24"
		b$CHR <- suppressWarnings(as.integer(b$CHR)); b$START <- suppressWarnings(as.numeric(b$START)); b$END <- suppressWarnings(as.numeric(b$END)); b$GENE <- as.character(b$GENE)
		b[is.finite(b$CHR) & is.finite(b$START) & is.finite(b$END) & nzchar(b$GENE), , drop = FALSE]
	}
	bed <- read_gene_bed(cis_bed, "cis_bed")
	primary_bed <- read_gene_bed(mh_plot_bed, "mh_plot_bed")
	if (!is.null(cis_gene)) {
		if (is.null(bed)) stop("cis_bed is required when cis_gene is specified")
		cis <- bed[bed$GENE == as.character(cis_gene), , drop = FALSE]
		if (!nrow(cis)) stop("No cis_bed entry for cis_gene: ", cis_gene)
		cis_center <- (min(cis$START) + max(cis$END)) / 2
		cis_start <- cis_center - locus_size / 2; cis_end <- cis_center + locus_size / 2
		group[dat$CHR == cis$CHR[1] & dat$POS >= cis_start & dat$POS <= cis_end] <- as.character(cis_gene)
	}
	max_per_chr <- as.integer(other_top_max_per_chr)
	if (!is.finite(max_per_chr) || max_per_chr < 0) stop("other_top_max_per_chr must be a non-negative integer")
	eligible <- which(group == "Background" & dat$P <= genomewide_p)
	other_groups <- other_labels <- character()
	if (length(eligible) && max_per_chr > 0) {
		# Define significant loci from SNP coordinates first, independently of gene coverage.
		loci <- do.call(rbind, lapply(split(eligible, dat$CHR[eligible]), function(ii) {
			ii <- ii[order(dat$POS[ii])]
			loc <- cumsum(c(TRUE, diff(dat$POS[ii]) > locus_size))
			do.call(rbind, lapply(split(ii, loc), function(jj) {
				lead <- jj[which.min(dat$P[jj])]
				data.frame(CHR = dat$CHR[lead], POS = dat$POS[lead], P = dat$P[lead],
					LOCUS_START = min(dat$POS[jj]) - locus_size / 2,
					LOCUS_END = max(dat$POS[jj]) + locus_size / 2)
			}))
		}))
		loci <- loci[order(loci$CHR, loci$P, loci$POS), , drop = FALSE]
		keep <- unlist(lapply(split(seq_len(nrow(loci)), loci$CHR), head, max_per_chr), use.names = FALSE)
		loci <- loci[keep, , drop = FALSE]

		near_genes <- function(b, chr, pos, radius = 1e6) {
			if (is.null(b) || !nrow(b)) return(NULL)
			x <- b[b$CHR == chr, , drop = FALSE]
			if (!nrow(x)) return(NULL)
			x$DIST <- pmax(x$START - pos, pos - x$END, 0)
			x <- x[x$DIST <= radius, , drop = FALSE]
			x[order(x$DIST, x$START, x$GENE), , drop = FALSE]
		}
		gene_label <- function(chr, pos) {
			# Exact overlap, primary first and cis BED only as fallback.
			x <- near_genes(primary_bed, chr, pos, 0)
			if (is.null(x) || !nrow(x)) x <- near_genes(bed, chr, pos, 0)
			if (!is.null(x) && nrow(x)) return(paste(head(unique(x$GENE), 3), collapse = ","))
			# No overlap: genes within +/-1 Mb, again primary then fallback.
			x <- near_genes(primary_bed, chr, pos, 1e6)
			if (is.null(x) || !nrow(x)) x <- near_genes(bed, chr, pos, 1e6)
			if (!is.null(x) && nrow(x)) return(paste(head(unique(x$GENE), 3), collapse = ","))
			# No nearby gene: report nearest upstream and downstream genes with
			# signed boundary distance (upstream negative, downstream positive).
			x <- if (!is.null(primary_bed) && any(primary_bed$CHR == chr)) primary_bed[primary_bed$CHR == chr, , drop = FALSE] else bed[bed$CHR == chr, , drop = FALSE]
			if (is.null(x) || !nrow(x)) return(paste0("chr", chr, ":", format(pos, scientific = FALSE)))
			up <- x[x$END < pos, , drop = FALSE]; down <- x[x$START > pos, , drop = FALSE]
			labs <- character()
			if (nrow(up)) { u <- up[which.max(up$END), ]; labs <- c(labs, sprintf("%s (%+.1fMb)", u$GENE, - (pos - u$END) / 1e6)) }
			if (nrow(down)) { d <- down[which.min(down$START), ]; labs <- c(labs, sprintf("%s (%+.1fMb)", d$GENE, (d$START - pos) / 1e6)) }
			paste(labs, collapse = ",")
		}
		for (i in seq_len(nrow(loci))) {
			key <- paste0("Other locus ", i)
			inside <- group == "Background" & dat$CHR == loci$CHR[i] & dat$POS >= loci$LOCUS_START[i] & dat$POS <= loci$LOCUS_END[i]
			group[inside] <- key
			other_groups <- c(other_groups, key)
			other_labels <- c(other_labels, gene_label(loci$CHR[i], loci$POS[i]))
		}
	}
	group <- droplevels(factor(group, levels = unique(c("Background", as.character(cis_gene), other_groups))))
	ann <- list(group)
	if (!is.null(cis_gene) && as.character(cis_gene) %in% levels(group)) ann[[as.character(cis_gene)]] <- list(col = cis_color, fill = cis_color, cex = 1.05, label = list(show = TRUE))
	for (i in seq_along(other_groups)) if (other_groups[i] %in% levels(group)) ann[[other_groups[i]]] <- list(col = other_top_color, fill = other_top_color, cex = 1.05, label = list(show = TRUE, text = other_labels[i]))
	if (is.null(main)) main <- sub("\\.(gz|bgz|tsv|txt)$", "", basename(smalled_gwas_file), ignore.case = TRUE)
	.mh_plot_vectors(dat$CHR, dat$POS, dat$P, sig.level = genomewide_p, annotate = ann, col = col, main = main, ...)
}
