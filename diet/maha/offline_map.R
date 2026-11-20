# Local boundary-only version of the shared plot_map helper. No tile services.
plot_map_offline <- function(dat1, centers = NULL, separate = TRUE, width = 50, nrow = 1,
			dist_outlier = 0.98, transparency = "tdi", center.lst = center.lst,
			cell_size = 1000, alpha_range = c(0.18, 0.95),
			basemap = FALSE) {

	req <- c("center", "home_east", "home_north", transparency)
	if (length(setdiff(req, names(dat1))) > 0) stop("Missing columns: ", paste(setdiff(req, names(dat1)), collapse = ", "))
	if (isTRUE(basemap)) stop("Online basemaps are disabled.")
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
	# Use only boundary data shipped with the preinstalled package.
	if (!requireNamespace("rnaturalearthdata", quietly = TRUE))
		stop("Preinstalled rnaturalearthdata is required for offline maps.")
	boundaries <- new.env()
	utils::data("countries50", package = "rnaturalearthdata", envir = boundaries)
	uk <- sf::st_as_sf(boundaries$countries50)
	uk <- sf::st_transform(uk[uk$admin == "United Kingdom", ], 4326)
	if (!nrow(uk)) stop("United Kingdom boundary missing from local countries50 data.")
	adm <- NULL

	draw <- function(z, bb, ttl, fill_var, cols, legend = TRUE) {
		sb <- if (separate) scale_bar(bb, 10) else NULL
		p <- ggplot()

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
