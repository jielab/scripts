# EMS120 v2 helper functions (2026-08-09): statistical summaries, plotting, phone rules, and policy-window utilities.
# Disease-model cross-validation and keyword learning live in f/ems120.py; execution flow lives in f/ems120.R.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Packages, paths, and shared inputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (!exists("fig_theme")) {
	fig_theme <- function(base_size = 11) {
		theme_minimal(base_size = base_size) +
			theme(
				axis.line = element_line(linewidth = 0.35),
				axis.title = element_text(face = "bold"),
				axis.text = element_text(face = "bold"),
				plot.title = element_text(face = "bold"),
				strip.text = element_text(face = "bold"),
				panel.grid.minor = element_blank()
			)
	}
}

if (!exists("save_plot")) {
	save_plot <- function(p, file, width = 8.6, height = 6, dpi = 600, ...) {
		ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi, limitsize = FALSE, ...)
	}
}

roll3 <- function(x) zoo::rollmean(x, 3, fill = NA, align = "center")
roll7 <- function(x) zoo::rollmean(x, 7, fill = NA, align = "center")

sig_star <- function(p) dplyr::case_when(
	is.na(p) ~ "",
	p < 0.0001 ~ "***",
	p < 0.001 ~ "**",
	p < 0.01 ~ "*",
	TRUE ~ ""
)

# Conventional significance labels used when the figure is meant to show
# formal hypothesis tests (* <0.05, ** <0.01, *** <0.001).
sig_star05 <- function(p) dplyr::case_when(
	is.na(p) ~ "",
	p < 0.001 ~ "***",
	p < 0.01 ~ "**",
	p < 0.05 ~ "*",
	TRUE ~ ""
)

fmt_p <- function(p, digits = 3) {
	p <- suppressWarnings(as.numeric(p))
	ifelse(!is.finite(p), "NA", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = digits)))
}

clamp <- function(x, lo = -Inf, hi = Inf) pmin(pmax(x, lo), hi)
div1 <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b != 0, a / b, NA_real_)

scale01 <- function(x) {
	x <- as.numeric(x); ok <- is.finite(x); out <- rep(NA_real_, length(x))
	if (!any(ok)) return(out)
	rg <- range(x[ok], na.rm = TRUE)
	out[ok] <- if (diff(rg) == 0) 0.5 else (x[ok] - rg[1]) / diff(rg)
	out
}

scale_0_10 <- function(x) 10 * scale01(x)

z_score <- function(x) {
	x <- as.numeric(x); ok <- is.finite(x); out <- rep(NA_real_, length(x))
	if (!any(ok)) return(out)
	s <- sd(x[ok], na.rm = TRUE)
	out[ok] <- if (!is.finite(s) || s == 0) 0 else (x[ok] - mean(x[ok], na.rm = TRUE)) / s
	out
}

calc_smd <- function(x, g) {
	x <- as.numeric(x); g <- as.character(g)
	x1 <- x[g == "high" & is.finite(x)]
	x0 <- x[g == "low" & is.finite(x)]
	if (length(x1) < 20 || length(x0) < 20) return(NA_real_)
	sp <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
	if (!is.finite(sp) || sp == 0) return(NA_real_)
	(mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / sp
}

get_or <- function(fit, target_term) {
	x <- tryCatch(broom::tidy(fit), error = function(e) tibble()) %>%
		filter(.data$term == .env$target_term)
	if (!nrow(x)) {
		return(tibble(term = target_term, OR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
	}
	tibble(
		term = x$term,
		OR = exp(x$estimate),
		lo = exp(x$estimate - 1.96 * x$std.error),
		hi = exp(x$estimate + 1.96 * x$std.error),
		p = x$p.value
	)
}

q_group_by_year <- function(score, year, probs = c(.2, .8)) {
	grp <- rep(NA_character_, length(score))
	for (yy in sort(unique(year))) {
		idx <- which(year == yy & is.finite(score))
		if (length(idx) < 100) next
		r <- rank(score[idx], ties.method = "average", na.last = "keep")
		pct <- (r - 0.5) / length(idx)
		grp[idx[pct <= probs[1]]] <- "low"
		grp[idx[pct >= probs[2]]] <- "high"
	}
	factor(grp, levels = c("low", "high"))
}

score_decile_by_year <- function(score, year) {
	out <- rep(NA_integer_, length(score))
	for (yy in sort(unique(year))) {
		idx <- which(year == yy & is.finite(score))
		if (length(idx) < 100) next
		r <- rank(score[idx], ties.method = "average", na.last = "keep")
		pct <- (r - 0.5) / length(idx)
		out[idx] <- pmin(10L, pmax(1L, ceiling(10 * pct)))
	}
	out
}

read_center_image <- function(path, circle_alpha = TRUE, radius = 0.492) {
	ext <- tolower(tools::file_ext(path))
	if (ext %in% c("jpg", "jpeg")) {
		img <- jpeg::readJPEG(path)
		nr <- dim(img)[1]; nc <- dim(img)[2]
		img <- array(c(img, rep(1, nr * nc)), dim = c(nr, nc, 4))
	} else if (ext %in% c("png")) {
		img <- png::readPNG(path)
		if (length(dim(img)) == 2) {
			nr <- dim(img)[1]; nc <- dim(img)[2]
			img <- array(c(img, img, img, rep(1, nr * nc)), dim = c(nr, nc, 4))
		}
		if (dim(img)[3] == 3) {
			nr <- dim(img)[1]; nc <- dim(img)[2]
			img <- array(c(img, rep(1, nr * nc)), dim = c(nr, nc, 4))
		}
	} else {
		stop("Unsupported image format: ", path)
	}

	nr <- dim(img)[1]; nc <- dim(img)[2]
	yy <- (seq_len(nr) - 0.5 - nr / 2) / nr
	xx <- (seq_len(nc) - 0.5 - nc / 2) / nc
	rr <- outer(yy, xx, function(y, x) sqrt(x^2 + y^2))
	if (circle_alpha) img[, , 4][rr > radius] <- 0

	# Some downloaded PNGs contain a checkerboard background instead of true alpha.
	# Remove near-white / pale-grey checkerboard pixels outside the circular medallion.
	edge_bg <- rr > min(0.46, radius - 0.02)
	near_checker <- (
		(img[, , 1] > 0.90 & img[, , 2] > 0.90 & img[, , 3] > 0.90) |
		(abs(img[, , 1] - img[, , 2]) < 0.02 & abs(img[, , 1] - img[, , 3]) < 0.02 & img[, , 1] > 0.75)
	)
	img[, , 4][edge_bg & near_checker] <- 0
	img
}

close_hourly_line <- function(dat, x_col = "hour", y_col = "pct_smooth") {
	dat <- dat %>% arrange(.data[[x_col]])
	if (!nrow(dat)) return(dat)
	bind_rows(
		dat,
		dat %>% slice(1) %>% mutate(!!x_col := 24L)
	)
}

score_bin_label <- function(x) {
	x <- suppressWarnings(as.numeric(x))
	out <- rep(NA_character_, length(x))
	ok <- is.finite(x)
	out[!ok] <- "NA"
	out[ok & x <= 0] <- "0"
	k <- ceiling(x[ok & x > 0])
	out[ok & x > 0] <- sprintf("(%d,%d]", k - 1L, k)
	out
}

print_phone_score_bins <- function(d, year, score_vars = c("phone.sco0", "phone.sco1", "phone.sco2", "phone.sco3")) {
	cat(sprintf("\nPhone score summary, year %s:\n", year))
	for (v in score_vars) {
		cat(sprintf("\n%s\n", v))
		if (!v %in% names(d)) {
			cat("  missing\n")
			next
		}
		x <- suppressWarnings(as.numeric(d[[v]]))
		print(table(score_bin_label(x), useNA = "always"))
	}
}

clean_xia_colnames <- function(x) {
	x <- trimws(as.character(x))
	x <- gsub("（米）|\\(米\\)", "", x)
	x <- gsub("\\s+", "", x)
	x
}

merge_xia_geo_year <- function(dat, year, dir_geo_xia, vars.basic.ems = NULL) {
	if (is.null(dat) || !nrow(dat)) {
		return(list(data = dat, log = tibble(year = as.integer(year), records_raw = 0L, records_geo_xlsx = NA_integer_, matched_records = 0L, status = "empty_raw")))
	}
	geo_file <- file.path(dir_geo_xia, paste0(year, ".geo.xlsx"))
	if (!file.exists(geo_file)) {
		warning("Missing xia geo file: ", geo_file)
		return(list(data = dat, log = tibble(year = as.integer(year), records_raw = nrow(dat), records_geo_xlsx = NA_integer_, matched_records = 0L, status = "missing_geo_xlsx")))
	}

	xia <- suppressMessages(readxl::read_excel(geo_file))
	names(xia) <- make.unique(clean_xia_colnames(names(xia)), sep = "_")
	if (!"ID" %in% names(xia)) xia <- xia %>% mutate(ID = row_number(), .before = 1)
	if (!"ID" %in% names(dat)) dat <- dat %>% mutate(ID = row_number(), .before = 1)

	xia <- xia %>% mutate(.merge_key = as.character(ID), .xia_row = row_number()) %>% distinct(.merge_key, .keep_all = TRUE)
	dat2 <- dat %>% dplyr::select(-matches("\\.(x|y)$")) %>% mutate(.merge_key = as.character(ID))
	payload_cols <- setdiff(names(xia), c("ID", "地址", ".merge_key", ".xia_row"))
	drop_from_main <- intersect(names(dat2), payload_cols)
	dat_join <- dat2 %>% dplyr::select(-any_of(drop_from_main))
	xia_join <- xia %>% dplyr::select(.merge_key, all_of(payload_cols))
	merged <- dat_join %>% left_join(xia_join, by = ".merge_key")

	check_cols <- c("地址", "地址类型", "接车地址经度", "接车地址纬度", "房价指数")
	log_row <- tibble(
		year = as.integer(year),
		records_raw = nrow(dat),
		records_geo_xlsx = nrow(xia),
		matched_records = sum(!is.na(dat2$.merge_key) & dat2$.merge_key %in% xia$.merge_key),
		dropped_main_columns = paste(drop_from_main, collapse = ";"),
		status = "merged"
	)
	for (cc in check_cols) {
		log_row[[paste0("missing_", cc)]] <- if (cc %in% names(merged)) sum(is.na(merged[[cc]]) | trimws(as.character(merged[[cc]])) == "", na.rm = TRUE) else NA_integer_
	}

	merged <- merged %>% dplyr::select(-.merge_key)
	cat(sprintf("Year %s, %s records from %s; unmatched/missing: 地址 %s, 地址类型 %s, 接车地址经度 %s, 接车地址纬度 %s\n",
		year, nrow(xia), basename(geo_file), log_row$missing_地址, log_row$missing_地址类型, log_row$missing_接车地址经度, log_row$missing_接车地址纬度))
	list(data = merged, log = log_row)
}

empty_panel <- function(title = NULL, msg = "No analyzable data") {
	ggplot() + annotate("text", x = 0, y = 0, label = msg, fontface = "bold", size = 4) +
		labs(title = title) + theme_void(base_size = 10) +
		theme(plot.title = element_text(face = "bold", size = 10))
}

clean_phone_value <- function(x) {
	x <- as.character(x)
	x <- trimws(x)
	x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
	x <- gsub("[^0-9]", "", x)
	x <- sub("^0+", "", x)
	ifelse(nchar(x) == 11, x, NA_character_)
}

score_phone_simple <- function(phone_tail) {
	p <- as.character(phone_tail)
	p <- trimws(p)
	p[is.na(p)] <- ""
	out_score <- rep(NA_real_, length(p))
	out_reason <- rep("", length(p))
	ok <- grepl("^[0-9]{8}$", p)
	has4 <- ok & grepl("4", p, fixed = TRUE)
	out_score[has4] <- 0
	out_reason[has4] <- "has4 => 0"
	no4 <- ok & !has4
	if (any(no4)) {
		count_digit <- function(z, d) lengths(regmatches(z, gregexpr(d, z, fixed = TRUE)))
		n8 <- count_digit(p[no4], "8")
		n9 <- count_digit(p[no4], "9")
		n6 <- count_digit(p[no4], "6")
		n1 <- count_digit(p[no4], "1")
		raw <- n8 * 1 + (n9 + n6 + n1) * 0.5
		s <- pmin(10, raw / 8 * 10)
		out_score[no4] <- s
		out_reason[no4] <- sprintf("raw=8(%s*1)+9(%s*0.5)+6(%s*0.5)+1(%s*0.5)=%.1f; rescaled to %.1f/10", n8, n9, n6, n1, raw, s)
	}
	tibble(phone.sco = out_score, phone.sco.reason = out_reason)
}

audit_phone_vars <- function(path, year, phone_var_patterns) {
	dat <- suppressMessages(readxl::read_excel(path))
	names(dat) <- trimws(names(dat))
	phone_vars <- unique(unlist(lapply(phone_var_patterns, function(p) grep(p, names(dat), value = TRUE))))
	if (!length(phone_vars)) {
		cat(sprintf("year %s: no matched phone variables\n", year))
		return(list(
			missing = tibble(year = as.integer(year), variable = NA_character_, n = nrow(dat), missing_n = nrow(dat), missing_pct = 1),
			concordance = tibble(year = integer(), var1 = character(), var2 = character(), n_both = integer(), concordant_n = integer(), concordance = numeric())
		))
	}
	phone_clean <- setNames(lapply(phone_vars, function(v) clean_phone_value(dat[[v]])), phone_vars)
	missing <- tibble(
		year = as.integer(year),
		variable = phone_vars,
		n = nrow(dat),
		missing_n = vapply(phone_clean, function(x) sum(is.na(x)), integer(1)),
		missing_pct = missing_n / n
	)
	concordance <- tibble(year = integer(), var1 = character(), var2 = character(), n_both = integer(), concordant_n = integer(), concordance = numeric())
	if (length(phone_vars) >= 2) {
		pairs <- utils::combn(phone_vars, 2, simplify = FALSE)
		concordance <- purrr::map_dfr(pairs, function(pp) {
			x <- phone_clean[[pp[1]]]
			y <- phone_clean[[pp[2]]]
			ok <- !is.na(x) & !is.na(y)
			tibble(
				year = as.integer(year),
				var1 = pp[1],
				var2 = pp[2],
				n_both = sum(ok),
				concordant_n = sum(x[ok] == y[ok]),
				concordance = ifelse(sum(ok) > 0, sum(x[ok] == y[ok]) / sum(ok), NA_real_)
			)
		})
	}
	cat(sprintf("year %s:\n", year))
	for (i in seq_len(nrow(missing))) {
		cat(sprintf(
			"  %s, missing %s/%s (%.1f%%)\n",
			missing$variable[i], missing$missing_n[i], missing$n[i], 100 * missing$missing_pct[i]
		))
	}
	if (nrow(concordance)) {
		for (i in seq_len(nrow(concordance))) {
			cat(sprintf(
				"  concordance: %s vs %s = %.3f (%s/%s)\n",
				concordance$var1[i], concordance$var2[i], concordance$concordance[i],
				concordance$concordant_n[i], concordance$n_both[i]
			))
		}
	}
	flush.console()
	list(missing = missing, concordance = concordance)
}

max_run_fun <- function(z) max(rle(strsplit(z, "", fixed = TRUE)[[1]])$lengths)

run_digit_fun <- function(z, target) {
	r <- rle(strsplit(z, "", fixed = TRUE)[[1]])
	ifelse(any(r$values == target), max(r$lengths[r$values == target]), 0)
}

seq_max_fun <- function(z) {
	x <- as.integer(strsplit(z, "", fixed = TRUE)[[1]])
	dx <- diff(x)
	best <- cur <- 1
	last <- 99L
	for (v in dx) {
		if (v %in% c(1, -1) && v == last) cur <- cur + 1
		else if (v %in% c(1, -1)) cur <- 2
		else cur <- 1
		last <- ifelse(v %in% c(1, -1), v, 99L)
		best <- max(best, cur)
	}
	best
}

build_daily_dx <- function(start_date, end_date) {
	yrs <- sort(unique(lubridate::year(seq.Date(start_date, end_date, by = "day"))))
	raw <- bind_rows(lapply(yrs, function(y) {
		d <- dat1.list[[as.character(y)]]
		if (is.null(d) || !nrow(d)) return(tibble())
		d %>% transmute(
			date = as.Date(日期),
			dx_grp = factor(as.character(dx_grp), levels = dxs.all),
			phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high"))
		) %>% filter(date >= start_date, date <= end_date, phone.luck %in% c("low", "high"), dx_grp %in% dxs.all)
	}))
	daily_total <- raw %>% count(date, phone.luck, name = "total_group_day") %>% complete(date = seq.Date(start_date, end_date, by = "day"), phone.luck = c("low", "high"), fill = list(total_group_day = 0))
	daily_dx <- raw %>% count(date, phone.luck, dx_grp, name = "count") %>%
		complete(date = seq.Date(start_date, end_date, by = "day"), phone.luck = c("low", "high"), dx_grp = dxs.all, fill = list(count = 0)) %>%
		left_join(daily_total, by = c("date", "phone.luck")) %>%
		mutate(rate_per_1000 = div1(count * 1000, total_group_day), high = as.integer(phone.luck == "high"), dow = factor(lubridate::wday(date, label = TRUE, week_start = 1)))
	list(raw = raw, daily_total = daily_total, daily_dx = daily_dx)
}

fit_period_did <- function(dat, period_var = "period", ref_period, contrast_periods) {
	dat <- dat %>% mutate(periodF = stats::relevel(factor(.data[[period_var]]), ref = ref_period))
	map_dfr(dxs.all, function(dx) {
		d <- dat %>% filter(dx_grp == dx, total_group_day > 0)
		if (nrow(d) < 10 || sum(d$count) < 20) return(tibble(dx_grp = dx, period = contrast_periods, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
		fit <- tryCatch(glm(count ~ high * periodF + dow + offset(log(total_group_day)), family = poisson(), data = d), error = function(e) NULL)
		if (is.null(fit)) return(tibble(dx_grp = dx, period = contrast_periods, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
		td <- broom::tidy(fit)
		map_dfr(contrast_periods, function(pp) {
			term <- paste0("high:periodF", pp)
			x <- td %>% filter(term == .env$term)
			if (!nrow(x)) return(tibble(dx_grp = dx, period = pp, RR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
			tibble(dx_grp = dx, period = pp, RR = exp(x$estimate), lo = exp(x$estimate - 1.96 * x$std.error), hi = exp(x$estimate + 1.96 * x$std.error), p = x$p.value)
		})
	}) %>% group_by(period) %>% mutate(p_adj = p.adjust(p, "BH"), sig = sig_star(p_adj), sig05 = sig_star05(p_adj), p_label = fmt_p(p_adj)) %>% ungroup()
}

fit_total_period_did <- function(dat, period_var = "period", ref_period, contrast_periods) {
	dat <- dat %>% filter(phone.luck %in% c("low", "high"), total_group_day >= 0) %>%
		mutate(high = as.integer(phone.luck == "high"), periodF = stats::relevel(factor(.data[[period_var]]), ref = ref_period))
	if (!"dow" %in% names(dat)) dat <- dat %>% mutate(dow = factor(lubridate::wday(date, label = TRUE, week_start = 1)))
	empty <- tibble(period = contrast_periods, RR = NA_real_, lo = NA_real_, hi = NA_real_, p_poisson = NA_real_, overdispersion = NA_real_, p = NA_real_, p_adj = NA_real_, sig05 = "", p_label = "NA")
	if (nrow(dat) < 10 || length(unique(dat$high)) < 2) return(empty)
	fit <- tryCatch(glm(total_group_day ~ high * periodF + dow, family = poisson(), data = dat), error = function(e) NULL)
	if (is.null(fit)) return(empty)
	td <- broom::tidy(fit)
	disp <- suppressWarnings(sum(stats::residuals(fit, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(fit))
	if (!is.finite(disp) || disp < 1) disp <- 1
	out <- purrr::map_dfr(contrast_periods, function(pp) {
		term <- paste0("high:periodF", pp)
		x <- td %>% filter(term == .env$term)
		if (!nrow(x)) return(tibble(period = pp, RR = NA_real_, lo = NA_real_, hi = NA_real_, p_poisson = NA_real_, overdispersion = disp, p = NA_real_))
		se <- x$std.error * sqrt(disp)
		tibble(period = pp, RR = exp(x$estimate), lo = exp(x$estimate - 1.96 * se), hi = exp(x$estimate + 1.96 * se), p_poisson = x$p.value, overdispersion = disp, p = 2 * stats::pnorm(abs(x$estimate / se), lower.tail = FALSE))
	})
	out %>% mutate(p_adj = p.adjust(p, "BH"), sig05 = sig_star05(p_adj), p_label = fmt_p(p_adj))
}


fit_period_did_counts <- function(dat, period_var = "period", ref_period, contrast_periods, min_cell_n = 20) {
	did <- fit_period_did(dat, period_var = period_var, ref_period = ref_period, contrast_periods = contrast_periods)
	if (!nrow(did)) return(did)

	cnt0 <- dat %>%
		filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high")) %>%
		mutate(period_chr = as.character(.data[[period_var]])) %>%
		filter(period_chr %in% c(ref_period, contrast_periods)) %>%
		group_by(dx_grp, period_chr, phone.luck) %>%
		summarise(
			days = n_distinct(date),
			calls = sum(count, na.rm = TRUE),
			group_calls = sum(total_group_day, na.rm = TRUE),
			rate = ifelse(group_calls > 0, calls / group_calls, NA_real_),
			.groups = "drop"
		)

	cnt <- purrr::map_dfr(contrast_periods, function(pp) {
		x <- cnt0 %>%
			filter(period_chr %in% c(ref_period, pp)) %>%
			mutate(period_role = ifelse(period_chr == ref_period, "ref", "event"),
				cell = paste(period_role, phone.luck, sep = "_")) %>%
			dplyr::select(dx_grp, cell, calls, group_calls, rate) %>%
			pivot_wider(names_from = cell, values_from = c(calls, group_calls, rate), values_fill = 0)
		for (cc in c("calls_ref_low", "calls_ref_high", "calls_event_low", "calls_event_high", "group_calls_ref_low", "group_calls_ref_high", "group_calls_event_low", "group_calls_event_high")) {
			if (!cc %in% names(x)) x[[cc]] <- 0
		}
		x %>% mutate(
			period = pp,
			min_cell_calls = pmin(calls_ref_low, calls_ref_high, calls_event_low, calls_event_high, na.rm = TRUE),
			low_count = is.finite(min_cell_calls) & min_cell_calls < min_cell_n,
			n_label = sprintf("event n=%s/%s; ref n=%s/%s", calls_event_high, calls_event_low, calls_ref_high, calls_ref_low)
		)
	})

	did %>%
		left_join(cnt, by = c("dx_grp", "period")) %>%
		mutate(
			low_count = ifelse(is.na(low_count), TRUE, low_count),
			label = ifelse(is.finite(RR), sprintf("%.2f%s%s\n(%.2f–%.2f)", RR, sig05, ifelse(low_count, "\u2020", ""), lo, hi), ""),
			short_label = ifelse(is.finite(RR), sprintf("%.2f%s%s", RR, sig05, ifelse(low_count, "\u2020", "")), "")
		)
}

fit_mix_global_lrt <- function(dat, period_var = "period", ref_period, contrast_periods) {
	purrr::map_dfr(contrast_periods, function(pp) {
		d <- dat %>%
			filter(dx_grp %in% dxs.all, phone.luck %in% c("low", "high"), total_group_day > 0) %>%
			mutate(period_chr = as.character(.data[[period_var]])) %>%
			filter(period_chr %in% c(ref_period, pp)) %>%
			mutate(
				high = as.integer(phone.luck == "high"),
				periodF = stats::relevel(factor(period_chr), ref = ref_period),
				dxF = stats::relevel(factor(as.character(dx_grp), levels = dxs.all), ref = dxs.all[1])
			)
		if (!"dow" %in% names(d)) d <- d %>% mutate(dow = factor(lubridate::wday(date, label = TRUE, week_start = 1)))
		d <- d %>% filter(!is.na(dxF), !is.na(periodF), !is.na(dow))
		if (nrow(d) < 24 || length(unique(d$high)) < 2 || length(unique(d$periodF)) < 2) {
			return(tibble(period = pp, df = NA_real_, deviance = NA_real_, p_chisq = NA_real_, overdispersion = NA_real_, p_overdisp = NA_real_, p_label = "NA"))
		}
		f0 <- tryCatch(glm(count ~ dxF + high + periodF + dxF:periodF + dxF:high + high:periodF + dow + offset(log(total_group_day)), family = poisson(), data = d), error = function(e) NULL)
		f1 <- tryCatch(glm(count ~ dxF * high * periodF + dow + offset(log(total_group_day)), family = poisson(), data = d), error = function(e) NULL)
		if (is.null(f0) || is.null(f1)) {
			return(tibble(period = pp, df = NA_real_, deviance = NA_real_, p_chisq = NA_real_, overdispersion = NA_real_, p_overdisp = NA_real_, p_label = "NA"))
		}
		a <- tryCatch(anova(f0, f1, test = "Chisq"), error = function(e) NULL)
		if (is.null(a) || nrow(a) < 2) {
			return(tibble(period = pp, df = NA_real_, deviance = NA_real_, p_chisq = NA_real_, overdispersion = NA_real_, p_overdisp = NA_real_, p_label = "NA"))
		}
		dev <- suppressWarnings(as.numeric(a$Deviance[2]))
		df <- suppressWarnings(as.numeric(a$Df[2]))
		p0 <- suppressWarnings(as.numeric(a$`Pr(>Chi)`[2]))
		disp <- suppressWarnings(sum(stats::residuals(f1, type = "pearson")^2, na.rm = TRUE) / stats::df.residual(f1))
		p_over <- ifelse(is.finite(dev) & is.finite(df) & df > 0, stats::pchisq(dev / max(1, disp), df = df, lower.tail = FALSE), NA_real_)
		tibble(
			period = pp,
			df = df,
			deviance = dev,
			p_chisq = p0,
			overdispersion = disp,
			p_overdisp = p_over,
			p_label = fmt_p(p_over)
		)
	}) %>%
		mutate(test = "global disease-specific high-score × period interaction; overdispersion-adjusted P used for figure title")
}


plot_dx_raw_trend <- function(dat_list, years, dx_var, dxs.cn, level = c("raw", "group"), show = c("percent", "count")) {
	level <- match.arg(level)
	show <- match.arg(show)
	dxs_raw <- trimws(as.character(unlist(dxs.cn, use.names = FALSE)))
	map_grp2 <- utils::stack(dxs.cn) %>% setNames(c("dx_raw", "group")) %>% mutate(across(everything(), function(x) trimws(as.character(x))))
	dat <- purrr::map_dfr(years, function(y) {
		d <- dat_list[[as.character(y)]]
		if (is.null(d) || !nrow(d)) return(tibble())
		d0 <- d %>% transmute(dx_raw = trimws(as.character(.data[[dx_var]]))) %>% filter(!is.na(dx_raw), dx_raw != "", dx_raw != "NA", dx_raw != "NaN")
		out <- if (level == "raw") {
			d0 %>% filter(dx_raw %in% dxs_raw) %>% count(group = dx_raw, name = "count")
		} else {
			d0 %>% left_join(map_grp2, by = "dx_raw") %>% count(group, name = "count")
		}
		out %>% mutate(group = trimws(as.character(group))) %>% filter(!is.na(group), group != "", group != "NA", group != "NaN") %>% mutate(year = y, pct = count / sum(count))
	})
	lev <- dat %>% filter(year == max(years), !is.na(group), group != "", group != "NA") %>% arrange(desc(pct)) %>% pull(group) %>% unique()
	dat <- dat %>% filter(!is.na(group), group != "", group != "NA") %>% mutate(group = factor(group, levels = lev), y = if (show == "percent") pct else count) %>% filter(!is.na(group))
	p <- ggplot(dat, aes(year, y, color = group, group = group)) +
		geom_line(linewidth = .95) + geom_point(size = 1.9) +
		scale_x_continuous(breaks = years) +
		scale_color_discrete(drop = TRUE, na.translate = FALSE) +
		{if (show == "percent") scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0)) else scale_y_continuous(expand = c(0, 0))} +
		labs(title = NULL, x = "Year", y = if (show == "percent") "Percentage" else "Count", color = NULL) +
		fig_theme(base_size = 9.5) +
		theme(
			legend.position = "right",
			legend.text = element_text(size = 8, face = "plain"),
			legend.title = element_text(size = 8, face = "plain"),
			legend.key.size = grid::unit(0.35, "cm")
		)
	list(plot = p, data = dat)
}

hourly_label_map <- c(
	Traffic = "Traffic",
	Intoxication = "Intoxication",
	Trauma = "Trauma",
	CVD = "CVD",
	Respiratory = "Respiratory",
	Psychiatric = "Psychiatric",
	`NCD-Other` = "NCD-Other",
	Death = "Death"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Phone/ML pipeline helpers from ems120_main.f.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Auto-extracted helper functions for ems120.R
# Keep analysis settings and execution flow in the main script.
quasi_quarter <- function(x, low_p = 0.25, middle_p = 0.5, high_p = 0.25) { out <- rep(NA_character_, length(x)); ok <- is.finite(x); if (!any(ok)) return(factor(out, levels = c("low", "middle", "high"))); tab <- table(x[ok]); tab <- tab[order(as.numeric(names(tab)))]; low_cut <- as.numeric(names(tab)[which(cumsum(tab) >= low_p * sum(tab))[1]]); high_cut <- as.numeric(names(rev(tab))[which(cumsum(rev(tab)) >= high_p * sum(tab))[1]]); out[ok] <- ifelse(x[ok] <= low_cut, "low", ifelse(x[ok] >= high_cut, "high", "middle")); factor(out, levels = c("low", "middle", "high")) }

normalize_phone_tail_r <- function(x) { x <- ifelse(is.na(x), "", as.character(x)); x <- sub("\\.0$", "", trimws(x)); x <- gsub("[^0-9]", "", x); ifelse(nchar(x) == 11, substr(x, 4, 11), x) }

seq_matches_sco0 <- function(p, seq_min = 3) { out <- list(); start <- 1; cur_step <- NA_integer_; cur_len <- 1; flush <- function(s, l, st) if (l >= seq_min && !is.na(st) && st %in% c(1L, -1L)) out[[length(out) + 1]] <<- list(sub = substr(p, s, s + l - 1), len = l, step = st); for (i in 2:8) { diff <- as.integer(substr(p, i, i)) - as.integer(substr(p, i - 1, i - 1)); if (abs(diff) == 1) { if (is.na(cur_step)) { cur_step <- diff; cur_len <- 2; start <- i - 1 } else if (diff == cur_step) cur_len <- cur_len + 1 else { flush(start, cur_len, cur_step); cur_step <- diff; cur_len <- 2; start <- i - 1 } } else { flush(start, cur_len, cur_step); cur_step <- NA_integer_; cur_len <- 1 } }; flush(start, cur_len, cur_step); out }

score_phone_sco0_one <- function(p) { p <- normalize_phone_tail_r(p); if (!(nchar(p) == 8 && grepl("^[0-9]{8}$", p))) return(c(score = NA_character_, reason = "")); terms <- c("base(4)"); raw <- 4; pat <- 0; for (d in names(c("8" = 1, "9" = 1, "1" = .5, "6" = .5, "4" = -4))) { w <- c("8" = 1, "9" = 1, "1" = .5, "6" = .5, "4" = -4)[[d]]; cc <- lengths(regmatches(p, gregexpr(d, p, fixed = TRUE))); if (cc > 0) { raw <- raw + cc * w; terms <- c(terms, sprintf("%s(%d*%g)", d, cc, w)) } }; rr <- rle(strsplit(p, "", fixed = TRUE)[[1]]); for (i in seq_along(rr$lengths)) if (rr$lengths[i] >= 3) { add <- rr$lengths[i] * .8; pat <- pat + add; terms <- c(terms, sprintf("run %s(%d*0.8)", paste(rep(rr$values[i], rr$lengths[i]), collapse = ""), rr$lengths[i])) }; for (z in seq_matches_sco0(p, 3)) { add <- z$len * .6; pat <- pat + add; terms <- c(terms, sprintf("seq %s(%d*0.6)", z$sub, z$len)) }; if (substr(p, 1, 4) == substr(p, 5, 8) && length(unique(strsplit(substr(p, 1, 4), "", fixed = TRUE)[[1]])) > 1) { pat <- pat + 3.5; terms <- c(terms, sprintf("rep4 %sx2(3.5)", substr(p, 1, 4))) }; if (paste(rep(substr(p, 1, 2), 4), collapse = "") == p && length(unique(strsplit(substr(p, 1, 2), "", fixed = TRUE)[[1]])) > 1) { pat <- pat + 2; terms <- c(terms, sprintf("rep2 %sx4(2)", substr(p, 1, 2))) }; if (substr(p, 5, 8) == paste(rev(strsplit(substr(p, 1, 4), "", fixed = TRUE)[[1]]), collapse = "") && length(unique(strsplit(substr(p, 1, 4), "", fixed = TRUE)[[1]])) > 1) { pat <- pat + 2.5; terms <- c(terms, sprintf("mirror %s|%s(2.5)", substr(p, 1, 4), substr(p, 5, 8))) }; if (pat > 6) { terms <- c(terms, "pat_cap(6)"); pat <- 6 }; final <- max(0, min(10, round(raw + pat, 1))); c(score = as.character(final), reason = paste0(paste(terms, collapse = " + "), sprintf(" = %g", final))) }

normalize_phone_score_names <- function(d) { d }

apply_phone_luck_use <- function(d) { if (is.null(d) || !nrow(d)) return(d); d <- normalize_phone_score_names(d); if (!"phone.luck.sco0" %in% names(d)) d$phone.luck.sco0 <- with(d, case_when(is.na(phone.sco0) ~ NA_character_, phone.sco0 <= 2 ~ "low", phone.sco0 <= 7 ~ "middle", phone.sco0 <= 10 ~ "high", TRUE ~ NA_character_)); d$phone.luck.sco0 <- factor(as.character(d$phone.luck.sco0), levels = c("low", "middle", "high")); for (i in 1:3) { score_i <- paste0("phone.sco", i); luck_i <- paste0("phone.luck.sco", i); if (score_i %in% names(d)) d[[luck_i]] <- quasi_quarter(as.numeric(d[[score_i]]), 0.25, 0.5, 0.25) }; score_col <- paste0("phone.sco", sub("^sco", "", phone.score.use)); reason_col <- paste0(score_col, ".reason"); luck_col <- if (phone.luck.use %in% names(d)) phone.luck.use else "phone.luck.quarter"; d %>% mutate(phone.sco = as.numeric(.data[[score_col]]), phone.sco.reason = as.character(.data[[reason_col]]), phone.luck.quarter = quasi_quarter(phone.sco, 0.25, 0.5, 0.25), phone.luck = if (phone.luck.use == "phone.luck.quarter") phone.luck.quarter else .data[[luck_col]], phone.luck = factor(as.character(phone.luck), levels = c("low", "middle", "high"))) }

plot_design_flow <- function(nodes, edges = NULL, title = NULL, subtitle = NULL, file = NULL, width = 8.6, height = 4.8) {
	stopifnot(all(c("id", "x", "y", "label", "type") %in% names(nodes)))
	if (is.null(edges)) edges <- tibble::tibble(from = head(nodes$id, -1), to = tail(nodes$id, -1))
	nodes <- nodes %>% mutate(label = stringr::str_replace_all(label, "\\n", "\n"), type = factor(type, levels = unique(type)))
	edges <- edges %>% left_join(nodes %>% dplyr::select(from = id, x0 = x, y0 = y), by = "from") %>% left_join(nodes %>% dplyr::select(to = id, x1 = x, y1 = y), by = "to")
	p <- ggplot() +
		geom_segment(data = edges, aes(x = x0 + 0.42 * sign(x1 - x0), y = y0, xend = x1 - 0.42 * sign(x1 - x0), yend = y1),
			arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"), linewidth = 0.65, color = "grey35") +
		geom_label(data = nodes, aes(x, y, label = label, fill = type), linewidth = 0.25, label.r = grid::unit(0.13, "lines"),
			fontface = "bold", size = 3.1, lineheight = 0.92, color = "black") +
		scale_fill_brewer(palette = "Set2", guide = "none") +
		labs(title = title, subtitle = subtitle) +
		coord_cartesian(xlim = range(nodes$x) + c(-0.85, 0.85), ylim = range(nodes$y) + c(-0.65, 0.65), expand = FALSE) +
		theme_void(base_size = 11) +
		theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5), plot.subtitle = element_text(size = 9.5, hjust = 0.5, color = "grey25"), plot.margin = margin(12, 12, 12, 12))
	if (!is.null(file)) save_plot(p, file, width = width, height = height, dpi = 600)
	p
}

make_ml_design_figures <- function() {
	ml_phone_nodes <- tibble::tribble(
		~id, ~x, ~y, ~label, ~type,
		"input", 1, 2, "Phone-number tail\n8-digit cleaning", "Input",
		"sco0", 2.4, 2, "sco0\npattern rule\nkept for comparison", "Rule",
		"sco1", 3.8, 2, "sco1\nsimple rule\ncontains 4 => 0", "Rule",
		"sco2", 5.2, 2, "sco2\nadvanced rule\npattern rescue", "Rule",
		"sco3", 6.6, 2, "sco3\nexpert labels\nML learning", "ML",
		"output", 8.0, 2, "low/middle/high\nphone.sco.summary.xlsx\nFigS1", "Output"
	)
	plot_design_flow(ml_phone_nodes, title = "ml_phone design: phone-number luckiness score", subtitle = "From transparent rules to pattern rescue and expert-label-based machine learning", file = "design_ml_phone.png")

	ml_geo_nodes <- tibble::tribble(
		~id, ~x, ~y, ~label, ~type,
		"input", 1, 2, "Pickup address\naddress type", "Input",
		"clean", 2.6, 2, "Cleaning/truncation\nmissingness handling", "Rule",
		"model", 4.2, 2, "HFL transformer\naddress classification", "ML",
		"validate", 5.8, 2, "Residential records\nyearly housing linkage", "Validation",
		"output", 7.4, 2, "geo.type1\nFig4 housing\ngeo summary", "Output"
	)
	plot_design_flow(ml_geo_nodes, title = "ml_geo design: geolocation classification and validation", subtitle = "Address classification followed by yearly housing and environmental linkage", file = "design_ml_geo.png")

	ml_dx_nodes <- tibble::tribble(
		~id, ~x, ~y, ~label, ~type,
		"input", 1, 2, "Call reason/chief complaint/history\ninitial and supplemental diagnosis", "Input",
		"kw", 2.7, 2, "dx.type0\nkeyword rules\ninterpretable QC", "Rule",
		"hfl", 4.4, 2, "dx.type1\nHFL transformer\nmain classifier", "ML",
		"map", 6.1, 2, "raw disease labels\nmap to six phenotypes", "Mapping",
		"output", 7.8, 2, "dx_grp\nconcordance\nFigS2", "Output"
	)
	plot_design_flow(ml_dx_nodes, title = "ml_dx design: disease-text classification", subtitle = "Keyword rules support interpretability/QC; HFL transformer provides main disease labels", file = "design_ml_dx.png")
}
