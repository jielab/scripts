dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
indir <- file.path(dir0, "data", "ukb", "phe")
maha_outdir <- Sys.getenv("MAHA_OUTDIR", unset = file.path(dir0, "analysis", "maha"))
.this_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_file <- if (length(.this_arg)) sub("^--file=", "", .this_arg[[1]]) else file.path(dir0, "scripts", "maha", "f", "MAHA_ukb.R")
.this_dir <- dirname(normalizePath(.this_file, winslash = "/", mustWork = FALSE))
.helper_dir <- Sys.getenv("MAHA_HELPER_DIR", unset = file.path(dir0, "scripts", "0f"))

.bootstrap_files <- c(
  file.path(.helper_dir, "0phe.f.R"),
  file.path(.helper_dir, "assoc.f.R"),
  file.path(.helper_dir, "pred.f.R"),
  file.path(.helper_dir, "plot.f.R"),
  file.path(.this_dir, "comm.f.R")
)
.bootstrap_missing <- .bootstrap_files[!file.exists(.bootstrap_files)]
if (length(.bootstrap_missing)) {
  cat("[UKB INPUT CHECK] Required code/helper file(s) do not exist:
", paste0("  MISSING: ", .bootstrap_missing, collapse = "
"), "
", sep = "")
  stop("UKB bootstrap input check failed.", call. = FALSE)
}

pacman::p_load(writexl, tidyverse, survival, flextable, patchwork, gtsummary, broom, cowplot)
invisible(lapply(c("0phe.f.R", "assoc.f.R", "pred.f.R", "plot.f.R"), function(f) source(file.path(.helper_dir, f))))
source(file.path(.this_dir, "comm.f.R"))

cohort_prefix <- "ukb"
maha_auxdir <- Sys.getenv("MAHA_AUXDIR", unset = file.path(maha_outdir, cohort_prefix))
dir.create(maha_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(maha_auxdir, recursive = TRUE, showWarnings = FALSE)
cohort_file <- function(path) file.path(maha_auxdir, basename(path))
pub_file <- function(path) file.path(maha_outdir, basename(path))
.save_plot_unprefixed <- save_plot
save_plot <- function(plot, filename, ...) .save_plot_unprefixed(plot, cohort_file(filename), ...)
save_pub_plot <- function(plot, filename, ...) .save_plot_unprefixed(plot, pub_file(filename), ...)
write_xlsx <- function(x, path, ...) writexl::write_xlsx(x, cohort_file(path), ...)
write_pub_xlsx <- function(x, path, ...) writexl::write_xlsx(x, pub_file(path), ...)
maha_map <- c(
	maha = "MAHA",
	maha_bal = "MAHA-balanced",
	maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy",
	maha_noprotein = "MAHA-no protein"
)

diet.lst <- c(
	medito = "MEDI-Touch", medi24 = "MEDI", dash = "DASH", mind = "MIND",
	hpdi = "hPDI", ahei = "AHEI", phdi = "PHDI", modern = "MODERN", digm = "DIGM",
	maha = "MAHA", maha_bal = "MAHA", maha_strict = "MAHA", maha_nodairy = "MAHA", maha_noprotein = "MAHA"
)

maha_ref <- "maha"
maha_label <- unname(maha_map[maha_ref])
maha_sum <- paste0("diet.", maha_ref, ".sum")
maha_pts <- paste0("diet.", maha_ref, ".pts")
maha_s100 <- paste0("diet.", maha_ref, ".s100")
maha_q5 <- paste0("diet.", maha_ref, ".q5")
maha_s100_q5 <- paste0("diet.", maha_ref, ".s100_q5")
maha_3c <- paste0("diet.", maha_ref, ".3c")
maha_hml <- paste0("diet.", maha_ref, ".hml")

diet.sum.cols <- paste0("diet.", names(diet.lst), ".sum")
diet.pts.cols <- paste0("diet.", names(diet.lst), ".pts")
diet.s100.cols <- paste0("diet.", names(diet.lst), ".s100")
diet.q5.cols <- paste0("diet.", names(diet.lst), ".q5")
diet.s100_q5.cols <- paste0("diet.", names(diet.lst), ".s100_q5")
diet.order <- c("MEDI-Touch", "MEDI", "DASH", "MIND", "hPDI", "AHEI", "PHDI", "MODERN", "DIGM", maha_label)

covs <- c("age", "sex.f", "center", "tdi", "PC1", "PC2", "bmi.pts", "bp.pts", "nonhdl.pts", "smoke.pts", "pa.pts", "sleep.pts")
diet.inc <- c(maha_ref, "dash", "mind", "medi24")
diet.inc.pts <- paste0("diet.", diet.inc, ".pts")
diet.inc.s100 <- paste0("diet.", diet.inc, ".s100")
diet.inc.q5 <- paste0("diet.", diet.inc, ".q5")
diet.inc.s100_q5 <- paste0("diet.", diet.inc, ".s100_q5")

Y.inc <- c("cvd_cad", "cvd_stroke_i", "cvd_hfail", "t2dm", "ckd", "death")
Y.cols <- setNames(c("#B97AF7", "#F26D60", "#D98C2B", "#7CAE00", "#16B9C0", "#1F78B4"), dx.lst[Y.inc])

group_pct_th <- 0.40
mk_hml <- function(x, p = 0.40) { q1 <- quantile(x, probs = p, na.rm = TRUE, names = FALSE); q2 <- quantile(x, probs = 1 - p, na.rm = TRUE, names = FALSE); case_when(is.na(x) ~ NA_character_, x <= q1 ~ "low", x >= q2 ~ "high", TRUE ~ "middle") }

get_mode <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA else names(sort(table(x), decreasing = TRUE))[1] }

make_typical_value <- function(x) { x2 <- x[!is.na(x)]; if (length(x2) == 0) NA else if (is.numeric(x)) mean(x2) else if (is.factor(x)) factor(get_mode(x2), levels = levels(x)) else if (is.character(x)) get_mode(x2) else x2[1] }

var2lab <- function(x) {
	x0 <- gsub("^diet\\.|\\.(s100_q5|s100|sum|pts|q5|3c|hml)$", "", x)
	out <- ifelse(x0 %in% names(maha_map), unname(maha_map[x0]), unname(diet.lst[x0]))
	ifelse(is.na(out), x, out)
}

zstd <- function(x) { x <- as.numeric(x); if (is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x))); as.numeric(scale(x)) }

std_10_90 <- function(x) { x <- as.numeric(x); q <- quantile(x, c(0.1, 0.9), na.rm = TRUE, names = FALSE); if (!all(is.finite(q)) || q[1] == q[2]) rep(NA_real_, length(x)) else (x - q[1]) / (q[2] - q[1]) }

score_q5 <- function(x) { x <- as.numeric(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x); out[ok] <- c(0, 25, 50, 75, 100)[dplyr::ntile(x[ok], 5)]; out }

score_0_100 <- function(x) { x <- as.numeric(x); ok <- is.finite(x); if (sum(ok) < 2 || diff(range(x[ok])) == 0) rep(NA_real_, length(x)) else { out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100; out } }

score_0_100_q5 <- function(x) { x <- as.numeric(x); ok <- is.finite(x); if (sum(ok) < 2 || diff(range(x[ok])) == 0) rep(NA_real_, length(x)) else { z <- rep(NA_real_, length(x)); z[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100; as.numeric(as.character(cut(z, c(-Inf, 20, 40, 60, 80, Inf), labels = c(0, 25, 50, 75, 100), right = FALSE))) } }

make_phewas_panel <- function(diet_name, show_x = TRUE) {
	d <- res %>% filter(Diet == diet_name) %>% left_join(phewas_xmap, by = c("category", "phenotype"))
	lab_d <- d %>% filter(is.finite(logp)) %>% slice_max(logp, n = 12, with_ties = FALSE)
	panel_title <- recode(diet_name, MAHA = "a. MAHA PheWAS", DASH = "b. DASH PheWAS", .default = diet_name)
	ggplot(d, aes(phe_x, logp, color = category)) +
		geom_point(size = 0.85, alpha = 0.85) +
		geom_hline(yintercept = bonf_line, color = "#F2A3A3", linewidth = 0.35) +
		ggrepel::geom_text_repel(
			data = lab_d, aes(label = description),
			size = 2.0, fontface = "bold", color = "black", segment.color = "black", segment.size = 0.25,
			box.padding = 0.18, point.padding = 0.10, max.overlaps = Inf, min.segment.length = 0,
			show.legend = FALSE
		) +
		scale_color_manual(values = phewas_cols, guide = "none") +
		scale_x_continuous(
			breaks = phewas_centers$phe_x,
			labels = phewas_centers$category,
			expand = expansion(mult = c(0.01, 0.01))
		) +
		labs(title = panel_title, x = if (show_x) "Phenotypes" else NULL, y = expression(-log[10](italic(p)))) +
		theme_bw(base_size = 9) +
		theme(
			plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
			axis.title = element_text(face = "bold"),
			axis.text.y = element_text(face = "bold"),
			axis.title.x = if (show_x) element_text(size = 9, margin = margin(t = 4)) else element_blank(),
			axis.text.x = if (show_x) element_text(size = 5.5, angle = -30, vjust = 1, hjust = 1, face = "bold") else element_blank(),
			axis.ticks.x = if (show_x) element_line(linewidth = 0.25) else element_blank(),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			plot.margin = margin(4, 4, 4, 4)
		)
}

fit_disc_reftrad <- function(dat, Y, trad_nm, trad_lab, covs, maha_nm = "maha", maha_lab = "MAHA") {
	trad_var <- paste0("diet.", trad_nm, ".hml"); maha_var <- paste0("diet.", maha_nm, ".hml")
	d0 <- dat %>%
		transmute(time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]], across(all_of(covs)), trad = .data[[trad_var]], maha = .data[[maha_var]]) %>%
		filter(trad %in% c("low", "high"), maha %in% c("low", "high")) %>%
		drop_na()
	g.ref <- paste0(trad_lab, " high + ", maha_lab, " low")
	g.alt <- paste0(trad_lab, " low + ", maha_lab, " high")
	d0 <- d0 %>% mutate(discord = factor(case_when(trad == "high" & maha == "low" ~ g.ref, trad == "low" & maha == "high" ~ g.alt, TRUE ~ NA_character_), levels = c(g.ref, g.alt)))
	coxph(Surv(time, event) ~ discord + ., data = d0 %>% filter(!is.na(discord)) %>% dplyr::select(time, event, all_of(covs), discord)) %>%
		broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
		filter(grepl("^discord", term)) %>%
		transmute(MAHA = maha_lab, Pattern = trad_lab, Outcome = factor(unname(dx.lst[Y]), levels = outcome_levels), contrast = paste0(g.alt, " vs ", g.ref), estimate, conf.low, conf.high, p.value)
}

fit_aic_reftrad <- function(dat, Y, trad_nm, trad_lab, covs, maha_nm = "maha", maha_lab = "MAHA") {
	d0 <- dat %>%
		transmute(
			time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]], across(all_of(covs)),
			maha = as.numeric(scale(.data[[paste0("diet.", maha_nm, ".pts")]])),
			trad = as.numeric(scale(.data[[paste0("diet.", trad_nm, ".pts")]]))
		) %>%
		drop_na()
	get_aic <- function(v) AIC(coxph(Surv(time, event) ~ ., data = d0 %>% dplyr::select(time, event, all_of(covs), all_of(v))))
	a0 <- get_aic(character(0)); am <- get_aic("maha"); at <- get_aic("trad"); ab <- get_aic(c("maha", "trad"))
	tibble(
		MAHA = maha_lab, Pattern = trad_lab, Outcome = unname(dx.lst[Y]),
		Base = a0, MAHA_model = am, Traditional = at, Both = ab,
		add_MAHA_given_trad = at - ab, add_trad_given_MAHA = am - ab,
		maha_only_information = pmax(at - ab, 0), traditional_only_information = pmax(am - ab, 0),
		net_traditional_gain = (am - ab) - (at - ab)
	)
}

plot_cuminc_joint <- function(dat1, Y_time, Y_event, Xfacet, Xline, covs, t0 = 10, Xfacetlab = NULL, Xlinelab = NULL, title = NULL, leg_position = "right", leg_direction = "vertical") {
	d <- dat1[!is.na(dat1[[Xfacet]]) & !is.na(dat1[[Xline]]), ]
	lvf <- levels(factor(d[[Xfacet]])); lvl <- levels(factor(d[[Xline]]))
	d$xf <- factor(d[[Xfacet]], levels = lvf); d$xl <- factor(d[[Xline]], levels = lvl)
	covs2 <- covs[covs %in% names(d)]
	d <- d[complete.cases(d[, c("xf", "xl", Y_time, Y_event, covs2), drop = FALSE]), ]
	form <- as.formula(paste0("survival::Surv(", Y_time, ", ", Y_event, ") ~ xf + xl", if (length(covs2)) paste0(" + ", paste(covs2, collapse = " + ")) else ""))
	fit <- survival::coxph(form, data = d)
	nd <- expand.grid(xf = factor(lvf, levels = lvf), xl = factor(lvl, levels = lvl))
	if (length(covs2)) for (v in covs2) nd[[v]] <- if (is.numeric(d[[v]])) mean(d[[v]], na.rm = TRUE) else names(which.max(table(d[[v]])))
	tt <- seq(0, t0, by = 0.1)
	res <- bind_rows(lapply(seq_len(nrow(nd)), function(i) { s <- summary(survival::survfit(fit, newdata = nd[i, , drop = FALSE]), times = tt, extend = TRUE); tibble(xf = nd$xf[i], xl = nd$xl[i], time = s$time, y = 1 - s$surv, cl = 1 - s$upper, cu = 1 - s$lower) }))
	ggplot(res, aes(time, y, color = xl, fill = xl)) +
		geom_ribbon(aes(ymin = cl, ymax = cu), alpha = 0.15, linewidth = 0) +
		geom_line(linewidth = 1.1) +
		facet_wrap(~xf, nrow = 1, strip.position = "bottom") +
		scale_x_continuous(limits = c(0, t0), breaks = c(0, 5, 10), labels = c("0", "5", "10")) +
		scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
		scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") + theme_minimal() +
		labs(title = title, x = if (is.null(Xfacetlab)) "Follow-up time (years)" else Xfacetlab, y = "Cumulative incidence", color = if (is.null(Xlinelab)) Xline else Xlinelab, fill = if (is.null(Xlinelab)) Xline else Xlinelab) +
		guides(color = guide_legend(nrow = 1, byrow = TRUE), fill = guide_legend(nrow = 1, byrow = TRUE)) +
		theme(
			plot.title = element_text(hjust = 0.5, face = "bold"),
			axis.title = element_text(face = "bold"),
			axis.title.x = element_text(size = 12),
			axis.text = element_text(size = 11, face = "bold"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			legend.text = element_text(size = 10, face = "bold", color = "grey30", margin = margin(l = 4, r = 12)),
			legend.title = element_text(size = 12, face = "bold", color = "grey30", margin = margin(r = 10)),
			legend.key.width = grid::unit(0.42, "cm"),
			legend.spacing.x = grid::unit(0.16, "cm"),
			legend.position = leg_position,
			legend.direction = leg_direction,
			legend.justification = "center",
			strip.placement = "outside",
			strip.background = element_blank(),
			strip.text = element_text(size = 10, face = "bold", color = "grey30")
		)
}

extract_joint_risk <- function(dat1, Y_time, Y_event, Xfacet, Xline, covs, times = 10, Xfacetlab = NULL, Xlinelab = NULL, panel = NULL) {
	d <- dat1[!is.na(dat1[[Xfacet]]) & !is.na(dat1[[Xline]]), ]
	lvf <- levels(factor(d[[Xfacet]])); lvl <- levels(factor(d[[Xline]]))
	d$xf <- factor(d[[Xfacet]], levels = lvf); d$xl <- factor(d[[Xline]], levels = lvl)
	covs2 <- covs[covs %in% names(d)]
	d <- d[complete.cases(d[, c("xf", "xl", Y_time, Y_event, covs2), drop = FALSE]), ]
	form <- as.formula(paste0("survival::Surv(", Y_time, ", ", Y_event, ") ~ xf + xl", if (length(covs2)) paste0(" + ", paste(covs2, collapse = " + ")) else ""))
	fit <- survival::coxph(form, data = d)
	nd <- expand.grid(xf = factor(lvf, levels = lvf), xl = factor(lvl, levels = lvl))
	if (length(covs2)) for (v in covs2) {
		x <- d[[v]]
		if (is.numeric(x)) nd[[v]] <- mean(x, na.rm = TRUE)
		else if (is.factor(x)) nd[[v]] <- factor(names(which.max(table(x))), levels = levels(x))
		else nd[[v]] <- names(which.max(table(x)))
	}
	cnt <- d %>%
		group_by(xf, xl) %>%
		summarise(N = n(), Events = sum(.data[[Y_event]], na.rm = TRUE), .groups = "drop")
	res <- bind_rows(lapply(seq_len(nrow(nd)), function(i) {
		s <- summary(survival::survfit(fit, newdata = nd[i, , drop = FALSE]), times = times, extend = TRUE)
		tibble(
			Facet_var = ifelse(is.null(Xfacetlab), Xfacet, Xfacetlab),
			Facet_group = as.character(nd$xf[i]),
			Line_var = ifelse(is.null(Xlinelab), Xline, Xlinelab),
			Line_group = as.character(nd$xl[i]),
			time = s$time,
			risk = 1 - s$surv,
			risk_lower = 1 - s$upper,
			risk_upper = 1 - s$lower
		)
	})) %>%
		left_join(cnt %>% mutate(Facet_group = as.character(xf), Line_group = as.character(xl)) %>% dplyr::select(Facet_group, Line_group, N, Events), by = c("Facet_group", "Line_group")) %>%
		mutate(
			Panel = panel,
			Outcome = dx.lst[Y],
			risk_pct = 100 * risk,
			risk_lower_pct = 100 * risk_lower,
			risk_upper_pct = 100 * risk_upper
		) %>%
		dplyr::select(Panel, Outcome, Facet_var, Facet_group, Line_var, Line_group, time, N, Events, risk, risk_lower, risk_upper, risk_pct, risk_lower_pct, risk_upper_pct)
	res
}

ap_effect_from_assoc <- function(assoc_tbl, source_label, margins = c(1.05, 1.075, 1.10), priors = c(1.05, 1.10, 1.20)) {
	num_col <- function(d, nm) if (!nm %in% names(d)) rep(NA_real_, nrow(d)) else suppressWarnings(as.numeric(as.character(d[[nm]])))
	chr_col <- function(d, nm) if (!nm %in% names(d)) rep(NA_character_, nrow(d)) else as.character(d[[nm]])

	hr0 <- num_col(assoc_tbl, "estimate"); lo0 <- num_col(assoc_tbl, "conf.low"); hi0 <- num_col(assoc_tbl, "conf.high")
	b0 <- num_col(assoc_tbl, "beta"); se0 <- num_col(assoc_tbl, "se"); p0 <- num_col(assoc_tbl, "p.value")

	base <- tibble(
		Source = source_label,
		Outcome = chr_col(assoc_tbl, "Outcome"),
		Diet = chr_col(assoc_tbl, "Diet"),
		beta = ifelse(is.finite(hr0) & hr0 > 0, log(hr0), b0),
		se = ifelse(is.finite(lo0) & lo0 > 0 & is.finite(hi0) & hi0 > lo0, (log(hi0) - log(lo0)) / (2 * 1.96), se0),
		P = p0,
		N_total = num_col(assoc_tbl, "N_total"),
		N_event = num_col(assoc_tbl, "N_event")
	) %>%
		mutate(HR = exp(beta), LCI95 = exp(beta - 1.96 * se), UCI95 = exp(beta + 1.96 * se)) %>%
		filter(is.finite(beta), is.finite(se), se > 0, is.finite(HR), HR > 0)

	purrr::map_dfr(seq_len(nrow(base)), function(i) {
		d <- base[i, ]
		tost <- purrr::map_dfr(margins, \(m) ap_tost_one(d$beta, d$se, m))
		pow <- purrr::map_dfr(margins, \(m) ap_power_one(d$se, m)) %>%
			dplyr::select(margin_hr, MDE_logHR, MDE_HR_lower, MDE_HR_upper, powered_for_margin)
		bf <- tibble(
			prior_95_hr = priors,
			BF01 = vapply(priors, \(p) ap_bf01_one(d$beta, d$se, p), numeric(1)),
			BF10 = 1 / BF01,
			BF_interpretation = case_when(
				BF01 >= 3 ~ "evidence_for_point_null",
				BF01 <= 1 / 3 ~ "evidence_against_point_null",
				TRUE ~ "inconclusive"
			)
		)

		tidyr::crossing(tost %>% left_join(pow, by = "margin_hr"), bf) %>%
			mutate(
				Source = d$Source, Outcome = d$Outcome, Diet = d$Diet,
				N_total = d$N_total, N_event = d$N_event,
				beta = d$beta, se = d$se, HR = d$HR, LCI95 = d$LCI95, UCI95 = d$UCI95, P = d$P
			) %>%
			dplyr::select(Source, Outcome, Diet, N_total, N_event, beta, se, HR, LCI95, UCI95, P, everything())
	})
}

ap_fit_head2head_ukb <- function(dat, Y, trad_nm, trad_lab, covs, maha_nm = "maha", maha_lab = "MAHA") {
	maha_var <- paste0("diet.", maha_nm, ".pts")
	trad_var <- paste0("diet.", trad_nm, ".pts")
	need <- c(paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), maha_var, trad_var, covs)
	if (!all(need %in% names(dat))) return(tibble())

	d0 <- dat %>%
		transmute(
			time = .data[[paste0(Y, ".t2e")]],
			event = .data[[paste0(Y, ".Yt2e")]],
			maha = as.numeric(scale(.data[[maha_var]])),
			trad = as.numeric(scale(.data[[trad_var]])),
			across(all_of(covs))
		) %>%
		drop_na()

	if (nrow(d0) < 200 || sum(d0$event == 1) < 20) return(tibble())

	f0 <- as.formula(paste0("Surv(time, event) ~ ", paste(covs, collapse = " + ")))
	fm <- as.formula(paste0("Surv(time, event) ~ maha + ", paste(covs, collapse = " + ")))
	ft <- as.formula(paste0("Surv(time, event) ~ trad + ", paste(covs, collapse = " + ")))
	fb <- as.formula(paste0("Surv(time, event) ~ maha + trad + ", paste(covs, collapse = " + ")))

	fit0 <- coxph(f0, data = d0, ties = "efron")
	fitm <- coxph(fm, data = d0, ties = "efron")
	fitt <- coxph(ft, data = d0, ties = "efron")
	fitb <- coxph(fb, data = d0, ties = "efron")

	co <- coef(fitb); vv <- vcov(fitb)
	if (!all(c("maha", "trad") %in% names(co))) return(tibble())

	diff <- unname(co["maha"] - co["trad"])
	se_diff <- sqrt(vv["maha", "maha"] + vv["trad", "trad"] - 2 * vv["maha", "trad"])
	lr_maha_given_trad <- as.numeric(2 * (logLik(fitb) - logLik(fitt)))
	lr_trad_given_maha <- as.numeric(2 * (logLik(fitb) - logLik(fitm)))

	tibble(
		Outcome = unname(dx.lst[Y]),
		Comparison = paste0(maha_lab, " vs ", trad_lab),
		N_total = nrow(d0),
		N_event = sum(d0$event == 1),
		HR_maha_mutual = exp(unname(co["maha"])),
		HR_trad_mutual = exp(unname(co["trad"])),
		logHR_difference_maha_minus_trad = diff,
		HR_ratio_maha_vs_trad = exp(diff),
		SE_difference = se_diff,
		P_equal_coefficients = 2 * pnorm(-abs(diff / se_diff)),
		AIC_base = AIC(fit0),
		AIC_maha = AIC(fitm),
		AIC_trad = AIC(fitt),
		AIC_both = AIC(fitb),
		LRT_P_add_MAHA_given_trad = pchisq(lr_maha_given_trad, df = 1, lower.tail = FALSE),
		LRT_P_add_trad_given_MAHA = pchisq(lr_trad_given_maha, df = 1, lower.tail = FALSE)
	)
}

ap_clean_equiv <- function(x, dataset_label, margin = ap_margin_main, prior = ap_prior_main) {
	x %>%
		filter(abs(margin_hr - margin) < 1e-8, abs(prior_95_hr - prior) < 1e-8) %>%
		mutate(
			Dataset = dataset_label,
			Outcome = factor(as.character(Outcome), levels = outcome_levels),
			Diet = factor(as.character(Diet), levels = diet_levels.ab),
			HR_95CI = ap_ci_lab(HR, LCI95, UCI95),
			HR_90CI = ap_ci_lab(HR, CI90_low_HR, CI90_high_HR),
			TOST_P_label = ap_p_lab(TOST_p),
			BF01_label = ap_bf_lab(BF01),
			MDE80_HR_label = sprintf("%.2f-%.2f", MDE_HR_lower, MDE_HR_upper),
			Equivalence_5pct = factor(
				ifelse(equivalent, "Equivalent within HR 0.95-1.05", "Not equivalent/inconclusive"),
				levels = c("Equivalent within HR 0.95-1.05", "Not equivalent/inconclusive")
			),
			Powered_5pct = factor(
				ifelse(powered_for_margin, "Powered for HR 1.05", "Not powered for HR 1.05"),
				levels = c("Powered for HR 1.05", "Not powered for HR 1.05")
			),
			BF_summary = case_when(
				BF01 >= 3 ~ "BF supports point null",
				BF01 <= 1/3 ~ "BF against point null",
				TRUE ~ "BF inconclusive"
			)
		) %>%
		arrange(Outcome, Diet)
}

ap_clean_head2head <- function(x, dataset_label) {
	x %>%
		mutate(
			Dataset = dataset_label,
			Outcome = factor(as.character(Outcome), levels = outcome_levels),
			Comparison = factor(as.character(Comparison), levels = c("MAHA vs DASH", "MAHA vs MIND", "MAHA vs MEDI")),
			Comp2 = factor(gsub("^MAHA vs ", "", as.character(Comparison)), levels = c("DASH", "MIND", "MEDI")),
			LCI_ratio = exp(logHR_difference_maha_minus_trad - 1.96 * SE_difference),
			UCI_ratio = exp(logHR_difference_maha_minus_trad + 1.96 * SE_difference),
			Ratio_95CI = ap_ci_lab(HR_ratio_maha_vs_trad, LCI_ratio, UCI_ratio),
			P_equal_label = ap_p_lab(P_equal_coefficients)
		) %>%
		arrange(Outcome, Comparison)
}

ap_plot_main <- function(d, title, margin = ap_margin_main) {
	shift <- c(MAHA = -0.27, DASH = -0.09, MIND = 0.09, MEDI = 0.27)
	lab_layout <- tibble(
		Diet = factor(c("MAHA", "DASH", "MIND", "MEDI"), levels = diet_levels.ab),
		lab_col = c(1, 2, 1, 2),
		lab_y_shift = c(0.16, 0.16, -0.16, -0.16)
	)
	d <- d %>%
		mutate(y = as.numeric(Outcome) + unname(shift[as.character(Diet)]), lab = HR_95CI) %>%
		left_join(lab_layout, by = "Diet") %>%
		mutate(y_lab = as.numeric(Outcome) + lab_y_shift, x_lab = if_else(lab_col == 1, 1.095, 1.225))

	ggplot(d, aes(color = Diet)) +
		annotate("rect", xmin = 1 / margin, xmax = margin, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "grey70") +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
		geom_vline(xintercept = c(1 / margin, margin), linetype = 3, color = "grey55", linewidth = 0.6) +
		geom_segment(aes(x = LCI95, xend = UCI95, y = y, yend = y), linewidth = 1) +
		geom_point(aes(x = HR, y = y), size = 2.75) +
		geom_text(aes(x = x_lab, y = y_lab, label = lab, color = Diet), hjust = 0, size = 3.05, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = diet_cols.ab, breaks = diet_levels.ab) +
		scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
		coord_cartesian(xlim = c(min(0.77, d$LCI95, 1 / margin, na.rm = TRUE), 1.32), clip = "off") +
		labs(title = title, x = "Hazard ratio", y = NULL, color = NULL) +
		ap_theme(14) +
		theme(
			legend.position = "bottom",
			axis.text.y = element_text(face = "bold", size = 13),
			plot.margin = margin(5.5, 34, 5.5, 36)
		)
}

ap_plot_head2head <- function(d, title) {
	shift <- c(DASH = -0.22, MIND = 0, MEDI = 0.22)
	d <- d %>% mutate(y = as.numeric(Outcome) + unname(shift[as.character(Comp2)]), lab = paste0(Ratio_95CI, "; Pdiff=", P_equal_label))
	comp_cols <- c(DASH = diet_cols.ab[["DASH"]], MIND = diet_cols.ab[["MIND"]], MEDI = diet_cols.ab[["MEDI"]])

	ggplot(d, aes(color = Comp2)) +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
		geom_segment(aes(x = LCI_ratio, xend = UCI_ratio, y = y, yend = y), linewidth = 1) +
		geom_point(aes(x = HR_ratio_maha_vs_trad, y = y), size = 2.75) +
		geom_text(aes(x = 1.645, y = y, label = lab, color = Comp2), hjust = 1, size = 3.45, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = comp_cols, breaks = c("DASH", "MIND", "MEDI")) +
		scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
		coord_cartesian(xlim = c(min(0.98, d$LCI_ratio, na.rm = TRUE), 1.66), clip = "off") +
		labs(title = title, x = "HR ratio: MAHA coefficient / established-score coefficient", y = NULL, color = NULL) +
		ap_theme(14) +
		theme(
			legend.position = "bottom",
			axis.text.y = element_blank(),
			axis.ticks.y = element_line(linewidth = 0.45, color = "black"),
			axis.ticks.length.y = grid::unit(4.2, "pt"),
			plot.margin = margin(5.5, 28, 5.5, 8)
		)
}

ap_plot_maha_equiv_bf <- function(d, title, margin = ap_margin_main) {
	d <- d %>%
		filter(as.character(Diet) == "MAHA") %>%
		mutate(
			y = as.numeric(Outcome),
			lab = paste0("TOST=", TOST_P_label, "; BF=", BF01_label),
			Outcome_chr = as.character(Outcome),
			Equivalence_5pct = factor(as.character(Equivalence_5pct), levels = c("Equivalent within HR 0.95-1.05", "Not equivalent/inconclusive"))
		)

	ggplot(d, aes(color = Outcome_chr)) +
		annotate("rect", xmin = 1 / margin, xmax = margin, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "grey70") +
		geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
		geom_vline(xintercept = c(1 / margin, margin), linetype = 3, color = "grey55", linewidth = 0.6) +
		geom_segment(aes(x = CI90_low_HR, xend = CI90_high_HR, y = y, yend = y), linewidth = 1.2) +
		geom_point(aes(x = HR, y = y, shape = Equivalence_5pct), fill = "white", size = 3.4, stroke = 1.2) +
		geom_text(aes(x = 1.165, y = y, label = lab), hjust = 1, size = 3.3, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = fig3_outcome_cols, guide = "none") +
		scale_shape_manual(
			values = c("Equivalent within HR 0.95-1.05" = 16, "Not equivalent/inconclusive" = 1),
			breaks = c("Equivalent within HR 0.95-1.05", "Not equivalent/inconclusive"),
			drop = FALSE
		) +
		ap_shape_guide +
		scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
		coord_cartesian(xlim = c(min(0.88, d$CI90_low_HR, 1 / margin, na.rm = TRUE), 1.18), clip = "off") +
		labs(title = title, x = "MAHA hazard ratio with 90% CI", y = NULL, shape = NULL) +
		ap_theme(14) +
		ap_leg_theme +
		theme(
			axis.text.y = element_text(face = "bold", size = 13),
			plot.margin = margin(5.5, 28, 5.5, 36)
		)
}

ap_plot_maha_power <- function(d, title, margin = ap_margin_main) {
	d <- d %>%
		filter(as.character(Diet) == "MAHA") %>%
		mutate(
			y = as.numeric(Outcome),
			lab = MDE80_HR_label,
			Outcome_chr = as.character(Outcome),
			Powered_5pct = factor(as.character(Powered_5pct), levels = c("Powered for HR 1.05", "Not powered for HR 1.05"))
		)

	ggplot(d, aes(color = Outcome_chr)) +
		geom_vline(xintercept = margin, linetype = 3, color = "grey55", linewidth = 0.7) +
		geom_segment(aes(x = MDE_HR_lower, xend = MDE_HR_upper, y = y, yend = y), linewidth = 1.2) +
		geom_point(aes(x = MDE_HR_upper, y = y, shape = Powered_5pct), fill = "white", size = 3.4, stroke = 1.2) +
		geom_text(aes(x = 1.112, y = y, label = lab), hjust = 1, size = 3.3, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = fig3_outcome_cols, guide = "none") +
		scale_shape_manual(
			values = c("Powered for HR 1.05" = 16, "Not powered for HR 1.05" = 1),
			breaks = c("Powered for HR 1.05", "Not powered for HR 1.05"),
			drop = FALSE
		) +
		ap_shape_guide +
		scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
		coord_cartesian(xlim = c(0.90, 1.12), clip = "off") +
		labs(title = title, x = "Minimum detectable HR at 80% power", y = NULL, shape = NULL) +
		ap_theme(14) +
		ap_leg_theme +
		theme(
			axis.text.y = element_blank(),
			axis.ticks.y = element_line(linewidth = 0.45, color = "black"),
			axis.ticks.length.y = grid::unit(4.2, "pt"),
			plot.margin = margin(5.5, 28, 5.5, 8)
		)
}

ap_make_fig3_outputs <- function(equiv_main_tbl, head_tbl, dataset_label, out_prefix = "Fig3.precision_attenuation", equiv_alt_tbl = NULL) {
	ed <- ap_clean_equiv(equiv_main_tbl, dataset_label)
	hd <- ap_clean_head2head(head_tbl, dataset_label)

	p1 <- ap_plot_main(ed, "a. Main associations with small-effect zone")
	p2 <- ap_plot_head2head(hd, "b. Direct attenuation versus established scores")
	p3 <- ap_plot_maha_equiv_bf(ed, "c. MAHA equivalence and Bayes evidence")
	p4 <- ap_plot_maha_power(ed, "d. MAHA detectable effect")

	fig_gap_x <- ggplot() + theme_void()
	fig_gap_y <- ggplot() + theme_void()

	fig <- cowplot::plot_grid(
		cowplot::plot_grid(p1, fig_gap_x, p2, nrow = 1, rel_widths = c(1.02, 0.055, 1.00), align = "none"),
		fig_gap_y,
		cowplot::plot_grid(p3, fig_gap_x, p4, nrow = 1, rel_widths = c(1.02, 0.055, 1.00), align = "none"),
		ncol = 1,
		rel_heights = c(1.00, 0.055, 1.00),
		align = "none"
	)

	save_plot(fig, paste0(out_prefix, ".png"), width = 16.0, height = 11.8, dpi = 320, bg = "white")

	write_xlsx(
		list(
			Fig3_main_summary = ed %>%
				dplyr::select(Dataset, Outcome, Diet, N_total, N_event, HR_95CI, HR_90CI, P, margin_hr, TOST_p, Equivalence_5pct, BF01, BF_summary, MDE80_HR_label, Powered_5pct),
			Fig3_MAHA_summary = ed %>%
				filter(as.character(Diet) == "MAHA") %>%
				dplyr::select(Dataset, Outcome, Diet, N_total, N_event, HR_95CI, HR_90CI, P, margin_hr, TOST_p, Equivalence_5pct, BF01, BF_summary, MDE80_HR_label, Powered_5pct),
			Fig3_head2head = hd %>%
				dplyr::select(Dataset, Outcome, Comparison, N_total, N_event, HR_maha_mutual, HR_trad_mutual, HR_ratio_maha_vs_trad, LCI_ratio, UCI_ratio, P_equal_coefficients, Ratio_95CI),
			Fig3_raw_main = equiv_main_tbl,
			Fig3_raw_altMAHA = if (is.null(equiv_alt_tbl) || nrow(equiv_alt_tbl) == 0) tibble() else equiv_alt_tbl,
			Fig3_raw_head2head = head_tbl
		),
		paste0(out_prefix, ".out.xlsx")
	)

	invisible(list(equiv = ed, head = hd, fig = fig))
}

z1 <- function(x) as.numeric(scale(as.numeric(x)))
r2 <- function(fit, adj = FALSE) if (adj) summary(fit)$adj.r.squared else summary(fit)$r.squared
cor1 <- function(x, y) tryCatch({
	d <- complete.cases(x, y)
	if (sum(d) < 3 || sd(x[d]) == 0 || sd(y[d]) == 0) return(NA_real_)
	cor(x[d], y[d], method = "spearman")
}, error = \(e) NA_real_)

cor_p <- function(x, y) tryCatch({
	d <- complete.cases(x, y)
	if (sum(d) < 3 || sd(x[d]) == 0 || sd(y[d]) == 0) return(NA_real_)
	cor.test(x[d], y[d], method = "spearman", exact = FALSE)$p.value
}, error = \(e) NA_real_)

decomp_one <- function(d, diet_lab) {
	d0 <- d %>% filter(Diet == diet_lab) %>% mutate(score_z = z1(score))
	m1 <- lm(score_z ~ tdi_z, d0); m2 <- lm(score_z ~ center_name, d0); m3 <- lm(score_z ~ spatial_tile, d0)
	m4 <- lm(score_z ~ tdi_z + center_name, d0); m5 <- lm(score_z ~ tdi_z + spatial_tile, d0)
	tibble(
		Diet = diet_lab, N = nrow(d0),
		R2_TDI_only = r2(m1), R2_center_only = r2(m2), R2_grid_only = r2(m3),
		adjR2_grid_only = r2(m3, TRUE),
		R2_TDI_plus_center = r2(m4), R2_TDI_plus_grid = r2(m5),
		unique_TDI_given_center = r2(m4) - r2(m2),
		unique_center_given_TDI = r2(m4) - r2(m1),
		unique_TDI_given_grid = r2(m5) - r2(m3),
		unique_grid_given_TDI = r2(m5) - r2(m1)
	)
}

get_wkappa <- function(x, y) {
	d0 <- data.frame(
		x = factor(x, levels = c(0, 25, 50, 75, 100), ordered = TRUE),
		y = factor(y, levels = c(0, 25, 50, 75, 100), ordered = TRUE)
	) %>% drop_na()
	if (nrow(d0) < 50) return(NA_real_)
	irr::kappa2(d0, weight = "squared")$value
}

agree_stat <- function(x, y) {
	d0 <- data.frame(x = as.numeric(x), y = as.numeric(y)) %>% drop_na()
	if (nrow(d0) < 50) return(c(n = nrow(d0), exact = NA, within1 = NA))
	c(n = nrow(d0), exact = mean(d0$x == d0$y), within1 = mean(abs(d0$x - d0$y) <= 25))
}

run_joint_risk_one <- function(dat, Y, covs, t0 = 10) {
	d2 <- dat %>% transmute(time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]], across(all_of(covs)), dash = factor(diet.dash.3c, levels = c("low", "middle", "high")), maha = factor(.data[[maha_3c]], levels = c("low", "middle", "high"))) %>% drop_na()
	fit <- coxph(Surv(time, event) ~ dash * maha + ., data = d2 %>% dplyr::select(time, event, dash, maha, all_of(covs)))
	nd <- crossing(dash = factor(c("low", "middle", "high"), levels = c("low", "middle", "high")), maha = factor(c("low", "middle", "high"), levels = c("low", "middle", "high"))); for (v in covs) nd[[v]] <- make_typical_value(d2[[v]]); stopifnot(!anyNA(nd))
	bind_rows(lapply(seq_len(nrow(nd)), function(i) { ss <- summary(survfit(fit, newdata = nd[i, , drop = FALSE]), times = t0, extend = TRUE); data.frame(Y = Y, dash = as.character(nd$dash[i]), maha = as.character(nd$maha[i]), risk10 = round(100 * (1 - ss$surv), 2)) }))
}

plot_joint_heat <- function(d, ttl) {
	d <- d %>% mutate(dash = factor(dash, levels = c("low", "middle", "high")), maha = factor(maha, levels = c("low", "middle", "high")))
	ggplot(d, aes(x = maha, y = dash, fill = risk10)) +
		geom_tile(color = "white", linewidth = 0.8) +
		geom_text(aes(label = sprintf("%.2f%%", risk10)), size = 4.0) +
		scale_x_discrete(labels = c(low = "Low", middle = "Middle", high = "High")) +
		scale_y_discrete(labels = c(low = "Low", middle = "Middle", high = "High")) +
		scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
		labs(title = ttl, x = "MAHA adherence", y = "DASH adherence", fill = "10-year risk") +
		theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.margin = margin(10, 5.5, 10, 5.5))
}

collapse1 <- function(dat, v, per = 0.8) {
	cc <- grep(paste0("^", v, "_i[0-4]$"), names(dat), value = TRUE)
	if (length(cc) == 0) rep(NA_real_, nrow(dat)) else rowMeans1(dat[, cc, drop = FALSE], per = per)
}

mk_ci <- function(b, se) c(est = exp(b), low = exp(b - 1.96 * se), high = exp(b + 1.96 * se))

fit_disc <- function(v, dat.item, covs) {
	d0 <- dat.item %>%
		transmute(y_disc, x0 = suppressWarnings(as.numeric(.data[[v]])), across(any_of(covs))) %>%
		filter(!is.na(y_disc), is.finite(x0)) %>%
		drop_na(any_of(covs)) %>%
		mutate(across(where(is.factor), droplevels))

	if (nrow(d0) < 200 || dplyr::n_distinct(d0$x0) < 5 || min(table(d0$y_disc)) < 50) return(NULL)

	xbin <- d0$x0 > 0
	tab <- table(d0$y_disc, xbin)
	if (length(tab) < 4 || sum(xbin) < 20 || min(tab) < 10) return(NULL)

	d0$x <- zstd(d0$x0)
	if (!any(is.finite(d0$x)) || is.na(sd(d0$x, na.rm = TRUE)) || sd(d0$x, na.rm = TRUE) == 0) return(NULL)

	fit <- tryCatch(
		suppressWarnings(glm(as.formula(paste0("y_disc ~ x + ", paste(covs, collapse = " + "))), family = binomial(), data = d0)),
		error = function(e) NULL
	)
	if (is.null(fit)) return(NULL)

	tb <- broom::tidy(fit) %>% filter(term == "x")
	if (nrow(tb) == 0 || any(!is.finite(c(tb$estimate, tb$std.error, tb$p.value)))) return(NULL)

	ci <- mk_ci(tb$estimate, tb$std.error)
	tibble(item = v, OR = ci["est"], conf.low = ci["low"], conf.high = ci["high"], p.value = tb$p.value, beta = tb$estimate, n = nrow(d0))
}

fit_cox <- function(v, Y, dat.item0, covs) {
	d0 <- dat %>%
		transmute(eid = as.character(eid), time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]], across(any_of(covs))) %>%
		left_join(dat.item0 %>% dplyr::select(eid, all_of(v)), by = "eid") %>%
		transmute(time, event, x0 = suppressWarnings(as.numeric(.data[[v]])), across(any_of(covs))) %>%
		filter(!is.na(time), !is.na(event), is.finite(x0)) %>%
		drop_na(any_of(covs)) %>%
		mutate(across(where(is.factor), droplevels))

	if (nrow(d0) < 200 || sum(d0$event == 1, na.rm = TRUE) < 20 || dplyr::n_distinct(d0$x0) < 5) return(NULL)

	xbin <- d0$x0 > 0
	if (sum(xbin) < 20 || sum(d0$event == 1 & xbin) < 5 || sum(d0$event == 1 & !xbin) < 5) return(NULL)

	d0$x <- zstd(d0$x0)
	if (!any(is.finite(d0$x)) || is.na(sd(d0$x, na.rm = TRUE)) || sd(d0$x, na.rm = TRUE) == 0) return(NULL)

	fit <- tryCatch(
		suppressWarnings(coxph(as.formula(paste0("Surv(time, event) ~ x + ", paste(covs, collapse = " + "))), data = d0, ties = "efron", singular.ok = TRUE)),
		error = function(e) NULL
	)
	if (is.null(fit)) return(NULL)

	tb <- broom::tidy(fit) %>% filter(term == "x")
	if (nrow(tb) == 0 || any(!is.finite(c(tb$estimate, tb$std.error, tb$p.value)))) return(NULL)

	ci <- mk_ci(tb$estimate, tb$std.error)
	tibble(item = v, Y = Y, Outcome = dx.lst[Y], HR = ci["est"], conf.low = ci["low"], conf.high = ci["high"], p.value = tb$p.value, beta = tb$estimate)
}

ukb_component_keys <- c("protein", "dairy", "veg", "fruit", "wholegrain", "fat", "upf", "alcohol", "sodium")
ukb_component_cols <- paste0(ukb_component_keys, ".maha")
ukb_component_labels <- c(
  protein = "Protein", dairy = "Dairy", veg = "Vegetables", fruit = "Fruit",
  wholegrain = "Whole grains", fat = "Healthy fats", upf = "Low UPF/refined foods",
  alcohol = "Low alcohol", sodium = "Low sodium"
)
ukb_jobs <- suppressWarnings(as.integer(Sys.getenv("MAHA_JOBS", unset = ifelse(.Platform$OS.type == "windows", "1", "4"))))
if (!is.finite(ukb_jobs) || ukb_jobs < 1) ukb_jobs <- 1L
ukb_null_max <- suppressWarnings(as.integer(Sys.getenv("MAHA_NULL_MAX", unset = "0")))
if (!is.finite(ukb_null_max) || ukb_null_max < 0) ukb_null_max <- 0L

ukb_make_sign_grid <- function(p, max_null = 0L) {
  g <- expand.grid(rep(list(c(-1, 1)), p), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  names(g) <- ukb_component_keys[seq_len(p)]
  obs <- which(rowSums(g == 1) == p)
  stopifnot(length(obs) == 1L)
  nul <- setdiff(seq_len(nrow(g)), obs)
  if (max_null > 0L && length(nul) > max_null) {
    pick <- unique(round(seq(1, length(nul), length.out = max_null)))
    nul <- nul[pick]
  }
  g[c(obs, nul), , drop = FALSE]
}

ukb_parallel_lapply <- function(X, FUN) {
  if (.Platform$OS.type != "windows" && ukb_jobs > 1L) parallel::mclapply(X, FUN, mc.cores = ukb_jobs, mc.preschedule = TRUE) else lapply(X, FUN)
}

ukb_validate_maha <- function(dat) {
  step_header("Reviewer validation before UKB main analysis: MAHA is not an arbitrary score")
  need <- c("eid", "death.t2e", "death.Yt2e", covs, ukb_component_cols,
            "diet.dash.sum", "diet.mind.sum", "diet.medi24.sum")
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("UKB validation missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  d0 <- dat %>%
    dplyr::select(all_of(need)) %>%
    filter(is.finite(death.t2e), death.t2e > 0) %>%
    drop_na()
  if (nrow(d0) < 1000 || sum(d0$death.Yt2e == 1) < 100) {
    stop("UKB validation sample is unexpectedly small: N=", nrow(d0), "; deaths=", sum(d0$death.Yt2e == 1), call. = FALSE)
  }
  message("[UKB VALIDATE] common complete-case N=", nrow(d0), "; deaths=", sum(d0$death.Yt2e == 1))

  X <- as.matrix(d0[, ukb_component_cols, drop = FALSE]); storage.mode(X) <- "double"
  primary <- rowMeans(X)
  estz <- sapply(c("diet.dash.sum", "diet.mind.sum", "diet.medi24.sum"), function(v) zstd(d0[[v]]))
  anchor <- rowMeans(estz)
  patterns <- ukb_make_sign_grid(ncol(X), ukb_null_max)

  fit_pattern <- function(i) {
    sg <- as.numeric(patterns[i, ])
    raw_score <- as.numeric(X %*% sg) / ncol(X)
    dd <- d0; dd$score <- zstd(raw_score)
    fm <- as.formula(paste0("Surv(death.t2e, death.Yt2e) ~ score + ", paste(covs, collapse = " + ")))
    fit <- tryCatch(coxph(fm, data = dd, ties = "efron", singular.ok = TRUE, model = FALSE, x = FALSE, y = FALSE), error = function(e) NULL)
    beta <- se <- p <- NA_real_
    if (!is.null(fit)) {
      tt <- broom::tidy(fit) %>% filter(term == "score")
      if (nrow(tt)) { beta <- tt$estimate[1]; se <- tt$std.error[1]; p <- tt$p.value[1] }
    }
    tibble(
      pattern_id = i, pattern = paste(ifelse(sg > 0, "+", "-"), collapse = ""),
      is_prespecified_MAHA = all(sg == 1), N = nrow(dd), deaths = sum(dd$death.Yt2e == 1),
      beta = beta, se = se, HR = exp(beta), CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se), P = p,
      protective_Z = -beta / se,
      rho_established_anchor = suppressWarnings(cor(raw_score, anchor, method = "spearman", use = "complete.obs"))
    )
  }
  message("[UKB VALIDATE] fitting ", nrow(patterns), " matched direction patterns; workers=", ukb_jobs)
  pattern_res <- bind_rows(ukb_parallel_lapply(seq_len(nrow(patterns)), fit_pattern))
  obs <- pattern_res %>% filter(is_prespecified_MAHA)
  nul <- pattern_res %>% filter(!is_prespecified_MAHA)
  summary_tbl <- tibble(
    Cohort = "UK Biobank", N = obs$N, Deaths = obs$deaths,
    MAHA_HR_per_SD = obs$HR, MAHA_CI_low = obs$CI_low, MAHA_CI_high = obs$CI_high,
    Direction_null_n = nrow(nul), MAHA_protective_Z = obs$protective_Z,
    Mortality_percentile_vs_direction_null = mean(nul$protective_Z <= obs$protective_Z, na.rm = TRUE),
    Mortality_empirical_one_sided_P = (1 + sum(nul$protective_Z >= obs$protective_Z, na.rm = TRUE)) / (1 + sum(is.finite(nul$protective_Z))),
    MAHA_anchor_Spearman = obs$rho_established_anchor,
    Anchor_percentile_vs_direction_null = mean(nul$rho_established_anchor <= obs$rho_established_anchor, na.rm = TRUE),
    Anchor_empirical_one_sided_P = (1 + sum(nul$rho_established_anchor >= obs$rho_established_anchor, na.rm = TRUE)) / (1 + sum(is.finite(nul$rho_established_anchor)))
  )

  loo_specs <- c("Primary MAHA", paste0("Without ", unname(ukb_component_labels)))
  loo_drop <- c(NA_character_, ukb_component_keys)
  loo <- purrr::map2_dfr(loo_specs, loo_drop, function(lab, drop_key) {
    keep <- if (is.na(drop_key)) ukb_component_keys else setdiff(ukb_component_keys, drop_key)
    cols <- paste0(keep, ".maha")
    raw_score <- rowMeans(as.matrix(d0[, cols, drop = FALSE]))
    dd <- d0; dd$score <- zstd(raw_score)
    fm <- as.formula(paste0("Surv(death.t2e, death.Yt2e) ~ score + ", paste(covs, collapse = " + ")))
    fit <- coxph(fm, data = dd, ties = "efron", singular.ok = TRUE, model = FALSE, x = FALSE, y = FALSE)
    tt <- broom::tidy(fit) %>% filter(term == "score")
    tibble(label = lab, omitted = ifelse(is.na(drop_key), "None", drop_key), N = nrow(dd), deaths = sum(dd$death.Yt2e == 1),
           beta = tt$estimate[1], se = tt$std.error[1], HR = exp(tt$estimate[1]),
           CI_low = exp(tt$estimate[1] - 1.96 * tt$std.error[1]), CI_high = exp(tt$estimate[1] + 1.96 * tt$std.error[1]),
           P = tt$p.value[1], rho_primary = suppressWarnings(cor(raw_score, primary, method = "spearman")))
  })

  p1 <- ggplot(nul, aes(protective_Z)) +
    geom_histogram(bins = 35, boundary = 0) + geom_vline(xintercept = obs$protective_Z, linewidth = 0.9) +
    labs(title = "a. Mortality criterion validity", subtitle = "511 matched scores use the same 9 components; only directions differ",
         x = "Protective Z statistic (-beta/SE)", y = "Matched-null scores") + theme_classic(base_size = 11)
  p2 <- ggplot(nul, aes(rho_established_anchor)) +
    geom_histogram(bins = 35, boundary = 0) + geom_vline(xintercept = obs$rho_established_anchor, linewidth = 0.9) +
    labs(title = "b. Convergent construct validity", subtitle = "Anchor = mean z(DASH, MIND, Mediterranean diet)",
         x = "Spearman rho with established-diet anchor", y = "Matched-null scores") + theme_classic(base_size = 11)
  loo_plot <- loo %>% mutate(label = factor(label, levels = rev(label)))
  p3 <- ggplot(loo_plot, aes(HR, label)) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45) +
    geom_segment(aes(x = CI_low, xend = CI_high, yend = label), linewidth = 0.6) + geom_point(size = 2) +
    labs(title = "c. Leave-one-component-out mortality", x = "HR per 1-SD higher score", y = NULL) + theme_classic(base_size = 11)
  p4 <- ggplot(loo_plot, aes(rho_primary, label)) + geom_segment(aes(x = 0, xend = rho_primary, yend = label), linewidth = 0.45) +
    geom_point(size = 2) + coord_cartesian(xlim = c(0, 1)) +
    labs(title = "d. Rank stability after omitting one domain", x = "Spearman rho with primary MAHA", y = NULL) + theme_classic(base_size = 11)
  fig <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "UK Biobank: prespecified MAHA versus matched arbitrary scores")
  save_plot(fig, "FigS1.validate.png", width = 14, height = 10.5, dpi = 500, bg = "white")
  write_xlsx(list(summary = summary_tbl, direction_null = pattern_res, leave_one_component_out = loo), "FigS1.validate.out.xlsx")
  print(summary_tbl)
  invisible(list(summary = summary_tbl, pattern = pattern_res, loo = loo))
}

ukb_construct_profile <- function(dat) {
  step_header("Reviewer validation before UKB main analysis: MAHA construct profile")
  need <- c(ukb_component_cols, "diet.maha.sum", "diet.dash.sum", "diet.mind.sum", "diet.medi24.sum")
  miss <- setdiff(need, names(dat)); if (length(miss)) stop("UKB construct profile missing: ", paste(miss, collapse = ", "))
  d <- dat %>% dplyr::select(all_of(need)) %>% drop_na()
  X <- as.data.frame(d[, ukb_component_cols, drop = FALSE]); names(X) <- ukb_component_keys
  X$maha <- d$diet.maha.sum
  X$decile <- dplyr::ntile(X$maha, 10)
  profile <- X %>% group_by(decile) %>% summarise(across(all_of(ukb_component_keys), mean), .groups = "drop") %>%
    pivot_longer(all_of(ukb_component_keys), names_to = "component", values_to = "mean_score") %>%
    mutate(Component = unname(ukb_component_labels[component]))
  metrics <- purrr::map_dfr(ukb_component_keys, function(k) {
    others <- setdiff(ukb_component_keys, k)
    rest <- rowMeans(as.matrix(X[, others, drop = FALSE]))
    d1 <- mean(X[[k]][X$decile == 1], na.rm = TRUE); d10 <- mean(X[[k]][X$decile == 10], na.rm = TRUE)
    tibble(component = k, Component = ukb_component_labels[[k]],
           item_rest_rho = suppressWarnings(cor(X[[k]], rest, method = "spearman", use = "complete.obs")),
           rho_primary_MAHA = suppressWarnings(cor(X[[k]], X$maha, method = "spearman", use = "complete.obs")),
           mean_decile1 = d1, mean_decile10 = d10, D10_minus_D1 = d10 - d1)
  })
  conv <- tibble(
    Comparator = c("DASH", "MIND", "Mediterranean diet", "Established-diet anchor"),
    Spearman_rho = c(
      cor(d$diet.maha.sum, d$diet.dash.sum, method = "spearman"),
      cor(d$diet.maha.sum, d$diet.mind.sum, method = "spearman"),
      cor(d$diet.maha.sum, d$diet.medi24.sum, method = "spearman"),
      cor(zstd(d$diet.maha.sum), rowMeans(cbind(zstd(d$diet.dash.sum), zstd(d$diet.mind.sum), zstd(d$diet.medi24.sum))), method = "spearman")
    )
  )
  p1 <- ggplot(profile, aes(decile, mean_score, group = Component, color = Component)) + geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    scale_x_continuous(breaks = 1:10) + labs(title = "a. All nine domains across MAHA deciles", x = "Primary MAHA decile", y = "Mean component score (0-10)", color = NULL) + theme_classic(base_size = 11)
  mplot <- metrics %>% mutate(Component = factor(Component, levels = rev(Component)))
  p2 <- ggplot(mplot, aes(D10_minus_D1, Component)) + geom_vline(xintercept = 0, linetype = 2) + geom_segment(aes(x = 0, xend = D10_minus_D1, yend = Component)) + geom_point(size = 2) +
    labs(title = "b. Component separation", x = "Mean score difference: decile 10 - decile 1", y = NULL) + theme_classic(base_size = 11)
  p3 <- ggplot(mplot, aes(item_rest_rho, Component)) + geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 2) +
    labs(title = "c. Corrected item-rest correlation", subtitle = "Descriptive only: MAHA is a formative index, not a psychometric scale", x = "Spearman rho", y = NULL) + theme_classic(base_size = 11)
  p4 <- ggplot(conv, aes(Spearman_rho, reorder(Comparator, Spearman_rho))) + geom_segment(aes(x = 0, xend = Spearman_rho, yend = reorder(Comparator, Spearman_rho))) + geom_point(size = 2.3) +
    coord_cartesian(xlim = c(0, 1)) + labs(title = "d. Convergent validity", x = "Spearman rho", y = NULL) + theme_classic(base_size = 11)
  fig <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "UK Biobank: MAHA construct profile before outcome benchmarking")
  save_plot(fig, "FigS2.construct_profile.png", width = 14, height = 10.5, dpi = 500, bg = "white")
  write_xlsx(list(decile_profile = profile, component_metrics = metrics, convergent_validity = conv), "FigS2.construct_profile.out.xlsx")
  invisible(list(profile = profile, metrics = metrics, convergent = conv))
}

setwd2(maha_outdir)

log_file <- cohort_file("maha.log")
if (file.exists(log_file)) invisible(file.remove(log_file))
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
on.exit({ while (sink.number() > 0) sink(); close(log_con) }, add = TRUE)

options(width = 200, warn = 1)

ukb_required_inputs <- c(
  all_rds = file.path(indir, "Rdata", "all.rds"),
  maha_components_rds = file.path(indir, "Rdata", "diet.maha.rds"),
  foods_dictionary = file.path(dir0, "files", "foods.xlsx"),
  diet_dictionary = file.path(indir, "common", "diet.lst"),
  diet_raw_rds = file.path(indir, "Rdata", "diet0.rds"),
  diet_fudan_rds = file.path(indir, "Rdata", "diet.fudan.rds")
)
cat("\n[UKB INPUT CHECK] Required analysis inputs:\n")
for (nm in names(ukb_required_inputs)) {
  f <- ukb_required_inputs[[nm]]
  cat(sprintf("  %-24s %s : %s%s\n", nm, if (file.exists(f)) "OK" else "MISSING", f, if (file.exists(f)) "" else " 文件不存在"))
}
ukb_missing_inputs <- ukb_required_inputs[!file.exists(ukb_required_inputs)]
if (length(ukb_missing_inputs)) {
  cat("[UKB INPUT CHECK] STOP. The following required file(s) do not exist:\n",
      paste0("  MISSING: ", unname(ukb_missing_inputs), " 文件不存在", collapse = "\n"), "\n", sep = "")
  stop("UKB required input file(s) missing. See ", log_file, call. = FALSE)
}
cat("[UKB INPUT CHECK] PASS\n\n")

maha_configure_steps(
	c("data_qc", "table1", "fig1", "fig2", "fig3", "fig5", "figS1", "figS2", "figS3", "figS4", "survdiag"),
	dependencies = list(
		table1 = "data_qc",
		fig1 = "data_qc",
		fig2 = "data_qc",
		fig3 = c("data_qc", "fig2"),
		fig5 = "data_qc",
		figS1 = "data_qc",
		figS2 = c("data_qc", "fig1"),
		figS3 = "data_qc",
		figS4 = "data_qc",
		survdiag = "data_qc"
	)
)

maha_run_step("data_qc", {
step_header("Data QC")
dat0 <- readRDS(paste0(indir, "/Rdata/all.rds"))

dat <- dat0 %>%
	filter(ethnic.c == "White") %>%
	mutate(
		across(all_of(diet.sum.cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.pts"),
		across(all_of(diet.sum.cols), score_0_100, .names = "{sub('\\\\.sum$', '', .col)}.s100"),
		across(all_of(diet.sum.cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5"),
		across(all_of(diet.sum.cols), score_0_100_q5, .names = "{sub('\\\\.sum$', '', .col)}.s100_q5"),
		across(all_of(diet.sum.cols), f3c, .names = "{sub('\\\\.sum$', '', .col)}.3c")
	)

for (Y in names(dx.lst)) {
	dat[grep(paste0("^", Y, "\\.([Y]?(t2e|r2e))$"), names(dat))] <- NULL
	dat <- t2e(dat, "cvd", paste0("fod_icd10_", Y), "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
}

dat <- dat %>% mutate(across(any_of(diet.inc.pts), ~ mk_hml(., group_pct_th), .names = "{sub('\\\\.pts$', '', .col)}.hml"))

if (!all(ukb_component_cols %in% names(dat))) {
  comp <- readRDS(ukb_required_inputs[["maha_components_rds"]]) %>% as.data.frame()
  if (!"eid" %in% names(comp)) stop("UKB diet.maha.rds does not contain eid.", call. = FALSE)
  miss_comp <- setdiff(ukb_component_cols, names(comp))
  if (length(miss_comp)) stop("UKB diet.maha.rds missing MAHA component columns: ", paste(miss_comp, collapse = ", "), call. = FALSE)
  comp <- comp %>% dplyr::select(eid, all_of(ukb_component_cols))
  dat <- dat %>% dplyr::select(-any_of(ukb_component_cols)) %>% left_join(comp, by = "eid")
}

ukb_validation <- ukb_validate_maha(dat)
ukb_construct <- ukb_construct_profile(dat)
})

maha_run_step("table1", {
step_header("Table 1")
score.cols <- c("bmi.pts","bp.pts","nonhdl.pts","smoke.pts","pa.pts","sleep.pts",diet.s100.cols)
lbls <- list(
	age ~ "Age", female ~ "Female", tdi ~ "TDI", bmi ~ "BMI", energy_kcal ~ "Total Energy Intake, kcal/d", bmi.pts ~ "BMI score", bp.pts ~ "Blood pressure score", nonhdl.pts ~ "Non-HDL cholesterol score", smoke.pts ~ "Smoking score", pa.pts ~ "Physical activity score", sleep.pts ~ "Sleep score"
)
for (n in names(diet.lst)) lbls[[paste0("diet.", n, ".s100")]] <- diet.lst[[n]]
reset_gtsummary_theme(); theme_gtsummary_compact()
table1 <- dat %>% mutate(
	MAHA_Group = factor(.data[[maha_3c]], c("low","middle","high"), c("Low Adherence","Moderate Adherence","High Adherence")),
	energy_kcal = energy.kJ * 0.239, female = if_else(tolower(as.character(sex.f)) == "female", 1L, 0L, missing = NA_integer_)
	) %>% dplyr::select(MAHA_Group, age, female, tdi, bmi, energy_kcal, all_of(score.cols)) %>% drop_na(MAHA_Group) %>%
	tbl_summary(by = MAHA_Group, type = list(female ~ "dichotomous", all_of(c("age","tdi","bmi","energy_kcal",score.cols)) ~ "continuous"), value = list(female ~ 1), statistic = list(all_continuous() ~ "{mean} ({sd})", female ~ "{n} ({p}%)"), digits = list(all_continuous() ~ 1, female ~ c(0,1)), label = lbls, missing = "no") %>%
	add_overall() %>% modify_header(label = "**Baseline Characteristics**") %>% modify_footnote(all_stat_cols() ~ "Continuous variables are presented as mean (SD)") %>% bold_labels()
table1$table_styling$footnote <- NULL
table1 %>% as_flex_table() %>% flextable::set_table_properties(layout = "autofit") %>% flextable::save_as_docx(path = cohort_file("Table1.docx"))
write_xlsx(table1$table_body %>% as.data.frame(), "Table1.out.xlsx")
})

maha_run_step("fig1", {
step_header("Figure 1: PheWAS")
pacman::p_load(PheWAS, writexl, ggrepel)

dat1 <- dat
Xs <- diet.inc.s100
phewas.rds <- cohort_file("phewas.res.rds")
phewas.res <- if (file.exists(phewas.rds)) {
	readRDS(phewas.rds)
} else {
	x <- plot_phewas(dat1, phecode = NA, Xs = Xs, varX = vars.basic)
	saveRDS(x, phewas.rds)
	x
}
phewas.res$plots <- Map(\(p, x)
	p + labs(title = var2lab(x)) + theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.title = element_text(face = "bold"),
		axis.title.x = element_text(size = 11, margin = margin(t = 6)),
		axis.text.x = element_text(size = 8, angle = -30, vjust = 1, hjust = 1, face = "bold")
	), phewas.res$plots, Xs
)
names(phewas.res$plots) <- Xs; phewas.res$plots

lab <- setNames(unname(diet.lst[gsub("^diet\\.|\\.(sum|pts|q5|s100|s100_q5)$", "", Xs)]), Xs)

res <- phewas.res$res %>%
	as_tibble() %>%
	mutate(
		Diet = recode(snp, !!!lab),
		category = tolower(category),
		description = stringr::str_squish(description),
		logp = -log10(p),
		FDR = p.adjust(p, "BH"),
		sig_bonf = bonferroni %in% TRUE
	) %>%
	filter(!is.na(Diet), is.finite(p), p > 0)

write_xlsx(res %>% filter(logp >= 10), "Fig1.phewas.top.xlsx")

sum_score <- res %>% group_by(Diet) %>% summarise(
	n_bonf = sum(sig_bonf), n_fdr = sum(FDR < 0.05), mean_logp = mean(logp), median_logp = median(logp),
	n_neg_bonf = sum(sig_bonf & beta < 0), n_pos_bonf = sum(sig_bonf & beta > 0), .groups = "drop"
) %>% arrange(desc(n_bonf))

sum_cat <- res %>% filter(Diet %in% c("DASH", "MAHA")) %>%
	group_by(Diet, category) %>%
	summarise(n_bonf = sum(sig_bonf), prop_bonf = mean(sig_bonf), mean_logp = mean(logp), .groups = "drop")

dm <- res %>% filter(Diet %in% c("DASH", "MAHA")) %>%
	dplyr::select(Diet, phenotype, description, category, beta, p, logp, sig_bonf) %>%
	pivot_wider(names_from = Diet, values_from = c(beta, p, logp, sig_bonf), names_sep = ".") %>%
	mutate(
		same_direction = sign(beta.DASH) == sign(beta.MAHA),
		delta_logp = logp.DASH - logp.MAHA,
		delta_beta = beta.DASH - beta.MAHA,
		sig_pattern = case_when(
			sig_bonf.DASH & sig_bonf.MAHA ~ "Both",
			sig_bonf.DASH & !sig_bonf.MAHA ~ "DASH_only",
			!sig_bonf.DASH & sig_bonf.MAHA ~ "MAHA_only", TRUE ~ "Neither"
		)
	)

dash_only <- dm %>% filter(sig_pattern == "DASH_only") %>% arrange(desc(logp.DASH))
maha_only <- dm %>% filter(sig_pattern == "MAHA_only") %>% arrange(desc(logp.MAHA))

cat_diff <- dm %>%
	group_by(category) %>%
	summarise(
		n_both = sum(sig_pattern == "Both"),
		n_dash_only = sum(sig_pattern == "DASH_only"),
		n_maha_only = sum(sig_pattern == "MAHA_only"),
		mean_delta_logp = mean(delta_logp),
		prop_same_direction = mean(same_direction),
		.groups = "drop"
	) %>%
	arrange(desc(mean_delta_logp))

death_tbl <- res %>% filter(description == dx.lst["death"]) %>% arrange(p) %>% dplyr::select(Diet, beta, SE, p, logp, sig_bonf, n_total)

phewas_xmap <- phewas.res$plots[[Xs[1]]]$data %>%
	transmute(category = tolower(category), phenotype = as.character(phenotype), phe_x = seq) %>%
	distinct(category, phenotype, .keep_all = TRUE)
phewas_centers <- phewas_xmap %>%
	group_by(category) %>%
	summarise(phe_x = mean(phe_x), first_x = min(phe_x), .groups = "drop") %>%
	arrange(first_x)
phewas_cols <- setNames(scales::hue_pal()(n_distinct(phewas_xmap$category)), sort(unique(phewas_xmap$category)))
bonf_line <- -log10(0.05 / n_distinct(phewas_xmap$phenotype))

	cat_span <- phewas_xmap %>%
		group_by(category) %>%
		summarise(first_x = min(phe_x), last_x = max(phe_x), phe_x = mean(phe_x), .groups = "drop") %>%
		mutate(cat_w = pmax(5, (last_x - first_x) * 0.18))
	bar_w0 <- (cat_span %>% filter(category == "neoplasms") %>% pull(cat_w))[1] * 0.42
	if (!is.finite(bar_w0)) bar_w0 <- median(cat_span$cat_w, na.rm = TRUE) * 0.42
	cat_span <- cat_span %>% mutate(bar_w = bar_w0, offset = bar_w0 * 0.72)
	cat_ymax <- max(sum_cat$n_bonf, na.rm = TRUE) * 1.14
	cat_labs <- cat_span %>% mutate(label_y = cat_ymax, text_y = cat_ymax * 0.972, label = paste0(" ", stringr::str_trim(category)))
	cat_bar <- sum_cat %>%
		left_join(cat_span, by = "category") %>%
		mutate(x = phe_x + if_else(Diet == "DASH", -offset, offset))
	cat_missing <- cat_span %>%
		left_join(sum_cat %>% group_by(category) %>% summarise(any_bar = any(n_bonf > 0), .groups = "drop"), by = "category") %>%
		filter(!any_bar %in% TRUE)
	dash_lines <- cat_bar %>%
		filter(Diet == "DASH", n_bonf > 2) %>%
		group_by(category, x, bar_w, n_bonf) %>%
		summarise(y = list(seq(2, n_bonf - 0.5, by = 4)), .groups = "drop") %>%
		unnest(y) %>%
		mutate(x0 = x - bar_w * 0.38, x1 = x + bar_w * 0.38)
	p_cat <- ggplot(cat_bar) +
		geom_rect(aes(xmin = x - bar_w / 2, xmax = x + bar_w / 2, ymin = 0, ymax = n_bonf, fill = category), color = "white", linewidth = 0.15, show.legend = FALSE) +
		geom_segment(data = dash_lines, aes(x = x0, xend = x1, y = y, yend = y), inherit.aes = FALSE, color = "white", linewidth = 0.25) +
		geom_point(data = cat_labs, aes(x = phe_x, y = label_y), inherit.aes = FALSE,
			shape = 25, size = 2.7, stroke = .25, fill = "black", color = "black") +
		geom_text(data = cat_labs, aes(x = phe_x, y = text_y, label = label), inherit.aes = FALSE, angle = -90, hjust = 0, vjust = 0.5, size = 2.35, fontface = "bold", color = "black", lineheight = 0.8) +
		geom_text(data = cat_missing, aes(x = phe_x, y = cat_ymax * 0.045, label = "x", color = category), inherit.aes = FALSE, size = 3.0, fontface = "bold", vjust = 0.5, show.legend = FALSE) +
		scale_fill_manual(values = phewas_cols, guide = "none") +
		scale_color_manual(values = phewas_cols, guide = "none") +
		scale_x_continuous(limits = range(phewas_xmap$phe_x), expand = expansion(mult = c(0.01, 0.01))) +
		scale_y_continuous(limits = c(0, cat_ymax), expand = expansion(mult = c(0, 0))) +
		coord_cartesian(clip = "on") +
		labs(title = "c. Significant associations", x = NULL, y = NULL, fill = NULL) +
		theme_bw(base_size = 9) +
		theme(
			plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
			legend.position = "none",
			axis.title = element_text(face = "bold"),
			axis.title.y = element_text(size = 8, margin = margin(r = 4)),
			axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		axis.text.y = element_text(face = "bold"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
	plot.margin = margin(8, 4, 8, 4)
	)

lab_dat <- dm %>% filter(pmin(p.DASH, p.MAHA) <= 1e-10) %>% slice_max(abs(delta_logp), n = 12) %>% mutate(dx = ifelse(logp.DASH >= logp.MAHA, 0.7, -0.7), dy = ifelse(logp.DASH >= logp.MAHA, -0.5, 0.5), hjust = ifelse(dx > 0, 0, 1))
p_cmp_p <- ggplot(dm, aes(logp.DASH, logp.MAHA)) +
	geom_point(alpha = 0.45) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
	geom_text(data = lab_dat, aes(x = logp.DASH + dx, y = logp.MAHA + dy, label = description, hjust = hjust), size = 3, fontface = "bold", check_overlap = TRUE) +
	labs(x = "DASH: -log10(P)", y = "MAHA: -log10(P)") +
	theme_bw(base_size = 14) +
	theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"))

lab_dat_b <- dm %>% filter(pmin(p.DASH, p.MAHA) <= 1e-10) %>% slice_max(abs(delta_beta), n = 12) %>% mutate(dx = ifelse(beta.DASH >= beta.MAHA, 0.0006, -0.0006), dy = ifelse(beta.DASH >= beta.MAHA, -0.0006, 0.0006), hjust = ifelse(dx > 0, 0, 1))
p_cmp_beta <- ggplot(dm, aes(beta.DASH, beta.MAHA)) +
	geom_point(alpha = 0.45) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
	geom_text(data = lab_dat_b, aes(x = beta.DASH + dx, y = beta.MAHA + dy, label = description, hjust = hjust), size = 3, fontface = "bold", check_overlap = TRUE) +
	labs(x = "DASH: beta", y = "MAHA: beta") +
	theme_bw(base_size = 14) +
	theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"))

phewas_maha <- make_phewas_panel("MAHA", show_x = FALSE)
phewas_dash <- make_phewas_panel("DASH", show_x = FALSE)
Fig1 <- phewas_maha / phewas_dash / p_cat
Fig1 <- Fig1 + plot_layout(heights = c(0.8, 0.8, 0.70))
save_plot(Fig1, "Fig1.phewas.png", width = 8.6, height = 8.2, dpi = 320)
write_xlsx(
	list(
		phewas_full = res,
		score_summary = sum_score,
		category_summary = sum_cat,
		dash_only = dash_only,
		maha_only = maha_only,
		category_difference = cat_diff,
		death = death_tbl
	),
	"Fig1.phewas.out.xlsx"
)
print(sum_score)
print(sum_cat, n = 50)
print(head(dash_only, 50))
print(maha_only, n = 20)
print(cat_diff, n = 30)
print(death_tbl)
})

maha_run_step("fig2", {
step_header("Figure 2: main and alternative MAHA specifications")

maha.extra <- setdiff(names(maha_map), "maha")
maha.extra <- maha.extra[paste0("diet.", maha.extra, ".sum") %in% names(dat)]
for (nm in maha.extra) {
	sumv <- paste0("diet.", nm, ".sum"); ptsv <- paste0("diet.", nm, ".pts")
	s100v <- paste0("diet.", nm, ".s100"); q5v <- paste0("diet.", nm, ".q5")
	s100q5v <- paste0("diet.", nm, ".s100_q5"); v3c <- paste0("diet.", nm, ".3c"); hmlv <- paste0("diet.", nm, ".hml")
	if (!ptsv %in% names(dat)) dat[[ptsv]] <- std_10_90(dat[[sumv]])
	if (!s100v %in% names(dat)) dat[[s100v]] <- score_0_100(dat[[sumv]])
	if (!q5v %in% names(dat)) dat[[q5v]] <- score_q5(dat[[sumv]])
	if (!s100q5v %in% names(dat)) dat[[s100q5v]] <- score_0_100_q5(dat[[sumv]])
	if (!v3c %in% names(dat)) dat[[v3c]] <- f3c(dat[[sumv]])
	if (!hmlv %in% names(dat)) dat[[hmlv]] <- mk_hml(dat[[ptsv]], group_pct_th)
}

outcome_levels <- rev(unname(dx.lst[Y.inc]))

diet_cols.ab <- c(MAHA = "#F26D60", DASH = "#16B9C0", MIND = "#B97AF7", MEDI = "#7CAE00")

fig3_outcome_cols <- c(
	"Coronary artery disease" = "#B97AF7",
	"Ischemic stroke" = "#F26D60",
	"Heart failure" = "#E68613",
	"Type 2 diabetes" = "#7CAE00",
	"Chronic kidney disease" = "#16B9C0",
	"All-cause mortality" = "#2C7FB8"
)
diet_levels.ab <- c("MAHA", "DASH", "MIND", "MEDI")
dat.ab <- dat
for (x in paste0("diet.", c("maha", "dash", "mind", "medi24"), ".pts")) dat.ab[[x]] <- as.numeric(scale(dat.ab[[x]]))

assoc.ab <- assoc_reg(
	dat.ab, paste0("diet.", c("maha", "dash", "mind", "medi24"), ".pts"),
	covs, names(dx.lst), type = "t2e"
) %>%
	mutate(
		Outcome = factor(unname(dx.lst[Outcome]), levels = outcome_levels),
		Diet = factor(
			c(maha = "MAHA", dash = "DASH", mind = "MIND", medi24 = "MEDI")[gsub("^diet\\.|\\.(pts|q5|qt)$", "", Exposure)],
			levels = diet_levels.ab
		)
	) %>%
	filter(Outcome %in% outcome_levels, Diet %in% diet_levels.ab)

diet.pairs <- c(dash = "DASH", medi24 = "MEDI", mind = "MIND")
disc.ab <- purrr::imap_dfr(diet.pairs, ~ purrr::map_dfr(Y.inc, fit_disc_reftrad, dat = dat, trad_nm = .y, trad_lab = .x, covs = covs, maha_nm = "maha", maha_lab = "MAHA")) %>%
	mutate(Pattern = factor(Pattern, levels = c("DASH", "MEDI", "MIND")))
aic.ab <- purrr::imap_dfr(diet.pairs, ~ purrr::map_dfr(Y.inc, fit_aic_reftrad, dat = dat, trad_nm = .y, trad_lab = .x, covs = covs, maha_nm = "maha", maha_lab = "MAHA"))

shift.ab <- c(MAHA = -0.27, DASH = -0.09, MIND = 0.09, MEDI = 0.27)
pA.dat <- assoc.ab %>% arrange(Outcome, Diet) %>% mutate(y = as.numeric(Outcome) + unname(shift.ab[as.character(Diet)]), lab = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high))
xA <- range(c(1, pA.dat$conf.low, pA.dat$conf.high), na.rm = TRUE)
forest_x_half <- max(abs(c(pA.dat$conf.low, pA.dat$conf.high) - 1), na.rm = TRUE)
forest_x_half <- max(0.18, forest_x_half * 1.10)
forest_xlim <- c(1 - forest_x_half, 1 + forest_x_half)
forest_x_txt <- 1.10

pA <- ggplot(pA.dat, aes(color = Diet)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
	geom_segment(aes(x = conf.low, xend = conf.high, y = y, yend = y), linewidth = 1) +
	geom_point(aes(x = estimate, y = y), size = 2.55) +
	geom_text(aes(x = forest_x_txt, y = y, label = lab, color = Diet), hjust = 0, size = 3.05, fontface = "bold", show.legend = FALSE) +
	scale_color_manual(values = diet_cols.ab, breaks = diet_levels.ab) +
	scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
	coord_cartesian(xlim = forest_xlim, clip = "off") +
	labs(title = "a. Main associations", x = "Hazard ratio", y = NULL) +
	theme_classic(base_size = 16) +
	theme(plot.title = element_text(face = "bold", hjust = 0, size = 16), legend.position = "none", axis.text.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"), plot.margin = margin(5.5, 34, 5.5, 42))

shift.b <- c(DASH = -0.22, MEDI = 0, MIND = 0.22)
pB.dat <- disc.ab %>% arrange(Outcome, Pattern) %>% mutate(y = as.numeric(Outcome) + unname(shift.b[as.character(Pattern)]))
xB <- range(c(1, pB.dat$conf.low, pB.dat$conf.high), na.rm = TRUE)

pB <- ggplot(pB.dat, aes(color = Pattern)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
	geom_segment(aes(x = conf.low, xend = conf.high, y = y, yend = y), linewidth = 1) +
	geom_point(aes(x = estimate, y = y), size = 3.2) +
	scale_color_manual(values = diet_cols.ab, breaks = diet_levels.ab) +
	scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
	coord_cartesian(xlim = c(xB[1] - 0.03, xB[2] + 0.05), clip = "off") +
	labs(title = "b. Discordant contrasts", x = "Hazard ratio", y = NULL) +
	theme_classic(base_size = 16) +
	theme(plot.title = element_text(face = "bold", hjust = 0, size = 16), legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(face = "bold"))

print(assoc.ab); print(disc.ab); print(aic.ab)

assoc.cd <- NULL
alt_lab <- c(maha_bal = "balanced", maha_strict = "strict", maha_nodairy = "no dairy", maha_noprotein = "no protein")
alt_lab <- alt_lab[maha.extra]
alt_levels <- unname(alt_lab)
alt_cols <- c("balanced" = "#E76F51", "strict" = "#C77DFF", "no dairy" = "#2A9D8F", "no protein" = "#577590")
dat.cd <- dat
for (x in paste0("diet.", maha.extra, ".pts")) dat.cd[[x]] <- as.numeric(scale(dat.cd[[x]]))

assoc.cd <- assoc_reg(dat.cd, paste0("diet.", maha.extra, ".pts"), covs, names(dx.lst), type = "t2e") %>%
	mutate(
		Outcome = factor(unname(dx.lst[Outcome]), levels = outcome_levels),
		Diet = factor(unname(alt_lab[gsub("^diet\\.|\\.(pts|q5|qt)$", "", Exposure)]), levels = alt_levels)
	) %>%
	filter(Outcome %in% outcome_levels, Diet %in% alt_levels)

shift.cd <- stats::setNames(seq(-0.27, 0.27, length.out = length(alt_levels)), alt_levels)
pS5B.dat <- assoc.cd %>% arrange(Outcome, Diet) %>% mutate(y = as.numeric(Outcome) + unname(shift.cd[as.character(Diet)]), lab = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high))
xS5B <- range(c(1, pS5B.dat$conf.low, pS5B.dat$conf.high), na.rm = TRUE)
forest_x_half <- max(abs(c(pA.dat$conf.low, pA.dat$conf.high, pS5B.dat$conf.low, pS5B.dat$conf.high) - 1), na.rm = TRUE)
forest_x_half <- max(0.18, forest_x_half * 1.10)
forest_xlim <- c(1 - forest_x_half, 1 + forest_x_half)
forest_x_txt <- 1.10

pS5A <- pA + coord_cartesian(xlim = forest_xlim, clip = "off") + labs(title = "a. Main associations") + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 4)), plot.margin = margin(5.5, 34, 0, 42))
pS5B <- ggplot(pS5B.dat, aes(color = Diet)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.8) +
	geom_segment(aes(x = conf.low, xend = conf.high, y = y, yend = y), linewidth = 1) +
	geom_point(aes(x = estimate, y = y), shape = 21, fill = "white", stroke = 1.0, size = 2.55) +
	geom_text(aes(x = forest_x_txt, y = y, label = lab, color = Diet), hjust = 0, size = 3.05, fontface = "bold", show.legend = FALSE) +
	scale_color_manual(values = alt_cols[alt_levels], breaks = alt_levels) +
	scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
	coord_cartesian(xlim = forest_xlim, clip = "off") +
	labs(title = "d. Alternative MAHA specifications", x = "Hazard ratio", y = NULL) +
	theme_classic(base_size = 16) +
	theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 4)), legend.position = "none", axis.text.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"), plot.margin = margin(5.5, 34, 0, 42))

g_leg.s5.main <- ggplot(data.frame(Diet = factor(diet_levels.ab, levels = diet_levels.ab)), aes(1, 1, color = Diet)) +
	geom_point(size = 4) +
	scale_color_manual(values = diet_cols.ab, breaks = diet_levels.ab) +
	guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 5))) +
	theme_void() +
	theme(
		legend.position = "bottom",
		legend.justification = "center",
		legend.title = element_blank(),
		legend.text = element_text(size = 12, face = "bold", margin = margin(l = 4, r = 14)),
		legend.key.width = grid::unit(0.50, "cm"),
		legend.spacing.x = grid::unit(0.10, "cm"),
		legend.margin = margin(0, 0, 0, 0),
		legend.box.margin = margin(0, 0, 0, 0)
	)
leg.s5.main0 <- cowplot::get_legend(g_leg.s5.main)
leg.s5.main <- cowplot::plot_grid(
	ggplot() + theme_void(), leg.s5.main0,
	nrow = 1, rel_widths = c(0.24, 0.76)
)

g_leg.s5.alt <- ggplot(data.frame(Diet = factor(alt_levels, levels = alt_levels)), aes(1, 1, color = Diet)) +
	geom_point(size = 4) +
	scale_color_manual(values = alt_cols[alt_levels], breaks = alt_levels) +
	guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(shape = 21, fill = "white", stroke = 1.1, size = 5))) +
	theme_void() +
	theme(
		legend.position = "bottom",
		legend.justification = "center",
		legend.title = element_blank(),
		legend.text = element_text(size = 12, face = "bold", margin = margin(l = 5, r = 16)),
		legend.key.width = grid::unit(0.56, "cm"),
		legend.spacing.x = grid::unit(0.12, "cm"),
		legend.margin = margin(0, 0, 0, 0),
		legend.box.margin = margin(0, 0, 0, 0)
	)
leg.s5.alt0 <- cowplot::get_legend(g_leg.s5.alt)
leg.s5.alt <- cowplot::plot_grid(
	ggplot() + theme_void(), leg.s5.alt0,
	nrow = 1, rel_widths = c(0.24, 0.76)
)

alt_note <- tibble(
	note = "Panel A copies Figure 2A main associations. Panel B uses the same association style for alternative MAHA specifications; hollow circles distinguish alternative MAHA estimates from the solid-circle main dietary patterns."
)

fig2_out <- list(
	note = alt_note,
	Fig2_main_assoc = assoc.ab %>% transmute(Outcome = as.character(Outcome), Diet = as.character(Diet), HR = estimate, LCI = conf.low, UCI = conf.high, P = p.value, N_total, N_event),
	Fig2_main_discordant = disc.ab %>% transmute(Outcome = as.character(Outcome), MAHA = "MAHA", Pattern = as.character(Pattern), contrast, HR = estimate, LCI = conf.low, UCI = conf.high, P = p.value),
	Fig2_main_aic = aic.ab,
	Alternative_MAHA_assoc = assoc.cd %>% transmute(Outcome = as.character(Outcome), Diet = as.character(Diet), HR = estimate, LCI = conf.low, UCI = conf.high, P = p.value, N_total, N_event)
)

print(assoc.cd)

step_header("Figure 2 right panels: 10-year all-cause mortality risk")
dat1 <- dat %>% mutate(diet.dash.3c = factor(diet.dash.3c, levels = c("low", "middle", "high")), maha3c = factor(.data[[maha_3c]], levels = c("low", "middle", "high")))

Y <- "death"

risk10_by_dash <- extract_joint_risk(dat1, paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), "diet.dash.3c", "maha3c", vars.basic, times = 10, Xfacetlab = "DASH", Xlinelab = "MAHA", panel = "b. 10-year risk by DASH strata")
risk10_by_maha <- extract_joint_risk(dat1, paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), "maha3c", "diet.dash.3c", vars.basic, times = 10, Xfacetlab = "MAHA", Xlinelab = "DASH", panel = "c. 10-year risk by MAHA strata")
incidence_by_dash <- extract_joint_risk(dat1, paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), "diet.dash.3c", "maha3c", vars.basic, times = seq(0, 10, by = 0.1), Xfacetlab = "DASH", Xlinelab = "MAHA", panel = "e. Incidence by DASH strata")
incidence_by_maha <- extract_joint_risk(dat1, paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), "maha3c", "diet.dash.3c", vars.basic, times = seq(0, 10, by = 0.1), Xfacetlab = "MAHA", Xlinelab = "DASH", panel = "f. Incidence by MAHA strata")

pa <- plot_risk(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), X1 = "diet.dash.3c", X2 = "maha3c", covs = vars.basic, method = "10years", t0 = 10, X1lab = "DASH", X2lab = "MAHA", group_bgcolor = TRUE, title = NULL, leg_position = "top", leg_direction = "horizontal", tab = TRUE)
pb <- plot_risk(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), X1 = "maha3c", X2 = "diet.dash.3c", covs = vars.basic, method = "10years", t0 = 10, X1lab = "MAHA", X2lab = "DASH", group_bgcolor = TRUE, title = NULL, leg_position = "top", leg_direction = "horizontal", tab = TRUE)
fig2_title_theme <- theme(
	plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 4)),
	plot.title.position = "plot"
)

pa <- pa + labs(title = "b. 10-year risk by DASH strata") + fig2_title_theme +
	theme(
		legend.justification = "center",
		legend.text = element_text(color = "grey30", face = "bold", size = 12, margin = margin(l = 4, r = 10)),
		legend.title = element_text(color = "grey30", face = "bold", size = 12, margin = margin(r = 10)),
		legend.key.width = grid::unit(0.44, "cm"),
		legend.spacing.x = grid::unit(0.14, "cm"),
		axis.title.x = element_text(face = "bold", size = 12),
		axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.margin = margin(5.5, 14, 5.5, 5.5)
	)
pb <- pb + labs(title = "c. 10-year risk by MAHA strata") + fig2_title_theme +
	theme(
		legend.justification = "center",
		legend.text = element_text(color = "grey30", face = "bold", size = 12, margin = margin(l = 4, r = 10)),
		legend.title = element_text(color = "grey30", face = "bold", size = 12, margin = margin(r = 10)),
		legend.key.width = grid::unit(0.44, "cm"),
		legend.spacing.x = grid::unit(0.14, "cm"),
		axis.title.x = element_text(face = "bold", size = 12),
		axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.margin = margin(5.5, 14, 5.5, 5.5)
	)

pc <- plot_cuminc_joint(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), Xfacet = "diet.dash.3c", Xline = "maha3c", covs = vars.basic, t0 = 10, Xfacetlab = "DASH", Xlinelab = "MAHA", title = NULL, leg_position = "top", leg_direction = "horizontal")
pd <- plot_cuminc_joint(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), Xfacet = "maha3c", Xline = "diet.dash.3c", covs = vars.basic, t0 = 10, Xfacetlab = "MAHA", Xlinelab = "DASH", title = NULL, leg_position = "top", leg_direction = "horizontal")
pc <- pc + labs(title = "e. Incidence by DASH strata") + fig2_title_theme +
	theme(
		legend.text = element_text(color = "grey30", face = "bold", size = 10, margin = margin(l = 4, r = 8)),
		legend.title = element_text(color = "grey30", face = "bold", size = 12, margin = margin(r = 10)),
		legend.key.width = grid::unit(0.44, "cm"),
		legend.spacing.x = grid::unit(0.12, "cm"),
		axis.title.x = element_text(face = "bold", size = 12),
		axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.margin = margin(5.5, 18, 5.5, 5.5)
	)
pd <- pd + labs(title = "f. Incidence by MAHA strata") + fig2_title_theme +
	theme(
		legend.text = element_text(color = "grey30", face = "bold", size = 10, margin = margin(l = 4, r = 8)),
		legend.title = element_text(color = "grey30", face = "bold", size = 12, margin = margin(r = 10)),
		legend.key.width = grid::unit(0.44, "cm"),
		legend.spacing.x = grid::unit(0.12, "cm"),
		axis.title.x = element_text(face = "bold", size = 12),
		axis.title.y = element_text(face = "bold", margin = margin(r = 4)),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		plot.margin = margin(5.5, 18, 5.5, 5.5)
	)

pS5A <- pS5A + fig2_title_theme
pS5B <- pS5B + fig2_title_theme

forest_top <- cowplot::plot_grid(
	pS5A, leg.s5.main,
	ncol = 1,
	rel_heights = c(1.08, 0.075),
	align = "none"
)
forest_bottom <- cowplot::plot_grid(
	pS5B, leg.s5.alt,
	ncol = 1,
	rel_heights = c(1.08, 0.075),
	align = "none"
)
fig2_row_gap <- ggplot() + theme_void()

fig2_row1 <- cowplot::plot_grid(
	forest_top, pa, pb,
	nrow = 1,
	rel_widths = c(1.22, 0.78, 0.78),
	align = "none"
)
fig2_row2 <- cowplot::plot_grid(
	forest_bottom, pc, pd,
	nrow = 1,
	rel_widths = c(1.22, 0.78, 0.78),
	align = "none"
)

Fig2 <- cowplot::plot_grid(
	fig2_row1, fig2_row_gap, fig2_row2,
	ncol = 1,
	rel_heights = c(1.00, 0.045, 1.00),
	align = "none"
)
save_plot(Fig2, "Fig2.prospective.png", width = 16.0, height = 11.0, dpi = 320)
fig2_tab <- dat1 %>% count(diet.dash.3c, maha3c, name = "N") %>% rename(diet.maha.3c = maha3c) %>% mutate(Outcome = dx.lst[Y])
fig2_out <- c(fig2_out, list(
	Fig2_10y_risk_by_DASH = risk10_by_dash,
	Fig2_10y_risk_by_MAHA = risk10_by_maha,
	Fig2_incidence_by_DASH = incidence_by_dash,
	Fig2_incidence_by_MAHA = incidence_by_maha,
	Fig2_strata_counts = fig2_tab
))
write_xlsx(fig2_out, "Fig2.prospective.out.xlsx")
print(risk10_by_dash)
print(risk10_by_maha)
print(fig2_tab)
})

maha_run_step("fig3", {
ap_ukb_equiv_main <- ap_effect_from_assoc(
	assoc.ab %>% filter(as.character(Diet) %in% c("MAHA", "DASH", "MIND", "MEDI")),
	"UK Biobank main outcomes"
)

ap_ukb_equiv_alt <- if (exists("assoc.cd")) {
	ap_effect_from_assoc(assoc.cd, "UK Biobank alternative MAHA")
} else {
	tibble()
}

ap_ukb_head2head <- purrr::map_dfr(Y.inc, function(y) {
	purrr::map2_dfr(c("dash", "mind", "medi24"), c("DASH", "MIND", "MEDI"), function(nm, lab) {
		ap_fit_head2head_ukb(dat, y, nm, lab, covs, maha_nm = "maha", maha_lab = "MAHA")
	})
})

ap_margin_main <- 1.05
ap_prior_main  <- 1.10

if (!exists("diet_levels.ab")) diet_levels.ab <- c("MAHA", "DASH", "MIND", "MEDI")
if (!exists("outcome_levels")) outcome_levels <- rev(unname(dx.lst[Y.inc]))
if (!exists("fig3_outcome_cols")) {
	fig3_outcome_cols <- c(
		"Coronary artery disease" = "#B97AF7",
		"Ischemic stroke" = "#F26D60",
		"Heart failure" = "#E68613",
		"Type 2 diabetes" = "#7CAE00",
		"Chronic kidney disease" = "#16B9C0",
		"All-cause mortality" = "#2C7FB8"
	)
}

ap_shape_guide <- guides(
	shape = guide_legend(
		nrow = 1,
		byrow = TRUE,
		keywidth = grid::unit(0.56, "cm"),
		override.aes = list(size = 3.4, stroke = 1.2)
	)
)

ap_leg_theme <- theme(
	legend.position = "bottom",
	legend.spacing.x = grid::unit(0.95, "cm"),
	legend.key.width = grid::unit(0.52, "cm"),
	legend.text = element_text(face = "bold", margin = margin(l = 4, r = 26))
)

ap_fig3_ukb <- ap_make_fig3_outputs(
	ap_ukb_equiv_main,
	ap_ukb_head2head,
	"UK Biobank",
	"Fig3.precision_attenuation",
	equiv_alt_tbl = if (exists("ap_ukb_equiv_alt")) ap_ukb_equiv_alt else NULL
)
})

maha_run_step_cached("fig5", c(
	cohort_file("Fig5.geography.png"),
	cohort_file("Fig5.geography.out.xlsx")
), {
step_header("Figure 5: geography vs TDI clustering")
pacman::p_load(dplyr, tidyr, ggplot2, patchwork, sf, rnaturalearth, ggspatial, rosm, raster, scales)

fig5_tile_cache <- file.path(dir0, "data", "ukb", "map", "maha", "rosm.cache")
fig5_tile_files <- list.files(file.path(fig5_tile_cache, "osm"), pattern = "[.]png$", full.names = TRUE)
if (length(fig5_tile_files) < 90L) {
	stop("Fig5 persistent OSM cache is incomplete: expected at least 90 PNG tiles in ",
		file.path(fig5_tile_cache, "osm"), "; found ", length(fig5_tile_files), call. = FALSE)
}

if (TRUE) {
center.lst.inc <- c("c11020" = "Croydon", "c11018" = "Hounslow", "c11011" = "Bristol", "c11005" = "Edinburgh")
score_map <- c(dash = "DASH", medi24 = "MEDI", mind = "MIND", maha = "MAHA")
score_cols <- paste0("diet.", names(score_map), ".s100")
fs_title <- 11; fs_text <- 9.5; fs_small <- 9; tile_size <- 2000
blue_axis <- "#2171B5"

miss <- setdiff(c("home_east", "home_north", "tdi", "center", score_cols), names(dat))
if (length(miss)) stop("Missing required columns for Fig5: ", paste(miss, collapse = ", "))

fig5_wide <- dat %>%
	filter(as.character(center) %in% names(center.lst.inc)) %>%
	transmute(
		center_code = as.character(center),
		center_name = factor(center.lst.inc[as.character(center)], levels = unname(center.lst.inc)),
		home_east = as.numeric(home_east), home_north = as.numeric(home_north), tdi = as.numeric(tdi),
		MAHA = as.numeric(diet.maha.s100), DASH = as.numeric(diet.dash.s100),
		MIND = as.numeric(diet.mind.s100), MEDI = as.numeric(diet.medi24.s100)
	) %>%
	filter(if_all(c(home_east, home_north, tdi, MAHA, DASH, MIND, MEDI), is.finite)) %>%
	group_by(center_name) %>%
	mutate(dist = sqrt((home_east - median(home_east))^2 + (home_north - median(home_north))^2)) %>%
	filter(dist <= quantile(dist, 0.98, na.rm = TRUE)) %>%
	ungroup() %>% dplyr::select(-dist) %>%
	mutate(spatial_tile = interaction(center_name, floor(home_east / tile_size), floor(home_north / tile_size), drop = TRUE), tdi_z = z1(tdi))

valid_tiles <- fig5_wide %>% count(spatial_tile) %>% filter(n >= 30) %>% pull(spatial_tile)
fig5_long <- fig5_wide %>% pivot_longer(all_of(unname(score_map)), names_to = "Diet", values_to = "score") %>%
	mutate(Diet = factor(Diet, levels = unname(score_map)))
fig5_model_long <- fig5_wide %>% filter(spatial_tile %in% valid_tiles) %>%
	mutate(spatial_tile = droplevels(spatial_tile)) %>%
	pivot_longer(all_of(unname(score_map)), names_to = "Diet", values_to = "score") %>%
	mutate(Diet = factor(Diet, levels = unname(score_map)))
variance_decomp <- bind_rows(lapply(levels(fig5_model_long$Diet), function(x) decomp_one(fig5_model_long, x)))
variance_plot_dat <- variance_decomp %>%
	dplyr::select(Diet, TDI = R2_TDI_only, Centers = R2_center_only, `Spatial tiles` = R2_grid_only,
		`TDI | tiles` = unique_TDI_given_grid, `Tiles | TDI` = unique_grid_given_TDI) %>%
	pivot_longer(-Diet, names_to = "Component", values_to = "R2") %>% mutate(R2_pct = pmax(0, 100*R2))
cluster_ratio <- variance_decomp %>% transmute(Diet,
	Geography_to_TDI_ratio = R2_grid_only / pmax(R2_TDI_only, 1e-6),
	TDI_R2_pct = 100*R2_TDI_only, Spatial_grid_R2_pct = 100*R2_grid_only)

fig5_map_dat <- fig5_wide %>% mutate(center = as.character(center_name), MAHA_DASH_discordance = abs(MAHA - DASH))
p_map_body <- plot_map(fig5_map_dat,
	centers = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), separate = TRUE, width = 50, nrow = 2,
	transparency = "MAHA_DASH_discordance", center.lst = center.lst,
	basemap = TRUE, basemap_type = "osm", tile_zoom = 11, tile_cachedir = fig5_tile_cache)
p_map_body[[2]] <- p_map_body[[2]] + scale_fill_manual(values = c(Hounslow = "#F28E2B"), guide = "none")
p_map <- wrap_plots(
	ggplot() + annotate("text", .5, .5, label = "", fontface = "bold", size = fs_title / 3) +
		xlim(0, 1) + ylim(0, 1) + theme_void(),
	wrap_elements(full = p_map_body), ncol = 1, heights = c(.05, 1)
) & theme(plot.margin = margin(2, 0, 2, 4))

geo_plot <- variance_decomp %>% transmute(Diet=factor(Diet, levels=rev(unname(score_map))),
	TDI=100*R2_TDI_only, Tiles=100*R2_grid_only, TDI_after_tiles=100*pmax(0,unique_TDI_given_grid), Tiles_after_TDI=100*pmax(0,unique_grid_given_TDI))
p_decomp <- ggplot(geo_plot, aes(y=Diet)) +
	geom_segment(aes(x=TDI, xend=Tiles, yend=Diet), color="grey72", linewidth=1.2) +
	geom_point(aes(x=TDI, color="TDI"), size=3) + geom_point(aes(x=Tiles, color="Spatial tiles"), size=3) +
	geom_text(aes(x=TDI, label=sprintf("%.2f",TDI)), nudge_y=.21, size=fs_small/3.2, fontface="bold", color="#E76F51") +
	geom_text(aes(x=Tiles, label=sprintf("%.2f",Tiles)), nudge_y=-.21, size=fs_small/3.2, fontface="bold", color=blue_axis) +
	scale_color_manual(values=c(TDI="#E76F51", `Spatial tiles`=blue_axis)) +
	labs(title=" ", x="Variance explained (R2, %)", y=NULL, color=NULL) + theme_classic(base_size=fs_text) +
	theme(plot.title=element_text(face="bold",hjust=.5,size=fs_title), axis.text=element_text(face="bold",size=fs_text), axis.title=element_text(face="bold"), legend.position="top", legend.text=element_text(face="bold"))
p_ratio <- ggplot(geo_plot, aes(TDI_after_tiles, Tiles_after_TDI, color=Diet, label=Diet)) +
	geom_abline(slope=1,intercept=0,linetype=2,color="grey60") + geom_segment(aes(x=0,y=0,xend=TDI_after_tiles,yend=Tiles_after_TDI), linewidth=.55, alpha=.45) +
	geom_point(size=3.2) + ggrepel::geom_text_repel(fontface="bold",size=fs_small/3.1,show.legend=FALSE,box.padding=.25) +
	scale_color_manual(values=c(MAHA="#F8766D",DASH="#00BFC4",MIND="#B77CFF",MEDI="#7CAE00"),guide="none") +
	scale_x_continuous(labels=scales::label_number(accuracy=.01)) + scale_y_continuous(labels=scales::label_number(accuracy=.01)) +
	labs(title="c. Geographic contribution", x="TDI after spatial tiles (R2, %)", y="Spatial tiles after TDI (R2, %)") + theme_classic(base_size=fs_text) +
	theme(plot.title=element_text(face="bold",hjust=.5,size=fs_title), axis.text=element_text(face="bold",size=fs_text), axis.title=element_text(face="bold",size=fs_text))

right_col <- p_decomp / plot_spacer() / p_ratio + plot_layout(heights=c(1, .10, 1)) & theme(plot.margin=margin(2,4,2,0))
Fig5 <- cowplot::ggdraw() +
	cowplot::draw_plot(p_map, x=0, y=0, width=.69, height=1) +
	cowplot::draw_plot(right_col, x=.64, y=0, width=.34, height=1) +
	cowplot::draw_label("a. MAHA-DASH discordance", x=.35, y=.985, size=fs_title, fontface="bold") +
	cowplot::draw_label("b. Socioeconomic vs spatial variation", x=.831, y=.985, size=fs_title, fontface="bold")
fig5_print_with_retry <- function(plot, attempts = 8L) {
	for (attempt in seq_len(attempts)) {
		ok <- tryCatch({ print(plot); TRUE }, error = function(e) {
			message("Fig5 OSM render attempt ", attempt, "/", attempts,
				" failed; retrying missing tiles. Cause: ", conditionMessage(e))
			FALSE
		})
		if (ok) return(invisible(plot))
		Sys.sleep(min(2 * attempt, 10))
	}
	stop("Fig5 OSM street tiles still incomplete after ", attempts,
		" render attempts. No completion cache was recorded.", call. = FALSE)
}
fig5_print_with_retry(Fig5)
save_plot(Fig5, "Fig5.geography.png", width = 12.5, height = 8, dpi = 320, bg = "white")
save_pub_plot(Fig5, "Fig5.ukb.geography.png", width = 12.5, height = 8, dpi = 320, bg = "white")

center_score_summary <- fig5_wide %>% group_by(center_name) %>% summarise(N = n(), TDI_mean = mean(tdi), TDI_sd = sd(tdi),
	across(all_of(unname(score_map)), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"), MAHA_minus_DASH_mean = mean(MAHA - DASH), .groups = "drop")
tdi_cor_by_center <- fig5_long %>% group_by(center_name, Diet) %>% summarise(N = n(), r_spearman = cor1(score, tdi), P = cor_p(score, tdi), .groups = "drop")
tile_summary <- fig5_wide %>% group_by(center_name, spatial_tile) %>% summarise(N = n(), TDI_mean = mean(tdi),
	MAHA_mean = mean(MAHA), DASH_mean = mean(DASH), MIND_mean = mean(MIND), MEDI_mean = mean(MEDI), MAHA_minus_DASH = mean(MAHA - DASH), .groups = "drop") %>%
	arrange(center_name, desc(abs(MAHA_minus_DASH)))
fig5_out <- list(center_score_summary = center_score_summary, tdi_cor_by_center = tdi_cor_by_center,
	variance_decomposition = variance_decomp, variance_plot_data = variance_plot_dat,
	geography_to_tdi_ratio = cluster_ratio, tile_summary = tile_summary)
write_xlsx(fig5_out, "Fig5.geography.out.xlsx")
write_pub_xlsx(fig5_out, "Fig5.ukb.geography.out.xlsx")
print(center_score_summary); print(variance_decomp); print(cluster_ratio)
} else {

center.lst.inc <- c("c11020" = "Croydon", "c11018" = "Hounslow", "c11011" = "Bristol", "c11005" = "Edinburgh")
score_map <- c(MAHA = "diet.maha.s100", DASH = "diet.dash.s100", MIND = "diet.mind.s100", MEDI = "diet.medi24.s100")
tile_size <- 2000
min_tile_n <- 30

need_fig5 <- unique(c(
  "home_east", "home_north", "tdi", "center", "age", "sex.f", "death.t2e", "death.Yt2e",
  unname(score_map), setdiff(covs, "center")
))
miss_fig5 <- setdiff(need_fig5, names(dat))
if (length(miss_fig5)) stop("UKB Fig5 missing columns: ", paste(miss_fig5, collapse = ", "), call. = FALSE)

geo <- dat %>%
  filter(as.character(center) %in% names(center.lst.inc)) %>%
  transmute(
    eid = eid,
    center_code = as.character(center),
    center_name = factor(unname(center.lst.inc[as.character(center)]), levels = unname(center.lst.inc)),
    home_east = as.numeric(home_east), home_north = as.numeric(home_north),
    age = as.numeric(age), sex.f = sex.f, tdi = as.numeric(tdi),
    death.t2e = as.numeric(death.t2e), death.Yt2e = as.integer(death.Yt2e),
    across(all_of(setdiff(covs, c("age", "sex.f", "center", "tdi")))),
    MAHA = as.numeric(diet.maha.s100),
    DASH = as.numeric(diet.dash.s100),
    MIND = as.numeric(diet.mind.s100),
    MEDI = as.numeric(diet.medi24.s100)
  ) %>%
  filter(if_all(c(home_east, home_north, age, tdi, MAHA, DASH, MIND, MEDI), is.finite)) %>%
  group_by(center_name) %>%
  mutate(radial_distance = sqrt((home_east - median(home_east))^2 + (home_north - median(home_north))^2)) %>%
  filter(radial_distance <= quantile(radial_distance, 0.98, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(-radial_distance) %>%
  mutate(
    anchor = rowMeans(cbind(DASH, MIND, MEDI), na.rm = FALSE),
    discordance = MAHA - anchor,
    tile_east = (floor(home_east / tile_size) + 0.5) * tile_size,
    tile_north = (floor(home_north / tile_size) + 0.5) * tile_size,
    spatial_tile = interaction(center_name, floor(home_east / tile_size), floor(home_north / tile_size), drop = TRUE),
    age_z = zstd(age), tdi_z = zstd(tdi)
  )

geo_cc <- geo %>% drop_na(discordance, age_z, sex.f, tdi_z, center_name)
fit_disc_geo <- lm(discordance ~ age_z + sex.f + tdi_z + center_name, data = geo_cc)
geo_cc$discordance_resid <- residuals(fit_disc_geo)

tile_map <- geo_cc %>%
  group_by(center_name, spatial_tile, tile_east, tile_north) %>%
  summarise(
    N = n(), residual_discordance = mean(discordance_resid),
    raw_MAHA_minus_anchor = mean(discordance), TDI = mean(tdi), .groups = "drop"
  ) %>%
  filter(N >= min_tile_n)
if (nrow(tile_map) < 20) stop("UKB Fig5 has too few 2-km tiles after N>=30 privacy/stability filter.", call. = FALSE)

lim <- quantile(abs(tile_map$residual_discordance), 0.98, na.rm = TRUE, names = FALSE)
if (!is.finite(lim) || lim == 0) lim <- max(abs(tile_map$residual_discordance), na.rm = TRUE)

p_map <- ggplot(tile_map, aes(tile_east, tile_north, fill = residual_discordance)) +
  geom_tile(width = tile_size * 0.94, height = tile_size * 0.94) +
  facet_wrap(~center_name, scales = "free", ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), oob = scales::squish,
                       name = "MAHA − established\nanchor (residual points)") +
  labs(
    title = "a. Neighborhoods where MAHA and established diets disagree",
    subtitle = "2-km cells; N ≥ 30; residualized for age, sex, deprivation, and assessment centre",
    x = "British National Grid easting", y = "Northing"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        panel.grid = element_blank(), axis.text = element_text(size = 7), legend.title = element_text(face = "bold"))

partial_r2 <- function(y, d) {
  dd <- d %>% mutate(y = zstd(.data[[y]])) %>% drop_na(y, age_z, sex.f, center_name, tdi_z)
  m0 <- lm(y ~ age_z + sex.f + center_name, data = dd)
  m1 <- lm(y ~ age_z + sex.f + center_name + tdi_z, data = dd)
  sse0 <- sum(residuals(m0)^2); sse1 <- sum(residuals(m1)^2)
  pmax(0, (sse0 - sse1) / sse0)
}
spatial_icc <- function(y, d) {
  dd <- d %>% filter(spatial_tile %in% tile_map$spatial_tile) %>% mutate(y = zstd(.data[[y]])) %>%
    drop_na(y, age_z, sex.f, center_name, tdi_z, spatial_tile)
  if (nrow(dd) < 100) return(NA_real_)
  fit <- lm(y ~ age_z + sex.f + center_name + tdi_z, data = dd)
  rr <- residuals(fit)
  g <- droplevels(factor(dd$spatial_tile)); n <- length(rr); k <- nlevels(g)
  if (k < 2 || n <= k) return(NA_real_)
  gm <- tapply(rr, g, mean); gn <- table(g); grand <- mean(rr)
  ssb <- sum(as.numeric(gn) * (gm - grand)^2); ssw <- sum((rr - gm[g])^2)
  msb <- ssb / (k - 1); msw <- ssw / (n - k); n0 <- (n - sum(as.numeric(gn)^2) / n) / (k - 1)
  pmax(0, (msb - msw) / (msb + (n0 - 1) * msw))
}
structure_tbl <- bind_rows(lapply(names(score_map), function(lbl) tibble(
  Diet = lbl,
  TDI_partial_R2 = partial_r2(lbl, geo_cc),
  Spatial_ICC = spatial_icc(lbl, geo_cc)
))) %>% mutate(Diet = factor(Diet, levels = names(score_map)))
structure_long <- structure_tbl %>%
  pivot_longer(c(TDI_partial_R2, Spatial_ICC), names_to = "Metric", values_to = "value") %>%
  mutate(Metric = recode(Metric, TDI_partial_R2 = "Deprivation: partial R²", Spatial_ICC = "Neighborhood: residual ICC"))

p_structure <- ggplot(structure_long, aes(Diet, 100 * value, group = Metric, shape = Metric)) +
  geom_point(size = 3.2) + geom_line(aes(linetype = Metric), linewidth = 0.7) +
  labs(title = "b. Does each score track real population structure?", x = NULL, y = "Independent structure explained (%)", shape = NULL, linetype = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(face = "bold"), legend.position = "bottom")

valid_person_tiles <- tile_map$spatial_tile
mort_covs <- setdiff(covs, c("center", "age", "sex.f", "tdi"))
fit_within_tile <- function(lbl) {
  dd <- geo %>% filter(spatial_tile %in% valid_person_tiles) %>%
    mutate(score_z = zstd(.data[[lbl]]), spatial_tile = droplevels(factor(spatial_tile))) %>%
    dplyr::select(death.t2e, death.Yt2e, score_z, age, sex.f, tdi, all_of(mort_covs), spatial_tile) %>% drop_na()
  fm <- as.formula(paste0("Surv(death.t2e, death.Yt2e) ~ score_z + age + sex.f + tdi",
                         if (length(mort_covs)) paste0(" + ", paste(mort_covs, collapse = " + ")) else "",
                         " + strata(spatial_tile)"))
  fit <- coxph(fm, data = dd, ties = "efron")
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "score_z")
  tibble(Diet = lbl, N = nrow(dd), Deaths = sum(dd$death.Yt2e == 1), HR = tt$estimate, CI_low = tt$conf.low, CI_high = tt$conf.high, P = tt$p.value)
}
within_tile <- bind_rows(lapply(names(score_map), fit_within_tile)) %>% mutate(Diet = factor(Diet, levels = rev(names(score_map))))

p_mort <- ggplot(within_tile, aes(HR, Diet)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_segment(aes(x = CI_low, xend = CI_high, yend = Diet), linewidth = 0.7) +
  geom_point(size = 3) +
  scale_x_log10() +
  labs(title = "c. Mortality association within the same neighborhood", subtitle = "Cox model stratified by 2-km cell", x = "HR per 1-SD healthier diet score (log scale)", y = NULL) +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))

Fig5 <- p_map / (p_structure | p_mort) + plot_layout(heights = c(1.45, 1))
print(Fig5)
save_plot(Fig5, "Fig5.geography.png", width = 13.2, height = 10.2, dpi = 360, bg = "white")

geo_sensitivity_one <- function(tile_m, min_n) {
  dd <- geo %>% mutate(
    tile_x = floor(home_east / tile_m), tile_y = floor(home_north / tile_m),
    tile = interaction(center_name, tile_x, tile_y, drop = TRUE)
  )
  ok <- dd %>% count(tile) %>% filter(n >= min_n) %>% pull(tile)
  if (length(ok) < 5) return(tibble())
  bind_rows(lapply(names(score_map), function(lbl) {
    dx <- dd %>% filter(tile %in% ok) %>% mutate(score_z = zstd(.data[[lbl]]), age_z2 = zstd(age), tdi_z2 = zstd(tdi), tile = droplevels(factor(tile))) %>%
      drop_na(score_z, age_z2, sex.f, center_name, tdi_z2, tile, death.t2e, death.Yt2e)
    fit0 <- lm(score_z ~ age_z2 + sex.f + center_name + tdi_z2, data = dx)
    rr <- residuals(fit0); g <- dx$tile; n <- length(rr); k <- nlevels(g)
    icc <- NA_real_
    if (k >= 2 && n > k) {
      gm <- tapply(rr, g, mean); gn <- table(g); grand <- mean(rr)
      ssb <- sum(as.numeric(gn) * (gm - grand)^2); ssw <- sum((rr - gm[g])^2)
      msb <- ssb/(k-1); msw <- ssw/(n-k); n0 <- (n - sum(as.numeric(gn)^2)/n)/(k-1)
      icc <- pmax(0, (msb-msw)/(msb+(n0-1)*msw))
    }
    mc <- setdiff(covs, c("center", "age", "sex.f", "tdi")); dm <- dx %>% dplyr::select(death.t2e, death.Yt2e, score_z, age, sex.f, tdi, all_of(mc), tile) %>% drop_na()
    fm <- as.formula(paste0("Surv(death.t2e, death.Yt2e) ~ score_z + age + sex.f + tdi",
                           if(length(mc)) paste0(" + ", paste(mc, collapse=" + ")) else "", " + strata(tile)"))
    cf <- tryCatch(broom::tidy(coxph(fm, data=dm, ties="efron"), exponentiate=TRUE, conf.int=TRUE) %>% filter(term=="score_z"), error=function(e) tibble())
    tibble(tile_m = tile_m, min_n = min_n, Diet = lbl, N = nrow(dm), Tiles = nlevels(dm$tile), Spatial_ICC = icc,
           HR = if(nrow(cf)) cf$estimate else NA_real_, CI_low = if(nrow(cf)) cf$conf.low else NA_real_, CI_high = if(nrow(cf)) cf$conf.high else NA_real_)
  }))
}
geo_sens <- bind_rows(lapply(c(1000, 2000, 5000), function(ts) bind_rows(lapply(c(20, 30, 50), function(nn) geo_sensitivity_one(ts, nn)))))

p_sens1 <- ggplot(geo_sens, aes(factor(tile_m/1000), 100*Spatial_ICC, group=Diet, linetype=Diet, shape=Diet)) +
  geom_point() + geom_line() + facet_wrap(~min_n, labeller=labeller(min_n=function(x) paste0("minimum N = ", x))) +
  labs(title="a. Spatial ICC across neighborhood definitions", x="Grid width (km)", y="Residual ICC (%)", linetype=NULL, shape=NULL) + theme_classic(base_size=10)
p_sens2 <- ggplot(geo_sens, aes(factor(tile_m/1000), HR, group=Diet, linetype=Diet, shape=Diet)) +
  geom_hline(yintercept=1, linetype=3) + geom_point() + geom_line() + facet_wrap(~min_n, labeller=labeller(min_n=function(x) paste0("minimum N = ", x))) +
  labs(title="b. Within-neighborhood mortality HR across definitions", x="Grid width (km)", y="HR per 1 SD", linetype=NULL, shape=NULL) + theme_classic(base_size=10)
FigS7 <- p_sens1 / p_sens2 + plot_layout(guides="collect") & theme(legend.position="bottom")
save_plot(FigS7, "FigS7.geography_sensitivity.png", width=12.5, height=9, dpi=360, bg="white")
write_xlsx(list(tile_map = tile_map, score_structure = structure_tbl, within_neighborhood_mortality = within_tile, sensitivity = geo_sens), "Fig5.geography.out.xlsx")
write_xlsx(list(sensitivity = geo_sens), "FigS7.geography_sensitivity.out.xlsx")
print(structure_tbl); print(within_tile)
}
})

maha_run_step("survdiag", {
step_header("UKB follow-up and proportional-hazards diagnostics")
pacman::p_load(tidyverse, survival, patchwork, broom, writexl, scales)

diag_score_map <- c(
  DASH = "diet.dash.pts",
  MEDI = "diet.medi24.pts",
  MIND = "diet.mind.pts",
  MAHA = "diet.maha.pts"
)

diag_outcome_levels <- unname(dx.lst[Y.inc])
diag_outcome_cols <- setNames(
  c("#D24DC2", "#4A82E3", "#00A2A6", "#009C22", "#9A8500", "#D46057"),
  diag_outcome_levels
)

miss_diag_scores <- setdiff(unname(diag_score_map), names(dat))
if (length(miss_diag_scores)) {
  stop("Survival diagnostics missing dietary-score column(s): ",
       paste(miss_diag_scores, collapse = ", "), call. = FALSE)
}

reverse_km_quantiles <- function(time, event) {
  dd <- tibble(time = suppressWarnings(as.numeric(time)),
               event = suppressWarnings(as.integer(event))) %>%
    filter(is.finite(time), time > 0, event %in% c(0L, 1L))
  if (nrow(dd) < 10) return(c(Q25 = NA_real_, Median = NA_real_, Q75 = NA_real_))

  sf <- survival::survfit(survival::Surv(time, 1L - event) ~ 1, data = dd)

  qq <- tryCatch({
    qobj <- stats::quantile(sf, probs = c(.25, .50, .75), conf.int = FALSE)
    qraw <- if (is.list(qobj) && !is.null(qobj$quantile)) qobj$quantile else qobj
    qnum <- as.numeric(qraw)
    if (length(qnum) >= 3) qnum[seq_len(3)] else rep(NA_real_, 3)
  }, error = function(e) {
    message("[FOLLOW-UP WARNING] reverse-KM quantiles failed: ", conditionMessage(e))
    rep(NA_real_, 3)
  })

  if (length(qq) != 3) qq <- rep(NA_real_, 3)
  names(qq) <- c("Q25", "Median", "Q75")
  qq
}

followup_common <- purrr::map_dfr(Y.inc, function(Y) {
  time_var <- paste0(Y, ".t2e")
  event_var <- paste0(Y, ".Yt2e")
  need <- unique(c(time_var, event_var, unname(diag_score_map), covs))
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("Follow-up summary missing columns for ", Y, ": ", paste(miss, collapse = ", "), call. = FALSE)
  dd <- dat %>% dplyr::select(all_of(need)) %>% drop_na() %>%
    filter(is.finite(.data[[time_var]]), .data[[time_var]] > 0,
           .data[[event_var]] %in% c(0, 1))
  q <- reverse_km_quantiles(dd[[time_var]], dd[[event_var]])
  py <- sum(dd[[time_var]], na.rm = TRUE)
  nev <- sum(dd[[event_var]] == 1, na.rm = TRUE)
  tibble(
    Outcome_key = Y,
    Outcome = unname(dx.lst[Y]),
    N = nrow(dd),
    Events = nev,
    Event_percent = 100 * nev / nrow(dd),
    Person_years = py,
    Event_rate_per_1000PY = 1000 * nev / py,
    ReverseKM_Q25 = unname(q["Q25"]),
    ReverseKM_median = unname(q["Median"]),
    ReverseKM_Q75 = unname(q["Q75"]),
    Observed_followup_median = median(dd[[time_var]], na.rm = TRUE),
    Observed_followup_min = min(dd[[time_var]], na.rm = TRUE),
    Observed_followup_max = max(dd[[time_var]], na.rm = TRUE)
  )
}) %>%
  mutate(
    Outcome = factor(Outcome, levels = diag_outcome_levels),
    Followup_label = ifelse(
      is.finite(ReverseKM_median),
      sprintf("%.1f (%.1f–%.1f)", ReverseKM_median, ReverseKM_Q25, ReverseKM_Q75),
      sprintf("observed %.1f", Observed_followup_median)
    )
  )

followup_by_diet <- purrr::map_dfr(Y.inc, function(Y) {
  time_var <- paste0(Y, ".t2e")
  event_var <- paste0(Y, ".Yt2e")
  purrr::imap_dfr(diag_score_map, function(score_var, diet_lab) {
    need <- unique(c(time_var, event_var, score_var, covs))
    dd <- dat %>% dplyr::select(all_of(need)) %>% drop_na() %>%
      filter(is.finite(.data[[time_var]]), .data[[time_var]] > 0,
             .data[[event_var]] %in% c(0, 1))
    q <- reverse_km_quantiles(dd[[time_var]], dd[[event_var]])
    py <- sum(dd[[time_var]], na.rm = TRUE)
    nev <- sum(dd[[event_var]] == 1, na.rm = TRUE)
    tibble(
      Outcome_key = Y,
      Outcome = unname(dx.lst[Y]),
      Diet = diet_lab,
      N = nrow(dd), Events = nev,
      Person_years = py,
      ReverseKM_Q25 = unname(q["Q25"]),
      ReverseKM_median = unname(q["Median"]),
      ReverseKM_Q75 = unname(q["Q75"]),
      Observed_followup_median = median(dd[[time_var]], na.rm = TRUE),
      Observed_followup_max = max(dd[[time_var]], na.rm = TRUE)
    )
  })
})

safe_zph <- function(fit) {
  attempts <- list(
    list(terms = TRUE,  singledf = TRUE,  label = "terms_singledf"),
    list(terms = TRUE,  singledf = FALSE, label = "terms"),
    list(terms = FALSE, singledf = FALSE, label = "coefficients")
  )
  errs <- character()
  for (a in attempts) {
    z <- tryCatch(
      survival::cox.zph(fit, transform = "km", terms = a$terms,
                        singledf = a$singledf, global = TRUE),
      error = function(e) {
        errs <<- c(errs, paste0(a$label, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(z)) return(list(zph = z, method = a$label, error = paste(errs, collapse = " | ")))
  }
  message("[PH WARNING] cox.zph failed after all fallbacks: ", paste(errs, collapse = " | "))
  list(zph = NULL, method = NA_character_, error = paste(errs, collapse = " | "))
}

fit_ph_one <- function(Y, diet_label, score_var) {
  time_var <- paste0(Y, ".t2e")
  event_var <- paste0(Y, ".Yt2e")

  need <- unique(c(time_var, event_var, score_var, covs))
  dd <- dat %>%
    transmute(
      time = .data[[time_var]],
      event = .data[[event_var]],
      score = zstd(.data[[score_var]]),
      across(all_of(covs))
    ) %>%
    drop_na() %>%
    filter(is.finite(time), time > 0, event %in% c(0, 1)) %>%
    mutate(across(where(is.factor), droplevels))

  n_event <- sum(dd$event == 1)

  empty_summary <- function(reason = NA_character_) {
    tibble(
      Outcome_key = Y,
      Outcome = unname(dx.lst[Y]),
      Diet = diet_label,
      N = nrow(dd),
      Events = n_event,
      Exposure_PH_chisq = NA_real_,
      Exposure_PH_df = NA_real_,
      Exposure_PH_P = NA_real_,
      Global_PH_chisq = NA_real_,
      Global_PH_df = NA_real_,
      Global_PH_P = NA_real_,
      PH_method = NA_character_,
      PH_error = reason,
      PH_center_removed = NA,
      PH_covariate_set = NA_character_,
      PH_fallback_reason = reason
    )
  }

  if (nrow(dd) < 200 || n_event < 20) {
    return(list(
      summary = empty_summary("Too few observations/events for PH diagnostic."),
      detail = tibble()
    ))
  }

  parse_zph <- function(zph, covariate_set, center_removed, method, error_text,
                        fallback_reason = NA_character_) {
    if (is.null(zph)) return(NULL)

    ztab <- as.data.frame(zph$table) %>% rownames_to_column("Term")
    chisq_col <- intersect(c("chisq", "Chisq"), names(ztab))[1]
    p_col <- intersect(c("p", "P", "p.value"), names(ztab))[1]
    df_col <- intersect(c("df", "Df"), names(ztab))[1]

    getv <- function(row, col) {
      if (length(col) == 0 || is.na(col[1]) || nrow(row) == 0) return(NA_real_)
      suppressWarnings(as.numeric(row[[col[1]]][1]))
    }

    er <- ztab %>% filter(Term == "score")
    gr <- ztab %>% filter(Term == "GLOBAL")
    ep <- getv(er, p_col)
    gp <- getv(gr, p_col)

    detail <- ztab %>%
      mutate(
        Outcome_key = Y,
        Outcome = unname(dx.lst[Y]),
        Diet = diet_label,
        N = nrow(dd),
        Events = n_event,
        PH_center_removed = center_removed,
        PH_covariate_set = covariate_set,
        PH_method = method,
        PH_fallback_reason = fallback_reason
      ) %>%
      relocate(
        Outcome_key, Outcome, Diet, N, Events,
        PH_center_removed, PH_covariate_set, PH_method,
        PH_fallback_reason, Term
      )

    list(
      valid = is.finite(ep) && is.finite(gp),
      summary = tibble(
        Outcome_key = Y,
        Outcome = unname(dx.lst[Y]),
        Diet = diet_label,
        N = nrow(dd),
        Events = n_event,
        Exposure_PH_chisq = getv(er, chisq_col),
        Exposure_PH_df = getv(er, df_col),
        Exposure_PH_P = ep,
        Global_PH_chisq = getv(gr, chisq_col),
        Global_PH_df = getv(gr, df_col),
        Global_PH_P = gp,
        PH_method = method,
        PH_error = error_text,
        PH_center_removed = center_removed,
        PH_covariate_set = covariate_set,
        PH_fallback_reason = fallback_reason
      ),
      detail = detail
    )
  }

  run_one_ph <- function(rhs_covs, center_removed = FALSE,
                         fallback_reason = NA_character_) {
    rhs <- c("score", rhs_covs)
    fm <- as.formula(
      paste("Surv(time, event) ~", paste(rhs, collapse = " + "))
    )

    fit <- tryCatch(
      survival::coxph(
        fm, data = dd, ties = "efron",
        x = TRUE, model = TRUE, singular.ok = TRUE
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(
        valid = FALSE,
        summary = NULL,
        detail = tibble(),
        error = paste0("coxph: ", conditionMessage(fit))
      ))
    }

    zres <- safe_zph(fit)
    if (is.null(zres$zph)) {
      return(list(
        valid = FALSE,
        summary = NULL,
        detail = tibble(),
        error = zres$error
      ))
    }

    parsed <- parse_zph(
      zres$zph,
      covariate_set = if (center_removed) "Primary covariates minus center" else "Primary covariates",
      center_removed = center_removed,
      method = zres$method,
      error_text = zres$error,
      fallback_reason = fallback_reason
    )

    if (is.null(parsed)) {
      return(list(
        valid = FALSE,
        summary = NULL,
        detail = tibble(),
        error = "cox.zph returned no usable table."
      ))
    }

    if (!parsed$valid) {
      parsed$error <- paste0(
        "cox.zph returned non-finite dietary-score and/or GLOBAL P value. ",
        "Exposure P=", parsed$summary$Exposure_PH_P,
        "; Global P=", parsed$summary$Global_PH_P
      )
    } else {
      parsed$error <- NA_character_
    }
    parsed
  }

  primary <- run_one_ph(covs, center_removed = FALSE)

  if (isTRUE(primary$valid)) {
    return(list(summary = primary$summary, detail = primary$detail))
  }

  fallback_covs <- setdiff(covs, "center")
  fallback_reason <- paste0(
    "Primary PH diagnostic not estimable; reran without center. Primary error: ",
    ifelse(is.null(primary$error) || is.na(primary$error), "unknown", primary$error)
  )

  message(
    "[PH FALLBACK] ", unname(dx.lst[Y]), " : ", diet_label,
    " -- retrying without center."
  )

  fallback <- run_one_ph(
    fallback_covs,
    center_removed = TRUE,
    fallback_reason = fallback_reason
  )

  if (isTRUE(fallback$valid)) {
    return(list(summary = fallback$summary, detail = fallback$detail))
  }

  both_errors <- paste0(
    fallback_reason,
    " | No-center fallback error: ",
    ifelse(is.null(fallback$error) || is.na(fallback$error), "unknown", fallback$error)
  )

  message(
    "[PH WARNING] ", unname(dx.lst[Y]), " : ", diet_label,
    " -- still not estimable after removing center."
  )

  list(
    summary = empty_summary(both_errors) %>%
      mutate(
        PH_center_removed = TRUE,
        PH_covariate_set = "Primary covariates minus center",
        PH_fallback_reason = fallback_reason
      ),
    detail = tibble()
  )
}

ph_results <- list()
for (Y in Y.inc) {
  for (diet_lab in names(diag_score_map)) {
    message("[PH] ", unname(dx.lst[Y]), " : ", diet_lab)
    ph_results[[paste(Y, diet_lab, sep = "__")]] <-
      fit_ph_one(Y, diet_lab, diag_score_map[[diet_lab]])
  }
}

PH_summary <- bind_rows(lapply(ph_results, `[[`, "summary")) %>%
  mutate(
    Outcome = factor(Outcome, levels = diag_outcome_levels),
    Diet = factor(Diet, levels = names(diag_score_map)),
    Exposure_PH_flag = ifelse(is.na(Exposure_PH_P), NA, Exposure_PH_P < .05),
    Global_PH_flag = ifelse(is.na(Global_PH_P), NA, Global_PH_P < .05)
  )
PH_all_terms <- bind_rows(lapply(ph_results, `[[`, "detail"))
PH_flagged <- PH_summary %>%
  filter((!is.na(Exposure_PH_P) & Exposure_PH_P < .05) |
         (!is.na(Global_PH_P) & Global_PH_P < .05)) %>%
  arrange(Exposure_PH_P, Global_PH_P)

fit_time_split_one <- function(Y, diet_label, score_var, cut_year = 5) {
  time_var <- paste0(Y, ".t2e")
  event_var <- paste0(Y, ".Yt2e")
  need <- unique(c(time_var, event_var, score_var, covs))
  dd <- dat %>%
    transmute(
      time = .data[[time_var]], event = .data[[event_var]],
      score = zstd(.data[[score_var]]), across(all_of(covs))
    ) %>%
    drop_na() %>%
    filter(is.finite(time), time > 0, event %in% c(0, 1)) %>%
    mutate(across(where(is.factor), droplevels), id = row_number())

  if (nrow(dd) < 200 || sum(dd$event == 1) < 20) return(tibble())

  ds <- survival::survSplit(
    Surv(time, event) ~ ., data = dd, cut = cut_year,
    start = "tstart", end = "time", event = "event", episode = "period"
  ) %>%
    mutate(
      score_early = ifelse(period == 1, score, 0),
      score_late  = ifelse(period >= 2, score, 0)
    )

  fm <- as.formula(paste0(
    "Surv(tstart, time, event) ~ score_early + score_late + ",
    paste(covs, collapse = " + "), " + strata(period)"
  ))
  fit <- tryCatch(
    survival::coxph(fm, data = ds, ties = "efron", singular.ok = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble())

  b <- stats::coef(fit); V <- stats::vcov(fit)
  if (!all(c("score_early", "score_late") %in% names(b))) return(tibble())
  se1 <- sqrt(V["score_early", "score_early"])
  se2 <- sqrt(V["score_late", "score_late"])
  dlt <- b["score_late"] - b["score_early"]
  sed <- sqrt(V["score_late", "score_late"] + V["score_early", "score_early"] -
              2 * V["score_late", "score_early"])
  pint <- if (is.finite(sed) && sed > 0) 2 * pnorm(-abs(dlt / sed)) else NA_real_

  tibble(
    Outcome_key = Y, Outcome = unname(dx.lst[Y]), Diet = diet_label,
    N = nrow(dd), Events = sum(dd$event == 1), Cut_year = cut_year,
    HR_0_5 = exp(b["score_early"]),
    LCI_0_5 = exp(b["score_early"] - 1.96 * se1),
    UCI_0_5 = exp(b["score_early"] + 1.96 * se1),
    HR_after5 = exp(b["score_late"]),
    LCI_after5 = exp(b["score_late"] - 1.96 * se2),
    UCI_after5 = exp(b["score_late"] + 1.96 * se2),
    P_time_interaction = pint
  )
}

time_split_sensitivity <- purrr::map_dfr(Y.inc, function(Y) {
  purrr::imap_dfr(diag_score_map, function(score_var, diet_lab) {
    fit_time_split_one(Y, diet_lab, score_var, cut_year = 5)
  })
})

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < .001, "<0.001", sprintf("%.3f", x)))
}
ph_class <- function(x) {
  factor(case_when(
    is.na(x) ~ "Not estimable",
    x < .01 ~ "P < 0.01",
    x < .05 ~ "P = 0.01–0.049",
    TRUE ~ "P ≥ 0.05"
  ), levels = c("P ≥ 0.05", "P = 0.01–0.049", "P < 0.01", "Not estimable"))
}

p_fu <- ggplot(followup_common, aes(x = ReverseKM_median, y = Outcome, color = Outcome)) +
  geom_segment(aes(x = ReverseKM_Q25, xend = ReverseKM_Q75, yend = Outcome), linewidth = .9) +
  geom_point(size = 2.8) +
  geom_text(aes(label = Followup_label), hjust = -0.08, vjust = 0,
            nudge_y = .20, size = 3.1, fontface = "bold") +
  scale_color_manual(values = diag_outcome_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.02, .28))) +
  scale_y_discrete(expand = expansion(add = c(.55, .78))) +
  labs(title = "a. Follow-up duration", subtitle = "Reverse Kaplan–Meier median (IQR); common complete-case sample",
       x = "Follow-up time (years)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))

maha_split_plot <- time_split_sensitivity %>%
  filter(Diet == "MAHA") %>%
  mutate(Outcome = factor(Outcome, levels = diag_outcome_levels)
  ) %>%
  pivot_longer(
    cols = c(HR_0_5, HR_after5), names_to = "Period", values_to = "HR"
  ) %>%
  mutate(
    LCI = ifelse(Period == "HR_0_5", LCI_0_5, LCI_after5),
    UCI = ifelse(Period == "HR_0_5", UCI_0_5, UCI_after5),
    Period = recode(Period, HR_0_5 = "0–5 years", HR_after5 = ">5 years")
  )

p_split <- ggplot(maha_split_plot, aes(x = HR, y = Outcome, color = Period)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = .16,
                 position = position_dodge(width = .42), linewidth = .7) +
  geom_point(position = position_dodge(width = .42), size = 2.5) +
  labs(title = "b. MAHA time-period sensitivity",
       subtitle = "Separate HRs during years 0–5 and >5",
       x = "Hazard ratio per 1 SD higher MAHA", y = NULL, color = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"),
        legend.position = "bottom")

p_ph_exp <- PH_summary %>%
  mutate(
    Status = ph_class(Exposure_PH_P),
    P_label = paste0(fmt_p(Exposure_PH_P), ifelse(PH_center_removed %in% TRUE, "†", ""))
  ) %>%
  ggplot(aes(x = Diet, y = Outcome, fill = Status)) +
  geom_tile(color = "white", linewidth = .45) +
  geom_text(aes(label = P_label), size = 3.0, fontface = "bold") +
  labs(title = "c. Dietary-score PH test", subtitle = "Schoenfeld test; † center omitted only when full-model test was not estimable",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "bottom")

p_ph_global <- PH_summary %>%
  mutate(
    Status = ph_class(Global_PH_P),
    P_label = paste0(fmt_p(Global_PH_P), ifelse(PH_center_removed %in% TRUE, "†", ""))
  ) %>%
  ggplot(aes(x = Diet, y = Outcome, fill = Status)) +
  geom_tile(color = "white", linewidth = .45) +
  geom_text(aes(label = P_label), size = 3.0, fontface = "bold") +
  labs(title = "d. Global PH test", subtitle = "Schoenfeld test; † center omitted only when full-model test was not estimable",
       x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), legend.position = "bottom")

Fig_survdiag <- (p_fu | p_split) / (p_ph_exp | p_ph_global)
save_plot(Fig_survdiag, "FigS5.survival_diagnostics.png", width = 15.5, height = 11.0, dpi = 500, bg = "white")
write_xlsx(
  list(
    followup_summary = followup_common %>% mutate(Outcome = as.character(Outcome)),
    followup_by_diet = followup_by_diet,
    PH_summary = PH_summary %>% mutate(Outcome = as.character(Outcome), Diet = as.character(Diet)),
    PH_flagged = PH_flagged %>% mutate(Outcome = as.character(Outcome), Diet = as.character(Diet)),
    PH_all_terms = PH_all_terms,
    time_split_sensitivity = time_split_sensitivity
  ),
  "FigS5.survival_diagnostics.out.xlsx"
)

print(followup_common)
print(PH_summary)
if (nrow(PH_flagged)) {
  message("[PH] Models with exposure-specific or global P < 0.05:")
  print(PH_flagged)
} else {
  message("[PH] No exposure-specific or global PH test had P < 0.05.")
}
})

maha_run_step("figS1", {
step_header("Figure S3: diet score distributions and correlation")
pacman::p_load(tidyverse, patchwork, cowplot)
dat1 <- dat

figS1_score_map <- c(
	maha = "MAHA", dash = "DASH", mind = "MIND", medi24 = "MEDI",
	maha_bal = "MAHA-balanced", maha_strict = "MAHA-strict",
	maha_nodairy = "MAHA-no dairy", maha_noprotein = "MAHA-no protein"
)
figS1_cols <- c(
	"MAHA" = "#F26D60", "DASH" = "#16B9C0", "MIND" = "#B97AF7", "MEDI" = "#7CAE00",
	"MAHA-balanced" = "#E76F51", "MAHA-strict" = "#C77DFF",
	"MAHA-no dairy" = "#2A9D8F", "MAHA-no protein" = "#577590"
)
figS1_dist_axis_text_size <- 8.2
figS1_top_height <- 0.76 * 1.30
figS1_gap_height <- 0.10
figS1_bottom_height <- 1.24
figS1_plot_height <- 9.6 * (figS1_top_height + figS1_gap_height + figS1_bottom_height) / (0.76 + 1.24)

for (nm in names(figS1_score_map)) {
	sumv <- paste0("diet.", nm, ".sum"); s100v <- paste0("diet.", nm, ".s100")
	if (!s100v %in% names(dat1) && sumv %in% names(dat1)) dat1[[s100v]] <- score_0_100(dat1[[sumv]])
}

figS1_vars <- paste0("diet.", names(figS1_score_map), ".s100")
missing_figS1 <- setdiff(figS1_vars, names(dat1))
if (length(missing_figS1)) stop("Missing FigS1 score variables: ", paste(missing_figS1, collapse = ", "))
figS1_var_labs <- setNames(unname(figS1_score_map), figS1_vars)

plot_dat <- dat1 %>%
	dplyr::select(all_of(figS1_vars)) %>%
	pivot_longer(everything(), names_to = "var", values_to = "score") %>%
	mutate(score = as.numeric(score), Diet_score = factor(figS1_var_labs[var], levels = rev(unname(figS1_score_map)))) %>%
	filter(is.finite(score), !is.na(Diet_score))

plot_sum <- plot_dat %>%
	group_by(Diet_score) %>%
	summarise(
		N = n(), Mean = mean(score), SD = sd(score), Median = median(score),
		P25 = quantile(score, 0.25) %>% as.numeric(), P75 = quantile(score, 0.75) %>% as.numeric(),
		Min = min(score), Max = max(score), .groups = "drop"
	) %>%
	mutate(across(c(Mean, SD, Median, P25, P75, Min, Max), ~ round(.x, 1))) %>%
	as.data.frame()
density_dat <- plot_dat %>% mutate(Diet_score = as.character(Diet_score)) %>%
	group_by(Diet_score) %>% group_modify(~ {
		dn <- density(.x$score, from = 0, to = 100, n = 256, na.rm = TRUE)
		tibble(score = dn$x, density = dn$y)
	}) %>% ungroup()

p_dist_ukb <- ggplot(plot_dat, aes(x = score, y = Diet_score, fill = Diet_score)) +
	geom_violin(width = 0.82, scale = "width", trim = TRUE, alpha = 0.72, color = NA) +
	geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.92, color = "grey25") +
	scale_fill_manual(values = figS1_cols, guide = "none") +
	scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0.01, 0.01))) +
	labs(x = "Harmonized dietary score (0-100)", y = NULL, title = "a. Score distribution (UKB)") +
	theme_classic(base_size = 10) +
	theme(
		legend.position = "none",
		plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
		axis.text.y = element_text(face = "bold", size = figS1_dist_axis_text_size),
		axis.text.x = element_text(face = "bold", size = figS1_dist_axis_text_size),
		axis.title.x = element_text(face = "bold", size = figS1_dist_axis_text_size),
		plot.margin = margin(5.5, 12, 5.5, 5.5)
	)

score_mat <- dat1 %>% dplyr::select(all_of(figS1_vars)) %>% mutate(across(everything(), as.numeric))
cor_mat <- cor(score_mat, use = "pairwise.complete.obs", method = "spearman")
dimnames(cor_mat) <- rep(list(unname(figS1_score_map)), 2)
cor_dat <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE) %>%
	mutate(Var1 = factor(Var1, levels = rev(unname(figS1_score_map))), Var2 = factor(Var2, levels = unname(figS1_score_map)), lab = sprintf("%.2f", Freq))

p_cor_ukb <- ggplot(cor_dat, aes(Var2, Var1, fill = Freq)) +
	geom_tile(color = "white", linewidth = 0.35) +
	geom_text(aes(label = lab), size = 2.55, fontface = "bold") +
	scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.5, 1), oob = scales::squish, name = NULL) +
	coord_fixed() +
	labs(x = NULL, y = NULL, title = "c. Spearman correlation (UKB)") +
	theme_minimal(base_size = 10) +
	theme(
		plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
		axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, face = "bold", size = 8.2),
		axis.text.y = element_text(face = "bold", size = 8.2),
		panel.grid = element_blank(),
		legend.position = "right",
		legend.key.height = unit(0.55, "cm"),
		plot.margin = margin(5.5, 5.5, 5.5, 5.5)
	)

FigS3 <- p_dist_ukb / plot_spacer() / p_cor_ukb +
  plot_layout(heights = c(figS1_top_height, figS1_gap_height, figS1_bottom_height), guides = "collect") & theme(legend.position = "top")
save_plot(FigS3, "FigS3.score_distribution_concordance.png", width = 7.2, height = figS1_plot_height, dpi = 500, bg = "white")
write_xlsx(list(score_distribution = plot_sum, density_curve = density_dat,
		spearman_matrix = as.data.frame(cor_mat, check.names = FALSE)),
           "FigS3.score_distribution_concordance.out.xlsx")
print(plot_sum)
print(round(cor_mat, 2))
})

maha_run_step("figS2", {
step_header("eFigure 2: head-to-head comparison of PheWAS")
pacman::p_load(tidyverse, ggrepel, patchwork)

dm2 <- dm %>% transmute(phenotype, description, beta.MAHA = as.numeric(beta.MAHA), beta.DASH = as.numeric(beta.DASH), logp.MAHA = as.numeric(logp.MAHA), logp.DASH = as.numeric(logp.DASH), delta_logp = as.numeric(delta_logp), delta_beta = as.numeric(delta_beta), sig_pattern = recode(as.character(sig_pattern), DASH_only = "DASH only", MAHA_only = "MAHA only")) %>% filter(if_all(c(beta.MAHA, beta.DASH, logp.MAHA, logp.DASH), is.finite)) %>% mutate(sig_pattern = ifelse(sig_pattern %in% c("Both", "DASH only", "MAHA only"), sig_pattern, "Neither"), logp.MAHA.cap = pmin(logp.MAHA, 15), logp.DASH.cap = pmin(logp.DASH, 15))
bonf <- -log10(0.05 / nrow(dm2)); cols <- c("Both" = "#3B82F6", "DASH only" = "#F97316", "MAHA only" = "#10B981")
lab_p <- bind_rows(dm2 %>% filter(sig_pattern == "DASH only") %>% slice_max(delta_logp, n = 10), dm2 %>% filter(sig_pattern == "MAHA only") %>% slice_min(delta_logp, n = 4), dm2 %>% filter(sig_pattern == "Both") %>% slice_max(abs(delta_logp), n = 4)) %>% distinct(phenotype, .keep_all = TRUE)
lab_b <- bind_rows(dm2 %>% filter(sig_pattern == "DASH only") %>% slice_max(abs(delta_beta), n = 3), dm2 %>% filter(sig_pattern == "MAHA only") %>% slice_max(abs(delta_beta), n = 2), dm2 %>% filter(sig_pattern == "Both") %>% slice_max(abs(delta_beta), n = 2)) %>% distinct(phenotype, .keep_all = TRUE)
brange <- quantile(c(dm2$beta.DASH, dm2$beta.MAHA), c(0.003, 0.997), na.rm = TRUE, names = FALSE)
brange <- brange + c(-1, 1) * diff(brange) * 0.06
beta_axis_lab <- function(x) ifelse(abs(x) < 1e-12, "0", formatC(x, format = "e", digits = 0))
thm <- theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), legend.position = "top", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 14))

p1 <- ggplot() +
	geom_point(data = dm2 %>% filter(sig_pattern == "Neither"), aes(logp.DASH.cap, logp.MAHA.cap), color = "grey70", alpha = 0.35, size = 1) +
	geom_point(data = dm2 %>% filter(sig_pattern != "Neither"), aes(logp.DASH.cap, logp.MAHA.cap, color = sig_pattern, shape = sig_pattern), alpha = 0.9, size = 2.3) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") + geom_vline(xintercept = bonf, linetype = 3, color = "grey60") + geom_hline(yintercept = bonf, linetype = 3, color = "grey60") +
	geom_label_repel(data = lab_p, aes(logp.DASH.cap, logp.MAHA.cap, label = description, color = sig_pattern),
		seed = 20260830, size = 3.05, fontface = "bold", fill = scales::alpha("white", .84),
		label.size = .12, label.padding = unit(.08, "lines"), box.padding = .55,
		point.padding = .24, force = 3.0, force_pull = .08, max.time = 12,
		max.iter = 200000, min.segment.length = 0, segment.alpha = .55,
		max.overlaps = Inf, show.legend = FALSE) +
	scale_color_manual(values = cols) + scale_shape_manual(values = c("Both" = 16, "DASH only" = 17, "MAHA only" = 15)) +
	coord_cartesian(xlim = c(0, 16), ylim = c(0, 16), expand = FALSE, clip = "off") + scale_x_continuous(breaks = seq(0, 15, 3)) + scale_y_continuous(breaks = seq(0, 15, 3)) +
	labs(x = "DASH: -log10(P)", y = "MAHA: -log10(P)", title = "a. Association strength") + thm

p2 <- ggplot() +
	geom_point(data = dm2 %>% filter(sig_pattern == "Neither"), aes(beta.DASH, beta.MAHA), color = "grey70", alpha = 0.35, size = 1) +
	geom_point(data = dm2 %>% filter(sig_pattern != "Neither"), aes(beta.DASH, beta.MAHA, color = sig_pattern, shape = sig_pattern), alpha = 0.9, size = 2.3) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") + geom_hline(yintercept = 0, linetype = 3, color = "grey60") + geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
	geom_label_repel(data = lab_b, aes(beta.DASH, beta.MAHA, label = description, color = sig_pattern),
		seed = 20260830, size = 2.75, fontface = "bold", fill = scales::alpha("white", .84),
		label.size = .10, label.padding = unit(.07, "lines"), box.padding = .58,
		point.padding = .28, force = 2.6, force_pull = .10, max.time = 12,
		max.iter = 200000, min.segment.length = 0, segment.alpha = .55,
		max.overlaps = Inf, show.legend = FALSE) +
	scale_color_manual(values = cols) + scale_shape_manual(values = c("Both" = 16, "DASH only" = 17, "MAHA only" = 15)) +
	coord_cartesian(xlim = brange, ylim = brange, expand = FALSE) +
	scale_x_continuous(labels = beta_axis_lab) + scale_y_continuous(labels = beta_axis_lab) +
	labs(x = "DASH: beta", y = "MAHA: beta", title = "b. Effect size") + thm

FigS2 <- (p1 / plot_spacer() / p2) + plot_layout(heights = c(1, 0.08, 1), guides = "collect") & theme(legend.position = "top"); FigS2
save_plot(FigS2, "FigS2.phewas_head2head.png", width = 16, height = 16, dpi = 500, bg = "white")
write_xlsx(list(phewas_compare = dm2, label_p = lab_p, label_beta = lab_b), "FigS2.phewas_head2head.out.xlsx")
print(dm2, n = 30)
print(lab_p)
print(lab_b)
})

maha_run_step("figS3", {
step_header("eFigure 4: 10-year joint DASH and MAHA categories")

tab2.risk <- bind_rows(lapply(Y.inc, function(Y) run_joint_risk_one(dat, Y, covs)))
tab2.risk.out <- tab2.risk %>% mutate(Outcome = dx.lst[Y]) %>% dplyr::select(Outcome, dash, maha, risk10)
write_xlsx(tab2.risk.out, "FigS4.joint_risk.out.xlsx")

plist.risk <- lapply(Y.inc, function(Y) plot_joint_heat(tab2.risk %>% filter(Y == !!Y), dx.lst[Y]))
FigS3 <- wrap_plots(plist.risk, ncol = 2) + plot_layout(ncol = 2)
save_plot(FigS3, "FigS4.joint_risk.png", width = 8, height = 8, dpi = 320)
print(tab2.risk.out)
})

maha_run_step("figS4", {
step_header("Figure S6: MAHA and DASH components")

if (!exists("field")) field <- readxl::read_excel(paste0(dir0, "/files/foods.xlsx"), sheet = "Sheet2")
food.item.labels <- readxl::read_excel(paste0(dir0, "/files/foods.xlsx"), sheet = "Sheet1") %>%
	transmute(item = paste0("p", as.character(FieldID)), food_item_label = trimws(as.character(Food.Item))) %>%
	filter(!is.na(food_item_label), food_item_label != "") %>% distinct(item, .keep_all = TRUE)

diet.dict <- read.delim(
	ifelse(file.exists(paste0(indir, "/common/diet.lst")), paste0(indir, "/common/diet.lst"), "diet.lst"),
	sep = "\t", header = FALSE, fill = TRUE, quote = "", comment.char = ""
) %>%
	as_tibble() %>%
	transmute(item = trimws(V1), item_label = gsub("_", " ", trimws(V2), fixed = TRUE)) %>%
	distinct(item, .keep_all = TRUE)

dairy.vars <- c("p26150", "p26096", "p102820", "p102830", "p102840", "p102860", "p102880", "p102890", "p102900", "p102910")
veg.vars <- c("p104090", "p104240", "p104370", "p104140", "p104160", "p104180", "p104300", "p104310", "p104220", "p104230", "p104260", "p104340", "p104350", "p104150", "p104170", "p104290", "p104330", "p104130", "p104190", "p104270", "p104360", "p104060", "p104070", "p104200", "p104210", "p104250", "p104320", "p104380")
fruit.vars <- c("p104450", "p104560", "p104470", "p104480", "p104490", "p104530", "p104540", "p104430", "p104420", "p104550", "p104410", "p104580", "p104500", "p104460", "p104510", "p104570", "p104520", "p104440", "p104590")
fat.vars <- c("p26110", "oil", "p26062", "p26063", "p102440", "p102450", "p103160", "p104100")
wholegrain.vars <- c("p100770", "p100810", "p100800", "p100840", "p100850", "bread.wholemeal", "baguette.wholemeal", "bap.wholemeal", "roll.wholemeal", "bread.seeded", "baguette.seeded", "bap.seeded", "roll.seeded", "p102720", "p102740", "p102780")
refined.vars <- c("p100820", "p101230", "p101240", "p101250", "p101260", "p101270", "p102710", "p102730", "p102750", "p102760", "p102770", "bread.white", "baguette.white", "bap.white", "roll.white", "p100830", "p100860")
ssb.vars <- c("p100170", "p100180", "p100550", "p100160", "p100540", "p100530", "p100230", "p100190", "p100200", "p100210", "p100220")
sweet.vars <- c("p26064", "p102120", "p102140", "p102150", "p102170", "p102180", "p102190", "p102200", "p102210", "p102220", "p102230", "p102260", "p102270", "p102280", "p102290", "p102300", "p102310", "p102320", "p102330", "p102340", "p102350", "p102360", "p102370", "p102380", "p102060", "p102010", "p102020", "p102050", "p102070", "p101970", "p101980", "p101990", "p102030")
hp.vars <- c("p102460", "p102470", "p102480", "p102500", "p102040", "p102000", "p104020", "p103050", "p103170", "p103180", "p103010", "p103070", "p103080")
alcohol.vars <- c("p26067", "p26138", "p26151", "p26152", "p26153")
meat.vars <- c("p103010", "p103050", "p103070", "p103080")

raw.items <- unique(c(dairy.vars, veg.vars, fruit.vars, fat.vars, wholegrain.vars, refined.vars, ssb.vars, sweet.vars, hp.vars, alcohol.vars))
raw.items <- raw.items[grepl("^p\\d+$", raw.items)]
derived.items <- unique(c(setdiff(c(dairy.vars, veg.vars, fruit.vars, fat.vars, wholegrain.vars, refined.vars, ssb.vars, sweet.vars, hp.vars, alcohol.vars), raw.items), "w.p26005", "w.p26052"))

diet.raw <- readRDS(paste0(indir, "/Rdata/diet0.rds")) %>% mutate(eid = as.character(eid))
dat.food <- readRDS(paste0(indir, "/Rdata/diet.fudan.rds")) %>% mutate(eid = as.character(eid))
need.raw <- c("eid", grep(paste0("^(", paste(raw.items, collapse = "|"), ")_i[0-4]$"), names(diet.raw), value = TRUE))

dat.raw <- diet.raw[, need.raw, drop = FALSE] %>% dplyr::select(eid)
for (v in raw.items) dat.raw[[v]] <- collapse1(diet.raw, v)

dat.item0 <- dat.raw %>% left_join(dat.food %>% dplyr::select(eid, any_of(derived.items), energy.kJ), by = "eid")

dash.items <- unique(c(field %>% filter(DASH > 0) %>% pull(name), "w.p26052"))
maha.items <- unique(c(dairy.vars, veg.vars, fruit.vars, fat.vars, wholegrain.vars, refined.vars, ssb.vars, sweet.vars, hp.vars, alcohol.vars, "w.p26005", "w.p26052"))

extra.dict <- tribble(
	~item, ~item_label,
	"oil", "Healthy oil", "bread.white", "White bread", "bread.wholemeal", "Wholemeal bread",
	"bread.seeded", "Seeded bread", "baguette.white", "White baguette", "baguette.wholemeal", "Wholemeal baguette",
	"bap.white", "White bap", "bap.wholemeal", "Wholemeal bap", "roll.white", "White roll",
	"roll.wholemeal", "Wholemeal roll", "w.p26005", "Protein intake", "w.p26052", "Sodium"
)

item.map <- tibble(item = union(raw.items, derived.items)) %>%
	filter(item %in% names(dat.item0)) %>%
	mutate(
		in_dash = item %in% dash.items,
		in_maha = item %in% maha.items,
		origin = case_when(in_dash & in_maha ~ "Both", in_dash ~ "DASH only", in_maha ~ "MAHA only", TRUE ~ "Other"),
		origin2 = case_when(
			item %in% alcohol.vars ~ "Alcohol",
			item %in% dairy.vars ~ "Dairy",
			item %in% meat.vars ~ "Meat",
			item %in% ssb.vars ~ "Drinks",
			item %in% refined.vars ~ "Refined grains",
			item %in% c(sweet.vars, setdiff(hp.vars, meat.vars)) ~ "Sweets/processed foods",
			item %in% fruit.vars ~ "Fruit",
			item %in% veg.vars ~ "Vegetables",
			item %in% wholegrain.vars ~ "Whole grains",
			item %in% fat.vars ~ "Fats/oils",
			item == "w.p26005" ~ "Protein",
			item == "w.p26052" ~ "Sodium",
			TRUE ~ "Other"
		)
	) %>%
	left_join(diet.dict, by = "item") %>%
	left_join(food.item.labels, by = "item") %>%
	left_join(extra.dict, by = "item", suffix = c("", ".x")) %>%
	mutate(item_label = coalesce(item_label.x, food_item_label, item_label, item)) %>%
	dplyr::select(-item_label.x, -food_item_label)

dat.item <- dat %>%
	transmute(eid = as.character(eid), maha = diet.maha.sum, dash = diet.dash.sum, across(any_of(covs))) %>%
	mutate(
		maha.hml = mk_hml(maha, group_pct_th),
		dash.hml = mk_hml(dash, group_pct_th),
		disc.grp = case_when(maha.hml == "high" & dash.hml == "low" ~ "MH_DL", maha.hml == "low" & dash.hml == "high" ~ "ML_DH", TRUE ~ NA_character_)
	) %>%
	filter(!is.na(disc.grp)) %>%
	left_join(dat.item0, by = "eid") %>%
	mutate(y_disc = as.integer(disc.grp == "MH_DL"))

res.item.disc <- purrr::map_dfr(item.map$item, fit_disc, dat.item = dat.item, covs = covs) %>%
	left_join(item.map, by = "item") %>%
	filter(is.finite(beta), is.finite(conf.low), is.finite(conf.high)) %>%
	arrange(desc(abs(beta)))

top.items <- res.item.disc %>%
	filter(n >= 500, abs(beta) < log(5)) %>%
	arrange(desc(abs(beta))) %>%
	slice_head(n = 20) %>%
	pull(item)

if (length(top.items) == 0) top.items <- res.item.disc %>% arrange(desc(abs(beta))) %>% slice_head(n = 20) %>% pull(item)

res.item.cox <- purrr::map_dfr(top.items, \(v) purrr::map_dfr(Y.inc, \(Y) fit_cox(v, Y, dat.item0, covs))) %>%
	left_join(item.map, by = "item") %>%
	filter(is.finite(beta), is.finite(conf.low), is.finite(conf.high))

origin2.levels <- c("Alcohol", "Dairy", "Drinks", "Fats/oils", "Fruit", "Vegetables", "Meat", "Protein", "Refined grains", "Whole grains", "Sweets/processed foods", "Sodium", "Other")

item.class <- item.map %>%
	filter(item %in% top.items) %>%
	dplyr::select(item, origin2) %>%
	mutate(origin2 = factor(origin2, levels = origin2.levels))

ord.tbl <- res.item.disc %>%
	filter(item %in% top.items) %>%
	dplyr::select(item, item_label, beta.disc = beta) %>%
	distinct() %>%
	left_join(item.class, by = "item") %>%
	arrange(origin2, desc(beta.disc)) %>%
	mutate(item_label = factor(item_label, levels = rev(unique(item_label))))

res.item.cox2 <- res.item.cox %>%
	filter(item %in% top.items) %>%
	dplyr::select(item, Outcome, HR, p.value, beta, item_label) %>%
	left_join(ord.tbl %>% dplyr::select(item, item_label.ord = item_label, origin2), by = "item") %>%
	mutate(
		item_label = item_label.ord,
		origin2 = factor(as.character(origin2), levels = origin2.levels),
		Outcome = factor(Outcome, levels = dx.lst[Y.inc]),
		lab = ifelse(p.value < 0.05, sprintf("%.2f", HR), "")
	)

p4 <- ggplot(res.item.cox2, aes(Outcome, item_label, fill = beta)) +
	geom_tile(color = "white", linewidth = 0.5) +
	geom_text(aes(label = lab), size = 3.0) +
	scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
	facet_grid(origin2 ~ ., scales = "free_y", space = "free_y", switch = "y", drop = TRUE) +
	labs(title = "Dietary components underlying MAHA-DASH discordance", x = NULL, y = NULL, fill = "log(HR)") +
	theme_bw(base_size = 13) +
	theme(
		strip.placement = "outside",
		strip.background = element_rect(fill = "grey95", color = NA),
		strip.text.y.left = element_text(angle = 0, face = "bold"),
		axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
		plot.title = element_text(face = "bold", hjust = .5, size = 15, margin = margin(b = 8)),
		panel.spacing.y = unit(0.15, "lines")
	)

save_plot(p4, "FigS6.component_analysis.png", width = 10.5, height = 8.4, dpi = 320)

write_xlsx(
	list(
		item_map = item.map,
		item_discordance = res.item.disc,
		top_items = data.frame(item = top.items),
		ordering_table = ord.tbl %>% as.data.frame(),
		top_item_outcome = res.item.cox,
		top_item_outcome_heatmap = res.item.cox2 %>% as.data.frame()
	),
	"FigS6.component_analysis.out.xlsx"
)

print(head(res.item.disc, 50))
print(res.item.cox, n = 50)
})
