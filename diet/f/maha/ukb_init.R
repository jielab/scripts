# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", os_path("DIR_DAT"), os_path("DIR_DAT"))
indir <- Sys.getenv("DIET_UKB", os_path("DIR_UKB"))
maha_outdir <- Sys.getenv("MAHA_OUTDIR", unset = file.path(Sys.getenv("DIR_DAT", if (.Platform$OS.type == "windows") os_path("DIR_DAT") else os_path("DIR_DAT")), "analysis", "maha"))
.this_file <- maha_runner
.this_dir <- maha_code
.helper_dir <- Sys.getenv("MAHA_HELPER_DIR", unset = file.path(dir0, "scripts", "0f"))

.bootstrap_files <- c(
	file.path(.helper_dir, "phenotype.R"),
	file.path(.helper_dir, "association.R"),
	file.path(.helper_dir, "prediction.R"),
	file.path(.helper_dir, "plotting.R"),
	file.path(maha_code, "common.R")
)
.bootstrap_missing <- .bootstrap_files[!file.exists(.bootstrap_files)]
if (length(.bootstrap_missing)) {
	cat("[UKB INPUT CHECK] Required code/helper file(s) do not exist:
", paste0("  MISSING: ", .bootstrap_missing, collapse = "
"), "
", sep = "")
	stop("UKB bootstrap input check failed.", call. = FALSE)
}

maha_load_packages(writexl, tidyverse, survival, flextable, patchwork, gtsummary, broom, cowplot)
invisible(lapply(c("phenotype.R", "association.R", "prediction.R", "plotting.R"), function(f) source(file.path(.helper_dir, f))))
source(file.path(maha_code, "common.R"))


# 🚩 ukb_followup
# MAHA-specific timing contract. No outcome-dependent choice of landmark.
maha_followup_version <- "2026-09-06-webq-window-v1"
maha_date <- function(x) {
	if (inherits(x, "Date")) return(x)
	if (is.numeric(x)) return(as.Date(x, origin = "1970-01-01"))
	as.Date(as.character(x))
}

maha_webq_window <- function(raw) {
	date_cols <- paste0("p105010_i", 0:4)
	stopifnot(all(c("eid", date_cols) %in% names(raw)), !anyDuplicated(raw$eid))
	dates <- lapply(raw[date_cols], maha_date)
	last <- do.call(pmax, c(dates, list(na.rm = TRUE)))
	last[!is.finite(as.numeric(last))] <- as.Date(NA)
	# A recorded duration/energy without a date makes its temporal ordering unknown.
	undated <- rep(FALSE, nrow(raw))
	for (k in 0:4) {
		present <- rep(FALSE, nrow(raw))
		for (v in paste0(c("p20082_i", "p26002_i"), k)) {
			if (v %in% names(raw)) present <- present | !is.na(raw[[v]])
		}
		undated <- undated | (present & is.na(dates[[k + 1L]]))
	}
	data.frame(eid = as.character(raw$eid), webq_window_end = last,
				webq_undated = undated, webq_n_dates = rowSums(!is.na(as.data.frame(dates))))
}

maha_set_landmark <- function(dat, window, score_cols) {
	i <- match(as.character(dat$eid), window$eid)
	dat$webq_window_end <- window$webq_window_end[i]
	dat$webq_undated <- window$webq_undated[i]
	scored <- rowSums(!is.na(dat[score_cols])) > 0
	invalid <- scored & (is.na(dat$webq_window_end) | is.na(dat$webq_undated) | dat$webq_undated)
	dat$webq_timing_excluded <- invalid
	if (any(invalid)) message("[LANDMARK] Excluding ", sum(invalid), " scored participants with missing/undated WebQ dates from survival analyses.")
	attend <- maha_date(dat$date_attend)
	dat$diet_landmark <- pmax(attend, dat$webq_window_end)
	dat$diet_landmark[invalid] <- as.Date(NA)
	if ("start_date" %in% names(dat)) {
		inconsistent <- !is.na(dat$start_date) & !is.na(dat$diet_landmark) & maha_date(dat$start_date) > dat$diet_landmark
		if (any(inconsistent)) stop("Stored score start_date is later than the raw WebQ window.")
	}
	dat$age.landmark <- as.numeric(dat$age) + as.numeric(dat$diet_landmark - attend) / 365.25
	dat
}

maha_outcome_followup <- function(dat, outcome, end_date) {
	required <- c("diet_landmark", "date_lost", "date_death")
	event_col <- if (outcome == "death") "date_death" else paste0("fod_icd10_", outcome)
	if (!all(c(required, event_col) %in% names(dat))) stop("Missing follow-up input for ", outcome)
	start <- maha_date(dat$diet_landmark)
	event_date <- maha_date(dat[[event_col]])
	dead <- maha_date(dat$date_death)
	censor <- pmin(maha_date(dat$date_lost), dead, as.Date(end_date), na.rm = TRUE)
	prevalent <- !is.na(start) & !is.na(event_date) & event_date <= start
	unknown <- rep(FALSE, nrow(dat))
	# Use only evidence for THIS disease, never a general CVD flag for every outcome.
	if (outcome != "death") for (v in c(paste0("cnt_icd10_", outcome), paste0("icd10Ct_", outcome))) {
		if (v %in% names(dat)) unknown <- unknown | (is.na(event_date) & !is.na(dat[[v]]) & dat[[v]] > 0)
	}
	eligible <- !is.na(start) & !is.na(censor) & censor > start & !prevalent & !unknown
	incident <- eligible & !is.na(event_date) & event_date > start & event_date <= censor
	stop_date <- censor
	stop_date[incident] <- event_date[incident]
	dat[[paste0(outcome, ".t2e")]] <- ifelse(eligible, as.numeric(stop_date - start) / 365.25, NA_real_)
	dat[[paste0(outcome, ".Yt2e")]] <- ifelse(eligible, as.integer(incident), NA_integer_)
	attr(dat, "maha_outcome_audit") <- data.frame(outcome = outcome, N_input = nrow(dat),
		N_eligible = sum(eligible), events = sum(incident), prevalent_through_landmark = sum(prevalent),
		unknown_event_date = sum(unknown), no_positive_followup = sum(is.na(start) | is.na(censor) | censor <= start),
		timing_excluded = sum(dat$webq_timing_excluded %in% TRUE),
		time_origin = "max(enrolment date, last of all recorded WebQ dates)", admin_end = as.character(end_date))
	dat
}

# Raw P values are immutable; numerical floor and common display cap are separate.
maha_phewas_p <- function(res) {
	stopifnot(all(c("p", "snp") %in% names(res)))
	if (any(!is.na(res$p) & (res$p < 0 | res$p > 1))) stop("Invalid PheWAS P value")
	res$p_raw <- res$p
	res$p_underflow <- !is.na(res$p) & res$p == 0
	res$logp <-  - log10(pmax(res$p, .Machine$double.xmin))
	res$logp_display <- pmin(res$logp, 50)
	res$display_capped <- !is.na(res$logp) & res$logp > 50
	res$FDR <- ave(res$p, res$snp, FUN = function(p) p.adjust(p, "BH"))
	res
}

maha_raw_phewas <- function(dat, Xs, covs, phenotype_file, cache_dir, replace = FALSE) {
	# This package's addPhecodeInfo expects its pheinfo dataset on the search path.
	suppressPackageStartupMessages(library(PheWAS))
	dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
	phenos <- readRDS(phenotype_file)
	if (!"id" %in% names(phenos)) stop("PheWAS phenotype input has no id column")
	phenos$id <- as.character(phenos$id)
	dat$eid <- as.character(dat$eid)
	# Full-content fingerprint prevents stale resumption after input or model changes.
	tf <- tempfile(fileext = ".rds"); on.exit(unlink(tf), add = TRUE)
	saveRDS(list(version = "raw-p-v1", dat = dat[, c("eid", Xs, covs)],
					phenotypes = phenos, package = as.character(packageVersion("PheWAS"))), tf, version = 2)
	signature <- unname(tools::md5sum(tf))
	cores <- suppressWarnings(as.integer(Sys.getenv("PHEWAS_CORES", "1")))
	if (is.na(cores) || cores < 1) cores <- 1L
	out <- lapply(Xs, function(x) {
		path <- file.path(cache_dir, paste0(x, ".raw.rds"))
		cached <- if (!replace && file.exists(path)) readRDS(path) else NULL
		if (is.list(cached) && identical(cached$signature, signature)) return(cached$result)
		message("[RAW PHEWAS] Fitting ", x)
		g <- dat[, c("eid", x)]; names(g)[1] <- "id"
		cv <- dat[, c("eid", covs)]; names(cv)[1] <- "id"
		r <- PheWAS::phewas(phenotypes = phenos, genotypes = g, covariates = cv,
						cores = cores, significance.threshold = "bonferroni")
		r <- as.data.frame(PheWAS::addPhecodeInfo(r))
		names(r)[names(r) == "group"] <- "category"
		r <- maha_phewas_p(r)
		saveRDS(list(signature = signature, result = r), path)
		r
	})
	list(res = dplyr::bind_rows(out), signature = signature)
}



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

covs <- c("age.landmark", "sex.f", "center", "tdi", "PC1", "PC2", "bmi.pts", "bp.pts", "nonhdl.pts", "smoke.pts", "pa.pts", "sleep.pts")
vars.risk <- c("age.landmark", "sex", "tdi", "PC1", "PC2")

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

score_0_100_q5 <- function(x) { x <- as.numeric(x); ok <- is.finite(x); if (sum(ok) < 2 || diff(range(x[ok])) == 0) rep(NA_real_, length(x)) else { z <- rep(NA_real_, length(x)); z[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100; as.numeric(as.character(cut(z, c( - Inf, 20, 40, 60, 80, Inf), labels = c(0, 25, 50, 75, 100), right = FALSE))) } }

make_phewas_panel <- function(diet_name, show_x = TRUE) {
    d <- res %>% filter(Diet == diet_name)
    exposure <- names(lab)[unname(lab) == diet_name][1]
    info <- phewas_sample_info(dat1, exposure, vars.basic, phewas_plot_ids, type = "continuous")
    plot_phewas_results(d, title = paste0(diet_name, " PheWAS"), sample_info = info,
        significance_threshold = 10^(-bonf_line), xmap = phewas_xmap,
        colors = phewas_cols, show_x = show_x, label_n = 12, base_size = 9)
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
		facet_wrap( ~ xf, nrow = 1, strip.position = "bottom") +
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
	res$N_model <- nrow(d)
	res$Events_model <- sum(d[[Y_event]])
	res$adjustment <- paste(covs2, collapse = "; ")
	res$model <- "Cox; DASH + MAHA (no interaction)"
	res$prediction <- "Fixed profile: numeric means; categorical modes in model sample"
	res$covariate_profile <- paste(paste0(covs2, "=", vapply(nd[1, covs2, drop = FALSE], as.character, character(1))), collapse = "; ")
	res$time_origin <- "End of all recorded WebQ assessments (or enrolment if later)"
	res$N_complete9 <- sum(d$maha_n_components == 9)
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
		P_equal_coefficients = 2 * pnorm( - abs(diff / se_diff)),
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
				BF01 <= 1 / 3 ~ "BF against point null",
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
	shift <- c(MAHA =  - 0.27, DASH =  - 0.09, MIND = 0.09, MEDI = 0.27)
	lab_layout <- tibble(
		Diet = factor(c("MAHA", "DASH", "MIND", "MEDI"), levels = diet_levels.ab),
		lab_col = c(1, 2, 1, 2),
		lab_y_shift = c(0.16, 0.16, - 0.16, - 0.16)
	)
	d <- d %>%
		mutate(y = as.numeric(Outcome) + unname(shift[as.character(Diet)]), lab = HR_95CI) %>%
		left_join(lab_layout, by = "Diet") %>%
		mutate(y_lab = as.numeric(Outcome) + lab_y_shift, x_lab = if_else(lab_col == 1, 1.095, 1.225))

	ggplot(d, aes(color = Diet)) +
		annotate("rect", xmin = 1 / margin, xmax = margin, ymin =  - Inf, ymax = Inf, alpha = 0.10, fill = "grey70") +
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
	shift <- c(DASH =  - 0.22, MIND = 0, MEDI = 0.22)
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
		annotate("rect", xmin = 1 / margin, xmax = margin, ymin =  - Inf, ymax = Inf, alpha = 0.10, fill = "grey70") +
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
	bind_rows(lapply(seq_len(nrow(nd)), function(i) { ss <- summary(survfit(fit, newdata = nd[i, , drop = FALSE]), times = t0, extend = TRUE); data.frame(Y = Y, dash = as.character(nd$dash[i]), maha = as.character(nd$maha[i]), risk10 = 100 * (1 - ss$surv),
			N_model = nrow(d2), Events_model = sum(d2$event),
			N_group = sum(d2$dash == nd$dash[i] & d2$maha == nd$maha[i]),
			Events_group = sum(d2$event[d2$dash == nd$dash[i] & d2$maha == nd$maha[i]]),
			adjustment = paste(covs, collapse = "; "), model = "Cox; DASH * MAHA interaction",
			prediction = "Fixed profile: numeric means; categorical modes in model sample",
			covariate_profile = paste(paste0(covs, "=", vapply(nd[i, covs, drop = FALSE], as.character, character(1))), collapse = "; "),
			time_origin = "End of all recorded WebQ assessments (or enrolment if later)") }))
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
	g <- expand.grid(rep(list(c( - 1, 1)), p), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
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
			protective_Z =  - beta / se,
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
