# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
# 🚩 CHNS: MAHA replication / validation pipeline
# Single-cohort replication with explicit CHNS inputs and CFCT FOODCODE maps.

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", os_path("DIR_DAT"), os_path("DIR_DAT"))
project_dir <- Sys.getenv("DIET_CHNS", Sys.getenv("CHNS_PROJECT_DIR", file.path(Sys.getenv("DIR_DAT", dir0), "chns")))
raw_dir <- file.path(project_dir, "raw")
taobao_dir <- file.path(project_dir, "CHNS2.0-TaoB", "CHNS原始数据")
id2019_dir <- file.path(taobao_dir, "Master_ID_201908")
foodcode_dir <- file.path(project_dir, "MAHA_foodcode")
food_map_default <- file.path(foodcode_dir, "chns_foodcode_map.tsv")
pexam_public_default <- file.path(raw_dir, "individual", "pexampub12.sav")
nutrient_map_default <- file.path(foodcode_dir, "chns_foodcode_nutrients.tsv")

outdir <- Sys.getenv("MAHA_OUTDIR", unset = file.path(Sys.getenv("DIR_DAT", if (.Platform$OS.type == "windows") os_path("DIR_DAT") else os_path("DIR_DAT")), "analysis", "maha"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

maha_load_packages(tidyverse, haven, survival, broom, writexl, patchwork, cowplot, readxl, scales)
.this_file <- maha_runner
.this_dir <- maha_code
source(file.path(maha_code, "common.R"))

cohort_prefix <- "chns"
maha_auxdir <- Sys.getenv("MAHA_AUXDIR", unset = file.path(outdir, cohort_prefix))
dir.create(maha_auxdir, recursive = TRUE, showWarnings = FALSE)
cohort_file <- function(path) file.path(maha_auxdir, basename(path))
pub_file <- function(path) file.path(outdir, basename(path))
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 360, bg = "white") {
	bounded_ggsave(cohort_file(filename), plot = plot, width = width, height = height, dpi = dpi, bg = bg)
}
save_pub_plot <- function(plot, filename, width = 8, height = 6, dpi = 360, bg = "white") {
	bounded_ggsave(pub_file(filename), plot = plot, width = width, height = height, dpi = dpi, bg = bg)
}
write_xlsx <- function(x, path, ...) writexl::write_xlsx(x, cohort_file(path), ...)
write_pub_xlsx <- function(x, path, ...) writexl::write_xlsx(x, pub_file(path), ...)
step_header <- function(x) cat("\n", strrep("=", 88), "\n", x, "\n", strrep("=", 88), "\n", sep = "")

analysis_spec <- "primary"
if (!identical(analysis_spec, "primary")) stop("CHNS primary model must remain primary in this publication pipeline.", call. = FALSE)

baseline_req <- trimws(Sys.getenv("CHNS_BASELINE_WAVE", unset = "2011"))
followup_req <- trimws(Sys.getenv("CHNS_FOLLOWUP_END", unset = "2015"))
min_followup_years <- suppressWarnings(as.integer(Sys.getenv("CHNS_MIN_FOLLOWUP_YEARS", unset = "4")))
min_age <- suppressWarnings(as.integer(Sys.getenv("CHNS_MIN_AGE", unset = "40")))
min_diet_days <- suppressWarnings(as.integer(Sys.getenv("CHNS_MIN_DIET_DAYS", unset = "2")))
maha_jobs <- suppressWarnings(as.integer(Sys.getenv("MAHA_JOBS", unset = ifelse(.Platform$OS.type == "windows", "1", "4"))))
null_max <- suppressWarnings(as.integer(Sys.getenv("MAHA_NULL_MAX", unset = "0")))
full_map_min_cov <- suppressWarnings(as.numeric(Sys.getenv("CHNS_MIN_MAP_COVERAGE", unset = "0.90")))
require_nutrient_map <- suppressWarnings(as.integer(Sys.getenv("CHNS_REQUIRE_NUTRIENT_MAP", unset = "1")))
require_enhanced_map <- suppressWarnings(as.integer(Sys.getenv("CHNS_REQUIRE_ENHANCED_MAP", unset = "1")))
run_mortality_sweep <- suppressWarnings(as.integer(Sys.getenv("CHNS_RUN_MORTALITY_SWEEP", unset = "1")))
mortality_sweep_child <- identical(Sys.getenv("CHNS_MORTALITY_SWEEP_CHILD", unset = "0"), "1")
mortality_sweep_wave_req <- trimws(Sys.getenv("CHNS_MORTALITY_SWEEP_WAVES", unset = "2004,2006,2009,2011"))
mortality_model_age <- suppressWarnings(as.integer(Sys.getenv("CHNS_MORTALITY_MODEL_AGE", unset = "40")))
mortality_event_target <- suppressWarnings(as.integer(Sys.getenv("CHNS_MORTALITY_EVENT_TARGET", unset = "500")))
if (!require_nutrient_map %in% c(0L, 1L)) require_nutrient_map <- 1L
if (!require_enhanced_map %in% c(0L, 1L)) require_enhanced_map <- 1L
if (!run_mortality_sweep %in% c(0L, 1L)) run_mortality_sweep <- 1L
if (!is.finite(mortality_model_age) || mortality_model_age < 18L) mortality_model_age <- 40L
if (!is.finite(mortality_event_target) || mortality_event_target < 1L) mortality_event_target <- 500L

if (!is.finite(min_followup_years) || min_followup_years < 1) min_followup_years <- 4L
if (!is.finite(min_age) || min_age < 18) min_age <- 40L
if (!is.finite(min_diet_days) || min_diet_days < 2) min_diet_days <- 2L
if (!is.finite(maha_jobs) || maha_jobs < 1) maha_jobs <- 1L
if (!is.finite(null_max) || null_max < 0) null_max <- 0L
if (!is.finite(full_map_min_cov) || full_map_min_cov <= 0 || full_map_min_cov > 1) full_map_min_cov <- .90

# 🚩 Utilities
chns_refresh_scales <- function(dat, variables) {
	# Native plot/group representations are derived from the published raw sums.
	for (v in variables) if (v %in% names(dat)) {
		stem <- sub("\\.sum$", "", v)
		dat[[paste0(stem, ".s100")]] <- score100(dat[[v]])
		dat[[paste0(stem, ".3c")]] <- f3c(dat[[v]])
	}
	dat
}

num <- function(x) {
	if (inherits(x, "haven_labelled") || inherits(x, "labelled")) x <- haven::zap_labels(x)
	suppressWarnings(as.numeric(x))
}
recode_chns_school_years <- function(x) {
	z <- num(x)
	# Estimated years: 6 primary + 3 lower secondary + 3 upper secondary;
	# technical school follows lower secondary; college 6+ years is capped at 18.
	dplyr::case_when(
		z == 0 ~ 0,
		z %in% 11:16 ~ z - 10,
		z %in% 21:26 ~ z - 14,
		z %in% 27:29 ~ z - 17,
		z %in% 31:36 ~ z - 18,
		TRUE ~ NA_real_
	)
}
cache_identical_inputs <- function(fn) {
	# Run-local only: reuse intervals only when all input data and the definition
	# are identical. Downstream models still fit independently on these intervals.
	entries <- list()
	function(d, mortality_definition) {
		for(entry in entries) {
			if(identical(entry$definition, mortality_definition) && identical(entry$data, d))
				return(entry$result)
		}
		message("[TIME UPDATED] Building intervals: ", mortality_definition, "; input rows=", nrow(d))
		result <- fn(d, mortality_definition)
		entries[[length(entries) + 1L]] <<- list(data = d, definition = mortality_definition, result = result)
		message("[TIME UPDATED] Cached intervals: ", nrow(result), " rows")
		result
	}
}
zstd <- function(x) {
	x <- num(x); s <- sd(x, na.rm = TRUE)
	if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
	as.numeric(scale(x))
}
first_num <- function(x) {
	x <- num(x); x <- x[is.finite(x)]
	if (length(x)) x[1] else NA_real_
}
qscore <- function(x, reverse = FALSE) {
	x <- num(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
	if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(out)
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- as.numeric(as.character(cut(r, c( - Inf, .2, .4, .6, .8, Inf),
						labels = c(0, 2.5, 5, 7.5, 10), right = TRUE)))
	if (reverse) out[ok] <- 10 - out[ok]
	out
}
qscore5 <- function(x, reverse = FALSE) {
	x <- num(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
	if (sum(ok) < 5 || length(unique(x[ok])) < 2) return(out)
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- as.numeric(as.character(cut(r, c( - Inf, .2, .4, .6, .8, Inf),
						labels = 1:5, right = TRUE)))
	if (reverse) out[ok] <- 6 - out[ok]
	out
}
qscore5_grouped <- function(x, group, reverse = FALSE) {
	x <- num(x); g <- as.character(group); out <- rep(NA_real_, length(x))
	for (lev in unique(g[!is.na(g)])) {
		ii <- which(g == lev & is.finite(x))
		if (length(ii)) out[ii] <- qscore5(x[ii], reverse = reverse)
	}
	out
}
energy_residual_grouped <- function(x, kcal, group) {
	x <- num(x); e <- num(kcal); g <- as.character(group); out <- rep(NA_real_, length(x))
	for (lev in unique(g[!is.na(g)])) {
		ii <- which(g == lev & is.finite(x) & is.finite(e))
		if (length(ii) < 20L || stats::sd(e[ii]) == 0) next
		fit <- stats::lm(x[ii] ~ e[ii])
		out[ii] <- stats::residuals(fit) + mean(x[ii], na.rm = TRUE)
	}
	out
}
energy_density <- function(x, kcal, multiplier = 1000) {
	x <- num(x); e <- num(kcal)
	ifelse(is.finite(x) & is.finite(e) & e > 0, x / e * multiplier, NA_real_)
}
score100 <- function(x) {
	x <- num(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x)
	if (sum(ok) < 2 || diff(range(x[ok])) == 0) return(out)
	out[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100
	out
}
f3c <- function(x) {
	x <- num(x); out <- rep(NA_character_, length(x)); ok <- is.finite(x)
	if (sum(ok) < 3 || length(unique(x[ok])) < 2) return(factor(out, levels = c("low", "middle", "high")))
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- as.character(cut(r, c( - Inf, 1 / 3, 2 / 3, Inf), labels = c("low", "middle", "high")))
	factor(out, levels = c("low", "middle", "high"))
}
rowmean_min <- function(d, vars, min_n = 1L) {
	vars <- vars[vars %in% names(d)]
	if (!length(vars)) return(rep(NA_real_, nrow(d)))
	m <- as.matrix(d[, vars, drop = FALSE])
	nok <- rowSums(is.finite(m))
	o <- rowMeans(m, na.rm = TRUE)
	o[nok < min_n] <- NA_real_
	o
}
rowmean_vecs <- function(..., min_n = 1L) {
	m <- cbind(...)
	nok <- rowSums(is.finite(m))
	o <- rowMeans(m, na.rm = TRUE)
	o[nok < min_n] <- NA_real_
	o
}
sum_known <- function(x) {
	x <- num(x)
	if (!any(is.finite(x))) return(NA_real_)
	sum(x, na.rm = TRUE)
}
weighted_sum_known <- function(g, flag) {
	g <- num(g); flag <- num(flag)
	ok <- is.finite(g) & is.finite(flag)
	if (!any(ok)) return(NA_real_)
	sum(g[ok] * flag[ok], na.rm = TRUE)
}
safe_cor <- function(x, y) {
	x <- num(x); y <- num(y)
	ok <- is.finite(x) & is.finite(y)
	if (sum(ok) < 2L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(NA_real_)
	suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}
yes1 <- function(x) {
	z <- num(x)
	ifelse(z == 1, 1L, ifelse(z %in% c(0, 2), 0L, NA_integer_))
}
pref15 <- function(x) {
	# Prefer the embedded CHNS value labels rather than assuming numeric orientation.
	# This protects against releases in which the raw 1-5 coding is reversed.
	labs <- attr(x, "labels")
	z <- num(x); out <- rep(NA_real_, length(z))
	if(!is.null(labs) && length(labs)) {
		n <- tolower(names(labs)); val <- suppressWarnings(as.numeric(labs))
		rec <- rep(NA_real_, length(val))
		rec[grepl("dislike.*very|very.*dislike", n)] <- 1
		rec[grepl("dislike.*some|some.*dislike", n)] <- 2
		rec[grepl("neutral", n)] <- 3
		rec[grepl("like.*some|some.*like", n)] <- 4
		rec[grepl("like.*very|very.*like", n)] <- 5
		for(i in which(is.finite(rec))) out[z == val[i]] <- rec[i]
		# Explicit 'not eat' / non-consumption categories remain missing.
		if(any(is.finite(out))) return(out)
	}
	# Fallback consistent with the published CHNS food-preference analysis.
	out[z %in% 1:5] <- z[z %in% 1:5]
	out
}
moderate_alcohol_score <- function(g, female) {
	g <- num(g); female <- as.integer(female)
	out <- rep(NA_real_, length(g)); ok <- is.finite(g) & female %in% c(0L, 1L); out[ok] <- 0
	out[ok & female == 1L & g >= 5 & g <= 15] <- 10
	out[ok & female == 0L & g >= 10 & g <= 25] <- 10
	out
}
good_covariates <- function(d, covs, min_complete = .50) {
	covs <- intersect(covs, names(d))
	covs[vapply(covs, function(v) {
		x <- d[[v]]
		mean(!is.na(x)) >= min_complete && length(unique(x[!is.na(x)])) > 1
	}, logical(1))]
}
parse_year_request <- function(x, label) {
	if (!nzchar(x) || tolower(x) == "auto") return(NA_integer_)
	z <- suppressWarnings(as.integer(x))
	if (!is.finite(z)) stop(label, " must be an integer year or 'auto': ", x, call. = FALSE)
	z
}
valid_years <- function(x) {
	z <- num(x)
	sort(unique(z[z >= 1989 & z <= 2035]))
}
pick <- function(d, cands, required = TRUE) {
	cands <- toupper(cands)
	nm <- cands[cands %in% names(d)]
	if (length(nm)) return(nm[1])
	if (required) stop("Required variable missing: ", paste(cands, collapse = " | "), call. = FALSE)
	NA_character_
}
colnum <- function(d, nm) {
	if (is.na(nm) || !nm %in% names(d)) rep(NA_real_, nrow(d)) else num(d[[nm]])
}
read_any <- function(f, col_select = NULL) {
	ext <- tolower(tools::file_ext(f))
	if (ext == "sas7bdat") {
		# CHNS SAS names use mixed case (for example, IDind), while downstream
		# processing standardizes all names to upper case. Select requested fields
		# case-insensitively so the pre-read projection does not silently drop them.
		d <- if (is.null(col_select)) haven::read_sas(f) else haven::read_sas(
			f, col_select = matches(paste0("^(?:", paste(col_select, collapse = "|"), ")$"), ignore.case = TRUE)
		)
	} else if (ext == "sav") {
		d <- haven::read_sav(f)
		if (!is.null(col_select)) d <- d %>% dplyr::select(any_of(col_select))
	} else if (ext == "xpt") {
		d <- haven::read_xpt(f)
		if (!is.null(col_select)) d <- d %>% dplyr::select(any_of(col_select))
	} else if (ext == "csv") {
		d <- readr::read_csv(f, show_col_types = FALSE)
	} else if (ext %in% c("tsv", "txt")) {
		d <- readr::read_tsv(f, show_col_types = FALSE)
	} else if (ext %in% c("xlsx", "xls")) {
		d <- readxl::read_excel(f)
	} else stop("Unsupported file type: ", f)
	d <- as_tibble(d); names(d) <- toupper(names(d)); d
}
print_source <- function(nm, x, required = TRUE) {
	ok <- is.character(x[["path"]]) && length(x[["path"]]) == 1 && !is.na(x[["path"]]) && file.exists(x[["path"]])
	cat(sprintf("  %-12s %-7s %-14s : %s%s\n", nm, if(ok)"OK" else if(required)"MISSING" else "OPTIONAL",
					x[["source"]], ifelse(is.na(x[["path"]]), "<not found>", x[["path"]]),
					if(ok)"" else if(required)" 文件不存在" else ""))
	invisible(ok)
}
linear_adequacy <- function(x, target, max_score) {
	x <- num(x); out <- pmin(pmax(x / target, 0), 1) * max_score
	out[!is.finite(x)] <- NA_real_; out
}
linear_types_score <- function(n) {
	n <- num(n)
	out <- ifelse(!is.finite(n), NA_real_, ifelse(n <= 5, 0, ifelse(n >= 12, 10, (n - 5) / 7 * 10)))
	out
}




fit_bin_score <- function(d, outcome, score, covs = c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "bmi", "kcal")) {
	dd <- d; dd$score <- zstd(score)
	covs <- good_covariates(dd, covs)
	use <- c(outcome, "score", covs); dd <- dd %>% dplyr::select(all_of(use)) %>% drop_na()
	if(nrow(dd) < 100 || sum(dd[[outcome]] == 1) < 10 || length(unique(dd$score)) < 2)
		return(tibble(N = nrow(dd), Events = sum(dd[[outcome]] == 1), beta = NA_real_, se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, Covariates = paste(covs, collapse = " + ")))
	fm <- as.formula(paste(outcome, "~ score", if(length(covs))paste0(" + ", paste(covs, collapse = " + "))else""))
	fit <- tryCatch(glm(fm, data = dd, family = binomial()), error = function(e)NULL)
	if(is.null(fit)) return(tibble(N = nrow(dd), Events = sum(dd[[outcome]] == 1), beta = NA_real_, se = NA_real_, OR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, Covariates = paste(covs, collapse = " + ")))
	tt <- broom::tidy(fit) %>% filter(term == "score")
	tibble(N = nrow(dd), Events = sum(dd[[outcome]] == 1), beta = tt$estimate[1], se = tt$std.error[1],
				 OR = exp(tt$estimate[1]), CI_low = exp(tt$estimate[1] - 1.96 * tt$std.error[1]),
				 CI_high = exp(tt$estimate[1] + 1.96 * tt$std.error[1]), P = tt$p.value[1], Covariates = paste(covs, collapse = " + "))
}
fit_linear_pref <- function(d, outcome, score, covs = c("age", "female", "province", "urban", "log_income", "education_years")) {
	dd <- d; dd$score <- zstd(score)
	covs <- good_covariates(dd, covs); dd <- dd %>% dplyr::select(all_of(c(outcome, "score", covs))) %>% drop_na()
	if(nrow(dd) < 100) return(tibble(N = nrow(dd), beta = NA_real_, se = NA_real_, P = NA_real_))
	fm <- as.formula(paste(outcome, "~ score", if(length(covs))paste0(" + ", paste(covs, collapse = " + "))else""))
	fit <- lm(fm, data = dd); tt <- broom::tidy(fit) %>% filter(term == "score")
	tibble(N = nrow(dd), beta = tt$estimate[1], se = tt$std.error[1], P = tt$p.value[1])
}
fit_linear_score <- function(d, outcome, score, covs = c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "bmi", "kcal")) {
	dd <- d; dd$score <- zstd(score)
	covs <- good_covariates(dd, covs); dd <- dd %>% dplyr::select(all_of(c(outcome, "score", covs))) %>% drop_na()
	if(nrow(dd) < 100 || length(unique(dd$score)) < 2) return(tibble(N = nrow(dd), beta = NA_real_, se = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, Covariates = paste(covs, collapse = " + ")))
	fm <- as.formula(paste(outcome, "~ score", if(length(covs))paste0(" + ", paste(covs, collapse = " + "))else""))
	fit <- tryCatch(lm(fm, data = dd), error = function(e)NULL)
	if(is.null(fit)) return(tibble(N = nrow(dd), beta = NA_real_, se = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, Covariates = paste(covs, collapse = " + ")))
	tt <- broom::tidy(fit) %>% filter(term == "score")
	tibble(N = nrow(dd), beta = tt$estimate[1], se = tt$std.error[1], CI_low = tt$estimate[1] - 1.96 * tt$std.error[1], CI_high = tt$estimate[1] + 1.96 * tt$std.error[1], P = tt$p.value[1], Covariates = paste(covs, collapse = " + "))
}
# Publication helper: direct MAHA-versus-established coefficient comparison in a
# mutually adjusted logistic model. The reported OR ratio is exp(beta_MAHA-beta_trad):
# values >1 indicate that MAHA is less protective than the established comparator.
fit_head2head_bin <- function(d, trad_var, trad_lab, maha_var = "diet.maha.sum",
			outcome = "death",
			covs = c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "bmi", "kcal")) {
	dd <- d
	dd$maha_z <- zstd(dd[[maha_var]])
	dd$trad_z <- zstd(dd[[trad_var]])
	covs <- good_covariates(dd, covs)
	dd <- dd %>% dplyr::select(all_of(c(outcome, "maha_z", "trad_z", covs))) %>% drop_na()
	if (nrow(dd) < 100 || sum(dd[[outcome]] == 1) < 10) return(tibble())
	fm <- as.formula(paste(outcome, "~ maha_z + trad_z", if(length(covs)) paste0(" + ", paste(covs, collapse = " + ")) else ""))
	fit <- tryCatch(glm(fm, data = dd, family = binomial()), error = function(e) NULL)
	if (is.null(fit)) return(tibble())
	cf <- coef(fit); vc <- vcov(fit)
	if (!all(c("maha_z", "trad_z") %in% names(cf))) return(tibble())
	b1 <- unname(cf[["maha_z"]]); b2 <- unname(cf[["trad_z"]])
	dff <- b1 - b2
	se_diff <- sqrt(vc["maha_z", "maha_z"] + vc["trad_z", "trad_z"] - 2 * vc["maha_z", "trad_z"])
	tibble(
		Outcome = "All-cause mortality", Comparison = paste0("MAHA vs ", trad_lab), Comp2 = trad_lab,
		N_total = nrow(dd), N_event = sum(dd[[outcome]] == 1),
		OR_maha_mutual = exp(b1), OR_trad_mutual = exp(b2),
		logOR_difference_maha_minus_trad = dff, OR_ratio_maha_vs_trad = exp(dff),
		SE_difference = se_diff,
		LCI_ratio = exp(dff - 1.96 * se_diff), UCI_ratio = exp(dff + 1.96 * se_diff),
		P_equal_coefficients = 2 * pnorm( - abs(dff / se_diff))
	)
}

make_joint_logistic_risk <- function(d, facet_var, line_var, facet_lab, line_lab,
			outcome = "death",
			covs = c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "bmi", "kcal")) {
	dd <- d
	dd$facet <- factor(dd[[facet_var]], levels = c("low", "middle", "high"))
	dd$line <- factor(dd[[line_var]], levels = c("low", "middle", "high"))
	covs <- good_covariates(dd, covs)
	keep <- c(outcome, "facet", "line", covs)
	dd <- dd %>% dplyr::select(all_of(keep)) %>% drop_na()
	if (nrow(dd) < 100 || sum(dd[[outcome]] == 1) < 10) return(tibble())
	fm <- as.formula(paste(outcome, "~ facet + line", if(length(covs))paste0(" + ", paste(covs, collapse = " + "))else""))
	fit <- glm(fm, data = dd, family = binomial())
	nd <- expand.grid(facet = factor(c("low", "middle", "high"), levels = c("low", "middle", "high")),
				line = factor(c("low", "middle", "high"), levels = c("low", "middle", "high")))
	for(v in covs) {
		x <- dd[[v]]
		if(is.numeric(x)) nd[[v]] <- mean(x, na.rm = TRUE)
		else if(is.factor(x)) nd[[v]] <- factor(names(which.max(table(x))), levels = levels(x))
		else nd[[v]] <- names(which.max(table(x)))
	}
	pr <- predict(fit, newdata = nd, type = "link", se.fit = TRUE)
	nd$risk <- plogis(pr$fit); nd$risk_lower <- plogis(pr$fit - 1.96 * pr$se.fit); nd$risk_upper <- plogis(pr$fit + 1.96 * pr$se.fit)
	cells <- dd %>% count(facet, line, name = "N")
	nd <- left_join(nd, cells, by = c("facet", "line"))
	nd %>% transmute(
		Outcome = "All-cause mortality", Followup_years = followup_end - baseline_wave,
		Facet_var = facet_lab, Facet_group = as.character(facet), Line_var = line_lab, Line_group = as.character(line),
		N, risk, risk_lower, risk_upper,
		risk_pct = 100 * risk, risk_lower_pct = 100 * risk_lower, risk_upper_pct = 100 * risk_upper
	)
}

make_precision_table <- function(mort_tbl, margin_hr = 1.05, prior_95_hr = 1.10) {
	mort_tbl %>% filter(is.finite(beta), is.finite(se), se > 0) %>% rowwise() %>% do({
		r <- .
		tost <- ap_tost_one(r$beta, r$se, margin_hr = margin_hr)
		pw <- ap_power_one(r$se, margin_hr = margin_hr)
		bf <- ap_bf01_one(r$beta, r$se, prior_95_hr = prior_95_hr)
		tibble(
			Diet = r$Diet, N_total = r$N, N_event = r$Events, beta = r$beta, se = r$se, OR = r$OR, LCI95 = r$CI_low, UCI95 = r$CI_high, P = r$P,
			OR_95CI = ap_ci_lab(r$OR, r$CI_low, r$CI_high),
			OR_90CI = ap_ci_lab(r$OR, tost$CI90_low_HR, tost$CI90_high_HR),
			margin_hr = margin_hr, CI90_low_OR = tost$CI90_low_HR, CI90_high_OR = tost$CI90_high_HR,
			TOST_p = tost$TOST_p, Equivalence_5pct = ifelse(isTRUE(tost$equivalent), "Equivalent within OR 0.95-1.05", "Not equivalent/inconclusive"),
			prior_95_hr = prior_95_hr, BF01 = bf, BF_summary = case_when(bf >= 3 ~ "BF favors near-zero", bf <= 1 / 3 ~ "BF favors non-zero", TRUE ~ "BF inconclusive"),
			MDE_OR_lower = pw$MDE_HR_lower, MDE_OR_upper = pw$MDE_HR_upper,
			MDE80_OR_label = sprintf("%.2f-%.2f", pw$MDE_HR_lower, pw$MDE_HR_upper),
			Powered_5pct = ifelse(isTRUE(pw$powered_for_margin), "Powered for OR 1.05", "Not powered for OR 1.05")
		)
	}) %>% ungroup()
}

par_lapply <- function(X, FUN) {
	if(.Platform$OS.type != "windows" && maha_jobs > 1) parallel::mclapply(X, FUN, mc.cores = maha_jobs, mc.preschedule = TRUE) else lapply(X, FUN)
}


# 🚩 Primary CHNS mortality model helpers
primary_covariates <- c("age", "female", "province", "urban", "log_income", "education_years", "current_smoke", "bmi", "kcal")

add_primary_survival_fields <- function(d, mortality_definition = "end_wave") {
	d <- as_tibble(d)
	dod_y <- if ("DOD_Y" %in% names(d)) num(d$DOD_Y) else rep(NA_real_, nrow(d))
	dod_rpt <- if ("DOD_RPT" %in% names(d)) num(d$DOD_RPT) else rep(NA_real_, nrow(d))
	death_time <- ifelse(
		is.finite(dod_y) & dod_y > baseline_wave & dod_y <= followup_end,
		dod_y - baseline_wave,
		ifelse(is.finite(dod_rpt) & dod_rpt > baseline_wave & dod_rpt <= followup_end,
					 dod_rpt - baseline_wave, NA_real_)
	)
	death_reported <- if ("death" %in% names(d)) num(d$death) == 1 else rep(FALSE, nrow(d))
	has_postbaseline_death_date <- (is.finite(dod_y) & dod_y > baseline_wave) |
		(is.finite(dod_rpt) & dod_rpt > baseline_wave)
	invalid_reported_death <- !is.na(death_reported) & death_reported & !has_postbaseline_death_date
	death_event <- ifelse(invalid_reported_death, NA_integer_, as.integer(is.finite(death_time)))
	if (identical(mortality_definition, "end_wave")) {
		eligible <- (if ("followup_known" %in% names(d)) d$followup_known == 1 else is.finite(d$death)) & !invalid_reported_death
		event <- ifelse(eligible, death_event, NA_integer_)
		censor_time <- rep(followup_end - baseline_wave, nrow(d))
	} else if (identical(mortality_definition, "last_contact")) {
		last_contact <- if ("last_contact_year" %in% names(d)) num(d$last_contact_year) else rep(NA_real_, nrow(d))
		eligible <- !invalid_reported_death & (death_event == 1 | (is.finite(last_contact) & last_contact > baseline_wave))
		event <- ifelse(eligible, death_event, NA_integer_)
		censor_time <- pmin(last_contact, followup_end, na.rm = TRUE) - baseline_wave
		censor_time[!is.finite(last_contact)] <- NA_real_
	} else stop("Unknown mortality definition: ", mortality_definition, call. = FALSE)
	d$surv_event <- event
	d$surv_time <- ifelse(event == 1, pmax(death_time, 1 / 365), censor_time)
	d$surv_time[!is.finite(d$surv_time) | d$surv_time <= 0] <- NA_real_
	d$mortality_definition <- mortality_definition
	d$invalid_reported_death_date <- as.integer(invalid_reported_death)
	d
}

q4safe <- function(x) {
	x <- num(x)
	out <- rep(NA_character_, length(x))
	ok <- is.finite(x)
	if (sum(ok) < 8 || length(unique(x[ok])) < 4) return(factor(out, levels = c("Q1", "Q2", "Q3", "Q4")))
	r <- dplyr::percent_rank(x[ok])
	out[ok] <- as.character(cut(r, c( - Inf, .25, .50, .75, Inf), labels = c("Q1", "Q2", "Q3", "Q4"), right = TRUE))
	factor(out, levels = c("Q1", "Q2", "Q3", "Q4"))
}

primary_filter <- function(d, require_three_days = FALSE, exclude_first_year = FALSE, exclude_baseline_disease = FALSE,
			mortality_definition = "end_wave") {
	dd <- add_primary_survival_fields(d, mortality_definition = mortality_definition) %>% mutate(.row_id = seq_len(nrow(d)))
	if (require_three_days) dd <- dd %>% filter(n_diet_days >= 3)
	if (exclude_first_year) dd <- dd %>% filter(surv_event == 0 | surv_time > 1)
	if (exclude_baseline_disease) {
		disease_vars <- intersect(c("hypertension", "t2dm_self", "mi_self", "stroke_self", "cancer_self"), names(dd))
		if (length(disease_vars)) {
			known <- rowSums(is.finite(as.matrix(dd[, disease_vars, drop = FALSE])))
			cases <- rowSums(as.matrix(dd[, disease_vars, drop = FALSE]) == 1, na.rm = TRUE)
			dd$baseline_disease <- ifelse(known == 0, NA_integer_, as.integer(cases > 0))
			dd <- dd %>% filter(baseline_disease == 0)
		}
	}
	dd
}

fit_primary_cox <- function(d, score, covs = primary_covariates, require_three_days = FALSE,
			no_bmi = FALSE, exclude_first_year = FALSE, exclude_baseline_disease = FALSE,
			mortality_definition = "end_wave") {
	dd <- primary_filter(d, require_three_days, exclude_first_year, exclude_baseline_disease,
				mortality_definition = mortality_definition)
	score_map <- tibble(.row_id = seq_len(nrow(d)), score_raw = num(score))
	dd <- dd %>% left_join(score_map, by = ".row_id") %>% mutate(score_z = zstd(score_raw))
	covs_requested <- unique(as.character(covs))
	if(no_bmi) covs_requested <- setdiff(covs_requested, "bmi")
	covs_use <- intersect(covs_requested, names(dd))
	covs_use <- good_covariates(dd, covs_use)
	covs_dropped <- setdiff(covs_requested, covs_use)
	dd <- dd %>% select(all_of(c("surv_time", "surv_event", "score_z", covs_use))) %>% drop_na()
	empty <- function() tibble(N = nrow(dd), Events = sum(dd$surv_event == 1, na.rm = TRUE), Person_years = sum(dd$surv_time, na.rm = TRUE),
				beta = NA_real_, se = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
				Covariates = paste(covs_use, collapse = " + "),
				Requested_covariates = paste(covs_requested, collapse = " + "),
				Dropped_covariates = paste(covs_dropped, collapse = " + "))
	if(nrow(dd) < 100 || sum(dd$surv_event == 1, na.rm = TRUE) < 10 || length(unique(dd$score_z)) < 2) return(empty())
	fm <- as.formula(paste0("Surv(surv_time,surv_event) ~ score_z",
					if(length(covs_use)) paste0(" + ", paste(covs_use, collapse = " + ")) else ""))
	fit <- tryCatch(coxph(fm, data = dd, ties = "efron"), error = function(e)NULL)
	if(is.null(fit)) return(empty())
	tt <- broom::tidy(fit) %>% filter(term == "score_z")
	if(!nrow(tt)) return(empty())
	tibble(N = nrow(dd), Events = sum(dd$surv_event == 1, na.rm = TRUE), Person_years = sum(dd$surv_time, na.rm = TRUE),
				 beta = tt$estimate[1], se = tt$std.error[1], HR = exp(tt$estimate[1]),
				 CI_low = exp(tt$estimate[1] - 1.96 * tt$std.error[1]), CI_high = exp(tt$estimate[1] + 1.96 * tt$std.error[1]),
				 P = tt$p.value[1], Covariates = paste(covs_use, collapse = " + "),
				 Requested_covariates = paste(covs_requested, collapse = " + "),
				 Dropped_covariates = paste(covs_dropped, collapse = " + "))
}

fit_primary_quartile <- function(d, score, covs = primary_covariates) {
	dd <- primary_filter(d)
	qmap <- tibble(.row_id = seq_len(nrow(d)), q4 = q4safe(score))
	dd <- dd %>% left_join(qmap, by = ".row_id")
	covs_use <- good_covariates(dd, intersect(covs, names(dd)))
	dd <- dd %>% select(all_of(c("surv_time", "surv_event", "q4", covs_use))) %>% drop_na()
	if(nrow(dd) < 100 || sum(dd$surv_event == 1) < 10 || nlevels(droplevels(dd$q4)) < 4)
		return(tibble(N = nrow(dd), Events = sum(dd$surv_event == 1), Q4_vs_Q1 = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_))
	fm <- as.formula(paste0("Surv(surv_time,surv_event) ~ q4",
					if(length(covs_use)) paste0(" + ", paste(covs_use, collapse = " + ")) else ""))
	fit <- tryCatch(coxph(fm, data = dd, ties = "efron"), error = function(e)NULL)
	if(is.null(fit)) return(tibble(N = nrow(dd), Events = sum(dd$surv_event == 1), Q4_vs_Q1 = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_))
	tt <- broom::tidy(fit) %>% filter(term == "q4Q4")
	if(!nrow(tt)) return(tibble(N = nrow(dd), Events = sum(dd$surv_event == 1), Q4_vs_Q1 = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_))
	tibble(N = nrow(dd), Events = sum(dd$surv_event == 1), Q4_vs_Q1 = exp(tt$estimate[1]),
				 CI_low = exp(tt$estimate[1] - 1.96 * tt$std.error[1]), CI_high = exp(tt$estimate[1] + 1.96 * tt$std.error[1]), P = tt$p.value[1])
}

fit_primary_head2head <- function(d, trad_var, trad_lab, maha_var = "diet.maha.sum", covs = primary_covariates) {
	dd <- primary_filter(d)
	dd$maha_z <- zstd(dd[[maha_var]])
	dd$trad_z <- zstd(dd[[trad_var]])
	covs_use <- good_covariates(dd, intersect(covs, names(dd)))
	dd <- dd %>% select(all_of(c("surv_time", "surv_event", "maha_z", "trad_z", covs_use))) %>% drop_na()
	if(nrow(dd) < 100 || sum(dd$surv_event == 1) < 10) return(tibble())
	fm <- as.formula(paste0("Surv(surv_time,surv_event) ~ maha_z + trad_z",
					if(length(covs_use)) paste0(" + ", paste(covs_use, collapse = " + ")) else ""))
	fit <- tryCatch(coxph(fm, data = dd, ties = "efron"), error = function(e)NULL)
	if(is.null(fit)) return(tibble())
	cf <- coef(fit); vc <- vcov(fit)
	if(!all(c("maha_z", "trad_z") %in% names(cf))) return(tibble())
	b1 <- unname(cf[["maha_z"]]); b2 <- unname(cf[["trad_z"]]); dff <- b1 - b2
	se_diff <- sqrt(vc["maha_z", "maha_z"] + vc["trad_z", "trad_z"] - 2 * vc["maha_z", "trad_z"])
	tibble(Outcome = "All-cause mortality", Comparison = paste0("MAHA vs ", trad_lab), Comp2 = trad_lab,
				 N_total = nrow(dd), N_event = sum(dd$surv_event == 1), HR_maha_mutual = exp(b1), HR_trad_mutual = exp(b2),
				 logHR_difference_maha_minus_trad = dff, HR_ratio_maha_vs_trad = exp(dff), Ratio = exp(dff),
				 SE_difference = se_diff, LCI_ratio = exp(dff - 1.96 * se_diff), UCI_ratio = exp(dff + 1.96 * se_diff),
				 P_equal_coefficients = 2 * pnorm( - abs(dff / se_diff)), Effect_measure = "HR")
}

make_precision_table_hr <- function(mort_tbl, margin_hr = 1.05, prior_95_hr = 1.10) {
	mort_tbl %>% filter(is.finite(beta), is.finite(se), se > 0) %>% rowwise() %>% do({
		r <- .
		tost <- ap_tost_one(r$beta, r$se, margin_hr = margin_hr)
		pw <- ap_power_one(r$se, margin_hr = margin_hr)
		bf <- ap_bf01_one(r$beta, r$se, prior_95_hr = prior_95_hr)
		tibble(Diet = r$Diet, N_total = r$N, N_event = r$Events, beta = r$beta, se = r$se, Effect = r$HR, HR = r$HR,
					 LCI95 = r$CI_low, UCI95 = r$CI_high, P = r$P, Effect_measure = "HR",
					 Effect_95CI = ap_ci_lab(r$HR, r$CI_low, r$CI_high),
					 Effect_90CI = ap_ci_lab(r$HR, tost$CI90_low_HR, tost$CI90_high_HR),
					 margin_effect = margin_hr, CI90_low = tost$CI90_low_HR, CI90_high = tost$CI90_high_HR,
					 TOST_p = tost$TOST_p, Equivalence_5pct = ifelse(isTRUE(tost$equivalent), "Equivalent within HR 0.95-1.05", "Not equivalent/inconclusive"),
					 prior_95_effect = prior_95_hr, BF01 = bf,
					 BF_summary = case_when(bf >= 3 ~ "BF favors near-zero", bf <= 1 / 3 ~ "BF favors non-zero", TRUE ~ "BF inconclusive"),
					 MDE_lower = pw$MDE_HR_lower, MDE_upper = pw$MDE_HR_upper,
					 MDE80_label = sprintf("%.2f-%.2f", pw$MDE_HR_lower, pw$MDE_HR_upper),
					 Powered_5pct = ifelse(isTRUE(pw$powered_for_margin), "Powered for HR 1.05", "Not powered for HR 1.05"))
	}) %>% ungroup()
}
