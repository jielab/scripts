# C1 associations, inherited-score helpers and the dedicated PGS analysis.
# LE8_C1_MODE=helpers loads PGS functions only; pgs/correlate execute one analysis.
.c1_mode <- Sys.getenv("LE8_C1_MODE", unset = if (identical(Sys.getenv("LE8_JOB"), "pgs_focus")) "pgs" else if (sys.nframe() == 0L || identical(Sys.getenv("LE8_JOB"), "c1_correlate")) "correlate" else "helpers")
if (!.c1_mode %in% c("helpers", "pgs", "correlate")) stop("Unknown LE8_C1_MODE: ", .c1_mode)

# C1 inherited-score scan utilities.
#
# A COJO PGS is fixed at conception, but it is not the biomarker concentration
# "measured at birth".  These functions therefore keep PGS and adult measured
# omics in parallel columns and never substitute one for the other.

C1_PGS_SCAN_VERSION <- "2026-09-21.full-scan-provenance"

# c1 pgs helpers
# c1.correlate.R
# Statistical primitives for measured / biomarker-PGS comparisons.
# No individual records are exported. G is captured genetic prediction; R is
# the remainder, not an environmental, lifestyle, or disease-consequence fraction.
PGS_FOCUS_VERSION <- "2026-09-21.matched-locus-v1"
pgs_num <- function(key, default) {
	x <- suppressWarnings(as.numeric(Sys.getenv(key, as.character(default)))) ; if (length(x) != 1L || is.na(x)) stop("Invalid ", key) ; x
}
pgs_csv <- function(key, default = "") unique(Filter(nzchar, trimws(strsplit(Sys.getenv(key, default), ",", fixed = TRUE)[[1]])))
pgs_bind <- function(x) as.data.frame(data.table::rbindlist(x, fill = TRUE, use.names = TRUE))
pgs_hash <- function(x) {
	f <- tempfile() ; on.exit(unlink(f)) ; saveRDS(x, f, compress = FALSE, version = 2) ; unname(tools::md5sum(f))
}
pgs_stamp <- function(files) {
	files <- sort(unique(files[!is.na(files) & nzchar(files)])) ; data.frame(file = files, size = file.info(files)$size, mtime = as.numeric(file.info(files)$mtime))
}
pgs_ids <- function(d, label) {
	if (!"eid" %in% names(d)) stop(label, " lacks eid")
	d$eid <- as.character(d$eid)
	if (anyNA(d$eid) || any(!nzchar(d$eid)) || anyDuplicated(d$eid)) stop(label, " has missing/duplicate eid")
	d
}
pgs_complete <- function(d, cols) {
	if (length(setdiff(cols, names(d)))) return(rep(FALSE, nrow(d)))
	ok <- complete.cases(d[, cols, drop = FALSE])
	for (v in cols) if (is.numeric(d[[v]])) ok <- ok & is.finite(d[[v]])
	ok
}
pgs_folds <- function(group, k = 5L, seed = 2026L) {
	group <- as.character(group) ; if (anyNA(group) || any(!nzchar(group))) stop("Missing PGS fold group")
	u <- sort(unique(group)) ; if (length(u) < k) stop("Too few independent fold groups")
	set.seed(seed) ; f <- sample(rep(seq_len(k), length.out = length(u))) ; f[match(group, u)]
}
pgs_classify <- function(bm, bg, qm, qg, alpha = .05) {
	out <- rep("Neither supported", length(bm)) ; valid <- is.finite(bm) & is.finite(bg) & is.finite(qm) & is.finite(qg)
	m <- is.finite(qm) & qm < alpha ; g <- is.finite(qg) & qg < alpha
	out[m & !g] <- "Measured only" ; out[g & !m] <- "PGS only"
	out[which(valid & m & g & bm * bg >= 0)] <- "Both supported: concordant"
	out[which(valid & m & g & bm * bg < 0)] <- "Both supported: opposite"
	out[!valid] <- "Unavailable comparison" ; out
}
pgs_adjust <- function(d, p = "p", groups = character(), out = "FDR") {
	if (!nrow(d) || !p %in% names(d)) return(d)
	key <- if (length(groups)) do.call(paste, c(d[groups], sep = "\r")) else rep("all", nrow(d))
	d[[out]] <- NA_real_
	for (g in unique(key)) {
		ii <- which(key == g) ; d[[out]][ii] <- p.adjust(d[[p]][ii], "BH")
	} ; d
}
pgs_cox <- function(d, xs, covars, fold = FALSE, min_events = pgs_num("PGS_MIN_EVENTS", 20)) {
	covars <- unique(setdiff(covars, xs)) ; need <- c(".time", ".event", xs, covars)
	empty <- data.frame(
		term = xs, beta = NA_real_, se = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_,
		N = nrow(d), events = sum(d$.event == 1, na.rm = TRUE), status = "insufficient data", warning = "", dropped_constant_covariates = ""
	)
	result <- list(effects = empty, contrast = data.frame(), loglik = NA_real_, vcov = NULL)
	if (length(setdiff(need, names(d)))) {
		result$effects$status <- paste("missing covariates:", paste(setdiff(need, names(d)), collapse = ",")) ; return(result)
	}
	# Do not silently change the common sample between nested models.
	if (!all(pgs_complete(d, need)) || any(d$.time <= 0) || any(!d$.event %in% c(0, 1))) {
		result$effects$status <- "invalid common model sample" ; return(result)
	}
	if (nrow(d) < pgs_num("PGS_MIN_N", 200) || sum(d$.event) < min_events) return(result)
	if (any(vapply(d[xs], function(x) !is.finite(sd(x)) || sd(x) <= 1e-12, logical(1)))) {
		result$effects$status <- "constant exposure" ; return(result)
	}
	drop <- covars[vapply(d[covars], function(x) length(unique(x)) < 2L, logical(1))] ; covars <- setdiff(covars, drop)
	rhs <- paste(sprintf("`%s`", c(xs, covars)), collapse = " + ")
	if (fold && length(unique(d$.fold)) > 1) rhs <- paste(rhs, "strata(.fold)", sep = " + ")
	ff <- as.formula(paste("survival::Surv(.time,.event) ~", rhs)) ; environment(ff) <- environment()
	strata <- survival::strata ; warnings <- character()
	args <- list(
		formula = ff, data = d, ties = "efron", model = FALSE, x = FALSE, y = FALSE,
		control = survival::coxph.control(iter.max = 40)
	)
	if (".group" %in% names(d) && anyDuplicated(d$.group)) {
		args$cluster <- d$.group ; args$robust <- TRUE
	}
	fit <- tryCatch(withCallingHandlers(do.call(survival::coxph, args), warning = function(w) {
		warnings <<- c(warnings, conditionMessage(w)) ; invokeRestart("muffleWarning")
	}), error = function(e) e)
	if (inherits(fit, "error")) {
		result$effects$status <- paste("fit failed:", conditionMessage(fit)) ; return(result)
	}
	cc <- coef(fit) ; V <- vcov(fit) ; ii <- match(xs, names(cc)) ; b <- unname(cc[ii]) ; s <- sqrt(diag(V)[ii])
	valid <- is.finite(b) & is.finite(s) & s > 0
	if (any(grepl("infinite|converg", warnings, ignore.case = TRUE))) valid[] <- FALSE
	result$effects$beta <- ifelse(valid, b, NA_real_) ; result$effects$se <- ifelse(valid, s, NA_real_)
	result$effects$p <- ifelse(valid, 2 * pnorm( - abs(b / s)), NA_real_)
	result$effects$lo <- ifelse(valid, b - 1.96 * s, NA_real_) ; result$effects$hi <- ifelse(valid, b + 1.96 * s, NA_real_)
	result$effects$status <- ifelse(valid, "ok", "non-estimable or convergence warning")
	result$effects$warning <- paste(unique(warnings), collapse = "; ")
	result$effects$dropped_constant_covariates <- paste(drop, collapse = ";")
	result$loglik <- fit$loglik[2] ; result$vcov <- V
	if (all(c(".G", ".R") %in% xs) && all(valid)) {
		delta <- unname(cc[".G"] - cc[".R"]) ; v <- V[".G", ".G"] + V[".R", ".R"] - 2 * V[".G", ".R"]
		result$contrast <- data.frame(
			beta_difference = delta, se_difference = if (v > 0) sqrt(v) else NA_real_,
			covariance_GR = V[".G", ".R"], p = if (v > 0) 2 * pnorm( - abs(delta / sqrt(v))) else NA_real_,
			lo = if (v > 0) delta - 1.96 * sqrt(v) else NA_real_, hi = if (v > 0) delta + 1.96 * sqrt(v) else NA_real_,
			N = nrow(d), events = sum(d$.event), uncertainty = "Conditional on cross-fitted calibration; see refitted bootstrap"
		)
	}
	result
}
pgs_incident <- function(d, cols) d[pgs_complete(d, c(cols, ".time", ".event", ".prev")) &
	d$.prev %in% 0 & is.finite(d$.time) & d$.time > 0 & d$.event %in% c(0, 1), , drop = FALSE]
pgs_pair <- function(d, covars, feature, adjustment = "basic", scope = "existing_full") {
	z <- pgs_incident(d, c(".m", ".g", covars)) ; id <- pgs_hash(sort(z$eid))
	z$.M <- as.numeric(scale(z$.m)) ; z$.P <- as.numeric(scale(z$.g))
	fits <- list(measured = pgs_cox(z, ".M", covars), pgs = pgs_cox(z, ".P", covars), joint = pgs_cox(z, c(".M", ".P"), covars))
	rows <- lapply(names(fits), function(model) {
		a <- fits[[model]]$effects ; a$feature <- feature ; a$model <- model ; a$adjustment <- adjustment ; a$scope <- scope ; a$sample_hash <- id ; a$covariates <- paste(covars, collapse = ";") ; a$unit <- "log HR per own SD in the identical complete-case sample" ; a
	})
	a <- pgs_bind(rows) ; m <- fits$measured$effects ; g <- fits$pgs$effects ; j <- fits$joint$effects
	lr <- 2 * (fits$joint$loglik - fits$measured$loglik)
	s <- data.frame(feature, scope, adjustment,
		N = nrow(z), events = sum(z$.event), sample_hash = id,
		measured_beta = m$beta, measured_se = m$se, measured_p = m$p, pgs_beta = g$beta, pgs_se = g$se, pgs_p = g$p,
		measured_joint_beta = j$beta[match(".M", j$term)], pgs_joint_beta = j$beta[match(".P", j$term)],
		pgs_joint_p = j$p[match(".P", j$term)], pgs_increment_LRT_p = if (!anyDuplicated(z$.group) && is.finite(lr)) pchisq(max(0, lr), 1, lower.tail = FALSE) else NA_real_,
		LRT_status = if (anyDuplicated(z$.group)) "withheld for clustered observations; use robust joint Wald" else "iid partial-likelihood diagnostic",
		correlation = if (nrow(z) > 2 && all(is.finite(z$.M)) && all(is.finite(z$.P))) cor(z$.M, z$.P) else NA_real_,
		covariates = paste(covars, collapse = ";"), comparison = "identical people, follow-up and adjustment; own-SD associations, not MR"
	)
	list(effects = a, summary = s)
}
pgs_calibrate <- function(d, covars, k = 5L, seed = 2026L, diagnostics = TRUE) {
	d <- as.data.frame(d)
	if (!".group" %in% names(d)) d$.group <- d$eid
	if (!".fold" %in% names(d)) d$.fold <- pgs_folds(d$.group, k, seed)
	for (nm in c(".M", ".G", ".R", ".pred", ".pred0")) d[[nm]] <- NA_real_
	reports <- list() ; ok <- pgs_complete(d, c(".m", ".g", covars))
	for (f in sort(unique(d$.fold))) {
		train <- d[d$.fold != f & d$.prev %in% 0 & ok, , drop = FALSE] ; ii <- which(d$.fold == f & ok)
		if (nrow(train) < 100 || !length(ii) || sd(train$.m) <= 1e-12 || sd(train$.g) <= 1e-12) next
		om <- mean(train$.m) ; os <- sd(train$.m) ; gm <- mean(train$.g) ; gs <- sd(train$.g)
		train$.y <- (train$.m - om) / os ; train$.p <- (train$.g - gm) / gs
		cv <- covars[vapply(train[covars], function(x) length(unique(x)) > 1, logical(1))]
		fit <- tryCatch(lm(reformulate(c(".p", cv), ".y"), train), error = function(e) NULL)
		fit0 <- if (diagnostics) tryCatch(lm(reformulate(if (length(cv)) cv else "1", ".y"), train), error = function(e) NULL) else NULL
		if (is.null(fit) || (diagnostics && is.null(fit0)) || !is.finite(coef(fit)[".p"])) next
		test <- d[ii, , drop = FALSE] ; test$.p <- (test$.g - gm) / gs ; b <- unname(coef(fit)[".p"])
		p <- if (diagnostics) tryCatch(as.numeric(predict(fit, test)), error = function(e) rep(NA_real_, nrow(test))) else rep(NA_real_, nrow(test))
		p0 <- if (diagnostics) tryCatch(as.numeric(predict(fit0, test)), error = function(e) rep(NA_real_, nrow(test))) else rep(NA_real_, nrow(test))
		d$.M[ii] <- (test$.m - om) / os ; d$.G[ii] <- b * test$.p ; d$.R[ii] <- d$.M[ii] - d$.G[ii]
		d$.pred[ii] <- p ; d$.pred0[ii] <- p0
		reports[[length(reports) + 1L]] <- data.frame(
			fold = f, N_train = nrow(train), N_test = length(ii), slope = b,
			omic_mean = om, omic_sd = os, pgs_mean = gm, pgs_sd = gs,
			calibration = "baseline disease-free training groups; future outcomes not used"
		)
	}
	ii <- d$.prev %in% 0 & pgs_complete(d, c(".M", ".pred", ".pred0"))
	mse0 <- sum((d$.M[ii] - d$.pred0[ii]) ^ 2) ; r2 <- if (mse0 > 0) 1 - sum((d$.M[ii] - d$.pred[ii]) ^ 2) / mse0 else NA_real_
	folds <- pgs_bind(reports) ; b <- folds$slope
	audit <- data.frame(
		N_calibration = sum(ii), partial_R2 = r2,
		slope_mean = if (length(b)) mean(b) else NA_real_, slope_min = if (length(b)) min(b) else NA_real_, slope_max = if (length(b)) max(b) else NA_real_,
		orientation = if (!length(b)) "unavailable" else if (all(b > 0)) "positive in all folds" else if (all(b < 0)) "negative in all folds" else "unstable across folds",
		weak_capture = !is.finite(r2) || r2 < pgs_num("PGS_MIN_PARTIAL_R2", .005),
		identity_max_error = if (any(is.finite(d$.M))) max(abs(d$.M - d$.G - d$.R), na.rm = TRUE) else NA_real_,
		calibration_folds = length(b), unit = "training-fold whole-biomarker SD; G/R not separately standardized"
	)
	list(data = d, folds = folds, audit = audit)
}
pgs_component_models <- function(d, covars, feature, scope = "existing_full", landmark = 0, end = Inf) {
	z <- pgs_incident(d, c(".M", ".G", ".R", covars)) ; z <- z[z$.time > landmark, , drop = FALSE]
	z$.event <- as.integer(z$.event == 1 & z$.time <= end) ; z$.time <- pmin(z$.time, end) - landmark
	fits <- list(
		measured = pgs_cox(z, ".M", covars, TRUE), genetic = pgs_cox(z, ".G", covars, TRUE),
		remaining = pgs_cox(z, ".R", covars, TRUE), joint = pgs_cox(z, c(".G", ".R"), covars, TRUE)
	)
	effects <- pgs_bind(lapply(names(fits), function(nm) {
		a <- fits[[nm]]$effects ; a$model <- nm ; a
	}))
	for (nm in c("effects", "contrast")) {
		a <- if (nm == "effects") effects else fits$joint$contrast
		if (nrow(a)) {
			a$feature <- feature ; a$scope <- scope ; a$landmark <- landmark ; a$end <- end ; a$sample_hash <- pgs_hash(sort(z$eid)) ; a$covariates <- paste(covars, collapse = ";")
		}
		if (nm == "effects") effects <- a else contrast <- a
	}
	list(effects = effects, contrast = contrast)
}
pgs_bootstrap <- function(d, covars, feature, B = 100L, k = 5L, seed = 2026L) {
	if (B < 1) return(data.frame())
	if (!".group" %in% names(d)) d$.group <- d$eid
	if (!".fold" %in% names(d)) d$.fold <- pgs_folds(d$.group, k, seed)
	groups <- split(seq_len(nrow(d)), d$.group) ; rows <- vector("list", B)
	for (b in seq_len(B)) {
		set.seed(seed + b) ; ii <- unlist(groups[sample(seq_along(groups), length(groups), replace = TRUE)], use.names = FALSE)
		# Repeated copies of an individual/family always remain in the same fold.
		cal <- pgs_calibrate(d[ii, , drop = FALSE], covars, k, seed, diagnostics = FALSE)
		fitdata <- pgs_incident(cal$data, c(".G", ".R", covars))
		res <- pgs_cox(fitdata, c(".G", ".R"), covars, TRUE)
		z <- res$effects ; delta <- res$contrast$beta_difference
		rows[[b]] <- data.frame(
			replicate = b, genetic = if (nrow(z)) z$beta[match(".G", z$term)] else NA_real_,
			remaining = if (nrow(z)) z$beta[match(".R", z$term)] else NA_real_, difference = if (length(delta)) delta[1] else NA_real_
		)
	}
	draws <- pgs_bind(rows)
	pgs_bind(lapply(c("genetic", "remaining", "difference"), function(term) {
		x <- draws[[term]] ; x <- x[is.finite(x)] ; n <- length(x) ; valid <- n >= max(30, ceiling(.8 * B))
		data.frame(feature, term,
			requested = B, successful = n, lo = if (valid) unname(quantile(x, .025)) else NA_real_,
			hi = if (valid) unname(quantile(x, .975)) else NA_real_,
			sign_tail = if (valid) min(1, 2 * min((sum(x <= 0) + 1) / (n + 1), (sum(x >= 0) + 1) / (n + 1))) else NA_real_,
			status = if (valid) "ok" else "too few successful refits",
			uncertainty = "Participant/family resampling with calibration refit; conditional on selected candidate and original GWAS weights"
		)
	}))
}
pgs_lifestyle <- function(d, covars, components, feature) {
	pgs_bind(lapply(components, function(component) pgs_bind(lapply(c(".M", ".G", ".R"), function(part) {
		cols <- unique(c(part, component, covars)) ; z <- d[d$.prev %in% 0 & pgs_complete(d, cols), , drop = FALSE]
		ans <- data.frame(feature, component, part, N = nrow(z), beta = NA_real_, se = NA_real_, p = NA_real_, status = "unavailable")
		if (nrow(z) < 200 || !component %in% names(z) || sd(z[[component]]) <= 0 || sd(z[[part]]) <= 0) return(ans)
		z$.x <- as.numeric(scale(z[[component]])) ; z$.y <- z[[part]]
		cv <- setdiff(covars, component) ; cv <- cv[vapply(z[cv], function(x) length(unique(x)) > 1, logical(1))]
		f <- tryCatch(lm(reformulate(c(".x", cv), ".y"), z), error = function(e) NULL)
		if (is.null(f)) return(ans) ; s <- coef(summary(f)) ; if (!".x" %in% rownames(s)) return(ans)
		ans$beta <- s[".x", 1] ; ans$se <- s[".x", 2] ; ans$p <- s[".x", 4]
		if (anyDuplicated(z$.group)) {
			X <- model.matrix(f) ; B <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
			if (!is.null(B)) {
				u <- rowsum(X * as.numeric(residuals(f)), z$.group, reorder = FALSE) ; ng <- nrow(u)
				V <- B %*% crossprod(u) %*% B
				ans$se <- sqrt(V[".x", ".x"] * ng / (ng - 1) * (nrow(z) - 1) / (nrow(z) - ncol(X)))
				ans$p <- 2 * pt( - abs(ans$beta / ans$se), df = ng - 1)
			} else {
				ans$se <- NA_real_ ; ans$p <- NA_real_
			}
		}
		ans$status <- "association only; not mediation or environmental fraction" ; ans
	}))))
}
pgs_prevalent <- function(d, covars, feature) {
	# Baseline case-control association is kept separate from incident hazards.
	z <- d[pgs_complete(d, c(".M", ".G", ".R", ".prev", covars)) & d$.prev %in% c(0, 1), , drop = FALSE]
	if (nrow(z) < pgs_num("PGS_MIN_N", 200) || min(table(factor(z$.prev, levels = 0 : 1))) < pgs_num("PGS_MIN_EVENTS", 20)) return(data.frame())
	cv <- covars[vapply(z[covars], function(x) length(unique(x)) > 1, logical(1))]
	pgs_bind(lapply(list(measured = ".M", joint = c(".G", ".R")), function(xs) {
		if (any(vapply(z[xs], function(x) sd(x) <= 1e-12, logical(1)))) return(data.frame())
		warns <- character() ; f <- tryCatch(withCallingHandlers(glm(reformulate(c(xs, cv), ".prev"), data = z, family = binomial()),
			warning = function(w) {
				warns <<- c(warns, conditionMessage(w)) ; invokeRestart("muffleWarning")
			}
		), error = function(e) NULL)
		if (is.null(f) || !f$converged || length(warns)) return(data.frame(feature, status = "prevalent logistic non-estimable"))
		V <- vcov(f)
		if (anyDuplicated(z$.group)) {
			X <- model.matrix(f) ; score <- X * as.numeric(z$.prev - fitted(f)) ; u <- rowsum(score, z$.group, reorder = FALSE)
			V <- V %*% crossprod(u) %*% V
		}
		b <- coef(f)[xs] ; se <- sqrt(diag(V)[xs]) ; data.frame(feature,
			term = xs, model = if (length(xs) == 1) "measured" else "joint",
			beta = unname(b), se = unname(se), p = 2 * pnorm( - abs(b / se)), lo = b - 1.96 * se, hi = b + 1.96 * se,
			N = nrow(z), cases = sum(z$.prev), status = "ok", unit = "log OR per training-fold biomarker SD; not comparable in magnitude with log HR"
		)
	}))
}


# c1.correlate.R
# Called by c1.correlate.R; independent pgs_focus entry avoids rerunning MR/ML.
pgs_focus_anchors <- function(layer) {
	default <- if (layer == "protein") "PCSK9,LPA,GDF15,NTPROBNP,MMP12,CCL19,CCL21,CXCL13,AGER,IL15,TNFRSF4" else
		"L_VLDL_TG.pct,L_VLDL_TG,Total_TG,ApoB,VLDL_size,GlycA,Albumin,Phe,Lactate,DHA,LA.pct"
	pgs_csv(if (layer == "protein") "PGS_PROT_ANCHORS" else "PGS_MET_ANCHORS", default)
}
pgs_focus_selection <- function(summary, layer) {
	d <- summary ; ranked <- if (nrow(d)) d$feature[order(!(d$evidence_pattern == "Both supported: opposite"),
		pmax(d$measured_FDR, d$pgs_FDR), d$measured_p,
		na.last = TRUE
	)] else character()
	unique(c(intersect(pgs_focus_anchors(layer), d$feature), head(ranked, pgs_num("PGS_DEEP_MAX", 24))))
}
pgs_focus_export <- function(tables, root, prefix = "c1.pgs_focus") {
	dir.create(root, recursive = TRUE, showWarnings = FALSE)
	for (nm in names(tables)) {
		d <- tables[[nm]] ; if (!is.data.frame(d)) next
		if (!ncol(d)) d <- data.frame(status = "unavailable / not estimated")
		data.table::fwrite(d, file.path(root, paste0(prefix, ".", nm, ".csv")))
	}
	wb <- openxlsx::createWorkbook()
	for (nm in names(tables)) {
		d <- tables[[nm]] ; if (!is.data.frame(d)) next
		if (!ncol(d)) d <- data.frame(status = "unavailable / not estimated")
		# Excel has a hard row limit; CSV retains all rows if a large family exceeds it.
		if (nrow(d) > 1048500) d <- data.frame(status = "Full table exceeds Excel row limit; see companion CSV", rows = nrow(d))
		sh <- substr(nm, 1, 31) ; openxlsx::addWorksheet(wb, sh)
		if (nrow(d)) openxlsx::writeDataTable(wb, sh, d) else openxlsx::writeData(wb, sh, d)
		openxlsx::freezePane(wb, sh, firstRow = TRUE) ; openxlsx::setColWidths(wb, sh, seq_len(ncol(d)), 18)
	}
	openxlsx::saveWorkbook(wb, file.path(root, paste0(prefix, ".xlsx")), overwrite = TRUE)
}
pgs_focus_deep <- function(d, covars, feature, layer, boot = 0L) {
	k <- as.integer(pgs_num("PGS_FOLDS", 5)) ; seed <- as.integer(pgs_num("SEED", 2026))
	cal <- pgs_calibrate(d, covars, k, seed) ; cal$audit$feature <- feature ; cal$folds$feature <- rep(feature, nrow(cal$folds))
	main <- pgs_component_models(cal$data, covars, feature)
	cuts <- sort(unique(as.numeric(pgs_csv("PGS_TIME_CUTS", "0,2,5,10,Inf"))))
	if (length(cuts) < 2 || anyNA(cuts) || cuts[1] != 0 || any(diff(cuts) <= 0)) stop("Invalid PGS_TIME_CUTS")
	windows <- lapply(seq_len(length(cuts) - 1L), function(i) pgs_component_models(cal$data, covars, feature, landmark = cuts[i], end = cuts[i + 1]))
	landmarks <- lapply(c(2, 5), function(L) pgs_component_models(cal$data, covars, feature, landmark = L))
	sens <- list() ; sensitivity_sets <- list(
		LE4 = unique(c(covars, pgs_csv("PGS_LE4_COVARS", "diet.pts,pa.pts,smoke.pts,sleep.pts"))),
		LE8 = unique(c(covars, get0("vars.le8", ifnotfound = character()))),
		treatment = unique(c(covars, pgs_csv("PGS_TREATMENT_COVARS", "drug.lipid,drug.htn,drug.dm")))
	)
	for (nm in names(sensitivity_sets)) {
		cv <- sensitivity_sets[[nm]] ; z <- d[pgs_complete(d, c(".m", ".g", cv)), , drop = FALSE]
		if (identical(cv, covars)) next
		# Both models use the SAME restricted sample; no silent covariate dropping.
		a <- pgs_pair(z, covars, feature, paste0(nm, "_restricted_basic"))
		b <- pgs_pair(z, cv, feature, paste0(nm, "_adjusted"))
		sens[[nm]] <- pgs_bind(list(a$summary, b$summary))
	}
	# Outlier influence sensitivity is separate from the primary preprocessing.
	w <- d ; ii <- w$.prev %in% 0 & is.finite(w$.m)
	if (sum(ii) > 200) {
		lim <- quantile(w$.m[ii], c(.01, .99)) ; w$.m <- pmin(pmax(w$.m, lim[1]), lim[2]) ; sens$winsor <- pgs_pair(w, covars, feature, "measured_winsor_1_99")$summary
	}
	life <- pgs_lifestyle(cal$data, covars, intersect(get0("vars.le8", ifnotfound = character()), names(d)), feature)
	list(
		calibration = cal$audit, folds = cal$folds, components = main$effects, contrasts = main$contrast,
		prevalent = pgs_prevalent(cal$data, covars, feature),
		windows = pgs_bind(lapply(windows, `[[`, "effects")), window_contrasts = pgs_bind(lapply(windows, `[[`, "contrast")),
		landmarks = pgs_bind(lapply(landmarks, `[[`, "effects")), landmark_contrasts = pgs_bind(lapply(landmarks, `[[`, "contrast")),
		sensitivity = pgs_bind(sens), lifestyle = life,
		bootstrap = pgs_bootstrap(d, covars, feature, B = boot, k = k, seed = seed)
	)
}
run_c1_pgs_focus <- function(layer, outdir) {
	rd <- file.path(outdir, "c1_correlate") ; dir.create(rd, recursive = TRUE, showWarnings = FALSE)
	sf <- find_c1_pgs_file(layer)
	if (is.na(sf)) {
		out <- list(status = data.frame(status = "unavailable", detail = "Biomarker PGS input missing"))
		pgs_focus_export(out, rd) ; saveRDS(out, file.path(rd, "c1.pgs_focus.rds")) ; stop("Biomarker PGS input missing for ", layer, "; focused analysis stopped")
	}
	primary <- pgs_csv("PGS_BASIC_COVARS", paste(if (length(get0("le8_custom_covars", ifnotfound = character()))) le8_custom_covars else vars.basic, collapse = ","))
	# A supplied genetic score is not a disease PRS, even when its name is similar.
	scores <- pgs_ids(read_c1_pgs(sf), "PGS") ; biom <- pgs_ids(if (layer == "protein") read_prot() else read_met(), "Measured omics")
	feats <- setdiff(names(biom), "eid") ; mp <- map_c1_pgs_columns(feats, names(scores))
	only <- pgs_csv("PGS_FEATURES") ; if (length(only)) mp <- mp[intersect(names(mp), only)]
	if (!length(mp)) stop("No exact measured/PGS feature matches for ", layer)
	ph <- pgs_ids(read_all(), "Phenotypes") |>
		filter_analysis_cohort() |>
		make_outcome(Y)
	if (length(setdiff(primary, names(ph)))) stop("Missing primary PGS covariates: ", paste(setdiff(primary, names(ph)), collapse = ","))
	ph$.prev <- make_prevalent_status(ph, Y) ; ph$.time <- ph[[paste0(Y, ".t2e")]] ; ph$.event <- ph[[paste0(Y, ".Yt2e")]]
	# A diagnosis count without an onset date is an unresolved baseline state,
	# not a healthy control for calibration or prevalent analysis.
	dc <- if (Y %in% names(ph)) Y else paste0("icd10Ct_", Y)
	yd <- paste0("fod_icd10_", Y)
	if (all(c(dc, yd) %in% names(ph))) {
		positive <- if (inherits(ph[[dc]], "Date")) !is.na(ph[[dc]]) else suppressWarnings(as.numeric(ph[[dc]])) > 0
		unknown <- is.na(ph[[yd]]) & positive %in% TRUE ; ph$.prev[unknown] <- NA_real_
	}
	groupcol <- Sys.getenv("PGS_GROUP_COLUMN", "")
	if (nzchar(groupcol) && !groupcol %in% names(ph)) stop("PGS_GROUP_COLUMN missing: ", groupcol)
	ph$.group <- if (nzchar(groupcol)) as.character(ph[[groupcol]]) else ph$eid
	if (anyNA(ph$.group) || any(!nzchar(ph$.group))) stop("Missing PGS group IDs; do not silently split related participants")
	need <- unique(c(
		"eid", ".prev", ".time", ".event", ".group", primary, vars.le8,
		pgs_csv("PGS_TREATMENT_COVARS", "drug.lipid,drug.htn,drug.dm")
	))
	base <- ph[, intersect(need, names(ph)), drop = FALSE] ; rm(ph) ; invisible(gc())
	ids <- intersect(intersect(base$eid, biom$eid), scores$eid)
	base <- base[match(ids, base$eid), , drop = FALSE]
	im <- match(ids, biom$eid) ; ig <- match(ids, scores$eid)
	base$.fold <- pgs_folds(base$.group, as.integer(pgs_num("PGS_FOLDS", 5)), as.integer(pgs_num("SEED", 2026)))
	inputs <- c(
		sf, file.path(indir, "Rdata", c("all.rds", if (layer == "protein") "prot.rds" else "met.rds")),
		Sys.getenv("LE8_BASELINE_MED_FILE", file.path(indir, "rap/vip.tab.gz"))
	)
	sig <- pgs_hash(list(
		version = PGS_FOCUS_VERSION, files = pgs_stamp(inputs), Y = Y, primary = primary,
		baseline = LE8_BASELINE_VERSION, options = le8_analysis_options(), features = names(mp), group = groupcol,
		end = as.character(date_follow_end), folds = base$.fold, settings = Sys.getenv(c("PGS_FOLDS", "PGS_DEEP_MAX", "PGS_MIN_N", "PGS_MIN_EVENTS", "PGS_BOOT", "PGS_BOOT_FEATURES", "PGS_TIME_CUTS", "PGS_LE4_COVARS", "PGS_TREATMENT_COVARS", "LE8_BASELINE_MED_COLUMNS", "LE8_BASELINE_MAP", "PGS_PROT_ANCHORS", "PGS_MET_ANCHORS", "PGS_MIN_PARTIAL_R2", "PGS_GWAS_SOURCE", "PGS_GWAS_OVERLAP"))
	))
	cache <- le8_cache_dir("pgs_focus", layer, sig) ; dir.create(cache, recursive = TRUE, showWarnings = FALSE)
	getdata <- function(f) {
		d <- base ; d$.m <- suppressWarnings(as.numeric(biom[[f]][im])) ; d$.g <- suppressWarnings(as.numeric(scores[[mp[[f]]]][ig])) ; d
	}
	message("[LE8] START PGS/", layer, " matched scan | ", length(mp), " features; ", nrow(base), " matched participants")
	scanfile <- file.path(cache, "scan.rds") ; pairs <- if (file.exists(scanfile) && !LE8_REPLACE) readRDS(scanfile) else {
		ans <- lapply(names(mp), function(f) pgs_pair(getdata(f), primary, f)) ; saveRDS(ans, scanfile, compress = FALSE) ; ans
	}
	paired <- pgs_bind(lapply(pairs, `[[`, "summary")) ; effects <- pgs_bind(lapply(pairs, `[[`, "effects")) ; rm(pairs)
	paired$measured_FDR <- p.adjust(paired$measured_p, "BH") ; paired$pgs_FDR <- p.adjust(paired$pgs_p, "BH")
	paired$pgs_joint_FDR <- p.adjust(paired$pgs_joint_p, "BH")
	# Intersection-union test for opposite directions, with two possible
	# orientations corrected before BH over the full matched scan.
	zm <- paired$measured_beta / paired$measured_se ; zg <- paired$pgs_beta / paired$pgs_se
	paired$opposite_conjunction_p <- pmin(1, 2 * pmin(pmax(pnorm( - zm), pnorm(zg)), pmax(pnorm(zm), pnorm( - zg))))
	paired$opposite_conjunction_FDR <- p.adjust(paired$opposite_conjunction_p, "BH")
	paired$opposite_conjunction_supported <- is.finite(paired$opposite_conjunction_FDR) & paired$opposite_conjunction_FDR < .05
	paired$evidence_pattern <- pgs_classify(paired$measured_beta, paired$pgs_beta, paired$measured_FDR, paired$pgs_FDR)
	effects <- pgs_adjust(effects, groups = c("model", "term"))
	candidates <- pgs_focus_selection(paired, layer) ; paired$deep_selected <- paired$feature %in% candidates
	selection <- data.frame(
		feature = candidates, selection = ifelse(candidates %in% pgs_focus_anchors(layer), "declared anchor", "exploratory same-sample ranking"),
		caveat = "Same data used for candidate selection; candidate-only FDR is not independent confirmation"
	)
	defaultboot <- if (layer == "protein") "PCSK9,LPA,CCL19,CCL21" else "L_VLDL_TG.pct,GlycA,Lactate"
	bootfeatures <- pgs_csv("PGS_BOOT_FEATURES", defaultboot) ; B <- as.integer(pgs_num("PGS_BOOT", 100))
	message("[LE8] DONE PGS/", layer, " matched scan | opposite=", sum(paired$evidence_pattern == "Both supported: opposite"))
	message("[LE8] START PGS/", layer, " component analysis | ", length(candidates), " candidates; bootstrap=", B, " for ", length(intersect(candidates, bootfeatures)), " anchors")
	deep <- lapply(candidates, function(f) {
		cf <- file.path(cache, paste0("feature_", pgs_hash(f), ".rds"))
		if (file.exists(cf) && !LE8_REPLACE) return(readRDS(cf))
		z <- pgs_focus_deep(getdata(f), primary, f, layer, if (f %in% bootfeatures) B else 0L)
		tmp <- paste0(cf, ".tmp") ; saveRDS(z, tmp, compress = FALSE) ; if (!file.rename(tmp, cf)) stop("Cannot save PGS feature checkpoint") ; z
	})
	nm <- unique(unlist(lapply(deep, names))) ; out <- setNames(lapply(nm, function(n) pgs_bind(lapply(deep, `[[`, n))), nm)
	out$paired <- paired ; out$matched_models <- effects ; out$selection <- selection
	reference_file <- file.path(rd, paste0(if (layer == "protein") "pwas" else "mwas", "_pgs_incident_full_genetic.csv"))
	if (file.exists(reference_file)) {
		ref <- as.data.frame(data.table::fread(reference_file, showProgress = FALSE))
		if (nrow(ref)) {
			ref$comparison_role <- "Historical full-genetic-cohort reference; not matched to current measured analysis; see original adjustment metadata"
			out$full_cohort_reference <- ref
		}
	}
	if (layer == "metabolite") out$composition <- pgs_composition_analysis(base, biom, im, mp, getdata, primary)
	for (n in intersect(c("components", "windows", "landmarks", "prevalent"), names(out))) out[[n]] <- pgs_adjust(out[[n]], groups = intersect(c("scope", "model", "term", "landmark", "end"), names(out[[n]])))
	for (n in intersect(c("contrasts", "window_contrasts", "landmark_contrasts"), names(out))) out[[n]] <- pgs_adjust(out[[n]], groups = intersect(c("scope", "landmark", "end"), names(out[[n]])))
	out$lifestyle <- pgs_adjust(out$lifestyle, groups = intersect(c("part"), names(out$lifestyle)))
	out$status <- data.frame(
		status = "ok", version = PGS_FOCUS_VERSION, signature = sig, layer, outcome = Y, N_matched = nrow(base), features = length(mp), estimable_pairs = sum(is.finite(paired$measured_p) & is.finite(paired$pgs_p)), deep_features = length(candidates),
		adjustment = paste(primary, collapse = ";"), PGS_source = sf, GWAS_source = Sys.getenv("PGS_GWAS_SOURCE", "unknown"),
		GWAS_overlap = Sys.getenv("PGS_GWAS_OVERLAP", "unknown"), group_column = if (nzchar(groupcol)) groupcol else "eid; family structure not supplied",
		interpretation = "Association decomposition, not causal partition or biomarker concentration at birth"
	)
	# Source-specific PGSs are built only with verified build/alleles/genotypes.
	if (exists("pgs_source_analysis", mode = "function")) {
		src <- pgs_source_analysis(layer, outdir, candidates, base, biom, im, primary, paired, getdata)
		out <- c(out, src)
	}
	pgs_focus_export(out, rd) ; saveRDS(out, file.path(rd, "c1.pgs_focus.rds"), compress = "gzip")
	message("[LE8] DONE PGS/", layer, " component analysis")
	out
}

pgs_composition_analysis <- function(base, biom, im, mp, getdata, covars) {
	targets <- intersect(pgs_csv("PGS_COMPOSITION_TARGETS", "L_VLDL_TG.pct"), names(mp))
	controls <- pgs_csv("PGS_COMPOSITION_CONTROLS", "L_VLDL_TG,Total_TG,ApoB,VLDL_size")
	rows <- list()
	for (f in targets) for (v in controls) for (kind in c("measured", "pgs")) {
		tag <- paste(f, v, kind)
		if (!v %in% names(biom) || (kind == "pgs" && !v %in% names(mp))) {
			rows[[tag]] <- data.frame(feature = f, control = v, control_type = kind, status = "control unavailable") ; next
		}
		d <- getdata(f) ; d$.burden <- if (kind == "measured") as.numeric(biom[[v]][im]) else getdata(v)$.g
		d <- d[pgs_complete(d, c(".m", ".g", ".burden", covars)), , drop = FALSE]
		a <- pgs_pair(d, covars, f, "burden_restricted_basic")$summary
		b <- pgs_pair(d, c(covars, ".burden"), f, "burden_conditioned")$summary
		z <- pgs_bind(list(a, b)) ; z$control <- v ; z$control_type <- kind
		z$status <- "Conditional sensitivity; potential mediator/collider adjustment, not a causal direct effect"
		rows[[tag]] <- z
	}
	pgs_bind(rows)
}


# c1.correlate.R
# Source-specific score sensitivity. Existing full PGSs cannot be algebraically
# relabelled as cis/non-MHC scores: retained alleles must actually be rescored.
pgs_shared_mhc <- function(build) {
	if (!build %in% c("37", "38")) stop("Verified build 37/38 required for MHC annotation")
	get <- function(n) {
		v <- Sys.getenv(n, "") ; if (!nzchar(v)) {
			x <- get0(n, ifnotfound = NULL) ; if (!is.null(x)) v <- as.character(x)
		} ; suppressWarnings(as.numeric(v))
	}
	a <- get(paste0("MHC_START_b", build)) ; b <- get(paste0("MHC_END_b", build))
	if (length(a) != 1 || length(b) != 1 || !is.finite(a) || !is.finite(b) || a >= b)
		stop("Shared MHC_START_b", build, "/MHC_END_b", build, " unavailable: run through le8.sh (sources PHE_F). No substitute interval is used.")
	c(start = a, end = b)
}
pgs_source_masks <- function(w, build, regions = data.frame()) {
	mhc <- pgs_shared_mhc(build) ; ch <- sub("^chr", "", as.character(w$CHR)) ; pos <- as.numeric(w$POS)
	if (any(!is.finite(pos)) || anyNA(ch)) stop("Missing source-score variant coordinates")
	in_mhc <- ch == "6" & pos >= mhc[1] & pos <= mhc[2]
	kind <- sub("^.* / ", "", w$component)
	masks <- list(full_recomputed = rep(TRUE, nrow(w)), no_MHC = !in_mhc)
	# Protein gene-cis and metabolite lead-local remain explicitly distinct.
	for (k in intersect(c("cis", "trans", "local", "distal"), unique(kind))) {
		masks[[k]] <- kind %in% k
		if (k %in% c("trans", "distal")) masks[[paste0(k, "_no_MHC")]] <- kind %in% k & !in_mhc
	}
	if (nrow(regions)) for (i in seq_len(nrow(regions))) {
		r <- regions[i, ] ; if (!all(c("name", "chr", "start", "end", "build") %in% names(r))) stop("Invalid PGS_EXCLUDE_REGIONS schema")
		if (as.character(r$build) != build) stop("Exclusion-region build differs from score weights")
		if (!is.finite(r$start) || !is.finite(r$end) || r$start > r$end) stop("Invalid exclusion interval")
		masks[[paste0("without_", make.names(r$name))]] <- !(ch == sub("^chr", "", as.character(r$chr)) & pos >= r$start & pos <= r$end)
	}
	masks
}
pgs_plink <- function(args, log) {
	exe <- Sys.getenv("PLINK2", Sys.which("plink2")) ; if (!nzchar(exe)) stop("plink2 unavailable")
	rc <- system2(exe, vapply(args, shQuote, character(1)), stdout = log, stderr = log)
	if (rc != 0) stop("PLINK2 failed; see ", log)
}
pgs_pfile_files <- function(prefix) c(paste0(prefix, c(".pgen", ".psam")), if (file.exists(paste0(prefix, ".pvar"))) paste0(prefix, ".pvar") else paste0(prefix, ".pvar.zst"))
pgs_pfile_args <- function(prefix) c("--pfile", prefix, if (!file.exists(paste0(prefix, ".pvar"))) "vzs")
pgs_prepare_genotype <- function(w, ids, pfiledir, cache) {
	allfiles <- unlist(lapply(unique(w$CHR), function(ch) pgs_pfile_files(file.path(pfiledir, paste0("chr", ch)))))
	if (any(!file.exists(allfiles))) stop("Source-score genotype files unavailable: ", paste(head(allfiles[!file.exists(allfiles)], 3), collapse = ";"))
	key <- pgs_hash(list(inputs = pgs_stamp(allfiles), variants = sort(unique(w$SNP)), ids = sort(ids)))
	rd <- file.path(cache, key) ; dir.create(rd, recursive = TRUE, showWarnings = FALSE)
	keep <- file.path(rd, "keep.txt") ; data.table::fwrite(data.frame(`#IID` = ids, check.names = FALSE), keep, sep = "\t")
	for (ch in unique(w$CHR)) {
		prefix <- file.path(rd, paste0("chr", ch)) ; done <- paste0(prefix, ".complete")
		if (file.exists(done) && all(file.exists(pgs_pfile_files(prefix)))) next
		extract <- paste0(prefix, ".variants.txt") ; writeLines(unique(w$SNP[w$CHR == ch]), extract)
		pgs_plink(c(
			pgs_pfile_args(file.path(pfiledir, paste0("chr", ch))), "--extract", extract,
			"--keep", keep, "--make-pgen", "--threads", Sys.getenv("N_CORES", "2"), "--out", prefix
		), paste0(prefix, ".console.log"))
		file.create(done)
	}
	rd
}
pgs_score_sources <- function(w, masks, pfiledir, cache) {
	# Exact ID + effect-allele matching; no strand guessing, no duplicate-ID skip.
	if (anyDuplicated(w$SNP) || any(!grepl("^[ACGT]+$", toupper(w$effect_allele))) || any(!is.finite(w$weight))) stop("Invalid/duplicate score alleles or weights")
	scopes <- names(masks) ; tot <- NULL ; used <- list() ; status <- list()
	for (ch in unique(w$CHR)) {
		ix <- which(as.character(w$CHR) == as.character(ch)) ; a <- w[ix, ] ; prefix <- file.path(cache, paste0("chr", ch)) ; dir.create(cache, recursive = TRUE, showWarnings = FALSE)
		weight <- data.frame(ID = a$SNP, EA = toupper(a$effect_allele))
		for (n in scopes) weight[[n]] <- a$weight * as.numeric(masks[[n]][ix])
		sf <- paste0(prefix, ".weights.tsv") ; data.table::fwrite(weight, sf, sep = "\t")
		pgs_plink(c(
			pgs_pfile_args(file.path(pfiledir, paste0("chr", ch))), "--score", sf, "1", "2", "header-read", "no-mean-imputation",
			"list-variants", "cols=+scoresums", "--score-col-nums", paste0("3-", ncol(weight)), "--threads", Sys.getenv("N_CORES", "2"), "--out", prefix
		), paste0(prefix, ".console.log"))
		vf <- paste0(prefix, ".sscore.vars") ; if (!file.exists(vf)) stop("PLINK did not retain scored-variant list")
		matched <- readLines(vf, warn = FALSE) ; hit <- a$SNP %in% matched
		aa <- a ; aa$matched_by_plink <- hit ; used[[length(used) + 1L]] <- aa
		z <- as.data.frame(data.table::fread(paste0(prefix, ".sscore"), showProgress = FALSE)) ; names(z) <- sub("^#", "", names(z))
		if (!all(c("IID", paste0(scopes, "_SUM")) %in% names(z))) stop("Unexpected PLINK multi-score output schema")
		z$IID <- as.character(z$IID) ; if (anyDuplicated(z$IID)) stop("Duplicate IID in genotype score output")
		sc <- data.frame(eid = z$IID) ; for (n in scopes) sc[[n]] <- z[[paste0(n, "_SUM")]]
		# With no-mean-imputation, flag participants with any missing dosage in
		# the union. The primary source comparison uses this exact common set.
		sc$.complete_genotype <- z$ALLELE_CT >= 2 * sum(hit)
		if (is.null(tot)) tot <- sc else {
			i <- match(tot$eid, sc$eid) ; if (anyNA(i) || nrow(sc) != nrow(tot)) stop("Genotype sample sets differ across chromosomes")
			for (n in scopes) tot[[n]] <- tot[[n]] + sc[[n]][i]
			tot$.complete_genotype <- tot$.complete_genotype & sc$.complete_genotype[i]
		}
	}
	u <- pgs_bind(used)
	for (n in scopes) {
		mask <- masks[[n]][match(u$SNP, w$SNP)]
		status[[n]] <- data.frame(
			scope = n, requested = sum(mask), matched = sum(mask & u$matched_by_plink),
			match_fraction = if (sum(mask)) mean(u$matched_by_plink[mask]) else NA_real_
		)
	}
	list(scores = tot, variants = u, status = pgs_bind(status))
}
pgs_source_analysis <- function(layer, outdir, candidates, base, biom, im, covars, paired, getdata) {
	rd <- file.path(outdir, "c1_correlate") ; wf <- Sys.getenv(
		if (layer == "protein") "PGS_PROT_WEIGHTS" else "PGS_MET_WEIGHTS",
		file.path(outdir, "c2_cause", "c2.genetic_score_weights.tsv")
	)
	status <- list() ; counts <- list() ; variants <- list() ; effects <- list() ; summaries <- list() ; components <- list() ; calibration <- list() ; region_audit <- list() ; contrasts <- list()
	fail <- function(detail) list(source_status = data.frame(status = "unavailable", detail = detail))
	if (toupper(Sys.getenv("PGS_BUILD_SOURCES", "AUTO")) %in% c("FALSE", "NO", "0")) return(fail("Explicitly disabled PGS_BUILD_SOURCES"))
	if (!file.exists(wf)) return(fail(paste("COJO weight manifest missing:", wf)))
	w <- as.data.frame(data.table::fread(wf, showProgress = FALSE))
	if (!all(c("feature", "SNP", "CHR", "POS", "effect_allele", "weight", "component") %in% names(w))) return(fail("Invalid COJO weight manifest schema"))
	w$CHR <- sub("^chr", "", as.character(w$CHR))
	source_default <- unique(c(
		intersect(candidates, c("PCSK9", "LPA", "CCL19", "CCL21", "CXCL13", "L_VLDL_TG.pct", "GlycA", "Lactate")),
		head(paired$feature[paired$evidence_pattern == "Both supported: opposite"], 3)
	))
	chosen <- intersect(pgs_csv("PGS_SOURCE_FEATURES", paste(source_default, collapse = ",")), candidates)
	if (!length(chosen)) return(fail("No source-score candidate among measured/PGS matches"))
	private_cache <- file.path(Sys.getenv("PGS_SOURCE_CACHE_DIR", file.path(indir, ".le8_pgs_source_cache")), layer)
	dir.create(private_cache, recursive = TRUE, showWarnings = FALSE)
	project <- if (layer == "protein") dir.X else dir.met.gwas
	# The C2 manifest retains QTL coordinates; verify against its GWAS QC build.
	build_for <- function(f) {
		declared <- Sys.getenv("PGS_WEIGHTS_BUILD", "")
		qcf <- file.path(project, "common", f, "qc", paste0(f, ".grch"))
		qc <- if (file.exists(qcf)) trimws(readLines(qcf, n = 1, warn = FALSE)) else ""
		qc <- sub("^(GRCh|b)", "", qc, ignore.case = TRUE)
		if (nzchar(declared) && nzchar(qc) && declared != qc) stop("PGS_WEIGHTS_BUILD conflicts with QTL QC metadata")
		b <- if (nzchar(declared)) declared else qc
		if (!b %in% c("37", "38")) stop("Weight build unverified; supply PGS_WEIGHTS_BUILD or QTL qc/<feature>.grch")
		b
	}
	regions_file <- Sys.getenv("PGS_EXCLUDE_REGIONS", "")
	regions <- if (nzchar(regions_file)) as.data.frame(data.table::fread(regions_file)) else data.frame()
	for (f in chosen) {
		tryCatch(
			{
				a <- w[w$feature == f, , drop = FALSE] ; if (!nrow(a)) stop("No COJO weights for candidate")
				if (any(!a$CHR %in% as.character(1 : 22))) stop("Autosomal source analysis requires autosomal weights")
				build <- build_for(f)
				annotation_status <- "Metabolite lead-local/distal; not gene cis"
				if (layer == "protein") {
					bed <- get0("prot_bed_file", ifnotfound = "")
					bbuild <- Sys.getenv("LE8_PROT_BED_BUILD", "")
					if (!nzchar(bbuild) && grepl("(b|[.])38[.]", basename(bed))) bbuild <- "38"
					if (!nzchar(bbuild) && grepl("(b|[.])37[.]", basename(bed))) bbuild <- "37"
					a$component <- "COJO PGS / unknown"
					annotation_status <- "cis/trans withheld: annotation build unverified or differs from QTL build"
					if (file.exists(bed) && identical(build, bbuild)) {
						ann <- as.data.frame(data.table::fread(bed, header = FALSE, showProgress = FALSE))
						hit <- which(as.character(ann[[4]]) == f)
						if (length(hit) == 1 && all(is.finite(as.numeric(ann[hit, 2 : 3])))) {
							pad <- pgs_num("C2_CIS_WINDOW_BP", 1e6)
							cis <- a$CHR == sub("^chr", "", as.character(ann[hit, 1])) & a$POS >= as.numeric(ann[hit, 2]) + 1 - pad & a$POS <= as.numeric(ann[hit, 3]) + pad
							a$component <- paste0("COJO PGS / ", ifelse(cis, "cis", "trans"))
							annotation_status <- paste("Gene annotation verified in build", build)
						}
					}
				}
				local_regions <- regions
				# A leave-hub-out sensitivity uses observed weight-table coordinates,
				# never a hard-coded SH2B3/MHC interval or another build's coordinates.
				for (tag in pgs_csv("PGS_HUB_VARIANTS", "rs3184504")) {
					hit <- w[w$SNP == tag, , drop = FALSE]
					if (nrow(hit)) {
						same <- vapply(hit$feature, function(ff) identical(tryCatch(build_for(ff), error = function(e) ""), build), logical(1))
						loc <- unique(hit[same, c("CHR", "POS"), drop = FALSE])
						if (nrow(loc) == 1) {
							pad <- pgs_num("PGS_HUB_WINDOW_BP", 1e6)
							r <- data.frame(name = paste0(tag, "_region"), chr = loc$CHR, start = max(1, loc$POS - pad), end = loc$POS + pad, build = build)
							overlap <- a$CHR == as.character(r$chr) & a$POS >= r$start & a$POS <= r$end
							if (any(overlap)) local_regions <- pgs_bind(list(local_regions, r))
						}
					}
				}
				if (nrow(local_regions)) {
					rr <- local_regions ; rr$feature <- f ; region_audit[[f]] <- rr
				}
				masks <- pgs_source_masks(a, build, local_regions)
				pd <- Sys.getenv("PGS_PFILE_DIR", file.path("/mnt/f/gen/ukb", build, "imp"))
				genotypebuild <- Sys.getenv("PGS_PFILE_BUILD", if (identical(pd, file.path("/mnt/f/gen/ukb", build, "imp"))) build else "")
				if (!identical(genotypebuild, build)) stop("Custom PGS_PFILE_DIR requires matching PGS_PFILE_BUILD; liftover is not guessed")
				inputfiles <- unlist(lapply(unique(a$CHR), function(ch) pgs_pfile_files(file.path(pd, paste0("chr", ch)))))
				sig <- pgs_hash(list(version = PGS_FOCUS_VERSION, w = a, masks = masks, build = build, pfiles = pgs_stamp(inputfiles), ids = sort(base$eid)))
				cache <- file.path(private_cache, "scores", sig) ; dir.create(cache, recursive = TRUE, showWarnings = FALSE)
				file <- file.path(cache, "scores.rds")
				res <- if (file.exists(file) && !LE8_REPLACE) readRDS(file) else {
					# Cached subset makes scoring multiple source columns cheap. Each
					# feature is isolated so a missing chromosome cannot erase other candidates.
					pp <- pgs_prepare_genotype(a, base$eid, pd, file.path(private_cache, "genotypes"))
					z <- pgs_score_sources(a, masks, pp, cache) ; saveRDS(z, file) ; z
				}
				res$status$feature <- f ; res$status$build <- build ; res$status$MHC_interval <- paste(pgs_shared_mhc(build), collapse = "-")
				counts[[f]] <- res$status ; res$variants$feature <- f ; res$variants$build <- build ; variants[[f]] <- res$variants
				dd <- base ; dd$.m <- as.numeric(biom[[f]][im]) ; j <- match(dd$eid, res$scores$eid)
				good <- !is.na(j) & res$scores$.complete_genotype[j] %in% TRUE
				dd <- dd[good, , drop = FALSE] ; j <- j[good]
				validscopes <- res$status$scope[res$status$matched > 0 & res$status$match_fraction >= pgs_num("PGS_MIN_VARIANT_MATCH", .9)]
				existing <- getdata(f) ; dd$.g <- existing$.g[match(dd$eid, existing$eid)]
				common <- pgs_complete(dd, c(".m", ".g", covars))
				for (n in validscopes) common <- common & is.finite(res$scores[[n]][j])
				dd <- dd[common, , drop = FALSE] ; j <- j[common]
				r0 <- pgs_pair(dd, covars, f, scope = "existing_full_common_genotype")
				summaries[[paste(f, "existing")]] <- r0$summary ; effects[[paste(f, "existing")]] <- r0$effects
				fullcor <- if (sum(is.finite(dd$.g)) > 2) cor(dd$.g, res$scores$full_recomputed[j], use = "complete.obs") else NA_real_
				for (n in validscopes) {
					d <- dd ; d$.g <- res$scores[[n]][j]
					r <- pgs_pair(d, covars, f, scope = n) ; effects[[paste(f, n)]] <- r$effects ; summaries[[paste(f, n)]] <- r$summary
					ca <- pgs_calibrate(d, covars) ; ca$audit$feature <- f ; ca$audit$scope <- n ; calibration[[paste(f, n)]] <- ca$audit
					co <- pgs_component_models(ca$data, covars, f, scope = n) ; components[[paste(f, n)]] <- co$effects ; contrasts[[paste(f, n)]] <- co$contrast
				}
				status[[f]] <- data.frame(
					feature = f, status = if (length(validscopes)) "ok" else "no adequately matched scopes", build = build,
					annotation_status = annotation_status, full_score_correlation = fullcor, N_source_common = nrow(dd), N_genotype_complete = sum(good), N_genotype_incomplete = sum(!good), detail = "Existing PGS untouched; common sample across recomputed scopes. Cis/trans labels inherit C2 annotation; check its build."
				)
			},
			error = function(e) {
				status[[f]] <<- data.frame(feature = f, status = "unavailable", detail = conditionMessage(e))
			}
		)
	}
	list(
		source_status = pgs_bind(status), source_regions = pgs_bind(region_audit), source_counts = pgs_bind(counts), source_variants = pgs_bind(variants),
		source_models = pgs_adjust(pgs_bind(effects), groups = c("scope", "model", "term")), source_paired = pgs_bind(summaries),
		source_components = pgs_adjust(pgs_bind(components), groups = c("scope", "model", "term")), source_calibration = pgs_bind(calibration),
		source_contrasts = pgs_adjust(pgs_bind(contrasts), groups = "scope")
	)
}


# c1.correlate.R
# Aggregate-only figures. Every rendered page has a companion workbook.
pgs_theme <- function() ggplot2::theme_classic(base_size = 11) + ggplot2::theme(
	plot.title = ggplot2::element_text(face = "bold", size = 12), plot.subtitle = ggplot2::element_text(size = 9, color = "#526172"),
	legend.position = "bottom", legend.title = ggplot2::element_blank(), plot.margin = ggplot2::margin(10, 12, 10, 10)
)
pgs_palette <- c(
	"Both supported: opposite" = "#B34359", "Both supported: concordant" = "#287D8E",
	"Measured only" = "#B39557", "PGS only" = "#7964A5", "Neither supported" = "#C1C8CE", "Unavailable comparison" = "#E0E3E5"
)
pgs_ok <- function(d, cols) is.data.frame(d) && nrow(d) > 0 && all(cols %in% names(d))
pgs_pick <- function(x, n = 8) {
	if (!pgs_ok(x$paired, c("feature", "evidence_pattern"))) return(character())
	d <- x$paired ; anchor <- if (pgs_ok(x$selection, c("feature", "selection"))) x$selection$feature[x$selection$selection == "declared anchor"] else character()
	ranked <- d$feature[order(d$evidence_pattern != "Both supported: opposite", pmax(d$measured_FDR, d$pgs_FDR), d$measured_p, na.last = TRUE)]
	# Ranked discordants plus declared biological controls; selection is exploratory.
	unique(head(c(head(ranked, n %/% 2), anchor, ranked), n))
}
pgs_forest_plot <- function(d, title, subtitle = "", unit = "Log HR", label = "feature", color = "series") {
	if (!pgs_ok(d, c(label, "beta", "lo", "hi", color))) return(NULL)
	d <- d[is.finite(d$beta) & is.finite(d$lo) & is.finite(d$hi), , drop = FALSE] ; if (!nrow(d)) return(NULL)
	d$.label <- factor(d[[label]], levels = rev(unique(d[[label]]))) ; d$.series <- d[[color]]
	ggplot2::ggplot(d, ggplot2::aes(beta, .label, color = .series)) +
		ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
		ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = .15, position = ggplot2::position_dodge(.55)) +
		ggplot2::geom_point(position = ggplot2::position_dodge(.55), size = 2) +
		ggplot2::scale_color_manual(values = rep(c("#287D8E", "#B34359", "#7964A5", "#B39557", "#526172", "#59A596", "#A57761", "#555555"), length.out = max(1, length(unique(d$.series))))) +
		ggplot2::labs(title = title, subtitle = subtitle, x = unit, y = NULL, color = NULL) +
		pgs_theme()
}
pgs_locus_panels <- function(features, loci) {
	pp <- list()
	if (pgs_ok(features, c("observed_beta", "pgs_beta", "triangulation_class"))) {
		a <- features[is.finite(features$observed_beta) & is.finite(features$pgs_beta), , drop = FALSE]
		if (nrow(a)) pp$comparison <- ggplot2::ggplot(a, ggplot2::aes(pgs_beta, observed_beta, color = triangulation_class)) +
			ggplot2::geom_hline(yintercept = 0, color = "grey80") +
			ggplot2::geom_vline(xintercept = 0, color = "grey80") +
			ggplot2::geom_point(alpha = .65) +
			ggplot2::scale_color_manual(values = pgs_palette) +
			ggplot2::labs(
				title = "Measured and PGS associations", subtitle = if (all(a$matched_comparison)) "Identical participants and adjustment" else "Contains legacy unmatched screening; rerun pgs_focus",
				x = "PGS log HR / own SD", y = "Measured log HR / own SD", color = NULL
			) +
			pgs_theme()
	}
	if (pgs_ok(loci, c("feature", "locus", "locus_class", "PP.H4", "PP.H4_robust_min", "PP.H3"))) {
		rank <- match(loci$feature, features$feature) ; a <- head(loci[order(rank, - loci$PP.H4, na.last = TRUE), ], 14)
		a$.label <- paste(a$feature, a$locus_class, a$locus, sep = " | ")
		z <- pgs_bind(lapply(c("PP.H3", "PP.H4", "PP.H4_robust_min"), function(n) data.frame(label = a$.label, metric = n, posterior = a[[n]])))
		z$label <- factor(z$label, levels = rev(unique(a$.label)))
		pp$loci <- ggplot2::ggplot(z, ggplot2::aes(metric, label, fill = posterior)) +
			ggplot2::geom_tile(color = "white") +
			ggplot2::geom_text(ggplot2::aes(label = ifelse(is.finite(posterior), sprintf("%.2f", posterior), "NA")), size = 2.7) +
			ggplot2::scale_fill_gradient(low = "#F0F4F5", high = "#247F8C", limits = c(0, 1), na.value = "#DEDEDE") +
			ggplot2::scale_x_discrete(labels = c("PP.H3" = "H3", "PP.H4" = "H4 default", "PP.H4_robust_min" = "H4 min prior")) +
			ggplot2::labs(title = "Locus support and prior sensitivity", subtitle = "All regions retained in workbook; NA means not tested", x = NULL, y = NULL, fill = NULL) +
			pgs_theme() +
			ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7), legend.position = "none")
		z <- unique(loci[, c("region_cluster", "cluster_features", "cluster_robust_features")]) ; z <- head(z[order( - z$cluster_robust_features, - z$cluster_features), ], 10)
		pp$clusters <- ggplot2::ggplot(z, ggplot2::aes(cluster_robust_features, reorder(region_cluster, cluster_robust_features))) +
			ggplot2::geom_col(fill = "#7964A5", width = .7) +
			ggplot2::labs(title = "Shared-region concentration", subtitle = "Overlapping windows, not independent causal signals", x = "Biomarkers with robust H4", y = NULL) +
			pgs_theme()
	}
	pp
}
pgs_main_panels <- function(x, tri = data.frame(), loci = data.frame()) {
	pp <- list() ; chosen <- pgs_pick(x) ; d <- x$paired
	if (pgs_ok(d, c("feature", "pgs_beta", "measured_beta", "evidence_pattern"))) {
		a <- d[is.finite(d$pgs_beta) & is.finite(d$measured_beta), , drop = FALSE] ; a$label <- ifelse(a$feature %in% chosen, a$feature, NA_character_)
		if (nrow(a)) pp$paired <- ggplot2::ggplot(a, ggplot2::aes(pgs_beta, measured_beta, color = evidence_pattern)) +
			ggplot2::geom_hline(yintercept = 0, color = "grey80") +
			ggplot2::geom_vline(xintercept = 0, color = "grey80") +
			ggplot2::geom_point(alpha = .65, size = 1.7) +
			ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 2.8, seed = 2026, max.overlaps = 20, na.rm = TRUE) +
			ggplot2::scale_color_manual(values = pgs_palette) +
			ggplot2::labs(
				title = "Measured level versus biomarker PGS",
				subtitle = "Identical people and covariates; BH FDR within each full scan", x = "PGS log HR / own SD", y = "Measured log HR / own SD", color = NULL
			) +
			pgs_theme()
	}
	a <- x$matched_models
	if (pgs_ok(a, c("feature", "model", "beta", "lo", "hi"))) {
		a <- a[a$feature %in% chosen & a$model %in% c("measured", "pgs"), ] ; a$series <- ifelse(a$model == "measured", "Measured", "Biomarker PGS")
		pp$effects <- pgs_forest_plot(a, "Direction and uncertainty", "Separate models, identical sample", "Log HR per own SD")
	}
	a <- x$components
	if (pgs_ok(a, c("feature", "term", "model", "beta", "lo", "hi"))) {
		a <- a[a$feature %in% chosen & a$model == "joint", ] ; a$series <- ifelse(a$term == ".G", "PGS-captured G", "Remaining R")
		pp$components <- pgs_forest_plot(a, "Captured and remaining components", "Joint model; conditional 95% CI; refitted bootstrap in workbook", "Log HR per whole-biomarker SD")
	}
	a <- x$calibration
	if (pgs_ok(a, c("feature", "partial_R2", "weak_capture"))) {
		a <- a[a$feature %in% chosen & is.finite(a$partial_R2), ]
		if (nrow(a)) pp$calibration <- ggplot2::ggplot(a, ggplot2::aes(partial_R2, reorder(feature, partial_R2), color = weak_capture)) +
			ggplot2::geom_vline(xintercept = 0, color = "grey75") +
			ggplot2::geom_point(size = 2.7) +
			ggplot2::scale_color_manual(values = c("FALSE" = "#287D8E", "TRUE" = "#B39557"), labels = c("Adequate capture", "Weak capture")) +
			ggplot2::labs(title = "Does PGS capture the measured biomarker?", subtitle = "Held-out incremental prediction; negative values retained", x = "Cross-fitted partial R²", y = NULL, color = NULL) +
			pgs_theme()
	}
	a <- x$window_contrasts
	if (pgs_ok(a, c("feature", "beta_difference", "lo", "hi", "landmark", "end"))) {
		a <- a[a$feature %in% head(chosen, 4) & is.finite(a$beta_difference), ]
		if (nrow(a)) {
			a$window <- paste0(a$landmark, "–", ifelse(is.finite(a$end), a$end, "end")) ; a$window <- factor(a$window, levels = unique(a$window[order(a$landmark)]))
			pp$time <- ggplot2::ggplot(a, ggplot2::aes(window, beta_difference, color = feature, group = feature)) +
				ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
				ggplot2::geom_line() +
				ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = .12, position = ggplot2::position_dodge(.25)) +
				ggplot2::geom_point(position = ggplot2::position_dodge(.25)) +
				ggplot2::labs(title = "Does discordance persist with follow-up?", subtitle = "Risk-set intervals; joint contrast includes G/R covariance", x = "Years since baseline", y = "βG − βR (conditional 95% CI)", color = NULL) +
				pgs_theme()
		}
	}
	if (nrow(loci)) {
		# Order by the same declared figure candidates, not a second best-H4 ranking.
		tr <- if (nrow(tri)) tri[order(match(tri$feature, chosen), na.last = TRUE), , drop = FALSE] else data.frame(feature = chosen)
		lp <- pgs_locus_panels(tr, loci) ; pp$loci <- lp$loci
	}
	Filter(Negate(is.null), pp)
}
pgs_save_panels <- function(pp, rd, stem, tables, title, caption = "") {
	# Publication and facet-expanded panels can arrive without list names.
	# Assign keys before filtering/pagination so they retain their input identity.
	keys <- names(pp) ; if (is.null(keys)) keys <- rep("", length(pp))
	missing <- is.na(keys) | !nzchar(keys)
	keys[missing] <- paste0("panel_", which(missing)) ; names(pp) <- keys
	pp <- Filter(Negate(is.null), pp) ; if (!length(pp)) return(invisible(NULL))
	dir.create(rd, recursive = TRUE, showWarnings = FALSE) ; manifest <- list()
	group <- le8_figure_rule(paste0(stem, '.png'))$group
	for (start in seq(1, length(pp), by = 6L)) {
		page <- pp[start : min(length(pp), start + 5L)] ; name <- if (length(pp) <= 6) stem else paste0(stem, ".page", ceiling(start / 6))
		p <- patchwork::wrap_plots(page, ncol = 2) + patchwork::plot_annotation(
			title = title,
			caption = stringr::str_wrap(caption, 140), tag_levels = "A", theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 17), plot.caption = ggplot2::element_text(hjust = 0, size = 9))
		)
		height <- ceiling(length(page) / 2) * 4.7 + 1
		ggplot2::ggsave(file.path(rd, paste0(name, ".png")), p, width = 17, height = height, dpi = as.integer(Sys.getenv("LE8_FINAL_DPI", "300")), bg = "white", limitsize = FALSE)
		le8_table_plot_exports(page, rd, name)
		provenance <- data.frame(panel = LETTERS[seq_along(page)], key = names(page), title = vapply(page, function(p) paste(as.character(p$labels$title), collapse = " "), character(1)))
		pgs_focus_export(c(list(panels = provenance, caption = data.frame(caption = caption)), tables), rd, name)
		manifest[[length(manifest) + 1]] <- data.frame(file = paste0(name, ".png"), group = group, panels = length(page), source = "PGS focus aggregate results")
	}
	mf <- file.path(rd, "figure_manifest.csv") ; old <- if (file.exists(mf)) as.data.frame(data.table::fread(mf)) else data.frame()
	new <- pgs_bind(manifest)
	if (nrow(old) && 'group' %in% names(old)) {
		old$group <- sub('^(c[1-5][.])?Fig[0-9]+[.]', '', old$group)
		stale <- old$file[old$group == group & !old$file %in% new$file]
		old <- old[old$group != group, , drop = FALSE]
		if (length(stale)) unlink(file.path(rd, stale))
	}
	data.table::fwrite(pgs_bind(list(old, new)), mf)
	le8_refresh_figure_files(rd)
	invisible(pp)
}
pgs_plot_focus <- function(x, outdir, tri = data.frame(), loci = data.frame()) {
	caption <- "Exploratory association decomposition. G denotes the part captured by the supplied PGS; R includes uncaptured genetics, lifestyle, disease and measurement. Opposite associations do not establish antagonistic pleiotropy. Conditional CIs do not include calibration uncertainty; anchor bootstrap refits calibration."
	rd <- file.path(outdir, "c1_correlate")
	pgs_save_panels(pgs_main_panels(x, tri, loci), rd, "c1.FigPGS.overview", x, "Measured levels and genetic prediction", caption)
	chosen <- pgs_pick(x, 12) ; extra <- list()
	a <- x$landmarks
	if (pgs_ok(a, c("feature", "model", "term", "landmark", "beta", "lo", "hi"))) {
		a <- a[a$feature %in% head(chosen, 6) & a$model == "joint", ] ; a$series <- paste(a$term, "after", a$landmark, "y")
		extra$landmarks <- pgs_forest_plot(a, "Minimum lead-time sensitivity", unit = "Joint log HR / whole-biomarker SD")
	}
	a <- x$bootstrap
	if (pgs_ok(a, c("feature", "term", "lo", "hi", "status"))) {
		a <- a[a$term == "difference" & a$status == "ok", ] ; b <- x$contrasts
		if (nrow(a) && pgs_ok(b, c("feature", "beta_difference"))) {
			a$beta <- b$beta_difference[match(a$feature, b$feature)] ; a$series <- "Refitted bootstrap" ;
			extra$bootstrap <- pgs_forest_plot(a, "Refitted bootstrap of βG − βR", "Participant / family resampling; selected anchors only", unit = "Joint coefficient difference")
		}
	}
	a <- x$lifestyle
	if (pgs_ok(a, c("feature", "component", "part", "beta"))) {
		a <- a[a$feature %in% head(chosen, 8) & a$part == ".R" & is.finite(a$beta), ]
		if (nrow(a)) extra$lifestyle <- ggplot2::ggplot(a, ggplot2::aes(component, feature, fill = beta)) +
			ggplot2::geom_tile(color = "white") +
			ggplot2::scale_fill_gradient2(low = "#B34359", mid = "white", high = "#287D8E") +
			ggplot2::labs(title = "LE8 connections with remaining component R", subtitle = "Adjusted contemporaneous associations; not intervention effects", x = NULL, y = NULL, fill = NULL) +
			pgs_theme() +
			ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
	}
	a <- x$source_models
	if (pgs_ok(a, c("feature", "scope", "model", "beta", "lo", "hi"))) {
		a <- a[a$model == "pgs", ] ; a$series <- a$scope
		extra$sources <- pgs_forest_plot(a, "Actual rescoring by genetic source", "Common genotype-complete participants; see allele audit", unit = "PGS log HR per own SD")
	}
	a <- x$prevalent
	if (pgs_ok(a, c("feature", "model", "term", "beta", "lo", "hi"))) {
		a <- a[a$feature %in% chosen & a$model == "joint", ] ; a$series <- ifelse(a$term == ".G", "PGS-captured G", "Remaining R")
		extra$prevalent <- pgs_forest_plot(a, "Baseline prevalent cases: G and R", "Descriptive case-control analysis; separate from incident hazards", unit = "Log OR / whole-biomarker SD")
	}
	pgs_save_panels(extra, rd, "c1.FigPGS.sensitivity", x, "Discordance sensitivity analyses", caption)
}
pgs_publication_joint <- function(trait, analysis_root, layers = c("prot", "met")) {
	if (!all(c("prot", "met") %in% layers)) return(invisible(NULL))
	root <- file.path(analysis_root, trait) ; objects <- list() ; tri <- list() ; loci <- list() ; tables <- list() ; pp <- list()
	for (layer in c("prot", "met")) {
		f <- file.path(root, layer, "c1_correlate", "c1.pgs_focus.rds")
		if (!file.exists(f)) return(invisible(NULL)) ; x <- readRDS(f)
		if (!pgs_ok(x$paired, c("feature", "evidence_pattern"))) return(invisible(NULL))
		cf <- file.path(root, layer, "c3_coloc", "c3.coloc_summary.csv")
		co <- if (file.exists(cf)) as.data.frame(data.table::fread(cf)) else data.frame()
		tr <- make_c3_pgs_integration(coloc_summary = co, focus = x) ; lo <- attr(tr, "loci")
		parts <- pgs_main_panels(x, tr, lo)
		if (!is.null(parts$paired)) pp[[paste0(layer, "_paired")]] <- parts$paired + ggplot2::labs(title = paste(toupper(layer), "| Measured versus PGS"))
		if (!is.null(parts$components)) pp[[paste0(layer, "_components")]] <- parts$components + ggplot2::labs(title = paste(toupper(layer), "| Captured and remaining"))
		for (n in names(x)) if (is.data.frame(x[[n]])) tables[[paste0(layer, "_", n)]] <- x[[n]]
		tables[[paste0(layer, "_loci")]] <- lo
		pick <- pgs_pick(x, 6) ; cc <- x$calibration
		if (pgs_ok(cc, c("feature", "partial_R2"))) {
			cc <- cc[cc$feature %in% pick, ] ; cc$feature <- paste(toupper(layer), cc$feature, sep = ": ") ; objects[[layer]] <- cc
		}
		tr$feature <- paste(toupper(layer), tr$feature, sep = ": ") ; tri[[layer]] <- tr
		if (nrow(lo)) {
			lo$feature <- paste(toupper(layer), lo$feature, sep = ": ") ; loci[[layer]] <- lo
		}
	}
	cc <- pgs_bind(objects)
	if (nrow(cc)) {
		# Use all selected cross-layer calibration rows, not a new discovery ranking.
		cc <- cc[is.finite(cc$partial_R2), ]
		if (nrow(cc)) pp$calibration <- ggplot2::ggplot(cc, ggplot2::aes(partial_R2, reorder(feature, partial_R2), color = weak_capture)) +
			ggplot2::geom_vline(xintercept = 0, color = "grey75") +
			ggplot2::geom_point(size = 2.5) +
			ggplot2::scale_color_manual(values = c("FALSE" = "#287D8E", "TRUE" = "#B39557")) +
			ggplot2::labs(title = "PGS capture across both omic layers", x = "Cross-fitted partial R²", y = NULL, color = "Weak capture") +
			pgs_theme()
	}
	loc <- pgs_locus_panels(pgs_bind(tri), pgs_bind(loci)) ; if (!is.null(loc$loci)) pp$loci <- loc$loci
	if (length(pp)) pgs_save_panels(
		pp, le8_final_dir(trait = trait, root = analysis_root), "Fig1.PGS_integrated", tables, paste(trait, "| Measured and genetically predicted biomarkers"),
		"Exploratory matched-cohort comparisons within each omic layer. G and R share the whole-biomarker scale; R is not pure lifestyle. The protein and metabolite cohorts need not contain the same people. Colocalization remains locus-specific. See each layer's Fig1 for follow-up contrasts and its workbook for refitted bootstrap and source-score sensitivity."
	)
	invisible(pp)
}


find_c1_pgs_file <- function(layer) {
	explicit <- Sys.getenv("C1_PGS_FILE", unset = "")
	automatic <- file.path(indir, "Rdata", if (layer == "protein") "prot.pgs.rds" else "met.pgs.rds")
	hit <- unique(c(explicit, automatic))
	hit <- hit[nzchar(hit) & file.exists(hit) & file.size(hit) > 0]
	if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = FALSE) else NA_character_
}

c1_pgs_signature <- function(layer) {
	f <- find_c1_pgs_file(layer)
	if (is.na(f)) return("missing")
	i <- file.info(f)
	paste(normalizePath(f, winslash = "/", mustWork = FALSE), i$size,
		format(i$mtime, "%Y-%m-%dT%H:%M:%S%z"),
		sep = "|"
	)
}

read_c1_pgs <- function(file) {
	x <- if (grepl("\\.rds$", file, ignore.case = TRUE)) readRDS(file) else
		data.table::fread(file, showProgress = FALSE, check.names = FALSE)
	x <- as_tibble(x)
	if (!"eid" %in% names(x)) stop("PGS file must contain eid: ", file, call. = FALSE)
	# Keep this conversion: UKB joins otherwise fail when one input stores eid
	# as integer64 and another stores it as character.
	pgs_ids(x, "Biomarker PGS input")
}

map_c1_pgs_columns <- function(features, nms) {
	one <- function(f) {
		candidates <- c(
			paste0(f, ".pgs"), paste0(f, "_pgs"), paste0(f, ".PGS"),
			paste0(f, "_PGS"), paste0(f, "_GRS"), paste0("GRS_", f), f
		)
		hit <- candidates[candidates %in% nms]
		if (length(hit)) hit[[1]] else NA_character_
	}
	ans <- setNames(vapply(features, one, character(1)), features)
	ans[!is.na(ans)]
}

empty_c1_pgs_assoc <- function(features, score_columns = character()) {
	tibble(
		term = features, score_column = unname(score_columns[features]),
		estimate = NA_real_, beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
		conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_,
		N_total = NA_integer_, N_event = NA_integer_, FDR = NA_real_
	)
}

run_c1_pgs_scan <- function(
	layer, features, covars, outcome = Y, rawdir = NULL,
	overlap_eids = NULL
) {
	score_file <- find_c1_pgs_file(layer)
	overlap_eids <- as.character(overlap_eids %||% character())
	overlap_signature <- if (!length(overlap_eids)) "none" else le8_hash_object(sort(unique(overlap_eids)))
	signature <- pgs_hash(list(C1_PGS_SCAN_VERSION, c1_pgs_signature(layer), overlap_signature,
		outcome = outcome, covars = covars, features = features, baseline = LE8_BASELINE_VERSION,
		options = le8_analysis_options(), phenotypes = pgs_stamp(file.path(indir, "Rdata/all.rds"))
	))
	cache <- if (is.null(rawdir)) NA_character_ else file.path(rawdir, "c1.pgs_scan.rds")
	if (!is.na(cache) && cache_valid(cache)) {
		old <- tryCatch(readRDS(cache), error = function(e) NULL)
		if (is.list(old) && identical(old$signature, signature) &&
			all(c("status", "incident", "prevalent", "attained_age") %in% names(old))) {
			message("C1/", layer, ": reuse inherited-omic PGS scans")
			return(old)
		}
	}
	empty <- list(
		signature = signature, score_file = score_file,
		status = tibble(status = "unavailable", detail = if (is.na(score_file))
			"prot.pgs.rds/met.pgs.rds was not found" else "No assayed feature matched a PGS column"),
		incident = empty_c1_pgs_assoc(features), prevalent = empty_c1_pgs_assoc(features),
		attained_age = empty_c1_pgs_assoc(features),
		incident_same_omic = empty_c1_pgs_assoc(features),
		prevalent_same_omic = empty_c1_pgs_assoc(features),
		attained_age_same_omic = empty_c1_pgs_assoc(features),
		vldl_conditional = tibble(), score_map = character()
	)
	if (is.na(score_file)) {
		warning("C1/", layer, ": inherited-score file is unavailable; PGS panels will be blank.", call. = FALSE)
		return(empty)
	}
	scores <- read_c1_pgs(score_file)
	score_map <- map_c1_pgs_columns(features, names(scores))
	max_scores <- suppressWarnings(as.integer(Sys.getenv("C1_PGS_MAX", unset = "0")))
	if (is.finite(max_scores) && max_scores > 0) score_map <- head(score_map, max_scores)
	if (!length(score_map)) return(empty)

	need <- unique(c(
		"eid", "ethnic.c", covars, "birth_date", "date_attend", "date_lost",
		"date_death", paste0("fod_icd10_", outcome)
	))
	ph <- read_all(need) |>
		filter_analysis_cohort() |>
		make_outcome(outcome) |>
		add_attained_age_time(outcome)
	ph$prevalent_status <- make_prevalent_status(ph, outcome)
	# Do not remove this conversion; see read_c1_pgs().
	ph$eid <- as.character(ph$eid)
	idx <- match(scores$eid, ph$eid)
	keep <- !is.na(idx)
	scores <- scores[keep, unique(c("eid", unname(score_map))), drop = FALSE]
	base <- ph[idx[keep], , drop = FALSE]
	base$.omic_overlap <- if (length(overlap_eids)) base$eid %in% overlap_eids else TRUE
	rm(ph) ; invisible(gc())
	tvar <- paste0(outcome, ".t2e") ; evar <- paste0(outcome, ".Yt2e")
	covars <- intersect(covars, names(base))

	message(
		"C1/", layer, ": scan ", length(score_map),
		" inherited omic scores in the genotyped cohort (sequential; full PGS matrix remains in memory)"
	)
	rows <- lapply(names(score_map), function(feature) {
		gcol <- score_map[[feature]]
		d <- base
		d$.pgs_score <- suppressWarnings(as.numeric(scores[[gcol]]))
		inc <- cox_scan(d, ".pgs_score", covars, outcome, time_var = tvar, event_var = evar) |>
			mutate(term = feature, score_column = gcol)
		prev <- logistic_scan(d, ".pgs_score", covars, "prevalent_status") |>
			mutate(term = feature, score_column = gcol)
		age <- cox_scan_delayed_entry(d, ".pgs_score", setdiff(
			covars,
			unique(c("age", grep("^age($|[._])", covars, value = TRUE, ignore.case = TRUE)))
		), outcome) |>
			mutate(term = feature, score_column = gcol)
		ds <- d[d$.omic_overlap %in% TRUE, , drop = FALSE]
		inc_same <- cox_scan(ds, ".pgs_score", covars, outcome, time_var = tvar, event_var = evar) |>
			mutate(term = feature, score_column = gcol)
		prev_same <- logistic_scan(ds, ".pgs_score", covars, "prevalent_status") |>
			mutate(term = feature, score_column = gcol)
		age_same <- cox_scan_delayed_entry(ds, ".pgs_score", setdiff(
			covars,
			unique(c("age", grep("^age($|[._])", covars, value = TRUE, ignore.case = TRUE)))
		), outcome) |>
			mutate(term = feature, score_column = gcol)
		list(
			incident = inc, prevalent = prev, attained_age = age,
			incident_same_omic = inc_same, prevalent_same_omic = prev_same,
			attained_age_same_omic = age_same
		)
	})
	finish <- function(kind) bind_rows(lapply(rows, `[[`, kind)) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		arrange(p.value)
	vldl_conditional <- tibble()
	if (layer == "metabolite" && "L_VLDL_TG.pct" %in% names(score_map)) {
		comparator_features <- intersect(
			c("L_VLDL_TG", "Total_TG", "ApoB", "VLDL_size"),
			names(score_map)
		)
		model_sets <- c(
			list(`Target PGS only` = character()),
			setNames(lapply(comparator_features, c), paste0("+ ", comparator_features))
		)
		if (length(comparator_features) > 1L)
			model_sets[["+ all burden/size PGS"]] <- comparator_features
		dd <- base
		target_name <- ".pgs_L_VLDL_TG_pct"
		dd[[target_name]] <- suppressWarnings(as.numeric(scores[[score_map[["L_VLDL_TG.pct"]]]]))
		zmap <- setNames(character(length(comparator_features)), comparator_features)
		for (feature in comparator_features) {
			zn <- paste0(".pgs_", make.names(feature)) ; zmap[[feature]] <- zn
			dd[[zn]] <- suppressWarnings(as.numeric(scores[[score_map[[feature]]]]))
		}
		same_n_need <- unique(c(tvar, evar, covars, target_name, unname(zmap)))
		dd <- dd[complete.cases(dd[, same_n_need, drop = FALSE]), , drop = FALSE]
		for (nm in c(target_name, unname(zmap))) dd[[nm]] <- as.numeric(scale(dd[[nm]]))
		vldl_conditional <- map_dfr(names(model_sets), function(model_name) {
			adjust_features <- model_sets[[model_name]] ; xs <- c(target_name, unname(zmap[adjust_features]))
			empty_row <- tibble(
				model = model_name, target = "L_VLDL_TG.pct PGS", adjusters =
					paste(adjust_features, collapse = ";"), beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
				conf.high = NA_real_, p.value = NA_real_, N_total = nrow(dd), N_event = sum(dd[[evar]] == 1),
				max_abs_score_correlation = NA_real_, exposure_condition_number = NA_real_
			)
			if (nrow(dd) < 500 || sum(dd[[evar]] == 1) < 20) return(empty_row)
			ff <- as.formula(paste0(
				"Surv(", bt(tvar), ",", bt(evar), ") ~ ",
				paste(bt(c(xs, covars)), collapse = " + ")
			))
			fit <- tryCatch(coxph(ff, dd, ties = "efron"), error = function(e) NULL)
			if (is.null(fit)) return(empty_row)
			sm <- coef(summary(fit)) ; if (!target_name %in% rownames(sm)) return(empty_row)
			b <- sm[target_name, "coef"] ; se <- sm[target_name, "se(coef)"]
			cm <- if (length(xs) > 1L) cor(dd[, xs, drop = FALSE]) else matrix(1, 1, 1)
			max_cor <- if (length(xs) > 1L) max(abs(cm[1, - 1]), na.rm = TRUE) else 0
			cond <- if (length(xs) > 1L) tryCatch(kappa(cm), error = function(e) NA_real_) else 1
			tibble(
				model = model_name, target = "L_VLDL_TG.pct PGS", adjusters = paste(adjust_features, collapse = ";"),
				beta = b, std.error = se, conf.low = b - 1.96 * se, conf.high = b + 1.96 * se,
				p.value = sm[target_name, "Pr(>|z|)"], N_total = nrow(dd), N_event = sum(dd[[evar]] == 1),
				max_abs_score_correlation = max_cor, exposure_condition_number = cond
			)
		}) |> mutate(
			FDR = p.adjust(p.value, "BH"),
			interpretation = "Exploratory multivariable PGS association; not multivariable MR"
		)
	}
	ans <- list(
		signature = signature, score_file = score_file,
		status = tibble(
			status = "ok", detail = paste(length(score_map), "matched PGS columns"),
			score_file = score_file, genotype_rows = nrow(scores), phenotype_matches = nrow(base),
			same_omic_rows = sum(base$.omic_overlap),
			adjustment = paste(covars, collapse = ";")
		),
		incident = finish("incident"), prevalent = finish("prevalent"),
		attained_age = finish("attained_age"),
		incident_same_omic = finish("incident_same_omic"),
		prevalent_same_omic = finish("prevalent_same_omic"),
		attained_age_same_omic = finish("attained_age_same_omic"),
		vldl_conditional = vldl_conditional, score_map = score_map
	)
	if (!is.na(cache)) saveRDS(ans, cache, compress = "xz")
	ans
}

build_pgs_actual_concordance <- function(pgs_full, observed, pgs_same = NULL) {
	analyses <- c(
		incident = "Incident", prevalent = "Baseline prevalent",
		attained_age = "Attained-age delayed entry"
	)
	map_dfr(names(analyses), function(nm) {
		g <- as_tibble(pgs_full[[nm]] %||% tibble())
		s <- as_tibble((pgs_same %||% list())[[nm]] %||% tibble())
		o <- as_tibble(observed[[nm]] %||% tibble())
		if (!nrow(g) && !nrow(s) && !nrow(o)) return(tibble())
		if (nrow(g) && !"score_column" %in% names(g)) g$score_column <- NA_character_
		if (nrow(s) && !"score_column" %in% names(s)) s$score_column <- NA_character_
		gg <- if (nrow(g)) g |> transmute(
			feature = term, score_column,
			pgs_full_beta = beta, pgs_full_se = std.error, pgs_full_p = p.value,
			pgs_full_FDR = FDR, pgs_full_N = N_total, pgs_full_events = N_event
		) else
			tibble(
				feature = character(), score_column = character(), pgs_full_beta = double(),
				pgs_full_se = double(), pgs_full_p = double(), pgs_full_FDR = double(),
				pgs_full_N = double(), pgs_full_events = double()
			)
		ss <- if (nrow(s)) s |> transmute(
			feature = term,
			pgs_same_beta = beta, pgs_same_se = std.error, pgs_same_p = p.value,
			pgs_same_FDR = FDR, pgs_same_N = N_total, pgs_same_events = N_event
		) else
			tibble(
				feature = character(), pgs_same_beta = double(), pgs_same_se = double(),
				pgs_same_p = double(), pgs_same_FDR = double(), pgs_same_N = double(),
				pgs_same_events = double()
			)
		oo <- if (nrow(o)) o |> transmute(
			feature = term, observed_beta = beta,
			observed_se = std.error, observed_p = p.value, observed_FDR = FDR,
			observed_N = N_total, observed_events = N_event
		) else
			tibble(
				feature = character(), observed_beta = double(), observed_se = double(),
				observed_p = double(), observed_FDR = double(), observed_N = double(),
				observed_events = double()
			)
		full_join(gg, ss, by = "feature") |>
			full_join(oo, by = "feature") |>
			mutate(
				analysis = analyses[[nm]],
				# Compatibility aliases now refer explicitly to the full genetic cohort.
				pgs_beta = pgs_full_beta, pgs_p = pgs_full_p, pgs_FDR = pgs_full_FDR,
				pgs_N = pgs_full_N, pgs_events = pgs_full_events,
				pgs_supported = is.finite(pgs_full_FDR) & pgs_full_FDR < .05,
				observed_supported = is.finite(observed_FDR) & observed_FDR < .05,
				sign_concordant = is.finite(pgs_full_beta) & is.finite(observed_beta) &
					sign(pgs_full_beta) == sign(observed_beta),
				same_omic_sign_concordant = is.finite(pgs_same_beta) & is.finite(observed_beta) &
					sign(pgs_same_beta) == sign(observed_beta),
				beta_difference_full = observed_beta - pgs_full_beta,
				beta_difference_same_omic = observed_beta - pgs_same_beta,
				evidence_pattern = case_when(
					pgs_supported & observed_supported & sign_concordant ~ "Both, same direction",
					pgs_supported & observed_supported & !sign_concordant ~ "Both, opposite direction",
					pgs_supported & !observed_supported ~ "PGS only",
					!pgs_supported & observed_supported ~ "Observed only",
					TRUE ~ "Neither at FDR 5%"
				)
			)
	})
}

plot_pgs_actual_concordance <- function(x, anchors = character()) {
	d <- as_tibble(x)
	fmt_n <- function(v) {
		n <- sort(unique(as.integer(v[is.finite(v)])))
		if (!length(n)) return("NA")
		z <- if (length(n) == 1L) n else range(n)
		paste(format(z, big.mark = ",", scientific = FALSE, trim = TRUE), collapse = "–")
	}
	one_panel <- function(xvar, nvar, panel_title, xlab) {
		z <- d |> filter(is.finite(.data[[xvar]]), is.finite(observed_beta))
		if (!nrow(z)) return(blank_plot(panel_title, "No matched PGS–observed effect pair was estimable"))
		facets <- z |>
			group_by(analysis) |>
			summarise(
				facet = paste0(
					first(analysis), "\nPGS N = ", fmt_n(.data[[nvar]]),
					"; measured N = ", fmt_n(observed_N)
				), .groups = "drop"
			)
		z <- z |>
			left_join(facets, by = "analysis") |>
			group_by(analysis) |>
			mutate(
				.discord = abs(observed_beta - .data[[xvar]]),
				label = ifelse(feature %in% anchors | min_rank(desc(.discord)) <= 6, feature, NA_character_),
				direction = ifelse(sign(observed_beta) == sign(.data[[xvar]]),
					"Same direction", "Opposite direction"
				)
			) |>
			ungroup()
		st <- z |>
			group_by(facet) |>
			summarise(
				rho = suppressWarnings(cor(.data[[xvar]], observed_beta, method = "spearman", use = "complete.obs")),
				stat_label = ifelse(is.finite(rho), sprintf("Spearman rho = %.2f", rho), "Spearman rho = NA"),
				.groups = "drop"
			)
		ggplot(z, aes(x = .data[[xvar]], y = observed_beta, color = direction)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey58") +
			geom_hline(yintercept = 0, color = "grey82") +
			geom_vline(xintercept = 0, color = "grey82") +
			geom_point(alpha = .62, size = 1.65) +
			geom_text(
				data = st, aes(x =  - Inf, y = Inf, label = stat_label), inherit.aes = FALSE,
				hjust =  - .04, vjust = 1.25, size = 2.5, fontface = "bold", color = "grey30"
			) +
			ggrepel::geom_text_repel(aes(label = label),
				size = 2.25, max.overlaps = 24,
				seed = 91, show.legend = FALSE
			) +
			facet_wrap( ~ facet, scales = "free", nrow = 1) +
			scale_color_manual(values = c(
				"Same direction" = "#1B9E77",
				"Opposite direction" = "#D7301F"
			)) +
			labs(
				title = panel_title,
				subtitle = "Effects, not P values, are compared; both predictors are standardized within their own model",
				x = xlab, y = "Measured adult omic: log HR/OR per 1 SD", color = NULL
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")
	}
	pa <- one_panel(
		"pgs_full_beta", "pgs_full_N",
		"a. Full genetic cohort: inherited versus measured effect size",
		"Inherited omic PGS: log HR/OR per 1 SD"
	)
	pb <- one_panel(
		"pgs_same_beta", "pgs_same_N",
		"b. Same-omic-subset sensitivity",
		"Inherited omic PGS in assay subset: log HR/OR per 1 SD"
	)
	pa / pb + plot_layout(heights = c(1, 1)) +
		plot_annotation(caption = paste0(
			"PGS is fixed at conception but is not a biomarker measured at birth. ",
			"Effect concordance is descriptive and is not an MR or mediation estimate."
		))
}

# Cross-module PGS triangulation and matched-PGS associations
# Locus-resolved triangulation. A genome-wide PGS is never an input to coloc.
# Keep every tested region; overlapping regions are not independent discoveries.
pgs_coloc_loci <- function(x, h4 = .70) {
	d <- as.data.frame(x) ; if (!nrow(d) || !all(c("feature", "PP.H4") %in% names(d))) return(data.frame())
	numeric_cols <- c("chr", "start", "end", "PP.H3", "PP.H4_robust_min", "credible_set_n")
	for (n in numeric_cols) if (!n %in% names(d)) d[[n]] <- NA_real_
	for (n in c("locus", "locus_class", "status", "lead_shared")) if (!n %in% names(d)) d[[n]] <- NA_character_
	d$chr <- sub("^chr", "", as.character(d$chr), ignore.case = TRUE)
	d$prior_tested <- is.finite(d$PP.H4_robust_min)
	d$robust_coloc <- d$status %in% "ok" & d$prior_tested & d$PP.H4_robust_min >= h4
	d$locus_evidence <- ifelse(!d$status %in% "ok", "Unavailable / failed",
		ifelse(d$robust_coloc, "H4 robust across tested priors", ifelse(d$PP.H4 >= h4 & !d$prior_tested,
			"H4 default; prior sensitivity unavailable", ifelse(d$PP.H4 >= h4, "H4 sensitive to priors",
				ifelse(is.finite(d$PP.H3) & d$PP.H3 >= h4, "H3: distinct signals favored", "No strong H3/H4 support")
			)
		))
	)
	d$region_cluster <- NA_character_
	for (ch in unique(d$chr[!is.na(d$chr)])) {
		ix <- which(d$chr == ch & is.finite(d$start) & is.finite(d$end) & d$start <= d$end)
		ix <- ix[order(d$start[ix], d$end[ix])] ; right <-  - Inf ; k <- 0L
		for (i in ix) {
			if (d$start[i] > right) {
				k <- k + 1L ; right <- d$end[i]
			} else right <- max(right, d$end[i]) ; d$region_cluster[i] <- paste0("chr", ch, "_region", k)
		}
	}
	# Unknown coordinates are not merged based only on a reused lead label.
	missing <- which(is.na(d$region_cluster)) ; d$region_cluster[missing] <- paste0("unmapped_", missing)
	d$cluster_features <- vapply(d$region_cluster, function(k) length(unique(d$feature[d$region_cluster == k])), integer(1))
	d$cluster_robust_features <- vapply(d$region_cluster, function(k) length(unique(d$feature[d$region_cluster == k & d$robust_coloc])), integer(1))
	d$interpretation <- "Overlapping windows define a region cluster, not one fine-mapped causal signal; H4 is locus-specific"
	d
}
make_c3_pgs_integration <- function(c1 = list(), c2 = list(), coloc_summary = data.frame(), h4 = .70, focus = NULL) {
	loci <- pgs_coloc_loci(coloc_summary, h4)
	paired <- if (is.list(focus)) focus$paired else NULL
	if (is.data.frame(paired) && nrow(paired)) {
		z <- as.data.frame(paired) ; names(z)[names(z) == "measured_beta"] <- "observed_beta"
		names(z)[names(z) == "measured_p"] <- "observed_p" ; names(z)[names(z) == "measured_FDR"] <- "observed_FDR"
		z$matched_comparison <- TRUE
	} else {
		a <- c1$association_adj2 ; if (is.null(a)) a <- c1$association
		p <- c1$pgs_incident
		conv <- function(d, pgs = FALSE) {
			if (is.null(d) || !nrow(d)) return(data.frame(feature = character()))
			pre <- if (pgs) "pgs" else "observed" ; o <- data.frame(feature = d$term)
			o[[paste0(pre, "_beta")]] <- d$beta ; o[[paste0(pre, "_p")]] <- d$p.value ; o[[paste0(pre, "_FDR")]] <- d$FDR ; o
		}
		z <- merge(conv(a), conv(p, TRUE), by = "feature", all = TRUE) ; z$matched_comparison <- rep(FALSE, nrow(z))
	}
	features <- unique(c(z$feature, loci$feature))
	if (!length(features)) {
		z <- data.frame(feature = character()) ; attr(z, "loci") <- loci ; return(z)
	}
	z <- merge(data.frame(feature = features), z, by = "feature", all.x = TRUE)
	if (!nrow(z)) {
		attr(z, "loci") <- loci ; return(z)
	}
	for (n in c("observed_beta", "observed_p", "observed_FDR", "pgs_beta", "pgs_p", "pgs_FDR")) if (!n %in% names(z)) z[[n]] <- NA_real_
	if (!"matched_comparison" %in% names(z)) z$matched_comparison <- FALSE
	z$matched_comparison[is.na(z$matched_comparison)] <- FALSE
	z$observed_supported <- is.finite(z$observed_FDR) & z$observed_FDR < .05
	z$pgs_supported <- is.finite(z$pgs_FDR) & z$pgs_FDR < .05
	z$pgs_observed_sign_match <- is.finite(z$observed_beta) & is.finite(z$pgs_beta) & z$observed_beta * z$pgs_beta >= 0
	z$triangulation_class <- pgs_classify(z$observed_beta, z$pgs_beta, z$observed_FDR, z$pgs_FDR)
	z$n_loci_tested <- z$n_robust_loci <- z$n_region_clusters <- 0L ; z$robust_PP4 <- NA_real_ ; z$coloc_supported <- FALSE
	z$coloc_summary <- "Not tested / unavailable"
	for (i in seq_len(nrow(z))) {
		a <- loci[loci$feature == z$feature[i], , drop = FALSE] ; if (!nrow(a)) next
		z$n_loci_tested[i] <- sum(a$status %in% "ok") ; z$n_robust_loci[i] <- sum(a$robust_coloc)
		z$n_region_clusters[i] <- length(unique(a$region_cluster[a$robust_coloc]))
		z$robust_PP4[i] <- if (any(is.finite(a$PP.H4_robust_min))) max(a$PP.H4_robust_min, na.rm = TRUE) else NA_real_
		z$coloc_supported[i] <- any(a$robust_coloc)
		z$coloc_summary[i] <- paste(unique(paste(a$locus_class, a$locus_evidence, sep = ": ")), collapse = "; ")
	}
	z$inherited_locus_pattern <- z$observed_supported & z$pgs_supported & z$pgs_observed_sign_match & z$coloc_supported
	z$reactive_compatible_pattern <- FALSE # Absence of genetic support is not evidence of consequence.
	z$interpretation <- ifelse(z$matched_comparison,
		"Identical-sample associations; opposite directions do not establish antagonistic pleiotropy or causal partition",
		"Legacy different-sample screening only; run pgs_focus before comparing direction or statistical strength"
	)
	if (is.list(focus) && is.data.frame(focus$calibration) && nrow(focus$calibration)) z <- merge(z, focus$calibration, by = "feature", all.x = TRUE)
	z <- z[order(!(z$triangulation_class == "Both supported: opposite"), z$pgs_p, z$observed_p), , drop = FALSE]
	attr(z, "loci") <- loci ; z
}
read_c3_pgs_integration <- function(layer, outdir, coloc_summary) {
	get <- function(module, file) {
		p <- file.path(outdir, module, file) ; if (file.exists(p)) tryCatch(readRDS(p), error = function(e) list()) else list()
	}
	x <- make_c3_pgs_integration(
		get("c1_correlate", "c1.res.rds"), get("c2_cause", "c2.res.rds"),
		coloc_summary, get0("H4_STRONG", ifnotfound = .70), get("c1_correlate", "c1.pgs_focus.rds")
	)
	rd <- file.path(outdir, "c3_coloc") ; dir.create(rd, recursive = TRUE, showWarnings = FALSE)
	loci <- attr(x, "loci") ; clusters <- if (nrow(loci)) unique(loci[, c("region_cluster", "chr", "cluster_features", "cluster_robust_features")]) else data.frame()
	pgs_focus_export(list(features = x, loci = loci, region_clusters = clusters), rd, "c3.pgs_focus")
	# Refresh the established CSV name as well; downstream readers must not see
	# an old "PGS weak" label after a focused rerun.
	data.table::fwrite(as.data.frame(x), file.path(rd, "c3.pgs_observed_coloc_triangulation.csv"))
	x
}
plot_c3_pgs_integration <- function(x, outdir) {
	if (!nrow(x)) return(invisible(NULL))
	rd <- file.path(outdir, "c3_coloc") ; dir.create(rd, recursive = TRUE, showWarnings = FALSE)
	pp <- pgs_locus_panels(x, attr(x, "loci"))
	if (length(pp)) pgs_save_panels(
		pp, rd, "c3.Fig5.pgs_coloc_triangulation",
		list(features = as.data.frame(x), loci = attr(x, "loci")), "Measured / PGS discordance and locus evidence"
	)
	invisible(x)
}

# Matched biomarker PGS -> measured omic. This is calibration/association,
# not a disease PRS -> all proteins scan and not identification of mediation.
le8_pgs_bridge <- function(dat, features, fold, covars, disease, membership, layer, rawdir) {
	empty <- list(scan = tibble(), bridges = tibble(), status = tibble(status = "PGS input unavailable"))
	emit <- function(out) {
		for (nm in names(out)) write_raw_csv(out[[nm]], paste0("c4.matched_PGS_", nm, ".csv"), rawdir)
		out
	}
	f <- Sys.getenv("C4_PGS_FILE", unset = find_c1_pgs_file(layer))
	if (is.na(f) || !nzchar(f) || !file.exists(f)) return(emit(empty))
	message("C4: load matched omic PGS: ", f)
	scores <- read_c1_pgs(f)
	if (anyDuplicated(scores$eid) || anyDuplicated(dat$eid)) stop("Duplicate eid in omic/PGS input")
	mapping <- map_c1_pgs_columns(features, setdiff(names(scores), "eid"))
	if (!length(mapping)) return(emit(empty))
	idx <- match(as.character(dat$eid), scores$eid)
	scores <- scores[idx, unname(mapping), drop = FALSE]
	covars <- intersect(covars, names(dat)) ; rows <- list()
	for (h in unique(fold)) for (feature in names(mapping)) {
		y <- suppressWarnings(as.numeric(dat[[feature]])) ; g <- as.numeric(scores[[mapping[[feature]]]])
		keep <- fold == h & is.finite(y) & is.finite(g) & complete.cases(dat[, covars, drop = FALSE])
		n <- sum(keep) ; r <- p <- NA_real_
		if (n >= 100) {
			# Subset covariates only: copying thousands of omic columns for every
			# QR fit made a genome-wide bridge needlessly slow and memory intensive.
			M <- model.matrix(reformulate(covars), dat[keep, covars, drop = FALSE]) ; q <- qr(M)
			yr <- qr.resid(q, y[keep]) ; gr <- qr.resid(q, g[keep]) ; den <- sqrt(sum(yr ^ 2) * sum(gr ^ 2)) ; df <- n - q$rank - 1
			if (is.finite(den) && den > 0 && df > 2) {
				r <- max( - .999999, min(.999999, sum(yr * gr) / den)) ; p <- 2 * pt(abs(r) * sqrt(df / (1 - r * r)), df, lower.tail = FALSE)
			}
		}
		rows[[length(rows) + 1L]] <- tibble(feature, score_column = mapping[[feature]], split = h, n, r, p.value = p)
		if (length(rows) %% 500L == 0L) message("C4 matched PGS: ", length(rows), " / ", length(unique(fold)) * length(mapping), " fits")
	}
	rm(scores) ; invisible(gc())
	scan <- bind_rows(rows) |>
		group_by(split) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		ungroup()
	a <- scan |>
		filter(split == "discovery") |>
		select(feature, score_column, n_disc = n, r_disc = r, p_disc = p.value, FDR_disc = FDR)
	b <- scan |>
		filter(split == "replication") |>
		select(feature, n_rep = n, r_rep = r, p_rep = p.value, FDR_rep = FDR)
	z <- left_join(a, b, by = "feature") |> mutate(
		replicated = coalesce(FDR_disc < .05 & FDR_rep < .05 & r_disc * r_rep > 0, FALSE),
		partial_R2_rep = r_rep ^ 2, edge = "biomarker PGS -> measured biomarker",
		causal_mediation_identified = FALSE, source_file = f
	)
	if (nrow(membership)) z <- left_join(z, membership |> select(any_of(c("feature", "primary_component", "set", "strict_YS"))), by = "feature")
	if (nrow(disease)) z <- left_join(z, disease |> transmute(feature = term, observed_disease_p = p.value), by = "feature")
	c1file <- file.path(dirname(rawdir), "c1_correlate", "c1.res.rds")
	if (file.exists(c1file)) {
		c1 <- readRDS(c1file) ; pg <- c1$pgs_incident %||% tibble()
		if (nrow(pg)) z <- left_join(z, pg |> transmute(feature = term, PGS_disease_beta = beta, PGS_disease_p = p.value, PGS_disease_FDR = FDR), by = "feature")
	}
	status <- tibble(
		status = "ok", matched_scores = length(mapping), replicated_edges = sum(z$replicated),
		note = "Split-replicated partial association; discovery GWAS overlap and score transportability require separate assessment"
	)
	for (nm in c("scan", "bridges", "status")) {
		value <- switch(nm,
			scan = scan,
			bridges = z,
			status = status
		)
		write_raw_csv(value, paste0("c4.matched_PGS_", nm, ".csv"), rawdir)
	}
	list(scan = scan, bridges = z, status = status)
}


if (.c1_mode == "pgs") {
	source(file.path(Sys.getenv("LE8_FDIR", "."), "0.common.R"))
for (layer in c(if (prot_DO) "protein", if (met_DO) "metabolite")) {
		outdir <- if (layer == "protein")
			out.prot else out.met
		x <- run_c1_pgs_focus(layer, outdir)
		cf <- file.path(outdir, "c3_coloc", "c3.coloc_summary.csv")
		co <- if (file.exists(cf))
			as.data.frame(data.table::fread(cf, showProgress = FALSE)) else data.frame()
		tri <- read_c3_pgs_integration(layer, outdir, co)
		pgs_plot_focus(x, outdir, tri, attr(tri, "loci"))
		plot_c3_pgs_integration(tri, outdir)
		rm(x, tri)
		invisible(gc())
	}
} else if (.c1_mode == "correlate") {
# C1: prospective/prevalent associations, trajectories, clusters, and enrichment.

suppressPackageStartupMessages({
	fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(fdir, "0.common.R"))
})
# Endpoint ascertainment is separate from biomarker lead time. In particular,
# zero early AMR events must not be interpreted as biological latency.
le8_endpoint_audit <- function(dat, features, covars, Y, root) {
	tv <- paste0(Y, ".t2e") ; ev <- paste0(Y, ".Yt2e") ; bv <- paste0(Y, ".b2e")
	cuts <- c(0, .5, 1, 2, 5, 10, Inf)
	windows <- bind_rows(lapply(seq_len(length(cuts) - 1L), function(i) {
		at_risk <- is.finite(dat[[tv]]) & dat[[tv]] > cuts[i]
		cases <- at_risk & dat[[ev]] == 1 & dat[[tv]] <= cuts[i + 1]
		data.frame(
			window_lo = cuts[i], window_hi = cuts[i + 1], N_at_risk = sum(at_risk),
			events = sum(cases, na.rm = TRUE), interpretation = "Observed event distribution; verify recording coverage"
		)
	}))
	write.csv(windows, file.path(root, "c1.endpoint_event_windows.csv"), row.names = FALSE)
	calendar <- if ("date_attend" %in% names(dat)) as.Date(dat$date_attend) else rep(as.Date(NA), nrow(dat))
	dy <- paste0("fod_icd10_", Y)
	if (dy %in% names(dat)) {
		dates <- as.Date(dat[[dy]]) ; z <- table(format(dates[!is.na(dates)], "%Y"))
		write.csv(data.frame(year = names(z), recorded_diagnoses = as.integer(z)), file.path(root, "c1.endpoint_calendar_years.csv"), row.names = FALSE)
	}
	mf <- Sys.getenv("LE8_ENDPOINT_MANIFEST", "")
	out <- data.frame(Y,
		prevalent = sum(is.finite(dat[[bv]]) & dat[[bv]] <= 0),
		early_2y_events = sum(dat[[ev]] == 1 & dat[[tv]] <= 2, na.rm = TRUE),
		registry_status = "unknown; absence of a diagnosis record is not confirmed absence of infection/resistance"
	)
	if (nzchar(mf)) {
		m <- read.csv(mf, stringsAsFactors = FALSE)
		required <- c("Y", "source", "definition", "coverage_start", "coverage_end", "prebaseline_capture")
		if (!all(required %in% names(m))) stop("Endpoint manifest missing fields: ", paste(setdiff(required, names(m)), collapse = ","))
		m <- m[m$Y == Y, , drop = FALSE]
		if (nrow(m) == 1) {
			write.csv(m, file.path(root, "c1.endpoint_manifest_used.csv"), row.names = FALSE)
			start <- as.Date(m$coverage_start) ; end <- as.Date(m$coverage_end)
			if (is.na(start) || is.na(end) || end <= start) stop("Registry dates must be valid, verified dates")
			d <- dat ; d$.registry_entry <- pmax(0, as.numeric(start - calendar) / 365.25)
			stop_time <- as.numeric(end - calendar) / 365.25
			d$.registry_event <- as.integer(d[[ev]] == 1 & d[[tv]] <= stop_time)
			d$.registry_exit <- pmin(d[[tv]], stop_time)
			# Explicitly exclude any already observed disease at/before entry.
			eligible <- is.finite(d$.registry_entry) & is.finite(d$.registry_exit) & d$.registry_exit > d$.registry_entry
			d <- d[eligible, , drop = FALSE]
			z <- cox_scan_delayed_entry(d, features, covars, Y, entry_var = ".registry_entry", exit_var = ".registry_exit", event_var = ".registry_event")
			z$time_scale <- "Years after blood draw, delayed entry at verified registry start"
			z$estimand <- "First recorded Y during covered follow-up; no claim of first lifetime disease if prebaseline capture incomplete"
			write.csv(z, file.path(root, "c1.registry_delayed_entry.csv"), row.names = FALSE)
			out$registry_status <- paste("Sensitivity completed;", m$source, "; prebaseline capture:", m$prebaseline_capture)
		} else out$registry_status <- "No unique manifest row for this Y; no registry correction assumed"
	}
	write.csv(out, file.path(root, "c1.endpoint_audit.csv"), row.names = FALSE)
	out
}

LE8_JOB <- "c1_correlate"
C1_CODE_VERSION <- "2026-09-15.final-systematic"
C1_SCAN_VERSION <- "2026-09-15.final-systematic"

TOP_N <- as.integer(Sys.getenv("C1_TOP_N", unset = "30"))
YY_TOP <- as.integer(Sys.getenv("C1_YY_TOP", unset = "6"))
GRADIENT_TOP <- as.integer(Sys.getenv("C1_GRADIENT_TOP", unset = "10"))
CLUSTER_TOP <- as.integer(Sys.getenv("C1_CLUSTER_TOP", unset = "100"))
ATTENUATION_TOP <- as.integer(Sys.getenv("C1_ATTENUATION_TOP", unset = "500"))
YY_BINS <- as.integer(Sys.getenv("C1_YY_BINS", unset = "28"))
YY_MAX_YEAR <- as.numeric(Sys.getenv("C1_YY_MAX_YEAR", unset = "16"))
GRADIENT_STEP <- as.numeric(Sys.getenv("C1_GRADIENT_STEP", unset = "1"))
CLUSTER_STEP <- as.numeric(Sys.getenv("C1_CLUSTER_STEP", unset = "0.5"))
MIN_BIN_N <- as.integer(Sys.getenv("C1_MIN_BIN_N", unset = "20"))
CLUSTER_K_MAX <- min(4L, as.integer(Sys.getenv("C1_CLUSTER_K_MAX", unset = "4")))
CLUSTER_GAP_B <- as.integer(Sys.getenv("C1_CLUSTER_GAP_B", unset = "30"))
CLUSTER_STABILITY_B <- as.integer(Sys.getenv("C1_CLUSTER_STABILITY_B", unset = "40"))
C1_MIN_EVENT <- 20L
C1_YY_SPAR <- as.numeric(Sys.getenv("C1_YY_SPAR", unset = "0.62"))
C1_DIRECTION_ANCHORS <- unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
	unset = "PCSK9,LPA,GDF15,NTPROBNP,MMP12,L_VLDL_TG.pct,L_VLDL_TG,Total_TG,ApoB"
), ",", fixed = TRUE)[[1]]))
C1_VLDL_DEEP_FEATURES <- c(
	paste0(c("XS", "S", "M", "L", "XL", "XXL"), "_VLDL_TG.pct"),
	paste0(c("XS", "S", "M", "L", "XL", "XXL"), "_VLDL_TG"),
	"L_VLDL_L", "L_VLDL_CE.pct", "L_VLDL_FC.pct", "L_VLDL_PL.pct",
	"Total_TG", "ApoB", "VLDL_size"
)
C1_LANDMARK_YEARS <- sort(unique(as.numeric(strsplit(Sys.getenv("C1_LANDMARK_YEARS", unset = "0.5,1,2,5,10"), ",", fixed = TRUE)[[1]])))
C1_LANDMARK_YEARS <- C1_LANDMARK_YEARS[is.finite(C1_LANDMARK_YEARS) & C1_LANDMARK_YEARS > 0]
C1_LANDMARK_TOP <- as.integer(Sys.getenv("C1_LANDMARK_TOP", unset = "500"))
C1_RISK_TOP <- as.integer(Sys.getenv("C1_RISK_TOP", unset = "16"))
C1_RISK_CUTS <- sort(unique(as.numeric(strsplit(Sys.getenv("C1_RISK_CUTS",
	unset = "0,0.5,1,2,5,10,16"
), ",", fixed = TRUE)[[1]])))
C1_RISK_CUTS <- C1_RISK_CUTS[is.finite(C1_RISK_CUTS) & C1_RISK_CUTS >= 0]
if (length(C1_RISK_CUTS) < 2L) C1_RISK_CUTS <- c(0, .5, 1, 2, 5, 10, 16)
split_env_names <- function(name, unset = "") {
	z <- trimws(strsplit(Sys.getenv(name, unset = unset), ",", fixed = TRUE)[[1]])
	unique(z[nzchar(z)])
}
# The primary omics model adjusts for the four behavioral LE8 components only.
# BMI, non-HDL-C, HbA1c and blood pressure can be biological intermediates or
# overlap the assayed biomarker; conditioning on them can overadjust the omic
# effect. The former full-LE8 model is retained as a sensitivity analysis.
C1_LE4_COVARS <- split_env_names(
	"C1_LE4_COVARS",
	"diet.pts,pa.pts,smoke.pts,sleep.pts"
)
C1_FULL_LE8_SENSITIVITY <- truthy(Sys.getenv("C1_FULL_LE8_SENSITIVITY",
	unset = Sys.getenv("C1_MET_FULL_LE8_SENSITIVITY", unset = "TRUE")
))
# Optional known baseline treatment variables, supplied by the local phenotype
# dictionary (for example a lipid-lowering-medication indicator).  No variable
# name is guessed silently.
C1_TREATMENT_VARS <- split_env_names("C1_TREATMENT_VARS")

# Publication typography. This overrides the shared theme only inside C1.
theme_5c <- function(base_size = 12) {
	theme_classic(base_size = base_size) +
		theme(
			plot.title = element_text(face = "bold", size = base_size * 1.14, hjust = 0),
			plot.subtitle = element_text(face = "bold", size = base_size * .94, color = "grey30"),
			axis.title = element_text(face = "bold", size = base_size * 1.08),
			axis.text = element_text(face = "bold", size = base_size, color = "black"),
			legend.title = element_text(face = "bold"),
			legend.text = element_text(face = "bold"),
			strip.background = element_blank(),
			strip.text = element_text(face = "bold", size = base_size * 1.02),
			panel.grid.major.y = element_line(color = "grey91", linewidth = .25),
			panel.grid.minor = element_blank(),
			plot.margin = margin(9, 13, 9, 13)
		)
}
forest_theme <- function(base_size = 10) theme_5c(base_size) + theme(panel.grid.major.y = element_blank())

assoc_has_results <- function(x) is.data.frame(x) && nrow(x) > 0 && any(is.finite(x$p.value))

assoc_empty_message <- function(x, case_label, min_event = C1_MIN_EVENT) {
	nt <- suppressWarnings(max(x$N_total, na.rm = TRUE)) ; ne <- suppressWarnings(max(x$N_event, na.rm = TRUE))
	if (!is.finite(nt)) nt <- NA_integer_ ; if (!is.finite(ne)) ne <- NA_integer_
	paste0(
		"Complete-case N = ", format(nt, big.mark = ","), "\n",
		case_label, " cases = ", format(ne, big.mark = ","),
		if (is.finite(ne) && ne < min_event) paste0(" (<", min_event, " required)") else ""
	)
}

assoc_blank_plot <- function(title, message) {
	ggplot() +
		annotate("text", x = 0, y = .16, label = title, fontface = "bold", size = 7) +
		annotate("text", x = 0, y =  - .12, label = message, fontface = "bold", color = "grey25", size = 5) +
		xlim( - 1, 1) +
		ylim( - 1, 1) +
		theme_void(base_size = 16)
}


# 🚩 Baseline prevalent association: logistic regression

logistic_scan <- function(dat, xs, covars, y, scale_x = TRUE, min_n = 500, min_case = 20) {
	xs <- intersect(xs, names(dat)) ; covars <- intersect(covars, names(dat))
	bind_rows(parallel_map(xs, function(x) {
		d <- dat[, unique(c(y, x, covars)), drop = FALSE]
		d <- d[complete.cases(d), , drop = FALSE]
		nc <- sum(d[[y]] == 1, na.rm = TRUE)
		empty <- tibble(
			term = x, estimate = NA_real_, beta = NA_real_, std.error = NA_real_,
			conf.low = NA_real_, conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_,
			N_total = nrow(d), N_event = nc
		)
		if (nrow(d) < min_n || nc < min_case || length(unique(d[[y]])) < 2) return(empty)
		d[[x]] <- suppressWarnings(as.numeric(d[[x]])) ; sx <- sd(d[[x]], na.rm = TRUE)
		if (!is.finite(sx) || sx <= 0) return(empty)
		if (scale_x) d[[x]] <- as.numeric(scale(d[[x]]))
		f <- reformulate(c(x, covars), response = y)
		fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
		if (is.null(fit)) return(empty)
		sm <- coef(summary(fit)) ; if (!x %in% rownames(sm)) return(empty)
		b <- sm[x, "Estimate"] ; se <- sm[x, "Std. Error"]
		tibble(
			term = x, estimate = exp(b), beta = b, std.error = se,
			conf.low = exp(b - 1.96 * se), conf.high = exp(b + 1.96 * se),
			statistic = b / se, p.value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
			N_total = nrow(d), N_event = nc
		)
	})) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		arrange(p.value)
}

# Among baseline-prevalent cases, estimate whether current protein abundance
# varies with elapsed time since diagnosis.  This explicitly tests the
# disease -> protein interpretation that a time-agnostic case/control model
# cannot address.
prevalent_duration_scan <- function(dat, xs, bvar, covars, min_n = 80L) {
	xs <- intersect(xs, names(dat)) ; covars <- intersect(covars, names(dat))
	bind_rows(parallel_map(xs, function(x) {
		d <- dat[, unique(c(bvar, x, covars)), drop = FALSE] |>
			filter(is.finite(.data[[bvar]]), .data[[bvar]] < 0) |>
			transmute(.protein = suppressWarnings(as.numeric(.data[[x]])), .duration = log1p( - .data[[bvar]]), across(all_of(covars)))
		d <- d[complete.cases(d), , drop = FALSE]
		empty <- tibble(term = x, beta = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = nrow(d))
		if (nrow(d) < min_n || sd(d$.protein) <= 0 || sd(d$.duration) <= 0) return(empty)
		d$.protein <- as.numeric(scale(d$.protein)) ; d$.duration <- as.numeric(scale(d$.duration))
		fit <- tryCatch(lm(reformulate(c(".duration", covars), response = ".protein"), d), error = function(e) NULL)
		if (is.null(fit)) return(empty) ; sm <- coef(summary(fit)) ; if (!".duration" %in% rownames(sm)) return(empty)
		tibble(
			term = x, beta = sm[".duration", "Estimate"], std.error = sm[".duration", "Std. Error"],
			statistic = sm[".duration", "t value"], p.value = sm[".duration", "Pr(>|t|)"], N_total = nrow(d)
		)
	})) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		arrange(p.value)
}


# Formal diagnosis-window associations. Incident cases in [lo, hi) are
# compared with people known to remain disease-free and observed through hi;
# later cases therefore serve as valid controls for earlier windows. Baseline-
# prevalent duration bands use all non-prevalent participants as the reference.

plot_risk_window_scan <- function(z, anchors = C1_DIRECTION_ANCHORS) {
	if (!nrow(z) || !any(is.finite(z$beta))) return(blank_plot(
		"Risk-set diagnosis-window associations",
		"No diagnosis window had enough cases and eligible controls"
	))
	chosen <- unique(c(anchors, z |> filter(is.finite(p.value)) |> group_by(term) |>
		summarise(best_p = min(p.value), .groups = "drop") |> arrange(best_p) |> slice_head(n = C1_RISK_TOP) |> pull(term)))
	d <- z |>
		filter(term %in% chosen, is.finite(beta), is.finite(std.error)) |>
		mutate(term = factor(term, levels = rev(chosen)), window = sprintf("%g–%g y", window_lo, window_hi))
	ggplot(d, aes(time, beta, color = side, fill = side, group = side)) +
		geom_hline(yintercept = 0, color = "grey72") +
		geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
		geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .10, color = NA) +
		geom_line(linewidth = .70) +
		geom_point(aes(size = N_case), alpha = .88) +
		facet_wrap( ~ term, scales = "free_y", ncol = 4) +
		scale_x_continuous(limits = c( - 16, 16), breaks = c( - 16, - 10, - 5, - 2, - 1, 0, 1, 2, 5, 10, 16)) +
		scale_color_manual(values = c(`Pre-baseline prevalent` = "#8C6BB1", `Post-baseline incident` = "#2C7FB8")) +
		scale_fill_manual(values = c(`Pre-baseline prevalent` = "#8C6BB1", `Post-baseline incident` = "#2C7FB8")) +
		labs(
			title = "Risk-set associations across recorded diagnosis time",
			subtitle = "Incident controls remain observed and disease-free through each window; coefficients are adjusted log odds ratios per 1-SD biomarker",
			x = "Years relative to baseline blood draw", y = "Adjusted beta (95% CI)", color = NULL, fill = NULL, size = "Cases"
		) +
		theme_5c(8) +
		theme(legend.position = "bottom")
}

build_directionality_table <- function(incident, prevalent, duration, landmark = tibble(), birthline = tibble(), reverse = tibble()) {
	slim <- function(z, prefix) {
		if (!nrow(z) || !all(c("term", "beta", "p.value", "FDR") %in% names(z)))
			z <- tibble(term = character(), beta = numeric(), p.value = numeric(), FDR = numeric())
		z |>
			select(term, beta, p.value, FDR) |>
			rename(! !paste0("beta_", prefix) := beta, ! !paste0("p_", prefix) := p.value, ! !paste0("FDR_", prefix) := FDR)
	}
	# A 2-year result must never silently become the column labelled 5-year.
	landmark_target <- if (nrow(landmark) && any(landmark$landmark_years == 5, na.rm = TRUE)) 5 else NA_real_
	lm5 <- if (nrow(landmark) && is.finite(landmark_target)) landmark |>
		filter(landmark_years == landmark_target) |>
		select(term, beta_landmark5 = beta, p_landmark5 = p.value, FDR_landmark5 = FDR) else
		tibble(term = character(), beta_landmark5 = numeric(), p_landmark5 = numeric(), FDR_landmark5 = numeric())
	full_join(slim(incident, "incident"), slim(prevalent, "prevalent"), by = "term") |>
		full_join(slim(duration, "duration"), by = "term") |>
		full_join(lm5, by = "term") |>
		full_join(slim(birthline, "birthline"), by = "term") |>
		full_join(slim(reverse, "reverse"), by = "term") |>
		mutate(
			incident_support = is.finite(FDR_incident) & FDR_incident < .05,
			prevalent_support = is.finite(FDR_prevalent) & FDR_prevalent < .05,
			duration_support = is.finite(FDR_duration) & FDR_duration < .05,
			distal_support = is.finite(FDR_landmark5) & FDR_landmark5 < .05,
			birthline_support = is.finite(FDR_birthline) & FDR_birthline < .05,
			landmark5_tested = is.finite(FDR_landmark5),
			# Missing or nonsignificant distal tests do not establish attenuation.
			# Duration association is a disease-state clue, not reverse causation.
			reactive_compatible = prevalent_support & duration_support,
			direction_class = case_when(
				distal_support & prevalent_support ~ "distal + disease-state (mixed)",
				distal_support & !prevalent_support ~ "distal antecedent supported",
				reactive_compatible ~ "near-diagnosis / reactive-compatible",
				!incident_support & prevalent_support ~ "established-disease associated",
				incident_support & !landmark5_tested ~ "incident association; distal test unavailable",
				incident_support ~ "incident association only",
				TRUE ~ "unresolved"
			),
			incident_score = sign(beta_incident) * pmin(12, - log10(pmax(p_incident, 1e-300))),
			prevalent_score = sign(beta_prevalent) * pmin(12, - log10(pmax(p_prevalent, 1e-300))),
			duration_score = sign(beta_duration) * pmin(12, - log10(pmax(p_duration, 1e-300))),
			landmark5_score = sign(beta_landmark5) * pmin(12, - log10(pmax(p_landmark5, 1e-300))),
			birthline_score = sign(beta_birthline) * pmin(12, - log10(pmax(p_birthline, 1e-300))),
			reverse_score = sign(beta_reverse) * pmin(12, - log10(pmax(p_reverse, 1e-300)))
		) |>
		arrange(match(direction_class, c(
			"near-diagnosis / reactive-compatible", "distal + disease-state (mixed)",
			"distal antecedent supported", "incident association only", "established-disease associated", "unresolved"
		)), p_incident)
}

plot_directionality_triage <- function(z, anchors = C1_DIRECTION_ANCHORS) {
	if (!nrow(z)) return(blank_plot("Temporal directionality triage", "No estimable protein results"))
	pal <- c(
		"near-diagnosis / reactive-compatible" = "#B45A4D", "distal + disease-state (mixed)" = "#8C6BB1",
		"distal antecedent supported" = "#2B8CBE", "incident association only" = "#74A9CF",
		"established-disease associated" = "#B38A3E", "unresolved" = "grey75"
	)
	lab <- z |>
		filter(term %in% anchors) |>
		bind_rows(z |> filter(reactive_compatible) |> slice_min(p_incident, n = 6)) |>
		distinct(term, .keep_all = TRUE)
	pa <- ggplot(z, aes(beta_incident, beta_prevalent, color = direction_class)) +
		geom_hline(yintercept = 0, color = "grey85") +
		geom_vline(xintercept = 0, color = "grey85") +
		geom_point(alpha = .65, size = 1.8) +
		geom_text_repel(data = lab, aes(label = term), size = 3, fontface = "bold", max.overlaps = Inf) +
		scale_color_manual(values = pal) +
		labs(
			title = "a. Incident versus baseline-prevalent association",
			subtitle = "Prevalent association is disease-state compatible, not proof that disease caused the biomarker",
			x = "Incident Cox beta", y = "Prevalent logistic beta", color = NULL
		) +
		theme_5c(9)
	pbdat <- z |> filter(is.finite(beta_landmark5), is.finite(beta_duration))
	pb <- if (!nrow(pbdat)) blank_plot("b. Distal versus established-disease evidence", "Landmark or case-only duration result was unavailable") else
		ggplot(pbdat, aes(beta_landmark5, beta_duration, color = direction_class)) +
			geom_hline(yintercept = 0, color = "grey85") +
			geom_vline(xintercept = 0, color = "grey85") +
			geom_point(alpha = .7, size = 1.8) +
			geom_text_repel(data = pbdat |> filter(term %in% anchors), aes(label = term), size = 3.2, fontface = "bold", max.overlaps = Inf) +
			scale_color_manual(values = pal) +
			labs(
				title = "b. Five-year landmark versus case-only duration",
				subtitle = "Far-horizon persistence argues against a purely near-diagnosis signal; duration remains survivor/treatment sensitive",
				x = "Incident Cox beta after 5-year landmark", y = "Years-since-diagnosis slope", color = NULL
			) +
			theme_5c(9)
	(pa | pb) + plot_layout(guides = "collect") &
		theme(legend.position = "bottom", legend.box = "vertical")
}

plot_directionality_supplement <- function(z, anchors = C1_DIRECTION_ANCHORS) {
	if (!nrow(z)) return(blank_plot("Direction-of-time evidence detail", "No estimable biomarker results"))
	pal <- c(
		"near-diagnosis / reactive-compatible" = "#B45A4D", "distal + disease-state (mixed)" = "#8C6BB1",
		"distal antecedent supported" = "#2B8CBE", "incident association only" = "#74A9CF",
		"established-disease associated" = "#B38A3E", "unresolved" = "grey75"
	)
	chosen <- unique(c(anchors, z |> filter(reactive_compatible | distal_support) |> slice_min(p_incident, n = 14) |> pull(term)))
	heat <- z |>
		filter(term %in% chosen) |>
		select(term, incident_score, birthline_score, landmark5_score, prevalent_score, duration_score) |>
		pivot_longer( - term, names_to = "analysis", values_to = "signed_evidence") |>
		mutate(analysis = factor(analysis,
			levels = c("incident_score", "birthline_score", "landmark5_score", "prevalent_score", "duration_score"),
			labels = c("Incident Cox", "Attained-age Cox", "5-y landmark", "Prevalent logistic", "Duration slope")
		))
	pc <- ggplot(heat, aes(analysis, factor(term, levels = rev(chosen)), fill = signed_evidence)) +
		geom_tile(color = "white") +
		scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c( - 12, 12), name = "signed\n-log10(P)") +
		labs(title = "a. Signed evidence by temporal model", x = NULL, y = NULL) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1))
	pd <- z |>
		count(direction_class, name = "biomarkers") |>
		mutate(direction_class = factor(direction_class, levels = names(pal))) |>
		ggplot(aes(biomarkers, direction_class, fill = direction_class)) +
		geom_col(width = .68) +
		geom_text(aes(label = biomarkers), hjust =  - .12, fontface = "bold") +
		scale_x_continuous(expand = expansion(mult = c(0, .18))) +
		scale_fill_manual(values = pal, guide = "none") +
		labs(title = "b. Full-screen evidence classes", x = "Biomarkers", y = NULL) +
		theme_5c(9)
	pc | pd
}

plot_landmark_and_attained_age <- function(landmark, incident, birthline, anchors = C1_DIRECTION_ANCHORS) {
	top_landmark <- if (nrow(landmark)) landmark |>
		filter(is.finite(p.value)) |>
		group_by(term) |>
		summarise(best_p = min(p.value), .groups = "drop") |>
		slice_min(best_p, n = 8, with_ties = FALSE) |>
		pull(term) else character()
	chosen <- unique(c(anchors, top_landmark)) ; lm <- landmark |>
		filter(term %in% chosen, is.finite(beta), is.finite(std.error)) |>
		mutate(term = factor(term, levels = rev(chosen)), lo = beta - 1.96 * std.error, hi = beta + 1.96 * std.error)
	pa <- if (!nrow(lm)) blank_plot("a. Landmark persistence", "No landmark estimate was available") else
		ggplot(lm, aes(landmark_years, beta, color = term, group = term)) +
			geom_hline(yintercept = 0, linetype = 3, color = "grey65") +
			geom_ribbon(aes(ymin = lo, ymax = hi, fill = term), alpha = .08, color = NA) +
			geom_line(linewidth = .85) +
			geom_point(size = 1.8) +
			scale_x_continuous(breaks = C1_LANDMARK_YEARS) +
			labs(
				title = "a. Landmark persistence",
				subtitle = "Risk sets begin 0.5, 1, 2, 5 or 10 years after baseline; estimates remain prospective",
				x = "Event-free landmark after baseline (years)", y = "Log HR per 1-SD baseline biomarker", color = NULL, fill = NULL
			) +
			theme_5c(9) +
			theme(legend.position = "bottom")
	z <- incident |>
		select(term, beta_baseline = beta, p_baseline = p.value) |>
		inner_join(birthline |> select(term, beta_attained = beta, p_attained = p.value), by = "term") |>
		filter(is.finite(beta_baseline), is.finite(beta_attained))
	rr <- if (nrow(z) > 2) suppressWarnings(cor(z$beta_baseline, z$beta_attained, use = "complete.obs")) else NA_real_
	lab <- z |>
		filter(term %in% anchors) |>
		bind_rows(z |> mutate(delta = abs(beta_attained - beta_baseline)) |> slice_max(delta, n = 6)) |>
		distinct(term, .keep_all = TRUE)
	pb <- if (!nrow(z)) blank_plot("b. Time-scale sensitivity", "Attained-age result was unavailable") else
		ggplot(z, aes(beta_baseline, beta_attained)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
			geom_hline(yintercept = 0, color = "grey88") +
			geom_vline(xintercept = 0, color = "grey88") +
			geom_point(alpha = .55, color = "#3F78A8") +
			ggrepel::geom_text_repel(data = lab, aes(label = term), size = 2.8, fontface = "bold", max.overlaps = Inf, seed = 91) +
			annotate("text", x = Inf, y =  - Inf, label = sprintf("r = %.3f", rr), hjust = 1.08, vjust =  - .5, fontface = "bold", size = 3.2) +
			labs(
				title = "b. Baseline-time versus attained-age Cox",
				subtitle = "Attained-age model uses delayed entry at blood draw; this is not a birth-cohort omics measurement",
				x = "Baseline-time Cox beta", y = "Attained-age delayed-entry Cox beta"
			) +
			theme_5c(9)
	pa | pb
}

plot_reverse_time_exploratory <- function(reverse, prevalent, duration, anchors = C1_DIRECTION_ANCHORS) {
	z <- reverse |>
		select(term, beta_reverse = beta, p_reverse = p.value) |>
		left_join(prevalent |> select(term, beta_prevalent = beta), by = "term") |>
		left_join(duration |> select(term, beta_duration = beta), by = "term")
	lab <- z |>
		filter(term %in% anchors) |>
		bind_rows(z |> filter(is.finite(p_reverse)) |> slice_min(p_reverse, n = 8)) |>
		distinct(term, .keep_all = TRUE)
	one <- function(y, title, ylab) {
		d <- z |> filter(is.finite(beta_reverse), is.finite(.data[[y]])) ; ll <- lab |> semi_join(d, by = "term")
		if (!nrow(d)) return(blank_plot(title, "No exploratory reverse-time estimate was available"))
		ggplot(d, aes(beta_reverse, .data[[y]])) +
			geom_hline(yintercept = 0, color = "grey88") +
			geom_vline(xintercept = 0, color = "grey88") +
			geom_point(alpha = .58, color = "#8C6BB1") +
			ggrepel::geom_text_repel(data = ll, aes(label = term), size = 2.8, fontface = "bold", max.overlaps = Inf, seed = 92) +
			labs(title = title, x = "Legacy reverse-time Cox beta", y = ylab) +
			theme_5c(9)
	}
	(one("beta_prevalent", "a. Reverse-time versus prevalence", "Prevalent logistic beta") |
		one("beta_duration", "b. Reverse-time versus case duration", "Years-since-diagnosis slope")) +
		plot_annotation(
			title = "Exploratory only: incompatible reverse-time risk origins",
			subtitle = "Cases use time since diagnosis, controls use age since birth; these panels are excluded from directionality grading"
		)
}

# Same-N attenuation: both Cox models use the complete adj2 sample for each biomarker.
same_sample_attenuation <- function(dat, features, tvar, evar, covs_basic, covs_adj2) {
	min_basic_beta <- as.numeric(Sys.getenv("C1_ATTENUATION_MIN_ABS_BETA", unset = "0.02"))
	if (!is.finite(min_basic_beta) || min_basic_beta <= 0) min_basic_beta <- 0.02
	features <- intersect(features, names(dat))
	bind_rows(parallel_map(features, function(x) {
		need <- unique(c(tvar, evar, x, covs_adj2)) ; d <- dat[, need, drop = FALSE]
		d <- d[complete.cases(d), , drop = FALSE]
		ne <- sum(d[[evar]] == 1, na.rm = TRUE)
		if (nrow(d) < 500 || ne < C1_MIN_EVENT) return(tibble(
			term = x, N = nrow(d), events = ne,
			beta_basic_sameN = NA_real_, p_basic_sameN = NA_real_, beta_adj2_sameN = NA_real_, p_adj2_sameN = NA_real_,
			effect_ratio = NA_real_, log2_effect_ratio = NA_real_, sign_flip = NA, attenuation_pct = NA_real_
		))
		d[[x]] <- as.numeric(scale(suppressWarnings(as.numeric(d[[x]]))))
		fit_one <- function(covars) {
			f <- as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(x, covars)), collapse = " + ")))
			fit <- tryCatch(coxph(f, d, ties = "efron"), error = function(e) NULL)
			if (is.null(fit)) return(c(beta = NA_real_, p = NA_real_))
			sm <- coef(summary(fit)) ; if (!x %in% rownames(sm)) return(c(beta = NA_real_, p = NA_real_))
			c(beta = sm[x, "coef"], p = sm[x, "Pr(>|z|)"])
		}
		a <- fit_one(covs_basic) ; b <- fit_one(covs_adj2)
		stable <- is.finite(a[["beta"]]) && abs(a[["beta"]]) >= min_basic_beta
		ratio <- if (stable) abs(b[["beta"]]) / abs(a[["beta"]]) else NA_real_
		tibble(
			term = x, N = nrow(d), events = ne, beta_basic_sameN = a[["beta"]], p_basic_sameN = a[["p"]],
			beta_adj2_sameN = b[["beta"]], p_adj2_sameN = b[["p"]],
			effect_ratio = ratio, log2_effect_ratio = ifelse(is.finite(ratio) & ratio > 0, log2(ratio), NA_real_),
			sign_flip = ifelse(is.finite(a[["beta"]]) & is.finite(b[["beta"]]), sign(a[["beta"]]) != sign(b[["beta"]]), NA),
			attenuation_pct = ifelse(is.finite(ratio), 100 * (1 - ratio), NA_real_)
		)
	})) |> arrange(p_basic_sameN)
}


# 🚩 Yin-Yang: baseline = 0, prevalent < 0, incident > 0
residualize_z <- function(d, x, covars) {
	covars <- intersect(covars, names(d)) ; z <- suppressWarnings(as.numeric(d[[x]]))
	if (!length(covars)) return(as.numeric(scale(z)))
	dd <- d[, c(x, covars), drop = FALSE] ; dd[[x]] <- z
	f <- reformulate(covars, response = x)
	fit <- tryCatch(lm(f, dd, na.action = na.exclude), error = function(e) NULL)
	if (is.null(fit)) return(as.numeric(scale(z)))
	as.numeric(scale(residuals(fit)))
}

global_trajectory_p <- function(time, z, df = 3) {
	ok <- is.finite(time) & is.finite(z)
	if (sum(ok) < 80 || length(unique(time[ok])) < 8) return(NA_real_)
	d <- data.frame(time = time[ok], z = z[ok])
	f0 <- tryCatch(lm(z ~ 1, d), error = function(e) NULL)
	f1 <- tryCatch(lm(z ~ splines::ns(time, df = df), d), error = function(e) NULL)
	if (is.null(f0) || is.null(f1)) return(NA_real_)
	a <- tryCatch(anova(f0, f1), error = function(e) NULL)
	if (is.null(a) || nrow(a) < 2) NA_real_ else as.numeric(a$`Pr(>F)`[2])
}

make_yy_panel <- function(dat, features, bvar, covars, title, ylab = "Mean biomarker z-score",
				bins = YY_BINS, max_year = YY_MAX_YEAR, min_bin_n = MIN_BIN_N,
				smooth_lines = TRUE) {
	features <- intersect(features, names(dat)) ; if (!length(features)) return(blank_plot(title, "No feature available"))
	d0 <- dat |> filter(is.finite(.data[[bvar]]))
	lim <- min(max_year, max(abs(d0[[bvar]]), na.rm = TRUE)) ; if (!is.finite(lim) || lim <= 0) lim <- max_year
	br <- seq( - lim, lim, length.out = bins + 1) ; mids <- (head(br, - 1) + tail(br, - 1)) / 2
	hist0 <- hist(d0[[bvar]], breaks = br, plot = FALSE)
	h <- tibble(
		mid = mids, xmin = head(br, - 1), xmax = tail(br, - 1), n = hist0$counts,
		side = ifelse(mids < 0, "Prevalent (Yang)", "Incident (Yin)")
	)
	rows <- map_dfr(features, function(x) {
		need <- unique(c(bvar, x, covars)) ; dd <- d0[, need, drop = FALSE] ; dd <- dd[complete.cases(dd), , drop = FALSE]
		if (nrow(dd) < 80) return(tibble())
		dd$z <- residualize_z(dd, x, covars)
		dd$bin <- cut(dd[[bvar]], br, include.lowest = TRUE, labels = FALSE)
		pg <- global_trajectory_p(dd[[bvar]], dd$z)
		sm <- dd |>
			filter(is.finite(z), !is.na(bin)) |>
			group_by(bin) |>
			summarise(mean = mean(z), sd = sd(z), se = sd / sqrt(n()), N = n(), .groups = "drop") |>
			filter(N >= min_bin_n) |>
			mutate(time = mids[bin], feature = x, p_global = pg, N_total = nrow(dd))
		# Smooth only the observed bin summaries.  A weighted smoothing spline is
		# deterministic, does not extrapolate, and prevents small extreme bins from
		# dominating the visual trajectory.  Raw bin means remain in the output.
		if (nrow(sm) >= 5) {
			fit <- tryCatch(stats::smooth.spline(sm$time, sm$mean, w = sqrt(sm$N), spar = C1_YY_SPAR),
				error = function(e) NULL
			)
			sm$mean_smooth <- if (is.null(fit)) sm$mean else
				as.numeric(predict(fit, x = sm$time)$y)
		} else sm$mean_smooth <- sm$mean
		sm
	})
	if (!nrow(rows)) return(blank_plot(title, "No bin reached the minimum sample size"))
	labs0 <- rows |>
		group_by(feature) |>
		summarise(p_global = first(p_global), N_total = first(N_total), .groups = "drop") |>
		mutate(label = sprintf(
			"%s (N=%s, Pglobal=%s)", feature, format(N_total, big.mark = ","),
			ifelse(is.finite(p_global), format.pval(p_global, digits = 2, eps = 1e-99), "NA")
		))
	rows <- rows |> left_join(labs0 |> select(feature, label), by = "feature")
	yr <- range(rows$mean, na.rm = TRUE) ; if (!all(is.finite(yr)) || diff(yr) < .2) yr <- c( - 1, 1)
	# Put the actual case-count histogram on the primary scale and biomarker
	# trajectories on the secondary scale, matching yy.R Fig1 panel a.
	pad <- .14 * diff(yr) ; yr <- yr + c( - pad, pad) ; hmax <- max(h$n, 1)
	plot_max <- hmax * 1.12
	scale_factor <- plot_max / diff(yr)
	z_to_count <- function(z) (z - yr[1]) * scale_factor
	rows <- rows |> mutate(
		line_value = if (smooth_lines) mean_smooth else mean,
		y_plot = z_to_count(line_value), y_raw = z_to_count(mean)
	)
	ggplot() +
		geom_col(data = h, aes(mid, n, fill = side), width = diff(br)[1] * .98, color = "white", alpha = .55, inherit.aes = FALSE) +
		geom_vline(xintercept = 0, linetype = 2, color = "grey38", linewidth = .75) +
		geom_hline(yintercept = z_to_count(0), linetype = 3, color = "grey62") +
		geom_line(data = rows, aes(time, y_plot, color = label, group = label), linewidth = 1.05) +
		geom_point(
			data = rows, aes(time, y_raw, color = label), size = 1.65,
			alpha = if (smooth_lines) .42 else .90
		) +
		scale_fill_manual(values = c(`Prevalent (Yang)` = "grey58", `Incident (Yin)` = "#78B7E5"), guide = "none") +
		scale_y_continuous(
			name = "Patient count", limits = c(0, plot_max),
			sec.axis = sec_axis( ~ . / scale_factor + yr[1], name = ylab)
		) +
		labs(
			title = title,
			subtitle = if (smooth_lines)
				"Weighted smoothing-spline lines; circles are observed bin means"
			else "Original display: solid circles are observed bin means joined without smoothing",
			x = "Years relative to baseline", color = NULL
		) +
		coord_cartesian(xlim = c( - lim, lim), clip = "off") +
		theme_5c(11) +
		guides(color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(linewidth = 1.1, size = 2.8))) +
		theme(
			legend.position = c(.985, .985), legend.justification = c(1, 1), legend.direction = "horizontal", legend.box = "horizontal",
			legend.background = element_rect(fill = scales::alpha("white", .82), color = "grey75", linewidth = .25),
			legend.text = element_text(size = 8.2, face = "bold"), legend.margin = margin(3, 4, 3, 4)
		)
}


# 🚩 Fine temporal gradient (no extrapolated smoothed curve)
make_gradient_panel <- function(dat, features, bvar, side = c("incident", "prevalent"), step = .25,
				min_bin_n = 20, title = NULL, covars = character(),
				max_year = YY_MAX_YEAR) {
	side <- match.arg(side) ; features <- intersect(features, names(dat))
	# Selection differs by panel (incident top ten versus prevalent top ten), but
	# both panels deliberately show the complete -16..+16 diagnosis-time domain.
	# This is a cross-sectional pseudo-trajectory from one baseline sample per
	# person, not longitudinal within-person change.
	dd <- dat |>
		filter(is.finite(.data[[bvar]]), abs(.data[[bvar]]) <= max_year) |>
		mutate(.plot_time = .data[[bvar]])
	domain_note <- if (side == "incident") "Features ranked by incident Cox P" else "Features ranked by baseline-prevalent logistic P"
	rows <- map_dfr(features, function(x) {
		x0 <- suppressWarnings(as.numeric(dd[[x]])) ; if (sum(is.finite(x0)) < 50) return(tibble())
		z <- residualize_z(dd, x, covars)
		tb <- pmax( - max_year + step / 2, pmin(
			max_year - step / 2,
			floor((dd$.plot_time + max_year) / step) * step - max_year + step / 2
		))
		tibble(time = tb, z = z) |>
			filter(is.finite(time), is.finite(z)) |>
			group_by(time) |>
			summarise(mean = mean(z), N = n(), .groups = "drop") |>
			mutate(feature = x, reliable = N >= min_bin_n)
	})
	if (!nrow(rows)) return(blank_plot(title %||% "Temporal gradient", "No sufficiently populated time bins"))
	ord <- rev(features[features %in% unique(rows$feature)]) ; rows$feature <- factor(rows$feature, levels = ord)
	ggplot(rows, aes(time, feature)) +
		geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
		geom_point(aes(size = pmin(N, 200), fill = mean, alpha = reliable), shape = 21, color = "grey55", stroke = .15) +
		scale_fill_gradient2(low = "#3B78A8", mid = "#F7F7F7", high = "#C86B4A", midpoint = 0, name = "Mean z") +
		scale_size_continuous(range = c(.6, 4), name = "Bin N") +
		scale_alpha_manual(values = c(`TRUE` = .92, `FALSE` = .35), guide = "none") +
		scale_x_continuous(limits = c( - max_year, max_year), breaks = seq( - max_year, max_year, by = 4), expand = expansion(mult = c(.01, .01))) +
		labs(
			title = title, subtitle = paste0(domain_note, "; behavioral-LE4-residualized; faint circles have N < ", min_bin_n),
			x = "Years of recorded diagnosis relative to baseline blood draw", y = NULL
		) +
		theme_5c(10) +
		theme(legend.position = "right", axis.text.y = element_text(size = 8.5, face = "bold"))
}


# 🚩 Publication-style YY circle matrix
# Used as the first row of Fig4. No line joins adjacent circles: each dot is an
# observed time-bin summary.
make_gradient_circle_panel <- function(dat, features, bvar, step = GRADIENT_STEP, min_bin_n = MIN_BIN_N,
				max_year = YY_MAX_YEAR, title = "a. Temporal biomarker gradients") {
	features <- intersect(features, names(dat)) ; d <- dat |> filter(is.finite(.data[[bvar]]), abs(.data[[bvar]]) <= max_year)
	rows <- map_dfr(features, function(x) {
		x0 <- suppressWarnings(as.numeric(d[[x]])) ; if (sum(is.finite(x0)) < 50) return(tibble())
		z <- as.numeric(scale(x0)) ; tb <- floor(d[[bvar]] / step) * step + step / 2
		tibble(time = tb, z = z) |>
			filter(is.finite(time), is.finite(z)) |>
			group_by(time) |>
			summarise(mean = mean(z), N = n(), .groups = "drop") |>
			filter(N >= min_bin_n) |>
			mutate(feature = x)
	})
	if (!nrow(rows)) return(blank_plot(title, "No sufficiently populated time bins"))
	# Keep strongest proteins at the top, matching the conventional circle-matrix style.
	ord <- rev(features[features %in% unique(rows$feature)]) ; rows$feature <- factor(rows$feature, levels = ord)
	ggplot(rows, aes(time, feature)) +
		geom_vline(xintercept = 0, linetype = 2, color = "grey38", linewidth = .7) +
		geom_point(aes(size = N, fill = mean), shape = 21, color = "grey72", stroke = .22, alpha = .95) +
		scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0, name = "Mean z") +
		scale_size_continuous(range = c(.8, 4.6), name = "Bin N") +
		scale_x_continuous(breaks = pretty(c( - max_year, max_year), n = 7)) +
		labs(title = title, x = "Years relative to baseline", y = NULL) +
		theme_5c(10) +
		theme(legend.position = "right", axis.text.y = element_text(size = 8.6, face = "bold"))
}


# 🚩 Data-driven clustering: silhouette + gap + temporal-bin bootstrap stability
adjusted_rand <- function(a, b) {
	tab <- table(a, b) ; n <- sum(tab) ; if (n < 2) return(NA_real_)
	c2 <- function(x) x * (x - 1) / 2
	nij <- sum(c2(tab)) ; ai <- sum(c2(rowSums(tab))) ; bj <- sum(c2(colSums(tab))) ; tot <- c2(n)
	expected <- ai * bj / tot ; denom <- .5 * (ai + bj) - expected
	if (!is.finite(denom) || denom == 0) return(NA_real_)
	(nij - expected) / denom
}

mean_silhouette <- function(distmat, cl) {
	n <- length(cl) ; if (n < 3 || length(unique(cl)) < 2) return(NA_real_)
	s <- vapply(seq_len(n), function(i) {
		same <- which(cl == cl[i] & seq_len(n) != i) ; a <- if (length(same)) mean(distmat[i, same]) else 0
		oth <- setdiff(unique(cl), cl[i]) ; b <- min(vapply(oth, function(g) mean(distmat[i, cl == g]), numeric(1)))
		if (max(a, b) == 0) 0 else (b - a) / max(a, b)
	}, numeric(1))
	mean(s, na.rm = TRUE)
}

within_ss <- function(x, cl) {
	sum(vapply(unique(cl), function(g) {
		z <- x[cl == g, , drop = FALSE] ; if (nrow(z) < 2) return(0) ; sum(rowSums((z - matrix(colMeans(z), nrow(z), ncol(z), byrow = TRUE)) ^ 2))
	}, numeric(1)))
}

trajectory_matrix <- function(dat, features, bvar, step = .5, min_bin_n = 15, max_year = YY_MAX_YEAR) {
	features <- intersect(features, names(dat)) ; d <- dat |> filter(is.finite(.data[[bvar]]), abs(.data[[bvar]]) <= max_year)
	grid <- seq( - max_year + step / 2, max_year - step / 2, by = step)
	raw <- map_dfr(features, function(x) {
		z <- as.numeric(scale(suppressWarnings(as.numeric(d[[x]])))) ; tb <- floor((d[[bvar]] + max_year) / step) * step - max_year + step / 2
		tibble(time = tb, z = z) |>
			filter(is.finite(z), is.finite(time)) |>
			group_by(time) |>
			summarise(mean = mean(z), N = n(), .groups = "drop") |>
			filter(N >= min_bin_n) |>
			mutate(feature = x)
	})
	if (!nrow(raw)) return(NULL)
	wide <- raw |>
		select(feature, time, mean) |>
		pivot_wider(names_from = time, values_from = mean)
	rn <- wide$feature ; X <- as.matrix(wide[, - 1, drop = FALSE]) ; rownames(X) <- rn
	# Internal interpolation only for clustering. Plotted raw trajectories are never extrapolated.
	for (i in seq_len(nrow(X))) {
		ok <- which(is.finite(X[i, ])) ; if (length(ok) >= 2) {
			miss <- which(!is.finite(X[i, ]) & seq_len(ncol(X)) >= min(ok) & seq_len(ncol(X)) <= max(ok))
			if (length(miss)) X[i, miss] <- approx(ok, X[i, ok], xout = miss, rule = 1)$y
		}
		m <- mean(X[i, ], na.rm = TRUE) ; if (!is.finite(m)) m <- 0 ; X[i, !is.finite(X[i, ])] <- m
	}
	X <- t(scale(t(X))) ; X[!is.finite(X)] <- 0
	list(X = X, raw = raw, grid = grid)
}

choose_cluster_k <- function(X, kmax = 4, gap_B = 30, stability_B = 40, seed = SEED) {
	n <- nrow(X) ; ks <- 2 : min(kmax, n - 1) ; if (!length(ks)) return(list(k = 1, metrics = tibble()))
	D <- as.matrix(dist(X)) ; hc <- hclust(as.dist(D), method = "ward.D2")
	sil <- vapply(ks, function(k) mean_silhouette(D, cutree(hc, k)), numeric(1))
	# Gap statistic on a uniform reference box with the same feature dimensions.
	set.seed(seed) ; obsW <- vapply(ks, function(k) within_ss(X, cutree(hc, k)), numeric(1)) ; obsW <- pmax(obsW, 1e-12)
	ref_logW <- matrix(NA_real_, nrow = gap_B, ncol = length(ks))
	lo <- apply(X, 2, min) ; hi <- apply(X, 2, max)
	for (b in seq_len(gap_B)) {
		XR <- sweep(matrix(runif(n * ncol(X)), nrow = n), 2, hi - lo, "*") ; XR <- sweep(XR, 2, lo, "+")
		hR <- hclust(dist(XR), method = "ward.D2")
		ref_logW[b, ] <- vapply(ks, function(k) log(pmax(within_ss(XR, cutree(hR, k)), 1e-12)), numeric(1))
	}
	gap <- colMeans(ref_logW, na.rm = TRUE) - log(obsW) ; gap_se <- apply(ref_logW, 2, sd, na.rm = TRUE) * sqrt(1 + 1 / gap_B)
	# Stability: resample time dimensions, cluster again, compare with original labels by ARI.
	stab <- vapply(ks, function(k) {
		orig <- cutree(hc, k) ; z <- rep(NA_real_, stability_B)
		for (b in seq_len(stability_B)) {
			# Sparse outcomes can leave only one or two populated temporal bins.
			# Never request more columns than exist when bootstrapping stability.
			n_take <- min(ncol(X), max(1L, ceiling(.8 * ncol(X))))
			cols <- sort(sample(seq_len(ncol(X)), n_take, replace = FALSE))
			cb <- cutree(hclust(dist(X[, cols, drop = FALSE]), method = "ward.D2"), k)
			z[b] <- adjusted_rand(orig, cb)
		}
		median(z, na.rm = TRUE)
	}, numeric(1))
	min_cluster_n <- vapply(ks, function(k) min(table(cutree(hc, k))), numeric(1))
	min_allowed <- max(3L, ceiling(.05 * n))
	met <- tibble(
		k = ks, silhouette = sil, gap = gap, gap_se = gap_se, stability = stab,
		min_cluster_n = min_cluster_n, passes_min_cluster = min_cluster_n >= min_allowed
	)
	rank_metric <- function(x) rank( - ifelse(is.finite(x), x, - Inf), ties.method = "min")
	met <- met |> mutate(rank_sum = rank_metric(silhouette) + rank_metric(gap) + rank_metric(stability))
	eligible <- met |> filter(passes_min_cluster)
	if (!nrow(eligible)) eligible <- met |> filter(k == min(k))
	best <- eligible |>
		arrange(rank_sum, desc(silhouette), desc(stability), desc(gap)) |>
		slice(1) |>
		pull(k)
	list(k = best, metrics = met, hc = hc)
}

make_cluster_figure <- function(dat, features, bvar, step = CLUSTER_STEP, max_year = YY_MAX_YEAR,
				side = c("incident", "prevalent"), panel_label = "b", fixed_k = NULL,
				saved_membership = NULL, saved_metrics = NULL) {
	side <- match.arg(side)
	train_dat <- dat |> filter(if (side == "incident") .data[[bvar]] > 0 else .data[[bvar]] < 0)
	tm <- trajectory_matrix(train_dat, features, bvar, step = step, min_bin_n = max(10L, floor(MIN_BIN_N * .75)), max_year = max_year)
	if (is.null(tm) || nrow(tm$X) < 3 || (!is.null(saved_membership) && nrow(saved_membership) == 0L)) {
		p0 <- blank_plot(
			paste0(panel_label, ". ", ifelse(side == "incident", "Incident", "Baseline-prevalent"), " trajectory clusters"),
			"Insufficient trajectories"
		)
		return(list(plot = p0, cluster_plot = p0, diagnostic_plot = p0, cluster = tibble(), metrics = tibble(), raw = tibble(), k = 1L))
	}
	if (!is.null(saved_membership)) {
		cl <- saved_membership |> select(feature, cluster)
		ck <- list(k = length(unique(cl$cluster)), metrics = saved_metrics %||% tibble())
	} else {
		ck <- choose_cluster_k(tm$X, CLUSTER_K_MAX, CLUSTER_GAP_B, CLUSTER_STABILITY_B)
		if (!is.null(fixed_k)) ck$k <- min(max(1L, as.integer(fixed_k)), nrow(tm$X))
		cl <- tibble(feature = rownames(tm$X), cluster = cutree(ck$hc, ck$k))
	}
	full_tm <- trajectory_matrix(dat, cl$feature, bvar, step = step, min_bin_n = max(10L, floor(MIN_BIN_N * .75)), max_year = max_year)
	lines <- (full_tm$raw %||% tm$raw) |> inner_join(cl, by = "feature")
	cnt <- cl |> count(cluster, name = "proteins")
	lines <- lines |>
		left_join(cnt, by = "cluster") |>
		mutate(panel = factor(cluster, levels = sort(unique(cluster)), labels = paste0("Cluster ", sort(unique(cluster)), " (N=", cnt$proteins[match(sort(unique(cluster)), cnt$cluster)], ")")))
	centroid <- lines |>
		group_by(cluster, panel, time) |>
		summarise(estimate = weighted.mean(mean, pmax(N, 1), na.rm = TRUE), lo = quantile(mean, .10, na.rm = TRUE), hi = quantile(mean, .90, na.rm = TRUE), N = sum(N), .groups = "drop") |>
		group_by(cluster, panel) |>
		group_modify( ~ {
			q <- .x |>
				filter(is.finite(time), is.finite(estimate)) |>
				arrange(time) ; if (nrow(q) < 4) return(q)
			fit <- tryCatch(smooth.spline(q$time, q$estimate, w = sqrt(pmax(q$N, 1)), spar = C1_YY_SPAR), error = function(e) NULL)
			if (!is.null(fit)) q$estimate <- predict(fit, q$time)$y ; q
		}) |>
		ungroup()
	pA <- ggplot(lines, aes(time, mean, group = feature)) +
		geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
		geom_line(color = "grey55", alpha = .16, linewidth = .36) +
		geom_ribbon(data = centroid, aes(x = time, ymin = lo, ymax = hi, group = cluster, fill = factor(cluster)), inherit.aes = FALSE, alpha = .16, color = NA) +
		geom_line(data = centroid, aes(time, estimate, group = cluster, color = factor(cluster)), inherit.aes = FALSE, linewidth = 1.45) +
		facet_wrap( ~ panel, nrow = 1, scales = "free_y") +
		scale_color_brewer(palette = "Dark2", guide = "none") +
		scale_fill_brewer(palette = "Dark2", guide = "none") +
		labs(
			title = paste0(
				panel_label, ". ", ifelse(side == "incident", "Incident", "Baseline-prevalent"),
				" trajectory clusters (top ", nrow(tm$X), " biomarkers; k = ", ck$k, ")"
			),
			subtitle = "Clusters are learned on the named side; both Yin and Yang data are displayed. Thick lines are weighted smoothing splines; ribbons are 10th-90th percentiles.",
			x = "Years relative to baseline", y = "Mean biomarker z-score"
		) +
		theme_5c(10)
	# Do NOT min-max rescale the three criteria together. A constant stability
	# series can legitimately equal 1, but rescaling formerly collapsed the
	# display into horizontal 0/1 lines and hid the silhouette/gap structure.
	md <- ck$metrics |>
		select(k, silhouette, gap, stability) |>
		pivot_longer( - k, names_to = "criterion", values_to = "value") |>
		mutate(criterion = factor(criterion,
			levels = c("silhouette", "gap", "stability"),
			labels = c("Silhouette", "Gap statistic", "Bootstrap stability (ARI)")
		))
	bestpts <- md |> filter(k == ck$k)
	pB <- ggplot(md, aes(k, value, group = criterion)) +
		geom_vline(xintercept = ck$k, linetype = 2, color = "grey40") +
		geom_line(linewidth = .85, color = "#3F6F9F") +
		geom_point(size = 2, color = "#3F6F9F") +
		geom_point(data = bestpts, size = 3.1, shape = 21, fill = "white", stroke = .9, color = "#B2182B") +
		facet_wrap( ~ criterion, nrow = 1, scales = "free_y") +
		scale_x_continuous(breaks = unique(md$k)) +
		labs(title = "c. Cluster-number diagnostics", x = "Number of clusters (k)", y = "Native criterion value") +
		theme_5c(9) +
		theme(legend.position = "none")
	list(
		plot = pA / pB + plot_layout(heights = c(3.2, 1)), cluster_plot = pA, diagnostic_plot = pB,
		cluster = cl, metrics = ck$metrics, raw = tm$raw, k = ck$k
	)
}


# 🚩 Adjusted quintile analysis
plot_quantile_top <- function(dat, features, tvar, evar, covars, title_prefix) {
	features <- intersect(features, names(dat)) ; covars <- intersect(covars, names(dat))
	map(features, function(x) {
		need <- unique(c(tvar, evar, x, covars)) ; d <- dat[, need, drop = FALSE] ; d <- d[complete.cases(d), , drop = FALSE]
		if (nrow(d) < 500 || sum(d[[evar]] == 1) < 20) return(blank_plot(x, "Insufficient complete-case follow-up"))
		d$value <- suppressWarnings(as.numeric(d[[x]])) ; d$quintile <- factor(ntile(d$value, 5), levels = 1 : 5, labels = paste0("Q", 1 : 5))
		sf <- tryCatch(survfit(Surv(d[[tvar]], d[[evar]]) ~ d$quintile), error = function(e) NULL)
		if (is.null(sf)) return(blank_plot(x, "Survival model failed"))
		ss <- summary(sf) ; sdf <- tibble(time = ss$time, surv = ss$surv, lo = ss$lower, hi = ss$upper, strata = as.character(ss$strata)) |>
			mutate(quintile = str_extract(strata, "Q[1-5]"))
		f <- as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ relevel(quintile,'Q1') + ", paste(bt(covars), collapse = " + ")))
		fit <- tryCatch(coxph(f, d, ties = "efron"), error = function(e) NULL) ; sm <- if (is.null(fit)) NULL else coef(summary(fit))
		rn <- if (is.null(sm)) character() else rownames(sm) ; i5 <- grep("Q5", rn, fixed = TRUE)[1]
		ann <- if (is.na(i5) || !length(i5)) "" else sprintf(
			"Adjusted Q5 vs Q1: HR %.2f (%.2f–%.2f), P=%s",
			exp(sm[i5, "coef"]), exp(sm[i5, "coef"] - 1.96 * sm[i5, "se(coef)"]), exp(sm[i5, "coef"] + 1.96 * sm[i5, "se(coef)"]),
			format.pval(sm[i5, "Pr(>|z|)"], digits = 2, eps = 1e-99)
		)
		p <- ggplot(sdf, aes(time, 1 - surv, color = quintile, fill = quintile)) +
			geom_step(linewidth = .7) +
			geom_ribbon(aes(ymin = 1 - hi, ymax = 1 - lo), alpha = .07, color = NA) +
			annotate("text", x = Inf, y = Inf, label = ann, hjust = 1.03, vjust = 1.5, size = 3, fontface = "bold") +
			scale_color_brewer(palette = "Spectral", direction =  - 1) +
			scale_fill_brewer(palette = "Spectral", direction =  - 1) +
			labs(title = x, x = "Years after baseline", y = "Cumulative incidence", color = NULL, fill = NULL) +
			theme_5c(10) +
			theme(legend.position = "top")
		curves <- as.data.frame(sdf)
		curves$feature <- x
		results <- list(curves = curves)
		if (length(i5) && !is.na(i5)) results$contrast <- data.frame(
			feature = x, contrast = 'Q5 versus Q1', beta = sm[i5, 'coef'], SE = sm[i5, 'se(coef)'],
			HR = exp(sm[i5, 'coef']), lower95 = exp(sm[i5, 'coef'] - 1.96 * sm[i5, 'se(coef)']),
			upper95 = exp(sm[i5, 'coef'] + 1.96 * sm[i5, 'se(coef)']), p = sm[i5, 'Pr(>|z|)'],
			N = nrow(d), events = sum(d[[evar]] == 1)
		)
		attr(p, 'le8_result_tables') <- results
		p
	})
}


# 🚩 Functional enrichment
# Protein enrichment prefers reproducible offline GO analysis with the assayed
# proteome as universe, then falls back to g:Profiler.
run_offline_go_enrichment <- function(sig, background) {
	# A direct hypergeometric ORA avoids making clusterProfiler a hard dependency.
	# AnnotationDbi + org.Hs.eg.db are sufficient; GO.db is used only for labels.
	need <- c("AnnotationDbi", "org.Hs.eg.db")
	if (!all(vapply(need, requireNamespace, logical(1), quietly = TRUE))) return(tibble())
	mp <- tryCatch(
		suppressMessages(AnnotationDbi::select(
			org.Hs.eg.db::org.Hs.eg.db,
			keys = background, keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")
		)),
		error = function(e) tibble()
	) |>
		as_tibble() |>
		filter(!is.na(SYMBOL), !is.na(GO), ONTOLOGY %in% c("BP", "MF", "CC")) |>
		distinct(SYMBOL, GO, ONTOLOGY)
	if (!nrow(mp)) return(tibble())
	nm <- if (requireNamespace("GO.db", quietly = TRUE)) tryCatch(suppressMessages(AnnotationDbi::select(
		GO.db::GO.db,
		keys = unique(mp$GO), keytype = "GOID", columns = "TERM"
	)), error = function(e) tibble()) |>
		as_tibble() |>
		transmute(GO = GOID, term_name = TERM) |>
		distinct(GO, .keep_all = TRUE) else
		tibble(GO = unique(mp$GO), term_name = unique(mp$GO))
	map_dfr(c("BP", "MF", "CC"), function(ont) {
		mo <- mp |> filter(ONTOLOGY == ont) ; if (!nrow(mo)) return(tibble())
		universe <- intersect(background, unique(mo$SYMBOL)) ; query <- intersect(sig, universe)
		if (length(universe) < 20 || length(query) < 2) return(tibble())
		mo |>
			filter(SYMBOL %in% universe) |>
			group_by(GO) |>
			summarise(term_size = n_distinct(SYMBOL), intersection_size = n_distinct(SYMBOL[SYMBOL %in% query]), .groups = "drop") |>
			filter(term_size >= 5, intersection_size >= 1) |>
			mutate(
				p_raw = phyper(intersection_size - 1, term_size, length(universe) - term_size, length(query), lower.tail = FALSE),
				adjusted_p = p.adjust(p_raw, "BH"), source = paste0("GO:", ont)
			) |>
			left_join(nm, by = "GO") |>
			transmute(source,
				term_id = GO, term_name = coalesce(term_name, GO), adjusted_p,
				intersection_size = as.integer(intersection_size), term_size = as.integer(term_size)
			)
	})
}

run_functional_enrichment <- function(assoc, layer, background = assoc$term) {
	sig <- assoc |>
		filter(is.finite(p.value), p.value * nrow(assoc) < .05) |>
		pull(term) |>
		unique()
	background <- intersect(unique(as.character(background)), unique(as.character(assoc$term)))
	empty <- tibble(
		source = character(), term_id = character(), term_name = character(), adjusted_p = numeric(),
		intersection_size = integer(), term_size = integer(), query_size = integer(), background_n = integer(), enrichment_ratio = numeric()
	)
	if (!length(sig)) return(empty)
	# Metabolites do not have a GO-style gene ontology.  Use the curated UKB
	# met.lst hierarchy as the assay-specific universe and perform the same
	# over-representation test at super-group, group and subgroup levels.
	if (layer == "metabolite") {
		ann <- read_met_annotation(background) |> filter(trait %in% background)
		if (!nrow(ann)) return(empty)
		one_level <- function(col, source_name) {
			z <- ann |>
				transmute(trait, term = coalesce(na_if(as.character(.data[[col]]), ""), "Other")) |>
				distinct()
			map_dfr(sort(unique(z$term)), function(tt) {
				members <- z |>
					filter(term == tt) |>
					pull(trait) |>
					unique() ; m <- length(intersect(members, background)) ; k <- length(intersect(members, sig))
				if (m < 2 || k < 1) return(tibble())
				tibble(
					source = source_name, term_id = paste(source_name, tt, sep = ":"), term_name = tt,
					p_raw = phyper(k - 1, m, length(background) - m, length(sig), lower.tail = FALSE), intersection_size = k, term_size = m
				)
			})
		}
		ans <- bind_rows(
			one_level("super_group", "Metabolite super-group"), one_level("group", "Metabolite group"),
			one_level("subgroup", "Metabolite subgroup")
		)
		if (!nrow(ans)) return(empty)
		return(ans |> group_by(source) |> mutate(adjusted_p = p.adjust(p_raw, "BH")) |> ungroup() |>
			mutate(
				query_size = length(sig), background_n = length(background),
				enrichment_ratio = (intersection_size / query_size) / (term_size / background_n)
			) |>
			select(source, term_id, term_name, adjusted_p, intersection_size, term_size, query_size, background_n, enrichment_ratio) |> arrange(adjusted_p))
	}
	if (layer != "protein") return(empty)
	ans <- run_offline_go_enrichment(sig, background) ; used_custom_bg <- nrow(ans) > 0
	if (!nrow(ans) && requireNamespace("gprofiler2", quietly = TRUE)) {
		gp <- tryCatch(gprofiler2::gost(sig, organism = "hsapiens", correction_method = "fdr", sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), custom_bg = background, domain_scope = "custom", user_threshold = 1), error = function(e) NULL)
		if (is.null(gp) || is.null(gp$result) || !nrow(gp$result)) {
			used_custom_bg <- FALSE
			gp <- tryCatch(gprofiler2::gost(sig, organism = "hsapiens", correction_method = "fdr", sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"), domain_scope = "annotated", user_threshold = 1), error = function(e) NULL)
		}
		if (!is.null(gp) && !is.null(gp$result) && nrow(gp$result)) ans <- as_tibble(gp$result) |>
			transmute(source,
				term_id = native, term_name = name, adjusted_p = p_value,
				intersection_size = as.integer(intersection_size), term_size = as.integer(term_size), query_size = as.integer(query_size)
			)
	}
	if (!nrow(ans)) {
		py <- Sys.which("python3") ; script <- file.path(Sys.getenv("LE8_FDIR"), "c1.abm.py")
		if (nzchar(py) && file.exists(script)) {
			sf <- tempfile(fileext = ".txt") ; bf <- tempfile(fileext = ".txt") ; of <- tempfile(fileext = ".csv") ; on.exit(unlink(c(sf, bf, of)), add = TRUE)
			writeLines(sig, sf) ; writeLines(background, bf) ; status <- suppressWarnings(system2(py, c(script, "annotations", "--gprofiler", sf, bf, of), stdout = FALSE, stderr = FALSE))
			if (status == 0 && file.exists(of)) {
				ans <- as_tibble(data.table::fread(of, showProgress = FALSE)) ; used_custom_bg <- TRUE
			}
		}
	}
	# A network-free, non-GO fallback prevents a blank figure when Bioconductor
	# annotation is unavailable.  It is explicitly labelled as protein assay-group
	# enrichment and therefore is not confused with pathway enrichment.
	if (!nrow(ans)) {
		ann <- tryCatch(read_prot_bed(background) |> filter(protein %in% background), error = function(e) tibble())
		if (nrow(ann)) {
			ans <- map_dfr(sort(unique(ann$group)), function(g) {
				members <- ann |>
					filter(group == g) |>
					pull(protein) |>
					unique() ; m <- length(members) ; k <- length(intersect(members, sig))
				if (m < 2 || k < 1) return(tibble())
				tibble(
					source = "protein assay group", term_id = paste0("PROT:", g), term_name = g,
					adjusted_p = phyper(k - 1, m, length(background) - m, length(sig), lower.tail = FALSE),
					intersection_size = k, term_size = m, query_size = length(sig), background_n = length(background)
				)
			})
			if (nrow(ans)) ans <- ans |> mutate(adjusted_p = p.adjust(adjusted_p, "BH")) ; used_custom_bg <- TRUE
		}
	}
	if (!nrow(ans)) return(empty)
	for (nm in c("term_size", "query_size")) if (!nm %in% names(ans)) ans[[nm]] <- NA_integer_
	ans |>
		mutate(
			query_size = coalesce(as.integer(query_size), length(sig)),
			background_n = if (used_custom_bg) length(background) else NA_integer_,
			enrichment_ratio = ifelse(is.finite(term_size) & term_size > 0 & is.finite(background_n) & background_n > 0,
				(intersection_size / query_size) / (term_size / background_n), NA_real_
			)
		) |>
		select(source, term_id, term_name, adjusted_p, intersection_size, term_size, query_size, background_n, enrichment_ratio) |>
		arrange(adjusted_p)
}

plot_functional_enrichment <- function(enrich, n_sig, assoc = NULL, layer = "protein", panel_prefix = NULL) {
	add_prefix <- function(x) if (is.null(panel_prefix) || !nzchar(panel_prefix)) x else paste0(panel_prefix, "\n", x)
	if (!nrow(enrich)) {
		return(blank_plot(
			add_prefix(paste0("Functional enrichment (", n_sig, " significant biomarkers)")),
			if (layer == "protein") "No GO/pathway result was available; biomarker P values are not shown as pathway enrichment"
			else "No metabolite hierarchy contained an eligible enriched class"
		))
	}
	levels0 <- c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "protein assay group", "Metabolite super-group", "Metabolite group", "Metabolite subgroup")
	labs0 <- c(
		`GO:BP` = "BP", `GO:MF` = "MF", `GO:CC` = "CC", KEGG = "KEGG", REAC = "Reactome", `protein assay group` = "protein group",
		`Metabolite super-group` = "Super-group", `Metabolite group` = "Group", `Metabolite subgroup` = "Subgroup"
	)
	is_met_enrichment <- any(str_detect(enrich$source, "^Metabolite"))
	levels0 <- levels0[levels0 %in% unique(enrich$source)]
	d <- enrich |>
		filter(is.finite(adjusted_p)) |>
		group_by(source) |>
		slice_min(adjusted_p, n = 6, with_ties = FALSE) |>
		ungroup() |>
		mutate(
			source = factor(source, levels = levels0, labels = unname(labs0[levels0])), score =  - log10(pmax(adjusted_p, .Machine$double.xmin)),
			term_label = str_wrap(term_name, 34), term_key = paste(source, term_id, sep = "___")
		) |>
		arrange(source, adjusted_p) |>
		mutate(term_key = factor(term_key, levels = rev(unique(term_key))))
	if (!nrow(d)) return(blank_plot(add_prefix("Functional enrichment"), "No finite enrichment result was available"))
	use_ratio <- all(is.finite(d$enrichment_ratio))
	d$x_value <- if (use_ratio) d$enrichment_ratio else d$score
	ggplot(d, aes(x_value, term_key, color = score, size = intersection_size)) +
		{
			if (use_ratio) geom_vline(xintercept = 1, linetype = 3, color = "grey70") else geom_blank()
		} +
		geom_point(alpha = .9) +
		facet_grid(source ~ ., scales = "free_y", space = "free_y") +
		scale_y_discrete(labels = setNames(d$term_label, as.character(d$term_key))) +
		scale_color_viridis_c(option = "C", direction =  - 1, name = expression( - log[10](FDR))) +
		scale_size_continuous(range = c(2.2, 7), name = "Overlapping biomarkers") +
		labs(
			title = add_prefix(paste0("Functional enrichment of ", n_sig, " Bonferroni-significant biomarkers")),
			subtitle = paste0(
				if (is_met_enrichment) "Assay-specific metabolite hierarchy" else if (any(d$source == "protein assay group")) "protein assay-group fallback (not GO/pathway evidence)" else if (any(is.na(d$background_n))) "Annotated-genome background" else "All assayed biomarkers as background",
				"; ", sum(d$adjusted_p < .05), " displayed terms have FDR < 0.05"
			),
			x = if (use_ratio) "Fold enrichment" else expression( - log[10](FDR)), y = NULL
		) +
		theme_5c(9) +
		theme(
			legend.position = "right", panel.grid.major.y = element_line(color = "grey93", linewidth = .25),
			strip.text.y = element_text(angle = 0), plot.title.position = "plot", plot.margin = margin(7, 13, 7, 13)
		)
}

plot_sameN_attenuation <- function(att, assoc_adj2, top_n = 30L) {
	all_d <- att |>
		filter(is.finite(beta_basic_sameN), is.finite(beta_adj2_sameN)) |>
		left_join(assoc_adj2 |> select(term, p_adj2 = p.value), by = "term") |>
		mutate(
			direction = ifelse(beta_adj2_sameN >= 0, "Positive", "Inverse"),
			log2_ratio_plot = pmax( - 2, pmin(2, log2_effect_ratio))
		)
	if (!nrow(all_d)) return(blank_plot("Same-sample covariate attenuation", "No same-sample comparison was estimable"))
	pA <- ggplot(all_d, aes(beta_basic_sameN, beta_adj2_sameN, color = log2_ratio_plot, shape = sign_flip)) +
		geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
		geom_hline(yintercept = 0, color = "grey88") +
		geom_vline(xintercept = 0, color = "grey88") +
		geom_point(alpha = .68, size = 1.8) +
		scale_color_gradient2(
			low = "#2C7FB8", mid = "grey78", high = "#D7301F", midpoint = 0,
			breaks = c( - 2, - 1, 0, 1, 2), labels = c("0.25x", "0.5x", "1x", "2x", "4x"), name = "Adj2/basic\n|beta| ratio"
		) +
		scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), na.value = 1, name = "Direction flip") +
		labs(
			title = "a. Same-participant basic versus behavioral-LE4 effects",
			subtitle = "Colour is omitted when |basic beta| < 0.02, where percentage attenuation is unstable",
			x = "Basic Cox beta", y = "Adj2 Cox beta"
		) +
		theme_5c(9)
	d <- all_d |>
		arrange(p_adj2) |>
		slice_head(n = min(top_n, 24L)) |>
		mutate(
			term = factor(term, levels = rev(term)),
			attenuation_label = case_when(
				is.na(effect_ratio) ~ "basic beta near 0",
				sign_flip ~ sprintf("%.2fx; flip", effect_ratio), TRUE ~ sprintf("%.2fx", effect_ratio)
			)
		)
	xr <- range(c(d$beta_basic_sameN, d$beta_adj2_sameN), na.rm = TRUE) ; span <- diff(xr) ; if (!is.finite(span) || span == 0) span <- 1
	ann_x <- xr[2] + .10 * span
	pB <- ggplot(d, aes(y = term)) +
		geom_vline(xintercept = 0, color = "grey82") +
		geom_segment(aes(x = beta_basic_sameN, xend = beta_adj2_sameN, yend = term), color = "grey72", linewidth = .8) +
		geom_point(aes(x = beta_basic_sameN, color = "Basic"), size = 2.5) +
		geom_point(aes(x = beta_adj2_sameN, color = "Behavioral LE4"), size = 2.8) +
		geom_text(aes(x = ann_x, label = attenuation_label), hjust = 0, size = 3, fontface = "bold", color = "grey35") +
		scale_color_manual(values = c(Basic = "grey55", `Behavioral LE4` = "#3F78A8"), name = NULL) +
		scale_x_continuous(limits = c(xr[1] - .04 * span, xr[2] + .25 * span)) +
		labs(
			title = "b. Strongest behavioral-LE4 associations",
			subtitle = "Basic and LE4 Cox models use the same complete-case sample; labels show |LE4 beta| / |basic beta|",
			x = "Log hazard ratio per 1-SD biomarker", y = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "top")
	pA | pB
}

c1_panel_n <- function(x) {
	n <- suppressWarnings(as.integer(x$N_total))
	n <- sort(unique(n[is.finite(n)]))
	if (!length(n)) return("NA")
	if (length(n) == 1L) return(format(n, big.mark = ",", scientific = FALSE, trim = TRUE))
	paste(format(range(n), big.mark = ",", scientific = FALSE, trim = TRUE), collapse = "–")
}

plot_c1_fig12 <- function(layer, features_all, pgs_incident, pgs_prevalent,
				pgs_attained_age, assoc_adj2, prevalent_adj2,
				birthline_adj2, outdir) {
	omic_name <- if (layer == "protein") "protein" else "metabolite"
	measured_model <- if (le8_custom_adjustment()) paste0(length(le8_custom_covars), " selected covariates") else "behavioral-LE4 adjusted"
	pgs_title <- function(panel, analysis, x, include_outcome = FALSE) {
		paste0(
			panel, ". ", analysis, if (include_outcome) paste0(" ", Y) else "",
			" — inherited ", omic_name, " PGS, N = ", c1_panel_n(x)
		)
	}
	circle_pgs_title <- function(panel, analysis, x) {
		str_replace(pgs_title(panel, analysis, x, TRUE), ", N = ", ",\nN = ")
	}
	observed_title <- function(panel, analysis, x, method = "") {
		paste0(
			panel, ". ", analysis, if (nzchar(method)) paste0(" ", method) else "",
			", N = ", c1_panel_n(x)
		)
	}

	# The inherited-score panels use everyone with a matched PGS and phenotype
	# record. They must not be restricted to participants with the corresponding
	# adult omics assay; that restriction discards most of the genetic cohort.
	if (layer == "protein") {
		p1a <- plot_pwas_manhattan(pgs_incident, features_all, pgs_title("a", "Incident", pgs_incident, TRUE))
		p1b <- plot_pwas_manhattan(pgs_prevalent, features_all, pgs_title("b", "Baseline-prevalent", pgs_prevalent, TRUE))
		p1c <- plot_pwas_manhattan(assoc_adj2, features_all, observed_title("c", paste0("Incident ", Y), assoc_adj2, paste0("— ", measured_model, " Cox")))
		p1d <- plot_pwas_manhattan(prevalent_adj2, features_all, observed_title("d", paste0("Baseline-prevalent ", Y), prevalent_adj2, paste0("— ", measured_model, " logistic")))
		p1e <- plot_pwas_manhattan(pgs_attained_age, features_all, pgs_title("e", "Attained-age", pgs_attained_age, TRUE))
		p1f <- plot_pwas_manhattan(birthline_adj2, features_all, observed_title("f", paste0("Attained-age ", Y), birthline_adj2, paste0("— measured adult level, ", measured_model)))
		save_plot((p1a | p1c) / plot_spacer() / (p1b | p1d) / plot_spacer() / (p1e | p1f) +
			plot_layout(heights = c(1, .06, 1, .06, 1)), "c1.Fig1.mh.png", 24, 19, outdir = outdir)
		p2a <- plot_volcano(pgs_incident, "beta", "term", pgs_title("a", "Incident", pgs_incident), .05 / max(1, nrow(pgs_incident)), TOP_N)
		p2b <- plot_volcano(pgs_prevalent, "beta", "term", pgs_title("b", "Prevalent", pgs_prevalent), .05 / max(1, nrow(pgs_prevalent)), TOP_N)
		p2c <- plot_volcano(assoc_adj2, "beta", "term", observed_title("c", "Incident ProtWAS", assoc_adj2, paste0("— ", measured_model)), .05 / max(1, nrow(assoc_adj2)), TOP_N)
		p2d <- plot_volcano(prevalent_adj2, "beta", "term", observed_title("d", "Prevalent ProtWAS", prevalent_adj2, paste0("— ", measured_model, " logistic")), .05 / max(1, nrow(prevalent_adj2)), TOP_N)
		p2e <- plot_volcano(pgs_attained_age, "beta", "term", pgs_title("e", "Attained-age", pgs_attained_age), .05 / max(1, nrow(pgs_attained_age)), TOP_N)
		p2f <- plot_volcano(birthline_adj2, "beta", "term", observed_title("f", "Attained-age ProtWAS", birthline_adj2, paste0("— measured adult level, ", measured_model)), .05 / max(1, nrow(birthline_adj2)), TOP_N)
		save_plot((p2a | p2c) / plot_spacer() / (p2b | p2d) / plot_spacer() / (p2e | p2f) +
			plot_layout(heights = c(1, .06, 1, .06, 1)), "c1.Fig2.vc.png", 19, 19, outdir = outdir)
	} else {
		ann <- read_met_annotation(features_all) ; annot <- function(z) z |>
			mutate(trait = term) |>
			left_join(ann, by = "trait") |>
			mutate(label = coalesce(label, term), group = coalesce(group, "Other"))
		circ_data <- list(annot(pgs_incident), annot(pgs_prevalent), annot(assoc_adj2), annot(prevalent_adj2), annot(pgs_attained_age), annot(birthline_adj2))
		shared_beta_lim <- as.numeric(Sys.getenv('C1_CIRCLE_BETA_CAP', unset = '.3'))
		shared_fdr_cap <- as.numeric(Sys.getenv('C1_CIRCLE_FDR_CAP', unset = '50'))
		panel_caps <- vapply(circ_data, function(d) {
			scores <- le8_met_circle_data(d)$neglog10_FDR
			min(shared_fdr_cap, max(5, ceiling(max(0, scores, na.rm = TRUE) / 5) * 5))
		}, numeric(1))
		circ <- function(d, title) {
			scores <- le8_met_circle_data(d)$neglog10_FDR
			cap <- min(shared_fdr_cap, max(5, ceiling(max(0, scores, na.rm = TRUE) / 5) * 5))
			plot_met_circle(d, title, beta_limit = shared_beta_lim, fdr_cap = cap)
		}
		p1a <- circ(circ_data[[1]], circle_pgs_title("a", "Incident", pgs_incident))
		p1b <- circ(circ_data[[2]], circle_pgs_title("b", "Prevalent", pgs_prevalent))
		p1c <- circ(circ_data[[3]], observed_title("c", paste0("Incident ", Y, " MWAS"), assoc_adj2, paste0("— ", measured_model)))
		p1d <- circ(circ_data[[4]], observed_title("d", paste0("Prevalent ", Y, " MWAS"), prevalent_adj2, paste0("— ", measured_model, " logistic")))
		p1e <- circ(circ_data[[5]], circle_pgs_title("e", "Attained-age", pgs_attained_age))
		p1f <- circ(circ_data[[6]], observed_title("f", paste0("Attained-age ", Y), birthline_adj2, paste0("— measured adult level, ", measured_model)))
		# Paper_t2dm Fig2 style: large radial panels and a dedicated shared-key row.
		# A guide_area prevents four independent legends from shrinking the circles.
		p_circle <- wrap_plots(list(p1a, p1c, p1b, p1d, p1e, p1f), ncol = 2, guides = 'collect') +
			plot_annotation(caption = paste0(
				'Category-wise radial associations. Each panel selects BH FDR < 0.05 from all metabolites tested in that analysis. ',
				'Radius: -log10(FDR), with each panel cap stated above its circle (maximum ', shared_fdr_cap, '); colour: log effect per SD, capped at +/- ', shared_beta_lim,
				'. Incident and attained-age: log(HR); prevalent: log(OR). Left: inherited PGS; right: baseline measured metabolites. ',
				'Labels use assay identifiers; full names and uncapped results are in c1.circular_associations.csv.'
			))
		circle_analyses <- c('PGS_incident', 'PGS_prevalent', 'measured_incident', 'measured_prevalent', 'PGS_attained_age', 'measured_attained_age')
		circle_table <- map2_dfr(
			circ_data, circle_analyses,
 ~ le8_met_circle_data(.x) |> mutate(analysis = .y, colour_cap = shared_beta_lim)
		)
		circle_table$bar_cap <- panel_caps[match(circle_table$analysis, circle_analyses)]
		write_raw_csv(circle_table, 'c1.circular_associations.csv', le8_job_dir(outdir, 'c1_correlate'))
		save_plot(p_circle, "c1.Fig1.circular.png", 24, 36.8, outdir = outdir)
		p2a <- plot_volcano(pgs_incident, "beta", "term", pgs_title("a", "Incident", pgs_incident), .05 / max(1, nrow(pgs_incident)), 12, x_quantile = .985)
		p2b <- plot_volcano(pgs_prevalent, "beta", "term", pgs_title("b", "Prevalent", pgs_prevalent), .05 / max(1, nrow(pgs_prevalent)), 12, x_quantile = .985)
		p2c <- plot_volcano(assoc_adj2, "beta", "term", observed_title("c", "Incident MWAS", assoc_adj2, paste0("— ", measured_model)), .05 / max(1, nrow(assoc_adj2)), 12, x_quantile = .985)
		p2d <- plot_volcano(prevalent_adj2, "beta", "term", observed_title("d", "Prevalent MWAS", prevalent_adj2, paste0("— ", measured_model, " logistic")), .05 / max(1, nrow(prevalent_adj2)), 12, x_quantile = .985)
		p2e <- plot_volcano(pgs_attained_age, "beta", "term", pgs_title("e", "Attained-age", pgs_attained_age), .05 / max(1, nrow(pgs_attained_age)), 12, x_quantile = .985)
		p2f <- plot_volcano(birthline_adj2, "beta", "term", observed_title("f", "Attained-age MWAS", birthline_adj2, paste0("— measured adult level, ", measured_model)), .05 / max(1, nrow(birthline_adj2)), 12, x_quantile = .985)
		save_plot((p2a | p2c) / plot_spacer() / (p2b | p2d) / plot_spacer() / (p2e | p2f) +
			plot_layout(heights = c(1, .06, 1, .06, 1)), "c1.Fig2.vc.png", 19, 19, outdir = outdir)
	}
}

run_vldl_tg_measured_models <- function(dat, covars, tvar, evar) {
	target <- "L_VLDL_TG.pct"
	comparators <- intersect(c("L_VLDL_TG", "L_VLDL_L", "Total_TG", "ApoB", "VLDL_size"), names(dat))
	if (!target %in% names(dat)) return(tibble())
	zname <- function(x) paste0(".z_", make.names(x))
	needed <- unique(c(tvar, evar, covars, target, comparators))
	d <- dat[, needed, drop = FALSE]
	for (v in c(target, comparators)) d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
	d <- d[complete.cases(d), , drop = FALSE]
	d <- d[is.finite(d[[tvar]]) & d[[tvar]] > 0 & d[[evar]] %in% c(0, 1), , drop = FALSE]
	for (v in c(target, comparators)) d[[zname(v)]] <- as.numeric(scale(d[[v]]))
	fit_one <- function(dd, exposure, adjusters = character(), model_group = "Target conditional", model = "") {
		xs <- c(exposure, adjusters) ; empty <- tibble(model_group, model, exposure,
			adjusters =
				paste(adjusters, collapse = ";"), beta = NA_real_, std.error = NA_real_, conf.low = NA_real_,
			conf.high = NA_real_, p.value = NA_real_, N_total = nrow(dd), N_event = sum(dd[[evar]] == 1),
			max_abs_exposure_correlation = NA_real_, exposure_condition_number = NA_real_
		)
		if (nrow(dd) < 500 || sum(dd[[evar]] == 1) < 20) return(empty)
		ff <- as.formula(paste0(
			"Surv(", bt(tvar), ",", bt(evar), ") ~ ",
			paste(bt(c(xs, covars)), collapse = " + ")
		))
		fit <- tryCatch(coxph(ff, dd, ties = "efron"), error = function(e) NULL)
		if (is.null(fit)) return(empty) ; sm <- coef(summary(fit)) ; if (!exposure %in% rownames(sm)) return(empty)
		b <- sm[exposure, "coef"] ; se <- sm[exposure, "se(coef)"]
		cm <- if (length(xs) > 1L) cor(dd[, xs, drop = FALSE]) else matrix(1, 1, 1)
		max_cor <- if (length(xs) > 1L) max(abs(cm[1, - 1]), na.rm = TRUE) else 0
		cond <- if (length(xs) > 1L) tryCatch(kappa(cm), error = function(e) NA_real_) else 1
		tibble(model_group, model, exposure,
			adjusters = paste(adjusters, collapse = ";"),
			beta = b, std.error = se, conf.low = b - 1.96 * se, conf.high = b + 1.96 * se,
			p.value = sm[exposure, "Pr(>|z|)"], N_total = nrow(dd), N_event = sum(dd[[evar]] == 1),
			max_abs_exposure_correlation = max_cor, exposure_condition_number = cond
		)
	}
	target_z <- zname(target)
	conditional <- bind_rows(
		fit_one(d, target_z, model = "Target only"),
		map_dfr(comparators, function(v) fit_one(d, target_z, zname(v), model = paste0("+ ", v))),
		if (length(comparators) > 1L) fit_one(d, target_z, zname(comparators),
			model = "+ all burden/size traits"
		) else tibble()
	)

	component_vars <- c("L_VLDL_TG.pct", "L_VLDL_CE.pct", "L_VLDL_FC.pct", "L_VLDL_PL.pct")
	balances <- tibble()
	if (all(component_vars %in% names(dat))) {
		bd <- dat[, unique(c(tvar, evar, covars, component_vars)), drop = FALSE]
		for (v in component_vars) bd[[v]] <- suppressWarnings(as.numeric(bd[[v]]))
		bd <- bd[complete.cases(bd) & apply(bd[, component_vars, drop = FALSE] > 0, 1, all), , drop = FALSE]
		bd <- bd[is.finite(bd[[tvar]]) & bd[[tvar]] > 0 & bd[[evar]] %in% c(0, 1), , drop = FALSE]
		bd$.TG_vs_nonTG_balance <- log(bd$L_VLDL_TG.pct) -
			rowMeans(log(as.matrix(bd[, setdiff(component_vars, "L_VLDL_TG.pct"), drop = FALSE])))
		for (v in setdiff(component_vars, "L_VLDL_TG.pct"))
			bd[[paste0(".TG_vs_", make.names(v))]] <- log(bd$L_VLDL_TG.pct) - log(bd[[v]])
		balance_vars <- c(".TG_vs_nonTG_balance", paste0(".TG_vs_", make.names(
			setdiff(component_vars, "L_VLDL_TG.pct")
		)))
		for (v in balance_vars) bd[[v]] <- as.numeric(scale(bd[[v]]))
		balances <- map_dfr(balance_vars, function(v) fit_one(bd, v,
			model_group = "Compositional log-ratio",
			model = str_remove(v, "^\\.")
		))
	}
	bind_rows(conditional, balances) |> mutate(
		FDR = p.adjust(p.value, "BH"),
		interpretation = ifelse(model_group == "Target conditional",
			"Same-N measured target association conditional on correlated lipid burden/size traits",
			"Log-ratio balance within large VLDL; exploratory compositional analysis"
		)
	)
}

build_vldl_tg_deep_dive <- function(
	pgs_full, pgs_same, measured_le4,
	measured_basic, measured_full_le8, risk_window_le4 = tibble(),
	pgs_conditional = tibble(), measured_conditional = tibble()
) {
	size_levels <- c("XS", "S", "M", "L", "XL", "XXL")
	pct_features <- paste0(size_levels, "_VLDL_TG.pct")
	absolute_features <- paste0(size_levels, "_VLDL_TG")
	key_features <- c("L_VLDL_TG.pct", "L_VLDL_TG", "Total_TG", "ApoB", "VLDL_size")
	wanted <- unique(c(pct_features, absolute_features, key_features))
	collect_set <- function(x, predictor, model) {
		labels <- c(
			incident = "Incident", prevalent = "Baseline prevalent",
			attained_age = "Attained age"
		)
		map_dfr(names(labels), function(nm) {
			z <- as_tibble(x[[nm]] %||% tibble())
			if (!nrow(z) || !"term" %in% names(z)) return(tibble())
			z |>
				filter(term %in% wanted) |>
				transmute(
					feature = term, analysis = labels[[nm]],
					predictor, model, beta = as.numeric(beta), se = as.numeric(std.error),
					p_value = as.numeric(p.value), FDR = as.numeric(FDR),
					N = as.numeric(N_total), events = as.numeric(N_event)
				)
		})
	}
	d <- bind_rows(
		collect_set(pgs_full, "Inherited", "PGS: full genetic cohort"),
		collect_set(pgs_same, "Inherited", "PGS: same-omic subset"),
		collect_set(measured_basic, "Measured", "Measured: basic"),
		collect_set(measured_le4, "Measured", "Measured: behavioral LE4"),
		collect_set(measured_full_le8, "Measured", "Measured: full LE8 sensitivity")
	) |>
		mutate(
			conf.low = beta - 1.96 * se, conf.high = beta + 1.96 * se,
			vldl_size = str_extract(feature, "^(XS|S|M|L|XL|XXL)"),
			metric = case_when(
				feature %in% pct_features ~ "TG / total lipids (%)",
				feature %in% absolute_features ~ "Absolute TG", TRUE ~ "Comparator"
			),
			target = feature == "L_VLDL_TG.pct"
		)
	cols <- c(
		"PGS: full genetic cohort" = "#7B3294", "PGS: same-omic subset" = "#C2A5CF",
		"Measured: basic" = "#7F7F7F", "Measured: behavioral LE4" = "#008837",
		"Measured: full LE8 sensitivity" = "#D95F02"
	)
	keep_models <- c("PGS: full genetic cohort", "Measured: behavioral LE4")

	a <- d |>
		filter(analysis == "Incident", feature %in% key_features, model %in% keep_models) |>
		mutate(feature = factor(feature, levels = rev(key_features)), model = factor(model, levels = keep_models))
	pa <- if (!nrow(a)) blank_plot("a. Target and burden comparators", "No estimable effects") else
		ggplot(a, aes(beta, feature, color = model)) +
			geom_vline(xintercept = 0, color = "grey78") +
			geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
				height = .11,
				position = position_dodge(width = .52)
			) +
			geom_point(size = 2.2, position = position_dodge(width = .52)) +
			scale_color_manual(values = cols, drop = FALSE) +
			labs(
				title = "a. Incident CAD: composition versus burden",
				subtitle = "L_VLDL_TG.pct is a composition ratio; absolute TG and ApoB are separate traits",
				x = "Log HR per 1 SD", y = NULL, color = NULL
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")

	size_panel <- function(metric_name, panel, title) {
		z <- d |>
			filter(analysis == "Incident", metric == metric_name, model %in% keep_models) |>
			mutate(vldl_size = factor(vldl_size, levels = size_levels), model = factor(model, levels = keep_models))
		if (!nrow(z)) return(blank_plot(paste0(panel, ". ", title), "No estimable effects"))
		ggplot(z, aes(vldl_size, beta, color = model, group = model)) +
			geom_hline(yintercept = 0, color = "grey78") +
			geom_ribbon(aes(
				ymin = conf.low, ymax = conf.high,
				fill = model
			), alpha = .10, color = NA) +
			geom_line(linewidth = .7) +
			geom_point(aes(shape = target), size = 2.2) +
			scale_color_manual(values = cols, drop = FALSE) +
			scale_fill_manual(values = cols, drop = FALSE) +
			scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
			labs(
				title = paste0(panel, ". ", title), x = "VLDL subclass size", y = "Incident log HR per 1 SD",
				color = NULL, fill = NULL
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")
	}
	pb <- size_panel("TG / total lipids (%)", "b", "TG percentage across VLDL sizes")
	pc <- size_panel("Absolute TG", "c", "Absolute TG across VLDL sizes")

	dd <- d |>
		filter(feature == "L_VLDL_TG.pct", model %in% keep_models) |>
		mutate(
			analysis = factor(analysis, levels = c("Incident", "Attained age", "Baseline prevalent")),
			model = factor(model, levels = keep_models)
		)
	pd <- if (!nrow(dd)) blank_plot("d. Target across disease-time analyses", "No estimable effects") else
		ggplot(dd, aes(beta, analysis, color = model)) +
			geom_vline(xintercept = 0, color = "grey78") +
			geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
				height = .10,
				position = position_dodge(width = .48)
			) +
			geom_point(size = 2.3, position = position_dodge(width = .48)) +
			scale_color_manual(values = cols, drop = FALSE) +
			labs(
				title = "d. L_VLDL_TG.pct across disease-time analyses",
				subtitle = "Incident/attained-age are prospective; prevalent is treatment- and survivor-sensitive",
				x = "Log HR or log OR per 1 SD", y = NULL, color = NULL
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")

	sensitivity_models <- c(
		"Measured: basic", "Measured: behavioral LE4",
		"Measured: full LE8 sensitivity"
	)
	ee <- d |>
		filter(
			analysis == "Incident", feature %in% c("L_VLDL_TG.pct", "L_VLDL_TG", "Total_TG", "ApoB"),
			model %in% sensitivity_models
		) |>
		mutate(
			feature = factor(feature, levels = rev(c("L_VLDL_TG.pct", "L_VLDL_TG", "Total_TG", "ApoB"))),
			model = factor(model, levels = sensitivity_models)
		)
	pe <- if (!nrow(ee)) blank_plot("e. Adjustment sensitivity", "No estimable effects") else
		ggplot(ee, aes(beta, feature, color = model)) +
			geom_vline(xintercept = 0, color = "grey78") +
			geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
				height = .08,
				position = position_dodge(width = .55)
			) +
			geom_point(size = 2, position = position_dodge(width = .55)) +
			scale_color_manual(values = cols, drop = FALSE) +
			labs(
				title = "e. Measured-trait adjustment sensitivity",
				subtitle = "Behavioral LE4 is primary; full LE8 can condition on biomarker intermediates",
				x = "Incident log HR per 1 SD", y = NULL, color = NULL
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")

	ff <- as_tibble(pgs_conditional)
	if (nrow(ff) && all(c("model", "beta", "conf.low", "conf.high") %in% names(ff)))
		ff <- ff |> mutate(model = factor(model, levels = rev(model))) else ff <- tibble()
	pf <- if (!nrow(ff)) blank_plot("f. Conditional PGS models", "Conditional PGS output unavailable") else
		ggplot(ff, aes(beta, model)) +
			geom_vline(xintercept = 0, color = "grey78") +
			geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .12, color = "#7B3294") +
			geom_point(color = "#7B3294", size = 2.3) +
			labs(
				title = "f. Does the target PGS survive burden/size PGS adjustment?",
				subtitle = "Exploratory multivariable score association; this is not multivariable MR",
				x = "L_VLDL_TG.pct PGS incident log HR", y = NULL
			) +
			theme_5c(8)

	gg <- as_tibble(risk_window_le4)
	if (nrow(gg) && all(c("term", "time", "beta", "conf.low", "conf.high", "side") %in% names(gg)))
		gg <- gg |>
			filter(term %in% c("L_VLDL_TG.pct", "L_VLDL_TG", "Total_TG", "ApoB")) |>
			mutate(term = factor(term, levels = c("L_VLDL_TG.pct", "L_VLDL_TG", "Total_TG", "ApoB"))) else
		gg <- tibble()
	pg <- if (!nrow(gg)) blank_plot(
		"g. Diagnosis-anchored risk-set trajectory",
		"Trajectory is generated on the behavioral-LE4 rerun"
	) else
		ggplot(gg, aes(time, beta, color = side, fill = side, group = side)) +
			geom_hline(yintercept = 0, color = "grey78") +
			geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
			geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .12, color = NA) +
			geom_line(linewidth = .75) +
			geom_point(aes(size = N_case), alpha = .85) +
			facet_wrap( ~ term, nrow = 1, scales = "free_y") +
			scale_color_manual(values = c(
				"Pre-baseline prevalent" = "#D95F02",
				"Post-baseline incident" = "#3F78A8"
			)) +
			scale_fill_manual(values = c(
				"Pre-baseline prevalent" = "#D95F02",
				"Post-baseline incident" = "#3F78A8"
			)) +
			labs(
				title = "g. Diagnosis-anchored measured-trait risk-set trajectory",
				subtitle = "One baseline measurement per person; this is not a within-person longitudinal trajectory",
				x = "Years from baseline diagnosis boundary", y = "LE4-adjusted log effect (OR / HR)",
				color = NULL, fill = NULL, size = "Cases"
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")

	hh <- as_tibble(measured_conditional)
	if (nrow(hh) && all(c("model_group", "model", "beta", "conf.low", "conf.high") %in% names(hh)))
		hh <- hh |> mutate(model = factor(model, levels = rev(unique(model)))) else hh <- tibble()
	ph <- if (!nrow(hh)) blank_plot(
		"h. Conditional and compositional measured models",
		"Measured deep-model output unavailable"
	) else
		ggplot(hh, aes(beta, model)) +
			geom_vline(xintercept = 0, color = "grey78") +
			geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .10, color = "#008837") +
			geom_point(color = "#008837", size = 2.2) +
			facet_wrap( ~ model_group, scales = "free_y", nrow = 1) +
			labs(
				title = "h. Same-N conditional and within-particle log-ratio analyses",
				subtitle = "Condition numbers and score correlations are written to the audit tables",
				x = "Incident log HR per 1 SD", y = NULL
			) +
			theme_5c(8)

	list(
		data = d, trajectory = gg, pgs_conditional = ff, measured_conditional = hh,
		figure = (pa | pb) / (pc | pd) / (pe | pf) / pg / ph +
			plot_layout(heights = c(1, 1, 1.05, 1.08, 1.05)) +
			plot_annotation(
				title = "L_VLDL_TG.pct deep dive: inherited composition signal versus measured burden",
				caption = paste0(
					"Percentage traits are compositional. Negative PGS association is not evidence that ",
					"higher absolute VLDL triglyceride burden is protective; prioritize conditional MR/colocalization."
				)
			)
	)
}


# 🚩 Main C1
# Temporal associations and fixed-baseline risk sets

risk_window_scan <- function(dat, xs, tvar, evar, bvar, covars, cuts = C1_RISK_CUTS) {
	xs <- intersect(xs, names(dat))
	cuts <- sort(unique(cuts))
	inc <- map_dfr(seq_len(length(cuts) - 1L), function(i) map_dfr(xs, function(x) le8_interval_cox(
		dat, x, tvar,
		evar, covars, cuts[i], cuts[i + 1]
	))) |>
		mutate(side = "Post-baseline incident")
	# Prevalence remains logistic, explicitly separate from prospective HR.
	prev <- map_dfr(seq_len(length(cuts) - 1L), function(i) {
		lo <- cuts[i]
		hi <- cuts[i + 1]
		pc <- dat$prevalent_status == 1 & is.finite(dat[[bvar]]) & ( - dat[[bvar]]) >= lo & ( - dat[[bvar]]) < hi
		cc <- dat$prevalent_status == 0
		keep <- replace_na(pc, FALSE) | replace_na(cc, FALSE)
		d <- dat[keep, , drop = FALSE]
		d$.window_case <- as.integer(replace_na(pc[keep], FALSE))
		z <- logistic_scan(d, xs, covars, ".window_case")
		z |>
			transmute(term, beta, std.error,
				conf.low = beta - 1.96 * std.error, conf.high = beta + 1.96 * std.error,
				p.value, N_total, N_event, N_case = N_event, N_control = N_total - N_event, person_years = NA_real_,
				window_lo = lo, window_hi = hi, time =  - (lo + hi) / 2, effect_measure = "log OR; baseline prevalence",
				status = ifelse(is.finite(p.value), "ok", "insufficient data"), side = "Pre-baseline prevalent"
			)
	})
	bind_rows(prev, inc) |>
		group_by(side) |>
		mutate(FDR = p.adjust(p.value, "BH")) |>
		ungroup()
}

le8_time_heterogeneity <- function(dat, xs, tvar, evar, covars, cuts = c(0, 2, 5, 10, 16)) {
	if (!length(xs))
		return(tibble(
			feature = character(), N = integer(), events = integer(), p_time_interaction = numeric(),
			status = character(), FDR = numeric()
		))
	map_dfr(xs, function(x) {
		cols <- unique(c(x, tvar, evar, covars))
		d <- as.data.frame(dat[, cols, drop = FALSE])
		d <- d[complete.cases(d), , drop = FALSE]
		d <- d[d[[tvar]] > 0 & d[[evar]] %in% c(0, 1), , drop = FALSE]
		ans <- tibble(feature = x, N = nrow(d), events = sum(d[[evar]] == 1), p_time_interaction = NA_real_, status = "insufficient data")
		if (nrow(d) < 500 || sum(d[[evar]] == 1) < 40 || sd(d[[x]]) <= 0)
			return(ans)
		d$.z <- as.numeric(scale(d[[x]]))
		pieces <- lapply(seq_len(length(cuts) - 1L), function(j) {
			z <- d[d[[tvar]] > cuts[j], , drop = FALSE]
			z$.start <- cuts[j]
			z$.end <- pmin(z[[tvar]], cuts[j + 1])
			z$.ev <- as.integer(z[[evar]] == 1 & z[[tvar]] <= cuts[j + 1])
			z$.window <- j
			z
		})
		long <- bind_rows(pieces)
		long$.window <- factor(long$.window)
		cv <- paste(bt(covars), collapse = " + ")
		cv <- if (nzchar(cv))
			paste0(" + ", cv) else ""
		f0 <- as.formula(paste0("survival::Surv(.start,.end,.ev) ~ strata(.window) + .z", cv))
		f1 <- as.formula(paste0("survival::Surv(.start,.end,.ev) ~ strata(.window) + .z:factor(.window)", cv))
		a <- tryCatch(survival::coxph(f0, long), error = function(e) NULL)
		b <- tryCatch(survival::coxph(f1, long), error = function(e) NULL)
		if (is.null(a) || is.null(b))
			return(ans)
		df <- sum(is.finite(coef(b))) - sum(is.finite(coef(a)))
		if (df > 0) {
			ans$p_time_interaction <- pchisq(max(0, 2 * (b$loglik[2] - a$loglik[2])), df, lower.tail = FALSE)
			ans$status <- "interval-slope likelihood-ratio test"
		}
		ans
	}) |>
		mutate(FDR = p.adjust(p_time_interaction, "BH"))
}

le8_c1_additions <- function(dat, layer, covars, tvar, evar) {
	anchors <- intersect(
		le8_csv_env("C1_DIRECTION_ANCHORS", "PCSK9,LPA,GDF15,NTPROBNP,MMP12,L_VLDL_TG.pct,L_VLDL_TG,Total_TG,ApoB"),
		names(dat)
	)
	heterogeneity <- le8_time_heterogeneity(dat, anchors, tvar, evar, covars)
	rawdir <- le8_job_dir(if (layer == "protein")
		out.prot else out.met, "c1_correlate")
	paired <- tibble()
	status <- tibble(status = "PGS input unavailable")
	sf <- find_c1_pgs_file(layer)
	if (length(sf) == 1L && !is.na(sf) && length(anchors)) {
		scores <- read_c1_pgs(sf)
		mp <- map_c1_pgs_columns(anchors, names(scores))
		if (length(mp)) {
			scores <- scores[, unique(c("eid", unname(mp))), drop = FALSE]
			d <- dat
			d$eid <- as.character(d$eid)
			d <- inner_join(d, scores, by = "eid")
			rm(scores)
			invisible(gc())
			paired <- map_dfr(names(mp), function(x) {
				z <- d[, unique(c(x, mp[[x]], tvar, evar, covars)), drop = FALSE]
				z <- z[complete.cases(z), , drop = FALSE]
				z <- z[z[[tvar]] > 0 & z[[evar]] %in% c(0, 1), , drop = FALSE]
				z$.measured <- as.numeric(z[[x]])
				z$.pgs <- as.numeric(z[[mp[[x]]]])
				singles <- cox_scan(z, c(".measured", ".pgs"), covars, Y, time_var = tvar, event_var = evar) |>
					mutate(feature = x, model = "separate; identical people and covariates", component = ifelse(term ==
						".measured", "Measured", "PGS"), unit = "per own SD; not interchangeable concentration units")
				if (nrow(z) >= 500 && sum(z[[evar]] == 1) >= 20 && sd(z$.measured) > 0 && sd(z$.pgs) > 0) {
					z$.measured <- as.numeric(scale(z$.measured))
					z$.pgs <- as.numeric(scale(z$.pgs))
					f <- as.formula(paste0("survival::Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(
						".measured",
						".pgs", covars
					)), collapse = " + ")))
					fit <- tryCatch(survival::coxph(f, z), error = function(e) NULL)
					if (!is.null(fit)) {
						sm <- coef(summary(fit))
						j <- intersect(c(".measured", ".pgs"), rownames(sm))
						joint <- tibble(
							term = j, beta = sm[j, "coef"], std.error = sm[j, "se(coef)"], p.value = sm[
								j,
								"Pr(>|z|)"
							], N_total = nrow(z), N_event = sum(z[[evar]] == 1), feature = x, model = "joint; identical people and covariates",
							component = ifelse(j == ".measured", "Measured", "PGS"), unit = "per own SD; not interchangeable concentration units"
						) |>
							mutate(
								conf.low = exp(beta - 1.96 * std.error), conf.high = exp(beta + 1.96 * std.error),
								estimate = exp(beta)
							)
						singles <- bind_rows(singles, joint)
					}
				}
				singles
			}) |>
				mutate(FDR = p.adjust(p.value, "BH"))
			status <- tibble(status = "ok", N_overlap = nrow(d), matched_anchors = length(mp), note = "This paired comparison is not MR; significance is not compared across unequal sample sizes")
		} else status <- tibble(status = "No configured anchor matched a PGS column")
	}
	write_raw_csv(paired, "c1.paired_pgs_measured.csv", rawdir)
	write_raw_csv(heterogeneity, "c1.temporal_heterogeneity.csv", rawdir)
	le8_plot_c1_temporal(paired, heterogeneity, status, anchors, layer)
	list(paired_associations = paired, time_heterogeneity = heterogeneity, paired_status = status)
}

le8_plot_c1_temporal <- function(paired, heterogeneity, status, anchors, layer) {
	if (nrow(paired)) {
		p1 <- paired |>
			mutate(lo = beta - 1.96 * std.error, hi = beta + 1.96 * std.error) |>
			ggplot(aes(beta, feature, color = component)) +
			geom_vline(xintercept = 0) +
			geom_errorbarh(aes(
				xmin = lo,
				xmax = hi
			), height = 0.1, position = position_dodge(width = 0.4)) +
			geom_point(position = position_dodge(width = 0.4)) +
			facet_wrap( ~ model) +
			labs(
				title = "a. Paired measured and PGS associations", x = "Log HR per own SD",
				y = NULL, color = NULL
			) +
			theme_5c(9)
	} else p1 <- blank_plot("a. Paired measured and PGS associations", if (!length(anchors))
		"No configured anchor is present in this omics layer" else status$status[1])
	p2 <- if (any(is.finite(heterogeneity$p_time_interaction)))
		heterogeneity |>
			ggplot(aes( - log10(pmax(p_time_interaction, 1e-300)), feature)) +
			geom_point() +
			labs(
				title = "b. Does the baseline association vary across follow-up intervals?",
				x = "-log10(time-interaction P)", y = NULL
			) +
			theme_5c(9) else blank_plot("b. Association across follow-up intervals", if (!length(anchors))
		"No configured anchor is present in this omics layer" else "No estimable time-interaction test; consult temporal_heterogeneity.csv")
	save_plot(p1 / p2, "c1.Fig16.paired_temporal_validation.png", 16, 10, outdir = if (layer == "protein")
		out.prot else out.met)
}

landmark_incident_scan <- function(dat, xs, tvar, evar, covars, landmarks = C1_LANDMARK_YEARS) {
	xs <- intersect(xs, names(dat))
	d0 <- dat
	for (x in xs) {
		z <- as.numeric(d0[[x]])
		sd0 <- sd(z, na.rm = TRUE)
		d0[[x]] <- if (is.finite(sd0) && sd0 > 0)
			(z - mean(z, na.rm = TRUE)) / sd0 else NA_real_
	}
	map_dfr(landmarks, function(h) {
		d <- d0 |>
			filter(is.finite(.data[[tvar]]), !is.na(.data[[evar]]), .data[[tvar]] > h) |>
			mutate(.landmark_time = .data[[tvar]] - h, .landmark_event = .data[[evar]])
		z <- cox_scan(d, xs, covars, Y, scale_x = FALSE, time_var = ".landmark_time", event_var = ".landmark_event")
		z |>
			mutate(landmark_years = h, N_risk = N_total, events_after_landmark = N_event, exposure_unit = "SD fixed in the baseline omics cohort, not restandardized by landmark")
	})
}

run_c1_layer <- function(layer = c("protein", "metabolite")) {
	if (LE8_REUSE_RESULTS) return(le8_restore_outputs(match.arg(layer), "c1_correlate"))
	# Check reusable results and initialize the analysis output directory.
	layer <- match.arg(layer)
	le8_begin_analysis(layer, "c1_correlate")
	.le8_analysis_env <- environment()
	on.exit(le8_finish_analysis(layer, "c1_correlate", .le8_analysis_env), add = TRUE)

	layer <- match.arg(layer) ; outdir <- if (layer == "protein") out.prot else out.met ; setwd2(outdir)
	rawdir <- le8_job_dir(outdir, LE8_JOB) ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	scan_cache <- file.path(rawdir, "c1.scan.rds") ; selected_cache <- file.path(rawdir, "c1.res.rds")
	# Reuse legacy version-named scans without making version changes a cache miss.
	if (!LE8_REPLACE && !file.exists(scan_cache)) {
		legacy <- list.files(rawdir, pattern = "^c1[.]scan[.].+[.]rds$", full.names = TRUE)
		legacy <- legacy[order(file.info(legacy)$mtime, decreasing = TRUE)]
		for (p in legacy) {
			z <- tryCatch(readRDS(p), error = function(e) NULL)
			if (is.list(z) && all(c("association_adj2", "prevalent_adj2") %in% names(z))) {
				scan_cache <- p ; break
			}
		}
	}
	if (truthy(Sys.getenv("C1_FIG12_ONLY", unset = "FALSE"))) {
		if (!file.exists(selected_cache) || file.size(selected_cache) <= 0)
			stop("C1_FIG12_ONLY requires an existing C1 result: ", selected_cache, call. = FALSE)
		old <- readRDS(selected_cache)
		if (!identical(old$meta$code_version %||% NA_character_, C1_CODE_VERSION))
			stop("C1_FIG12_ONLY cannot relabel a pre-update result; rerun C1 with the current code first.", call. = FALSE)
		needed <- c("pgs_incident", "pgs_prevalent", "pgs_attained_age", "association_adj2", "prevalent_adj2", "birthline_adj2")
		missing <- setdiff(needed, names(old)) ; if (length(missing))
			stop("Existing C1 result lacks fields required for Fig1/2: ", paste(missing, collapse = ", "), call. = FALSE)
		features_all <- old$input_feature_annotation_audit$feature %||%
			unique(c(old$pgs_incident$term, old$association_adj2$term))
		plot_c1_fig12(
			layer, features_all, old$pgs_incident, old$pgs_prevalent,
			old$pgs_attained_age, old$association_adj2, old$prevalent_adj2,
			old$birthline_adj2, outdir
		)
		message("Completed C1 Fig1/2-only redraw: ", outdir)
		return(invisible(old))
	}
	pgs_signature <- c1_pgs_signature(layer)
	if (cache_valid(selected_cache)) {
		old <- tryCatch(readRDS(selected_cache), error = function(e) NULL)
		if (is.list(old) && identical(old$meta$code_version, C1_CODE_VERSION) && identical(old$meta$pgs_signature, pgs_signature) &&
			all(c("pgs_incident", "pgs_prevalent", "pgs_attained_age") %in% names(old))) {
			cache_message(paste0("C1/", layer), selected_cache) ; return(le8_restore_outputs(layer, "c1_correlate"))
		}
		message("C1/", layer, ": cache incomplete; recomputing inherited-score and observed-level outputs")
	}

	biom0 <- if (layer == "protein") read_prot() else read_met() ; features_all <- setdiff(names(biom0), "eid")
	scan <- if (cache_valid(scan_cache)) tryCatch(readRDS(scan_cache), error = function(e) NULL) else NULL
	scan_contract <- list(
		version = "final-landmark-family", landmark_all = Sys.getenv("C1_LANDMARK_ALL", "TRUE"), endpoint_manifest = Sys.getenv("LE8_ENDPOINT_MANIFEST", ""), layer = layer, le4_covars = sort(C1_LE4_COVARS),
		full_le8_sensitivity = C1_FULL_LE8_SENSITIVITY,
		treatment_vars = sort(C1_TREATMENT_VARS)
	)
	reuse <- is.list(scan) && identical(scan$scan_contract, scan_contract) && all(c("association_adj2", "prevalent_adj2") %in% names(scan))
	if (!is.null(scan) && !reuse)
		message("C1/", layer, ": scan cache incomplete; recomputing association scans")

	# If the expensive scans exist, load only the features needed for figures + same-N table.
	if (reuse) {
		a0 <- scan$association_adj2 %||% scan$association_basic ; p0 <- scan$prevalent_adj2 %||% scan$prevalent_basic
		fig_features <- unique(c(
			C1_DIRECTION_ANCHORS, C1_VLDL_DEEP_FEATURES,
			a0$term[is.finite(a0$p.value) & a0$p.value < .05 / nrow(a0)],
			a0 |> filter(is.finite(p.value)) |> slice_min(p.value, n = max(CLUSTER_TOP, ATTENUATION_TOP, TOP_N), with_ties = FALSE) |> pull(term),
			p0 |> filter(is.finite(p.value)) |> slice_min(p.value, n = max(GRADIENT_TOP, TOP_N), with_ties = FALSE) |> pull(term)
		))
		biom <- biom0[, intersect(c("eid", fig_features), names(biom0)), drop = FALSE]
	} else biom <- biom0
	rm(biom0) ; invisible(gc())

	need <- unique(c(
		"eid", "ethnic.c", vars.basic, vars.le8, C1_TREATMENT_VARS,
		"birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_", Y)
	))
	# Restrict C1 to participants with the omics assay.
	# Filter the narrow phenotype table before materializing the wide join.
	dat <- read_all(need) |>
		filter_analysis_cohort() |>
		inner_join(biom, by = "eid") |>
		make_outcome(Y) |>
		add_attained_age_time(Y)
	rm(biom) ; invisible(gc())
	features <- intersect(features_all, names(dat)) ; covs_basic <- intersect(vars.basic, names(dat))
	covs_adj2_full <- intersect(unique(c(vars.adj2, C1_TREATMENT_VARS)), names(dat))
	covs_adj2 <- intersect(if (le8_custom_adjustment()) le8_custom_covars else unique(c(vars.basic, C1_LE4_COVARS, C1_TREATMENT_VARS)), names(dat))
	le8_covariates_removed <- setdiff(covs_adj2_full, covs_adj2)
	c1_covs_use_name <- if (covs_use_name == "basic") "basic" else
		"adj_le4"
	tvar <- paste0(Y, ".t2e") ; evar <- paste0(Y, ".Yt2e") ; bvar <- paste0(Y, ".b2e") ; bivar <- paste0(Y, ".bi2e")
	rvar <- paste0(Y, ".r2e") ; revar <- paste0(Y, ".Yr2e")
	# Attained age is already the analysis time, so an additional linear age
	# covariate would double-adjust age.  Sex and the remaining covariates stay.
	age_covs <- unique(c("age", grep("^age($|[._])", unique(c(covs_basic, covs_adj2)), value = TRUE, ignore.case = TRUE)))
	birth_covs_basic <- setdiff(covs_basic, age_covs) ; birth_covs_adj2 <- setdiff(covs_adj2, age_covs)
	birth_covs_adj2_full <- setdiff(covs_adj2_full, age_covs)
	dat$prevalent_status <- make_prevalent_status(dat, Y)
	vldl_measured_models <- if (layer == "metabolite")
		run_vldl_tg_measured_models(dat, covs_adj2, tvar, evar) else tibble()

	scan_label <- paste0("C1/", layer, if (layer == "protein") " PWAS" else " MWAS")
	scan_started <- le8_stage_start(scan_label, paste0("features=", length(features_all), " cache=", reuse))
	if (!reuse) {
		message("C1/", layer, ": incident Cox scans: basic + behavioral LE4")
		assoc_basic <- cox_scan(dat, features, covs_basic, Y, time_var = tvar, event_var = evar)
		assoc_adj2 <- cox_scan(dat, features, covs_adj2, Y, time_var = tvar, event_var = evar)
		message("C1/", layer, ": attained-age Cox scans with delayed entry at baseline")
		birthline_basic <- cox_scan_delayed_entry(dat, features, birth_covs_basic, Y)
		birthline_adj2 <- cox_scan_delayed_entry(dat, features, birth_covs_adj2, Y)
		message("C1/", layer, ": baseline-prevalent logistic scans: basic + behavioral LE4")
		prevalent_basic <- logistic_scan(dat, features, covs_basic, "prevalent_status")
		prevalent_adj2 <- logistic_scan(dat, features, covs_adj2, "prevalent_status")
		run_full_le8_sensitivity <- length(le8_covariates_removed) > 0 &&
			C1_FULL_LE8_SENSITIVITY
		if (run_full_le8_sensitivity) {
			message("C1/", layer, ": full-LE8 sensitivity (non-primary)")
			assoc_adj2_full_le8 <- cox_scan(dat, features, covs_adj2_full, Y, time_var = tvar, event_var = evar)
			birthline_adj2_full_le8 <- cox_scan_delayed_entry(dat, features, birth_covs_adj2_full, Y)
			prevalent_adj2_full_le8 <- logistic_scan(dat, features, covs_adj2_full, "prevalent_status")
		} else {
			assoc_adj2_full_le8 <- tibble() ; birthline_adj2_full_le8 <- tibble()
			prevalent_adj2_full_le8 <- tibble()
		}
		message("C1/", layer, ": case-only diagnosis-duration and landmark incident scans")
		# Retained at the user's request as an exploratory legacy analysis.  Its
		# cases and controls have incompatible time origins, so it is never used in
		# directionality classification or Final evidence grading.
		reverse_adj2 <- if (all(c(rvar, revar) %in% names(dat)))
			cox_scan(dat, features, covs_adj2, Y, time_var = rvar, event_var = revar) |>
				mutate(analysis_status = "exploratory only: incompatible time origins") else
			tibble(term = features, beta = NA_real_, p.value = NA_real_, FDR = NA_real_, analysis_status = "not estimable")
		duration_adj2 <- prevalent_duration_scan(dat, features, bvar, covs_adj2)
		landmark_features <- unique(c(
			C1_DIRECTION_ANCHORS,
			assoc_adj2 |> filter(is.finite(p.value)) |> slice_min(p.value, n = C1_LANDMARK_TOP, with_ties = FALSE) |> pull(term),
			prevalent_adj2 |> filter(is.finite(p.value)) |> slice_min(p.value, n = C1_LANDMARK_TOP, with_ties = FALSE) |> pull(term)
		))
		if (truthy(Sys.getenv("C1_LANDMARK_ALL", unset = "TRUE"))) landmark_features <- features
		landmark_adj2 <- landmark_incident_scan(dat, landmark_features, tvar, evar, covs_adj2, C1_LANDMARK_YEARS)
		landmark_adj2$FDR_full_family <- p.adjust(landmark_adj2$p.value, "BH", n = length(features) * length(C1_LANDMARK_YEARS))
		risk_features <- unique(c(
			C1_DIRECTION_ANCHORS,
			assoc_adj2 |> filter(is.finite(p.value)) |> slice_min(p.value, n = C1_RISK_TOP, with_ties = FALSE) |> pull(term),
			prevalent_adj2 |> filter(is.finite(p.value)) |> slice_min(p.value, n = C1_RISK_TOP, with_ties = FALSE) |> pull(term)
		))
		risk_window_adj2 <- risk_window_scan(dat, risk_features, tvar, evar, bvar, covs_adj2, C1_RISK_CUTS)
		# Same-N attenuation only needs the strongest features; set C1_ATTENUATION_TOP=0 for all.
		att_features <- assoc_basic |>
			filter(is.finite(p.value)) |>
			arrange(p.value) |>
			pull(term)
		if (ATTENUATION_TOP > 0) att_features <- head(att_features, ATTENUATION_TOP)
		attenuation_sameN <- same_sample_attenuation(dat, att_features, tvar, evar, covs_basic, covs_adj2)
		scan <- list(
			scan_contract = scan_contract,
			association_basic = assoc_basic, association_adj2 = assoc_adj2,
			birthline_basic = birthline_basic, birthline_adj2 = birthline_adj2,
			prevalent_basic = prevalent_basic, prevalent_adj2 = prevalent_adj2, reverse_adj2 = reverse_adj2,
			duration_adj2 = duration_adj2, landmark_adj2 = landmark_adj2, risk_window_adj2 = risk_window_adj2,
			attenuation_sameN = attenuation_sameN,
			association_adj2_full_le8 = assoc_adj2_full_le8,
			birthline_adj2_full_le8 = birthline_adj2_full_le8,
			prevalent_adj2_full_le8 = prevalent_adj2_full_le8
		)
		saveRDS(scan, scan_cache, compress = "xz")
	} else {
		assoc_basic <- scan$association_basic ; assoc_adj2 <- scan$association_adj2
		birthline_basic <- scan$birthline_basic %||% cox_scan_delayed_entry(dat, features, birth_covs_basic, Y)
		birthline_adj2 <- scan$birthline_adj2 %||% cox_scan_delayed_entry(dat, features, birth_covs_adj2, Y)
		prevalent_basic <- scan$prevalent_basic ; prevalent_adj2 <- scan$prevalent_adj2
		reverse_adj2 <- scan$reverse_adj2 %||% tibble(term = features, beta = NA_real_, p.value = NA_real_, FDR = NA_real_)
		duration_adj2 <- scan$duration_adj2 %||% tibble(term = features, beta = NA_real_, p.value = NA_real_, FDR = NA_real_)
		landmark_adj2 <- scan$landmark_adj2 %||% tibble()
		risk_features <- unique(c(
			C1_DIRECTION_ANCHORS,
			assoc_adj2 |> filter(is.finite(p.value)) |> slice_min(p.value, n = C1_RISK_TOP, with_ties = FALSE) |> pull(term),
			prevalent_adj2 |> filter(is.finite(p.value)) |> slice_min(p.value, n = C1_RISK_TOP, with_ties = FALSE) |> pull(term)
		))
		risk_window_adj2 <- scan$risk_window_adj2 %||% risk_window_scan(dat, risk_features, tvar, evar, bvar, covs_adj2, C1_RISK_CUTS)
		attenuation_sameN <- scan$attenuation_sameN %||% tibble()
		assoc_adj2_full_le8 <- scan$association_adj2_full_le8 %||% tibble()
		birthline_adj2_full_le8 <- scan$birthline_adj2_full_le8 %||% tibble()
		prevalent_adj2_full_le8 <- scan$prevalent_adj2_full_le8 %||% tibble()
		message("C1/", layer, ": reuse association scans and regenerate figures")
	}

	le8_stage_done(scan_label, scan_started)

	endpoint_audit <- le8_endpoint_audit(dat, features, covs_adj2, Y, rawdir)
	assoc <- if (c1_covs_use_name == "basic") assoc_basic else assoc_adj2
	assoc_prevalent <- if (c1_covs_use_name == "basic") prevalent_basic else prevalent_adj2
	pgs <- le8_stage(paste0("C1/", layer, " PGS"), run_c1_pgs_scan(layer, features_all, covs_basic, Y, rawdir, overlap_eids = dat$eid), paste0("features=", length(features_all)))
	pgs_incident <- pgs$incident ; pgs_prevalent <- pgs$prevalent ; pgs_attained_age <- pgs$attained_age
	pgs_incident_same <- pgs$incident_same_omic ; pgs_prevalent_same <- pgs$prevalent_same_omic
	pgs_attained_age_same <- pgs$attained_age_same_omic
	pgs_concordance <- build_pgs_actual_concordance(
		list(incident = pgs_incident, prevalent = pgs_prevalent, attained_age = pgs_attained_age),
		list(incident = assoc_adj2, prevalent = prevalent_adj2, attained_age = birthline_adj2),
		list(incident = pgs_incident_same, prevalent = pgs_prevalent_same, attained_age = pgs_attained_age_same)
	)
	vldl_deep <- if (layer == "metabolite") build_vldl_tg_deep_dive(
		pgs_full = list(incident = pgs_incident, prevalent = pgs_prevalent, attained_age = pgs_attained_age),
		pgs_same = list(
			incident = pgs_incident_same, prevalent = pgs_prevalent_same,
			attained_age = pgs_attained_age_same
		),
		measured_le4 = list(incident = assoc_adj2, prevalent = prevalent_adj2, attained_age = birthline_adj2),
		measured_basic = list(incident = assoc_basic, prevalent = prevalent_basic, attained_age = birthline_basic),
		measured_full_le8 = list(
			incident = assoc_adj2_full_le8, prevalent = prevalent_adj2_full_le8,
			attained_age = birthline_adj2_full_le8
		), risk_window_le4 = risk_window_adj2,
		pgs_conditional = pgs$vldl_conditional %||% tibble(),
		measured_conditional = vldl_measured_models
	) else
		list(
			data = tibble(), trajectory = tibble(), pgs_conditional = tibble(),
			measured_conditional = tibble(), figure = NULL
		)
	prefix <- ifelse(layer == "protein", "pwas", "mwas")
	input_feature_audit <- layer_annotation_audit(layer, features_all)
	write_raw_csv(input_feature_audit, "c1.input_feature_annotation_audit.csv", rawdir)
	write_raw_csv(assoc_basic, paste0(prefix, "_incident_basic.csv"), rawdir) ; write_raw_csv(assoc_adj2, paste0(prefix, "_incident_adj2.csv"), rawdir)
	write_raw_csv(birthline_basic, paste0(prefix, "_birthline_attained_age_basic.csv"), rawdir)
	write_raw_csv(birthline_adj2, paste0(prefix, "_birthline_attained_age_adj2.csv"), rawdir)
	write_raw_csv(prevalent_basic, paste0(prefix, "_prevalent_basic.csv"), rawdir) ; write_raw_csv(prevalent_adj2, paste0(prefix, "_prevalent_adj2.csv"), rawdir)
	write_raw_csv(reverse_adj2, paste0(prefix, "_prevalent_reverse_cox_adj2.csv"), rawdir)
	write_raw_csv(duration_adj2, paste0(prefix, "_prevalent_duration_adj2.csv"), rawdir)
	write_raw_csv(landmark_adj2, paste0(prefix, "_incident_landmark_adj2.csv"), rawdir)
	write_raw_csv(risk_window_adj2, paste0(prefix, "_diagnosis_window_riskset_adj2.csv"), rawdir)
	write_raw_csv(attenuation_sameN, "adj2_attenuation_sameN.csv", rawdir)
	write_raw_csv(pgs$status, "c1.pgs_status.csv", rawdir)
	write_raw_csv(pgs_incident, paste0(prefix, "_pgs_incident_full_genetic.csv"), rawdir)
	write_raw_csv(pgs_prevalent, paste0(prefix, "_pgs_prevalent_full_genetic.csv"), rawdir)
	write_raw_csv(pgs_attained_age, paste0(prefix, "_pgs_attained_age_full_genetic.csv"), rawdir)
	write_raw_csv(pgs_incident_same, paste0(prefix, "_pgs_incident_same_omic.csv"), rawdir)
	write_raw_csv(pgs_prevalent_same, paste0(prefix, "_pgs_prevalent_same_omic.csv"), rawdir)
	write_raw_csv(pgs_attained_age_same, paste0(prefix, "_pgs_attained_age_same_omic.csv"), rawdir)
	write_raw_csv(pgs_concordance, "c1.pgs_actual_concordance.csv", rawdir)
	if (layer == "metabolite")
		write_raw_csv(vldl_deep$data, "c1.L_VLDL_TG_pct_deep_dive.csv", rawdir)
	if (layer == "metabolite")
		write_raw_csv(vldl_deep$trajectory, "c1.L_VLDL_TG_pct_riskset_trajectory.csv", rawdir)
	if (layer == "metabolite")
		write_raw_csv(vldl_deep$pgs_conditional, "c1.L_VLDL_TG_pct_pgs_conditional.csv", rawdir)
	if (layer == "metabolite")
		write_raw_csv(vldl_deep$measured_conditional, "c1.L_VLDL_TG_pct_measured_conditional.csv", rawdir)
	write_raw_csv(assoc_adj2_full_le8, paste0(prefix, "_incident_adj2_full_le8_sensitivity.csv"), rawdir)
	write_raw_csv(prevalent_adj2_full_le8, paste0(prefix, "_prevalent_adj2_full_le8_sensitivity.csv"), rawdir)
	write_raw_csv(birthline_adj2_full_le8, paste0(prefix, "_birthline_adj2_full_le8_sensitivity.csv"), rawdir)

	cohort <- tibble(layer,
		N_omics = nrow(dat), incident_events = sum(dat[[evar]] == 1, na.rm = TRUE),
		prevalent_cases = sum(dat$prevalent_status == 1, na.rm = TRUE), features = length(features_all),
		annotation_matched = sum(input_feature_audit$annotation_matched, na.rm = TRUE),
		annotation_unmatched = sum(!input_feature_audit$annotation_matched, na.rm = TRUE),
		attained_age_N = sum(is.finite(dat$.attained_entry) & is.finite(dat$.attained_exit)),
		bi2e_available = sum(is.finite(dat[[bivar]])),
		PGS_matched = length(pgs$score_map), PGS_file_signature = pgs_signature,
		LE8_components = length(intersect(vars.le8, names(dat))), covs_use = c1_covs_use_name,
		LE4_components = length(intersect(C1_LE4_COVARS, names(dat))),
		le8_covariates_removed = paste(le8_covariates_removed, collapse = ";"),
		overlap_covariates_removed = paste(le8_covariates_removed, collapse = ";"),
		treatment_covariates = paste(intersect(C1_TREATMENT_VARS, names(dat)), collapse = ";")
	)
	write_raw_csv(cohort, "c1.cohort.csv", rawdir)

	top <- assoc |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = TOP_N, with_ties = FALSE) |>
		pull(term)
	top6 <- head(top, YY_TOP) ; top_inc <- assoc_adj2 |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = GRADIENT_TOP, with_ties = FALSE) |>
		pull(term)
	top_prev <- prevalent_adj2 |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = GRADIENT_TOP, with_ties = FALSE) |>
		pull(term)
	cluster_features <- assoc |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = CLUSTER_TOP, with_ties = FALSE) |>
		pull(term)

	# Fig1 + Fig2: the third row is an attained-age sensitivity analysis with
	# delayed entry at baseline. The left column is a fixed-at-conception omic
	# PGS in the full genetic cohort; the right column is the measured adult omic
	# level with behavioral-LE4 adjustment. Basic observed-level analyses remain in the
	# workbook and raw files but are no longer displayed here.
	plot_c1_fig12(
		layer, features_all, pgs_incident, pgs_prevalent, pgs_attained_age,
		assoc_adj2, prevalent_adj2, birthline_adj2, outdir
	)

	# Fig3: same participants in basic and behavioral-LE4 residualized YY panels.
	yy_need <- unique(c(bvar, top6, covs_adj2)) ; yy_dat <- dat[complete.cases(dat[, intersect(yy_need, names(dat)), drop = FALSE]) & is.finite(dat[[bvar]]), intersect(yy_need, names(dat)), drop = FALSE]
	omic_name <- ifelse(layer == "protein", "protein", "metabolite")
	p3a <- make_yy_panel(yy_dat, top6, bvar, covs_basic, paste0("a. Yin–Yang ", omic_name, " patterns — basic adjustment"),
		ifelse(layer == "protein", "Mean protein z-score", "Mean metabolite z-score"),
		smooth_lines = FALSE
	)
	adj2_label <- if (le8_custom_adjustment()) paste0("adjusted for ", paste(le8_custom_covars, collapse = ", ")) else "basic + behavioral LE4 adjustment"
	p3b <- make_yy_panel(yy_dat, top6, bvar, covs_adj2, paste0("b. Yin–Yang trajectories — ", adj2_label),
		ifelse(layer == "protein", "Mean protein z-score", "Mean metabolite z-score"),
		smooth_lines = TRUE
	)
	save_plot(p3a / p3b + plot_layout(heights = c(1, 1)), "c1.Fig3.yy_top.png", 16, 13.5, outdir = outdir)

	# Fig4 shows adjusted diagnosis-timed cross-sectional omics gradients.
	p4i <- make_gradient_panel(dat, top_inc, bvar, "incident", GRADIENT_STEP, MIN_BIN_N,
		paste0("a. Incident (Yin): top ", length(top_inc), if (le8_custom_adjustment()) " adjusted biomarkers" else " behavioral-LE4 biomarkers"),
		covars = covs_adj2
	)
	p4p <- make_gradient_panel(dat, top_prev, bvar, "prevalent", GRADIENT_STEP, MIN_BIN_N,
		paste0("b. Baseline-prevalent (Yang): top ", length(top_prev), if (le8_custom_adjustment()) " adjusted biomarkers" else " behavioral-LE4 biomarkers"),
		covars = covs_adj2
	)
	save_plot((p4i / p4p + plot_layout(guides = "collect")) & theme(legend.position = "right"),
		"c1.Fig4.gradient.png", 17, 9.6,
		outdir = outdir
	)
	gradient_rank <- bind_rows(
		assoc_adj2 |> filter(term %in% top_inc) |> transmute(analysis = "incident_behavioral_LE4", feature = term, p_value = p.value) |> arrange(p_value) |> mutate(rank = row_number()),
		prevalent_adj2 |> filter(term %in% top_prev) |> transmute(analysis = "prevalent_behavioral_LE4", feature = term, p_value = p.value) |> arrange(p_value) |> mutate(rank = row_number())
	)
	write_raw_csv(gradient_rank, "c1.gradient_top10_provenance.csv", rawdir)

	# Fig5: separately selected and clustered incident and prevalent adj2 sets.
	cluster_prev <- prevalent_adj2 |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = CLUSTER_TOP, with_ties = FALSE) |>
		pull(term)
	cli <- make_cluster_figure(dat, cluster_features, bvar, CLUSTER_STEP, YY_MAX_YEAR, "incident", "a", fixed_k = NULL)
	clp <- make_cluster_figure(dat, cluster_prev, bvar, CLUSTER_STEP, YY_MAX_YEAR, "prevalent", "b", fixed_k = NULL)
	# Keep the five top-level patchwork rows explicit.  Nesting the three-row
	# cluster_main object here is flattened by patchwork, so a three-entry outer
	# heights vector leaves room for only three of the resulting five panels.
	cluster_figure <- cli$cluster_plot / plot_spacer() / clp$cluster_plot / plot_spacer() /
		(cli$diagnostic_plot | clp$diagnostic_plot) +
		plot_layout(heights = c(1, .10, 1, .05, .38))
	save_plot(cluster_figure,
		"c1.Fig5.gradient_cluster.png", 18, 10.5,
		outdir = outdir
	)
	cl_members <- bind_rows(cli$cluster |> mutate(analysis = "incident_behavioral_LE4"), clp$cluster |> mutate(analysis = "prevalent_behavioral_LE4"))
	cl_metrics <- bind_rows(cli$metrics |> mutate(analysis = "incident_behavioral_LE4"), clp$metrics |> mutate(analysis = "prevalent_behavioral_LE4"))
	write_raw_csv(cl_members, "c1.cluster_membership.csv", rawdir) ; write_raw_csv(cl_metrics, "c1.cluster_selection.csv", rawdir)

	# Fig6: adjusted Q5-vs-Q1 HRs on the same adj2 complete-case sample; no subtitle.
	qplots <- plot_quantile_top(dat, top6, tvar, evar, covs_adj2, paste0("Adjusted for ", paste(covs_adj2, collapse = ", ")))
	save_plot(wrap_plots(qplots, ncol = 3), "c1.Fig6.quantile_top.png", 16, 10, outdir = outdir)

	# Fig8: functional coherence of incident and baseline-prevalent adj2 scans.
	enrich <- run_functional_enrichment(assoc_adj2, layer, features_all)
	enrich_prev <- run_functional_enrichment(prevalent_adj2, layer, features_all)
	n_sig <- sum(is.finite(assoc_adj2$p.value) & assoc_adj2$p.value * nrow(assoc_adj2) < .05)
	n_sig_prev <- sum(is.finite(prevalent_adj2$p.value) & prevalent_adj2$p.value * nrow(prevalent_adj2) < .05)
	write_raw_csv(enrich, "c1.enrichment_incident_sig.csv", rawdir) ; write_raw_csv(enrich_prev, "c1.enrichment_prevalent_sig.csv", rawdir)
	pe_i <- plot_functional_enrichment(enrich, n_sig, assoc_adj2, layer, if (le8_custom_adjustment()) "a. Incident (Yin), selected covariates" else "a. Incident (Yin), behavioral LE4")
	pe_p <- plot_functional_enrichment(enrich_prev, n_sig_prev, prevalent_adj2, layer, if (le8_custom_adjustment()) "b. Baseline-prevalent (Yang), selected covariates" else "b. Baseline-prevalent (Yang), behavioral LE4")
	if (layer != "protein") save_plot((pe_i / pe_p + plot_layout(guides = "collect")) & theme(legend.position = "right"),
		"c1.Fig8.enrich_sig.png", 18, 10.5,
		outdir = outdir
	)

	# Fig9 separates distal antecedent prediction from disease-state evidence.
	# GDF15 and PCSK9 are anchors by default, but the table is generated for all
	# assayed proteins/metabolites.
	directionality <- build_directionality_table(assoc_adj2, prevalent_adj2, duration_adj2, landmark_adj2, birthline_adj2, reverse_adj2)
	write_raw_csv(directionality, "c1.directionality_triage.csv", rawdir)
	save_plot(plot_directionality_triage(directionality, C1_DIRECTION_ANCHORS),
		"c1.Fig9.directionality_triage.png", 16, 9,
		outdir = outdir
	)
	save_plot(plot_landmark_and_attained_age(landmark_adj2, assoc_adj2, birthline_adj2, C1_DIRECTION_ANCHORS),
		"c1.Fig10.landmark_birthline_sensitivity.png", 16, 9,
		outdir = outdir
	)
	save_plot(plot_risk_window_scan(risk_window_adj2, C1_DIRECTION_ANCHORS),
		"c1.Fig11.diagnosis_window_riskset.png", 18, 11,
		outdir = outdir
	)
	save_plot(plot_directionality_supplement(directionality, C1_DIRECTION_ANCHORS),
		"c1.Fig12.directionality_detail.png", 14, 7.5,
		outdir = outdir
	)
	save_plot(plot_reverse_time_exploratory(reverse_adj2, prevalent_adj2, duration_adj2, C1_DIRECTION_ANCHORS),
		"c1.Fig13.reverse_time_exploratory.png", 15.5, 8,
		outdir = outdir
	)
	save_plot(plot_pgs_actual_concordance(pgs_concordance, C1_DIRECTION_ANCHORS),
		"c1.Fig14.pgs_actual_concordance.png", 18, 10.5,
		outdir = outdir
	)
	if (layer == "metabolite") save_plot(vldl_deep$figure,
		"c1.Fig15.L_VLDL_TG_pct_deep_dive.png", 19, 26,
		outdir = outdir
	)

	out <- list(
		meta = module_meta(layer, extra = list(
			N = nrow(dat), events = sum(dat[[evar]] == 1, na.rm = TRUE),
			covs_use = c1_covs_use_name, covariates = covs_adj2,
			le4_covariates = intersect(C1_LE4_COVARS, names(dat)),
			le8_covariates_removed = le8_covariates_removed,
			overlap_covariates_removed = le8_covariates_removed,
			treatment_covariates = intersect(C1_TREATMENT_VARS, names(dat)),
			code_version = C1_CODE_VERSION, pgs_signature = pgs_signature,
			pgs_scan_signature = pgs$signature,
			pgs_matched = length(pgs$score_map), pgs_interpretation = "fixed-at-conception score; not a biomarker measured at birth"
		)),
		cohort = cohort, association = assoc, association_basic = assoc_basic, association_adj2 = assoc_adj2,
		association_LE4 = assoc_adj2, association_LE8 = assoc_adj2_full_le8,
		association_adj2_full_le8_sensitivity = assoc_adj2_full_le8,
		birthline_basic = birthline_basic, birthline_adj2 = birthline_adj2,
		birthline_adj2_full_le8_sensitivity = birthline_adj2_full_le8,
		prevalent = assoc_prevalent, prevalent_basic = prevalent_basic, prevalent_adj2 = prevalent_adj2,
		prevalent_adj2_full_le8_sensitivity = prevalent_adj2_full_le8,
		reverse_prevalent = reverse_adj2, prevalent_duration = duration_adj2, landmark_incident = landmark_adj2,
		diagnosis_window_riskset = risk_window_adj2, directionality = directionality,
		attenuation = attenuation_sameN, enrichment = enrich, enrichment_prevalent = enrich_prev,
		pgs_status = pgs$status, pgs_incident = pgs_incident, pgs_prevalent = pgs_prevalent,
		pgs_attained_age = pgs_attained_age, pgs_incident_same_omic = pgs_incident_same,
		pgs_prevalent_same_omic = pgs_prevalent_same, pgs_attained_age_same_omic = pgs_attained_age_same,
		pgs_actual_concordance = pgs_concordance,
		L_VLDL_TG_pct_deep_dive = vldl_deep$data,
		L_VLDL_TG_pct_riskset_trajectory = vldl_deep$trajectory,
		L_VLDL_TG_pct_pgs_conditional = vldl_deep$pgs_conditional,
		L_VLDL_TG_pct_measured_conditional = vldl_deep$measured_conditional,
		input_feature_annotation_audit = input_feature_audit, top_features = top, gradient_top10_provenance = gradient_rank, clusters = cl_members, cluster_selection = cl_metrics
	)
	if (layer == "protein") {
		out$pwas_incident <- assoc ; out$pwas_prevalent <- assoc_prevalent ; out$top_proteins <- top
	} else {
		out$MWAS <- assoc ; out$MWAS_prevalent <- assoc_prevalent
	}
	write_xlsx2(list(
		cohort = cohort, input_feature_audit = input_feature_audit, association = assoc, prevalent = assoc_prevalent, incident_basic = assoc_basic, incident_adj2 = assoc_adj2,
		birthline_basic = birthline_basic, birthline_adj2 = birthline_adj2,
		prevalent_basic = prevalent_basic, prevalent_adj2 = prevalent_adj2,
		association_adj2_full_le8_sensitivity = assoc_adj2_full_le8,
		prevalent_adj2_full_le8_sensitivity = prevalent_adj2_full_le8,
		birthline_adj2_full_le8_sensitivity = birthline_adj2_full_le8,
		attenuation_sameN = attenuation_sameN,
		reverse_prevalent = reverse_adj2, prevalent_duration = duration_adj2, incident_landmark = landmark_adj2,
		diagnosis_window_riskset = risk_window_adj2, directionality = directionality,
		pgs_status = pgs$status, pgs_incident = pgs_incident, pgs_prevalent = pgs_prevalent,
		pgs_attained_age = pgs_attained_age, pgs_incident_same_omic = pgs_incident_same,
		pgs_prevalent_same_omic = pgs_prevalent_same, pgs_attained_same_omic = pgs_attained_age_same,
		pgs_actual_concordance = pgs_concordance,
		L_VLDL_TG_pct_deep_dive = vldl_deep$data,
		L_VLDL_TG_pct_riskset_trajectory = vldl_deep$trajectory,
		L_VLDL_TG_pct_pgs_conditional = vldl_deep$pgs_conditional,
		L_VLDL_TG_pct_measured_conditional = vldl_deep$measured_conditional,
		gradient_top10_provenance = gradient_rank,
		cluster_membership = cl_members, cluster_selection = cl_metrics,
		enrichment_incident = enrich, enrichment_prevalent = enrich_prev
	), "c1.out.xlsx")
	out$mock <- le8_mock_c1(dat, assoc_adj2, enrich, layer, covs_adj2, tvar, evar, outdir, enrich_prev)
	finalize_outputs(LE8_JOB, outdir) ; saveRDS(out, selected_cache, compress = "xz") ; out
}

if (prot_DO) {
	invisible(le8_stage("C1/protein", run_c1_layer("protein"))) ; gc(full = TRUE)
}
if (met_DO) {
	invisible(le8_stage("C1/metabolite", run_c1_layer("metabolite"))) ; gc(full = TRUE)
}
}
