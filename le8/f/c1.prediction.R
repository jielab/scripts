# C1 conventional prediction: training-only selection, clinical/omics comparisons.
# Statistical settings retain their existing FINAL_* names for run compatibility.
c1_eid_hash <- function(ids) digest::digest(paste0(paste0(sort(unique(as.character(ids))), "\n"), collapse=""), algo="sha256", serialize=FALSE)
c1_lockpaper_protocol <- function() {
    path <- Sys.getenv("LE8_OUTER_ROSTER", "")
    if (!nzchar(path)) stop("C1 LOCKPAPER requires frozen shared outer roster")
    record <- jsonlite::fromJSON(paste0(path, ".protocol.json"))
    if (!identical(record$release, "2026-10-17.c1-lockpaper-v1") ||
        !identical(record$roster_sha256, digest::digest(file=path, algo="sha256")))
        stop("C1 LOCKPAPER outer roster integrity/protocol mismatch")
    if (nzchar(record$group_file) &&
        !identical(record$group_file_sha256, digest::digest(file=record$group_file, algo="sha256")))
        stop("C1 LOCKPAPER group map SHA256 mismatch")
    if (!is.null(record$relatedness_source_file) &&
        !identical(record$relatedness_source_sha256, digest::digest(file=record$relatedness_source_file, algo="sha256")))
        stop("C1 LOCKPAPER relatedness source SHA256 mismatch")
    audit_fields <- c("relatedness_source_file", "relatedness_source_sha256", "group_map_sha256",
        "pair_rows", "related_participant_N", "singleton_participant_N", "component_N",
        "largest_component_size", "grouping_unit", "grouping_source", "family_disjoint_outer_split")
    c(list(release=record$release, outer_roster=normalizePath(path, mustWork=TRUE), outer_roster_sha256=record$roster_sha256),
      record[intersect(audit_fields, names(record))])
}
c1_export_benchmark <- function(pred, layer, rawdir, horizon) {
    local <- file.path(rawdir, "_LOCAL_ONLY")
    dir.create(local, recursive = TRUE, showWarnings = FALSE, mode = "0700")
    path <- file.path(local, "benchmark.json.gz")
    pred$eid <- as.character(pred$eid)
    connection <- gzfile(path, "wt")
    tryCatch(writeLines(jsonlite::toJSON(list(
        format = "le8-c1-benchmark-v1", trait = Y, layer = if (layer == "protein") "prot" else "met",
        horizon = horizon, score_kind = "risk_ranking", predictions = pred,
        lockpaper_protocol = c(c1_lockpaper_protocol(), list(validation_eid_hash=c1_eid_hash(pred$eid))),
        analysis_options = le8_analysis_options(), interpretation = "Frozen outer-test predictions; exploratory internal validation"
    ), dataframe = "rows", auto_unbox = TRUE, digits = NA, na = "null"), connection),
    finally = close(connection))
    Sys.chmod(path, "0600")
}
# C1 prediction: unified score phenotyping and prediction figures.

suppressPackageStartupMessages({
	fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(fdir, "0.common.R")) ;
	# C1 prediction methods
	# Independent prediction functions extracted from the retained LE8 implementation.
	# Shared preprocessing/survival helpers are loaded via 0.common.R.
	# Changing this method must not change the held-out sample or test-set selection rules.
	fit_gaussian_glmnet <- function(tr, te, vars, yvar, nfolds = 5, alpha = .5) {
		if (!requireNamespace("glmnet", quietly = TRUE)) stop("Prediction requires glmnet.", call. = FALSE)
		x <- impute_train_test(tr, te, vars) ; if (ncol(x$tr) < 1) return(NULL) ; y <- suppressWarnings(as.numeric(tr[[yvar]])) ; ok <- is.finite(y)
		if (sum(ok) < 200 || sd(y[ok]) == 0) return(NULL)
		fit <- tryCatch(glmnet::cv.glmnet(x$tr[ok, , drop = FALSE], y[ok], family = "gaussian", alpha = alpha, nfolds = min(nfolds, max(3, floor(sum(ok) / 50))), standardize = FALSE, keep = TRUE), error = function(e) NULL) ; if (is.null(fit)) return(NULL)
		b <- as.matrix(coef(fit, s = "lambda.1se")) ; sel <- rownames(b)[b[, 1] != 0] ; sel <- setdiff(sel, "(Intercept)")
		# Use inner-fold out-of-fold predictions for the outer-training rows. This prevents
		# the second-stage Cox model from being calibrated on in-sample omics scores.
		li <- which.min(abs(fit$lambda - fit$lambda.1se)) ; trp <- rep(NA_real_, nrow(tr))
		if (!is.null(fit$fit.preval) && length(li)) trp[ok] <- as.numeric(fit$fit.preval[, li])
		if (sum(is.finite(trp)) < sum(ok) * .8) stop("Insufficient OOF predictions; in-sample substitution is forbidden")
		list(tr = trp, te = as.numeric(predict(fit, x$te, s = "lambda.1se")), selected = sel)
	}

	fit_cox_glmnet <- function(tr, te, vars, tvar, evar, nfolds = 5, alpha = .5, penalty_factor = NULL) {
		if (!requireNamespace("glmnet", quietly = TRUE)) stop("Prediction requires glmnet.", call. = FALSE)
		x <- impute_train_test(tr, te, vars) ; if (ncol(x$tr) < 2) return(NULL) ; time <- tr[[tvar]] ; event <- tr[[evar]] ; ok <- is.finite(time) & is.finite(event) & time > 0
		if (sum(ok) < 300 || sum(event[ok] == 1) < 30) return(NULL)
		pf <- rep(1, ncol(x$tr)) ; names(pf) <- colnames(x$tr) ; if (!is.null(penalty_factor)) {
			hit <- intersect(names(penalty_factor), names(pf)) ; pf[hit] <- penalty_factor[hit]
		}
		fit <- tryCatch(glmnet::cv.glmnet(x$tr[ok, , drop = FALSE], Surv(time[ok], event[ok]), family = "cox", alpha = alpha, nfolds = min(nfolds, max(3, floor(sum(event[ok]) / 10))), standardize = FALSE, penalty.factor = pf, keep = TRUE), error = function(e) NULL) ; if (is.null(fit)) return(NULL)
		b <- as.matrix(coef(fit, s = "lambda.1se")) ; sel <- rownames(b)[b[, 1] != 0]
		li <- which.min(abs(fit$lambda - fit$lambda.1se)) ; trp <- rep(NA_real_, nrow(tr))
		if (!is.null(fit$fit.preval) && length(li)) trp[ok] <- as.numeric(fit$fit.preval[, li])
		if (sum(is.finite(trp)) < sum(ok) * .8) stop("Insufficient OOF predictions; in-sample substitution is forbidden")
		list(tr = trp, te = as.numeric(predict(fit, x$te, s = "lambda.1se", type = "link")), selected = sel, beta = b[sel, 1])
	}


	pradeep_population_scores <- function(train, target, biom_lists, tvar, evar,
					method = "binomial", inner = 10, seed = 2026) {
		if (!requireNamespace("glmnet", quietly = TRUE)) stop("Prediction requires glmnet.", call. = FALSE)
		out <- vector("list", length(biom_lists)) ; names(out) <- names(biom_lists)
		for (i in seq_along(biom_lists)) {
			nm <- names(biom_lists)[i] ; vv <- intersect(biom_lists[[i]], intersect(names(train), names(target)))
			vv <- vv[vapply(train[, vv, drop = FALSE], function(z) {
				z <- suppressWarnings(as.numeric(z)) ; is.finite(stats::var(z, na.rm = TRUE)) && stats::var(z, na.rm = TRUE) > 0
			}, logical(1))]
			if (!length(vv)) next
			xtr <- as.matrix(data.frame(lapply(train[, vv, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))
			xall <- as.matrix(data.frame(lapply(target[, vv, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))
			storage.mode(xtr) <- storage.mode(xall) <- "double"
			med <- apply(xtr, 2, function(z) {
				m <- median(z[is.finite(z)], na.rm = TRUE) ; ifelse(is.finite(m), m, 0)
			})
			for (j in seq_along(med)) {
				xtr[!is.finite(xtr[, j]), j] <- med[j] ; xall[!is.finite(xall[, j]), j] <- med[j]
			}
			y <- switch(method,
				gaussian = as.numeric(train[[evar]]),
				binomial = as.integer(train[[evar]] == 1),
				cox = survival::Surv(train[[tvar]], train[[evar]])
			)
			ok <- is.finite(train[[tvar]]) & !is.na(train[[evar]])
			if (method == "binomial" && length(unique(y[ok])) < 2L) next
			set.seed(seed + i)
			nf <- if (method == "binomial") max(3L, min(inner, min(table(y[ok])))) else inner
			fit <- glmnet::cv.glmnet(xtr[ok, , drop = FALSE], y[ok],
				family = method, alpha = 1, nfolds = nf,
				type.measure = ifelse(method == "binomial", "auc", "deviance")
			)
			bb <- as.matrix(coef(fit, s = "lambda.1se")) ; nsel <- sum(bb[ - 1, 1] != 0)
			score <- as.numeric(predict(fit, newx = xall, s = "lambda.1se", type = ifelse(method == "binomial", "link", "link")))
			out[[i]] <- tibble::tibble(eid = target$eid, biom_set = nm, model = "Biomarkers", score = score, n_selected = nsel)
		}
		dplyr::bind_rows(out)
	}


	# C1 prediction
	# Independent prediction functions extracted from the retained LE8 implementation.
	# Shared preprocessing/survival helpers are loaded via 0.common.R.
	# Changing this method must not change the held-out sample or test-set selection rules.
	split_yu_prediction <- function(dat, biom_vars, basic, tvar, evar, horizon = 10,
					test_frac = .20, seed = 2026,
					checkpoint_file = NULL, force = FALSE) {
		if (!requireNamespace("lightgbm", quietly = TRUE))
			stop("Yu-fair prediction requires the R package 'lightgbm' (included in environment.yml).", call. = FALSE)
		if (!force && !is.null(checkpoint_file) && file.exists(checkpoint_file) && file.size(checkpoint_file) > 0)
			return(readRDS(checkpoint_file))
		biom_vars <- intersect(biom_vars, names(dat)) ; basic <- intersect(basic, names(dat))
		if (!length(biom_vars)) return(tibble::tibble())
		set.seed(seed) ; test_id <- sample(seq_len(nrow(dat)), max(1L, round(nrow(dat) * test_frac)))
		tr_id <- setdiff(seq_len(nrow(dat)), test_id)
		known <- (dat[[evar]][tr_id] == 1 & dat[[tvar]][tr_id] <= horizon) | dat[[tvar]][tr_id] > horizon
		tr_id <- tr_id[which(known)] ; y <- as.integer(dat[[evar]][tr_id] == 1 & dat[[tvar]][tr_id] <= horizon)
		if (length(unique(y)) < 2L) stop("Yu-fair training set has fewer than two horizon-status classes.", call. = FALSE)

		numeric_matrix <- function(rows, vars) {
			x <- as.matrix(data.frame(lapply(dat[rows, vars, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE))
			storage.mode(x) <- "double" ; x
		}
		xb_tr <- numeric_matrix(tr_id, biom_vars) ; xb_te <- numeric_matrix(test_id, biom_vars)
		clinical_matrices <- function() {
			if (!length(basic)) return(list(tr = matrix(nrow = length(tr_id), ncol = 0), te = matrix(nrow = length(test_id), ncol = 0)))
			# Clinical predictors are numeric scores/covariates. Construct the two
			# matrices separately so missing values can never make model.matrix drop
			# participants and desynchronise rows from tr_id/test_id.
			a <- numeric_matrix(tr_id, basic) ; b <- numeric_matrix(test_id, basic)
			med <- apply(a, 2, function(z) {
				v <- median(z[is.finite(z)], na.rm = TRUE) ; ifelse(is.finite(v), v, 0)
			})
			for (j in seq_along(med)) {
				a[!is.finite(a[, j]), j] <- med[j] ; b[!is.finite(b[, j]), j] <- med[j]
			}
			if (nrow(a) != length(tr_id) || nrow(b) != length(test_id)) stop("Yu-fair clinical matrix row contract failed.", call. = FALSE)
			list(tr = a, te = b)
		}
		xc <- clinical_matrices()
		params <- list(
			objective = "binary", metric = "auc", max_depth = 15L, num_leaves = 10L,
			bagging_fraction = .70, bagging_freq = 1L, learning_rate = .01,
			feature_fraction = .70, max_bin = 63L, verbosity =  - 1L,
			num_threads = max(1L, as.integer(Sys.getenv("N_CORES", unset = "4"))), seed = as.integer(seed)
		)
		fit_one <- function(xtr, xte, model_name) {
			ds <- lightgbm::lgb.Dataset(data = xtr, label = y)
			fit <- lightgbm::lgb.train(params = params, data = ds, nrounds = 500L, verbose =  - 1)
			score <- as.numeric(predict(fit, xte)) ; imp <- tryCatch(lightgbm::lgb.importance(fit), error = function(e) NULL)
			nsel <- if (is.null(imp)) ncol(xtr) else sum(imp$Gain > 0, na.rm = TRUE)
			tibble::tibble(
				eid = dat$eid[test_id], biom_set = "Yu / LightGBM", model = model_name,
				time = dat[[tvar]][test_id], event = dat[[evar]][test_id], score = score,
				cindex = fcidx(Surv(dat[[tvar]][test_id], dat[[evar]][test_id]), score), n_selected = nsel
			)
		}
		ans <- dplyr::bind_rows(
			fit_one(xb_tr, xb_te, "Biomarkers"),
			fit_one(cbind(xb_tr, xc$tr), cbind(xb_te, xc$te), "Combined")
		)
		if (!is.null(checkpoint_file)) {
			dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE) ; saveRDS(ans, checkpoint_file, compress = FALSE)
		}
		ans
	}


	# final.R

	# Method-specific functions live in independent files for future development.
	.le8_method_dir <- Sys.getenv("LE8_FDIR", file.path(Sys.getenv("DIRSCRIPT"), "f"))
	# Final-specific nested cross-validation and omics prediction utilities.
	# Generic prediction helpers (including fcidx and RE_pred) come from
	# /mnt/d/scripts/0f/prediction.R, sourced by 0.common.R.

	impute_train_test <- function(tr, te, vars) {
		vars <- intersect(vars, intersect(names(tr), names(te))) ; if (!length(vars)) return(list(tr = matrix(nrow = nrow(tr), ncol = 0), te = matrix(nrow = nrow(te), ncol = 0), vars = character()))
		a <- tr[, vars, drop = FALSE] ; b <- te[, vars, drop = FALSE]
		a[] <- lapply(a, function(x) suppressWarnings(as.numeric(x))) ; b[] <- lapply(b, function(x) suppressWarnings(as.numeric(x)))
		keep <- vapply(a, function(x) sum(is.finite(x)) >= 100 && is.finite(sd(x, na.rm = TRUE)) && sd(x, na.rm = TRUE) > 0, logical(1)) ; a <- a[, keep, drop = FALSE] ; b <- b[, keep, drop = FALSE]
		if (!ncol(a)) return(list(tr = matrix(nrow = nrow(tr), ncol = 0), te = matrix(nrow = nrow(te), ncol = 0), vars = character()))
		cen <- vapply(a, function(x) median(x[is.finite(x)], na.rm = TRUE), numeric(1)) ; sca <- vapply(a, function(x) sd(x, na.rm = TRUE), numeric(1)) ; sca[!is.finite(sca) | sca == 0] <- 1
		# `[` keeps a one-column tibble as a list-backed data frame, so
		# `is.finite(a[, j])` fails when upstream joins return a tibble. Work on
		# the underlying vectors to support both base data.frames and tibbles.
		for (j in seq_along(cen)) {
			a[[j]][!is.finite(a[[j]])] <- cen[j]
			b[[j]][!is.finite(b[[j]])] <- cen[j]
		}
		list(
			tr = scale(as.matrix(a), center = cen, scale = sca), te = scale(as.matrix(b), center = cen, scale = sca),
			vars = colnames(a), center = cen, scale = sca
		)
	}

	auc_binary <- function(y, score) {
		id <- which(is.finite(y) & is.finite(score)) ; y <- y[id] ; score <- score[id] ; if (length(y) < 50 || length(unique(y)) < 2) return(NA_real_)
		r <- rank(score) ; n1 <- sum(y == 1) ; n0 <- sum(y == 0) ; (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
	}

	ipcw_at_horizon <- function(time, event, horizon) {
		ok <- is.finite(time) & is.finite(event) & time > 0 ; time <- time[ok] ; event <- event[ok]
		sf <- tryCatch(survival::survfit(Surv(time, 1 - event) ~ 1), error = function(e) NULL)
		if (is.null(sf)) return(list(ok = ok, weight = rep(NA_real_, length(time)), case = rep(FALSE, length(time)), control = rep(FALSE, length(time))))
		Gfun <- stats::stepfun(sf$time, c(1, sf$surv), right = TRUE)
		Gt <- pmax(Gfun(horizon), .02) ; case <- event == 1 & time <= horizon ; control <- time > horizon
		w <- rep(0, length(time)) ; w[case] <- 1 / pmax(Gfun(pmax(time[case] - 1e-8, 0)), .02) ; w[control] <- 1 / Gt
		list(ok = ok, weight = w, case = case, control = control)
	}
	weighted_time_auc <- function(time, event, score, horizon) {
		id <- which(is.finite(time) & is.finite(event) & is.finite(score) & time > 0) ; if (length(id) < 100) return(NA_real_)
		ip <- ipcw_at_horizon(time[id], event[id], horizon) ; d <- data.frame(score = score[id], cw = ifelse(ip$case, ip$weight, 0), tw = ifelse(ip$control, ip$weight, 0))
		d <- d[d$cw > 0 | d$tw > 0, , drop = FALSE] ; if (sum(d$cw) == 0 || sum(d$tw) == 0) return(NA_real_)
		g <- dplyr::as_tibble(d) |>
			dplyr::group_by(score) |>
			dplyr::summarise(cw = sum(cw), tw = sum(tw), .groups = "drop") |>
			dplyr::arrange(score) |>
			dplyr::mutate(tw_before = dplyr::lag(cumsum(tw), default = 0))
		sum(g$cw * (g$tw_before + .5 * g$tw)) / (sum(g$cw) * sum(g$tw))
	}
	ipcw_brier <- function(time, event, risk, horizon) {
		id <- which(is.finite(time) & is.finite(event) & is.finite(risk) & time > 0) ; if (length(id) < 100) return(NA_real_)
		ip <- ipcw_at_horizon(time[id], event[id], horizon) ; target <- as.numeric(ip$control) ; w <- ip$weight ; sum(w * (target - (1 - risk[id])) ^ 2) / sum(w)
	}

	calibration_metrics <- function(y, risk) {
		id <- which(is.finite(y) & is.finite(risk)) ; if (length(id) < 100 || length(unique(y[id])) < 2) return(c(intercept = NA_real_, slope = NA_real_))
		lp <- qlogis(pmin(pmax(risk[id], 1e-6), 1 - 1e-6)) ; slope <- tryCatch(coef(glm(y[id] ~ lp, family = binomial()))[[2]], error = function(e) NA_real_)
		intercept <- tryCatch(coef(glm(y[id] ~ offset(lp), family = binomial()))[[1]], error = function(e) NA_real_) ; c(intercept = intercept, slope = slope)
	}


	discover_ys_training <- function(tr,biom_vars,le8_vars,basic_vars,block=100,fdr=.05,specificity_cut=.35,max_per_component=40,seed=2026) {
		biom_vars<-intersect(biom_vars,names(tr));le8_vars<-intersect(le8_vars,names(tr));basic_vars<-intersect(basic_vars,names(tr))
		group<-if(".group" %in% names(tr)) tr$.group else if("eid" %in% names(tr)) le8_participant_groups(tr) else seq_len(nrow(tr))
		fold<-le8_group_folds(group,2,seed)
		a<-le8_proxy_map(tr[fold==1,,drop=FALSE],biom_vars,le8_vars,basic_vars)
		b<-le8_proxy_map(tr[fold==2,,drop=FALSE],biom_vars,le8_vars,basic_vars)
		if(!nrow(a)||!nrow(b)) return(list())
		p<-inner_join(a |> select(feature,component,r_disc=r,FDR_disc=FDR),b |> select(feature,component,r_rep=r,FDR_rep=FDR),by=c("feature","component")) |>
			mutate(strict=is.finite(FDR_disc)&is.finite(FDR_rep)&FDR_disc<fdr&FDR_rep<fdr&sign(r_disc)==sign(r_rep),strength=pmin(abs(r_disc),abs(r_rep)),specificity=NA_real_)
		chosen<-p |> filter(strict) |> group_by(component) |> arrange(desc(strength),feature,.by_group=TRUE) |> slice_head(n=max_per_component) |> ungroup() |>
			mutate(component=sub("[.]pts$","",component),proxy_level="YS_strict_fold")
		ans<-setNames(lapply(sub("[.]pts$","",le8_vars),function(cmp) chosen$feature[chosen$component==cmp]),sub("[.]pts$","",le8_vars))
		attr(ans,"selection_info")<-chosen;ans
	}

	cox_predict <- function(tr, te, vars, tvar, evar, horizon = 10) {
		vars <- intersect(vars, intersect(names(tr), names(te))) ; if (!length(vars)) return(list(lp = rep(NA_real_, nrow(te)), risk = rep(NA_real_, nrow(te))))
		f <- as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(vars), collapse = " + ")))
		dtr <- tr[, c(tvar, evar, vars), drop = FALSE] ; dtr <- dtr[complete.cases(dtr), , drop = FALSE] ; if (nrow(dtr) < 200 || sum(dtr[[evar]] == 1) < 20) return(list(lp = rep(NA_real_, nrow(te)), risk = rep(NA_real_, nrow(te))))
		fit <- tryCatch(coxph(f, dtr, x = TRUE), error = function(e) NULL) ; if (is.null(fit)) return(list(lp = rep(NA_real_, nrow(te)), risk = rep(NA_real_, nrow(te))))
		id <- complete.cases(te[, vars, drop = FALSE]) ; lp <- risk <- rep(NA_real_, nrow(te)) ; lp[id] <- tryCatch(as.numeric(predict(fit, newdata = te[id, , drop = FALSE], type = "lp")), error = function(e) NA_real_)
		bh <- basehaz(fit, centered = TRUE)
		H0 <- if (!nrow(bh) || horizon < min(bh$time, na.rm = TRUE)) 0 else bh$hazard[max(which(bh$time <= horizon))]
		if (!is.finite(H0)) H0 <- tail(bh$hazard[is.finite(bh$hazard)], 1)
		risk[id] <- 1 - exp( - H0 * exp(lp[id])) ; list(lp = lp, risk = pmin(pmax(risk, 0), 1), fit = fit)
	}

	nested_cv_omics <- function(dat, biom_vars, le8_vars, basic_vars, prs_vars, tvar, evar, k = 5, inner = 5, horizon = 10, seed = 2026) {
		dat$fold_id <- make_folds(dat, evar, k, seed) ; pred_rows <- list() ; sel_rows <- list() ; ii <- 0L
		for (fd in seq_len(k)) {
			tr <- dat[dat$fold_id != fd, , drop = FALSE] ; te <- dat[dat$fold_id == fd, , drop = FALSE] ; idx <- which(dat$fold_id == fd)
			ys_sets <- discover_ys_training(tr, biom_vars, le8_vars, basic_vars, seed = seed + fd) ; ys_info <- attr(ys_sets, "selection_info") %||% tibble::tibble() ; ys_selected <- unique(unlist(ys_sets)) ; ys_tr <- ys_te <- data.frame(row.names = seq_len(nrow(tr))) ; ys_te <- data.frame(row.names = seq_len(nrow(te)))
			for (cmp in names(ys_sets)) {
				v <- ys_sets[[cmp]] ; fit <- fit_gaussian_glmnet(tr, te, v, paste0(cmp, ".pts"), inner) ; vn <- paste0("YS_", cmp) ; ys_tr[[vn]] <- if (is.null(fit)) NA_real_ else fit$tr ; ys_te[[vn]] <- if (is.null(fit)) NA_real_ else fit$te
				if (!is.null(fit) && length(fit$selected)) {
					si <- tibble::tibble(fold = fd, set = "YS", component = cmp, feature = fit$selected) |> dplyr::left_join(ys_info |> dplyr::select(feature, component, proxy_level), by = c("feature", "component")) ; sel_rows[[length(sel_rows) + 1]] <- si
				}
			}
			tr2 <- dplyr::bind_cols(tr, ys_tr) ; te2 <- dplyr::bind_cols(te, ys_te) ; ys_score_vars <- grep("^YS_", names(tr2), value = TRUE) ; ys_score_vars <- ys_score_vars[vapply(tr2[ys_score_vars], function(x) sum(is.finite(x)) >= 200, logical(1))]
			ns <- fit_cox_glmnet(tr, te, biom_vars, tvar, evar, inner) ; tr2$NS_score <- if (is.null(ns)) NA_real_ else ns$tr ; te2$NS_score <- if (is.null(ns)) NA_real_ else ns$te ; if (!is.null(ns) && length(ns$selected)) sel_rows[[length(sel_rows) + 1]] <- data.frame(fold = fd, set = "NS", component = NA, feature = ns$selected)
			plus_pool <- setdiff(biom_vars, ys_selected) ; plus <- fit_cox_glmnet(tr, te, plus_pool, tvar, evar, inner) ; plus_sel <- if (is.null(plus)) character() else plus$selected
			union <- unique(c(ys_selected, plus_sel)) ; pf <- setNames(rep(1, length(union)), union) ; pf[ys_selected] <- .6 ; ysp <- fit_cox_glmnet(tr, te, union, tvar, evar, inner, penalty_factor = pf) ; tr2$YSP_score <- if (is.null(ysp)) NA_real_ else ysp$tr ; te2$YSP_score <- if (is.null(ysp)) NA_real_ else ysp$te
			ysp_plus_selected <- if (is.null(ysp)) character() else intersect(ysp$selected, plus_sel)
			if (length(ysp_plus_selected)) sel_rows[[length(sel_rows) + 1]] <- data.frame(fold = fd, set = "YSP_plus", component = NA, feature = ysp_plus_selected)
			models <- list(Base = basic_vars, PRS = c(basic_vars, prs_vars), LE8 = c(basic_vars, le8_vars), YS = c(basic_vars, ys_score_vars), LE8_YS = c(basic_vars, le8_vars, ys_score_vars), PRS_LE8 = c(basic_vars, prs_vars, le8_vars), PRS_LE8_YS = c(basic_vars, prs_vars, le8_vars, ys_score_vars), PRS_LE8_YSP = c(basic_vars, prs_vars, le8_vars, "YSP_score"), PRS_LE8_NS = c(basic_vars, prs_vars, le8_vars, "NS_score"))
			for (nm in names(models)) {
				z <- cox_predict(tr2, te2, models[[nm]], tvar, evar, horizon) ; ii <- ii + 1 ; pred_rows[[ii]] <- data.frame(row_id = idx, fold = fd, model = nm, lp = z$lp, risk = z$risk, time = te[[tvar]], event = te[[evar]])
			}
		}
		pred <- dplyr::bind_rows(pred_rows) ; selected <- dplyr::bind_rows(sel_rows)
		if (!nrow(selected)) selected <- tibble::tibble(fold = integer(), set = character(), component = character(), feature = character(), proxy_level = character())
		metrics <- pred |>
			dplyr::group_by(fold, model) |>
			dplyr::group_modify(function(d, key) {
				eligible <- d$time >= horizon | (d$event == 1 & d$time <= horizon) ; y10 <- as.integer(d$event == 1 & d$time <= horizon) ; cal <- calibration_metrics(y10[eligible], d$risk[eligible]) ; tibble::tibble(cindex = fcidx(Surv(d$time, d$event), d$lp), AUC10 = weighted_time_auc(d$time, d$event, d$risk, horizon), Brier10 = ipcw_brier(d$time, d$event, d$risk, horizon), cal_intercept = cal[["intercept"]], cal_slope = cal[["slope"]])
			}) |>
			dplyr::ungroup()
		summary <- metrics |>
			dplyr::group_by(model) |>
			dplyr::summarise(dplyr::across(c(cindex, AUC10, Brier10, cal_intercept, cal_slope), list(mean =  ~ mean(.x, na.rm = TRUE), sd =  ~ sd(.x, na.rm = TRUE))), .groups = "drop")
		stability <- selected |>
			dplyr::count(set, component, feature, name = "fold_count") |>
			dplyr::mutate(selection_frequency = fold_count / k) |>
			dplyr::arrange(set, component, dplyr::desc(selection_frequency))
		list(pred = pred, metrics = metrics, summary = summary, selected = selected, stability = stability, folds = dat$fold_id)
	}

	risk_deciles <- function(pred, model) {
		d <- pred |>
			dplyr::filter(.data$model == .env$model, is.finite(lp), is.finite(time), is.finite(event), time > 0) |>
			dplyr::mutate(decile = dplyr::ntile(lp, 10), decile = factor(decile))
		fit <- tryCatch(coxph(Surv(time, event) ~ decile, d), error = function(e) NULL) ; hr <- if (is.null(fit)) tibble::tibble() else {
			sm <- coef(summary(fit)) ; tibble::tibble(term = rownames(sm), HR = exp(sm[, "coef"]), lo = exp(sm[, "coef"] - 1.96 * sm[, "se(coef)"]), hi = exp(sm[, "coef"] + 1.96 * sm[, "se(coef)"]), p = sm[, "Pr(>|z|)"])
		}
		cum <- d |>
			dplyr::group_by(decile) |>
			dplyr::summarise(N = dplyr::n(), events = sum(event), risk10 = mean(risk, na.rm = TRUE), .groups = "drop") ; list(rows = d, HR = hr, cum = cum)
	}

	# Nature Aging-style sequential panel: C1-significant biomarkers are ordered by
	# Wald P value, added one at a time, and compared with a paired DeLong test on a
	# held-out subset. Stop after two consecutive non-significant AUC gains.
	sequential_forward_panel <- function(dat, ranked, tvar, evar, horizon = 10, seed = 2026,
					alpha = .05, max_extended = 30) {
		ranked <- intersect(unique(ranked), names(dat)) ; extended <- head(ranked, max_extended)
		empty_log <- tibble::tibble(step = integer(), feature = character(), auc = numeric(), delta_auc = numeric(), delong_p = numeric(), accepted = logical())
		if (!length(extended) || !requireNamespace("pROC", quietly = TRUE))
			return(list(parsimonious = head(extended, min(5, length(extended))), extended = extended, log = empty_log))
		set.seed(seed) ; test_id <- sample(seq_len(nrow(dat)), max(1L, round(.2 * nrow(dat))))
		eligible <- is.finite(dat[[tvar]]) & !is.na(dat[[evar]]) & (dat[[tvar]] >= horizon | dat[[evar]] == 1)
		tr <- dat[ - test_id, , drop = FALSE] ; te <- dat[test_id, , drop = FALSE]
		tr <- tr[eligible[ - test_id], , drop = FALSE] ; te <- te[eligible[test_id], , drop = FALSE]
		ytr <- as.integer(tr[[evar]] == 1 & tr[[tvar]] <= horizon) ; yte <- as.integer(te[[evar]] == 1 & te[[tvar]] <= horizon)
		if (length(unique(ytr)) < 2 || length(unique(yte)) < 2) return(list(parsimonious = head(extended, min(5, length(extended))), extended = extended, log = empty_log))
		med <- vapply(tr[, extended, drop = FALSE], function(x) median(as.numeric(x), na.rm = TRUE), numeric(1))
		prep <- function(d) {
			x <- as.data.frame(lapply(d[, extended, drop = FALSE], as.numeric)) ; names(x) <- extended ; for (v in extended) x[[v]][!is.finite(x[[v]])] <- med[[v]] ; x
		}
		xtr <- prep(tr) ; xte <- prep(te) ; old_roc <- NULL ; last_sig <- 0L ; nonsig <- 0L ; stop_at <- NA_integer_ ; rows <- list()
		for (i in seq_along(extended)) {
			use <- extended[seq_len(i)] ; dd <- data.frame(y = ytr, xtr[, use, drop = FALSE]) ; fit <- tryCatch(glm(y ~ ., dd, family = binomial()), error = function(e) NULL)
			score <- if (is.null(fit)) rep(NA_real_, nrow(xte)) else tryCatch(as.numeric(predict(fit, newdata = xte[, use, drop = FALSE], type = "response")), error = function(e) rep(NA_real_, nrow(xte)))
			roc <- tryCatch(pROC::roc(yte, score, quiet = TRUE, direction = "<"), error = function(e) NULL) ; auc <- if (is.null(roc)) NA_real_ else as.numeric(roc$auc)
			dp <- if (is.null(old_roc) || is.null(roc)) NA_real_ else tryCatch(as.numeric(pROC::roc.test(old_roc, roc, paired = TRUE, method = "delong")$p.value), error = function(e) NA_real_)
			gain <- if (is.null(old_roc)) TRUE else is.finite(dp) && dp < alpha && auc > as.numeric(old_roc$auc)
			if (gain) {
				last_sig <- i ; nonsig <- 0L
			} else nonsig <- nonsig + 1L
			rows[[i]] <- tibble::tibble(step = i, feature = extended[[i]], auc = auc, delta_auc = if (is.null(old_roc)) NA_real_ else auc - as.numeric(old_roc$auc), delong_p = dp, accepted = gain)
			old_roc <- roc
			if (nonsig >= 2L) {
				stop_at <- i ; break
			}
		}
		if (last_sig < 1L) last_sig <- min(1L, length(extended))
		# The evaluated prefix at the stopping rule is the parsimonious panel.
		panel_n <- if (is.finite(stop_at)) stop_at else last_sig
		list(parsimonious = extended[seq_len(panel_n)], extended = extended, log = dplyr::bind_rows(rows))
	}

	# Pradeep-style evidence-set prediction: outer 80/20 splits and 10-fold CV for lambda.
	evidence_biom_lists <- function(outdir, biom_vars, layer, pred.prot.use = NULL, pred.met.use = NULL) {
		rr <- function(job) {
		if (job != "c1_correlate") return(list(data = list(), available = FALSE, status = "outside C1 baseline"))
			f <- file.path(le8_job_dir(outdir, job), paste0(sub("_.*$", "", job), ".res.rds")) ; if (file.exists(f)) readRDS(f) else list()
		}
		r1 <- rr("c1_correlate") ; r2 <- rr("c2_cause")
		c1 <- r1[["association"]] %||% tibble::tibble()
		c2 <- r2[["MR"]] %||% tibble::tibble()
		pwas <- if (nrow(c1) && "term" %in% names(c1)) {
			q <- if ("FDR" %in% names(c1)) c1$FDR else p.adjust(c1$p.value, "BH")
			unique(as.character(c1$term[is.finite(q) & q < .05]))
		} else character()
		mr <- if (nrow(c2) && "exposure" %in% names(c2)) {
			q <- if ("FDR_all" %in% names(c2)) c2$FDR_all else p.adjust(c2$pval, "BH")
			unique(as.character(c2$exposure[is.finite(q) & q < .05]))
		} else character()
		user <- if (layer == "protein") pred.prot.use else pred.met.use
		z <- list(All = biom_vars, `PWAS/MWAS significant` = base::intersect(pwas, biom_vars), `MR significant` = base::intersect(mr, biom_vars))
		if (!is.null(user) && length(base::intersect(user, biom_vars))) z$`User specified` <- base::intersect(user, biom_vars)
		empty <- names(z)[lengths(z) == 0L]
		if (length(empty)) warning("No biomarkers found for: ", paste(empty, collapse = ", "), call. = FALSE)
		z
	}

	split_re_prediction <- function(dat, biom_lists, basic, tvar, evar, method = "cox", test_frac = .20, inner = 10, seed = 2026, checkpoint_dir = NULL, force = FALSE) {
		set.seed(seed)
		test <- sample(seq_len(nrow(dat)), max(1L, round(nrow(dat) * test_frac)))
		dat$fold_id <- 0L ; dat$fold_id[test] <- 1L
		configs <- list(`Model 0` = list(fx = "RE", method = method, varX = basic, vars.basic = character(), opt = list(inner_cv_n = inner, lambda_rule = "lambda.min", alpha = 1)))
		for (sn in names(biom_lists)) if (length(biom_lists[[sn]])) {
			configs[[paste(sn, "Biomarkers", sep = "||")]] <- list(fx = "RE", method = method, varX = biom_lists[[sn]], vars.basic = character(), opt = list(inner_cv_n = inner, lambda_rule = "lambda.1se", alpha = 1))
			configs[[paste(sn, "Combined", sep = "||")]] <- list(fx = "RE", method = method, varX = biom_lists[[sn]], vars.basic = basic, opt = list(inner_cv_n = inner, lambda_rule = "lambda.1se", alpha = 1))
		}
		if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
		zz <- vector("list", length(configs))
		predict_single <- function(nm, var) {
			tr <- dat[dat$fold_id != 1L, , drop = FALSE] ; te <- dat[dat$fold_id == 1L, , drop = FALSE]
			med <- median(suppressWarnings(as.numeric(tr[[var]])), na.rm = TRUE) ; if (!is.finite(med)) med <- 0
			xtr <- suppressWarnings(as.numeric(tr[[var]])) ; xte <- suppressWarnings(as.numeric(te[[var]]))
			xtr[!is.finite(xtr)] <- med ; xte[!is.finite(xte)] <- med
			dtr <- data.frame(time = tr[[tvar]], event = tr[[evar]], x = xtr)
			ok <- is.finite(dtr$time) & !is.na(dtr$event) & is.finite(dtr$x)
			fit <- switch(method,
				binomial = tryCatch(glm(I(event == 1) ~ x, data = dtr[ok, , drop = FALSE], family = binomial()), error = function(e) NULL),
				gaussian = tryCatch(lm(event ~ x, data = dtr[ok, , drop = FALSE]), error = function(e) NULL),
				cox = tryCatch(coxph(Surv(time, event) ~ x, data = dtr[ok, , drop = FALSE]), error = function(e) NULL)
			)
			if (is.null(fit)) return(tibble::tibble())
			score <- tryCatch(as.numeric(predict(fit, newdata = data.frame(x = xte), type = if (method == "binomial") "response" else if (method == "cox") "lp" else "response")), error = function(e) rep(NA_real_, nrow(te)))
			ci <- if (method == "binomial") tryCatch(as.numeric(pROC::auc(pROC::roc(te[[evar]], score, quiet = TRUE))), error = function(e) NA_real_) else fcidx(Surv(te[[tvar]], te[[evar]]), score)
			data.frame(eid = te$eid, fold = 1L, model = nm, fx = "RE", method = method, time = te[[tvar]], event = te[[evar]], risk = score, cidx = ci, n_sel = 1L)
		}
		for (i in seq_along(configs)) {
			nm <- names(configs)[i] ; cf <- configs[[i]]
			tag <- gsub("[^A-Za-z0-9]+", "_", nm)
			cov_tag <- paste0("c", length(basic), "_", sum(utf8ToInt(paste(sort(basic), collapse = "|"))))
			cp <- if (is.null(checkpoint_dir)) NA_character_ else file.path(checkpoint_dir, sprintf("%02d_%s_n%d_%s_%s.rds", i, tag, length(cf$varX), cov_tag, method))
			if (!force && !is.na(cp) && file.exists(cp) && file.size(cp) > 0) {
				message("Reuse prediction checkpoint: ", basename(cp)) ; zz[[i]] <- readRDS(cp)
			} else {
				zz[[i]] <- if (length(cf$varX) == 1L && !length(cf$vars.basic)) predict_single(nm, cf$varX[[1]]) else
					pred_risk_by_fold(1L, dat, tvar, evar, basic, stats::setNames(list(cf), nm), id.var = "eid")
				if (!is.na(cp) && nrow(zz[[i]])) saveRDS(zz[[i]], cp, compress = FALSE)
			}
			invisible(gc())
		}
		z <- dplyr::bind_rows(zz) ; rm(zz) ; invisible(gc())
		if (!nrow(z)) return(tibble::tibble())
		z |>
			tidyr::separate(model, c("biom_set", "model"), sep = "\\|\\|", fill = "left", extra = "merge") |>
			dplyr::mutate(biom_set = ifelse(model == "Model 0", NA_character_, biom_set)) |>
			dplyr::select(eid, biom_set, model, time, event, score = risk, cindex = cidx, n_selected = n_sel)
	}

	# Apply the final fixed model population-wide; use held-out predictions for metrics.

	# Yu-fair comparison: publication-aligned LightGBM using the same held-out rows
	# as split_re_prediction. Missing biomarker values are left for LightGBM's
	# native missing-value handling, as in the reference implementation.

	plot_evidence_prediction <- function(pred, horizon = 10) {
		if (!nrow(pred)) return(blank_plot("Evidence-based prediction", "No out-of-fold prediction was available"))
		sets <- unique(stats::na.omit(pred$biom_set)) ; clinical <- pred |> dplyr::filter(model == "Model 0")
		pred <- dplyr::bind_rows(pred |> dplyr::filter(model != "Model 0"), lapply(sets, function(s) clinical |> dplyr::mutate(biom_set = s))) |>
			dplyr::mutate(biom_set = factor(biom_set, levels = sets)) |>
			dplyr::group_by(biom_set, model) |>
			dplyr::mutate(z = as.numeric(scale(score)), quintile = dplyr::ntile(score, 5), decile = dplyr::ntile(score, 10)) |>
			dplyr::ungroup()
		pal <- c(`Model 0` = "#333333", Biomarkers = "#A6DDA0", Combined = "#32B43C")
		selected_n <- pred |>
			dplyr::filter(model == "Biomarkers") |>
			dplyr::group_by(biom_set) |>
			dplyr::summarise(n_biom = as.integer(stats::median(n_selected, na.rm = TRUE)), .groups = "drop")
		biom_legend <- paste0("Biom (N=", paste(selected_n$n_biom[match(sets, as.character(selected_n$biom_set))], collapse = "/"), ")")
		pd <- pred |>
			dplyr::filter(model == "Biomarkers") |>
			dplyr::mutate(status = factor(event, 0 : 1, c("Controls", "Cases"))) |>
			dplyr::group_by(biom_set) |>
			dplyr::mutate(threshold = quantile(z[event == 0], .95, na.rm = TRUE), DR = mean(z[event == 1] >= dplyr::first(threshold), na.rm = TRUE)) |>
			dplyr::ungroup()
		ann <- pd |> dplyr::distinct(biom_set, threshold, DR)
		a <- ggplot2::ggplot(pd, ggplot2::aes(z, fill = status)) +
			ggplot2::geom_density(alpha = .48) +
			ggplot2::geom_vline(data = ann, ggplot2::aes(xintercept = threshold)) +
			ggplot2::geom_text(data = ann, ggplot2::aes(x = Inf, y = Inf, label = sprintf("DR = %.1f%%\nFPR = 5.0%%", 100 * DR)), hjust = 1.05, vjust = 1.2, inherit.aes = FALSE, size = 2.7) +
			ggplot2::scale_fill_manual(values = c(Controls = "#F8E8D2", Cases = "#C89D6B")) +
			ggplot2::labs(x = "Biomarker score (s.d.)", y = "Density", fill = NULL) +
			theme_5c(9) +
			ggplot2::theme(legend.position = "top", legend.justification = "center", legend.box.just = "center")
		cum <- pred |>
			dplyr::group_by(biom_set, model, quintile) |>
			dplyr::group_modify( ~ {
				sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = .x) ; tibble::tibble(time = sf$time, cuminc = 1 - sf$surv)
			}) |>
			dplyr::ungroup()
		b <- ggplot2::ggplot(cum |> dplyr::filter(model == "Biomarkers"), ggplot2::aes(time, cuminc, color = factor(quintile))) +
			ggplot2::geom_step() +
			ggplot2::coord_cartesian(xlim = c(0, horizon)) +
			ggplot2::scale_color_brewer(palette = "YlOrBr", direction = 1) +
			ggplot2::labs(x = "Follow-up time (years)", y = "Cumulative incidence", color = "Score quintile") +
			theme_5c(9) +
			ggplot2::theme(legend.position = "top", legend.justification = "center", legend.box.just = "center")
		rate <- pred |>
			dplyr::filter(model == "Biomarkers") |>
			dplyr::group_by(biom_set, decile) |>
			dplyr::summarise(rate = 1000 * sum(event) / sum(time), .groups = "drop") |>
			dplyr::filter(rate > 0)
		c <- ggplot2::ggplot(rate, ggplot2::aes(decile * 10 - 5, rate, color = decile)) +
			ggplot2::geom_point(size = 2) +
			ggplot2::scale_color_gradient(low = "#FFD29B", high = "#A65E00", guide = "none") +
			ggplot2::scale_y_log10() +
			ggplot2::labs(x = "Biomarker score percentile", y = "Incidence rate per 1,000 years") +
			theme_5c(9)
		roc <- pred |>
			dplyr::group_by(biom_set, model) |>
			dplyr::group_modify( ~ {
				r <- pROC::roc(.x$event, .x$score, quiet = TRUE) ; tibble::tibble(fpr = 1 - r$specificities, tpr = r$sensitivities, auc = as.numeric(r$auc))
			}) |>
			dplyr::ungroup()
		labs <- roc |>
			dplyr::group_by(biom_set, model) |>
			dplyr::summarise(auc = dplyr::first(auc), .groups = "drop") |>
			dplyr::left_join(selected_n, by = "biom_set") |>
			dplyr::mutate(model_label = dplyr::case_when(model == "Model 0" ~ "Clinical", model == "Biomarkers" ~ paste0("Biom (N=", n_biom, ")"), TRUE ~ "Combined"), label = sprintf("%s: %.3f", model_label, auc), x = .98, y = c(.12, .06, .18)[match(model, c("Model 0", "Biomarkers", "Combined"))])
		d <- ggplot2::ggplot(roc, ggplot2::aes(fpr, tpr, color = model)) +
			ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey80") +
			ggplot2::geom_step() +
			ggplot2::geom_text(data = labs, ggplot2::aes(x = x, y = y, label = label, color = model), hjust = 1, inherit.aes = FALSE, size = 2.5) +
			ggplot2::scale_color_manual(values = pal, breaks = c("Biomarkers", "Combined", "Model 0"), labels = c(biom_legend, "Combined", "Clinical")) +
			ggplot2::labs(x = "FPR", y = "True positive rate", color = NULL) +
			theme_5c(9) +
			ggplot2::theme(legend.position = "top", legend.justification = "center", legend.box.just = "center")
		(a / b / c / d) + patchwork::plot_layout(guides = "keep") + patchwork::plot_annotation(title = "Prediction by biomarker evidence set: 80% training / 20% testing", theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5))) & ggplot2::facet_wrap( ~ biom_set, nrow = 1, scales = "free")
	}
})
LE8_JOB <- "c1_prediction"
FINAL_CODE_VERSION <- "2026-10-15.c1-prediction"
FINAL_MODEL_VERSION <- "2026-10-15.c1-prediction"

OUTER <- as.integer(Sys.getenv("FINAL_OUTER_FOLDS", unset = "5"))
INNER <- as.integer(Sys.getenv("FINAL_INNER_FOLDS", unset = "10"))
HORIZON <- as.numeric(Sys.getenv("FINAL_HORIZON", unset = "10"))
FINAL_NESTED_CV <- truthy(Sys.getenv("FINAL_NESTED_CV", unset = "FALSE"))
FINAL_TEST_FRAC <- as.numeric(Sys.getenv("FINAL_TEST_FRAC", unset = "0.20"))
FINAL_LGB_FOLDS <- as.integer(Sys.getenv("FINAL_LGB_FOLDS", unset = "5"))
FINAL_LGB_NROUND <- as.integer(Sys.getenv("FINAL_LGB_NROUND", unset = "500"))
FINAL_YU_PRESELECT_N <- as.integer(Sys.getenv("FINAL_YU_PRESELECT_N", unset = "257"))
FINAL_YU_LIST <- Sys.getenv("FINAL_YU_LIST", unset = "")
FINAL_EXTENDED_N <- as.integer(Sys.getenv("FINAL_EXTENDED_N", unset = "30"))
FINAL_YY_BINS <- as.integer(Sys.getenv("FINAL_YY_BINS", unset = "18"))
FINAL_TEMPORAL_STEP <- as.numeric(Sys.getenv("FINAL_TEMPORAL_STEP", unset = "1"))
FINAL_TIME_CAP <- as.numeric(Sys.getenv("FINAL_TIME_CAP", unset = "16"))
FINAL_BOOT <- as.integer(Sys.getenv("FINAL_AUC_BOOT", unset = "150"))
FINAL_MIN_BIN_N <- as.integer(Sys.getenv("FINAL_MIN_BIN_N", unset = "15"))
FINAL_TOPN_MAX <- as.integer(Sys.getenv("FINAL_TOPN_MAX", unset = "30"))
FINAL_TOPN_BOOT <- as.integer(Sys.getenv("FINAL_TOPN_BOOT", unset = "150"))
FINAL_LEAD_MIN_CASES <- as.integer(Sys.getenv("FINAL_LEAD_MIN_CASES", unset = "100"))
FINAL_MECHANISM_N <- as.integer(Sys.getenv("FINAL_MECHANISM_N", unset = "10"))
FINAL_BOOT_MAX_CASES <- as.integer(Sys.getenv("FINAL_BOOT_MAX_CASES", unset = "1500"))
FINAL_BOOT_MAX_CONTROLS <- as.integer(Sys.getenv("FINAL_BOOT_MAX_CONTROLS", unset = "6000"))
FINAL_PAIR_MAX_CASES <- as.integer(Sys.getenv("FINAL_PAIR_MAX_CASES", unset = "1500"))
FINAL_DISTAL_LANDMARK <- as.numeric(Sys.getenv("FINAL_DISTAL_LANDMARK", unset = "5"))
FINAL_DISTAL_TOP <- as.integer(Sys.getenv("FINAL_DISTAL_TOP", unset = "500"))
FINAL_INC_PROT <- Sys.getenv("FINAL_INC_PROT", unset = "")
FINAL_INC_MET <- Sys.getenv("FINAL_INC_MET", unset = "")
FINAL_LEAD_ANCHORS <- unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
	unset = "PCSK9,LPA,GDF15,NTPROBNP,MMP12"
), ",", fixed = TRUE)[[1]]))

# Publication typography for Final only.
theme_5c <- function(base_size = 12) {
	theme_classic(base_size = base_size) + theme(
		plot.title = element_text(face = "bold", size = base_size * 1.12, hjust = 0),
		plot.subtitle = element_text(face = "bold", size = base_size * .92, color = "grey30"),
		axis.title = element_text(face = "bold", size = base_size * 1.03), axis.text = element_text(face = "bold", size = base_size, color = "black"),
		legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), strip.background = element_blank(), strip.text = element_text(face = "bold"),
		panel.grid.major.y = element_line(color = "grey91", linewidth = .25), panel.grid.minor = element_blank(), plot.margin = margin(7, 10, 7, 10)
	)
}


# 🚩 Fixed outer split and score fitting
make_outer_split <- function(dat, evar = NULL, test_frac = FINAL_TEST_FRAC, seed = SEED) {
	c1_lockpaper_protocol()
	shared <- le8_outer_roles(dat)
	shared
}

fit_glmnet_score <- function(dat, vars, tvar, evar, split, label, inner = INNER, alpha = 1,
			lambda_rule = c("lambda.1se", "lambda.min"), penalty_factor = NULL) {
	lambda_rule <- match.arg(lambda_rule)
	if (!requireNamespace("glmnet", quietly = TRUE)) stop("C1 prediction requires glmnet.", call. = FALSE)
	vars <- intersect(unique(vars), names(dat)) ; trid <- which(split == "training") ; vaid <- which(split == "validation")
	if (!length(vars) || length(trid) < 500 || length(vaid) < 100) return(NULL)
	tr <- dat[trid, , drop = FALSE] ; va <- dat[vaid, , drop = FALSE] ; x <- impute_train_test(tr, va, vars)
	if (ncol(x$tr) < 1) return(NULL)
	ok <- is.finite(tr[[tvar]]) & !is.na(tr[[evar]]) & tr[[tvar]] > 0
	y <- Surv(tr[[tvar]][ok], tr[[evar]][ok]) ; events <- as.integer(tr[[evar]][ok] == 1)
	if (sum(ok) < 300 || sum(events) < 30) return(NULL)
	k <- min(inner, max(3L, min(table(events))))
	foldid <- make_folds(tr, evar, k, SEED + 17)
	pf <- setNames(rep(1, ncol(x$tr)), x$vars)
	if (!is.null(penalty_factor)) {
		supplied <- suppressWarnings(as.numeric(penalty_factor)) ; names(supplied) <- names(penalty_factor)
		if (!is.null(names(penalty_factor))) {
			hit <- intersect(names(penalty_factor), names(pf)) ; pf[hit] <- supplied[match(hit, names(penalty_factor))]
		} else if (length(supplied) == length(pf)) pf[] <- supplied
	}
	pf[!is.finite(pf) | pf <= 0] <- 1
	fit <- tryCatch(glmnet::cv.glmnet(x$tr[ok, , drop = FALSE], y,
		family = "cox", alpha = alpha, foldid = foldid[ok],
		type.measure = "C", standardize = FALSE, keep = TRUE, penalty.factor = unname(pf)
	), error = function(e) NULL)
	if (is.null(fit)) return(NULL)
	lam <- if (lambda_rule == "lambda.min") fit$lambda.min else fit$lambda.1se
	li <- which.min(abs(fit$lambda - lam)) ; oof <- rep(NA_real_, nrow(tr))
	if (!is.null(fit$fit.preval) && length(li)) oof[ok] <- as.numeric(fit$fit.preval[, li])
	val <- tryCatch(as.numeric(predict(fit, x$te, s = lambda_rule, type = "link")), error = function(e) rep(NA_real_, nrow(va)))
	b <- as.matrix(coef(fit, s = lambda_rule)) ; sel <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
	rows <- bind_rows(
		tibble(eid = tr$eid, row_id = trid, split = "training", time = tr[[tvar]], event = tr[[evar]], score = oof),
		tibble(eid = va$eid, row_id = vaid, split = "validation", time = va[[tvar]], event = va[[evar]], score = val)
	) |>
		mutate(method = label, n_features = length(x$vars), n_selected = length(sel), engine = "glmnet_cox", lambda_rule = lambda_rule)
	list(
		rows = rows, selected = sel, available = x$vars, fit = fit, engine = "glmnet_cox", lambda_rule = lambda_rule,
		penalty_factor = pf, preprocess = list(vars = x$vars, center = x$center, scale = x$scale)
	)
}

fit_lightgbm_score <- function(dat, vars, tvar, evar, split, label = "Yu-style / LightGBM", k = FINAL_LGB_FOLDS, nrounds = FINAL_LGB_NROUND) {
	if (!requireNamespace("lightgbm", quietly = TRUE)) {
		warning("R package 'lightgbm' is not installed; the Yu-style LightGBM score will be skipped.", call. = FALSE) ; return(NULL)
	}
	vars <- intersect(unique(vars), names(dat)) ; trid <- which(split == "training") ; vaid <- which(split == "validation")
	if (!length(vars) || length(trid) < 500 || length(vaid) < 100) return(NULL)
	tr <- dat[trid, , drop = FALSE] ; va <- dat[vaid, , drop = FALSE]
	# Yu-style LightGBM keeps missing values natively. Only remove invariant columns.
	keep <- vapply(tr[, vars, drop = FALSE], function(z) {
		z <- suppressWarnings(as.numeric(z)) ; sum(is.finite(z)) >= 100 && is.finite(sd(z, na.rm = TRUE)) && sd(z, na.rm = TRUE) > 0
	}, logical(1))
	vars <- vars[keep] ; if (length(vars) < 2) return(NULL)
	to_mat <- function(d) {
		x <- as.matrix(data.frame(lapply(d[, vars, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE)) ; storage.mode(x) <- "double" ; x
	}
	xtr <- to_mat(tr) ; xva <- to_mat(va)
	eligible <- is.finite(tr[[tvar]]) & !is.na(tr[[evar]]) & tr[[tvar]] > 0 &
		((tr[[evar]] == 1 & tr[[tvar]] <= HORIZON) | tr[[tvar]] > HORIZON)
	y <- as.integer(tr[[evar]][eligible] == 1 & tr[[tvar]][eligible] <= HORIZON) ; if (sum(eligible) < 300 || length(unique(y)) < 2) return(NULL)
	nfold <- min(k, max(3L, min(table(y)))) ; set.seed(SEED + 23) ; foldid <- integer(length(y))
	for (value in 0 : 1) {
		ii <- which(y == value) ; foldid[ii] <- sample(rep(seq_len(nfold), length.out = length(ii)))
	}
	oof <- rep(NA_real_, nrow(tr))
	params <- list(
		objective = "binary", metric = "auc", max_depth = 15L, num_leaves = 10L,
		bagging_fraction = .70, bagging_freq = 1L, learning_rate = .01, feature_fraction = .70, max_bin = 63L,
		lambda_l2 = 1, verbosity =  - 1L, num_threads = max(1L, N_CORES), seed = SEED
	)
	for (fd in sort(unique(foldid))) {
		it <- which(foldid != fd) ; iv <- which(foldid == fd) ; if (length(unique(y[it])) < 2) next
		ds <- lightgbm::lgb.Dataset(data = xtr[eligible, , drop = FALSE][it, , drop = FALSE], label = y[it])
		ff <- tryCatch(lightgbm::lgb.train(params = params, data = ds, nrounds = nrounds, verbose =  - 1), error = function(e) NULL)
		if (!is.null(ff)) oof[which(eligible)[iv]] <- as.numeric(predict(ff, xtr[eligible, , drop = FALSE][iv, , drop = FALSE]))
	}
	ds <- lightgbm::lgb.Dataset(data = xtr[eligible, , drop = FALSE], label = y) ; fit <- tryCatch(lightgbm::lgb.train(params = params, data = ds, nrounds = nrounds, verbose =  - 1), error = function(e) NULL)
	if (is.null(fit)) return(NULL) ; val <- as.numeric(predict(fit, xva))
	imp <- tryCatch(lightgbm::lgb.importance(fit), error = function(e) NULL) ; sel <- if (is.null(imp)) vars else as.character(imp$Feature[imp$Gain > 0])
	rows <- bind_rows(
		tibble(eid = tr$eid, row_id = trid, split = "training", time = tr[[tvar]], event = tr[[evar]], score = oof),
		tibble(eid = va$eid, row_id = vaid, split = "validation", time = va[[tvar]], event = va[[evar]], score = val)
	) |>
		mutate(method = label, n_features = length(vars), n_selected = length(sel), engine = "lightgbm", lambda_rule = NA_character_)
	list(rows = rows, selected = sel, available = vars, fit = fit, importance = imp, engine = "lightgbm", horizon = HORIZON, preprocess = list(vars = vars))
}

empty_score_rows <- function() {
	tibble(
		eid = numeric(), row_id = integer(), split = character(), time = numeric(), event = numeric(),
		score = numeric(), method = character(), n_features = integer(), n_selected = integer(), engine = character(),
		lambda_rule = character(), score_z = numeric(), score_native = numeric(), score_scale = character()
	)
}

empty_prediction_rows <- function() {
	tibble(
		eid = numeric(), row_id = integer(), time = numeric(), event = numeric(), score = numeric(),
		biom_set = character(), model = character(), n_selected = integer(), cindex = numeric()
	)
}

standardize_scores <- function(rows) {
	if (is.null(rows) || !nrow(rows) || !all(c("method", "score", "split") %in% names(rows))) return(empty_score_rows())
	rows |>
		group_by(method) |>
		group_modify(function(d, key) {
			m <- mean(d$score[d$split == "training"], na.rm = TRUE) ; s <- sd(d$score[d$split == "training"], na.rm = TRUE) ; if (!is.finite(s) || s == 0) s <- 1
			# A Cox linear predictor is not a probability; plogis(lp) was therefore a
			# misleading "native 0-1" scale. Preserve raw model output and use the
			# training-standardized score for cross-model temporal displays.
			d |> mutate(
				score_z = (score - m) / s, score_native = score,
				score_scale = ifelse(engine == "lightgbm", "binary model output", "Cox linear predictor")
			)
		}) |>
		ungroup()
}

predict_prevalent_rows <- function(obj, prevalent, method, bvar, training_rows) {
	if (is.null(obj) || !nrow(prevalent) || is.null(obj$fit)) return(tibble())
	expected <- obj$preprocess$vars %||% obj$available
	missing <- setdiff(expected, names(prevalent)) ; if (length(missing)) {
		warning(method, ": prevalent scoring missing ", length(missing), " model variables; no partial-column prediction was attempted.", call. = FALSE)
		return(tibble())
	}
	vars <- expected ; if (!length(vars)) return(tibble())
	if (identical(obj$engine, "lightgbm")) {
		x <- as.matrix(data.frame(lapply(prevalent[, vars, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE)) ; storage.mode(x) <- "double"
		score <- tryCatch(as.numeric(predict(obj$fit, x)), error = function(e) rep(NA_real_, nrow(prevalent)))
	} else {
		cen <- obj$preprocess$center ; sca <- obj$preprocess$scale
		if (is.null(cen) || is.null(sca) || !all(vars %in% names(cen)) || !all(vars %in% names(sca))) return(tibble())
		cen <- cen[vars] ; sca <- sca[vars] ; x <- as.data.frame(prevalent[, vars, drop = FALSE]) ; x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
		for (j in seq_along(vars)) x[!is.finite(x[, j]), j] <- cen[[j]]
		xm <- scale(as.matrix(x), center = cen, scale = sca)
		score <- tryCatch(as.numeric(predict(obj$fit, newx = xm, s = obj$lambda_rule, type = "link")), error = function(e) rep(NA_real_, nrow(prevalent)))
	}
	tr <- training_rows |> filter(method == .env$method, split == "training", is.finite(score)) ; m <- mean(tr$score, na.rm = TRUE) ; s <- sd(tr$score, na.rm = TRUE) ; if (!is.finite(s) || s == 0) s <- 1
	tibble(
		eid = prevalent$eid, row_id = NA_integer_, split = "prevalent", time = prevalent[[bvar]], event = 1, score = score,
		method = method, n_features = length(vars), n_selected = length(obj$selected %||% character()), engine = obj$engine,
		lambda_rule = obj$lambda_rule %||% NA_character_, score_z = (score - m) / s,
		score_native = score, score_scale = ifelse(obj$engine == "lightgbm", "binary model output", "Cox linear predictor")
	)
}

score_coverage_audit <- function(score_rows, expected_methods, prevalent_n) {
	expected <- tidyr::crossing(method = expected_methods, split = c("training", "validation", "prevalent")) |>
		mutate(expected_n = ifelse(split == "prevalent", prevalent_n, NA_integer_))
	got <- score_rows |>
		group_by(method, split) |>
		summarise(rows = n(), finite_score = sum(is.finite(score)), finite_z = sum(is.finite(score_z)), .groups = "drop")
	expected |>
		left_join(got, by = c("method", "split")) |>
		mutate(across(c(rows, finite_score, finite_z), ~ replace_na(.x, 0L)),
			coverage = ifelse(rows > 0, finite_z / rows, 0), status = case_when(rows == 0 ~ "missing rows", finite_z == 0 ~ "all predictions non-finite", coverage < .95 ~ "partial coverage", TRUE ~ "ok")
		)
}

fit_combined_score <- function(dat, score_obj, clinical, tvar, evar, split, label) {
	if (is.null(score_obj)) return(NULL)
	vars <- unique(c(score_obj$available, clinical))
	if (identical(score_obj$engine, "lightgbm")) fit_lightgbm_score(dat, vars, tvar, evar, split, paste0(label, " + clinical"))
	else {
		pf <- score_obj$penalty_factor %||% setNames(rep(1, length(score_obj$available)), score_obj$available)
		pf <- c(pf, setNames(rep(1, length(setdiff(clinical, names(pf)))), setdiff(clinical, names(pf))))
		fit_glmnet_score(dat, vars, tvar, evar, split, paste0(label, " + clinical"),
			lambda_rule = "lambda.1se", penalty_factor = pf
		)
	}
}

make_validation_bundle <- function(dat, score_obj, clinical_obj, clinical, tvar, evar, split, label) {
	if (is.null(score_obj)) return(tibble())
	bio <- score_obj$rows |>
		filter(split == "validation") |>
		transmute(eid, row_id, time, event, score, biom_set = label, model = "Biomarkers", n_selected = first(score_obj$rows$n_selected))
	comb <- fit_combined_score(dat, score_obj, clinical, tvar, evar, split, label)
	cc <- if (is.null(comb)) tibble() else comb$rows |>
		filter(split == "validation") |>
		transmute(eid, row_id, time, event, score, biom_set = label, model = "Combined", n_selected = first(score_obj$rows$n_selected))
	cl <- if (is.null(clinical_obj)) tibble() else clinical_obj$rows |>
		filter(split == "validation") |>
		transmute(eid, row_id, time, event, score, biom_set = label, model = "Model 0", n_selected = first(score_obj$rows$n_selected))
	bind_rows(bio, cc, cl) |>
		group_by(biom_set, model) |>
		mutate(cindex = fcidx(Surv(time, event), score)) |>
		ungroup()
}


# 🚩 Common temporal helpers
yy_score_summary <- function(rows, method, bins = FINAL_YY_BINS, min_bin_n = FINAL_MIN_BIN_N) {
	d <- rows |> filter(method == .env$method, event == 1, is.finite(time), is.finite(score_z)) ; if (!nrow(d)) return(list(lines = tibble(), hist = tibble()))
	lim <- max(d$time, na.rm = TRUE) ; br <- seq(0, lim, length.out = bins + 1) ; mids <- (head(br, - 1) + tail(br, - 1)) / 2
	lines <- d |>
		mutate(bin = cut(time, br, include.lowest = TRUE, labels = FALSE)) |>
		filter(!is.na(bin)) |>
		group_by(split, bin) |>
		summarise(mean = mean(score_z), sd = sd(score_z), se = sd / sqrt(n()), N = n(), N_total = n_distinct(eid), .groups = "drop") |>
		filter(N >= min_bin_n) |>
		mutate(year = mids[bin])
	hh <- hist(d$time, breaks = br, plot = FALSE) ; histdf <- tibble(xmin = head(br, - 1), xmax = tail(br, - 1), mid = mids, n = hh$counts)
	list(lines = lines, hist = histdf)
}

risk_set_pairs <- function(rows, dat, method, score_col = c("score_z", "score_native"), seed = SEED, max_cases = FINAL_PAIR_MAX_CASES) {
	score_col <- match.arg(score_col)
	keepcov <- intersect(c("eid", "age", "sex", "ethnic.c"), names(dat))
	d <- rows |>
		filter(method == .env$method, is.finite(.data[[score_col]]), is.finite(time)) |>
		left_join(dat |> select(all_of(keepcov)), by = "eid")
	out <- list() ; ii <- 0L ; set.seed(seed)
	for (sp in intersect(c("training", "validation"), unique(d$split))) {
		ca <- d |>
			filter(split == sp, event == 1) |>
			arrange(time) ; co <- d |> filter(split == sp, event == 0)
		if (!nrow(ca) || !nrow(co)) next
		if (nrow(ca) > max_cases) ca <- ca |>
			slice_sample(n = max_cases) |>
			arrange(time)
		for (i in seq_len(nrow(ca))) {
			pool <- which(is.finite(co$time) & co$time >= ca$time[i]) ; if (!length(pool)) next
			if ("sex" %in% names(d) && !is.na(ca$sex[i])) {
				same <- pool[co$sex[pool] == ca$sex[i]] ; if (length(same)) pool <- same
			}
			if ("ethnic.c" %in% names(d) && !is.na(ca$ethnic.c[i])) {
				same <- pool[as.character(co$ethnic.c[pool]) == as.character(ca$ethnic.c[i])] ; if (length(same)) pool <- same
			}
			if ("age" %in% names(d) && is.finite(ca$age[i])) {
				ad <- abs(co$age[pool] - ca$age[i]) ; pool <- pool[order(ad, na.last = TRUE)][seq_len(min(20, length(pool)))]
			}
			j <- sample(pool, 1) ; ii <- ii + 1L
			out[[ii]] <- tibble(
				pair_id = paste(sp, i, sep = "_"), split = sp, index_time = ca$time[i], case_eid = ca$eid[i], control_eid = co$eid[j],
				case_score = ca[[score_col]][i], control_score = co[[score_col]][j]
			)
		}
	}
	bind_rows(out)
}

pair_long <- function(pairs) {
	if (!nrow(pairs)) return(tibble())
	bind_rows(
		pairs |> transmute(pair_id, split, index_time, year =  - index_time, eid = case_eid, group = ifelse(split == "training", "Cases: training", "Cases: validation"), score = case_score),
		pairs |> transmute(pair_id, split, index_time, year =  - index_time, eid = control_eid, group = "Controls", score = control_score)
	)
}

trapz_mean <- function(x, y) {
	ok <- is.finite(x) & is.finite(y) ; x <- x[ok] ; y <- y[ok] ; if (length(x) < 2) return(NA_real_) ; o <- order(x) ; x <- x[o] ; y <- y[o] ; span <- max(x) - min(x) ; if (span <= 0) return(NA_real_) ; sum(diff(x) * (head(y, - 1) + tail(y, - 1)) / 2) / span
}
temporal_implied_auc <- function(pairs, step = FINAL_TEMPORAL_STEP) {
	if (!nrow(pairs)) return(c(dbar = NA_real_, implied_auc = NA_real_))
	z <- pairs |>
		mutate(bin = floor(index_time / step) * step + step / 2, diff = case_score - control_score) |>
		group_by(bin) |>
		summarise(d = mean(diff, na.rm = TRUE), N = n(), .groups = "drop") |>
		filter(N >= 5)
	dbar <- trapz_mean(z$bin, z$d) ; c(dbar = dbar, implied_auc = ifelse(is.finite(dbar), pnorm(dbar / sqrt(2)), NA_real_))
}
match_bootstrap_pairs <- function(v, seed, max_cases = FINAL_PAIR_MAX_CASES) {
	set.seed(seed) ; ca <- v |> filter(event == 1) ; co <- v |> filter(event == 0) ; if (!nrow(ca) || !nrow(co)) return(tibble())
	if (nrow(ca) > max_cases) ca <- ca |> slice_sample(n = max_cases)
	out <- list() ; for (i in seq_len(nrow(ca))) {
		pool <- which(is.finite(co$time) & co$time >= ca$time[i]) ; if (!length(pool)) next
		if ("sex" %in% names(v) && !is.na(ca$sex[i])) {
			same <- pool[co$sex[pool] == ca$sex[i]] ; if (length(same)) pool <- same
		}
		if ("ethnic.c" %in% names(v) && !is.na(ca$ethnic.c[i])) {
			same <- pool[as.character(co$ethnic.c[pool]) == as.character(ca$ethnic.c[i])] ; if (length(same)) pool <- same
		}
		if ("age" %in% names(v) && is.finite(ca$age[i])) {
			ad <- abs(co$age[pool] - ca$age[i]) ; pool <- pool[order(ad, na.last = TRUE)][seq_len(min(20, length(pool)))]
		}
		j <- sample(pool, 1) ; out[[length(out) + 1]] <- tibble(index_time = ca$time[i], case_score = ca$score_z[i], control_score = co$score_z[j])
	}
	bind_rows(out)
}
bootstrap_auc_link <- function(method_rows, dat, B = FINAL_BOOT) {
	keepcov <- intersect(c("eid", "age", "sex", "ethnic.c"), names(dat))
	v <- method_rows |>
		filter(split == "validation") |>
		left_join(dat |> select(all_of(keepcov)), by = "eid") |>
		filter(is.finite(score_z), !is.na(event), is.finite(time))
	if (nrow(v) < 200) return(tibble())
	ca <- v |> filter(event == 1) ; co <- v |> filter(event == 0) ; if (!nrow(ca) || !nrow(co)) return(tibble())
	base_pairs <- match_bootstrap_pairs(v, SEED + 1999, FINAL_PAIR_MAX_CASES)
	map_dfr(seq_len(B), function(b) {
		set.seed(SEED + 1000 + b)
		vb <- bind_rows(
			ca[sample(seq_len(nrow(ca)), min(nrow(ca), FINAL_BOOT_MAX_CASES), replace = TRUE), , drop = FALSE],
			co[sample(seq_len(nrow(co)), min(nrow(co), FINAL_BOOT_MAX_CONTROLS), replace = TRUE), , drop = FALSE]
		)
		roc_auc <- weighted_time_auc(vb$time, vb$event, vb$score_z, HORIZON)
		pairs <- if (nrow(base_pairs)) base_pairs[sample(seq_len(nrow(base_pairs)), nrow(base_pairs), replace = TRUE), , drop = FALSE] else tibble() ; ta <- temporal_implied_auc(pairs)
		tibble(bootstrap = b, roc_auc = roc_auc, temporal_d = ta[["dbar"]], trajectory_implied_auc = ta[["implied_auc"]])
	})
}

# Held-out mechanism benchmark.  Feature ranking, signs, Cox weights, means and
# scales are learned only in the outer training split.  The comparators are a
# locked training-selected top-1 biomarker, an unweighted top-N mean and a
# training marginal-Cox-weighted top-N score.  We also report the distribution of the N
# individual biomarkers; no "best" biomarker is selected on validation data.
topn_score_benchmark <- function(dat, ranked, screen, tvar, evar, split, clinical = character(),
			max_n = FINAL_TOPN_MAX, B = FINAL_TOPN_BOOT) {
	empty <- list(rows = tibble(), summary = tibble(), bootstrap = tibble(), individuals = tibble(), individual_range = tibble())
	ranked <- head(intersect(ranked, names(dat)), max(1L, max_n)) ; if (!length(ranked)) return(empty)
	it <- which(split == "training") ; iv <- which(split == "validation") ; if (length(it) < 300 || length(iv) < 100) return(empty)
	center <- vapply(dat[it, ranked, drop = FALSE], function(x) mean(suppressWarnings(as.numeric(x)), na.rm = TRUE), numeric(1))
	scale0 <- vapply(dat[it, ranked, drop = FALSE], function(x) sd(suppressWarnings(as.numeric(x)), na.rm = TRUE), numeric(1))
	keep <- is.finite(center) & is.finite(scale0) & scale0 > 0 ; ranked <- ranked[keep] ; center <- center[keep] ; scale0 <- scale0[keep]
	if (!length(ranked)) return(empty)
	make_x <- function(ii) {
		x <- as.matrix(data.frame(lapply(dat[ii, ranked, drop = FALSE], function(z) suppressWarnings(as.numeric(z))), check.names = FALSE)) ; storage.mode(x) <- "double"
		for (j in seq_along(ranked)) x[!is.finite(x[, j]), j] <- center[j]
		sweep(sweep(x, 2, center, "-"), 2, scale0, "/")
	}
	xtr <- make_x(it) ; xva <- make_x(iv) ; w <- screen$beta[match(ranked, screen$term)] ; w[!is.finite(w) | w == 0] <- 1
	clinical <- intersect(clinical, names(dat)) ; base <- tibble(eid = dat$eid[iv], time = dat[[tvar]][iv], event = dat[[evar]][iv], row_id = iv)
	standardize_pair <- function(tr, va) {
		m <- mean(tr, na.rm = TRUE) ; s <- sd(tr, na.rm = TRUE) ; if (!is.finite(s) || s == 0) s <- 1 ; list(tr = (tr - m) / s, va = (va - m) / s)
	}
	metric_vec <- function(score) {
		d <- base |> mutate(score = score) ; auc <- weighted_time_auc(d$time, d$event, d$score, HORIZON)
		hz <- d |>
			filter(is.finite(score), is.finite(time), !is.na(event)) |>
			mutate(class = case_when(event == 1 & time <= HORIZON ~ "Case", time > HORIZON ~ "Control", TRUE ~ NA_character_)) |>
			filter(!is.na(class))
		ca <- hz$score[hz$class == "Case"] ; co <- hz$score[hz$class == "Control"]
		delta <- if (length(ca) && length(co)) mean(ca) - mean(co) else NA_real_
		pooled <- if (length(ca) > 1 && length(co) > 1) sqrt(((length(ca) - 1) * var(ca) + (length(co) - 1) * var(co)) / (length(ca) + length(co) - 2)) else NA_real_
		dd <- dat[iv, unique(c(tvar, evar, clinical)), drop = FALSE] ; dd$.score <- score ; dd <- dd[complete.cases(dd), , drop = FALSE] ; dd <- dd[dd[[tvar]] > 0, , drop = FALSE]
		rhs <- c(".score", clinical) ; fit <- if (nrow(dd) >= 100 && sum(dd[[evar]] == 1) >= 20) tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(rhs), collapse = " + "))), dd, ties = "efron"), error = function(e) NULL) else NULL
		sm <- if (is.null(fit)) NULL else coef(summary(fit)) ; beta <- if (!is.null(sm) && ".score" %in% rownames(sm)) sm[".score", "coef"] else NA_real_ ; se <- if (!is.null(sm) && ".score" %in% rownames(sm)) sm[".score", "se(coef)"] else NA_real_
		inc <- d |> filter(event == 1, is.finite(time), is.finite(score)) ; slope <- if (nrow(inc) >= 20) tryCatch(unname(coef(lm(score ~ time, data = inc))[2]), error = function(e) NA_real_) else NA_real_
		tibble(
			AUC = auc, beta = beta, beta_se = se, beta_lo = beta - 1.96 * se, beta_hi = beta + 1.96 * se,
			mean_separation = delta, within_group_SD = pooled, cohen_d = delta / pooled, score_SD = sd(score, na.rm = TRUE),
			prediagnostic_steepness =  - slope, cases = sum(d$event == 1, na.rm = TRUE), N_validation = nrow(d)
		)
	}
	boot_vec <- function(score, n, model_index) {
		if (B <= 0) return(tibble()) ; ca <- which(base$event == 1 & is.finite(score)) ; co <- which(base$event == 0 & is.finite(score))
		if (length(ca) < 20 || length(co) < 30) return(tibble())
		map_dfr(seq_len(B), function(b) {
			set.seed(SEED + 50000 + 1000 * n + 100 * model_index + b)
			ii <- c(sample(ca, min(length(ca), FINAL_BOOT_MAX_CASES), replace = TRUE), sample(co, min(length(co), FINAL_BOOT_MAX_CONTROLS), replace = TRUE))
			d <- base[ii, , drop = FALSE] |> mutate(score = score[ii]) ; auc <- weighted_time_auc(d$time, d$event, d$score, HORIZON)
			hz <- d |>
				mutate(class = case_when(event == 1 & time <= HORIZON ~ "Case", time > HORIZON ~ "Control", TRUE ~ NA_character_)) |>
				filter(!is.na(class))
			ca0 <- hz$score[hz$class == "Case"] ; co0 <- hz$score[hz$class == "Control"]
			delta <- if (length(ca0) && length(co0)) mean(ca0) - mean(co0) else NA_real_ ; pooled <- if (length(ca0) > 1 && length(co0) > 1) sqrt(((length(ca0) - 1) * var(ca0) + (length(co0) - 1) * var(co0)) / (length(ca0) + length(co0) - 2)) else NA_real_
			inc <- d |> filter(event == 1, is.finite(time), is.finite(score)) ; slope <- if (nrow(inc) >= 20) tryCatch(unname(coef(lm(score ~ time, data = inc))[2]), error = function(e) NA_real_) else NA_real_
			tibble(bootstrap = b, AUC = auc, mean_separation = delta, within_group_SD = pooled, cohen_d = delta / pooled, prediagnostic_steepness =  - slope)
		})
	}
	extra_n <- if (length(ranked) >= 15L) seq(15L, length(ranked), by = 5L) else integer() ; ns <- unique(pmin(length(ranked), c(seq_len(min(10L, length(ranked))), extra_n, length(ranked))))
	individual <- map_dfr(seq_along(ranked), function(j) metric_vec(xva[, j] * sign(w[j])) |> mutate(rank = j, feature = ranked[j]))
	summaries <- list() ; boots <- list() ; trajectory <- list() ; kk <- 0L ; target_n <- ns[which.min(abs(ns - min(FINAL_MECHANISM_N, length(ranked))))]
	for (n in ns) {
		ww <- w[seq_len(n)] ; sgn <- sign(ww) ; sgn[sgn == 0] <- 1
		weighted <- standardize_pair(as.numeric(xtr[, seq_len(n), drop = FALSE] %*% ww), as.numeric(xva[, seq_len(n), drop = FALSE] %*% ww))
		unweighted <- standardize_pair(rowMeans(sweep(xtr[, seq_len(n), drop = FALSE], 2, sgn, "*")), rowMeans(sweep(xva[, seq_len(n), drop = FALSE], 2, sgn, "*")))
		top1 <- list(tr = xtr[, 1] * sign(w[1]), va = xva[, 1] * sign(w[1]))
		models <- list(
			`Top-1 biomarker (locked)` = top1, `Unweighted top-N mean` = unweighted,
			`Marginal-Cox-weighted top-N` = weighted
		)
		cm <- suppressWarnings(cor(xtr[, seq_len(n), drop = FALSE], use = "pairwise.complete.obs")) ; avgcor <- if (n > 1) mean(abs(cm[row(cm) != col(cm)]), na.rm = TRUE) else 0
		for (mi in seq_along(models)) {
			nm <- names(models)[mi] ; sc <- models[[mi]]$va ; weights0 <- if (mi == 1) c(1, rep(0, n - 1)) else if (mi == 2) sgn else ww
			effn <- if (sum(abs(weights0)) > 0) 1 / sum((abs(weights0) / sum(abs(weights0))) ^ 2) else NA_real_
			kk <- kk + 1L ; summaries[[kk]] <- metric_vec(sc) |> mutate(N = n, model = nm, feature = if (nm == "Top-1 biomarker (locked)") ranked[1] else paste0("top ", n), average_abs_correlation = avgcor, effective_N = effn)
			boots[[kk]] <- boot_vec(sc, n, mi) |> mutate(N = n, model = nm)
			if (n == target_n) trajectory[[kk]] <- base |> transmute(eid, time, event, N = n, model = nm, score = sc)
		}
	}
	summary <- bind_rows(summaries) ; boot <- bind_rows(boots)
	if (nrow(boot)) summary <- summary |> left_join(boot |> group_by(N, model) |> summarise(
		AUC_boot_SD = sd(AUC, na.rm = TRUE), separation_boot_SD = sd(mean_separation, na.rm = TRUE),
		slope_boot_SD = sd(prediagnostic_steepness, na.rm = TRUE), AUC_lo = quantile(AUC, .025, na.rm = TRUE), AUC_hi = quantile(AUC, .975, na.rm = TRUE), .groups = "drop"
	), by = c("N", "model"))
	for (nm in c("AUC_boot_SD", "separation_boot_SD", "slope_boot_SD", "AUC_lo", "AUC_hi")) if (!nm %in% names(summary)) summary[[nm]] <- NA_real_
	irange <- map_dfr(ns, function(n) individual |>
		filter(rank <= n) |>
		summarise(
			N = n, AUC_q25 = quantile(AUC, .25, na.rm = TRUE), AUC_median = median(AUC, na.rm = TRUE), AUC_q75 = quantile(AUC, .75, na.rm = TRUE), AUC_max = max(AUC, na.rm = TRUE),
			beta_q25 = quantile(beta, .25, na.rm = TRUE), beta_median = median(beta, na.rm = TRUE), beta_q75 = quantile(beta, .75, na.rm = TRUE), cohen_d_median = median(cohen_d, na.rm = TRUE)
		))
	list(rows = bind_rows(trajectory), summary = summary, bootstrap = boot, individuals = individual, individual_range = irange)
}


# 🚩 The five reusable row functions.  All Final score figures call these directly.
final_row1_score_profile <- function(score_rows, label, n_selected = NA_integer_) {
	d <- score_rows |> filter(method == .env$label, split == "validation", is.finite(score_z), is.finite(time), time > 0, !is.na(event))
	ttl <- paste0(label, " (N=", ifelse(is.finite(n_selected), format(as.integer(n_selected), big.mark = ","), "NA"), ")")
	if (nrow(d) < 100) return(blank_plot(ttl, "Insufficient validation scores"))
	status <- factor(d$event, levels = 0 : 1, labels = c("Controls", "Cases")) ; xr <- range(d$score_z, na.rm = TRUE) ; if (diff(xr) <= 0) xr <- xr + c( - 1, 1)
	den <- map_dfr(levels(status), function(g) {
		x <- d$score_z[status == g] ; if (length(x) < 5) return(tibble()) ; z <- density(x, from = xr[1], to = xr[2], n = 256, na.rm = TRUE) ; tibble(x = z$x, density = z$y, status = g)
	})
	rate <- d |>
		mutate(bin = ntile(score_z, 10)) |>
		group_by(bin) |>
		summarise(x = mean(score_z), rate = 1000 * sum(event) / sum(time), N = n(), .groups = "drop") |>
		filter(is.finite(rate))
	dmax <- max(den$density, na.rm = TRUE) ; if (!is.finite(dmax) || dmax <= 0) dmax <- 1
	rr <- range(rate$rate, na.rm = TRUE) ; if (!all(is.finite(rr))) rr <- c(0, 1) ; if (diff(rr) < 1e-9) rr <- rr + c( - .5, .5)
	a <- .78 * dmax / diff(rr) ; offset <- .08 * dmax - a * rr[1] ; rate <- rate |> mutate(y = a * rate + offset)
	ggplot() +
		geom_area(data = den, aes(x, density, fill = status, group = status), alpha = .38, position = "identity") +
		geom_line(data = den, aes(x, density, color = status, group = status), linewidth = .55) +
		geom_line(data = rate, aes(x, y), color = "black", linewidth = .85) +
		geom_point(data = rate, aes(x, y), shape = 21, fill = "white", size = 2.2, stroke = .75) +
		scale_y_continuous(name = "Density", sec.axis = sec_axis( ~ (.x - offset) / a, name = "Incidence / 1,000 person-years")) +
		scale_fill_manual(values = c(Controls = "grey78", Cases = "#D9A066")) +
		scale_color_manual(values = c(Controls = "grey48", Cases = "#A8612B")) +
		coord_cartesian(xlim = xr, clip = "off") +
		labs(title = ttl, x = "Biomarker score (s.d.)", fill = NULL, color = NULL) +
		theme_5c(9) +
		theme(legend.position = "top", legend.justification = "center")
}

ipcw_roc_curve <- function(time, event, score, horizon = HORIZON) {
	ok <- is.finite(time) & !is.na(event) & is.finite(score) & time > 0
	time <- time[ok] ; event <- event[ok] ; score <- score[ok] ; if (length(score) < 100) return(tibble())
	ip <- ipcw_at_horizon(time, event, horizon) ; cw <- ifelse(ip$case, ip$weight, 0) ; tw <- ifelse(ip$control, ip$weight, 0)
	keep <- cw > 0 | tw > 0 ; if (sum(cw[keep]) <= 0 || sum(tw[keep]) <= 0) return(tibble())
	d <- tibble(score = score[keep], cw = cw[keep], tw = tw[keep]) |>
		group_by(score) |>
		summarise(cw = sum(cw), tw = sum(tw), .groups = "drop") |>
		arrange(desc(score)) |>
		mutate(tpr = cumsum(cw) / sum(cw), fpr = cumsum(tw) / sum(tw))
	bind_rows(tibble(score = Inf, tpr = 0, fpr = 0), d |> select(score, tpr, fpr), tibble(score =  - Inf, tpr = 1, fpr = 1))
}

final_row2_roc <- function(pred, label, boot = tibble()) {
	d <- pred |> filter(biom_set == .env$label) ; if (!nrow(d)) return(blank_plot("Validation ROC", "No validation predictions"))
	roc <- d |>
		group_by(model) |>
		group_modify( ~ {
			z <- ipcw_roc_curve(.x$time, .x$event, .x$score, HORIZON) ; if (!nrow(z)) return(tibble())
			z |> mutate(auc = weighted_time_auc(.x$time, .x$event, .x$score, HORIZON))
		}) |>
		ungroup()
	if (!nrow(roc)) return(blank_plot("Validation ROC", "ROC unavailable"))
	al <- roc |>
		group_by(model) |>
		summarise(auc = first(auc), .groups = "drop") |>
		mutate(txt = sprintf("%s: %.3f", case_when(model == "Model 0" ~ "Clinical", model == "Biomarkers" ~ "Biomarkers", TRUE ~ "Combined"), auc), x = .03, y = c(.97, .90, .83)[match(model, c("Model 0", "Biomarkers", "Combined"))])
	p <- ggplot(roc, aes(fpr, tpr, color = model)) +
		geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey78") +
		geom_step(linewidth = .9) +
		geom_text(data = al, aes(x = x, y = y, label = txt, color = model), hjust = 0, vjust = 1, inherit.aes = FALSE, show.legend = FALSE, size = 2.8, fontface = "bold") +
		scale_color_manual(values = c(`Model 0` = "grey35", Biomarkers = "#C77732", Combined = "#287A58"), breaks = c("Biomarkers", "Combined", "Model 0"), labels = c("Biomarkers", "Combined", "Clinical")) +
		labs(x = "False positive rate", y = "True positive rate", color = NULL, subtitle = paste0(HORIZON, "-year IPCW cumulative/dynamic ROC")) +
		theme_5c(9) +
		theme(legend.position = "top")
	if (nrow(boot) >= 20) {
		rr <- suppressWarnings(cor(boot$trajectory_implied_auc, boot$roc_auc, use = "complete.obs"))
		ins <- ggplot(boot, aes(trajectory_implied_auc, roc_auc)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey65") +
			geom_point(alpha = .28, size = .7) +
			geom_smooth(method = "lm", se = FALSE, linewidth = .45, color = "grey25") +
			annotate("text", x = Inf, y =  - Inf, label = sprintf("r=%.2f", rr), hjust = 1.08, vjust =  - .4, size = 2.4, fontface = "bold") +
			labs(x = "Temporal AUC", y = "ROC AUC") +
			theme_classic(base_size = 6) +
			theme(axis.title = element_text(face = "bold"), axis.text = element_text(size = 5), plot.background = element_rect(fill = "white", color = "grey75"))
		p <- p + patchwork::inset_element(ins, left = .38, bottom = .10, right = .84, top = .55, align_to = "panel", on_top = TRUE)
	}
	p + labs(title = NULL)
}

final_row3_followup <- function(score_rows, label, horizon = HORIZON) {
	allm <- score_rows |> filter(method == .env$label, is.finite(score_z)) ; tr <- allm |> filter(split == "training")
	if (nrow(tr) < 100) return(blank_plot("Diagnosis burden across baseline", "No training score distribution"))
	br <- unique(as.numeric(quantile(tr$score_z, probs = seq(0, 1, .2), na.rm = TRUE))) ; if (length(br) < 6) br <- seq(min(tr$score_z), max(tr$score_z), length.out = 6)
	addq <- function(x) x |> mutate(quintile = factor(cut(score_z, breaks = br, include.lowest = TRUE, labels = 1 : 5), levels = 1 : 5, labels = paste0("Q", 1 : 5)))
	pos <- addq(allm |> filter(split == "validation", is.finite(time), time > 0, !is.na(event)))
	cp <- if (!nrow(pos)) tibble() else pos |>
		filter(!is.na(quintile)) |>
		group_by(quintile) |>
		group_modify( ~ {
			sf <- survfit(Surv(time, event) ~ 1, data = .x) ; tibble(time = sf$time, cuminc = 1 - sf$surv)
		}) |>
		ungroup() |>
		mutate(phase = "Incident follow-up")
	p_inc <- if (!nrow(cp)) blank_plot("Incident follow-up", "No validation survival curve") else ggplot(cp, aes(time, cuminc, color = quintile)) +
		geom_line(linewidth = .78) +
		coord_cartesian(xlim = c(0, FINAL_TIME_CAP)) +
		scale_color_brewer(palette = "YlOrBr", direction = 1) +
		labs(title = "Incident follow-up", x = "Years after baseline", y = "Cumulative incidence", color = "Score quintile") +
		theme_5c(8) +
		theme(legend.position = "top")
	# Do not draw a pseudo-survival curve before baseline. Prevalent cases are a
	# cross-sectional disease-state comparison, so report their proportion in
	# each score quintile using the full scored cohort.
	base <- addq(allm |> filter(split %in% c("training", "validation", "prevalent"))) |> filter(!is.na(quintile))
	prev <- base |>
		group_by(quintile) |>
		summarise(N = n(), prevalent = sum(split == "prevalent"), proportion = prevalent / N, .groups = "drop")
	if (!any(base$split == "prevalent")) prev <- prev[0, , drop = FALSE]
	p_prev <- if (!nrow(prev)) blank_plot("Baseline prevalence", "No prevalent score") else ggplot(prev, aes(quintile, proportion, fill = quintile)) +
		geom_col(width = .72, show.legend = FALSE) +
		geom_text(aes(label = scales::percent(proportion, accuracy = .1)), vjust =  - .25, size = 2.6, fontface = "bold") +
		scale_fill_brewer(palette = "YlOrBr", direction = 1) +
		scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .16))) +
		labs(title = "Baseline prevalence", subtitle = "Cross-sectional; not a reverse-time survival curve", x = "Score quintile", y = "Prevalent proportion") +
		theme_5c(8)
	p_prev | p_inc
}

final_row4_yy <- function(score_rows, label) {
	d <- score_rows |>
		filter(method == .env$label, event == 1, is.finite(time), time >=  - FINAL_TIME_CAP, time <= FINAL_TIME_CAP, is.finite(score_z)) |>
		mutate(
			group = case_when(split == "prevalent" ~ "Prevalent (Yang)", split == "training" ~ "Incident: training", TRUE ~ "Incident: validation"),
			bin = floor(time / FINAL_TEMPORAL_STEP) * FINAL_TEMPORAL_STEP + FINAL_TEMPORAL_STEP / 2
		)
	ln <- d |>
		group_by(group, bin) |>
		summarise(mean = mean(score_z), N = n(), N_total = n_distinct(eid), .groups = "drop") |>
		mutate(reliable = N >= FINAL_MIN_BIN_N)
	h <- d |>
		count(bin, name = "n") |>
		transmute(bin, xmin = bin - FINAL_TEMPORAL_STEP / 2, xmax = bin + FINAL_TEMPORAL_STEP / 2, n)
	if (!nrow(ln)) return(blank_plot("Temporal score trajectory", "No sufficiently populated bins"))
	ln <- ln |>
		group_by(group) |>
		mutate(group_label = sprintf("%s (N=%s)", first(group), format(first(N_total), big.mark = ","))) |>
		ungroup()
	yr <- range(ln$mean, na.rm = TRUE) ; if (diff(yr) < .2) yr <- yr + c( - .5, .5) ; pad <- .13 * diff(yr) ; yr <- yr + c( - pad, pad)
	count_max <- max(h$n, 1, na.rm = TRUE) ; plot_max <- count_max * 1.15 ; scale_factor <- plot_max / diff(yr)
	ln <- ln |> mutate(y_plot = (mean - yr[1]) * scale_factor)
	lev <- unique(ln |> arrange(factor(group, levels = c("Prevalent (Yang)", "Incident: training", "Incident: validation"))) |> pull(group_label)) ; base <- sub(" \\(N=.*$", "", lev) ; pal0 <- c(`Prevalent (Yang)` = "grey45", `Incident: training` = "#2C7FB8", `Incident: validation` = "#D95F02") ; pal <- setNames(unname(pal0[base]), lev)
	missing_side <- tibble(x = c( - FINAL_TIME_CAP / 2, FINAL_TIME_CAP / 2), label = c(if (any(d$split == "prevalent")) NA_character_ else "Prevalence trajectory unavailable", if (any(d$split != "prevalent")) NA_character_ else "Incidence trajectory unavailable")) |> filter(!is.na(label))
	ggplot() +
		geom_rect(data = h, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = n, fill = bin >= 0), alpha = .31, inherit.aes = FALSE) +
		scale_fill_manual(values = c(`TRUE` = "#78B7E5", `FALSE` = "grey60"), guide = "none") +
		geom_vline(xintercept = 0, linetype = 2, color = "grey38") +
		geom_hline(yintercept = (0 - yr[1]) * scale_factor, linetype = 3, color = "grey62") +
		geom_smooth(data = ln |> filter(reliable), aes(bin, y_plot, color = group_label, group = group_label), method = "loess", se = FALSE, span = .65, linewidth = 1.0) +
		geom_point(data = ln, aes(bin, y_plot, color = group_label, alpha = reliable), size = 1.45) +
		geom_text(data = missing_side, aes(x = x, y = plot_max * .52, label = label), inherit.aes = FALSE, color = "grey45", fontface = "bold", size = 3) +
		scale_alpha_manual(values = c(`TRUE` = .78, `FALSE` = .25), guide = "none") +
		scale_y_continuous(name = "Patient count", limits = c(0, plot_max), expand = expansion(mult = c(0, .02)), sec.axis = sec_axis( ~ .x / scale_factor + yr[1], name = "Biomarker score (s.d.)")) +
		scale_color_manual(values = pal) +
		coord_cartesian(xlim = c( - FINAL_TIME_CAP, FINAL_TIME_CAP)) +
		labs(x = "Recorded diagnosis relative to baseline blood draw (years)", color = NULL) +
		theme_5c(9) +
		guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
		theme(legend.position = "top")
}

final_row5_preclinical <- function(score_rows, label) {
	d <- score_rows |> filter(method == .env$label, is.finite(score_z))
	cases <- d |>
		filter(event == 1, is.finite(time), time >=  - FINAL_TIME_CAP, time <= FINAL_TIME_CAP) |>
		mutate(
			group = case_when(split == "prevalent" ~ "Cases: prevalent", split == "training" ~ "Cases: training", TRUE ~ "Cases: validation"),
			bin = floor(time / FINAL_TEMPORAL_STEP) * FINAL_TEMPORAL_STEP + FINAL_TEMPORAL_STEP / 2
		)
	sm <- cases |>
		group_by(group, bin) |>
		summarise(mean = mean(score_z), se = sd(score_z) / sqrt(n()), N = n(), .groups = "drop") |>
		filter(N >= FINAL_MIN_BIN_N)
	if (!nrow(sm)) return(blank_plot("Baseline-anchored native score", "No sufficiently populated case bins"))
	ctrl <- d |>
		filter(event == 0, split == "validation") |>
		summarise(mean = mean(score_z, na.rm = TRUE), se = sd(score_z, na.rm = TRUE) / sqrt(sum(is.finite(score_z))))
	nlab <- cases |>
		group_by(group) |>
		summarise(N = n_distinct(eid), .groups = "drop") |>
		mutate(group_label = paste0(group, " (N=", format(N, big.mark = ","), ")"))
	sm <- sm |> left_join(nlab, by = "group")
	# Standardized scores are comparable across Cox and binary engines. They are
	# deliberately not labelled as calibrated probabilities.
	ylo <- min(sm$mean - sm$se, na.rm = TRUE) ; yhi <- max(sm$mean + sm$se, na.rm = TRUE)
	if (!all(is.finite(c(ylo, yhi)))) c(ylo, yhi) <- c(0, 1)
	if (yhi - ylo < .02) {
		mid <- (ylo + yhi) / 2 ; ylo <- mid - .01 ; yhi <- mid + .01
	}
	ypad <- .12 * (yhi - ylo) ; ylim_native <- c(ylo - ypad, yhi + ypad)
	gl <- unique(nlab$group_label) ; base_group <- sub(" \\(N=.*$", "", gl) ; pal0 <- c(`Cases: prevalent` = "grey45", `Cases: training` = "#2C7FB8", `Cases: validation` = "#D95F02") ; pal <- setNames(unname(pal0[base_group]), gl)
	missing_side <- tibble(x = c( - FINAL_TIME_CAP / 2, FINAL_TIME_CAP / 2), label = c(if (any(cases$split == "prevalent")) NA_character_ else "Prevalence score unavailable", if (any(cases$split != "prevalent")) NA_character_ else "Incidence score unavailable")) |> filter(!is.na(label))
	ggplot(sm, aes(bin, mean, color = group_label, group = group_label)) +
		geom_hline(data = ctrl, aes(yintercept = mean), inherit.aes = FALSE, linetype = 3, color = "grey45") +
		geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .14, linewidth = .35, alpha = .55) +
		geom_point(size = 1.8) +
		geom_smooth(method = "loess", se = FALSE, span = .75, linewidth = 1.0) +
		geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
		geom_text(data = missing_side, aes(x = x, y = mean(ylim_native), label = label), inherit.aes = FALSE, color = "grey45", fontface = "bold", size = 3) +
		scale_color_manual(values = pal) +
		coord_cartesian(xlim = c( - FINAL_TIME_CAP, FINAL_TIME_CAP), ylim = ylim_native) +
		labs(
			x = "Recorded diagnosis relative to baseline blood draw (years)", y = "Model score (training s.d.)", color = NULL,
			subtitle = "Dotted line: validation-control mean"
		) +
		theme_5c(9) +
		guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
		theme(legend.position = "top", legend.text = element_text(size = 7.4, face = "bold"))
}

final_make_5row_grid <- function(labels, pred, score_rows, pairs_native, boot_by_method, title = NULL) {
	if (!length(labels)) return(blank_plot(title %||% "Final prediction", "No score was available"))
	# Preserve the requested grid even when a model is unavailable; otherwise an
	# empty C4 cache silently collapses Fig3 into one completely blank canvas.
	nsel <- score_rows |>
		group_by(method) |>
		summarise(
			n_selected = as.integer(median(n_selected, na.rm = TRUE)),
			n_features = as.integer(median(n_features, na.rm = TRUE)), .groups = "drop"
		)
	getn <- function(lb) {
		col <- if (lb == "User specified") "n_features" else "n_selected" ; z <- nsel[[col]][match(lb, nsel$method)] ; if (length(z) && is.finite(z)) z else NA_integer_
	}
	method_rows <- lapply(labels, function(lb) wrap_plots(list(
		final_row1_score_profile(score_rows, lb, getn(lb)),
		final_row2_roc(pred, lb, boot_by_method[[lb]] %||% tibble()),
		final_row3_followup(score_rows, lb),
		final_row4_yy(score_rows, lb),
		final_row5_preclinical(score_rows, lb)
	), nrow = 1, widths = c(1.15, 1.05, 1, 1, 1.08)))
	spaced <- list() ; for (i in seq_along(method_rows)) {
		spaced[[length(spaced) + 1]] <- method_rows[[i]] ; if (i < length(method_rows)) spaced[[length(spaced) + 1]] <- plot_spacer()
	}
	p <- wrap_plots(spaced, ncol = 1, heights = rep(c(1, .08), length.out = length(spaced)))
	if (!is.null(title)) p <- p + plot_annotation(title = title, theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = .5)))
	p
}


# 🚩 Lead-time and consolidation figures (Figs 4-10)
make_individual_lead_rows <- function(dat, features, screen, split, tvar, evar) {
	features <- intersect(unique(features), names(dat)) ; it <- which(split == "training") ; iv <- which(split == "validation")
	map_dfr(features, function(x) {
		tr <- suppressWarnings(as.numeric(dat[[x]][it])) ; va <- suppressWarnings(as.numeric(dat[[x]][iv])) ; m <- mean(tr, na.rm = TRUE) ; s <- sd(tr, na.rm = TRUE)
		if (!is.finite(s) || s == 0) return(tibble()) ; sgn <- sign(screen$beta[match(x, screen$term)]) ; if (!is.finite(sgn) || sgn == 0) sgn <- 1
		tibble(
			eid = dat$eid[iv], split = "validation", time = dat[[tvar]][iv], event = dat[[evar]][iv], score_z = (va - m) / s * sgn,
			method = x, kind = "Single biomarker", n_features = 1L, n_selected = 1L
		)
	})
}

leadtime_auc_table <- function(score_rows, methods, max_h = FINAL_TIME_CAP) {
	if (!requireNamespace("pROC", quietly = TRUE)) return(tibble())
	if (!"kind" %in% names(score_rows)) score_rows$kind <- "Omic score"
	hs <- sort(unique(c(.5, seq(1, max(1, floor(max_h)), by = 1))))
	map_dfr(intersect(methods, unique(score_rows$method)), function(mm) {
		d0 <- score_rows |> filter(method == mm, split == "validation", is.finite(time), is.finite(score_z), !is.na(event))
		kind0 <- first(d0$kind) %||% "Omic score"
		map_dfr(hs, function(h) {
			# Landmark-style discrimination: incident cases diagnosed at least h years
			# after baseline versus controls observed event-free for at least h years.
			d <- d0 |> filter(time >= h)
			nc <- sum(d$event == 1) ; nn <- sum(d$event == 0)
			if (nc < 20 || nn < 50) return(tibble(
				method = mm, kind = kind0, horizon = h, cases = nc, controls = nn, AUC = NA_real_, lo = NA_real_, hi = NA_real_,
				AUC_definition = "minimum-lead-time case/control AUC"
			))
			roc <- tryCatch(pROC::roc(d$event, d$score_z, quiet = TRUE, direction = "<"), error = function(e) NULL)
			if (is.null(roc)) return(tibble(
				method = mm, kind = kind0, horizon = h, cases = nc, controls = nn, AUC = NA_real_, lo = NA_real_, hi = NA_real_,
				AUC_definition = "minimum-lead-time case/control AUC"
			))
			ci <- tryCatch(as.numeric(pROC::ci.auc(roc, method = "delong")), error = function(e) c(NA_real_, NA_real_, NA_real_))
			tibble(
				method = mm, kind = kind0, horizon = h, cases = nc, controls = nn, AUC = as.numeric(roc$auc), lo = ci[1], hi = ci[3],
				AUC_definition = "minimum-lead-time case/control AUC"
			)
		})
	})
}

leadtime_headline <- function(lead, min_cases = FINAL_LEAD_MIN_CASES) {
	if (!nrow(lead) || !all(c("method", "kind", "horizon", "AUC", "lo", "hi", "cases") %in% names(lead)))
		return(tibble(
			method = character(), kind = character(), last_supported_horizon = numeric(),
			AUC_at_horizon = numeric(), lo_at_horizon = numeric(), hi_at_horizon = numeric(), cases_at_horizon = integer()
		))
	lead |>
		group_by(method, kind) |>
		group_modify(function(z, key) {
			z <- z |> arrange(horizon) ; supported <- is.finite(z$AUC) & is.finite(z$lo) & z$lo > .5 & z$cases >= min_cases
			first_fail <- match(FALSE, supported, nomatch = length(supported) + 1L) ; last <- first_fail - 1L
			if (last < 1) return(tibble()) ; tibble(last_supported_horizon = z$horizon[last], AUC_at_horizon = z$AUC[last], lo_at_horizon = z$lo[last], hi_at_horizon = z$hi[last], cases_at_horizon = z$cases[last])
		}) |>
		ungroup()
}

prediction_window_auc_table <- function(score_rows, methods) {
	if (!requireNamespace("pROC", quietly = TRUE)) return(tibble())
	if (!"kind" %in% names(score_rows)) score_rows$kind <- "Omic score"
	windows <- tibble(window = c("Within 5 y", "Within 10 y", "Over 10 y", "Over 12 y"), type = c("within", "within", "over", "over"), cut = c(5, 10, 10, 12))
	ans <- map_dfr(intersect(methods, unique(score_rows$method)), function(mm) {
		d0 <- score_rows |> filter(method == mm, split == "validation", is.finite(time), is.finite(score_z), !is.na(event)) ; kind0 <- first(d0$kind) %||% "Omic score"
		map_dfr(seq_len(nrow(windows)), function(i) {
			w <- windows[i, ] ; type0 <- w$type[[1]] ; cut0 <- w$cut[[1]] ; window0 <- w$window[[1]]
			d <- if (type0 == "within") d0 |>
				mutate(class = case_when(event == 1 & time <= cut0 ~ 1, time > cut0 ~ 0, TRUE ~ NA_real_)) |>
				filter(!is.na(class)) else d0 |>
				filter(time > cut0) |>
				mutate(class = event)
			nc <- sum(d$class == 1) ; nn <- sum(d$class == 0) ; if (nc < 20 || nn < 50) return(tibble(
				method = mm, kind = kind0, window = window0, cases = nc, controls = nn, AUC = NA_real_, lo = NA_real_, hi = NA_real_,
				AUC_definition = "prespecified binary prediction-window AUC"
			))
			roc <- tryCatch(pROC::roc(d$class, d$score_z, quiet = TRUE, direction = "<"), error = function(e) NULL) ; if (is.null(roc)) return(tibble())
			ci <- tryCatch(as.numeric(pROC::ci.auc(roc, method = "delong")), error = function(e) c(NA_real_, NA_real_, NA_real_))
			tibble(
				method = mm, kind = kind0, window = window0, cases = nc, controls = nn, AUC = as.numeric(roc$auc), lo = ci[1], hi = ci[3],
				AUC_definition = "prespecified binary prediction-window AUC"
			)
		})
	})
	if (!nrow(ans)) return(tibble(
		method = character(), kind = character(), window = factor(levels = windows$window),
		cases = integer(), controls = integer(), AUC = numeric(), lo = numeric(), hi = numeric(), AUC_definition = character()
	))
	ans |> mutate(window = factor(window, levels = windows$window))
}

plot_leadtime_prediction <- function(lead, windows = tibble()) {
	if (!nrow(lead) || !any(is.finite(lead$AUC))) return(blank_plot("Prediction before onset", "No lead-time estimate met the case/control requirement"))
	ok <- lead |> filter(is.finite(AUC)) ; hmax <- leadtime_headline(lead)
	curve_panel <- function(kind0, title) {
		d <- ok |> filter(kind == kind0) ; hm <- hmax |> filter(kind == kind0)
		if (!nrow(d)) return(blank_plot(title, "No eligible method"))
		ggplot(d, aes(horizon, AUC)) +
			geom_hline(yintercept = .5, linetype = 2, color = "grey55") +
			geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#6BAED6", alpha = .20) +
			geom_line(color = "#2C7FB8", linewidth = .85) +
			geom_point(color = "#2C7FB8", size = 1.25) +
			geom_vline(data = hm, aes(xintercept = last_supported_horizon), color = "#D7301F", linetype = 2, linewidth = .8) +
			geom_point(data = hm, aes(last_supported_horizon, AUC_at_horizon), color = "#D7301F", size = 3) +
			geom_label(
				data = hm, aes(last_supported_horizon, AUC_at_horizon, label = paste0("up to ", format(last_supported_horizon, trim = TRUE), " y\n", cases_at_horizon, " cases")),
				color = "#D7301F", fill = "white", label.size = .18, size = 2.5, fontface = "bold", nudge_y = .075
			) +
			facet_wrap( ~ method, ncol = 2) +
			scale_x_continuous(breaks = c(.5, 2, 5, 8, 10, 12, 14, 16)) +
			scale_y_continuous(limits = c(.45, 1)) +
			labs(
				title = title, subtitle = "Red marker = farthest consecutive horizon passing the prespecified criterion",
				x = "Minimum years from blood draw to recorded diagnosis", y = "Held-out AUC (95% CI)"
			) +
			theme_5c(8)
	}
	pA <- curve_panel("Omic score", "A. How many years ahead do scores discriminate?")
	pB <- curve_panel("Single biomarker", "B. How many years ahead do individual biomarkers discriminate?")
	pC <- if (!nrow(hmax)) blank_plot("C. Supported lead time", paste0("No method had lower 95% AUC CI > 0.50 with at least ", FINAL_LEAD_MIN_CASES, " cases")) else
		ggplot(hmax, aes(last_supported_horizon, fct_reorder(method, last_supported_horizon))) +
			geom_segment(aes(x = 0, xend = last_supported_horizon, yend = fct_reorder(method, last_supported_horizon)), color = "grey75", linewidth = 2.2) +
			geom_point(color = "#D7301F", size = 3.2) +
			geom_text(aes(label = paste0(" ", format(last_supported_horizon, trim = TRUE), " y; AUC ", sprintf("%.2f", AUC_at_horizon), "; n=", cases_at_horizon)), hjust = 0, size = 2.8, fontface = "bold") +
			facet_grid(kind ~ ., scales = "free_y", space = "free_y") +
			coord_cartesian(xlim = c(0, max(hmax$last_supported_horizon) * 1.55), clip = "off") +
			labs(title = "C. Headline lead time (red = supported, not maximum follow-up)", x = "Years before recorded diagnosis", y = NULL) +
			theme_5c(9)
	counts <- lead |>
		group_by(horizon) |>
		summarise(cases = max(cases, na.rm = TRUE), controls = max(controls, na.rm = TRUE), .groups = "drop") |>
		pivot_longer(c(cases, controls), names_to = "sample", values_to = "N")
	pD <- ggplot(counts, aes(horizon, N, color = sample)) +
		geom_vline(xintercept = c(.5, 1, 2, 5, 10), linetype = 3, color = "grey82") +
		geom_line(linewidth = .9) +
		geom_point(size = 1.5) +
		scale_color_manual(values = c(cases = "#D7301F", controls = "#2C7FB8"), labels = c(cases = "Incident cases", controls = "Eligible controls")) +
		labs(title = "E. Information remaining at each lead time", subtitle = paste0("Headline requires lower 95% AUC CI > 0.50, at least ", FINAL_LEAD_MIN_CASES, " cases, and no earlier failure"), x = "Minimum lead time (years)", y = "Participants", color = NULL) +
		theme_5c(9) +
		theme(legend.position = "top")
	ww <- windows |> filter(kind == "Omic score", is.finite(AUC))
	pW <- if (!nrow(ww)) blank_plot("D. Prespecified prediction windows", "Window-specific AUC was unavailable") else ggplot(ww, aes(AUC, fct_reorder(method, AUC), color = window)) +
		geom_errorbarh(aes(xmin = lo, xmax = hi), height = .15, position = position_dodge(width = .55)) +
		geom_point(position = position_dodge(width = .55), size = 2) +
		geom_vline(xintercept = .5, linetype = 3, color = "grey65") +
		scale_x_continuous(limits = c(.45, 1)) +
		labs(title = "D. MASLD-style prespecified windows", subtitle = "Within 5/10 years and diagnosis more than 10/12 years after baseline", x = "Held-out AUC (95% CI)", y = NULL, color = NULL) +
		theme_5c(8) +
		theme(legend.position = "top")
	(pA | pB) / plot_spacer() / (pC | pW | pD) + plot_layout(heights = c(1.35, .07, 1))
}

topn_trajectory_summary <- function(rows, max_cases = FINAL_PAIR_MAX_CASES, step = FINAL_TEMPORAL_STEP) {
	if (!nrow(rows)) return(tibble())
	map_dfr(unique(rows$model), function(mm) {
		d <- rows |> filter(model == mm, is.finite(score), is.finite(time), !is.na(event)) ; ca <- d |> filter(event == 1, time <= FINAL_TIME_CAP)
		if (!nrow(ca) || nrow(d) < 2) return(tibble()) ; set.seed(SEED + 810 + match(mm, unique(rows$model)))
		if (nrow(ca) > max_cases) ca <- ca |> slice_sample(n = max_cases)
		# Incidence-density sampling: a comparator may be event-free or may become
		# a case later, provided they are still under observation at the index
		# case's diagnosis time.  This estimates score separation within the same
		# risk set instead of contrasting only with never-cases.
		pairs <- map_dfr(seq_len(nrow(ca)), function(i) {
			pool <- which(d$eid != ca$eid[i] & d$time >= ca$time[i]) ; if (!length(pool)) return(tibble()) ; j <- sample(pool, 1)
			tibble(years_before = ca$time[i], difference = ca$score[i] - d$score[j])
		})
		pairs |>
			mutate(bin = floor(years_before / step) * step + step / 2) |>
			group_by(bin) |>
			summarise(mean_difference = mean(difference), se = sd(difference) / sqrt(n()), N = n(), .groups = "drop") |>
			filter(N >= 10) |>
			mutate(model = mm, lo = mean_difference - 1.96 * se, hi = mean_difference + 1.96 * se)
	})
}

plot_topn_mechanism <- function(topn) {
	z <- topn$summary ; ir <- topn$individual_range
	if (!nrow(z)) return(blank_plot("Why combine biomarkers?", "Top-N benchmark was unavailable"))
	pal <- c(`Top-1 biomarker (locked)` = "#595959", `Unweighted top-N mean` = "#2C7FB8", `Marginal-Cox-weighted top-N` = "#D95F02")
	pA <- ggplot(z, aes(N, AUC, color = model, fill = model)) +
		geom_hline(yintercept = .5, linetype = 3, color = "grey65") +
		{
			if (nrow(ir)) geom_ribbon(data = ir, aes(N, ymin = AUC_q25, ymax = AUC_q75), inherit.aes = FALSE, fill = "grey72", alpha = .35) else geom_blank()
		} +
		{
			if (nrow(ir)) geom_line(data = ir, aes(N, AUC_max), inherit.aes = FALSE, color = "grey45", linetype = 2, linewidth = .65) else geom_blank()
		} +
		geom_ribbon(aes(ymin = AUC_lo, ymax = AUC_hi), alpha = .08, color = NA) +
		geom_line(linewidth = .9) +
		geom_point(size = 1.45) +
		scale_color_manual(values = pal) +
		scale_fill_manual(values = pal) +
		labs(
			title = "A. Held-out discrimination",
			subtitle = "Grey band: IQR of the N individual biomarkers; dashed grey: their optimistic validation maximum",
			x = "Training-ranked biomarkers included (N)", y = "10-year IPCW AUC", color = NULL, fill = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "top")
	pB <- ggplot(z, aes(N, beta, color = model, group = model)) +
		geom_hline(yintercept = 0, linetype = 3, color = "grey65") +
		geom_ribbon(aes(ymin = beta_lo, ymax = beta_hi, fill = model), alpha = .08, color = NA) +
		geom_line(linewidth = .9) +
		geom_point(size = 1.45) +
		scale_color_manual(values = pal) +
		scale_fill_manual(values = pal) +
		labs(
			title = "B. Effect per training SD",
			subtitle = "Validation Cox beta adjusted for the same clinical covariates", x = "N", y = "Log HR per 1-SD score", color = NULL, fill = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "none")
	sep <- z |>
		select(N, model, mean_separation, within_group_SD, cohen_d) |>
		pivot_longer(c(mean_separation, within_group_SD, cohen_d), names_to = "metric", values_to = "value") |>
		mutate(metric = recode(metric, mean_separation = "Case-control mean separation", within_group_SD = "Pooled within-group SD", cohen_d = "Signal/noise (Cohen d)"))
	pC <- ggplot(sep, aes(N, value, color = model)) +
		geom_line(linewidth = .85) +
		geom_point(size = 1.3) +
		facet_wrap( ~ metric, ncol = 1, scales = "free_y") +
		scale_color_manual(values = pal) +
		labs(
			title = "C. What changes when signals are combined?",
			subtitle = "A score helps when separation grows relative to within-group noise; raw SD alone is not the mechanism",
			x = "N", y = NULL, color = NULL
		) +
		theme_5c(8) +
		theme(legend.position = "none")
	stab <- z |>
		select(N, model, AUC_boot_SD, separation_boot_SD) |>
		pivot_longer(c(AUC_boot_SD, separation_boot_SD), names_to = "metric", values_to = "value") |>
		filter(is.finite(value)) |>
		mutate(metric = recode(metric, AUC_boot_SD = "Bootstrap SD of AUC", separation_boot_SD = "Bootstrap SD of mean separation"))
	pD <- if (!nrow(stab)) blank_plot("D. Sampling stability", "Bootstrap estimates unavailable") else ggplot(stab, aes(N, value, color = model)) +
		geom_line(linewidth = .85) +
		geom_point(size = 1.3) +
		facet_wrap( ~ metric, ncol = 1, scales = "free_y") +
		scale_color_manual(values = pal) +
		labs(title = "D. Sampling stability", x = "N", y = "Bootstrap SD", color = NULL) +
		theme_5c(8) +
		theme(legend.position = "none")
	tr <- topn_trajectory_summary(topn$rows)
	pE <- if (!nrow(tr)) blank_plot("E. Diagnosis-anchored separation", "Trajectory rows unavailable") else ggplot(tr, aes(bin, mean_difference, color = model, fill = model)) +
		geom_hline(yintercept = 0, linetype = 3, color = "grey60") +
		geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .10, color = NA) +
		geom_line(linewidth = .95) +
		geom_point(size = 1.3) +
		scale_color_manual(values = pal) +
		scale_fill_manual(values = pal) +
		scale_x_reverse(limits = c(FINAL_TIME_CAP, 0), breaks = c(16, 12, 10, 5, 2, 0)) +
		labs(
			title = paste0("E. Score separation up to diagnosis (N=", unique(topn$rows$N)[1], ")"),
			subtitle = "Incidence-density risk-set matching; cross-sectional baselines in different people, not within-person trajectories",
			x = "Years before the recorded diagnosis", y = "Risk-set case − comparator score (95% CI)", color = NULL, fill = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "top")
	div <- z |>
		filter(model == "Marginal-Cox-weighted top-N") |>
		select(N, average_abs_correlation, effective_N) |>
		pivot_longer(c(average_abs_correlation, effective_N), names_to = "metric", values_to = "value") |>
		mutate(metric = recode(metric, average_abs_correlation = "Average |correlation|", effective_N = "Weight-based effective N"))
	pF <- ggplot(div, aes(N, value)) +
		geom_line(color = "#6A3D9A", linewidth = .9) +
		geom_point(color = "#6A3D9A", size = 1.5) +
		facet_wrap( ~ metric, ncol = 1, scales = "free_y") +
		labs(title = "F. Redundancy versus diversification", subtitle = "Highly correlated biomarkers add less independent information", x = "N", y = NULL) +
		theme_5c(8)
	(pA | pB | pC) / plot_spacer() / (pD | pE | pF) + plot_layout(heights = c(1, .07, 1))
}

attained_age_score_sensitivity <- function(score_rows, dat, methods, clinical, tvar, evar) {
	empty <- tibble(
		method = character(), time_scale = character(), beta = numeric(), se = numeric(),
		N = integer(), events = integer(), lo = numeric(), hi = numeric()
	)
	eligible_methods <- intersect(methods, unique(score_rows$method))
	if (!length(eligible_methods)) return(empty)
	age_covs <- unique(c("age", grep("^age($|[._])", clinical, value = TRUE, ignore.case = TRUE))) ; att_covs <- setdiff(clinical, age_covs)
	# Join by participant identifier rather than cached row position: an older
	# reusable score-model cache remains valid even if phenotype columns change.
	base_dat <- dat |> select(eid, all_of(unique(c(tvar, evar, ".attained_entry", ".attained_exit", clinical))))
	map_dfr(eligible_methods, function(mm) {
		d <- score_rows |>
			filter(method == mm, split == "validation", is.finite(score_z)) |>
			select(eid, score_z) |>
			inner_join(base_dat, by = "eid")
		fit_one <- function(kind) {
			covs <- if (kind == "Baseline time") clinical else att_covs ; surv <- if (kind == "Baseline time") paste0("Surv(", bt(tvar), ",", bt(evar), ")") else paste0("Surv(", bt(".attained_entry"), ",", bt(".attained_exit"), ",", bt(evar), ")")
			time_cols <- if (kind == "Baseline time") tvar else c(".attained_entry", ".attained_exit")
			rhs <- c("score_z", covs) ; dd <- d[, unique(c(rhs, time_cols, evar)), drop = FALSE] ; dd <- dd[complete.cases(dd), , drop = FALSE]
			if (kind != "Baseline time") dd <- dd[dd$.attained_exit > dd$.attained_entry, , drop = FALSE]
			if (nrow(dd) < 100 || sum(dd[[evar]] == 1) < 20) {
				message(
					"C1 prediction: score time-scale sensitivity unavailable for ", mm, " / ", kind,
					"; complete-case N=", nrow(dd), ", events=", sum(dd[[evar]] == 1), " (requires N >= 100 and events >= 20)"
				)
				return(empty)
			}
			fit <- tryCatch(coxph(as.formula(paste0(surv, " ~ ", paste(bt(rhs), collapse = " + "))), dd, ties = "efron"),
				error = function(e) {
					message("C1 prediction: score time-scale sensitivity fit failed for ", mm, " / ", kind, ": ", conditionMessage(e)) ; NULL
				}
			)
			sm <- if (is.null(fit)) NULL else coef(summary(fit)) ; if (is.null(sm) || !"score_z" %in% rownames(sm)) return(empty)
			tibble(method = mm, time_scale = kind, beta = sm["score_z", "coef"], se = sm["score_z", "se(coef)"], N = nrow(dd), events = sum(dd[[evar]] == 1))
		}
		bind_rows(fit_one("Baseline time"), fit_one("Attained age, delayed entry"))
	}) |> mutate(lo = beta - 1.96 * se, hi = beta + 1.96 * se)
}

plot_attained_age_score_sensitivity <- function(z) {
	if (!nrow(z)) return(blank_plot("Score time-scale sensitivity", "No attained-age score model was estimable"))
	ggplot(z, aes(beta, fct_reorder(method, beta), color = time_scale)) +
		geom_vline(xintercept = 0, linetype = 3, color = "grey65") +
		geom_errorbarh(aes(xmin = lo, xmax = hi), height = .16, position = position_dodge(width = .45)) +
		geom_point(position = position_dodge(width = .45), size = 2.2) +
		scale_color_manual(values = c(`Baseline time` = "#2C7FB8", `Attained age, delayed entry` = "#D95F02")) +
		labs(
			title = "Prediction-score sensitivity to the survival time scale",
			subtitle = "Attained-age follow-up starts at baseline age; adult omics are not treated as birth measurements",
			x = "Validation log HR per training-SD score", y = NULL, color = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "top")
}

plot_evidence_matrix <- function(ev) {
	if (!nrow(ev)) return(blank_plot("5C evidence matrix", "No upstream evidence table was available"))
	grade_levels <- c("A", "B", "C", "D")
	keep <- ev |>
		mutate(.grade = match(evidence_grade, grade_levels)) |>
		arrange(.grade, desc(evidence_count), obs_p) |>
		slice_head(n = 40) |>
		pull(feature)
	z0 <- ev |>
		filter(feature %in% keep) |>
		transmute(feature,
			`C1 correlation` = C1, `Distal support` = as.integer(distal_support), `C2 MR` = C2_MR, `Reverse MR` = C2_REVERSE_MR, `C2 DPG` = C2_DANDELION,
			`C3 coloc` = C3, `C4 connection` = C4, `Final stability` = ifelse(FINAL_available, C5, NA_integer_),
			`Reactive signal` =  - as.integer(coalesce(reactive_compatible, FALSE))
		) |>
		pivot_longer( - feature, names_to = "stage", values_to = "support")
	grade <- ev |>
		filter(feature %in% keep) |>
		transmute(feature, stage = "Grade", support = NA_integer_, grade_label = evidence_grade)
	z <- bind_rows(z0 |> mutate(grade_label = NA_character_), grade) |>
		mutate(
			feature = factor(feature, levels = rev(keep)),
			stage = factor(stage, levels = c("C1 correlation", "Distal support", "C2 MR", "Reverse MR", "C2 DPG", "C3 coloc", "C4 connection", "Final stability", "Reactive signal", "Grade"))
		)
	ggplot(z, aes(stage, feature, fill = factor(support))) +
		geom_tile(color = "white", linewidth = .35) +
		geom_text(data = z |> filter(stage == "Grade"), aes(label = grade_label), fontface = "bold", size = 3.1) +
		scale_fill_manual(
			values = c(`-1` = "#D9A441", `0` = "grey92", `1` = "#3B8064"), na.value = "white",
			name = "Evidence", labels = c(`-1` = "Reactive-compatible", `0` = "No", `1` = "Yes")
		) +
		labs(
			title = "5C evidence consolidation with causal safeguards",
			subtitle = "Reactive evidence changes the biomarker role; it no longer vetoes concordant MR + robust colocalization",
			x = NULL, y = NULL
		) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(), legend.position = "top")
}

plot_causal_reactive_map <- function(ev) {
	if (!nrow(ev)) return(blank_plot("Causal–reactive map", "No evidence table"))
	z <- ev |> mutate(
		genetic_strength = coalesce(pmin(12, - log10(pmax(MR_p, 1e-300))), 0) + 4 * pmax(0, pmin(1, coalesce(PP.H4_robust_min, 0))),
		reactive_strength = pmax(0, abs(coalesce(prevalent_score, 0))) + pmax(0, abs(coalesce(duration_score, 0))),
		label = ifelse(feature %in% c("GDF15", "NTPROBNP", "NPPB", "PCSK9") | evidence_grade %in% c("A", "B") & min_rank(obs_p) <= 18, feature, NA_character_)
	)
	ggplot(z, aes(genetic_strength, reactive_strength, color = consolidation_role)) +
		geom_point(aes(size = pmin(12, - log10(pmax(obs_p, 1e-300)))), alpha = .72) +
		geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 25, seed = SEED, fontface = "bold") +
		labs(
			title = "Causal–reactive evidence map", subtitle = "Upper-right biomarkers may be both causal and disease-responsive; neither axis alone proves mechanism",
			x = "Genetic causal/locus evidence", y = "Established-disease/reactive evidence", color = NULL, size = "Incident\nevidence"
		) +
		theme_5c(9) +
		theme(legend.position = "bottom")
}

plot_performance_benchmark <- function(psum) {
	if (!nrow(psum)) return(blank_plot("Held-out prediction benchmark", "No validation summary was available"))
	z <- psum |>
		select(biom_set, model, AUC, C_index) |>
		pivot_longer(c(AUC, C_index), names_to = "metric", values_to = "value") |>
		filter(is.finite(value))
	ggplot(z, aes(value, fct_reorder(biom_set, value, max), color = model, shape = model)) +
		geom_vline(xintercept = .5, linetype = 3, color = "grey65") +
		geom_point(size = 2.4, position = position_dodge(width = .38)) +
		facet_wrap( ~ metric, scales = "free_x") +
		scale_color_manual(values = c(`Model 0` = "grey35", Biomarkers = "#C77732", Combined = "#287A58"), labels = c(`Model 0` = "Clinical")) +
		labs(title = "Held-out prediction benchmark", subtitle = "All comparisons use the same outer validation split", x = "Performance", y = NULL, color = NULL, shape = NULL) +
		theme_5c(9) +
		theme(legend.position = "top")
}

plot_incremental_performance <- function(psum) {
	if (!nrow(psum)) return(blank_plot("Incremental discrimination", "No validation summary was available"))
	z <- psum |>
		group_by(biom_set) |>
		summarise(
			clinical_AUC = AUC[match("Model 0", model)], biom_AUC = AUC[match("Biomarkers", model)], combined_AUC = AUC[match("Combined", model)],
			clinical_C = C_index[match("Model 0", model)], biom_C = C_index[match("Biomarkers", model)], combined_C = C_index[match("Combined", model)], .groups = "drop"
		) |>
		mutate(
			delta_AUC_combined = combined_AUC - clinical_AUC, delta_C_combined = combined_C - clinical_C,
			delta_AUC_biom = biom_AUC - clinical_AUC, delta_C_biom = biom_C - clinical_C
		) |>
		select(biom_set, starts_with("delta_")) |>
		pivot_longer( - biom_set, names_to = "contrast", values_to = "delta") |>
		filter(is.finite(delta)) |>
		mutate(metric = ifelse(str_detect(contrast, "AUC"), "ΔAUC", "ΔC-index"), model = ifelse(str_detect(contrast, "combined"), "Clinical + omics", "Omics only"))
	if (!nrow(z)) return(blank_plot("Incremental discrimination", "Could not align clinical, omics, and combined models"))
	ggplot(z, aes(delta, fct_reorder(biom_set, delta, max), color = model, shape = model)) +
		geom_vline(xintercept = 0, color = "grey60") +
		geom_point(size = 2.5, position = position_dodge(width = .35)) +
		facet_wrap( ~ metric, scales = "free_x") +
		labs(title = "Incremental discrimination over Model 0", x = "Performance difference", y = NULL, color = NULL, shape = NULL) +
		theme_5c(9) +
		theme(legend.position = "top")
}

plot_complexity_performance <- function(psum) {
	z <- psum |> filter(model %in% c("Biomarkers", "Combined"), is.finite(AUC), is.finite(n_selected), n_selected > 0)
	if (!nrow(z)) return(blank_plot("Parsimony versus discrimination", "No finite feature-count/performance pairs"))
	z <- z |> mutate(label = ifelse(min_rank(desc(AUC)) <= 6 | n_selected <= 5, biom_set, NA_character_))
	ggplot(z, aes(n_selected, AUC, color = model)) +
		geom_hline(yintercept = .5, linetype = 3, color = "grey65") +
		geom_point(size = 2.6) +
		ggrepel::geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 15, seed = SEED, show.legend = FALSE) +
		scale_x_log10() +
		labs(title = "Model parsimony versus validation AUC", x = "Selected biomarkers (log scale)", y = "AUC", color = NULL) +
		theme_5c(10) +
		theme(legend.position = "top")
}

score_correlation_table <- function(score_rows) {
	w <- score_rows |>
		filter(split == "validation", is.finite(score_z)) |>
		select(eid, method, score_z) |>
		distinct() |>
		pivot_wider(names_from = method, values_from = score_z)
	if (nrow(w) < 20 || ncol(w) < 3) return(tibble())
	m <- suppressWarnings(cor(as.data.frame(w[, - 1, drop = FALSE]), use = "pairwise.complete.obs", method = "spearman"))
	as.data.frame(as.table(m), stringsAsFactors = FALSE) |>
		as_tibble() |>
		rename(method1 = Var1, method2 = Var2, rho = Freq)
}

plot_score_correlation <- function(corr) {
	if (!nrow(corr)) return(blank_plot("Prediction-score concordance", "Fewer than two validation scores were available"))
	ggplot(corr, aes(method1, method2, fill = rho)) +
		geom_tile(color = "white") +
		geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4) +
		scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0, limits = c( - 1, 1), name = "Spearman ρ") +
		labs(title = "Concordance among prediction paradigms", x = NULL, y = NULL) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

subgroup_auc_table <- function(score_rows, dat, methods) {
	keep <- intersect(c("eid", "age", "sex"), names(dat)) ; if (!all(c("eid", "age", "sex") %in% keep)) return(tibble())
	if (!length(intersect(methods, unique(score_rows$method)))) return(tibble())
	dd <- score_rows |>
		filter(method %in% methods, split == "validation", is.finite(score_z)) |>
		left_join(dat |> select(all_of(keep)), by = "eid")
	agecuts <- quantile(dd$age, c(0, 1 / 3, 2 / 3, 1), na.rm = TRUE) ; if (anyDuplicated(agecuts)) return(tibble())
	dd <- dd |> mutate(age_group = cut(age, agecuts, include.lowest = TRUE, dig.lab = 3), sex_group = paste0("Sex category ", sex))
	bind_rows(dd |> mutate(domain = "Age tertile", subgroup = as.character(age_group)), dd |> mutate(domain = "Reported sex", subgroup = sex_group)) |>
		group_by(method, domain, subgroup) |>
		group_modify( ~ {
			if (nrow(.x) < 100 || sum(.x$event == 1) < 15 || sum(.x$event == 0) < 30) return(tibble())
			auc <- weighted_time_auc(.x$time, .x$event, .x$score_z, HORIZON) ; if (!is.finite(auc)) return(tibble())
			set.seed(SEED + nrow(.x)) ; bs <- replicate(100, {
				ii <- sample(seq_len(nrow(.x)), nrow(.x), replace = TRUE) ; weighted_time_auc(.x$time[ii], .x$event[ii], .x$score_z[ii], HORIZON)
			})
			ci <- quantile(bs, c(.025, .975), na.rm = TRUE) ; tibble(N = nrow(.x), events = sum(.x$event == 1), AUC = auc, lo = ci[[1]], hi = ci[[2]])
		}) |>
		ungroup()
}

plot_subgroup_auc <- function(sg) {
	if (!nrow(sg)) return(blank_plot("Subgroup discrimination", "Age and sex subgroup estimates were unavailable"))
	ggplot(sg, aes(AUC, interaction(method, subgroup, lex.order = TRUE), color = method)) +
		geom_vline(xintercept = .5, linetype = 3, color = "grey65") +
		geom_errorbarh(aes(xmin = lo, xmax = hi), height = .12) +
		geom_point(size = 2) +
		facet_wrap( ~ domain, scales = "free_y", ncol = 1) +
		labs(title = "Validation discrimination across prespecified subgroups", subtitle = paste0(HORIZON, "-year IPCW AUC; reported sex is retained in its source coding"), x = "AUC (bootstrap 95% CI)", y = NULL, color = NULL) +
		theme_5c(8) +
		theme(legend.position = "top")
}


# 🚩 Evidence table: C2 includes MR and DANDELION, C3 coloc, C4 connection.
make_evidence_table <- function(layer, outdir, stability = tibble()) {
	rr <- function(job) {
		if (job != "c1_correlate") return(list(data = list(), available = FALSE, status = "outside C1 baseline"))
		f <- file.path(le8_job_dir(outdir, job), paste0(sub("_.*$", "", job), ".res.rds"))
		if (!file.exists(f)) return(list(data = list(), available = FALSE, status = "missing"))
		z <- tryCatch(readRDS(f), error = function(e) e)
		if (inherits(z, "condition")) return(list(data = list(), available = FALSE, status = paste0("unreadable: ", conditionMessage(z))))
		le8_check_options(z)
		declared <- str_to_lower(as.character(z$meta$status %||% "ok")[[1]])
		usable <- !declared %in% c("unavailable", "no qtl locus constructed", "failed", "not run")
		list(
			data = if (usable) z else list(), available = usable,
			status = if (usable) "available" else paste0("available file; ", declared)
		)
	}
	u1 <- rr("c1_correlate") ; u2 <- rr("c2_cause") ; u3 <- rr("c3_coloc") ; u4 <- rr("c4_connect")
	c1 <- u1$data ; c2 <- u2$data ; c3 <- u3$data ; c4 <- u4$data
	a0 <- c1$association %||% tibble() ; a <- if (nrow(a0)) a0 |> transmute(feature = term, obs_beta = beta, obs_p = p.value, obs_FDR = FDR) else tibble(feature = character())
	g0 <- c1$pgs_incident %||% tibble() ; g <- if (nrow(g0)) g0 |> transmute(
		feature = term,
		pgs_beta = beta, pgs_p = p.value, pgs_FDR = FDR, pgs_N = N_total
	) else tibble(feature = character())
	di0 <- c1$directionality %||% tibble() ; di <- if (nrow(di0) && all(c("term", "direction_class") %in% names(di0))) di0 |>
		transmute(
			feature = term, temporal_class = direction_class, reactive_compatible = coalesce(reactive_compatible, FALSE),
			distal_support = coalesce(distal_support, FALSE), FDR_landmark5 = FDR_landmark5,
			prevalent_score = prevalent_score, duration_score = duration_score, landmark5_score = landmark5_score
		) else
		tibble(feature = character(), temporal_class = character(), reactive_compatible = logical(), distal_support = logical(), FDR_landmark5 = numeric(), prevalent_score = numeric(), duration_score = numeric(), landmark5_score = numeric())
	m0 <- c2$MR %||% tibble() ; m <- if(nrow(m0)) le8_select_mr(m0,layer) |>
		transmute(feature=exposure,MR_analysis=analysis,MR_beta=b,MR_p=pval,MR_FDR=FDR_all) else tibble(feature=character())
	r0 <- c2$MR_reverse %||% tibble() ; rmr <- if (nrow(r0)) r0 |>
		filter(is.finite(pval)) |>
		transmute(feature, reverse_MR_beta = b, reverse_MR_p = pval, reverse_MR_FDR = FDR_reverse) else tibble(feature = character())
	dan <- c2$DANDELION %||% list() ; d0 <- as_tibble(dan$targets %||% tibble())
	if (nrow(d0) && !"consolidation_eligible" %in% names(d0)) {
		primary_input <- if ("gene_evidence_type" %in% names(d0))
			!str_detect(str_to_lower(coalesce(d0$gene_evidence_type, "")), "magma|adapted") else FALSE
		tested_n <- suppressWarnings(as.numeric((dan$ptrans_dimensions %||% c(NA_real_))[[1]]))
		frac <- if (is.finite(tested_n) && tested_n > 0) nrow(d0) / tested_n else NA_real_
		d0$consolidation_eligible <- primary_input & !(is.finite(frac) && frac > .25)
	}
	dd <- if (nrow(d0)) d0 |> transmute(
		feature = gene2, DANDELION_p, DANDELION_loci = n_distal_loci,
		DANDELION_sensitivity = TRUE, DANDELION_primary = as.logical(consolidation_eligible)
	) else
		tibble(feature = character(), DANDELION_p = numeric(), DANDELION_loci = integer(), DANDELION_sensitivity = logical(), DANDELION_primary = logical())
	co0 <- as_tibble(c3$summary %||% tibble())
	if (nrow(co0) && !"PP.H4_robust_min" %in% names(co0)) co0$PP.H4_robust_min <- NA_real_
	same_region <- le8_same_locus_evidence(m0,co0,layer)
	region <- if(nrow(same_region)) same_region |> group_by(feature) |> summarise(same_region_eligible=any(eligible %in% TRUE),.groups="drop") else tibble(feature=character(),same_region_eligible=logical())
	co <- if (nrow(co0) && all(c("feature", "PP.H4") %in% names(co0))) co0 |>
		filter(status == "ok") |>
		mutate(.robust = PP.H4_robust_min,policy_coloc_pass=le8_coloc_pass(pick(everything()))) |>
		group_by(feature) |>
		slice_max(.robust, n = 1, with_ties = FALSE) |>
		ungroup() |>
		select(feature, PP.H4, PP.H4_robust_min, everything()) else tibble(feature = character())
	cx0 <- c4$membership %||% tibble() ; cx <- if (nrow(cx0) && "feature" %in% names(cx0)) cx0 else tibble(feature = character())
	st <- if (nrow(stability) && all(c("feature", "selection_frequency") %in% names(stability))) stability |>
		group_by(feature) |>
		summarise(max_selection_frequency = max(selection_frequency, na.rm = TRUE), .groups = "drop") else tibble(feature = character(), max_selection_frequency = numeric())
	z <- Reduce(function(x, y) full_join(x, y, by = "feature"), list(a, g, di, m, rmr, dd, co, region, cx, st)) ; if (!nrow(z)) return(tibble())
	defaults <- list(
		same_region_eligible=FALSE, policy_coloc_pass=FALSE, obs_FDR = NA_real_, obs_p = NA_real_, MR_p = NA_real_, MR_FDR = NA_real_, DANDELION_sensitivity = FALSE,
		DANDELION_primary = FALSE, PP.H4 = NA_real_, PP.H4_robust_min = NA_real_, strict_YS = FALSE,
		max_selection_frequency = NA_real_, reactive_compatible = FALSE, distal_support = FALSE, reverse_MR_FDR = NA_real_,
		prevalent_score = NA_real_, duration_score = NA_real_, landmark5_score = NA_real_,
		pgs_beta = NA_real_, pgs_p = NA_real_, pgs_FDR = NA_real_, pgs_N = NA_real_
	)
	for (nm in names(defaults)) if (!nm %in% names(z)) z[[nm]] <- defaults[[nm]]
	z |>
		mutate(
			C1 = as.integer(is.finite(obs_FDR) & obs_FDR < .05),
			PGS_tested = is.finite(pgs_p), C1_PGS = ifelse(PGS_tested, as.integer(pgs_FDR < .05), NA_integer_),
			PGS_observed_sign_match = PGS_tested & is.finite(obs_beta) & sign(pgs_beta) == sign(obs_beta),
			C2_MR = if (u2$available) as.integer(is.finite(MR_FDR) & MR_FDR < le8_evidence_policy()$mr_fdr) else NA_integer_,
			C2_REVERSE_MR = if (u2$available) as.integer(is.finite(reverse_MR_FDR) & reverse_MR_FDR < le8_evidence_policy()$mr_fdr) else NA_integer_,
			C2_DANDELION = if (u2$available) as.integer(coalesce(DANDELION_primary, FALSE)) else NA_integer_,
			C2_DANDELION_sensitivity = if (u2$available) as.integer(coalesce(DANDELION_sensitivity, FALSE)) else NA_integer_,
			C2 = if (u2$available) as.integer(coalesce(C2_MR, 0L) == 1L | coalesce(C2_DANDELION, 0L) == 1L) else NA_integer_,
			C3 = if (u3$available) as.integer(coalesce(policy_coloc_pass,FALSE)) else NA_integer_,
			C4 = if (u4$available) as.integer(coalesce(strict_YS, FALSE)) else NA_integer_,
			C5 = as.integer(is.finite(max_selection_frequency) & max_selection_frequency >= .6),
			FINAL_available = is.finite(max_selection_frequency), reactive_compatible = coalesce(reactive_compatible, FALSE),
			distal_support = coalesce(distal_support, FALSE), temporal_warning = reactive_compatible,
			C2_status = u2$status, C3_status = u3$status, C4_status = u4$status,
			evidence_count = rowSums(cbind(C1, C2, C3, C4, ifelse(FINAL_available, C5, NA_integer_)), na.rm = TRUE),
			causal_locus_evidence = coalesce(C2_MR, 0L) == 1L | coalesce(C2_DANDELION, 0L) == 1L | coalesce(C3, 0L) == 1L,
			evidence_grade = case_when(
				same_region_eligible %in% TRUE & C1 == 1 ~ "A",
				causal_locus_evidence & (C1 == 1 | distal_support) ~ "B",
				evidence_count >= 2 ~ "C", TRUE ~ "D"
			),
			consolidation_role = case_when(
				evidence_grade %in% c("A", "B") & reactive_compatible ~ "Causal + reactive mixed biomarker",
				evidence_grade == "A" ~ "Cis/local region-supported candidate",
				coalesce(C2_REVERSE_MR, 0L) == 1L & !causal_locus_evidence ~ "Disease-liability-responsive biomarker",
				distal_support & !causal_locus_evidence ~ "Distal predictive biomarker",
				reactive_compatible & !causal_locus_evidence ~ "Reactive/diagnostic-compatible biomarker",
				evidence_grade == "B" ~ "Multi-domain causal/locus candidate",
				evidence_grade == "C" ~ "Multi-domain supported candidate", TRUE ~ "Single-domain or insufficient"
			),
			# This is an auditable prior for transparent penalized prediction, not a
			# fitted causal effect. Reactive-compatible features without locus evidence
			# receive the smallest weight; a name such as GDF15 is never hard-coded.
			evidence_weight_prior = case_when(
				causal_locus_evidence & distal_support & !reactive_compatible ~ 1.00,
				causal_locus_evidence & reactive_compatible ~ .75,
				causal_locus_evidence ~ .90,
				coalesce(C1_PGS, 0L) == 1L & distal_support ~ .80,
				distal_support ~ .65,
				reactive_compatible ~ .25,
				C1 == 1 ~ .50, TRUE ~ .35
			),
			prediction_role = case_when(
				causal_locus_evidence ~ "Mechanism-supported candidate",
				distal_support ~ "Distal prediction candidate",
				reactive_compatible ~ "Downweight; diagnostic/reactive-compatible",
				coalesce(C1_PGS, 0L) == 1L ~ "Inherited-score candidate; audit pleiotropy",
				TRUE ~ "Exploratory"
			)
		) |>
		arrange(factor(evidence_grade, levels = c("A", "B", "C", "D")), desc(evidence_count), obs_p)
}


# 🚩 Main C1 prediction
run_c1_prediction_layer <- function(layer = c("protein", "metabolite")) {
	if (LE8_REUSE_RESULTS) return(le8_restore_outputs(match.arg(layer), "c1_prediction"))
	# Check reusable results and initialize the analysis output directory.
	layer <- match.arg(layer)
	le8_begin_analysis(layer, "c1_prediction")
	.le8_analysis_env <- environment()
	on.exit(le8_finish_analysis(layer, "c1_prediction", .le8_analysis_env), add = TRUE)

	layer <- match.arg(layer) ; outdir <- if (layer == "protein") out.prot else out.met ; setwd2(outdir) ; rawdir <- le8_job_dir(outdir, LE8_JOB) ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE) ; setwd2(rawdir)
	module_file <- function(job) file.path(le8_job_dir(outdir, job), paste0(sub("_.*$", "", job), ".res.rds"))
	safe_module <- function(job, required = FALSE) {
		f <- module_file(job)
		if (job != "c1_correlate") return(list(data = list(), available = FALSE, status = "outside C1 baseline", file = f))
		if (!file.exists(f)) {
			if (required) stop("C1 prediction requires C1 for the same layer; missing: ", f, call. = FALSE)
			return(list(data = list(), available = FALSE, status = "missing", file = f))
		}
		z <- tryCatch(readRDS(f), error = function(e) e)
		if (inherits(z, "condition")) {
			if (required) stop("C1 prediction could not read required C1 result: ", conditionMessage(z), call. = FALSE)
			warning("C1 prediction: optional ", job, " result is unreadable; corresponding panels will be unavailable: ", conditionMessage(z), call. = FALSE)
			return(list(data = list(), available = FALSE, status = "unreadable", file = f))
		}
		le8_check_options(z)
		declared <- str_to_lower(as.character(z$meta$status %||% "ok")[[1]])
		usable <- !declared %in% c("unavailable", "no qtl locus constructed", "failed", "not run")
		if (!usable) return(list(data = list(), available = FALSE, status = paste0("available file; ", declared), file = f))
		list(data = z, available = TRUE, status = "available", file = f)
	}
	u1 <- safe_module("c1_correlate", TRUE) ; u2 <- safe_module("c2_cause") ; u3 <- safe_module("c3_coloc") ; u4 <- safe_module("c4_connect")
	upstream_signature <- paste0("c1=", u1$status, ";c2=", u2$status, ";c3=", u3$status, ";c4=", u4$status)
	upstream_audit <- tibble(
		module = c("C1", "C2", "C3", "C4"),
		required = c(TRUE, FALSE, FALSE, FALSE), available = c(u1$available, u2$available, u3$available, u4$available),
		status = c(u1$status, u2$status, u3$status, u4$status), file = c(u1$file, u2$file, u3$file, u4$file)
	)
	write_raw_csv(upstream_audit, "upstream_availability.csv", rawdir)
	message("C1 prediction/", layer, ": ", upstream_signature)
	final_cache <- file.path(rawdir, "res.rds")
	if (cache_valid(final_cache)) {
		old <- tryCatch(readRDS(final_cache), error = function(e) NULL)
		if (is.list(old) && all(c("meta", "scores", "prediction", "summary") %in% names(old))) {
			le8_check_final_request(old)
			cache_message(paste0("C1 prediction/", layer), final_cache) ; return(le8_restore_outputs(layer, "c1_prediction"))
		}
		message("C1 prediction/", layer, ": cache incomplete; recomputing score coverage and IPCW figures")
	}
	all_glmnet_label <- if (layer == "protein") "Pradeep-style / glmnet" else "All-metabolite / glmnet"
	lightgbm_label <- if (layer == "protein") "Yu-style / LightGBM" else "MWAS-ranked / LightGBM"
	biom <- if (layer == "protein") read_prot() else read_met() ; biom_vars <- setdiff(names(biom), "eid")
	incfile <- if (layer == "protein") FINAL_INC_PROT else FINAL_INC_MET
	if (nzchar(incfile) && file.exists(incfile)) {
		req <- unique(scan(incfile, what = "character", quiet = TRUE)) ; biom_vars <- intersect(biom_vars, req) ; biom <- biom[, c("eid", biom_vars), drop = FALSE]
	}
	all0 <- read_all() ; prs_vars <- character() # Disease PRS requires an explicit manifest in joint Final; no ambiguous name matching
	need <- unique(c("eid", "ethnic.c", le8_custom_covars, vars.basic, vars.le8, prs_vars, "birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_",Y), le8_y_date(Y), Sys.getenv("LE8_GROUP_COLUMN",Sys.getenv("PGS_GROUP_COLUMN",""))))
	dat_all <- all0[, intersect(need, names(all0)), drop = FALSE] |>
		filter_analysis_cohort() |>
		make_outcome(Y) |>
		add_attained_age_time(Y) |>
		inner_join(biom, by = "eid") ; rm(all0, biom) ; gc()
	tvar <- paste0(Y, ".t2e") ; evar <- paste0(Y, ".Yt2e") ; bvar <- paste0(Y, ".b2e") ; bivar <- paste0(Y, ".bi2e")
	clinical <- intersect(covs_use, names(dat_all)) ; biom_vars <- intersect(biom_vars, names(dat_all))
	prevalent_dat <- dat_all |> filter(is.finite(.data[[bvar]]), .data[[bvar]] <= 0)
	dat <- dat_all[complete.cases(dat_all[, c(tvar, evar), drop = FALSE]) & dat_all[[tvar]] > 0, , drop = FALSE] ; rm(dat_all) ; invisible(gc())
	split <- make_outer_split(dat, evar = evar) ; dat$final_split <- split
	# C1_LOCKPAPER::2026-10-17.c1-lockpaper-v1
	protocol <- c1_lockpaper_protocol()
	write_raw_csv(tibble::as_tibble(protocol), "c1.grouping_audit.csv", rawdir)
	write_raw_csv(tibble::tibble(protocol="2026-10-17.c1-lockpaper-v1",outer_roster=normalizePath(Sys.getenv("LE8_OUTER_ROSTER"),mustWork=TRUE),
		outer_roster_sha256=protocol$outer_roster_sha256,
		N=nrow(dat),training=sum(split=="training"),validation=sum(split=="validation"),events=sum(dat[[evar]]==1),
		training_eid_hash=c1_eid_hash(dat$eid[split=="training"]),
		validation_eid_hash=c1_eid_hash(dat$eid[split=="validation"])),
		"c1.lockpaper_protocol.csv",rawdir)
	message("C1 prediction/", layer, ": N=", nrow(dat), ", events=", sum(dat[[evar]] == 1), ", training=", sum(split == "training"), ", validation=", sum(split == "validation"), ", biomarkers=", length(biom_vars))

	c1f <- u1$file ; c1 <- u1$data ; a1 <- c1$association %||% tibble()
	if (!nrow(a1)) stop("C1 prediction requires a non-empty C1 association table for the same layer: ", c1f, call. = FALSE)
	# Re-screen C1 associations inside the outer training set.  The full-cohort
	# C1 table remains part of 5C evidence consolidation but is not allowed to
	# choose features for held-out prediction.
	screen_cache <- file.path(rawdir, "training_association_screen.rds") ; training_screen <- read_stage_cache(screen_cache)
	if (is.null(training_screen)) {
		message("C1 prediction/", layer, ": training-only PWAS/MWAS screen")
		training_screen <- cox_scan(dat[split == "training", , drop = FALSE], biom_vars, clinical, Y, time_var = tvar, event_var = evar)
		write_stage_cache(training_screen, screen_cache)
	} else message("C1 prediction/", layer, ": reuse training-only association screen")
	write_raw_csv(training_screen, "training_association_screen.csv", rawdir)
	ranked <- if (nrow(training_screen)) training_screen |>
		filter(is.finite(p.value)) |>
		arrange(p.value) |>
		pull(term) |>
		intersect(biom_vars) else biom_vars
	pwas <- if (nrow(training_screen)) training_screen |>
		filter(is.finite(FDR), FDR < .05) |>
		arrange(p.value) |>
		pull(term) |>
		intersect(biom_vars) else character()
	c2f <- u2$file ; c2 <- u2$data ; m2 <- c2$MR %||% tibble()
	mrset <- if (nrow(m2) && truthy(Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE"))) m2 |>
		filter(is.finite(FDR_all), FDR_all < le8_evidence_policy()$mr_fdr) |>
		arrange(pval) |>
		pull(exposure) |>
		unique() |>
		intersect(biom_vars) else character()
	local_class <- if (layer == "protein") "cis" else "local"
	local_mrset <- if (nrow(m2) && truthy(Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE"))) m2 |>
		filter(analysis == local_class, is.finite(FDR_all), FDR_all < le8_evidence_policy()$mr_fdr) |>
		arrange(pval) |>
		pull(exposure) |>
		unique() |>
		intersect(biom_vars) else character()
	# Training-only five-year landmark screen: this is the predictive component
	# least compatible with a purely imminent-diagnosis signal.
	distal_cache <- file.path(rawdir, paste0("distal_landmark_screen.", FINAL_CODE_VERSION, ".rds")) ; distal_screen <- read_stage_cache(distal_cache)
	if (is.null(distal_screen)) {
		distal_candidates <- head(ranked, FINAL_DISTAL_TOP) ; td <- dat[split == "training" & is.finite(dat[[tvar]]) & dat[[tvar]] > FINAL_DISTAL_LANDMARK, , drop = FALSE]
		if (nrow(td) >= 300 && sum(td[[evar]] == 1, na.rm = TRUE) >= 30 && length(distal_candidates)) {
			td$.distal_time <- td[[tvar]] - FINAL_DISTAL_LANDMARK ; td$.distal_event <- td[[evar]]
			distal_screen <- cox_scan(td, distal_candidates, clinical, Y, time_var = ".distal_time", event_var = ".distal_event") |> mutate(landmark_years = FINAL_DISTAL_LANDMARK)
		} else distal_screen <- tibble(term = character(), beta = numeric(), p.value = numeric(), FDR = numeric(), landmark_years = numeric())
		write_stage_cache(distal_screen, distal_cache)
	}
	write_raw_csv(distal_screen, "training_distal_landmark_screen.csv", rawdir)
	distalset <- distal_screen |>
		filter(is.finite(FDR), FDR < .05) |>
		arrange(p.value) |>
		pull(term) |>
		intersect(biom_vars)
	c3f <- u3$file ; c3 <- u3$data ; co3 <- as_tibble(c3$summary %||% tibble())
	if (nrow(co3) && !"PP.H4_robust_min" %in% names(co3)) co3$PP.H4_robust_min <- NA_real_
	if (nrow(co3) && !"status" %in% names(co3)) co3$status <- "ok"
	colocset <- if (nrow(co3) && truthy(Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE"))) co3 |>
		filter(le8_coloc_pass(pick(everything()))) |>
		pull(feature) |>
		unique() |>
		intersect(biom_vars) else character()
	# The causal score is deliberately stricter than a union of all MR and all
	# coloc hits: require local/cis MR and robust colocalization for the same
	# feature. Distal/trans-only MR remains in the separate MR score.
	geneticset <- if (u2$available && u3$available && truthy(Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE"))) le8_same_locus_candidates(m2, co3, layer) else character()
	hybridset <- if (u2$available && u3$available && length(geneticset) && length(distalset)) unique(c(distalset, geneticset)) else character()
	# A compact, transparent alternative to a black-box all-omic score. Candidate
	# selection uses only the outer training screen plus genetic/locus evidence.
	# Lower glmnet penalty factors favor distal and causal/locus-supported features;
	# training-ranked features lacking either support are retained but downweighted.
	mechanismset <- head(intersect(unique(c(geneticset, distalset, head(ranked, FINAL_MECHANISM_N))), biom_vars), FINAL_MECHANISM_N)
	mechanism_prior <- tibble(
		feature = mechanismset,
		distal_training = feature %in% distalset, genetic_locus = feature %in% geneticset
	) |>
		mutate(
			prior_class = case_when(
				genetic_locus & distal_training ~ "genetic + distal",
				genetic_locus ~ "genetic/locus", distal_training ~ "distal training signal", TRUE ~ "training-ranked only"
			),
			penalty_factor = 1.0,
			interpretation = "Uniform glmnet penalty; evidence defines candidates, not hard-coded biomarker weights"
		)
	mechanism_pf <- setNames(mechanism_prior$penalty_factor, mechanism_prior$feature)
	write_raw_csv(mechanism_prior, "mechanism_weight_prior.csv", rawdir)
	# For this analysis, "User specified" means the ten strongest training-screen biomarkers.
	user <- head(ranked, 10L)
	yuvars <- if (nzchar(FINAL_YU_LIST) && file.exists(FINAL_YU_LIST)) intersect(unique(scan(FINAL_YU_LIST, what = "character", quiet = TRUE)), biom_vars) else head(ranked, FINAL_YU_PRESELECT_N)
	if (!nzchar(FINAL_YU_LIST) && layer == "protein") message("C1 prediction: FINAL_YU_LIST not supplied; Yu column is a Yu-style 257-protein proxy using the C1 ranking, not an exact reproduction of the published 257-protein panel.")
	if (!nzchar(FINAL_YU_LIST) && layer == "metabolite") message("C1 prediction: metabolite LightGBM uses the C1 MWAS ranking (up to FINAL_YU_PRESELECT_N available metabolites).")

	# Nature Aging-style sequential screen, learned only inside the outer training sample.
	seq0 <- sequential_forward_panel(dat[split == "training", , drop = FALSE], pwas, tvar, evar, HORIZON, SEED + 71, max_extended = FINAL_EXTENDED_N)
	pars <- intersect(seq0$parsimonious, biom_vars) ; extended <- intersect(seq0$extended, biom_vars)
	write_raw_csv(seq0$log, "sequential_forward_training.csv", rawdir)

	# C4 connection-defined candidate sets for Fig3.
	c4f <- u4$file ; c4 <- u4$data ; l4 <- c4$lists %||% list()
	required_c4_lists <- c("NS", "YS_all", "YSP_all")
	c4_contract_ok <- u4$available && all(required_c4_lists %in% names(l4))
	if (u4$available && !c4_contract_ok) warning("C1 prediction: C4 result lacks the final list contract; C4 panels will be unavailable.", call. = FALSE)
	if (!c4_contract_ok) l4 <- list(NS = character(), YS_all = character(), YSP_all = character())
	ys4 <- if (c4_contract_ok) le8_training_connection_set(dat[split == "training", , drop = FALSE], biom_vars, intersect(vars.le8, names(dat)), intersect(vars.basic, names(dat)), rawdir) else character()
	# YS is selected without using disease outcome. Rebuild the outcome-ranked
	# NS and YSP-plus portions from the outer training screen to avoid leakage.
	ns4 <- if (c4_contract_ok) head(ranked, as.integer(le8_num_env("FINAL_CONNECTION_NS_N", 160))) else character()
	plus_n <- if (c4_contract_ok) as.integer(le8_num_env("FINAL_CONNECTION_PLUS_N", 80)) else 0L
	plus4 <- if (plus_n > 0L) head(setdiff(ranked, ys4), plus_n) else character()
	c4sets <- list(`C4 NS` = ns4, `C4 YS` = ys4, `C4 YSplus` = unique(c(ys4, plus4)))
	write_raw_csv(bind_rows(imap(c4sets, ~ tibble(set = .y, feature = .x))), "connection_sets_training.csv", rawdir)
	empty_c4 <- names(c4sets)[lengths(c4sets) == 0L]
	if (length(empty_c4)) message("C1 prediction: optional C4 score panels unavailable for: ", paste(empty_c4, collapse = ", "))

	module_tag <- paste0("c2", as.integer(u2$available), "c3", as.integer(u3$available), "c4", as.integer(c4_contract_ok))
	cache <- file.path(rawdir, paste0("score_models.", FINAL_MODEL_VERSION, ".", module_tag, ".rds")) ; sc <- if (cache_valid(cache)) tryCatch(readRDS(cache), error = function(e) NULL) else NULL
	if (is.null(sc)) {
		message("C1 prediction: fit ", all_glmnet_label, " all-biomarker score")
		pr <- fit_glmnet_score(dat, biom_vars, tvar, evar, split, all_glmnet_label, lambda_rule = "lambda.1se")
		message("C1 prediction: fit ", lightgbm_label, " on ", length(yuvars), " C1-ranked biomarkers")
		yu <- fit_lightgbm_score(dat, yuvars, tvar, evar, split, lightgbm_label)
		us <- fit_glmnet_score(dat, user, tvar, evar, split, "User specified", lambda_rule = "lambda.1se")
		pw <- fit_glmnet_score(dat, pwas, tvar, evar, split, "PWAS/MWAS significant", lambda_rule = "lambda.1se")
		mr <- fit_glmnet_score(dat, mrset, tvar, evar, split, "MR significant", lambda_rule = "lambda.1se")
		ds <- fit_glmnet_score(dat, distalset, tvar, evar, split, "Distal antecedent", lambda_rule = "lambda.1se")
		gc <- fit_glmnet_score(dat, geneticset, tvar, evar, split, "Genetic-region evidence", lambda_rule = "lambda.1se")
		hy <- fit_glmnet_score(dat, hybridset, tvar, evar, split, "Hybrid triangulated", lambda_rule = "lambda.1se")
		mw <- fit_glmnet_score(dat, mechanismset, tvar, evar, split, "Evidence-selected compact",
			lambda_rule = "lambda.1se", penalty_factor = mechanism_pf
		)
		pa <- fit_glmnet_score(dat, pars, tvar, evar, split, "Parsimonious", lambda_rule = "lambda.1se")
		ex <- fit_glmnet_score(dat, extended, tvar, evar, split, "Extended", lambda_rule = "lambda.1se")
		ns <- fit_glmnet_score(dat, c4sets[["C4 NS"]], tvar, evar, split, "C4 NS", lambda_rule = "lambda.1se")
		ys <- fit_glmnet_score(dat, c4sets[["C4 YS"]], tvar, evar, split, "C4 YS", lambda_rule = "lambda.1se")
		yp <- fit_glmnet_score(dat, c4sets[["C4 YSplus"]], tvar, evar, split, "C4 YSplus", lambda_rule = "lambda.1se")
		gp <- fit_glmnet_score(dat, prs_vars, tvar, evar, split, "Genetic PRS", lambda_rule = "lambda.min")
		clinical_obj <- fit_glmnet_score(dat, clinical, tvar, evar, split, "Clinical", lambda_rule = "lambda.min")
		objs <- list(
			Pradeep = pr, Yu = yu, User = us, PWAS = pw, MR = mr, Distal = ds, Genetic = gc, Hybrid = hy,
			Mechanism = mw, Parsimonious = pa, Extended = ex, C4_NS = ns, C4_YS = ys, C4_YSplus = yp, PRS = gp
		)
		score_rows <- standardize_scores(bind_rows(lapply(Filter(Negate(is.null), objs), `[[`, "rows")))
		sc <- c(objs, list(
			Clinical = clinical_obj, rows = score_rows, split = split,
			sets = list(
				PWAS = pwas, MR = mrset, Distal = distalset, Genetic = geneticset, Hybrid = hybridset,
				Mechanism = mechanismset, Parsimonious = pars, Extended = extended, Yu = yuvars, User = user, C4 = c4sets
			)
		))
		saveRDS(sc, cache, compress = "xz")
	} else {
		message("C1 prediction: reuse score-model cache") ; score_rows <- sc$rows ; clinical_obj <- sc$Clinical
	}

	object_map <- setNames(
		list(
			sc$Pradeep, sc$Yu, sc$User, sc$PWAS, sc$MR, sc$Distal, sc$Genetic, sc$Hybrid, sc$Mechanism, sc$Parsimonious, sc$Extended,
			sc$C4_NS, sc$C4_YS, sc$C4_YSplus, sc$PRS
		),
		c(
			all_glmnet_label, lightgbm_label, "User specified", "PWAS/MWAS significant", "MR significant", "Distal antecedent", "Genetic-region evidence", "Hybrid triangulated", "Evidence-selected compact", "Parsimonious", "Extended",
			"C4 NS", "C4 YS", "C4 YSplus", "Genetic PRS"
		)
	)
	base_score_rows <- score_rows |> filter(split != "prevalent")
	prev_rows <- imap_dfr(object_map, function(obj, nm) predict_prevalent_rows(obj, prevalent_dat, nm, bvar, base_score_rows))
	score_rows <- bind_rows(base_score_rows, prev_rows)
	score_rows$kind <- "Omic score"
	score_audit <- score_coverage_audit(score_rows, names(object_map), nrow(prevalent_dat))
	bad_prev <- score_audit |> filter(split == "prevalent", status != "ok")
	if (nrow(bad_prev)) warning("C1 prediction prevalent scoring coverage issue for: ", paste(paste0(bad_prev$method, " [", bad_prev$status, "]"), collapse = "; "), call. = FALSE)
	write_raw_csv(score_audit, "score_coverage_audit.csv", rawdir)
	write_raw_csv(score_rows, "person_scores.csv", rawdir)
	prior <- NULL
	pred <- prior$prediction %||% tibble()
	if (!nrow(pred)) {
		bundles <- imap(object_map, function(obj, nm) if (is.null(obj)) tibble() else make_validation_bundle(dat, obj, clinical_obj, clinical, tvar, evar, split, nm))
		pred <- bind_rows(bundles)
	}
	if (!nrow(pred)) pred <- empty_prediction_rows()
	write_raw_csv(pred, "prediction_validation_rows.csv", rawdir)
	c1_export_benchmark(pred, layer, rawdir, HORIZON)
	psum <- if (nrow(pred)) pred |>
		group_by(biom_set, model) |>
		summarise(
			N = n(), events = sum(event), C_index = first(cindex),
			AUC = weighted_time_auc(time, event, score, HORIZON), AUC_horizon_years = HORIZON,
			AUC_definition = "IPCW cumulative/dynamic", n_selected = median(n_selected, na.rm = TRUE), .groups = "drop"
		) else
		tibble(
			biom_set = character(), model = character(), N = integer(), events = integer(), C_index = numeric(),
			AUC = numeric(), AUC_horizon_years = numeric(), AUC_definition = character(), n_selected = numeric()
		)
	write_raw_csv(psum, "prediction_summary.csv", rawdir)

	method_names <- unique(score_rows$method)
	fig1_order <- c(all_glmnet_label, lightgbm_label, "User specified")
	fig2_order <- c("Distal antecedent", "Genetic-region evidence", "Hybrid triangulated", "Evidence-selected compact")
	fig3_order <- c("C4 NS", "C4 YS", "C4 YSplus")
	# Write the primary panels before optional matching/bootstrap diagnostics.
	# If a very large metabolite run is interrupted later, the principal C1 prediction
	# figures still exist and the output log identifies the unfinished step.
	save_plot(final_make_5row_grid(fig1_order, pred, score_rows, list(), list(), "Prediction paradigms"), "Fig1.simple_pred.png", 24, 13.5, outdir = outdir)
	save_plot(final_make_5row_grid(fig2_order, pred, score_rows, list(), list(), "Evidence-screened prediction"), "Fig2.screened_pred.png", 24, 13.5, outdir = outdir)
	save_plot(final_make_5row_grid(fig3_order, pred, score_rows, list(), list(), "C4 connection-guided prediction"), "Fig3.connection_pred.png", 24, 13.5, outdir = outdir)
	lead_methods <- intersect(c(all_glmnet_label, lightgbm_label, "Distal antecedent", "Genetic-region evidence", "Hybrid triangulated", "Evidence-selected compact", "C4 YSplus"), method_names)
	individual_features <- unique(c(head(ranked, 3L), intersect(FINAL_LEAD_ANCHORS, biom_vars))) |> head(6L)
	individual_lead_rows <- make_individual_lead_rows(dat, individual_features, training_screen, split, tvar, evar)
	lead_input <- bind_rows(score_rows |> filter(method %in% lead_methods), individual_lead_rows)
	lead <- leadtime_auc_table(lead_input, unique(lead_input$method), FINAL_TIME_CAP)
	lead_headline <- leadtime_headline(lead, FINAL_LEAD_MIN_CASES)
	lead_windows <- prediction_window_auc_table(lead_input, unique(lead_input$method))
	write_raw_csv(lead, "leadtime_discrimination.csv", rawdir) ; write_raw_csv(lead_headline, "leadtime_headline.csv", rawdir) ; write_raw_csv(lead_windows, "leadtime_prediction_windows.csv", rawdir)
	save_plot(plot_leadtime_prediction(lead, lead_windows), "Fig4.leadtime_prediction.png", 20, 11.25, outdir = outdir)
	topn <- topn_score_benchmark(dat, ranked, training_screen, tvar, evar, split, clinical, FINAL_TOPN_MAX, FINAL_TOPN_BOOT)
	write_raw_csv(topn$rows, "score_vs_topN_trajectory_rows.csv", rawdir)
	write_raw_csv(topn$summary, "score_vs_topN.csv", rawdir) ; write_raw_csv(topn$bootstrap, "score_vs_topN_bootstrap.csv", rawdir)
	write_raw_csv(topn$individuals, "score_vs_topN_individuals.csv", rawdir) ; write_raw_csv(topn$individual_range, "score_vs_topN_individual_range.csv", rawdir)
	save_plot(plot_topn_mechanism(topn), "Fig5.score_vs_topN_mechanism.png", 19, 10.7, outdir = outdir)
	attained_score <- attained_age_score_sensitivity(score_rows, dat, lead_methods, clinical, tvar, evar)
	write_raw_csv(attained_score, "attained_age_score_sensitivity.csv", rawdir)
	save_plot(plot_attained_age_score_sensitivity(attained_score), "Fig12.attained_age_sensitivity.png", 14, 8.5, outdir = outdir)
	# Matching/trajectory bootstrap is diagnostic and is not consumed by the
	# five-row grid.  Limit it to the main paradigm scores so the 88k-person
	# metabolite branch cannot stall before any figure is written.
	pair_methods <- if (nrow(dat) > 100000L) character() else intersect(fig1_order, method_names)
	if (!length(pair_methods) && nrow(dat) > 100000L) message("C1 prediction/", layer, ": skip optional matched-pair bootstrap for N > 100,000; lead-time and top-N uncertainty are still estimated")
	pairs_native <- prior$preclinical_pairs %||% setNames(vector("list", length(pair_methods)), pair_methods)
	boot_by_method <- prior$yy_auc_bootstrap %||% setNames(vector("list", length(pair_methods)), pair_methods)
	for (i in seq_along(pair_methods)) {
		nm <- pair_methods[[i]]
		if (is.null(pairs_native[[nm]]) || !nrow(pairs_native[[nm]])) pairs_native[[nm]] <- tryCatch(risk_set_pairs(score_rows, dat, nm, "score_z", SEED + 100 + i), error = function(e) {
			warning(nm, ": pair construction failed: ", conditionMessage(e), call. = FALSE) ; tibble()
		})
		if (FINAL_BOOT > 0 && (is.null(boot_by_method[[nm]]) || !nrow(boot_by_method[[nm]]))) boot_by_method[[nm]] <- tryCatch(bootstrap_auc_link(score_rows |> filter(method == nm), dat, FINAL_BOOT), error = function(e) {
			warning(nm, ": bootstrap diagnostic failed: ", conditionMessage(e), call. = FALSE) ; tibble()
		})
		if (is.null(boot_by_method[[nm]])) boot_by_method[[nm]] <- tibble()
	}
	write_raw_csv(bind_rows(imap(pairs_native, ~ .x |> mutate(method = .y))), "preclinical_pairs_all.csv", rawdir)
	write_raw_csv(bind_rows(imap(boot_by_method, ~ .x |> mutate(method = .y))), "yy_auc_bootstrap_all.csv", rawdir)

	f1 <- final_make_5row_grid(fig1_order, pred, score_rows, pairs_native, boot_by_method, "Prediction paradigms")
	f2 <- final_make_5row_grid(fig2_order, pred, score_rows, pairs_native, boot_by_method, "Evidence-screened prediction")
	f3 <- final_make_5row_grid(fig3_order, pred, score_rows, pairs_native, boot_by_method, "C4 connection-guided prediction")
	save_plot(f1, "Fig1.simple_pred.png", 24, 13.5, outdir = outdir)
	save_plot(f2, "Fig2.screened_pred.png", 24, 13.5, outdir = outdir)
	save_plot(f3, "Fig3.connection_pred.png", 24, 13.5, outdir = outdir)

	# Optional nested-CV consolidation remains available; it is not required for Figs 1–3.
	cv <- list(pred = tibble(), metrics = tibble(), summary = tibble(), selected = tibble(), stability = tibble())
	if (FINAL_NESTED_CV) {
		core <- unique(c(tvar, evar, vars.basic, vars.le8)) ; dn <- dat[complete.cases(dat[, intersect(core, names(dat)), drop = FALSE]), , drop = FALSE]
		cv <- nested_cv_omics(dn, biom_vars, intersect(vars.le8, names(dn)), intersect(vars.basic, names(dn)), character(), tvar, evar, OUTER, INNER, HORIZON, SEED)
		write_raw_csv(cv$pred, "prediction_rows_nested.csv", rawdir) ; write_raw_csv(cv$metrics, "metrics_by_fold.csv", rawdir) ; write_raw_csv(cv$summary, "metrics_summary.csv", rawdir) ; write_raw_csv(cv$selected, "selected_by_fold.csv", rawdir) ; write_raw_csv(cv$stability, "selection_stability.csv", rawdir)
	}
	ev <- make_evidence_table(layer, outdir, cv$stability) ; write_raw_csv(ev, "evidence_consolidation.csv", rawdir)
	set_sizes <- tibble(
		set = c(paste0(all_glmnet_label, " (all)"), paste0(lightgbm_label, " (preselected)"), "User", "PWAS/MWAS significant", "MR significant", paste0(local_class, " MR significant"), "Robust colocalization", "Distal antecedent", "Genetic-region evidence", "Hybrid triangulated", "Evidence-selected compact", "Parsimonious", "Extended", "C4 NS", "C4 YS", "C4 YSplus", "Genetic PRS"),
		N = c(length(biom_vars), length(yuvars), length(user), length(pwas), length(mrset), length(local_mrset), length(colocset), length(distalset), length(geneticset), length(hybridset), length(mechanismset), length(pars), length(extended), length(c4sets[[1]]), length(c4sets[[2]]), length(c4sets[[3]]), length(prs_vars))
	)
	candidate_sets <- bind_rows(
		tibble(set = "PWAS/MWAS significant", feature = pwas), tibble(set = "MR significant", feature = mrset),
		tibble(set = paste0(local_class, " MR significant"), feature = local_mrset), tibble(set = "Robust colocalization", feature = colocset),
		tibble(set = "Distal antecedent", feature = distalset), tibble(set = "Genetic-region evidence", feature = geneticset),
		tibble(set = "Hybrid triangulated", feature = hybridset), tibble(set = "Evidence-selected compact", feature = mechanismset),
		tibble(set = "Parsimonious", feature = pars),
		tibble(set = "Extended", feature = extended), bind_rows(imap(c4sets, ~ tibble(set = .y, feature = .x)))
	) |> distinct()
	write_raw_csv(set_sizes, "input_set_sizes.csv", rawdir)
	write_raw_csv(candidate_sets, "candidate_sets.csv", rawdir)
	corr <- score_correlation_table(score_rows) ; write_raw_csv(corr, "score_correlations.csv", rawdir)
	subgroup_methods <- intersect(c(all_glmnet_label, "Parsimonious", "C4 YSplus"), method_names)
	sg <- subgroup_auc_table(score_rows, dat, subgroup_methods) ; write_raw_csv(sg, "subgroup_auc.csv", rawdir)
	save_plot(plot_evidence_matrix(ev) | plot_causal_reactive_map(ev), "Fig6.evidence_matrix.png", 18, 10.5, outdir = outdir)
	save_plot(plot_performance_benchmark(psum), "Fig7.performance_benchmark.png", 13, 9, outdir = outdir)
	save_plot(plot_incremental_performance(psum), "Fig8.incremental_performance.png", 13, 9, outdir = outdir)
	save_plot(plot_complexity_performance(psum), "Fig9.parsimony_performance.png", 10.5, 8, outdir = outdir)
	save_plot(plot_score_correlation(corr), "Fig10.score_concordance.png", 11, 9.5, outdir = outdir)
	save_plot(plot_subgroup_auc(sg), "Fig11.subgroup_discrimination.png", 12, 10, outdir = outdir)
	write_xlsx2(list(
		upstream_availability = upstream_audit, prediction_summary = psum, input_set_sizes = set_sizes, candidate_sets = candidate_sets, mechanism_weight_prior = mechanism_prior, training_screen = training_screen,
		distal_landmark_screen = distal_screen, score_coverage_audit = score_audit, sequential_forward = seq0$log,
		score_method_summary = score_rows |> group_by(method, kind, split) |> summarise(rows = n(), finite_scores = sum(is.finite(score_z)), mean_score = mean(score_z, na.rm = TRUE), sd_score = sd(score_z, na.rm = TRUE), .groups = "drop"),
		preclinical_pairs = bind_rows(imap(pairs_native, ~ .x |> mutate(method = .y))), yy_auc_bootstrap = bind_rows(imap(boot_by_method, ~ .x |> mutate(method = .y))),
		leadtime_discrimination = lead, leadtime_headline = lead_headline, leadtime_windows = lead_windows,
		score_vs_topN = topn$summary, topN_individuals = topn$individuals, topN_individual_range = topn$individual_range,
		attained_age_sensitivity = attained_score, evidence_consolidation = ev, score_correlations = corr, subgroup_AUC = sg, nested_cv_summary = cv$summary
	), "prediction_panels.xlsx")
	out <- list(
		meta = module_meta(layer, extra = list(N = nrow(dat), events = sum(dat[[evar]] == 1), training = sum(split == "training"), validation = sum(split == "validation"), code_version = FINAL_CODE_VERSION, upstream_signature = upstream_signature)),
		upstream_availability = upstream_audit,
		scores = sc, prediction = pred, summary = psum, set_sizes = set_sizes, sequential = seq0, preclinical_pairs = pairs_native, yy_auc_bootstrap = boot_by_method,
		score_coverage_audit = score_audit, candidate_sets = candidate_sets, mechanism_weight_prior = mechanism_prior, distal_screen = distal_screen, leadtime = lead, leadtime_headline = lead_headline, leadtime_windows = lead_windows,
		individual_lead_rows = individual_lead_rows, score_vs_topN = topn, attained_age_sensitivity = attained_score,
		evidence = ev, score_correlations = corr, subgroup_AUC = sg, nested_cv = cv
	)
	out$training_request <- le8_final_request()
	out$mock <- le8_mock_final(out, outdir)
	saveRDS(out, final_cache, compress = "xz") ; finalize_outputs(LE8_JOB, outdir) ; out
}

if (prot_DO) {
	invisible(run_c1_prediction_layer("protein")) ; gc(full = TRUE)
}
if (met_DO) {
	invisible(run_c1_prediction_layer("metabolite")) ; gc(full = TRUE)
}
