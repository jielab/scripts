# Final reports, per-layer prediction, joint validation and pipeline orchestration.
# Each mode runs in a separate process and reads the same external analysis tree.
.final_mode <- Sys.getenv("LE8_FINAL_MODE", unset = if (Sys.getenv("LE8_JOB") == "final_prediction") "pipeline" else "report")
if (!.final_mode %in% c("report", "reference", "joint", "pipeline")) stop("Unknown LE8_FINAL_MODE: ", .final_mode)
if (.final_mode == "reference") {
	# Final: unified score phenotyping and prediction figures.

	suppressPackageStartupMessages({
		fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
		source(file.path(fdir, "0.common.R")) ;
		# final methods
		# final.R
		# Independent prediction functions extracted from the retained LE8 implementation.
		# Shared preprocessing/survival helpers are loaded by final.R via 0.common.R.
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


		# final.R
		# Independent prediction functions extracted from the retained LE8 implementation.
		# Shared preprocessing/survival helpers are loaded by final.R via 0.common.R.
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
				# Final clinical predictors are numeric scores/covariates. Construct the two
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
	LE8_JOB <- "final_prediction"
	dir.create(le8_final_dir(), recursive = TRUE, showWarnings = FALSE)
	FINAL_CODE_VERSION <- "2026-09-15.final-reference"
	FINAL_MODEL_VERSION <- "2026-09-15.final-reference"
	writeLines("Retained per-layer reference: training-only measured selection; different omic cohorts; not joint prot/met/PGS validation.", file.path(le8_final_dir(), "reference_scope.txt"))
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
		set.seed(seed)
		if (!is.null(evar) && evar %in% names(dat)) {
			strata <- split(seq_len(nrow(dat)), as.character(dat[[evar]]), drop = TRUE)
			val <- unlist(lapply(strata, function(ii) sample(ii, max(1L, round(length(ii) * test_frac)))), use.names = FALSE)
		} else val <- sample(seq_len(nrow(dat)), max(1L, round(nrow(dat) * test_frac)))
		ifelse(seq_len(nrow(dat)) %in% val, "validation", "training")
	}

	fit_glmnet_score <- function(dat, vars, tvar, evar, split, label, inner = INNER, alpha = 1,
				lambda_rule = c("lambda.1se", "lambda.min"), penalty_factor = NULL) {
		lambda_rule <- match.arg(lambda_rule)
		if (!requireNamespace("glmnet", quietly = TRUE)) stop("Final requires glmnet.", call. = FALSE)
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
						"Final: score time-scale sensitivity unavailable for ", mm, " / ", kind,
						"; complete-case N=", nrow(dd), ", events=", sum(dd[[evar]] == 1), " (requires N >= 100 and events >= 20)"
					)
					return(empty)
				}
				fit <- tryCatch(coxph(as.formula(paste0(surv, " ~ ", paste(bt(rhs), collapse = " + "))), dd, ties = "efron"),
					error = function(e) {
						message("Final: score time-scale sensitivity fit failed for ", mm, " / ", kind, ": ", conditionMessage(e)) ; NULL
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
		m0 <- c2$MR %||% tibble() ; m <- if (nrow(m0)) m0 |>
			filter(is.finite(pval)) |>
			group_by(exposure) |>
			slice_min(pval, n = 1, with_ties = FALSE) |>
			ungroup() |>
			transmute(feature = exposure, MR_analysis = analysis, MR_beta = b, MR_p = pval, MR_FDR = FDR_all) else tibble(feature = character())
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
		co <- if (nrow(co0) && all(c("feature", "PP.H4") %in% names(co0))) co0 |>
			filter(status == "ok") |>
			mutate(.robust = PP.H4_robust_min) |>
			group_by(feature) |>
			slice_max(.robust, n = 1, with_ties = FALSE) |>
			ungroup() |>
			select(feature, PP.H4, PP.H4_robust_min, everything()) else tibble(feature = character())
		cx0 <- c4$membership %||% tibble() ; cx <- if (nrow(cx0) && "feature" %in% names(cx0)) cx0 else tibble(feature = character())
		st <- if (nrow(stability) && all(c("feature", "selection_frequency") %in% names(stability))) stability |>
			group_by(feature) |>
			summarise(max_selection_frequency = max(selection_frequency, na.rm = TRUE), .groups = "drop") else tibble(feature = character(), max_selection_frequency = numeric())
		z <- Reduce(function(x, y) full_join(x, y, by = "feature"), list(a, g, di, m, rmr, dd, co, cx, st)) ; if (!nrow(z)) return(tibble())
		defaults <- list(
			obs_FDR = NA_real_, obs_p = NA_real_, MR_p = NA_real_, MR_FDR = NA_real_, DANDELION_sensitivity = FALSE,
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
				C2_MR = if (u2$available) as.integer(is.finite(MR_FDR) & MR_FDR < .05) else NA_integer_,
				C2_REVERSE_MR = if (u2$available) as.integer(is.finite(reverse_MR_FDR) & reverse_MR_FDR < .05) else NA_integer_,
				C2_DANDELION = if (u2$available) as.integer(coalesce(DANDELION_primary, FALSE)) else NA_integer_,
				C2_DANDELION_sensitivity = if (u2$available) as.integer(coalesce(DANDELION_sensitivity, FALSE)) else NA_integer_,
				C2 = if (u2$available) as.integer(coalesce(C2_MR, 0L) == 1L | coalesce(C2_DANDELION, 0L) == 1L) else NA_integer_,
				C3 = if (u3$available) as.integer(is.finite(PP.H4_robust_min) & PP.H4_robust_min >= .7) else NA_integer_,
				C4 = if (u4$available) as.integer(coalesce(strict_YS, FALSE)) else NA_integer_,
				C5 = as.integer(is.finite(max_selection_frequency) & max_selection_frequency >= .6),
				FINAL_available = is.finite(max_selection_frequency), reactive_compatible = coalesce(reactive_compatible, FALSE),
				distal_support = coalesce(distal_support, FALSE), temporal_warning = reactive_compatible,
				C2_status = u2$status, C3_status = u3$status, C4_status = u4$status,
				evidence_count = rowSums(cbind(C1, C2, C3, C4, ifelse(FINAL_available, C5, NA_integer_)), na.rm = TRUE),
				causal_locus_evidence = coalesce(C2_MR, 0L) == 1L | coalesce(C2_DANDELION, 0L) == 1L | coalesce(C3, 0L) == 1L,
				evidence_grade = case_when(
					coalesce(C2_MR, 0L) == 1L & coalesce(C3, 0L) == 1L & C1 == 1 ~ "A",
					causal_locus_evidence & (C1 == 1 | distal_support) ~ "B",
					evidence_count >= 2 ~ "C", TRUE ~ "D"
				),
				consolidation_role = case_when(
					evidence_grade %in% c("A", "B") & reactive_compatible ~ "Causal + reactive mixed biomarker",
					evidence_grade == "A" ~ "High-priority causal candidate",
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


	# 🚩 Main Final
	run_final_layer <- function(layer = c("protein", "metabolite")) {
		if (LE8_REUSE_RESULTS) return(le8_restore_outputs(match.arg(layer), "final_prediction"))
		# Check reusable results and initialize the analysis output directory.
		layer <- match.arg(layer)
		le8_begin_analysis(layer, "final_prediction")
		.le8_analysis_env <- environment()
		on.exit(le8_finish_analysis(layer, "final_prediction", .le8_analysis_env), add = TRUE)

		layer <- match.arg(layer) ; outdir <- if (layer == "protein") out.prot else out.met ; setwd2(outdir) ; rawdir <- le8_job_dir(outdir, LE8_JOB) ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE) ; setwd2(rawdir)
		module_file <- function(job) file.path(le8_job_dir(outdir, job), paste0(sub("_.*$", "", job), ".res.rds"))
		safe_module <- function(job, required = FALSE) {
			f <- module_file(job)
			if (!file.exists(f)) {
				if (required) stop("Final requires C1 for the same layer; missing: ", f, call. = FALSE)
				return(list(data = list(), available = FALSE, status = "missing", file = f))
			}
			z <- tryCatch(readRDS(f), error = function(e) e)
			if (inherits(z, "condition")) {
				if (required) stop("Final could not read required C1 result: ", conditionMessage(z), call. = FALSE)
				warning("Final: optional ", job, " result is unreadable; corresponding panels will be unavailable: ", conditionMessage(z), call. = FALSE)
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
		message("Final/", layer, ": ", upstream_signature)
		final_cache <- file.path(rawdir, "res.rds")
		if (cache_valid(final_cache)) {
			old <- tryCatch(readRDS(final_cache), error = function(e) NULL)
			if (is.list(old) && all(c("meta", "scores", "prediction", "summary") %in% names(old))) {
				le8_check_final_request(old)
				cache_message(paste0("Final/", layer), final_cache) ; return(le8_restore_outputs(layer, "final_prediction"))
			}
			message("Final/", layer, ": cache incomplete; recomputing score coverage and IPCW figures")
		}
		all_glmnet_label <- if (layer == "protein") "Pradeep-style / glmnet" else "All-metabolite / glmnet"
		lightgbm_label <- if (layer == "protein") "Yu-style / LightGBM" else "MWAS-ranked / LightGBM"
		biom <- if (layer == "protein") read_prot() else read_met() ; biom_vars <- setdiff(names(biom), "eid")
		incfile <- if (layer == "protein") FINAL_INC_PROT else FINAL_INC_MET
		if (nzchar(incfile) && file.exists(incfile)) {
			req <- unique(scan(incfile, what = "character", quiet = TRUE)) ; biom_vars <- intersect(biom_vars, req) ; biom <- biom[, c("eid", biom_vars), drop = FALSE]
		}
		all0 <- read_all() ; prs_vars <- character() # Disease PRS requires an explicit manifest in joint Final; no ambiguous name matching
		need <- unique(c("eid", "ethnic.c", le8_custom_covars, vars.basic, vars.le8, prs_vars, "birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_", Y)))
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
		message("Final/", layer, ": N=", nrow(dat), ", events=", sum(dat[[evar]] == 1), ", training=", sum(split == "training"), ", validation=", sum(split == "validation"), ", biomarkers=", length(biom_vars))

		c1f <- u1$file ; c1 <- u1$data ; a1 <- c1$association %||% tibble()
		if (!nrow(a1)) stop("Final requires a non-empty C1 association table for the same layer: ", c1f, call. = FALSE)
		# Re-screen C1 associations inside the outer training set.  The full-cohort
		# C1 table remains part of 5C evidence consolidation but is not allowed to
		# choose features for held-out prediction.
		screen_cache <- file.path(rawdir, "training_association_screen.rds") ; training_screen <- read_stage_cache(screen_cache)
		if (is.null(training_screen)) {
			message("Final/", layer, ": training-only PWAS/MWAS screen")
			training_screen <- cox_scan(dat[split == "training", , drop = FALSE], biom_vars, clinical, Y, time_var = tvar, event_var = evar)
			write_stage_cache(training_screen, screen_cache)
		} else message("Final/", layer, ": reuse training-only association screen")
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
			filter(is.finite(FDR_all), FDR_all < .05) |>
			arrange(pval) |>
			pull(exposure) |>
			unique() |>
			intersect(biom_vars) else character()
		local_class <- if (layer == "protein") "cis" else "local"
		local_mrset <- if (nrow(m2) && truthy(Sys.getenv("FINAL_GENETIC_EVIDENCE_INDEPENDENT", unset = "FALSE"))) m2 |>
			filter(analysis == local_class, is.finite(FDR_all), FDR_all < .05) |>
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
			filter(status == "ok", PP.H4_robust_min >= .7) |>
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
		if (!nzchar(FINAL_YU_LIST) && layer == "protein") message("Final: FINAL_YU_LIST not supplied; Yu column is a Yu-style 257-protein proxy using the C1 ranking, not an exact reproduction of the published 257-protein panel.")
		if (!nzchar(FINAL_YU_LIST) && layer == "metabolite") message("Final: metabolite LightGBM uses the C1 MWAS ranking (up to FINAL_YU_PRESELECT_N available metabolites).")

		# Nature Aging-style sequential screen, learned only inside the outer training sample.
		seq0 <- sequential_forward_panel(dat[split == "training", , drop = FALSE], pwas, tvar, evar, HORIZON, SEED + 71, max_extended = FINAL_EXTENDED_N)
		pars <- intersect(seq0$parsimonious, biom_vars) ; extended <- intersect(seq0$extended, biom_vars)
		write_raw_csv(seq0$log, "sequential_forward_training.csv", rawdir)

		# C4 connection-defined candidate sets for Fig3.
		c4f <- u4$file ; c4 <- u4$data ; l4 <- c4$lists %||% list()
		required_c4_lists <- c("NS", "YS_all", "YSP_all")
		c4_contract_ok <- u4$available && all(required_c4_lists %in% names(l4))
		if (u4$available && !c4_contract_ok) warning("Final: C4 result lacks the final list contract; C4 panels will be unavailable.", call. = FALSE)
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
		if (length(empty_c4)) message("Final: optional C4 score panels unavailable for: ", paste(empty_c4, collapse = ", "))

		module_tag <- paste0("c2", as.integer(u2$available), "c3", as.integer(u3$available), "c4", as.integer(c4_contract_ok))
		cache <- file.path(rawdir, paste0("score_models.", FINAL_MODEL_VERSION, ".", module_tag, ".rds")) ; sc <- if (cache_valid(cache)) tryCatch(readRDS(cache), error = function(e) NULL) else NULL
		if (is.null(sc)) {
			message("Final: fit ", all_glmnet_label, " all-biomarker score")
			pr <- fit_glmnet_score(dat, biom_vars, tvar, evar, split, all_glmnet_label, lambda_rule = "lambda.1se")
			message("Final: fit ", lightgbm_label, " on ", length(yuvars), " C1-ranked biomarkers")
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
			message("Final: reuse score-model cache") ; score_rows <- sc$rows ; clinical_obj <- sc$Clinical
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
		if (nrow(bad_prev)) warning("Final prevalent scoring coverage issue for: ", paste(paste0(bad_prev$method, " [", bad_prev$status, "]"), collapse = "; "), call. = FALSE)
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
		# If a very large metabolite run is interrupted later, the principal C5
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
		if (!length(pair_methods) && nrow(dat) > 100000L) message("Final/", layer, ": skip optional matched-pair bootstrap for N > 100,000; lead-time and top-N uncertainty are still estimated")
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
		invisible(run_final_layer("protein")) ; gc(full = TRUE)
	}
	if (met_DO) {
		invisible(run_final_layer("metabolite")) ; gc(full = TRUE)
	}
} else if (.final_mode == "joint") {
	# Final: one joint incident cohort, one validation scheme, any prespecified Y.
	fdir <- Sys.getenv("LE8_FDIR", file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(fdir, "0.common.R"))
	final.out <- le8_final_dir(kind = "joint")
	dir.create(final.out, recursive = TRUE, showWarnings = FALSE)
	# final joint helpers final.R Common-cohort prediction for any Y. No C1/C2/C3 selection enters
	# validation.
	final_csv <- function(name, default) Filter(nzchar, trimws(strsplit(Sys.getenv(name, default), ",", fixed = TRUE)[[1]]))
	final_assert_ids <- function(x) {
		if (!"eid" %in% names(x) || anyNA(x$eid) || anyDuplicated(x$eid))
			stop("Missing/duplicate participant IDs")
		x$eid <- as.character(x$eid)
		x
	}
	final_group_folds <- function(group, k, seed) {
		if (anyNA(group) || any(!nzchar(as.character(group))))
			stop("Missing split group")
		g <- sort(unique(as.character(group)))
		if (length(g) < k)
			stop("Too few independent groups")
		set.seed(seed)
		f <- sample(rep(seq_len(k), length.out = length(g)))
		f[match(group, g)]
	}
	final_landmark <- function(d, L, end) {
		covered <- if (".entry" %in% names(d))
			is.finite(d$.entry) & d$.entry <= L else rep(TRUE, nrow(d))
		d <- d[covered & is.finite(d$.time) & d$.time > L & d$.event %in% 0 : 1, , drop = FALSE]
		# Preserve held-out follow-up beyond the evaluation horizon: dynamic AUC controls require time > horizon,
		# not administratively truncated equality.
		d$.time <- d$.time - L
		d
	}
	final_prs_manifest_input <- function(Y) {
		f <- Sys.getenv("FINAL_PRS_MANIFEST", "")
		if (!nzchar(f))
			stop("Final requires FINAL_PRS_MANIFEST: explicit disease PRS for each Y; omic PGS is not a disease PRS")
		m <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
		required <- c("Y", "file", "column", "source", "build", "ancestry", "sample_overlap")
		if (!all(required %in% names(m)))
			stop("PRS manifest requires: ", paste(required, collapse = ","))
		m <- m[m$Y == Y, , drop = FALSE]
		if (nrow(m) != 1 || anyNA(m[, required]) || any(!nzchar(as.matrix(m[, required]))))
			stop("Exactly one complete PRS manifest row is required for ", Y)
		if (m$sample_overlap != "none")
			stop("Use an external disease PRS whose discovery excludes this cohort. Internally derived PRS requires end-to-end outer-fold GWAS/scoring and is not supported by a single precomputed score column.")
		if (!grepl("^(/|[A-Za-z]:)", m$file))
			m$file <- file.path(dirname(f), m$file)
		d <- if (grepl("[.]rds$", m$file, ignore.case = TRUE))
			readRDS(m$file) else data.table::fread(m$file, data.table = FALSE)
		d <- final_assert_ids(as.data.frame(d))
		if (!m$column %in% names(d) || !is.numeric(d[[m$column]]))
			stop("PRS column missing or nonnumeric")
		d <- data.frame(eid = d$eid, disease_PRS = d[[m$column]])
		list(data = d[is.finite(d$disease_PRS), , drop = FALSE], manifest = m)
	}
	# Exact endpoint naming, never a fuzzy match to biomarker PGS or another trait.  An explicit manifest takes
	# precedence. Unknown discovery overlap is labelled, not silently treated as independent; no participant's
	# missing PRS is imputed.
	final_prs_input <- function(Y, ph) {
		column <- paste0(Y, ".pgs")
		unavailable <- function(status, m = NULL) {
			if (is.null(m))
				m <- data.frame(
					Y = Y, file = "all.rds", column = column, source = "not supplied", build = "unknown",
					ancestry = "unknown", sample_overlap = "unknown"
				)
			m$status <- status
			m$enabled <- FALSE
			list(data = data.frame(eid = as.character(ph$eid)), manifest = m, enabled = FALSE)
		}
		if (nzchar(Sys.getenv("FINAL_PRS_MANIFEST", ""))) {
			obj <- final_prs_manifest_input(Y)
			obj$manifest$status <- "explicit manifest; discovery independence declared by user"
		} else {
			if (!column %in% names(ph))
				return(unavailable(paste("absent:", column, "; non-PRS comparisons continue")))
			if (!is.numeric(ph[[column]]))
				return(unavailable(paste("unusable nonnumeric column:", column)))
			overlap <- Sys.getenv("FINAL_DISEASE_PRS_OVERLAP", "unknown")
			if (!overlap %in% c("none", "unknown", "yes"))
				stop("FINAL_DISEASE_PRS_OVERLAP must be none, unknown or yes")
			m <- data.frame(Y = Y, file = "all.rds", column = column, source = Sys.getenv(
				"FINAL_DISEASE_PRS_SOURCE",
				"all.rds exact endpoint column; provenance not supplied"
			), build = Sys.getenv(
				"FINAL_DISEASE_PRS_BUILD",
				"unknown"
			), ancestry = Sys.getenv("FINAL_DISEASE_PRS_ANCESTRY", "unknown"), sample_overlap = overlap)
			if (overlap == "yes")
				return(unavailable("known discovery overlap; precomputed PRS excluded from validation", m))
			g <- ph[[column]]
			obj <- list(
				data = data.frame(eid = as.character(ph$eid[is.finite(g)]), disease_PRS = g[is.finite(g)]),
				manifest = m
			)
			obj$manifest$status <- if (overlap == "none")
				"exact all.rds endpoint PRS; independence declared by user" else "exploratory endpoint PRS; discovery overlap unknown; outer CV does not resolve PRS discovery leakage"
		}
		z <- obj$data$disease_PRS
		if (length(z) < 2L || !is.finite(sd(z)) || sd(z) <= 0)
			return(unavailable("PRS has insufficient finite values or zero variance", obj$manifest))
		obj$enabled <- TRUE
		obj$manifest$enabled <- TRUE
		obj$manifest$N_available <- length(z)
		obj
	}
	final_join_omics <- function(ph, prot, met, prs) {
		ph <- final_assert_ids(ph)
		map <- list()
		for (layer in c("prot", "met")) {
			b <- final_assert_ids(if (layer == "prot")
				prot else met)
			ff <- setdiff(names(b), "eid")
			if (any(grepl("([.]pgs$|[.]prs$|^fod_|^date_|[.]t2e$|[.]Yt2e$)", ff, ignore.case = TRUE)))
				stop("Non-assay columns in ", layer, " matrix")
			if (!all(vapply(b[ff], is.numeric, logical(1))))
				stop("Omics matrices require numeric assay columns only")
			map[[layer]] <- data.frame(feature = paste0(layer, "__", ff), assay = ff, layer = layer)
			names(b)[match(ff, names(b))] <- map[[layer]]$feature
			ph <- merge(ph, b, by = "eid", sort = FALSE)
		}
		ph <- merge(ph, final_assert_ids(prs), by = "eid", sort = FALSE)
		list(data = ph[order(ph$eid), , drop = FALSE], map = bind_rows(map))
	}
	final_attach_omic_pgs <- function(dat, map) {
		mf <- Sys.getenv("FINAL_OMIC_PGS_MANIFEST", "")
		if (!nzchar(mf))
			return(list(data = dat, map = map, inputs = character(), status = "not supplied; disease PRS and biomarker PGS are distinct"))
		m <- read.csv(mf, stringsAsFactors = FALSE)
		req <- c("layer", "file", "map_file", "source", "sample_overlap")
		if (!all(req %in% names(m)) || anyNA(m[, req]) || anyDuplicated(m$layer) || !setequal(m$layer, c("prot", "met")) ||
			any(m$sample_overlap != "none"))
			stop("Omic-PGS manifest requires prot/met rows and independent discovery (sample_overlap=none)")
		inputs <- mf
		map$pgs_feature <- NA_character_
		for (i in seq_len(nrow(m))) {
			path <- function(p) if (grepl("^(/|[A-Za-z]:)", p))
				p else file.path(dirname(mf), p)
			f <- path(m$file[i])
			mp <- path(m$map_file[i])
			inputs <- c(inputs, f, mp)
			lookup <- read.csv(mp, stringsAsFactors = FALSE)
			if (!all(c("assay", "column") %in% names(lookup)) || anyDuplicated(lookup$assay) || anyDuplicated(lookup$column))
				stop("Omic-PGS map needs one unique column per assay")
			g <- final_assert_ids(as.data.frame(if (grepl("[.]rds$", f, ignore.case = TRUE))
				readRDS(f) else data.table::fread(f, data.table = FALSE)))
			hit <- which(map$layer == m$layer[i] & map$assay %in% lookup$assay)
			cols <- lookup$column[match(map$assay[hit], lookup$assay)]
			if (!all(cols %in% names(g)) || !all(vapply(g[cols], is.numeric, logical(1))))
				stop("PGS mapped columns missing or nonnumeric")
			names_new <- paste0("pgs__", map$feature[hit])
			map$pgs_feature[hit] <- names_new
			g <- g[, c("eid", cols), drop = FALSE]
			names(g) <- c("eid", names_new)
			dat <- merge(dat, g, by = "eid", sort = FALSE)
		}
		list(data = dat[order(dat$eid), , drop = FALSE], map = map, inputs = inputs, status = "independent matched biomarker PGS available")
	}
	# Outcome-free, covariate-adjusted domain ranking. This is association strength, not a mediation estimate. A
	# feature may be relevant to both domains.
	final_domain_rank <- function(train, features, target, basic) {
		d <- as.data.frame(train)
		ok <- is.finite(d[[target]])
		if (sum(ok) < 100)
			stop("Insufficient observed training target: ", target)
		d <- d[ok, , drop = FALSE]
		z <- le8_prepare_prediction_matrix(d, d, basic)
		q <- qr(cbind(1, z$train))
		y <- qr.resid(q, d[[target]])
		sy <- sqrt(sum(y * y))
		if (!is.finite(sy) || sy <= 0)
			stop("Invariant domain target")
		blocks <- split(features, ceiling(seq_along(features) / 64))
		out <- map_dfr(blocks, function(bb) {
			X <- as.matrix(d[, bb, drop = FALSE])
			storage.mode(X) <- "double"
			observed <- colSums(is.finite(X))
			for (j in seq_len(ncol(X))) {
				v <- X[, j]
				med <- median(v[is.finite(v)], na.rm = TRUE)
				if (!is.finite(med))
					med <- 0
				v[!is.finite(v)] <- med
				X[, j] <- v
			}
			X <- qr.resid(q, X)
			den <- sqrt(colSums(X * X)) * sy
			r <- as.numeric(crossprod(X, y)) / den
			r[observed < 100 | !is.finite(r)] <- NA_real_
			tibble(feature = bb, domain = target, r = r, N_observed = observed)
		})
		out |>
			filter(is.finite(r)) |>
			arrange(desc(abs(r)), feature)
	}
	final_balance <- function(queues, k) {
		ans <- character()
		while (length(ans) < k && any(lengths(queues) > 0)) for (j in seq_along(queues)) {
			queues[[j]] <- setdiff(queues[[j]], ans)
			if (length(queues[[j]]) && length(ans) < k) {
				ans <- c(ans, queues[[j]][1])
				queues[[j]] <- queues[[j]][ - 1]
			}
		}
		ans
	}
	final_ys_panel <- function(queues, k, ns = NULL) {
		# Match the NS allocation as well as total assay budget.
		panel <- character()
		for (layer in c("prot", "met")) {
			kl <- if (layer == "prot")
				ceiling(k / 2) else floor(k / 2)
			q <- lapply(queues, function(v) v[startsWith(v, paste0(layer, "__"))])
			required <- if (is.null(ns))
				kl else max(1, floor(0.8 * kl))
			keep <- final_balance(q, required)
			if (length(keep) < required)
				return(character())
			if (!is.null(ns))
				keep <- head(c(keep, setdiff(ns[startsWith(ns, paste0(layer, "__"))], keep)), kl)
			panel <- c(panel, keep)
		}
		panel
	}
	final_gaussian_proxy <- function(train, test, target, features, seed) {
		ok <- is.finite(train[[target]])
		tr <- train[ok, , drop = FALSE]
		if (sum(ok) < 100)
			stop("Too few observed proxy targets")
		X <- le8_prepare_prediction_matrix(tr, test, features)
		if (!ncol(X$train))
			stop("No usable proxy assays")
		if (ncol(X$train) == 1) {
			fit <- lm.fit(cbind(1, X$train), tr[[target]])
			pred <- as.numeric(cbind(1, X$test) %*% fit$coefficients)
			b <- fit$coefficients
		} else {
			fd <- final_group_folds(tr$.group, 5, seed)
			fit <- glmnet::cv.glmnet(X$train, tr[[target]],
				family = "gaussian", alpha = 0, foldid = fd, standardize = TRUE,
				type.measure = "mse"
			)
			pred <- as.numeric(predict(fit, X$test, s = "lambda.1se"))
			b <- as.numeric(coef(fit, s = "lambda.1se"))
		}
		if (any(!is.finite(pred)))
			stop("Non-finite proxy prediction; no in-sample fallback")
		list(pred = pred, features = features, preprocess = X$audit, coefficient = data.frame(variable = c(
			"(Intercept)",
			colnames(X$train)
		), beta = b))
	}
	# Every training proxy prediction excludes that participant's entire group from ranking, proxy fitting and
	# lambda tuning. Yang never includes test groups.
	final_crossfit_proxy <- function(train, test, yang, features, target, k, basic, seed) {
		fd <- final_group_folds(train$.group, 5, seed)
		oof <- rep(NA_real_, nrow(train))
		fitone <- function(a, b, extra, s) {
			learn <- bind_rows(a, extra)
			rank <- final_domain_rank(learn, features, target, basic)
			panel <- head(rank$feature, k)
			if (length(panel) != k)
				stop("Insufficient eligible proxy panel")
			final_gaussian_proxy(learn, b, target, panel, s)
		}
		for (j in seq_len(5)) {
			hold <- train$.group[fd == j]
			extra <- yang[!yang$.group %in% hold, , drop = FALSE]
			oof[fd == j] <- fitone(train[fd != j, , drop = FALSE], train[fd == j, , drop = FALSE], extra, seed + j)$pred
		}
		obj <- fitone(train, test, yang, seed + 20)
		if (any(!is.finite(oof)))
			stop("Incomplete cross-fitted proxy; cannot fit disease model")
		list(train = oof, test = obj$pred, final = obj)
	}
	final_inflammation_score <- function(train, test, map) {
		a <- toupper(map$assay)
		hit <- map$feature[map$layer == "prot" & a %in% c("GDF15", "IL6", "IL-6")]
		if (length(hit) < 2)
			return(NULL)
		X <- le8_prepare_prediction_matrix(train, test, hit)
		if (ncol(X$train) != length(hit))
			return(NULL)
		ztr <- rowMeans(X$train)
		zte <- rowMeans(X$test)
		threshold <- unname(quantile(ztr, 0.75))
		list(
			group = ifelse(zte >= threshold, "higher inflammatory-marker burden", "lower inflammatory-marker burden"),
			threshold = threshold
		)
	}
	final_paired_delta <- function(a, b, horizon, B, seed) {
		z <- merge(a, b, by = c("eid", "fold", "time", "event"), suffixes = c(".a", ".b"))
		if (nrow(z) != nrow(a) || nrow(z) != nrow(b))
			return(tibble(status = "incomplete common predictions"))
		stat <- function(ix) {
			w <- le8_ipcw(z$time[ix], z$event[ix], horizon)
			if (w$status != "ok")
				return(c(AUC = NA_real_, Brier = NA_real_))
			c(AUC = le8_weighted_auc(z$risk.a[ix], w$y, w$w) - le8_weighted_auc(z$risk.b[ix], w$y, w$w), Brier = mean(w$w *
				((w$y - z$risk.a[ix]) ^ 2 - (w$y - z$risk.b[ix]) ^ 2)))
		}
		est <- stat(seq_len(nrow(z)))
		set.seed(seed)
		# Resample independent groups, retaining repeated relatives together.
		groups <- split(seq_len(nrow(z)), z$group.a)
		ng <- length(groups)
		bs <- if (B > 0)
			replicate(B, stat(unlist(groups[sample.int(ng, ng, replace = TRUE)], use.names = FALSE))) else matrix(NA_real_, 2, 0)
		ci <- function(v, p) {
			v <- v[is.finite(v)]
			if (length(v) >= max(20, 0.8 * B))
				unname(quantile(v, p)) else NA_real_
		}
		tibble(status = if (all(is.finite(est)))
			"ok" else "insufficient supported outcome follow-up", delta_AUC = est[1], AUC_lo = ci(bs[1, ], 0.025), AUC_hi = ci(bs[1, ], 0.975), delta_Brier = est[2], Brier_lo = ci(bs[2, ], 0.025), Brier_hi = ci(bs[2, ], 0.975), uncertainty = "paired group bootstrap of frozen out-of-fold predictions; training uncertainty excluded")
	}


	# final.R
	final_joint_outputs <- function(
		results, yin, map, root, K, B, seed, primary_budget, primary_landmark, signature,
		cached_validation = NULL
	) {
		private <- file.path(root, "_private")
		final_cleanup_stale_summaries(private)
		directory <- final_summary_directory(private)
		on.exit(unlink(directory, recursive = TRUE), add = TRUE)
		store <- final_stage_joint_results(results, signature, directory)
		if (is.null(yin)) {
			if (is.null(cached_validation))
				stop("Cached Final summaries require saved cohort validation")
			yin <- do.call(final_validate_cached_cohorts, c(list(cohorts = store$cohorts), cached_validation))
		}
		combine <- function(nm) final_store_table(store, nm)
		members <- combine("members") |>
			filter(feature %in% map$feature) |>
			distinct()
		diag <- combine("diagnostics")
		keys <- c("model", "budget", "arm", "landmark", "horizon")
		metrics <- cal <- dec <- foldmetrics <- pairs <- strata <- burden <- heterogeneity <- list()
		marker_groups <- combine("marker_groups")
		design <- read.csv(file.path(root, "factorial_design.csv"), stringsAsFactors = FALSE)
		factorial <- genetic_strata <- list()
		blocks <- unique(store$blocks[c("landmark", "arm")])
		prediction_export <- file.path(directory, "out_of_fold_predictions.csv.gz")
		for (block in seq_len(nrow(blocks))) {
			L0 <- blocks$landmark[block]
			a0 <- blocks$arm[block]
			pred <- final_store_predictions(store, L0, a0)
			message(
				"Final summary: block ", block, "/", nrow(blocks), "; landmark=", L0, "; arm=", a0, "; prediction rows=",
				nrow(pred)
			)
			groups <- split(seq_len(nrow(pred)), do.call(interaction, c(pred[keys], list(drop = TRUE, lex.order = TRUE))))
			for (ii in groups) {
				d <- pred[ii, , drop = FALSE]
				id <- d[1, keys, drop = FALSE]
				expected <- length(final_output_ids(yin, id$landmark))
				if (nrow(d) != expected || anyDuplicated(d$eid) || length(unique(d$fold)) != K) {
					metrics[[length(metrics) + 1L]] <- cbind(id,
						status = "incomplete outer predictions; do not compare",
						N = nrow(d)
					)
					next
				}
				z <- le8_evaluate_risk(d$time, d$event, d$risk, id$horizon, id$model, id$budget, "nested outer validation",
					B = 0
				)
				z$metrics$arm <- id$arm
				z$metrics$landmark <- id$landmark
				metrics[[length(metrics) + 1L]] <- z$metrics
				cal[[length(cal) + 1L]] <- z$calibration |>
					mutate(arm = id$arm, landmark = id$landmark)
				dec[[length(dec) + 1L]] <- z$decision |>
					mutate(arm = id$arm, landmark = id$landmark, budget = id$budget)
				for (fd in seq_len(K)) {
					q <- d[d$fold == fd, , drop = FALSE]
					q$.c_time <- pmin(q$time, id$horizon)
					q$.c_event <- as.integer(q$event == 1 & q$time <= id$horizon)
					cc <- tryCatch(survival::concordance(survival::Surv(.c_time, .c_event) ~ lp, data = q, reverse = TRUE),
						error = function(e) NULL
					)
					foldmetrics[[length(foldmetrics) + 1L]] <- cbind(id,
						fold = fd, N = nrow(q), events = sum(q$event),
						Harrell_C = if (is.null(cc))
							NA_real_ else cc$concordance, variance = if (is.null(cc))
							NA_real_ else cc$var, scope = "within-fold Harrell C restricted to the stated prediction horizon"
					)
				}
				if (id$budget == primary_budget && id$model %in% c(
					"Clinical_ProtMet_NS", "Clinical_ProtMet_YS_YinYang",
					"Clinical_ProtMet_PRS_NS"
				)) {
					for (bg in unique(d$inflammatory_burden)) {
						q <- d[d$inflammatory_burden == bg, , drop = FALSE]
						zz <- le8_evaluate_risk(q$time, q$event, q$risk, id$horizon, id$model, id$budget, "exploratory marker-burden subgroup",
							B = 0
						)
						burden[[length(burden) + 1L]] <- zz$metrics |>
							mutate(burden = bg, arm = id$arm, landmark = id$landmark)
					}
				}
			}
			# Prespecified contrasts. No selection of the best validation budget/model.
			contrasts <- data.frame(model = c(
				"Clinical_ProtMet_PRS_NS", "Clinical_ProtMet_PRS_NS", "Clinical_ProtMet_YS_YinYang",
				"Clinical_ProtMet_YS_YinYang", "Clinical_ProtMet_YSplus_YinYang", "Clinical_ProtMet_PRS_YS_YinYang",
				"Clinical_AddProxy_joint_YinYang", "Clinical_ReplaceProxy_joint_YinYang"
			), reference = c(
				"Clinical_ProtMet_NS",
				"Clinical_PRS", "Clinical_ProtMet_NS", "Clinical_ProtMet_YS_Yin", "Clinical_ProtMet_NS", "Clinical_ProtMet_PRS_NS",
				"Clinical", "Clinical"
			))
			contrasts <- bind_rows(contrasts, data.frame(
				model = c(
					"Clinical_MatchedMeasured_PGS_NS", "Clinical_MatchedMeasured_PGS_NS",
					"Clinical_MatchedMeasured_PGS_PRS_NS", "Clinical_ProtMet_NS", "Clinical_ProtMet_NS", "Clinical_ProtMet_YSplus_YinYang"
				),
				reference = c(
					"Clinical_MatchedMeasured_NS", "Clinical_MatchedPGS_NS", "Clinical_MatchedMeasured_PGS_NS",
					"Clinical_Protein_sharedBudget_NS", "Clinical_Metabolite_sharedBudget_NS", "Clinical_ProtMet_YSplus_Yin"
				)
			))
			contrasts <- bind_rows(contrasts, data.frame(model = c(
				"Clinical_ProtMet_PGSselected_NS", "Clinical_ProtMet_PGSselected_NS",
				"Clinical_ProtMet_PGSselected_PRS_NS", "Clinical_ProtMet_YS_PGSselected_YinYang"
			), reference = c(
				"Clinical_ProtMet_NS",
				"Clinical_PGSselected_NS", "Clinical_ProtMet_PGSselected_NS", "Clinical_ProtMet_PGSselected_NS"
			)))
			for (L in sort(unique(pred$landmark))) for (a in unique(pred$arm)) for (k in sort(unique(pred$budget[pred$budget >
				0]))) for (i in seq_len(nrow(contrasts))) {
				co <- contrasts[i, ]
				x <- pred[pred$model == co$model & pred$landmark == L & pred$arm == a & pred$budget == k, , drop = FALSE]
				kb <- if (co$reference %in% c("Clinical", "Clinical_PRS"))
					0 else k
				y <- pred[pred$model == co$reference & pred$landmark == L & pred$arm == a & pred$budget == kb, , drop = FALSE]
				expected <- length(final_output_ids(yin, L))
				if (nrow(x) != expected || nrow(y) != expected || anyDuplicated(x$eid) || anyDuplicated(y$eid))
					next
				h <- x$horizon[1]
				delta <- final_paired_delta(x, y, h, B, seed + as.integer(L * 100) + k)
				if (k == primary_budget && co$reference == "Clinical_ProtMet_NS" && co$model %in% c(
					"Clinical_ProtMet_YS_YinYang",
					"Clinical_ProtMet_YSplus_YinYang"
				)) {
					mg <- marker_groups[marker_groups$landmark == L, , drop = FALSE]
					mg <- mg[, setdiff(names(mg), c("fold", "landmark")), drop = FALSE]
					he <- final_marker_heterogeneity(x, y, mg, h, B, seed + as.integer(L * 100) + k)
					if (nrow(he))
						heterogeneity[[length(heterogeneity) + 1L]] <- he |>
							mutate(model = co$model, reference = co$reference, landmark = L, arm = a, budget = k)
				}
				pairs[[length(pairs) + 1L]] <- delta |>
					mutate(model = co$model, reference = co$reference, budget = k, landmark = L, arm = a, primary = (k ==
						primary_budget & L == primary_landmark & a == "all_assays" & co$model == "Clinical_ProtMet_YS_YinYang" &
						co$reference == "Clinical_ProtMet_NS"))
			}
			# Frozen training 75th-percentile thresholds in each fold. Show observed risk in the four PRS x omic
			# strata; no cutpoints optimized on held-out outcomes.
			for (L in sort(unique(pred$landmark))) for (a in unique(pred$arm)) {
				p <- pred[pred$model == "ProtMet_only_NS" & pred$budget == primary_budget & pred$landmark == L & pred$arm ==
					a, , drop = FALSE]
				g <- pred[pred$model == "PRS_only" & pred$landmark == L & pred$arm == a, , drop = FALSE]
				if (!nrow(p) || !nrow(g))
					next
				d <- merge(p, g[, c("eid", "high_training_q75")], by = "eid", suffixes = c(".omics", ".PRS"))
				d$stratum <- paste0("omics ", ifelse(d$high_training_q75.omics, "high", "lower"), " / PRS ", ifelse(d$high_training_q75.PRS,
					"high", "lower"
				))
				iw <- le8_ipcw(d$time, d$event, d$horizon[1])
				if (iw$status != "ok")
					next
				d$w <- iw$w
				d$y <- iw$y
				strata[[length(strata) + 1L]] <- d |>
					group_by(stratum) |>
					summarise(N = n(), events = sum(y), observed_risk_ipcw = sum(w * y) / sum(w), .groups = "drop") |>
					mutate(landmark = L, arm = a, interpretation = "Cause-specific net-risk display; high defined from training 75th percentile, not validated clinical cutoffs")
			}
			genetic_strata[[block]] <- final_genetic_measured_strata(pred, primary_budget)
			factorial[[block]] <- final_factorial_outputs(pred, yin, design, K, B, seed, primary_budget, primary_landmark)
			data.table::fwrite(pred, prediction_export, compress = "gzip", append = block > 1L)
			final_release_predictions(store, L0, a0)
			rm(pred, groups)
			invisible(gc())
		}
		overlap <- members |>
			select(model, budget, arm, fold, landmark, feature)
		overlap_rows <- list()
		for (L in unique(overlap$landmark)) for (a in unique(overlap$arm)) for (k in unique(overlap$budget)) for (fd in seq_len(K)) {
			z <- overlap[overlap$landmark == L & overlap$arm == a & overlap$budget == k & overlap$fold == fd, , drop = FALSE]
			for (nm in c("Clinical_ProtMet_YS_YinYang", "Clinical_ProtMet_YSplus_YinYang")) {
				u <- z$feature[z$model == nm]
				v <- z$feature[z$model == "Clinical_ProtMet_YS_Yin"]
				if (length(u) && length(v))
					overlap_rows[[length(overlap_rows) + 1L]] <- tibble(
						model = nm, reference = "Clinical_ProtMet_YS_Yin",
						fold = fd, landmark = L, arm = a, budget = k, identical_panel = setequal(u, v), Jaccard = length(intersect(
							u,
							v
						)) / length(union(u, v))
					)
			}
		}
		het <- bind_rows(heterogeneity)
		if (nrow(het))
			het$FDR <- p.adjust(het$p, "BH")
		tables <- list(
			biomarker_PGS_omics_strata = bind_rows(genetic_strata), marker_heterogeneity = het, marker_definitions = combine("marker_audit"),
			domain_support = combine("domain_status"), pgs_reconstruction = combine("pgs_reconstruction"), pgs_decomposition_coefficients = combine("pgs_coefficients"),
			cross_omic_links = combine("cross_omic"), predictor_inventory = combine("members"), metrics = bind_rows(metrics),
			paired_contrasts = bind_rows(pairs), calibration = bind_rows(cal), decision_curves = bind_rows(dec), fold_C_index = bind_rows(foldmetrics),
			panel_members = left_join(members, map, by = "feature"), panel_overlap = bind_rows(overlap_rows), panel_counts = members |>
				left_join(map, by = "feature") |>
				group_by(model, budget, arm, fold, landmark) |>
				summarise(
					total_assays = n_distinct(feature), protein_assays = n_distinct(feature[layer == "prot"]),
					metabolite_assays = n_distinct(feature[layer == "met"]), .groups = "drop"
				), panel_stability = members |>
				group_by(model, budget, arm, landmark, feature) |>
				summarise(selected_folds = n_distinct(fold), selection_frequency = n_distinct(fold) / K, .groups = "drop"),
			fit_diagnostics = diag, coefficients = combine("coefficients"), model_preprocessing = combine("risk_preprocessing"),
			baseline_hazards = combine("baseline_hazards"), proxy_validation = combine("proxy_validation"), proxy_coefficients = combine("proxy_weights"),
			proxy_preprocessing = combine("proxy_preprocessing"), domain_associations = combine("domain_associations"),
			training_screens = combine("screen"), genetic_training_screens = combine("genetic_screen"), PRS_omics_strata = bind_rows(strata),
			inflammatory_burden = bind_rows(burden)
		)
		factorial_tables <- setNames(
			lapply(names(factorial[[1]]), function(nm) bind_rows(lapply(factorial, `[[`, nm))),
			names(factorial[[1]])
		)
		tables <- c(tables, factorial_tables)
		provenance <- read.csv(file.path(root, "prs_provenance.csv"), stringsAsFactors = FALSE)
		prs_scope <- paste(provenance$status, collapse = "; ")
		for (nm in names(tables)[startsWith(names(tables), "factorial_")]) if (nrow(tables[[nm]]))
			tables[[nm]]$disease_PRS_scope <- prs_scope
		for (nm in names(tables)) write_raw_csv(tables[[nm]], paste0("", nm, ".csv"), root)
		primary_rows <- tables$paired_contrasts
		primary_ok <- nrow(primary_rows) > 0 && any(primary_rows$primary %in% TRUE & primary_rows$status == "ok")
		write.csv(data.frame(
			budget = primary_budget, landmark = primary_landmark, comparison = "Clinical_ProtMet_YS_YinYang vs Clinical_ProtMet_NS",
			status = if (primary_ok)
				"estimated; inspect effect and uncertainty" else "not estimable; no alternative primary selected"
		), file.path(root, "primary_status.csv"), row.names = FALSE)
		if (!file.rename(prediction_export, file.path(root, "_private", "out_of_fold_predictions.csv.gz")))
			stop("Cannot publish Final out-of-fold predictions")
		write.csv(
			data.frame(
				signature = signature, mode = if (is.null(cached_validation))
					"current run" else "frozen-checkpoint recovery", folds = K, bootstrap = B, seed = seed, primary_budget = primary_budget,
				primary_landmark = primary_landmark, prediction_blocks = nrow(blocks)
			), file.path(root, "summary_provenance.csv"),
			row.names = FALSE
		)
		writeLines(
			c(
				"Final FINAL — internal, common-cohort, grouped outer validation.", paste(
					"Primary: YS YinYang vs NS; total budget",
					primary_budget, "; landmark", primary_landmark
				), "Fixed end at baseline year 10 by default: landmark 5 evaluates years 5–10 among those still event-free and observed at year 5.",
				"Each landmark retrains and reselects using eligible training participants. This does not identify causal biomarkers.",
				"Clinical includes continuous BMI and non-HDL. Proxy replacement and incremental addition are separate tests.",
				"Proxy preprocessing/selection/fit are cross-fitted; missing OOF values never become in-sample predictions.",
				"Natriuretic/GDF15 omission excludes candidates before equal-budget re-selection; it is not an estimate of causal status.",
				"No automatic inflammatory/non-inflammatory disease subtype is inferred from these markers.", "IPCW assumes independent censoring; death is censored. Risks are net risks, not real-world competing-risk CIFs.",
				"Paired bootstrap conditions on trained models; it does not include full training/selection uncertainty.",
				"Secondary contrasts, strata, domain reconstructions and cell enrichments are exploratory.", "Clinical deployment claims require external validation, recalibration and competing-risk assessment."
			),
			file.path(root, "interpretation.txt")
		)
		tables$metrics$validation_scope <- "Internal measured-omics validation; no external confirmation"
		pgsmodels <- grepl("PGS", tables$metrics$model)
		tables$metrics$validation_scope[pgsmodels] <- paste("Biomarker PGS:", paste(readLines(file.path(root, "omic_PGS_status.txt")),
			collapse = " "
		))
		prsmodels <- grepl("PRS", tables$metrics$model) | tables$metrics$model %in% design$model[design$G]
		tables$metrics$validation_scope[prsmodels] <- paste(
			tables$metrics$validation_scope[prsmodels], "Disease PRS:",
			prs_scope
		)
		write_raw_csv(tables$metrics, "metrics.csv", root)
		final_factorial_plots(tables, design, root, primary_budget, primary_landmark)
		final_joint_plots(tables, root, primary_budget)
		c5_cell_output(tables$panel_members, map, root)
		final_final_panels(tables, root, primary_budget)
		final_evidence_output(tables$panel_members, root)
		request <- if (is.null(cached_validation)) le8_final_request() else {
			previous <- file.path(root, 'res.rds')
			if (file.exists(previous)) readRDS(previous)$training_request else NULL
		}
		final_write_joint_checkpoint(list(
			signature = signature, training_request = request, tables = tables, complete = all(diag$status == "ok"),
			summary_complete = TRUE, summary_settings = list(
				bootstrap = B, seed = seed, primary_budget = primary_budget,
				primary_landmark = primary_landmark
			)
		), file.path(root, "res.rds"))
		final_retire_fold_checkpoints(store, root, signature)
		invisible(tables)
	}
	final_evidence_output <- function(panels, root) {
		# Post-validation annotation only. Never reuse these full-cohort statistics for screening, tuning or
		# claiming independent replication of prediction.
		ledger <- evidence <- list()
		for (layer in c("prot", "met")) {
			candidates <- unique(panels$assay[panels$layer == layer])
			sources <- data.frame(module = c("c1_correlate", "c2_cause", "c2_cause", "c3_coloc", "c4_connect"), file = c(
				if (layer ==
					"prot") "pwas_incident_adj2.csv" else "mwas_incident_adj2.csv", "c2.MR_all.csv", "c2.dandelion_targets_all.csv",
				"c3.coloc_summary.csv", "c4.mediation_all.csv"
			), key = c(
				"term", "exposure", "feature", "feature",
				"feature"
			))
			for (i in seq_len(nrow(sources))) {
				s <- sources[i, ]
				path <- file.path(out.base, layer, s$module, s$file)
				resultfile <- file.path(dirname(path), paste0(substr(s$module, 1, 2), ".res.rds"))
				settings <- if (file.exists(resultfile)) readRDS(resultfile)$meta$analysis_options else NULL
				current <- identical(settings, le8_analysis_options())
				status <- if (!file.exists(path))
					"unavailable" else if (current)
					"current baseline contract" else "legacy or unverified baseline contract; rerun before inference"
				ledger[[length(ledger) + 1L]] <- tibble(layer, module = s$module, file = s$file, status, used_for_prediction = FALSE)
				if (!file.exists(path))
					next
				z <- as.data.frame(data.table::fread(path))
				if (!s$key %in% names(z))
					next
				z <- z[as.character(z[[s$key]]) %in% candidates, , drop = FALSE]
				if (!nrow(z))
					next
				z$assay <- as.character(z[[s$key]])
				z$record <- seq_len(nrow(z))
				names(z)[names(z) == "layer"] <- "source_layer"
				fields <- setdiff(names(z), c("assay", "record"))
				z[fields] <- lapply(z[fields], as.character)
				evidence[[length(evidence) + 1L]] <- tidyr::pivot_longer(z,
					cols = all_of(fields), names_to = "statistic",
					values_to = "value"
				) |>
					mutate(layer, module = s$module, source_file = s$file, baseline_contract = status, interpretation = "Descriptive triangulation; correlated evidence, no causal voting score")
			}
		}
		write_raw_csv(bind_rows(ledger), "evidence_availability.csv", root)
		write_raw_csv(bind_rows(evidence), "evidence_long.csv", root)
	}
	final_joint_plots <- function(tables, root, k) {
		savefig <- function(p, n, w = 11, h = 7) ggplot2::ggsave(file.path(root, n), p, width = w, height = h, dpi = 200)
		m <- tables$metrics |>
			filter(status == "ok", budget %in% c(0, k), arm == "all_assays")
		keep <- c(
			"Clinical", "Clinical_PRS", "Clinical_Protein_NS", "Clinical_Metabolite_NS", "Clinical_ProtMet_NS",
			"Clinical_ProtMet_PRS_NS", "Clinical_ProtMet_YS_YinYang"
		)
		if (nrow(m))
			savefig(ggplot(m |>
				filter(model %in% keep), aes(AUC, reorder(model, AUC), color = factor(landmark))) +
				geom_point(size = 2) +
				labs(x = "Out-of-fold IPCW AUC", y = NULL, color = "Landmark (years)", title = "Joint modality comparisons on the same participants") +
				theme_bw(), "Fig1.joint_modalities.png")
		p <- tables$paired_contrasts
		if (nrow(p))
			savefig(ggplot(p |>
				filter(budget == k, arm == "all_assays"), aes(delta_AUC, paste(model, "vs", reference), color = factor(landmark))) +
				geom_vline(xintercept = 0, linetype = 2) +
				geom_errorbar(aes(xmin = AUC_lo, xmax = AUC_hi),
					orientation = "y",
					width = 0.15
				) +
				geom_point() +
				labs(
					x = "Paired delta AUC (conditional 95% bootstrap CI)", y = NULL,
					color = "Landmark"
				) +
				theme_bw(), "Fig2.paired_increment.png", 13, 8)
		q <- tables$proxy_validation
		if (nrow(q))
			savefig(ggplot(q |>
				filter(arm == "all_assays"), aes(modality, R2_vs_training_mean, color = cohort)) +
				geom_point(position = position_jitter(
					width = 0.08,
					height = 0
				)) +
				facet_grid(target ~ landmark) +
				labs(y = "Held-out reconstruction R²", x = NULL, title = "BMI and lipid proxies: reconstruction is distinct from disease benefit") +
				theme_bw(), "Fig3.domain_proxies.png")
		q <- tables$PRS_omics_strata
		if (nrow(q))
			savefig(ggplot(q |>
				filter(arm == "all_assays"), aes(stratum, observed_risk_ipcw, fill = stratum)) +
				geom_col() +
				facet_wrap( ~ landmark) +
				labs(y = "Observed IPCW net risk", x = NULL, title = "Complementary PRS and omic information") +
				theme_bw() +
				theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "none"), "Fig4.PRS_omics_strata.png")
		q <- tables$calibration |>
			filter(model %in% keep, budget %in% c(0, k), arm == "all_assays")
		if (nrow(q))
			savefig(ggplot(q, aes(predicted, observed_ipcw, color = model)) +
				geom_abline(
					slope = 1, intercept = 0,
					linetype = 2
				) +
				geom_line() +
				geom_point() +
				facet_wrap( ~ landmark) +
				theme_bw() +
				labs(
					x = "Predicted net risk",
					y = "Observed IPCW risk"
				), "Fig5.calibration.png")
		q <- tables$metrics |>
			filter(status == "ok", budget == k, model %in% c("Clinical_ProtMet_NS", "Clinical_ProtMet_YS_YinYang"))
		if (nrow(q))
			savefig(
				ggplot(q, aes(landmark, AUC, color = arm, linetype = model)) +
					geom_line() +
					geom_point() +
					theme_bw() +
					labs(x = "Retraining landmark (years)", y = "IPCW AUC", title = "Longer lead-time and prespecified marker omission"),
				"Fig6.leadtime_omission.png"
			)
	}
	c5_cell_output <- function(panels, map, root) {
		prot <- map[map$layer == "prot", , drop = FALSE]
		# Exact symbols only; NTPROBNP is a fragment of NPPB, not a separate gene.
		prot$gene <- toupper(prot$assay)
		prot$gene[prot$gene == "NTPROBNP"] <- "NPPB"
		mf <- Sys.getenv("FINAL_ASSAY_GENE_MAP", "")
		if (nzchar(mf)) {
			g <- read.csv(mf, stringsAsFactors = FALSE)
			if (!all(c("assay", "gene") %in% names(g)) || anyDuplicated(g$assay))
				stop("Assay-gene map requires unique assay rows")
			prot$gene <- g$gene[match(prot$assay, g$assay)]
		}
		universe <- data.frame(assay = prot$feature, gene = prot$gene)
		p <- panels |>
			filter(layer == "prot") |>
			mutate(model = paste(model, budget, arm, landmark, fold, sep = "|")) |>
			select(model, feature) |>
			distinct()
		write.csv(universe, file.path(root, "c5.cell_assay_universe.csv"), row.names = FALSE)
		write.csv(p, file.path(root, "c5.cell_panels.csv"), row.names = FALSE)
		atlas <- Sys.getenv("C5_CELL_ATLAS", unset = if (Sys.info()[["sysname"]] == "Windows") "F:/annot/cellage/atlas.csv" else "/mnt/f/annot/cellage/atlas.csv")
		code <- file.path(fdir, "c5.cellulation.py")
		status <- system2(Sys.getenv("PYTHON_BIN", "python3"), vapply(c(
			code, "--annotate", "--universe", file.path(
				root,
				"c5.cell_assay_universe.csv"
			), "--atlas", atlas, "--panels", file.path(root, "c5.cell_panels.csv"), "--outdir",
			root, "--prefix", "c5.cell"
		), shQuote, character(1)))
		if (status != 0)
			stop("Cell annotation failed; inspect Python dependencies and assay map")
	}


	# final.R FINAL additions to joint validation. No full-cohort C1 statistics enter here.

	final_endpoint_coverage <- function(dat, Y) {
		dat$.entry <- 0
		mf <- Sys.getenv("LE8_ENDPOINT_MANIFEST", "")
		status <- "No endpoint coverage manifest supplied; baseline coverage must be verified"
		if (nzchar(mf)) {
			m <- read.csv(mf, stringsAsFactors = FALSE)
			if (!all(c("Y", "coverage_start", "coverage_end", "prebaseline_capture") %in% names(m)))
				stop("Incomplete endpoint manifest")
			m <- m[m$Y == Y, , drop = FALSE]
			if (nrow(m) > 1)
				stop("Duplicate endpoint manifest rows")
			if (nrow(m) == 1) {
				a <- as.Date(m$coverage_start)
				b <- as.Date(m$coverage_end)
				if (is.na(a) || is.na(b) || b <= a)
					stop("Verified registry start/end dates are required")
				baseline <- as.Date(dat$date_attend)
				dat$.entry <- pmax(0, as.numeric(a - baseline) / 365.25)
				end <- as.numeric(b - baseline) / 365.25
				dat$.event <- as.integer(dat$.event == 1 & dat$.time <= end)
				dat$.time <- pmin(dat$.time, end)
				status <- paste(
					"Only participants already covered at each landmark enter that risk set; prebaseline capture",
					m$prebaseline_capture
				)
			}
		}
		write.csv(data.frame(Y, status), file.path(final.out, "endpoint_coverage.csv"), row.names = FALSE)
		dat
	}

	final_load_biomarker_pgs <- function(dat, map) {
		mf <- Sys.getenv("FINAL_OMIC_PGS_MANIFEST", "")
		if (nzchar(mf)) {
			# Explicit provenance is preferred; keep independent inputs on their original strict path. Automatic
			# existing UKB PGS is exploratory below.
			result <- final_attach_omic_pgs(dat, map)
			write.csv(read.csv(mf, stringsAsFactors = FALSE), file.path(final.out, "biomarker_PGS_provenance.csv"),
				row.names = FALSE
			)
			return(result)
		}
		source(file.path(fdir, "c1.correlate.R"))
		inputs <- character()
		audit <- list()
		map$pgs_feature <- NA_character_
		for (layer in c("prot", "met")) {
			f <- file.path(indir, "Rdata", paste0(layer, ".pgs.rds"))
			if (!file.exists(f)) {
				audit[[layer]] <- data.frame(layer, status = "not supplied", mapped = 0)
				next
			}
			g <- final_assert_ids(as.data.frame(read_c1_pgs(f)))
			g <- g[g$eid %in% dat$eid, , drop = FALSE]
			mi <- which(map$layer == layer)
			lookup <- map_c1_pgs_columns(map$assay[mi], names(g))
			hit <- mi[map$assay[mi] %in% names(lookup)]
			cols <- unname(lookup[map$assay[hit]])
			if (!length(cols)) {
				audit[[layer]] <- data.frame(layer, status = "no matched scores", mapped = 0)
				next
			}
			if (anyDuplicated(cols) || !all(vapply(g[cols], is.numeric, logical(1))))
				stop("Ambiguous or nonnumeric matched PGS")
			nf <- paste0("pgs__", map$feature[hit])
			map$pgs_feature[hit] <- nf
			g <- g[, c("eid", cols), drop = FALSE]
			names(g) <- c("eid", nf)
			# Retain the measured common cohort. PGS missingness never becomes zero.
			dat <- merge(dat, g, by = "eid", all.x = TRUE, sort = FALSE)
			inputs <- c(inputs, f)
			audit[[layer]] <- data.frame(layer,
				status = "existing biomarker PGS; discovery overlap unverified; exploratory",
				mapped = length(hit)
			)
		}
		if (all(is.na(map$pgs_feature)))
			map$pgs_feature <- NULL
		write.csv(bind_rows(audit), file.path(final.out, "biomarker_PGS_provenance.csv"), row.names = FALSE)
		list(data = dat[order(dat$eid), , drop = FALSE], map = map, inputs = inputs, status = if (length(inputs)) "exploratory: automatic biomarker PGS; discovery independence unverified" else "biomarker PGS unavailable")
	}

	final_replicated_domains <- function(train, features, targets, basic, seed) {
		# Stable donor hashing: adding Yang must not move existing Yin donors between the discovery and replication
		# halves.
		stable_half <- function(g) {
			v <- as.double(seed)
			for (x in utf8ToInt(as.character(g))) v <- (v * 131 + x) %% 2147483629
			as.integer(v %% 2) + 1L
		}
		h <- vapply(train$.group, stable_half, integer(1))
		ranks <- list()
		audit <- list()
		minr <- as.numeric(Sys.getenv("FINAL_DOMAIN_MIN_R", "0.10"))
		if (!is.finite(minr) || minr < 0 || minr > 1)
			stop("FINAL_DOMAIN_MIN_R must be in [0,1]")
		for (tg in targets) {
			if (!tg %in% names(train)) {
				audit[[tg]] <- tibble(domain = tg, status = "missing target", eligible = 0L)
				next
			}
			one <- function(hh) {
				z <- tryCatch(final_domain_rank(train[h == hh, , drop = FALSE], features, tg, basic), error = function(e) tibble())
				if (!nrow(z))
					return(z)
				df <- pmax(1, z$N_observed - length(basic) - 2)
				z$p <- 2 * pt(abs(z$r) * sqrt(df / pmax(1e-12, 1 - z$r ^ 2)), df = df, lower.tail = FALSE)
				z$q <- p.adjust(z$p, "BH", n = length(features))
				z
			}
			a <- one(1)
			b <- one(2)
			if (!nrow(a) || !nrow(b)) {
				audit[[tg]] <- tibble(domain = tg, status = "insufficient training target", eligible = 0L)
				next
			}
			z <- inner_join(a, b, by = c("feature", "domain"), suffix = c("1", "2")) |>
				mutate(eligible = q1 < 0.05 & q2 < 0.05 & sign(r1) == sign(r2) & pmin(abs(r1), abs(r2)) >= minr, strength = pmin(
					abs(r1),
					abs(r2)
				)) |>
				arrange(desc(strength), feature)
			ranks[[tg]] <- z
			audit[[tg]] <- tibble(domain = tg, status = if (any(z$eligible))
				"supported in both training halves" else "no eligible proxies", eligible = sum(z$eligible))
		}
		list(ranks = lapply(ranks, function(z) z[z$eligible, , drop = FALSE]), all = bind_rows(ranks), audit = bind_rows(audit))
	}

	# Orthogonalization is a statistical decomposition, not causal partitioning.  The genetic-predicted part uses
	# only that biomarker's own PGS. The remainder includes acquired biology, untagged genetic effects,
	# confounding and error.
	final_pgs_decompose <- function(train, test, features, map, seed) {
		tr <- train
		te <- test
		gp <- resid <- character()
		aud <- coef <- list()
		folds <- final_group_folds(train$.group, 5, seed)
		for (f in features) {
			g <- map$pgs_feature[match(f, map$feature)]
			if (length(g) != 1 || is.na(g) || !g %in% names(train))
				next
			observed <- is.finite(train[[f]]) & is.finite(train[[g]])
			if (sum(observed) < 100 || sd(train[[g]][observed]) <= 0)
				next
			fit_predict <- function(ii, dd) {
				ok <- observed & ii
				if (sum(ok) < 50)
					return(NULL)
				fit <- lm.fit(cbind(1, train[[g]][ok]), train[[f]][ok])
				b <- fit$coefficients
				if (any(!is.finite(b)))
					return(NULL)
				list(pred = b[1] + b[2] * dd[[g]], coef = b)
			}
			oof <- rep(NA_real_, nrow(train))
			for (h in seq_len(5)) {
				z <- fit_predict(folds != h, train[folds == h, , drop = FALSE])
				if (!is.null(z))
					oof[folds == h] <- z$pred
			}
			obj <- fit_predict(rep(TRUE, nrow(train)), test)
			if (is.null(obj) || any(!is.finite(oof[observed])))
				next
			gn <- paste0("genpart__", f)
			rn <- paste0("remainder__", f)
			tr[[gn]] <- oof
			te[[gn]] <- obj$pred
			tr[[rn]] <- train[[f]] - oof
			te[[rn]] <- test[[f]] - obj$pred
			gp <- c(gp, gn)
			resid <- c(resid, rn)
			ok <- is.finite(test[[f]]) & is.finite(obj$pred)
			obs <- test[[f]][ok]
			pr <- obj$pred[ok]
			den <- sum((obs - mean(train[[f]][observed])) ^ 2)
			aud[[f]] <- tibble(feature = f, pgs = g, N_training = sum(observed), N_validation = sum(ok), R2 = if (sum(ok) >
				2 && den > 0)
				1 - sum((obs - pr) ^ 2) / den else NA_real_, measured_PGS_r = if (sum(ok) > 2)
				cor(obs, pr) else NA_real_, interpretation = "Prediction of adult measurement from matched PGS; remainder is not an environmental or causal effect")
			coef[[f]] <- tibble(feature = f, pgs = g, intercept = obj$coef[1], slope = obj$coef[2])
		}
		list(train = tr, test = te, genetic = gp, remainder = resid, audit = bind_rows(aud), coefficient = bind_rows(coef))
	}

	final_cross_omic_links <- function(train, test, features, basic) {
		pp <- features[startsWith(features, "prot__")]
		mm <- features[startsWith(features, "met__")]
		if (!length(pp) || !length(mm))
			return(tibble())
		getr <- function(d) {
			cv <- le8_prepare_prediction_matrix(d, d, basic)$train
			q <- qr(cbind(1, cv))
			# Pairwise observed assays, with basic-covariate residuals; training and held-out cohorts are kept
			# separate and all requested pairs are retained.
			residual <- function(v) {
				ok <- is.finite(v)
				ans <- rep(NA_real_, length(v))
				if (sum(ok) > length(basic) + 5)
					ans[ok] <- lm.fit(cbind(1, cv[ok, , drop = FALSE]), v[ok])$residuals
				ans
			}
			a <- sapply(d[, pp, drop = FALSE], residual)
			b <- sapply(d[, mm, drop = FALSE], residual)
			if (is.null(dim(a)))
				a <- matrix(a, ncol = 1, dimnames = list(NULL, pp))
			if (is.null(dim(b)))
				b <- matrix(b, ncol = 1, dimnames = list(NULL, mm))
			map_dfr(pp, function(p) map_dfr(mm, function(m) {
				ok <- is.finite(a[, p]) & is.finite(b[, m])
				n <- sum(ok)
				r <- if (n > 20)
					cor(a[ok, p], b[ok, m]) else NA_real_
				df <- max(1, n - length(basic) - 2)
				pv <- if (is.finite(r))
					2 * pt(abs(r) * sqrt(df / max(1e-12, 1 - r * r)), df, lower.tail = FALSE) else NA_real_
				tibble(protein = p, metabolite = m, r = r, N = n, p = pv)
			}))
		}
		z <- inner_join(getr(train), getr(test), by = c("protein", "metabolite"), suffix = c("_training", "_validation"))
		z$q_validation <- p.adjust(z$p_validation, "BH")
		z$replicated_direction <- sign(z$r_training) == sign(z$r_validation)
		z$interpretation <- "Partial correlation; no directionality or protein-to-metabolite mediation established"
		z
	}

	final_marker_strata <- function(train, test, map) {
		definitions <- list(GDF15_IL6_median = list(assay = c("GDF15", "IL6"), q = 0.5), GDF15_IL6_q75 = list(assay = c(
			"GDF15",
			"IL6"
		), q = 0.75), GlycA_median = list(assay = "GLYCA", q = 0.5))
		ans <- data.frame(eid = test$eid)
		audit <- list()
		for (nm in names(definitions)) {
			spec <- definitions[[nm]]
			hit <- map$feature[match(spec$assay, toupper(map$assay))]
			if (anyNA(hit)) {
				ans[[nm]] <- "unavailable"
				audit[[nm]] <- tibble(definition = nm, status = "required marker absent")
				next
			}
			mu <- vapply(train[hit], mean, numeric(1), na.rm = TRUE)
			ss <- vapply(train[hit], sd, numeric(1), na.rm = TRUE)
			if (any(!is.finite(ss) | ss <= 0)) {
				ans[[nm]] <- "unavailable"
				next
			}
			sc <- function(d) rowMeans(sweep(sweep(as.matrix(d[hit]), 2, mu), 2, ss, "/"))
			tr <- sc(train)
			te <- sc(test)
			cut <- unname(quantile(tr, spec$q, na.rm = TRUE))
			ans[[nm]] <- ifelse(!is.finite(te), "unavailable", ifelse(te <= cut, "lower", "higher"))
			audit[[nm]] <- tibble(definition = nm, assay = hit, center = mu, sd = ss, cutoff = cut, status = "training threshold; complete-marker score")
		}
		list(groups = ans, audit = bind_rows(audit))
	}

	final_marker_heterogeneity <- function(a, b, groups, horizon, B, seed) {
		z <- merge(a, b, by = c("eid", "fold", "time", "event"), suffixes = c(".a", ".b"))
		if (nrow(z) != nrow(a) || nrow(z) != nrow(b))
			return(tibble())
		z <- merge(z, groups, by = "eid")
		out <- list()
		for (nm in setdiff(names(groups), "eid")) {
			d <- z[z[[nm]] %in% c("lower", "higher"), , drop = FALSE]
			if (!all(c("lower", "higher") %in% d[[nm]]))
				next
			stat <- function(ix) {
				sapply(c("lower", "higher"), function(g) {
					ii <- ix[d[[nm]][ix] == g]
					w <- le8_ipcw(d$time[ii], d$event[ii], horizon)
					if (w$status != "ok")
						return(NA_real_)
					le8_weighted_auc(d$risk.a[ii], w$y, w$w) - le8_weighted_auc(d$risk.b[ii], w$y, w$w)
				})
			}
			est <- stat(seq_len(nrow(d)))
			set.seed(seed)
			gg <- split(seq_len(nrow(d)), d$group.a)
			bs <- if (B > 0)
				replicate(B, stat(unlist(gg[sample.int(length(gg), length(gg), replace = TRUE)], use.names = FALSE))) else matrix(NA_real_, 2, 0)
			dd <- bs[2, ] - bs[1, ]
			dd <- dd[is.finite(dd)]
			diff <- unname(est[2] - est[1])
			se <- sd(dd)
			out[[nm]] <- tibble(
				definition = nm, N = nrow(d), delta_lower = est[1], delta_higher = est[2], difference = diff,
				lo = if (length(dd) >= 20)
					unname(quantile(dd, 0.025)) else NA_real_, hi = if (length(dd) >= 20)
					unname(quantile(dd, 0.975)) else NA_real_, p = if (length(dd) >= 20 && is.finite(se) && se > 0)
					2 * pnorm(abs(diff / se), lower.tail = FALSE) else NA_real_, interpretation = "Exploratory difference in delta AUC, group bootstrap of frozen fits; marker burden is not a causal subtype"
			)
		}
		bind_rows(out)
	}

	final_genetic_measured_strata <- function(pred, budget) {
		rows <- list()
		for (L in unique(pred$landmark)) for (arm in unique(pred$arm)) {
			a <- pred[pred$model == "MatchedMeasured_only_NS" & pred$budget == budget & pred$landmark == L & pred$arm ==
				arm, , drop = FALSE]
			b <- pred[pred$model == "MatchedPGS_only_NS" & pred$budget == budget & pred$landmark == L & pred$arm ==
				arm, , drop = FALSE]
			if (!nrow(a) || !nrow(b) || !setequal(a$eid, b$eid))
				next
			z <- merge(a, b[, c("eid", "high_training_q75")], by = "eid", suffixes = c(".measured", ".PGS"))
			z$stratum <- paste0("measured ", ifelse(z$high_training_q75.measured, "high", "lower"), " / PGS ", ifelse(z$high_training_q75.PGS,
				"high", "lower"
			))
			for (g in unique(z$stratum)) {
				d <- z[z$stratum == g, , drop = FALSE]
				w <- le8_ipcw(d$time, d$event, d$horizon[1])
				rows[[length(rows) + 1L]] <- tibble(
					stratum = g, landmark = L, arm = arm, N = nrow(d), events = w$N_case,
					observed_risk = if (w$status == "ok")
						sum(w$w * w$y) / sum(w$w) else NA_real_, status = w$status, interpretation = "Training upper-quartile thresholds; same matched biomarkers; acquired or genetic causation is not identified"
				)
			}
		}
		bind_rows(rows)
	}


	# final.R Fixed-panel comparisons of Clinical (C), proteome (P), metabolome (M), and the endpoint's
	# disease PRS (G). Biomarker PGS are a separate analysis.
	final_factorial_design <- function(has_prs = TRUE) {
		z <- expand.grid(C = c(FALSE, TRUE), P = c(FALSE, TRUE), M = c(FALSE, TRUE), G = c(FALSE, TRUE))
		z <- z[rowSums(z) > 0, , drop = FALSE]
		z$model <- apply(z, 1, function(r) paste0("F_", paste(c("C", "P", "M", "G")[as.logical(r)], collapse = "")))
		z$label <- apply(z[c("C", "P", "M", "G")], 1, function(r) paste(c("Clinical", "Protein", "Metabolite", "Disease PRS")[as.logical(r)],
			collapse = " + "
		))
		z$available <- !z$G | has_prs
		z$status <- ifelse(z$available, "scheduled", "endpoint PRS unavailable; model not fitted")
		z$budget_definition <- "k assays PER included measured layer; P+M uses 2k, G is one precomputed endpoint score"
		z[order(rowSums(z[c("C", "P", "M", "G")]), z$model), , drop = FALSE]
	}
	final_factorial_edges <- function(design) {
		z <- design[design$available, , drop = FALSE]
		ans <- list()
		for (i in seq_len(nrow(z))) for (j in seq_len(nrow(z))) {
			a <- as.logical(unlist(z[i, c("C", "P", "M", "G")], use.names = FALSE))
			b <- as.logical(unlist(z[j, c("C", "P", "M", "G")], use.names = FALSE))
			if (sum(a & !b) == 1 && !any(b & !a))
				ans[[length(ans) + 1L]] <- data.frame(model = z$model[i], reference = z$model[j], added_layer = c(
					"C",
					"P", "M", "G"
				)[a & !b], clinical_background = b[1], focus = z$model[i] == "F_CPMG" && z$model[j] %in%
					c("F_CPM", "F_CPG", "F_CMG"))
		}
		do.call(rbind, ans)
	}
	final_factorial_cutpoints <- function() {
		s <- Sys.getenv("FINAL_RISK_CUTS", "")
		if (!nzchar(trimws(s)))
			return(numeric())
		v <- suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1]]))
		if (!length(v) || any(!is.finite(v) | v <= 0 | v >= 1) || anyDuplicated(v))
			stop("FINAL_RISK_CUTS requires unique probabilities strictly between 0 and 1")
		sort(v)
	}
	# Compare ONLY within-fold pairs, because separate Cox models can have different LP scales. Administratively
	# restrict all C-index calculations to the same prediction horizon as AUC/Brier.
	final_within_fold_c <- function(time, event, risk, fold, horizon) {
		count <- c(0, 0, 0)
		for (f in unique(fold)) {
			ii <- which(fold == f)
			if (length(ii) < 2)
				next
			d <- data.frame(tt = pmin(time[ii], horizon), ee = as.integer(event[ii] == 1 & time[ii] <= horizon), rr = risk[ii])
			obj <- tryCatch(survival::concordance(survival::Surv(tt, ee) ~ rr, data = d, reverse = TRUE), error = function(e) NULL)
			if (!is.null(obj))
				count <- count + as.numeric(obj$count[c("concordant", "discordant", "tied.x")])
		}
		if (any(!is.finite(count)) || sum(count) <= 0)
			return(NA_real_)
		(count[1] + 0.5 * count[3]) / sum(count)
	}
	final_factorial_reclassification <- function(new, old, y, w, cuts = numeric()) {
		meanpart <- function(x, case) {
			ii <- y == case & w > 0
			if (!any(ii) || sum(w[ii]) <= 0)
				return(NA_real_)
			sum(w[ii] * x[ii]) / sum(w[ii])
		}
		direction <- sign(new - old)
		ne <- meanpart(direction, 1)
		nn <-  - meanpart(direction, 0)
		out <- c(IDI = meanpart(new - old, 1) - meanpart(new - old, 0), continuous_NRI = ne + nn, NRI_event = ne, NRI_nonevent = nn)
		if (length(cuts)) {
			change <- sign(findInterval(new, cuts) - findInterval(old, cuts))
			ce <- meanpart(change, 1)
			cn <-  - meanpart(change, 0)
			out <- c(out, categorical_NRI = ce + cn, categorical_NRI_event = ce, categorical_NRI_nonevent = cn)
		}
		out
	}
	final_factorial_shapley <- function(values) {
		# All 3! orderings of P/M/G added to C; exact layer-subset decomposition.
		players <- c("P", "M", "G")
		permutations <- list(
			c("P", "M", "G"), c("P", "G", "M"), c("M", "P", "G"), c("M", "G", "P"), c("G", "P", "M"),
			c("G", "M", "P")
		)
		needed <- c("F_C", "F_CP", "F_CM", "F_CG", "F_CPM", "F_CPG", "F_CMG", "F_CPMG")
		if (!all(needed %in% names(values)) || any(!is.finite(values[needed])))
			return(setNames(rep(NA_real_, 3), players))
		out <- setNames(numeric(3), players)
		for (order in permutations) {
			included <- character()
			previous <- "F_C"
			for (p in order) {
				included <- c(included, p)
				current <- paste0("F_C", paste(players[players %in% included], collapse = ""))
				out[p] <- out[p] + (values[[current]] - values[[previous]]) / length(permutations)
				previous <- current
			}
		}
		out
	}
	final_factorial_stat <- function(d, risk, edges, horizon, cuts, lp = risk) {
		iw <- le8_ipcw(d$time, d$event, horizon)
		model <- matrix(NA_real_, ncol(risk), 3, dimnames = list(colnames(risk), c("AUC", "C_index", "Brier")))
		for (j in seq_len(ncol(risk))) {
			r <- risk[, j]
			if (iw$status == "ok") {
				model[j, "AUC"] <- le8_weighted_auc(r, iw$y, iw$w)
				model[j, "Brier"] <- mean(iw$w * (iw$y - r) ^ 2)
			}
			model[j, "C_index"] <- final_within_fold_c(d$time, d$event, lp[, j], d$fold, horizon)
		}
		out <- numeric()
		for (nm in rownames(model)) for (mt in colnames(model)) out[paste("model", nm, mt, sep = "|")] <- model[
			nm,
			mt
		]
		reclass_names <- names(final_factorial_reclassification(c(0.1, 0.2), c(0.1, 0.2), c(0, 1), c(1, 1), cuts))
		for (i in seq_len(nrow(edges))) {
			a <- edges$model[i]
			b <- edges$reference[i]
			for (mt in colnames(model)) out[paste("delta", a, b, mt, sep = "|")] <- model[a, mt] - model[b, mt]
			re <- if (iw$status == "ok")
				final_factorial_reclassification(risk[, a], risk[, b], iw$y, iw$w, cuts) else setNames(rep(NA_real_, length(reclass_names)), reclass_names)
			for (mt in names(re)) out[paste("delta", a, b, mt, sep = "|")] <- re[[mt]]
		}
		for (mt in colnames(model)) {
			value <- model[, mt]
			if (mt == "Brier")
				value <-  - value
			sh <- final_factorial_shapley(value)
			for (p in names(sh)) out[paste("shapley", p, mt, sep = "|")] <- sh[[p]]
		}
		out
	}
	final_factorial_data <- function(pred, design, k, L, arm, expected_ids, K) {
		models <- design$model[design$available]
		lst <- list()
		audit <- list()
		for (nm in models) {
			ds <- design[design$model == nm, ]
			budget <- if (ds$P || ds$M)
				k else 0
			x <- pred[pred$model == nm & pred$budget == budget & pred$landmark == L & pred$arm == arm, , drop = FALSE]
			ok <- nrow(x) == length(expected_ids) && !anyDuplicated(x$eid) && setequal(x$eid, expected_ids) && length(unique(x$fold)) ==
				K && all(is.finite(x$risk)) && all(is.finite(x$lp))
			audit[[nm]] <- data.frame(model = nm, status = if (ok)
				"complete" else "missing/failed outer-fold predictions; excluded from comparisons", N = nrow(x))
			if (ok)
				lst[[nm]] <- x[match(expected_ids, x$eid), , drop = FALSE]
		}
		if (!length(lst))
			return(list(audit = do.call(rbind, audit), data = NULL, risk = NULL))
		d <- lst[[1]]
		for (x in lst) if (!all(vapply(
			c("eid", "fold", "time", "event", "group"), function(v) identical(x[[v]], d[[v]]),
			logical(1)
		)))
			stop("Factorial models do not share the same outcomes/folds/participants")
		risk <- do.call(cbind, lapply(lst, `[[`, "risk"))
		colnames(risk) <- names(lst)
		lp <- do.call(cbind, lapply(lst, `[[`, "lp"))
		colnames(lp) <- names(lst)
		list(data = d, risk = risk, lp = lp, models = lst, audit = do.call(rbind, audit))
	}


	# final.R Paired layer increments and uncertainty from frozen out-of-fold predictions.
	final_factorial_outputs <- function(pred, yin, design, K, B, seed, primary_budget, primary_landmark) {
		cuts <- final_factorial_cutpoints()
		edges_all <- final_factorial_edges(design)
		metrics <- contrasts <- shapley <- status <- strata <- correlations <- roc <- list()
		all_boot <- tolower(Sys.getenv("FINAL_FACTORIAL_BOOT_ALL", "FALSE")) %in% c("true", "1", "yes")
		budgets <- sort(unique(pred$budget[pred$budget > 0 & startsWith(pred$model, "F_")]))
		for (L in sort(unique(pred$landmark))) for (arm in unique(pred$arm)) for (k in budgets) {
			ids <- if (inherits(yin, "final_cached_cohorts"))
				as.character(yin$eid[yin$landmark == L]) else final_landmark(yin, L, Inf)$eid
			z <- final_factorial_data(pred, design, k, L, arm, ids, K)
			id <- data.frame(budget_per_layer = k, landmark = L, arm = arm)
			status[[length(status) + 1L]] <- cbind(id, z$audit)
			if (is.null(z$data))
				next
			d <- z$data
			r <- z$risk
			h <- d$horizon[1]
			edges <- edges_all[edges_all$model %in% colnames(r) & edges_all$reference %in% colnames(r), , drop = FALSE]
			est <- final_factorial_stat(d, r, edges, h, cuts, z$lp)
			nboot <- if (all_boot || (k == primary_budget && arm == "all_assays"))
				B else 0L
			message(
				"Final layer comparisons: L=", L, "; ", arm, "; k/layer=", k, "; N=", nrow(d), "; models=", ncol(r),
				"; paired bootstrap=", nboot
			)
			bs <- matrix(NA_real_, length(est), nboot, dimnames = list(names(est), NULL))
			if (nboot > 0) {
				g <- split(seq_len(nrow(d)), d$group)
				ng <- length(g)
				set.seed(seed + as.integer(100 * L) + k)
				for (b in seq_len(nboot)) {
					ix <- unlist(g[sample.int(ng, ng, replace = TRUE)], use.names = FALSE)
					bs[, b] <- final_factorial_stat(d[ix, , drop = FALSE], r[ix, , drop = FALSE], edges, h, cuts, z$lp[ix, ,
						drop = FALSE
					])[names(est)]
				}
			}
			intervals <- function(nm) {
				v <- bs[nm, ]
				v <- v[is.finite(v)]
				ci <- if (length(v) >= max(20, ceiling(0.8 * nboot)))
					quantile(v, c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
				data.frame(
					estimate = unname(est[nm]), lo = ci[1], hi = ci[2], bootstrap_success = length(v), bootstrap_requested = nboot,
					status = if (is.finite(est[nm]))
						"estimated" else "not estimable"
				)
			}
			for (nm in names(est)) {
				pieces <- strsplit(nm, "|", fixed = TRUE)[[1]]
				row <- cbind(id, horizon = h, intervals(nm))
				row$N <- nrow(d)
				row$events_by_horizon <- sum(d$event == 1 & d$time <= h)
				row$uncertainty <- "paired group bootstrap; fixed trained models; screening/training uncertainty excluded"
				if (pieces[1] == "model") {
					ds <- design[design$model == pieces[2], ]
					row$model <- pieces[2]
					row$metric <- pieces[3]
					row$label <- ds$label
					row$protein_assays <- as.integer(ds$P) * k
					row$metabolite_assays <- as.integer(ds$M) * k
					row$total_assays <- row$protein_assays + row$metabolite_assays
					row$disease_PRS_scores <- as.integer(ds$G)
					row$clinical <- ds$C
					metrics[[length(metrics) + 1L]] <- row
				} else if (pieces[1] == "delta") {
					ee <- edges[edges$model == pieces[2] & edges$reference == pieces[3], ]
					row$model <- pieces[2]
					row$reference <- pieces[3]
					row$metric <- pieces[4]
					row$added_layer <- ee$added_layer
					row$clinical_background <- ee$clinical_background
					row$focus <- ee$focus & k == primary_budget & L == primary_landmark & arm == "all_assays"
					row$cutpoints <- paste(cuts, collapse = ",")
					contrasts[[length(contrasts) + 1L]] <- row
				} else {
					row$layer <- pieces[2]
					row$metric <- if (pieces[3] == "Brier")
						"Brier_reduction" else pieces[3]
					row$interpretation <- "mean predictive increment across all layer-addition orders conditional on Clinical; not causal contribution"
					shapley[[length(shapley) + 1L]] <- row
				}
			}
			if (k != primary_budget)
				next
			iw <- le8_ipcw(d$time, d$event, h)
			if (all(c("F_P", "F_M", "F_G") %in% names(z$models)) && iw$status == "ok") {
				flags <- data.frame(P = z$models$F_P$high_training_q75, M = z$models$F_M$high_training_q75, G = z$models$F_G$high_training_q75)
				labels <- apply(flags, 1, function(x) paste0("P", as.integer(x[1]), " M", as.integer(x[2]), " G", as.integer(x[3])))
				for (s in sort(unique(labels))) {
					ii <- which(labels == s)
					obs <- sum(iw$w[ii] * iw$y[ii]) / sum(iw$w[ii])
					q <- tryCatch(summary(survival::survfit(survival::Surv(time, event) ~ 1, data = d[ii, , drop = FALSE]),
						times = h, extend = FALSE
					), error = function(e) NULL)
					km <- lo <- hi <- NA_real_
					if (!is.null(q) && length(q$surv) == 1) {
						km <- 1 - q$surv
						lo <- 1 - q$upper
						hi <- 1 - q$lower
					}
					strata[[length(strata) + 1L]] <- cbind(id,
						horizon = h, stratum = s, N = length(ii), events = sum(iw$y[ii]),
						observed_ipcw = if (is.finite(obs))
							obs else NA_real_, KM_net_risk = km, KM_lo = lo, KM_hi = hi, interpretation = "1 = above fold-training 75th percentile of the layer-only predicted LP; net risk, no causal subtype"
					)
				}
			}
			# Correlations within folds avoid mixing differently scaled risk scores.
			single <- intersect(c("F_P", "F_M", "F_G"), colnames(r))
			if (length(single) > 1)
				for (a in seq_len(length(single) - 1L)) for (b in seq.int(a + 1, length(single))) for (fd in sort(unique(d$fold))) {
					ii <- d$fold == fd
					rho <- suppressWarnings(cor(z$lp[ii, single[a]], z$lp[ii, single[b]], method = "spearman"))
					correlations[[length(correlations) + 1L]] <- cbind(id,
						fold = fd, score_a = single[a], score_b = single[b],
						N = sum(ii), Spearman_r = rho
					)
				}
			if (iw$status == "ok" && arm == "all_assays")
				for (nm in colnames(r)) {
					p <- r[, nm]
					th <- c(
						Inf, sort(unique(quantile(p, seq(0, 1, length.out = 101), names = FALSE)), decreasing = TRUE),
 - Inf
					)
					aa <- vapply(th, function(t) sum(iw$w * iw$y * (p >= t)) / sum(iw$w * iw$y), numeric(1))
					bb <- vapply(th, function(t) sum(iw$w * (1 - iw$y) * (p >= t)) / sum(iw$w * (1 - iw$y)), numeric(1))
					roc[[length(roc) + 1L]] <- cbind(id, horizon = h, model = nm, data.frame(
						threshold = th, sensitivity = aa,
						false_positive_rate = bb
					))
				}
		}
		list(
			factorial_metrics = bind_rows(metrics), factorial_contrasts = bind_rows(contrasts), factorial_shapley = bind_rows(shapley),
			factorial_status = bind_rows(status), factorial_strata = bind_rows(strata), factorial_score_correlations = bind_rows(correlations),
			factorial_ROC = bind_rows(roc)
		)
	}


	# final.R Two compact six-panel figures; all numerical outputs remain in root CSVs.
	final_factorial_plots <- function(t, design, root, k, L) {
		blank <- function(title) ggplot() +
			annotate("text", x = 0, y = 0, label = "Not estimable / input unavailable") +
			theme_void() +
			labs(title = title)
		style <- theme_bw(base_size = 10) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
		short <- function(x) gsub("([CPMG])(?=[CPMG])", "\\1 + ", sub("^F_", "", x), perl = TRUE)
		forest <- function(d, title, xlab) {
			if (!nrow(d))
				return(blank(title))
			ggplot(d, aes(estimate, reorder(label, estimate))) +
				geom_errorbar(aes(xmin = lo, xmax = hi),
					orientation = "y",
					width = 0.15, na.rm = TRUE
				) +
				geom_point() +
				labs(title = title, x = xlab, y = NULL) +
				style
		}
		save <- function(p, nm) {
			ggsave(file.path(root, paste0(nm, ".png")), p,
				width = 18, height = 12,
				dpi = 240, limitsize = FALSE
			)
		}
		m <- t$factorial_metrics
		d <- t$factorial_contrasts
		s <- t$factorial_shapley
		if (!nrow(m)) {
			save(blank("Factorial comparison unavailable"), "Fig23.PRS_prot_met_prediction")
			return(invisible(NULL))
		}
		m$label <- short(m$model)
		q <- m |>
			filter(budget_per_layer == k, landmark == L, arm == "all_assays", status == "estimated")
		a <- forest(q |>
			filter(metric == "C_index"), "A. All available model combinations", "Within-fold, horizon-restricted Harrell C")
		b <- forest(q |>
			filter(metric == "AUC"), "B. Fixed-horizon discrimination", "Out-of-fold IPCW AUC")
		x <- if (nrow(d))
			d |>
				filter(budget_per_layer == k, landmark == L, arm == "all_assays", model == "F_CPMG", reference %in%
					c("F_CPM", "F_CPG", "F_CMG"), metric == "C_index") |>
				mutate(label = paste0("Add ", added_layer, " to ", short(reference))) else data.frame()
		c <- forest(x, "C. Each layer added after the other two", "Paired delta C (conditional 95% CI)")
		if (nrow(x))
			c <- c + geom_vline(xintercept = 0, linetype = 2)
		clinical_models <- design$model[design$C & design$available]
		x <- t$calibration
		if (nrow(x))
			x <- x |>
				filter(model %in% clinical_models, arm == "all_assays", landmark == L, budget %in% c(0, k))
		e <- if (!nrow(x))
			blank("D. Calibration") else ggplot(x, aes(predicted, observed_ipcw, color = short(model))) +
			geom_abline(
				slope = 1, intercept = 0,
				linetype = 2
			) +
			geom_line() +
			geom_point(size = 1) +
			labs(
				title = "D. Calibration on held-out predictions",
				x = "Predicted net risk", y = "Observed IPCW risk", color = NULL
			) +
			style
		x <- t$decision_curves
		if (nrow(x))
			x <- x |>
				filter(model %in% clinical_models, arm == "all_assays", landmark == L, budget %in% c(0, k))
		f <- if (!nrow(x))
			blank("E. Decision curves") else ggplot(x, aes(threshold, net_benefit, color = short(model))) +
			geom_line() +
			geom_line(aes(y = treat_all),
				color = "grey50", linetype = 2
			) +
			geom_hline(yintercept = 0, linetype = 3) +
			labs(
				title = "E. Exploratory net benefit",
				x = "Risk threshold", y = "Net benefit", color = NULL
			) +
			style
		x <- if (nrow(s))
			s |>
				filter(budget_per_layer == k, landmark == L, arm == "all_assays", metric == "AUC", status == "estimated") |>
				mutate(label = layer) else data.frame()
		g <- forest(x, "F. Average contribution across addition orders", "AUC gain allocated across P / M / G")
		if (nrow(x))
			g <- g + geom_vline(xintercept = 0, linetype = 2)
		caption <- paste0(
			"C=clinical; P=proteins; M=metabolites; G=disease PRS, distinct from biomarker PGS. k=",
			k, " per measured layer; P+M uses ", 2 * k, " assays. Shared participants and folds. Landmark ", L, ". Conditional bootstrap CI excludes training uncertainty. G discovery overlap: see prs_provenance.csv. Death censored; net risks, not competing-risk CIFs."
		)
		caption <- paste(strwrap(caption, width = 160), collapse = "\n")
		annotation_theme <- theme(plot.caption = element_text(hjust = 0, size = 9), plot.title = element_text(
			size = 15,
			face = "bold"
		))
		save((a | b | c) / (e | f | g) + plot_annotation(
			title = paste(Y, "— disease PRS, proteome and metabolome"),
			caption = caption, theme = annotation_theme
		), "Fig23.PRS_prot_met_prediction")
		x <- t$factorial_ROC
		a <- blank("A. Time-dependent ROC")
		if (nrow(x)) {
			x <- x |>
				filter(landmark == L, model %in% clinical_models)
			if (nrow(x))
				a <- ggplot(x, aes(false_positive_rate, sensitivity, color = short(model))) +
					geom_line() +
					geom_abline(
						slope = 1,
						intercept = 0, linetype = 2
					) +
					coord_equal() +
					labs(
						title = "A. IPCW ROC on the common cohort",
						x = "False positive rate", y = "Sensitivity", color = NULL
					) +
					style
		}
		x <- t$factorial_score_correlations
		b <- blank("B. Layer-score correlations")
		if (nrow(x)) {
			x <- x |>
				filter(landmark == L, arm == "all_assays") |>
				group_by(score_a, score_b) |>
				summarise(r = mean(Spearman_r, na.rm = TRUE), .groups = "drop")
			if (nrow(x))
				b <- ggplot(x, aes(short(score_a), short(score_b), fill = r)) +
					geom_tile() +
					geom_text(aes(label = sprintf(
						"%.2f",
						r
					))) +
					scale_fill_gradient2(low = "#326da8", mid = "white", high = "#ba493e", limits = c( - 1, 1)) +
					labs(title = "B. Mean within-fold Spearman correlation", x = NULL, y = NULL, fill = "r") +
					style
		}
		x <- t$factorial_strata
		c <- blank("C. Three-layer risk strata")
		if (nrow(x)) {
			x <- x |>
				filter(landmark == L, arm == "all_assays")
			if (nrow(x))
				c <- ggplot(x, aes(stratum, KM_net_risk)) +
					geom_col(fill = "#437a9c") +
					geom_errorbar(aes(
						ymin = KM_lo,
						ymax = KM_hi
					), width = 0.15, na.rm = TRUE) +
					geom_text(aes(label = paste0(events, "/", N)),
						vjust =  - 0.4,
						size = 2.5
					) +
					labs(
						title = "C. High/lower layer-score combinations", x = "1 = above training 75th percentile; labels = events / N",
						y = "KM net risk (95% CI)"
					) +
					style +
					theme(axis.text.x = element_text(angle = 45, hjust = 1))
		}
		x <- m |>
			filter(budget_per_layer == k, arm == "all_assays", metric == "C_index", model %in% clinical_models, status ==
				"estimated")
		e <- if (!nrow(x))
			blank("D. Lead time") else ggplot(x, aes(landmark, estimate, color = short(model))) +
			geom_line() +
			geom_point() +
			geom_errorbar(aes(
				ymin = lo,
				ymax = hi
			), width = 0.1, na.rm = TRUE) +
			labs(
				title = "D. Refit at each event-free landmark", x = "Years since sampling",
				y = "Horizon-restricted C-index", color = NULL
			) +
			style
		x <- if (nrow(d))
			d |>
				filter(budget_per_layer == k, landmark == L, arm == "all_assays", model == "F_CPMG", reference %in%
					c("F_CPM", "F_CPG", "F_CMG"), metric %in% c("IDI", "continuous_NRI", "categorical_NRI")) |>
				mutate(label = paste("Add", added_layer)) else data.frame()
		f <- if (!nrow(x))
			blank("E. Reclassification") else ggplot(x, aes(estimate, label)) +
			geom_vline(xintercept = 0, linetype = 2) +
			geom_errorbar(aes(
				xmin = lo,
				xmax = hi
			), orientation = "y", width = 0.1, na.rm = TRUE) +
			geom_point() +
			facet_wrap( ~ metric, scales = "free_x") +
			labs(title = "E. Paired reclassification (exploratory)", x = "Estimate and conditional 95% CI", y = NULL) +
			style
		x <- m |>
			filter(
				landmark == L, arm == "all_assays", metric == "C_index", model %in% c("F_CP", "F_CM", "F_CPM", "F_CPMG"),
				status == "estimated"
			)
		g <- if (!nrow(x))
			blank("F. Assay requirements") else ggplot(x, aes(total_assays, estimate, color = short(model))) +
			geom_line() +
			geom_point() +
			labs(
				title = "F. Performance versus actual assay count",
				x = "Total measured assays (PGS reported separately)", y = "Horizon-restricted C-index", color = NULL
			) +
			style
		save((a | b | c) / (e | f | g) + plot_annotation(
			title = paste(Y, "— complementarity, lead time and assay requirements"),
			caption = caption, theme = annotation_theme
		), "Fig24.PRS_prot_met_complementarity")
	}


	# final.R Compact multi-panel figures supplement (and retain) the existing figures.
	final_final_panels <- function(t, root, k) {
		blank <- function(title) ggplot() +
			annotate("text", x = 0, y = 0, label = "Not estimable / input unavailable") +
			theme_void() +
			labs(title = title)
		style <- theme_bw(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
		save <- function(p, name) {
			ggsave(file.path(root, paste0(name, ".png")), p,
				width = 18, height = 12,
				dpi = 240, limitsize = FALSE
			)
		}
		short <- function(x) gsub("Clinical_", "C + ", gsub("ProtMet", "P+M", gsub("Matched", "", x)))
		m <- t$metrics |>
			filter(status == "ok", arm == "all_assays", budget %in% c(0, k))
		keep <- c(
			"Clinical", "Clinical_PRS", "Clinical_Protein_NS", "Clinical_Metabolite_NS", "Clinical_ProtMet_NS",
			"Clinical_ProtMet_PRS_NS"
		)
		a <- if (!nrow(m))
			blank("A. Joint modality comparisons") else ggplot(m |>
			filter(model %in% keep), aes(AUC, short(model), color = factor(landmark))) +
			geom_point() +
			labs(
				title = "A. Same participants; declared assay budgets",
				x = "Out-of-fold IPCW AUC", y = NULL, color = "Landmark"
			) +
			style
		d <- t$paired_contrasts |>
			filter(budget == k, arm == "all_assays", grepl("PGS|sharedBudget", paste(model, reference)))
		b <- if (!nrow(d))
			blank("B. Incremental PGS / modality information") else ggplot(d, aes(delta_AUC, paste(short(model), "vs", short(reference)), color = factor(landmark))) +
			geom_vline(
				xintercept = 0,
				linetype = 2
			) +
			geom_errorbar(aes(xmin = AUC_lo, xmax = AUC_hi), orientation = "y", width = 0.15) +
			geom_point() +
			labs(
				title = "B. Paired increments; measured reference matches PGS coverage", x = "Delta AUC with conditional 95% CI",
				y = NULL, color = "Landmark"
			) +
			style
		d <- t$pgs_reconstruction
		c <- if (!nrow(d))
			blank("C. Matched PGS capture of measured omics") else ggplot(d |>
			filter(arm == "all_assays"), aes(sub("__.*$", "", feature), R2, color = factor(landmark))) +
			geom_hline(
				yintercept = 0,
				linetype = 2
			) +
			geom_boxplot(outlier.shape = NA) +
			geom_point(position = position_jitter(
				width = 0.1, height = 0,
				seed = 2026
			), alpha = 0.45, size = 0.7) +
			labs(
				title = "C. Adult measurement predicted from its PGS", x = NULL,
				y = "Held-out reconstruction R²", color = "Landmark"
			) +
			style
		d <- m |>
			filter(grepl("PGS_captured|PGS_remainder|MatchedMeasured|MatchedPGS", model))
		e <- if (!nrow(d))
			blank("D. Genetic-score captured / remaining information") else ggplot(d, aes(AUC, short(model), color = factor(landmark))) +
			geom_point() +
			labs(
				title = "D. Statistical decomposition; not causal partition",
				x = "Out-of-fold IPCW AUC", y = NULL, color = "Landmark"
			) +
			style
		d <- t$cross_omic_links
		f <- if (!nrow(d))
			blank("E. Protein–metabolite connections") else {
			z <- d |>
				filter(arm == "all_assays") |>
				group_by(protein, metabolite) |>
				summarise(r = mean(r_validation, na.rm = TRUE), .groups = "drop")
			ggplot(z, aes(sub("met__", "", metabolite), sub("prot__", "", protein), fill = r)) +
				geom_tile() +
				scale_fill_gradient2(
					low = "#4374a9",
					mid = "white", high = "#ba514e", limits = c( - 1, 1)
				) +
				labs(
					title = "E. Held-out partial correlations, selected pairs",
					x = NULL, y = NULL, fill = "Mean r"
				) +
				style +
				theme(axis.text.x = element_text(
					angle = 60, hjust = 1,
					size = 6
				), axis.text.y = element_text(size = 6))
		}
		d <- t$marker_heterogeneity
		g <- if (!nrow(d))
			blank("F. Inflammation-stratum heterogeneity") else ggplot(d |>
			filter(arm == "all_assays"), aes(difference, paste(definition, short(model)), color = factor(landmark))) +
			geom_vline(xintercept = 0, linetype = 2) +
			geom_errorbar(aes(xmin = lo, xmax = hi),
				orientation = "y",
				width = 0.15
			) +
			geom_point() +
			labs(
				title = "F. Difference in improvement: higher minus lower burden",
				x = "Difference in delta AUC", y = NULL, color = "Landmark"
			) +
			style
		save(
			(a | b | c) / (e | f | g) + plot_annotation(
				title = paste(Y, "— measured proteome, metabolome and matched biomarker PGS"),
				caption = "P=proteins; M=metabolites; C=clinical. Disease PRS is separate from biomarker PGS. Automatic UKB-derived PGS are exploratory until discovery overlap is resolved. All contrasts and failures are retained in CSV outputs."
			),
			"Fig20.integrated_omics_PGS"
		)
		d <- t$domain_support
		a <- if (!nrow(d))
			blank("A. Supported domains") else {
			z <- d |>
				filter(arm == "all_assays") |>
				group_by(domain, cohort) |>
				summarise(n = median(eligible), .groups = "drop")
			ggplot(z, aes(cohort, domain, fill = log1p(n))) +
				geom_tile() +
				geom_text(aes(label = n)) +
				labs(
					title = "A. Training-replicated domain candidates",
					x = NULL, y = NULL, fill = "log(1+n)"
				) +
				style
		}
		d <- t$proxy_validation
		b <- if (!nrow(d))
			blank("B. Proxy reconstruction") else ggplot(d |>
			filter(arm == "all_assays"), aes(modality, R2_vs_training_mean, color = cohort)) +
			geom_point(position = position_jitter(
				width = 0.12,
				height = 0, seed = 2026
			)) +
			facet_wrap( ~ target) +
			labs(
				title = "B. Reconstruct measured domains in held-out participants",
				x = NULL, y = "R² vs training mean"
			) +
			style
		d <- m |>
			filter(grepl("Proxy|^Clinical$", model))
		c <- if (!nrow(d))
			blank("C. Proxy addition versus replacement") else ggplot(d, aes(AUC, short(model), color = factor(landmark))) +
			geom_point() +
			labs(
				title = "C. Disease benefit is tested separately",
				x = "Out-of-fold IPCW AUC", y = NULL, color = "Landmark"
			) +
			style
		d <- t$panel_counts |>
			filter(budget == k, arm == "all_assays", model %in% c(
				"Clinical_ProtMet_NS", "Clinical_ProtMet_YS_YinYang",
				"Clinical_ProtMet_YSplus_YinYang"
			))
		e <- if (!nrow(d))
			blank("D. Actual measured assay requirements") else ggplot(d, aes(short(model), total_assays, color = factor(landmark))) +
			geom_point(position = position_jitter(
				width = 0.1,
				height = 0, seed = 2026
			)) +
			labs(
				title = "D. Score count is not assay count", x = NULL, y = "Measured assays",
				color = "Landmark"
			) +
			style +
			theme(axis.text.x = element_text(angle = 20, hjust = 1))
		d <- t$coefficients |>
			filter(budget == k, arm == "all_assays", fold == 1, model == "Clinical_ProtMet_YS_YinYang")
		f <- if (!nrow(d))
			blank("E. Transparent model coefficients") else ggplot(d, aes(beta, variable, color = factor(landmark))) +
			geom_vline(xintercept = 0, linetype = 2) +
			geom_point() +
			labs(
				title = "E. Standardized coefficients, outer fold 1", x = "Log-hazard coefficient (training scale)",
				y = NULL, color = "Landmark"
			) +
			style
		cf <- file.path(root, "c5.cell.coverage.csv")
		d <- if (file.exists(cf))
			read.csv(cf) else data.frame()
		g <- blank("F. External cell-type interpretation")
		# Coverage is displayed only after inspecting the annotation output schema.
		if (nrow(d) && all(c("model", "cell_labelled_genes", "unique_genes") %in% names(d))) {
			d <- d |>
				filter(grepl(paste0("[|]", k, "[|]all_assays[|]"), model)) |>
				mutate(method = sub("[|].*$", "", model)) |>
				filter(method %in% c("Clinical_ProtMet_NS", "Clinical_ProtMet_YS_YinYang", "Clinical_ProtMet_YSplus_YinYang"))
			if (nrow(d))
				g <- ggplot(d, aes(short(method), cell_labelled_genes / pmax(1, unique_genes))) +
					geom_boxplot() +
					geom_point(position = position_jitter(
						width = 0.1,
						height = 0, seed = 2026
					)) +
					labs(
						title = "F. External cell-label coverage across folds", x = NULL,
						y = "Labelled / unique panel genes"
					) +
					style +
					theme(axis.text.x = element_text(angle = 20, hjust = 1))
		}
		save(
			(a | b | c) / (e | f | g) + plot_annotation(
				title = paste(Y, "— supported LE8 domains and interpretable panels"),
				caption = "One to eight candidate domains; unsupported domains contribute no forced proxies. BMI/non-HDL reconstruction does not demonstrate a better clinical measure. No inflammatory-driven subtype is inferred from marker thresholds."
			),
			"Fig21.domains_interpretation"
		)
		keep <- c("Clinical", "Clinical_ProtMet_NS", "Clinical_ProtMet_YS_YinYang")
		d <- t$metrics |>
			filter(status == "ok", model %in% keep, budget %in% c(0, k))
		a <- if (!nrow(d))
			blank("A. Lead time and marker omission") else ggplot(d, aes(landmark, AUC, color = short(model), linetype = arm)) +
			geom_line() +
			geom_point() +
			labs(
				title = "A. Retraining at each landmark",
				x = "Event-free landmark after sampling (years)", y = "IPCW AUC", color = NULL, linetype = NULL
			) +
			style
		d <- t$calibration
		b <- if (!nrow(d))
			blank("B. Calibration") else ggplot(d |>
			filter(model %in% keep, budget %in% c(0, k), arm == "all_assays"), aes(predicted, observed_ipcw, color = short(model))) +
			geom_abline(slope = 1, intercept = 0, linetype = 2) +
			geom_line() +
			geom_point() +
			facet_wrap( ~ landmark) +
			labs(title = "B. Out-of-fold calibration", x = "Predicted net risk", y = "Observed IPCW net risk", color = NULL) +
			style
		d <- t$decision_curves
		c <- if (!nrow(d))
			blank("C. Decision curves") else ggplot(d |>
			filter(model %in% keep, budget %in% c(0, k), arm == "all_assays"), aes(threshold, net_benefit, color = short(model))) +
			geom_hline(yintercept = 0, linetype = 2) +
			geom_line() +
			facet_wrap( ~ landmark) +
			labs(
				title = "C. Exploratory net benefit",
				x = "Risk threshold", y = "Net benefit", color = NULL
			) +
			style
		d <- t$biomarker_PGS_omics_strata
		e <- if (!nrow(d))
			blank("D. Measured score × biomarker PGS") else ggplot(d |>
			filter(arm == "all_assays"), aes(stratum, observed_risk, fill = stratum)) +
			geom_col() +
			geom_text(aes(label = paste0(
				"N=",
				N
			)), vjust =  - 0.3, size = 2.7) +
			facet_wrap( ~ landmark) +
			labs(
				title = "D. Complementary measured and matched genetic scores",
				x = NULL, y = "Observed IPCW net risk"
			) +
			style +
			theme(axis.text.x = element_text(
				angle = 35, hjust = 1,
				size = 7
			), legend.position = "none")
		save(
			(a | b) / (c | e) + plot_annotation(
				title = paste(Y, "— temporal validation and genetic-score strata"),
				caption = "Each landmark uses a different eligible risk set and a shorter horizon to the fixed baseline-year endpoint. Cutoffs are learned in training. Death is censored: these displays estimate net risk, not competing-risk incidence."
			),
			"Fig22.temporal_calibration_genetic_strata"
		)
	}


	# final.R Disk-backed aggregation: retain only one fold while staging and one landmark/ablation block
	# while evaluating. No participant rows are dropped.  Only known model-copy filenames for a durable, readable
	# fold are eligible.  Never remove fold checkpoints, roles, predictions, results, or newer files from a
	# concurrently recomputing fold. Unknown files and symlinks stay intact.
	final_prune_fold_models <- function(checkpoint, signature = NULL, object = NULL, dry_run = FALSE) {
		empty <- data.frame(file = character(), bytes = numeric(), action = character())
		if (!file.exists(checkpoint))
			return(empty)
		if (basename(dirname(checkpoint)) != "_cache" || nzchar(Sys.readlink(checkpoint)))
			return(empty)
		info <- file.info(checkpoint)
		if (is.na(info$mtime) || is.na(info$size) || info$isdir)
			return(empty)
		z <- if (is.null(object))
			tryCatch(readRDS(checkpoint), error = function(e) NULL) else object
		if (!is.list(z) || !is.character(z$signature) || length(z$signature) != 1L || is.na(z$signature) || !nzchar(z$signature) ||
			(!is.null(signature) && !identical(z$signature, signature)))
			return(empty)
		r <- z$result
		if (!is.list(r) || !is.data.frame(r$predictions) || !nrow(r$predictions) || !is.data.frame(r$coefficients) ||
			!nrow(r$coefficients) || !is.data.frame(r$marker_groups) || !nrow(r$marker_groups))
			return(empty)
		mg <- r$marker_groups
		if (!all(c("eid", "landmark", "fold") %in% names(mg)))
			return(empty)
		L <- unique(mg$landmark)
		fd <- unique(mg$fold)
		if (length(L) != 1L || length(fd) != 1L || anyNA(c(L, fd)) || anyNA(mg$eid) || anyDuplicated(mg$eid) || !all(c(
			"eid",
			"landmark", "fold"
		) %in% names(r$predictions)) || anyNA(r$predictions[c("eid", "landmark", "fold")]) ||
			!all(r$predictions$eid %in% mg$eid) || !all(r$predictions$fold == fd) || !all(r$predictions$landmark ==
			L) || basename(checkpoint) != paste0("L", L, ".fold", fd, ".rds"))
			return(empty)
		prefix <- paste0("L", L, ".fold", fd, ".")
		names <- character()
		d <- r$diagnostics
		if (is.data.frame(d) && all(c("model", "arm", "budget", "status") %in% names(d))) {
			d <- d[d$status %in% "ok" & is.finite(d$budget) & !is.na(d$model) & !is.na(d$arm), , drop = FALSE]
			names <- paste0(prefix, d$arm, ".", d$model, ".k", d$budget, ".rds")
		}
		d <- r$proxy_validation
		if (is.data.frame(d) && nrow(d) && all(c("arm", "cohort", "modality", "target") %in% names(d)))
			names <- c(names, paste0(prefix, d$arm, ".proxy.", d$cohort, ".", d$modality, ".", d$target, ".rds"))
		names <- unique(names[!grepl("[/\\\\]", names)])
		root <- dirname(dirname(checkpoint))
		private <- file.path(root, "_private")
		if (nzchar(Sys.readlink(private)))
			return(empty)
		files <- file.path(private, names)
		fi <- file.info(files)
		keep <- !is.na(fi$size) & !is.na(fi$mtime) & !fi$isdir & fi$mtime <= info$mtime & !nzchar(Sys.readlink(files))
		files <- files[keep]
		fi <- fi[keep, , drop = FALSE]
		if (!length(files))
			return(empty)
		ans <- data.frame(file = files, bytes = fi$size, action = if (dry_run)
			"eligible" else "removed")
		if (dry_run)
			return(ans)
		# Recheck at deletion time so recently replaced/open-for-writing copies stay.
		deleted <- vapply(seq_along(files), function(i) {
			now <- file.info(files[i])
			durable <- file.info(checkpoint)
			unchanged <- isTRUE(now$size == fi$size[i] && now$mtime == fi$mtime[i] && durable$size == info$size &&
				durable$mtime == info$mtime && !nzchar(Sys.readlink(files[i])))
			unchanged && unlink(files[i]) == 0L && !file.exists(files[i])
		}, logical(1))
		ans$action[!deleted] <- "kept: changed or removal failed"
		audit <- transform(ans, checkpoint = basename(checkpoint), time = format(Sys.time(), "%F %T %z"))
		log <- file.path(root, "temporary_cleanup.csv")
		write.table(audit, log, sep = ",", row.names = FALSE, col.names = !file.exists(log), append = file.exists(log))
		message(
			"Final cleanup: removed ", sum(deleted), " redundant model copies; ", sprintf("%.2f", sum(ans$bytes[deleted]) / 1024 ^ 3),
			" GiB; retained ", basename(checkpoint)
		)
		ans
	}

	final_process_token <- function(pid = Sys.getpid()) {
		path <- file.path("/proc", as.character(pid), "stat")
		if (!file.exists(path))
			return(NA_character_)
		tryCatch(
			{
				fields <- strsplit(sub("^.*[)] ", "", suppressWarnings(readLines(path, warn = FALSE))[1]), " +")[[1]]
				if (length(fields) < 20L)
					NA_character_ else fields[20]
			},
			error = function(e) NA_character_
		)
	}
	final_boot_id <- function() {
		tryCatch(suppressWarnings(readLines("/proc/sys/kernel/random/boot_id", warn = FALSE)[1]), error = function(e) NA_character_)
	}
	final_summary_directory <- function(private) {
		dir.create(private, recursive = TRUE, showWarnings = FALSE)
		d <- tempfile(paste0(".summary-", Sys.getpid(), "-"), tmpdir = private)
		dir.create(d, mode = "0700")
		saveRDS(
			list(pid = Sys.getpid(), host = Sys.info()[["nodename"]], token = final_process_token(), boot = final_boot_id()),
			file.path(d, ".owner.rds")
		)
		d
	}
	final_cleanup_stale_summaries <- function(private) {
		if (!dir.exists("/proc") || !dir.exists(private))
			return(invisible(NULL))
		dirs <- list.files(private, pattern = "^[.]summary-", full.names = TRUE, all.files = TRUE)
		for (d in dirs) {
			if (!dir.exists(d) || nzchar(Sys.readlink(d)))
				next
			o <- tryCatch(suppressWarnings(readRDS(file.path(d, ".owner.rds"))), error = function(e) NULL)
			if (!is.list(o) || !identical(o$host, Sys.info()[["nodename"]]) || length(o$pid) != 1L || !is.numeric(o$pid) ||
				!is.finite(o$pid) || o$pid < 1 || length(o$token) != 1L || !is.character(o$token) || is.na(o$token))
				next
			current <- final_process_token(o$pid)
			# An unreadable live /proc entry is ambiguous; keep it. A changed start token identifies PID reuse and
			# must not preserve a dead run's scratch.
			dead <- !dir.exists(file.path("/proc", as.character(o$pid)))
			boot <- final_boot_id()
			reboot <- is.character(o$boot) && length(o$boot) == 1L && !is.na(o$boot) && !is.na(boot) && !identical(
				boot,
				o$boot
			)
			if (dead || reboot || (!is.na(current) && !identical(current, o$token))) {
				unlink(d, recursive = TRUE)
				message("Final cleanup: removed abandoned summary scratch ", basename(d))
			}
		}
		invisible(NULL)
	}

	final_write_joint_checkpoint <- function(object, path) {
		tmp <- tempfile(".fold-", tmpdir = dirname(path))
		on.exit(unlink(tmp), add = TRUE)
		saveRDS(object, tmp)
		if (!file.rename(tmp, path))
			stop("Cannot publish Final fold checkpoint: ", path)
	}

	# A completed summary, including plots/annotations and the published prediction export, replaces fold
	# checkpoints. Failed summaries never reach this cleanup.
	final_retire_fold_checkpoints <- function(store, root, signature) {
		empty <- data.frame(file = character(), bytes = numeric(), action = character())
		final <- tryCatch(readRDS(file.path(root, "res.rds")), error = function(e) NULL)
		if (!isTRUE(final$summary_complete) || !identical(final$signature, signature) || !file.exists(file.path(
			root,
			"_private", "out_of_fold_predictions.csv.gz"
		)))
			return(empty)
		inventory <- store$checkpoints
		if (!is.data.frame(inventory) || !nrow(inventory))
			return(empty)
		cache <- normalizePath(file.path(root, "_cache"), mustWork = FALSE)
		rows <- lapply(seq_len(nrow(inventory)), function(i) {
			p <- inventory$file[i]
			now <- file.info(p)
			safe <- identical(normalizePath(dirname(p), mustWork = FALSE), cache) && grepl(
				"^L[0-9.]+[.]fold[0-9]+[.]rds$",
				basename(p)
			) && !nzchar(Sys.readlink(p)) && !nzchar(Sys.readlink(dirname(p))) && isTRUE(now$size ==
				inventory$size[i] && now$mtime == inventory$mtime[i])
			removed <- safe && unlink(p) == 0L && !file.exists(p)
			data.frame(file = p, bytes = inventory$size[i], action = if (removed)
				"removed" else "kept: changed or unavailable")
		})
		ans <- do.call(rbind, rows)
		log <- file.path(root, "fold_cache_cleanup.csv")
		write.table(transform(ans, time = format(Sys.time(), "%F %T %z")), log,
			sep = ",", row.names = FALSE, col.names = !file.exists(log),
			append = file.exists(log)
		)
		message("Final cleanup: summary outputs published; removed ", sum(ans$action == "removed"), " consumed fold checkpoints")
		ans
	}

	final_stage_joint_results <- function(paths, signature, directory) {
		dir.create(directory, recursive = TRUE, showWarnings = FALSE, mode = "0700")
		checkpoint_info <- file.info(paths)
		checkpoint_info$file <- paths
		tables <- blocks <- cohorts <- list()
		seen <- character()
		for (i in seq_along(paths)) {
			message("Final summary: read fold ", i, "/", length(paths), " (", basename(paths[i]), ")")
			z <- tryCatch(readRDS(paths[i]), error = function(e) stop(
				"Unreadable Final checkpoint: ", paths[i], ": ",
				conditionMessage(e)
			))
			if (!identical(z$signature, signature))
				stop("Final checkpoint signature mismatch: ", paths[i])
			r <- z$result
			p <- r$predictions
			mg <- r$marker_groups
			if (!is.data.frame(p) || !nrow(p) || !is.data.frame(mg) || !nrow(mg) || !all(c("eid", "fold", "landmark") %in%
				names(mg)))
				stop("Incomplete Final checkpoint: ", paths[i])
			L <- unique(mg$landmark)
			fd <- unique(mg$fold)
			if (length(L) != 1L || length(fd) != 1L || anyNA(mg$eid) || anyDuplicated(mg$eid) || !all(p$fold == fd) ||
				!all(p$landmark == L) || !all(p$eid %in% mg$eid))
				stop("Inconsistent Final checkpoint participant/fold inventory: ", paths[i])
			key <- paste(L, fd, sep = "/")
			if (key %in% seen)
				stop("Duplicate Final landmark/fold: ", key)
			seen <- c(seen, key)
			cohorts[[i]] <- mg[, c("eid", "fold", "landmark"), drop = FALSE]
			for (nm in setdiff(names(r), "predictions")) {
				file <- file.path(directory, paste0("fold-", i, "-", nm, ".rds"))
				saveRDS(r[[nm]], file)
				tables[[nm]] <- c(tables[[nm]], file)
			}
			for (a in unique(p$arm)) {
				file <- file.path(directory, paste0("prediction-", length(blocks) + 1L, ".rds"))
				saveRDS(p[p$arm == a, , drop = FALSE], file)
				blocks[[length(blocks) + 1L]] <- data.frame(landmark = L, arm = a, fold = fd, file = file)
			}
			final_prune_fold_models(paths[i], signature, object = z)
			rm(z, r, p, mg)
			invisible(gc())
		}
		if (!length(blocks))
			stop("No Final predictions to summarize")
		list(tables = tables, blocks = bind_rows(blocks), cohorts = bind_rows(cohorts), checkpoints = checkpoint_info[,
			c("file", "size", "mtime"),
			drop = FALSE
		])
	}

	final_store_table <- function(store, name) {
		bind_rows(lapply(store$tables[[name]], readRDS))
	}
	final_store_predictions <- function(store, landmark, arm) {
		b <- store$blocks
		bind_rows(lapply(b$file[b$landmark == landmark & b$arm == arm], readRDS))
	}
	final_release_predictions <- function(store, landmark, arm) {
		b <- store$blocks
		unlink(b$file[b$landmark == landmark & b$arm == arm])
		invisible(NULL)
	}
	final_output_ids <- function(yin, L) {
		if (inherits(yin, "final_cached_cohorts"))
			return(as.character(yin$eid[yin$landmark == L]))
		as.character(final_landmark(yin, L, Inf)$eid)
	}

	final_validate_cached_cohorts <- function(cohorts, availability, roles, K) {
		if (!all(c("eid", ".fold") %in% names(roles)) || anyNA(roles$eid) || anyDuplicated(roles$eid))
			stop("Invalid saved Final outer-fold roles")
		if (!setequal(unique(cohorts$landmark), availability$landmark))
			stop("Missing Final landmark checkpoints")
		for (i in seq_len(nrow(availability))) {
			a <- availability[i, ]
			d <- cohorts[cohorts$landmark == a$landmark, , drop = FALSE]
			j <- match(d$eid, roles$eid)
			if (nrow(d) != a$N || anyDuplicated(d$eid) || anyNA(j) || !setequal(unique(d$fold), seq_len(K)) || any(d$fold !=
				roles$.fold[j]))
				stop("Incomplete or inconsistent saved Final cohort at landmark ", a$landmark)
		}
		cohorts <- cohorts[order(match(cohorts$landmark, availability$landmark), match(cohorts$eid, roles$eid)), ,
			drop = FALSE
		]
		class(cohorts) <- c("final_cached_cohorts", class(cohorts))
		cohorts
	}

	# Completed numerical results outlive the disposable fold checkpoints. Normal runs and explicit recovery use
	# the same check before loading raw omics.
	final_completed_joint_outputs <- function(root, trait, B, seed, primary_budget, primary_landmark, replace = FALSE) {
		if (isTRUE(replace))
			return(NULL)
		required <- file.path(root, c(
			"cohort.csv", "res.rds", "metrics.csv", "summary_provenance.csv",
			"_private/out_of_fold_predictions.csv.gz"
		))
		fi <- file.info(required)
		if (anyNA(fi$size) || any(fi$isdir) || any(fi$size <= 0))
			return(NULL)
		cohort <- tryCatch(read.csv(required[1], stringsAsFactors = FALSE), error = function(e) NULL)
		if (!is.data.frame(cohort) || nrow(cohort) != 1L || !all(c("Y", "signature") %in% names(cohort)) || !identical(
			cohort$Y,
			trait
		) || is.na(cohort$signature) || !nzchar(cohort$signature))
			return(NULL)
		final <- tryCatch(readRDS(required[2]), error = function(e) NULL)
		settings <- list(bootstrap = B, seed = seed, primary_budget = primary_budget, primary_landmark = primary_landmark)
		if (!is.list(final) || !isTRUE(final$summary_complete) || !identical(final$signature, cohort$signature) ||
			!isTRUE(all.equal(final$summary_settings, settings)) || !is.list(final$tables) || !is.data.frame(final$tables$metrics) ||
			!nrow(final$tables$metrics))
			return(NULL)
		if (!truthy(Sys.getenv('FINAL_JOINT_SUMMARY_ONLY', 'FALSE'))) le8_check_final_request(final)
		message("Final: reuse completed joint outputs; no raw-data loading or model retraining. ", "Use --replace TRUE or a new output directory to recompute after changing inputs/model settings.")
		final$tables
	}

	# Explicit recovery uses the frozen cache's cohort, folds and signature. It never relabels old fits using
	# newly loaded all.rds/omics or new training code.
	final_resume_joint_outputs <- function(root, trait, B, seed, primary_budget, primary_landmark) {
		completed <- final_completed_joint_outputs(root, trait, B, seed, primary_budget, primary_landmark)
		if (!is.null(completed))
			return(invisible(completed))
		required <- c(
			"cohort.csv", "landmark_availability.csv", "assay_inventory.csv", "factorial_design.csv",
			"prs_provenance.csv", "omic_PGS_status.txt", "_private/roles.csv.gz"
		)
		missing <- required[!file.exists(file.path(root, required))]
		if (length(missing))
			stop("Final summary-only recovery needs saved metadata: ", paste(missing, collapse = ", "))
		cohort <- read.csv(file.path(root, "cohort.csv"), stringsAsFactors = FALSE)
		if (nrow(cohort) != 1L || cohort$Y != trait || is.na(cohort$signature) || !nzchar(cohort$signature))
			stop("Invalid saved Final cohort provenance")
		K <- as.integer(cohort$outer_folds)
		if (!is.finite(K) || K < 3L)
			stop("Invalid saved Final fold count")
		a <- read.csv(file.path(root, "landmark_availability.csv"))
		a <- a[a$available %in% TRUE, , drop = FALSE]
		if (!nrow(a) || !primary_landmark %in% a$landmark)
			stop("Primary landmark unavailable in saved Final run")
		paths <- unlist(lapply(a$landmark, function(L) file.path(root, "_cache", paste0(
			"L", L, ".fold", seq_len(K),
			".rds"
		))), use.names = FALSE)
		if (any(!file.exists(paths)))
			stop("Final summary-only requires every completed fold; missing: ", paste(basename(paths[!file.exists(paths)]),
				collapse = ", "
			))
		roles <- as.data.frame(data.table::fread(file.path(root, "_private", "roles.csv.gz"), colClasses = c(eid = "character")))
		map <- read.csv(file.path(root, "assay_inventory.csv"), stringsAsFactors = FALSE)
		message("Final: summary-only recovery of ", length(paths), " frozen fold checkpoints; no model retraining")
		validation <- list(availability = a, roles = roles, K = K)
		final_joint_outputs(paths, NULL, map, root, K, B, seed, primary_budget, primary_landmark, cohort$signature,
			cached_validation = validation
		)
	}


	LE8_JOB <- "final_prediction"
	dir.create(le8_final_dir(), recursive = TRUE, showWarnings = FALSE)
	if (!requireNamespace("glmnet", quietly = TRUE)) stop("Final requires glmnet for cross-fitted proxy learning")
	risk_solver <- Sys.getenv("FINAL_RISK_SOLVER", "ridge")
	if (!risk_solver %in% c("cox", "ridge")) stop("FINAL_RISK_SOLVER must be cox or ridge")
	if (!prot_DO || !met_DO) stop("Joint Final requires --biom prot,met; all comparisons use their common cohort")
	if (nzchar(Sys.getenv("FINAL_INC_PROT", "")) || nzchar(Sys.getenv("FINAL_INC_MET", ""))) stop("Outcome-derived include lists cannot be reused for joint validation. Supply full assay matrices.")
	K <- as.integer(Sys.getenv("FINAL_OUTER_FOLDS", "5"))
	B <- as.integer(Sys.getenv("FINAL_BOOT", "200"))
	budgets <- as.numeric(final_csv("FINAL_ASSAY_BUDGETS", "5,10,50"))
	landmarks <- as.numeric(final_csv("FINAL_LANDMARKS", "0,2,5"))
	end <- as.numeric(Sys.getenv("FINAL_END_YEARS", "10"))
	if (!is.finite(K) || K < 3 || !is.finite(B) || B < 0 || any(!is.finite(budgets) | budgets < 2 | budgets > 100 |
		budgets != floor(budgets)) || any(!is.finite(landmarks) | landmarks < 0 | landmarks >= end)) stop("Invalid fold, budget or landmark settings")
	budgets <- sort(unique(budgets))
	landmarks <- sort(unique(landmarks))
	primary_budget <- as.integer(Sys.getenv("FINAL_PRIMARY_BUDGET", "10"))
	primary_landmark <- as.numeric(Sys.getenv("FINAL_PRIMARY_LANDMARK", "5"))
	if (!primary_budget %in% budgets || !primary_landmark %in% landmarks) stop("Prespecified primary budget/landmark must be included")
	if (truthy(Sys.getenv("FINAL_JOINT_SUMMARY_ONLY", "FALSE"))) {
		if (LE8_REPLACE)
			stop("FINAL_JOINT_SUMMARY_ONLY cannot be combined with --replace TRUE")
		final_resume_joint_outputs(final.out, Y, B, SEED, primary_budget, primary_landmark)
		quit(save = "no", status = 0)
	}
	completed <- final_completed_joint_outputs(final.out, Y, B, SEED, primary_budget, primary_landmark, replace = LE8_REPLACE)
	if (!is.null(completed)) quit(save = "no", status = 0)
	rm(completed)
	clinical <- final_csv("FINAL_CLINICAL_VARS", paste(unique(c(
		vars.basic, "smoke.pts", "sbp", "hba1c_ngsp", "bmi",
		"nonhdl"
	)), collapse = ","))
	basic <- vars.basic
	targets <- final_csv("FINAL_DOMAINS", "bmi,nonhdl,hba1c_ngsp,sbp,smoke.pts,diet.pts,pa.pts,sleep.pts")
	proxy_targets <- final_csv("FINAL_PROXY_TARGETS", "bmi,nonhdl")
	if (length(targets) < 1 || length(targets) > 8 || length(proxy_targets) < 1 || length(proxy_targets) > 8) stop("Choose 1–8 candidate domains / proxy targets")
	if (primary_budget < length(proxy_targets)) stop("Primary assay budget must cover requested proxy targets")
	if (!all(c("bmi", "nonhdl") %in% clinical)) stop("Continuous BMI and non-HDL are required clinical comparators")
	if (any(grepl("(^fod_|^date_|[.]Yt2e$|[.]t2e$|[.]b2e$|^drug[.]|^dm[.]|^htn[.])", clinical))) stop("Unapproved diagnosis/medication/outcome-derived predictor; use explicit baseline fields")
	if (any(grepl("[.](pgs|prs)$|^disease_PRS$", c(clinical, basic, targets, proxy_targets), ignore.case = TRUE))) stop("Keep disease PRS / biomarker PGS out of the clinical and LE8-domain comparator fields")
	groupcol <- Sys.getenv("FINAL_GROUP_COLUMN", "")
	ph <- read_all(unique(c(
		"eid", "ethnic.c", clinical, basic, targets, proxy_targets, "birth_date", "date_attend",
		"date_lost", "date_death", le8_y_date(), groupcol, paste0(Y, ".pgs")
	))) |>
		filter_analysis_cohort() |>
		make_outcome(Y)
	if (length(setdiff(c(clinical, basic), names(ph)))) stop("Missing clinical baseline fields")
	prs <- final_prs_input(Y, ph)
	has_prs <- isTRUE(prs$enabled)
	prs$manifest$file[prs$manifest$file == "all.rds"] <- file.path(indir, "Rdata/all.rds")
	# First count the common measured cohort. Never impute a missing disease PRS.
	joined <- final_join_omics(ph, read_prot(), read_met(), data.frame(eid = as.character(ph$eid)))
	prs_flow <- data.frame(stage = "clinical_prot_met_intersection", N = nrow(joined$data))
	if (has_prs) {
		joined$data <- merge(joined$data, prs$data, by = "eid", all.x = TRUE, sort = FALSE)
		prs_flow <- rbind(prs_flow, data.frame(stage = "finite_endpoint_PRS_in_intersection", N = sum(is.finite(joined$data$disease_PRS))))
	}
	pgs <- final_load_biomarker_pgs(joined$data, joined$map)
	dat <- pgs$data
	map <- pgs$map
	pgs_inputs <- pgs$inputs
	pgs_status <- pgs$status
	rm(ph, joined, pgs)
	invisible(gc())
	dat$.time <- dat[[paste0(Y, ".t2e")]]
	dat$.event <- dat[[paste0(Y, ".Yt2e")]]
	dat <- final_endpoint_coverage(dat, Y)
	dat$.baseline_disease <- is.finite(dat[[paste0(Y, ".b2e")]]) & dat[[paste0(Y, ".b2e")]] <= 0
	if (nzchar(groupcol) && !groupcol %in% names(dat)) stop("Group column missing")
	dat$.group <- if (nzchar(groupcol)) as.character(dat[[groupcol]]) else dat$eid
	yang <- dat[dat$.baseline_disease, , drop = FALSE]
	yin <- dat[!dat$.baseline_disease & is.finite(dat$.time) & dat$.time > dat$.entry & dat$.event %in% 0 : 1, , drop = FALSE]
	prs_flow <- rbind(prs_flow, data.frame(stage = "incident_measured_cohort_before_PRS", N = nrow(yin)))
	if (has_prs) {
		eligible <- is.finite(yin$disease_PRS)
		if (sum(eligible) >= 1000 && sum(yin$.event[eligible]) >= 100 && sd(yin$disease_PRS[eligible]) > 0) {
			yin <- yin[eligible, , drop = FALSE]
			yang <- yang[is.finite(yang$disease_PRS), , drop = FALSE]
			prs_flow <- rbind(prs_flow, data.frame(
				stage = "common_incident_cohort_for_ALL_models_with_and_without_PRS",
				N = nrow(yin)
			))
		} else {
			has_prs <- FALSE
			prs$enabled <- FALSE
			prs$manifest$enabled <- FALSE
			prs$manifest$status <- paste(prs$manifest$status, "; insufficient common incident PRS cohort, non-PRS suite continues")
		}
	}
	write.csv(prs_flow, file.path(final.out, "prs_cohort_flow.csv"), row.names = FALSE)
	if (length(intersect(yin$eid, yang$eid))) stop("Yin/Yang overlap")
	if (nrow(yin) < 1000 || sum(yin$.event) < 100) stop("Insufficient common incident cohort; do not substitute a different cohort per modality")
	yin$.fold <- final_group_folds(yin$.group, K, SEED)
	private <- file.path(final.out, "_private")
	dir.create(private, recursive = TRUE, showWarnings = FALSE)
	cache <- file.path(final.out, "_cache")
	dir.create(cache, recursive = TRUE, showWarnings = FALSE)
	writeLines(c("*", "!.gitignore"), file.path(private, ".gitignore"))
	writeLines(c("*", "!.gitignore"), file.path(cache, ".gitignore"))
	input_files <- c(
		file.path(indir, "Rdata", c("all.rds", "prot.rds", "met.rds")), prs$manifest$file, Sys.getenv("FINAL_PRS_MANIFEST"),
		Sys.getenv("LE8_ENDPOINT_MANIFEST"), pgs_inputs
	)
	fi <- file.info(input_files)
	signature <- le8_hash_object(list(
		version = "2026-09-16.factorial-prs", Y = Y, inputs = data.frame(
			path = input_files,
			size = fi$size, mtime = as.character(fi$mtime)
		), code = tools::md5sum(file.path(fdir, c(
			"final.R", "c1.correlate.R",
			"0.common.R"
		))), options = le8_analysis_options(), settings = Sys.getenv()[grepl(
			"^(FINAL_|DATE_FOLLOW_END)",
			names(Sys.getenv())
		)], ids = yin$eid, fold = yin$.fold, clinical = clinical, budgets = budgets, landmarks = landmarks,
		end = end, seed = SEED
	))
	write.csv(prs$manifest, file.path(final.out, "prs_provenance.csv"), row.names = FALSE)
	factorial_design <- final_factorial_design(has_prs)
	write.csv(factorial_design, file.path(final.out, "factorial_design.csv"), row.names = FALSE)
	write.csv(
		data.frame(
			Y = Y, N_joint = nrow(dat), N_Yin = nrow(yin), N_Yang = nrow(yang), incident_events = sum(yin$.event),
			outer_folds = K, grouping = if (nzchar(groupcol)) groupcol else "eid; relatedness not controlled", signature = signature
		),
		file.path(final.out, "cohort.csv"),
		row.names = FALSE
	)
	data.table::fwrite(yin[, c("eid", ".fold", ".group")], file.path(private, "roles.csv.gz"), compress = "gzip")
	write.csv(map, file.path(final.out, "assay_inventory.csv"), row.names = FALSE)
	writeLines(pgs_status, file.path(final.out, "omic_PGS_status.txt"))
	# yin/yang already own the needed rows; the full joined matrix is no longer used.
	rm(dat)
	invisible(gc())
	drop_assays <- toupper(final_csv("FINAL_EXCLUDE_ASSAYS", "GDF15,NTPROBNP,NPPB"))
	features <- map$feature
	arms <- list(all_assays = features, omit_GDF15_natriuretic = map$feature[!toupper(map$assay) %in% drop_assays])
	if (truthy(Sys.getenv("FINAL_INCLUDE_LIPID_ABLATION", "TRUE"))) {
		direct_lipid <- toupper(final_csv("FINAL_DIRECT_LIPID_ASSAYS", "Non_HDL_C,NonHDL_C,Total_C,HDL_C"))
		arms$omit_direct_lipid_reconstruction <- map$feature[!(map$layer == "met" & toupper(map$assay) %in% direct_lipid)]
	}
	write.csv(data.frame(arm = names(arms), available_assays = lengths(arms)), file.path(final.out, "ablation_inventory.csv"),
		row.names = FALSE
	)
	all_results <- character()
	landmark_availability <- map_dfr(landmarks, function(L) {
		d <- final_landmark(yin, L, end)
		ne <- vapply(seq_len(K), function(fd) sum(d$.fold == fd & d$.event == 1 & d$.time <= end - L), numeric(1))
		ok <- nrow(d) >= 1000 && all(ne >= 10) && all(sum(ne) - ne >= 50)
		tibble(landmark = L, N = nrow(d), events = sum(ne), min_test_events = min(ne), available = ok, status = if (ok)
			"estimable" else "insufficient covered event-free participants/events; no substitute landmark chosen")
	})
	write.csv(landmark_availability, file.path(final.out, "landmark_availability.csv"), row.names = FALSE)
	landmarks_to_run <- landmark_availability$landmark[landmark_availability$available]
	if (!length(landmarks_to_run)) stop("No prespecified landmark is estimable; see landmark_availability.csv")
	for (L in landmarks_to_run) for (fd in seq_len(K)) {
		checkpoint <- file.path(cache, paste0("L", L, ".fold", fd, ".rds"))
		old <- if (file.exists(checkpoint) && !LE8_REPLACE)
			tryCatch(readRDS(checkpoint), error = function(e) NULL) else NULL
		reusable <- identical(old$signature, signature) && is.list(old$result)
		if (reusable)
			final_prune_fold_models(checkpoint, signature, object = old)
		rm(old)
		invisible(gc())
		if (reusable) {
			message("Final: reuse landmark ", L, ", outer fold ", fd, "/", K)
			all_results <- c(all_results, checkpoint)
			next
		}
		tr <- final_landmark(yin[yin$.fold != fd, , drop = FALSE], L, end)
		te <- final_landmark(yin[yin$.fold == fd, , drop = FALSE], L, end)
		tr$.event <- as.integer(tr$.event == 1 & tr$.time <= end - L)
		tr$.time <- pmin(tr$.time, end - L)
		yy <- yang[!yang$.group %in% yin$.group[yin$.fold == fd], , drop = FALSE]
		if (length(intersect(tr$.group, te$.group)) || length(intersect(yy$.group, te$.group)))
			stop("Group leakage")
		if (sum(tr$.event) < 50 || sum(te$.event) < 10)
			stop("Too few events at landmark/fold; reduce prespecified model complexity or obtain more data")
		message(Y, ": joint Final landmark ", L, ", outer fold ", fd, "/", K, "; events ", sum(tr$.event), " / ", sum(te$.event))
		# Identical screening reference for all panels; screening fits use outer train only.
		screen <- cox_scan(tr, features, clinical, Y, time_var = ".time", event_var = ".event")
		ranked <- screen |>
			filter(is.finite(p.value)) |>
			arrange(p.value, term) |>
			pull(term)
		genetic_screen <- tibble()
		genetic_ranked <- character()
		if ("pgs_feature" %in% names(map)) {
			gs <- map$pgs_feature[!is.na(map$pgs_feature)]
			if (length(gs)) {
				genetic_screen <- cox_scan(tr, gs, clinical, Y, time_var = ".time", event_var = ".event")
				genetic_ranked <- genetic_screen |>
					filter(is.finite(p.value)) |>
					arrange(p.value, term) |>
					pull(term)
			}
		}
		preds <- members <- diagnostics <- coefficients <- proxy_validation <- proxy_weights <- proxy_preprocessing <- domain_tables <- list()
		risk_preprocessing <- baseline_hazards <- list()
		pgs_reconstruction <- pgs_coefficients <- cross_omic <- domain_status <- list()
		marker <- final_marker_strata(tr, te, map)
		inflam <- final_inflammation_score(tr, te, map)
		fit_cache <- new.env(parent = emptyenv())
		addfit <- function(model, k, arm, cv, ff, train = tr, test = te) {
			if (!has_prs && "disease_PRS" %in% c(cv, ff))
				return(invisible(NULL))
			cacheable <- missing(train) && missing(test)
			ck <- paste(paste(cv, collapse = "|"), paste(ff, collapse = "|"), sep = "::")
			obj <- if (cacheable && exists(ck, fit_cache, inherits = FALSE))
				get(ck, fit_cache) else le8_fit_budget_model(train, test, cv, ff, ".time", ".event", solver = risk_solver, seed = SEED + fd)
			if (cacheable)
				assign(ck, obj, fit_cache)
			transformed <- ff[grepl("^(genpart__|remainder__)", ff)]
			inherited_dependencies <- if (length(transformed))
				paste0("pgs__", sub("^(genpart__|remainder__)", "", transformed)) else character()
			measured_dependencies <- unique(c(ff[ff %in% map$feature], sub("^remainder__", "", ff[startsWith(ff, "remainder__")])))
			id <- length(diagnostics) + 1L
			diagnostics[[id]] <<- tibble(model,
				budget = k, arm, fold = fd, landmark = L, status = obj$status, requested_predictors = length(ff),
				fitted_predictors = obj$N_selected %||% NA_integer_, direct_measured_assays = length(measured_dependencies),
				disease_PRS_predictors = as.integer("disease_PRS" %in% c(cv, ff)), genetic_score_predictors = length(unique(c(ff[startsWith(
					ff,
					"pgs__"
				)], inherited_dependencies)))
			)
			if (length(ff))
				members[[length(members) + 1L]] <<- tibble(model, budget = k, arm, fold = fd, landmark = L, feature = ff)
			if ("disease_PRS" %in% cv)
				members[[length(members) + 1L]] <<- tibble(model, budget = k, arm, fold = fd, landmark = L, feature = "disease_PRS")
			dependencies <- sub("^remainder__", "", ff[startsWith(ff, "remainder__")])
			if (length(dependencies))
				members[[length(members) + 1L]] <<- tibble(model, budget = k, arm, fold = fd, landmark = L, feature = dependencies)
			if (obj$status != "ok")
				return(invisible(NULL))
			risk <- le8_risk_at(obj, end - L)
			tx <- le8_prepare_prediction_matrix(train, train, unique(c(cv, ff)))
			lp_training <- drop(tx$train %*% obj$coefficient$beta) - obj$lp_center
			q75 <- unname(quantile(lp_training, 0.75))
			preds[[length(preds) + 1L]] <<- tibble(
				eid = test$eid, group = test$.group, fold = fd, time = test$.time,
				event = test$.event, model, budget = k, arm, landmark = L, horizon = end - L, risk = risk, lp = obj$lp,
				high_training_q75 = obj$lp >= q75, inflammatory_burden = if (is.null(inflam))
					"unavailable" else inflam$group
			)
			coefficients[[length(coefficients) + 1L]] <<- obj$coefficient |>
				mutate(model, budget = k, arm, fold = fd, landmark = L, lambda = obj$lambda, condition_number = obj$condition_number)
			risk_preprocessing[[length(risk_preprocessing) + 1L]] <<- obj$preprocess |>
				mutate(model, budget = k, arm, fold = fd, landmark = L)
			bh <- obj$baseline_hazard
			baseline_hazards[[length(baseline_hazards) + 1L]] <<- bh[!duplicated(bh$hazard), , drop = FALSE] |>
				mutate(model, budget = k, arm, fold = fd, landmark = L, lp_center = obj$lp_center)
			# Fold tables already retain these values; never persist per-model RDS copies.
			invisible(obj)
		}
		for (arm in names(arms)) {
			pool <- arms[[arm]]
			rr <- ranked[ranked %in% pool]
			addfit("Clinical", 0, arm, clinical, character())
			addfit("PRS_only", 0, arm, "disease_PRS", character())
			addfit("Clinical_PRS", 0, arm, c(clinical, "disease_PRS"), character())
			for (i in which(!factorial_design$P & !factorial_design$M & factorial_design$available)) {
				ds <- factorial_design[i, ]
				cv <- c(if (ds$C) clinical, if (ds$G) "disease_PRS")
				addfit(ds$model, 0, arm, cv, character())
			}
			domains <- list()
			for (cohort in c("Yin", "YinYang")) {
				learn <- if (cohort == "Yin")
					tr else bind_rows(tr, yy)
				dom <- final_replicated_domains(learn, pool, targets, basic, SEED + fd + 411)
				domains[[cohort]] <- dom$ranks
				domain_tables[[length(domain_tables) + 1L]] <- dom$all |>
					mutate(cohort, arm, fold = fd, landmark = L)
				domain_status[[length(domain_status) + 1L]] <- dom$audit |>
					mutate(cohort, arm, fold = fd, landmark = L)
			}
			for (k in budgets) {
				pp <- head(rr[startsWith(rr, "prot__")], k)
				mm <- head(rr[startsWith(rr, "met__")], k)
				# Joint NS budget is TOTAL assays: ceil(k/2) protein, floor(k/2) metabolite.
				joint <- c(head(pp, ceiling(k / 2)), head(mm, floor(k / 2)))
				# Fixed layer panels: adding M retains ALL k proteins (k+k assays).  The legacy equal-total-budget
				# comparison below remains unchanged.
				for (i in which((factorial_design$P | factorial_design$M) & factorial_design$available)) {
					ds <- factorial_design[i, ]
					ff <- c(if (ds$P) pp, if (ds$M) mm)
					cv <- c(if (ds$C) clinical, if (ds$G) "disease_PRS")
					if ((ds$P && length(pp) != k) || (ds$M && length(mm) != k))
						next
					addfit(ds$model, k, arm, cv, ff)
				}
				for (mode in c("Protein", "Metabolite", "ProtMet")) {
					ff <- switch(mode,
						Protein = pp,
						Metabolite = mm,
						ProtMet = joint
					)
					if (length(ff) != k)
						stop("Insufficient measured assays for requested budget")
					addfit(paste0(mode, "_only_NS"), k, arm, character(), ff)
					addfit(paste0("Clinical_", mode, "_NS"), k, arm, clinical, ff)
				}
				addfit("Clinical_Protein_sharedBudget_NS", k, arm, clinical, joint[startsWith(joint, "prot__")])
				addfit("Clinical_Metabolite_sharedBudget_NS", k, arm, clinical, joint[startsWith(joint, "met__")])
				addfit("Clinical_ProtMet_PRS_NS", k, arm, c(clinical, "disease_PRS"), joint)
				gp <- character()
				if (length(genetic_ranked)) {
					eligible <- map$pgs_feature[map$feature %in% pool & !is.na(map$pgs_feature)]
					gr <- genetic_ranked[genetic_ranked %in% eligible]
					gp <- c(head(gr[startsWith(gr, "pgs__prot__")], ceiling(k / 2)), head(
						gr[startsWith(gr, "pgs__met__")],
						floor(k / 2)
					))
					if (length(gp) == k) {
						addfit("Clinical_PGSselected_NS", k, arm, clinical, gp)
						addfit("Clinical_ProtMet_PGSselected_NS", k, arm, clinical, c(joint, gp))
						addfit("Clinical_ProtMet_PGSselected_PRS_NS", k, arm, c(clinical, "disease_PRS"), c(joint, gp))
						genetically_ranked_assays <- map$feature[match(gp, map$pgs_feature)]
						addfit("Clinical_GeneticSelectedMeasured_NS", k, arm, clinical, genetically_ranked_assays)
					}
				}
				if ("pgs_feature" %in% names(map)) {
					gg <- map$pgs_feature[match(joint, map$feature)]
					valid <- vapply(
						gg, function(g) !is.na(g) && g %in% names(tr) && sum(is.finite(tr[[g]])) >= 100,
						logical(1)
					)
					matched <- joint[valid]
					gg <- gg[valid]
					if (length(gg) >= 2) {
						# The reduced matched subset gets its OWN measured reference.
						addfit("MatchedPGS_only_NS", k, arm, character(), gg)
						addfit("MatchedMeasured_only_NS", k, arm, character(), matched)
						addfit("Clinical_MatchedMeasured_NS", k, arm, clinical, matched)
						addfit("Clinical_MatchedPGS_NS", k, arm, clinical, gg)
						addfit("Clinical_MatchedMeasured_PGS_NS", k, arm, clinical, c(matched, gg))
						addfit("Clinical_MatchedMeasured_PGS_PRS_NS", k, arm, c(clinical, "disease_PRS"), c(
							matched,
							gg
						))
						if (k == primary_budget) {
							decomp <- final_pgs_decompose(tr, te, matched, map, SEED + fd + 312)
							pgs_reconstruction[[length(pgs_reconstruction) + 1L]] <- decomp$audit |>
								mutate(fold = fd, landmark = L, arm)
							pgs_coefficients[[length(pgs_coefficients) + 1L]] <- decomp$coefficient |>
								mutate(fold = fd, landmark = L, arm)
							if (length(decomp$genetic) >= 2) {
								addfit("Clinical_PGS_captured", k, arm, clinical, decomp$genetic, decomp$train, decomp$test)
								addfit("Clinical_PGS_remainder", k, arm, clinical, decomp$remainder, decomp$train, decomp$test)
								addfit(
									"Clinical_PGS_captured_remainder", k, arm, clinical, c(decomp$genetic, decomp$remainder),
									decomp$train, decomp$test
								)
							}
						}
					}
				}
				for (cohort in c("Yin", "YinYang")) {
					queues <- lapply(domains[[cohort]], function(z) z$feature)
					ys <- final_ys_panel(queues, k)
					plus <- final_ys_panel(queues, k, rr)
					if (length(ys) != k || length(plus) != k) {
						diagnostics[[length(diagnostics) + 1L]] <- tibble(
							model = paste0("YS_", cohort), budget = k,
							arm, fold = fd, landmark = L, status = "insufficient replicated domain proxies; no forced pillar or NS substitution"
						)
						next
					}
					addfit(paste0("Clinical_ProtMet_YS_", cohort), k, arm, clinical, ys)
					addfit(paste0("Clinical_ProtMet_YSplus_", cohort), k, arm, clinical, plus)
					addfit(paste0("Clinical_ProtMet_PRS_YS_", cohort), k, arm, c(clinical, "disease_PRS"), ys)
					if (length(gp) == k)
						addfit(paste0("Clinical_ProtMet_YS_PGSselected_", cohort), k, arm, clinical, c(ys, gp))
					if (k == primary_budget && cohort == "YinYang")
						cross_omic[[length(cross_omic) + 1L]] <- final_cross_omic_links(tr, te, union(joint, ys), basic) |>
							mutate(fold = fd, landmark = L, arm)
				}
			}
			# Domain proxy models at one prespecified budget. Report union assay count; two predicted targets do not
			# mean a two-assay clinical panel.
			pk <- primary_budget
			for (cohort in c("Yin", "YinYang")) {
				extra <- if (cohort == "Yin")
					yy[FALSE, , drop = FALSE] else yy
				for (modality in c("prot", "met", "joint")) {
					candidates <- if (modality == "joint")
						pool else pool[startsWith(pool, paste0(modality, "__"))]
					trp <- tr
					tep <- te
					union_panel <- character()
					proxy_ok <- TRUE
					for (tg in proxy_targets) {
						kp <- pk %/% length(proxy_targets) + as.integer(match(tg, proxy_targets) <= pk %% length(proxy_targets))
						obj <- tryCatch(final_crossfit_proxy(tr, te, extra, candidates, tg, kp, basic, SEED + fd + 111),
							error = function(e) {
								diagnostics[[length(diagnostics) + 1L]] <<- tibble(model = paste0(
									"Proxy_", tg, "_", modality,
									"_", cohort
								), budget = pk, arm, fold = fd, landmark = L, status = conditionMessage(e))
								NULL
							}
						)
						if (is.null(obj)) {
							proxy_ok <- FALSE
							next
						}
						nv <- paste0("proxy_", tg)
						trp[[nv]] <- obj$train
						tep[[nv]] <- obj$test
						union_panel <- union(union_panel, obj$final$features)
						ok <- is.finite(te[[tg]])
						pred <- obj$test[ok]
						obs <- te[[tg]][ok]
						denom <- sum((obs - mean(tr[[tg]], na.rm = TRUE)) ^ 2)
						proxy_validation[[length(proxy_validation) + 1L]] <- tibble(
							fold = fd, landmark = L, arm, cohort,
							modality, target = tg, N = sum(ok), R2_vs_training_mean = 1 - sum((obs - pred) ^ 2) / denom, RMSE = sqrt(mean((obs -
								pred) ^ 2)), correlation = cor(obs, pred), scope = "Held-out reconstruction of continuous measurement; not evidence of intervention responsiveness"
						)
						proxy_weights[[length(proxy_weights) + 1L]] <- obj$final$coefficient |>
							mutate(fold = fd, landmark = L, arm, cohort, modality, target = tg)
						proxy_preprocessing[[length(proxy_preprocessing) + 1L]] <- obj$final$preprocess |>
							mutate(fold = fd, landmark = L, arm, cohort, modality, target = tg)
					}
					if (!proxy_ok)
						next
					for (action in c("Replace", "Add")) {
						if (action == "Replace" && !all(proxy_targets %in% clinical))
							next
						cv <- if (action == "Replace")
							setdiff(clinical, proxy_targets) else clinical
						label <- paste0("Clinical_", action, "Proxy_", modality, "_", cohort)
						addfit(label, pk, arm, cv, paste0("proxy_", proxy_targets), trp, tep)
						# Correct the feature inventory: measured inputs, not two score names.
						members[[length(members) + 1L]] <- tibble(
							model = label, budget = pk, arm, fold = fd, landmark = L,
							feature = union_panel
						)
					}
				}
			}
		}
		result <- list(
			predictions = bind_rows(preds), members = bind_rows(members), diagnostics = bind_rows(diagnostics),
			risk_preprocessing = bind_rows(risk_preprocessing), baseline_hazards = bind_rows(baseline_hazards), coefficients = bind_rows(coefficients),
			proxy_validation = bind_rows(proxy_validation), proxy_weights = bind_rows(proxy_weights), proxy_preprocessing = bind_rows(proxy_preprocessing),
			domain_associations = bind_rows(domain_tables), domain_status = bind_rows(domain_status), pgs_reconstruction = bind_rows(pgs_reconstruction),
			pgs_coefficients = bind_rows(pgs_coefficients), cross_omic = bind_rows(cross_omic), marker_groups = marker$groups |>
				mutate(fold = fd, landmark = L), marker_audit = marker$audit |>
				mutate(fold = fd, landmark = L), genetic_screen = genetic_screen |>
				mutate(fold = fd, landmark = L), screen = screen |>
				mutate(fold = fd, landmark = L)
		)
		final_write_joint_checkpoint(list(signature = signature, result = result), checkpoint)
		final_prune_fold_models(checkpoint, signature, object = list(signature = signature, result = result))
		all_results <- c(all_results, checkpoint)
		rm(
			result, preds, members, diagnostics, coefficients, proxy_validation, proxy_weights, proxy_preprocessing,
			domain_tables, pgs_reconstruction, pgs_coefficients, cross_omic, domain_status, fit_cache, risk_preprocessing,
			baseline_hazards
		)
		rm(list = intersect(c(
			"tr", "te", "yy", "learn", "trp", "tep", "decomp", "obj", "extra", "marker", "dom", "domains",
			"screen", "genetic_screen"
		), ls()))
		invisible(gc())
	}
	# Summaries need only eligibility/IDs; release the assay matrices and final fold's working copies before
	# reading any saved predictions.
	yin <- yin[, intersect(c("eid", ".time", ".event", ".entry"), names(yin)), drop = FALSE]
	rm(list = intersect(c(
		"dat", "yang", "tr", "te", "yy", "learn", "trp", "tep", "decomp", "obj", "extra", "marker",
		"dom", "domains", "screen", "genetic_screen"
	), ls()))
	invisible(gc())
	final_joint_outputs(all_results, yin, map, final.out, K, B, SEED, primary_budget, primary_landmark, signature)
} else if (.final_mode == "pipeline") {
	# Final orchestrator: retain per-layer analyses, add joint validation, and always consolidate available
	# C1–Final evidence at analysis/le8/final/[Y]/.
	fdir <- Sys.getenv("LE8_FDIR", file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(fdir, "0.common.R"))
	LE8_JOB <- "final_prediction"
	dir.create(le8_final_dir(), recursive = TRUE, showWarnings = FALSE)
	stage_log <- list()
	run_stage <- function(name, mode, env = character()) {
		started <- le8_stage_start(paste0("Final/", name))
		logfile <- file.path(le8_final_dir(), paste0(name, ".log"))
		status <- tryCatch(system2(file.path(R.home("bin"), "Rscript"), shQuote(file.path(fdir, "final.R")),
			stdout = logfile,
			stderr = logfile, env = c(paste0("LE8_FINAL_MODE=", mode), env)
		), error = function(e) {
			writeLines(conditionMessage(e), logfile)
			1L
		})
		if (status == 0)
			le8_stage_done(paste0("Final/", name), started) else message("[LE8] FAIL Final/", name, " exit=", status, " log=", logfile)
		stage_log[[length(stage_log) + 1L]] <<- data.frame(stage = name, status = if (status == 0)
			"complete" else "failed", detail = paste("Exit", status, ";", basename(logfile)))
	}
	skip <- function(stage, detail) stage_log[[length(stage_log) + 1L]] <<- data.frame(stage,
		status = "unavailable",
		detail
	)
	if (truthy(Sys.getenv("FINAL_RUN_REFERENCE", "TRUE"))) {
		for (layer in c(if (prot_DO) "prot", if (met_DO) "met")) {
			d <- if (layer == "prot")
				out.prot else out.met
			if (file.exists(file.path(le8_job_dir(d, "c1_correlate"), "c1.res.rds")))
				run_stage(paste0("reference_", layer), "reference", env = c(paste0(
					"BIOM=",
					layer
				), paste0("PROT_DO=", if (layer == "prot") "TRUE" else "FALSE"), paste0("MET_DO=", if (layer ==
					"met") "TRUE" else "FALSE"))) else skip(paste0("reference_", layer), "C1 RDS unavailable; other layers/stages continue")
		}
	}
	if (truthy(Sys.getenv("FINAL_RUN_JOINT", "TRUE"))) {
		if (prot_DO && met_DO)
			run_stage("joint", "joint") else skip("joint", "Select --biom prot,met for a common-cohort joint comparison; per-layer modules remain available")
	}
	code <- file.path(fdir, "final.py")
	atlas_started <- le8_stage_start("Final/systematic_atlas")
	status <- tryCatch(system2(Sys.getenv("PYTHON_BIN", "python3"), vapply(c(code, "tables", "--root", out.base), shQuote, character(1))),
		error = function(e) 1L
	)
	if (status == 0) le8_stage_done("Final/systematic_atlas", atlas_started) else message(
		"[LE8] FAIL Final/systematic_atlas exit=",
		status
	)
	stage_log[[length(stage_log) + 1L]] <- data.frame(
		stage = "systematic_atlas", status = if (status == 0) "complete" else "failed",
		detail = paste("Exit", status)
	)
	log <- bind_rows(stage_log)
	write.csv(log, file.path(le8_final_dir(), "pipeline_status.csv"), row.names = FALSE)
	if (any(log$status == "failed")) stop("One or more Final stages failed. Completed stages and upstream analyses are preserved; see pipeline_status.csv")
} else if (.final_mode == "report") {
	# Publication assembly from aggregate module outputs only. No participant data, model fitting, validation-set
	# selection, or internet access is required.
	suppressPackageStartupMessages({
		library(data.table)
		library(dplyr)
		library(tidyr)
		library(ggplot2)
		library(patchwork)
		library(openxlsx)
		library(purrr)
	})
	`%||%` <- function(x, y) if (is.null(x)) y else x
	source(file.path(Sys.getenv("LE8_FDIR", unset = "."), "0.common.R"))
	source(file.path(Sys.getenv("LE8_FDIR", unset = "."), "c1.correlate.R"))

	# final panels All plot data are exported with figure/panel provenance. Missing panels are omitted, and each
	# page contains no more than six actual axes.
	pub_theme <- function() theme_classic(base_size = 11) + theme(
		plot.title = element_text(face = "bold", size = 12),
		plot.subtitle = element_text(size = 9, color = "#526172"), axis.text = element_text(color = "#26364A"), legend.position = "bottom",
		legend.title = element_blank(), panel.grid.major.y = element_line(color = "#E9EDF1", linewidth = 0.25), plot.margin = margin(
			9,
			12, 9, 9
		)
	)
	pub_cols <- c("#286B8B", "#BC5366", "#6D629A", "#C58D37", "#388B78", "#687784", "#9C674E", "#6493B0")
	pub_ok <- function(d, cols) is.data.frame(d) && nrow(d) > 0 && all(cols %in% names(d))
	pub_num <- function(x) suppressWarnings(as.numeric(x))
	pub_forest <- function(d, label, beta, lo, hi, title, xlab = "Effect estimate", null = 0, n = 12) {
		if (!pub_ok(d, c(label, beta, lo, hi)))
			return(NULL)
		d <- d[is.finite(d[[beta]]) & is.finite(d[[lo]]) & is.finite(d[[hi]]), , drop = FALSE]
		d <- head(d, n)
		if (!nrow(d))
			return(NULL)
		d$.label <- factor(d[[label]], levels = rev(unique(d[[label]])))
		ggplot(d, aes(x = .data[[beta]], y = .label)) +
			geom_vline(xintercept = null, color = "grey65", linetype = 2) +
			geom_segment(aes(x = .data[[lo]], xend = .data[[hi]], yend = .label), color = pub_cols[1]) +
			geom_point(
				color = pub_cols[1],
				size = 2.2
			) +
			labs(title = title, x = xlab, y = NULL) +
			pub_theme()
	}
	pub_count <- function(d, col, title) {
		if (!pub_ok(d, col))
			return(NULL)
		z <- d |>
			count(.data[[col]], name = "n")
		z$.label <- stringr::str_wrap(as.character(z[[col]]), 32)
		ggplot(z, aes(n, reorder(.label, n))) +
			geom_col(fill = pub_cols[1], width = 0.68) +
			geom_text(aes(label = n),
				hjust =  - 0.15, size = 3
			) +
			scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
			labs(
				title = title,
				x = "Biomarkers", y = NULL
			) +
			pub_theme()
	}
	pub_enrich <- function(d, title, asset_dir) {
		if (!pub_ok(d, c("term_name", "adjusted_p", "source")))
			return(NULL)
		d <- d |>
			filter(is.finite(adjusted_p), adjusted_p < 0.05) |>
			arrange(adjusted_p) |>
			group_by(source) |>
			slice_head(n = 3) |>
			ungroup() |>
			arrange(adjusted_p) |>
			slice_head(n = 10)
		if (!nrow(d))
			return(NULL)
		f <- file.path(asset_dir, "go_terms.tsv")
		if (file.exists(f)) {
			labs <- fread(f)
			i <- match(d$term_name, labs$term)
			hit <- !is.na(i)
			d$term_name[hit] <- labs$TERM[i[hit]]
		}
		d$.label <- factor(stringr::str_wrap(d$term_name, 40), levels = rev(unique(stringr::str_wrap(d$term_name, 40))))
		ggplot(d, aes( - log10(pmax(adjusted_p, 1e-300)), .label, color = source)) +
			geom_point(size = 3) +
			scale_color_manual(values = pub_cols) +
			labs(title = title, x = expression( - log[10](FDR)), y = NULL) +
			pub_theme()
	}
	publication_run <- function(trait, layer, analysis_root) {
		stopifnot(layer %in% c("prot", "met"))
		source_root <- file.path(analysis_root, trait, layer)
		root <- le8_final_dir(layer, "report", trait, analysis_root)
		if (!dir.exists(source_root)) {
			message("FINAL: no analysis directory: ", source_root)
			return(invisible(NULL))
		}
		modules <- c(
			c1 = "c1_correlate", c2 = "c2_cause", c3 = "c3_coloc", c4 = "c4_connect", focus = "c4_connect",
			prediction = "final_prediction", extra = "final_supp", cell = "c5_cellulation"
		)
		# Regenerate this report at the same location; remove stale prior pages.
		if (dir.exists(root)) unlink(root, recursive = TRUE)
		dir.create(root, recursive = TRUE, showWarnings = FALSE)
		module_dirs <- setNames(file.path(source_root, modules), names(modules))
		module_dirs["prediction"] <- le8_final_dir(layer, "prediction", trait, analysis_root)
		module_dirs["extra"] <- file.path(root, "supplement")
		audit <- list()
		tables <- list()
		provenance <- list()
		figures <- list()
		captions <- character()
		read <- function(key, module, file) {
			if (layer == "met")
				file <- sub("^pwas_", "mwas_", file)
			path <- file.path(module_dirs[[module]], file)
			d <- if (file.exists(path))
				tryCatch(as_tibble(fread(path, showProgress = FALSE)), error = function(e) tibble()) else tibble()
			audit[[key]] <<- tibble(table = key, source = path, rows = nrow(d), bytes = if (file.exists(path))
				file.info(path)$size else NA_real_, modified = if (file.exists(path))
				as.character(file.info(path)$mtime) else NA_character_, status = if (nrow(d))
				"available" else "absent or no rows")
			tables[[key]] <<- d
			d
		}
		cohort <- read("cohort", "c1", "c1.cohort.csv")
		a <- read("incident", "c1", "pwas_incident_adj2.csv")
		prev <- read("prevalent", "c1", "pwas_prevalent_adj2.csv")
		temporal <- read("temporal", "c1", "c1.directionality_triage.csv")
		pgs <- read("pgs", "c1", "c1.pgs_actual_concordance.csv")
		paired <- read("paired", "c1", "c1.paired_pgs_measured.csv")
		enrich <- read("enrich", "c1", "c1.enrichment_incident_sig.csv")
		enrich_prev <- read("enrich_prev", "c1", "c1.enrichment_prevalent_sig.csv")
		mr <- read("mr", "c2", "c2.MR_all.csv")
		grades <- read("grades", "c2", "c2.evidence_grades.csv")
		dan <- read("dandelion_audit", "c2", "c2.dandelion_input_audit.csv")
		co <- read("coloc", "c3", "c3.coloc_summary.csv")
		# Reuse the project's strict overlap helper without sourcing any model worker.
		helper <- parse(file.path(Sys.getenv("LE8_FDIR", unset = "."), "0.common.R"))
		for (ex in helper) if (is.call(ex) && identical(ex[[1]], as.name("<-")) && identical(ex[[2]], as.name("le8_same_locus_evidence")))
			eval(ex)
		same <- le8_same_locus_evidence(mr, co, if (layer == "prot")
			"protein" else "metabolite")
		tables$same_region <- same
		membership <- read("membership", "c4", "c4.proxy_membership_YS_YSP_NS.csv")
		med <- read("mediation", "c4", "c4.mediation_all.csv")
		bridge <- read("bridges", "c4", "c4.genetic_omic_disease_bridges.csv")
		matched_pgs <- read("matched_pgs", "c4", "c4.matched_PGS_bridges.csv")
		cigma <- read("cigma_status", "cell", "c5.cellulation_status.csv")
		read("cigma_annotation", "cell", "c5.CIGMA_annotation.csv")
		fm <- read("focus_metrics", "focus", "c4.focus.metrics.csv")
		fc <- read("focus_contrasts", "focus", "c4.focus.contrasts.csv")
		fp <- read("focus_proxy_accuracy", "focus", "c4.focus.proxy_accuracy.csv")
		fpc <- read("focus_pillars", "focus", "c4.focus.pillar_counts.csv")
		read("focus_design", "focus", "c4.focus.design.csv")
		read("focus_members", "focus", "c4.focus.panel_members.csv")
		imaging <- read("imaging_status", "c4", "c4.imaging_status.csv")
		perf <- read("performance", "prediction", "prediction_summary.csv")
		budget <- read("budget", "prediction", "review_budget_metrics.csv")
		cal <- read("calibration", "prediction", "review_calibration.csv")
		delta <- read("paired_delta", "prediction", "review_paired_delta_CI.csv")
		design <- read("prediction_design", "prediction", "review_design.csv")
		lead <- read("leadtime", "prediction", "leadtime_discrimination.csv")
		ev <- read("prediction_evidence", "prediction", "evidence_consolidation.csv")
		panel_members <- read("panel_members", "prediction", "review_panel_members.csv")
		go_dir <- le8_go_dir()
		add <- function(fig, p, keys, label) {
			if (is.null(p))
				return(invisible(NULL))
			figures[[fig]][[length(figures[[fig]]) + 1L]] <<- p
			provenance[[length(provenance) + 1L]] <<- tibble(
				figure = fig, panel = LETTERS[length(figures[[fig]])],
				title = label, tables = paste(keys, collapse = ";"), selection = "Prespecified panel recipe; display ranking only, not prediction selection"
			)
		}
		# Figure 1: four panels; counts, volcano, anchors and prevalent concordance.
		if (pub_ok(a, c("term", "beta", "p.value", "estimate", "conf.low", "conf.high"))) {
			z <- a |>
				filter(is.finite(beta), is.finite(p.value)) |>
				mutate(support = ifelse(p.value < 0.05 / nrow(a), "Bonferroni", "Other"))
			p <- ggplot(z, aes(beta, - log10(pmax(p.value, 1e-300)), color = support)) +
				geom_point(alpha = 0.65, size = 1.1) +
				geom_hline(yintercept =  - log10(0.05 / nrow(a)), linetype = 2, color = "grey55") +
				scale_color_manual(values = c(
					Bonferroni = pub_cols[2],
					Other = "#B6C1CC"
				)) +
				labs(title = "Incident association scan", x = "Log hazard ratio per SD", y = expression( - log[10](P))) +
				pub_theme()
			add("Fig1", p, "incident", "Incident association scan")
			counts <- tibble(
				analysis = c("Incident: tested", "Incident: Bonferroni", "Prevalent: tested", "Prevalent: Bonferroni"),
				n = c(nrow(a), sum(a$p.value < 0.05 / nrow(a), na.rm = TRUE), nrow(prev), if (nrow(prev)) sum(prev$p.value <
					0.05 / nrow(prev), na.rm = TRUE) else 0)
			)
			tables$association_counts <- counts
			add(
				"Fig1", ggplot(counts, aes(n, reorder(analysis, n))) +
					geom_col(fill = pub_cols[1], width = 0.65) +
					geom_text(aes(label = n), hjust =  - 0.12, size = 3.5) +
					scale_x_continuous(expand = expansion(mult = c(
						0,
						0.18
					))) +
					labs(title = "Association yield", x = "Biomarkers", y = NULL) +
					pub_theme(), "association_counts",
				"Association yield"
			)
			aa <- a |>
				arrange(p.value)
			add("Fig1", pub_forest(
				aa, "term", "estimate", "conf.low", "conf.high", "Leading incident associations",
				"Hazard ratio per SD", 1
			), "incident", "Leading incident associations")
			if (pub_ok(prev, c("term", "beta"))) {
				z <- inner_join(a |>
					select(term, incident = beta), prev |>
					select(term, prevalent = beta), by = "term")
				add("Fig1", ggplot(z, aes(incident, prevalent)) +
					geom_hline(yintercept = 0, color = "grey85") +
					geom_vline(
						xintercept = 0,
						color = "grey85"
					) +
					geom_point(color = pub_cols[1], alpha = 0.5, size = 1) +
					labs(
						title = "Incident versus prevalent associations",
						x = "Incident log HR", y = "Prevalent log OR"
					) +
					pub_theme(), c("incident", "prevalent"), "Incident versus prevalent")
			}
		}
		captions["Fig1"] <- "Baseline prevalent cases contribute to descriptive and prevalent analyses; incident models exclude baseline disease. Bonferroni correction is within each omics scan. Prevalent odds ratios and incident hazard ratios estimate different quantities."
		# Figure 2: inherited propensity, cis/local MR, colocalization, and diagnostic coverage.
		if (pub_ok(paired, c("feature", "model", "component", "beta", "std.error"))) {
			z <- paired |>
				filter(grepl("^separate;", model), is.finite(beta), is.finite(std.error))
			shown <- head(unique(z$feature), 10)
			z <- z |>
				filter(feature %in% shown) |>
				mutate(feature = factor(feature, levels = rev(shown)), lo = beta - 1.96 * std.error, hi = beta + 1.96 *
					std.error)
			if (nrow(z))
				add("Fig2", ggplot(z, aes(beta, feature, color = component)) +
					geom_vline(
						xintercept = 0, color = "grey70",
						linetype = 2
					) +
					geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.2, position = position_dodge(width = 0.5)) +
					geom_point(position = position_dodge(width = 0.5)) +
					scale_color_manual(values = pub_cols) +
					labs(
						title = "Paired measured and inherited propensity",
						subtitle = "Identical participants and covariates within each biomarker", x = "Log HR per own SD (95% CI)",
						y = NULL, color = NULL
					) +
					pub_theme(), "paired", "Paired measured and inherited propensity")
		}
		primary <- if (layer == "prot")
			"cis" else "local"
		if (pub_ok(mr, c("analysis", "b", "se", "pval", "exposure"))) {
			mm <- mr |>
				filter(analysis == primary) |>
				arrange(pval) |>
				mutate(lo = b - 1.96 * se, hi = b + 1.96 * se)
			add(
				"Fig2", pub_forest(mm, "exposure", "b", "lo", "hi", paste(primary, "MR estimates"), "Effect on GWAS outcome scale"),
				"mr", "Primary MR estimates"
			)
		}
		if (pub_ok(co, c("PP.H4", "PP.H4_robust_min", "locus_class", "status"))) {
			z <- co |>
				filter(status == "ok", is.finite(PP.H4), is.finite(PP.H4_robust_min))
			if (nrow(z))
				add("Fig2", ggplot(z, aes(PP.H4, PP.H4_robust_min, color = locus_class)) +
					geom_point(alpha = 0.65) +
					geom_hline(yintercept = 0.7, linetype = 2) +
					scale_color_manual(values = pub_cols) +
					labs(
						title = "Colocalization and prior sensitivity",
						x = "PP.H4, default prior", y = "Minimum PP.H4 across priors"
					) +
					pub_theme(), "coloc", "Colocalization prior sensitivity")
		}
		add("Fig2", pub_count(grades, "evidence_grade", "MR diagnostic coverage"), "grades", "MR diagnostic coverage")
		captions["Fig2"] <- "PGS is an inherited genetic score, not a protein concentration measured at birth. Paired measured and PGS estimates use identical people and covariates, standardized to their own SD; their units are not interchangeable. Cis/local MR is primary; trans/distal associations remain exploratory. Colocalization supports a shared region, not proof of causal direction."
		# Figure 3: four panels capturing LE8 supervision, reproducibility and mediation.
		if (pub_ok(membership, c("in_YS", "in_YSP_plus", "in_NS", "primary_component", "r_disc", "r_rep"))) {
			z <- membership |>
				filter(in_YS %in% TRUE)
			pillars <- tibble(primary_component = c("diet", "pa", "smoke", "sleep", "bmi", "nonhdl", "hba1c", "bp")) |>
				left_join(z |>
					count(primary_component), by = "primary_component") |>
				mutate(n = coalesce(n, 0L))
			tables$pillar_counts <- pillars
			add(
				"Fig3", ggplot(pillars, aes(factor(primary_component, levels = primary_component), n)) +
					geom_col(
						fill = pub_cols[1],
						width = 0.65
					) +
					geom_text(aes(label = n), vjust =  - 0.4) +
					scale_y_continuous(expand = expansion(mult = c(
						0,
						0.15
					))) +
					labs(title = "LE8 pillars represented among YS proxies", x = NULL, y = "Biomarkers") +
					pub_theme(),
				"pillar_counts", "LE8 pillar coverage"
			)
			counts <- tibble(set = c("YS", "YSP plus", "NS"), n = c(sum(membership$in_YS %in% TRUE), sum(membership$in_YSP_plus %in%
				TRUE), sum(membership$in_NS %in% TRUE)))
			tables$set_counts <- counts
			add("Fig3", ggplot(counts, aes(set, n, fill = set)) +
				geom_col(width = 0.6, show.legend = FALSE) +
				geom_text(aes(label = n),
					vjust =  - 0.4
				) +
				scale_fill_manual(values = pub_cols) +
				scale_y_continuous(expand = expansion(mult = c(
					0,
					0.15
				))) +
				labs(title = "Proxy and comparator sets", x = NULL, y = "Biomarkers (sets may overlap)") +
				pub_theme(), "set_counts", "Proxy and comparator sets")
			if (nrow(z))
				add("Fig3", ggplot(z, aes(r_disc, r_rep, color = primary_component)) +
					geom_abline(
						slope = 1, intercept = 0,
						color = "grey65", linetype = 2
					) +
					geom_point(alpha = 0.6, size = 1.4) +
					scale_color_manual(values = pub_cols) +
					labs(title = "Discovery and replication profiles", x = "Discovery partial correlation", y = "Replication partial correlation") +
					pub_theme(), "membership", "Replication profiles")
		}
		if (pub_ok(med, c("component", "feature", "indirect_beta", "indirect_lo", "indirect_hi", "FDR_indirect"))) {
			z <- med |>
				arrange(FDR_indirect) |>
				mutate(label = paste(component, feature, sep = " / "))
			add("Fig3", pub_forest(
				z, "label", "indirect_beta", "indirect_lo", "indirect_hi", "Exploratory mediation paths",
				"Product-of-coefficients indirect effect", 0, 10
			), "mediation", "Mediation paths")
		}
		captions["Fig3"] <- "LE8 → omics → disease ← omics ← genetics is the organizing hypothesis. YS and YSP-plus are supervised and additional proxies; NS is an independently defined unsupervised comparator and may overlap them. Baseline mediation is exploratory and does not establish intervention effects. An empty genetic bridge table means that connection was not demonstrated."
		# Figure 4: held-out discrimination, assay budgets, calibration, paired uncertainty.
		if (pub_ok(perf, c("biom_set", "model", "AUC", "C_index"))) {
			z <- perf |>
				filter(is.finite(AUC))
			add(
				"Fig4", ggplot(z, aes(AUC, biom_set, color = model)) +
					geom_point(
						position = position_dodge(width = 0.5),
						size = 2
					) +
					scale_color_manual(values = pub_cols) +
					scale_y_discrete(labels = function(x) stringr::str_wrap(
						x,
						25
					)) +
					labs(title = "Held-out discrimination", x = "IPCW AUC (recorded horizon)", y = NULL) +
					pub_theme(),
				"performance", "Held-out discrimination"
			)
		}
		if (pub_ok(budget, c("ablation", "horizon", "actual_n_assays", "AUC", "paradigm", "status")) && all(LE8_ASSAY_BUDGETS %in%
			budget$actual_n_assays)) {
			z <- budget |>
				filter(ablation == "none", horizon == 10, status == "ok")
			if (nrow(z))
				add("Fig4", ggplot(z, aes(actual_n_assays, AUC, color = paradigm, group = paradigm)) +
					geom_line(na.rm = TRUE) +
					geom_point(size = 2) +
					scale_color_manual(values = pub_cols) +
					labs(
						title = "Performance at explicit assay budgets",
						x = "Actual fitted assays", y = "10-year IPCW AUC"
					) +
					pub_theme(), "budget", "Assay budget comparison")
		}
		if (pub_ok(cal, c("horizon", "ablation", "budget", "predicted", "observed_ipcw", "model"))) {
			z <- cal |>
				filter(horizon == 10, ablation == "none", budget %in% c(0, 10))
			if (nrow(z))
				add("Fig4", ggplot(z, aes(predicted, observed_ipcw, color = model)) +
					geom_abline(
						slope = 1, intercept = 0,
						linetype = 2, color = "grey65"
					) +
					geom_line() +
					geom_point(size = 1.4) +
					scale_color_manual(values = pub_cols) +
					labs(title = "Calibration at a fixed 10-assay budget", x = "Predicted 10-year risk", y = "Observed IPCW risk") +
					pub_theme(), "calibration", "Calibration at 10 assays")
		}
		if (pub_ok(delta, c("horizon", "model", "delta_AUC_lo", "delta_AUC_hi"))) {
			z <- delta |>
				filter(horizon == 10, is.finite(delta_AUC_lo), is.finite(delta_AUC_hi), !grepl("_[0-9]+$", model) |
					pub_num(sub(".*_", "", model)) %in% LE8_ASSAY_BUDGETS) |>
				slice_head(n = 12)
			if (nrow(z))
				add("Fig4", ggplot(z, aes(y = reorder(model, delta_AUC_lo))) +
					geom_vline(
						xintercept = 0, linetype = 2,
						color = "grey55"
					) +
					geom_segment(aes(x = delta_AUC_lo, xend = delta_AUC_hi, yend = reorder(
						model,
						delta_AUC_lo
					)), linewidth = 1.2, color = pub_cols[1]) +
					labs(
						title = "Paired improvement versus clinical model",
						x = "95% bootstrap interval for ΔAUC (10 years)", y = NULL
					) +
					pub_theme(), "paired_delta", "Paired AUC intervals")
		}
		captions["Fig4"] <- "Frozen models evaluated on the held-out sample. The clinical comparator is defined in the design sheet; it is not automatically SCORE2. Assay budgets and 10-year horizon are prespecified in this presentation. Bootstrap intervals describe frozen fits, not the uncertainty of the entire model-selection process. Deaths are censored; these risks are not competing-risk cumulative incidences."
		# Prefer the direct test of the article's hypothesis once it exists, including null/adverse results. Do not
		# choose panels by significance or best AUC.
		if (pub_ok(fpc, c("cohort", "component", "n")) && pub_ok(fp, c("component", "model", "delta_R2")) && length(figures[["Fig3"]]) >=
			4) {
			figures[["Fig3"]][[1]] <- ggplot(fpc, aes(component, n, fill = cohort)) +
				geom_col(position = "dodge") +
				labs(title = "LE8 proxy discovery: Yin and Yin + Yang", x = NULL, y = "Replicated proxies") +
				pub_theme() +
				theme(axis.text.x = element_text(angle = 35, hjust = 1))
			display_budget <- max(fm$budget, na.rm = TRUE)
			z <- fp |>
				filter(grepl(paste0("_", display_budget, "$"), model))
			figures[["Fig3"]][[2]] <- ggplot(z, aes(component, model, fill = delta_R2)) +
				geom_tile() +
				scale_fill_gradient2(
					low = pub_cols[2],
					mid = "white", high = pub_cols[1], midpoint = 0
				) +
				labs(title = paste(
					"Held-out LE8 reconstruction at",
					display_budget, "assays"
				), x = NULL, y = NULL, fill = "Incremental R²") +
				pub_theme() +
				theme(axis.text.x = element_text(
					angle = 35,
					hjust = 1
				))
			for (i in seq_along(provenance)) if (provenance[[i]]$figure == "Fig3" && provenance[[i]]$panel %in% c(
				"A",
				"B"
			)) {
				provenance[[i]]$tables <- if (provenance[[i]]$panel == "A")
					"focus_pillars" else "focus_proxy_accuracy"
				provenance[[i]]$title <- if (provenance[[i]]$panel == "A")
					"Yin/Yang proxy discovery" else "Held-out LE8 reconstruction"
			}
			captions["Fig3"] <- paste(captions["Fig3"], "Panels A/B use only outer-training donors for discovery and fit; Yang donors contribute to proxy learning. Incremental R² is measured beyond basic covariates on held-out incident participants, and is not intervention responsiveness.")
		}
		if (pub_ok(fm, c("stratum", "landmark", "budget", "AUC", "paradigm")) && pub_ok(fc, c(
			"stratum", "landmark",
			"model", "delta_AUC", "delta_lo", "delta_hi"
		))) {
			figures[["Fig4"]] <- list()
			provenance <- Filter(function(x) x$figure != "Fig4", provenance)
			z <- fm |>
				filter(stratum == "All", landmark == 0, budget > 0)
			clinical_auc <- fm$AUC[fm$model == "Clinical" & fm$stratum == "All" & fm$landmark == 0][1]
			add("Fig4", ggplot(z, aes(budget, AUC, color = paradigm)) +
				geom_line() +
				geom_point() +
				geom_hline(
					yintercept = clinical_auc,
					linetype = 2, color = "grey50"
				) +
				scale_x_continuous(breaks = sort(unique(z$budget))) +
				labs(
					title = "Matched assay budgets",
					subtitle = sprintf("Dashed line: clinical covariates alone (AUC %.3f)", clinical_auc), x = "Measured assays",
					y = "10-year IPCW AUC"
				) +
				pub_theme(), "focus_metrics", "Matched assay budgets")
			z <- fc |>
				filter(stratum == "All", landmark == 0)
			add("Fig4", ggplot(z, aes(delta_AUC, reorder(model, delta_AUC), color = contrast)) +
				geom_vline(
					xintercept = 0,
					linetype = 2
				) +
				geom_errorbar(aes(xmin = delta_lo, xmax = delta_hi), orientation = "y", width = 0.2) +
				geom_point() +
				labs(
					title = "Supervision and added-Yang contrasts", x = "Paired ΔAUC with 95% bootstrap interval",
					y = NULL
				) +
				pub_theme(), "focus_contrasts", "Supervision and Yang contrasts")
			z <- fm |>
				filter(stratum == "All", budget %in% c(0, max(fm$budget, na.rm = TRUE)))
			add(
				"Fig4", ggplot(z, aes(landmark, AUC, color = paradigm)) +
					geom_line() +
					geom_point() +
					labs(
						title = "Conditional prediction to baseline year 10",
						x = "Disease-free landmark (years)", y = "IPCW AUC after landmark"
					) +
					pub_theme(), "focus_metrics",
				"Landmark validation"
			)
			z <- fm |>
				filter(stratum != "All", landmark == 0, budget %in% c(0, max(fm$budget, na.rm = TRUE)))
			if (nrow(z))
				add(
					"Fig4", ggplot(z, aes(AUC, paradigm, color = stratum)) +
						geom_point(position = position_dodge(width = 0.4)) +
						labs(title = "Baseline inflammation strata", x = "10-year IPCW AUC", y = NULL) +
						pub_theme(), "focus_metrics",
					"Baseline inflammation strata"
				) else add(
				"Fig4", pub_count(tables$focus_members, "model", "Assays in each evaluated panel"), "focus_members",
				"Panel composition"
			)
			captions["Fig4"] <- "Same incident test cohort, same assay budget, same Cox model and clinical covariates. Yang donors enter only proxy discovery. YSplus reserves 80% of slots (rounded up) for YS by default; see focus_design. Low baseline inflammation is an outcome-independent marker stratum, not a validated CAD subtype. These exploratory comparisons were designed after inspecting prior results; bootstrap intervals condition on frozen fits and require external confirmation. Death is censored."
		}
		if (pub_ok(matched_pgs, c("feature", "r_disc", "r_rep", "PGS_disease_FDR")) && length(figures[["Fig3"]]) >=
			3) {
			disease_label <- paste0("PGS-", trait, " FDR < 0.05")
			z <- matched_pgs |>
				filter(is.finite(r_disc), is.finite(r_rep)) |>
				mutate(disease_support = case_when(!is.finite(PGS_disease_FDR) ~ "Unavailable", PGS_disease_FDR < 0.05 ~
					disease_label, TRUE ~ "Not detected"))
			figures[["Fig3"]][[3]] <- ggplot(z, aes(r_disc, r_rep, color = disease_support)) +
				geom_abline(
					slope = 1,
					intercept = 0, linetype = 2, color = "grey70"
				) +
				geom_point(alpha = 0.55, size = 1.2) +
				scale_color_manual(values = pub_cols) +
				labs(
					title = "Matched inherited propensity → measured omic", subtitle = "Each biomarker is paired with its own PGS",
					x = "Discovery partial correlation", y = "Replication partial correlation", color = NULL
				) +
				pub_theme()
			for (i in seq_along(provenance)) if (provenance[[i]]$figure == "Fig3" && provenance[[i]]$panel == "C") {
				provenance[[i]]$tables <- "matched_pgs"
				provenance[[i]]$title <- "Matched PGS calibration"
			}
			captions["Fig3"] <- sub("An empty genetic bridge table means that connection was not demonstrated.", "",
				captions["Fig3"],
				fixed = TRUE
			)
			captions["Fig3"] <- paste(captions["Fig3"], "Panel C uses two training halves. Source-GWAS overlap may inflate PGS calibration; shared associations do not identify a causal mediation chain.")
		}
		# Conservative evidence matrix built afresh, never endorses historical prediction grades.
		if (pub_ok(a, c("term", "p.value", "FDR"))) {
			strict <- a |>
				transmute(feature = term, observed = ifelse(is.finite(FDR), FDR < 0.05, NA), p = p.value)
			pg <- if (pub_ok(pgs, c("analysis", "feature", "pgs_FDR")))
				pgs |>
					filter(analysis == "Incident") |>
					distinct(feature, .keep_all = TRUE) |>
					transmute(feature, PGS = ifelse(is.finite(pgs_FDR), pgs_FDR < 0.05, NA)) else tibble(feature = character(), PGS = logical())
			me <- if (pub_ok(mr, c("analysis", "exposure", "FDR_all")))
				mr |>
					filter(analysis == primary) |>
					group_by(exposure) |>
					summarise(MR = if (all(!is.finite(FDR_all)))
						NA else any(FDR_all < 0.05, na.rm = TRUE), .groups = "drop") |>
					rename(feature = exposure) else tibble(feature = character(), MR = logical())
			ys <- if (pub_ok(membership, c("feature", "in_YS")))
				membership |>
					distinct(feature, .keep_all = TRUE) |>
					transmute(feature, YS = in_YS %in% TRUE) else tibble(feature = character(), YS = logical())
			sl <- if (nrow(same))
				same |>
					group_by(feature) |>
					summarise(same_region = any(eligible %in% TRUE), .groups = "drop") else tibble(feature = character(), same_region = logical())
			strict <- strict |>
				left_join(pg, by = "feature") |>
				left_join(me, by = "feature") |>
				left_join(sl, by = "feature") |>
				left_join(ys, by = "feature") |>
				mutate(role = case_when(same_region %in% TRUE ~ "Cis/local region-supported candidate", MR %in% TRUE ~
					"MR support; locus unresolved", YS %in% TRUE ~ "LE8 proxy", TRUE ~ "Association / unresolved")) |>
				arrange(desc(same_region %in% TRUE), p)
			tables$publication_evidence <- strict
			top <- head(strict, 18)
			hm <- top |>
				pivot_longer(c(observed, PGS, MR, same_region, YS), names_to = "domain", values_to = "support") |>
				mutate(domain = factor(domain, levels = c("observed", "PGS", "MR", "same_region", "YS"), labels = c(
					"Observed",
					"PGS", "MR", "Same region", "YS"
				)), feature = factor(feature, levels = rev(top$feature)), state = case_when(is.na(support) ~
					"Unavailable", support ~ "Supported", TRUE ~ "Not detected"))
			add("Fig5", ggplot(hm, aes(domain, feature, fill = state)) +
				geom_tile(color = "white", linewidth = 0.5) +
				scale_fill_manual(values = c(Supported = pub_cols[1], `Not detected` = "#E5EAF0", Unavailable = "#F4D9B3")) +
				labs(title = "Conservative cross-domain evidence", x = NULL, y = NULL) +
				pub_theme(), c(
				"publication_evidence",
				"same_region"
			), "Cross-domain evidence")
			add(
				"Fig5", pub_count(strict, "role", "Evidence roles, without causal proof"), "publication_evidence",
				"Evidence roles"
			)
		}
		add("Fig5", pub_enrich(enrich, "Functional context of incident associations", go_dir), "enrich", "Incident enrichment")
		if (pub_ok(lead, c("horizon", "AUC", "kind", "method"))) {
			z <- lead |>
				filter(kind == "Omic score", method %in% c(
					"C4 YS", "C4 YSplus", "C4 NS", "Pradeep-style / glmnet",
					"Yu-style / LightGBM"
				), is.finite(AUC))
			if (nrow(z))
				add("Fig5", ggplot(z, aes(horizon, AUC, color = method)) +
					geom_line() +
					geom_point(size = 1.5) +
					scale_color_manual(values = pub_cols) +
					labs(title = "Discrimination with minimum lead time", x = "Minimum years before diagnosis", y = "Case/control AUC") +
					pub_theme(), "leadtime", "Minimum lead-time discrimination")
		}
		captions["Fig5"] <- "Cis/local MR and conservative colocalization must overlap retained instruments in the same region. This is region-level corroboration, not signal-resolved causality; missing tests are distinct from negative tests. The minimum-lead-time case/control AUC is not the IPCW AUC in Fig4. Annotation enrichment is contextual evidence."
		# The generic association landscape is retained as a supplement. Main Figure 1 now joins C1 matched scores,
		# calibrated components and C3 loci.
		old_fig1 <- figures[["Fig1"]]
		if (length(old_fig1))
			pgs_save_panels(old_fig1, file.path(root, "supplement"), "association_landscape", tables[c(
				"incident",
				"prevalent", "cohort"
			)], paste(trait, toupper(layer), "association context"), captions[["Fig1"]])
		figures[["Fig1"]] <- list()
		provenance <- Filter(function(z) !identical(z$figure, "Fig1"), provenance)
		pf <- file.path(source_root, "c1_correlate", "c1.pgs_focus.rds")
		focus_pgs <- if (file.exists(pf))
			tryCatch(readRDS(pf), error = function(e) list()) else list()
		if (pgs_ok(focus_pgs$paired, c("feature", "evidence_pattern"))) {
			triangulation <- make_c3_pgs_integration(coloc_summary = co, focus = focus_pgs)
			loci <- attr(triangulation, "loci")
			panels <- pgs_main_panels(focus_pgs, triangulation, loci)
			for (nm in names(focus_pgs)) if (is.data.frame(focus_pgs[[nm]]))
				tables[[paste0("PGS_", nm)]] <- focus_pgs[[nm]]
			tables$PGS_loci <- loci
			tables$PGS_triangulation <- as.data.frame(triangulation)
			keys <- grep("^PGS_", names(tables), value = TRUE)
			for (nm in names(panels)) add("Fig1", panels[[nm]], keys, paste("PGS discordance:", nm))
			audit$PGS_focus <- tibble(
				table = "PGS_focus", source = pf, rows = nrow(focus_pgs$paired), bytes = file.info(pf)$size,
				modified = as.character(file.info(pf)$mtime), status = "available"
			)
		} else audit$PGS_focus <- tibble(
			table = "PGS_focus", source = pf, rows = 0L, bytes = NA_real_, modified = NA_character_,
			status = "Run pgs_focus; legacy unequal-sample directions are not used as main evidence"
		)
		captions["Fig1"] <- "Identical-sample measured/PGS associations and cross-fitted captured (G) versus remaining (R) components. Joint contrasts include covariance; conditional CIs and refitted-bootstrap intervals are distinguished. R is not a pure lifestyle fraction. Colocalization is locus-specific, with prior sensitivity retained. Candidate selection is exploratory; opposite directions do not prove antagonistic pleiotropy."
		manifest <- if (length(provenance))
			bind_rows(provenance) else tibble(figure = character(), panel = character(), title = character(), tables = character(), selection = character())
		manifest$analysis <- manifest$figure
		source_audit <- bind_rows(audit)
		# Preserve old main/supplement versions before assembling the new edition.
		old <- list.files(root, pattern = "^(Fig[0-9]+|FigS[0-9]+)\\.(png|xlsx)$", full.names = TRUE)
		le8_archive_figure_files(old, root)
		status <- list()
		figure_number <- 0L
		rendered_manifest <- list()
		dpi <- as.integer(Sys.getenv("LE8_FINAL_DPI", unset = "320"))
		for (fig in paste0("Fig", 1 : 5)) {
			expanded <- lapply(figures[[fig]], le8_expand_facets)
			pp <- unlist(expanded, recursive = FALSE)
			n <- length(pp)
			if (!n) {
				status[[fig]] <- tibble(figure = NA_character_, analysis = fig, panels = 0, status = "withheld: no usable panels")
				next
			}
			keys <- unique(unlist(strsplit(manifest$tables[manifest$analysis == fig], ";", fixed = TRUE)))
			source_rows <- manifest[manifest$analysis == fig, ]
			panel_rows <- map_dfr(seq_along(expanded), function(i) {
				if (!length(expanded[[i]]))
					return(tibble())
				map_dfr(expanded[[i]], function(p) source_rows[i, ] |>
					mutate(source_panel = panel, title = paste(as.character(p$labels$title), collapse = " ")))
			})
			if (n > 6L) {
				pgs_save_panels(
					pp[7 : n], file.path(root, "supplement"), paste0("overflow.", fig), tables[keys],
					paste(trait, toupper(layer), fig, "additional panels"), captions[[fig]]
				)
				pp <- pp[1 : 6]
				panel_rows <- panel_rows[1 : 6, ]
				n <- 6L
			}
			actual <- character()
			for (start in seq(1L, n, by = 6L)) {
				page <- pp[start : min(n, start + 5L)]
				figure_number <- figure_number + 1L
				dest <- fig
				actual <- c(actual, dest)
				page_rows <- panel_rows[start : min(n, start + 5L), ] |>
					mutate(figure = dest, panel = LETTERS[seq_along(page)])
				rendered_manifest[[length(rendered_manifest) + 1L]] <- page_rows
				sheets <- c(list(provenance = page_rows, caption = tibble(caption = captions[[fig]])), tables[keys])
				p <- wrap_plots(page, ncol = if (length(page) == 1)
					1 else 2, widths = c(1, 1)) + plot_annotation(title = paste(trait, toupper(layer), "|", switch(fig,
					Fig1 = "Measured versus genetically predicted biomarkers",
					Fig2 = "Genetic triangulation",
					Fig3 = "LE8 connections",
					Fig4 = "Prediction and assay budget",
					Fig5 = "Evidence consolidation"
				)), caption = stringr::str_wrap(gsub("Fig4", "the prediction figure",
					captions[[fig]],
					fixed = TRUE
				), 155), tag_levels = "A", theme = theme(plot.title = element_text(
					size = 17,
					face = "bold"
				), plot.caption = element_text(size = 9, hjust = 0), plot.tag = element_text(face = "bold")))
				ggsave(file.path(root, paste0(dest, ".png")), p, width = 19, height = ceiling(length(page) / 2) * 6 +
					0.9, dpi = dpi, bg = "white", limitsize = FALSE)
				pub_workbook(sheets, file.path(root, paste0(dest, ".xlsx")))
			}
			status[[fig]] <- tibble(figure = paste(actual, collapse = ";"), analysis = fig, panels = n, status = "rendered")
		}
		manifest <- if (length(rendered_manifest))
			bind_rows(rendered_manifest) else tibble(
			figure = character(), panel = character(), title = character(), tables = character(), selection = character(),
			analysis = character(), source_panel = character()
		)
		# Retain every available module figure in final supplements. Manuscript selection must not silently delete
		# a view from the reviewable output set.  No raster tiling or stretching: preserve original resolution and
		# aspect ratio.
		supp <- list()
		sn <- 0L
		for (code in names(modules)) {
			module <- modules[[code]]
			mf <- file.path(module_dirs[[code]], "figure_manifest.csv")
			if (!file.exists(mf))
				next
			fm <- as_tibble(fread(mf))
			files <- file.path(module_dirs[[code]], fm$file[fm$panels <= 6])
			for (src in files[file.exists(files)]) {
				sn <- sn + 1L
				dest <- paste0("FigS", sn, ".png")
				if (!file.copy(src, file.path(root, dest), overwrite = TRUE))
					stop("Cannot copy supplement: ", src)
				sx <- sub("[.]png$", ".xlsx", src)
				if (file.exists(sx) && !file.copy(sx, file.path(root, sub("[.]png$", ".xlsx", dest)), overwrite = TRUE))
					stop("Cannot copy supplementary workbook: ", sx)
				supp[[sn]] <- tibble(figure = dest, source = src, module = module, review = "Retained module figure; original source layout and previous editions are preserved")
			}
		}
		pub_workbook(c(
			list(sources = source_audit, figure_status = bind_rows(status), panels = manifest, supplements = bind_rows(supp)),
			tables[c(
				"cohort", "prediction_design", "dandelion_audit", "imaging_status", "bridges", "matched_pgs",
				"cigma_status", "cigma_annotation", "focus_metrics", "focus_contrasts", "focus_proxy_accuracy", "focus_pillars",
				"focus_design", "focus_members", "paired", "same_region", "publication_evidence", "panel_members",
				"enrich_prev"
			)]
		), file.path(root, "TablesS.xlsx"))
		fwrite(manifest, file.path(root, "figure_manifest.csv"))
		supp_table <- if (length(supp))
			bind_rows(supp) else tibble(figure = character(), source = character(), module = character(), review = character())
		fwrite(supp_table, file.path(root, "supplement_manifest.csv"))
		fwrite(bind_rows(status), file.path(root, "status.csv"))
		pub_report(root, trait, layer, tables, source_audit, bind_rows(status), captions)
		message(
			"FINAL: ", root, "; ", sum(vapply(status, function(z) z$status == "rendered", logical(1))), " main figures; ",
			sn, " supplements"
		)
		invisible(status)
	}
	pub_workbook <- function(sheets, file) {
		wb <- createWorkbook()
		sheets <- sheets[!vapply(sheets, is.null, logical(1))]
		names(sheets) <- make.unique(substr(names(sheets), 1, 28))
		for (nm in names(sheets)) {
			z <- as.data.frame(sheets[[nm]])
			if (!ncol(z))
				z <- data.frame(status = "No available result")
			addWorksheet(wb, nm)
			if (nrow(z))
				writeDataTable(wb, nm, z, tableStyle = "TableStyleMedium2") else writeData(wb, nm, z)
			freezePane(wb, nm, firstRow = TRUE)
			setColWidths(wb, nm, cols = seq_len(ncol(z)), widths = 18)
		}
		saveWorkbook(wb, file, overwrite = TRUE)
	}
	pub_report <- function(root, trait, layer, t, source_audit, status, captions) {
		lines <- c(
			paste("#", trait, "/", layer, "— 结果与写作清单"), "", paste("生成时间：", Sys.time()),
			"", "本摘要读取模块内的聚合结果；缺失结果不会写成阴性发现。论文数值以本次输出的来源表和实际运行设计为准。",
			""
		)
		if (pub_ok(t$cohort, c("N_omics", "incident_events", "prevalent_cases", "features"))) {
			d <- t$cohort[1, ]
			lines <- c(lines, sprintf(
				"Omics cohort：%s 人；incident events %s；baseline prevalent cases %s；%s 个 biomarker。各模型 complete-case N 需单独报告。",
				d$N_omics, d$incident_events, d$prevalent_cases, d$features
			), "")
		}
		for (key in c("incident", "prevalent")) {
			d <- t[[key]]
			if (!pub_ok(d, c("p.value", "FDR")))
				next
			lines <- c(lines, sprintf(
				"%s：测试 %d 项，Bonferroni 显著 %d 项，BH-FDR <0.05 为 %d 项。",
				key, nrow(d), sum(d$p.value < 0.05 / nrow(d), na.rm = TRUE), sum(d$FDR < 0.05, na.rm = TRUE)
			))
		}
		if (pub_ok(t$PGS_paired, c("evidence_pattern", "feature"))) {
			d <- t$PGS_paired
			opposite <- d$feature[d$evidence_pattern == "Both supported: opposite"]
			lines <- c(
				lines, "", sprintf(
					"同样本 measured/PGS：%d 项可匹配指标；两边均 BH-FDR <0.05 且方向相反 %d 项。",
					nrow(d), length(opposite)
				), paste("方向相反候选：", if (length(opposite)) paste(opposite, collapse = ", ") else "无达到该证据标准的候选"),
				"G/R 的效应比较、校准能力、随访窗口及位点证据见 Fig1.xlsx；R 不等于纯生活方式来源。"
			)
		}
		if (pub_ok(t$performance, c("model", "biom_set", "AUC", "C_index"))) {
			d <- t$performance |>
				filter(model == "Combined", is.finite(AUC))
			lines <- c(
				lines, "", "验证集 Combined 模型（逐模型报告，不按验证集结果挑选最终模型）：",
				""
			)
			for (i in seq_len(nrow(d))) lines <- c(lines, sprintf(
				"- %s：C-index %.3f；%s 年 IPCW AUC %.3f；%s 个 biomarker。",
				d$biom_set[i], d$C_index[i], d$AUC_horizon_years[i], d$AUC[i], d$n_selected[i]
			))
		}
		if (pub_ok(t$same_region, "eligible"))
			lines <- c(lines, "", sprintf(
				"同一 cis/local 区域 MR + 保守共定位支持：%d 个 biomarker（区域层面的候选，不等于因果证明）。",
				n_distinct(t$same_region$feature[t$same_region$eligible %in% TRUE])
			))
		if (pub_ok(t$mediation, "FDR_indirect"))
			lines <- c(lines, sprintf(
				"探索性 mediation：%d 条路径，FDR <0.05 为 %d 条。", nrow(t$mediation),
				sum(t$mediation$FDR_indirect < 0.05, na.rm = TRUE)
			))
		lines <- c(
			lines, "", "## 尚未建立的结论与写作边界", "", "- PGS 反映遗传倾向，不能当作出生时实测 omics；需交代来源 GWAS、样本重叠、权重和预测能力。",
			"- 诊断前时间模式保留用于文章比较，方法中注明为不同人的一次基线采样。",
			"- 预测有效性与病因作用分开评价；无 MR 支持不等于没有临床预测价值。", "- 当前 clinical model 和论文中的 SCORE2/协变量不一定一致，不能仅凭 AUC 接近宣称严格复现。",
			"- 历史预测证据 的最小 P cis/trans 混选、不同区域合并和启发式权重不作为最终因果分级；Fig5/同区域表采用独立的保守规则。",
			"- 同期 LE8—omics mediation 不足以证明 lifestyle intervention 的因果效果。", "- 外部验证、锁定小 panel 后的完整选择流程验证、竞争风险及绝对风险可迁移性仍需另行建立。"
		)
		if (pub_ok(t$matched_pgs, "replicated"))
			lines <- c(lines, sprintf(
				"- 匹配 biomarker PGS 的遗传—实测 omic bridge：%d 项检验，%d 项两半样本复制；这不等同于因果中介链。",
				nrow(t$matched_pgs), sum(t$matched_pgs$replicated %in% TRUE)
			)) else if (is.null(t$bridges) || !nrow(t$bridges))
			lines <- c(lines, "- 遗传—omics—疾病 bridge 表没有结果：四维 connection 尚未由此模块贯通。")
		if (pub_ok(t$focus_metrics, c("stratum", "landmark", "model", "AUC"))) {
			z <- t$focus_metrics |>
				filter(stratum == "All", landmark == 0)
			lines <- c(lines, "", "## 同等 assay 预算的监督比较", "", paste0("- ", z$model, ": AUC=", sprintf(
				"%.4f",
				z$AUC
			)), "", "这些比较属于已查看旧结果之后的探索性分析；不能将最高 AUC 当作预先锁定的成功模型。Yin/Yang 的差异和 NS 对照见 focus_contrasts。")
		}
		if (pub_ok(t$coloc, "susie_status") && all(t$coloc$susie_status != "ok", na.rm = TRUE))
			lines <- c(lines, "- 当前没有成功的 SuSiE 信号级 fine-mapping；ABF 的 H4 条件 SNP 后验不能当作 trait fine-mapping credible set。")
		if (pub_ok(t$dandelion_audit, c("metric", "value")) && any(t$dandelion_audit$metric == "primary_eligible" &
			t$dandelion_audit$value == "FALSE"))
			lines <- c(lines, "- DANDELION 当前标记 primary_eligible=FALSE，仅作为已记录的敏感性分析，不宣称完整复现原始方法。")
		missing <- source_audit |>
			filter(status != "available")
		if (nrow(missing))
			lines <- c(lines, "", "缺失或空表：", paste0("- ", missing$table, "：", missing$source))
		lines <- c(
			lines, "", "## 图表与文章结构", "", paste0(
				"- ", status$figure, "：", status$status, "（",
				status$panels, " panels）"
			), "", "正文可围绕“LE8-supervised、遗传证据分层且限制 assay 数量的疾病风险评估”组织；遗传、可干预解释和预测收益必须分别有对应结果，不能由统计关联推导 intervention promise。",
			"", "Fig1 实测与遗传预测的分歧 → Fig2 遗传位点证据 → Fig3 LE8 connections → Fig4 预测与 assay budget → Fig5 证据整合。",
			"Fig*.xlsx 保存来源数据和 panel 映射；TablesS.xlsx 保存设计、缺口、同区域证据和入选附图路径。",
			""
		)
		writeLines(lines, file.path(root, "results.md"), useBytes = TRUE)
		writeLines(c("# Figure legends", "", unlist(lapply(names(captions), function(n) c(
			paste("##", n), "", captions[[n]],
			""
		)))), file.path(root, "legends.md"))
	}


	traits <- Sys.getenv("Y", unset = "cvd_cad")
	layers <- strsplit(Sys.getenv("BIOM", unset = "prot,met"), ",", fixed = TRUE)[[1]]
	for (layer in layers) publication_run(traits, layer, Sys.getenv("LE8_ANALYSIS_ROOT", unset = "/mnt/d/analysis/le8"))
	le8_layer_figure_correspondence(file.path(Sys.getenv("LE8_ANALYSIS_ROOT", unset = "/mnt/d/analysis/le8"), traits))
	pgs_publication_joint(traits, Sys.getenv("LE8_ANALYSIS_ROOT", unset = "/mnt/d/analysis/le8"), layers)
}
