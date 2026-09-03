# 🚩 Packages, paths, and shared inputs
.pred_var_ok <- function(x) is.character(x) || is.factor(x) || (is.numeric(x) && is.finite(var(x, na.rm = TRUE)) && var(x, na.rm = TRUE) > 0)
.pred_pick_vars <- function(dat, vars) vars[vars %in% names(dat) & vapply(dat[, vars, drop = FALSE], .pred_var_ok, logical(1))]
.pred_rebuild_test <- function(te, need) { need <- intersect(unique(need), names(te)); if (!length(need)) te[0, , drop = FALSE] else te[complete.cases(te[, need, drop = FALSE]), , drop = FALSE] }

fcidx <- function(y, pred) {
	cc <- survival::concordance(y ~ pred)$concordance
	ifelse(is.na(cc), NA_real_, max(cc, 1 - cc))
}

RE_pred <- function(datTR, datTE, t2e.var, event.var, varX, vars.basic = NULL,
                    method = c("gaussian", "binomial", "cox"), opt = list()) {
	method <- match.arg(method)
	vb <- .pred_pick_vars(datTR, vars.basic %||% character())
	vx <- .pred_pick_vars(datTR, setdiff(varX, vb)); vv <- unique(c(vb, vx))
	need <- c(t2e.var, event.var, vv)
	# Outcomes must be observed. Predictors are imputed below using training-only
	# statistics, matching the held-out design used in the Pradeep reproduction.
	itr <- which(complete.cases(datTR[, c(t2e.var,event.var), drop = FALSE]))
	ite <- which(complete.cases(datTE[, c(t2e.var,event.var), drop = FALSE]))
	dtr <- datTR[itr, need, drop = FALSE]; dte <- datTE[ite, need, drop = FALSE]
	if (nrow(dtr) < 30 || nrow(dte) < 5 || !length(vv)) return(NULL)
	for (v in vv) {
		if (is.numeric(dtr[[v]]) || is.integer(dtr[[v]])) {
			m <- stats::median(as.numeric(dtr[[v]]),na.rm=TRUE); if(!is.finite(m)) m <- 0
			dtr[[v]][!is.finite(as.numeric(dtr[[v]]))] <- m; dte[[v]][!is.finite(as.numeric(dte[[v]]))] <- m
		} else {
			a <- as.character(dtr[[v]]); b <- as.character(dte[[v]])
			tb <- table(a,useNA="no"); mode0 <- if(length(tb)) names(which.max(tb)) else "Missing"
			a[is.na(a)|!nzchar(a)] <- mode0; lv <- unique(a); b[is.na(b)|!b%in%lv] <- mode0
			dtr[[v]] <- factor(a,levels=lv); dte[[v]] <- factor(b,levels=lv)
		}
	}
	ff <- as.formula(paste0("~", paste(paste0("`", vv, "`"), collapse = "+")))
	xtr <- model.matrix(ff, dtr)[, -1, drop = FALSE]; xte <- model.matrix(ff, dte)[, -1, drop = FALSE]
	keep <- apply(xtr, 2, function(z) is.finite(var(z)) && var(z) > 0)
	xtr <- xtr[, keep, drop = FALSE]; xte <- xte[, colnames(xtr), drop = FALSE]
	pf <- ifelse(colnames(xtr) %in% vb, 0, 1)
	lambda_rule <- opt$lambda_rule %||% "lambda.min"; inner_cv_n <- opt$inner_cv_n %||% 10
	type.measure <- opt$type.measure %||% switch(method, binomial = "auc", cox = "deviance", gaussian = "deviance")
	ytr <- switch(method, gaussian = as.numeric(dtr[[event.var]]), binomial = as.integer(dtr[[event.var]] == 1),
	               cox = survival::Surv(dtr[[t2e.var]], dtr[[event.var]]))
	if (method == "binomial" && length(unique(ytr)) < 2) return(NULL)
	nfolds0 <- if (method == "binomial") max(3, min(inner_cv_n, min(table(ytr)))) else inner_cv_n
	fit <- glmnet::cv.glmnet(xtr, ytr, family = method, alpha = opt$alpha %||% ifelse(method == "gaussian", .5, 1),
	                         nfolds = nfolds0, type.measure = type.measure, penalty.factor = pf)
	pred <- as.numeric(predict(fit, newx = xte, s = lambda_rule,
	                           type = ifelse(method == "binomial", "response", "link")))
	cidx <- if (method == "cox") fcidx(survival::Surv(dte[[t2e.var]], dte[[event.var]]), pred) else if (method == "binomial")
		as.numeric(pROC::auc(pROC::roc(as.integer(dte[[event.var]] == 1), pred, quiet = TRUE))) else
		suppressWarnings(stats::cor(as.numeric(dte[[event.var]]), pred, use = "complete.obs"))
	bb <- as.matrix(coef(fit, s = lambda_rule))
	list(pred = pred, cidx = cidx, n_sel = sum(bb[, 1] != 0), lambda = fit[[lambda_rule]],
	     test_index = ite, fit = if (isTRUE(opt$keep_fit)) fit else NULL,
	     coef = tibble::tibble(term = rownames(bb), beta = as.numeric(bb[, 1])) |> dplyr::filter(beta != 0))
}

ML_pred <- function(datTR, datTE, t2e.var, event.var, varX, vars.basic = NULL, method = c("xgb", "rfsrc"), opt = list()) {
	method <- match.arg(method)
	vv <- varX[varX %in% names(datTR) & vapply(datTR[, varX, drop = FALSE], \(x) is.numeric(x) && is.finite(var(x, na.rm = TRUE)) && var(x, na.rm = TRUE) > 0, logical(1))]
	need <- c(t2e.var, event.var, vv)
	dtr <- datTR[complete.cases(datTR[, need, drop = FALSE]), need, drop = FALSE]
	dte <- datTE[complete.cases(datTE[, need, drop = FALSE]), need, drop = FALSE]
	n_sel <- length(vv)
	if (!n_sel) return(list(pred = rep(0, nrow(dte)), cidx = 0.5, n_sel = 0))
	if (nrow(dtr) < 30 || nrow(dte) < 5) return(NULL)
	if (method == "xgb") {
		lab <- ifelse(as.numeric(dtr[[event.var]]) == 1, pmax(as.numeric(dtr[[t2e.var]]), 0.01), -pmax(as.numeric(dtr[[t2e.var]]), 0.01))
		fit <- xgboost::xgb.train(params = modifyList(list(objective = "survival:cox", eval_metric = "cox-nloglik", eta = 0.01, max_depth = 4, subsample = 0.7, colsample_bytree = 0.7, nthread = 1), opt$params %||% list()), data = xgboost::xgb.DMatrix(as.matrix(dtr[, vv, drop = FALSE]), label = lab), nrounds = opt$nrounds %||% 1000, verbose = 0)
		pred <- -as.numeric(predict(fit, xgboost::xgb.DMatrix(as.matrix(dte[, vv, drop = FALSE]))))
	} else {
		fit <- randomForestSRC::rfsrc(as.formula(paste0("Surv(", t2e.var, ",", event.var, ")~.")), data = dtr[, c(t2e.var, event.var, vv), drop = FALSE], ntree = opt$ntree %||% 500, nodesize = opt$nodesize %||% 10, mtry = opt$mtry %||% max(1, floor(sqrt(n_sel))), nsplit = opt$nsplit %||% 0, samptype = opt$samptype %||% "swor", sampsize = min(opt$sampsize %||% nrow(dtr), nrow(dtr)), block.size = 1)
		pred <- -as.numeric(predict(fit, dte)$predicted)
	}
	cc <- tryCatch(concordance(Surv(dte[[t2e.var]], dte[[event.var]]) ~ pred)$concordance, error = function(e) NA_real_)
	list(pred = pred, cidx = ifelse(!is.na(cc) && cc < .5, 1 - cc, cc), n_sel = n_sel)
}

DL_pred <- function(datTR, datTE, t2e.var, event.var, varX, vars.basic = NULL, method = c("tabpfn", "ft_transformer"), opt = list()) {
	method <- match.arg(method)
	prep10y <- function(datTR, datTE, t2e.var, event.var, varX, id_col = "eid", H = 10, nmax = 100000) {
		if (!id_col %in% names(datTR) || !id_col %in% names(datTE)) id_col <- names(datTR)[1]
		tt1 <- as.numeric(datTR[[t2e.var]]); ev1 <- as.integer(datTR[[event.var]] %in% c(1, "1", TRUE, "TRUE"))
		tt2 <- as.numeric(datTE[[t2e.var]]); ev2 <- as.integer(datTE[[event.var]] %in% c(1, "1", TRUE, "TRUE"))
		datTR10 <- datTR[, c(id_col, varX), drop = FALSE]; datTR10$y10 <- as.integer(ev1 == 1 & tt1 <= H)
		keep_te <- which(is.finite(tt2) & is.finite(ev2) & (tt2 >= H | ev2 == 1))
		datTE10 <- datTE[keep_te, c(id_col, varX), drop = FALSE]; datTE10$y10 <- if (length(keep_te)) as.integer(ev2[keep_te] == 1 & tt2[keep_te] <= H) else integer(0)
		if (nrow(datTR10) < 200 || nrow(datTE10) < 50 || length(unique(datTR10$y10)) < 2 || length(unique(datTE10$y10)) < 2) return(NULL)
		var_cat <- varX[sapply(datTR10[, varX, drop = FALSE], \(x) is.factor(x) || is.character(x))]; var_num <- setdiff(varX, var_cat)
		if (length(var_num)) { datTR10[var_num] <- lapply(datTR10[var_num], as.numeric); datTE10[var_num] <- lapply(datTE10[var_num], as.numeric) }
		if (length(var_cat)) { datTR10[var_cat] <- lapply(datTR10[var_cat], as.character); datTE10[var_cat] <- lapply(datTE10[var_cat], as.character) }
		idx1 <- seq_len(nrow(datTR10))
		if (nrow(datTR10) > nmax) {
			i1 <- which(datTR10$y10 == 1L); i0 <- which(datTR10$y10 == 0L); if (length(i1) < 50 || length(i0) < 50) return(NULL)
			n_each <- floor(nmax / 2); idx1 <- c(sample(i1, min(length(i1), n_each)), sample(i0, min(length(i0), n_each)))
		}
		list(datTR10 = datTR10, datTE10 = datTE10, keep_te = keep_te, var_cat = var_cat, var_num = var_num, idx1 = idx1)
	}
	vv <- .pred_pick_vars(datTR, varX); rec <- prep10y(datTR, datTE, t2e.var, event.var, vv, id_col = opt$id_col %||% "eid", H = opt$H %||% 10, nmax = opt$nmax %||% 100000)
	if (is.null(rec)) return(list(pred = rep(NA_real_, nrow(datTE)), auc = NA_real_, cidx = tryCatch(fcidx(Surv(datTE[[t2e.var]], datTE[[event.var]]), rep(NA_real_, nrow(datTE))), error = function(e) NA_real_), n_sel = length(vv)))
	dtr <- rec$datTR10; dte <- rec$datTE10; keep_te <- rec$keep_te; idx1 <- rec$idx1; var_cat <- rec$var_cat; var_num <- rec$var_num
	if (method == "tabpfn") {
		TabPFN$fit(as.matrix(dtr[idx1, vv, drop = FALSE]), as.integer(dtr$y10[idx1])); p2 <- TabPFN$predict_proba(as.matrix(dte[, vv, drop = FALSE])); p_hat <- if (is.matrix(p2) && ncol(p2) >= 2) as.numeric(p2[, 2]) else as.numeric(p2[, 1])
		if (exists("torch", inherits = TRUE)) try(torch$cuda$empty_cache(), silent = TRUE)
	} else {
		ntr <- nrow(dtr[idx1, , drop = FALSE]); vidx <- sample.int(ntr, max(1L, floor((opt$valid_frac %||% 0.15) * ntr)))
		pred_res <- ft_transformer_predict(train_in = dtr[idx1, c(vv, "y10"), drop = FALSE][-vidx, , drop = FALSE], valid_in = dtr[idx1, c(vv, "y10"), drop = FALSE][vidx, , drop = FALSE], test_df = dte[, vv, drop = FALSE], event_col = "y10", var_num = var_num, var_cat = var_cat)
		pc <- grep("_1_probability$", names(pred_res), value = TRUE); p_hat <- if (length(pc)) as.numeric(pred_res[[pc[1]]]) else as.numeric(pred_res[[ncol(pred_res)]])
		if (exists("torch", inherits = TRUE)) try(torch$cuda$empty_cache(), silent = TRUE)
	}
	pred <- rep(NA_real_, nrow(datTE)); pred[keep_te] <- p_hat
	auc <- tryCatch(as.numeric(pROC::auc(pROC::roc(dte$y10, p_hat, quiet = TRUE))), error = function(e) NA_real_)
	cidx <- tryCatch(fcidx(Surv(datTE[[t2e.var]], datTE[[event.var]]), pred), error = function(e) NA_real_)
	list(pred = pred, auc = auc, cidx = as.numeric(cidx), n_sel = length(vv))
}

pred_risk_by_fold <- function(fold, dat1, t2e.var, event.var, vars.basic.use0, models_list, id.var = "eid", verbose = TRUE) {
	tr <- dat1[dat1$fold_id != fold, , drop = FALSE]; te <- dat1[dat1$fold_id == fold, , drop = FALSE]
	if (sum(tr[[event.var]] == 1, na.rm = TRUE) < 30 || sum(te[[event.var]] == 1, na.rm = TRUE) < 5 || length(unique(te[[event.var]])) < 2) return(NULL)
	vb <- intersect(vars.basic.use0, names(tr))
	bind_rows(lapply(names(models_list), function(M) {
		c <- models_list[[M]]; if (verbose) message(" Fold ", fold, " | Model: ", M, " | ", Sys.time())
		vx <- .pred_pick_vars(tr, intersect(c$varX, names(tr))); if (!length(vx)) { if (verbose) message(" Skip: ", M); return(NULL) }
		pred_fun <- switch(c$fx, RE = RE_pred, ML = ML_pred, DL = DL_pred, NULL); if (is.null(pred_fun)) return(NULL)
		vb_model <- if (!is.null(c$vars.basic)) c$vars.basic else vb
		out <- tryCatch(pred_fun(tr, te, t2e.var, event.var, vx, vb_model, c$method %||% "cox", c$opt %||% list()), error = function(e) { if (verbose) message(" ", M, ": ", e$message); NULL })
		if (is.null(out) || is.null(out$pred)) return(NULL)
		dte <- if (!is.null(out$test_index)) te[out$test_index, , drop = FALSE] else .pred_rebuild_test(te, c(t2e.var, event.var, vx, vb_model))
		n0 <- min(nrow(dte), length(out$pred)); if (n0 == 0) return(NULL)
		res <- data.frame(eid = if (id.var %in% names(dte)) dte[[id.var]][1:n0] else seq_len(n0), fold = fold, model = M, fx = c$fx, method = c$method, time = dte[[t2e.var]][1:n0], event = dte[[event.var]][1:n0], risk = as.numeric(out$pred)[1:n0], cidx = out$cidx %||% NA_real_, n_sel = out$n_sel %||% NA_real_)
		rm(out, dte); invisible(gc()); res
	}))
}


# 🚩 Prediction summary + figures
.pred_lab <- function(x) dplyr::recode(x, clinical = "Clinical", protein = "Proteins", combined = "Combined", .default = x)
.std01 <- function(x) as.numeric(scale(x))
get_risk_df <- function(dat, risk_col, time_col, event_col) dat %>% transmute(risk = .data[[risk_col]], time = .data[[time_col]], event = as.integer(.data[[event_col]])) %>% filter(is.finite(risk), is.finite(time), !is.na(event))

summ_pred_cv <- function(pred.cv.long) pred.cv.long %>%
	group_by(model, fold, fx, method) %>% summarise(cidx = first(cidx), n_sel = first(n_sel), .groups = "drop") %>%
	group_by(model, fx, method) %>% summarise(cidx_median = median(cidx, na.rm = TRUE), cidx_mean = mean(cidx, na.rm = TRUE), cidx_sd = sd(cidx, na.rm = TRUE), cidx_se = cidx_sd / sqrt(sum(!is.na(cidx))), n_sel_median = median(n_sel, na.rm = TRUE), .groups = "drop") %>%
	mutate(model_lab = factor(.pred_lab(model), levels = .pred_lab(c("clinical", "protein", "combined")))) %>% arrange(desc(cidx_median))

pred_make_risk_dat <- function(dat1, pred.cv.long, Y, id.var = "eid") dat1 %>%
	dplyr::select(all_of(c(id.var, paste0(Y, c(".t2e", ".Yt2e"))))) %>%
	left_join(pred.cv.long %>% dplyr::select(all_of(id.var), model, risk) %>% mutate(model = paste0(Y, ".", model, ".cv.cox.risk")) %>% pivot_wider(names_from = model, values_from = risk), by = id.var)

pred_metric_tab <- function(dat, Y, models = c("clinical", "protein", "combined"), fprs = c(.01, .05, .10)) {
	time_col <- paste0(Y, ".t2e"); event_col <- paste0(Y, ".Yt2e")
	purrr::map_df(models, function(m) {
		risk_col <- paste0(Y, ".", m, ".cv.cox.risk"); df <- get_risk_df(dat, risk_col, time_col, event_col)
		roc0 <- pROC::roc(df$event, df$risk, levels = c(0, 1), quiet = TRUE, ci = TRUE)
		dr <- purrr::map_df(fprs, \(fp) { thr <- quantile(df$risk[df$event == 0], 1 - fp, na.rm = TRUE); tibble(FPR = fp, DR = mean(df$risk[df$event == 1] >= thr, na.rm = TRUE), PPV = mean(df$event[df$risk >= thr] == 1, na.rm = TRUE)) })
		fit1 <- tryCatch(coxph(Surv(time, event) ~ scale(risk), data = df), error = function(e) NULL)
		hr1 <- if (is.null(fit1)) tibble(HR = NA_real_, HR_l = NA_real_, HR_u = NA_real_, P = NA_real_) else broom::tidy(fit1, exponentiate = TRUE, conf.int = TRUE) %>% transmute(HR = estimate, HR_l = conf.low, HR_u = conf.high, P = p.value)
		tibble(model = m, model_lab = .pred_lab(m), AUC = as.numeric(roc0$auc), AUC_l = as.numeric(roc0$ci)[1], AUC_u = as.numeric(roc0$ci)[3], N = nrow(df), Events = sum(df$event == 1)) %>% bind_cols(hr1) %>% tidyr::crossing(dr)
	})
}

pred_vars_risk <- function(dat1, t2e.var, event.var, models_list, id.var = "eid", s = NULL) {
	purrr::map_df(names(models_list), function(M) {
		c <- models_list[[M]]
		vx <- .pred_pick_vars(dat1, intersect(c$varX, names(dat1)))
		if(!length(vx) || c$fx != "RE") return(NULL)

		need <- c(t2e.var, event.var, vx)
		d <- dat1[complete.cases(dat1[, need, drop = FALSE]), need, drop = FALSE]
		if(nrow(d) < 30 || sum(d[[event.var]] == 1, na.rm = TRUE) < 10) return(NULL)

		for(v in vx) if(is.character(d[[v]]) || is.factor(d[[v]])) d[[v]] <- factor(d[[v]])
		x <- model.matrix(as.formula(paste0("~", paste(paste0("`", vx, "`"), collapse = "+"))), d)[, -1, drop = FALSE]
		x <- x[, apply(x, 2, var) > 0, drop = FALSE]

		family <- c$opt$family %||% "cox"
		lambda_rule <- s %||% c$opt$lambda_rule %||% "lambda.min"
		inner_cv_n <- c$opt$inner_cv_n %||% 10
		type.measure <- c$opt$type.measure %||% ifelse(family == "binomial", "auc", "deviance")

		y <- if(family == "binomial") as.integer(d[[event.var]] == 1) else survival::Surv(d[[t2e.var]], d[[event.var]])
		if(family == "binomial" && length(unique(y)) < 2) return(NULL)
		nfolds0 <- if(family == "binomial") max(3, min(inner_cv_n, as.numeric(table(y)))) else inner_cv_n

		fit <- glmnet::cv.glmnet(x, y, family = family, alpha = c$opt$alpha %||% 1,
		                          nfolds = nfolds0, type.measure = type.measure)

		lambda0 <- fit[[lambda_rule]]
		bb <- as.matrix(coef(fit, s = lambda0))
		out <- tibble(term = rownames(bb), beta = as.numeric(bb[, 1])) %>% filter(beta != 0)
		if(!nrow(out)) return(NULL)

		n_selected0 <- sum(out$term != "(Intercept)")

		out %>%
			mutate(model = M, model_lab = .pred_lab(M), fx = c$fx, method = c$method,
			       family = family, lambda_rule = lambda_rule, lambda = lambda0,
			       N = nrow(d), Events = sum(d[[event.var]] == 1, na.rm = TRUE),
			       n_input = ncol(x), n_selected = n_selected0,
			       exp_beta = exp(beta), abs_beta = abs(beta),
			       rank = rank(-abs_beta, ties.method = "first")) %>%
			arrange(model, rank) %>%
			dplyr::select(model, model_lab, fx, method, family, lambda_rule, rank, term, beta,
			       exp_beta, abs_beta, lambda, N, Events, n_input, n_selected)
	})
}

pred_vars_risk_cv <- function(dat1, t2e.var, event.var, models_list, fold.var = "fold_id", s = NULL) {
	purrr::map_df(sort(unique(dat1[[fold.var]])), function(fold) {
		tr <- dat1[dat1[[fold.var]] != fold, , drop = FALSE]
		pred_vars_risk(tr, t2e.var, event.var, models_list, s = s) %>% mutate(fold = fold, .before = model)
	})
}

summ_pred_vars_cv <- function(pred.vars.cv) pred.vars.cv %>%
	group_by(model, model_lab, fx, method, term) %>%
	summarise(n_folds = n_distinct(fold), beta_mean = mean(beta, na.rm = TRUE), beta_sd = sd(beta, na.rm = TRUE), beta_median = median(beta, na.rm = TRUE), beta_min = min(beta, na.rm = TRUE), beta_max = max(beta, na.rm = TRUE), select_freq = n_folds / max(pred.vars.cv$fold, na.rm = TRUE), HR_mean = exp(beta_mean), abs_beta_mean = abs(beta_mean), .groups = "drop") %>%
	group_by(model) %>% mutate(rank = rank(-abs_beta_mean, ties.method = "first")) %>% ungroup() %>% arrange(model, rank)

pred_delong_tab <- function(dat, Y, ref = "clinical", models = c("protein", "combined")) {
	time_col <- paste0(Y, ".t2e"); event_col <- paste0(Y, ".Yt2e")
	roc_ref <- pROC::roc(dat[[event_col]], dat[[paste0(Y, ".", ref, ".cv.cox.risk")]], levels = c(0, 1), quiet = TRUE)
	purrr::map_df(models, function(m) {
		roc_m <- pROC::roc(dat[[event_col]], dat[[paste0(Y, ".", m, ".cv.cox.risk")]], levels = c(0, 1), quiet = TRUE)
		tibble(compare = paste0(m, " vs ", ref), AUC_ref = as.numeric(roc_ref$auc), AUC_model = as.numeric(roc_m$auc), delta_AUC = AUC_model - AUC_ref, P_delong = tryCatch(pROC::roc.test(roc_ref, roc_m, method = "delong")$p.value, error = function(e) NA_real_))
	})
}

# Plot below.
plot_risk <- function(dat1, Y_time, Y_event, X1, X2, covs = vars.basic, method = "10years", t0 = 10, X1lab = NULL, X2lab = NULL, group_bgcolor = FALSE, title = NULL, leg_position = "topleft", leg_direction = "vertical", tab = FALSE) {
	need_vars <- c(Y_time, Y_event, X1, X2)
	miss_vars <- setdiff(need_vars, names(dat1))
	if (length(miss_vars)) stop("plot_risk: missing column(s): ", paste(miss_vars, collapse = ", "), call. = FALSE)
	d <- dat1[!is.na(dat1[[X1]]) & !is.na(dat1[[X2]]), ]; lv1 <- levels(factor(d[[X1]])); lv2 <- levels(factor(d[[X2]]))
	if (!nrow(d)) stop("plot_risk: no rows after filtering non-missing ", X1, " and ", X2, call. = FALSE)
	d$x1 <- factor(d[[X1]], levels = lv1); d$x2 <- factor(d[[X2]], levels = lv2)
	if (method == "1000PYs") {
		res <- d %>% group_by(x1, x2) %>% summarise(e = sum(.data[[Y_event]], na.rm = TRUE), p = sum(.data[[Y_time]], na.rm = TRUE), .groups = "drop") %>% mutate(y = e / p * 1000, cl = pmax(0, (e - 1.96 * sqrt(e)) / p * 1000), cu = (e + 1.96 * sqrt(e)) / p * 1000, lbl = sprintf("%.1f", y), y_tit = "Rate (per 1000 PY)")
		if (group_bgcolor) bg <- d %>% dplyr::count(x1, x2, name = "n") %>% left_join(res %>% dplyr::select(x1, x2, y), by = c("x1", "x2")) %>% group_by(x1) %>% summarise(y = weighted.mean(y, n), .groups = "drop") %>% mutate(bg_lbl = sprintf("%.1f", y))
	} else {
		covs <- covs[covs %in% names(d)]; d <- d[complete.cases(d[, c("x1", "x2", Y_time, Y_event, covs), drop = FALSE]), ]
		if (!nrow(d)) stop("plot_risk: no complete rows for ", paste(c(X1, X2, Y_time, Y_event, covs), collapse = ", "), call. = FALSE)
		if (!any(d[[Y_event]] == 1, na.rm = TRUE)) stop("plot_risk: no events in ", Y_event, " after filtering", call. = FALSE)
		form <- as.formula(paste0("survival::Surv(", Y_time, ", ", Y_event, ") ~ x1 + x2", if (length(covs)) paste0(" + ", paste(covs, collapse = " + ")) else ""))
		fit <- survival::coxph(form, data = d); res <- expand.grid(x1 = factor(lv1, levels = lv1), x2 = factor(lv2, levels = lv2))
		if (length(covs)) for (v in covs) res[[v]] <- if (is.numeric(d[[v]])) mean(d[[v]], na.rm = TRUE) else names(which.max(table(d[[v]])))
		s <- summary(survival::survfit(fit, newdata = res), times = t0, extend = TRUE)
		res <- res %>% mutate(y = 1 - as.numeric(s$surv), cl = 1 - as.numeric(s$upper), cu = 1 - as.numeric(s$lower), lbl = sprintf("%.1f%%", y * 100), y_tit = sprintf("%g-year Risk", t0))
		if (group_bgcolor) bg <- d %>% dplyr::count(x1, x2, name = "n") %>% left_join(res %>% dplyr::select(x1, x2, y), by = c("x1", "x2")) %>% group_by(x1) %>% summarise(y = weighted.mean(y, n), .groups = "drop") %>% mutate(bg_lbl = sprintf("%.1f%%", y * 100))
	}
	if (tab) {
		count_tab <- d %>% group_by(x1, x2) %>% summarise(N = n(), Events = sum(.data[[Y_event]], na.rm = TRUE), .groups = "drop")
		cat("\n--- Group Statistics and Calculated Risks ---\n")
		print(left_join(count_tab, res %>% dplyr::select(x1, x2, Risk = lbl), by = c("x1", "x2"))); if (group_bgcolor) print(bg)
	}
	ymax <- max(res$cu, if (group_bgcolor) bg$y else 0, na.rm = TRUE); if (group_bgcolor) bg$y_bg <- ymax * 0.01; lp <- if (is.character(leg_position) && leg_position == "topleft") c(0.02, 0.98) else leg_position
	ggplot(res, aes(x1, y, fill = x2)) +
		{if (group_bgcolor) geom_col(data = bg, aes(x1, y), inherit.aes = FALSE, fill = "grey85", width = 0.92)} +
		geom_col(position = position_dodge(0.8), width = 0.7) +
		geom_errorbar(aes(ymin = cl, ymax = cu), width = 0.2, position = position_dodge(0.8)) +
		geom_text(aes(label = lbl, y = cu), vjust = -0.5, position = position_dodge(0.8), size = 3, fontface = "bold") +
		{if (group_bgcolor) geom_text(data = bg, aes(x1, y_bg, label = bg_lbl), inherit.aes = FALSE, fontface = "bold", size = 3, vjust = 0)} +
		scale_y_continuous(labels = if (method == "10years") scales::label_percent() else scales::label_comma(), expand = expansion(mult = c(0, 0.2))) +
		scale_fill_brewer(palette = "Set2") + theme_minimal() +
		labs(title = title, x = if (is.null(X1lab)) X1 else X1lab, y = res$y_tit[1], fill = if (is.null(X2lab)) X2 else X2lab) +
		theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(face = "bold"), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.position = lp, legend.justification = c(0, 1), legend.direction = leg_direction, legend.background = element_rect(fill = alpha("white", 0.5), color = NA)) +
		guides(fill = guide_legend(ncol = if (leg_direction == "vertical") 1 else NULL))
}

plot_cindex_bar <- function(pred.cv.sum, title = NULL, base_size = 16) {
	ggplot(pred.cv.sum, aes(model_lab, cidx_mean)) +
		geom_col(width = .62, fill = "grey35") +
		geom_errorbar(aes(ymin = cidx_mean - cidx_se, ymax = cidx_mean + cidx_se), width = .16, linewidth = .8) +
		geom_text(aes(label = sprintf("%.3f", cidx_mean), y = cidx_mean + cidx_se), vjust = -.7, fontface = "bold", size = 4.5) +
		coord_flip() + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .12))) +
		labs(title = title, x = NULL, y = "Cross-fitted AUC / C-index") + fig_theme(base_size)
}

plot_roc_multi <- function(dat, Y, models = c("clinical", "protein", "combined"), title = NULL, base_size = 16) {
	time_col <- paste0(Y, ".t2e"); event_col <- paste0(Y, ".Yt2e")
	df <- purrr::map_df(models, function(m) {
		d <- get_risk_df(dat, paste0(Y, ".", m, ".cv.cox.risk"), time_col, event_col); ro <- pROC::roc(d$event, d$risk, levels = c(0, 1), quiet = TRUE, ci = TRUE)
		tibble(FPR = 1 - ro$specificities, TPR = ro$sensitivities, model = m, model_lab = paste0(.pred_lab(m), "  AUC=", sprintf("%.3f", as.numeric(ro$auc))))
	})
	ggplot(df, aes(FPR, TPR, colour = model_lab)) + geom_line(linewidth = 1.15) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") + coord_equal() +
		labs(title = title, x = "False positive rate", y = "True positive rate", colour = NULL) + fig_theme(base_size) + theme(legend.position = c(.68, .18), legend.background = element_rect(fill = scales::alpha("white", .75), color = NA))
}

plot_density <- function(dat, risk_col, time_col, event_col, thr_prob = .95, title = NULL, xlab = "Risk score", base_size = 16) {
	df <- get_risk_df(dat, risk_col, time_col, event_col); thr <- quantile(df$risk[df$event == 0], thr_prob, na.rm = TRUE)
	DR <- mean(df$risk[df$event == 1] >= thr, na.rm = TRUE); FPR <- mean(df$risk[df$event == 0] >= thr, na.rm = TRUE)
	n1 <- sum(df$event == 1); n0 <- sum(df$event == 0); ymax <- max(density(df$risk)$y, na.rm = TRUE)
	ggplot(df, aes(risk, after_stat(density), fill = factor(event, c(0,1), c(sprintf("Controls, N=%s", format(n0, big.mark=",")), sprintf("Cases, N=%s", format(n1, big.mark=",")))))) +
		geom_density(alpha = .55, colour = NA) + geom_vline(xintercept = thr, colour = "grey25", linewidth = .7) +
		annotate("label", x = thr, y = ymax * .90, label = sprintf("DR = %.1f%%\nFPR = %.1f%%", 100 * DR, 100 * FPR), hjust = -.05, vjust = 1, size = 4, label.size = NA) +
		labs(title = title, x = xlab, y = "Density", fill = NULL) + fig_theme(base_size) + theme(legend.position = c(.78, .88), legend.background = element_rect(fill = scales::alpha("white", .75), color = NA))
}

plot_density_multi <- function(dat, Y, models = c("clinical", "protein", "combined"), title = NULL, base_size = 16) {
	time_col <- paste0(Y, ".t2e"); event_col <- paste0(Y, ".Yt2e")
	df <- purrr::map_df(models, \(m) get_risk_df(dat, paste0(Y, ".", m, ".cv.cox.risk"), time_col, event_col) %>% mutate(model_lab = factor(.pred_lab(m), .pred_lab(models)), risk = .std01(risk)))
	ggplot(df, aes(risk, after_stat(density), fill = factor(event, c(0,1), c("Controls", "Cases")))) + geom_density(alpha = .55, colour = NA) + facet_wrap(~model_lab, nrow = 1) + labs(title = title, x = "Standardized risk score", y = "Density", fill = NULL) + fig_theme(base_size) + theme(legend.position = "top")
}

plot_surv_quantile <- function(dat, risk_col, time_col, event_col, q = 5, title = NULL, base_size = 16) {
	qn <- q
	labs_q <- paste0("Q", seq_len(qn))
	df <- get_risk_df(dat, risk_col, time_col, event_col) %>%
		mutate(risk_q = ntile(risk, qn),
		       risk_q = factor(risk_q, levels = seq_len(qn), labels = labs_q))

	sf <- broom::tidy(survival::survfit(survival::Surv(time, event) ~ risk_q, data = df)) %>%
		mutate(risk_q = gsub("^risk_q=", "", strata),
		       risk_q = factor(risk_q, levels = labs_q),
		       cuminc = 1 - estimate)

	ggplot(sf, aes(time, cuminc, colour = risk_q)) +
		geom_step(linewidth = 1.05) +
		scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
		labs(title = title, x = "Follow-up time (years)", y = "Cumulative incidence", colour = "Risk quintile") +
		fig_theme(base_size) +
		theme(
			legend.position = c(.03, .97),
			legend.justification = c(0, 1),
			legend.background = element_rect(fill = scales::alpha("white", .75), color = NA)
		)
}

plot_hr_quantile <- function(dat, risk_col, time_col, event_col, q = 5, title = NULL, base_size = 16) {
	qn <- q
	df <- get_risk_df(dat, risk_col, time_col, event_col) %>%
		mutate(risk_q = ntile(risk, qn),
		       risk_q = factor(risk_q, levels = seq_len(qn), labels = paste0("Q", seq_len(qn))))

	fit <- survival::coxph(survival::Surv(time, event) ~ risk_q, data = df)
	d <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
		filter(grepl("^risk_q", term)) %>%
		mutate(lab = gsub("^risk_q", "", term),
		       lab = factor(lab, levels = paste0("Q", 2:qn)))

	ggplot(d, aes(estimate, lab)) +
		geom_vline(xintercept = 1, linetype = "dashed", colour = "grey65") +
		geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .18, linewidth = .8) +
		geom_point(size = 2.8) +
		scale_x_log10(
			breaks = c(1, 3, 10),
			expand = expansion(mult = c(.12, .08))
		) +
		labs(title = title, subtitle = "Reference = lowest risk quintile", x = "Hazard ratio", y = NULL) +
		fig_theme(base_size)
}

plot_inci_decile <- function(dat, risk_col, time_col, event_col, title = NULL, log_y = TRUE, base_size = 16) {
	df <- get_risk_df(dat, risk_col, time_col, event_col)
	d0 <- df %>% mutate(dec = ntile(risk, 10)) %>% group_by(dec) %>% summarise(N = n(), Events = sum(event), PY = sum(time), ir = Events / PY * 1000, .groups = "drop") %>% mutate(perc = (dec - .5) * 10, lab = sprintf("%d", Events))
	ggplot(d0, aes(perc, ir)) + geom_line(linewidth = 1) + geom_point(size = 2.8) + geom_text(aes(label = lab), vjust = -1, size = 3.8, fontface = "bold") + scale_x_continuous(breaks = seq(10, 90, 20), limits = c(0, 100)) + {if (log_y) scale_y_log10() else scale_y_continuous()} + labs(title = title, subtitle = "Numbers above points are events", x = "Risk score percentile", y = "Incidence rate per 1,000 PY") + fig_theme(base_size)
}

plot_dr_fpr <- function(dat, Y, models = c("clinical", "protein", "combined"), fprs = c(.01, .05, .10), title = NULL, base_size = 16) {
	d <- pred_metric_tab(dat, Y, models, fprs) %>% mutate(model_lab = factor(model_lab, .pred_lab(models)), FPR_lab = paste0("FPR ", 100 * FPR, "%"))
	ggplot(d, aes(model_lab, DR, fill = FPR_lab)) + geom_col(position = position_dodge(.75), width = .65) + geom_text(aes(label = sprintf("%.1f", 100 * DR)), position = position_dodge(.75), vjust = -.45, size = 3.8, fontface = "bold") + scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, .18))) + labs(title = title, x = NULL, y = "Detection rate", fill = NULL) + fig_theme(base_size) + theme(legend.position = "top")
}

plot_calibration_decile <- function(dat, risk_col, time_col, event_col, title = NULL, base_size = 16) {
	df <- get_risk_df(dat, risk_col, time_col, event_col) %>% mutate(dec = ntile(risk, 10)) %>% group_by(dec) %>% summarise(pred = mean(.std01(risk), na.rm = TRUE), obs = mean(event == 1), n = n(), e = sum(event), .groups = "drop")
	ggplot(df, aes(pred, obs)) + geom_line(linewidth = 1) + geom_point(size = 2.8) + geom_text(aes(label = e), vjust = -1, fontface = "bold", size = 3.8) + scale_y_continuous(labels = scales::label_percent(accuracy = 1)) + labs(title = title, subtitle = "Observed crude event proportion by predicted-risk decile; labels are events", x = "Mean standardized risk score", y = "Observed event proportion") + fig_theme(base_size)
}
