# 🚩 Packages, paths, and shared inputs
fcidx <- function(surv_obj, pred_v) {
	c_val <- concordance(surv_obj ~ pred_v)$concordance
	ifelse(c_val < 0.5, 1 - c_val, c_val)
}

vvar <- function(fit) tryCatch(stats::vcov(fit), error = function(e) vcov(fit))
.assoc_empty <- function(x, y_item, type, n_total = NA_integer_, msg = NA_character_) data.frame(term = x, estimate = NA_real_, std.error = NA_real_, conf.low = NA_real_, conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = n_total, N_event = NA_integer_, Outcome = y_item, Exposure = x, Model_Type = type, error = msg)

assoc_reg <- function(dat1, Xs, varX, Y, type = "bt", scale_X = TRUE) {
	if (length(Y) > 1) return(do.call(rbind, lapply(Y, \(yy) assoc_reg(dat1, Xs, varX, yy, type, scale_X))))
	y_item <- as.character(Y); vars <- varX
	do.call(rbind, lapply(Xs, function(x) {
		if (type == "t2e") need <- c(paste0(y_item, c(".t2e", ".Yt2e")), x, vars) else if (type == "t2e.tdc") {
			Yother <- sub("\\.t2e$", "", grep("\\.t2e$", vars, value = TRUE)); covars <- setdiff(vars, paste0(Yother, ".t2e")); need <- c("eid", paste0(y_item, c(".t2e", ".Yt2e")), x, covars, unlist(lapply(Yother, \(o) paste0(o, c(".t2e", ".Yt2e", ".Yr2e")))))
		} else need <- c(y_item, x, vars)
		need <- unique(need); d0 <- dat1[, need, drop = FALSE]
		if (type == "t2e.tdc") d0 <- d0[complete.cases(d0[, setdiff(names(d0), grep("\\.t2e$", names(d0), value = TRUE)), drop = FALSE]), , drop = FALSE] else d0 <- d0[complete.cases(d0), , drop = FALSE]
		if (nrow(d0) < 50) return(.assoc_empty(x, y_item, type, n_total = nrow(d0), msg = "too few complete cases"))
		if (is.numeric(d0[[x]])) { d0 <- d0[is.finite(d0[[x]]), , drop = FALSE]; if (scale_X) d0[[x]] <- as.numeric(scale(d0[[x]])) }
		for (v in vars) if (v %in% names(d0) && is.numeric(d0[[v]]) && !grepl("\\.t2e$", v)) d0 <- d0[is.finite(d0[[v]]), , drop = FALSE]
		formu <- as.formula(paste0(y_item, " ~ ", x, " + ", paste(vars, collapse = " + ")))
		tryCatch({
			if (type == "qt") { fit <- lm(formu, d0); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- cf[x, 1]; se <- cf[x, 2]; est <- b; cl <- b - 1.96 * se; ch <- b + 1.96 * se
			} else if (type == "bt") { fit <- glm(formu, d0, family = binomial()); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- cf[x, 1]; se <- cf[x, 2]; est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else if (type == "t2e") { f2 <- as.formula(paste0("survival::Surv(", y_item, ".t2e, ", y_item, ".Yt2e) ~ ", paste(c(x, vars), collapse = " + "))); fit <- survival::coxph(f2, d0, ties = "efron"); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- unname(cf[x, "coef"]); se <- unname(cf[x, "se(coef)"]); est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else if (type == "t2e.tdc") {
				Yother <- sub("\\.t2e$", "", grep("\\.t2e$", vars, value = TRUE)); covars <- setdiff(vars, paste0(Yother, ".t2e")); if (!"eid" %in% names(d0)) stop("eid not found")
				for (o in Yother) d0[[paste0(o, ".tdc")]] <- ifelse(d0[[paste0(o, ".Yr2e")]] == 1, 0, ifelse(d0[[paste0(o, ".Yt2e")]] == 1, d0[[paste0(o, ".t2e")]], NA))
				tstop <- d0[[paste0(y_item, ".t2e")]]; evt <- d0[[paste0(y_item, ".Yt2e")]]
				cmd_str <- paste0("survival::tmerge(data1 = d0, data2 = d0, id = eid, tstop = tstop, event = survival::event(tstop, evt)", if (length(Yother)) paste0(", ", paste0(paste0(Yother, ".tdc"), " = survival::tdc(", paste0(Yother, ".tdc"), ")", collapse = ", ")) else "", ")")
				tdc_dat <- eval(parse(text = cmd_str))
				f2 <- as.formula(paste0("survival::Surv(tstart, tstop, event) ~ ", paste(c(x, covars, paste0(Yother, ".tdc")), collapse = " + "))); fit <- survival::coxph(f2, tdc_dat, ties = "efron"); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- unname(cf[x, "coef"]); se <- unname(cf[x, "se(coef)"]); est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else if (type == "ordinal") { fit <- MASS::polr(formu, d0, Hess = TRUE, method = "logistic"); b0 <- coef(fit); j <- match(x, names(b0)); if (is.na(j)) stop("term dropped"); b <- unname(b0[j]); V <- vvar(fit); se <- sqrt(V[j, j]); est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else stop("unknown type")
			if (!is.finite(b) || !is.finite(se) || se <= 0) stop("non-finite coef/se")
			z <- b / se; p <- 2 * pnorm(abs(z), lower.tail = FALSE); Ntot <- if (type == "t2e.tdc") nrow(tdc_dat) else nrow(d0)
			Nev <- if (type == "bt") sum(d0[[y_item]] == 1, na.rm = TRUE) else if (type %in% c("t2e", "t2e.tdc")) sum(d0[[paste0(y_item, ".Yt2e")]] == 1, na.rm = TRUE) else NA_integer_
			data.frame(term = x, estimate = est, std.error = se, conf.low = cl, conf.high = ch, statistic = z, p.value = p, N_total = Ntot, N_event = Nev, Outcome = y_item, Exposure = x, Model_Type = type, error = NA_character_)
		}, error = function(e) .assoc_empty(x, y_item, type, n_total = nrow(d0), msg = conditionMessage(e)))
	}))
}

assoc_me <- function(dat1, M, varX, X, Y, type, sims = 100, boot = FALSE) {
	Mok <- paste0("`", M, "`")
	Xok <- paste0("`", X, "`")
	formu_m <- reformulate(c(Xok, varX), response = Mok)
	fit.X2M <- lm(formu_m, data = dat1, na.action = na.exclude)
	fit.X2M$call$formula <- formu_m
	if (type == "qt") {
		formu_y <- reformulate(c(Mok, Xok, varX), response = paste0("`", Y, "`"))
		fit.M2Y <- lm(formu_y, data = dat1, na.action = na.exclude)
		fit.M2Y$call$formula <- formu_y
	} else if (type == "bt") {
		formu_y <- reformulate(c(Mok, Xok, varX), response = paste0("`", Y, ".Yt2e`"))
		fit.M2Y <- glm(formu_y, data = dat1, family = binomial(), na.action = na.exclude)
		fit.M2Y$call$formula <- formu_y
	} else if (type == "cox") {
		rhs <- paste(c(Mok, Xok, varX), collapse = " + ")
		formu_y <- as.formula(paste0("survival::Surv(`", Y, ".t2e`, `", Y, ".Yt2e`) ~ ", rhs))
		fit.M2Y <- coxph(formu_y, data = dat1, na.action = na.exclude, x = TRUE, y = TRUE)
		fit.M2Y$call$formula <- formu_y
	}
	me <- mediate(model.m = fit.X2M, model.y = fit.M2Y, treat = X, mediator = M, boot = boot, sims = sims)
	return(data.frame(M = M, Y = Y, ACME = me$d0, ACME_p = me$d0.p, ADE = me$z0, ADE_p = me$z0.p,
		Total = me$tau.coef, Total_p = me$tau.p, Prop_med = me$n0, Prop_Med_p = me$n0.p, stringsAsFactors = FALSE
	))
}
