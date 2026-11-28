# 🚩 Association helpers
fcidx <- function(surv_obj, pred_v) {
	# Same risk-score contract as prediction.R, independent of source order.
	survival::concordance(surv_obj ~ pred_v, reverse = TRUE)$concordance
}

vvar <- function(fit) tryCatch(stats::vcov(fit), error = function(e) vcov(fit))
.assoc_empty <- function(x, y_item, type, n_total = NA_integer_, msg = NA_character_) data.frame(term = x, estimate = NA_real_, std.error = NA_real_, conf.low = NA_real_, conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_, N_total = n_total, N_event = NA_integer_, Outcome = y_item, Exposure = x, Model_Type = type, error = msg)


# The two disease status columns deliberately use different missingness rules:
# prevalent = (Yt2e=NA, Yr2e=1), incident = (Yt2e=1, Yr2e=NA),
# never diagnosed during observed follow-up = (Yt2e=0, Yr2e=0).
# Requiring both status columns to be complete would remove every case.
.assoc_tdc_prepare <- function(dat, y_item, x, covars, other) {
	target <- paste0(y_item, c(".t2e", ".Yt2e"))
	base <- unique(c("eid", target, x, covars))
	need <- unique(c(base, unlist(lapply(other, function(o) paste0(o, c(".t2e", ".Yt2e", ".Yr2e"))))))
	d <- dat[, need, drop = FALSE]
	reserved <- c("tstart", "tstop", "event", if (length(other)) c(paste0(other, ".tdc"), paste0(other, ".tdc_time")))
	if (any(reserved %in% names(d))) stop("Time-dependent model input uses a reserved interval column name")
	tt <- d[[target[1]]]; ev <- d[[target[2]]]
	if (!is.numeric(tt)) stop("Target follow-up must be numeric")
	keep <- complete.cases(d[, base, drop = FALSE]) & is.finite(tt) & tt > 0 & ev %in% c(0, 1)
	for (o in other) {
		yt <- d[[paste0(o, ".Yt2e")]]; yr <- d[[paste0(o, ".Yr2e")]]
		ot <- d[[paste0(o, ".t2e")]]
		if (!is.numeric(ot)) stop("Time-dependent disease follow-up must be numeric: ", o)
		prevalent <- yr %in% 1 & (is.na(yt) | yt %in% 0)
		incident <- yt %in% 1 & (is.na(yr) | yr %in% 0) & is.finite(ot) & ot > 0
		# A noncase must be observed through the target follow-up; a shorter
		# disease observation period is not evidence of remaining disease-free.
		never <- yt %in% 0 & yr %in% 0 & is.finite(ot) & ot >= tt
		known <- prevalent | incident | never
		known[is.na(known)] <- FALSE
		keep <- keep & known
		transition <- rep(NA_real_, nrow(d))
		transition[prevalent] <- 0
		transition[incident] <- ot[incident]
		d[[paste0(o, ".tdc_time")]] <- transition
	}
	keep[is.na(keep)] <- FALSE
	d[keep, , drop = FALSE]
}

.assoc_tdc_intervals <- function(d, y_item, other) {
	if (!nrow(d)) stop("No valid time-dependent survival observations")
	if (anyDuplicated(d$eid)) stop("Time-dependent survival input requires one row per eid")
	updates <- data.frame(eid = d$eid, .target_time = d[[paste0(y_item, ".t2e")]],
		.target_event = d[[paste0(y_item, ".Yt2e")]])
	# event()/tdc() are special operators evaluated by tmerge, not exported
	# functions in the survival namespace. Construct a call without text parsing.
	merge_call <- quote(survival::tmerge(data1 = d, data2 = updates, id = eid,
		tstop = .target_time, event = event(.target_time, .target_event)))
	for (o in other) {
		time_name <- paste0(o, ".tdc_time")
		updates[[time_name]] <- d[[time_name]]
		merge_call[[paste0(o, ".tdc")]] <- call("tdc", as.name(time_name))
	}
	eval(merge_call)
}


# 🚩 Regression associations
assoc_reg <- function(dat1, Xs, varX, Y, type = "bt", scale_X = TRUE) {
	if (length(Y) > 1) return(do.call(rbind, lapply(Y, \(yy) assoc_reg(dat1, Xs, varX, yy, type, scale_X))))
	y_item <- as.character(Y); vars <- varX
	do.call(rbind, lapply(Xs, function(x) {
		if (type == "t2e") need <- c(paste0(y_item, c(".t2e", ".Yt2e")), x, vars) else if (type == "t2e.tdc") {
			Yother <- sub("\\.t2e$", "", grep("\\.t2e$", vars, value = TRUE)); covars <- setdiff(vars, paste0(Yother, ".t2e")); need <- c("eid", paste0(y_item, c(".t2e", ".Yt2e")), x, covars, unlist(lapply(Yother, \(o) paste0(o, c(".t2e", ".Yt2e", ".Yr2e")))))
		} else need <- c(y_item, x, vars)
		need <- unique(need); d0 <- dat1[, need, drop = FALSE]
		if (type == "t2e.tdc") d0 <- .assoc_tdc_prepare(d0, y_item, x, covars, Yother) else d0 <- d0[complete.cases(d0), , drop = FALSE]
		if (nrow(d0) < 50) return(.assoc_empty(x, y_item, type, n_total = nrow(d0), msg = "too few complete cases"))
		if (is.numeric(d0[[x]])) { d0 <- d0[is.finite(d0[[x]]), , drop = FALSE]; if (scale_X) d0[[x]] <- as.numeric(scale(d0[[x]])) }
		for (v in vars) if (v %in% names(d0) && is.numeric(d0[[v]]) && !grepl("\\.t2e$", v)) d0 <- d0[is.finite(d0[[v]]), , drop = FALSE]
		formu <- as.formula(paste0(y_item, " ~ ", paste(c(x, vars), collapse = " + ")))
		tryCatch({
			if (type == "qt") { fit <- lm(formu, d0); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- cf[x, 1]; se <- cf[x, 2]; est <- b; cl <- b - 1.96 * se; ch <- b + 1.96 * se
			} else if (type == "bt") { fit <- glm(formu, d0, family = binomial()); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- cf[x, 1]; se <- cf[x, 2]; est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else if (type == "t2e") { f2 <- as.formula(paste0("survival::Surv(", y_item, ".t2e, ", y_item, ".Yt2e) ~ ", paste(c(x, vars), collapse = " + "))); fit <- survival::coxph(f2, d0, ties = "efron"); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- unname(cf[x, "coef"]); se <- unname(cf[x, "se(coef)"]); est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else if (type == "t2e.tdc") {
				tdc_dat <- .assoc_tdc_intervals(d0, y_item, Yother)
				f2 <- as.formula(paste0("survival::Surv(tstart, tstop, event) ~ ", paste(c(x, covars, if (length(Yother)) paste0(Yother, ".tdc")), collapse = " + "))); fit <- survival::coxph(f2, tdc_dat, ties = "efron"); cf <- coef(summary(fit)); if (!x %in% rownames(cf)) stop("term dropped"); b <- unname(cf[x, "coef"]); se <- unname(cf[x, "se(coef)"]); est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else if (type == "ordinal") { fit <- MASS::polr(formu, d0, Hess = TRUE, method = "logistic"); b0 <- coef(fit); j <- match(x, names(b0)); if (is.na(j)) stop("term dropped"); b <- unname(b0[j]); V <- vvar(fit); se <- sqrt(V[j, j]); est <- exp(b); cl <- exp(b - 1.96 * se); ch <- exp(b + 1.96 * se)
			} else stop("unknown type")
			if (!is.finite(b) || !is.finite(se) || se <= 0) stop("non-finite coef/se")
			z <- b / se; p <- 2 * pnorm(abs(z), lower.tail = FALSE); Ntot <- if (type == "t2e.tdc") length(unique(tdc_dat$eid)) else nrow(d0)
			Nev <- if (type == "bt") sum(d0[[y_item]] == 1, na.rm = TRUE) else if (type %in% c("t2e", "t2e.tdc")) sum(d0[[paste0(y_item, ".Yt2e")]] == 1, na.rm = TRUE) else NA_integer_
			data.frame(term = x, estimate = est, std.error = se, conf.low = cl, conf.high = ch, statistic = z, p.value = p, N_total = Ntot, N_event = Nev, Outcome = y_item, Exposure = x, Model_Type = type, error = NA_character_)
		}, error = function(e) .assoc_empty(x, y_item, type, n_total = nrow(d0), msg = conditionMessage(e)))
	}))
}
