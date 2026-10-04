# 🚩 c4.connect
# C4 connections, matched-budget panel validation and reconstruction.
# LE8_C4_MODE selects connect/panel/reconstruct; the dispatcher sets the job/action.
.c4_mode <- Sys.getenv("LE8_C4_MODE", unset = if (identical(Sys.getenv("LE8_PANEL_ACTION"), "reconstruct")) "reconstruct" else if (identical(Sys.getenv("LE8_JOB"), "c4_panel_validation") || toupper(Sys.getenv("LE8_FOCUS_FUNCTIONS_ONLY")) %in% c("TRUE", "1", "YES")) "panel" else "connect")
if (!.c4_mode %in% c("connect", "panel", "reconstruct")) stop("Unknown LE8_C4_MODE: ", .c4_mode)
if (.c4_mode %in% c("panel", "reconstruct")) {
# Focused, prospective test of LE8 supervision and use of prevalent (Yang) donors.  All risk models use
# incident training participants only. Yang enters proxy learning.
suppressPackageStartupMessages({
	fdir <- Sys.getenv("LE8_FDIR", unset = "/mnt/d/scripts/le8/f")
	source(file.path(fdir, "0.common.R"))
	source(file.path(fdir, "c1.correlate.R"))

	# c4 focus extensions LE8-first selection and transparent comparison of distinct assay panels.  Neither
	# function below accesses disease outcomes or validation measurements.
	focus_heterogeneity <- function(contrasts, boot) {
		if (!nrow(boot) || !nrow(contrasts))
			return(tibble())
		lo <- "Low baseline inflammation"
		hi <- "High baseline inflammation"
		specs <- contrasts |>
			filter(stratum %in% c(lo, hi)) |>
			distinct(model, reference, landmark, horizon)
		out <- map_dfr(seq_len(nrow(specs)), function(i) {
			s <- specs[i, ]
			q <- contrasts |>
				filter(model == s$model, reference == s$reference, landmark == s$landmark, stratum %in% c(lo, hi))
			if (!all(c(lo, hi) %in% q$stratum))
				return(tibble())
			a <- boot |>
				filter(model == s$model, landmark == s$landmark, stratum %in% c(lo, hi))
			b <- boot |>
				filter(model == s$reference, landmark == s$landmark, stratum %in% c(lo, hi))
			z <- inner_join(a, b, by = c("replicate", "stratum", "landmark", "horizon"), suffix = c("_a", "_b")) |>
				mutate(delta = AUC_a - AUC_b) |>
				select(replicate, stratum, delta) |>
				pivot_wider(names_from = stratum, values_from = delta)
			if (!all(c(lo, hi) %in% names(z)))
				return(tibble())
			v <- z[[hi]] - z[[lo]]
			v <- v[is.finite(v)]
			if (length(v) < 20)
				return(tibble())
			est <- q$delta_AUC[match(hi, q$stratum)] - q$delta_AUC[match(lo, q$stratum)]
			se <- sd(v)
			p <- if (is.finite(se) && se > 0)
				2 * pnorm(abs(est / se), lower.tail = FALSE) else NA_real_
			bind_cols(s, tibble(delta_high_minus_low = est, lo = unname(quantile(v, 0.025)), hi = unname(quantile(
				v,
				0.975
			)), p_heterogeneity = p, inference = "Frozen-fit exploratory difference in delta AUC; independent within-stratum bootstrap; not inflammatory causation"))
		})
		if (nrow(out))
			out$FDR <- p.adjust(out$p_heterogeneity, "BH")
		out
	}

	focus_balanced_panel <- function(membership, k, components) {
		z <- membership |>
			filter(selected %in% TRUE, is.finite(r1), is.finite(r2)) |>
			mutate(strength = pmin(abs(r1), abs(r2))) |>
			arrange(desc(strength), feature) |>
			distinct(feature, .keep_all = TRUE)
		queues <- lapply(components, function(cmp) z$feature[z$component == cmp])
		names(queues) <- components
		panel <- character()
		while (length(panel) < k && any(lengths(queues) > 0)) {
			for (cmp in components) if (length(queues[[cmp]]) && length(panel) < k) {
				panel <- c(panel, queues[[cmp]][1])
				queues[[cmp]] <- queues[[cmp]][ - 1]
			}
		}
		unique(panel)
	}

	focus_panel_audit <- function(designs, membership, components) {
		coverage <- map_dfr(designs, function(ds) {
			mm <- membership |>
				filter(cohort == ds$cohort, selected %in% TRUE, feature %in% ds$features)
			counts <- table(factor(mm$component, levels = components))
			tibble(
				model = ds$name, component = components, assays = as.integer(counts), assay_budget = length(ds$features),
				selection_cohort = ds$cohort, scope = "Replicated partial association labels, not proof of cell origin or intervention responsiveness"
			)
		})
		overlap <- map_dfr(seq_along(designs), function(i) {
			if (i == 1)
				return(tibble())
			map_dfr(seq_len(i - 1L), function(j) {
				a <- designs[[i]]
				b <- designs[[j]]
				if (a$budget != b$budget || a$budget == 0)
					return(tibble())
				tibble(model = a$name, reference = b$name, budget = a$budget, intersection = length(intersect(
					a$features,
					b$features
				)), union = length(union(a$features, b$features)), Jaccard = length(intersect(
					a$features,
					b$features
				)) / length(union(a$features, b$features)), identical_panel = setequal(a$features, b$features))
			})
		})
		list(panel_coverage = coverage, panel_overlap = overlap)
	}
})
# Stable result key used by existing saved panels and reports.
LE8_JOB <- "c4_panel_validation"
# Edit f/0.common.R; C4_FOCUS_BUDGETS=5,10,50 overrides one run.
C4_FOCUS_BUDGETS <- LE8_ASSAY_BUDGETS

# Screening and proxy discovery do not depend on assay budgets/bootstrap count.  Keep their cache keys
# separate so changing c(5,10,50) only refits the panels.
focus_stage_signature <- function(name, layer, inputs) {
	fi <- file.info(inputs)
	funcs <- c(
		"cox_scan", "focus_split", "make_outcome", "t2e", "filter_analysis_cohort", "le8_training_connection_set",
		"le8_select_phenotypes"
	)
	le8_hash_object(list(
		stage = name, layer = layer, outcome = Y, inputs = inputs, size = fi$size, mtime = fi$mtime,
		options = le8_analysis_options(), clinical = covs_use, basic = vars.basic, components = vars.le8, follow_end = date_follow_end,
		seed = SEED, proxy_assignment = "fixed Yin halves v1", requested_components = Sys.getenv(
			"C4_FOCUS_COMPONENTS",
			paste(vars.le8, collapse = ",")
		), clinical_default = unique(c(
			vars.basic, "smoke.pts", "sbp", "hba1c_ngsp",
			"bmi", "nonhdl"
		)), algorithms = lapply(funcs, function(n) deparse(get(n, mode = "function"))), cap = Sys.getenv("FINAL_CONNECTION_MAX_N",
			unset = "60000"
		)
	))
}

focus_split <- function(dat, evar, seed = SEED) {
	set.seed(seed)
	strata <- split(seq_len(nrow(dat)), as.character(dat[[evar]]), drop = TRUE)
	val <- unlist(lapply(strata, function(ii) sample(ii, max(1L, round(length(ii) * 0.2)))), use.names = FALSE)
	ifelse(seq_len(nrow(dat)) %in% val, "validation", "training")
}

focus_proxy_accuracy <- function(train, test, features, components, covars, model) {
	map_dfr(components, function(cmp) {
		tr <- train[is.finite(train[[cmp]]), , drop = FALSE]
		te <- test[is.finite(test[[cmp]]), , drop = FALSE]
		if (nrow(tr) < 100 || nrow(te) < 50)
			return(tibble())
		predict_one <- function(vars) {
			x <- le8_prepare_prediction_matrix(tr, te, vars)
			xx <- cbind(Intercept = 1, x$train)
			xt <- cbind(Intercept = 1, x$test)
			b <- lm.fit(xx, tr[[cmp]])$coefficients
			if (anyNA(b))
				return(rep(NA_real_, nrow(te)))
			drop(xt %*% b)
		}
		y <- te[[cmp]]
		mu <- mean(tr[[cmp]])
		den <- sum((y - mu) ^ 2)
		base <- predict_one(covars)
		full <- predict_one(unique(c(covars, features)))
		only <- predict_one(features)
		tibble(model, component = cmp, N = nrow(te), R2_omics = 1 - sum((y - only) ^ 2) / den, R2_basic = 1 - sum((y -
			base) ^ 2) / den, R2_basic_omics = 1 - sum((y - full) ^ 2) / den, delta_R2 = (sum((y - base) ^ 2) - sum((y -
			full) ^ 2)) / den, note = "Common Yin-trained held-out LE8 reconstruction; not intervention response")
	})
}

focus_inflammation <- function(train, test, features) {
	markers <- intersect(le8_csv_env("C4_FOCUS_INFLAMMATION", "GDF15,IL6,CRP"), features)
	score <- function(d, ref) {
		x <- sapply(markers, function(v) (as.numeric(d[[v]]) - ref[v, "center"]) / ref[v, "sd"])
		if (is.null(dim(x)))
			x <- matrix(x, ncol = length(markers))
		ans <- rowMeans(x)
		ans[rowSums(is.finite(x)) != length(markers)] <- NA_real_
		ans
	}
	if (length(markers) < 2L)
		return(list(group = rep("unavailable", nrow(test)), audit = tibble(status = "fewer than two prespecified inflammatory markers")))
	ref <- t(vapply(
		markers, function(v) c(center = mean(train[[v]], na.rm = TRUE), sd = sd(train[[v]], na.rm = TRUE)),
		numeric(2)
	))
	if (any(!is.finite(ref)) || any(ref[, "sd"] <= 0))
		return(list(group = rep("unavailable", nrow(test)), audit = tibble(status = "invalid training marker distribution")))
	tr <- score(train, ref)
	te <- score(test, ref)
	cut <- median(tr[is.finite(tr)])
	group <- ifelse(!is.finite(te), "unavailable", ifelse(te <= cut, "Low baseline inflammation", "High baseline inflammation"))
	list(group = group, audit = tibble(
		marker = markers, center = ref[, "center"], sd = ref[, "sd"], cutoff = cut,
		definition = "Complete-marker mean z score; training median threshold; NOT a non-inflammatory CAD subtype"
	))
}

focus_panel <- function(ranked, ys, k, kind, ys_fraction = 0.8) {
	r <- ranked[ranked %in% ys]
	if (kind == "NS")
		return(head(ranked, k))
	if (kind == "YS")
		return(head(r, k))
	# Fix the supervised allocation before seeing validation performance.
	nys <- ceiling(k * ys_fraction)
	if (length(r) < nys)
		return(character())
	unique(c(head(r, nys), head(setdiff(ranked, ys), k - nys)))
}

focus_contrasts <- function(metrics, boot, budgets) {
	empty <- tibble(
		stratum = character(), landmark = numeric(), horizon = numeric(), AUC_a = numeric(), AUC_b = numeric(),
		model = character(), reference = character(), contrast = character(), delta_AUC = numeric(), inference = character()
	)
	if (!nrow(metrics))
		return(empty)
	specs <- list()
	for (k in budgets) {
		for (cohort in c("Yin", "YinYang")) for (kind in c("YS", "YSplus", "YSbalanced")) specs[[length(specs) +
			1L]] <- c(paste(kind, cohort, k, sep = "_"), paste("NS", k, sep = "_"), "Supervision vs NS")
		for (kind in c("YS", "YSplus", "YSbalanced")) specs[[length(specs) + 1L]] <- c(paste(kind, "YinYang", k,
			sep = "_"
		), paste(kind, "Yin", k, sep = "_"), "Added Yang for proxy learning")
	}
	q <- function(x, p) if (sum(is.finite(x)) >= 20)
		as.numeric(quantile(x, p, na.rm = TRUE)) else NA_real_
	result <- map_dfr(specs, function(s) {
		a <- metrics |>
			filter(model == s[1])
		b <- metrics |>
			filter(model == s[2])
		if (!nrow(a) || !nrow(b))
			return(tibble())
		z <- inner_join(a |>
			select(stratum, landmark, horizon, AUC_a = AUC), b |>
			select(stratum, landmark, horizon, AUC_b = AUC), by = c("stratum", "landmark", "horizon"))
		if (nrow(boot)) {
			ba <- boot |>
				filter(model == s[1])
			bb <- boot |>
				filter(model == s[2])
			zz <- inner_join(ba, bb, by = c("replicate", "stratum", "landmark", "horizon"), suffix = c("_a", "_b")) |>
				mutate(delta = AUC_a - AUC_b) |>
				group_by(stratum, landmark, horizon) |>
				summarise(
					delta_lo = q(delta, 0.025), delta_hi = q(delta, 0.975), valid = sum(is.finite(delta)),
					.groups = "drop"
				)
			z <- left_join(z, zz, by = c("stratum", "landmark", "horizon"))
		}
		z |>
			mutate(model = s[1], reference = s[2], contrast = s[3], delta_AUC = AUC_a - AUC_b, inference = "Exploratory paired validation bootstrap; multiple budgets/strata, no external confirmation")
	})
	bind_rows(empty, result)
}

focus_plot <- function(metrics, contrasts, proxy, pillars, rd) {
	if (!nrow(metrics))
		return(invisible(NULL))
	a <- metrics |>
		filter(stratum == "All", landmark == 0, budget > 0) |>
		ggplot(aes(budget, AUC, color = paradigm)) +
		geom_line() +
		geom_point() +
		scale_x_continuous(breaks = sort(unique(metrics$budget[metrics$budget >
			0]))) +
		labs(
			title = "a. Same assay budget and incident validation cohort", x = "Measured proteins/metabolites",
			y = "10-year IPCW AUC", color = NULL
		) +
		theme_5c(9)
	clinical_auc <- metrics$AUC[metrics$model == "Clinical" & metrics$stratum == "All" & metrics$landmark == 0]
	if (length(clinical_auc))
		a <- a + geom_hline(yintercept = clinical_auc[1], linetype = 2, color = "grey45") + labs(subtitle = sprintf(
			"Dashed line: clinical covariates alone (AUC %.3f)",
			clinical_auc[1]
		))
	# Older saved results can contain a zero-column tibble when no pair exists.
	z <- if (nrow(contrasts))
		contrasts |>
			filter(stratum == "All", landmark == 0) else contrasts
	if (nrow(z)) {
		b <- ggplot(z, aes(delta_AUC, reorder(model, delta_AUC), color = contrast)) +
			geom_vline(
				xintercept = 0,
				linetype = 2
			) +
			geom_point() +
			labs(
				title = "b. Paired comparison with NS and with Yin", x = "Change in AUC",
				y = NULL, color = NULL
			) +
			theme_5c(9)
		if (all(c("delta_lo", "delta_hi") %in% names(z)))
			b <- b + geom_errorbar(aes(xmin = delta_lo, xmax = delta_hi), orientation = "y", width = 0.2)
	} else {
		b <- ggplot() +
			annotate("text", x = 0, y = 0, label = "No eligible paired comparisons\nSee design_status and fit_diagnostics") +
			labs(title = "b. Paired comparison with NS and with Yin") +
			theme_void()
	}
	c <- proxy |>
		ggplot(aes(component, model, fill = delta_R2)) +
		geom_tile() +
		scale_fill_gradient2(
			low = "#A64E59", mid = "white",
			high = "#24748A", midpoint = 0
		) +
		labs(
			title = "c. LE8 reconstruction in held-out participants", x = NULL,
			y = NULL, fill = "Incremental R²"
		) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 40, hjust = 1))
	d <- pillars |>
		ggplot(aes(component, n, fill = cohort)) +
		geom_col(position = "dodge") +
		labs(
			title = "d. Replicated proxies by pillar, learned in training",
			x = NULL, y = "Eligible features", fill = NULL
		) +
		theme_5c(9) +
		theme(axis.text.x = element_text(
			angle = 40,
			hjust = 1
		))
	p <- (a | b) / (c | d) + plot_annotation(caption = "Yang is used for LE8-proxy discovery only. All risk coefficients are fitted in incident training participants. No mechanistic subtype is inferred.")
	le8_queue_figure(p, file.path(rd, "c4.Fig16.equal_budget.png"), 19, 13, 220)
	le8_flush_figures(rd)
}

focus_write_outputs <- function(tables, rd) {
	for (nm in setdiff(names(tables), "PGS")) write_raw_csv(tables[[nm]], paste0("c4.focus.", nm, ".csv"), rd)
	focus_plot(tables$metrics, tables$contrasts, tables$proxy_accuracy, tables$pillar_counts, rd)
	sheets <- tables[setdiff(names(tables), "PGS")]
	if (is.list(tables$PGS))
		for (nm in names(tables$PGS)) sheets[[paste0("PGS_", nm)]] <- tables$PGS[[nm]]
	openxlsx::write.xlsx(sheets, file.path(rd, "c4.focus.xlsx"), overwrite = TRUE)
}

run_c4_panel_validation <- function(layer) {
	outdir <- if (layer == "protein")
		out.prot else out.met
	rd <- le8_job_dir(outdir, "c4_panel_validation")
	dir.create(rd, recursive = TRUE, showWarnings = FALSE)
	budget_values <- suppressWarnings(as.numeric(le8_csv_env("C4_FOCUS_BUDGETS", paste(C4_FOCUS_BUDGETS, collapse = ","))))
	if (any(!is.finite(budget_values)) || any(budget_values != floor(budget_values)))
		stop("Assay budgets must be whole numbers")
	budgets <- sort(unique(as.integer(budget_values)))
	solver <- Sys.getenv("C4_FOCUS_SOLVER", unset = "ridge")
	if (!solver %in% c("cox", "ridge"))
		stop("C4_FOCUS_SOLVER must be cox or ridge")
	B <- as.integer(le8_num_env("C4_FOCUS_BOOT", 200))
	fraction <- le8_num_env("C4_FOCUS_YS_FRACTION", 0.8)
	if (anyNA(budgets) || !length(budgets) || any(budgets < 1 | budgets > 500) || B < 0 || fraction <= 0 || fraction >
		1)
		stop("Invalid C4_FOCUS settings (budgets: 1..500)")
	biomfile <- file.path(indir, "Rdata", if (layer == "protein")
		"prot.rds" else "met.rds")
	inputs <- c(file.path(indir, "Rdata/all.rds"), biomfile)
	fi <- file.info(inputs)
	signature <- le8_hash_object(list(
		inputs = inputs, size = fi$size, mtime = fi$mtime, options = le8_analysis_options(),
		code = tools::md5sum(file.path(fdir, c("c4.connect.R", "c1.correlate.R", "0.common.R"))), settings = Sys.getenv()[grepl(
			"^(C4_FOCUS|FINAL_CONNECTION|C4_PGS|C1_PGS|DATE_FOLLOW_END)",
			names(Sys.getenv())
		)], budgets = budgets, B = B, seed = SEED, solver = solver
	))
	cache <- file.path(rd, "c4.validation.rds")
	if (!LE8_REPLACE && file.exists(cache)) {
		old <- readRDS(cache)
		if (LE8_REUSE_RESULTS) {
			le8_check_options(old)
			if (!all(c("metrics", "panel_members", "design") %in% names(old$tables)))
				stop("C4 completed result lacks presentation tables: ", cache)
		}
		if (LE8_REUSE_RESULTS || identical(old$signature, signature)) {
			message("C4 focus: matching completed result reused; regenerate presentation files")
			focus_write_outputs(old$tables, rd)
			return(invisible(old))
		}
	}
	stage <- function(name, expr) {
		f <- file.path(rd, paste0(name, ".rds"))
		old <- if (file.exists(f))
			readRDS(f) else NULL
		key <- focus_stage_signature(name, layer, inputs)
		if (!LE8_REPLACE && identical(old$signature, key))
			return(old$value)
		value <- force(expr)
		saveRDS(list(signature = key, value = value), f)
		value
	}
	message("C4 focus: loading current phenotypes and ", layer)
	biom <- if (layer == "protein")
		read_prot() else read_met()
	features <- setdiff(names(biom), "eid")
	clinical <- if (le8_custom_adjustment())
		le8_custom_covars else unique(c(vars.basic, "smoke.pts", "sbp", "hba1c_ngsp", "bmi", "nonhdl"))
	requested_components <- le8_csv_env("C4_FOCUS_COMPONENTS", paste(vars.le8, collapse = ","))
	if (length(requested_components) < 1 || length(requested_components) > 8)
		stop("Request 1–8 candidate LE8 domains; unsupported domains need not contribute assays")
	need <- unique(c(
		"eid", "ethnic.c", clinical, vars.basic, vars.le8, requested_components, "birth_date", "date_attend",
		"date_lost", "date_death", le8_y_date()
	))
	ph <- read_all(need) |>
		filter_analysis_cohort() |>
		make_outcome(Y)
	ph$eid <- as.character(ph$eid)
	biom$eid <- as.character(biom$eid)
	if (anyDuplicated(ph$eid) || anyDuplicated(biom$eid))
		stop("Duplicate phenotype/omic eid")
	dat <- inner_join(ph, biom, by = "eid")
	rm(ph, biom)
	invisible(gc())
	tvar <- paste0(Y, ".t2e")
	evar <- paste0(Y, ".Yt2e")
	bvar <- paste0(Y, ".b2e")
	yang <- dat[is.finite(dat[[bvar]]) & dat[[bvar]] <= 0, , drop = FALSE]
	yin <- dat[is.finite(dat[[tvar]]) & dat[[tvar]] > 0 & dat[[evar]] %in% c(0, 1), , drop = FALSE]
	if (length(intersect(yin$eid, yang$eid)))
		stop("Yin/Yang sets overlap")
	fold <- focus_split(yin, evar)
	train <- yin[fold == "training", , drop = FALSE]
	test <- yin[fold == "validation", , drop = FALSE]
	# The same Yin donor remains in the same discovery/replication half when Yang donors are added; resplitting
	# Yin would confound the comparison.
	set.seed(SEED + 411)
	train$.le8_proxy_half <- sample(rep(1 : 2, length.out = nrow(train)))
	yang$.le8_proxy_half <- sample(rep(1 : 2, length.out = nrow(yang)))
	rm(dat, yin)
	invisible(gc())
	if (nrow(train) < 1000 || nrow(test) < 100)
		stop("Insufficient omic participants")
	if (length(setdiff(clinical, names(train))))
		stop("Clinical covariates missing")
	missing_domains <- setdiff(requested_components, names(train))
	write_raw_csv(tibble(component = requested_components, status = ifelse(requested_components %in% missing_domains,
		"unavailable", "candidate; support determined within training"
	)), "c4.focus.domain_availability.csv", rd)
	components <- intersect(requested_components, names(train))
	basic <- intersect(vars.basic, names(train))
	message("C4 focus: training=", nrow(train), ", validation=", nrow(test), ", Yang=", nrow(yang))
	data.table::fwrite(bind_rows(
		tibble(eid = train$eid, role = "incident_training"), tibble(eid = test$eid, role = "incident_validation"),
		tibble(eid = yang$eid, role = "Yang_proxy_training")
	), file.path(rd, "c4.focus.roles.csv.gz"), compress = "gzip")
	screen <- stage("training_screen", cox_scan(train, features, clinical, Y, time_var = tvar, event_var = evar))
	ranked <- screen |>
		filter(is.finite(p.value)) |>
		arrange(p.value, term) |>
		pull(term)
	memberships <- list()
	sets <- list()
	for (cohort in c("Yin", "YinYang")) {
		message("C4 focus: learn ", cohort, " LE8 proxies")
		subdir <- file.path(rd, cohort)
		dir.create(subdir, showWarnings = FALSE)
		learn <- if (cohort == "Yin")
			train else bind_rows(train, yang)
		sets[[cohort]] <- stage(paste0("proxy_", cohort), le8_training_connection_set(
			learn, features, components,
			basic, subdir
		))
		f <- file.path(subdir, "final.connection_membership_training_only.csv")
		if (file.exists(f))
			memberships[[cohort]] <- as_tibble(fread(f)) |>
				mutate(cohort = cohort)
	}
	membership <- bind_rows(memberships)
	if (!ncol(membership))
		membership <- tibble(
			feature = character(), component = character(), cohort = character(), selected = logical(),
			r1 = numeric(), r2 = numeric(), FDR1 = numeric(), FDR2 = numeric(), specificity = numeric()
		)
	pillars <- membership |>
		filter(selected) |>
		count(cohort, component) |>
		complete(cohort = c("Yin", "YinYang"), component = components, fill = list(n = 0L))
	designs <- list(list(name = "Clinical", features = character(), budget = 0, paradigm = "Clinical", cohort = "Yin"))
	design_status <- list()
	for (k in budgets) for (kind in c("NS", "YS", "YSplus")) for (cohort in if (kind == "NS")
		"Yin" else c("Yin", "YinYang")) {
		fs <- focus_panel(ranked, sets[[cohort]], k, kind, fraction)
		nm <- if (kind == "NS")
			paste(kind, k, sep = "_") else paste(kind, cohort, k, sep = "_")
		design_status[[length(design_status) + 1L]] <- tibble(model = nm, budget = k, available = length(fs), status = if (length(fs) ==
			k)
			"ready" else "insufficient eligible assays")
		if (length(fs) == k)
			designs[[length(designs) + 1L]] <- list(name = nm, features = fs, budget = k, paradigm = if (kind ==
				"NS") kind else paste(kind, cohort), cohort = cohort)
	}
	# A separate LE8-first arm: balance available pillars, rank by replicated proxy strength only. Do not force
	# unvalidated sleep/activity proxies.
	for (k in budgets) for (cohort in c("Yin", "YinYang")) {
		fs <- focus_balanced_panel(membership |>
			filter(.data$cohort == .env$cohort), k, components)
		nm <- paste("YSbalanced", cohort, k, sep = "_")
		design_status[[length(design_status) + 1L]] <- tibble(model = nm, budget = k, available = length(fs), status = if (length(fs) ==
			k)
			"ready" else "insufficient replicated proxies")
		if (length(fs) == k)
			designs[[length(designs) + 1L]] <- list(name = nm, features = fs, budget = k, paradigm = paste(
				"YSbalanced",
				cohort
			), cohort = cohort)
	}
	panel_audit <- focus_panel_audit(designs, membership, components)
	inflammatory <- focus_inflammation(train, test, features)
	stratum <- list(All = seq_len(nrow(test)))
	if (layer == "protein")
		for (g in c("Low baseline inflammation", "High baseline inflammation")) stratum[[g]] <- which(inflammatory$group ==
			g)
	checkpoint_file <- file.path(rd, "c4.focus.prediction_checkpoint.rds")
	checkpoint <- if (file.exists(checkpoint_file))
		readRDS(checkpoint_file) else NULL
	if (!LE8_REPLACE && identical(checkpoint$signature, signature)) {
		message("C4 focus: reuse completed prediction checkpoint")
		tables <- checkpoint$tables
		met <- tables$metrics
		contrasts <- tables$contrasts
		pr <- tables$proxy_accuracy
	} else {
		metrics <- boot <- proxy <- members <- calibration <- decision <- coefficients <- preprocessing <- diagnostics <- baseline_hazards <- list()
		mi <- 0L
		for (ds in designs) {
			message("C4 focus: fit ", ds$name)
			obj <- le8_fit_budget_model(train, test, clinical, ds$features, tvar, evar, solver = solver, seed = SEED +
				27)
			diagnostics[[length(diagnostics) + 1L]] <- tibble(
				model = ds$name, status = obj$status, fit_method = obj$fit_method %||%
					solver, condition_number = obj$condition_number %||% NA_real_, lambda = obj$lambda %||% NA_real_,
				warnings = obj$warnings %||% ""
			)
			members[[length(members) + 1L]] <- tibble(model = ds$name, feature = if (length(ds$features))
				ds$features else NA_character_, budget = ds$budget, status = obj$status)
			if (obj$status != "ok")
				next
			coefficients[[length(coefficients) + 1L]] <- obj$coefficient |>
				mutate(model = ds$name, lp_center = obj$lp_center)
			preprocessing[[length(preprocessing) + 1L]] <- obj$preprocess |>
				mutate(model = ds$name)
			baseline_hazards[[length(baseline_hazards) + 1L]] <- as_tibble(obj$baseline_hazard) |>
				mutate(model = ds$name, lp_center = obj$lp_center)
			# Common reconstruction training cohort isolates panel selection from the effect of fitting
			# reconstruction coefficients in a different population.
			proxy[[length(proxy) + 1L]] <- focus_proxy_accuracy(train, test, ds$features, components, basic, ds$name)
			for (g in names(stratum)) for (L in c(0, 2, 5)) {
				ii <- stratum[[g]]
				ii <- ii[test[[tvar]][ii] > L]
				if (length(ii) < 100)
					next
				# Conditional risk from landmark L to baseline year 10, with a frozen fit.
				bh <- obj$baseline_hazard
				H <- function(t) {
					i <- which(bh$time <= t)
					if (length(i))
						bh$hazard[max(i)] else 0
				}
				risk <-  - expm1( - max(0, H(10) - H(L)) * exp(pmin(30, obj$lp[ii])))
				ev <- le8_evaluate_risk(test[[tvar]][ii] - L, test[[evar]][ii], risk, 10 - L, ds$name, ds$budget,
					ds$paradigm,
					B = B, seed = SEED + 1000 * match(g, names(stratum)) + as.integer(L * 10)
				)
				mi <- mi + 1L
				metrics[[mi]] <- ev$metrics |>
					mutate(stratum = g, landmark = L, actual_assays = obj$N_selected, fit_method = obj$fit_method)
				if (nrow(ev$bootstrap))
					boot[[mi]] <- ev$bootstrap |>
						mutate(stratum = g, landmark = L)
				if (nrow(ev$calibration))
					calibration[[mi]] <- ev$calibration |>
						mutate(stratum = g, landmark = L)
				if (nrow(ev$decision))
					decision[[mi]] <- ev$decision |>
						mutate(stratum = g, landmark = L)
			}
		}
		met <- bind_rows(metrics)
		bo <- bind_rows(boot)
		pr <- bind_rows(proxy)
		contrasts <- focus_contrasts(met, bo, budgets)
		heterogeneity <- focus_heterogeneity(contrasts, bo)
		tables <- list(
			metrics = met, contrasts = contrasts, heterogeneity = heterogeneity, proxy_accuracy = pr,
			panel_members = bind_rows(members), pillar_counts = pillars, membership = membership, inflammation_definition = inflammatory$audit,
			design_status = bind_rows(design_status), training_screen = screen, calibration = bind_rows(calibration),
			decision = bind_rows(decision), model_coefficients = bind_rows(coefficients), preprocessing = bind_rows(preprocessing),
			baseline_hazards = bind_rows(baseline_hazards), fit_diagnostics = bind_rows(diagnostics), panel_overlap = panel_audit$panel_overlap,
			panel_coverage = panel_audit$panel_coverage, design = tibble(item = c(
				"training", "YinYang", "test",
				"YSP allocation", "inference"
			), value = c(
				paste(nrow(train), "incident participants"), paste(
					nrow(yang),
					"additional prevalent donors for proxy learning only"
				), paste(nrow(test), "incident participants; shared 80/20 split algorithm with Final"),
				paste0("ceil(budget * ", fraction, ") YS; remainder disease-ranked outside YS"), "Post hoc hypothesis development after inspecting earlier results; frozen-fit internal bootstrap, requires new external validation"
			))
		)
		for (nm in names(tables)) write_raw_csv(tables[[nm]], paste0("c4.focus.", nm, ".csv"), rd)
		data.table::fwrite(bo, file.path(rd, "c4.focus.bootstrap.csv.gz"), compress = "gzip")
		saveRDS(list(signature = signature, tables = tables), checkpoint_file)
	}
	if (truthy(Sys.getenv("C4_FOCUS_PGS", unset = "TRUE"))) {
		# Save bridge results under the established C4 directory as well as this run.
		pd <- bind_rows(train, yang)
		set.seed(SEED + 811)
		half <- sample(rep(c("discovery", "replication"), length.out = nrow(pd)))
		mm <- membership |>
			filter(cohort == "YinYang") |>
			transmute(feature, primary_component = component, strict_YS = selected)
		prd <- le8_job_dir(outdir, "c4_connect")
		dir.create(prd, recursive = TRUE, showWarnings = FALSE)
		tables$PGS <- le8_pgs_bridge(pd, features, half, basic, screen, mm, layer, prd)
	}
	out <- list(signature = signature, meta = module_meta(layer, extra = list(status = "ok")), tables = tables)
	saveRDS(out, cache)
	focus_write_outputs(tables, rd)
	invisible(out)
}
run_c4_reconstruction <- function(layer) {
	B <- as.integer(Sys.getenv("C4_EXPLAIN_BOOT", "500"))
	if (!is.finite(B) || B < 50) stop("C4_EXPLAIN_BOOT must be at least 50")
	rd <- file.path(out.base, layer, "c4_connect")
	roles_file <- file.path(rd, "c4.focus.roles.csv.gz")
	panel_file <- file.path(rd, "c4.focus.panel_members.csv")
	if (!file.exists(roles_file) || !file.exists(panel_file))
		stop("Frozen C4 roles/panels missing for ", layer, "; run c4_panel_validation first. No new split is invented.")
	roles <- as.data.frame(fread(roles_file))
	panels <- as.data.frame(fread(panel_file))
	if (!all(c("eid", "role") %in% names(roles)) || anyDuplicated(roles$eid))
		stop("Invalid frozen roles")
	if (!all(c("model", "feature", "budget") %in% names(panels)))
		stop("Invalid frozen panel table")
	if ("status" %in% names(panels))
		panels <- panels[panels$status == "ok", , drop = FALSE]
	panels <- panels[!is.na(panels$feature) & nzchar(panels$feature), , drop = FALSE]
	wanted <- unique(c("eid", vars.basic, vars.le8, Sys.getenv("C4_EXPLAIN_GROUP", "")))
	ph <- as.data.frame(read_all(wanted))
	omics <- as.data.frame(if (layer == "prot")
		read_prot() else read_met())
	ph$eid <- as.character(ph$eid)
	omics$eid <- as.character(omics$eid)
	roles$eid <- as.character(roles$eid)
	if (anyDuplicated(ph$eid) || anyDuplicated(omics$eid))
		stop("Duplicate participant IDs")
	features <- unique(panels$feature)
	if (length(setdiff(features, names(omics))))
		stop("Frozen assay missing from current input: ", paste(setdiff(features, names(omics)), collapse = ","))
	dat <- merge(ph, omics[, c("eid", features), drop = FALSE], by = "eid", sort = FALSE)
	dat <- merge(dat, roles, by = "eid", sort = FALSE)
	expected <- roles$eid[roles$role %in% c("incident_training", "incident_validation")]
	if (!setequal(expected, dat$eid[dat$role %in% c("incident_training", "incident_validation")]))
		stop("Current inputs do not cover the original frozen train/test roster")
	train <- dat[dat$role == "incident_training", , drop = FALSE]
	test <- dat[dat$role == "incident_validation", , drop = FALSE]
	if (length(intersect(train$eid, test$eid)))
		stop("Training/validation overlap")
	train <- train[order(train$eid), , drop = FALSE]
	test <- test[order(test$eid), , drop = FALSE]
	groupcol <- Sys.getenv("C4_EXPLAIN_GROUP", "")
	if (nzchar(groupcol) && !groupcol %in% names(test))
		stop("Configured group column missing")
	components <- intersect(vars.le8, names(train))
	basic <- intersect(vars.basic, names(train))
	metrics <- contrasts <- boots <- list()
	for (cmp in components) {
		tr <- train[is.finite(train[[cmp]]), , drop = FALSE]
		te <- test[is.finite(test[[cmp]]), , drop = FALSE]
		if (nrow(tr) < 100 || nrow(te) < 50)
			next
		yy <- te[[cmp]]
		mu <- mean(tr[[cmp]])
		den <- sum((yy - mu) ^ 2)
		if (!is.finite(den) || den <= 0)
			next
		fit_once <- function(fs) {
			z <- le8_prepare_prediction_matrix(tr, te, unique(c(basic, fs)))
			a <- cbind(Intercept = 1, z$train)
			b <- cbind(Intercept = 1, z$test)
			cf <- lm.fit(a, tr[[cmp]])$coefficients
			if (any(!is.finite(cf)))
				stop("Singular reconstruction fit; ", cmp)
			drop(b %*% cf)
		}
		pp <- list(Clinical = fit_once(character()))
		for (model in unique(panels$model)) pp[[model]] <- fit_once(panels$feature[panels$model == model])
		err <- sapply(pp, function(p) (yy - p) ^ 2)
		baseerr <- (yy - mu) ^ 2
		for (model in colnames(err)) metrics[[length(metrics) + 1L]] <- tibble(model,
			component = cmp, N = nrow(te),
			budget = if (model == "Clinical")
				0L else unique(panels$budget[panels$model == model])[1], R2_basic_omics = 1 - sum(err[, model]) / den, RMSE = sqrt(mean(err[
,
				model
			])), scope = "Original frozen C4 validation roster; common Yin reconstruction training"
		)
		# Every comparison uses exactly the same resampled participants/groups.
		groups <- if (nzchar(groupcol))
			as.character(te[[groupcol]]) else te$eid
		if (anyNA(groups) || any(!nzchar(groups)))
			stop("Missing bootstrap group")
		ug <- unique(groups)
		blocks <- split(seq_len(nrow(te)), factor(groups, levels = ug))
		set.seed(SEED + 10000 + match(cmp, components))
		boot_r2 <- matrix(NA_real_, B, ncol(err), dimnames = list(NULL, colnames(err)))
		for (j in seq_len(B)) {
			ix <- unlist(blocks[sample.int(length(blocks), length(blocks), replace = TRUE)], use.names = FALSE)
			dn <- sum(baseerr[ix])
			if (dn > 0)
				boot_r2[j, ] <- 1 - colSums(err[ix, , drop = FALSE]) / dn
		}
		for (model in colnames(err)[grepl("^YS(plus|balanced)?_(Yin|YinYang)_[0-9]+$", colnames(err))]) {
			k <- as.integer(sub(".*_", "", model))
			ref <- paste0("NS_", k)
			if (!ref %in% colnames(err))
				next
			estimate <- (sum(err[, ref]) - sum(err[, model])) / den
			samples <- boot_r2[, model] - boot_r2[, ref]
			se <- sd(samples, na.rm = TRUE)
			p <- if (is.finite(se) && se > 0)
				2 * pnorm(abs(estimate / se), lower.tail = FALSE) else if (estimate == 0)
				1 else NA_real_
			qs <- quantile(samples, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
			contrasts[[length(contrasts) + 1L]] <- tibble(model,
				reference = ref, budget = k, component = cmp,
				N = nrow(te), delta_R2_vs_NS = estimate, lo = qs[1], hi = qs[2], p, valid_boot = sum(is.finite(samples)),
				uncertainty = "Paired group bootstrap; reconstructed frozen-panel fits; selection/training uncertainty excluded",
				grouping = if (nzchar(groupcol))
					groupcol else "participant; original split relatedness control unverified", interpretation = "Reconstruction fidelity, not causal intervention responsiveness"
			)
			boots[[length(boots) + 1L]] <- tibble(model,
				reference = ref, budget = k, component = cmp, replicate = seq_len(B),
				delta_R2 = samples
			)
		}
	}
	ct <- bind_rows(contrasts)
	if (nrow(ct))
		ct$FDR_all_proxy_contrasts <- p.adjust(ct$p, "BH")
	for (x in list(c("metrics", "c4.explain.metrics.csv"), c("contrasts", "c4.explain.contrasts.csv"), c(
		"bootstrap",
		"c4.explain.bootstrap.csv"
	))) {
		d <- switch(x[1],
			metrics = bind_rows(metrics),
			contrasts = ct,
			bootstrap = bind_rows(boots)
		)
		fwrite(d, file.path(rd, x[2]))
	}
	fwrite(tibble(role = c("roles", "panels"), path = c(roles_file, panel_file), md5 = unname(tools::md5sum(c(
		roles_file,
		panel_file
	))), note = "Reuses original roles and panels; no outcome or validation-based selection"), file.path(
		rd,
		"c4.explain.provenance.csv"
	))
	message("C4 explanation contrasts saved: ", nrow(ct), " / ", layer)
}

if (.c4_mode == "reconstruct") {
	for (layer in strsplit(BIOM, ",", fixed = TRUE)[[1]]) run_c4_reconstruction(layer)
} else {
	if (!truthy(Sys.getenv("LE8_FOCUS_FUNCTIONS_ONLY", unset = "FALSE"))) {
		if (prot_DO)
			le8_stage("C4/focus/protein", run_c4_panel_validation("protein"))
		if (met_DO)
			le8_stage("C4/focus/metabolite", run_c4_panel_validation("metabolite"))
	}
}
} else {
# C4: Connection — LE8-supervised proxy discovery, module assignment, associational mediation, interactions, and nonlinearity.
suppressPackageStartupMessages({
	fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(fdir, "0.common.R"))

	# c4 helpers
	# c4.connect.R
	# State-dependent omic co-abundance remodeling.
	#
	# This adapts the experimental paper's "reference network -> differential
	# edges -> convergent hubs" logic without claiming that plasma covariance is a
	# physical protein-protein interaction. A known-PPI file can annotate edges,
	# but AP-MS/variant perturbation evidence is required before calling them PPI
	# rewiring.

	C4_NETWORK_TOP <- as.integer(Sys.getenv("C4_NETWORK_TOP", unset = "80"))
	C4_NETWORK_DELTA <- as.numeric(Sys.getenv("C4_NETWORK_DELTA", unset = "0.20"))
	C4_NETWORK_MIN_N <- as.integer(Sys.getenv("C4_NETWORK_MIN_N", unset = "250"))
	C4_NETWORK_ANCHORS <- unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
		unset = "PCSK9,LPA,GDF15,NTPROBNP,MMP12"
	), ",", fixed = TRUE)[[1]]))

	read_c4_ppi_reference <- function() {
		f <- Sys.getenv("C4_PPI_FILE", unset = "")
		if (!nzchar(f) || !file.exists(f) || file.size(f) <= 0) return(tibble())
		x <- tryCatch(as_tibble(data.table::fread(f,
			showProgress = FALSE,
			check.names = FALSE
		)), error = function(e) tibble())
		if (ncol(x) < 2) return(tibble())
		names(x)[1 : 2] <- c("feature1", "feature2")
		x |>
			transmute(
				feature1 = as.character(feature1), feature2 = as.character(feature2),
				ppi_key = ifelse(feature1 < feature2, paste(feature1, feature2, sep = "||"),
					paste(feature2, feature1, sep = "||")
				)
			) |>
			distinct(ppi_key, .keep_all = TRUE)
	}

	select_c4_network_features <- function(
		disease, sets, modules, available,
		max_n = C4_NETWORK_TOP
	) {
		ranked <- if (nrow(disease) && all(c("term", "p.value") %in% names(disease))) disease |>
			filter(is.finite(p.value)) |>
			arrange(p.value) |>
			pull(term) else character()
		supervised <- unique(c(
			sets$YS_strict %||% character(), sets$YS %||% character(),
			modules$YS_core %||% character(), as.character(sets$membership$feature %||% character())
		))
		preferred <- unique(c(C4_NETWORK_ANCHORS, supervised, ranked))
		head(preferred[preferred %in% available], max_n)
	}

	c4_residual_matrix <- function(d, features, covars) {
		covars <- intersect(covars, names(d)) ; features <- intersect(features, names(d))
		if (!length(features)) return(list(x = matrix(numeric(), 0, 0), n = 0L))
		if (length(covars)) d <- d[complete.cases(d[, covars, drop = FALSE]), , drop = FALSE]
		if (nrow(d) < C4_NETWORK_MIN_N) return(list(x = matrix(numeric(), 0, 0), n = nrow(d)))
		x <- as.matrix(data.frame(lapply(d[, features, drop = FALSE], function(z) {
			z <- suppressWarnings(as.numeric(z)) ; med <- median(z, na.rm = TRUE)
			if (!is.finite(med)) med <- 0
			z[!is.finite(z)] <- med ; z
		}), check.names = FALSE))
		colnames(x) <- features ; keep <- apply(x, 2, function(z) is.finite(sd(z)) && sd(z) > 0)
		x <- x[, keep, drop = FALSE]
		if (!ncol(x)) return(list(x = x, n = nrow(d)))
		if (length(covars)) {
			mm <- model.matrix(reformulate(covars), d)
			x <- qr.resid(qr(mm), x)
		}
		x <- scale(x) ; x[!is.finite(x)] <- 0
		list(x = x, n = nrow(x))
	}

	c4_corr_edges <- function(ref, alt, comparison, ppi = tibble()) {
		common <- intersect(colnames(ref$x), colnames(alt$x))
		if (length(common) < 2 || ref$n < 4 || alt$n < 4) return(tibble())
		r0 <- cor(ref$x[, common, drop = FALSE], use = "pairwise.complete.obs")
		r1 <- cor(alt$x[, common, drop = FALSE], use = "pairwise.complete.obs")
		ij <- which(upper.tri(r0), arr.ind = TRUE)
		ans <- tibble(
			feature1 = common[ij[, 1]], feature2 = common[ij[, 2]],
			reference_r = r0[ij], state_r = r1[ij], delta_r = r1[ij] - r0[ij],
			reference_N = ref$n, state_N = alt$n, comparison = comparison
		) |>
			mutate(
				se_delta_z = sqrt(1 / pmax(reference_N - 3, 1) + 1 / pmax(state_N - 3, 1)),
				fisher_delta = atanh(cap(state_r, .999)) - atanh(cap(reference_r, .999)),
				z_diff = fisher_delta / se_delta_z,
				p.value = 2 * pnorm(abs(z_diff), lower.tail = FALSE), FDR = p.adjust(p.value, "BH"),
				ci_low = fisher_delta - 1.96 * se_delta_z, ci_high = fisher_delta + 1.96 * se_delta_z,
				edge_key = ifelse(feature1 < feature2, paste(feature1, feature2, sep = "||"),
					paste(feature2, feature1, sep = "||")
				),
				remodeling_class = case_when(
					sign(reference_r) != sign(state_r) & abs(reference_r) >= .15 & abs(state_r) >= .15 ~ "sign reversal",
					abs(reference_r) < .15 & abs(state_r) >= .30 ~ "gained co-abundance",
					abs(reference_r) >= .30 & abs(state_r) < .15 ~ "lost co-abundance",
					abs(state_r) > abs(reference_r) ~ "strengthened",
					TRUE ~ "weakened"
				),
				remodeled = FDR < .05 & abs(delta_r) >= C4_NETWORK_DELTA
			)
		if (nrow(ppi)) ans <- ans |>
			left_join(ppi |> transmute(
				edge_key = ppi_key,
				known_PPI_backbone = TRUE
			), by = "edge_key") |>
			mutate(known_PPI_backbone = coalesce(known_PPI_backbone, FALSE))
		else ans$known_PPI_backbone <- FALSE
		ans |> arrange(FDR, desc(abs(delta_r)))
	}

	make_c4_state_hubs <- function(edges, sets, modules) {
		if (!nrow(edges)) return(tibble())
		hubs <- bind_rows(
			edges |> transmute(comparison,
				feature = feature1, remodeled,
				delta_strength = abs(delta_r), known_PPI_backbone
			),
			edges |> transmute(comparison,
				feature = feature2, remodeled,
				delta_strength = abs(delta_r), known_PPI_backbone
			)
		) |>
			group_by(comparison, feature) |>
			summarise(
				remodeled_edges = sum(remodeled, na.rm = TRUE),
				remodeling_burden = sum(delta_strength[remodeled], na.rm = TRUE),
				known_PPI_edges = sum(remodeled & known_PPI_backbone, na.rm = TRUE), .groups = "drop"
			)
		sm <- as_tibble(sets$membership %||% tibble())
		if (nrow(sm)) sm <- sm |> select(any_of(c("feature", "set", "primary_component", "strict_YS")))
		mm <- as_tibble(modules$membership %||% tibble())
		if (nrow(mm)) mm <- mm |> select(any_of(c("feature", "module", "module_stability", "YS_core")))
		if (nrow(sm)) hubs <- hubs |> left_join(sm, by = "feature")
		if (nrow(mm)) hubs <- hubs |> left_join(mm, by = "feature")
		hubs |> arrange(comparison, desc(remodeled_edges), desc(remodeling_burden))
	}

	plot_c4_state_network <- function(net, layer, outdir) {
		edges <- as_tibble(net$edges %||% tibble()) ; hubs <- as_tibble(net$hubs %||% tibble())
		if (!nrow(edges)) {
			msg <- as.character((net$status$detail %||% "No estimable state comparison")[[1]])
			save_plot(blank_plot("State-dependent omic co-abundance remodeling", msg),
				"c4.Fig10.state_network_remodeling.png", 11, 7,
				outdir = outdir
			)
			save_plot(blank_plot("Differential co-abundance edge audit", msg),
				"c4.Fig11.state_network_edges.png", 11, 7,
				outdir = outdir
			)
			return(invisible(NULL))
		}
		top_edges <- edges |>
			filter(remodeled) |>
			group_by(comparison) |>
			slice_min(FDR, n = 70, with_ties = FALSE) |>
			ungroup()
		pa <- if (!nrow(top_edges)) blank_plot(
			"a. Differential co-abundance edges",
			"No edge passed both FDR and effect-size thresholds"
		) else
			ggplot(top_edges, aes(feature1, feature2, fill = delta_r, size =  - log10(pmax(FDR, 1e-30)))) +
				geom_point(shape = 21, color = "grey45") +
				facet_wrap( ~ comparison, scales = "free") +
				scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
				labs(
					title = "a. State-dependent co-abundance remodeling",
					subtitle = "Residual correlation changes relative to participants disease-free through year 10",
					x = NULL, y = NULL, fill = expression(Delta * r), size = expression( - log[10](FDR))
				) +
				theme_5c(7) +
				theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
		ht <- hubs |>
			group_by(comparison) |>
			slice_max(remodeled_edges, n = 15, with_ties = FALSE) |>
			ungroup() |>
			filter(remodeled_edges > 0) |>
			mutate(label = paste0(feature, " (", remodeled_edges, ")"))
		pb <- if (!nrow(ht)) blank_plot("b. Remodeling hubs") else
			ggplot(ht, aes(remodeling_burden, fct_reorder(label, remodeling_burden), fill = comparison)) +
				geom_col() +
				facet_wrap( ~ comparison, scales = "free_y") +
				labs(
					title = "b. Convergent remodeling hubs", x = expression(Sigma * abs(Delta * r)), y = NULL,
					fill = NULL
				) +
				theme_5c(8) +
				theme(legend.position = "none")
		save_plot(pa | pb, "c4.Fig10.state_network_remodeling.png", 18, 10, outdir = outdir)

		forest <- edges |>
			filter(remodeled) |>
			group_by(comparison) |>
			slice_min(FDR, n = 25, with_ties = FALSE) |>
			ungroup() |>
			mutate(
				edge = paste(feature1, feature2, sep = " — "),
				edge = fct_reorder(edge, delta_r)
			)
		pf <- if (!nrow(forest)) blank_plot("Differential co-abundance edge audit") else
			ggplot(forest, aes(fisher_delta, edge, color = remodeling_class, shape = known_PPI_backbone)) +
				geom_vline(xintercept = 0, color = "grey70") +
				geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = .08) +
				geom_point(size = 2.1) +
				facet_wrap( ~ comparison, scales = "free_y") +
				scale_shape_manual(
					values = c(`TRUE` = 17, `FALSE` = 16),
					labels = c(`TRUE` = "Known-PPI annotation", `FALSE` = "No supplied PPI annotation")
				) +
				labs(
					title = "Differential residual-correlation edge audit",
					subtitle = paste0(
						ifelse(layer == "protein", "Protein", "Metabolite"),
						" co-abundance is not a physical interaction assay; known-PPI status is annotation only"
					),
					x = expression(Delta * " Fisher z-transformed residual correlation"), y = NULL, color = NULL, shape = NULL
				) +
				theme_5c(8) +
				theme(legend.position = "bottom")
		save_plot(pf, "c4.Fig11.state_network_edges.png", 15, 11, outdir = outdir)
	}

	run_c4_state_network <- function(
		dat, disease, sets, modules, features, covars,
		tvar, evar, bvar, layer, rawdir, outdir
	) {
		cache <- file.path(rawdir, "c4.state_network.rds")
		selected <- select_c4_network_features(disease, sets, modules, features)
		old <- read_stage_cache(cache)
		if (is.list(old) && all(c("features", "edges", "hubs", "state_counts", "status") %in% names(old))) {
			plot_c4_state_network(old, layer, outdir) ; return(old)
		}
		distal <- is.finite(dat[[tvar]]) & dat[[tvar]] >= 10
		near <- dat[[evar]] == 1 & is.finite(dat[[tvar]]) & dat[[tvar]] > 0 & dat[[tvar]] <= 2
		prevalent <- is.finite(dat[[bvar]]) & dat[[bvar]] < 0
		groups <- list(
			`Near-diagnosis incident vs distal` = near,
			`Baseline-prevalent vs distal` = prevalent
		)
		ppi <- read_c4_ppi_reference()
		ref <- c4_residual_matrix(dat[distal, , drop = FALSE], selected, covars)
		alts <- lapply(groups, function(ii) c4_residual_matrix(dat[ii %in% TRUE, , drop = FALSE], selected, covars))
		counts <- tibble(
			state = c("Distal reference", names(groups)),
			N = c(ref$n, vapply(alts, `[[`, integer(1), "n"))
		)
		if (ref$n < C4_NETWORK_MIN_N || any(vapply(alts, `[[`, integer(1), "n") < C4_NETWORK_MIN_N)) {
			ans <- list(
				code_version = C4_CODE_VERSION, features = selected, edges = tibble(), hubs = tibble(),
				state_counts = counts, ppi_reference = ppi,
				status = tibble(status = "unavailable", detail = paste(
					"At least one state has fewer than",
					C4_NETWORK_MIN_N, "complete-covariate participants"
				))
			)
		} else {
			edges <- bind_rows(Map(function(alt, nm) c4_corr_edges(ref, alt, nm, ppi), alts, names(alts)))
			hubs <- make_c4_state_hubs(edges, sets, modules)
			ans <- list(
				code_version = C4_CODE_VERSION, features = selected, edges = edges, hubs = hubs,
				state_counts = counts, ppi_reference = ppi,
				status = tibble(status = "ok", detail = paste(
					length(selected),
					"features; residual co-abundance remodeling, not physical PPI rewiring"
				))
			)
		}
		write_stage_cache(ans, cache) ; plot_c4_state_network(ans, layer, outdir) ; ans
	}


	# c4.connect.R
	# Paper_t2dm Fig6f / design slide20: descriptive association ribbons.
	# Proxy display selection uses LE8 replication only, never disease P values.
	le8_c4_pathway_data <- function(sets, med, n_per_pillar = 2L) {
		d <- sets$YS_edges
		needed <- c('feature', 'component', 'r_rep', 'FDR_rep')
		if (!is.data.frame(d) || !all(needed %in% names(d)) || !nrow(d)) return(tibble())
		d <- d |>
			filter(is.finite(r_rep), is.finite(FDR_rep), FDR_rep < .05) |>
			arrange(component, desc(abs(r_rep)), feature) |>
			distinct(component, feature, .keep_all = TRUE) |>
			group_by(component) |>
			slice_head(n = n_per_pillar) |>
			ungroup()
		m <- sets$membership |> select(feature, disease_beta, disease_p, disease_FDR)
		if (anyDuplicated(m$feature)) stop('Duplicate biomarker membership in pathway plot')
		d <- d |>
			left_join(m, by = 'feature') |>
			mutate(
				risk = case_when(
					!is.finite(disease_p) ~ 'Disease association unavailable',
					!is.finite(disease_FDR) | disease_FDR >= .05 ~ 'No FDR disease association',
					disease_beta > 0 ~ 'Higher disease risk', TRUE ~ 'Lower disease risk'
				),
				left_sign = ifelse(r_rep > 0, 'Positive association', 'Inverse association'),
				right_sign = case_when(risk == 'Higher disease risk' ~ 'Positive association', risk == 'Lower disease risk' ~ 'Inverse association', TRUE ~ 'No supported association')
			)
		d$mediation_supported <- FALSE
		if (is.data.frame(med) && all(c('component', 'feature', 'FDR_indirect') %in% names(med))) {
			hits <- med |> filter(is.finite(FDR_indirect), FDR_indirect < .05)
			d$mediation_supported <- paste(d$component, d$feature) %in% paste(hits$component, hits$feature)
		}
		d
	}

	le8_plot_c4_pathways <- function(sets, med, layer, outdir) {
		n <- as.integer(Sys.getenv('C4_FLOW_PER_PILLAR', unset = '2'))
		if (is.na(n) || n < 1 || n > 10) stop('C4_FLOW_PER_PILLAR must be 1..10')
		d <- le8_c4_pathway_data(sets, med, n)
		rd <- le8_job_dir(outdir, 'c4_connect')
		write_raw_csv(d, 'c4.lifestyle_omics_risk_display.csv', rd)
		if (!nrow(d)) return(save_plot(blank_plot('LE8–omics–disease paths', 'No replicated LE8 proxy available'),
			'c4.Fig15.lifestyle_omics_risk.png', 18, 13,
			outdir = outdir
		))
		pillar_names <- c(
			diet = 'Diet', pa = 'Physical activity', smoke = 'Nicotine exposure', sleep = 'Sleep',
			bmi = 'Body mass index', nonhdl = 'Blood lipids', hba1c = 'Blood glucose', bp = 'Blood pressure'
		)
		d <- d |> mutate(pillar = coalesce(unname(pillar_names[sub('[.]pts$', '', component)]), component))
		d$.id <- seq_len(nrow(d))
		orders <- list(unique(d$pillar), unique(d$feature), intersect(c(
			'Higher disease risk', 'Lower disease risk',
			'No FDR disease association', 'Disease association unavailable'
		), d$risk))
		cols <- c('pillar', 'feature', 'risk') ; max_height <- nrow(d) + .45 * (max(lengths(orders)) - 1)
		nodes <- list() ; lanes <- list()
		for (stage in 1 : 3) {
			cursor <- (max_height - (nrow(d) + .45 * (length(orders[[stage]]) - 1))) / 2
			for (name in orders[[stage]]) {
				ids <- d$.id[d[[cols[stage]]] == name] ; h <- length(ids)
				nodes[[length(nodes) + 1L]] <- tibble(stage = stage, name = name, x = stage, ymin = cursor, ymax = cursor + h)
				lanes[[length(lanes) + 1L]] <- tibble(.id = ids, stage = stage, lo = cursor + seq_along(ids) - 1, hi = cursor + seq_along(ids))
				cursor <- cursor + h + .45
			}
		}
		nodes <- bind_rows(nodes) ; lanes <- bind_rows(lanes)
		ribbons <- list()
		for (stage in 1 : 2) for (id in d$.id) {
			a <- lanes |> filter(.id == id, stage == ! !stage) ; b <- lanes |> filter(.id == id, stage == ! !(stage + 1))
			t <- seq(0, 1, length.out = 60) ; s <- 3 * t ^ 2 - 2 * t ^ 3
			ribbons[[length(ribbons) + 1L]] <- tibble(
				path = paste(stage, id), x = c(stage + .06 + t * .88, rev(stage + .06 + t * .88)),
				y = c(a$lo + (b$lo - a$lo) * s, rev(a$hi + (b$hi - a$hi) * s)),
				sign = if (stage == 1) d$left_sign[id] else d$right_sign[id]
			)
		}
		nodes$fill <- '#C7CDCF'
		for (i in seq_len(nrow(nodes))) if (nodes$stage[i] == 1) {
			component <- d$component[match(nodes$name[i], d$pillar)]
			nodes$fill[i] <- if (component %in% names(cols_le8)) cols_le8[[component]] else '#5086A1'
		}
		nodes$fill[nodes$name == 'Higher disease risk'] <- '#DB726D'
		nodes$fill[nodes$name == 'Lower disease risk'] <- '#6D9ABD'
		nodes$label <- nodes$name
		supported <- unique(d$feature[d$mediation_supported]) ; nodes$label[nodes$stage == 2 & nodes$name %in% supported] <-
			paste0(nodes$label[nodes$stage == 2 & nodes$name %in% supported], ' *')
		p <- ggplot(bind_rows(ribbons), aes(x, y, group = path)) +
			geom_polygon(aes(fill = sign), alpha = .42, color = NA) +
			scale_fill_manual(values = c('Positive association' = '#E99194', 'Inverse association' = '#89B8D8', 'No supported association' = '#D5DADB'), drop = FALSE) +
			geom_rect(
				data = nodes, aes(xmin = x - .06, xmax = x + .06, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
				fill = nodes$fill, color = 'grey55', linewidth = .4
			) +
			geom_label(
				data = nodes |> filter(stage == 2), aes(x = x + .09, y = (ymin + ymax) / 2, label = label), inherit.aes = FALSE,
				hjust = 0, size = 3, linewidth = 0, fill = 'white', alpha = .93
			) +
			geom_text(data = nodes |> filter(stage == 1), aes(x = x - .10, y = (ymin + ymax) / 2, label = label), inherit.aes = FALSE, hjust = 1, size = 3.5) +
			geom_text(data = nodes |> filter(stage == 3), aes(x = x + .10, y = (ymin + ymax) / 2, label = label), inherit.aes = FALSE, hjust = 0, size = 3.5) +
			scale_y_reverse(breaks = NULL) +
			scale_x_continuous(breaks = NULL) +
			coord_cartesian(xlim = c(.45, 3.65), clip = 'off') +
			theme_void() +
			labs(
				title = paste('LE8-supervised', if (layer == 'protein') 'protein' else 'metabolite', 'paths to', Y),
				subtitle = paste0('Up to ', n, ' strongest replicated proxies per pillar; display selection does not use disease results'), fill = NULL, x = NULL, y = NULL,
				caption = paste(
					'LE8 score → omics → disease association. Higher LE8 scores indicate healthier status.',
					'Ribbon colour is the sign of each association; width counts displayed links, not effect size or proportion mediated.',
					'* At least one displayed component–biomarker path has mediation FDR < 0.05. These are associational paths, not proven causal mediation.',
					'All eligible proxies remain in the module tables.'
				)
			) +
			theme(
				plot.title = element_text(face = 'bold', size = 15), plot.subtitle = element_text(size = 11), legend.position = 'bottom',
				plot.caption = element_text(size = 9, hjust = 0), plot.margin = margin(20, 40, 20, 40),
				axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()
			)
		save_plot(p, 'c4.Fig15.lifestyle_omics_risk.png', 18, max(10, min(20, n_distinct(d$feature) * .55)), outdir = outdir)
		invisible(d)
	}


	# c4.connect.R
	# C4 brain structure connections: baseline omics -> imaging measurements.
	# Complete-case standardized linear models; no image imputation or causal claim.
	c4_img_manifest <- function(indir) {
		p <- file.path(indir, 'Rdata/img.dictionary.csv')
		if (!file.exists(p)) return(tibble())
		m <- as_tibble(data.table::fread(p))
		m <- m |>
			mutate(family = case_when(
				grepl('^aparc-Desikan_.*_area_', label) & !grepl('TotalSurface', label) ~ 'Cortical area',
				grepl('^aparc-Desikan_.*_volume_', label) & !grepl('Total', label) ~ 'Cortical volume',
				grepl('^aseg_[lr]h_volume_(Thalamus-Proper|Caudate|Putamen|Pallidum|Hippocampus|Amygdala|Accumbens-area|VentralDC)$', label) ~ 'Subcortical volume',
				grepl('dMRI_ProbtrackX_(FA|MD)_', label) ~ 'White matter FA/MD',
				field == 'p24486' ~ 'White matter hyperintensity',
				field == 'p26517' ~ 'Global grey matter', TRUE ~ NA_character_
			)) |>
			filter(!is.na(family))
		m |> mutate(adjust_TIV = !family %in% c('Cortical area', 'White matter FA/MD'))
	}
	c4_img_fit <- function(dat, feature, measure, covars, min_n = 100L) {
		columns <- unique(c(feature, measure, covars)) ; z <- as.data.frame(dat[, columns, drop = FALSE])
		ok <- complete.cases(z) & is.finite(z[[feature]]) & is.finite(z[[measure]])
		z <- z[ok, , drop = FALSE] ; n <- nrow(z)
		empty <- tibble(feature = feature, measure = measure, N = n, beta = NA_real_, SE = NA_real_, p = NA_real_, status = 'insufficient observations or variance')
		if (n < min_n || !is.finite(sd(z[[feature]])) || sd(z[[feature]]) == 0 || sd(z[[measure]]) == 0) return(empty)
		cv <- covars[vapply(z[covars], function(x) length(unique(x)) > 1, logical(1))]
		z[[feature]] <- as.numeric(scale(z[[feature]])) ; z[[measure]] <- as.numeric(scale(z[[measure]]))
		tryCatch(
			{
				fit <- lm(reformulate(c(feature, cv), response = measure), z)
				s <- coef(summary(fit)) ; if (!feature %in% rownames(s)) return(empty)
				tibble(feature = feature, measure = measure, N = n, beta = s[feature, 1], SE = s[feature, 2], p = s[feature, 4], status = 'ok')
			},
			error = function(e) {
				empty$status <- conditionMessage(e) ; empty
			}
		)
	}
	run_c4_imaging <- function(layer, disease = NULL, outdir = if (layer == 'protein') out.prot else out.met) {
		rawdir <- le8_job_dir(outdir, 'c4_connect') ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
		enabled <- truthy(Sys.getenv('C4_IMG_ENABLED', unset = 'TRUE'))
		missing_result <- function(reason) {
			status <- tibble(status = 'not run', detail = reason)
			write_raw_csv(status, 'c4.imaging_status.csv', rawdir)
			if (enabled) save_plot(blank_plot('Omics and brain structure', reason), 'c4.Fig14.imaging_atlas.png', 14, 9, outdir = outdir)
			list(status = status, associations = tibble(), fields = tibble())
		}
		if (!enabled) return(missing_result('C4_IMG_ENABLED=FALSE'))
		imgfile <- Sys.getenv('C4_IMG_FILE', unset = file.path(indir, 'Rdata/img.raw.rds'))
		if (!file.exists(imgfile)) return(missing_result('Prepare unimputed img.raw.rds with ukb/f/phe.R (phe.sh biom step)'))
		if (is.null(disease)) {
			f <- file.path(le8_job_dir(outdir, 'c1_correlate'), 'c1.res.rds')
			if (!file.exists(f)) return(missing_result('C1 outcome association table unavailable'))
			c1 <- readRDS(f) ; le8_check_options(c1) ; disease <- as_tibble(c1$association %||% c1$pwas_incident %||% c1$MWAS)
		}
		max_features <- as.integer(Sys.getenv('C4_IMG_MAX_FEATURES', unset = '100'))
		explicit <- trimws(strsplit(Sys.getenv('C4_IMG_FEATURES', unset = ''), ',', fixed = TRUE)[[1]])
		explicit <- explicit[nzchar(explicit)]
		if (length(explicit)) features <- unique(explicit) else {
			stopifnot(all(c('term', 'p.value') %in% names(disease)))
			features <- disease |>
				filter(is.finite(p.value), p.value < .05 / nrow(disease)) |>
				arrange(p.value) |>
				pull(term) |>
				unique() |>
				head(max_features)
		}
		if (!length(features)) return(missing_result('No C1 Bonferroni-significant feature; no data-driven threshold relaxation'))
		covars <- trimws(strsplit(Sys.getenv('C4_IMG_COVARS', unset = paste(if (le8_custom_adjustment()) le8_custom_covars else vars.basic, collapse = ',')), ',', fixed = TRUE)[[1]])
		datecol <- Sys.getenv('C4_IMG_OUTCOME_DATE', unset = le8_y_date())
		files <- c(
			imgfile, file.path(indir, 'Rdata/img.dictionary.csv'), file.path(indir, 'Rdata/img.visits.rds'),
			file.path(indir, 'Rdata/all.rds'), if (layer == 'protein') file.path(indir, 'Rdata/prot.rds') else file.path(indir, 'Rdata/met.rds')
		)
		signature <- list(
			version = 3L, features = features, covars = covars, datecol = datecol, white = Sys.getenv('LE8_WHITE_ONLY', unset = 'TRUE'),
			files = file.info(files)[, c('size', 'mtime')], paths = files
		)
		cache <- file.path(rawdir, 'c4.imaging.rds')
		if (file.exists(cache) && !truthy(Sys.getenv('C4_IMG_REPLACE', unset = 'FALSE'))) {
			old <- readRDS(cache) ; if (identical(old$signature, signature)) return(old)
		}
		message('C4 imaging: ', length(features), ' ', layer, ' features; complete-case models')
		all <- read_all(unique(c('eid', covars, 'center', 'ethnic.c', 'date_attend', datecol))) |> filter_analysis_cohort()
		required <- unique(c('eid', covars, 'center', 'date_attend', datecol))
		if (length(setdiff(required, names(all)))) return(missing_result(paste('Missing covariates/date:', paste(setdiff(required, names(all)), collapse = ','))))
		# Keep the outcome-free baseline cohort without requiring complete LE8 scores.
		all <- all[!is.na(all$date_attend) & (is.na(all[[datecol]]) | all[[datecol]] > all$date_attend), , drop = FALSE]
		all$eid <- as.character(all$eid) ; all$center <- factor(all$center)
		if (layer == 'protein') {
			bio <- read_prot() ; omic_input <- 'Existing prepared protein matrix (same input as C1)'
		} else {
			bio <- read_met() ; omic_input <- 'Existing prepared metabolite matrix'
		}
		bio$eid <- as.character(bio$eid) ; features <- intersect(features, names(bio))
		if (!length(features)) return(missing_result('Selected features absent from omics input'))
		img <- readRDS(imgfile) ; img$eid <- as.character(img$eid)
		stopifnot(!anyDuplicated(all$eid), !anyDuplicated(bio$eid), !anyDuplicated(img$eid))
		fields <- c4_img_manifest(indir)
		if (all(c('img_p26552', 'img_p26583') %in% names(img))) {
			img$img_total_cortical_gm <- img$img_p26552 + img$img_p26583
			fields <- bind_rows(fields, tibble(field = 'derived', label = 'Total cortical grey matter volume', variable = 'img_total_cortical_gm', family = 'Global grey matter', adjust_TIV = TRUE))
		}
		fields <- fields |> filter(variable %in% names(img))
		if (!nrow(fields)) return(missing_result('No supported brain measurements found in dictionary'))
		vfile <- file.path(indir, 'Rdata/img.visits.rds')
		if (file.exists(vfile)) {
			visits <- readRDS(vfile) ; visits$eid <- as.character(visits$eid)
			all <- merge(all, visits, by = 'eid', all.x = TRUE, sort = FALSE)
			if ('p54_i2' %in% names(all)) all$center <- factor(all$p54_i2)
			if ('p53_i2' %in% names(all)) all <- all[!is.na(all$p53_i2) & all$p53_i2 >= all$date_attend, , drop = FALSE]
		}
		covars <- unique(c(covars, 'center'))
		measures <- unique(c(fields$variable, 'img_p26521')) ; measures <- intersect(measures, names(img))
		dat <- merge(merge(all, bio[, c('eid', features), drop = FALSE], by = 'eid'), img[, c('eid', measures), drop = FALSE], by = 'eid')
		rm(all, bio, img) ; invisible(gc())
		if (!'img_p26521' %in% names(dat) && any(fields$adjust_TIV)) return(missing_result('TIV (26521) missing; volume models not fit without required adjustment'))
		rows <- parallel_map(seq_len(nrow(fields)), function(i) {
			f <- fields[i, ] ; cv <- unique(c(covars, if (f$adjust_TIV) 'img_p26521'))
			bind_rows(lapply(features, function(x) c4_img_fit(dat, x, f$variable, cv)))
		})
		associations <- bind_rows(rows) |>
			left_join(fields, by = c('measure' = 'variable')) |>
			mutate(FDR_all = p.adjust(p, 'BH'), bonferroni = pmin(1, p * n())) |>
			group_by(family) |>
			mutate(FDR_family = p.adjust(p, 'BH')) |>
			ungroup()
		status <- tibble(
			status = 'ok', participants = nrow(dat), features = length(features), measurements = nrow(fields), tested = sum(is.finite(associations$p)),
			input = omic_input, covariates = paste(covars, collapse = ','), outcome_date = datecol,
			interpretation = 'Baseline omics associated with later brain structure; no causal or within-person change estimate'
		)
		result <- list(signature = signature, status = status, associations = associations, fields = fields)
		saveRDS(result, cache) ; write_raw_csv(status, 'c4.imaging_status.csv', rawdir)
		write_raw_csv(associations, 'c4.imaging_associations.csv', rawdir) ; write_raw_csv(fields, 'c4.imaging_fields.csv', rawdir)
		plot_c4_imaging(associations, features, outdir)
		result
	}

	plot_c4_imaging <- function(associations, features, outdir) {
		le8_mock_c4(associations, outdir)
		top <- associations |>
			filter(is.finite(p)) |>
			arrange(p) |>
			slice_head(n = 20) |>
			mutate(pair = stringr::str_wrap(paste(feature, label, sep = ' → '), 48), pair = factor(pair, levels = rev(unique(pair))))
		if (nrow(top)) {
			g <- ggplot(top, aes(beta, pair, color = family)) +
				geom_vline(xintercept = 0, color = 'grey70') +
				geom_errorbarh(aes(xmin = beta - 1.96 * SE, xmax = beta + 1.96 * SE), height = .15) +
				geom_point(aes(shape = FDR_all < .05)) +
				labs(title = 'Omics and brain structure', subtitle = 'Top 20 associations; complete-case coefficients; exploratory brain context', x = 'Standardized beta (95% CI)', y = NULL, shape = 'Global FDR < 0.05', color = NULL) +
				theme_5c(9)
			save_plot(g, 'c4.Fig12.imaging_associations.png', 15, 13, outdir = outdir)
			summary <- associations |>
				group_by(feature, family) |>
				summarise(tested = sum(is.finite(p)), significant = sum(FDR_all < .05, na.rm = TRUE), .groups = 'drop')
			show <- summary |>
				group_by(feature) |>
				summarise(n = sum(significant), .groups = 'drop') |>
				arrange(desc(n), feature) |>
				slice_head(n = 30) |>
				pull(feature)
			summary <- summary |> filter(feature %in% show)
			g <- ggplot(summary, aes(family, feature, fill = significant)) +
				geom_tile(color = 'white') +
				scale_fill_viridis_c() +
				labs(title = 'Brain associations across measurement families', subtitle = 'Top 30 features by significant-pair count; FDR across all tested pairs', x = NULL, y = NULL, fill = 'Significant') +
				theme_5c(8) +
				theme(axis.text.x = element_text(angle = 25, hjust = 1))
			save_plot(g, 'c4.Fig13.imaging_overview.png', 13, max(7, length(features) * .15), outdir = outdir)
		}
	}

	attach_c4_imaging <- function(obj, layer, outdir, disease = NULL) {
		obj$imaging <- run_c4_imaging(layer, disease, outdir)
		wbfile <- le8_artifact_path('c4.out.xlsx', outdir)
		if (file.exists(wbfile)) {
			wb <- openxlsx::loadWorkbook(wbfile)
			for (nm in c('status', 'associations', 'fields')) {
				sheet <- paste0('imaging_', nm) ; if (sheet %in% names(wb)) openxlsx::removeWorksheet(wb, sheet)
				openxlsx::addWorksheet(wb, sheet) ; x <- obj$imaging[[nm]]
				if (nrow(x)) openxlsx::writeData(wb, sheet, x)
			}
			openxlsx::saveWorkbook(wb, wbfile, overwrite = TRUE)
		}
		obj
	}


	source(file.path(fdir, "c1.correlate.R"))
})
LE8_JOB <- "c4_connect"
C4_CODE_VERSION <- "2026-09-15.measured-cohort-first"
MAX_N <- as.integer(Sys.getenv("C4_MAX_N", unset = "60000"))
BLOCK <- as.integer(Sys.getenv("C4_BLOCK", unset = "80"))
FDR_CUT <- as.numeric(Sys.getenv("C4_FDR", unset = "0.05"))
SPEC_CUT <- as.numeric(Sys.getenv("C4_SPECIFICITY", unset = "0.35"))
YSP_PLUS_N <- as.integer(Sys.getenv("C4_YSP_PLUS", unset = "80"))
NS_N <- as.integer(Sys.getenv("C4_NS_TOP", unset = "160"))
MED_TOP <- as.integer(Sys.getenv("C4_MED_TOP", unset = "30"))
MED_BOOT <- as.integer(Sys.getenv("C4_MED_BOOT", unset = "100"))
MED_MAX_N <- as.integer(Sys.getenv("C4_MED_MAX_N", unset = "30000"))
C4_MODULE_MAX <- as.integer(Sys.getenv("C4_MODULE_MAX", unset = "1200"))
C4_MODULE_BOOT <- as.integer(Sys.getenv("C4_MODULE_BOOT", unset = "100"))
C4_MODULE_STABILITY <- as.numeric(Sys.getenv("C4_MODULE_STABILITY", unset = "0.70"))
C4_MODULE_K_MAX <- as.integer(Sys.getenv("C4_MODULE_K_MAX", unset = "8"))

stratified_sample <- function(d, event, max_n) {
	if (max_n <= 0 || nrow(d) <= max_n) return(d)
	if (!event %in% names(d)) return(slice_sample(d, n = max_n))
	n_all <- nrow(d)
	z <- d |>
		group_by(.data[[event]]) |>
		mutate(
			.rand = runif(n()), .rank = rank(.rand, ties.method = "first"),
			# Use doubles for the proportional allocation: max_n and
			# n() are integers, and their product can exceed 2^31 - 1.
			.take = max(1, ceiling(as.double(max_n) * n() / as.double(n_all)))
		) |>
		filter(.rank <= .take) |>
		ungroup() |>
		select( - .rand, - .rank, - .take)
	if (nrow(z) > max_n) z <- slice_sample(z, n = max_n)
	z
}

stratified_split <- function(d, event) {
	out <- rep(NA_character_, nrow(d)) ; set.seed(SEED)
	if (!event %in% names(d)) return(sample(rep(c("discovery", "replication"), length.out = nrow(d))))
	strata <- ifelse(is.na(d[[event]]), "missing", as.character(d[[event]]))
	for (v in unique(strata)) {
		id <- which(strata == v)
		out[id] <- sample(rep(c("discovery", "replication"), length.out = length(id)))
	}
	out
}

proxy_scan <- function(d, features, components, covars, split_name) {
	d <- d[complete.cases(d[, unique(c(components, covars)), drop = FALSE]), , drop = FALSE]
	features <- intersect(features, names(d)) ; components <- intersect(components, names(d)) ; covars <- intersect(covars, names(d))
	if (!length(features) || !length(components) || nrow(d) < 500) return(tibble())
	blocks <- split(features, ceiling(seq_along(features) / BLOCK)) ; out <- list() ; k <- 0L
	for (cmp in components) {
		other <- setdiff(components, cmp) ; mm <- model.matrix(reformulate(unique(c(covars, other))), d) ; q <- qr(mm)
		y <- std_num(d[[cmp]]) ; yr <- qr.resid(q, y) ; sy <- sqrt(sum(yr ^ 2)) ; df <- max(3, nrow(d) - ncol(mm) - 2)
		for (bb in blocks) {
			x <- as.matrix(d[, bb, drop = FALSE]) ; storage.mode(x) <- "double"
			for (j in seq_len(ncol(x))) {
				m <- median(x[, j], na.rm = TRUE) ; if (!is.finite(m)) m <- 0 ; x[!is.finite(x[, j]), j] <- m ; x[, j] <- std_num(x[, j]) ; x[!is.finite(x[, j]), j] <- 0
			}
			xr <- qr.resid(q, x) ; den <- sqrt(colSums(xr ^ 2)) * sy ; r <- as.numeric(crossprod(xr, yr) / pmax(den, 1e-12)) ; r <- cap(r, .999)
			z <- r * sqrt(df / pmax(1 - r ^ 2, 1e-9)) ; p <- 2 * pt(abs(z), df = df, lower.tail = FALSE) ; k <- k + 1L
			out[[k]] <- tibble(feature = bb, component = sub("\\.pts$", "", cmp), component_var = cmp, r = r, z = z, p.value = p, n = nrow(d), split = split_name)
		}
	}
	bind_rows(out) |>
		group_by(component) |>
		mutate(FDR_component = p.adjust(p.value, "BH")) |>
		ungroup() |>
		mutate(FDR_global = p.adjust(p.value, "BH"))
}

make_proxy_sets <- function(disc, rep, disease_assoc, annotation) {
	x <- disc |>
		select(feature, component, component_var, r_disc = r, z_disc = z, p_disc = p.value, FDR_disc = FDR_component) |>
		left_join(rep |> select(feature, component, r_rep = r, z_rep = z, p_rep = p.value, FDR_rep = FDR_component), by = c("feature", "component"))
	spec <- disc |>
		group_by(feature) |>
		mutate(absz = abs(z), specificity = absz / sum(absz, na.rm = TRUE)) |>
		slice_max(absz, n = 1, with_ties = FALSE) |>
		select(feature, primary_component = component, primary_component_var = component_var, specificity)
	x_primary <- x |> left_join(spec, by = "feature") ; primary <- x_primary |>
		filter(component == primary_component) |>
		mutate(
			same_direction = is.finite(r_rep) & sign(r_disc) == sign(r_rep),
			strict_YS = coalesce(FDR_disc < FDR_CUT & FDR_rep < FDR_CUT & same_direction & specificity >= SPEC_CUT, FALSE),
			YS_model = strict_YS, selection_rule = ifelse(strict_YS, "YS_strict: replicated and pillar-specific", "not supervised")
		)
	if (sum(primary$strict_YS, na.rm = TRUE) < 8) {
		fallback <- primary |>
			filter(!strict_YS) |>
			group_by(primary_component) |>
			arrange(FDR_disc, FDR_rep, desc(specificity)) |>
			slice_head(n = 5) |>
			ungroup() |>
			pull(feature)
		primary <- primary |> mutate(
			YS_model = strict_YS | feature %in% fallback,
			selection_rule = case_when(strict_YS ~ "YS_strict: replicated and pillar-specific", feature %in% fallback ~ "YS_exploratory: LE8-ranked fallback", TRUE ~ "not supervised")
		)
	}
	primary <- primary |> mutate(YS = YS_model, proxy_level = case_when(strict_YS ~ "YS_strict", YS_model ~ "YS_exploratory", TRUE ~ "not supervised"))
	ys_strict <- primary |>
		filter(strict_YS) |>
		pull(feature) |>
		unique() ; ys <- primary |>
		filter(YS_model) |>
		pull(feature) |>
		unique()
	dis <- disease_assoc
	if (!"beta" %in% names(dis)) dis$beta <- if ("estimate" %in% names(dis)) safe_log(dis$estimate) else NA_real_
	if (!"FDR" %in% names(dis)) dis$FDR <- p.adjust(dis$p.value, "BH")
	dis <- dis |> arrange(p.value)
	plus <- dis |>
		filter(!term %in% ys) |>
		slice_head(n = YSP_PLUS_N) |>
		pull(term) |>
		unique()
	ns <- dis |>
		slice_head(n = NS_N) |>
		pull(term) |>
		unique()
	membership <- tibble(feature = unique(c(ys, plus, ns))) |>
		mutate(
			in_YS = feature %in% ys, in_YSP_plus = feature %in% plus, in_NS = feature %in% ns,
			set = case_when(in_YS ~ "YS", in_YSP_plus ~ "YSP_plus", TRUE ~ "NS")
		) |>
		left_join(primary |> select(feature, primary_component, specificity, strict_YS, YS_model, proxy_level, selection_rule, r_disc, r_rep, FDR_disc, FDR_rep), by = "feature") |>
		left_join(dis |> select(feature = term, disease_beta = beta, disease_p = p.value, disease_FDR = FDR), by = "feature") |>
		left_join(annotation, by = "feature")
	list(
		primary = primary, YS = ys, YS_strict = ys_strict, YSP_plus = plus, YSP = unique(c(ys, plus)), NS = ns, membership = membership,
		YS_edges = x_primary |> filter(feature %in% ys, component == primary_component) |>
			left_join(primary |> select(feature, strict_YS, YS_model, proxy_level, selection_rule), by = "feature") |> left_join(annotation, by = "feature")
	)
}

mediation_one <- function(dat, component_var, feature, covars, tvar, evar, B = 100) {
	need <- unique(c(component_var, feature, covars, tvar, evar)) ; d <- dat[, intersect(need, names(dat)), drop = FALSE] ; d <- d[complete.cases(d), , drop = FALSE]
	if (nrow(d) < 1000 || sum(d[[evar]] == 1) < 50) return(tibble())
	d[[component_var]] <- std_num(d[[component_var]]) ; d[[feature]] <- std_num(d[[feature]])
	fa <- tryCatch(lm(reformulate(c(component_var, covars), response = feature), d), error = function(e) NULL)
	ft <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(component_var, covars)), collapse = " + "))), d), error = function(e) NULL)
	fd <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(component_var, feature, covars)), collapse = " + "))), d), error = function(e) NULL)
	if (any(vapply(list(fa, ft, fd), is.null, logical(1)))) return(tibble())
	sma <- coef(summary(fa)) ; smt <- coef(summary(ft)) ; smd <- coef(summary(fd))
	if (!component_var %in% rownames(sma) || !component_var %in% rownames(smt) ||
		!component_var %in% rownames(smd) || !feature %in% rownames(smd)) return(tibble())
	a <- sma[component_var, "Estimate"] ; ase <- sma[component_var, "Std. Error"] ; b <- smd[feature, "coef"] ; bse <- smd[feature, "se(coef)"]
	total <- smt[component_var, "coef"] ; direct <- smd[component_var, "coef"] ; ind <- a * b ; seind <- sqrt(b ^ 2 * ase ^ 2 + a ^ 2 * bse ^ 2) ; p <- ifelse(is.finite(seind) && seind > 0, 2 * pnorm(abs(ind / seind), lower.tail = FALSE), NA_real_)
	boot <- numeric() ; if (B > 0) {
		for (i in seq_len(B)) {
			id <- sample.int(nrow(d), replace = TRUE) ; db <- d[id, , drop = FALSE] ; za <- tryCatch(coef(lm(reformulate(c(component_var, covars), response = feature), db))[[component_var]], error = function(e) NA_real_) ; zb <- tryCatch(coef(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(component_var, feature, covars)), collapse = " + "))), db))[[feature]], error = function(e) NA_real_) ; boot[i] <- za * zb
		}
	}
	boot <- boot[is.finite(boot)] ; blo <- if (length(boot) >= 20) as.numeric(quantile(boot, .025, names = FALSE)) else NA_real_ ; bhi <- if (length(boot) >= 20) as.numeric(quantile(boot, .975, names = FALSE)) else NA_real_
	z <- tibble(
		component = sub("\\.pts$", "", component_var), component_var, feature, n = nrow(d), events = sum(d[[evar]] == 1), a_beta = a, b_beta = b, total_beta = total, direct_beta = direct, indirect_beta = ind, indirect_se = seind, indirect_p = p,
		indirect_lo = blo, indirect_hi = bhi, bootstrap_valid = length(boot), prop_mediated = ifelse(total != 0, ind / total, NA_real_),
		estimator = "Baseline linear-by-Cox coefficient product; associational decomposition",
		causal_identified = FALSE,
		interpretation = "prop_mediated is a descriptive coefficient ratio, not an identified natural indirect-effect fraction"
	)
	if (nrow(z)) {
		z$associational_product_ratio <- z$prop_mediated
		z$causal_mediation_identified <- FALSE
		z$modifiability_demonstrated <- FALSE
		z$interpretation <- "Same-baseline exposure and omic; coefficient-product ratio is descriptive, not an identified causal mediation proportion"
	}
	z
}


make_genetic_edges <- function(disc, rep, disease, membership) {
	# Preserve the edge schema when PRS inputs are absent or a scan is not estimable.
	# Plotting and list export must handle unavailable evidence without inventing edges.
	if (!nrow(disc) || !nrow(rep)) {
		message(
			"C4: genetic bridges unavailable; discovery rows=", nrow(disc),
			"; replication rows=", nrow(rep)
		)
		return(tibble(
			feature = character(), prs = character(), prs_var = character(),
			r_disc = double(), z_disc = double(), p_disc = double(), FDR_disc = double(),
			r_rep = double(), z_rep = double(), p_rep = double(), FDR_rep = double(),
			same_direction = logical(), replicated = logical(), disease_beta = double(),
			disease_p = double(), set = character(), primary_component = character(), bridge_score = double()
		))
	}
	dis <- disease ; if (!"beta" %in% names(dis)) dis$beta <- safe_log(dis$estimate)
	disc |>
		select(feature, prs = component, prs_var = component_var, r_disc = r, z_disc = z, p_disc = p.value, FDR_disc = FDR_component) |>
		left_join(rep |> select(feature, prs = component, r_rep = r, z_rep = z, p_rep = p.value, FDR_rep = FDR_component), by = c("feature", "prs")) |>
		mutate(same_direction = is.finite(r_rep) & sign(r_disc) == sign(r_rep), replicated = coalesce(FDR_disc < FDR_CUT & FDR_rep < FDR_CUT & same_direction, FALSE)) |>
		left_join(dis |> select(feature = term, disease_beta = beta, disease_p = p.value), by = "feature") |>
		left_join(membership |> select(feature, set, primary_component), by = "feature") |>
		mutate(bridge_score =  - log10(pmax(FDR_disc, 1e-300)) +  - log10(pmax(FDR_rep, 1e-300)) +  - log10(pmax(disease_p, 1e-300))) |>
		arrange(desc(replicated), desc(bridge_score))
}

supervised_modules <- function(disc, rep, sets) {
	p <- disc |>
		select(feature, component, z_disc = z) |>
		inner_join(rep |> select(feature, component, z_rep = z), by = c("feature", "component")) |>
		mutate(same_direction = is.finite(z_rep) & sign(z_disc) == sign(z_rep), z_joint = ifelse(same_direction, (z_disc + z_rep) / sqrt(2), .25 * (z_disc + z_rep)), strength = pmax(abs(z_disc), abs(z_rep), na.rm = TRUE))
	if (!nrow(p)) return(list(membership = tibble(), metrics = tibble(), profile = tibble(), YS_core = character()))
	keep <- p |>
		group_by(feature) |>
		summarise(strength = max(strength, na.rm = TRUE), .groups = "drop") |>
		arrange(desc(strength)) |>
		slice_head(n = C4_MODULE_MAX) |>
		pull(feature)
	pw <- p |>
		filter(feature %in% keep) |>
		select(feature, component, z_joint) |>
		pivot_wider(names_from = component, values_from = z_joint, values_fill = 0)
	rn <- pw$feature ; X <- as.matrix(pw[, - 1, drop = FALSE]) ; rownames(X) <- rn ; X[!is.finite(X)] <- 0
	X <- t(scale(t(X))) ; X[!is.finite(X)] <- 0
	n <- nrow(X) ; ks <- if (n >= 4) 2 : min(C4_MODULE_K_MAX, n - 1) else integer() ; D <- as.matrix(dist(X))
	sil_one <- function(cl) {
		mean(vapply(seq_len(n), function(i) {
			same <- which(cl == cl[i] & seq_len(n) != i) ; a <- if (length(same)) mean(D[i, same]) else 0 ; other <- setdiff(unique(cl), cl[i]) ; if (!length(other)) return(0) ; b <- min(vapply(other, function(g) mean(D[i, cl == g]), numeric(1))) ; if (max(a, b) == 0) 0 else (b - a) / max(a, b)
		}, numeric(1)), na.rm = TRUE)
	}
	hc <- hclust(as.dist(D), method = "ward.D2")
	metrics <- if (length(ks)) map_dfr(ks, function(k) {
		cl <- cutree(hc, k) ; tibble(k, silhouette = sil_one(cl), min_module_n = min(table(cl)))
	}) else tibble(k = 1L, silhouette = NA_real_, min_module_n = n)
	eligible <- metrics |> filter(min_module_n >= max(3L, ceiling(.03 * n))) ; if (!nrow(eligible)) eligible <- metrics |> slice(1)
	k <- eligible |>
		arrange(desc(silhouette), k) |>
		slice(1) |>
		pull(k) ; orig <- if (k == 1) rep(1L, n) else cutree(hc, k)
	stable <- rep(1, n)
	if (k > 1 && C4_MODULE_BOOT > 0 && ncol(X) > 1) {
		set.seed(SEED + 404) ; hit <- numeric(n)
		for (b in seq_len(C4_MODULE_BOOT)) {
			cols <- sample(seq_len(ncol(X)), ncol(X), replace = TRUE) ; cb <- cutree(hclust(dist(X[, cols, drop = FALSE]), method = "ward.D2"), k)
			tab <- table(cb, orig) ; mp <- apply(tab, 1, function(v) as.integer(colnames(tab)[which.max(v)])) ; aligned <- mp[as.character(cb)]
			hit <- hit + (aligned == orig)
		}
		stable <- hit / C4_MODULE_BOOT
	}
	membership <- tibble(feature = rn, module = orig, module_stability = stable) |>
		left_join(sets$primary |> select(feature, primary_component, strict_YS, YS_model, FDR_disc, FDR_rep), by = "feature") |>
		mutate(YS_core = coalesce(YS_model, FALSE) & module_stability >= C4_MODULE_STABILITY)
	prof <- p |>
		filter(feature %in% rn) |>
		left_join(membership |> select(feature, module, module_stability, YS_core, primary_component), by = "feature")
	list(membership = membership, metrics = metrics |> mutate(selected_k = .env$k), profile = prof, YS_core = membership |> filter(YS_core) |> pull(feature))
}

plot_supervised_atlas <- function(modules, sets, med, layer, outdir) {
	mem <- modules$membership ; pr <- modules$profile
	if (!nrow(mem) || !nrow(pr)) {
		save_plot(blank_plot("LE8-supervised biomarker atlas", "No replicated profile was clusterable"), "c4.Fig7.supervised_atlas.png", 10, 6, outdir = outdir)
		return(invisible(NULL))
	}
	top <- mem |>
		arrange(desc(YS_core), desc(module_stability), FDR_disc) |>
		slice_head(n = 40) |>
		pull(feature)
	hm <- pr |>
		filter(feature %in% top) |>
		mutate(feature = factor(feature, levels = rev(top)), component = factor(component, levels = names.le8), z = cap(z_joint, 8))
	pa <- ggplot(hm, aes(component, feature, fill = z)) +
		geom_tile(color = "white", linewidth = .15) +
		facet_grid(module ~ ., scales = "free_y", space = "free_y") +
		scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "replicated z") +
		labs(title = paste0("a. LE8-supervised ", layer, " modules"), subtitle = "Top replicated eight-pillar profiles; clustering is limited to the strongest 1,200 features", x = NULL, y = NULL) +
		theme_5c(7) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1))
	pb <- ggplot(mem, aes(module_stability, factor(module), color = YS_core)) +
		geom_jitter(height = .16, width = 0, alpha = .45, size = 1.4) +
		geom_vline(xintercept = C4_MODULE_STABILITY, linetype = 2) +
		scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D73027")) +
		labs(title = "b. Perturbation-bootstrap module stability", x = "Feature stability", y = "Module", color = "YS core") +
		theme_5c(9)
	funnel <- tibble(
		stage = factor(c("All profiled", "YS model", "YS strict", "YS core"), levels = rev(c("All profiled", "YS model", "YS strict", "YS core"))),
		n = c(nrow(mem), sum(mem$YS_model, na.rm = TRUE), sum(mem$strict_YS, na.rm = TRUE), sum(mem$YS_core, na.rm = TRUE))
	)
	pc <- ggplot(funnel, aes(n, stage, fill = stage)) +
		geom_col(width = .72) +
		geom_text(aes(label = n), hjust =  - .15, fontface = "bold") +
		scale_x_continuous(expand = expansion(mult = c(0, .18))) +
		guides(fill = "none") +
		labs(title = "c. Supervised selection funnel", x = "Features", y = NULL) +
		theme_5c(9)
	md <- med |>
		filter(is.finite(indirect_beta)) |>
		arrange(indirect_p) |>
		slice_head(n = 20) |>
		mutate(label = factor(paste(feature, component, sep = " | "), levels = rev(paste(feature, component, sep = " | "))))
	pd <- if (!nrow(md)) blank_plot("d. LE8–omic–disease paths", "No estimable path") else ggplot(md, aes(indirect_beta, label, color = component)) +
		geom_vline(xintercept = 0, color = "grey60") +
		geom_errorbar(aes(xmin = coalesce(indirect_lo, indirect_beta - 1.96 * indirect_se), xmax = coalesce(indirect_hi, indirect_beta + 1.96 * indirect_se)), orientation = "y", width = .05) +
		geom_point() +
		scale_color_manual(values = cols_le8) +
		labs(title = "d. Associational bridge estimates", x = "Indirect log-hazard effect", y = NULL, color = NULL) +
		forest_theme(8)
	save_plot((pa | pb) / (pc | pd) + plot_layout(widths = c(1.35, 1)), "c4.Fig7.supervised_atlas.png", 17, 14, outdir = outdir)
	save_plot(pc | pd, "c4.Fig9.selection_mediation.png", 14, 8, outdir = outdir)
}

plot_module_globe <- function(modules, outdir) {
	d0 <- as_tibble(modules$profile %||% tibble())
	if (!nrow(d0) || !all(c("YS_core", "feature", "module", "z_joint", "module_stability", "primary_component") %in% names(d0))) {
		save_plot(blank_plot("Stable LE8 module globe", "No YS-core feature passed stability"), "c4.Fig8.network_globe.png", 9, 6, outdir = outdir) ; return(invisible(NULL))
	}
	d <- d0 |>
		filter(YS_core) |>
		group_by(feature, module) |>
		slice_max(abs(z_joint), n = 1, with_ties = FALSE) |>
		ungroup()
	if (!nrow(d)) {
		save_plot(blank_plot("Stable LE8 module globe", "No YS-core feature passed stability"), "c4.Fig8.network_globe.png", 9, 6, outdir = outdir) ; return(invisible(NULL))
	}
	mods <- sort(unique(d$module)) ; cent <- tibble(module = mods, a = head(seq(0, 2 * pi, length.out = length(mods) + 1), - 1), x = .28 * cos(a), y = .28 * sin(a))
	d <- d |>
		group_by(module) |>
		mutate(a = head(seq(0, 2 * pi, length.out = n() + 1), - 1) + first(module) * .31, x = cos(a), y = sin(a)) |>
		ungroup() |>
		left_join(cent |> select(module, mx = x, my = y), by = "module")
	p <- ggplot(d) +
		geom_curve(aes(x = mx, y = my, xend = x, yend = y, color = primary_component, alpha = module_stability), curvature = .08, linewidth = .5) +
		geom_point(aes(x, y, size = abs(z_joint), fill = primary_component), shape = 21, color = "grey30") +
		geom_point(data = cent, aes(x, y), shape = 23, size = 5, fill = "#FFD92F") +
		geom_text(data = cent, aes(x, y, label = paste0("M", module)), fontface = "bold") +
		geom_text_repel(aes(x, y, label = feature), size = 2.3, max.overlaps = 35, seed = 14) +
		scale_color_manual(values = cols_le8) +
		scale_fill_manual(values = cols_le8) +
		coord_equal(clip = "off") +
		theme_void() +
		labs(title = "Stable LE8-supervised biomarker globe", subtitle = "Yellow diamonds are learned modules; outer nodes are stable YS-core biomarkers", color = "Primary pillar", fill = "Primary pillar", size = "Profile strength", alpha = "Stability") +
		theme(plot.title = element_text(face = "bold", size = 15), legend.position = "bottom", plot.margin = margin(20, 45, 20, 45))
	save_plot(p, "c4.Fig8.network_globe.png", 13, 11, outdir = outdir)
}

plot_c4 <- function(scan, sets, med, genetics, layer, outdir) {
	le8_plot_c4_pathways(sets, med, layer, outdir)
	topf <- sets$primary |>
		arrange(desc(strict_YS), desc(YS_model), FDR_disc) |>
		slice_head(n = 55) |>
		pull(feature)
	hm <- scan |>
		filter(split == "discovery", feature %in% topf) |>
		mutate(
			component = factor(component, levels = names.le8),
			feature = factor(feature, levels = rev(topf)), zcap = cap(z, 8)
		)
	p1 <- if (!nrow(hm)) blank_plot(paste0("LE8-supervised ", layer, " proxy map")) else ggplot(hm, aes(component, feature, fill = zcap)) +
		geom_tile(color = "white", linewidth = .18) +
		scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0, name = "Partial z") +
		labs(
			title = paste0("LE8-supervised ", layer, " proxy map"), subtitle = "Discovery associations; YS requires directionally concordant replication and pillar specificity",
			x = NULL, y = NULL
		) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1))
	save_plot(p1, "c4.Fig1.proxy_heatmap.png", 10.5, 11, outdir = outdir)

	flow <- sets$membership |>
		filter(set %in% c("YS", "YSP_plus")) |>
		count(group = coalesce(group, "Other"), primary_component = coalesce(primary_component, "Disease-selected"), set)
	if (requireNamespace("ggalluvial", quietly = TRUE) && nrow(flow)) {
		p2 <- ggplot(flow, aes(axis1 = group, axis2 = primary_component, axis3 = set, y = n)) +
			ggalluvial::geom_alluvium(aes(fill = primary_component), width = 1 / 12, alpha = .75) +
			ggalluvial::geom_stratum(width = 1 / 8, fill = "grey96", color = "grey55") +
			ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 2.5) +
			scale_x_discrete(limits = c("Omic group", "LE8 pillar", "Proxy set"), expand = c(.08, .08)) +
			scale_fill_manual(values = cols_le8, na.value = "grey65", guide = "none") +
			labs(title = "Omic group → LE8 pillar → supervised set", subtitle = "YS is LE8-supervised; YSP_plus adds disease-ranked features", y = "Features", x = NULL) +
			theme_5c(10)
	} else if (nrow(flow)) p2 <- ggplot(flow, aes(primary_component, fct_reorder(group, n, sum), size = n, color = set)) +
		geom_point(alpha = .8) +
		scale_size_continuous(range = c(2, 10)) +
		scale_color_manual(values = c(YS = "#D95F02", YSP_plus = "#1B9E77")) +
		labs(title = "Omic groups assigned to LE8 pillars", x = NULL, y = NULL, color = NULL, size = "Features") +
		theme_5c(10)
	else p2 <- blank_plot("Omic groups assigned to LE8 pillars")
	save_plot(p2, "c4.Fig2.group_pillar_flow.png", 13, 7.5, outdir = outdir)

	le <- sets$YS_edges |>
		arrange(FDR_disc) |>
		slice_head(n = 18) |>
		mutate(side = "LE8")
	ge <- genetics |>
		filter(replicated) |>
		arrange(desc(bridge_score)) |>
		slice_head(n = 18) |>
		mutate(side = "Genetics")
	if (nrow(le) || nrow(ge)) {
		lfeat <- unique(le$feature) ; rfeat <- unique(ge$feature) ; lpil <- unique(le$component) ; rprs <- unique(ge$prs)
		make_nodes <- function(nms, type, x, lo = .05, hi = .95) {
			nms <- unique(na.omit(as.character(nms)))
			if (!length(nms)) return(tibble(name = character(), type = character(), x = numeric(), y = numeric()))
			tibble(name = nms, type = type, x = x, y = seq(lo, hi, length.out = length(nms)))
		}
		node <- bind_rows(
			make_nodes(lpil, "LE8", - 2, .08, .92), make_nodes(lfeat, "omic_left", - 1),
			tibble(name = Y, type = "disease", x = 0, y = .5), make_nodes(rfeat, "omic_right", 1), make_nodes(rprs, "PRS", 2, .08, .92)
		)
		el <- le |>
			left_join(node |> filter(type == "LE8") |> select(component = name, x1 = x, y1 = y), by = "component") |>
			left_join(node |> filter(type == "omic_left") |> select(feature = name, x2 = x, y2 = y), by = "feature")
		er <- ge |>
			left_join(node |> filter(type == "PRS") |> select(prs = name, x1 = x, y1 = y), by = "prs") |>
			left_join(node |> filter(type == "omic_right") |> select(feature = name, x2 = x, y2 = y), by = "feature")
		toY <- bind_rows(node |> filter(type == "omic_left") |> transmute(x1 = x, y1 = y, x2 = 0, y2 = .5), node |> filter(type == "omic_right") |> transmute(x1 = x, y1 = y, x2 = 0, y2 = .5))
		p3 <- ggplot() +
			geom_curve(
				data = el, aes(x = x1, y = y1, xend = x2, yend = y2, color = component), curvature = .08, linewidth = .65,
				arrow = grid::arrow(length = grid::unit(1.6, "mm"))
			) +
			geom_curve(
				data = er, aes(x = x1, y = y1, xend = x2, yend = y2), curvature =  - .08, color = "#6A3D9A", linewidth = .65,
				arrow = grid::arrow(length = grid::unit(1.6, "mm"))
			) +
			geom_curve(
				data = toY, aes(x = x1, y = y1, xend = x2, yend = y2), curvature = .06, color = "grey68", alpha = .55,
				arrow = grid::arrow(length = grid::unit(1.5, "mm"))
			) +
			geom_point(data = node, aes(x, y, shape = type), size = 2.5) +
			geom_text(data = node, aes(x, y, label = name), size = 2.4, nudge_y = .022, check_overlap = TRUE) +
			scale_color_manual(values = cols_le8, guide = "none") +
			coord_cartesian(xlim = c( - 2.25, 2.25), ylim = c(0, 1), clip = "off") +
			theme_void() +
			labs(title = "LE8 → omics → disease ← omics ← genetics", subtitle = "Top replicated bridges; PRS associations are genetic anchors, not proof of mediation") +
			theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(color = "grey35"), plot.margin = margin(12, 25, 12, 25))
	} else p3 <- blank_plot("LE8 → omics → disease ← omics ← genetics", "No replicated bridge was available")
	save_plot(p3, "c4.Fig3.connection_bridge.png", 15, 8.5, outdir = outdir)

	evf <- unique(c(head(le$feature, 25), head(ge$feature, 25)))
	mat1 <- sets$primary |>
		filter(feature %in% evf) |>
		transmute(feature, metric = "LE8 proxy", value = z_disc, set = ifelse(strict_YS, "YS strict", ifelse(YS_model, "YS exploratory", "Other")))
	mat2 <- sets$membership |>
		filter(feature %in% evf) |>
		transmute(feature, metric = "Disease", value = sign(disease_beta) *  - log10(pmax(disease_p, 1e-300)), set = coalesce(set, "Other"))
	mat3 <- genetics |>
		filter(feature %in% evf, replicated) |>
		group_by(feature) |>
		slice_max(abs(z_disc), n = 1, with_ties = FALSE) |>
		ungroup() |>
		transmute(feature, metric = "Genetic anchor", value = z_disc, set = coalesce(set, "Other"))
	em <- bind_rows(mat1, mat2, mat3) |> mutate(value = cap(value, 8), feature = factor(feature, levels = rev(evf)), metric = factor(metric, levels = c("LE8 proxy", "Disease", "Genetic anchor")))
	p4 <- if (!nrow(em)) blank_plot("Connection evidence matrix") else ggplot(em, aes(metric, feature)) +
		geom_point(aes(size = abs(value), fill = value, shape = set), color = "grey40") +
		scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0) +
		scale_size_continuous(range = c(1.5, 7)) +
		labs(
			title = "Connection evidence matrix", subtitle = "Signed evidence across modifiable, proximal disease, and genetic dimensions",
			x = NULL, y = NULL, fill = "Signed evidence", size = "|Evidence|", shape = "Set"
		) +
		theme_5c(9) +
		theme(legend.position = "bottom")
	save_plot(p4, "c4.Fig4.connection_evidence.png", 11.5, 9.5, outdir = outdir)

	md <- med |>
		filter(is.finite(indirect_beta)) |>
		arrange(indirect_p) |>
		slice_head(n = 35) |>
		mutate(label = paste(feature, component, sep = " | "), label = factor(label, levels = rev(label)))
	p5 <- if (!nrow(md)) blank_plot("Associational LE8–omic–disease paths") else ggplot(md, aes(indirect_beta, label, color = component)) +
		geom_vline(xintercept = 0, color = "grey60") +
		geom_errorbar(aes(xmin = coalesce(indirect_lo, indirect_beta - 1.96 * indirect_se), xmax = coalesce(indirect_hi, indirect_beta + 1.96 * indirect_se)), orientation = "y", width = .08) +
		geom_point(aes(size = pmin(abs(prop_mediated), 1)), alpha = .85) +
		scale_color_manual(values = cols_le8) +
		labs(
			title = "Associational LE8–omic–disease paths", subtitle = "Same-baseline LE8 and omics: indirect effects do not establish temporal mediation",
			x = "Indirect log-hazard effect", y = NULL, color = "LE8 pillar", size = "|Proportion mediated|"
		) +
		forest_theme(8) +
		theme(legend.position = "bottom")
	save_plot(p5, "c4.Fig5.mediation_forest.png", 12, 10, outdir = outdir)

	if (nrow(md)) {
		pb <- ggplot(md, aes(total_beta, direct_beta, color = component)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
			geom_point(aes(size = abs(prop_mediated)), alpha = .8) +
			scale_color_manual(values = cols_le8) +
			labs(title = "a. Total versus direct effects", x = "Total effect", y = "Direct effect", color = NULL, size = "|Proportion mediated|") +
			theme_5c(9)
		pc <- ggplot(md, aes(100 * cap(prop_mediated, 1), - log10(pmax(indirect_p, 1e-30)), color = component)) +
			geom_vline(xintercept = 0, color = "grey60") +
			geom_point(alpha = .8) +
			scale_color_manual(values = cols_le8) +
			labs(title = "b. Proportion and indirect-effect evidence", x = "Percent mediated (capped)", y = expression( - log[10](P[indirect])), color = NULL) +
			theme_5c(9)
		p6 <- pb | pc
	} else p6 <- blank_plot("Mediation diagnostics")
	save_plot(p6, "c4.Fig6.mediation_diagnostics.png", 14, 7.5, outdir = outdir)
}

# Retained only to reproduce the pre-genetics C4 figure set.  Keeping the same
# name used to override the six-argument publication function above and made
# every current C4 run fail with an unused `genetics` argument.
plot_c4_legacy <- function(scan, sets, med, layer, outdir) {
	topf <- sets$primary |>
		arrange(desc(YS), FDR_disc) |>
		slice_head(n = 80) |>
		pull(feature)
	hm <- scan |>
		filter(split == "discovery", feature %in% topf) |>
		mutate(component = factor(component, levels = names.le8), feature = factor(feature, levels = rev(topf)), zcap = cap(z, 8))
	p1 <- if (!nrow(hm)) blank_plot(paste0("LE8-supervised ", layer, " proxy map"), "No discovery association was available") else ggplot(hm, aes(component, feature, fill = zcap)) +
		geom_tile(color = "white", linewidth = .15) +
		scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0, name = "z") +
		labs(title = paste0("LE8-supervised ", layer, " proxy map"), subtitle = "Discovery associations; proxy status requires replication and pillar specificity", x = NULL, y = NULL) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1))
	save_plot(p1, "c4.Fig1.proxy_heatmap.png", 10, 11, outdir = outdir)

	flow <- sets$membership |>
		filter(set %in% c("YS", "YSP_plus")) |>
		count(group = coalesce(group, "Other"), primary_component = coalesce(primary_component, "Disease-selected"), set)
	if (requireNamespace("ggalluvial", quietly = TRUE) && nrow(flow)) {
		p2 <- ggplot(flow, aes(axis1 = group, axis2 = primary_component, axis3 = set, y = n)) +
			ggalluvial::geom_alluvium(aes(fill = primary_component), width = 1 / 12, alpha = .75) +
			ggalluvial::geom_stratum(width = 1 / 8, fill = "grey94", color = "grey55") +
			ggalluvial::stat_stratum(geom = "text", aes(label = after_stat(stratum)), size = 2.6) +
			scale_x_discrete(limits = c("Omic group", "LE8 pillar", "Proxy set"), expand = c(.08, .08)) +
			scale_fill_manual(values = cols_le8, na.value = "grey65", guide = "none") +
			labs(title = "Omic group → LE8 pillar → supervised set", y = "Number of features", x = NULL) +
			theme_5c(10)
	} else if (nrow(flow)) {
		p2 <- ggplot(flow, aes(primary_component, fct_reorder(group, n, sum), size = n, color = set)) +
			geom_point(alpha = .8) +
			scale_size_continuous(range = c(2, 10)) +
			scale_color_manual(values = c(YS = "#D95F02", YSP_plus = "#1B9E77")) +
			labs(title = "Omic groups assigned to LE8 pillars", x = NULL, y = NULL, color = NULL, size = "Features") +
			theme_5c(10) +
			theme(axis.text.x = element_text(angle = 35, hjust = 1))
	} else p2 <- blank_plot("Omic groups assigned to LE8 pillars", "No proxy-set flow could be formed")
	save_plot(p2, "c4.Fig2.group_pillar_flow.png", 12, 7, outdir = outdir)

	edges <- sets$YS_edges |>
		arrange(FDR_disc) |>
		slice_head(n = 100)
	if (nrow(edges)) {
		pillars <- tibble(name = unique(edges$component), type = "pillar", angle = seq(0, 2 * pi, length.out = n_distinct(edges$component) + 1)[ - (n_distinct(edges$component) + 1)], r = .35)
		feats <- tibble(name = unique(edges$feature), type = "feature", angle = seq(0, 2 * pi, length.out = n_distinct(edges$feature) + 1)[ - (n_distinct(edges$feature) + 1)] + .12, r = 1)
		nodes <- bind_rows(pillars, feats) |> mutate(x = r * cos(angle), y = r * sin(angle)) ; ee <- edges |>
			left_join(nodes |> select(component = name, x1 = x, y1 = y), by = "component") |>
			left_join(nodes |> select(feature = name, x2 = x, y2 = y), by = "feature")
		p3 <- ggplot() +
			geom_curve(data = ee, aes(x = x1, y = y1, xend = x2, yend = y2, color = component, alpha = pmin(abs(z_disc) / 8, 1)), curvature = .12, linewidth = .55) +
			geom_point(data = nodes, aes(x, y, shape = type), size = 2.5) +
			ggrepel::geom_text_repel(data = nodes, aes(x, y, label = name), size = 2.3, max.overlaps = 35, seed = 11) +
			scale_color_manual(values = cols_le8) +
			scale_alpha(range = c(.2, .9), guide = "none") +
			coord_equal() +
			theme_void() +
			labs(title = "LE8 pillar–proxy module network", color = "Pillar") +
			theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
	} else p3 <- blank_plot("LE8 pillar–proxy module network")
	save_plot(p3, "c4.Fig3.module_network.png", 11, 9, outdir = outdir)

	# Globe/orbit view: the radius encodes pillar specificity and the angle encodes the primary LE8 component.
	orbit <- sets$membership |>
		filter(set %in% c("YS", "YSP_plus")) |>
		mutate(comp = factor(coalesce(primary_component, "disease"), levels = c(names.le8, "disease")), angle = 2 * pi * (as.numeric(comp) - 1) / nlevels(comp) + runif(n(), - .18, .18), radius = ifelse(set == "YS", .65 + .35 * coalesce(specificity, 0), 1.18), x = radius * cos(angle), y = radius * sin(angle), lab = ifelse(min_rank(disease_p) <= 25 | set == "YS" & min_rank(FDR_disc) <= 20, feature, NA_character_))
	p4 <- if (!nrow(orbit)) blank_plot("LE8–omics relationship globe", "No supervised or plus feature was available") else ggplot(orbit, aes(x, y)) +
		annotate("path", x = cos(seq(0, 2 * pi, length.out = 300)), y = sin(seq(0, 2 * pi, length.out = 300)), color = "grey80") +
		geom_segment(aes(x = 0, y = 0, xend = x, yend = y, color = primary_component), alpha = .12) +
		geom_point(aes(color = primary_component, shape = set, size =  - log10(pmax(disease_p, 1e-30))), alpha = .8) +
		ggrepel::geom_text_repel(aes(label = lab), size = 2.4, max.overlaps = 30, seed = 12, na.rm = TRUE) +
		scale_color_manual(values = c(cols_le8, disease = "grey35"), na.value = "grey60") +
		scale_size_continuous(range = c(1.5, 5), name = "Disease evidence") +
		coord_equal() +
		theme_void() +
		labs(title = "LE8–omics relationship globe", color = "Primary pillar", shape = "Set") +
		theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
	save_plot(p4, "c4.Fig4.relationship_globe.png", 11, 9, outdir = outdir)

	if (nrow(med)) {
		wheel <- med |>
			filter(is.finite(prop_mediated)) |>
			arrange(indirect_p) |>
			slice_head(n = 60) |>
			mutate(wheel_label = paste(feature, component, sep = " | "), wheel_label = factor(wheel_label, levels = wheel_label), id = row_number(), angle = 90 - 360 * (id - .5) / n(), hjust = ifelse(angle <  - 90, 1, 0), angle = ifelse(angle <  - 90, angle + 180, angle), value = 100 * cap(prop_mediated, 1))
		p5 <- if (!nrow(wheel)) blank_plot("Associational mediation wheel", "Mediation effects were estimable, but no finite proportion mediated was available") else ggplot(wheel, aes(id, abs(value), fill = component)) +
			geom_col(width = .8) +
			geom_hline(yintercept = 0) +
			geom_text(aes(y = abs(value) + 3, label = wheel_label, angle = angle, hjust = hjust), size = 2, na.rm = TRUE) +
			coord_polar(clip = "off") +
			scale_fill_manual(values = cols_le8) +
			ylim( - 10, max(abs(wheel$value), na.rm = TRUE) + 15) +
			labs(title = "Associational mediation wheel", subtitle = "Percent mediated is shown for screened LE8 → omic → disease paths", fill = "LE8 pillar") +
			theme_void(10) +
			theme(plot.title = element_text(face = "bold", hjust = .5), plot.subtitle = element_text(hjust = .5), legend.position = "bottom", plot.margin = margin(30, 55, 30, 55))
		save_plot(p5, "c4.Fig5.mediation_wheel.png", 11, 10, outdir = outdir)
		md <- med |>
			filter(is.finite(indirect_beta), is.finite(total_beta), is.finite(direct_beta)) |>
			arrange(indirect_p) |>
			slice_head(n = 50) |>
			mutate(label = paste(feature, component, sep = " | "), label = factor(label, levels = rev(label)))
		if (!nrow(md)) {
			save_plot(blank_plot("Mediation diagnostics", "No finite mediation diagnostic row was available"), "c4.Fig6.mediation_diagnostics.png", 9, 5, outdir = outdir)
		} else {
			pa <- ggplot(md, aes(indirect_beta, label)) +
				geom_vline(xintercept = 0, color = "grey60") +
				geom_errorbar(aes(xmin = coalesce(indirect_lo, indirect_beta - 1.96 * indirect_se), xmax = coalesce(indirect_hi, indirect_beta + 1.96 * indirect_se)), orientation = "y", width = 0) +
				geom_point(aes(color = component), size = 2) +
				scale_color_manual(values = cols_le8) +
				labs(title = "A. Indirect effects", x = "Indirect log-hazard effect", y = NULL, color = NULL) +
				forest_theme(9)
			pb <- ggplot(md, aes(total_beta, direct_beta, color = component)) +
				geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
				geom_point(aes(size = abs(prop_mediated)), alpha = .8) +
				scale_color_manual(values = cols_le8) +
				labs(title = "B. Total versus direct effects", x = "Total effect", y = "Direct effect", color = NULL, size = "|Proportion mediated|") +
				theme_5c(10)
			pc <- ggplot(md, aes(100 * cap(prop_mediated, 1), - log10(pmax(indirect_p, 1e-30)), color = component)) +
				geom_vline(xintercept = 0, color = "grey60") +
				geom_point(alpha = .8) +
				scale_color_manual(values = cols_le8) +
				labs(title = "C. Proportion and evidence", x = "Percent mediated (plot capped)", y = expression( - log[10](P[indirect])), color = NULL) +
				theme_5c(10)
			save_plot(pa | (pb / pc), "c4.Fig6.mediation_diagnostics.png", 14, 10, outdir = outdir)
		}
	} else {
		save_plot(blank_plot("Associational mediation wheel", "No screened YS path could be estimated"), "c4.Fig5.mediation_wheel.png", 9, 5, outdir = outdir)
		save_plot(blank_plot("Mediation diagnostics", "No screened YS path could be estimated"), "c4.Fig6.mediation_diagnostics.png", 9, 5, outdir = outdir)
	}
}

# Connection-domain summaries
le8_c4_additions <- function(out, outdir) {
	scan <- as_tibble(out$scan %||% tibble())
	if (!nrow(scan) || !"component" %in% names(scan))
		return(list(status = tibble(status = "connection scan unavailable")))
	scan <- scan |>
		mutate(
			connection_domain = ifelse(component %in% c("diet", "pa", "smoke", "sleep"), "Behavioral LE4", "Biological LE4"),
			causal_interpretation = "not established"
		)
	count <- scan |>
		group_by(connection_domain, component, split) |>
		summarise(tested = n(), FDR_supported = sum(is.finite(FDR_component) & FDR_component < 0.05), .groups = "drop")
	rd <- le8_job_dir(outdir, "c4_connect")
	write_raw_csv(scan, "c4.behavior_biology_associations.csv", rd)
	write_raw_csv(count, "c4.behavior_biology_summary.csv", rd)
	list(behavior_biology = count, association_scope = scan, interpretation = tibble(statement = c(
		"LE8 association does not prove that changing behavior will change this protein",
		"A non-PGS residual still includes uncaptured inherited variation", "Same-baseline path products are not an identified interventional mediation proportion"
	)))
}

run_c4_layer <- function(layer = c("protein", "metabolite")) {
	if (LE8_REUSE_RESULTS) return(le8_restore_outputs(match.arg(layer), "c4_connect"))
	# Check reusable results and initialize the analysis output directory.
	layer <- match.arg(layer)
	le8_begin_analysis(layer, "c4_connect")
	.le8_analysis_env <- environment()
	on.exit(le8_finish_analysis(layer, "c4_connect", .le8_analysis_env), add = TRUE)

	layer <- match.arg(layer) ; outdir <- if (layer == "protein") out.prot else out.met ; setwd2(outdir) ; rawdir <- le8_job_dir(outdir, LE8_JOB) ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE) ; cache <- file.path(rawdir, "c4.res.rds")
	if (cache_valid(cache)) {
		old <- tryCatch(readRDS(cache), error = function(e) NULL)
		if (is.list(old) && all(c("meta", "scan", "membership", "modules", "mediation") %in% names(old))) {
			cache_message(paste0("C4/", layer), cache)
			old <- attach_c4_imaging(old, layer, outdir)
			saveRDS(old, cache, compress = "xz") ; finalize_outputs(LE8_JOB, outdir)
			return(le8_restore_outputs(layer, "c4_connect"))
		}
		message("C4/", layer, ": cached result incomplete; recomputing")
	}
	c1f <- file.path(le8_job_dir(outdir, "c1_correlate"), "c1.res.rds") ; if (!file.exists(c1f)) stop("Run C1 first.", call. = FALSE) ; c1 <- readRDS(c1f) ; le8_check_options(c1) ; disease <- (c1$association %||% c1$pwas_incident %||% c1$MWAS) |> as_tibble() ; if (!"beta" %in% names(disease)) disease <- disease |> mutate(beta = safe_log(estimate))
	biom <- if (layer == "protein") read_prot() else read_met() ; features <- setdiff(names(biom), "eid") ; ann <- layer_annotation(layer, features)
	all0 <- read_all() ; prs_vars <- find_prs_vars(all0, Y, 4)
	need <- unique(c("eid", "ethnic.c", le8_custom_covars, vars.basic, vars.le8, prs_vars, "birth_date", "date_attend", "date_lost", "date_death", paste0("fod_icd10_", Y)))
	dat0 <- all0[, intersect(need, names(all0)), drop = FALSE] |>
		filter_analysis_cohort() |>
		make_outcome(Y) ; rm(all0) ; invisible(gc())
	tvar <- paste0(Y, ".t2e") ; evar <- paste0(Y, ".Yt2e") ; bvar <- paste0(Y, ".b2e")
	comps <- intersect(vars.le8, names(dat0)) ; covs <- intersect(if (le8_custom_adjustment()) le8_custom_covars else vars.basic, names(dat0))
	dat0$eid <- as.character(dat0$eid) ; biom$eid <- as.character(biom$eid)
	if (anyDuplicated(dat0$eid) || anyDuplicated(biom$eid)) stop("Duplicate phenotype/omic eid")
	# Restrict to measured donors BEFORE the cap. Sampling all UKB first discards
	# most proteomic donors because only a subset has the assay panel.
	cohort_audit <- tibble(stage = "eligible phenotype", N = nrow(dat0))
	dat0 <- dat0[dat0$eid %in% biom$eid & complete.cases(dat0[, comps, drop = FALSE]), , drop = FALSE]
	cohort_audit <- bind_rows(cohort_audit, tibble(stage = "measured omics and complete LE8", N = nrow(dat0)))
	set.seed(SEED) ; dat0 <- stratified_sample(dat0, evar, MAX_N)
	cohort_audit <- bind_rows(cohort_audit, tibble(stage = "after cap within measured cohort", N = nrow(dat0)))
	write_raw_csv(cohort_audit, "c4.cohort_sampling_audit.csv", rawdir)
	# Join the wide matrix only after eligibility/capping to keep peak memory low.
	dat <- inner_join(dat0, biom |> filter(eid %in% dat0$eid), by = "eid")
	rm(biom, dat0) ; invisible(gc())
	fold <- stratified_split(dat, evar)
	scan_cache <- file.path(rawdir, "c4.LE8_proxy_scan.rds") ; scan_stage <- read_stage_cache(scan_cache)
	if (is.null(scan_stage)) {
		disc <- proxy_scan(dat[fold == "discovery", ], features, comps, covs, "discovery")
		rep <- proxy_scan(dat[fold == "replication", ], features, comps, covs, "replication")
		scan_stage <- list(discovery = disc, replication = rep) ; write_stage_cache(scan_stage, scan_cache)
	} else {
		disc <- scan_stage$discovery ; rep <- scan_stage$replication ; message("C4/", layer, ": reuse LE8 proxy-scan stage")
	}
	if (!nrow(disc) || !nrow(rep)) {
		reason <- paste0(
			"No discovery/replication proxy scan was estimable (discovery rows=", nrow(disc),
			"; replication rows=", nrow(rep), "). Check LE8 variables and complete-case sample size."
		)
		message("C4/", layer, ": ", reason, " Writing blank panels and continuing.")
		suffix <- if (layer == "protein") c(
			"proxy_heatmap", "group_pillar_flow", "connection_bridge",
			"connection_evidence", "mediation_forest", "mediation_diagnostics", "supervised_atlas",
			"network_globe", "selection_mediation", "state_network_remodeling", "state_network_edges"
		) else c(
			"proxy_heatmap", "group_pillar_flow",
			"module_network", "relationship_globe", "mediation_wheel", "mediation_diagnostics",
			"supervised_atlas", "network_globe", "selection_mediation", "state_network_remodeling", "state_network_edges"
		)
		for (i in seq_along(suffix)) save_plot(blank_plot(paste0("C4 Figure ", i), reason),
			paste0("c4.Fig", i, ".", suffix[[i]], ".png"), 10, 6,
			outdir = outdir
		)
		lists <- list(
			YS_by_component = list(), YS_all = character(), YS_strict = character(), YS_core = character(),
			YSP_plus = character(), YSP_all = character(), NS = character(), index = tibble(),
			genetic_anchors = character(), PRS_variables = prs_vars
		)
		status <- tibble(status = "unavailable", detail = reason, discovery_rows = nrow(disc), replication_rows = nrow(rep))
		out <- list(
			meta = module_meta(layer, extra = list(code_version = C4_CODE_VERSION, status = "unavailable")),
			status = status, scan = bind_rows(disc, rep), primary = tibble(), membership = tibble(), YS_edges = tibble(),
			modules = list(membership = tibble(), metrics = tibble(), YS_core = character()), genetic_scan = tibble(),
			genetic_edges = tibble(), PRS_variables = prs_vars, mediation = tibble(), state_network = list(), lists = lists
		)
		out <- attach_c4_imaging(out, layer, outdir, disease)
		write_raw_csv(status, "c4.status.csv", rawdir) ; saveRDS(out, cache, compress = "xz")
		write_xlsx2(list(status = status, all_LE8_associations = out$scan), "c4.out.xlsx")
		finalize_outputs(LE8_JOB, outdir) ; return(out)
	}
	scan <- bind_rows(disc, rep) ; sets <- make_proxy_sets(disc, rep, disease, ann)
	modules <- supervised_modules(disc, rep, sets)
	write_raw_csv(scan, "c4.LE8_feature_associations.csv", rawdir) ; write_raw_csv(sets$primary, "c4.primary_pillar_assignment.csv", rawdir) ; write_raw_csv(sets$membership, "c4.proxy_membership_YS_YSP_NS.csv", rawdir) ; write_raw_csv(sets$YS_edges, "c4.YS_edges.csv", rawdir)
	write_raw_csv(modules$membership, "c4.supervised_module_membership.csv", rawdir) ; write_raw_csv(modules$metrics, "c4.supervised_module_selection.csv", rawdir)
	genetic_cache <- file.path(rawdir, "c4.PRS_proxy_scan.rds") ; genetic_stage <- read_stage_cache(genetic_cache)
	if (is.null(genetic_stage)) {
		gdisc <- if (length(prs_vars)) proxy_scan(dat[fold == "discovery", ], features, prs_vars, covs, "discovery") else tibble()
		grep0 <- if (length(prs_vars)) proxy_scan(dat[fold == "replication", ], features, prs_vars, covs, "replication") else tibble()
		genetic_stage <- list(discovery = gdisc, replication = grep0) ; write_stage_cache(genetic_stage, genetic_cache)
	} else {
		gdisc <- genetic_stage$discovery ; grep0 <- genetic_stage$replication ; message("C4/", layer, ": reuse PRS proxy-scan stage")
	}
	gscan <- bind_rows(gdisc, grep0) ; gedges <- make_genetic_edges(gdisc, grep0, disease, sets$membership)
	write_raw_csv(gscan, "c4.PRS_feature_associations.csv", rawdir) ; write_raw_csv(gedges, "c4.genetic_omic_disease_bridges.csv", rawdir)
	matched_pgs <- le8_pgs_bridge(dat, features, fold, covs, disease, sets$membership, layer, rawdir)
	pairs <- sets$YS_edges |>
		left_join(disease |> select(feature = term, disease_p = p.value), by = "feature") |>
		arrange(desc(strict_YS), FDR_disc, disease_p) |>
		slice_head(n = MED_TOP)
	med_need <- intersect(unique(c(tvar, evar, covs, comps, pairs$feature)), names(dat))
	med_dat <- stratified_sample(dat[, med_need, drop = FALSE], evar, MED_MAX_N)
	mediation_cache <- file.path(rawdir, "c4.mediation_stage.rds") ; med <- read_stage_cache(mediation_cache)
	if (is.null(med)) {
		med <- if (!nrow(pairs)) tibble() else map_dfr(seq_len(nrow(pairs)), function(i) {
			rw <- pairs[i, , drop = FALSE]
			mediation_one(med_dat, rw$component_var[[1]], rw$feature[[1]], unique(c(covs, setdiff(comps, rw$component_var[[1]]))), tvar, evar, MED_BOOT)
		})
		if (nrow(med)) med <- med |>
			mutate(FDR_indirect = p.adjust(indirect_p, "BH")) |>
			arrange(indirect_p)
		write_stage_cache(med, mediation_cache)
	} else message("C4/", layer, ": reuse mediation stage")
	write_raw_csv(med, "c4.mediation_all.csv", rawdir) ; plot_c4(scan, sets, med, gedges, layer, outdir)
	plot_supervised_atlas(modules, sets, med, layer, outdir) ; plot_module_globe(modules, outdir)
	state_network <- run_c4_state_network(
		dat, disease, sets, modules, features, covs, tvar, evar, bvar,
		layer, rawdir, outdir
	)
	write_raw_csv(state_network$status %||% tibble(), "c4.state_network_status.csv", rawdir)
	write_raw_csv(state_network$state_counts %||% tibble(), "c4.state_network_counts.csv", rawdir)
	write_raw_csv(state_network$edges %||% tibble(), "c4.state_network_edges.csv", rawdir)
	write_raw_csv(state_network$hubs %||% tibble(), "c4.state_network_hubs.csv", rawdir)
	lists <- list(
		YS_by_component = split(sets$YS_edges$feature, sets$YS_edges$component), YS_all = sets$YS, YS_strict = sets$YS_strict, YS_core = modules$YS_core, YSP_plus = sets$YSP_plus, YSP_all = sets$YSP, NS = sets$NS, index = sets$membership,
		genetic_anchors = gedges |> filter(replicated) |> pull(feature) |> unique(), PRS_variables = prs_vars
	)
	saveRDS(lists, file.path(rawdir, paste0("c4.", layer, "_lists.rds")), compress = "xz") ; if (layer == "protein") saveRDS(lists, file.path(rawdir, "c4.protein_lists.rds"), compress = "xz")
	out <- list(
		meta = module_meta(layer, extra = list(code_version = C4_CODE_VERSION, status = "ok", scan_N = nrow(dat), discovery_N = sum(fold == "discovery"), replication_N = sum(fold == "replication"), mediation_N = nrow(med_dat), YS_strict_n = length(sets$YS_strict), YS_model_n = length(sets$YS), PRS_n = length(prs_vars))),
		scan = scan, primary = sets$primary, membership = sets$membership, YS_edges = sets$YS_edges, modules = modules,
		genetic_scan = gscan, genetic_edges = gedges, matched_PGS = matched_pgs, PRS_variables = prs_vars, mediation = med,
		state_network = state_network, lists = lists
	)
	saveRDS(out, cache, compress = "xz") ; write_xlsx2(
		list(
			primary_assignment = sets$primary, proxy_membership = sets$membership, YS_edges = sets$YS_edges,
			supervised_modules = modules$membership, module_selection = modules$metrics, all_LE8_associations = scan, PRS_associations = gscan, genetic_omic_bridges = gedges, mediation = med,
			state_network_status = state_network$status %||% tibble(), state_counts = state_network$state_counts %||% tibble(),
			state_network_edges = state_network$edges %||% tibble(), state_network_hubs = state_network$hubs %||% tibble()
		),
		"c4.out.xlsx"
	)
	out <- attach_c4_imaging(out, layer, outdir, disease)
	saveRDS(out, cache, compress = "xz") ; finalize_outputs(LE8_JOB, outdir) ; out
}


# Sex-specific omics associations and pairwise LE8 component interactions.
C4_SEX_MAX <- as.integer(Sys.getenv("C4_SEX_MAX", unset = "120"))

run_c4_sex_interactions <- function(layer = c("protein", "metabolite")) {
	layer <- match.arg(layer)
	outdir <- if (layer == "protein")
		out.prot else out.met
	setwd2(outdir)
	rawdir <- le8_job_dir(outdir, LE8_JOB)
	dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	sex_cache <- file.path(rawdir, "c4.sex_interaction.csv")
	if (LE8_REUSE_RESULTS || cache_valid(sex_cache)) {
		res <- le8_figure_csv(rawdir, "c4.sex_interaction.csv")
		plot_c4_sex_interactions(res, layer, outdir)
		finalize_outputs(LE8_JOB, outdir)
		return(invisible(res))
	}
	c1f <- file.path(le8_job_dir(outdir, "c1_correlate"), "c1.res.rds")
	if (!file.exists(c1f))
		stop("Run C1 first.", call. = FALSE)
	c1 <- readRDS(c1f)
	a <- (c1$association %||% c1$pwas_incident %||% c1$MWAS) |>
		as_tibble()
	features <- a |>
		arrange(p.value) |>
		slice_head(n = C4_SEX_MAX) |>
		pull(term)
	biom <- if (layer == "protein")
		read_prot() else read_met()
	features <- intersect(features, setdiff(names(biom), "eid"))
	biom <- biom[, c("eid", features), drop = FALSE]
	need <- unique(c(
		"eid", "ethnic.c", vars.basic, "sex", "birth_date", "date_attend", "date_lost", "date_death",
		paste0("fod_icd10_", Y)
	))
	d <- read_all(need) |>
		left_join(biom, by = "eid") |>
		filter_analysis_cohort() |>
		make_outcome(Y)
	tvar <- paste0(Y, ".t2e")
	evar <- paste0(Y, ".Yt2e")
	covs <- setdiff(intersect(vars.basic, names(d)), "sex")
	d$sex_f <- factor(d$sex)
	features <- intersect(features, names(d))
	sex_levels <- levels(d$sex_f)
	sex_name <- function(x) if (identical(sex_levels, c("0", "1")))
		c(`0` = "Female", `1` = "Male")[[as.character(x)]] else paste("Sex category", x)
	res <- map_dfr(features, function(f) {
		z <- d[, unique(c(tvar, evar, f, "sex_f", covs)), drop = FALSE]
		z <- z[complete.cases(z), , drop = FALSE]
		z$sex_f <- droplevels(z$sex_f)
		if (nrow(z) < 500 || sum(z[[evar]] == 1) < 30 || nlevels(z$sex_f) < 2)
			return(tibble())
		z[[f]] <- std_num(z[[f]])
		rhs <- c(paste0(bt(f), "*sex_f"), bt(covs))
		fit <- tryCatch(coxph(
			as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(rhs, collapse = " + "))),
			z
		), error = function(e) NULL)
		if (is.null(fit))
			return(tibble())
		sm <- coef(summary(fit))
		ir <- grep(":", rownames(sm), value = TRUE)[1]
		if (!length(ir))
			ir <- NA_character_
		st <- map_dfr(levels(z$sex_f), function(sx) {
			zz <- droplevels(z[z$sex_f == sx, , drop = FALSE])
			rhs0 <- c(bt(f), bt(covs))
			ff <- tryCatch(coxph(
				as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(rhs0, collapse = " + "))),
				zz
			), error = function(e) NULL)
			if (is.null(ff))
				return(tibble(sex = sx, beta = NA_real_, se = NA_real_, p = NA_real_))
			ss <- coef(summary(ff))
			if (!f %in% rownames(ss))
				return(tibble(sex = sx, beta = NA_real_, se = NA_real_, p = NA_real_))
			tibble(sex = sx, beta = ss[f, "coef"], se = ss[f, "se(coef)"], p = ss[f, "Pr(>|z|)"])
		}) |>
			pivot_wider(names_from = sex, values_from = c(beta, se, p), names_glue = "{.value}_sex{sex}")
		bind_cols(tibble(feature = f, n = nrow(z), events = sum(z[[evar]] == 1), beta_interaction = if (!is.na(ir))
			sm[ir, "coef"] else NA_real_, se_interaction = if (!is.na(ir))
			sm[ir, "se(coef)"] else NA_real_, p_interaction = if (!is.na(ir))
			sm[ir, "Pr(>|z|)"] else NA_real_), st)
	})
	if (!nrow(res))
		res <- tibble(
			feature = character(), n = integer(), events = integer(), beta_interaction = numeric(), se_interaction = numeric(),
			p_interaction = numeric()
		)
	res <- res |>
		mutate(FDR_interaction = p.adjust(p_interaction, "BH"), interaction_contrast = if (length(sex_levels) >=
			2)
			paste0(sex_name(sex_levels[2]), " minus ", sex_name(sex_levels[1])) else NA_character_) |>
		arrange(p_interaction)
	write_raw_csv(res, "c4.sex_interaction.csv", rawdir)
	plot_c4_sex_interactions(res, layer, outdir)
	finalize_outputs(LE8_JOB, outdir)
	res
}
run_c4_le8_interactions <- function(layer) {
	outdir <- if (layer == "prot") out.prot else out.met
	setwd2(outdir)
	rawdir <- le8_job_dir(outdir, LE8_JOB)
	dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	interaction_cache <- file.path(rawdir, "c4.interactions.res.rds")
	if (LE8_REUSE_RESULTS || cache_valid(interaction_cache)) {
		obj <- le8_figure_result(interaction_cache, c("results", "surfaces"))
		plot_c4_le8_interactions(obj$results, obj$surfaces, outdir)
		finalize_outputs(LE8_JOB, outdir)
		return(invisible(obj))
	}
	need <- unique(c(
		"eid", "ethnic.c", vars.basic, vars.le8, "birth_date", "date_attend", "date_lost", "date_death",
		paste0("fod_icd10_", Y)
	))
	dat <- read_all(need) |>
		filter_analysis_cohort() |>
		make_outcome(Y)
	tvar <- paste0(Y, ".t2e")
	evar <- paste0(Y, ".Yt2e")
	comps <- intersect(vars.le8, names(dat))
	covs <- intersect(vars.basic, names(dat))
	if (length(comps) < 2)
		stop("C4 interactions require at least two LE8 component variables.", call. = FALSE)
	interaction_formula <- function(x, z, cv) as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(c(paste0(
		bt(x),
		"*", bt(z)
	), bt(cv)), collapse = " + ")))
	res <- map_dfr(combn(comps, 2, simplify = FALSE), function(v) {
		x <- v[[1]]
		z <- v[[2]]
		cv <- setdiff(covs, c(x, z))
		d <- dat[, unique(c(tvar, evar, x, z, cv)), drop = FALSE]
		d <- d[complete.cases(d), , drop = FALSE]
		if (nrow(d) < 1000 || sum(d[[evar]] == 1) < 50)
			return(tibble())
		d[[x]] <- std_num(d[[x]])
		d[[z]] <- std_num(d[[z]])
		fit <- tryCatch(coxph(interaction_formula(x, z, cv), d), error = function(e) NULL)
		if (is.null(fit))
			return(tibble())
		sm <- coef(summary(fit))
		ir <- grep(":", rownames(sm), value = TRUE)[1]
		if (is.na(ir) || !ir %in% rownames(sm))
			return(tibble())
		tibble(xvar = x, zvar = z, x = sub("\\.pts$", "", x), z = sub("\\.pts$", "", z), n = nrow(d), events = sum(d[[evar]] ==
			1), beta_interaction = sm[ir, "coef"], se_interaction = sm[ir, "se(coef)"], p_interaction = sm[
			ir,
			"Pr(>|z|)"
		])
	})
	if (!nrow(res))
		res <- tibble(
			x = character(), z = character(), xvar = character(), zvar = character(), n = integer(),
			events = integer(), beta_interaction = numeric(), se_interaction = numeric(), p_interaction = numeric()
		)
	res <- res |>
		mutate(FDR_interaction = p.adjust(p_interaction, "BH"), pair = paste(coalesce(LE8_LABS[x], x), "×", coalesce(
			LE8_LABS[z],
			z
		))) |>
		arrange(p_interaction)
	write_raw_csv(res, "c4.LE8_pairwise_interactions.csv", rawdir)
	top <- res |>
		slice_head(n = 4)
	grid <- if (!nrow(top))
		tibble() else map_dfr(seq_len(nrow(top)), function(i) {
		rw <- top[i, , drop = FALSE]
		x <- rw$x[[1]]
		z <- rw$z[[1]]
		xvar <- rw$xvar[[1]]
		zvar <- rw$zvar[[1]]
		cv <- setdiff(covs, c(xvar, zvar))
		d <- dat[, unique(c(tvar, evar, xvar, zvar, cv)), drop = FALSE]
		d <- d[complete.cases(d), , drop = FALSE]
		if (nrow(d) < 1000 || sum(d[[evar]] == 1) < 50)
			return(tibble())
		d[[xvar]] <- std_num(d[[xvar]])
		d[[zvar]] <- std_num(d[[zvar]])
		fit <- tryCatch(coxph(interaction_formula(xvar, zvar, cv), d), error = function(e) NULL)
		if (is.null(fit))
			return(tibble())
		g <- expand_grid(xvalue = seq( - 2, 2, length.out = 45), zvalue = seq( - 2, 2, length.out = 45))
		nd <- g
		names(nd)[1 : 2] <- c(xvar, zvar)
		for (v in cv) if (is.numeric(d[[v]]))
			nd[[v]] <- median(d[[v]], na.rm = TRUE) else {
			mode_value <- names(sort(table(d[[v]]), decreasing = TRUE))[1]
			nd[[v]] <- if (is.factor(d[[v]]))
				factor(mode_value, levels = levels(d[[v]])) else mode_value
		}
		lp <- tryCatch(as.numeric(predict(fit, newdata = nd, type = "lp")), error = function(e) rep(NA_real_, nrow(nd)))
		g |>
			mutate(pair = rw$pair[[1]], relative_hazard = exp(lp - median(lp, na.rm = TRUE)))
	}) |>
		filter(is.finite(relative_hazard))
	plot_c4_le8_interactions(res, grid, outdir)
	saveRDS(list(results = res, surfaces = grid), interaction_cache, compress = "xz")
	finalize_outputs(LE8_JOB, outdir)
}

plot_c4_sex_interactions <- function(res, layer, outdir) {
	sex_levels <- sub("^beta_sex", "", grep("^beta_sex", names(res), value = TRUE))
	sex_name <- function(x) if (identical(sex_levels, c("0", "1")))
		c(`0` = "Female", `1` = "Male")[[as.character(x)]] else paste("Sex category", x)
	bcols <- grep("^beta_sex", names(res), value = TRUE)
	if (length(bcols) >= 2) {
		x <- bcols[1]
		y <- bcols[2]
		res$label <- ifelse(min_rank(res$p_interaction) <= 20, res$feature, NA_character_)
		xl <- sex_name(sub("beta_sex", "", x))
		yl <- sex_name(sub("beta_sex", "", y))
		pA <- ggplot(res, aes(.data[[x]], .data[[y]])) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
			geom_hline(yintercept = 0, color = "grey75") +
			geom_vline(xintercept = 0, color = "grey75") +
			geom_point(aes(color = FDR_interaction <
				0.05), size = 2, alpha = 0.8) +
			ggrepel::geom_text_repel(aes(label = label),
				size = 2.7, seed = 20,
				max.overlaps = 25, na.rm = TRUE
			) +
			scale_color_manual(
				values = c(`TRUE` = "#D95F02", `FALSE` = "grey70"),
				guide = "none"
			) +
			labs(
				title = paste0("A. Sex-specific ", layer, " associations"), subtitle = "Orange: interaction FDR < 0.05",
				x = paste0("Log-HR per SD in ", xl), y = paste0("Log-HR per SD in ", yl)
			) +
			theme_5c(11)
		top <- res |>
			slice_head(n = 40) |>
			mutate(
				feature = factor(feature, levels = rev(feature)), lo = beta_interaction - 1.96 * se_interaction,
				hi = beta_interaction + 1.96 * se_interaction
			)
		contrast_subtitle <- unique(na.omit(res$interaction_contrast))
		contrast_subtitle <- if (length(contrast_subtitle))
			contrast_subtitle[[1]] else "Second category minus reference"
		pB <- ggplot(top, aes(beta_interaction, feature)) +
			geom_vline(xintercept = 0, color = "grey60") +
			geom_errorbarh(aes(
				xmin = lo,
				xmax = hi
			), height = 0.12) +
			geom_point(aes(color = FDR_interaction < 0.05), size = 2) +
			scale_color_manual(values = c(
				`TRUE` = "#D95F02",
				`FALSE` = "#4C78A8"
			), guide = "none") +
			labs(
				title = "B. Omic-by-sex interaction", subtitle = contrast_subtitle,
				x = "Interaction log-hazard coefficient", y = NULL
			) +
			forest_theme(9)
		save_plot(pA / pB + plot_layout(heights = c(0.8, 1.2)), "c4.Fig20.sex_interaction.png", 10, 11, outdir = outdir)
	} else save_plot(blank_plot("Sex interaction", "Sex-specific estimates could not be formed"), "c4.Fig20.sex_interaction.png",
		9, 5,
		outdir = outdir
	)
	write_xlsx2(list(sex_interaction = res), "c4.Fig20.sex_interaction.out.xlsx")
	invisible(res)
}

plot_c4_le8_interactions <- function(res, grid, outdir) {
	mat <- bind_rows(res |>
		transmute(row = x, col = z, beta = beta_interaction, p = p_interaction, FDR = FDR_interaction), res |>
		transmute(row = z, col = x, beta = beta_interaction, p = p_interaction, FDR = FDR_interaction)) |>
		mutate(row = factor(row, levels = rev(names.le8)), col = factor(col, levels = names.le8), label = ifelse(FDR <
			0.05, "*", ifelse(p < 0.05, "·", "")))
	pA <- if (!nrow(mat))
		blank_plot("A. All pairwise LE8 interactions", "No interaction model met the sample/event requirements") else ggplot(mat, aes(col, row, fill = cap(beta, 0.25))) +
		geom_tile(color = "white") +
		geom_text(aes(label = label),
			fontface = "bold", size = 5
		) +
		scale_fill_gradient2(
			low = "#2C7FB8", mid = "white", high = "#D7301F", midpoint = 0,
			name = "Interaction β"
		) +
		labs(
			title = "A. All pairwise LE8 interactions", subtitle = "* FDR < 0.05; · nominal P < 0.05",
			x = NULL, y = NULL
		) +
		theme_5c(10) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1))
	pB <- if (!nrow(grid))
		blank_plot("B. Interaction surfaces", "No finite adjusted predictions were available") else ggplot(grid, aes(xvalue, zvalue, fill = pmin(relative_hazard, 3))) +
		geom_raster() +
		geom_contour(aes(z = relative_hazard),
			color = "white", alpha = 0.55, bins = 6
		) +
		facet_wrap( ~ pair, nrow = 1) +
		scale_fill_viridis_c(
			name = "Relative hazard",
			option = "C"
		) +
		labs(
			title = "B. Joint LE8 risk surfaces", subtitle = "Axes are standardized LE8 component scores; surfaces are adjusted associations",
			x = "First LE8 component (SD)", y = "Second LE8 component (SD)"
		) +
		theme_5c(9)
	save_plot(pA / pB + plot_layout(heights = c(0.9, 1.1)), "c4.Fig21.LE8_component_interactions.png", 13, 10,
		outdir = outdir
	)
	write_xlsx2(list(interaction_tests = res, prediction_surfaces = grid), "c4.Fig21.LE8_component_interactions.out.xlsx")
}


# Nonlinear dose-response and nested-CV LE8 penalty scoring.
plot_c4_nonlin <- function(res, curves, nadir, outdir) {
	pA <- if (!nrow(res))
		blank_plot("A. Evidence for non-linearity") else ggplot(res, aes( - log10(pmax(p_nonlin, 1e-300)), fct_reorder(variable, p_nonlin, .desc = TRUE))) +
		geom_vline(
			xintercept =  - log10(0.05 / nrow(res)),
			linetype = 2, color = "grey55"
		) +
		geom_segment(aes(x = 0, xend =  - log10(pmax(p_nonlin, 1e-300)), yend = variable),
			color = "grey75"
		) +
		geom_point(aes(color = FDR_nonlin < 0.05), size = 2) +
		scale_color_manual(values = c(
			`TRUE` = "#D95F02",
			`FALSE` = "#4C78A8"
		), guide = "none") +
		labs(
			title = "A. Likelihood-ratio evidence for non-linearity",
			subtitle = "Spline (4 df) versus linear Cox term; dashed line is Bonferroni 0.05/N", x = expression( - log[10](P[nonlinear])),
			y = NULL
		) +
		theme_5c(9)
	pB <- if (!nrow(curves))
		blank_plot("B. Non-linear risk chart", "No variable met the sample/event requirements") else ggplot(curves, aes(x, HR)) +
		geom_hline(yintercept = 1, linetype = 2, color = "grey55") +
		geom_ribbon(aes(
			ymin = lo,
			ymax = hi
		), fill = "grey80", alpha = 0.55) +
		geom_line(color = "#2C7FB8", linewidth = 1) +
		geom_point(
			data = nadir,
			aes(nadir_x, nadir_HR), inherit.aes = FALSE, color = "#D95F02", size = 1.7
		) +
		facet_wrap( ~ variable,
			scales = "free_x",
			ncol = 3
		) +
		scale_y_log10() +
		labs(
			title = "B. Non-linear LE8 and omics-score risk chart", subtitle = "Natural-spline HRs (95% CI), median reference; orange point is the fitted minimum within the central 96%",
			x = NULL, y = "Hazard ratio (log scale)"
		) +
		theme_5c(9)
	p <- pA | pB
	save_plot(p, "c4.Fig22.spline_patterns.png", 12, 9, outdir = outdir)
	write_xlsx2(list(nonlin_tests = res, spline_curves = curves, nadir = nadir), "c4.Fig22.spline_patterns.out.xlsx")
}

plot_c4_penalty <- function(perf, choose, inner_perf, summ, hr, full, best_gamma, outdir) {
	# Legacy caches included adjustment coefficients with an NA quintile label.  Display only the four intended
	# risk-group contrasts; retain the source cache.
	if (nrow(hr))
		hr <- hr |>
			filter(str_detect(term, "^risk_groupRisk Q[2-5]$"), is.finite(HR), is.finite(lo), is.finite(hi), lo >
				0) |>
			mutate(group = factor(str_extract(term, "[2-5]$"), levels = 2 : 5))
	pA <- if (!nrow(perf))
		blank_plot("Cross-validated LE8 aggregation algorithms") else ggplot(perf, aes(cindex, fct_reorder(rule, cindex, mean), color = rule)) +
		geom_point(
			position = position_jitter(height = 0.08),
			alpha = 0.6
		) +
		stat_summary(fun = mean, geom = "point", size = 3) +
		labs(
			title = "A. Nested-CV LE8 aggregation algorithms",
			x = "Held-out C-index", y = NULL, color = NULL
		) +
		theme_5c(10) +
		theme(legend.position = "none")
	pB <- if (!nrow(choose))
		blank_plot("B. Penalty selected within each outer fold") else ggplot(choose, aes(chosen_gamma, factor(fold))) +
		geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
		geom_point(size = 3, color = "#D95F02") +
		scale_x_continuous(limits = range(c(0, inner_perf$gamma, choose$chosen_gamma),
			na.rm = TRUE
		)) +
		labs(
			title = paste0("B. Inner-CV penalty choices; modal gamma = ", best_gamma), subtitle = "Gamma = 0 retains the arithmetic mean without an additional penalty",
			x = "Selected gamma", y = "Outer validation fold"
		) +
		theme_5c(10)
	pC <- if (!nrow(hr))
		blank_plot("C. Disease risk across penalized-score quintiles") else ggplot(hr, aes(HR, group)) +
		geom_vline(xintercept = 1, color = "grey60") +
		geom_errorbarh(aes(
			xmin = lo,
			xmax = hi
		), height = 0.15) +
		geom_point(size = 2, color = "#4C78A8") +
		scale_x_log10() +
		labs(
			title = "C. Disease risk across penalized-score quintiles",
			x = "Hazard ratio versus lowest-risk quintile", y = "Higher-risk quintile"
		) +
		theme_5c(10)
	save_plot(pA / (pB | pC) + plot_layout(heights = c(0.8, 1.1)), "c4.Fig23.pass_fail_penalty.png", 13, 9,
		outdir = outdir
	)
	write_xlsx2(
		list(CV_by_fold = perf, selected_penalty = choose, inner_CV = inner_perf, summary = summ, risk_group_HR = hr),
		"c4.Fig23.pass_fail_penalty.out.xlsx"
	)
}

run_c4_nonlin <- function(layer) {
	outdir <- if (layer == "prot") out.prot else out.met
	setwd2(outdir)
	rawdir <- le8_job_dir(outdir, LE8_JOB)
	dir.create(rawdir, recursive = TRUE, showWarnings = FALSE)
	nonlin_cache <- file.path(rawdir, "c4.nonlin.res.rds")
	penalty_cache <- file.path(rawdir, "c4.penalty.res.rds")
	if (LE8_REUSE_RESULTS || (cache_valid(nonlin_cache) && cache_valid(penalty_cache))) {
		nonlin <- le8_figure_result(nonlin_cache, c("tests", "curves", "nadir"))
		penalty <- le8_figure_result(penalty_cache, c("performance", "selected", "inner", "summary", "risk_HR"))
		best_gamma <- if (nrow(penalty$selected))
			as.numeric(names(sort(table(penalty$selected$chosen_gamma), decreasing = TRUE))[1]) else 0
		full <- penalty$plot_scores
		# Selection decisions are already cached. A loess over a few repeated integer failure counts was
		# numerically unstable and did not evaluate the penalty.  Plot the actual nested-CV choices; no cohort
		# reload is needed for this view.
		plot_c4_nonlin(nonlin$tests, nonlin$curves, nonlin$nadir, outdir)
		plot_c4_penalty(
			penalty$performance, penalty$selected, penalty$inner, penalty$summary, penalty$risk_HR, full,
			best_gamma, outdir
		)
		finalize_outputs(LE8_JOB, outdir)
	} else {
		need <- unique(c(
			"eid", "ethnic.c", vars.basic, vars.le8, "sleep_duration", "sleep.hours", "sleep", "birth_date",
			"date_attend", "date_lost", "date_death", paste0("fod_icd10_", Y)
		))
		dat <- read_all(need) |>
			filter_analysis_cohort() |>
			make_outcome(Y)
		tvar <- paste0(Y, ".t2e")
		evar <- paste0(Y, ".Yt2e")
		covs <- intersect(vars.basic, names(dat))

		load_oof_score <- function(layer) {
			score_name <- paste0(layer, "_oof_score")
			empty_score <- function() tibble::tibble(eid = dat$eid[0])
			od <- if (layer == "prot")
				out.prot else out.met
			f <- file.path(le8_job_dir(od, "final_prediction"), "res.rds")
			if (!file.exists(f))
				return(empty_score())
			z <- readRDS(f)
			pr <- z$scores$rows %||% tibble()
			if (!all(c("method", "eid", "score_z", "split") %in% names(pr)))
				return(empty_score())
			preferred <- c("C4 YSplus", "Parsimonious", "Pradeep / glmnet")
			best <- preferred[preferred %in% unique(pr$method)][1]
			if (!length(best) || is.na(best))
				return(empty_score())
			# Training predictions are out-of-fold and validation predictions are fully held out; prevalent rows are
			# excluded from this incident-risk analysis.
			sc <- pr |>
				filter(method == best, split %in% c("training", "validation"), is.finite(score_z)) |>
				transmute(eid, score = score_z) |>
				distinct(eid, .keep_all = TRUE)
			names(sc)[2] <- score_name
			sc
		}
		dat <- left_join(dat, load_oof_score(layer),
			by = "eid"
		)

		vars0 <- unique(c(
			intersect(c("sleep_duration", "sleep.hours", "sleep"), names(dat)), intersect(vars.le8, names(dat)),
			grep("_oof_score$", names(dat), value = TRUE)
		))
		fit_nonlin <- function(v, return_curve = FALSE) {
			cv <- setdiff(covs, v)
			d <- dat[, unique(c(tvar, evar, v, cv)), drop = FALSE]
			d <- d[complete.cases(d), , drop = FALSE]
			d[[v]] <- suppressWarnings(as.numeric(d[[v]]))
			if (nrow(d) < 1000 || sum(d[[evar]] == 1) < 50 || !is.finite(sd(d[[v]])) || sd(d[[v]]) == 0)
				return(tibble())
			rhs_lin <- c(paste0("scale(", bt(v), ")"), bt(cv))
			rhs_spl <- c(paste0("splines::ns(", bt(v), ",df=4)"), bt(cv))
			fl <- tryCatch(coxph(
				as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(rhs_lin, collapse = " + "))),
				d
			), error = function(e) NULL)
			fs <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(rhs_spl, collapse = " + "))),
				d,
				x = TRUE
			), error = function(e) NULL)
			if (is.null(fl) || is.null(fs))
				return(tibble())
			lr <- tryCatch(anova(fl, fs, test = "LRT"), error = function(e) NULL)
			pnon <- if (!is.null(lr) && nrow(lr) >= 2)
				lr$`Pr(>|Chi|)`[2] else NA_real_
			if (!return_curve)
				return(tibble(
					variable = v, n = nrow(d), events = sum(d[[evar]] == 1), p_nonlin = pnon, AIC_linear = AIC(fl),
					AIC_spline = AIC(fs)
				))
			xg <- seq(quantile(d[[v]], 0.02), quantile(d[[v]], 0.98), length.out = 120)
			nd <- data.frame(xg)
			names(nd) <- v
			for (cc in cv) {
				if (is.numeric(d[[cc]]))
					nd[[cc]] <- median(d[[cc]], na.rm = TRUE) else nd[[cc]] <- factor(names(sort(table(d[[cc]]), decreasing = TRUE))[1], levels = levels(d[[cc]]))
			}
			nd_ref <- nd[1, , drop = FALSE]
			nd_ref[[v]] <- median(d[[v]], na.rm = TRUE)
			tt <- delete.response(terms(fs))
			X <- model.matrix(tt, data = nd, xlev = fs$xlevels)
			X0 <- model.matrix(tt, data = nd_ref, xlev = fs$xlevels)
			cn <- names(coef(fs))
			X <- X[, intersect(cn, colnames(X)), drop = FALSE]
			X0 <- X0[, colnames(X), drop = FALSE]
			b <- coef(fs)[colnames(X)]
			V <- vcov(fs)[colnames(X), colnames(X), drop = FALSE]
			XD <- sweep(X, 2, as.numeric(X0[1, ]), "-")
			lp <- as.numeric(XD %*% b)
			se <- sqrt(pmax(0, rowSums((XD %*% V) * XD)))
			tibble(variable = v, x = xg, HR = exp(lp), lo = exp(lp - 1.96 * se), hi = exp(lp + 1.96 * se), p_nonlin = pnon)
		}
		res <- map_dfr(vars0, fit_nonlin, return_curve = FALSE)
		if (!nrow(res))
			res <- tibble(
				variable = character(), n = integer(), events = integer(), p_nonlin = numeric(), AIC_linear = numeric(),
				AIC_spline = numeric()
			)
		res <- res |>
			mutate(FDR_nonlin = p.adjust(p_nonlin, "BH")) |>
			arrange(p_nonlin)
		curves <- map_dfr(vars0, fit_nonlin, return_curve = TRUE)
		if (!nrow(curves))
			curves <- tibble(
				variable = character(), x = numeric(), HR = numeric(), lo = numeric(), hi = numeric(),
				p_nonlin = numeric()
			)
		curves <- curves |>
			left_join(res |>
				select(variable, FDR_nonlin), by = "variable")
		nadir <- if (!nrow(curves))
			tibble() else curves |>
			group_by(variable) |>
			slice_min(HR, n = 1, with_ties = FALSE) |>
			ungroup() |>
			transmute(variable, nadir_x = x, nadir_HR = HR, p_nonlin, FDR_nonlin)
		res <- res |>
			left_join(nadir |>
				select(variable, nadir_x, nadir_HR), by = "variable")
		write_raw_csv(res, "c4.nonlin_tests.csv", rawdir)
		write_raw_csv(curves, "c4.nonlin_curves.csv", rawdir)
		plot_c4_nonlin(res, curves, nadir, outdir)
		saveRDS(list(tests = res, curves = curves, nadir = nadir), nonlin_cache, compress = "xz")

		# Nested-CV alternatives to the simple LE8 sum, including bottleneck and cross-validated penalty rules.
		# This is part of the same C4 job and output tree.
		penalty_dat <- dat
		comps <- intersect(vars.le8, names(penalty_dat))
		covs <- intersect(vars.basic, names(penalty_dat))
		if (length(comps) < 2)
			stop("C4 penalty analysis requires at least two LE8 component variables.", call. = FALSE)
		penalty_dat <- penalty_dat[complete.cases(penalty_dat[, unique(c(tvar, evar, comps, covs)), drop = FALSE]) &
			penalty_dat[[tvar]] > 0, , drop = FALSE]
		K <- as.integer(Sys.getenv("C4_PENALTY_FOLDS", unset = "5"))
		K_INNER <- as.integer(Sys.getenv("C4_PENALTY_INNER_FOLDS", unset = "4"))
		gammas <- c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3)
		gamma_name <- function(g) paste0("penalty_", gsub("\\.", "p", format(g, scientific = FALSE, trim = TRUE)))
		fit_scaler <- function(x) {
			lo <- vapply(x, min, numeric(1), na.rm = TRUE)
			hi <- vapply(x, max, numeric(1), na.rm = TRUE)
			list(lo = lo, span = pmax(hi - lo, 1e-08))
		}
		apply_scaler <- function(x, scaler) {
			z <- sweep(as.matrix(x), 2, scaler$lo, "-")
			z <- sweep(z, 2, scaler$span, "/")
			pmin(pmax(z, 0), 1)
		}
		add_scores <- function(d, z, gamma_values = gammas) {
			d$LE8_mean <- rowMeans(z)
			d$LE8_min <- apply(z, 1, min)
			d$LE8_geomean <- exp(rowMeans(log(pmax(z, 0.01))))
			d$n_fail <- rowSums(z < 0.4)
			d$LE8_passfail <- d$LE8_mean - 0.08 * pmax(d$n_fail - 1, 0)
			for (g in gamma_values) d[[gamma_name(g)]] <- d$LE8_mean - g * pmax(d$n_fail - 1, 0) ^ 2
			d
		}
		prepare_pair <- function(tr, te) {
			sc <- fit_scaler(tr[, comps, drop = FALSE])
			list(tr = add_scores(tr, apply_scaler(tr[, comps, drop = FALSE], sc)), te = add_scores(te, apply_scaler(te[,
				comps,
				drop = FALSE
			], sc)))
		}
		fit_test_score <- function(tr, te, score) {
			m <- mean(tr[[score]], na.rm = TRUE)
			ss <- sd(tr[[score]], na.rm = TRUE)
			if (!is.finite(ss) || ss == 0)
				ss <- 1
			tr[[score]] <- (tr[[score]] - m) / ss
			te[[score]] <- (te[[score]] - m) / ss
			fit <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(score, covs)),
				collapse = " + "
			))), tr), error = function(e) NULL)
			if (is.null(fit))
				return(tibble(score = score, beta = NA_real_, HR = NA_real_, p = NA_real_, cindex = NA_real_))
			lp <- tryCatch(as.numeric(predict(fit, newdata = te, type = "lp")), error = function(e) rep(NA_real_, nrow(te)))
			sm <- coef(summary(fit))
			tibble(
				score = score, beta = sm[score, "coef"], HR = exp(sm[score, "coef"]), p = sm[score, "Pr(>|z|)"],
				cindex = fcidx(Surv(te[[tvar]], te[[evar]]), lp)
			)
		}
		select_gamma_inner <- function(tr, seed) {
			folds <- make_folds(tr, evar, K_INNER, seed)
			rows <- lapply(seq_len(K_INNER), function(fd) {
				a <- tr[folds != fd, , drop = FALSE]
				b <- tr[folds == fd, , drop = FALSE]
				pp <- prepare_pair(a, b)
				map_dfr(gammas, function(g) fit_test_score(pp$tr, pp$te, gamma_name(g)) |>
					mutate(gamma = g, inner_fold = fd))
			})
			perf <- bind_rows(rows) |>
				group_by(gamma) |>
				summarise(mean_cindex = mean(cindex, na.rm = TRUE), valid_folds = sum(is.finite(cindex)), .groups = "drop") |>
				arrange(desc(mean_cindex), gamma)
			if (!nrow(perf) || !any(is.finite(perf$mean_cindex)))
				return(list(gamma = 0, performance = perf))
			list(gamma = perf$gamma[which.max(perf$mean_cindex)], performance = perf)
		}
		fold <- make_folds(penalty_dat, evar, K, SEED)
		rows <- list()
		chosen <- list()
		inner_all <- list()
		rr <- 0L
		base_rules <- c("LE8_mean", "LE8_min", "LE8_geomean", "LE8_passfail")
		for (fd in seq_len(K)) {
			tr0 <- penalty_dat[fold != fd, , drop = FALSE]
			te0 <- penalty_dat[fold == fd, , drop = FALSE]
			sel <- select_gamma_inner(tr0, SEED + fd)
			inner_best <- if (nrow(sel$performance) && any(is.finite(sel$performance$mean_cindex)))
				max(sel$performance$mean_cindex, na.rm = TRUE) else NA_real_
			chosen[[fd]] <- tibble(fold = fd, chosen_gamma = sel$gamma, chosen_penalty = gamma_name(sel$gamma), inner_mean_cindex = inner_best)
			inner_all[[fd]] <- sel$performance |>
				mutate(outer_fold = fd)
			pp <- prepare_pair(tr0, te0)
			for (score in c(base_rules, gamma_name(sel$gamma))) {
				rr <- rr + 1L
				rows[[rr]] <- fit_test_score(pp$tr, pp$te, score) |>
					mutate(fold = fd, gamma = ifelse(str_detect(score, "^penalty_"), sel$gamma, NA_real_), rule = ifelse(str_detect(
						score,
						"^penalty_"
					), "CV-selected big penalty", score))
			}
		}
		perf <- bind_rows(rows)
		choose <- bind_rows(chosen)
		inner_perf <- bind_rows(inner_all)
		summ <- perf |>
			group_by(rule) |>
			summarise(mean_cindex = mean(cindex, na.rm = TRUE), sd_cindex = sd(cindex, na.rm = TRUE), mean_HR = mean(HR,
				na.rm = TRUE
			), .groups = "drop") |>
			arrange(desc(mean_cindex))
		write_raw_csv(perf, "c4.penalty_CV_by_fold.csv", rawdir)
		write_raw_csv(choose, "c4.penalty_selected_gamma.csv", rawdir)
		write_raw_csv(inner_perf, "c4.penalty_inner_CV.csv", rawdir)
		write_raw_csv(summ, "c4.penalty_summary.csv", rawdir)
		best_gamma <- if (nrow(choose))
			as.numeric(names(sort(table(choose$chosen_gamma), decreasing = TRUE))[1]) else 0
		full <- prepare_pair(penalty_dat, penalty_dat)$tr
		bestp <- gamma_name(best_gamma)
		full$best_penalty <- full[[bestp]]
		full$risk_group <- factor(ntile( - full$best_penalty, 5), levels = 1 : 5, labels = paste0("Risk Q", 1 : 5))
		fit <- tryCatch(coxph(as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ risk_group + ", paste(bt(covs),
			collapse = " + "
		))), full), error = function(e) NULL)
		hr <- if (is.null(fit))
			tibble() else {
			sm <- coef(summary(fit))
			tibble(term = rownames(sm), HR = exp(sm[, "coef"]), lo = exp(sm[, "coef"] - 1.96 * sm[, "se(coef)"]), hi = exp(sm[
,
				"coef"
			] + 1.96 * sm[, "se(coef)"]), p = sm[, "Pr(>|z|)"]) |>
				filter(str_detect(term, "^risk_groupRisk Q[2-5]$")) |>
				mutate(group = str_extract(term, "[2-5]$"), group = factor(group, levels = 2 : 5))
		}
		plot_c4_penalty(perf, choose, inner_perf, summ, hr, full, best_gamma, outdir)
		saveRDS(list(performance = perf, selected = choose, inner = inner_perf, summary = summ, risk_HR = hr, plot_scores = full[
,
			c("n_fail", "best_penalty")
		]), penalty_cache, compress = "xz")
		finalize_outputs(LE8_JOB, outdir)
	}
}

if (prot_DO) {
	invisible(le8_stage("C4/connect/protein", run_c4_layer("protein")))
	invisible(le8_stage("C4/sex-interactions/protein", run_c4_sex_interactions("protein")))
	gc(full = TRUE)
}
if (met_DO) {
	invisible(le8_stage("C4/connect/metabolite", run_c4_layer("metabolite")))
	invisible(le8_stage("C4/sex-interactions/metabolite", run_c4_sex_interactions("metabolite")))
	gc(full = TRUE)
}
for (layer in c(if (prot_DO) "prot", if (met_DO) "met")) {
	invisible(le8_stage(paste0("C4/LE8-interactions/", layer), run_c4_le8_interactions(layer)))
	invisible(le8_stage(paste0("C4/nonlinear-and-penalty/", layer), run_c4_nonlin(layer)))
}
}
