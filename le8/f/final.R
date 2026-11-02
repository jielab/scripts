# Final reports and optional joint validation. Basic prediction belongs to C1.
# Each mode runs in a separate process and reads the same external analysis tree.
.final_mode <- Sys.getenv("LE8_FINAL_MODE", unset = if (Sys.getenv("LE8_JOB") == "final_prediction") "pipeline" else "report")
if (!.final_mode %in% c("report", "joint", "pipeline")) stop("Unknown LE8_FINAL_MODE: ", .final_mode)
if (.final_mode == "joint") {
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
		savefig <- function(p, n, w = 11, h = 7) bounded_ggsave(file.path(root, n), p, width = w, height = h, dpi = 200)
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
		source(file.path(fdir, "c1.correlate.R"))
		inputs <- character()
		audit <- list()
		map$pgs_feature <- NA_character_
		for (layer in c("prot", "met")) {
			f <- find_c1_pgs_file(layer)
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
			bounded_ggsave(file.path(root, paste0(nm, ".png")), p,
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
			bounded_ggsave(file.path(root, paste0(name, ".png")), p,
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
			prediction = "c1_prediction", extra = "final_supp", cell = "c5_cellulation"
		)
		# Regenerate this report at the same location; remove stale prior pages.
		if (dir.exists(root)) unlink(root, recursive = TRUE)
		dir.create(root, recursive = TRUE, showWarnings = FALSE)
		module_dirs <- setNames(file.path(source_root, modules), names(modules))
		module_dirs["prediction"] <- file.path(analysis_root, trait, layer, "c1_correlate", "prediction")
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
		fd <- read("deployed_concept_fidelity", "focus", "c4.focus.deployed_concept_fidelity.csv")
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
			captions["Fig3"] <- paste(captions["Fig3"], "Panels A/B use only outer-training donors for discovery and fit; Yang donors contribute to proxy learning. Post-hoc OLS incremental R² describes panel reconstruction beyond basic covariates. Fidelity of the frozen deployed concepts is exported separately with R², RMSE and calibration; neither estimates intervention responsiveness.")
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
						NA else any(FDR_all < le8_evidence_policy()$mr_fdr, na.rm = TRUE), .groups = "drop") |>
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
				bounded_ggsave(file.path(root, paste0(dest, ".png")), p, width = 19, height = ceiling(length(page) / 2) * 6 +
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
				"cigma_status", "cigma_annotation", "focus_metrics", "focus_contrasts", "focus_proxy_accuracy", "deployed_concept_fidelity", "focus_pillars",
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
