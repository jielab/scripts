# C2: cis/trans MR and DANDELION driver prioritization using gene-level disease evidence.

suppressPackageStartupMessages({
	.c2_fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
	source(file.path(.c2_fdir, "0.common.R"))

	# c2 helpers
	# c2.cause.R
	# C2 individual-level decomposition of an observed omic trait into a simple
	# COJO-weighted PGS component and a non-genetic residual.
	#
	# Auto-discovered inputs:
	#   <UKB_PHE>/Rdata/prot.pgs.rds  (eid, FEATURE.pgs, ...)
	#   <UKB_PHE>/Rdata/met.pgs.rds   (eid, FEATURE.pgs, ...)
	# C2_GENETIC_SCORE_FILE is an optional explicit override. This is an
	# exploratory in-sample decomposition, not external/cross-fitted prediction.

	find_c2_score_file <- function(layer) {
		explicit <- Sys.getenv("C2_GENETIC_SCORE_FILE", unset = "")
		automatic <- file.path(indir, "Rdata", if (layer == "protein") "prot.pgs.rds" else "met.pgs.rds")
		z <- unique(c(explicit, automatic)) ; z <- z[nzchar(z) & file.exists(z) & file.size(z) > 0]
		if (length(z)) normalizePath(z[[1]], winslash = "/", mustWork = FALSE) else NA_character_
	}

	read_c2_scores <- function(file) {
		x <- if (grepl("\\.rds$", file, ignore.case = TRUE)) readRDS(file) else
			data.table::fread(file, showProgress = FALSE, check.names = FALSE)
		x <- as_tibble(x)
		if (!"eid" %in% names(x)) stop("PGS file must contain eid: ", file, call. = FALSE)
		x
	}

	map_c2_score_columns <- function(features, nms) {
		one <- function(f) {
			z <- c(
				paste0(f, ".pgs"), paste0(f, "_pgs"), paste0(f, ".PGS"),
				paste0(f, "_PGS"), paste0(f, "_GRS"), paste0("GRS_", f), f
			)
			hit <- z[z %in% nms] ; if (length(hit)) hit[[1]] else NA_character_
		}
		x <- setNames(vapply(features, one, character(1)), features) ; x[!is.na(x)]
	}

	find_c2_heritability_file <- function(layer) {
		layer_explicit <- Sys.getenv(if (layer == "protein") "C2_PROT_HERITABILITY_FILE" else
			"C2_MET_HERITABILITY_FILE", unset = "")
		explicit <- Sys.getenv("C2_HERITABILITY_FILE", unset = "")
		stem <- if (layer == "protein") "prot" else "met"
		automatic <- c(
			file.path(indir, "Rdata", paste0(stem, ".heritability.csv")),
			file.path(indir, "Rdata", paste0(stem, ".heritability.rds"))
		)
		z <- unique(c(layer_explicit, explicit, automatic)) ; z <- z[nzchar(z) & file.exists(z) & file.size(z) > 0]
		if (length(z)) normalizePath(z[[1]], winslash = "/", mustWork = FALSE) else NA_character_
	}

	read_c2_heritability <- function(layer) {
		f <- find_c2_heritability_file(layer)
		empty <- list(
			data = tibble(feature = character(), snp_h2 = numeric(), snp_h2_se = numeric()),
			status = tibble(status = "unavailable", file = NA_character_, detail = "Optional SNP-heritability file was not found")
		)
		if (is.na(f)) return(empty)
		x <- tryCatch(if (grepl("\\.rds$", f, ignore.case = TRUE)) as_tibble(readRDS(f)) else
			as_tibble(data.table::fread(f, showProgress = FALSE, check.names = FALSE)), error = function(e) e)
		if (inherits(x, "condition")) return(le8_result_update(empty, list(status = tibble(status = "invalid", file = f, detail = conditionMessage(x)))))
		pick <- function(pattern) {
			i <- grep(pattern, names(x), ignore.case = TRUE)[1] ; if (is.na(i)) NA_character_ else names(x)[[i]]
		}
		fcol <- pick("^(feature|term|trait|exposure|protein|metabolite)$")
		hcol <- pick("^(snp_?h2|h2_?snp|h2|heritability)$")
		scol <- pick("^(snp_?h2_?se|h2_?se|heritability_?se|se)$")
		if (is.na(fcol) || is.na(hcol)) return(le8_result_update(empty, list(status = tibble(
			status = "invalid", file = f,
			detail = "Need feature/term/trait and snp_h2/h2/heritability columns"
		))))
		d <- tibble(
			feature = as.character(x[[fcol]]), snp_h2 = suppressWarnings(as.numeric(x[[hcol]])),
			snp_h2_se = if (is.na(scol)) NA_real_ else suppressWarnings(as.numeric(x[[scol]]))
		) |>
			filter(!is.na(feature), nzchar(feature), is.finite(snp_h2)) |>
			mutate(snp_h2 = pmax(0, pmin(1, snp_h2))) |>
			distinct(feature, .keep_all = TRUE)
		list(data = d, status = tibble(
			status = if (nrow(d)) "ok" else "invalid", file = f,
			detail = paste(nrow(d), "features with SNP-heritability estimates")
		))
	}

	component_association <- function(dd, x, covars, tvar, evar, prevalent = FALSE) {
		covars <- intersect(covars, names(dd))
		if (prevalent) {
			need <- unique(c(".prevalent", x, covars)) ; z <- dd[, need, drop = FALSE]
			z <- z[complete.cases(z), , drop = FALSE] ; events <- sum(z$.prevalent == 1)
			if (nrow(z) < 500 || events < 20 || length(unique(z$.prevalent)) < 2)
				return(c(beta = NA, se = NA, p = NA, N = nrow(z), events = events))
			fit <- tryCatch(glm(reformulate(c(x, covars), ".prevalent"), z, family = binomial()), error = function(e) NULL)
		} else {
			need <- unique(c(tvar, evar, x, covars)) ; z <- dd[, need, drop = FALSE]
			z <- z[complete.cases(z), , drop = FALSE]
			z <- z[is.finite(z[[tvar]]) & z[[tvar]] > 0 & z[[evar]] %in% c(0, 1), , drop = FALSE]
			events <- sum(z[[evar]] == 1)
			if (nrow(z) < 500 || events < 20) return(c(beta = NA, se = NA, p = NA, N = nrow(z), events = events))
			ff <- as.formula(paste0(
				"Surv(", bt(tvar), ",", bt(evar), ") ~ ",
				paste(bt(c(x, covars)), collapse = " + ")
			))
			fit <- tryCatch(coxph(ff, z, ties = "efron"), error = function(e) NULL)
		}
		if (is.null(fit)) return(c(beta = NA, se = NA, p = NA, N = nrow(z), events = events))
		sm <- coef(summary(fit)) ; if (!x %in% rownames(sm))
			return(c(beta = NA, se = NA, p = NA, N = nrow(z), events = events))
		if (prevalent) c(
			beta = sm[x, "Estimate"], se = sm[x, "Std. Error"], p = sm[x, "Pr(>|z|)"],
			N = nrow(z), events = events
		) else
			c(
				beta = sm[x, "coef"], se = sm[x, "se(coef)"], p = sm[x, "Pr(>|z|)"],
				N = nrow(z), events = events
			)
	}

	joint_component_association <- function(dd, covars, tvar, evar, prevalent = FALSE) {
		covars <- intersect(covars, names(dd)) ; xs <- c(".genetic_z", ".residual_z")
		if (prevalent) {
			need <- unique(c(".prevalent", xs, covars)) ; z <- dd[, need, drop = FALSE]
			z <- z[complete.cases(z), , drop = FALSE]
			fit <- if (nrow(z) >= 500 && sum(z$.prevalent == 1) >= 20)
				tryCatch(glm(reformulate(c(xs, covars), ".prevalent"), z, family = binomial()), error = function(e) NULL) else NULL
		} else {
			need <- unique(c(tvar, evar, xs, covars)) ; z <- dd[, need, drop = FALSE]
			z <- z[complete.cases(z), , drop = FALSE]
			z <- z[is.finite(z[[tvar]]) & z[[tvar]] > 0 & z[[evar]] %in% c(0, 1), , drop = FALSE]
			ff <- as.formula(paste0(
				"Surv(", bt(tvar), ",", bt(evar), ") ~ ",
				paste(bt(c(xs, covars)), collapse = " + ")
			))
			fit <- if (nrow(z) >= 500 && sum(z[[evar]] == 1) >= 20)
				tryCatch(coxph(ff, z, ties = "efron"), error = function(e) NULL) else NULL
		}
		if (is.null(fit)) return(tibble(component = xs, beta = NA_real_, se = NA_real_, p = NA_real_))
		sm <- coef(summary(fit))
		tibble(
			component = xs,
			beta = vapply(xs, function(x) if (x %in% rownames(sm)) sm[x, if (prevalent) "Estimate" else "coef"] else NA_real_, numeric(1)),
			se = vapply(xs, function(x) if (x %in% rownames(sm)) sm[x, if (prevalent) "Std. Error" else "se(coef)"] else NA_real_, numeric(1)),
			p = vapply(xs, function(x) if (x %in% rownames(sm)) sm[x, "Pr(>|z|)"] else NA_real_, numeric(1))
		)
	}


	plot_component_leadtime <- function(x) {
		if (is.null(x) || !is.data.frame(x) || !nrow(x) || !all(c("beta", "se") %in% names(x)))
			return(blank_plot(
				"PGS and residual lead-time associations",
				"PGS input unavailable; individual genetic decomposition was skipped"
			))
		d <- as_tibble(x) |> filter(is.finite(beta), is.finite(se))
		if (!nrow(d)) return(blank_plot(
			"PGS and residual lead-time associations",
			"No risk-set window estimate was available"
		))
		ord <- d |>
			group_by(feature) |>
			summarise(best_p = safe_min_finite(p), .groups = "drop") |>
			arrange(best_p) |>
			slice_head(n = 12) |>
			pull(feature)
		d <- d |>
			filter(feature %in% ord) |>
			mutate(feature = factor(feature, levels = rev(ord)))
		ggplot(d, aes(lead_mid, beta, color = component, group = component)) +
			geom_hline(yintercept = 0, color = "grey72") +
			geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = component), alpha = .10, color = NA) +
			geom_line(linewidth = .75) +
			geom_point(aes(size = N_case), alpha = .88) +
			facet_wrap( ~ feature, scales = "free_y", ncol = 3) +
			scale_x_continuous(breaks = c(.5, 1, 2, 5, 10, 16), limits = c(0, 16)) +
			labs(
				title = "Risk-set lead-time associations of observed, PGS-predicted and residual components",
				subtitle = "Cases are diagnosed in each window; controls remain observed and disease-free through its upper boundary",
				x = "Years from baseline to diagnosis window", y = "Log odds ratio per 1-SD",
				color = NULL, fill = NULL, size = "Cases"
			) +
			theme_5c(8) +
			theme(legend.position = "bottom")
	}


	# c2.cause.R
	# Effect/CI comparison inspired by Koprulu et al., Cell 2026, Figure 5.
	# Estimates are this run's MR results; the display does not re-estimate MR.
	le8_c2_paired_estimates <- function(mr, layer) {
		classes <- if (layer == 'protein') c('cis', 'trans') else c('local', 'distal')
		d <- mr |> filter(analysis %in% classes)
		if (anyDuplicated(d[, c('exposure', 'analysis')])) stop('Multiple primary MR rows for one exposure/class; specify the primary estimator before plotting')
		for (nm in c('n_IV', 'FDR_analysis', 'LD_status', 'instrument_snps')) if (!nm %in% names(d)) d[[nm]] <- NA
		d <- d |>
			group_by(analysis) |>
			mutate(FDR_analysis = coalesce(FDR_analysis, p.adjust(pval, 'BH'))) |>
			ungroup()
		wide <- d |>
			select(exposure, analysis, b, se, pval, n_IV, FDR_analysis, LD_status, instrument_snps) |>
			pivot_wider(names_from = analysis, values_from =  - c(exposure, analysis))
		for (field in c('b', 'se', 'pval', 'n_IV', 'FDR_analysis', 'LD_status', 'instrument_snps')) {
			for (i in 1 : 2) {
				nm <- paste0(field, '_', classes[i]) ; if (!nm %in% names(wide)) wide[[nm]] <- NA
				wide[[paste0(field, if (i == 1) '_primary' else '_distal')]] <- wide[[nm]]
			}
		}
		wide <- wide |> mutate(
			paired = is.finite(b_primary) & is.finite(b_distal) & is.finite(se_primary) &
				is.finite(se_distal) & se_primary > 0 & se_distal > 0,
			effect_difference = b_primary - b_distal, difference_se = sqrt(se_primary ^ 2 + se_distal ^ 2),
			Q_difference = ifelse(paired, (effect_difference / difference_se) ^ 2, NA_real_),
			I2_difference = ifelse(paired, ifelse(Q_difference > 0, pmax(0, (Q_difference - 1) / Q_difference) * 100, 0), NA_real_),
			p_heterogeneity = ifelse(paired, 2 * pnorm( - abs(effect_difference / difference_se)), NA_real_),
			FDR_heterogeneity = p.adjust(p_heterogeneity, 'BH'),
			heterogeneity = case_when(
				!paired ~ 'Unavailable', p_heterogeneity >= .01 ~ 'No detected difference',
				sign(b_primary) != sign(b_distal) ~ 'Different effects; opposite signs', TRUE ~ 'Different effects; same sign'
			),
			primary_supported = FDR_analysis_primary < .05,
			lo_primary = b_primary - 1.96 * se_primary, hi_primary = b_primary + 1.96 * se_primary,
			lo_distal = b_distal - 1.96 * se_distal, hi_distal = b_distal + 1.96 * se_distal
		)
		wide
	}

	le8_plot_cis_trans_comparison <- function(mr, assoc, layer) {
		cls <- if (layer == 'protein') c('cis', 'trans') else c('local', 'distal')
		wide <- le8_c2_paired_estimates(mr, layer)
		wide <- wide |> left_join(assoc |> select(exposure = term, obs_beta = beta, obs_se = std.error, obs_p = p.value), by = 'exposure')
		focus <- wide |> filter(paired, primary_supported)
		selection <- paste0(cls[1], ' MR FDR < 0.05 with an estimable ', cls[2], ' effect')
		if (!nrow(focus)) {
			focus <- wide |> filter(paired) ; selection <- 'All estimable pairs; no primary-class FDR discovery'
		}
		if (!nrow(focus)) return(list(plot = blank_plot('Paired genetic effect comparison', 'No biomarker has both instrument classes'), wide = wide))
		anchors <- if (layer == 'protein') c('PCSK9', 'LPA', 'GDF15', 'NTPROBNP', 'MMP12') else
			c('L_VLDL_TG.pct', 'L_VLDL_TG', 'Total_TG', 'ApoB', 'Glucose')
		labels <- head(unique(c(intersect(anchors, focus$exposure), focus$exposure[order(focus$pval_primary)])), 10)
		focus$label <- ifelse(focus$exposure %in% labels, focus$exposure, NA_character_)
		pal <- c(
			'Different effects; opposite signs' = '#BB4560', 'Different effects; same sign' = '#DA953F',
			'No detected difference' = '#3B7D9B'
		)
		pa <- ggplot(focus, aes(b_primary, b_distal)) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = 'grey55', linewidth = .5) +
			geom_hline(yintercept = 0, color = 'grey85', linewidth = .4) +
			geom_vline(xintercept = 0, color = 'grey85', linewidth = .4) +
			geom_segment(aes(x = lo_primary, xend = hi_primary, yend = b_distal, color = heterogeneity), alpha = .32, linewidth = .4) +
			geom_segment(aes(y = lo_distal, yend = hi_distal, xend = b_primary, color = heterogeneity), alpha = .32, linewidth = .4) +
			geom_point(aes(color = heterogeneity), size = 2.8, alpha = .95) +
			ggrepel::geom_text_repel(aes(label = label),
				size = 3.1, seed = SEED, max.overlaps = Inf, na.rm = TRUE,
				box.padding = .5, min.segment.length = 0, color = 'grey20'
			) +
			scale_color_manual(values = pal, drop = FALSE) +
			labs(
				title = paste0('A. ', cls[1], ' and ', cls[2], ' effects on ', Y),
				subtitle = paste0(nrow(focus), ' biomarkers | ', selection),
				x = paste0(cls[1], ' MR effect'), y = paste0(cls[2], ' MR effect'), color = NULL
			) +
			theme_5c(11) +
			theme(panel.grid = element_blank())
		# Fixed anchors first, followed by the strongest primary-class signals.
		candidates <- wide |>
			filter(paired) |>
			arrange(pval_primary, exposure)
		examples <- head(unique(c(intersect(anchors, candidates$exposure), candidates$exposure)), 12)
		forest <- bind_rows(lapply(1 : 2, function(i) {
			suffix <- if (i == 1) 'primary' else 'distal'
			candidates |>
				filter(exposure %in% examples) |>
				transmute(exposure,
					class = cls[i],
					b = .data[[paste0('b_', suffix)]], se = .data[[paste0('se_', suffix)]], n_IV = .data[[paste0('n_IV_', suffix)]],
					pval = .data[[paste0('pval_', suffix)]],
					y = match(exposure, rev(examples)) + if (i == 1) .16 else - .16
				)
		})) |> mutate(lo = b - 1.96 * se, hi = b + 1.96 * se)
		row_labels <- candidates[match(rev(examples), candidates$exposure), ] |>
			transmute(label = paste0(exposure, '  [', n_IV_primary, ' / ', n_IV_distal, ']')) |>
			pull(label)
		pb <- ggplot(forest, aes(b, y, color = class)) +
			geom_hline(yintercept = seq(.5, length(examples) + .5, by = 1), color = 'grey94', linewidth = .35) +
			geom_vline(xintercept = 0, color = 'grey65', linetype = 2, linewidth = .5) +
			geom_segment(aes(x = lo, xend = hi, yend = y), linewidth = .8) +
			geom_point(size = 2.6) +
			scale_color_manual(values = setNames(c('#D88D2D', '#337EAB'), cls), breaks = cls) +
			scale_y_continuous(breaks = seq_along(examples), labels = row_labels, expand = expansion(add = .5)) +
			labs(
				title = 'B. Paired estimates and 95% confidence intervals',
				subtitle = paste0('Fixed anchors + strongest ', cls[1], ' signals; [', cls[1], ' / ', cls[2], ' IV counts]'),
				x = 'MR effect (disease GWAS scale)', y = NULL, color = 'Instrument class'
			) +
			theme_5c(11) +
			theme(panel.grid = element_blank())
		wide$shown_scatter <- wide$exposure %in% focus$exposure ; wide$shown_forest <- wide$exposure %in% examples
		wide$heterogeneity_assumption <- 'Approximate independent-estimate comparison; covariance / cross-class LD not modelled'
		caption <- paste0(
			'Points and intervals are the saved primary MR estimates and 95% CIs. Dashed diagonal: equal effects. ',
			'Colour uses exploratory P(diff) < 0.01, with Var(diff) = SE1^2 + SE2^2; cross-class covariance is unavailable. ',
			'Q (1 df), I² and BH-adjusted P values are in c2.cis_trans_comparison.csv. ',
			sum(focus$n_IV_primary == 1 & focus$n_IV_distal == 1), ' / ', nrow(focus), ' scatter pairs use one IV in both classes. ',
			'The current instrument set has no protein-specific / phenome pleiotropy filter. ',
			'Primary-class MR significance alone does not establish colocalization or causal validity.'
		)
		list(plot = (pa | pb) + plot_layout(widths = c(1.05, 1)) + plot_annotation(caption = caption), wide = wide, forest = forest)
	}

	le8_plot_instrument_comparison <- function(mr, layer) {
		cls <- if (layer == 'protein') c('cis', 'trans') else c('local', 'distal')
		d <- mr |> filter(analysis %in% cls)
		paired <- d |>
			select(exposure, analysis, r2_median) |>
			pivot_wider(names_from = analysis, values_from = r2_median)
		if (!all(cls %in% names(paired))) return(blank_plot('Instrument architecture', 'One instrument class is unavailable'))
		paired <- paired |> filter(is.finite(.data[[cls[1]]]), is.finite(.data[[cls[2]]]), .data[[cls[1]]] > 0, .data[[cls[2]]] > 0)
		anchors <- if (layer == 'protein') c('PCSK9', 'LPA', 'GDF15', 'NTPROBNP', 'MMP12') else c('ApoB', 'Glucose', 'Total_TG', 'L_VLDL_TG.pct')
		paired$label <- ifelse(paired$exposure %in% anchors, paired$exposure, NA_character_)
		pa <- ggplot(paired, aes(100 * .data[[cls[1]]], 100 * .data[[cls[2]]])) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = 'grey65') +
			geom_point(color = '#397F9D', alpha = .55, size = 2) +
			ggrepel::geom_text_repel(aes(label = label), size = 3, seed = SEED, na.rm = TRUE, max.overlaps = Inf) +
			scale_x_log10() +
			scale_y_log10() +
			labs(
				title = 'A. Instrument strength across the same biomarkers',
				subtitle = paste0(nrow(paired), ' paired sets; median per-IV variance, without summing across SNPs'),
				x = paste0(cls[1], ' median per-IV partial R² (%)'), y = paste0(cls[2], ' median per-IV partial R² (%)')
			) +
			theme_5c(11) +
			theme(panel.grid = element_blank())
		d <- d |> mutate(status = case_when(
			!is.finite(b) | !is.finite(pval) ~ 'No estimate',
			grepl('fallback', LD_status, fixed = TRUE) ~ 'Single IV: LD unavailable', n_IV == 1 ~ 'Single eligible IV', TRUE ~ 'Multiple retained IVs'
		))
		counts <- d |>
			count(analysis, status) |>
			mutate(analysis = factor(analysis, levels = rev(cls)))
		pb <- ggplot(counts, aes(n, analysis, fill = status)) +
			geom_col(width = .55) +
			geom_text(aes(label = ifelse(n >= 10, n, '')), position = position_stack(vjust = .5), size = 3.5, color = 'white') +
			scale_fill_manual(values = c(
				'Single IV: LD unavailable' = '#BC6575', 'Single eligible IV' = '#6897AD',
				'Multiple retained IVs' = '#44836F', 'No estimate' = '#B7BDC4'
			)) +
			labs(
				title = 'B. What the saved MR models could use', subtitle = 'Retained instrument sets after harmonization and LD handling',
				x = 'Biomarker–instrument-class records', y = NULL, fill = NULL
			) +
			theme_5c(11) +
			theme(panel.grid = element_blank())
		(pa | pb) + plot_annotation(caption = 'Per-IV partial R² is not total genetic variance explained. LD-unavailable fallback is shown explicitly; it is not evidence that the biological architecture is monogenic.')
	}


	# c2.cause.R
	# Restored verbatim plotting functions from le8_review_backup_20260913.
	# Kept alongside the new Figure 5/6 comparisons; no estimator changes.
	le8_plot_legacy_mr_overview <- function(mr, assoc, layer) {
		local_name <- if (layer == "protein") "cis" else "local" ; distal_name <- if (layer == "protein") "trans" else "distal"
		wide <- mr |>
			select(exposure, analysis, b, se, pval) |>
			pivot_wider(names_from = analysis, values_from = c(b, se, pval)) |>
			left_join(assoc |> select(exposure = term, obs_beta = beta, obs_se = std.error, obs_p = p.value), by = "exposure")
		bl <- paste0("b_", local_name) ; sl <- paste0("se_", local_name) ; pl <- paste0("pval_", local_name) ; bd <- paste0("b_", distal_name) ; sd <- paste0("se_", distal_name) ; pd <- paste0("pval_", distal_name)
		a <- wide |>
			filter(is.finite(.data[[bl]]), is.finite(obs_beta)) |>
			mutate(label = ifelse(min_rank(pmin(obs_p, .data[[pl]], na.rm = TRUE)) <= 12, exposure, NA_character_))
		pA <- if (!nrow(a)) blank_plot(paste0("a. Observational versus ", local_name, " MR"), "No paired estimate") else ggplot(a, aes(obs_beta, .data[[bl]])) +
			geom_hline(yintercept = 0, color = "grey70") +
			geom_vline(xintercept = 0, color = "grey70") +
			geom_smooth(method = "lm", se = TRUE, color = "grey35", fill = "grey85", linewidth = .7) +
			geom_point(aes(color = .data[[pl]] < .05), size = 2) +
			ggrepel::geom_text_repel(aes(label = label), size = 2.8, fontface = "bold", seed = 1, max.overlaps = 20, na.rm = TRUE) +
			scale_color_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "grey70"), guide = "none") +
			labs(title = paste0("a. Observational versus ", local_name, "-MR effects"), x = "Observational effect", y = paste0(local_name, " MR effect")) +
			theme_5c(11)
		b <- wide |>
			filter(is.finite(.data[[bl]]), is.finite(.data[[bd]])) |>
			mutate(label = ifelse(min_rank(pmin(.data[[pl]], .data[[pd]], na.rm = TRUE)) <= 12, exposure, NA_character_))
		pB <- if (!nrow(b)) blank_plot(paste0("b. ", local_name, " versus ", distal_name), "No trait had both instrument classes") else ggplot(b, aes(.data[[bl]], .data[[bd]])) +
			geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
			geom_hline(yintercept = 0, color = "grey70") +
			geom_vline(xintercept = 0, color = "grey70") +
			geom_smooth(method = "lm", se = TRUE, color = "grey35", fill = "grey85", linewidth = .7) +
			geom_point(aes(color = sign(.data[[bl]]) == sign(.data[[bd]])), size = 2) +
			ggrepel::geom_text_repel(aes(label = label), size = 2.8, fontface = "bold", seed = 2, max.overlaps = 20, na.rm = TRUE) +
			scale_color_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "#D95F02"), guide = "none") +
			labs(title = paste0("b. ", local_name, " versus ", distal_name, " MR"), x = paste0(local_name, " effect"), y = paste0(distal_name, " effect")) +
			theme_5c(11)
		ev <- bind_rows(
			assoc |> transmute(exposure = term, evidence = "Observational", effect = beta, p = p.value),
			mr |> filter(analysis %in% c(local_name, distal_name)) |>
				transmute(exposure, evidence = ifelse(analysis == local_name, paste0("MR: ", local_name), paste0("MR: ", distal_name)), effect = b, p = pval)
		) |>
			filter(is.finite(p), is.finite(effect)) |>
			group_by(evidence) |>
			mutate(q = p.adjust(p, "BH")) |>
			ungroup()
		top_ev <- ev |>
			group_by(exposure) |>
			summarise(best_p = min(p), .groups = "drop") |>
			slice_min(best_p, n = 28, with_ties = FALSE) |>
			pull(exposure)
		ev <- ev |>
			filter(exposure %in% top_ev) |>
			mutate(
				score = stable_neglog10_p(p),
				signed_score = sign(effect) * pmin(score, 12), significant = q < .05,
				evidence = factor(evidence, levels = c("Observational", paste0("MR: ", local_name), paste0("MR: ", distal_name))),
				exposure = factor(exposure, levels = rev(top_ev))
			)
		pC <- if (!nrow(ev)) blank_plot("c. Evidence matrix") else ggplot(ev, aes(evidence, exposure)) +
			geom_tile(aes(fill = signed_score), color = "white", linewidth = .35) +
			geom_point(data = ev |> filter(significant), shape = 8, size = 2.1, color = "black") +
			scale_fill_gradient2(
				low = "#3F78A8", mid = "white", high = "#C86B4A", midpoint = 0, limits = c( - 12, 12),
				name = "signed\n-log10(P)"
			) +
			labs(
				title = "c. Signed cross-evidence matrix",
				subtitle = "Blue = protective; orange = risk-increasing; asterisk = within-analysis FDR < 0.05; scale capped at 12",
				x = NULL, y = NULL
			) +
			theme_5c(8) +
			theme(legend.position = "right")
		list(plot = (pA | pB) / pC + plot_layout(heights = c(1, .9)), wide = wide)
	}

	le8_plot_legacy_qtl_variance <- function(mr, layer) {
		cls <- if (layer == "protein") c("cis", "trans") else c("local", "distal")
		ps <- map(cls, function(cc) {
			d <- mr |>
				filter(analysis == cc, is.finite(r2_median), n_IV > 0) |>
				mutate(r2_med = 100 * r2_median, r2_lo = 100 * r2_q25, r2_hi = 100 * r2_q75, r2_p90_plot = 100 * r2_p90) |>
				slice_max(r2_p90_plot, n = 20, with_ties = FALSE) |>
				arrange(r2_med) |>
				mutate(exposure = factor(exposure, levels = exposure))
			if (!nrow(d)) return(blank_plot(paste0(toupper(substr(cc, 1, 1)), substr(cc, 2, nchar(cc)), " instruments"), "No valid instrument set"))
			ggplot(d, aes(r2_med, exposure)) +
				geom_segment(aes(x = r2_med, xend = r2_p90_plot, yend = exposure), color = "grey72", linewidth = .7) +
				geom_errorbarh(aes(xmin = r2_lo, xmax = r2_hi), height = .16, color = ifelse(cc %in% c("cis", "local"), "#D95F02", "#12AEB5"), linewidth = .8) +
				geom_point(aes(size = pmin(n_IV, 500)), color = ifelse(cc %in% c("cis", "local"), "#D95F02", "#12AEB5")) +
				geom_text(aes(x = r2_p90_plot, label = sprintf("median %.3f%%; K=%d", r2_med, n_IV)), hjust =  - .05, size = 2.35, fontface = "bold") +
				scale_x_continuous(expand = expansion(mult = c(.02, .34))) +
				scale_size_continuous(range = c(1.6, 4.8), name = "IV count") +
				labs(
					title = paste0(cc, " QTL instruments"), subtitle = "Per-IV partial R²: point = median; interval = IQR; grey tail = 90th percentile",
					x = "Per-instrument partial phenotypic R² (%)", y = NULL
				) +
				theme_5c(8)
		})
		(ps[[1]] | ps[[2]]) +
			plot_annotation(
				title = "QTL instrument-level variance architecture",
				subtitle = "This figure does not sum R² and does not use 1 - product(1 - R²); values cannot saturate merely because K is large",
				caption = "Instrument-level R² distribution",
				theme = theme(
					plot.title = element_text(face = "bold", size = 15),
					plot.subtitle = element_text(size = 10, color = "grey30"),
					plot.caption = element_text(size = 8, color = "grey40", hjust = 1)
				)
			)
	}


	# c2.cause.R
	# Fig6-style triangulation: independent cis/local MR, incident Cox and prevalent
	# logistic results. This adaptation uses within-analysis FDR, not the paper's
	# high-confidence MR + coloc gate. No choice of minimum P across cis/trans.
	le8_c2_observational_data <- function(mr, incident, prevalent, layer) {
		primary <- if (layer == 'protein') 'cis' else 'local'
		m <- mr |> filter(analysis == primary)
		if (anyDuplicated(m$exposure)) stop('Duplicate primary MR estimates')
		if (!'FDR_analysis' %in% names(m)) m$FDR_analysis <- p.adjust(m$pval, 'BH')
		m <- m |> transmute(exposure, b_mr = b, se_mr = se, p_mr = pval, q_mr = FDR_analysis)
		obs <- function(d, suffix) {
			if (anyDuplicated(d$term)) stop('Duplicate observational estimates')
			if (!'FDR' %in% names(d)) d$FDR <- p.adjust(d$p.value, 'BH')
			z <- d |> transmute(exposure = term, b = beta, se = std.error, p = p.value, q = FDR)
			names(z)[ - 1] <- paste0(names(z)[ - 1], '_', suffix) ; z
		}
		z <- full_join(m, obs(incident, 'incident'), by = 'exposure') |> full_join(obs(prevalent, 'prevalent'), by = 'exposure')
		for (k in c('mr', 'incident', 'prevalent')) {
			z[[paste0('available_', k)]] <- is.finite(z[[paste0('b_', k)]]) & is.finite(z[[paste0('se_', k)]]) &
				z[[paste0('se_', k)]] > 0 & is.finite(z[[paste0('q_', k)]])
			z[[paste0('supported_', k)]] <- z[[paste0('available_', k)]] & z[[paste0('q_', k)]] < .05
		}
		z
	}

	le8_c2_support_counts <- function(z) {
		levels <- c('Concordant', 'Discordant', 'No FDR support', 'Unavailable')
		rows <- list()
		for (endpoint in c('incident', 'prevalent')) for (direction in c('MR to observational', 'Observational to MR')) {
			source <- if (direction == 'MR to observational') 'mr' else endpoint
			target <- if (direction == 'MR to observational') endpoint else 'mr'
			d <- z[z[[paste0('supported_', source)]], , drop = FALSE]
			state <- ifelse(!d[[paste0('available_', target)]], 'Unavailable',
				ifelse(!d[[paste0('supported_', target)]], 'No FDR support',
					ifelse(sign(d[[paste0('b_', source)]]) == sign(d[[paste0('b_', target)]]), 'Concordant', 'Discordant')
				)
			)
			counts <- table(factor(state, levels = levels))
			rows[[length(rows) + 1L]] <- tibble(
				endpoint = endpoint, direction = direction, status = levels,
				n = as.integer(counts), denominator = nrow(d), fraction = if (nrow(d)) as.integer(counts) / nrow(d) else NA_real_
			)
		}
		bind_rows(rows)
	}

	le8_plot_mr_incident_prevalent <- function(mr, incident, prevalent, layer) {
		z <- le8_c2_observational_data(mr, incident, prevalent, layer)
		counts <- le8_c2_support_counts(z)
		if (!any(z$available_mr) || !any(z$available_incident) || !any(z$available_prevalent))
			return(list(plot = blank_plot('MR, incident and prevalent comparison', 'One or more result sources unavailable'), data = z, counts = counts))
		pal <- c('Concordant' = '#E37768', 'Discordant' = '#587DAB', 'No FDR support' = '#BDC5C4', 'Unavailable' = '#ECECEC')
		bars <- lapply(c('incident', 'prevalent'), function(ep) {
			d <- counts |>
				filter(endpoint == ep) |>
				mutate(
					status = factor(status, levels = names(pal)),
					row = paste0(ifelse(direction == 'MR to observational', 'MR → observed', 'Observed → MR'), '  (n=', denominator, ')')
				)
			ggplot(d, aes(fraction, row, fill = status)) +
				geom_col(width = .56, position = position_stack(reverse = TRUE), color = 'white', linewidth = .3) +
				geom_text(aes(label = ifelse(is.finite(fraction) & fraction >= .065, n, '')),
					position = position_stack(vjust = .5, reverse = TRUE), size = 3.5, color = 'grey20'
				) +
				scale_fill_manual(values = pal, drop = FALSE) +
				scale_x_continuous(labels = scales::label_percent(), limits = c(0, 1), expand = expansion(mult = c(0, .01))) +
				labs(
					title = paste0(if (ep == 'incident') 'A. Incident' else 'B. Prevalent', ' evidence convergence'),
					subtitle = 'Source set: FDR < 0.05; opposite evidence tested in the same biomarkers', x = 'Fraction of source discoveries', y = NULL, fill = NULL
				) +
				theme_5c(11)
		})
		anchors <- if (layer == 'protein') c('PCSK9', 'LPA', 'GDF15', 'NTPROBNP', 'MMP12') else
			c('L_VLDL_TG.pct', 'ApoB', 'Glucose', 'GlycA', 'Total_TG')
		eligible <- z |>
			filter(available_incident | available_prevalent) |>
			arrange(q_mr, p_mr, exposure)
		observed <- eligible |>
			arrange(p_incident) |>
			slice_head(n = 3) |>
			pull(exposure)
		chosen <- head(unique(c(intersect(anchors, eligible$exposure), observed, eligible$exposure)), 16)
		z$shown_forest <- z$exposure %in% chosen
		forests <- lapply(c('mr', 'incident', 'prevalent'), function(k) {
			d <- z |>
				filter(shown_forest) |>
				transmute(
					exposure = factor(exposure, levels = rev(chosen)),
					b = .data[[paste0('b_', k)]], se = .data[[paste0('se_', k)]], available = .data[[paste0('available_', k)]],
					supported = .data[[paste0('supported_', k)]], b_mr, supported_mr
				)
			d <- d |> mutate(status = case_when(
				!available ~ 'Unavailable', !supported ~ 'No FDR support',
				k == 'mr' ~ 'MR FDR support', !supported_mr ~ 'Observed FDR support',
				sign(b) == sign(b_mr) ~ 'Concordant', TRUE ~ 'Discordant'
			), lo = b - 1.96 * se, hi = b + 1.96 * se)
			title <- switch(k,
				mr = paste0('C. ', if (layer == 'protein') 'cis' else 'local', ' Mendelian randomization'),
				incident = 'D. Incident disease · Cox',
				prevalent = 'E. Prevalent disease · logistic'
			)
			ggplot(d, aes(y = exposure)) +
				geom_hline(yintercept = seq_along(chosen), color = 'grey94', linewidth = .4) +
				geom_vline(xintercept = 0, color = 'grey65', linetype = 2, linewidth = .4) +
				geom_segment(data = filter(d, available), aes(x = lo, xend = hi, yend = exposure, color = status), linewidth = .8) +
				geom_point(data = filter(d, available), aes(x = b, color = status, shape = supported), size = 2.7) +
				geom_point(data = filter(d, !available), aes(x = 0), shape = 4, color = 'grey65', size = 2) +
				scale_y_discrete(drop = FALSE) +
				scale_color_manual(values = c(pal, 'MR FDR support' = '#866AA3', 'Observed FDR support' = '#B38D43')) +
				scale_shape_manual(values = c('TRUE' = 16, 'FALSE' = 1)) +
				guides(color = 'none', shape = 'none') +
				labs(
					title = title, subtitle = if (k == 'mr') 'Genetic effect; 95% CI' else 'Measured baseline omic level; 95% CI',
					x = if (k == 'incident') 'Log hazard ratio' else 'Log odds ratio', y = NULL
				) +
				theme_5c(11) +
				theme(panel.grid = element_blank())
		})
		caption <- paste(
			'Fixed anchors + top incident signals + primary MR ranking; the same rows are retained in all three forests.',
			'Filled points: within-analysis FDR < 0.05; open points: no FDR support; ×: unavailable.',
			'Purple: MR support; gold: observed support without MR support; coral/blue: concordant/discordant supported effects.',
			'Cox HRs and logistic/MR ORs are different estimands; separate horizontal scales are used.',
			'MR significance alone is not causal validation. Unlike the reference Figure 6, this descriptive adaptation does not impose a colocalization gate.'
		)
		list(plot = wrap_plots(c(bars, forests), design = 'AAABBB\nCCDDEE', heights = c(.55, 1.5)) +
			plot_annotation(caption = caption), data = z, counts = counts)
	}

	le8_emit_restored_c2 <- function(mr, assoc, layer, outdir) {
		rd <- le8_job_dir(outdir, 'c2_cause')
		old <- le8_plot_legacy_mr_overview(mr, assoc, layer)
		save_plot(old$plot, 'c2.Fig20.observational_mr_overview.png', 20, 15, outdir = outdir)
		save_plot(le8_plot_legacy_qtl_variance(mr, layer), 'c2.Fig21.qtl_variance_ranked.png', 19, 12, outdir = outdir)
		prefix <- if (layer == 'protein') 'pwas' else 'mwas'
		files <- file.path(outdir, 'c1_correlate', paste0(prefix, c('_incident_adj2.csv', '_prevalent_adj2.csv')))
		if (all(file.exists(files))) {
			z <- le8_plot_mr_incident_prevalent(mr, as_tibble(data.table::fread(files[1])), as_tibble(data.table::fread(files[2])), layer)
			save_plot(z$plot, 'c2.Fig22.mr_incident_prevalent.png', 22, 14, outdir = outdir)
			write_raw_csv(z$data, 'c2.mr_incident_prevalent.csv', rd)
			write_raw_csv(z$counts, 'c2.mr_incident_prevalent_support.csv', rd)
		} else save_plot(blank_plot('MR, incident and prevalent comparison', 'C1 incident/prevalent tables unavailable'),
			'c2.Fig22.mr_incident_prevalent.png', 22, 14,
			outdir = outdir
		)
	}
})
# Result records contain data frames: recursive modifyList corrupts their rows.
le8_result_update <- function(base, updates) {
	base[names(updates)] <- updates
	base
}
LE8_JOB <- "c2_cause"
C2_CODE_VERSION <- "2026-09-05.5c-audit-v1"
C2_FEATURE_SCOPE <- Sys.getenv("C2_FEATURE_SCOPE", unset = "all_qtl")
if (!C2_FEATURE_SCOPE %in% c("all_qtl", "observational_top")) stop("C2_FEATURE_SCOPE must be all_qtl or observational_top")
MAX_FEATURES <- as.integer(Sys.getenv("C2_MAX_FEATURES", unset = if (C2_FEATURE_SCOPE == "all_qtl") "0" else "500"))
if (!is.finite(MAX_FEATURES) || MAX_FEATURES < 0) stop("C2_MAX_FEATURES must be nonnegative; 0 means all mapped QTL traits")
C2_STAGE_VERSION <- paste0("2026-10-02.ungated-", C2_FEATURE_SCOPE, "-", MAX_FEATURES)
TOP_FOREST <- as.integer(Sys.getenv("C2_TOP_FOREST", unset = "32"))
N_DEFAULT <- as.numeric(Sys.getenv("C2_SUMSTAT_N", unset = "100000"))
normalize_c2_mode <- function(x, name) {
	x <- str_to_lower(str_trim(as.character(x)[1])) ; x <- case_when(
		x == "top" ~ "Top", x %in% c("all", "true", "1", "yes", "y") ~ "All",
		x %in% c("none", "false", "0", "no", "n") ~ "None", TRUE ~ NA_character_
	)
	if (is.na(x)) stop(name, " must be Top, All, or None.", call. = FALSE)
	x
}
# These names are intentionally case-sensitive.  RUN_MRLINK2 and
# RUN_DANDELION are retired so logs and batch scripts expose one contract.
RUN_MRlink2 <- normalize_c2_mode(Sys.getenv("RUN_MRlink2", unset = "Top"), "RUN_MRlink2")
RUN_Dandelion <- normalize_c2_mode(Sys.getenv("RUN_Dandelion", unset = "Top"), "RUN_Dandelion")
C2_FIXED_TOP <- unique(trimws(strsplit(Sys.getenv("C2_TOP_CANDIDATES",
	unset = "PCSK9,LPA,GDF15,NTPROBNP,MMP12"
), ",", fixed = TRUE)[[1]]))
C2_TOP_MAX <- suppressWarnings(as.integer(Sys.getenv("C2_TOP_MAX", unset = "50")))
if (!is.finite(C2_TOP_MAX) || C2_TOP_MAX < 5L) C2_TOP_MAX <- 50L
DANDELION_FDR <- as.numeric(Sys.getenv("C2_DANDELION_FDR", unset = "0.10"))
DANDELION_CIS_BP <- as.numeric(Sys.getenv("C2_DANDELION_CIS_BP", unset = "5000000"))
DANDELION_GWS <- as.numeric(Sys.getenv("C2_DANDELION_GWS", unset = "5e-8"))
DANDELION_LEAD_BP <- as.numeric(Sys.getenv("C2_DANDELION_LEAD_BP", unset = "5000000"))
DANDELION_MAX_SNPS <- as.integer(Sys.getenv("C2_DANDELION_MAX_SNPS", unset = "100"))
DANDELION_MAX_GENE2 <- as.integer(Sys.getenv("C2_DANDELION_MAX_GENE2", unset = "0")) # 0 = all available
DANDELION_ALLOW_MAGMA <- truthy(Sys.getenv("C2_DANDELION_ALLOW_MAGMA", unset = "TRUE"))
DANDELION_MAX_TARGET_FRACTION <- as.numeric(Sys.getenv("C2_DANDELION_MAX_TARGET_FRACTION", unset = "0.25"))

c2_method_scope_audit <- function(layer) {
	tibble(
		method = c("cis/local MR", "trans/distal MR", "MR-link-2", "DANDELION", "CIGMA"),
		role = c(
			"causal estimation", "pleiotropy-sensitive causal estimation", "regional LD-aware causal estimation",
			"distal regulatory driver prioritization", "cell-type-shared/specific eQTL variance decomposition"
		),
		eligible_with_current_inputs = c(TRUE, TRUE, layer == "protein", layer == "protein", FALSE),
		decision = c(
			"core C2", "core C2 with stricter interpretation", "optional C2 regional sensitivity",
			ifelse(layer == "protein", "primary only with valid disease gene/rare-variant evidence; MAGMA fallback is sensitivity", "not defined for metabolites"),
			"do not add to C2 causation core; optional cell-type annotation after external CIGMA analysis"
		),
		required_input = c(
			"omic QTL + disease GWAS", "omic QTL + disease GWAS", "dense cis/local QTL + disease GWAS + LD reference",
			"distal QTL plus valid disease-gene evidence",
			"population-scale single-cell RNA-seq h5ad, donor/cell metadata, pseudobulk and kinship"
		),
		reference = c(
			NA_character_, NA_character_, NA_character_, NA_character_,
			"Nature 2026, doi:10.1038/s41586-026-10577-6"
		)
	)
}
if (!is.finite(DANDELION_MAX_TARGET_FRACTION) || DANDELION_MAX_TARGET_FRACTION <= 0 || DANDELION_MAX_TARGET_FRACTION > 1)
	DANDELION_MAX_TARGET_FRACTION <- .25

# Publication typography for C2 only.
theme_5c <- function(base_size = 12) {
	theme_classic(base_size = base_size) + theme(
		plot.title = element_text(face = "bold", size = base_size * 1.14, hjust = 0),
		plot.subtitle = element_text(size = base_size * .94, color = "grey30"),
		axis.title = element_text(size = base_size * 1.08),
		axis.text = element_text(size = base_size, color = "black"),
		legend.title = element_text(), legend.text = element_text(),
		strip.background = element_blank(), strip.text = element_text(face = "bold"),
		panel.grid.major.y = element_line(color = "grey91", linewidth = .25), panel.grid.minor = element_blank(),
		plot.margin = margin(9, 13, 9, 13)
	)
}
forest_theme <- function(base_size = 10) theme_5c(base_size) + theme(panel.grid.major.y = element_blank())

select_c2_top_candidates <- function(layer, assoc, mr, universe) {
	primary <- if (layer == "protein") "cis" else "local"
	fixed <- C2_FIXED_TOP[C2_FIXED_TOP %in% universe]
	obs <- as_tibble(assoc) |>
		filter(term %in% universe, is.finite(p.value)) |>
		arrange(p.value) |>
		slice_head(n = 20) |>
		pull(term)
	gen <- as_tibble(mr) |>
		filter(exposure %in% universe, analysis == primary, is.finite(pval)) |>
		arrange(pval) |>
		slice_head(n = 25) |>
		pull(exposure)
	x <- head(unique(c(fixed, gen, obs)), C2_TOP_MAX)
	tibble(feature = x, top_rank = seq_along(x), top_reason = case_when(
		feature %in% fixed ~ "prespecified anchor",
		feature %in% gen ~ paste0("top ", primary, " MR"), TRUE ~ "top C1 association"
	))
}

le8_execute_mrlink2 <- function(layer, rawdir, jobs, ygfile, mode = RUN_MRlink2) {
	linkdir <- le8_cache_dir("mrlink2", basename(dirname(rawdir)), str_to_lower(mode)) ; dir.create(linkdir, recursive = TRUE, showWarnings = FALSE)
	complete_file <- file.path(linkdir, "mrlink2.complete")
	if (mode == "None") {
		message("C2/", layer, ": MR-link-2 not run (RUN_MRlink2=None)")
		return(invisible(0L))
	}
	if (cache_valid(complete_file)) {
		cache_message(paste0("MR-link-2/", layer), complete_file)
		return(invisible(0L))
	}
	if (!nrow(jobs)) return(invisible(0L))
	jobs_file <- file.path(linkdir, "c2.mrlink2.jobs.tsv") ; write_raw_tsv(jobs, basename(jobs_file), linkdir)
	ref_dir <- Sys.getenv("MRLINK2_REF_PFILE_DIR", unset = "") ; ref_bed <- Sys.getenv("MRLINK2_REF_BED", unset = "")
	ref_pop <- toupper(Sys.getenv("MRLINK2_REF_POP", unset = "EUR"))
	ref_id_dir <- Sys.getenv("MRLINK2_REF_ID_DIR", unset = "")
	ref_samples <- Sys.getenv("MRLINK2_REF_SAMPLES", unset = "")
	if (!nzchar(ref_id_dir) && nzchar(ref_dir)) ref_id_dir <- file.path(dirname(ref_dir), "id")
	if (!nzchar(ref_samples) && nzchar(ref_dir)) ref_samples <- file.path(dirname(ref_dir), "samples.txt")
	ref_keep <- if (!nzchar(ref_bed) && nzchar(ref_id_dir) && ref_pop != "ALL") file.path(ref_id_dir, paste0(ref_pop, ".id.2col")) else ""
	ref_bfile_dir <- if (!nzchar(ref_bed) && nzchar(ref_dir)) file.path(dirname(ref_dir), "bfile", ref_pop) else ""
	if (nzchar(ref_bed)) {
		ref_keep <- "" ; ref_samples <- ""
	}
	sh <- file.path(Sys.getenv("LE8_FDIR"), "c2.cause.sh")
	status <- system2("bash", c(sh, "mr-link2", "--jobs", jobs_file, "--cad-gwas", ygfile, "--outdir", linkdir))
	if (status != 0) warning("MR-link-2 runner returned status ", status) ; invisible(status)
}


# MR-link-2 needs dense marginal summary statistics from one biologically
# defined cis/local region.  The cleaning pipeline's *.cis.gz file is therefore
# authoritative when it exists; its observed bounds also preserve the exact
# flank used by gwas_post.sh instead of silently replacing it with a 20-Mb
# lead-SNP window.
mrlink2_file_bounds <- function(file, preferred_chr = NA_character_) {
	empty <- tibble(chr = character(), start = numeric(), end = numeric(), n_region_variants = integer())
	if (length(file) != 1L || is.na(file) || !file.exists(file) || file.size(file) <= 0) return(empty)
	con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
	hdr <- tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) "")
	try(close(con), silent = TRUE)
	if (!length(hdr) || !nzchar(hdr)) return(empty)
	nms <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
	cchr <- .pick_col(nms, c("CHR", "CHROM", "CHROMOSOME")) ; cpos <- .pick_col(nms, c("POS", "BP", "POSITION", "BASE_PAIR_LOCATION"))
	if (is.na(cchr) || is.na(cpos)) return(empty)
	d <- tryCatch(data.table::fread(file, select = unique(c(cchr, cpos)), showProgress = FALSE, check.names = FALSE), error = function(e) NULL)
	if (is.null(d) || !nrow(d)) return(empty)
	z <- tibble(
		chr = str_remove(as.character(d[[cchr]]), regex("^chr", ignore_case = TRUE)),
		pos = suppressWarnings(as.numeric(d[[cpos]]))
	) |>
		filter(!is.na(chr), nzchar(chr), is.finite(pos), pos >= 1)
	if (!nrow(z)) return(empty)
	pref <- str_remove(as.character(preferred_chr)[1], regex("^chr", ignore_case = TRUE))
	if (!is.na(pref) && nzchar(pref) && any(z$chr == pref)) z <- z |> filter(chr == pref)
	z |>
		group_by(chr) |>
		summarise(start = floor(min(pos)), end = ceiling(max(pos)), n_region_variants = n(), .groups = "drop") |>
		arrange(desc(n_region_variants), start) |>
		slice(1)
}

build_mrlink2_jobs <- function(iv_list, annotation, layer, ygfile) {
	primary <- if (layer == "protein") "cis" else "local"
	fallback_bp <- as.numeric(Sys.getenv(if (layer == "protein") "C2_CIS_WINDOW_BP" else "C2_LOCAL_WINDOW_BP", unset = "1000000"))
	jobs <- imap_dfr(iv_list, function(obj, feature) {
		iv <- obj$instruments
		if (!nrow(iv)) return(tibble())
		lead_pool <- iv |> filter(analysis == primary, is.finite(POS), !is.na(CHR), is.finite(P))
		# MR-link-2 is a cis/local method.  Do not relabel a trans/distal-only
		# exposure as cis merely to force a job through the runner.
		if (!nrow(lead_pool)) return(tibble())
		lead <- lead_pool |>
			arrange(P) |>
			slice(1)
		cis_file <- obj$files$cis
		use_cis <- length(cis_file) == 1L && !is.na(cis_file) && file.exists(cis_file) && file.size(cis_file) > 0
		bounds <- if (use_cis) mrlink2_file_bounds(cis_file, lead$CHR) else tibble()
		if (nrow(bounds)) {
			qf <- cis_file ; chr <- bounds$chr[[1]] ; start <- bounds$start[[1]] ; end <- bounds$end[[1]]
			source <- "cis.gz observed bounds" ; n_region <- bounds$n_region_variants[[1]]
		} else {
			qf <- obj$files$full
			if (length(qf) != 1L || is.na(qf) || !file.exists(qf)) return(tibble())
			chr <- as.character(lead$CHR[[1]]) ; start <- max(1, lead$POS[[1]] - fallback_bp) ; end <- lead$POS[[1]] + fallback_bp
			source <- paste0("full summary; ", primary, " lead +/-", format(fallback_bp, scientific = FALSE), " bp fallback")
			n_region <- NA_integer_
			if (layer == "protein" && !is.null(annotation) && nrow(annotation)) {
				a <- annotation |>
					filter(.data$feature == .env$feature) |>
					slice(1)
				if (nrow(a) && !is.na(a$chr[[1]]) && is.finite(a$start[[1]]) && is.finite(a$end[[1]])) {
					chr <- as.character(a$chr[[1]]) ; start <- max(1, a$start[[1]] - fallback_bp) ; end <- a$end[[1]] + fallback_bp
					source <- paste0("full summary; gene +/-", format(fallback_bp, scientific = FALSE), " bp fallback")
				}
			}
		}
		tibble(
			omics = layer, trait = feature, exposure = qf, outcome = ygfile,
			region = paste0(chr, ":", floor(start), "-", ceiling(end)), region_source = source,
			region_variants = n_region, lead_snp = lead$SNP[[1]], lead_p = lead$P[[1]]
		)
	})
	if (!nrow(jobs) || !"exposure" %in% names(jobs)) return(tibble())
	jobs |> filter(!is.na(exposure), file.exists(exposure))
}

venn_membership_counts <- function(a, b, c) {
	u <- sort(unique(c(a, b, c))) ; z <- tibble(id = u, A = u %in% a, B = u %in% b, C = u %in% c)
	tibble(
		region = c("A", "B", "C", "AB", "AC", "BC", "ABC"),
		n = c(
			sum(z$A & !z$B & !z$C), sum(!z$A & z$B & !z$C), sum(!z$A & !z$B & z$C),
			sum(z$A & z$B & !z$C), sum(z$A & !z$B & z$C), sum(!z$A & z$B & z$C), sum(z$A & z$B & z$C)
		)
	)
}

plot_venn_evidence <- function(a, b, c, names3, title, highlight = "A") {
	# Circle identity is fixed across panels: local/cis = left, distal/trans =
	# top, observational = right.  Only the panel's reference set changes fill.
	theta <- seq(0, 2 * pi, length.out = 240)
	cent <- tibble(set = c("A", "B", "C"), cx = c( - .42, 0, .42), cy = c( - .15, .33, - .15))
	circ <- cent |>
		group_by(set) |>
		group_modify( ~ tibble(x = .x$cx + .66 * cos(theta), y = .x$cy + .66 * sin(theta)))
	pos <- tibble(
		region = c("A", "B", "C", "AB", "AC", "BC", "ABC"),
		x = c( - .72, 0, .72, - .25, 0, .25, 0), y = c( - .30, .70, - .30, .19, - .40, .19, .02)
	) |>
		left_join(venn_membership_counts(a, b, c), by = "region")
	lab <- tibble(
		x = c( - .88, 0, .88), y = c( - .94, 1.08, - .94),
		label = paste0(names3, "\n(n=", c(length(unique(a)), length(unique(b)), length(unique(c))), ")")
	)
	ggplot(circ, aes(x, y, group = set, color = set)) +
		geom_polygon(
			data = circ |> filter(set == highlight), aes(x, y, group = set), inherit.aes = FALSE,
			fill = "#FFD54F", color = NA, alpha = .88
		) +
		geom_path(linewidth = 1.25) +
		geom_text(data = pos, aes(x, y, label = n), inherit.aes = FALSE, fontface = "bold", size = 4) +
		geom_text(data = lab, aes(x, y, label = label), inherit.aes = FALSE, fontface = "bold", size = 3) +
		scale_color_manual(values = c(A = "#D95F02", B = "#1B9E77", C = "#4C78A8"), guide = "none") +
		coord_equal(xlim = c( - 1.30, 1.30), ylim = c( - 1.10, 1.20), clip = "off") +
		theme_void() +
		labs(title = title, subtitle = "Yellow circle = reference set") +
		theme(
			plot.title = element_text(face = "bold", size = 11, hjust = .5),
			plot.subtitle = element_text(face = "bold", size = 8.5, hjust = .5, color = "#8A5A00"),
			plot.margin = margin(5, 16, 5, 16)
		)
}

plot_c2_fig1 <- function(mr, assoc, layer) {
	local_name <- if (layer == "protein") "cis" else "local" ; distal_name <- if (layer == "protein") "trans" else "distal"
	mr0 <- mr |> filter(is.finite(pval)) ; obs_sig <- assoc |>
		filter(is.finite(p.value), p.value < .05 / max(1, nrow(assoc))) |>
		pull(term)
	local <- mr0 |> filter(analysis == local_name) ; distal <- mr0 |> filter(analysis == distal_name)
	local_sig <- local |>
		filter(FDR_analysis < .05) |>
		pull(exposure) ; distal_sig <- distal |>
		filter(FDR_analysis < .05) |>
		pull(exposure)
	topn <- min(TOP_FOREST, max(10L, floor(sqrt(max(1, nrow(mr0))) * 2)))
	top_local <- local |>
		slice_min(pval, n = topn, with_ties = FALSE) |>
		pull(exposure)
	top_distal <- distal |>
		slice_min(pval, n = topn, with_ties = FALSE) |>
		pull(exposure)
	top_obs <- assoc |>
		filter(is.finite(p.value)) |>
		slice_min(p.value, n = topn, with_ties = FALSE) |>
		pull(term)
	va <- plot_venn_evidence(
		top_local, distal_sig, obs_sig,
		c(paste0("Top ", local_name), paste0(distal_name, " FDR"), "Observational Bonf."), paste0("a. Anchored on ", local_name, " MR"), "A"
	)
	vb <- plot_venn_evidence(
		local_sig, top_distal, obs_sig,
		c(paste0(local_name, " FDR"), paste0("Top ", distal_name), "Observational Bonf."), paste0("b. Anchored on ", distal_name, " MR"), "B"
	)
	vc <- plot_venn_evidence(
		local_sig, distal_sig, top_obs,
		c(paste0(local_name, " FDR"), paste0(distal_name, " FDR"), "Top observational"), "c. Anchored on observational evidence", "C"
	)
	top_ids <- unique(c(head(top_local, 12), head(top_distal, 12), head(top_obs, 12)))
	ev <- bind_rows(
		assoc |> filter(term %in% top_ids) |> transmute(exposure = term, evidence = "Observational", b = beta, se = std.error, p = p.value),
		mr0 |> filter(exposure %in% top_ids, analysis %in% c(local_name, distal_name)) |> transmute(exposure, evidence = paste0("MR (", analysis, ")"), b, se, p = pval)
	) |>
		complete(exposure = top_ids, evidence = c("Observational", paste0("MR (", local_name, ")"), paste0("MR (", distal_name, ")"))) |>
		mutate(
			available = is.finite(b) & is.finite(se), lo = b - 1.96 * se, hi = b + 1.96 * se, b_plot = ifelse(available, b, 0),
			exposure = factor(exposure, levels = rev(top_ids))
		)
	pforest <- ggplot(ev, aes(b_plot, exposure)) +
		geom_vline(xintercept = 0, color = "grey60") +
		geom_errorbarh(data = ev |> filter(available), aes(xmin = lo, xmax = hi), height = 0, color = "grey55") +
		geom_point(data = ev |> filter(available), aes(color = p < .05), size = 2) +
		geom_point(data = ev |> filter(!available), shape = 4, size = 2.2, color = "grey72") +
		facet_grid(. ~ evidence, scales = "free_x") +
		scale_color_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#4C78A8"), guide = "none") +
		labs(title = "d. Effect estimates; × denotes no valid instrument result", x = "Effect per 1-SD higher omic trait", y = NULL) +
		forest_theme(8)
	best <- mr0 |>
		group_by(exposure) |>
		slice_min(pval, n = 1, with_ties = FALSE) |>
		ungroup() |>
		left_join(assoc |> select(exposure = term, obs_beta = beta, obs_se = std.error, obs_p = p.value), by = "exposure") |>
		mutate(concordance = case_when(
			FDR_all < .05 & sign(b) == sign(obs_beta) ~ "FDR-significant, concordant",
			FDR_all < .05 & sign(b) != sign(obs_beta) ~ "FDR-significant, discordant", TRUE ~ "No FDR MR support"
		))
	counts <- bind_rows(
		tibble(anchor = local_name, venn_membership_counts(top_local, distal_sig, obs_sig)),
		tibble(anchor = distal_name, venn_membership_counts(local_sig, top_distal, obs_sig)),
		tibble(anchor = "observational", venn_membership_counts(local_sig, distal_sig, top_obs))
	)
	composed <- (va | vb | vc) / pforest + plot_layout(heights = c(.82, 1.18)) +
		plot_annotation(
			caption = paste0(
				toupper(substr(local_name, 1, 1)), substr(local_name, 2, nchar(local_name)),
				" is always left; ", distal_name, " is always top; observational is always right. Yellow = panel reference set."
			),
			theme = theme(plot.caption = element_text(size = 8, color = "grey40", hjust = 1))
		)
	list(plot = composed, bars = counts, forest = ev, best = best)
}

plot_c2_fig2 <- function(mr, assoc, layer) {
	le8_plot_cis_trans_comparison(mr, assoc, layer)
}

plot_qtl_variance <- function(mr, layer) {
	le8_plot_instrument_comparison(mr, layer)
}

# C2 has one output contract.  Standard figure names are overwritten whenever
# figures are regenerated from a valid numerical cache.
save_c2_plot <- function(p, file, width, height, outdir) {
	save_plot(p, file, width, height, outdir = outdir)
}

plot_c2_fig4 <- function(mr, layer) {
	for (nm in c("egger_intercept_p", "steiger_support", "steiger_support_fraction", "steiger_n", "median_F", "r2_median")) if (!nm %in% names(mr)) mr[[nm]] <- NA
	d <- mr |>
		filter(is.finite(pval), is.finite(median_F), is.finite(r2_median)) |>
		mutate(
			class = analysis, significant = FDR_analysis < .05, label = ifelse(min_rank(pval) <= 8, exposure, NA_character_),
			q_score = stable_neglog10_p(Q_p), egger_score = stable_neglog10_p(egger_intercept_p),
			q_plot = compress_extreme_tail(q_score, 10, 4), egger_plot = compress_extreme_tail(egger_score, 10, 4)
		)
	if (!nrow(d)) return(blank_plot("MR sensitivity and instrument architecture", "No complete instrument metrics"))
	pa <- ggplot(d, aes(median_F, 100 * r2_median, color = class, size = pmin(n_IV, 500))) +
		geom_vline(xintercept = 10, linetype = 2, color = "grey55") +
		geom_point(alpha = .78) +
		ggrepel::geom_text_repel(aes(label = label), size = 2.35, seed = 31, max.overlaps = 10, na.rm = TRUE) +
		scale_x_log10() +
		labs(
			title = "a. Instrument strength and typical variance", subtitle = "Median per-IV metrics avoid the many-IV product-to-one artefact",
			x = "Median F statistic (log scale)", y = "Median per-IV partial R² (%)", color = NULL, size = "IV count"
		) +
		theme_5c(9)
	ds <- d |> filter(is.finite(q_plot), is.finite(egger_plot))
	pb <- if (!nrow(ds)) blank_plot("b. Heterogeneity and directional pleiotropy", "Egger diagnostics require at least three valid IVs") else {
		sx <- compressed_tail_scale(ds$q_score, 10, 4) ; sy <- compressed_tail_scale(ds$egger_score, 10, 4)
		ggplot(ds, aes(q_plot, egger_plot, color = class)) +
			geom_vline(xintercept =  - log10(.05), linetype = 2, color = "grey60") +
			geom_hline(yintercept =  - log10(.05), linetype = 2, color = "grey60") +
			geom_point(aes(size = pmin(n_IV, 500)), alpha = .78) +
			scale_x_continuous(breaks = sx$breaks, labels = sx$labels) +
			scale_y_continuous(breaks = sy$breaks, labels = sy$labels) +
			labs(
				title = "b. Heterogeneity and directional pleiotropy", subtitle = "Axes show -log10(P); values above 10 are monotonically compressed",
				x = expression( - log[10](P[Cochran ~ Q])), y = expression( - log[10](P[Egger ~ intercept])), color = NULL, size = "IV count"
			) +
			theme_5c(9)
	}
	dz <- d |> filter(is.finite(steiger_support_fraction))
	pc <- if (!nrow(dz)) blank_plot("c. Per-IV variance-direction heuristic", "No per-IV directionality comparison was estimable") else
		ggplot(dz, aes(class, steiger_support_fraction, color = class)) +
			geom_hline(yintercept = .5, linetype = 2, color = "grey55") +
			geom_violin(aes(fill = class), alpha = .12, color = NA, trim = FALSE) +
			geom_boxplot(width = .18, outlier.shape = NA, alpha = .55) +
			geom_jitter(aes(size = pmin(steiger_n, 500)), width = .08, alpha = .48) +
			scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
			labs(
				title = "c. Per-IV variance-direction heuristic", subtitle = "Fraction of harmonized IVs with R²(exposure) > R²(outcome); >50% is heuristic only; not a liability-scale direction test",
				x = NULL, y = "Supporting IV fraction", color = NULL, fill = NULL, size = "IV count"
			) +
			theme_5c(9) +
			theme(legend.position = "none")
	(pa | pb) / pc + plot_layout(heights = c(1, .72))
}

read_mrlink2_results <- function(rawdir, mode = RUN_MRlink2) {
	linkdir <- le8_cache_dir("mrlink2", basename(dirname(rawdir)), str_to_lower(mode))
	f <- file.path(linkdir, "mrlink2.all.tsv") ; s <- file.path(linkdir, "mrlink2.status.tsv")
	ans <- list(results = tibble(), status = tibble(), mode = mode)
	if (file.exists(f) && file.size(f) > 0) ans$results <- tryCatch(as_tibble(data.table::fread(f, showProgress = FALSE, check.names = FALSE)), error = function(e) tibble())
	if (file.exists(s) && file.size(s) > 0) ans$status <- tryCatch(as_tibble(data.table::fread(s, showProgress = FALSE, check.names = FALSE)), error = function(e) tibble())
	ans
}

plot_mrlink2_results <- function(x, outdir) {
	d <- x$results ; st <- x$status
	if (nrow(d)) {
		ac <- pick_col_local(names(d), c("^ALPHA$")) ; sc <- pick_col_local(names(d), c("^SE\\(ALPHA\\)$", "^SE_ALPHA$")) ; pc <- pick_col_local(names(d), c("^P\\(ALPHA\\)$", "^P_ALPHA$"))
		yc <- pick_col_local(names(d), c("^SIGMA_Y$", "^SIGMAY$")) ; ypc <- pick_col_local(names(d), c("^P\\(SIGMA_Y\\)$", "^P_SIGMA_Y$", "^PSIGMAY$"))
		nc <- pick_col_local(names(d), c("^M_SNPS_OVERLAP$", "^MOVERLAP$", "^N_SNPS$")) ; rc <- pick_col_local(names(d), c("^REGION$"))
		if (all(!is.na(c(ac, sc, pc)))) {
			num_or_na <- function(nm) if (is.na(nm)) rep(NA_real_, nrow(d)) else suppressWarnings(as.numeric(d[[nm]]))
			char_or_na <- function(nm) if (is.na(nm)) rep(NA_character_, nrow(d)) else as.character(d[[nm]])
			trait_col <- pick_col_local(names(d), c("^TRAIT$", "^EXPOSURE$"))
			z_all <- tibble(
				trait = char_or_na(trait_col), region = char_or_na(rc), alpha = num_or_na(ac), se = num_or_na(sc), p = num_or_na(pc),
				sigma_y = num_or_na(yc), p_sigma_y = num_or_na(ypc), n_overlap = num_or_na(nc)
			) |>
				filter(!is.na(trait), trait != "", is.finite(alpha), is.finite(se), se > 0, is.finite(p)) |>
				group_by(trait) |>
				slice_min(p, n = 1, with_ties = FALSE) |>
				ungroup() |>
				arrange(p)
			if (!nrow(z_all)) {
				pa <- blank_plot("a. MR-link-2 regional causal effects", "Mapped columns contained no finite estimates")
				pb <- blank_plot("b. MR-link-2 causal evidence", "Mapped columns contained no finite estimates")
				pc <- blank_plot("c. Regional pleiotropy audit", "Mapped columns contained no finite estimates")
			} else {
				z <- z_all |>
					slice_head(n = 32) |>
					mutate(lo = alpha - 1.96 * se, hi = alpha + 1.96 * se)
				robust_lim <- as.numeric(quantile(abs(c(z$lo, z$hi)), .98, na.rm = TRUE, names = FALSE))
				if (!is.finite(robust_lim) || robust_lim <= 0) robust_lim <- max(abs(c(z$lo, z$hi)), na.rm = TRUE)
				z <- z |> mutate(
					lo_plot = pmax(lo, - robust_lim), hi_plot = pmin(hi, robust_lim), clipped = lo <  - robust_lim | hi > robust_lim,
					trait = factor(trait, levels = rev(trait))
				)
				pa <- ggplot(z, aes(alpha, trait)) +
					geom_vline(xintercept = 0, color = "grey60") +
					geom_errorbarh(aes(xmin = lo_plot, xmax = hi_plot), height = .08) +
					geom_point(aes(color = p < .05 / max(1, nrow(z))), size = 2) +
					geom_point(data = z |> filter(clipped), aes(x = ifelse(alpha >= 0, robust_lim, - robust_lim), y = trait), shape = 17, size = 2.1, color = "grey35", inherit.aes = FALSE) +
					scale_color_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#4C78A8"), guide = "none") +
					coord_cartesian(xlim = c( - robust_lim, robust_lim)) +
					labs(
						title = "a. MR-link-2 regional causal effects", subtitle = "LD-aware estimate per trait; triangles mark CIs clipped at the robust display limit",
						x = expression(alpha), y = NULL
					) +
					forest_theme(8)
				zb <- z_all |> mutate(score = stable_neglog10_p(p), score_plot = compress_extreme_tail(score, 10, 4), label = ifelse(min_rank(p) <= 14, trait, NA_character_))
				sy <- compressed_tail_scale(zb$score, 10, 4)
				pb <- ggplot(zb, aes(alpha, score_plot)) +
					geom_vline(xintercept = 0, color = "grey60") +
					geom_hline(yintercept =  - log10(.05 / max(1, nrow(zb))), linetype = 2, color = "grey55") +
					geom_point(color = "#6A3D9A", size = 2, alpha = .78) +
					ggrepel::geom_text_repel(aes(label = label), size = 2.45, seed = 41, max.overlaps = 15, na.rm = TRUE) +
					scale_y_continuous(breaks = sy$breaks, labels = sy$labels) +
					labs(
						title = "b. MR-link-2 causal evidence", subtitle = "-log10(P) above 10 is monotonically compressed",
						x = expression(alpha), y = expression( - log[10](P[alpha]))
					) +
					theme_5c(9)
				audit_txt <- if (nrow(st) && "status" %in% names(st)) paste(paste0(names(table(st$status)), "=", as.integer(table(st$status))), collapse = "; ") else "status file unavailable"
				zc <- z_all |>
					filter(is.finite(sigma_y)) |>
					mutate(
						label = ifelse(min_rank(p) <= 12, trait, NA_character_),
						pleiotropy = ifelse(is.finite(p_sigma_y) & p_sigma_y < .05, "P(sigma_y) < 0.05", "No sigma_y support"),
						n_overlap = ifelse(is.finite(n_overlap), pmax(1, n_overlap), 1)
					)
				pc <- if (!nrow(zc)) blank_plot("c. Regional pleiotropy audit", "sigma_y fields were not returned by this MR-link-2 build") else
					ggplot(zc, aes(alpha, sigma_y, color = pleiotropy, size = n_overlap)) +
						geom_hline(yintercept = 0, color = "grey70") +
						geom_vline(xintercept = 0, color = "grey70") +
						geom_point(alpha = .75) +
						ggrepel::geom_text_repel(aes(label = label), size = 2.4, seed = 42, max.overlaps = 12, na.rm = TRUE) +
						scale_color_manual(values = c(`P(sigma_y) < 0.05` = "#D95F02", `No sigma_y support` = "#4C78A8")) +
						scale_size_continuous(range = c(1.5, 5), name = "Overlapping SNPs") +
						labs(
							title = "c. Regional causal effect and residual pleiotropy", subtitle = paste0("Execution audit: ", audit_txt),
							x = expression(alpha), y = expression(sigma[y]), color = NULL
						) +
						theme_5c(9)
			}
		} else {
			pa <- blank_plot("a. MR-link-2 regional causal effects", "Output columns could not be mapped")
			pb <- blank_plot("b. MR-link-2 causal evidence", "Output columns could not be mapped")
			pc <- blank_plot("c. Regional pleiotropy audit", "Output columns could not be mapped")
		}
	} else {
		why <- if (identical(x$mode, "None")) "RUN_MRlink2=None; no regional jobs were requested" else
			paste0("RUN_MRlink2=", x$mode %||% "unknown", "; no completed regional estimate")
		pa <- blank_plot("MR-link-2 execution status", why)
		pb <- blank_plot("Candidate selection", why)
		pc <- blank_plot("Regional pleiotropy audit", why)
	}
	save_plot(pa / (pb | pc) + plot_layout(heights = c(1.2, .8)), "c2.Fig5.mrlink2.png", 14.5, 11, outdir = outdir)
}


# 🚩 DANDELION helpers
pick_col_local <- function(nms, patterns) {
	nl <- toupper(nms)
	for (p in patterns) {
		i <- grep(p, nl, perl = TRUE)[1] ; if (!is.na(i)) return(nms[[i]])
	}
	NA_character_
}

find_magma_gene_file <- function(ygfile) {
	explicit <- Sys.getenv("C2_MAGMA_GENE_FILE", unset = "") ; if (nzchar(explicit) && file.exists(explicit)) return(normalizePath(explicit, winslash = "/", mustWork = FALSE))
	trait <- sub("\\.gz$", "", basename(ygfile), ignore.case = TRUE)
	# gwas_post contract: <project>/common/<trait>/{gwas,magma}.
	magma_dir <- gwas_magma_dir_from_clean_file(ygfile)
	explicit_candidates <- file.path(magma_dir, paste0(trait, ".genes.out"))
	discovered <- if (!is.na(magma_dir) && dir.exists(magma_dir)) list.files(magma_dir, pattern = "(genes\\.out|gene\\.out)$", full.names = TRUE, ignore.case = TRUE) else character()
	fs <- unique(c(explicit_candidates, discovered)) ; fs <- fs[file.exists(fs) & !grepl("(gsa|sets)\\.out$", fs, ignore.case = TRUE)]
	if (length(fs)) return(normalizePath(fs[[which.max(file.info(fs)$mtime)]], winslash = "/", mustWork = FALSE))
	NA_character_
}


map_magma_gene_ids <- function(ids, target_symbols) {
	ids <- as.character(ids) ; target_symbols <- unique(as.character(target_symbols)) ; over <- mean(ids %in% target_symbols, na.rm = TRUE)
	if (is.finite(over) && over > .05) return(ids)
	mapfile <- Sys.getenv("C2_MAGMA_GENE_MAP", unset = "")
	if (nzchar(mapfile) && file.exists(mapfile)) {
		m <- data.table::fread(mapfile, showProgress = FALSE) ; if (ncol(m) >= 2) {
			mp <- setNames(as.character(m[[2]]), as.character(m[[1]])) ; z <- unname(mp[ids]) ; z[is.na(z)] <- ids[is.na(z)] ; return(z)
		}
	}
	if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
		z <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")) ; z <- as.character(z) ; z[is.na(z)] <- ids[is.na(z)] ; return(z)
	}
	ids
}


read_lead_file <- function(file, ygfile) {
	if (is.na(file) || !file.exists(file)) return(tibble())
	d <- if (grepl("jma\\.cojo$", file, ignore.case = TRUE)) match_GRCH_table(file, ygfile) else data.table::fread(file, showProgress = FALSE, check.names = FALSE, fill = TRUE) ; nms <- names(d)
	sc <- pick_col_local(nms, c("^SNP$", "^RSID$", "^ID$")) ; cc <- pick_col_local(nms, c("^CHR$", "^CHROM$")) ; bc <- pick_col_local(nms, c("^BP$", "^POS$", "^POSITION$")) ; pc <- pick_col_local(nms, c("^P$", "^PVAL$", "^P_VALUE$", "^P_BOLT_LMM$"))
	if (is.na(sc)) return(tibble())
	z <- tibble(SNP = as.character(d[[sc]]), CHR = if (!is.na(cc)) as.character(d[[cc]]) else NA_character_, POS = if (!is.na(bc)) as.numeric(d[[bc]]) else NA_real_, P = if (!is.na(pc)) as.numeric(d[[pc]]) else NA_real_)
	miss <- z$SNP[!is.finite(z$P) | !is.finite(z$POS) | is.na(z$CHR)] ; if (length(miss)) {
		full <- read_sumstat_snps(ygfile, miss, N_DEFAULT) ; if (nrow(full)) z <- z |>
			left_join(full |> select(SNP, CHR2 = CHR, POS2 = POS, P2 = P), by = "SNP") |>
			mutate(CHR = coalesce(CHR, CHR2), POS = coalesce(POS, POS2), P = coalesce(P, P2)) |>
			select( - ends_with("2"))
	}
	z |>
		filter(!is.na(SNP), SNP != "", !is.na(CHR), CHR != "", is.finite(POS)) |>
		distinct(SNP, .keep_all = TRUE)
}

stream_gws_hits <- function(ygfile, threshold = DANDELION_GWS) {
	hdr <- tryCatch(readLines(if (grepl("\\.gz$", ygfile)) gzfile(ygfile, "rt") else file(ygfile, "rt"), n = 1), error = function(e) "")
	if (!nzchar(hdr)) return(tibble()) ; sep <- if (grepl("\\t", hdr)) "\t" else " " ; nms <- strsplit(trimws(hdr), "[[:space:]]+")[[1]]
	sc <- pick_col_local(nms, c("^SNP$", "^RSID$", "^ID$", "^VARIANT_ID$")) ; cc <- pick_col_local(nms, c("^CHR$", "^CHROM$")) ; bc <- pick_col_local(nms, c("^BP$", "^POS$", "^POSITION$")) ; pc <- pick_col_local(nms, c("^P$", "^PVAL$", "^P_VALUE$", "^P_BOLT_LMM$"))
	idx <- match(c(sc, cc, bc, pc), nms) ; if (any(!is.finite(idx))) return(tibble())
	if (.Platform$OS.type != "windows" && nzchar(Sys.which("awk"))) {
		dec <- if (grepl("\\.gz$", ygfile)) paste("gzip -cd", shQuote(ygfile)) else paste("cat", shQuote(ygfile))
		awk <- sprintf("awk 'BEGIN{OFS=\"\\t\"} NR==1{next} $%d<=%g {print $%d,$%d,$%d,$%d}'", idx[4], threshold, idx[1], idx[2], idx[3], idx[4])
		z <- tryCatch(data.table::fread(cmd = paste(dec, "|", awk), col.names = c("SNP", "CHR", "POS", "P"), showProgress = FALSE), error = function(e) NULL)
		if (!is.null(z)) return(as_tibble(z) |> mutate(SNP = as.character(SNP), CHR = as.character(CHR), POS = as.numeric(POS), P = as.numeric(P)) |>
			filter(!is.na(SNP), SNP != "", !is.na(CHR), CHR != "", is.finite(POS), is.finite(P), P <= threshold))
	}
	d <- data.table::fread(ygfile, select = unique(idx), showProgress = FALSE, check.names = FALSE)
	tibble(SNP = as.character(d[[sc]]), CHR = as.character(d[[cc]]), POS = as.numeric(d[[bc]]), P = as.numeric(d[[pc]])) |>
		filter(!is.na(SNP), SNP != "", !is.na(CHR), CHR != "", is.finite(POS), is.finite(P), P <= threshold)
}

distance_prune_leads <- function(z, window = DANDELION_LEAD_BP, max_n = DANDELION_MAX_SNPS) {
	if (!nrow(z)) return(z)
	window <- suppressWarnings(as.numeric(window)) ; if (length(window) != 1L || !is.finite(window) || window < 0) stop("DANDELION lead-pruning window must be a finite non-negative number.", call. = FALSE)
	max_n <- suppressWarnings(as.integer(max_n)) ; if (length(max_n) != 1L || is.na(max_n) || max_n < 0) stop("DANDELION_MAX_SNPS must be a non-negative integer.", call. = FALSE)
	z <- z |>
		mutate(CHR = str_remove(as.character(CHR), "^chr"), POS = suppressWarnings(as.numeric(POS)), P = suppressWarnings(as.numeric(P))) |>
		filter(!is.na(SNP), SNP != "", !is.na(CHR), CHR != "", is.finite(POS), is.finite(P)) |>
		arrange(P)
	if (!nrow(z)) return(z)
	keep <- logical(nrow(z)) ; picked <- list()
	for (i in seq_len(nrow(z))) {
		chr <- z$CHR[[i]] ; pos <- z$POS[[i]] ; old <- picked[[chr]] %||% numeric()
		if (!length(old) || all(abs(pos - old) > window)) {
			keep[[i]] <- TRUE ; picked[[chr]] <- c(old, pos)
		}
	}
	z <- z[keep, , drop = FALSE] ; if (max_n > 0) z <- head(z, max_n) ; z
}

find_disease_lead_candidates <- function(ygfile) {
	explicit <- Sys.getenv("C2_DANDELION_SNP_FILE", unset = "") ; cand <- character()
	if (nzchar(explicit)) cand <- c(cand, explicit)
	d <- dirname(ygfile) ; cand <- c(cand, list.files(d, pattern = "(jma\\.cojo$|\\.clumped$|lead.*\\.(txt|tsv|csv)$|indep.*\\.(txt|tsv|csv)$)", full.names = TRUE, ignore.case = TRUE))
	unique(cand[file.exists(cand)])
}

find_disease_leads <- function(ygfile) {
	cand <- find_disease_lead_candidates(ygfile)
	if (length(cand)) {
		for (f in cand) {
			z <- read_lead_file(f, ygfile) ; if (nrow(z) >= 2) {
				z <- z |> arrange(P) ; if (DANDELION_MAX_SNPS > 0) z <- head(z, DANDELION_MAX_SNPS) ; attr(z, "source") <- f ; return(z)
			}
		}
	}
	warning("DANDELION: no independent/lead-SNP file found beside outcome GWAS; using distance-pruned genome-wide significant SNPs. For the primary analysis, provide an LD-independent file via C2_DANDELION_SNP_FILE.", call. = FALSE)
	z <- distance_prune_leads(stream_gws_hits(ygfile, DANDELION_GWS)) ; attr(z, "source") <- "distance-pruned GWS from outcome GWAS" ; z
}

build_trans_p_matrix <- function(gene2, lead, base) {
	snps <- lead$SNP
	rows <- parallel_map(gene2, function(g) {
		fs <- find_qtl_files(g, base, "protein") ; f <- fs$full
		if (is.na(f) || !file.exists(f)) return(list(p = setNames(rep(NA_real_, length(snps)), snps), file = NA_character_, n = 0L, error = "full QTL missing"))
		err <- NA_character_ ; z <- tryCatch(read_sumstat_snps(f, snps, N_DEFAULT), error = function(e) {
			err <<- conditionMessage(e) ; tibble()
		})
		p <- rep(NA_real_, length(snps)) ; names(p) <- snps
		if (nrow(z)) p[z$SNP] <- z$P
		list(p = p, file = f, n = sum(is.finite(p)), error = err)
	})
	M <- do.call(rbind, lapply(rows, `[[`, "p")) ; rownames(M) <- gene2 ; colnames(M) <- snps ; storage.mode(M) <- "double"
	attr(M, "qtl_audit") <- tibble(
		gene2 = gene2, full_qtl = vapply(rows, function(x) x$file %||% NA_character_, character(1)),
		n_disease_snps_found = vapply(rows, function(x) as.integer(x$n %||% 0L), integer(1)),
		read_error = vapply(rows, function(x) x$error %||% NA_character_, character(1))
	)
	M
}

# Fall back to assayed-protein annotations when no gene annotation is supplied.
read_dandelion_gene_ref <- function(assayed_ref) {
	f <- Sys.getenv("C2_DANDELION_GENE_ANNOTATION", unset = "")
	if (!nzchar(f) || !file.exists(f)) {
		warning("DANDELION: C2_DANDELION_GENE_ANNOTATION was not supplied; SNP-to-cis-gene labels will be inferred only from assayed protein genes. Supply a genome-wide gene annotation table for the primary disease-driver analysis.", call. = FALSE)
		return(assayed_ref)
	}
	d <- data.table::fread(f, showProgress = FALSE, check.names = FALSE) ; nms <- names(d)
	gc <- pick_col_local(nms, c("^GENE_NAME$", "^GENE$", "^SYMBOL$", "^GENESYMBOL$", "^GENE_SYMBOL$"))
	tc <- pick_col_local(nms, c("^TYPE$", "^GENE_TYPE$", "^BIOTYPE$"))
	cc <- pick_col_local(nms, c("^CHROMOSOME$", "^CHR$", "^CHROM$"))
	sc <- pick_col_local(nms, c("^START$", "^GENE_START$", "^START_POS$"))
	ec <- pick_col_local(nms, c("^END$", "^GENE_END$", "^END_POS$"))
	if (any(is.na(c(gc, cc, sc, ec)))) stop("C2_DANDELION_GENE_ANNOTATION needs gene/symbol, chromosome, start and end columns: ", f, call. = FALSE)
	z <- tibble(
		gene_name = as.character(d[[gc]]), type = if (!is.na(tc)) as.character(d[[tc]]) else "protein_coding",
		Chromosome = paste0("chr", str_remove(as.character(d[[cc]]), regex("^chr", ignore_case = TRUE))),
		start = suppressWarnings(as.numeric(d[[sc]])), end = suppressWarnings(as.numeric(d[[ec]]))
	) |>
		filter(!is.na(gene_name), gene_name != "", is.finite(start), is.finite(end)) |>
		mutate(type = coalesce(na_if(type, ""), "protein_coding"), type = case_when(str_detect(str_to_lower(type), "protein") ~ "protein_coding", str_detect(str_to_lower(type), "linc|lnc") ~ "lincRNA", TRUE ~ type)) |>
		distinct(gene_name, .keep_all = TRUE)
	# Always retain exact coordinates for the assayed proteins used as gene2.
	bind_rows(assayed_ref, z |> filter(!gene_name %in% assayed_ref$gene_name)) |> distinct(gene_name, .keep_all = TRUE)
}


safe_min_finite <- function(x) {
	x <- x[is.finite(x)] ; if (length(x)) min(x) else NA_real_
}


# 🚩 Main C2
le8_count_available_instruments <- function(iv_list, ygwas, layer) {
	classes <- if (layer == "protein") c("cis", "trans") else c("local", "distal")
	imap_dfr(iv_list, function(obj, feature) {
		iv <- obj$instruments %||% tibble() ; fs <- obj$files %||% list(joint = NA_character_, full = NA_character_)
		map_dfr(classes, function(cc) {
			z <- if (nrow(iv)) iv |> filter(analysis == cc) else tibble() ; ov <- if (nrow(z)) sum(z$SNP %in% ygwas$SNP) else 0L
			hz <- if (nrow(z) && ov > 0) tryCatch(nrow(harmonize_sumstats(z, ygwas)), error = function(e) 0L) else 0L
			reason <- case_when(
				nrow(z) > 0 && hz > 0 ~ "Available",
				is.na(fs$joint) || !file.exists(fs$joint) ~ "No independent clump/COJO instrument file",
				!nrow(iv) ~ "No allele-complete independent instrument",
				!nrow(z) ~ paste0("No ", cc, " instrument after genomic classification"),
				ov == 0 ~ "No SNP overlap with outcome GWAS", hz == 0 ~ "Alleles could not be harmonized", TRUE ~ "Unavailable"
			)
			tibble(
				exposure = feature, analysis = cc, joint_file = fs$joint %||% NA_character_, full_file = fs$full %||% NA_character_,
				n_independent = nrow(iv), n_class = nrow(z), n_outcome_overlap = ov, n_harmonized = hz, status = reason
			)
		})
	})
}


restore_c2_figures <- function(cached, layer, outdir, rawdir = NULL) {
	mr <- as_tibble(cached$MR %||% tibble()) ; assoc <- as_tibble(cached$observational %||% tibble())
	if (!nrow(mr) || !nrow(assoc)) return(invisible(FALSE))
	fig1 <- plot_c2_fig1(mr, assoc, layer) ; save_c2_plot(fig1$plot, "c2.Fig1.prots.top.png", 17, 13.5, outdir = outdir)
	save_c2_plot(plot_qtl_variance(mr, layer), "c2.Fig2.pQTL_R2.png", 15.5, 10.5, outdir = outdir)
	fig3 <- plot_c2_fig2(mr, assoc, layer) ; save_plot(fig3$plot, "c2.Fig3.effect_concordance.png", 16, 12.5, outdir = outdir)
	write_raw_csv(fig3$wide, "c2.cis_trans_comparison.csv", le8_job_dir(outdir, "c2_cause"))
	le8_emit_restored_c2(mr, assoc, layer, outdir)
	save_plot(plot_c2_fig4(mr, layer), "c2.Fig4.sensitivity_architecture.png", 15.5, 11.5, outdir = outdir)
	if (!is.null(rawdir)) plot_mrlink2_results(read_mrlink2_results(rawdir), outdir)
	dan <- cached$DANDELION %||% list()
	if (layer == "protein") {
		plot_dandelion_results(dan$targets %||% tibble(), dan$pairs %||% tibble(), dan$gene_pairs %||% tibble(), outdir)
		plot_dandelion_landscape(dan$evidence_plot %||% tibble(), dan$targets %||% tibble(), outdir, dan$gene_evidence_type %||% "gene-level disease association")
		integration <- plot_dandelion_mr_integration(dan, mr, assoc, outdir)
		if (!is.null(rawdir) && nrow(dan$gene_pairs %||% tibble())) {
			p_gene <- tryCatch(read_gene_level_p(dan$gene_level_file %||% NA_character_, unique(as.character(dan$gene_pairs$gene2))), error = function(e) numeric())
			if (length(p_gene)) write_dandelion_package_network(dan$gene_pairs, p_gene, rawdir, outdir)
		}
		if (nrow(integration) && !is.null(rawdir)) write_raw_csv(integration, "c2.dandelion_mr_integration.csv", rawdir)
	}
	dirint <- cached$directionality_integration %||% tibble()
	save_plot(plot_c2_directionality(dirint), "c2.Fig10.directionality_causal.png", 16, 12, outdir = outdir)
	if (nrow(dirint) && !is.null(rawdir)) write_raw_csv(dirint, "c2.directionality_causal.csv", rawdir)
	invisible(TRUE)
}

read_c1_cache_for_c2 <- function(outdir) {
	c1file <- file.path(le8_job_dir(outdir, "c1_correlate"), "c1.res.rds")
	wait_seconds <- suppressWarnings(as.numeric(Sys.getenv("C2_C1_WAIT_SECONDS", unset = "60")))
	if (!is.finite(wait_seconds) || wait_seconds < 0) wait_seconds <- 0
	deadline <- Sys.time() + wait_seconds
	announced <- FALSE
	repeat{
		if (file.exists(c1file) && isTRUE(file.size(c1file) > 0)) {
			c1 <- tryCatch(readRDS(c1file), error = function(e) e)
			if (!inherits(c1, "condition")) return(c1)
			last_error <- conditionMessage(c1)
		} else last_error <- "file does not exist or is empty"
		if (Sys.time() >= deadline) {
			stop("C2 requires an existing, readable C1 result: ", c1file,
				" (", last_error, "). C2 does not run C1 automatically.",
				call. = FALSE
			)
		}
		if (!announced) {
			message("C2: waiting up to ", wait_seconds, " seconds for existing C1 result: ", c1file)
			announced <- TRUE
		}
		Sys.sleep(1)
	}
}

run_reverse_mr_stage <- function(iv_list, ygfile, layer) {
	# Disease-liability -> omic MR uses independent disease GWAS instruments as
	# exposure and each full pQTL/mQTL summary as outcome. It complements (and
	# does not replace) clinical lead-time sensitivity because liability is not
	# the same estimand as established disease or treatment.
	leads <- find_disease_leads(ygfile) ; if (!nrow(leads)) return(list(MR = tibble(), audit = tibble()))
	yiv <- read_sumstat_snps(ygfile, leads$SNP, N_DEFAULT) |> filter(SNP %in% leads$SNP)
	if (!nrow(yiv)) return(list(MR = tibble(), audit = tibble(feature = names(iv_list), status = "No disease-IV beta/SE in outcome GWAS")))
	rows <- parallel_map(names(iv_list), function(feature) {
		fs <- iv_list[[feature]]$files
		f <- fs$full %||% NA_character_
		if (length(f) != 1L || is.na(f) || !file.exists(f))
			return(list(result = tibble(), audit = tibble(feature = feature, status = "Full QTL summary missing", n_disease_IV = nrow(yiv), n_QTL_overlap = 0L)))
		qtl <- tryCatch(read_sumstat_snps(f, yiv$SNP, N_DEFAULT), error = function(e) tibble())
		if (!nrow(qtl)) return(list(result = tibble(), audit = tibble(feature = feature, status = "No disease-IV overlap in QTL", n_disease_IV = nrow(yiv), n_QTL_overlap = 0L)))
		rr <- run_mr(yiv, qtl, Y, "disease_liability_to_omic") |> mutate(feature = feature, .before = 1)
		list(result = rr, audit = tibble(feature = feature, status = ifelse(rr$n_IV > 0, "Available", "Harmonization failed"), n_disease_IV = nrow(yiv), n_QTL_overlap = nrow(qtl)))
	})
	mr <- bind_rows(lapply(rows, `[[`, "result")) ; if (nrow(mr)) mr <- mr |>
		mutate(FDR_reverse = p.adjust(pval, "BH")) |>
		arrange(pval)
	list(MR = mr, audit = bind_rows(lapply(rows, `[[`, "audit")), disease_instruments = yiv)
}

plot_bidirectional_mr <- function(forward, reverse, layer) {
	fw <- forward |>
		filter(is.finite(pval)) |>
		group_by(exposure) |>
		slice_min(pval, n = 1, with_ties = FALSE) |>
		ungroup() |>
		transmute(feature = exposure, forward_beta = b, forward_p = pval, forward_FDR = FDR_all)
	rv <- if (nrow(reverse) && all(c("feature", "b", "pval", "FDR_reverse", "n_IV") %in% names(reverse))) reverse |>
		filter(is.finite(pval)) |>
		transmute(feature, reverse_beta = b, reverse_p = pval, reverse_FDR = FDR_reverse, n_reverse_IV = n_IV) else
		tibble(feature = character(), reverse_beta = numeric(), reverse_p = numeric(), reverse_FDR = numeric(), n_reverse_IV = integer())
	z <- full_join(fw, rv, by = "feature") |> mutate(
		class = case_when(
			is.finite(forward_FDR) & forward_FDR < .05 & is.finite(reverse_FDR) & reverse_FDR < .05 ~ "Bidirectional genetic support",
			is.finite(forward_FDR) & forward_FDR < .05 ~ "Omic -> disease only",
			is.finite(reverse_FDR) & reverse_FDR < .05 ~ "Disease liability -> omic only", TRUE ~ "No FDR support"
		),
		x = sign(forward_beta) *  - log10(pmax(forward_p, 1e-300)), y = sign(reverse_beta) *  - log10(pmax(reverse_p, 1e-300)),
		label_p = pmin(forward_p, reverse_p, na.rm = TRUE), label_p = ifelse(is.finite(label_p), label_p, NA_real_)
	)
	pa <- if (!nrow(z)) blank_plot("Bidirectional MR", "No reverse-direction estimate") else ggplot(z, aes(x, y, color = class)) +
		geom_hline(yintercept = 0, color = "grey80") +
		geom_vline(xintercept = 0, color = "grey80") +
		geom_point(alpha = .72, size = 2) +
		geom_text_repel(data = z |> filter(class != "No FDR support", is.finite(label_p)) |> slice_min(label_p, n = 20), aes(label = feature), size = 2.8, max.overlaps = Inf) +
		labs(
			title = "a. Forward and reverse-direction MR", subtitle = "Axes are signed -log10(P); reverse MR estimates disease liability, not post-diagnosis treatment effects",
			x = "Omic -> disease MR", y = "Disease liability -> omic MR", color = NULL
		) +
		theme_5c(9) +
		theme(legend.position = "top")
	cnt <- z |> count(class, name = "biomarkers")
	pb <- if (!nrow(cnt)) blank_plot("b. Direction classes") else ggplot(cnt, aes(biomarkers, fct_reorder(class, biomarkers), fill = class)) +
		geom_col(width = .68, show.legend = FALSE) +
		geom_text(aes(label = biomarkers), hjust =  - .12, fontface = "bold") +
		scale_x_continuous(expand = expansion(mult = c(0, .20))) +
		labs(title = "b. Genetic direction classes", x = "Biomarkers", y = NULL) +
		theme_5c(9)
	pa | pb
}


le8_combine_directionality <- function(mr, c1, dan = list(), reverse_mr = tibble()) {
	d <- as_tibble(c1$directionality %||% tibble()) ; if (!nrow(d)) return(tibble())
	mb <- mr |>
		filter(is.finite(pval)) |>
		group_by(exposure) |>
		slice_min(pval, n = 1, with_ties = FALSE) |>
		ungroup() |>
		transmute(term = exposure, MR_analysis = analysis, MR_beta = b, MR_p = pval, MR_FDR = FDR_all)
	tg <- as_tibble(dan$targets %||% tibble())
	if (nrow(tg) && !"dpg_class" %in% names(tg)) tg$dpg_class <- "DANDELION prioritized target"
	if (nrow(tg) && !"consolidation_eligible" %in% names(tg)) {
		primary_input <- if ("gene_evidence_type" %in% names(tg))
			!str_detect(str_to_lower(coalesce(tg$gene_evidence_type, "")), "magma|adapted") else FALSE
		tested_n <- suppressWarnings(as.numeric((dan$ptrans_dimensions %||% c(NA_real_))[[1]]))
		frac <- if (is.finite(tested_n) && tested_n > 0) nrow(tg) / tested_n else NA_real_
		tg$consolidation_eligible <- primary_input & !(is.finite(frac) && frac > DANDELION_MAX_TARGET_FRACTION)
	}
	dg <- if (nrow(tg)) tg |> transmute(
		term = as.character(gene2), DANDELION_p,
		DANDELION_class = dpg_class %||% "DANDELION prioritized target",
		DANDELION_primary = as.logical(consolidation_eligible)
	) else
		tibble(term = character(), DANDELION_p = numeric(), DANDELION_class = character(), DANDELION_primary = logical())
	rv <- if (nrow(reverse_mr)) reverse_mr |> transmute(term = feature, reverse_MR_beta = b, reverse_MR_p = pval, reverse_MR_FDR = FDR_reverse) else
		tibble(term = character(), reverse_MR_beta = numeric(), reverse_MR_p = numeric(), reverse_MR_FDR = numeric())
	d |>
		left_join(mb, by = "term") |>
		left_join(rv, by = "term") |>
		left_join(dg, by = "term") |>
		mutate(
			MR_support = is.finite(MR_FDR) & MR_FDR < .05,
			reverse_MR_support = is.finite(reverse_MR_FDR) & reverse_MR_FDR < .05,
			DANDELION_sensitivity = is.finite(DANDELION_p),
			DANDELION_support = DANDELION_sensitivity & replace_na(DANDELION_primary, FALSE),
			reactive_warning = coalesce(reactive_compatible, FALSE),
			distal_antecedent = coalesce(distal_support, FALSE),
			causal_interpretation = case_when(
				MR_support & reactive_warning ~ "Forward MR + reactive-compatible (mixed)",
				MR_support & distal_antecedent ~ "Forward MR + distal antecedent",
				reverse_MR_support & !MR_support ~ "Disease liability -> omic MR supported",
				MR_support ~ "Forward MR supported",
				!MR_support & DANDELION_support & reactive_warning ~ "DANDELION + reactive-compatible",
				!MR_support & DANDELION_support ~ "DANDELION-prioritized",
				!MR_support & DANDELION_sensitivity ~ "DANDELION sensitivity only",
				reactive_warning ~ "Observational reactive-compatible",
				distal_antecedent ~ "Observational distal antecedent",
				TRUE ~ "no causal support"
			)
		) |>
		arrange(desc(MR_support), desc(DANDELION_support), p_incident)
}


plot_c2_directionality <- function(z, anchors = unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
							unset = "PCSK9,LPA,GDF15,NTPROBNP,MMP12"
						), ",", fixed = TRUE)[[1]]))) {
	if (!nrow(z)) return(blank_plot("Causal evidence with temporal directionality", "C1 directionality table was unavailable"))
	if (!"DANDELION_sensitivity" %in% names(z)) z$DANDELION_sensitivity <- z$DANDELION_support
	chosen <- unique(c(anchors, z |> filter(MR_support | DANDELION_sensitivity) |>
		arrange(pmin(MR_p, DANDELION_p, na.rm = TRUE)) |> slice_head(n = 24) |> pull(term)))
	long <- z |>
		filter(term %in% chosen) |>
		transmute(term,
			`Incident Cox` = incident_score, `5-y landmark` = landmark5_score,
			`Prevalent logistic` = prevalent_score, `Duration slope` = duration_score,
			`Forward MR` = sign(MR_beta) * pmin(12, - log10(pmax(MR_p, 1e-300))),
			`Reverse MR` = sign(reverse_MR_beta) * pmin(12, - log10(pmax(reverse_MR_p, 1e-300))),
			`DANDELION` = pmin(12, - log10(pmax(DANDELION_p, 1e-300)))
		) |>
		pivot_longer( - term, names_to = "evidence", values_to = "signed_score")
	pa <- ggplot(long, aes(evidence, factor(term, levels = rev(chosen)), fill = signed_score)) +
		geom_tile(color = "white") +
		scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c( - 12, 12), name = "signed\n-log10(P)") +
		labs(title = "a. Genetic evidence beside Yin-Yang temporal evidence", x = NULL, y = NULL) +
		theme_5c(8) +
		theme(axis.text.x = element_text(angle = 35, hjust = 1))
	sc <- z |> filter(is.finite(MR_beta), is.finite(beta_incident))
	pb <- if (!nrow(sc)) blank_plot("b. MR and incident association", "No overlapping estimate") else
		ggplot(sc, aes(MR_beta, beta_incident, color = causal_interpretation)) +
			geom_hline(yintercept = 0, color = "grey85") +
			geom_vline(xintercept = 0, color = "grey85") +
			geom_point(size = 2, alpha = .75) +
			geom_text_repel(data = sc |> filter(term %in% anchors | reactive_warning), aes(label = term), size = 3, max.overlaps = 20, fontface = "bold") +
			labs(title = "b. Genetic causality and reactive signal may coexist", subtitle = "Prevalent/duration evidence is a disease-state warning, not a veto on valid cis-MR", x = "MR beta", y = "Incident Cox beta", color = NULL) +
			theme_5c(9)
	pc <- z |>
		count(causal_interpretation, name = "proteins") |>
		ggplot(aes(proteins, reorder(causal_interpretation, proteins), fill = causal_interpretation)) +
		geom_col(width = .72) +
		guides(fill = "none") +
		labs(title = "c. Integrated interpretation", x = "Proteins", y = NULL) +
		theme_5c(9)
	pa / (pb | pc) + plot_layout(heights = c(1.25, 1))
}


# MR sensitivity, calibrated PGS decomposition and DANDELION analysis

le8_oof_decompose <- function(d, feature, gcol, covars, k = 5, seed = 2026) {
	d <- as.data.frame(d)
	d$.omic <- as.numeric(d[[feature]])
	d$.grs <- as.numeric(d[[gcol]])
	set.seed(seed)
	ord <- order(as.character(d$eid))
	fold <- integer(nrow(d))
	fold[ord] <- sample(rep(seq_len(k), length.out = nrow(d)))
	for (nm in c(".omic_z", ".genetic_z", ".residual_z", ".full_prediction", ".reduced_prediction")) d[[nm]] <- NA_real_
	reports <- list()
	for (f in seq_len(k)) {
		train <- d[fold != f & !is.na(d$.prevalent) & d$.prevalent == 0 & complete.cases(d[, c(
			".omic", ".grs",
			covars
		), drop = FALSE]), , drop = FALSE]
		test <- which(fold == f & complete.cases(d[, c(".omic", ".grs", covars), drop = FALSE]))
		if (nrow(train) < 200 || !length(test) || sd(train$.omic) <= 0 || sd(train$.grs) <= 0)
			next
		om <- mean(train$.omic)
		os <- sd(train$.omic)
		gm <- mean(train$.grs)
		gs <- sd(train$.grs)
		train$.y <- (train$.omic - om) / os
		train$.g <- (train$.grs - gm) / gs
		full <- tryCatch(lm(reformulate(c(".g", covars), ".y"), train), error = function(e) NULL)
		reduced <- tryCatch(lm(reformulate(covars, ".y"), train), error = function(e) NULL)
		if (is.null(full) || is.null(reduced) || !is.finite(coef(full)[".g"]))
			next
		te <- d[test, , drop = FALSE]
		te$.g <- (te$.grs - gm) / gs
		bg <- unname(coef(full)[".g"])
		d$.omic_z[test] <- (te$.omic - om) / os
		d$.genetic_z[test] <- bg * te$.g
		d$.residual_z[test] <- d$.omic_z[test] - d$.genetic_z[test]
		d$.full_prediction[test] <- tryCatch(as.numeric(predict(full, te)), error = function(e) rep(NA_real_, length(test)))
		d$.reduced_prediction[test] <- tryCatch(as.numeric(predict(reduced, te)), error = function(e) rep(
			NA_real_,
			length(test)
		))
		reports[[length(reports) + 1L]] <- tibble(feature,
			fold = f, N_train = nrow(train), N_test = length(test),
			genetic_beta = bg, omic_mean = om, omic_sd = os, pgs_mean = gm, pgs_sd = gs, calibration_rule = "baseline disease-free training fold; future outcomes not consulted"
		)
	}
	d$.calibration_fold <- factor(fold)
	ii <- !is.na(d$.prevalent) & d$.prevalent == 0 & is.finite(d$.full_prediction) & is.finite(d$.reduced_prediction) &
		is.finite(d$.omic_z)
	mse0 <- sum((d$.omic_z[ii] - d$.reduced_prediction[ii]) ^ 2)
	r2 <- if (mse0 > 0)
		1 - sum((d$.omic_z[ii] - d$.full_prediction[ii]) ^ 2) / mse0 else NA_real_
	list(data = d, folds = bind_rows(reports), partial_R2 = r2, N_calibration = sum(ii))
}

risk_set_component_scan <- function(dd, feature, covars, tvar, evar, cuts = c(0, 0.5, 1, 2, 5, 10, 16)) {
	components <- c(Observed = ".omic_z", `PGS-predicted` = ".genetic_z", Residual = ".residual_z")
	map_dfr(seq_len(length(cuts) - 1L), function(i) imap_dfr(components, function(x, label) {
		z <- le8_interval_cox(dd, x, tvar, evar, covars, cuts[i], cuts[i + 1], scale_x = FALSE)
		z |>
			transmute(feature,
				component = label, lead_lo = window_lo, lead_hi = window_hi, lead_mid = time, beta,
				se = std.error, conf.low, conf.high, p = p.value, N_case, N_control, N_risk = N_total, person_years,
				status, effect_measure, unit = "common training-fold omic scale"
			)
	}))
}

run_individual_genetic_decomposition <- function(layer, outdir, top_candidates) {
	sf <- find_c2_score_file(layer)
	empty <- list(status = tibble(status = "PGS unavailable"), summary = tibble(), trajectory = tibble(), folds = tibble())
	if (length(sf) != 1L || is.na(sf))
		return(empty)
	scores <- read_c2_scores(sf)
	biom <- if (layer == "protein")
		read_prot() else read_met()
	features <- unique(c(C2_FIXED_TOP, as.character(top_candidates$feature)))
	features <- head(intersect(features, names(biom)), le8_num_env("C2_DECOMP_MAX", 100))
	mp <- map_c2_score_columns(features, names(scores))
	if (!length(mp))
		return(empty)
	covars <- unique(c(vars.basic, le8_csv_env("C2_LE4_COVARS", "diet.pts,pa.pts,smoke.pts,sleep.pts"), le8_csv_env("C2_TREATMENT_VARS")))
	ph <- read_all(unique(c(
		"eid", "ethnic.c", covars, "birth_date", "date_attend", "date_lost", "date_death",
		paste0("fod_icd10_", Y)
	))) |>
		filter_analysis_cohort() |>
		make_outcome(Y)
	ph$.prevalent <- make_prevalent_status(ph, Y)
	for (nm in c("ph", "scores", "biom")) {
		z <- get(nm)
		z$eid <- as.character(z$eid)
		assign(nm, z)
	}
	d <- inner_join(ph, biom[, c("eid", names(mp)), drop = FALSE], by = "eid") |>
		inner_join(scores[, unique(c("eid", unname(mp))), drop = FALSE], by = "eid")
	rm(ph, scores, biom)
	invisible(gc())
	covars <- intersect(covars, names(d))
	tvar <- paste0(Y, ".t2e")
	evar <- paste0(Y, ".Yt2e")
	results <- lapply(names(mp), function(feature) {
		z <- le8_oof_decompose(d[, unique(c("eid", feature, mp[[feature]], covars, tvar, evar, ".prevalent")),
			drop = FALSE
		], feature, mp[[feature]], covars, k = as.integer(le8_num_env("C2_DECOMP_FOLDS", 5)), seed = SEED)
		dd <- z$data
		cv <- c(covars, ".calibration_fold")
		co <- component_association(dd, ".omic_z", cv, tvar, evar)
		cg <- component_association(dd, ".genetic_z", cv, tvar, evar)
		cr <- component_association(dd, ".residual_z", cv, tvar, evar)
		po <- component_association(dd, ".omic_z", cv, tvar, evar, TRUE)
		pg <- component_association(dd, ".genetic_z", cv, tvar, evar, TRUE)
		pr <- component_association(dd, ".residual_z", cv, tvar, evar, TRUE)
		ji <- joint_component_association(dd, cv, tvar, evar)
		jp <- joint_component_association(dd, cv, tvar, evar, TRUE)
		take <- function(x, term, col) {
			i <- match(term, x$component)
			if (is.na(i))
				NA_real_ else x[[col]][i]
		}
		out <- tibble(feature,
			score_column = mp[[feature]], status = if (any(is.finite(dd$.genetic_z)))
				"ok" else "calibration unavailable", N_calibration = z$N_calibration, N_incident = co[["N"]], incident_events = co[["events"]],
			genetic_beta = mean(z$folds$genetic_beta), genetic_partial_R2 = z$partial_R2, pgs_partial_R2 = z$partial_R2,
			incident_beta_observed = co[["beta"]], incident_se_observed = co[["se"]], incident_p_observed = co[["p"]],
			incident_beta_genetic = cg[["beta"]], incident_se_genetic = cg[["se"]], incident_p_genetic = cg[["p"]],
			incident_beta_residual = cr[["beta"]], incident_se_residual = cr[["se"]], incident_p_residual = cr[["p"]],
			incident_joint_beta_genetic = take(ji, ".genetic_z", "beta"), incident_joint_p_genetic = take(
				ji, ".genetic_z",
				"p"
			), incident_joint_beta_residual = take(ji, ".residual_z", "beta"), incident_joint_p_residual = take(
				ji,
				".residual_z", "p"
			), prevalent_beta_observed = po[["beta"]], prevalent_p_observed = po[["p"]],
			prevalent_beta_genetic = pg[["beta"]], prevalent_p_genetic = pg[["p"]], prevalent_beta_residual = pr[["beta"]],
			prevalent_p_residual = pr[["p"]], prevalent_joint_beta_genetic = take(jp, ".genetic_z", "beta"), prevalent_joint_p_genetic = take(
				jp,
				".genetic_z", "p"
			), prevalent_joint_beta_residual = take(jp, ".residual_z", "beta"), prevalent_joint_p_residual = take(
				jp,
				".residual_z", "p"
			), cross_fitted = TRUE, pgs_scope = "existing external/overlapping COJO score; pQTL GWAS overlap must be audited separately",
			causal_fraction_identifiable = FALSE, incident_genetic_signal_fraction_abs = NA_real_, interpretation = "Residual includes uncaptured genetics, environment, treatment, assay error and disease processes; not a disease-consequence fraction",
			unit = "common training-fold omic scale; components not separately standardized"
		)
		trajectory <- if (feature %in% head(names(mp), le8_num_env("C2_DECOMP_TRAJECTORY_MAX", 20)))
			risk_set_component_scan(dd, feature, cv, tvar, evar) else tibble()
		list(summary = out, folds = z$folds, trajectory = trajectory)
	})
	summary <- bind_rows(lapply(results, `[[`, "summary"))
	folds <- bind_rows(lapply(results, `[[`, "folds"))
	trajectory <- bind_rows(lapply(results, `[[`, "trajectory"))
	h2 <- read_c2_heritability(layer)
	summary <- left_join(summary, h2$data, by = "feature")
	summary$pgs_h2_coverage <- NA_real_ # partial R2 and external h2 are not a captured-causal-fraction estimator.
	for (p in grep("^(incident|prevalent)(_joint)?_p_", names(summary), value = TRUE)) summary[[paste0(
		"FDR_",
		p
	)]] <- p.adjust(summary[[p]], "BH")
	rd <- le8_job_dir(outdir, "c2_cause")
	write_raw_csv(summary, "c2.individual_genetic_decomposition.csv", rd)
	write_raw_csv(trajectory, "c2.genetic_component_leadtime.csv", rd)
	write_raw_csv(folds, "c2.decomposition_calibration_folds.csv", rd)
	list(
		status = tibble(status = "ok", cross_fitted = TRUE, calibration = "baseline disease-free training folds only"),
		summary = summary, trajectory = trajectory, folds = folds, heritability_status = h2$status
	)
}

plot_individual_genetic_decomposition <- function(summary) {
	z <- as_tibble(summary)
	if (!nrow(z) || !"pgs_partial_R2" %in% names(z))
		return(blank_plot("PGS decomposition", "No calibrated scores"))
	z <- z |>
		filter(status == "ok") |>
		slice_head(n = 30)
	a <- ggplot(z, aes(100 * pgs_partial_R2, reorder(feature, pgs_partial_R2))) +
		geom_col() +
		labs(
			title = "a. Out-of-fold partial R² (negative values retained)",
			x = "PGS incremental explained variance (%)", y = NULL
		) +
		theme_5c(9)
	eff <- bind_rows(
		z |>
			transmute(feature, component = "Observed", beta = incident_beta_observed, se = incident_se_observed), z |>
			transmute(feature, component = "PGS-predicted", beta = incident_beta_genetic, se = incident_se_genetic),
		z |>
			transmute(feature, component = "Residual", beta = incident_beta_residual, se = incident_se_residual)
	)
	b <- ggplot(eff, aes(beta, feature, color = component)) +
		geom_vline(xintercept = 0) +
		geom_errorbarh(aes(xmin = beta -
			1.96 * se, xmax = beta + 1.96 * se), height = 0.1, position = position_dodge(width = 0.5)) +
		geom_point(position = position_dodge(width = 0.5)) +
		labs(
			title = "b. Associations on a common omic scale", x = "Log HR per training-fold omic SD", y = NULL,
			color = NULL
		) +
		theme_5c(9)
	c <- ggplot(z, aes(N_incident, feature)) +
		geom_point() +
		labs(
			title = "c. Disease-model sample sizes (not calibration N)",
			x = "Incident-model N", y = NULL
		) +
		theme_5c(9)
	(a | b) / c + plot_annotation(caption = "No causal/consequence percentage is calculated. Cross-fitting calibration does not remove overlap in the original pQTL GWAS.")
}

le8_safe_call <- function(fun, args) {
	f <- names(formals(fun))
	if (!"..." %in% f)
		args <- args[names(args) %in% f]
	do.call(fun, args)
}

find_dandelion_gene_file <- function(ygfile) {
	f <- Sys.getenv("C2_DANDELION_GENE_P_FILE", unset = "")
	type <- Sys.getenv("C2_DANDELION_GENE_EVIDENCE", unset = "unverified")
	allowed <- c("WES_disease", "MAGMA_disease", "other_disease", "unverified")
	if (!type %in% allowed)
		stop("C2_DANDELION_GENE_EVIDENCE must be one of: ", paste(allowed, collapse = ", "))
	if (nzchar(f)) {
		if (!file.exists(f))
			return(list(
				path = NA_character_, evidence_type = "configured gene file missing", primary_eligible = FALSE,
				analysis_class = "unavailable"
			))
		matches <- identical(Sys.getenv("C2_DANDELION_GENE_OUTCOME", unset = ""), Y)
		unfiltered <- identical(Sys.getenv("C2_DANDELION_GENE_SCOPE", unset = "unknown"), "unfiltered")
		independent <- identical(Sys.getenv("C2_DANDELION_SAMPLE_OVERLAP", unset = "unknown"), "none")
		primary <- type == "WES_disease" && matches && unfiltered && independent
		return(list(path = f, evidence_type = type, primary_eligible = primary, analysis_class = if (primary) "WES-supported trans-pQTL adaptation" else "declared/unknown-source sensitivity"))
	}
	f <- if (DANDELION_ALLOW_MAGMA)
		find_magma_gene_file(ygfile) else NA_character_
	list(path = f, evidence_type = "MAGMA_disease", primary_eligible = FALSE, analysis_class = "MAGMA trans-pQTL sensitivity; common-variant evidence is not independent WES evidence")
}

read_gene_level_p <- function(file, target_symbols) {
	if (is.na(file) || !file.exists(file))
		return(numeric())
	d <- data.table::fread(file, showProgress = FALSE, check.names = FALSE)
	gc <- pick_col_local(names(d), c("^GENE$", "^GENE_ID$", "^SYMBOL$", "^GENE_SYMBOL$", "^GENEID$"))
	pc <- Sys.getenv("C2_DANDELION_GENE_P_COLUMN", unset = "")
	if (!nzchar(pc))
		pc <- pick_col_local(names(d), c("^P$", "^PVAL$", "^P_VALUE$", "^BURDEN_P$", "^P_BURDEN$", "^PVAL_BURDEN$"))
	if (is.na(gc) || is.na(pc) || !pc %in% names(d))
		stop("Gene-level disease input needs gene and explicit P columns")
	# Optional prespecified test/mask choice is applied before multiplicity handling.
	test <- Sys.getenv("C2_DANDELION_GENE_TEST", unset = "")
	if (nzchar(test)) {
		tc <- pick_col_local(names(d), c("^TEST$", "^MASK$", "^ANNOTATION$"))
		if (is.na(tc))
			stop("GENE_TEST configured but no TEST/MASK column")
		d <- d[as.character(d[[tc]]) == test, , drop = FALSE]
	}
	z <- tibble(gene = map_magma_gene_ids(d[[gc]], target_symbols), p = as.numeric(d[[pc]])) |>
		filter(!is.na(gene), nzchar(gene), is.finite(p), p >= 0, p <= 1) |>
		mutate(p = pmax(p, .Machine$double.xmin)) |>
		group_by(gene) |>
		summarise(n_tests = n(), p_raw_min = min(p), p = pmin(1, min(p) * n()), .groups = "drop")
	# Bonferroni across masks/tests, not an unadjusted minimum.
	if (!is.null(.le8_analysis_state$rawdir))
		write_raw_csv(z, "c2.dandelion_gene_p_aggregation.csv", .le8_analysis_state$rawdir)
	ans <- setNames(z$p, z$gene)
	attr(ans, "n_gene_tests") <- length(unique(as.character(d[[gc]])[!is.na(d[[gc]])]))
	ans
}

le8_dandelion_family <- function(evidence, universe, alpha = 0.1) {
	d <- evidence
	valid <- is.finite(d$DANDELION_p) & d$DANDELION_p >= 0 & d$DANDELION_p <= 1
	d$global_pair_BH <- d$global_pair_BY <- NA_real_
	d$global_pair_BH[valid] <- p.adjust(d$DANDELION_p[valid], "BH", n = nrow(d))
	d$global_pair_BY[valid] <- p.adjust(d$DANDELION_p[valid], "BY", n = nrow(d))
	d$maxP <- pmax(d$trans_p, d$gene_level_p)
	d$maxP_pair_BH <- NA_real_
	ii <- is.finite(d$maxP)
	d$maxP_pair_BH[ii] <- p.adjust(d$maxP[ii], "BH", n = nrow(d))
	rule <- toupper(Sys.getenv("C2_DANDELION_MULTIPLICITY", unset = "BY"))
	if (!rule %in% c("BH", "BY"))
		stop("C2_DANDELION_MULTIPLICITY must be BH or BY")
	qp <- if (rule == "BY")
		d$global_pair_BY else d$global_pair_BH
	d$significant <- is.finite(qp) & qp <= alpha
	# DACT unavailable is NA, not a non-significant result. MaxP still uses all eligible edges and remains
	# explicitly a different sensitivity method.
	tg <- map_dfr(universe, function(g) {
		all <- d[d$gene2 == g, , drop = FALSE]
		z <- all[is.finite(all$DANDELION_p), , drop = FALSE]
		m <- nrow(z)
		ps <- sort(z$DANDELION_p)
		tibble(
			gene2 = g, n_tested_pairs = m, n_eligible_trans_pairs = nrow(all), DANDELION_p = if (m)
				min(1, min(ps) * nrow(all)) else NA_real_, loo_p = if (m > 1)
				min(1, ps[2] * (nrow(all) - 1)) else if (m == 1)
				1 else NA_real_, n_selected_pairs = sum(z$significant), n_distal_loci = n_distinct(z$locus_id[z$significant]),
			gene_level_p = if (nrow(all))
				all$gene_level_p[1] else NA_real_, best_trans_p = if (nrow(all))
				min(all$trans_p) else NA_real_, maxP_target_p = if (nrow(all))
				min(1, min(all$maxP) * nrow(all)) else NA_real_
		)
	}) |>
		mutate(
			target_BH = p.adjust(DANDELION_p, "BH", n = length(universe)), target_BY = p.adjust(DANDELION_p,
				"BY",
				n = length(universe)
			), maxP_target_BH = p.adjust(maxP_target_p, "BH", n = length(universe)),
			leave_best_edge_BH = p.adjust(loo_p, "BH", n = length(universe)), selected = if (rule == "BY")
				is.finite(target_BY) & target_BY <= alpha else is.finite(target_BH) & target_BH <= alpha, selection_rule = rule, note = "Target P is Bonferroni-min across ALL valid edges; target family includes every target; leave-best-edge is descriptive sensitivity"
		)
	list(pairs = d, targets = tg)
}

le8_dandelion_plot_bundle <- function(dan, outdir) {
	tg <- dan$targets_all %||% tibble()
	e <- dan$evidence_plot %||% tibble()
	caption <- paste(dan$analysis_class %||% "unavailable", "; no causal effect size or mediation fraction is estimated")
	if (!nrow(tg)) {
		for (pair in list(c("c2.Fig6.dandelion.png", "DANDELION targets"), c(
			"c2.Fig7.dandelion_evidence.png",
			"DANDELION evidence"
		), c("c2.Fig9.dandelion_network.png", "DANDELION network"))) save_plot(blank_plot(
			pair[2],
			dan$status %||% "unavailable"
		), pair[1], 12, 7, outdir = outdir)
		return(invisible(NULL))
	}
	z <- tg |>
		arrange(target_BH) |>
		slice_head(n = 24) |>
		mutate(gene2 = factor(gene2, levels = rev(gene2)))
	a <- ggplot(z, aes( - log10(pmax(target_BH, 1e-300)), gene2)) +
		geom_point(aes(size = n_distal_loci, shape = selected)) +
		labs(
			title = "a. Target-level evidence (all targets corrected)", x = "-log10(target BH)", y = NULL, size = "Distinct upstream loci",
			shape = "Selected"
		) +
		theme_5c(9)
	b <- ggplot(z, aes( - log10(pmax(gene_level_p, 1e-300)), gene2)) +
		geom_point() +
		labs(
			title = "b. Gene-level disease evidence",
			x = "-log10(gene disease P)", y = NULL
		) +
		theme_5c(9)
	c <- ggplot(tg, aes( - log10(pmax(target_BH, 1e-300)), - log10(pmax(maxP_target_BH, 1e-300)))) +
		geom_point(aes(shape = selected)) +
		geom_abline(slope = 1, intercept = 0, linetype = 2) +
		labs(
			title = "c. DACT and conservative MaxP sensitivity",
			x = "-log10(DACT target BH)", y = "-log10(MaxP target BH)"
		) +
		theme_5c(9)
	d <- ggplot(z, aes( - log10(pmax(leave_best_edge_BH, 1e-300)), gene2)) +
		geom_point() +
		labs(
			title = "d. Remove the strongest edge (aggregation sensitivity)",
			x = "-log10(leave-best-edge BH)", y = NULL
		) +
		theme_5c(9)
	save_plot((a | b) / (c | d) + plot_annotation(caption = caption), "c2.Fig6.dandelion.png", 17, 12, outdir = outdir)
	if (nrow(e)) {
		a <- ggplot(e, aes( - log10(pmax(trans_p, 1e-300)), - log10(pmax(gene_level_p, 1e-300)))) +
			geom_point(aes(color = significant),
				alpha = 0.4
			) +
			labs(
				title = "a. Regulatory and gene-disease evidence are separate", x = "-log10(trans-pQTL P)",
				y = "-log10(gene disease P)"
			) +
			theme_5c(9)
		b <- dan$exposure_qc |>
			ggplot(aes(n_targets_tested, n_selected)) +
			geom_point() +
			labs(
				title = "b. Broad-regulator audit (not an automatic exclusion)",
				x = "Tested targets per locus", y = "Selected edges"
			) +
			theme_5c(9)
		save_plot((a | b) + plot_annotation(caption = caption), "c2.Fig7.dandelion_evidence.png", 16, 8, outdir = outdir)
	}
	edges <- dan$pairs %||% tibble()
	edges <- edges |>
		filter(gene2 %in% as.character(z$gene2))
	if (nrow(edges)) {
		genes <- unique(edges$gene2)
		loci <- unique(edges$locus_id)
		edges$yg <- match(edges$gene2, genes) / max(1, length(genes))
		edges$yl <- match(edges$locus_id, loci) / max(1, length(loci))
		a <- ggplot(edges) +
			geom_segment(aes(x = 0, xend = 1, y = yl, yend = yg, alpha =  - log10(pmax(
				global_pair_BH,
				1e-300
			)))) +
			geom_text(
				data = unique(edges[, c("locus_id", "yl")]), aes(x = 0, y = yl, label = locus_id),
				hjust = 1, size = 2.6
			) +
			geom_text(
				data = unique(edges[, c("gene2", "yg")]), aes(x = 1, y = yg, label = gene2),
				hjust = 0, size = 2.8
			) +
			xlim( - 0.8, 1.8) +
			theme_void() +
			labs(
				title = "Disease loci → target genes",
				subtitle = "SNP-to-upstream-gene labels are annotations, not proven effector genes", alpha = "Edge evidence"
			)
	} else a <- blank_plot("DANDELION network", "No edge passed the configured full-family threshold")
	save_plot(a + plot_annotation(caption = caption), "c2.Fig9.dandelion_network.png", 15, 10, outdir = outdir)
}

.le8_dandelion_impl <- function(layer, assoc, ann, base, ygfile, rawdir, outdir, top_candidates = tibble(), mode = RUN_Dandelion) {
	empty <- list(
		status = "not_run", pairs = tibble(), targets = tibble(), targets_all = tibble(), gene_pairs = tibble(),
		lead_snps = tibble(), snp_gene_map = tibble(), evidence_plot = tibble(), exposure_qc = tibble(), qtl_audit = tibble(),
		input_audit = tibble(), integration = tibble(), sig_gene2 = character(), non_sig_gene2 = character(), primary_eligible = FALSE,
		consolidation_eligible = FALSE, broad_selection_warning = FALSE, target_fraction = NA_real_, gene_level_file = NA_character_,
		gene_evidence_type = "unavailable", network_file = NA_character_, analysis_class = "not run"
	)
	if (mode == "None" || layer != "protein")
		return(le8_result_update(empty, list(status = if (layer != "protein") "protein-only adaptation" else "RUN_Dandelion=None")))
	input <- find_dandelion_gene_file(ygfile)
	if (is.na(input$path))
		return(le8_result_update(empty, list(status = "gene-level disease input unavailable", analysis_class = input$analysis_class)))
	gene_map <- le8_assay_genes(ann$feature)
	ref <- ann |>
		mutate(gene_name = unname(gene_map[feature])) |>
		filter(!is.na(chr), is.finite(start), is.finite(end)) |>
		arrange(feature) |>
		distinct(gene_name, .keep_all = TRUE) # prespecified lexicographic representative, never minimum P
	p_gene <- read_gene_level_p(input$path, ref$gene_name)
	ref <- ref |>
		filter(gene_name %in% names(p_gene))
	if (mode == "Top")
		ref <- ref |>
			filter(feature %in% top_candidates$feature)
	if (DANDELION_MAX_GENE2 > 0)
		ref <- head(ref, DANDELION_MAX_GENE2) # not ranked by gene-disease P
	if (nrow(ref) < 5)
		return(le8_result_update(empty, list(status = "fewer than five matched annotated targets", gene_level_file = input$path)))
	lead <- find_disease_leads(ygfile)
	src <- attr(lead, "source") %||% "unknown"
	lead <- lead |>
		filter(is.finite(P), P <= DANDELION_GWS, is.finite(POS), !is.na(CHR)) |>
		distinct(SNP, .keep_all = TRUE)
	# Attempt explicit LD validation of the selected disease variants.
	fake <- lead
	fake$source_file <- ygfile
	fake$BETA <- 1
	fake$SE <- 1
	ld <- le8_ld_for_iv(fake, sub("\\.gz$", "", basename(ygfile)))
	ld_verified <- !is.null(ld$R) && all(lead$SNP %in% rownames(ld$R))
	if (ld_verified)
		lead <- le8_ld_greedy(lead, ld$R, le8_num_env("C2_MR_LD_R2", 0.001))
	if (nrow(lead) < 2)
		return(le8_result_update(empty, list(status = "fewer than two eligible disease loci", gene_level_file = input$path)))
	lead$locus_id <- paste0("chr", sub("^chr", "", lead$CHR), ":", lead$POS)
	if (!ld_verified) {
		# Distance groups provide descriptive counts only; not called LD-independent.
		for (ch in unique(lead$CHR)) {
			ix <- which(lead$CHR == ch)
			ix <- ix[order(lead$POS[ix])]
			g <- cumsum(c(TRUE, diff(lead$POS[ix]) > DANDELION_LEAD_BP))
			for (k in unique(g)) {
				jj <- ix[g == k]
				lead$locus_id[jj] <- paste0("chr", sub("^chr", "", ch), ":", min(lead$POS[jj]), "-", max(lead$POS[jj]))
			}
		}
	}
	pt <- build_trans_p_matrix(ref$feature, lead, base)
	qa <- attr(pt, "qtl_audit")
	rownames(pt) <- ref$gene_name[match(rownames(pt), ref$feature)]
	refm <- ref |>
		transmute(gene_name,
			type = "protein_coding", Chromosome = paste0("chr", sub("^chr", "", chr)), start,
			end
		)
	pw <- p_gene[rownames(pt)]
	# Explicit local exclusion before all tests; missing tests remain NA, not P=1.
	for (i in seq_len(nrow(pt))) {
		r <- refm[match(rownames(pt)[i], refm$gene_name), ]
		local <- paste0("chr", sub("^chr", "", lead$CHR)) == r$Chromosome & lead$POS >= r$start - DANDELION_CIS_BP &
			lead$POS <= r$end + DANDELION_CIS_BP
		pt[i, local] <- NA_real_
	}
	med <- NULL
	package_status <- "package unavailable; MaxP sensitivity only"
	if (requireNamespace("DANDELION", quietly = TRUE)) {
		med <- tryCatch(le8_safe_call(DANDELION::med_gene, list(
			p.trans = pt, p.wes = pw, ref.table = as.data.frame(refm),
			gene1.list = colnames(pt), gene1.type = "SNP", SNP.ref = as.data.frame(lead |>
				transmute(SNP, SNPPos = POS, SNPChr = CHR)), target.fdr = DANDELION_FDR, dist = DANDELION_CIS_BP,
			n.cores = N_CORES, verbose = FALSE
		)), error = function(e) e)
		package_status <- if (inherits(med, "condition"))
			conditionMessage(med) else "ok"
		if (inherits(med, "condition"))
			med <- NULL
	}
	dact <- matrix(NA_real_, nrow(pt), ncol(pt), dimnames = dimnames(pt))
	native <- matrix(FALSE, nrow(pt), ncol(pt), dimnames = dimnames(pt))
	if (!is.null(med)) {
		rr <- intersect(rownames(pt), rownames(med$mat.p))
		cc <- intersect(colnames(pt), colnames(med$mat.p))
		dact[rr, cc] <- med$mat.p[rr, cc]
		native[rr, cc] <- med$mat.sig[rr, cc] != 0
	}
	evidence <- tibble(
		gene2 = rep(rownames(pt), times = ncol(pt)), rsid = rep(colnames(pt), each = nrow(pt)),
		trans_p = as.vector(pt), gene_level_p = rep(as.numeric(pw), times = ncol(pt)), DANDELION_p = as.vector(dact),
		package_selected = as.vector(native)
	) |>
		filter(is.finite(trans_p), is.finite(gene_level_p)) |>
		left_join(lead |>
			select(rsid = SNP, locus_id, CHR, POS, outcome_P = P), by = "rsid")
	fam <- le8_dandelion_family(evidence, rownames(pt), DANDELION_FDR)
	ev <- fam$pairs
	alltg <- fam$targets
	# Truncating the target family affects the empirical DACT null mixture.  BY correction cannot repair
	# selected inputs or dependent component P values.
	uncapped <- DANDELION_MAX_GENE2 == 0 && DANDELION_MAX_SNPS == 0
	qtl_unfiltered <- identical(Sys.getenv("C2_DANDELION_QTL_SCOPE", unset = "unknown"), "unfiltered")
	primary <- input$primary_eligible && mode == "All" && uncapped && qtl_unfiltered && ld_verified && package_status ==
		"ok"
	class <- paste(input$analysis_class, if (mode == "Top")
		"; selected Top family is exploratory" else "; all eligible target genes", if (!ld_verified)
		"; disease-locus LD unverified" else "; disease-locus LD verified")
	alltg <- alltg |>
		left_join(ref |>
			select(gene2 = gene_name, feature), by = "gene2") |>
		mutate(
			gene_evidence_type = input$evidence_type, primary_eligible = primary, analysis_class = class, consolidation_eligible = primary &
				selected, gene_level_bonferroni = gene_level_p < 0.05 / max(1, attr(p_gene, "n_gene_tests") %||% length(p_gene)),
			dpg_class = case_when(!is.finite(DANDELION_p) ~ "DACT not tested / unavailable", selected & gene_level_bonferroni ~
				"DACT + gene-disease evidence", selected ~ "DACT-prioritized hypothesis", gene_level_bonferroni ~
				"Gene-disease evidence only", TRUE ~ "Not selected")
		)
	tg <- alltg |>
		filter(selected)
	fraction <- nrow(tg) / max(1, nrow(alltg))
	broad <- fraction > DANDELION_MAX_TARGET_FRACTION
	# Broad selection triggers an audit flag, not a data-dependent veto threshold.
	for (nm in c("alltg", "tg")) {
		z <- get(nm)
		z$target_fraction <- fraction
		z$broad_selection_warning <- broad
		assign(nm, z)
	}
	pd <- ev |>
		filter(significant) |>
		left_join(
			alltg |>
				select(gene2, feature, primary_eligible, consolidation_eligible, analysis_class, gene_evidence_type),
			by = "gene2"
		)
	mp <- read_dandelion_snp_gene_map(lead)
	if (nrow(mp))
		pd <- left_join(pd, as_tibble(mp) |>
			group_by(SNP) |>
			summarise(gene1 = paste(unique(GeneSymbol), collapse = ";"), .groups = "drop"), by = c(rsid = "SNP"))
	if (!"gene1" %in% names(pd))
		pd$gene1 <- pd$rsid
	pd$gene1[is.na(pd$gene1)] <- pd$rsid[is.na(pd$gene1)]
	gp <- pd |>
		mutate(region = locus_id) # locus is the upstream unit; labels do not establish a causal gene.
	qc <- ev |>
		group_by(rsid, locus_id) |>
		summarise(
			n_targets_tested = n(), n_selected = sum(significant), selected_fraction = mean(significant),
			.groups = "drop"
		)
	audit <- tibble(metric = c(
		"mode", "gene_file", "gene_evidence", "analysis_class", "primary_eligible", "disease_LD_verified",
		"DACT_package_status", "target_universe", "valid_trans_tests", "selected_targets", "target_fraction", "multiplicity",
		"scope_boundary"
	), value = as.character(c(
		mode, input$path, input$evidence_type, class, primary, ld_verified,
		package_status, nrow(alltg), nrow(ev), nrow(tg), fraction, paste0(
			"All eligible DACT pairs, missing results retained in family; target Bonferroni-min then BH/BY over all targets; selection=",
			Sys.getenv("C2_DANDELION_MULTIPLICITY", unset = "BY")
		), "trans-pQTL adaptation; neither causal direction nor protein mediation fraction is identified"
	)))
	audit <- bind_rows(audit, tibble(metric = c(
		"package_version", "target_exposure_uncapped", "QTL_unfiltered",
		"sample_overlap", "native_selection", "extension_selection"
	), value = c(
		if (requireNamespace("DANDELION",
			quietly = TRUE
		)) as.character(utils::packageVersion("DANDELION")) else "unavailable", as.character(uncapped),
		as.character(qtl_unfiltered), Sys.getenv("C2_DANDELION_SAMPLE_OVERLAP", unset = "unknown"), "package_selected: native per-exposure q-value rule",
		"selected: project-specific target aggregation and BH/BY; not the native DANDELION target rule"
	)))
	ans <- le8_result_update(empty, list(
		status = if (package_status == "ok") "ok" else "DACT unavailable; see MaxP sensitivity",
		mode = mode, result = med, pairs = pd, targets = tg, targets_all = alltg, gene_pairs = gp, lead_snps = lead,
		snp_gene_map = mp, evidence_plot = ev, exposure_qc = qc, qtl_audit = qa, input_audit = audit, sig_gene2 = tg$gene2[tg$gene_level_bonferroni %in%
			TRUE], non_sig_gene2 = tg$gene2[!tg$gene_level_bonferroni %in% TRUE], gene_level_file = input$path,
		gene_evidence_type = input$evidence_type, analysis_class = class, primary_eligible = primary, consolidation_eligible = primary,
		target_fraction = fraction, broad_selection_warning = broad, ptrans_dimensions = dim(pt)
	))
	tabs <- list(
		pairs = pd, targets = tg, targets_all = alltg, gene_pairs = gp, lead_snps = lead, snp_to_cis_gene = mp,
		exposure_qc = qc, qtl_coverage = qa, input_audit = audit
	)
	for (nm in names(tabs)) write_raw_csv(tabs[[nm]], paste0("c2.dandelion_", nm, ".csv"), rawdir)
	data.table::fwrite(ev, file.path(rawdir, "c2.dandelion_tested_pairs.csv.gz"), compress = "gzip")
	le8_dandelion_plot_bundle(ans, outdir)
	invisible(tabs) # appended after the legacy workbook writer by the analysis completion handler
	ans
}

plot_dandelion_results <- function(tg, pd, gene_pair, outdir) invisible(NULL)

plot_dandelion_landscape <- function(evidence_plot, tg, outdir, gene_evidence_type) invisible(NULL)

write_dandelion_package_network <- function(gene_pair, p_gene, rawdir, outdir) NA_character_

plot_dandelion_mr_integration <- function(dandelion, mr, assoc, outdir) {
	tg <- dandelion$targets_all %||% tibble()
	if (!nrow(tg)) {
		save_plot(blank_plot("DANDELION and MR", "DANDELION target tests unavailable"), "c2.Fig8.dandelion_mr_integration.png",
			14, 8,
			outdir = outdir
		)
		return(tibble())
	}
	m <- mr |>
		filter(analysis == "cis") |>
		select(feature = exposure, MR_beta = b, MR_p = pval, MR_FDR = FDR_all)
	z <- tg |>
		left_join(m, by = "feature") |>
		left_join(assoc |>
			select(feature = term, observed_beta = beta, observed_p = p.value), by = "feature")
	a <- ggplot(z, aes( - log10(pmax(target_BH, 1e-300)), - log10(pmax(MR_FDR, 1e-300)))) +
		geom_point(aes(shape = primary_eligible)) +
		labs(
			title = "a. Orthogonal evidence shown side by side", x = "-log10(DACT target BH)", y = "-log10(cis-MR BH)",
			shape = "WES/All/LD-eligible"
		) +
		theme_5c(9)
	b <- z |>
		arrange(target_BH) |>
		slice_head(n = 20) |>
		ggplot(aes(observed_beta, reorder(gene2, observed_beta))) +
		geom_point() +
		geom_vline(xintercept = 0) +
		labs(
			title = "b. Observed association does not determine regulatory direction", x = "Observed incident log HR",
			y = NULL
		) +
		theme_5c(9)
	save_plot(a | b, "c2.Fig8.dandelion_mr_integration.png", 16, 9, outdir = outdir)
	z
}

read_dandelion_snp_gene_map <- function(lead) {
	f <- Sys.getenv("C2_DANDELION_SNP_GENE_MAP", unset = "")
	if (!nzchar(f) || !file.exists(f))
		return(tibble(SNP = character(), GeneSymbol = character()))
	d <- data.table::fread(f, showProgress = FALSE)
	sc <- pick_col_local(names(d), c("^SNP$", "^RSID$"))
	gc <- pick_col_local(names(d), c("^GENESYMBOL$", "^GENE_SYMBOL$", "^GENE$", "^SYMBOL$"))
	if (is.na(sc) || is.na(gc))
		stop("SNP gene map needs SNP and GeneSymbol columns")
	tibble(SNP = as.character(d[[sc]]), GeneSymbol = as.character(d[[gc]])) |>
		filter(SNP %in% lead$SNP, !is.na(GeneSymbol), nzchar(GeneSymbol)) |>
		distinct()
}

run_dandelion_step <- function(layer, assoc, ann, base, ygfile, rawdir, outdir, top_candidates = tibble(), mode = RUN_Dandelion) {
	tryCatch(
		{
			ans <- .le8_dandelion_impl(layer, assoc, ann, base, ygfile, rawdir, outdir, top_candidates, mode)
			failure <- file.path(rawdir, "c2.dandelion_failure.csv")
			if (file.exists(failure))
				unlink(failure)
			ans
		},
		error = function(e) {
			msg <- conditionMessage(e)
			warning("DANDELION input/method failure: ", msg, call. = FALSE)
			au <- tibble(status = "failed", message = msg)
			write_raw_csv(au, "c2.dandelion_failure.csv", rawdir)
			list(
				status = paste("failed:", msg), pairs = tibble(), targets = tibble(), targets_all = tibble(), gene_pairs = tibble(),
				evidence_plot = tibble(), exposure_qc = tibble(), qtl_audit = tibble(), input_audit = au, lead_snps = tibble(),
				snp_gene_map = tibble(), integration = tibble(), gene_level_file = NA_character_, gene_evidence_type = "unavailable",
				network_file = NA_character_, primary_eligible = FALSE, analysis_class = "Failed optional analysis; not a negative finding"
			)
		}
	)
}

grade_c2_evidence <- function(mr, layer, mrlink2 = list(results = tibble())) {
	primary <- if (layer == "protein")
		"cis" else "local"
	z <- mr |>
		filter(analysis == primary)
	if (!nrow(z))
		return(z)
	z <- z |>
		mutate(heterogeneous = ifelse(n_IV > 1 & is.finite(Q_p), Q_p < 0.05, NA), heterogeneity_tested = n_IV >
			1 & is.finite(Q_p), egger_tested = n_IV >= 3 & is.finite(egger_intercept_p), egger_warning = ifelse(egger_tested,
			egger_intercept_p < 0.05, NA
		), weighted_median_concordant = ifelse(n_IV >= 3 & is.finite(tsmr_weighted_median_b),
			sign(tsmr_weighted_median_b) == sign(b), NA
		), steiger_ok = NA, many_IV_warning = n_IV > 10, instrument_architecture = case_when(
			n_IV ==
				0 ~ "unavailable", n_IV == 1 ~ "single IV", n_IV <= 5 ~ "oligogenic (2-5 IVs)", n_IV <= 10 ~ "multi-IV (6-10)",
			TRUE ~ "many-IV (>10)"
		), evidence_grade = case_when(
			!is.finite(pval) ~ "U: unavailable", FDR_all <
				0.05 & n_IV == 1 ~ "S: single-IV support", FDR_all < 0.05 & heterogeneity_tested & egger_tested & !heterogeneous &
				!egger_warning ~ "A: diagnostics-compatible", FDR_all < 0.05 & coalesce(heterogeneous, FALSE) & coalesce(
				weighted_median_concordant,
				FALSE
			) ~ "B: heterogeneous", FDR_all < 0.05 ~ "C: sensitivity unresolved", pval < 0.05 ~ "D: nominal only",
			TRUE ~ "E: not detected"
		), grade_reason = case_when(
			n_IV == 1 ~ "Wald estimate; Egger, weighted-median and heterogeneity tests unavailable",
			evidence_grade == "A: diagnostics-compatible" ~ "FDR support; no detected Q/Egger warning, not proof of instrument validity",
			evidence_grade == "B: heterogeneous" ~ "FDR support with heterogeneity; weighted-median direction agrees",
			evidence_grade == "U: unavailable" ~ "Not tested / missing valid input", TRUE ~ "Provisional statistical evidence; requires same-region/signal corroboration"
		))
	rr <- as_tibble(mrlink2$results %||% tibble())
	z$MRlink2_p <- NA_real_
	z$MRlink2_FDR <- NA_real_
	if (nrow(rr)) {
		tc <- pick_col_local(names(rr), c("^TRAIT$", "^EXPOSURE$"))
		pc <- pick_col_local(names(rr), c("^P\\(ALPHA\\)$", "^P_ALPHA$"))
		if (!is.na(tc) && !is.na(pc)) {
			v <- tibble(exposure = as.character(rr[[tc]]), p = as.numeric(rr[[pc]])) |>
				filter(is.finite(p)) |>
				group_by(exposure) |>
				summarise(MRlink2_p = pmin(1, min(p) * n()), .groups = "drop") |>
				mutate(MRlink2_FDR = p.adjust(MRlink2_p, "BH"))
			i <- match(z$exposure, v$exposure)
			z$MRlink2_p <- v$MRlink2_p[i]
			z$MRlink2_FDR <- v$MRlink2_FDR[i]
		}
	}
	arrange(z, pval)
}

plot_c2_evidence_grades <- function(z, layer) {
	if (!nrow(z))
		return(blank_plot("C2 evidence grades", "No estimable result"))
	chosen <- unique(c(C2_FIXED_TOP, head(z$exposure[order(z$pval)], 24)))
	d <- z |>
		filter(exposure %in% chosen) |>
		mutate(exposure = factor(exposure, levels = rev(chosen)), lo = b - 1.96 * se, hi = b + 1.96 * se)
	pa <- ggplot(d, aes(b, exposure, color = evidence_grade)) +
		geom_vline(xintercept = 0) +
		geom_errorbarh(aes(
			xmin = lo,
			xmax = hi
		), height = 0.12) +
		geom_point(aes(shape = instrument_architecture), size = 2) +
		labs(
			title = "a. Marginal-effect MR and uncertainty",
			x = "MR effect (GWAS outcome scale)", y = NULL, color = NULL, shape = NULL
		) +
		theme_5c(9)
	pb <- z |>
		count(evidence_grade) |>
		ggplot(aes(n, reorder(evidence_grade, n), fill = evidence_grade)) +
		geom_col(show.legend = FALSE) +
		labs(
			title = "b. Missing diagnostics are not passed tests",
			x = "Biomarkers", y = NULL
		) +
		theme_5c(9)
	pa | pb
}

build_genetic_score_manifest <- function(iv_list, layer) {
	imap_dfr(iv_list, function(obj, feature) {
		z <- as_tibble(obj$score_instruments %||% tibble())
		if (!nrow(z))
			return(tibble())
		z |>
			transmute(feature, SNP, CHR, POS,
				effect_allele = EA, other_allele = NEA, weight = BETA, SE, P, EAF,
				N, component = paste0("COJO PGS / ", analysis), effect_type = "joint PGS weight; not MR input",
				interpretation = "Existing prot.pgs.rds/met.pgs.rds unchanged. Not a protein concentration measured at birth."
			)
	})
}

run_mrlink2_step <- function(...) {
	tryCatch(le8_execute_mrlink2(...), error = function(e) {
		warning("Optional MR-link-2 failed: ", conditionMessage(e), call. = FALSE)
		invisible(1L)
	})
}

integrate_c2_directionality <- function(mr, c1, dan = list(), reverse_mr = tibble()) {
	# Do not pick whichever of cis or trans has the smaller P as causal evidence.
	z <- le8_combine_directionality(mr[mr$analysis %in% c("cis", "local"), , drop = FALSE], c1, dan, reverse_mr)
	if (nrow(z)) {
		if ("causal_interpretation" %in% names(z))
			z$causal_interpretation <- gsub("no causal support", "not established / not detected", z$causal_interpretation,
				fixed = TRUE
			)
		z$MR_scope <- "cis/local primary; trans/distal remains in its own evidence table"
	}
	z
}

instrument_availability_audit <- function(...) {
	z <- le8_count_available_instruments(...)
	if (nrow(z) && "n_independent" %in% names(z)) {
		z$n_COJO_selected <- z$n_independent
		z$n_independent <- NA_integer_
		z$LD_interpretation <- "Input availability is not LD verification; final pruned counts and fallback status are in c2.MR_all.csv"
	}
	z
}

run_c2_layer <- function(layer = c("protein", "metabolite")) {
	if (LE8_REUSE_RESULTS) return(le8_restore_outputs(match.arg(layer), "c2_cause"))
	# Check reusable results and initialize the analysis output directory.
	layer <- match.arg(layer)
	le8_begin_analysis(layer, "c2_cause")
	.le8_analysis_env <- environment()
	on.exit(le8_finish_analysis(layer, "c2_cause", .le8_analysis_env), add = TRUE)

	layer <- match.arg(layer) ; outdir <- if (layer == "protein") out.prot else out.met ; setwd2(outdir)
	rawdir <- le8_job_dir(outdir, LE8_JOB) ; dir.create(rawdir, recursive = TRUE, showWarnings = FALSE) ; cache <- file.path(rawdir, "c2.res.rds")
	method_scope <- c2_method_scope_audit(layer) ; write_raw_csv(method_scope, "c2.method_scope_audit.csv", rawdir)
	h2_file0 <- find_c2_heritability_file(layer)
	h2_signature <- if (is.na(h2_file0)) "h2=missing" else paste0(
		"h2=", h2_file0,
		";size=", file.size(h2_file0), ";mtime=", format(file.info(h2_file0)$mtime, tz = "UTC", usetz = TRUE)
	)
	component_covariate_signature <- paste0(
		"LE4=", Sys.getenv("C2_LE4_COVARS",
			unset = "diet.pts,pa.pts,smoke.pts,sleep.pts"
		), ";treatment=",
		Sys.getenv("C2_TREATMENT_VARS", unset = Sys.getenv("C1_TREATMENT_VARS", unset = ""))
	)
	mode_signature <- paste0(
		"MRlink2=", RUN_MRlink2, ";Dandelion=", RUN_Dandelion,
		";", h2_signature, ";", component_covariate_signature
	)
	pgs_file0 <- find_c2_score_file(layer)
	pgs_check <- length(pgs_file0) == 1L && !is.na(pgs_file0) && nzchar(pgs_file0) &&
		file.exists(pgs_file0) && file.size(pgs_file0) > 0
	pgs_expected <- file.path(indir, "Rdata", if (layer == "protein") "prot.pgs.rds" else "met.pgs.rds")
	pgs_signature <- if (is.na(pgs_file0)) "PGS=missing" else paste0(
		"PGS=", pgs_file0, ";size=", file.size(pgs_file0),
		";mtime=", format(file.info(pgs_file0)$mtime, tz = "UTC", usetz = TRUE)
	)
	if (!pgs_check) message(
		"C2/", layer, ": PGS check: missing ", pgs_expected,
		"; skip individual genetic decomposition and lead-time analysis"
	)
	if (cache_valid(cache)) {
		old <- tryCatch(readRDS(cache), error = function(e) NULL)
		if (is.list(old) && all(c("meta", "MR", "MR_reverse", "MR_best") %in% names(old)) &&
			!grepl("^failed", old$DANDELION$status %||% "") && file.exists(file.path(rawdir, "c2.genetic_discovery_scope.csv")) && identical(as.character(data.table::fread(file.path(rawdir, "c2.genetic_discovery_scope.csv"))$scope[1]), C2_FEATURE_SCOPE)) {
			cache_message(paste0("C2/", layer), cache) ; return(le8_restore_outputs(layer, "c2_cause"))
		}
		message("C2/", layer, ": cache incomplete; reusing stage caches and rebuilding final outputs")
	}
	base0 <- if (layer == "protein") dir.X else dir.met.gwas
	c1 <- read_c1_cache_for_c2(outdir) ; assoc <- (c1$association %||% c1$pwas_incident %||% c1$MWAS) |> as_tibble() ; if (!"beta" %in% names(assoc)) assoc <- assoc |> mutate(beta = safe_log(estimate))
	base <- base0 ; ann <- layer_annotation(layer, unique(assoc$term)) ; ranked <- assoc |>
		filter(is.finite(p.value)) |>
		arrange(p.value) |>
		pull(term)
	# MR requires an independent COJO set plus the corresponding full summary
	# statistics for allele recovery. A lone cis extract or LD (.ldr.cojo) file
	# must not consume a candidate slot.
	mapped <- ranked[vapply(ranked, function(x) {
		fs <- find_qtl_files(x, base, layer) ; !is.na(fs$joint) && !is.na(fs$full)
	}, logical(1))]
	fixed_assayed <- intersect(C2_FIXED_TOP, unique(assoc$term))
	fixed_mapped <- fixed_assayed[vapply(fixed_assayed, function(x) {
		fs <- find_qtl_files(x, base, layer) ; !is.na(fs$joint) && !is.na(fs$full)
	}, logical(1))]
	if (C2_FEATURE_SCOPE == "all_qtl") {
		# Disease-association P values neither rank nor exclude assays in this arm.
		universe <- sort(unique(assoc$term))
		mapped <- universe[vapply(universe, function(x) {
			fs <- find_qtl_files(x, base, layer) ; !is.na(fs$joint) && !is.na(fs$full)
		}, logical(1))]
		selected <- if (MAX_FEATURES == 0) mapped else head(mapped, MAX_FEATURES)
		mr_candidates <- unique(c(fixed_mapped, selected))
		audit_features <- unique(c(fixed_assayed, mr_candidates))
	} else {
		selected <- if (MAX_FEATURES == 0) mapped else head(mapped, MAX_FEATURES)
		mr_candidates <- unique(c(fixed_mapped, selected))
		audit_features <- unique(c(fixed_assayed, head(ranked, if (MAX_FEATURES == 0) length(ranked) else MAX_FEATURES), mr_candidates))
	}
	write_raw_csv(tibble(
		scope = C2_FEATURE_SCOPE, maximum = MAX_FEATURES, mapped_traits = length(mapped),
		analyzed_traits = length(mr_candidates), observational_significance_gate = FALSE,
		observational_ranking = C2_FEATURE_SCOPE == "observational_top",
		caveat = if (C2_FEATURE_SCOPE == "all_qtl" && MAX_FEATURES > 0) "Lexicographic resource cap: not an exhaustive scan" else "See full instrument/QC coverage; PGS significance is not MR"
	), "c2.genetic_discovery_scope.csv", rawdir)
	if (!length(mr_candidates)) stop("C2/", layer, ": no QTL files mapped under ", base, call. = FALSE)
	input_table <- tibble(feature = audit_features, top_observational = feature %in% head(ranked, if (MAX_FEATURES == 0) length(ranked) else MAX_FEATURES), MR_candidate = feature %in% mr_candidates)
	write_raw_csv(input_table, "c2.input_features.csv", rawdir) ; message("C2/", layer, ": auditing ", length(audit_features), " traits; ", length(mr_candidates), " have mapped QTL files")
	stage_started <- le8_stage_start(paste0("C2/", layer, " instruments"))
	iv_cache <- file.path(rawdir, "c2.instruments.rds") ; iv_list <- read_stage_cache(iv_cache, C2_STAGE_VERSION)
	if (is.null(iv_list)) {
		message("C2/", layer, ": build independent instrument stage")
		iv_list <- setNames(map(audit_features, ~ read_qtl_instruments(.x, base, layer, ann)), audit_features)
		write_stage_cache(iv_list, iv_cache, C2_STAGE_VERSION)
	} else message("C2/", layer, ": reuse independent instrument stage")
	ygfile <- get_y_gwas_file(Y, TRUE)
	le8_stage_done(paste0("C2/", layer, " instruments"), stage_started)
	stage_started <- le8_stage_start(paste0("C2/", layer, " MR"))
	mr_cache <- file.path(rawdir, "c2.mr.rds") ; mr_stage <- read_stage_cache(mr_cache, C2_STAGE_VERSION)
	if (is.null(mr_stage)) {
		all_snps <- unique(unlist(map(iv_list, ~ .$instruments$SNP)))
		ygwas <- read_sumstat_snps(ygfile, all_snps, N_DEFAULT) ; if (!nrow(ygwas)) stop("No outcome-GWAS variants overlapped QTL instruments.", call. = FALSE)
		availability <- instrument_availability_audit(iv_list, ygwas, layer)
		mr <- imap_dfr(iv_list, function(obj, feature) {
			iv <- obj$instruments ; if (!nrow(iv)) return(tibble()) ; map_dfr(unique(iv$analysis), ~ run_mr(iv |> filter(analysis == .x), ygwas, feature, .x))
		})
		if (!nrow(mr)) stop("C2/", layer, ": no harmonized MR estimate.", call. = FALSE)
		mr <- mr |>
			group_by(analysis) |>
			mutate(FDR_analysis = p.adjust(pval, "BH")) |>
			ungroup() |>
			mutate(FDR_all = p.adjust(pval, "BH"))
		mr_stage <- list(MR = mr, availability = availability) ; write_stage_cache(mr_stage, mr_cache, C2_STAGE_VERSION)
	} else {
		mr <- mr_stage$MR ; availability <- mr_stage$availability ; message("C2/", layer, ": reuse harmonization/MR stage")
	}
	write_raw_csv(availability, "c2.instrument_availability.csv", rawdir) ; write_raw_csv(mr, "c2.MR_all.csv", rawdir)
	le8_stage_done(paste0("C2/", layer, " MR"), stage_started)
	stage_started <- le8_stage_start(paste0("C2/", layer, " reverse-MR"))
	reverse_cache <- file.path(rawdir, "c2.reverse_mr.rds") ; reverse_stage <- read_stage_cache(reverse_cache, C2_STAGE_VERSION)
	if (is.null(reverse_stage)) {
		message("C2/", layer, ": disease-liability -> omic reverse MR")
		reverse_stage <- run_reverse_mr_stage(iv_list, ygfile, layer) ; write_stage_cache(reverse_stage, reverse_cache, C2_STAGE_VERSION)
	} else message("C2/", layer, ": reuse reverse-direction MR stage")
	le8_stage_done(paste0("C2/", layer, " reverse-MR"), stage_started)
	reverse_mr <- reverse_stage$MR %||% tibble() ; reverse_audit <- reverse_stage$audit %||% tibble()
	write_raw_csv(reverse_mr, "c2.reverse_MR_all.csv", rawdir) ; write_raw_csv(reverse_audit, "c2.reverse_MR_audit.csv", rawdir)
	top_candidates <- select_c2_top_candidates(layer, assoc, mr, names(iv_list))
	write_raw_csv(top_candidates, "c2.top_candidates.csv", rawdir)
	genetic_score_manifest <- build_genetic_score_manifest(iv_list, layer)
	write_raw_tsv(genetic_score_manifest, "c2.genetic_score_weights.tsv", rawdir)
	individual_decomposition <- if (pgs_check)
		run_individual_genetic_decomposition(layer, outdir, top_candidates) else
		list(
			status = tibble(status = "not run", detail = paste0("PGS check failed; score file not found: ", pgs_expected)),
			summary = tibble(), trajectory = tibble()
		)
	write_raw_csv(individual_decomposition$status %||% tibble(), "c2.individual_genetic_decomposition_status.csv", rawdir)
	write_raw_csv(individual_decomposition$heritability_status %||% tibble(), "c2.heritability_status.csv", rawdir)
	save_plot(plot_individual_genetic_decomposition(individual_decomposition$summary %||% tibble()),
		"c2.Fig12.genetic_decomposition.png", 17, 13,
		outdir = outdir
	)
	save_plot(plot_component_leadtime(individual_decomposition$trajectory %||% tibble()),
		"c2.Fig13.genetic_component_leadtime.png", 17, 11,
		outdir = outdir
	)
	save_plot(plot_bidirectional_mr(mr, reverse_mr, layer), "c2.Fig11.bidirectional_mr.png", 16, 8.5, outdir = outdir)
	fig1 <- plot_c2_fig1(mr, assoc, layer) ; save_c2_plot(fig1$plot, "c2.Fig1.prots.top.png", 17, 13.5, outdir = outdir)
	save_c2_plot(plot_qtl_variance(mr, layer), "c2.Fig2.pQTL_R2.png", 15.5, 10.5, outdir = outdir)
	fig3 <- plot_c2_fig2(mr, assoc, layer) ; save_plot(fig3$plot, "c2.Fig3.effect_concordance.png", 16, 12.5, outdir = outdir)
	write_raw_csv(fig3$wide, "c2.cis_trans_comparison.csv", le8_job_dir(outdir, "c2_cause"))
	le8_emit_restored_c2(mr, assoc, layer, outdir)
	save_plot(plot_c2_fig4(mr, layer), "c2.Fig4.sensitivity_architecture.png", 15.5, 11.5, outdir = outdir)
	best <- fig1$best |> mutate(y =  - log10(pmax(pval, 1e-300)), direction = case_when(FDR_all < .05 & b > 0 ~ "Positive", FDR_all < .05 & b < 0 ~ "Inverse", TRUE ~ "NS"), label = ifelse(min_rank(pval) <= 25, exposure, NA_character_))
	top_r2 <- mr |>
		filter(is.finite(r2_median), n_IV > 0) |>
		mutate(
			r2_median_pct = 100 * r2_median, r2_q25_pct = 100 * r2_q25, r2_q75_pct = 100 * r2_q75,
			r2_p90_pct = 100 * r2_p90, r2_max_pct = 100 * r2_max
		) |>
		arrange(desc(r2_p90_pct))
	arch <- mr |>
		filter(is.finite(median_F), is.finite(r2_median)) |>
		mutate(
			evidence = case_when(FDR_all < .05 ~ "MR FDR < 0.05", Q_p < .05 ~ "Heterogeneous", TRUE ~ "Other"),
			label = ifelse(min_rank(pval) <= 15, exposure, NA_character_)
		)

	# DANDELION is complementary to MR: it asks whether a distal disease SNP has
	# trans-regulatory evidence to a gene/protein that itself has disease-genetic evidence.
	dan_cache <- file.path(rawdir, "c2.dandelion.rds") ; dandelion <- read_stage_cache(dan_cache, paste(RUN_Dandelion, C2_STAGE_VERSION, sep = "/"))
	if (!is.null(dandelion) && grepl("^failed", dandelion$status %||% "")) dandelion <- NULL
	if (is.null(dandelion)) {
		dandelion <- le8_stage(paste0("C2/", layer, " Dandelion"), run_dandelion_step(
			layer, assoc, ann, base, ygfile, rawdir, outdir,
			top_candidates, RUN_Dandelion
		), paste0("mode=", RUN_Dandelion))
		write_stage_cache(dandelion, dan_cache, paste(RUN_Dandelion, C2_STAGE_VERSION, sep = "/"))
	} else le8_stage(paste0("C2/", layer, " Dandelion"), invisible(NULL), "cache=reused")

	if (layer == "protein") {
		plot_dandelion_results(dandelion$targets %||% tibble(), dandelion$pairs %||% tibble(), dandelion$gene_pairs %||% tibble(), outdir)
		plot_dandelion_landscape(dandelion$evidence_plot %||% tibble(), dandelion$targets %||% tibble(), outdir, dandelion$gene_evidence_type %||% "gene-level disease association")
		dandelion$integration <- plot_dandelion_mr_integration(dandelion, mr, assoc, outdir)
		write_raw_csv(dandelion$integration %||% tibble(), "c2.dandelion_mr_integration.csv", rawdir)
		network_file <- dandelion$network_file %||% NA_character_
		network_missing <- length(network_file) != 1L || is.na(network_file) || !nzchar(network_file) || !file.exists(network_file)
		if (nrow(dandelion$gene_pairs %||% tibble()) && network_missing) {
			p_gene <- tryCatch(read_gene_level_p(dandelion$gene_level_file %||% NA_character_, unique(as.character(dandelion$gene_pairs$gene2))), error = function(e) numeric())
			if (length(p_gene)) dandelion$network_file <- write_dandelion_package_network(dandelion$gene_pairs, p_gene, rawdir, outdir)
		}
		if (!identical(dandelion$status %||% "", "ok")) {
			why <- paste0("DANDELION unavailable: ", dandelion$status %||% "no usable result")
			for (i in 6 : 9) save_plot(blank_plot(paste0("C2 Figure ", i), why),
				c(
					"c2.Fig6.dandelion.png", "c2.Fig7.dandelion_evidence.png",
					"c2.Fig8.dandelion_mr_integration.png", "c2.Fig9.dandelion_network.png"
				)[[i - 5]], 10, 6,
				outdir = outdir
			)
		} else if (!nrow(dandelion$gene_pairs %||% tibble())) {
			save_plot(blank_plot("DANDELION network", "No significant gene-pair network was available"),
				"c2.Fig9.dandelion_network.png", 10, 6,
				outdir = outdir
			)
		}
	} else {
		why <- "DANDELION is a gene/protein regulatory-network analysis and is not defined for metabolites"
		for (i in 6 : 9) save_plot(blank_plot(paste0("C2 Figure ", i), why),
			c(
				"c2.Fig6.dandelion.png", "c2.Fig7.dandelion_evidence.png",
				"c2.Fig8.dandelion_mr_integration.png", "c2.Fig9.dandelion_network.png"
			)[[i - 5]], 10, 6,
			outdir = outdir
		)
	}

	directionality_integration <- integrate_c2_directionality(mr, c1, dandelion, reverse_mr)
	write_raw_csv(directionality_integration, "c2.directionality_causal.csv", rawdir)
	save_plot(plot_c2_directionality(directionality_integration), "c2.Fig10.directionality_causal.png", 16, 12, outdir = outdir)

	jobs_all <- build_mrlink2_jobs(iv_list, ann, layer, ygfile)
	jobs_audit <- if (nrow(jobs_all)) jobs_all |>
		left_join(top_candidates |> select(trait = feature, top_rank, top_reason), by = "trait") |>
		mutate(
			mode = RUN_MRlink2, selected = case_when(
				RUN_MRlink2 == "All" ~ TRUE,
				RUN_MRlink2 == "Top" ~ is.finite(top_rank), TRUE ~ FALSE
			),
			selection_reason = case_when(
				selected & RUN_MRlink2 == "All" ~ "RUN_MRlink2=All",
				selected ~ top_reason, RUN_MRlink2 == "None" ~ "RUN_MRlink2=None", TRUE ~ "outside Top set"
			)
		) else
		tibble(trait = character(), mode = character(), selected = logical(), selection_reason = character())
	write_raw_csv(jobs_audit, "c2.mrlink2_job_audit.csv", rawdir)
	jobs <- if (nrow(jobs_all)) jobs_audit |>
		filter(selected) |>
		select(all_of(names(jobs_all))) else jobs_all
	le8_stage(paste0("C2/", layer, " MR-link-2"), run_mrlink2_step(layer, rawdir, jobs, ygfile, RUN_MRlink2), paste0("mode=", RUN_MRlink2, " jobs=", nrow(jobs)))
	mrlink2 <- read_mrlink2_results(rawdir, RUN_MRlink2) ; plot_mrlink2_results(mrlink2, outdir)
	# Older completed runs sometimes retained the combined estimates but not the
	# generated job manifest. Recover a minimal provenance audit instead of
	# publishing a non-empty result figure beside an empty audit.
	if (!nrow(jobs_audit) && nrow(mrlink2$results %||% tibble())) {
		trait_col <- pick_col_local(names(mrlink2$results), c("^TRAIT$", "^EXPOSURE$"))
		region_col <- pick_col_local(names(mrlink2$results), c("^REGION$"))
		trait <- if (is.na(trait_col)) rep(NA_character_, nrow(mrlink2$results)) else as.character(mrlink2$results[[trait_col]])
		region <- if (is.na(region_col)) rep(NA_character_, nrow(mrlink2$results)) else as.character(mrlink2$results[[region_col]])
		jobs_audit <- tibble(
			trait = trait, region = region, mode = RUN_MRlink2, selected = TRUE,
			selection_reason = "Recovered from completed MR-link-2 output; original job manifest absent"
		) |>
			filter(!is.na(trait), nzchar(trait)) |>
			distinct()
		write_raw_csv(jobs_audit, "c2.mrlink2_job_audit.csv", rawdir)
	}
	evidence_grades <- grade_c2_evidence(mr, layer, mrlink2)
	write_raw_csv(evidence_grades, "c2.evidence_grades.csv", rawdir)
	save_plot(plot_c2_evidence_grades(evidence_grades, layer),
		"c2.Fig14.evidence_grades.png", 17, 9,
		outdir = outdir
	)

	out <- list(
		meta = module_meta(layer, extra = list(
			code_version = C2_CODE_VERSION, mode_signature = mode_signature,
			pgs_signature = pgs_signature, h2_signature = h2_signature,
			RUN_MRlink2 = RUN_MRlink2, RUN_Dandelion = RUN_Dandelion
		)), MR = mr, MR_reverse = reverse_mr, MR_reverse_audit = reverse_audit, genetic_score_manifest = genetic_score_manifest, individual_decomposition = individual_decomposition, MR_best = best, observational = assoc, instrument_availability = availability,
		Fig1_summary = fig1$bars, Fig1_forest = fig1$forest, Fig3 = fig3$wide, R2_QTL = top_r2,
		architecture = arch, top_candidates = top_candidates, evidence_grades = evidence_grades,
		DANDELION = dandelion, directionality_integration = directionality_integration,
		MRLink2_job_audit = jobs_audit, MRLink2_jobs = jobs, MRLink2 = mrlink2, method_scope = method_scope
	)
	saveRDS(out, cache, compress = "xz")
	write_xlsx2(list(
		MR_all = mr, MR_reverse = reverse_mr, MR_reverse_audit = reverse_audit,
		method_scope = method_scope, evidence_grades = evidence_grades, top_candidates = top_candidates,
		genetic_score_manifest = genetic_score_manifest,
		genetic_decomp_status = individual_decomposition$status %||% tibble(),
		heritability_status = individual_decomposition$heritability_status %||% tibble(),
		genetic_decomp_summary = individual_decomposition$summary %||% tibble(),
		genetic_leadtime = individual_decomposition$trajectory %||% tibble(),
		MR_best = best, instrument_availability = availability, evidence_overlap = fig1$bars, effect_forest = fig1$forest, cis_local_vs_trans_distal = fig3$wide, QTL_R2 = top_r2, instrument_architecture = arch,
		DANDELION_input_audit = dandelion$input_audit %||% tibble(), DANDELION_QTL_coverage = dandelion$qtl_audit %||% tibble(), DANDELION_exposure_QC = dandelion$exposure_qc %||% tibble(),
		DANDELION_lead_snps = dandelion$lead_snps %||% tibble(), DANDELION_snp_gene_map = dandelion$snp_gene_map %||% tibble(), DANDELION_pairs = dandelion$pairs %||% tibble(), DANDELION_gene_pairs = dandelion$gene_pairs %||% tibble(), DANDELION_targets = dandelion$targets %||% tibble(), DANDELION_MR_integration = dandelion$integration %||% tibble(),
		directionality_causal = directionality_integration,
		MRLink2_job_audit = jobs_audit, MRLink2_jobs = jobs,
		MRLink2_results = mrlink2$results %||% tibble(),
		MRLink2_status = mrlink2$status %||% tibble()
	), "c2.out.xlsx")
	finalize_outputs(LE8_JOB, outdir) ; out
}

if (prot_DO) invisible(le8_stage("C2/protein", run_c2_layer("protein")))
if (met_DO) invisible(le8_stage("C2/metabolite", run_c2_layer("metabolite")))
