# CKM shared and stage-dependent associations.
# This file defines functions only. Incident analyses use cause-specific Cox
# models: non-CVD death ends follow-up and is not an incident CVD event.
# Association labels are observational candidates, never causal declarations.

ckm_a_default <- function(x, y) if (is.null(x)) y else x
ckm_a_bind <- function(x) data.table::rbindlist(x, use.names = TRUE, fill = TRUE)
ckm_a_csv <- function(x) {
    if (is.null(x) || !length(x) || all(is.na(x))) return(character())
    unique(Filter(nzchar, trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))))
}

ckm_a_options <- function(cfg) {
    a <- utils::modifyList(list(
        stages = 1:3, stage_column = "stage", covariates = c("age", "sex", "tdi", "PC1", "PC2"),
        covariates_by_source = list(), min_n = 100L, min_events = 20L,
        fdr = 0.05, equivalence_hr = 1.2, lag_years = c(0, 2, 5),
        workers = 1L, ph_test = TRUE, yang = TRUE, yang_death = TRUE
    ), ckm_a_default(cfg$association, list()))
    a$stages <- sort(unique(as.integer(unlist(a$stages))))
    a$covariates <- unique(as.character(unlist(a$covariates)))
    a$lag_years <- sort(unique(c(0, as.numeric(unlist(a$lag_years)))))
    a$workers <- as.integer(a$workers)
    if (length(a$stages) < 2L || anyNA(a$stages) || any(!a$stages %in% 0:3)) stop("association.stages must contain at least two stages from 0:3")
    if (anyNA(a$lag_years) || any(!is.finite(a$lag_years)) || any(a$lag_years < 0)) stop("association.lag_years must be finite and nonnegative")
    if (length(a$equivalence_hr) != 1L || !is.finite(a$equivalence_hr) || a$equivalence_hr <= 1) stop("association.equivalence_hr must exceed 1")
    if (length(a$fdr) != 1L || !is.finite(a$fdr) || a$fdr <= 0 || a$fdr >= 1) stop("association.fdr must be between 0 and 1")
    if (length(a$workers) != 1L || is.na(a$workers) || a$workers < 1L) stop("association.workers must be positive")
    if (!is.finite(a$min_n) || a$min_n < 10 || !is.finite(a$min_events) || a$min_events < 5) stop("Invalid association minimum sample/event counts")
    a
}

ckm_a_manifest <- function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!all(c("feature", "source", "column") %in% names(x))) stop("Exposure manifest requires feature, source, and column")
    if (!"reference" %in% names(x) && "ref" %in% names(x)) x$reference <- x$ref
    defaults <- list(type = "continuous", transform = "none", unit = "", reference = "", levels = "",
                     scale = "sd", adjustment_exclude = "", lipid_block = FALSE, definition_component = FALSE)
    for (v in names(defaults)) if (!v %in% names(x)) x[[v]] <- defaults[[v]]
    for (v in c("feature", "source", "column", "type", "transform", "unit", "reference", "levels", "scale", "adjustment_exclude")) {
        x[[v]] <- as.character(x[[v]])
        if (v %in% c("reference", "levels", "unit", "adjustment_exclude")) x[[v]][is.na(x[[v]])] <- ""
    }
    if (!nrow(x)) stop("Exposure manifest is empty")
    if (anyNA(x[, c("feature", "source", "column")]) || any(!nzchar(x$feature)) || any(!nzchar(x$source)) || any(!nzchar(x$column))) stop("Exposure manifest has blank identifiers")
    if (anyDuplicated(paste(x$source, x$feature, sep = "\r"))) stop("Exposure manifest repeats a source/feature pair")
    if (any(!x$type %in% c("continuous", "binary", "categorical"))) stop("Exposure type must be continuous, binary, or categorical")
    if (any(!x$transform %in% c("none", "log1p", "log2"))) stop("Exposure transform must be none, log1p, or log2")
    if (any(x$type != "continuous" & x$transform != "none")) stop("Only continuous exposures can be log transformed")
    if (any(x$type == "continuous" & !x$scale %in% c("sd", "unit"))) stop("Continuous exposure scale must be sd or unit")
    x
}

ckm_a_meta <- function(m, lag = 0, endpoint = "incident_stage4") {
    data.frame(feature = m$feature, source = m$source, column = m$column, type = m$type,
               transform = m$transform, unit = m$unit, reference = m$reference,
               lipid_block = m$lipid_block, definition_component = m$definition_component,
               endpoint = endpoint, lag_years = lag, stringsAsFactors = FALSE)
}

ckm_a_prepare_exposure <- function(x, m, reference_rows) {
    invalid <- 0L
    if (m$type == "continuous") {
        if (is.factor(x)) x <- as.character(x)
        z <- suppressWarnings(as.numeric(x))
        z[!is.finite(z)] <- NA_real_
        if (m$transform == "log1p") {
            invalid <- sum(!is.na(z) & z < 0); z[!is.na(z) & z < 0] <- NA_real_; z <- log1p(z)
        }
        if (m$transform == "log2") {
            invalid <- sum(!is.na(z) & z <= 0); z[!is.na(z) & z <= 0] <- NA_real_; z <- log2(z)
        }
        ref <- z[reference_rows & is.finite(z)]
        center <- if (length(ref)) mean(ref) else NA_real_
        scale <- if (m$scale == "sd") stats::sd(ref) else 1
        if (!is.finite(scale) || scale <= 0 || !length(ref)) {
            return(list(ok = FALSE, status = "constant_or_unavailable_exposure", terms = "per_SD", center = center, scale = scale))
        }
        val <- if (m$scale == "sd") (z - center) / scale else z
        mat <- matrix(val, ncol = 1L, dimnames = list(NULL, ".x1"))
        return(list(ok = TRUE, x = mat, terms = if (m$scale == "sd") "per_SD" else "per_unit",
                    center = center, scale = scale, reference = "", invalid = invalid))
    }
    if (m$type == "binary") {
        if (is.factor(x)) x <- as.character(x)
        if (is.logical(x)) x <- as.integer(x)
        z <- suppressWarnings(as.numeric(x))
        bad <- !is.na(x) & (!is.finite(z) | !z %in% 0:1)
        if (any(bad)) return(list(ok = FALSE, status = "binary_exposure_not_coded_0_1", terms = "1_vs_0", center = 0, scale = 1))
        return(list(ok = TRUE, x = matrix(z, ncol = 1L, dimnames = list(NULL, ".x1")),
                    terms = "1_vs_0", center = 0, scale = 1, reference = "0", invalid = 0L))
    }
    z <- as.character(x); z[is.na(x) | !nzchar(trimws(z))] <- NA_character_
    lev <- if (nzchar(m$levels)) strsplit(m$levels, "|", fixed = TRUE)[[1L]] else sort(unique(z[!is.na(z)]))
    lev <- unique(lev[nzchar(lev)])
    if (length(lev) < 2L) return(list(ok = FALSE, status = "fewer_than_two_exposure_levels", terms = "not_estimable", center = 0, scale = 1))
    if (any(!is.na(z) & !z %in% lev)) return(list(ok = FALSE, status = "unlisted_exposure_level", terms = "not_estimable", center = 0, scale = 1))
    reference <- if (nzchar(m$reference)) m$reference else lev[1L]
    if (!reference %in% lev) stop("Manifest reference is not among the levels for ", m$feature)
    lev <- c(reference, setdiff(lev, reference))
    mat <- vapply(lev[-1L], function(v) as.numeric(z == v), numeric(length(z)))
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1L)
    colnames(mat) <- paste0(".x", seq_len(ncol(mat)))
    list(ok = TRUE, x = mat, terms = paste(lev[-1L], "vs", reference), center = 0, scale = 1,
         reference = reference, invalid = 0L)
}

ckm_a_covariates <- function(m, a) {
    source_cov <- a$covariates_by_source[[m$source]]
    setdiff(unique(c(a$covariates, as.character(unlist(source_cov)))),
            unique(c(m$column, ckm_a_csv(m$adjustment_exclude))))
}

ckm_a_design <- function(d, prep, covariates, rows, stage_column, time_col = "time", event_col = "event", lag = 0, logistic = FALSE) {
    z <- data.frame(.t = if (logistic) rep(1, nrow(d)) else as.numeric(d[[time_col]]) - lag,
                    .e = as.integer(d[[event_col]]), .s = as.character(d[[stage_column]]), stringsAsFactors = FALSE)
    z <- cbind(z, prep$x)
    cmap <- stats::setNames(paste0(".c", seq_along(covariates)), covariates)
    for (v in covariates) {
        x <- d[[v]]
        if (is.character(x) || is.factor(x)) {
            x <- as.character(x); x[is.na(x) | !nzchar(trimws(x))] <- NA_character_; x <- factor(x)
        } else if (is.logical(x)) x <- as.integer(x)
        if (is.numeric(x)) x[!is.finite(x)] <- NA_real_
        z[[cmap[[v]]]] <- x
    }
    ok <- rows & stats::complete.cases(z)
    ok <- ok & is.finite(z$.t) & z$.t > 0 & z$.e %in% 0:1
    z <- droplevels(z[which(ok), , drop = FALSE])
    constants <- names(cmap)[vapply(cmap, function(v) length(unique(z[[v]])) < 2L, logical(1))]
    keep <- unname(cmap[setdiff(names(cmap), constants)])
    z$.s <- factor(z$.s)
    list(data = z, covariates = keep, covariate_map = cmap,
         dropped_constant_covariates = paste(constants, collapse = ","), rows = which(ok))
}

ckm_a_formula <- function(xterms, covariates, stratify = TRUE, logistic = FALSE) {
    strata <- survival::strata
    rhs <- c(xterms, covariates, if (stratify) "strata(.s)")
    stats::as.formula(paste(if (logistic) ".e" else "survival::Surv(.t,.e)", "~", paste(rhs, collapse = " + ")), env = environment())
}

ckm_a_fit <- function(formula, d, logistic = FALSE, keep_x = FALSE) {
    warnings <- character()
    fit <- tryCatch(withCallingHandlers({
        if (logistic) stats::glm(formula, data = d, family = stats::binomial(), control = stats::glm.control(maxit = 50)) else
            survival::coxph(formula, data = d, ties = "efron", singular.ok = FALSE, x = keep_x, y = keep_x,
                            model = FALSE, control = survival::coxph.control(iter.max = 50))
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
    warning_text <- paste(unique(warnings), collapse = "; ")
    if (inherits(fit, "error")) return(list(ok = FALSE, status = paste0("fit_error: ", conditionMessage(fit)), warning = warning_text))
    beta <- stats::coef(fit); V <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    bad_warning <- any(grepl("infinite|did not converge|ran out of iterations|Loglik converged before|fitted probabilities numerically 0 or 1", warnings, ignore.case = TRUE))
    if (bad_warning || (logistic && !isTRUE(fit$converged)) || any(!is.finite(beta)) || is.null(V) || any(!is.finite(V)) || any(diag(V) <= 0)) {
        return(list(ok = FALSE, status = "nonconverged_or_nonestimable_model", warning = warning_text))
    }
    list(ok = TRUE, fit = fit, beta = beta, V = V, status = "ok", warning = warning_text)
}

ckm_a_wald <- function(beta, V) {
    if (!length(beta) || any(!is.finite(beta)) || any(!is.finite(V))) return(NA_real_)
    w <- tryCatch(as.numeric(crossprod(beta, solve(V, beta))), error = function(e) NA_real_)
    if (!is.finite(w) || w < -1e-8) return(NA_real_)
    stats::pchisq(max(0, w), df = length(beta), lower.tail = FALSE)
}

ckm_a_effects <- function(fit, C, terms) {
    beta <- as.numeric(C %*% fit$beta)
    V <- C %*% fit$V %*% t(C)
    se <- sqrt(pmax(0, diag(V)))
    data.frame(term = terms, beta = beta, se = se, effect = exp(beta), ci_low = exp(beta - 1.96 * se),
               ci_high = exp(beta + 1.96 * se), p = 2 * stats::pnorm(-abs(beta / se)), stringsAsFactors = FALSE)
}

ckm_a_empty <- function(m, a, lag, status, terms = "not_estimable", endpoint = "incident_stage4") {
    meta <- ckm_a_meta(m, lag, endpoint)
    summary <- cbind(meta, data.frame(n = 0L, events = 0L, n_before_complete = 0L, n_missing = 0L,
        n_transform_invalid = 0L, scaling_mean = NA_real_, scaling_sd = NA_real_,
        beta = NA_real_, se = NA_real_, hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        p_shared = NA_real_, p_interaction = NA_real_, interaction_df = NA_integer_, p_equivalence = NA_real_,
        p_ph = NA_real_, p_ph_global = NA_real_, equivalence_hr = a$equivalence_hr,
        status = status, interaction_status = status, ph_status = "not_run",
        dropped_constant_covariates = "", warning = "", stringsAsFactors = FALSE))
    stage <- ckm_a_bind(lapply(a$stages, function(s) {
        cbind(meta[rep(1L, length(terms)), , drop = FALSE], data.frame(stage = s, term = terms, n = 0L, events = 0L,
            beta = NA_real_, se = NA_real_, hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
            p_stage = NA_real_, status = status, estimation = "joint_interaction_model",
            dropped_constant_covariates = "", stringsAsFactors = FALSE))
    }))
    pooled <- cbind(meta[rep(1L, length(terms)), , drop = FALSE], data.frame(term = terms, beta = NA_real_, se = NA_real_,
        hr = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_, status = status, stringsAsFactors = FALSE))
    pairs <- utils::combn(a$stages, 2L)
    contrasts <- ckm_a_bind(lapply(seq_len(ncol(pairs)), function(i) {
        cbind(meta[rep(1L, length(terms)), , drop = FALSE], data.frame(stage_a = pairs[1L, i], stage_b = pairs[2L, i], term = terms,
            beta_difference = NA_real_, se = NA_real_, ratio_hr = NA_real_, ci90_low = NA_real_, ci90_high = NA_real_,
            p_difference = NA_real_, p_equivalence = NA_real_, status = status, stringsAsFactors = FALSE))
    }))
    list(summary = summary, pooled = pooled, stage = stage, contrasts = contrasts)
}

# A sparse stage must not erase estimable effects in the other stages. These
# descriptive fallback estimates do not rescue an unavailable interaction test.
ckm_a_stage_fallback <- function(ans, z, design, xnames, terms, a) {
    for (s in a$stages) {
        ii <- which(ans$stage$stage == s)
        zs <- droplevels(z[as.character(z$.s) == as.character(s), , drop = FALSE])
        data.table::set(ans$stage, i = ii, j = "estimation", value = "separate_stage_model")
        if (nrow(zs) < a$min_n || sum(zs$.e) < a$min_events) {
            data.table::set(ans$stage, i = ii, j = "status", value = "insufficient_stage_sample_or_events"); next
        }
        if (any(vapply(zs[xnames], function(x) length(unique(x)) < 2L, logical(1)))) {
            data.table::set(ans$stage, i = ii, j = "status", value = "constant_or_absent_exposure_level"); next
        }
        constants <- design$covariates[vapply(design$covariates, function(v) length(unique(zs[[v]])) < 2L, logical(1))]
        cname <- names(design$covariate_map)[design$covariate_map %in% constants]
        data.table::set(ans$stage, i = ii, j = "dropped_constant_covariates", value = paste(cname, collapse = ","))
        fit <- ckm_a_fit(ckm_a_formula(xnames, setdiff(design$covariates, constants), stratify = FALSE), zs)
        if (!fit$ok) {data.table::set(ans$stage, i = ii, j = "status", value = fit$status); next}
        C <- matrix(0, nrow = length(xnames), ncol = length(fit$beta), dimnames = list(NULL, names(fit$beta)))
        C[cbind(seq_along(xnames), match(xnames, names(fit$beta)))] <- 1
        es <- ckm_a_effects(fit, C, terms)
        ps <- ckm_a_wald(fit$beta[xnames], fit$V[xnames, xnames, drop = FALSE])
        for (v in c("beta", "se", "hr", "ci_low", "ci_high", "p", "p_stage", "status")) {
            value <- switch(v, beta = es$beta, se = es$se, hr = es$effect, ci_low = es$ci_low,
                            ci_high = es$ci_high, p = es$p, p_stage = ps, status = "ok")
            data.table::set(ans$stage, i = ii, j = v, value = value)
        }
    }
    ans
}

ckm_a_incident <- function(d, m, prep, a, lag) {
    ans <- ckm_a_empty(m, a, lag, if (isTRUE(prep$ok)) "insufficient_data" else prep$status, prep$terms)
    if (!isTRUE(prep$ok)) return(ans)
    ans$summary$reference <- prep$reference; ans$stage$reference <- prep$reference
    ans$pooled$reference <- prep$reference; ans$contrasts$reference <- prep$reference
    ans$summary$scaling_mean <- prep$center; ans$summary$scaling_sd <- prep$scale
    ans$summary$n_transform_invalid <- prep$invalid
    rows <- !is.na(d[[a$stage_column]]) & d[[a$stage_column]] %in% a$stages &
        is.finite(d$time) & d$time > lag & !is.na(d$status) & d$status %in% 0:2
    ans$summary$n_before_complete <- sum(rows)
    covariates <- ckm_a_covariates(m, a)
    design <- ckm_a_design(d, prep, covariates, rows, a$stage_column, lag = lag)
    z <- design$data; ans$summary$n <- nrow(z); ans$summary$events <- sum(z$.e)
    ans$summary$n_missing <- sum(rows) - nrow(z)
    ans$summary$dropped_constant_covariates <- design$dropped_constant_covariates
    for (s in a$stages) {
        ii <- ans$stage$stage == s
        ans$stage$n[ii] <- sum(as.character(z$.s) == as.character(s))
        ans$stage$events[ii] <- sum(z$.e[as.character(z$.s) == as.character(s)])
    }
    if (nrow(z) < a$min_n || sum(z$.e) < a$min_events) return(ans)
    xnames <- colnames(prep$x)
    if (any(vapply(z[xnames], function(x) length(unique(x)) < 2L, logical(1)))) {
        ans$summary$status <- "constant_or_absent_exposure_level"; return(ans)
    }
    shared <- ckm_a_fit(ckm_a_formula(xnames, design$covariates), z, keep_x = isTRUE(a$ph_test))
    if (!shared$ok) {
        ans$summary$status <- shared$status; ans$summary$warning <- shared$warning; return(ans)
    }
    ans$summary$status <- "ok"; ans$summary$warning <- shared$warning
    ans$summary$p_shared <- ckm_a_wald(shared$beta[xnames], shared$V[xnames, xnames, drop = FALSE])
    C <- matrix(0, nrow = length(xnames), ncol = length(shared$beta), dimnames = list(NULL, names(shared$beta)))
    C[cbind(seq_along(xnames), match(xnames, names(shared$beta)))] <- 1
    e <- ckm_a_effects(shared, C, prep$terms)
    ans$pooled[, c("beta", "se", "hr", "ci_low", "ci_high", "p", "status")] <-
        data.frame(beta = e$beta, se = e$se, hr = e$effect, ci_low = e$ci_low, ci_high = e$ci_high, p = e$p, status = "ok")
    if (nrow(e) == 1L) {
        ans$summary[, c("beta", "se", "hr", "ci_low", "ci_high")] <-
            data.frame(beta = e$beta, se = e$se, hr = e$effect, ci_low = e$ci_low, ci_high = e$ci_high)
    }
    if (isTRUE(a$ph_test)) {
        ph <- tryCatch(suppressWarnings(survival::cox.zph(shared$fit, transform = "km", terms = FALSE, global = TRUE)), error = function(e) e)
        if (inherits(ph, "error")) ans$summary$ph_status <- paste0("error: ", conditionMessage(ph)) else {
            vals <- ph$table[intersect(xnames, rownames(ph$table)), "p"]
            if (length(vals) == length(xnames) && all(is.finite(vals))) {
                # Each dummy is its own model term. Bonferroni gives one valid,
                # conservative exposure-level PH test for a categorical feature.
                ans$summary$p_ph <- min(1, length(vals) * min(vals)); ans$summary$ph_status <- "ok"
            } else ans$summary$ph_status <- "nonestimable_exposure_PH_test"
            if ("GLOBAL" %in% rownames(ph$table)) ans$summary$p_ph_global <- ph$table["GLOBAL", "p"]
        }
    }
    counts <- unique(ans$stage[, c("stage", "n", "events")])
    if (any(counts$n < a$min_n | counts$events < a$min_events)) {
        ans$summary$interaction_status <- "insufficient_stage_sample_or_events"
        ans$stage$status <- "insufficient_stage_sample_or_events"
        ans$contrasts$status <- "insufficient_stage_sample_or_events"
        return(ckm_a_stage_fallback(ans, z, design, xnames, prep$terms, a))
    }
    iterms <- character()
    for (s in a$stages[-1L]) for (j in seq_along(xnames)) {
        nm <- paste0(".i", j, "s", s)
        z[[nm]] <- z[[xnames[j]]] * as.integer(as.character(z$.s) == as.character(s))
        iterms <- c(iterms, nm)
    }
    full <- ckm_a_fit(ckm_a_formula(c(xnames, iterms), design$covariates), z)
    if (!full$ok) {
        ans$summary$interaction_status <- full$status; ans$stage$status <- full$status; ans$contrasts$status <- full$status
        ans$summary$warning <- paste(Filter(nzchar, c(shared$warning, full$warning)), collapse = "; ")
        return(ckm_a_stage_fallback(ans, z, design, xnames, prep$terms, a))
    }
    lr <- 2 * (full$fit$loglik[2L] - shared$fit$loglik[2L])
    if (!is.finite(lr) || lr < -1e-5) {
        ans$summary$interaction_status <- "invalid_nested_likelihood"
        return(ckm_a_stage_fallback(ans, z, design, xnames, prep$terms, a))
    }
    ans$summary$interaction_status <- "ok"; ans$summary$interaction_df <- length(iterms)
    ans$summary$p_interaction <- stats::pchisq(max(0, lr), df = length(iterms), lower.tail = FALSE)
    slopes <- list()
    for (s in a$stages) {
        Cs <- matrix(0, nrow = length(xnames), ncol = length(full$beta), dimnames = list(NULL, names(full$beta)))
        Cs[cbind(seq_along(xnames), match(xnames, names(full$beta)))] <- 1
        if (s != a$stages[1L]) {
            nm <- paste0(".i", seq_along(xnames), "s", s)
            Cs[cbind(seq_along(xnames), match(nm, names(full$beta)))] <- 1
        }
        slopes[[as.character(s)]] <- Cs
        es <- ckm_a_effects(full, Cs, prep$terms)
        ps <- ckm_a_wald(as.numeric(Cs %*% full$beta), Cs %*% full$V %*% t(Cs))
        ii <- which(ans$stage$stage == s)
        data.table::set(ans$stage, i = ii, j = "beta", value = es$beta)
        data.table::set(ans$stage, i = ii, j = "se", value = es$se)
        data.table::set(ans$stage, i = ii, j = "hr", value = es$effect)
        data.table::set(ans$stage, i = ii, j = "ci_low", value = es$ci_low)
        data.table::set(ans$stage, i = ii, j = "ci_high", value = es$ci_high)
        data.table::set(ans$stage, i = ii, j = "p", value = es$p)
        data.table::set(ans$stage, i = ii, j = "p_stage", value = ps)
        data.table::set(ans$stage, i = ii, j = "status", value = "ok")
    }
    pairs <- utils::combn(a$stages, 2L); margin <- log(a$equivalence_hr)
    for (i in seq_len(ncol(pairs))) {
        sa <- pairs[1L, i]; sb <- pairs[2L, i]
        D <- slopes[[as.character(sa)]] - slopes[[as.character(sb)]]
        ed <- ckm_a_effects(full, D, prep$terms)
        # Intersection-union TOST: equality is not inferred from P_interaction>0.05.
        peq <- pmax(stats::pnorm(-(ed$beta + margin) / ed$se), stats::pnorm((ed$beta - margin) / ed$se))
        ii <- which(ans$contrasts$stage_a == sa & ans$contrasts$stage_b == sb)
        for (v in c("beta_difference", "se", "ratio_hr", "ci90_low", "ci90_high", "p_difference", "p_equivalence", "status")) {
            value <- switch(v, beta_difference = ed$beta, se = ed$se, ratio_hr = ed$effect,
                ci90_low = exp(ed$beta - stats::qnorm(.95) * ed$se), ci90_high = exp(ed$beta + stats::qnorm(.95) * ed$se),
                p_difference = ed$p, p_equivalence = peq, status = "ok")
            data.table::set(ans$contrasts, i = ii, j = v, value = value)
        }
    }
    # Equivalence must hold for every stage pair and every prespecified contrast.
    if (all(is.finite(ans$contrasts$p_equivalence))) ans$summary$p_equivalence <- max(ans$contrasts$p_equivalence)
    ans
}

ckm_a_yang <- function(d, m, prep, a, endpoint) {
    terms <- prep$terms
    meta <- ckm_a_meta(m, 0, endpoint)
    out <- cbind(meta[rep(1L, length(terms)), , drop = FALSE], data.frame(term = terms, n = 0L, events = 0L,
        beta = NA_real_, se = NA_real_, effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
        p_global = NA_real_, effect_measure = if (endpoint == "baseline_stage4_marker") "OR" else "HR",
        scaling_mean = ckm_a_default(prep$center, NA_real_), scaling_sd = ckm_a_default(prep$scale, NA_real_),
        status = "insufficient_data", dropped_constant_covariates = "", warning = "", stringsAsFactors = FALSE))
    if (!isTRUE(prep$ok)) {out$status <- prep$status; return(out)}
    logistic <- endpoint == "baseline_stage4_marker"
    dd <- as.data.frame(d)
    if (logistic) {
        dd$.yang_event <- as.integer(dd[[a$stage_column]] == 4L)
        rows <- !is.na(dd[[a$stage_column]]) & dd[[a$stage_column]] %in% c(a$stages, 4L)
        event_col <- ".yang_event"; time_col <- "time"
    } else {
        if (!all(c("time_death", "event_death") %in% names(dd))) {out$status <- "mortality_endpoint_unavailable"; return(out)}
        rows <- !is.na(dd[[a$stage_column]]) & dd[[a$stage_column]] == 4L &
            is.finite(dd$time_death) & dd$time_death > 0 & !is.na(dd$event_death) & dd$event_death %in% 0:1
        event_col <- "event_death"; time_col <- "time_death"
    }
    design <- ckm_a_design(dd, prep, ckm_a_covariates(m, a), rows, a$stage_column,
                           time_col = time_col, event_col = event_col, logistic = logistic)
    z <- design$data; out$n <- nrow(z); out$events <- sum(z$.e)
    out$reference <- prep$reference; out$dropped_constant_covariates <- design$dropped_constant_covariates
    if (nrow(z) < a$min_n || sum(z$.e) < a$min_events || (logistic && sum(z$.e == 0) < a$min_events)) return(out)
    xnames <- colnames(prep$x)
    if (any(vapply(z[xnames], function(x) length(unique(x)) < 2L, logical(1)))) {out$status <- "constant_or_absent_exposure_level"; return(out)}
    fit <- ckm_a_fit(ckm_a_formula(xnames, design$covariates, stratify = FALSE, logistic = logistic), z, logistic = logistic)
    out$warning <- fit$warning
    if (!fit$ok) {out$status <- fit$status; return(out)}
    C <- matrix(0, nrow = length(xnames), ncol = length(fit$beta), dimnames = list(NULL, names(fit$beta)))
    C[cbind(seq_along(xnames), match(xnames, names(fit$beta)))] <- 1
    es <- ckm_a_effects(fit, C, prep$terms)
    out[, c("beta", "se", "effect", "ci_low", "ci_high", "p")] <- es[, c("beta", "se", "effect", "ci_low", "ci_high", "p")]
    out$p_global <- ckm_a_wald(fit$beta[xnames], fit$V[xnames, xnames, drop = FALSE]); out$status <- "ok"
    out
}

ckm_a_adjust <- function(tab, p, out, groups, unique_keys = NULL) {
    tab <- data.table::as.data.table(tab)
    if (!nrow(tab)) {tab[[out]] <- numeric(); return(tab)}
    tab[[out]] <- NA_real_
    ncol <- paste0(out, "_planned_tests"); tab[[ncol]] <- NA_integer_
    key <- do.call(paste, c(as.data.frame(tab[, ..groups]), sep = "\r"))
    for (g in unique(key)) {
        ii <- which(key == g)
        if (is.null(unique_keys)) {
            pval <- tab[[p]][ii]
            # Failed/planned tests retain their place in the correction family.
            qval <- stats::p.adjust(pval, method = "BH", n = length(pval))
            tab[[out]][ii] <- qval; tab[[ncol]][ii] <- length(pval)
        } else {
            id <- do.call(paste, c(as.data.frame(tab[ii, ..unique_keys]), sep = "\r"))
            first <- !duplicated(id); pval <- tab[[p]][ii[first]]
            qval <- stats::p.adjust(pval, method = "BH", n = sum(first))
            tab[[out]][ii] <- qval[match(id, id[first])]; tab[[ncol]][ii] <- sum(first)
        }
    }
    tab
}

ckm_a_classify <- function(assoc, stages, a) {
    assoc <- data.table::as.data.table(assoc)
    assoc[, `:=`(pattern = "unresolved", direction = "unresolved", interpretation = "observational_association")]
    for (i in seq_len(nrow(assoc))) {
        r <- assoc[i]
        ss <- stages[feature == r$feature & source == r$source & endpoint == r$endpoint & lag_years == r$lag_years]
        if (r$status != "ok") {assoc$pattern[i] <- "not_estimable"; next}
        main <- is.finite(r$q_shared) && r$q_shared < a$fdr
        interaction <- is.finite(r$q_interaction) && r$q_interaction < a$fdr
        equivalence <- is.finite(r$q_equivalence) && r$q_equivalence < a$fdr
        all_stages <- nrow(ss) >= length(a$stages) && all(ss$status == "ok") && all(is.finite(ss$beta))
        scalar <- r$type %in% c("continuous", "binary") || length(unique(ss$term)) == 1L
        same <- all_stages && scalar && (all(ss$beta > 0) || all(ss$beta < 0))
        if (scalar && is.finite(r$beta)) assoc$direction[i] <- if (r$beta > 0) "risk_associated" else "protective_associated"
        if (main && equivalence && same) {
            # A detectable small difference may still be equivalent within the
            # prespecified margin. Keep q_interaction alongside this label.
            assoc$pattern[i] <- "shared_equivalent_within_margin"
        } else if (interaction) {
            sig <- ss[is.finite(q_stage) & q_stage < a$fdr & is.finite(p) & p < a$fdr]
            positive <- nrow(sig) && any(sig$beta > 0 & sig$ci_low > 1)
            negative <- nrow(sig) && any(sig$beta < 0 & sig$ci_high < 1)
            assoc$pattern[i] <- if (scalar && positive && negative) "stage_dependent_opposite_supported" else "stage_dependent"
        } else if (main && same && is.finite(r$q_interaction)) assoc$pattern[i] <- "consistent_no_detected_heterogeneity" else
            if (main) assoc$pattern[i] <- "pooled_association_stage_pattern_uncertain" else assoc$pattern[i] <- "no_FDR_supported_association"
        if (!scalar) assoc$direction[i] <- "categorical_contrasts_see_stage_table"
    }
    assoc
}

ckm_run_associations <- function(d, cfg, exposures = NULL) {
    if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("survival", quietly = TRUE)) stop("Install data.table and survival before associations")
    a <- ckm_a_options(cfg)
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    required <- c("eid", a$stage_column, "time", "status")
    if (length(setdiff(required, names(d)))) stop("CKM data lack: ", paste(setdiff(required, names(d)), collapse = ", "))
    d$eid <- as.character(d$eid)
    if (anyNA(d$eid) || any(!nzchar(d$eid)) || anyDuplicated(d$eid)) stop("CKM data require one nonmissing row per eid")
    d$event <- as.integer(d$status == 1L)
    if (any(!is.na(d$status) & !d$status %in% 0:2)) stop("status must be 0 (censored), 1 (stage4), or 2 (competing death)")
    if (is.null(exposures)) {
        if (is.null(cfg$exposure_manifest)) stop("Pass exposures=list(manifest=...,data=...) or cfg$exposure_manifest")
        exposures <- cfg$exposure_manifest
    }
    if (is.data.frame(exposures)) {
        manifest <- ckm_a_manifest(exposures); xdata <- d
    } else if (is.list(exposures) && all(c("manifest", "data") %in% names(exposures))) {
        manifest <- ckm_a_manifest(exposures$manifest); xdata <- as.data.frame(exposures$data, stringsAsFactors = FALSE)
    } else stop("exposures must be a manifest data.frame or list(manifest,data)")
    if (!"eid" %in% names(xdata)) stop("Exposure data lack eid")
    xdata$eid <- as.character(xdata$eid)
    if (anyNA(xdata$eid) || any(!nzchar(xdata$eid)) || anyDuplicated(xdata$eid)) stop("Exposure data require one nonmissing row per eid")
    match_index <- match(d$eid, xdata$eid)
    all_covariates <- unique(unlist(lapply(seq_len(nrow(manifest)), function(i) ckm_a_covariates(manifest[i, ], a))))
    missing_cov <- setdiff(all_covariates, names(d))
    if (length(missing_cov)) stop("Configured adjustment columns are missing: ", paste(missing_cov, collapse = ", "), ". Correct the data/configuration; covariates are never silently removed.")
    for (v in all_covariates) if (all(is.na(d[[v]]))) stop("Configured adjustment column is entirely missing: ", v)
    scale_rows <- !is.na(d[[a$stage_column]]) & d[[a$stage_column]] %in% a$stages
    worker <- function(i) {
        m <- manifest[i, , drop = FALSE]
        if (!m$column %in% names(xdata)) {
            prep <- list(ok = FALSE, status = "missing_exposure_column", terms = "not_estimable", center = NA_real_, scale = NA_real_)
        } else prep <- ckm_a_prepare_exposure(xdata[[m$column]][match_index], m, scale_rows)
        inc <- lapply(a$lag_years, function(lag) ckm_a_incident(d, m, prep, a, lag))
        yang <- list()
        if (isTRUE(a$yang)) yang[[length(yang) + 1L]] <- ckm_a_yang(d, m, prep, a, "baseline_stage4_marker")
        if (isTRUE(a$yang_death)) yang[[length(yang) + 1L]] <- ckm_a_yang(d, m, prep, a, "stage4_all_cause_mortality")
        list(summary = ckm_a_bind(lapply(inc, `[[`, "summary")), pooled = ckm_a_bind(lapply(inc, `[[`, "pooled")),
             stage = ckm_a_bind(lapply(inc, `[[`, "stage")), contrasts = ckm_a_bind(lapply(inc, `[[`, "contrasts")), yang = ckm_a_bind(yang))
    }
    # Fork workers share read-only matrices on Unix; Windows uses a safe serial loop.
    if (a$workers > 1L && .Platform$OS.type != "windows") {
        values <- parallel::mclapply(seq_len(nrow(manifest)), worker, mc.cores = min(a$workers, nrow(manifest)), mc.preschedule = TRUE, mc.set.seed = FALSE)
    } else values <- lapply(seq_len(nrow(manifest)), worker)
    failures <- vapply(values, inherits, logical(1), what = "try-error")
    if (any(failures)) stop("Association worker failed; no incomplete result family was published: ", paste(values[failures], collapse = "; "))
    summary <- ckm_a_bind(lapply(values, `[[`, "summary")); stage <- ckm_a_bind(lapply(values, `[[`, "stage"))
    pooled <- ckm_a_bind(lapply(values, `[[`, "pooled")); contrasts <- ckm_a_bind(lapply(values, `[[`, "contrasts"))
    yang <- ckm_a_bind(lapply(values, `[[`, "yang"))
    groups <- c("source", "endpoint", "lag_years")
    for (test in c("shared", "interaction", "equivalence", "ph")) summary <- ckm_a_adjust(summary, paste0("p_", test), paste0("q_", test), groups)
    stage <- ckm_a_adjust(stage, "p_stage", "q_stage", c(groups, "stage"), unique_keys = "feature")
    pooled <- ckm_a_adjust(pooled, "p", "q_contrast", groups)
    contrasts <- ckm_a_adjust(contrasts, "p_difference", "q_difference", groups)
    if (nrow(yang)) yang <- ckm_a_adjust(yang, "p_global", "q_global", c("source", "endpoint"), unique_keys = "feature")
    summary <- ckm_a_classify(summary, stage, a)
    out <- list(assoc = summary[lag_years == 0], pooled = pooled[lag_years == 0], stage = stage[lag_years == 0],
                contrasts = contrasts[lag_years == 0], sensitivity = summary[lag_years > 0],
                sensitivity_stage = stage[lag_years > 0], yang = yang, settings = a)
    if (!is.null(cfg$out_dir) && nzchar(cfg$out_dir)) {
        dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
        for (name in setdiff(names(out), "settings")) {
            tab <- out[[name]]
            if (is.data.frame(tab) && ncol(tab)) data.table::fwrite(tab, file.path(cfg$out_dir, paste0(name, ".tsv")), sep = "\t", na = "NA")
        }
    }
    out
}
