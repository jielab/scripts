# C4 connection helpers. No UKB I/O, disease fitting, or package startup here.
# Functions are deliberately usable in small synthetic tests with base R.
C4_CONNECTION_VERSION <- "20261118-canonical-outcome-v2"

c4_core_names <- function(x) {
  x <- as.character(x)
  unique(x[!is.na(x) & nzchar(x)])
}

c4_core_env_number <- function(name, default, lower = 0, upper = Inf) {
  z <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(z) != 1L || !is.finite(z) || z < lower || z > upper)
    stop("Invalid ", name)
  z
}

c4_concept_lambda <- function() {
  # Old C4_CONCEPT_RIDGE penalized the unnormalized sum of squares. Reusing
  # its typical value 10 as a normalized penalty would silently change the
  # model by orders of magnitude. Require an explicit new-scale choice.
  if (nzchar(Sys.getenv("C4_CONCEPT_RIDGE", "")) && !nzchar(Sys.getenv("C4_CONCEPT_LAMBDA", "")))
    stop("C4_CONCEPT_RIDGE is retired (old sum-of-squares penalty). Unset it to use the new normalized default 0.001, or explicitly set C4_CONCEPT_LAMBDA. The numerical values are not interchangeable.")
  c4_core_env_number("C4_CONCEPT_LAMBDA", .001, 0)
}

c4_core_group_cap <- function(groups, max_n, seed) {
  groups <- as.character(groups)
  if (anyNA(groups) || any(!nzchar(groups))) stop("Missing C4 selection group")
  if (length(groups) <= max_n) return(seq_along(groups))
  block <- split(seq_along(groups), groups)
  set.seed(seed)
  order <- sample(sort(names(block)))
  sizes <- lengths(block[order])
  chosen <- order[cumsum(sizes) <= max_n]
  if (!length(chosen)) stop("Selection cap is smaller than a complete family")
  sort(unlist(block[chosen], use.names = FALSE))
}

c4_core_covariates <- function(dat, covars) {
  # A full-rank basic-covariate matrix. No phenotype target enters this matrix.
  ans <- matrix(1, nrow(dat), 1L, dimnames = list(NULL, "Intercept"))
  for (nm in c4_core_names(covars)) {
    x <- dat[[nm]]
    if (is.numeric(x)) {
      if (any(!is.finite(x))) stop("Non-finite complete-case covariate: ", nm)
      if (stats::sd(x) > 0) ans <- cbind(ans, setNames(data.frame(x), nm))
    } else {
      x <- factor(x)
      if (nlevels(x) > 1L) {
        z <- stats::model.matrix(~ x)[, -1L, drop = FALSE]
        colnames(z) <- paste0(nm, "__", colnames(z))
        ans <- cbind(ans, z)
      }
    }
  }
  ans <- as.matrix(ans)
  q <- qr(ans)
  ans[, q$pivot[seq_len(q$rank)], drop = FALSE]
}

c4_joint_candidates <- function(membership, eligible, components, max_candidates = 600L) {
  eligible <- c4_core_names(eligible)
  components <- c4_core_names(components)
  if (!length(eligible)) return(character())
  max_candidates <- as.integer(max_candidates)
  if (max_candidates < 1L) stop("Invalid joint candidate cap")
  if (!is.data.frame(membership) || !all(c("feature", "component", "r1", "r2", "selected") %in% names(membership)))
    stop("YSjoint requires training replicated feature-domain membership")
  m <- membership[membership$feature %in% eligible & membership$component %in% components &
                    membership$selected %in% TRUE & is.finite(membership$r1) & is.finite(membership$r2), , drop = FALSE]
  m$strength <- pmin(abs(m$r1), abs(m$r2))
  m <- m[order(-m$strength, m$feature, m$component), , drop = FALSE]
  # Round-robin candidate admission protects small domains before joint ranking.
  # The final panel is selected by incremental multivariate reconstruction gain.
  queues <- lapply(components, function(cmp) c4_core_names(m$feature[m$component == cmp]))
  names(queues) <- components
  out <- character()
  while (length(out) < max_candidates && any(lengths(queues) > 0L)) {
    for (cmp in components) {
      queues[[cmp]] <- setdiff(queues[[cmp]], out)
      if (length(queues[[cmp]]) && length(out) < max_candidates) {
        out <- c(out, queues[[cmp]][1L])
        queues[[cmp]] <- queues[[cmp]][-1L]
      }
    }
  }
  out
}

c4_joint_rank <- function(train, candidates, components, covars, max_budget,
                          max_n = 30000L, seed = 2026L, variance_tolerance = 1e-10) {
  candidates <- c4_core_names(candidates)
  components <- c4_core_names(components)
  covars <- c4_core_names(covars)
  max_budget <- as.integer(max_budget)
  if (max_budget < 1L || max_n < 100L) stop("Invalid YSjoint budget/sample cap")
  if (length(intersect(candidates, c(components, covars))))
    stop("YSjoint assay, target and covariate namespaces overlap")
  required <- unique(c(candidates, components, covars))
  if (length(setdiff(required, names(train)))) stop("YSjoint input columns missing")
  if (!all(vapply(train[c(candidates, components)], is.numeric, logical(1))))
    stop("YSjoint assays and LE8 domain points must be numeric")
  # Selection uses only baseline assays, LE8 labels and declared basic covariates.
  # Disease endpoints may exist in train but are never read by this function.
  complete <- stats::complete.cases(train[c(components, covars)])
  for (nm in c(components, covars)) if (is.numeric(train[[nm]]))
    complete <- complete & is.finite(train[[nm]])
  rows <- which(complete)
  empty <- list(features = character(), path = data.frame(), domain_gain = data.frame(),
                feature_audit = data.frame(), status = "insufficient common LE8 labels",
                N_complete = length(rows), N_selection = 0L, components = components,
                objective = "equal-domain mean incremental LE8 R2 conditional on basic covariates")
  if (length(rows) < 100L || !length(candidates)) return(empty)
  groups <- if (".group" %in% names(train)) train$.group else if ("eid" %in% names(train)) train$eid else seq_len(nrow(train))
  rows <- rows[c4_core_group_cap(groups[rows], max_n, seed)]
  d <- train[rows, required, drop = FALSE]
  x <- as.matrix(d[candidates])
  center <- scale <- rep(NA_real_, ncol(x))
  for (j in seq_len(ncol(x))) {
    ok <- is.finite(x[, j])
    if (!any(ok)) next
    center[j] <- stats::median(x[ok, j])
    x[!ok, j] <- center[j]
    scale[j] <- stats::sd(x[, j])
  }
  usable <- is.finite(scale) & scale > 0
  audit <- data.frame(feature = candidates, median = center, sd = scale, usable = usable,
                      stringsAsFactors = FALSE)
  x <- x[, usable, drop = FALSE]
  if (!ncol(x)) { empty$status <- "no varying assays"; return(empty) }
  x <- sweep(sweep(x, 2L, center[usable], "-"), 2L, scale[usable], "/")
  y <- as.matrix(d[components])
  ysd <- apply(y, 2L, stats::sd)
  good_y <- is.finite(ysd) & ysd > 0
  if (!any(good_y)) { empty$status <- "no varying LE8 domains"; return(empty) }
  y <- sweep(sweep(y[, good_y, drop = FALSE], 2L, colMeans(y[, good_y, drop = FALSE]), "-"), 2L, ysd[good_y], "/")
  q <- qr(c4_core_covariates(d, covars))
  x <- qr.resid(q, x)
  y <- qr.resid(q, y)
  # One Gram matrix keeps the greedy joint procedure feasible for proteomics.
  # A rank-one Schur update removes what prior assays already explain.
  gram <- crossprod(x) / (nrow(x) - 1)
  cross <- crossprod(x, y) / (nrow(x) - 1)
  weights <- rep(1 / ncol(y), ncol(y))
  remaining <- colnames(x)
  selected <- character()
  steps <- domain <- list()
  cumulative <- rep(0, ncol(y))
  for (step in seq_len(min(max_budget, ncol(x)))) {
    v <- diag(gram)
    gains <- rowSums(sweep(cross * cross, 2L, weights, "*")) / v
    gains[!is.finite(gains) | v <= variance_tolerance | remaining %in% selected] <- -Inf
    order <- order(-gains, remaining)
    j <- order[1L]
    if (!is.finite(gains[j])) break
    individual <- as.numeric(cross[j, ] ^ 2 / v[j])
    cumulative <- cumulative + individual
    selected <- c(selected, remaining[j])
    steps[[step]] <- data.frame(step = step, feature = remaining[j], mean_delta_R2 = gains[j],
                               mean_cumulative_R2 = sum(weights * cumulative), residual_assay_variance = v[j])
    domain[[step]] <- data.frame(step = step, feature = remaining[j], component = colnames(y),
                                delta_R2 = individual, cumulative_R2 = cumulative, domain_weight = weights)
    g <- gram[, j]
    cc <- cross[j, ]
    gram <- gram - tcrossprod(g) / v[j]
    cross <- cross - outer(g, cc) / v[j]
    gram <- (gram + t(gram)) / 2
  }
  list(features = selected, path = if (length(steps)) do.call(rbind, steps) else data.frame(),
       domain_gain = if (length(domain)) do.call(rbind, domain) else data.frame(), feature_audit = audit,
       status = if (length(selected) >= max_budget) "ok" else "insufficient independent eligible assays",
       N_complete = sum(complete), N_selection = length(rows), families_selection = length(unique(groups[rows])),
       selected_ids = if ("eid" %in% names(train)) as.character(train$eid[rows]) else as.character(rows),
       components = colnames(y), requested_components = components,
       objective = "equal-domain mean incremental LE8 R2 conditional on basic covariates",
       interpretation = "training selection objective; held-out reconstruction is evaluated separately")
}

c4_rescue_panel <- function(ranked, supervised, k, fraction = 0.8,
                            exclusion = c("selected_panel", "entire_eligible_pool")) {
  exclusion <- match.arg(exclusion)
  ranked <- c4_core_names(ranked)
  supervised <- c4_core_names(supervised)
  n <- ceiling(k * fraction)
  if (length(supervised) < n) return(character())
  chosen <- head(supervised, n)
  forbidden <- if (exclusion == "selected_panel") chosen else supervised
  out <- unique(c(chosen, head(setdiff(ranked, forbidden), k - n)))
  if (length(out) == k) out else character()
}

c4_centered_ridge <- function(x, y, newx, lambda = 0.001) {
  x <- as.matrix(x); newx <- as.matrix(newx); y <- as.numeric(y)
  if (!identical(colnames(x), colnames(newx)) || ncol(x) < 1L || nrow(x) != length(y))
    stop("Invalid ridge prediction matrices")
  if (any(!is.finite(x)) || any(!is.finite(newx)) || any(!is.finite(y)) || !is.finite(lambda) || lambda < 0)
    stop("Non-finite ridge input")
  xmean <- colMeans(x); ymean <- mean(y)
  xc <- sweep(x, 2L, xmean, "-")
  target <- y - ymean
  # Mean squared loss + lambda * sum(beta^2); penalty does not change with donor N.
  gram <- crossprod(xc) / nrow(xc)
  rhs <- crossprod(xc, target) / nrow(xc)
  beta <- if (lambda > 0) solve(gram + diag(lambda, ncol(xc)), rhs) else qr.solve(xc, target)
  beta <- as.numeric(beta)
  intercept <- ymean - sum(xmean * beta)
  list(prediction = intercept + drop(newx %*% beta), train_prediction = intercept + drop(x %*% beta),
       beta = beta, intercept = intercept, xmean = xmean, ymean = ymean, lambda = lambda,
       loss = "mean squared error plus lambda times squared coefficients; unpenalized intercept")
}

c4_oof_support <- function(observed, predicted, reference, min_n = 100L, min_r2 = 0) {
  observed <- as.matrix(observed); predicted <- as.matrix(predicted); reference <- as.matrix(reference)
  if (!identical(dim(observed), dim(predicted)) || !identical(dim(observed), dim(reference)))
    stop("OOF support matrices do not align")
  rows <- lapply(seq_len(ncol(observed)), function(j) {
    ok <- is.finite(observed[, j]) & is.finite(predicted[, j]) & is.finite(reference[, j])
    n <- sum(ok); den <- sum((observed[ok, j] - reference[ok, j]) ^ 2)
    r2 <- if (n > 0 && den > 0) 1 - sum((observed[ok, j] - predicted[ok, j]) ^ 2) / den else NA_real_
    supported <- n >= min_n && is.finite(r2) && r2 > min_r2
    data.frame(component = colnames(observed)[j], N_OOF = n, OOF_R2 = r2, min_R2 = min_r2,
               supported = supported, support_source = "training OOF predictions only; held-out labels unused")
  })
  do.call(rbind, rows)
}

c4_resolve_replacement <- function(requested, available_components, support, measure_map,
                                   basic, measured, prs = character()) {
  requested <- c4_core_names(requested)
  if (any(!requested %in% unique(measure_map$component))) stop("Unknown replacement domain")
  supported <- support$component[support$supported %in% TRUE]
  replaced <- intersect(requested, intersect(available_components, supported))
  unsupported <- setdiff(requested, replaced)
  retained <- setdiff(available_components, replaced)
  forbidden <- unique(c(measure_map$variable[measure_map$component %in% replaced],
                        if (length(replaced)) c("le8", "le8.pts", "le8_score", "le4", "le4.pts", "le4_score")))
  list(replaced = replaced, retained = retained, unsupported_requested = unsupported,
       background = setdiff(unique(c(basic, measured, prs)), forbidden), forbidden = forbidden,
       status = if (length(replaced)) "ok" else "unavailable: no requested domain has training OOF support")
}
