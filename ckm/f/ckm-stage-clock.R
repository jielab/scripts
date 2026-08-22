#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival and regression helpers
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cox_one <- function(dat, x, covars = character(), lag = 0, min_n = 300, min_event = 15) {
  need <- unique(c("time", "event", x, covars))
  if (!all(need %in% names(dat))) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = 0, events = 0))
  dd <- dat[, need, drop = FALSE]
  tt <- safe_num(dd$time); ee <- as.integer(safe_num(dd$event) == 1)
  keep <- is.finite(tt) & tt > lag & !(ee == 1 & tt <= lag)
  dd <- dd[keep, , drop = FALSE]
  dd$time_lag <- safe_num(dd$time) - lag
  dd$event_lag <- as.integer(safe_num(dd$event) == 1)
  dd <- dd[complete.cases(dd[, unique(c("time_lag", "event_lag", x, covars)), drop = FALSE]), , drop = FALSE]
  if (nrow(dd) < min_n || sum(dd$event_lag) < min_event || sd(safe_num(dd[[x]]), na.rm = TRUE) == 0) {
    return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(dd), events = sum(dd$event_lag)))
  }
  fm <- as.formula(paste0("Surv(time_lag, event_lag) ~ ", paste(bt(c(x, covars)), collapse = " + ")))
  fit <- try(coxph(fm, data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, HR = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(dd), events = sum(dd$event_lag)))
  sm <- summary(fit)
  co <- sm$coef[1, ]; ci <- sm$conf.int[1, ]
  data.table(beta = co[["coef"]], se = co[["se(coef)"]], p = co[["Pr(>|z|)"]],
             HR = ci[["exp(coef)"]], lo = ci[["lower .95"]], hi = ci[["upper .95"]],
             n = nrow(dd), events = sum(dd$event_lag))
}

lm_one <- function(dat, y, x, covars = character(), min_n = 200) {
  need <- unique(c(y, x, covars))
  if (!all(need %in% names(dat))) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = 0))
  dd <- dat[, need, drop = FALSE]
  dd <- dd[complete.cases(dd), , drop = FALSE]
  if (nrow(dd) < min_n || length(unique(dd[[x]])) < 2) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd)))
  fit <- try(lm(as.formula(paste0(bt(y), " ~ ", paste(bt(c(x, covars)), collapse = " + "))), data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd)))
  co <- summary(fit)$coef
  rn <- rownames(co)
  hit <- which(rn == x | rn == bt(x))
  if (!length(hit)) return(data.table(beta = NA_real_, se = NA_real_, p = NA_real_, n = nrow(dd)))
  data.table(beta = co[hit[1], "Estimate"], se = co[hit[1], "Std. Error"], p = co[hit[1], "Pr(>|t|)"], n = nrow(dd))
}

calc_cindex <- function(time, event, lp) {
  ok <- is.finite(time) & is.finite(event) & is.finite(lp) & time > 0
  if (sum(ok) < 20 || sum(event[ok] == 1) < 2) return(NA_real_)
  fit <- try(concordance(Surv(time[ok], event[ok]) ~ lp[ok], reverse = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) NA_real_ else as.numeric(fit$concordance)
}

calc_auc <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]; score <- score[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(score)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_calibration_slope <- function(time, event, lp) {
  ok <- is.finite(time) & time > 0 & is.finite(event) & is.finite(lp)
  if (sum(ok) < 100 || sum(event[ok] == 1) < 10 || sd(lp[ok]) <= 1e-8) return(NA_real_)
  fit <- try(coxph(Surv(time[ok], event[ok]) ~ lp[ok]), silent = TRUE)
  if (inherits(fit, "try-error")) NA_real_ else as.numeric(coef(fit)[1])
}

make_folds <- function(strata, k, seed0 = seed) {
  set.seed(seed0)
  strata <- as.character(strata)
  strata[is.na(strata) | !nzchar(strata)] <- "Missing"
  fold <- integer(length(strata))
  for (s0 in sort(unique(strata))) {
    ii <- which(strata == s0)
    fold[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  fold
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Quantitative CKM stage clock
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

survfit_rmst <- function(sf, tau) {
  if (!is.finite(tau) || tau <= 0) return(NA_real_)
  keep <- is.finite(sf$time) & sf$time < tau
  tt <- c(0, sf$time[keep], tau)
  ss <- c(1, sf$surv[keep])
  if (length(ss) < length(tt)) ss <- c(ss, tail(ss, 1))
  sum(diff(tt) * head(ss, -1), na.rm = TRUE)
}

pava_monotone <- function(y, decreasing = TRUE, weights = NULL) {
  y <- safe_num(y)
  if (is.null(weights)) weights <- rep(1, length(y))
  weights <- safe_num(weights)
  weights[!is.finite(weights) | weights <= 0] <- 1
  if (!length(y)) return(y)
  z <- if (decreasing) -y else y

  means <- numeric()
  wts <- numeric()
  starts <- integer()
  ends <- integer()

  for (i in seq_along(z)) {
    means <- c(means, z[i])
    wts <- c(wts, weights[i])
    starts <- c(starts, i)
    ends <- c(ends, i)

    while (length(means) >= 2) {
      k <- length(means)
      if (!is.finite(means[k - 1]) || !is.finite(means[k]) || means[k - 1] <= means[k]) break
      new_w <- wts[k - 1] + wts[k]
      new_m <- (means[k - 1] * wts[k - 1] + means[k] * wts[k]) / new_w
      means[k - 1] <- new_m
      wts[k - 1] <- new_w
      ends[k - 1] <- ends[k]
      means <- means[-k]
      wts <- wts[-k]
      starts <- starts[-k]
      ends <- ends[-k]
    }
  }

  out <- rep(NA_real_, length(y))
  for (k in seq_along(means)) out[starts[k]:ends[k]] <- means[k]
  if (decreasing) -out else out
}

monotone_remaining <- function(x, tau = Inf) {
  x <- safe_num(x)
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(rep(NA_real_, length(x)))

  # Fill isolated missing values before the monotonic projection.
  for (i in seq_along(x)) if (is.na(x[i]) && i > 1) x[i] <- x[i - 1]
  for (i in rev(seq_along(x))) if (is.na(x[i]) && i < length(x)) x[i] <- x[i + 1]

  out <- pava_monotone(x, decreasing = TRUE)
  out <- pmax(out, 0)
  if (is.finite(tau)) out <- pmin(out, tau)
  out
}

weibull_rmst_eta <- function(eta, scale, tau) {
  eta <- safe_num(eta)
  out <- rep(NA_real_, length(eta))
  if (!is.finite(scale) || scale <= 0 || !is.finite(tau) || tau <= 0) return(out)

  ok <- is.finite(eta)
  if (!any(ok)) return(out)
  e <- eta[ok]
  log_u <- (log(tau) - e) / scale

  # When the predicted time scale is far beyond tau, RMST is numerically tau.
  far_right <- log_u < -25
  ans <- rep(NA_real_, length(e))
  ans[far_right] <- tau

  mid <- !far_right
  if (any(mid)) {
    u <- exp(pmin(log_u[mid], 700))
    log_p <- pgamma(u, shape = scale, lower.tail = TRUE, log.p = TRUE)
    log_rmst <- e[mid] + lgamma(1 + scale) + log_p
    ans[mid] <- exp(pmin(log_rmst, log(tau)))
  }

  ans[!is.finite(ans)] <- tau
  out[ok] <- pmin(pmax(ans, 0), tau)
  out
}

make_clock_reference <- function(dd, covars, n_ref = clock_ref_n, mode = clock_reference, seed_offset = 0L) {
  idx <- seq_len(nrow(dd))
  set.seed(seed + 917L + as.integer(seed_offset))
  if (length(idx) > n_ref) idx <- sample(idx, n_ref)
  ref <- dd[idx, unique(c("stage_f", covars)), drop = FALSE]

  if (mode == "healthy_le8") {
    healthy <- dd$ckm.stage == 0
    if (sum(healthy, na.rm = TRUE) < 500) healthy <- dd$ckm.stage %in% c(0, 1)
    for (v in intersect(clock_behavior_vars, names(ref))) {
      if (is.numeric(dd[[v]]) || is.integer(dd[[v]])) {
        q <- quantile(safe_num(dd[[v]][healthy]), 0.80, na.rm = TRUE, names = FALSE, type = 8)
        if (is.finite(q)) ref[[v]] <- q
      }
    }
  }
  ref
}

survreg_clock <- function(dat, covars, draws = clock_draws, tau = clock_tau, seed_offset = 0L) {
  dd <- dat[dat$yin & dat$ckm.stage %in% 0:3, , drop = FALSE]
  dd$stage_f <- factor(dd$ckm.stage, levels = 0:3)
  need <- unique(c("time", "event", "ckm.stage", "stage_f", covars))
  dd <- dd[
    complete.cases(dd[, need, drop = FALSE]) &
      is.finite(dd$time) & dd$time > 0,
    ,
    drop = FALSE
  ]

  if (nrow(dd) < 1000 || sum(dd$event) < 50) {
    stop("Too few stage 0-3 participants/events for the quantitative stage clock", call. = FALSE)
  }

  # Do not extrapolate the clock beyond adequately supported follow-up.
  support95 <- as.numeric(quantile(dd$time, 0.95, na.rm = TRUE, names = FALSE, type = 8))
  max_follow <- max(dd$time, na.rm = TRUE)
  tau_eff <- min(tau, support95, max_follow)
  if (!is.finite(tau_eff) || tau_eff <= 1) {
    stop("The effective RMST horizon is too short: ", tau_eff, " years", call. = FALSE)
  }

  fm <- as.formula(
    paste0(
      "Surv(time, event) ~ stage_f",
      if (length(covars)) paste0(" + ", paste(bt(covars), collapse = " + ")) else ""
    )
  )

  # x/model/y are deliberately not retained.  The previous x=TRUE fit objects,
  # repeated for crude/demog/le8 and then xz-compressed, caused a large memory
  # spike immediately after Fig1 was written.
  fit <- survreg(
    fm,
    data = dd,
    dist = "weibull",
    x = FALSE,
    y = FALSE,
    model = FALSE
  )

  ref0 <- make_clock_reference(dd, covars, n_ref = clock_ref_n, seed_offset = seed_offset)
  refs_df <- do.call(rbind, lapply(0:3, function(s) {
    z <- ref0
    z$stage_f <- factor(s, levels = 0:3)
    z$stage <- s
    z
  }))
  rownames(refs_df) <- NULL

  beta <- coef(fit)
  beta_names <- names(beta)
  trm <- delete.response(terms(fit))
  X <- model.matrix(trm, refs_df)
  X <- X[, beta_names, drop = FALSE]
  eta <- as.numeric(X %*% beta)
  rmst_i <- weibull_rmst_eta(eta, fit$scale, tau_eff)

  point <- data.table(stage = refs_df$stage, rmst = rmst_i)[
    ,
    .(remaining_years_raw = mean(rmst, na.rm = TRUE)),
    by = stage
  ][order(stage)]

  stage_counts <- as.data.table(dd)[
    ,
    .(clock_N = .N, clock_events = sum(event, na.rm = TRUE)),
    by = .(stage = ckm.stage)
  ]
  point <- merge(point, stage_counts, by = "stage", all.x = TRUE, sort = FALSE)
  point <- point[order(stage)]
  point[, remaining_years := monotone_remaining(remaining_years_raw, tau = tau_eff)]
  point[, `:=`(
    tau_requested = tau,
    tau_effective = tau_eff,
    lo = NA_real_,
    hi = NA_real_
  )]

  sim <- NULL
  if (draws > 1) {
    vv <- vcov(fit)
    vv_names <- rownames(vv)
    scale_name <- setdiff(vv_names, beta_names)
    scale_name <- scale_name[grepl("scale", scale_name, ignore.case = TRUE)][1]

    if (!is.na(scale_name) && nzchar(scale_name)) {
      mu <- c(beta, setNames(log(fit$scale), scale_name))
      vv_use <- vv[names(mu), names(mu), drop = FALSE]
    } else {
      mu <- beta
      vv_use <- vv[beta_names, beta_names, drop = FALSE]
    }

    bd <- try(MASS::mvrnorm(draws, mu = mu, Sigma = vv_use), silent = TRUE)
    if (!inherits(bd, "try-error")) {
      bd <- as.matrix(bd)
      colnames(bd) <- names(mu)

      set.seed(seed + 1907L + as.integer(seed_offset))
      draw_idx <- seq_len(nrow(ref0))
      if (length(draw_idx) > clock_draw_ref_n) draw_idx <- sample(draw_idx, clock_draw_ref_n)
      ref_draw0 <- ref0[draw_idx, , drop = FALSE]
      refs_draw <- do.call(rbind, lapply(0:3, function(s) {
        z <- ref_draw0
        z$stage_f <- factor(s, levels = 0:3)
        z$stage <- s
        z
      }))
      rownames(refs_draw) <- NULL
      Xd <- model.matrix(trm, refs_draw)
      Xd <- Xd[, beta_names, drop = FALSE]

      sim <- matrix(NA_real_, nrow = draws, ncol = 4)
      colnames(sim) <- as.character(0:3)
      progress_log(
        "stage clock: simulating coefficient draws=", fmt_int(draws),
        "; draw-reference N=", fmt_int(nrow(ref_draw0)),
        "; RMST horizon=", sprintf("%.2f", tau_eff), " y"
      )

      for (b in seq_len(draws)) {
        beta_b <- bd[b, beta_names]
        scale_b <- if (!is.na(scale_name) && scale_name %in% colnames(bd)) {
          exp(bd[b, scale_name])
        } else {
          fit$scale
        }
        eta_b <- as.numeric(Xd %*% beta_b)
        rmst_b <- weibull_rmst_eta(eta_b, scale_b, tau_eff)
        tmp <- data.table(stage = refs_draw$stage, rmst = rmst_b)[
          ,
          mean(rmst, na.rm = TRUE),
          by = stage
        ][order(stage)]$V1
        sim[b, ] <- monotone_remaining(tmp, tau = tau_eff)
      }

      ci <- rbindlist(lapply(0:3, function(s) {
        z <- sim[, as.character(s)]
        data.table(
          stage = s,
          lo = quantile(z, 0.025, na.rm = TRUE, names = FALSE),
          hi = quantile(z, 0.975, na.rm = TRUE, names = FALSE)
        )
      }))
      point[, c("lo", "hi") := NULL]
      point <- merge(point, ci, by = "stage", all.x = TRUE, sort = FALSE)
      point <- point[order(stage)]
    } else {
      warning("Stage-clock coefficient simulation failed: ", as.character(bd), call. = FALSE)
    }
  }

  vv <- vcov(fit)
  se <- sqrt(diag(vv))[beta_names]
  fit_summary <- data.table(
    term = beta_names,
    estimate = as.numeric(beta),
    se = as.numeric(se),
    z = as.numeric(beta / se),
    p = 2 * pnorm(-abs(beta / se))
  )
  fit_meta <- data.table(
    item = c(
      "method", "distribution", "estimand", "clock N", "events",
      "requested tau", "effective tau", "95th percentile follow-up",
      "maximum follow-up", "reference N", "draw-reference N", "draws"
    ),
    value = c(
      "standardized Weibull AFT restricted-mean clock",
      "Weibull",
      "population-equivalent RMST distance to incident CKM stage 4",
      nrow(dd), sum(dd$event), tau, tau_eff, support95,
      max_follow, nrow(ref0), min(nrow(ref0), clock_draw_ref_n), draws
    )
  )

  list(
    point = point,
    draws = sim,
    fit_summary = fit_summary,
    fit_meta = fit_meta
  )
}

estimate_stage4_to_death <- function(dat, tau = 20, origin = c("all", "prevalent", "incident")) {
  origin <- match.arg(origin)
  d <- dat
  stage4_date <- if ("date_first_clinical_cvd" %in% names(d)) {
    as_date2(d$date_first_clinical_cvd)
  } else {
    as_date2(d$date_ckm4)
  }
  baseline <- as_date2(d$date_attend)
  lost <- if ("date_lost" %in% names(d)) as_date2(d$date_lost) else as.Date(rep(NA, nrow(d)))
  death <- if ("date_death" %in% names(d)) as_date2(d$date_death) else as.Date(rep(NA, nrow(d)))
  admin <- rep(date_follow_end_global, nrow(d))
  censor_no_death <- pmin_date2(data.frame(lost = lost, admin = admin))
  censor_no_death[is.na(censor_no_death)] <- date_follow_end_global
  exit_date <- pmin_date2(data.frame(death = death, censor = censor_no_death))
  exit_date[is.na(exit_date)] <- censor_no_death[is.na(exit_date)]

  prevalent <- !is.na(stage4_date) & !is.na(baseline) & stage4_date <= baseline
  incident <- !is.na(stage4_date) & !is.na(baseline) & stage4_date > baseline & stage4_date <= exit_date
  origin_keep <- switch(origin, prevalent = prevalent, incident = incident, all = prevalent | incident)

  entry <- if (origin == "incident") rep(0, nrow(d)) else pmax(as.numeric(baseline - stage4_date) / 365.25, 0)
  exit <- as.numeric(exit_date - stage4_date) / 365.25
  ev <- as.integer(!is.na(death) & death <= censor_no_death & death > stage4_date)
  keep <- origin_keep & is.finite(entry) & is.finite(exit) & exit > entry

  if (sum(keep) < 100 || sum(ev[keep]) < 20) {
    return(data.table(
      cohort = origin, n = sum(keep), deaths = sum(ev[keep]), median = NA_real_,
      lo = NA_real_, hi = NA_real_, rmst = NA_real_
    ))
  }

  sf <- survfit(Surv(entry[keep], exit[keep], ev[keep]) ~ 1, conf.type = "log-log")
  tab <- summary(sf)$table
  med <- unname(tab["median"])
  lo <- unname(tab["0.95LCL"])
  hi <- unname(tab["0.95UCL"])
  data.table(
    cohort = origin,
    n = sum(keep),
    deaths = sum(ev[keep]),
    median = med,
    lo = lo,
    hi = hi,
    rmst = survfit_rmst(sf, tau)
  )
}

# Aalen-Johansen cumulative incidence of clinical stage-4 proxy, treating death
# before stage 4 as a competing event rather than non-informative censoring.
stage4_cif <- function(dat, horizon) {
  d <- dat
  baseline <- as_date2(d$date_attend)
  event_date <- if ("date_first_clinical_cvd" %in% names(d)) {
    as_date2(d$date_first_clinical_cvd)
  } else {
    as_date2(d$date_ckm4)
  }
  death <- if ("date_death" %in% names(d)) as_date2(d$date_death) else as.Date(rep(NA, nrow(d)))
  lost <- if ("date_lost" %in% names(d)) as_date2(d$date_lost) else as.Date(rep(NA, nrow(d)))
  admin <- rep(date_follow_end_global, nrow(d))
  censor_date <- pmin_date2(data.frame(lost = lost, admin = admin))
  censor_date[is.na(censor_date)] <- date_follow_end_global

  invalid_event <- is.na(event_date) | is.na(baseline) |
    (!is.na(event_date) & !is.na(baseline) & event_date <= baseline)
  invalid_death <- is.na(death) | is.na(baseline) |
    (!is.na(death) & !is.na(baseline) & death <= baseline)
  event_date[invalid_event] <- as.Date(NA)
  death[invalid_death] <- as.Date(NA)
  stop_date <- pmin_date2(data.frame(event = event_date, death = death, censor = censor_date))
  status <- rep(0L, nrow(d))
  event_first <- !is.na(event_date) & event_date == stop_date &
    (is.na(death) | event_date <= death) & event_date <= censor_date
  death_first <- !event_first & !is.na(death) & death == stop_date & death <= censor_date
  event_first[is.na(event_first)] <- FALSE
  death_first[is.na(death_first)] <- FALSE
  status[event_first] <- 1L
  status[death_first] <- 2L
  time <- as.numeric(stop_date - baseline) / 365.25
  keep <- is.finite(time) & time > 0
  if (sum(keep) < 100 || sum(status[keep] == 1L) < 5) return(NA_real_)

  status_f <- factor(status[keep], levels = 0:2, labels = c("censor", "stage4", "death"))
  fit <- try(survfit(Surv(time[keep], status_f) ~ 1), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  ss <- try(summary(fit, times = horizon, extend = TRUE), silent = TRUE)
  if (inherits(ss, "try-error") || is.null(ss$pstate)) return(NA_real_)
  ps <- ss$pstate
  if (is.null(dim(ps))) ps <- matrix(ps, nrow = 1)
  cn <- colnames(ps)
  j <- if (!is.null(cn) && any(grepl("stage4", cn, ignore.case = TRUE))) {
    which(grepl("stage4", cn, ignore.case = TRUE))[1]
  } else if (ncol(ps) >= 2) {
    2L
  } else {
    NA_integer_
  }
  if (!is.finite(j)) return(NA_real_)
  as.numeric(ps[nrow(ps), j])
}

spread_label_positions <- function(x, min_gap) {
  x <- safe_num(x)
  out <- x
  ok <- which(is.finite(x))
  if (length(ok) <= 1) return(out)
  ord <- ok[order(x[ok])]
  y <- x[ord]
  for (i in 2:length(y)) {
    if (y[i] - y[i - 1] < min_gap) y[i] <- y[i - 1] + min_gap
  }
  out[ord] <- y
  out
}


make_clock_reference_generic <- function(dd, covars, n_ref = clock_ref_n, mode = clock_reference, seed_offset = 0L) {
  idx <- seq_len(nrow(dd))
  set.seed(seed + 1117L + as.integer(seed_offset))
  if (length(idx) > n_ref) idx <- sample(idx, n_ref)
  ref <- dd[idx, covars, drop = FALSE]
  if (mode == "healthy_le8") {
    healthy <- dd$ckm.stage == 0
    if (sum(healthy, na.rm = TRUE) < 500) healthy <- dd$ckm.stage %in% c(0, 1)
    for (v in intersect(clock_behavior_vars, names(ref))) {
      if (is.numeric(dd[[v]]) || is.integer(dd[[v]])) {
        q <- quantile(safe_num(dd[[v]][healthy]), 0.80, na.rm = TRUE, names = FALSE, type = 8)
        if (is.finite(q)) ref[[v]] <- q
      }
    }
  }
  ref
}

make_latent_stage_road <- function(dat, covars, reference_mode = clock_reference) {
  dd <- dat[dat$ckm.stage %in% 0:4, , drop = FALSE]
  need <- unique(c("ckm.stage", covars))
  dd <- dd[complete.cases(dd[, need, drop = FALSE]), , drop = FALSE]
  if (nrow(dd) < 1000 || length(unique(dd$ckm.stage)) < 5) {
    stop("Too few complete participants or missing CKM categories for the latent stage road", call. = FALSE)
  }

  ref <- make_clock_reference_generic(
    dd,
    covars,
    n_ref = clock_ref_n,
    mode = reference_mode,
    seed_offset = 8300L
  )

  threshold_rows <- rbindlist(lapply(0:3, function(k) {
    y <- as.integer(dd$ckm.stage > k)
    empirical_above <- mean(y, na.rm = TRUE)
    fit_status <- "empirical"
    p_above <- empirical_above

    if (length(covars)) {
      fit_dat <- dd[, covars, drop = FALSE]
      fit_dat$.y <- y
      fm <- as.formula(paste0(".y ~ ", paste(bt(covars), collapse = " + ")))
      fit <- try(glm(fm, family = binomial(), data = fit_dat, model = FALSE, x = FALSE, y = FALSE), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        pr <- try(predict(fit, newdata = ref, type = "response"), silent = TRUE)
        if (!inherits(pr, "try-error") && any(is.finite(pr))) {
          p_above <- mean(pr[is.finite(pr)], na.rm = TRUE)
          fit_status <- "standardized logistic"
        }
      }
    }

    p_above <- pmin(pmax(p_above, 1e-6), 1 - 1e-6)
    p_le <- 1 - p_above
    data.table(
      threshold = k,
      contrast = paste0("s", k, " | s", k + 1L, "+"),
      N = nrow(dd),
      above_N = sum(y, na.rm = TRUE),
      empirical_p_above = empirical_above,
      reference_p_above = p_above,
      reference_p_le = p_le,
      cutpoint_raw = qlogis(p_le),
      fit_status = fit_status
    )
  }))

  # Numerical projection only; cumulative cutpoints must be strictly increasing.
  cutpoints <- safe_num(threshold_rows$cutpoint_raw)
  for (i in 2:length(cutpoints)) {
    if (!is.finite(cutpoints[i])) cutpoints[i] <- cutpoints[i - 1] + 1e-3
    if (cutpoints[i] <= cutpoints[i - 1]) cutpoints[i] <- cutpoints[i - 1] + 1e-3
  }
  threshold_rows[, cutpoint := cutpoints]

  positive_gaps <- diff(cutpoints)
  positive_gaps <- positive_gaps[is.finite(positive_gaps) & positive_gaps > 0]
  edge_gap <- if (length(positive_gaps)) median(positive_gaps) else 1
  centers <- c(
    cutpoints[1] - edge_gap / 2,
    (cutpoints[1] + cutpoints[2]) / 2,
    (cutpoints[2] + cutpoints[3]) / 2,
    (cutpoints[3] + cutpoints[4]) / 2,
    cutpoints[4] + edge_gap / 2
  )
  if (!all(is.finite(centers)) || diff(range(centers)) <= 0) centers <- 0:4
  latent_position <- 4 * (centers - min(centers)) / diff(range(centers))

  stage_counts <- as.data.table(dd)[, .(N = .N), by = .(stage = as.integer(ckm.stage))]
  road <- merge(
    data.table(
      stage = 0:4,
      latent_center = centers,
      position_from_s0 = latent_position
    ),
    stage_counts,
    by = "stage",
    all.x = TRUE,
    sort = TRUE
  )
  road[, `:=`(
    label = unname(stage_full_labels[as.character(stage)]),
    plot_label = unname(stage_short_labels[as.character(stage)]),
    stage_col = unname(stage_palette[as.character(stage)]),
    axis_type = "adjusted latent severity"
  )]

  intervals <- data.table(
    from_stage = 0:3,
    to_stage = 1:4,
    x_from = latent_position[1:4],
    x_to = latent_position[2:5]
  )
  intervals[, `:=`(
    latent_distance = x_to - x_from,
    x_mid = (x_from + x_to) / 2,
    segment = paste0("s", from_stage, "→s", to_stage),
    stage_col = unname(stage_palette[as.character(from_stage)])
  )]

  metadata <- data.table(
    item = c(
      "method", "outcome", "reference", "covariates", "scale",
      "interpretation", "limitation"
    ),
    value = c(
      "four standardized cumulative-logit models",
      "baseline ordered CKM stage s0-s4",
      reference_mode,
      paste(covars, collapse = ", "),
      "normalized latent severity position from 0 to 4",
      "relative adjusted separation of CKM categories",
      "not time and not an observed longitudinal transition barrier"
    )
  )

  list(road = road, intervals = intervals, thresholds = threshold_rows, metadata = metadata)
}


make_risk_calibrated_stage_road <- function(
  dat,
  covars,
  reference_mode = clock_reference,
  horizon = 10
) {
  dd <- dat[
    dat$yin & dat$ckm.stage %in% 0:3 &
      is.finite(dat$time) & dat$time > 0 & !is.na(dat$event),
    ,
    drop = FALSE
  ]
  dd$stage_f <- factor(dd$ckm.stage, levels = 0:3)
  need <- unique(c("time", "event", "ckm.stage", "stage_f", covars))
  dd <- dd[complete.cases(dd[, need, drop = FALSE]), , drop = FALSE]

  if (nrow(dd) < 1000 || sum(dd$event, na.rm = TRUE) < 50) {
    stop("Too few participants/events for the risk-calibrated CKM road", call. = FALSE)
  }

  fm <- as.formula(
    paste0(
      "Surv(time, event) ~ stage_f",
      if (length(covars)) paste0(" + ", paste(bt(covars), collapse = " + ")) else ""
    )
  )
  fit <- coxph(
    fm,
    data = dd,
    ties = "efron",
    model = FALSE,
    x = FALSE,
    y = FALSE
  )

  ref0 <- make_clock_reference_generic(
    dd,
    covars,
    n_ref = clock_ref_n,
    mode = reference_mode,
    seed_offset = 13100L
  )
  refs <- do.call(rbind, lapply(0:3, function(s0) {
    z <- ref0
    z$stage_f <- factor(s0, levels = 0:3)
    z$stage <- s0
    z
  }))
  rownames(refs) <- NULL

  bh <- basehaz(fit, centered = FALSE)
  ii <- which(is.finite(bh$time) & bh$time <= horizon)
  h0 <- if (length(ii)) bh$hazard[max(ii)] else 0

  lp <- try(predict(fit, newdata = refs, type = "lp", reference = "zero"), silent = TRUE)
  if (inherits(lp, "try-error")) {
    lp <- predict(fit, newdata = refs, type = "lp")
    lp <- lp - mean(lp[refs$stage == 0], na.rm = TRUE)
  }
  risk_i <- 1 - exp(-h0 * exp(safe_num(lp)))

  risk_tab <- data.table(stage = refs$stage, risk = risk_i)[
    ,
    .(standardized_risk = mean(risk, na.rm = TRUE)),
    by = stage
  ][order(stage)]
  risk_tab[, standardized_risk := pmin(pmax(standardized_risk, 1e-6), 1 - 1e-6)]
  risk_tab[, standardized_risk_monotone := pava_monotone(standardized_risk, decreasing = FALSE)]
  risk_tab[, standardized_risk_monotone := pmin(pmax(standardized_risk_monotone, 1e-6), 1 - 1e-6)]

  # Complementary-log-log transformation is the natural Cox-model risk scale.
  # It expands clinically meaningful differences among low absolute risks while
  # preserving their order. s0-s3 are normalized to 0-3; clinical s4 is anchored at 4.
  risk_tab[, cloglog_risk := log(-log1p(-standardized_risk_monotone))]
  rr <- range(risk_tab$cloglog_risk, na.rm = TRUE)
  if (!all(is.finite(rr)) || diff(rr) <= 1e-8) {
    risk_tab[, position_from_s0 := 0:3]
  } else {
    risk_tab[, position_from_s0 := 3 * (cloglog_risk - rr[1]) / diff(rr)]
  }

  stage_counts <- as.data.table(dat)[
    ckm.stage %in% 0:4,
    .(N = .N),
    by = .(stage = as.integer(ckm.stage))
  ]
  road <- rbind(
    risk_tab[, .(
      stage,
      standardized_risk,
      standardized_risk_monotone,
      cloglog_risk,
      position_from_s0
    )],
    data.table(
      stage = 4L,
      standardized_risk = NA_real_,
      standardized_risk_monotone = NA_real_,
      cloglog_risk = NA_real_,
      position_from_s0 = 4
    ),
    fill = TRUE
  )
  road <- merge(road, stage_counts, by = "stage", all.x = TRUE, sort = TRUE)
  road[, `:=`(
    position_from_stage0 = position_from_s0,
    label = unname(stage_full_labels[as.character(stage)]),
    plot_label = unname(stage_short_labels[as.character(stage)]),
    stage_col = unname(stage_palette[as.character(stage)]),
    axis_type = "risk-calibrated severity"
  )]

  intervals <- data.table(
    from_stage = 0:3,
    to_stage = 1:4,
    x_from = road[match(0:3, stage), position_from_s0],
    x_to = road[match(1:4, stage), position_from_s0]
  )
  intervals[, `:=`(
    severity_distance = fifelse(to_stage <= 3, x_to - x_from, NA_real_),
    distance_estimand = fifelse(to_stage <= 3, "prospective-risk separation", "clinical anchor; not estimated"),
    x_mid = (x_from + x_to) / 2,
    segment = paste0("s", from_stage, "→s", to_stage),
    stage_col = unname(stage_palette[as.character(from_stage)])
  )]

  model_table <- data.table(
    term = names(coef(fit)),
    estimate = as.numeric(coef(fit)),
    se = sqrt(diag(vcov(fit))),
    HR = exp(as.numeric(coef(fit)))
  )
  model_table[, `:=`(
    lo = exp(estimate - 1.96 * se),
    hi = exp(estimate + 1.96 * se),
    p = 2 * pnorm(-abs(estimate / se))
  )]

  metadata <- data.table(
    item = c(
      "method", "outcome", "risk horizon", "reference", "covariates",
      "risk transform", "scale", "interpretation", "limitation"
    ),
    value = c(
      "standardized Cox model",
      "incident CKM stage 4 among baseline s0-s3",
      horizon,
      reference_mode,
      paste(covars, collapse = ", "),
      "complementary log-log of standardized absolute risk",
      "s0-s3 normalized to 0-3; clinical s4 anchored at 4",
      "adjusted separation according to future stage-4 risk",
      "not time and not an observed adjacent-stage transition duration"
    )
  )

  list(
    road = road,
    intervals = intervals,
    standardized_risk = risk_tab,
    model = model_table,
    metadata = metadata,
    horizon = horizon
  )
}

# Risk-age is an interpretable quantitative CKM coordinate: the stage log-hazard
# contrast divided by the age log-hazard coefficient.  It answers how many
# chronological-age years carry the same prospective stage-4 hazard contrast.
# It is deliberately not described as time to the next CKM stage.
make_risk_age_stage_clock <- function(dat, covars, reference_age = risk_age_reference_age) {
  dd <- dat[
    dat$yin & dat$ckm.stage %in% 0:3 &
      is.finite(dat$time) & dat$time > 0 & !is.na(dat$event),
    , drop = FALSE
  ]
  dd$stage_f <- factor(dd$ckm.stage, levels = 0:3)
  other_covars <- setdiff(safe_covars(dd, covars), "age")
  need <- unique(c("time", "event", "age", "ckm.stage", "stage_f", other_covars))
  dd <- dd[complete.cases(dd[, need, drop = FALSE]), , drop = FALSE]
  if (nrow(dd) < 1000 || sum(dd$event) < 50) stop("Too few observations for risk-age CKM clock", call. = FALSE)

  fm <- as.formula(paste0(
    "Surv(time,event) ~ age + stage_f",
    if (length(other_covars)) paste0(" + ", paste(bt(other_covars), collapse = " + ")) else ""
  ))
  fit <- coxph(fm, data = dd, ties = "efron", model = FALSE, x = FALSE, y = FALSE)
  co <- coef(fit)
  vv <- vcov(fit)
  age_term <- names(co)[names(co) == "age"][1]
  if (is.na(age_term) || !is.finite(co[age_term]) || co[age_term] <= 0) {
    stop("Risk-age clock requires a positive estimable age coefficient", call. = FALSE)
  }

  out <- rbindlist(lapply(0:3, function(st) {
    if (st == 0) return(data.table(
      stage = 0L, logHR_vs_s0 = 0, risk_age_years = 0,
      lo = 0, hi = 0, equivalent_age = reference_age,
      N = sum(dd$ckm.stage == 0), events = sum(dd$event[dd$ckm.stage == 0])
    ))
    term <- paste0("stage_f", st)
    if (!term %in% names(co)) return(data.table(
      stage = st, logHR_vs_s0 = NA_real_, risk_age_years = NA_real_,
      lo = NA_real_, hi = NA_real_, equivalent_age = NA_real_,
      N = sum(dd$ckm.stage == st), events = sum(dd$event[dd$ckm.stage == st])
    ))
    a <- unname(co[term]); b <- unname(co[age_term])
    est <- a / b
    var_est <- vv[term, term] / b^2 + (a^2 * vv[age_term, age_term]) / b^4 -
      2 * a * vv[term, age_term] / b^3
    se_est <- sqrt(pmax(var_est, 0))
    data.table(
      stage = st, logHR_vs_s0 = a, risk_age_years = est,
      lo = est - 1.96 * se_est, hi = est + 1.96 * se_est,
      equivalent_age = reference_age + est,
      N = sum(dd$ckm.stage == st), events = sum(dd$event[dd$ckm.stage == st])
    )
  }))
  out[, `:=`(
    label = unname(stage_full_labels[as.character(stage)]),
    plot_label = unname(stage_short_labels[as.character(stage)]),
    axis_type = "prospective stage-4 risk-age",
    position_from_stage0 = risk_age_years,
    position_from_s0 = risk_age_years,
    reference_age = reference_age
  )]
  intervals <- data.table(
    from_stage = 0:2, to_stage = 1:3,
    from_years = out[match(0:2, stage), risk_age_years],
    to_years = out[match(1:3, stage), risk_age_years]
  )
  intervals[, `:=`(
    delta_risk_age_years = to_years - from_years,
    segment = paste0("s", from_stage, "→s", to_stage)
  )]
  model <- data.table(
    term = names(co), estimate = as.numeric(co), se = sqrt(diag(vv)),
    HR = exp(as.numeric(co))
  )
  model[, `:=`(
    lo = exp(estimate - 1.96 * se), hi = exp(estimate + 1.96 * se),
    p = 2 * pnorm(-abs(estimate / se))
  )]
  metadata <- data.table(
    item = c("method", "outcome", "reference stage", "reference age", "covariates", "interpretation", "limitation"),
    value = c(
      "Cox-equivalent risk-age: stage log-HR divided by age log-HR",
      "incident clinical-CVD CKM stage-4 proxy among baseline s0-s3",
      "s0", reference_age, paste(other_covars, collapse = ", "),
      "chronological-age years with an equivalent prospective stage-4 hazard contrast",
      "not elapsed time to the next CKM stage; assumes locally linear age effect on log hazard"
    )
  )
  list(road = out, intervals = intervals, model = model, metadata = metadata, fit = fit)
}

make_stage_clock <- function(dat, write_outputs = TRUE) {
  write_clock_tab <- if (write_outputs) {
    function(x, filename) write_tab(x, filename, dir = global_rawdir)
  } else function(...) invisible(NULL)
  write_clock_book <- if (write_outputs) {
    function(x, filename) write_book(x, filename, dir = global_outroot)
  } else function(...) invisible(NULL)
  save_clock_plot <- function(p, filename, width = 12, height = 8) {
    if (write_outputs) save_plot(p, filename, width = width, height = height, dir = global_outroot)
  }

  cov_sets <- get_clock_covariate_sets(dat)

  model_objects <- list()
  points <- data.table()
  fit_summaries <- data.table()
  fit_metadata <- data.table()

  for (i in seq_along(cov_sets)) {
    nm <- names(cov_sets)[i]
    t0 <- Sys.time()
    progress_log("Fitting stage-4-free RMST profile: ", nm)
    rr <- try(
      survreg_clock(
        dat,
        safe_covars(dat, cov_sets[[nm]]),
        draws = if (nm == clock_adjust_mode) clock_draws else 0,
        tau = clock_tau,
        seed_offset = i * 1000L
      ),
      silent = TRUE
    )
    if (inherits(rr, "try-error")) {
      warning("RMST model failed for ", nm, ": ", as.character(rr), call. = FALSE)
      next
    }

    model_objects[[nm]] <- rr
    pp <- copy(rr$point)
    pp[, model := nm]
    points <- rbind(points, pp, fill = TRUE)

    fs <- copy(rr$fit_summary)
    fs[, model := nm]
    fit_summaries <- rbind(fit_summaries, fs, fill = TRUE)

    fm <- copy(rr$fit_meta)
    fm[, model := nm]
    fit_metadata <- rbind(fit_metadata, fm, fill = TRUE)

    invisible(gc())
    progress_log("Fitting stage-4-free RMST profile: ", nm, " done in ", elapsed_since(t0))
  }

  if (!clock_adjust_mode %in% names(model_objects)) {
    stop(
      "The requested primary RMST adjustment model failed: ", clock_adjust_mode,
      ". Sensitivity models will not be substituted for the primary model.",
      call. = FALSE
    )
  }
  primary_name <- clock_adjust_mode

  primary_rmst <- copy(model_objects[[primary_name]]$point)[order(stage)]
  primary_rmst[, `:=`(
    label = unname(stage_full_labels[as.character(stage)]),
    plot_label = unname(stage_short_labels[as.character(stage)]),
    stage_col = unname(stage_palette[as.character(stage)]),
    model = primary_name,
    rmst_loss_from_s0 = remaining_years[stage == 0][1] - remaining_years
  )]
  tau_primary <- unique(primary_rmst$tau_effective)
  tau_primary <- tau_primary[is.finite(tau_primary)][1]
  if (!is.finite(tau_primary)) tau_primary <- clock_tau

  risk_by_stage <- rbindlist(lapply(0:3, function(s0) {
    dd <- dat[dat$yin & dat$ckm.stage == s0, , drop = FALSE]
    data.table(
      ckm.stage = s0,
      stage_label = unname(stage_short_labels[as.character(s0)]),
      N = nrow(dd),
      events = sum(dd$event, na.rm = TRUE),
      cif_5y = stage4_cif(dd, 5),
      cif_10y = stage4_cif(dd, 10)
    )
  }))

  death4 <- rbindlist(lapply(c("incident", "prevalent", "all"), function(origin) {
    estimate_stage4_to_death(dat, tau = clock_tau, origin = origin)
  }), fill = TRUE)

  # The internal quantitative road remains a standardized prospective-risk
  # coordinate for downstream models. Fig1a now displays adjusted absolute risk
  # directly and shows prevalent clinical stage 4 as a separate state, avoiding
  # an artificial estimated s3-to-s4 distance.
  risk_axis <- make_risk_calibrated_stage_road(
    dat,
    safe_covars(dat, cov_sets[[primary_name]]),
    reference_mode = clock_reference,
    horizon = min(10, tau_primary)
  )
  severity_road <- copy(risk_axis$road)
  severity_road[, `:=`(
    model = primary_name,
    clock_type = "risk-calibrated severity"
  )]
  severity_intervals <- copy(risk_axis$intervals)
  risk_age <- make_risk_age_stage_clock(
    dat,
    safe_covars(dat, cov_sets[[primary_name]]),
    reference_age = risk_age_reference_age
  )

  latent_sensitivity <- try(
    make_latent_stage_road(
      dat,
      safe_covars(dat, cov_sets[[primary_name]]),
      reference_mode = clock_reference
    ),
    silent = TRUE
  )

  # Panel a: adjusted absolute risk for s0-s3; s4 is a separate clinical state.
  stage_plot <- copy(severity_road)
  stage_plot <- merge(
    stage_plot,
    risk_age$road[, .(stage, risk_age_years)],
    by = "stage", all.x = TRUE, sort = FALSE
  )[order(stage)]
  max_risk_pct <- max(100 * stage_plot[stage <= 3, standardized_risk_monotone], na.rm = TRUE)
  min_risk_pct <- min(100 * stage_plot[stage <= 3, standardized_risk_monotone], na.rm = TRUE)
  clinical_gap <- max(1.5, 0.25 * max(max_risk_pct - min_risk_pct, 1))
  s4_x <- max_risk_pct + clinical_gap
  stage_plot[, plot_x := fifelse(stage <= 3, 100 * standardized_risk_monotone, s4_x)]
  stage_plot[, label_x := spread_label_positions(plot_x, max(0.5, 0.055 * diff(range(plot_x))))]
  stage_plot[, label_y := 1.12 + 0.055 * (stage %% 2)]
  # Requested presentation: sample size is part of each stage label.
  stage_plot[, stage_n_label := paste0(plot_label, " (", scales::comma(N), ")")]
  stage_plot[, risk_label := ifelse(
    stage <= 3 & is.finite(standardized_risk_monotone),
    paste0(
      "adj ", scales::percent(standardized_risk_monotone, accuracy = 0.1),
      ifelse(is.finite(risk_age_years), sprintf("; %.1f y", risk_age_years), "")
    ),
    "prevalent clinical CVD"
  )]
  stage_plot[, risk_label_y := ifelse(stage <= 3, 0.925 - 0.035 * (stage %% 2), 0.91)]

  interval_plot <- data.table(from_stage = 0:2, to_stage = 1:3)
  interval_plot[, `:=`(
    x_from = stage_plot[match(from_stage, stage), plot_x],
    x_to = stage_plot[match(to_stage, stage), plot_x]
  )]
  # data.table evaluates all right-hand sides in one `:=` call before
  # assigning the new columns. Create delta_risk_pp first so it is available
  # when distance_label is formatted in the following call.
  interval_plot[, delta_risk_pp := x_to - x_from]
  interval_plot[, `:=`(
    x_mid = (x_from + x_to) / 2,
    distance_label = sprintf("Δrisk=%.1f pp", delta_risk_pp),
    label_y = 0.82 - 0.045 * (from_stage %% 2),
    stage_col = unname(stage_palette[as.character(from_stage)])
  )]

  death_incident <- death4[cohort == "incident"]
  death_prevalent <- death4[cohort == "prevalent"]
  death_text <- if (nrow(death_incident) && is.finite(death_incident$median[1])) {
    sprintf("incident clinical CVD → death: median %.1f y\n(separate post-event time scale)", death_incident$median[1])
  } else if (nrow(death_incident) && is.finite(death_incident$rmst[1])) {
    sprintf("incident clinical CVD → death: RMST %.1f y\n(separate post-event time scale)", death_incident$rmst[1])
  } else if (nrow(death_prevalent) && is.finite(death_prevalent$median[1])) {
    sprintf("prevalent Yang CVD → death: median %.1f y\n(delayed entry; survivor-selected)", death_prevalent$median[1])
  } else {
    "clinical CVD → death: not estimable"
  }

  x_breaks <- stage_plot$plot_x
  x_labels <- ifelse(
    stage_plot$stage <= 3,
    sprintf("%.1f%%", stage_plot$plot_x),
    "clinical CVD"
  )
  x_range <- range(stage_plot$plot_x, na.rm = TRUE)
  x_pad <- max(0.8, 0.08 * diff(x_range))

  p1 <- ggplot() +
    geom_hline(yintercept = 1, linewidth = 0.9, color = "grey68") +
    geom_segment(
      data = interval_plot,
      aes(x = x_from, xend = x_to, y = 1, yend = 1, color = factor(from_stage)),
      linewidth = 2.4,
      lineend = "round",
      show.legend = FALSE
    ) +
    annotate(
      "segment",
      x = stage_plot[stage == 3, plot_x], xend = stage_plot[stage == 4, plot_x],
      y = 1, yend = 1,
      linetype = "dashed", linewidth = 0.8, color = "grey60"
    ) +
    geom_point(
      data = stage_plot,
      aes(x = plot_x, y = 1, fill = factor(stage)),
      shape = 21,
      color = "white",
      stroke = 1.1,
      size = 4.8,
      show.legend = FALSE
    ) +
    geom_segment(
      data = stage_plot,
      aes(x = plot_x, xend = label_x, y = 1.015, yend = label_y - 0.012),
      linewidth = 0.35,
      color = "grey55"
    ) +
    geom_text(
      data = stage_plot,
      aes(x = label_x, y = label_y, label = stage_n_label, color = factor(stage)),
      fontface = "bold",
      size = 3.5,
      show.legend = FALSE
    ) +
    geom_text(
      data = stage_plot,
      aes(x = plot_x, y = risk_label_y, label = risk_label),
      size = 2.45,
      fontface = "bold",
      color = "grey20"
    ) +
    geom_text(
      data = interval_plot,
      aes(x = x_mid, y = label_y, label = distance_label, color = factor(from_stage)),
      fontface = "bold",
      size = 2.8,
      show.legend = FALSE
    ) +
    annotate(
      "label",
      x = max(min_risk_pct, s4_x - 1.15 * clinical_gap),
      y = 1.27,
      label = death_text,
      hjust = 0,
      size = 2.8,
      fill = "white",
      color = "black",
      label.size = 0.25
    ) +
    scale_fill_manual(values = stage_palette, drop = FALSE) +
    scale_color_manual(values = stage_palette, drop = FALSE) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(x_range[1] - x_pad, x_range[2] + x_pad)) +
    coord_cartesian(ylim = c(0.70, 1.34), clip = "off") +
    labs(
      title = "a. CKM risk-calibrated stage profile",
      subtitle = paste0(
        risk_axis$horizon, "-year standardized stage-4 proxy risk; s4 is prevalent clinical CVD.\n",
        "adj = standardized risk; Δrisk = adjacent-stage percentage-point difference."
      ),
      x = "Risk-calibrated clinical-CVD burden",
      y = NULL
    ) +
    theme_ckm(11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(8, 20, 8, 8)
    )

  points[, model := factor(model, levels = names(cov_sets))]
  y_min <- min(points$remaining_years, na.rm = TRUE)
  y_max <- max(points$remaining_years, na.rm = TRUE)
  y_pad <- max(0.15, 0.08 * (y_max - y_min))
  p2 <- ggplot(points, aes(stage, remaining_years, group = model, color = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    scale_x_continuous(breaks = 0:3, labels = unname(stage_short_labels[as.character(0:3)])) +
    coord_cartesian(ylim = c(y_min - y_pad, min(clock_tau, y_max + y_pad))) +
    scale_color_manual(values = adjust_palette[names(cov_sets)], drop = FALSE) +
    labs(
      title = "c. Stage-4-free RMST profile",
      subtitle = paste0("Restricted mean stage-4-free time over ≤", sprintf("%.1f", tau_primary), " years."),
      x = "Baseline CKM stage",
      y = "Stage-4-free RMST (years)",
      color = "Adjustment"
    ) +
    theme_ckm(10)

  risk_long <- risk_by_stage %>%
    pivot_longer(c(cif_5y, cif_10y), names_to = "horizon", values_to = "risk") %>%
    mutate(
      stage_label = factor(stage_label, levels = unname(stage_short_labels[as.character(0:3)])),
      horizon = factor(
        horizon,
        levels = c("cif_5y", "cif_10y"),
        labels = c("5-year risk", "10-year risk")
      )
    )

  p3 <- ggplot(risk_long, aes(stage_label, risk, fill = horizon)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = risk_horizon_palette) +
    labs(
      title = "d. Observed cumulative incidence of the stage-4 proxy",
      x = "Baseline CKM stage",
      y = "Cumulative incidence (death competing)",
      fill = NULL
    ) +
    theme_ckm(10)

  save_clock_plot(
    p1 / p2 / p3,
    "Fig1.stage_clock.png",
    9,
    16
  )

  latent_road_sens <- if (!inherits(latent_sensitivity, "try-error")) copy(latent_sensitivity$road) else data.table()
  latent_intervals_sens <- if (!inherits(latent_sensitivity, "try-error")) copy(latent_sensitivity$intervals) else data.table()
  latent_thresholds_sens <- if (!inherits(latent_sensitivity, "try-error")) copy(latent_sensitivity$thresholds) else data.table()
  latent_metadata_sens <- if (!inherits(latent_sensitivity, "try-error")) copy(latent_sensitivity$metadata) else data.table(
    item = "status",
    value = paste0("latent sensitivity failed: ", as.character(latent_sensitivity))
  )

  clock_meta <- data.table(
    item = c(
      "clock cohort", "primary adjustment", "reference", "healthy-reference variables",
      "RMST estimand", "primary quantitative-stage estimand", "risk-calibrated sensitivity estimand", "primary risk horizon",
      "requested RMST horizon", "effective RMST horizon", "stage-4 endpoint",
      "stage 4-to-death method", "stage-4 incidence method", "latent sensitivity",
      "interpretation warning"
    ),
    value = c(
      "full UKB analysis cohort before omic restriction",
      primary_name,
      clock_reference,
      paste(clock_behavior_vars, collapse = ", "),
      "standardized restricted mean stage-4-free time for baseline s0-s3",
      "prospective CKM risk-age: stage log-HR divided by age log-HR; s0-s3 only",
      "standardized future-stage-4 risk on the complementary-log-log scale; s4 shown separately",
      risk_axis$horizon,
      clock_tau,
      tau_primary,
      paste0("earliest CHD/HF/stroke/PAD/AF date from ", stage4_source_primary, " source definition"),
      "incident and prevalent cohorts reported separately; prevalent analysis uses delayed entry",
      "Aalen-Johansen cumulative incidence with death as a competing event",
      "prevalence-derived cumulative-logit road retained as sensitivity only",
      "No stage-to-stage years are claimed because repeated s0-s3 transition dates are unavailable"
    )
  )

  write_clock_tab(primary_rmst, "stage_clock.rmst_profile.tsv")
  write_clock_tab(risk_age$road, "stage_clock.risk_age.tsv")
  write_clock_tab(risk_age$intervals, "stage_clock.risk_age_intervals.tsv")
  write_clock_tab(risk_age$model, "stage_clock.risk_age_model.tsv")
  write_clock_tab(risk_age$metadata, "stage_clock.risk_age_metadata.tsv")
  write_clock_tab(severity_road, "stage_clock.risk_calibrated_road.tsv")
  write_clock_tab(severity_intervals, "stage_clock.risk_calibrated_intervals.tsv")
  write_clock_tab(risk_axis$standardized_risk, "stage_clock.standardized_risk.tsv")
  write_clock_tab(risk_axis$model, "stage_clock.risk_model.tsv")
  write_clock_tab(risk_axis$metadata, "stage_clock.risk_axis_metadata.tsv")
  if (nrow(latent_road_sens)) write_clock_tab(latent_road_sens, "stage_clock.latent_road_sensitivity.tsv")
  if (nrow(latent_intervals_sens)) write_clock_tab(latent_intervals_sens, "stage_clock.latent_intervals_sensitivity.tsv")
  if (nrow(latent_thresholds_sens)) write_clock_tab(latent_thresholds_sens, "stage_clock.latent_thresholds_sensitivity.tsv")
  write_clock_tab(points, "stage_clock.models.tsv")
  write_clock_tab(risk_by_stage, "stage_clock.observed_risk.tsv")
  write_clock_tab(death4, "stage_clock.stage4_death.tsv")
  write_clock_tab(fit_summaries, "stage_clock.fit_coefficients.tsv")
  write_clock_tab(fit_metadata, "stage_clock.fit_metadata.tsv")
  write_clock_tab(clock_meta, "stage_clock.metadata.tsv")

  write_clock_book(
    list(
      clock_metadata = clock_meta,
      rmst_profile = primary_rmst,
      risk_age = risk_age$road,
      risk_age_intervals = risk_age$intervals,
      risk_age_model = risk_age$model,
      risk_age_metadata = risk_age$metadata,
      risk_calibrated_road = severity_road,
      risk_calibrated_intervals = severity_intervals,
      standardized_risk = risk_axis$standardized_risk,
      risk_model = risk_axis$model,
      risk_axis_metadata = risk_axis$metadata,
      latent_road_sensitivity = latent_road_sens,
      latent_intervals_sensitivity = latent_intervals_sens,
      latent_thresholds_sensitivity = latent_thresholds_sens,
      latent_metadata_sensitivity = latent_metadata_sens,
      models = points,
      observed_risk = risk_by_stage,
      stage4_death = death4,
      fit_coefficients = fit_summaries,
      fit_metadata = fit_metadata
    ),
    "Fig1.stage_clock.out.xlsx"
  )

  clock_save <- list(
    stage3_method = stage3_method,
    primary = risk_age$road,
    risk_age = risk_age$road,
    risk_age_intervals = risk_age$intervals,
    risk_age_model = risk_age$model,
    primary_risk = severity_road,
    primary_rmst = primary_rmst,
    intervals = severity_intervals,
    standardized_risk = risk_axis$standardized_risk,
    risk_model = risk_axis$model,
    latent_sensitivity = list(
      road = latent_road_sens,
      intervals = latent_intervals_sens,
      thresholds = latent_thresholds_sens,
      metadata = latent_metadata_sens
    ),
    models = points,
    observed_risk = risk_by_stage,
    stage4_death = death4,
    fit_summaries = fit_summaries,
    fit_metadata = fit_metadata,
    metadata = clock_meta,
    plots = list(profile = p1, rmst = p2, incidence = p3)
  )
  if (write_outputs) save_rds_safe(clock_save, file.path(global_rawdir, "stage_clock.rds"), compress = FALSE)

  invisible(gc())
  clock_save
}

get_clock_positions <- function(clock_obj, type = c("risk_age", "risk", "rmst_loss", "latent_sensitivity")) {
  type <- match.arg(type)
  if (identical(type, "rmst_loss")) {
    z <- as.data.table(clock_obj$primary_rmst)
    if (!"rmst_loss_from_s0" %in% names(z)) {
      z[, rmst_loss_from_s0 := remaining_years[stage == 0][1] - remaining_years]
    }
    out <- setNames(z[stage %in% 0:3, rmst_loss_from_s0], z[stage %in% 0:3, stage])
    out <- c(out, `4` = max(out, na.rm = TRUE))
    return(out)
  }
  if (identical(type, "latent_sensitivity")) {
    z <- as.data.table(clock_obj$latent_sensitivity$road)
    if (!nrow(z)) return(setNames(0:4, 0:4))
    return(setNames(z[stage %in% 0:4, position_from_s0], z[stage %in% 0:4, stage]))
  }
  if (identical(type, "risk")) {
    z <- as.data.table(clock_obj$primary_risk)
    return(setNames(z[stage %in% 0:3, position_from_stage0], z[stage %in% 0:3, stage]))
  }
  z <- as.data.table(clock_obj$risk_age %||% clock_obj$primary)
  setNames(z[stage %in% 0:3, position_from_stage0], z[stage %in% 0:3, stage])
}
