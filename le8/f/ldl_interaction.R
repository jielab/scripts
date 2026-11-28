# Single prespecified biomarker. Both components use the observed LDL training-fold SD.
ldl_interaction <- function(d, covars, tvar, evar, dx, grouping_unit="supplied_group") {
    G <- ".genetic_z"; R <- ".residual_z"
    cols <- unique(c(G, R, tvar, evar, ".group", ".prevalent", covars))
    stopifnot(all(cols %in% names(d)))
    d <- as.data.frame(d)
    d <- d[complete.cases(d[cols]) & d$.prevalent %in% 0 &
           is.finite(d[[G]]) & is.finite(d[[R]]) & is.finite(d[[tvar]]) &
           d[[tvar]] > 0 & d[[evar]] %in% 0:1, , drop = FALSE]
    ref <- if (nrow(d)) unname(quantile(d[[G]], c(.1,.5,.9))) else rep(NA_real_, 3)
    fit_window <- function(lo, hi, analysis, restricted = FALSE) {
        z <- d[d[[tvar]] > lo, , drop = FALSE]
        if (restricted) z <- if ("drug.lipid" %in% names(z))
            z[z$drug.lipid %in% 0, , drop = FALSE] else z[FALSE, , drop = FALSE]
        z$.stop <- pmin(z[[tvar]], hi) - lo
        z$.event <- as.integer(z[[evar]] == 1 & z[[tvar]] <= hi)
        n <- nrow(z); events <- sum(z$.event)
        out <- data.frame(dx = dx, feature = "Clinical_LDL_C", analysis = analysis,
            start_year = lo, end_year = if (is.infinite(hi)) NA_real_ else hi,
            N = n, events = events, clustering_unit = grouping_unit, cluster_N = length(unique(z$.group)),
            families = if (grouping_unit=="individual") NA_integer_ else length(unique(z$.group)), person_years = sum(z$.stop),
            beta_GxR = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
            p_interaction = NA_real_, status = "insufficient_events_or_variation")
        empty <- function(msg) { out$status <- msg; list(test = out, slopes = NULL) }
        if (n == 0L && restricted) return(empty("no_eligible_baseline_nonusers"))
        cv <- covars[vapply(z[covars], function(x) length(unique(x)) > 1L, logical(1))]
        if (events < 20L || length(unique(z[[G]])) < 2L || length(unique(z[[R]])) < 2L) return(list(test=out, slopes=NULL))
        ff <- reformulate(c(paste(G, R, sep=" * "), if (length(cv)) paste0("`", cv, "`")),
                          response="survival::Surv(.stop,.event)")
        fit <- tryCatch(withCallingHandlers(survival::coxph(ff, data=z, ties="efron",
                       cluster=.group, robust=TRUE, singular.ok=FALSE),
                       warning=function(w) stop(conditionMessage(w))), error=identity)
        if (inherits(fit, "error")) return(empty(paste("cox_fit_unavailable", conditionMessage(fit))))
        term <- paste(G, R, sep=":"); b <- stats::coef(fit); V <- stats::vcov(fit)
        if (!all(c(term, R) %in% names(b)) || any(!is.finite(b)) || any(!is.finite(V)) ||
                V[term,term] <= 0) return(empty("unestimable_interaction"))
        se <- sqrt(V[term,term]); beta <- unname(b[term]); critical <- qnorm(.975)
        out[c("beta_GxR","se","ci_low","ci_high","p_interaction")] <-
            list(beta, se, beta-critical*se, beta+critical*se, 2*pnorm(-abs(beta/se)))
        out$status <- "estimated"
        slope <- unname(b[R]) + ref*beta
        variance <- V[R,R] + 2*ref*V[R,term] + ref^2*V[term,term]
        err <- sqrt(pmax(0, variance))
        slopes <- data.frame(dx=dx, analysis=analysis, start_year=lo, end_year=out$end_year,
            genetic_percentile=c(10,50,90), genetic_value=ref,
            within_window_G_range=ref >= min(z[[G]]) & ref <= max(z[[G]]),
            residual_HR=exp(slope), low=exp(slope-critical*err), high=exp(slope+critical*err))
        list(test=out, slopes=slopes)
    }
    cuts <- c(0,.5,1,2,5,10,16)
    fits <- c(list(fit_window(0,Inf,"full_followup")),
              lapply(seq_len(length(cuts)-1L), function(i) fit_window(cuts[i],cuts[i+1L],"interval")),
              list(fit_window(0,Inf,"baseline_no_lipid_medication",TRUE)))
    tests <- do.call(rbind, lapply(fits, "[[", "test"))
    tests$p_interval_holm <- NA_real_
    ii <- which(tests$analysis == "interval")
    tests$p_interval_holm[ii] <- p.adjust(tests$p_interaction[ii], "holm", n=length(ii))
    list(tests=tests, slopes=do.call(rbind, lapply(fits, "[[", "slopes")),
         unit="One common training-fold observed LDL SD; Wald inference conditional on cross-fitted calibration, clustered by the reported grouping unit")
}
