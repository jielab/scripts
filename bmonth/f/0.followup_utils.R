# Focused, explicitly exploratory follow-up helpers. No UKB data are loaded here.
# Source 0.common.R and phewas_engine.R before model-fitting helpers.

fu_complete <- function(d, variables) {
  variables <- unique(variables)
  missing <- setdiff(variables, names(d))
  if (length(missing)) stop("Missing follow-up columns: ", paste(missing, collapse=", "))
  keep <- complete.cases(d[, variables, drop=FALSE])
  for (v in variables) if (is.numeric(d[[v]])) keep <- keep & is.finite(d[[v]])
  droplevels(d[keep, , drop=FALSE])
}

fu_sample_key <- function(d) {
  # Only a digest is published; participant identifiers remain local.
  bm_hash(sort(as.character(d$eid)))
}

fu_sibship_stratum <- function(size) {
  size <- suppressWarnings(as.numeric(as.character(size)))
  ans <- rep(NA_character_, length(size))
  valid <- is.finite(size) & size == floor(size) & size >= 2
  ans[valid & size <= 5] <- as.character(size[valid & size <= 5])
  ans[valid & size >= 6] <- "6+"
  ans
}

fu_harmonic_basis <- function(month) {
  if (any(!is.finite(month) | month < 1 | month > 12)) stop("Month must be between 1 and 12")
  cbind(bmonth_cos=-cos(2*pi*month/12), bmonth_sin=sin(2*pi*month/12))
}

fu_cycle_contrasts <- function(beta, covariance, reference=1L) {
  terms <- c("bmonth_cos", "bmonth_sin")
  if (!all(terms %in% names(beta)) || !all(terms %in% rownames(covariance)) ||
      !all(terms %in% colnames(covariance))) stop("A joint cos+sin coefficient vector and covariance are required")
  b <- beta[terms]; v <- covariance[terms, terms, drop=FALSE]
  if (any(!is.finite(c(b,v))) || max(abs(v-t(v))) > 1e-8 ||
      min(eigen(v, symmetric=TRUE, only.values=TRUE)$values) <= 0)
    stop("Invalid joint coefficient covariance")
  contrast <- sweep(fu_harmonic_basis(1:12), 2, fu_harmonic_basis(reference)[1,], "-")
  log_hr <- as.vector(contrast %*% b)
  se <- sqrt(pmax(0, rowSums((contrast %*% v)*contrast)))
  data.frame(month=1:12, reference_month=reference, log_HR=log_hr, SE=se,
    HR=exp(log_hr), conf.low=exp(log_hr-qnorm(.975)*se),
    conf.high=exp(log_hr+qnorm(.975)*se), interval="95% pointwise; joint model; reference January")
}

fu_adjust_family <- function(d, family_cols, p_col="p") {
  if (!nrow(d)) return(d)
  d$family_n_eligible <- 0L; d$family_n_valid <- 0L
  d$FDR_family <- NA_real_; d$bonferroni_family <- NA_real_
  keys <- do.call(interaction,c(d[,family_cols,drop=FALSE],list(drop=TRUE,lex.order=TRUE)))
  for (index in split(seq_len(nrow(d)), keys)) {
    eligible <- !is.na(d$eligible[index]) & d$eligible[index]
    valid <- eligible & is.finite(d[[p_col]][index]) & d[[p_col]][index] >= 0 & d[[p_col]][index] <= 1
    n <- sum(eligible)
    d$family_n_eligible[index] <- n; d$family_n_valid[index] <- sum(valid)
    if (!n) next
    # Failed eligible fits enter the planned family as p=1, but their raw p and
    # adjusted p remain NA. A smaller converged-only denominator is not used.
    padded <- rep(1, n)
    padded[which(valid[eligible])] <- d[[p_col]][index][valid]
    adj <- p.adjust(padded, "BH")
    d$FDR_family[index[valid]] <- adj[which(valid[eligible])]
    d$bonferroni_family[index[valid]] <- pmin(1, d[[p_col]][index][valid]*n)
  }
  d
}

fu_scalar <- function(d, exposure, phenotype, covariates, analysis, sample_set,
                      family_adjusted=FALSE, stratum="all", extra_note="") {
  covariates <- setdiff(unique(covariates), exposure)
  z <- fu_complete(d, c("eid", ".y", exposure, covariates))
  x <- if (nrow(z)) bm_formula_matrix(z, c(exposure, covariates), tested=exposure) else
    matrix(numeric(),nrow=0L,ncol=2L,dimnames=list(NULL,c("(Intercept)",exposure)))
  dropped_terms <- attr(x,"dropped_terms")
  # bm_scalar's second-column contract is retained explicitly, even if a
  # nuisance factor disappears in this particular phenotype/risk set.
  if (!exposure %in% colnames(x)) x <- cbind(x, setNames(data.frame(z[[exposure]]), exposure))
  x <- as.matrix(x[, c("(Intercept)", exposure, setdiff(colnames(x), c("(Intercept)", exposure))), drop=FALSE])
  spec <- list(name=exposure, rows=seq_len(nrow(z)), x=x, dropped_terms=dropped_terms)
  r <- bm_scalar(spec, as.integer(z$.y), phenotype)
  r$exposure <- exposure; r$analysis <- analysis; r$sample_set <- sample_set
  r$family_adjusted <- family_adjusted; r$sibship_stratum <- stratum
  r$sample_key <- fu_sample_key(z)
  r$covariates <- paste(covariates, collapse=" + ")
  r$selection_status <- "exploratory: candidates selected from the 2026-10-20 audit"
  r$interpretation_note <- extra_note
  if (is.finite(r$beta) && is.finite(r$SE)) {
    r$conf.low <- exp(r$beta-qnorm(.975)*r$SE)
    r$conf.high <- exp(r$beta+qnorm(.975)*r$SE)
  } else r$conf.low <- r$conf.high <- NA_real_
  r
}

# Retain descriptive deltas while making both sample and estimator changes
# explicit. Failed models cannot contribute a numeric difference.
fu_comparison <- function(a,b,comparison,expect_same_sample=TRUE) {
  valid <- function(r) isTRUE(grepl("^ok_",r$status)) && isTRUE(is.finite(r$beta))
  known <- function(x) length(x)==1L && !is.na(x) && nzchar(x)
  same_method <- if(known(a$method) && known(b$method)) a$method==b$method else NA
  same_sample <- if(known(a$sample_key) && known(b$sample_key)) a$sample_key==b$sample_key else NA
  notes <- "Descriptive log-odds difference; not a mediation proportion"
  if(!valid(a) || !valid(b)) notes <- c(notes,paste0("Invalid endpoint: ",a$status," -> ",b$status,"; delta unavailable"))
  if(isFALSE(same_method)) notes <- c(notes,"Estimation method also changed (估计方法同时改变); delta cannot be attributed solely to adjustment or sample selection")
  if(is.na(same_method)) notes <- c(notes,"Endpoint method unavailable")
  if(is.na(same_sample)) notes <- c(notes,"Endpoint sample identity unavailable")
  if(isFALSE(same_sample)) notes <- c(notes,if(expect_same_sample)
    "Expected common sample differs; delta also includes sample change" else "Different samples: sample-selection comparison")
  data.frame(phenotype=a$phenotype,exposure=a$exposure,comparison=comparison,
    from_analysis=a$analysis,to_analysis=b$analysis,from_method=a$method,to_method=b$method,
    from_status=a$status,to_status=b$status,from_sample_key=a$sample_key,to_sample_key=b$sample_key,
    from_beta=a$beta,to_beta=b$beta,delta=if(valid(a)&&valid(b)) b$beta-a$beta else NA_real_,
    same_method=same_method,same_sample=same_sample,comparison_note=paste(notes,collapse="; "),
    stringsAsFactors=FALSE)
}

fu_sample_adjustment <- function(g) {
  select <- function(a) g[g$analysis==a,,drop=FALSE]
  rows <- lapply(c("full_D2","common_D2_adultSES","full_D2_family_size",
    "common_D2_adultSES_family_size","common_D5_migration_family_size"),select)
  if(any(vapply(rows,nrow,integer(1))!=1L)) stop("Missing/duplicate targeted comparison endpoints")
  names(rows) <- c("full_D2","common_D2","full_D2_family","common_D2_family","common_D5_family")
  a<-rows[[1]]; b<-rows[[2]]; c<-rows[[3]]; f<-rows[[4]]; e<-rows[[5]]
  comparisons <- dplyr::bind_rows(
    fu_comparison(a,b,"sample_selection",FALSE),
    fu_comparison(a,c,"family_adjustment_full"),
    fu_comparison(f,e,"D2_to_D5_on_common_family_sample"))
  summary <- data.frame(phenotype=a$phenotype,exposure=a$exposure,N_full=a$n_total,N_common=b$n_total,
    cases_full=a$n_cases,cases_common=b$n_cases,stringsAsFactors=FALSE)
  for(name in names(rows)) {
    summary[[paste0("beta_",name)]] <- rows[[name]]$beta
    for(field in c("method","status","sample_key")) summary[[paste0(field,"_",name)]] <- rows[[name]][[field]]
  }
  for(i in seq_len(nrow(comparisons))) {
    label <- comparisons$comparison[i]
    summary[[paste0("delta_",label)]] <- comparisons$delta[i]
    for(field in c("same_method","same_sample","comparison_note"))
      summary[[paste0(label,"_",field)]] <- comparisons[[field]][i]
  }
  summary$note <- "Descriptive differences of log-odds coefficients; inspect endpoint methods/status/sample keys; not mediation proportions"
  list(summary=summary,comparisons=comparisons)
}

fu_age_basis <- function(d, df=4L) {
  # The source birth_date is month-15 imputation, not exact date of birth.
  # Use a single prespecified age source for this sensitivity sample; never
  # silently mix fractional age for some people and integer age for others.
  fractional <- "age_fractional" %in% names(d) && any(is.finite(d$age_fractional))
  source <- if (fractional) "age_fractional" else "age"
  age <- d[[source]]
  ok <- is.finite(age)
  if (sum(ok) < 20L || length(unique(age[ok])) <= df)
    return(list(data=d, terms=character(), source=source, status="insufficient_age_variation"))
  design <- splines::ns(age[ok], df=df)
  terms <- paste0(".age_ns", seq_len(ncol(design)))
  for (i in seq_along(terms)) {
    d[[terms[i]]] <- NA_real_; d[[terms[i]]][ok] <- design[,i]
  }
  list(data=d, terms=terms, source=source, status="ok")
}

fu_cox_design <- function(d, covariates, tested=c("bmonth_cos", "bmonth_sin"), tol=1e-9) {
  # Drop constant nuisance terms, but do not sacrifice adjustment variables
  # to make the exposure estimable. Test columns are checked after nuisance.
  cv <- setdiff(unique(covariates), tested)
  varying <- vapply(cv, function(v) length(unique(d[[v]])) > 1L, logical(1))
  dropped <- cv[!varying]; cv <- cv[varying]
  nuisance <- model.matrix(if (length(cv)) reformulate(cv) else ~1, d)
  qr0 <- qr(nuisance, tol=tol, LAPACK=FALSE)
  keep <- qr0$pivot[seq_len(qr0$rank)]
  dropped_columns <- setdiff(colnames(nuisance), colnames(nuisance)[keep])
  nuisance <- nuisance[,keep,drop=FALSE]
  test_x <- as.matrix(d[,tested,drop=FALSE])
  qr1 <- qr(cbind(nuisance,test_x), tol=tol, LAPACK=FALSE)
  if (qr1$rank - ncol(nuisance) != length(tested))
    return(list(ok=FALSE, note="joint cycle not estimable after nuisance adjustment"))
  nuisance <- nuisance[,setdiff(colnames(nuisance),"(Intercept)"),drop=FALSE]
  if (ncol(nuisance)) colnames(nuisance) <- sprintf(".Z%03d",seq_len(ncol(nuisance)))
  note <- c(if(length(dropped)) paste0("constant term: ",dropped),
            if(length(dropped_columns)) paste0("aliased nuisance: ",dropped_columns))
  list(ok=TRUE, x=cbind(test_x,nuisance),note=paste(note,collapse="; "))
}

fu_joint_cox <- function(dat, outcome, covariates, analysis="D4_joint_cycle") {
  dat <- ensure_outcome(dat, outcome)
  tv <- paste0(outcome,".t2e"); ev <- paste0(outcome,".Yt2e")
  tested <- c("bmonth_cos","bmonth_sin")
  d <- fu_complete(dat, c("eid",tv,ev,tested,covariates))
  r <- data.frame(outcome=outcome,analysis=analysis,N=nrow(d),events=sum(d[[ev]]==1),
    eligible=nrow(d)>=500L && sum(d[[ev]]==1)>=20L, status="skipped_insufficient_events",
    beta_cos=NA_real_, beta_sin=NA_real_, var_cos=NA_real_, cov_cos_sin=NA_real_, var_sin=NA_real_,
    chisq=NA_real_,df=2L,p=NA_real_,PH_global_p=NA_real_,PH_cos_p=NA_real_,PH_sin_p=NA_real_,
    method="joint Wald; unpenalized Cox; Efron ties",sample_key=fu_sample_key(d),
    covariates=paste(covariates,collapse=" + "),note="",stringsAsFactors=FALSE)
  empty <- function() list(summary=r,curve=data.frame())
  if (!r$eligible) return(empty())
  design <- fu_cox_design(d,covariates,tested)
  r$note <- design$note
  if (!design$ok) { r$status <- "failed_nonestimable"; return(empty()) }
  model_data <- data.frame(.time=d[[tv]],.event=d[[ev]],design$x,check.names=FALSE)
  warning_text <- character()
  f <- reformulate(colnames(design$x),response="survival::Surv(.time, .event)")
  fit <- tryCatch(withCallingHandlers(survival::coxph(f,data=model_data,ties="efron",x=TRUE,y=TRUE),
    warning=function(w) { warning_text <<- c(warning_text,conditionMessage(w)); invokeRestart("muffleWarning") }),
    error=function(e) { warning_text <<- c(warning_text,conditionMessage(e)); NULL })
  r$note <- paste(c(r$note,unique(warning_text)),collapse="; ")
  if (is.null(fit) || any(grepl("converg|infinite|singular|ran out",warning_text,ignore.case=TRUE))) {
    r$status <- "failed_cox_fit"; return(empty())
  }
  beta <- stats::coef(fit)[tested]; v <- stats::vcov(fit)[tested,tested,drop=FALSE]
  if (any(!is.finite(c(beta,v))) || min(eigen(v,symmetric=TRUE,only.values=TRUE)$values)<=0) {
    r$status <- "failed_cox_covariance"; return(empty())
  }
  chi <- tryCatch(as.numeric(t(beta) %*% solve(v,beta)),error=function(e)NA_real_)
  if (!is.finite(chi) || chi < -1e-7) {
    r$status <- "failed_cox_Wald"; return(empty())
  }
  chi <- max(0,chi)
  r$status <- "ok"; r$beta_cos <- beta[1]; r$beta_sin <- beta[2]
  r$var_cos <- v[1,1]; r$cov_cos_sin <- v[1,2]; r$var_sin <- v[2,2]
  r$chisq <- chi; r$p <- pchisq(chi,df=2,lower.tail=FALSE)
  ph <- tryCatch(survival::cox.zph(fit,transform="km"),error=function(e) e)
  if (inherits(ph,"error")) r$note <- paste(r$note,"PH diagnostic:",conditionMessage(ph)) else {
    tab <- ph$table
    if ("GLOBAL" %in% rownames(tab)) r$PH_global_p <- tab["GLOBAL","p"]
    if (tested[1] %in% rownames(tab)) r$PH_cos_p <- tab[tested[1],"p"]
    if (tested[2] %in% rownames(tab)) r$PH_sin_p <- tab[tested[2],"p"]
  }
  curve <- tryCatch(fu_cycle_contrasts(beta,v),error=function(e)e)
  if (inherits(curve,"error") || any(!is.finite(as.matrix(curve[,c("HR","conf.low","conf.high")])))) {
    r$status <- "failed_cycle_numeric_range"; r$p <- NA_real_
    r$note <- paste(r$note,"Monthly contrast or interval was not finite")
    return(empty())
  }
  curve$outcome <- outcome; curve$analysis <- analysis; curve$N <- nrow(d)
  curve$events <- r$events; curve$joint_p <- r$p
  list(summary=r,curve=curve)
}
