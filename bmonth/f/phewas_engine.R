# Engine v6: retain the packed-panel/memory-bounded execution design.
# A valid result requires numerical diagnostics, not just fit$converged.
# Each worker holds a design matrix and one phenotype, never the full panel.
BM_ENGINE_VERSION <- "6.0.0"
BM_PANEL_VERSION <- 2L
BM_FIT_SETTINGS <- list(mle_maxit=100L, br_maxit=100L, epsilon=1e-8,
                        qr_tol=1e-9, boundary_tol=10*.Machine$double.eps,
                        fallback="brglm2_AS_mean", joint_fallback="coefficient_Wald")
bm_atomic_rds <- function(value, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- tempfile('.writing-', tmpdir=dirname(path)); on.exit(unlink(tmp))
  saveRDS(value, tmp)
  if (!file.rename(tmp, path)) stop('Cannot publish checkpoint: ', path)
  invisible(path)
}
bm_hash <- function(value) {
  tmp <- tempfile(); on.exit(unlink(tmp)); saveRDS(value,tmp,compress=FALSE)
  unname(tools::md5sum(tmp))
}
bm_read_cache <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path),error=function(e) NULL)
}
bm_cached_panel <- function(directory, input_md5=NULL, verify_binary=TRUE) {
  manifest <- file.path(directory,'panel.rds'); binary <- file.path(directory,'panel.bin')
  meta <- bm_read_cache(manifest)
  if (is.list(meta) && identical(meta$format_version,BM_PANEL_VERSION) &&
      is.character(meta$panel_md5) && length(meta$panel_md5)==1L &&
      (is.null(input_md5) || identical(meta$panel_md5,input_md5)) &&
      is.numeric(meta$n) && length(meta$n)==1L && is.finite(meta$n) && meta$n>0 &&
      is.list(meta$counts) && identical(names(meta$counts),meta$phenotypes) &&
      all(vapply(meta$counts,function(x) is.numeric(x) && length(x)==3L &&
        identical(names(x),c("cases","controls","missing")) &&
        all(is.finite(x) & x>=0) && sum(x)==meta$n,logical(1))) &&
      is.character(meta$ids) &&
      length(meta$ids)==meta$n && is.character(meta$phenotypes) &&
      length(meta$phenotypes)>0L && file.exists(binary) &&
      isTRUE(file.info(binary)$size == as.double(meta$n)*length(meta$phenotypes)) &&
      (!verify_binary || identical(meta$packed_md5,unname(tools::md5sum(binary))))) {
    meta$path <- binary; return(meta)
  }
  NULL
}
bm_panel <- function(input, directory, input_md5=NULL) {
  if(is.null(input_md5)) input_md5 <- unname(tools::md5sum(input))
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  manifest <- file.path(directory,'panel.rds'); binary <- file.path(directory,'panel.bin')
  meta <- bm_cached_panel(directory,input_md5)
  if (!is.null(meta)) {
    message('[PheWAS] START Reusing packed phenotype panel'); return(meta)
  }
  ph <- readRDS(input)
  if (!is.data.frame(ph) || !'id' %in% names(ph) || ncol(ph)<2L || !nrow(ph))
    stop('PheWAS.rds must contain id and phenotype columns')
  ids <- as.character(ph$id)
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids) || anyDuplicated(names(ph)))
    stop('Invalid or duplicate phenotype IDs/column names')
  ys <- setdiff(names(ph),'id'); counts <- setNames(vector('list',length(ys)),ys)
  tmp <- tempfile('.panel-',tmpdir=directory)
  con <- file(tmp,'wb')
  on.exit({if (!is.null(con)) close(con); unlink(tmp)},add=TRUE)
  for (i in seq_along(ys)) {
    y <- ph[[ys[i]]]
    if (!is.logical(y)) {
      if (!is.numeric(y)) {
        original <- as.character(y); y <- suppressWarnings(as.numeric(original))
        if(any(!is.na(original) & is.na(y)))
          stop('Non-numeric phenotype value in column: ',ys[i])
      }
      if (any(!is.na(y) & !y %in% c(0,1)))
        stop("Non-binary phenotype column (expected 0/1/NA): ",ys[i])
      y <- y==1
    }
    counts[[i]] <- c(cases=sum(y==1,na.rm=TRUE),controls=sum(y==0,na.rm=TRUE),missing=sum(is.na(y)))
    packed <- as.integer(y); packed[is.na(packed)] <- 2L
    writeBin(as.raw(packed),con)
    if (i %% 200L==0L || i==length(ys))
      message('[PheWAS] DONE Packing phenotypes ',i,'/',length(ys))
  }
  close(con); con <- NULL
  if (!file.rename(tmp,binary)) stop('Cannot publish phenotype panel')
  meta <- list(format_version=BM_PANEL_VERSION,panel_md5=input_md5,
    packed_md5=unname(tools::md5sum(binary)),counts=counts,
    n=nrow(ph),ids=ids,phenotypes=ys,path=binary)
  bm_atomic_rds(meta,manifest); meta
}
bm_panel_isolated <- function(input, directory, engine, input_md5=NULL) {
  if(is.null(input_md5)) input_md5 <- unname(tools::md5sum(input))
  meta <- bm_cached_panel(directory,input_md5)
  if (!is.null(meta)) {
    message('[PheWAS] START Reusing packed phenotype panel'); return(meta)
  }
  # A cold read leaves a large R heap and GC trigger even after rm()/gc().
  # Pack in a fresh process which exits BEFORE any model workers are forked,
  # so they cannot inherit (and privately copy) that multi-gigabyte heap.
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  script <- tempfile('bmonth-pack-',fileext='.R'); on.exit(unlink(script))
  writeLines(c('args <- commandArgs(TRUE)', 'source(args[[1]])',
               'invisible(bm_panel(args[[2]], args[[3]], args[[4]]))'),script)
  message('[PheWAS] START Packing phenotype panel in a separate R process')
  status <- system2(file.path(R.home('bin'),'Rscript'),
    c('--vanilla',shQuote(script),shQuote(normalizePath(engine,mustWork=TRUE)),
      shQuote(normalizePath(input,mustWork=TRUE)),
      shQuote(normalizePath(directory,mustWork=TRUE)),shQuote(input_md5)))
  if (status!=0L) stop('Phenotype packing process failed (exit ',status,'); models were not started')
  meta <- bm_cached_panel(directory,input_md5)
  if (is.null(meta)) stop('Phenotype packing process did not produce a valid panel')
  message('[PheWAS] DONE Packing process exited; panel ready for model workers')
  meta
}
bm_panel_column <- function(panel, phenotype) {
  j <- match(phenotype,panel$phenotypes)
  if(is.na(j)) stop("Phenotype absent from packed panel: ",phenotype)
  con <- file(panel$path,"rb"); on.exit(close(con))
  seek(con,where=as.double(j-1L)*panel$n,origin="start")
  bytes <- readBin(con,what="raw",n=panel$n)
  if(length(bytes)!=panel$n) stop("Truncated packed phenotype panel")
  y <- as.integer(bytes); y[y==2L] <- NA_integer_
  data.frame(eid=panel$ids,.y=y,stringsAsFactors=FALSE)
}

bm_prepare <- function(dat, panel, covariates) {
  idx <- match(as.character(dat$eid),panel$ids)
  keep <- !is.na(idx) & complete.cases(dat[,covariates,drop=FALSE])
  for (v in covariates) if (is.numeric(dat[[v]])) keep <- keep & is.finite(dat[[v]])
  d <- droplevels(dat[keep,,drop=FALSE])
  if (!nrow(d)) stop('No participants with complete covariates and phenotype IDs')
  x <- bm_formula_matrix(d,covariates); rownames(x) <- NULL
  list(dat=d,rows=idx[keep],covariates=x,dropped_terms=attr(x,'dropped_terms'))
}
bm_context <- function(prepared, names, panel, kind=c('scalar','joint')) {
  kind <- match.arg(kind); dat <- prepared$dat; specs <- list()
  for (name in names) {
    extra <- if (kind=='scalar') name else switch(name,
      month11df='birth_month',
      cycle2df=c('bmonth_cos','bmonth_sin'),
      north_interaction2df=c('birth_north_100km','bmonth_cos','bmonth_sin',
                            'bmonth_cos_x_north','bmonth_sin_x_north'),
      east_interaction2df=c('birth_east_100km','bmonth_cos','bmonth_sin',
                           'bmonth_cos_x_east','bmonth_sin_x_east'),
      move_interaction2df=c('move_log1p','bmonth_cos','bmonth_sin',
                           'bmonth_cos_x_move','bmonth_sin_x_move'),
      stop('Unknown joint test: ',name))
    keep <- complete.cases(dat[,extra,drop=FALSE])
    for (v in extra) keep <- keep & is.finite(dat[[v]])
    base <- prepared$covariates[keep,,drop=FALSE]
    if (kind=='scalar') {
      base <- base[,setdiff(colnames(base),name),drop=FALSE]
      x <- cbind(base[,1,drop=FALSE],dat[[name]][keep],base[,-1,drop=FALSE])
      colnames(x)[2] <- name; null_cols <- integer()
    } else {
      if (name=='month11df') {
        m <- factor(dat$birth_month[keep],levels=1:12)
        extra_x <- model.matrix(~m)[,-1,drop=FALSE]
      } else extra_x <- as.matrix(dat[keep,extra,drop=FALSE])
      # Birth coordinates already occur in the primary adjustment matrix.
      # Add each main effect only once, then append the tested interactions.
      extra_x <- extra_x[,setdiff(colnames(extra_x),colnames(base)),drop=FALSE]
      x <- cbind(base,extra_x)
      null_cols <- if (grepl('interaction',name))
        which(!colnames(x) %in% tail(extra,2L)) else seq_len(ncol(base))
    }
    rownames(x) <- NULL
    specs[[name]] <- list(name=name,kind=kind,rows=prepared$rows[keep],x=x,null_cols=null_cols,
                         dropped_terms=prepared$dropped_terms)
  }
  for (i in seq_along(specs)) {
    specs[[i]]$null_key <- as.character(i)
    if (kind=='joint' && i>1L) for (j in seq_len(i-1L)) {
      a <- specs[[i]]; b <- specs[[j]]
      if (identical(a$rows,b$rows) &&
          identical(a$x[,a$null_cols,drop=FALSE],b$x[,b$null_cols,drop=FALSE])) {
        specs[[i]]$null_key <- specs[[j]]$null_key; break
      }
    }
  }
  list(specs=specs,panel=panel[c('n','path','phenotypes')])
}
# Only single-level nuisance variables are omitted before model.matrix().
# Aliasing is evaluated again AFTER outcome-specific missingness is applied.
bm_formula_matrix <- function(data, terms, tested=character()) {
  missing <- setdiff(terms,names(data))
  if (length(missing)) stop('Missing design variables: ',paste(missing,collapse=', '))
  constant <- terms[vapply(data[,terms,drop=FALSE],function(v)
    length(unique(v[!is.na(v)]))<2L,logical(1))]
  drop <- setdiff(constant,tested)
  use <- setdiff(terms,drop)
  # A constant tested numeric variable is deliberately retained and rejected
  # as non-estimable by bm_prune_design, rather than silently changing the test.
  d <- data
  for (v in intersect(constant,tested)) if (is.factor(d[[v]]) || is.character(d[[v]]))
    d[[v]] <- rep(0,nrow(d))
  x <- model.matrix(if(length(use)) reformulate(use) else ~1,data=d)
  attr(x,'dropped_terms') <- drop
  x
}

bm_rank <- function(x) {
  if (!ncol(x)) return(list(rank=0L,pivot=integer()))
  norms <- sqrt(colSums(x*x))
  norms[!is.finite(norms) | norms==0] <- 1
  qr(sweep(x,2,norms,'/'),tol=BM_FIT_SETTINGS$qr_tol,LAPACK=FALSE)
}

bm_prune_design <- function(x,tested=character()) {
  x <- as.matrix(x)
  fail <- function(reason) list(ok=FALSE,status=reason,x=x,dropped=character())
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x)))
    return(fail('failed_invalid_design'))
  if (is.null(colnames(x)) || anyDuplicated(colnames(x)))
    return(fail('failed_invalid_column_names'))
  if (length(setdiff(tested,colnames(x)))) return(fail('failed_missing_tested_term'))
  nuisance <- which(!colnames(x) %in% tested)
  # Keep an independent nuisance basis. Do not preferentially retain a tested
  # variable that is identical to a confounder: that effect is not identifiable.
  if (length(nuisance)) {
    q0 <- bm_rank(x[,nuisance,drop=FALSE])
    nuisance <- if(q0$rank) nuisance[q0$pivot[seq_len(q0$rank)]] else integer()
  }
  tested_cols <- match(tested,colnames(x))
  use <- sort(c(nuisance,tested_cols))
  if (!length(use)) return(fail('failed_empty_design'))
  z <- x[,use,drop=FALSE]
  if (bm_rank(z)$rank != ncol(z)) return(fail('failed_tested_term_not_estimable'))
  list(ok=TRUE,status='ok',x=z,dropped=setdiff(colnames(x),colnames(z)))
}

# These are sufficient (not exhaustive) separation certificates. A nuisance
# dummy containing only cases or only controls is included. Multivariable
# separation without such a certificate is caught when numerical diagnostics
# fail; absence of a certificate is never labelled "no separation".
bm_separation_hint <- function(x,y) {
  if (!all(c(0,1) %in% unique(y))) return('outcome_has_one_level')
  hits <- vapply(seq_len(ncol(x)),function(j) {
    a <- range(x[y==0,j]); b <- range(x[y==1,j])
    span <- max(a,b)-min(a,b)
    span>0 && (a[2]<=b[1] || b[2]<=a[1])
  },logical(1))
  paste(colnames(x)[hits],collapse='|')
}

bm_fit <- function(x, y, coefficient=NULL, strict=FALSE, tested=coefficient) {
  # strict is retained for compatibility. Every fit now uses the strict checks.
  tested <- unique(as.character(tested))
  out <- list(ok=FALSE,status='failed',method=NA_character_,note='',warnings='',
    converged=FALSE,boundary=NA,iterations=NA_integer_,fallback_reason='',
    mle_status=NA_character_,br_status=NA_character_,dropped_nuisance='',
    separation_hint='',min_fitted=NA_real_,max_fitted=NA_real_)
  prep <- bm_prune_design(x,tested)
  if (!prep$ok) {out$status<-prep$status;out$note<-prep$status;return(out)}
  x <- prep$x; out$dropped_nuisance <- paste(prep$dropped,collapse='|')
  if (length(y)!=nrow(x) || anyNA(y) || any(!y %in% c(0,1))) {
    out$status<-'failed_invalid_response';out$note<-out$status;return(out)
  }
  out$separation_hint <- bm_separation_hint(x,y)
  attempt <- function(method) {
    notes <- character()
    fit <- tryCatch(withCallingHandlers({
      if (method=='mle') stats::glm.fit(x=x,y=y,family=binomial(),
        control=glm.control(maxit=BM_FIT_SETTINGS$mle_maxit,epsilon=BM_FIT_SETTINGS$epsilon))
      else brglm2::brglmFit(x=x,y=y,family=binomial(),
        control=list(type='AS_mean',maxit=BM_FIT_SETTINGS$br_maxit,
                     epsilon=BM_FIT_SETTINGS$epsilon),singular.ok=FALSE)
    },warning=function(w) {notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')}),
      error=function(e) {notes<<-c(notes,conditionMessage(e));NULL})
    reasons <- character(); sm <- NULL
    if (is.null(fit)) reasons <- 'fit_error' else {
      if (!isTRUE(fit$converged)) reasons<-c(reasons,'nonconverged')
      if (isTRUE(fit$boundary)) reasons<-c(reasons,'boundary')
      if (any(grepl('converg|infinite|singular|numerically 0 or 1|boundary|NaN',notes,ignore.case=TRUE)))
        reasons<-c(reasons,'numerical_warning')
      if (method=='mle' && nzchar(out$separation_hint)) reasons<-c(reasons,'separation_certificate')
      mu<-fit$fitted.values
      if (length(mu)!=length(y) || any(!is.finite(mu)) ||
          any(mu<=BM_FIT_SETTINGS$boundary_tol | mu>=1-BM_FIT_SETTINGS$boundary_tol))
        reasons<-c(reasons,'fitted_probability_boundary')
      sm <- tryCatch(summary.glm(fit,dispersion=1),error=function(e) {
        notes<<-c(notes,conditionMessage(e));NULL
      })
      if (is.null(sm) || nrow(sm$coefficients)!=ncol(x) ||
          any(!is.finite(sm$coefficients)) || any(sm$coefficients[,2]<=0))
        reasons<-c(reasons,'invalid_coefficient_or_se')
      if (!is.null(sm)) {
        v<-sm$cov.unscaled
        if (!is.matrix(v) || any(!is.finite(v)) || any(diag(v)<=0) ||
            is.null(tryCatch(chol(v),error=function(e) NULL)))
          reasons<-c(reasons,'invalid_covariance')
        if(length(tested) && all(tested %in% rownames(sm$coefficients))) {
          b<-sm$coefficients[tested,1]; se<-sm$coefficients[tested,2]
          er<-exp(c(b,b-qnorm(.975)*se,b+qnorm(.975)*se))
          if(any(!is.finite(er)) || any(er<=0)) reasons<-c(reasons,'unrepresentable_OR_or_CI')
        }
      }
      if (!is.finite(fit$deviance) || fit$rank!=ncol(x)) reasons<-c(reasons,'invalid_deviance_or_rank')
    }
    list(ok=!length(reasons),fit=fit,summary=sm,
         status=if(length(reasons))paste(unique(reasons),collapse='|') else 'ok',
         warnings=paste(unique(notes),collapse='; '))
  }
  a<-attempt('mle');out$mle_status<-a$status
  chosen<-a; method<-'mle'
  if (!a$ok) {
    out$fallback_reason<-a$status
    if (!requireNamespace('brglm2',quietly=TRUE)) {
      out$status<-'failed_brglm2_unavailable';out$br_status<-'package_unavailable'
      out$warnings<-a$warnings;out$note<-paste(a$status,a$warnings,'Install brglm2',sep='; ')
      return(out)
    }
    chosen<-attempt('brglm2_AS_mean');method<-'brglm2_AS_mean';out$br_status<-chosen$status
  }
  out$method<-method;out$ok<-chosen$ok
  out$status<-if(chosen$ok)if(method=='mle')'ok_mle' else 'ok_brglm2' else 'failed_numerical'
  out$warnings<-paste(unique(c(a$warnings,chosen$warnings)[nzchar(c(a$warnings,chosen$warnings))]),collapse='; ')
  out$note<-paste(c(if(nzchar(out$fallback_reason))paste0('fallback: ',out$fallback_reason),
                   if(!chosen$ok)paste0('final: ',chosen$status),out$warnings),collapse='; ')
  f<-chosen$fit
  if(!is.null(f)) {
    out$converged<-isTRUE(f$converged);out$boundary<-isTRUE(f$boundary)
    out$iterations<-as.integer(f$iter)
    if(length(f$fitted.values) && all(is.finite(f$fitted.values))) {
      out$min_fitted<-min(f$fitted.values);out$max_fitted<-max(f$fitted.values)
    }
  }
  if(!chosen$ok) return(out)
  out$deviance<-f$deviance;out$rank<-f$rank
  out$coefficients<-f$coefficients;out$vcov<-chosen$summary$cov.unscaled
  if(!is.null(coefficient) && length(coefficient)==1L && coefficient %in% rownames(chosen$summary$coefficients))
    out$coef<-unname(chosen$summary$coefficients[coefficient,])
  out
}

bm_fit_diagnostics <- function(fit) {
  fit[c('status','method','converged','boundary','iterations','fallback_reason',
        'mle_status','br_status','warnings','dropped_nuisance','separation_hint',
        'min_fitted','max_fitted','note')]
}
bm_empty_diagnostics <- function() list(status='skipped',method=NA_character_,
  converged=FALSE,boundary=NA,iterations=NA_integer_,fallback_reason='',mle_status=NA_character_,
  br_status=NA_character_,warnings='',dropped_nuisance='',separation_hint='',
  min_fitted=NA_real_,max_fitted=NA_real_,note='')

bm_scalar <- function(spec, y, phenotype) {
  y <- y[spec$rows]; keep <- !is.na(y); y <- y[keep]
  n <- length(y); cases <- sum(y)
  r <- data.frame(phenotype=phenotype,snp=spec$name,adjustment=NA_character_,
    beta=NA_real_,SE=NA_real_,OR=NA_real_,conf.low=NA_real_,conf.high=NA_real_,p=NA_real_,type='logistic',
    n_total=n,n_cases=cases,n_controls=n-cases,allele_freq=NA_real_,HWE_p=NA_real_,
    n_no_snp=0L,eligible=FALSE,stringsAsFactors=FALSE)
  r<-cbind(r,as.data.frame(bm_empty_diagnostics(),stringsAsFactors=FALSE))
  if (n<20L || cases<20L || n-cases<20L) {
    r$status<-'skipped_low_case_control_count';r$note<-'Fewer than 20 cases or controls';return(r)
  }
  x <- if (all(keep)) spec$x else spec$x[keep,,drop=FALSE]
  if(!spec$name %in% colnames(x) || length(unique(x[,spec$name]))<2L) {
    r$status<-'skipped_nonvarying_exposure';r$note<-'Exposure does not vary';return(r)
  }
  # Eligibility is decided before model fitting; failed eligible tests remain
  # in the planned multiplicity family, with raw p=NA.
  r$eligible<-TRUE
  fit<-bm_fit(x,y,coefficient=spec$name)
  r[names(bm_fit_diagnostics(fit))]<-bm_fit_diagnostics(fit)
  r$dropped_nuisance<-paste(unique(c(spec$dropped_terms,
    strsplit(r$dropped_nuisance,'|',fixed=TRUE)[[1]])),collapse='|')
  if(fit$ok && !is.null(fit$coef)) {
    r$beta<-fit$coef[1];r$SE<-fit$coef[2];r$OR<-exp(r$beta);r$p<-fit$coef[4]
    r$conf.low<-exp(r$beta-qnorm(.975)*r$SE);r$conf.high<-exp(r$beta+qnorm(.975)*r$SE)
  }
  r
}

bm_joint <- function(spec, y, phenotype, null_fits) {
  y <- y[spec$rows]; keep <- !is.na(y); y <- y[keep]
  r<-data.frame(phenotype=phenotype,N=length(y),cases=sum(y),p=NA_real_,eligible=FALSE,
    statistic=NA_real_,df=NA_integer_,test_method=NA_character_,
    null_status=NA_character_,null_method=NA_character_,null_converged=FALSE,null_boundary=NA,
    null_warnings='',null_note='',
    beta_cos=NA_real_,beta_sin=NA_real_,var_cos=NA_real_,var_sin=NA_real_,cov_cos_sin=NA_real_,
    stringsAsFactors=FALSE)
  r<-cbind(r,as.data.frame(bm_empty_diagnostics(),stringsAsFactors=FALSE))
  if(r$N<500L || r$cases<20L || r$N-r$cases<20L) {
    r$status<-'skipped_low_case_control_count';return(r)
  }
  r$eligible<-TRUE
  x<-if(all(keep))spec$x else spec$x[keep,,drop=FALSE]
  tested<-setdiff(colnames(x),colnames(x)[spec$null_cols])
  if(!length(tested)) {r$status<-'failed_no_tested_terms';return(r)}
  # Full/null nuisance pruning must use the same reduced nuisance basis.
  prep<-bm_prune_design(x,tested)
  if(!prep$ok) {r$status<-prep$status;r$note<-prep$status;return(r)}
  x<-prep$x;null_names<-setdiff(colnames(x),tested)
  b<-bm_fit(x,y,strict=TRUE,tested=tested)
  r[names(bm_fit_diagnostics(b))]<-bm_fit_diagnostics(b)
  r$dropped_nuisance<-paste(unique(c(spec$dropped_terms,prep$dropped,strsplit(b$dropped_nuisance,'|',fixed=TRUE)[[1]])),collapse='|')
  if(!b$ok) return(r)
  # Outcome-specific pruning can differ by test. Include the actual nuisance
  # matrix identity in the null cache key to prevent sharing the wrong fit.
  key<-paste0(spec$null_key,'_',bm_hash(list(rows=which(keep),names=null_names)))
  if(!exists(key,envir=null_fits,inherits=FALSE))
    assign(key,bm_fit(x[,null_names,drop=FALSE],y,strict=TRUE),envir=null_fits)
  a<-get(key,envir=null_fits,inherits=FALSE)
  r$null_status<-a$status;r$null_method<-a$method;r$null_note<-a$note
  r$null_converged<-a$converged;r$null_boundary<-a$boundary;r$null_warnings<-a$warnings
  k<-length(tested)
  if(a$ok && a$method=='mle' && b$method=='mle') {
    dev<-a$deviance-b$deviance;df<-b$rank-a$rank
    if(df!=k || !is.finite(dev) || dev< -1e-7) {
      r$status<-'failed_invalid_nested_LRT';r$note<-'Full/null rank or likelihood difference invalid';return(r)
    }
    r$statistic<-max(0,dev);r$df<-as.integer(df);r$test_method<-'ML_LRT'
  } else {
    # A joint Wald statistic uses the coefficient vector and its full
    # covariance from ONE valid full model. It is not a penalized-deviance LRT.
    v<-b$vcov[tested,tested,drop=FALSE];bb<-b$coefficients[tested]
    w<-tryCatch(drop(crossprod(bb,solve(v,bb))),error=function(e) NA_real_)
    if(!is.finite(w) || w< -1e-7) {r$status<-'failed_invalid_joint_Wald';return(r)}
    r$statistic<-max(0,w);r$df<-as.integer(k)
    r$test_method<-if(b$method=='mle')'ML_coefficient_Wald_null_unusable' else 'BR_coefficient_Wald'
  }
  r$p<-pchisq(r$statistic,df=r$df,lower.tail=FALSE)
  if(all(c('bmonth_cos','bmonth_sin') %in% names(b$coefficients))) {
    r$beta_cos<-unname(b$coefficients['bmonth_cos']);r$beta_sin<-unname(b$coefficients['bmonth_sin'])
    r$var_cos<-b$vcov['bmonth_cos','bmonth_cos'];r$var_sin<-b$vcov['bmonth_sin','bmonth_sin']
    r$cov_cos_sin<-b$vcov['bmonth_cos','bmonth_sin']
  }
  r
}

bm_adjust_family <- function(rows) {
  r<-as.data.frame(rows)
  if(!all(c('p','eligible','status') %in% names(r)))
    stop('Multiplicity correction requires engine-v6 p, eligible and status columns')
  eligible<-!is.na(r$eligible) & r$eligible
  valid<-eligible & grepl('^ok_',r$status) & is.finite(r$p) & r$p>=0 & r$p<=1
  r$p_raw<-r$p;r$p[!valid]<-NA_real_
  r$p_for_correction<-ifelse(eligible,ifelse(valid,r$p,1),NA_real_)
  r$FDR<-rep(NA_real_,nrow(r));r$FDR[eligible]<-p.adjust(r$p_for_correction[eligible],'BH')
  r$n_tests_planned<-rep(sum(eligible),nrow(r));r$n_tests_valid<-rep(sum(valid),nrow(r))
  r$n_tests_unresolved<-rep(sum(eligible & !valid),nrow(r))
  r$bonferroni_threshold<-rep(.05/max(1,sum(eligible)),nrow(r))
  r$bonferroni<-valid & r$p<=r$bonferroni_threshold
  r$logp<--log10(pmax(r$p,.Machine$double.xmin))
  r$correction_family<-rep('eligible phenotypes within exposure and covariate model',nrow(r))
  r$unresolved_policy<-rep('p=1 for correction only; raw p remains NA; resolve before final claims',nrow(r))
  r
}

bm_batch <- function(task) {
  context <- get('bm_worker_context',envir=.GlobalEnv); panel <- context$panel
  con <- file(panel$path,'rb'); on.exit(close(con))
  results <- vector('list',length(task$indices))
  for (k in seq_along(task$indices)) {
    i <- task$indices[k]
    seek(con,where=as.double(i-1L)*panel$n,origin='start')
    bytes <- readBin(con,what='raw',n=panel$n)
    if (length(bytes)!=panel$n) stop('Truncated packed phenotype panel')
    y <- as.integer(bytes); y[y==2L] <- NA_integer_
    null_fits <- new.env(parent=emptyenv())
    results[[k]] <- lapply(context$specs,function(s) {
      r <- if (s$kind=='scalar') bm_scalar(s,y,panel$phenotypes[i]) else
        bm_joint(s,y,panel$phenotypes[i],null_fits)
      r$exposure <- s$name; r
    })
  }
  bm_atomic_rds(list(indices=task$indices,results=results),task$path)
  invisible(gc(FALSE))
  list(count=length(task$indices),pid=Sys.getpid())
}
bm_memory_available <- function() {
  gib <- 1024^3; available <- 24*gib
  if (file.exists('/proc/meminfo')) {
    line <- grep('^MemAvailable:',readLines('/proc/meminfo'),value=TRUE)
    if (length(line)) available <- min(available,as.numeric(gsub('[^0-9]','',line))*1024)
  }
  if (file.exists('/proc/self/cgroup')) {
    cg <- sub('^0::','',grep('^0::',readLines('/proc/self/cgroup'),value=TRUE))
    if (length(cg)) {
      directory <- paste0('/sys/fs/cgroup',cg)
      repeat {
        read_number <- function(name) suppressWarnings(as.numeric(tryCatch(
          readLines(file.path(directory,name),warn=FALSE),error=function(e) NA)))
        limit <- read_number('memory.max'); used <- read_number('memory.current')
        if (length(limit)==1L && length(used)==1L && is.finite(limit) && is.finite(used))
          available <- min(available,limit-used)
        if (directory=='/sys/fs/cgroup' || directory=='/') break
        directory <- dirname(directory)
      }
    }
  }
  max(0,available-2*gib)
}
bm_worker_count <- function(context, requested='auto') {
  cpus <- parallel::detectCores(); if (is.na(cpus)) cpus <- 1L
  target <- if (requested=='auto') min(8L,max(1L,cpus-2L)) else as.numeric(requested)
  if (!is.finite(target) || target<1) stop('Invalid PheWAS worker count')
  sizes <- vapply(context$specs,function(s) as.double(object.size(s$x)),numeric(1))
  n <- max(vapply(context$specs,function(s) nrow(s$x),integer(1)))
  per_worker <- sum(sizes)+10*max(sizes)+n*8*50+256*1024^2
  memory_workers <- floor(bm_memory_available()/per_worker)
  if (memory_workers<1L) stop('Insufficient free memory for one PheWAS model; retry after other jobs finish')
  as.integer(min(target,cpus,memory_workers,length(context$panel$phenotypes)))
}
# Reuse a prior cache identity only for an explicitly verified, execution-only
# revision. Exact source hashes fail closed on every subsequent source edit;
# data, model settings, package versions and other code hashes remain bound.
bm_checkpoint_signature <- function(signature, compatibility) {
  if (!file.exists(compatibility)) return(signature)
  aliases <- read.delim(compatibility,stringsAsFactors=FALSE,colClasses='character')
  columns <- c('engine_before','entry_before','engine_after','entry_after')
  if (!all(columns %in% names(aliases)) ||
      anyNA(aliases[,columns]) ||
      any(!grepl('^[0-9a-f]{32}$',unlist(aliases[,columns]))))
    stop('Invalid PheWAS runtime compatibility manifest: ',compatibility)
  if (length(signature$engine)<2L) stop('Missing PheWAS code signature')
  hit <- which(aliases$engine_after==signature$engine[1] &
               aliases$entry_after==signature$engine[2])
  if (length(hit)>1L) stop('Ambiguous PheWAS runtime compatibility manifest')
  if (length(hit)) {
    signature$engine[1:2] <- c(aliases$engine_before[hit],aliases$entry_before[hit])
    message('[PheWAS] Reusing model-compatible checkpoints with updated batch scheduling')
  }
  signature
}

# Keep a bounded number of short-lived workers, replacing each worker as soon
# as its checkpoint finishes. A slow batch no longer stalls the other slots.
# Results are assembled later in phenotype order, independently of completion
# order. Fresh forks retain the previous per-batch memory reclamation behavior.
bm_schedule_batches <- function(tasks, workers, on_complete) {
  accept <- function(value, task) {
    if (inherits(value,'try-error') || !is.list(value) ||
        !identical(value$count,length(task$indices)))
      stop('PheWAS worker failed; completed checkpoints are retained. ',
           task$path,': ',paste(value,collapse='; '))
    on_complete(value)
  }
  if (workers==1L) {
    for (task in tasks) accept(bm_batch(task),task)
    return(invisible(NULL))
  }
  active <- list(); next_task <- 1L
  on.exit({
    if (length(active)) {
      # Only children created by this scheduler are terminated on an error or
      # interrupt. Published checkpoints are atomic and remain resumable.
      for (job in active) try(tools::pskill(job$pid,signal=9L),silent=TRUE)
      suppressWarnings(try(parallel::mccollect(active,wait=TRUE),silent=TRUE))
    }
  },add=TRUE)
  invisible(gc(FALSE))
  while (next_task<=length(tasks) || length(active)) {
    while (next_task<=length(tasks) && length(active)<workers) {
      id <- as.character(next_task)
      active[[id]] <- parallel::mcparallel(bm_batch(tasks[[next_task]]),
        name=id,mc.set.seed=FALSE)
      next_task <- next_task+1L
    }
    completed <- parallel::mccollect(active,wait=FALSE,timeout=.2)
    active[names(completed)] <- NULL
    for (id in names(completed)) {
      accept(completed[[id]],tasks[[as.integer(id)]])
    }
  }
  invisible(NULL)
}

bm_run_group <- function(context, directory, requested='auto', replace=FALSE,
                         engine, batch_size=8L) {
  if (replace && dir.exists(directory)) unlink(directory,recursive=TRUE)
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  n <- length(context$panel$phenotypes)
  ids <- split(seq_len(n),ceiling(seq_len(n)/batch_size))
  tasks <- lapply(seq_along(ids),function(i)
    list(indices=unname(ids[[i]]),path=file.path(directory,sprintf('batch-%05d.rds',i))))
  valid <- function(task) {
    value <- bm_read_cache(task$path)
    is.list(value) && identical(value$indices,task$indices) && length(value$results)==length(task$indices)
  }
  cached <- vapply(tasks,valid,logical(1))
  done <- sum(lengths(ids[cached])); initial <- done
  label <- paste(names(context$specs),collapse='+')
  message('[PheWAS] START ',label,': ',done,'/',n,' phenotypes restored from checkpoints')
  pending <- tasks[!cached]
  if (length(pending)) {
    workers <- min(bm_worker_count(context,requested),length(pending))
    message('[PheWAS] START ',label,': ',workers,' workers; ',batch_size,
            ' phenotypes/checkpoint; dynamic scheduling')
    if (workers>1L && .Platform$OS.type=='windows')
      stop('Parallel PheWAS requires Linux/WSL; use one worker on Windows')
    assign('bm_worker_context',context,envir=.GlobalEnv)
    on.exit({if (exists('bm_worker_context',envir=.GlobalEnv,inherits=FALSE))
      rm('bm_worker_context',envir=.GlobalEnv)},add=TRUE)
    start <- proc.time()[['elapsed']]
    bm_schedule_batches(pending,workers,function(completed) {
      done <<- done+completed$count
      elapsed <- proc.time()[['elapsed']]-start
      eta <- elapsed/max(1,done-initial)*(n-done)/60
      message(sprintf('[PheWAS] DONE %s: %d/%d (%.1f%%), elapsed %.1f min, remaining ~%.1f min',
                      label,done,n,100*done/n,elapsed/60,eta))
    })
  }
  rows <- lapply(tasks,function(task) {
    if (!valid(task)) stop('Missing or invalid checkpoint: ',task$path)
    value <- readRDS(task$path)
    do.call(rbind,lapply(value$results,function(r) do.call(rbind,r)))
  })
  result <- do.call(rbind,rows); rownames(result) <- NULL; result
}
