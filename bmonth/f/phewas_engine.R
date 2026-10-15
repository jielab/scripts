# Each worker holds a design matrix and one phenotype, never the full panel.
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
bm_panel <- function(input, directory) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  manifest <- file.path(directory,'panel.rds'); binary <- file.path(directory,'panel.bin')
  meta <- bm_read_cache(manifest)
  if (!is.null(meta) && file.exists(binary) &&
      isTRUE(file.info(binary)$size == as.double(meta$n)*length(meta$phenotypes))) {
    meta$path <- binary; message('[PheWAS] START Reusing packed phenotype panel'); return(meta)
  }
  ph <- readRDS(input)
  if (!is.data.frame(ph) || !'id' %in% names(ph) || ncol(ph)<2L)
    stop('PheWAS.rds must contain id and phenotype columns')
  ids <- as.character(ph$id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || anyDuplicated(names(ph)))
    stop('Invalid or duplicate phenotype IDs/column names')
  ys <- setdiff(names(ph),'id'); tmp <- tempfile('.panel-',tmpdir=directory)
  con <- file(tmp,'wb')
  on.exit({if (!is.null(con)) close(con); unlink(tmp)},add=TRUE)
  for (i in seq_along(ys)) {
    y <- ph[[ys[i]]]
    if (!is.logical(y)) {
      if (!is.numeric(y)) y <- suppressWarnings(as.numeric(as.character(y)))
      y <- y==1
    }
    # Column-major bytes: 0=control, 1=case, 2=missing (never turn NA into control).
    packed <- as.integer(y); packed[is.na(packed)] <- 2L
    writeBin(as.raw(packed),con)
    if (i %% 200L==0L || i==length(ys))
      message('[PheWAS] DONE Packing phenotypes ',i,'/',length(ys))
  }
  close(con); con <- NULL
  if (!file.rename(tmp,binary)) stop('Cannot publish phenotype panel')
  meta <- list(n=nrow(ph),ids=ids,phenotypes=ys,path=binary)
  bm_atomic_rds(meta,manifest); meta
}
bm_prepare <- function(dat, panel, covariates) {
  idx <- match(as.character(dat$eid),panel$ids)
  keep <- !is.na(idx) & complete.cases(dat[,covariates,drop=FALSE])
  for (v in covariates) if (is.numeric(dat[[v]])) keep <- keep & is.finite(dat[[v]])
  d <- droplevels(dat[keep,,drop=FALSE])
  if (!nrow(d)) stop('No participants with complete covariates and phenotype IDs')
  x <- model.matrix(reformulate(covariates),d); rownames(x) <- NULL
  list(dat=d,rows=idx[keep],covariates=x)
}
bm_context <- function(prepared, names, panel, kind=c('scalar','joint')) {
  kind <- match.arg(kind); dat <- prepared$dat; specs <- list()
  for (name in names) {
    extra <- if (kind=='scalar') name else switch(name,
      month11df='birth_month', cycle2df=c('bmonth_cos','bmonth_sin'),
      north_interaction2df=c('birth_north_100km','bmonth_cos','bmonth_sin',
                            'bmonth_cos_x_north','bmonth_sin_x_north'),
      east_interaction2df=c('birth_east_100km','bmonth_cos','bmonth_sin',
                           'bmonth_cos_x_east','bmonth_sin_x_east'),
      stop('Unknown joint test: ',name))
    keep <- complete.cases(dat[,extra,drop=FALSE])
    for (v in extra) keep <- keep & is.finite(dat[[v]])
    base <- prepared$covariates[keep,,drop=FALSE]
    if (kind=='scalar') {
      x <- cbind(base[,1,drop=FALSE],dat[[name]][keep],base[,-1,drop=FALSE])
      colnames(x)[2] <- name; null_cols <- integer()
    } else {
      if (name=='month11df') {
        m <- factor(dat$birth_month[keep],levels=1:12)
        extra_x <- model.matrix(~m)[,-1,drop=FALSE]
      } else extra_x <- as.matrix(dat[keep,extra,drop=FALSE])
      x <- cbind(base,extra_x)
      null_cols <- seq_len(ncol(base)+if (grepl('interaction',name)) 3L else 0L)
    }
    rownames(x) <- NULL
    specs[[name]] <- list(name=name,kind=kind,rows=prepared$rows[keep],x=x,null_cols=null_cols)
  }
  # Reuse a null fit only when participants and covariate matrices are identical.
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
# Same IRLS/QR solver and defaults as stats::glm(), without rebuilding formulas
# and model matrices. Retain only summaries, not two large QR objects for an LRT.
bm_fit <- function(x, y, coefficient=NULL, strict=FALSE) {
  notes <- character()
  fit <- tryCatch(withCallingHandlers(glm.fit(x=x,y=y,family=binomial()),
    warning=function(w) {notes <<- c(notes,conditionMessage(w)); invokeRestart('muffleWarning')}),
    error=function(e) {notes <<- c(notes,conditionMessage(e)); NULL})
  ok <- !is.null(fit) && isTRUE(fit$converged)
  if (strict && any(grepl('converg|infinite|singular',notes,ignore.case=TRUE))) ok <- FALSE
  out <- list(ok=ok,note=paste(unique(notes),collapse='; '))
  if (!ok) return(out)
  out$deviance <- fit$deviance; out$rank <- fit$rank
  if (!is.null(coefficient)) {
    sm <- summary.glm(fit)$coefficients
    if (coefficient %in% rownames(sm)) out$coef <- unname(sm[coefficient,])
  }
  out
}
bm_scalar <- function(spec, y, phenotype) {
  y <- y[spec$rows]; keep <- !is.na(y); y <- y[keep]
  n <- length(y); cases <- sum(y)
  r <- data.frame(phenotype=phenotype,snp=spec$name,adjustment=NA_character_,
    beta=NA_real_,SE=NA_real_,OR=NA_real_,p=NA_real_,type='logistic',
    n_total=n,n_cases=cases,n_controls=n-cases,allele_freq=NA_real_,HWE_p=NA_real_,
    n_no_snp=0L,note='',stringsAsFactors=FALSE)
  if (n<20L || cases<20L || n-cases<20L) {
    r$note <- '[Error: < 20 cases or controls]'; return(r)
  }
  x <- if (all(keep)) spec$x else spec$x[keep,,drop=FALSE]
  if (length(unique(x[,2]))<2L) {
    r$note <- '[Error: non-varying phenotype or genotype]'; return(r)
  }
  fit <- bm_fit(x,y,coefficient=spec$name); r$note <- fit$note
  if (fit$ok && !is.null(fit$coef)) {
    r$beta <- fit$coef[1]; r$SE <- fit$coef[2]; r$OR <- exp(r$beta); r$p <- fit$coef[4]
  } else if (!nzchar(r$note)) r$note <- '[Error: model did not converge or coefficient is aliased]'
  r
}
bm_joint <- function(spec, y, phenotype, null_fits) {
  y <- y[spec$rows]; keep <- !is.na(y); y <- y[keep]
  r <- data.frame(phenotype=phenotype,N=length(y),cases=sum(y),p=NA_real_)
  if (r$N<500L || r$cases<20L || r$N-r$cases<20L) return(r)
  x <- if (all(keep)) spec$x else spec$x[keep,,drop=FALSE]
  key <- spec$null_key
  if (!exists(key,envir=null_fits,inherits=FALSE))
    assign(key,bm_fit(x[,spec$null_cols,drop=FALSE],y,strict=TRUE),envir=null_fits)
  a <- get(key,envir=null_fits,inherits=FALSE)
  if (!a$ok) return(r)
  b <- bm_fit(x,y,strict=TRUE)
  if (!b$ok) return(r)
  df <- b$rank-a$rank; dev <- a$deviance-b$deviance
  if (df>0 && is.finite(dev) && dev >= -1e-7)
    r$p <- pchisq(max(0,dev),df=df,lower.tail=FALSE)
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
  # Conservative allowance for glm.fit work arrays and the worker runtime.
  per_worker <- sum(sizes)+10*max(sizes)+n*8*50+256*1024^2
  memory_workers <- floor(bm_memory_available()/per_worker)
  if (memory_workers<1L) stop('Insufficient free memory for one PheWAS model; retry after other jobs finish')
  as.integer(min(target,cpus,memory_workers,length(context$panel$phenotypes)))
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
    message('[PheWAS] START ',label,': ',workers,' workers; ',batch_size,' phenotypes/checkpoint')
    if (workers>1L && .Platform$OS.type=='windows')
      stop('Parallel PheWAS requires Linux/WSL; use one worker on Windows')
    # Fork only after the full phenotype table has been released. Read-only
    # matrices share physical pages, and local pipes avoid TCP/WSL idle stalls.
    assign('bm_worker_context',context,envir=.GlobalEnv)
    on.exit({if (exists('bm_worker_context',envir=.GlobalEnv,inherits=FALSE))
      rm('bm_worker_context',envir=.GlobalEnv)},add=TRUE)
    start <- proc.time()[['elapsed']]
    waves <- split(seq_along(pending),ceiling(seq_along(pending)/workers))
    for (wave in waves) {
      invisible(gc(FALSE))
      completed <- if (workers==1L) lapply(pending[wave],bm_batch) else
        parallel::mclapply(pending[wave],bm_batch,mc.cores=workers,
          mc.preschedule=FALSE,mc.set.seed=FALSE,mc.cleanup=TRUE,mc.allow.recursive=FALSE)
      failed <- vapply(completed,function(x) inherits(x,'try-error') ||
        !is.list(x) || is.null(x$count),logical(1))
      if (any(failed)) stop('PheWAS worker failed; completed checkpoints are retained. ',
                            paste(completed[failed],collapse='; '))
      done <- done+sum(vapply(completed,`[[`,integer(1),'count'))
      elapsed <- proc.time()[['elapsed']]-start
      eta <- elapsed/max(1,done-initial)*(n-done)/60
      message(sprintf('[PheWAS] DONE %s: %d/%d (%.1f%%), elapsed %.1f min, remaining ~%.1f min',
                      label,done,n,100*done/n,elapsed/60,eta))
    }
  }
  rows <- lapply(tasks,function(task) {
    if (!valid(task)) stop('Missing or invalid checkpoint: ',task$path)
    value <- readRDS(task$path)
    do.call(rbind,lapply(value$results,function(r) do.call(rbind,r)))
  })
  result <- do.call(rbind,rows); rownames(result) <- NULL; result
}
