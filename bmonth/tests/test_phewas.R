# Compare optimized fits with the original PheWAS/glm algorithms, including
# sample mismatches, disease missingness, absent factor levels and resume.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1]); work <- normalizePath(args[2])
Sys.setenv(BMONTH_OUTDIR=file.path(work,'output'),BMONTH_INDIR=file.path(work,'input'))
source(file.path(root,'f','0.common.R'))
source(file.path(root,'f','phewas_engine.R'))
suppressPackageStartupMessages(library(PheWAS))
dat <- load_cohort(); ph <- readRDS(file.path(indir,'Rdata','PheWAS.rds'))
# Shuffle IDs, remove a cohort member, add a phenotype-only member, and drop a
# whole centre from one phenotype. No participant IDs are exposed in output.
ph <- ph[c(nrow(ph):2,1),]; ph <- ph[ph$id != dat$eid[100],]
ph <- rbind(ph,ph[1,]); ph$id[nrow(ph)] <- 'not-in-cohort'
ph[['714']][ph$id %in% dat$eid[dat$center=='B']] <- NA
# Numeric phenotype encoding must agree with the original conversion.
ph[['493']] <- as.integer(ph[['493']])
panel_input <- file.path(work,'test-panel.rds'); saveRDS(ph,panel_input)
panel <- bm_panel(panel_input,file.path(work,'test-panel'))
prepared <- bm_prepare(dat,panel,covars)
engine <- file.path(root,'f','phewas_engine.R')
passed <- 0L
check <- function(label,expr) {
  stopifnot(isTRUE(expr)); passed <<- passed+1L; cat('PASS:',label,'\n')
}
close_to <- function(a,b,tol=1e-8) isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol))
ph[-1] <- lapply(ph[-1],function(y) if(is.logical(y)) y else as_num(y)==1)
for (exposure in EXPOSURES) {
  context <- bm_context(prepared,exposure,panel,'scalar')
  folder <- file.path(work,'engine-test',exposure)
  # Small blocks exercise multiple processes on seven traits.
  actual <- bm_run_group(context,folder,'2',TRUE,engine,batch_size=2L)
  d <- complete_model_data(dat,c('eid',exposure,covars))
  g <- d[,c('eid',exposure)]; names(g)[1] <- 'id'
  cv <- d[,c('eid',covars)]; names(cv)[1] <- 'id'
  expected <- PheWAS::phewas(phenotypes=ph[ph$id %in% d$eid,],genotypes=g,covariates=cv,
    cores=1,additive.genotypes=FALSE,significance.threshold='bonferroni')
  expected <- expected[match(actual$phenotype,expected$phenotype),]
  check(paste(exposure,'beta/SE/OR/Wald P match original'),all(vapply(
    c('beta','SE','OR','p'),function(v)close_to(actual[[v]],expected[[v]]),logical(1))))
  # Original PheWAS leaves n_cases NA for constant phenotypes; verify tested ones.
  tested <- is.finite(expected$p)
  check(paste(exposure,'sample counts match'),all(actual$n_total==expected$n_total) &&
    all(actual$n_cases[tested]==expected$n_cases[tested]))
  files <- list.files(folder,full.names=TRUE); stamps <- file.info(files)$mtime
  resumed <- bm_run_group(context,folder,'1',FALSE,engine,batch_size=2L)
  check(paste(exposure,'resume preserves results and files'),identical(actual,resumed) &&
    identical(file.info(files)$mtime,stamps))
  if (exposure==EXPOSURES[1]) {
    serial <- bm_run_group(context,paste0(folder,'-serial'),'1',TRUE,engine,batch_size=2L)
    check('serial and two-worker results identical',identical(actual,serial))
    writeLines('interrupted write',files[1])
    repaired <- bm_run_group(context,folder,'2',FALSE,engine,batch_size=2L)
    check('corrupted checkpoint recomputed',identical(actual,repaired))
  }
}
reference_joint <- function(yname,test) {
  ytab <- ph[,c('id',yname)]; names(ytab) <- c('eid','.y')
  z <- inner_join(dat,ytab,by='eid'); terms0 <- covars
  add <- switch(test,month11df='.month',cycle2df=c('bmonth_cos','bmonth_sin'),
    north_interaction2df=c('bmonth_cos_x_north','bmonth_sin_x_north'),
    east_interaction2df=c('bmonth_cos_x_east','bmonth_sin_x_east'))
  if (test=='month11df') {
    z <- complete_model_data(z,c('.y','birth_month',covars)); z$.month <- factor(z$birth_month,levels=1:12)
  } else {
    if (test=='north_interaction2df') terms0 <- c(covars,'birth_north_100km','bmonth_cos','bmonth_sin')
    if (test=='east_interaction2df') terms0 <- c(covars,'birth_east_100km','bmonth_cos','bmonth_sin')
    z <- complete_model_data(z,c('.y',terms0,add))
  }
  r <- data.frame(phenotype=yname,N=nrow(z),cases=sum(z$.y),p=NA_real_)
  if (r$N<500 || r$cases<20 || r$N-r$cases<20) return(r)
  a <- model_try(glm(reformulate(terms0,response='.y'),z,family=binomial()))
  b <- model_try(glm(reformulate(c(terms0,add),response='.y'),z,family=binomial()))
  if (!is.null(a$fit) && !is.null(b$fit)) r$p <- anova(a$fit,b$fit,test='LRT')$`Pr(>Chi)`[2]
  r
}
groups <- list(c('month11df','cycle2df'),'north_interaction2df','east_interaction2df')
for (group in groups) {
  context <- bm_context(prepared,group,panel,'joint')
  if (length(group)==2L) check('calendar tests share identical null fits',
    identical(context$specs[[1]]$null_key,context$specs[[2]]$null_key))
  actual <- bm_run_group(context,file.path(work,'engine-test',paste(group,collapse='_')),
                         '2',TRUE,engine,batch_size=2L)
  for (test in group) {
    a <- actual[actual$exposure==test,]
    expected <- do.call(rbind,lapply(a$phenotype,reference_joint,test=test))
    check(paste(test,'LRT P and sample counts match glm/anova'),
      close_to(a$p,expected$p) && identical(a$N,expected$N) && close_to(a$cases,expected$cases))
  }
}
# An exposure with no eligible samples should produce skipped tests, not crash.
empty <- prepared; empty$dat$bmonth_cos <- NA_real_
ctx <- bm_context(empty,'bmonth_cos',panel,'scalar')
check('all-missing exposure skipped',is.na(bm_scalar(ctx$specs[[1]],rep(0,panel$n),'empty')$p))
# Worker I/O failures must reach the parent, rather than becoming successful NA rows.
bad <- bm_context(prepared,'bmonth_cos',panel,'scalar')
bad$panel$path <- file.path(work,'does-not-exist.bin')
failure <- tryCatch(bm_run_group(bad,file.path(work,'broken-worker'),'2',TRUE,engine,
                                batch_size=2L),error=identity)
check('worker failure propagates without publishing results',inherits(failure,'error') &&
      length(list.files(file.path(work,'broken-worker'),pattern='[.]rds$'))==0L)
# Changing a source or input fingerprint must select a distinct checkpoint key.
check('cache keys reflect input and model changes',bm_hash(list(input='a',code='a')) !=
      bm_hash(list(input='b',code='a')) && bm_hash(list(input='a',code='a')) != bm_hash(list(input='a',code='b')))
cat('ALL',passed,'PHEWAS ENGINE CHECKS PASSED\n')
