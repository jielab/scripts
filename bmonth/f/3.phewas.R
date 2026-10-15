source(file.path(Sys.getenv('BMONTH_FDIR'),'0.common.R'))
suppressPackageStartupMessages(library(PheWAS))
engine <- file.path(Sys.getenv('BMONTH_FDIR'),'phewas_engine.R')
source(engine)

run_bmonth_phewas <- function() {
  input <- file.path(indir,'Rdata','PheWAS.rds')
  cohort_file <- file.path(rawdir,'bmonth.cohort.rds')
  if (!file.exists(cohort_file)) stop('Run bmonth.sh prep first')
  if (!file.exists(input)) stop('Missing ',input)
  signature <- list(format=2L,inputs=unname(tools::md5sum(c(cohort_file,input))),
    engine=unname(tools::md5sum(c(engine,file.path(Sys.getenv('BMONTH_FDIR'),'3.phewas.R'),
                               file.path(Sys.getenv('BMONTH_FDIR'),'0.common.R')))),
    R=R.version.string,PheWAS=as.character(packageVersion('PheWAS')),
    covariates=covars,contrasts=getOption('contrasts'))
  cache <- file.path(rawdir,'phewas_cache_v2',bm_hash(signature))
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  bm_atomic_rds(signature,file.path(cache,'signature.rds'))
  panel <- bm_panel(input,cache)
  invisible(gc()) # Release the full input before creating workers.
  cohort <- load_cohort(); prepared <- bm_prepare(cohort,panel,covars)
  requested <- Sys.getenv('BMONTH_PHEWAS_CORES','auto')
  message('[PheWAS] START ',length(panel$phenotypes),' phenotypes; ',nrow(prepared$dat),
          ' matched participants with complete covariates; requested workers=',requested)
  message('[PheWAS] START Checkpoints: ',cache)
  source(file.path(Sys.getenv('BMONTH_FDIR'),'phewas_plot.R'))
  scalar <- joint <- NULL
  dictionary <- annotate_phecodes(data.frame(phenotype=panel$phenotypes))
  xmap <- make_phewas_xmap(dictionary)
  publish <- function(rows,kind) {
    r <- annotate_phecodes(rows); r$FDR <- p.adjust(r$p,'BH')
    r$bonferroni_threshold <- .05/max(1,sum(is.finite(r$p)))
    if (kind=='scalar') {
      r$p_raw <- r$p; r$logp <- -log10(pmax(r$p,.Machine$double.xmin))
      r$bonferroni <- r$p <= r$bonferroni_threshold
    } else r$beta <- 0
    bm_atomic_rds(r,file.path(rawdir,paste0('phewas_',unique(r$exposure),'.rds')))
    r
  }
  for (name in EXPOSURES) {
    context <- bm_context(prepared,name,panel,'scalar')
    rows <- bm_run_group(context,file.path(cache,paste0('scalar_',name)),requested,REPLACE,engine)
    r <- publish(rows,'scalar'); scalar <- bind_rows(scalar,r)
    bm_atomic_rds(scalar,file.path(rawdir,'phewas.rds')); write_csv2(scalar,'phewas.csv')
    plot_phewas_figures(scalar,cohort=cohort,phenotype_ids=panel$ids,xmap=xmap,only=name)
    message('[PheWAS] DONE Completed ',name,'; results and figures saved')
    rm(context,rows); invisible(gc())
  }
  groups <- list(c('month11df','cycle2df'),'north_interaction2df','east_interaction2df')
  for (names in groups) {
    context <- bm_context(prepared,names,panel,'joint')
    rows <- bm_run_group(context,file.path(cache,paste(names,collapse='_')),requested,REPLACE,engine)
    for (name in names) joint <- bind_rows(joint,publish(rows[rows$exposure==name,,drop=FALSE],'joint'))
    bm_atomic_rds(joint,file.path(rawdir,'phewas_joint.rds')); write_csv2(joint,'phewas_joint.csv')
    plot_phewas_figures(scalar,joint,cohort=cohort,phenotype_ids=panel$ids,xmap=xmap,only=names)
    message('[PheWAS] DONE Completed ',paste(names,collapse=', '),'; results and figures saved')
    rm(context,rows); invisible(gc())
  }
}
run_bmonth_phewas()
