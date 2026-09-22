# Explicit, trait-independent settings. Disk phenotype data remain unchanged.
le8_custom_covars<-unique(Filter(nzchar,trimws(strsplit(Sys.getenv('LE8_VARS_ADJ',''),'[,[:space:]]+')[[1]])))
le8_custom_adjustment<-function()length(le8_custom_covars)>0L
le8_y_date<-function(outcome=Y)Sys.getenv('LE8_Y_DATE',paste0('fod_icd10_',outcome))
le8_analysis_options<-function()list(y_date=le8_y_date(),vars_adj=le8_custom_covars,
  white_only=Sys.getenv('LE8_WHITE_ONLY','TRUE'),baseline_contract=get0('LE8_BASELINE_VERSION',ifnotfound='2026-09-21.baseline-med-source-v3'),
  baseline_med_columns=Sys.getenv('LE8_BASELINE_MED_COLUMNS',''),
  baseline_map=if(nzchar(Sys.getenv('LE8_BASELINE_MAP','')))tools::md5sum(Sys.getenv('LE8_BASELINE_MAP'))else '')
le8_check_options<-function(obj) {
  current<-le8_analysis_options();old<-obj$meta$analysis_options
  default<-list(y_date=paste0('fod_icd10_',Y),vars_adj=character(),white_only='TRUE')
  if(is.null(old))old<-default
  if(!identical(old,current))stop('Existing results use different analysis options. Use --replace TRUE or a new output directory; changing options must not silently relabel cached estimates.')
}
le8_select_phenotypes<-function(x) {
  source<-le8_y_date();required<-unique(c(source,le8_custom_covars))
  if(length(setdiff(required,names(x))))stop('Unknown --Y-date/--vars.adj fields: ',paste(setdiff(required,names(x)),collapse=', '))
  # Alias only inside this analysis, so all existing consumers use the selected date.
  x[[paste0('fod_icd10_',Y)]]<-x[[source]]
  x
}
if(le8_custom_adjustment()) {
  covs_use<-le8_custom_covars;covs_use_name<-'adj2'
}

le8_guard_stage_options<-function(layer,module) {
  outdir<-if(layer=='protein')out.prot else out.met
  d<-le8_job_dir(outdir,module);dir.create(d,recursive=TRUE,showWarnings=FALSE)
  f<-file.path(d,'analysis_options.rds');current<-le8_analysis_options()
  if(file.exists(f)&&!LE8_REPLACE&&!identical(readRDS(f),current))stop('Cached stage options differ. Use --replace TRUE or a new output directory.')
  if(!file.exists(f)&&length(list.files(d,pattern='[.]rds$'))&&!LE8_REPLACE){
    le8_check_options(list(meta=list()))
  }
  saveRDS(current,f)
}
