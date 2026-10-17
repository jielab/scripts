source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
suppressPackageStartupMessages(library(PheWAS))
engine <- file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R")
source(engine)

run_bmonth_phewas <- function() {
  input <- file.path(indir,"Rdata","PheWAS.rds")
  cohort_file <- file.path(rawdir,"bmonth.cohort.rds")
  if(!file.exists(cohort_file))stop("Run bmonth.sh prep first")
  if(!file.exists(input))stop("Missing ",input)

  signature <- list(
    format=4L,
    inputs=unname(tools::md5sum(c(cohort_file,input))),
    engine=unname(tools::md5sum(c(engine,file.path(Sys.getenv("BMONTH_FDIR"),"3.phewas.R"),
                                  file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R")))),
    R=R.version.string,PheWAS=as.character(packageVersion("PheWAS")),
    primary_exposures=EXPOSURES,
    context_exposures=CONTEXT_EXPOSURES,
    exploratory_traits=EXPLORATORY_TRAITS,
    primary_covariates=covars,
    context_covariates=context_covars,
    exploratory_covariates=exploratory_covars,
    contrasts=getOption("contrasts")
  )
  cache<-file.path(rawdir,"phewas_cache_v4",bm_hash(signature))
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  bm_atomic_rds(signature,file.path(cache,"signature.rds"))

  panel<-bm_panel(input,cache)
  invisible(gc())
  cohort<-load_cohort()

  # Prepare each covariate design only once. Context geography models are NOT
  # forced onto the smaller primary complete-case sample.
  prep_env <- new.env(parent=emptyenv())
  get_prepared <- function(cv){
    key <- paste(cv,collapse="|")
    if(!exists(key,envir=prep_env,inherits=FALSE))
      assign(key,bm_prepare(cohort,panel,cv),envir=prep_env)
    get(key,envir=prep_env,inherits=FALSE)
  }

  requested<-Sys.getenv("BMONTH_PHEWAS_CORES","auto")
  message("[PheWAS] primary exposures: ",paste(EXPOSURES,collapse=", "))
  message("[PheWAS] context exposures: ",paste(CONTEXT_EXPOSURES,collapse=", "))
  message("[PheWAS] exploratory birth-order traits: ",paste(EXPLORATORY_TRAITS,collapse=", "))
  message("[PheWAS] primary covariates: ",paste(covars,collapse=", "))
  message("[PheWAS] context covariates (OLD-compatible): ",paste(context_covars,collapse=", "))
  message("[PheWAS] checkpoints: ",cache)

  source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_plot.R"))
  scalar<-joint<-NULL
  dictionary<-annotate_phecodes(data.frame(phenotype=panel$phenotypes))
  xmap<-make_phewas_xmap(dictionary)

  publish<-function(rows,kind,model_label){
    r<-annotate_phecodes(rows)
    r$covariate_model <- model_label
    r$analysis_track <- analysis_track(r$exposure)
    if(kind=="scalar") r$exposure_unit <- exposure_unit(r$exposure)
    r$FDR<-p.adjust(r$p,"BH")
    r$bonferroni_threshold<-.05/max(1,sum(is.finite(r$p)))
    if(kind=="scalar"){
      r$p_raw<-r$p
      r$logp<--log10(pmax(r$p,.Machine$double.xmin))
      r$bonferroni<-r$p<=r$bonferroni_threshold
    } else r$beta<-0
    r
  }

  for(name in PHEWAS_EXPOSURES){
    cv <- phewas_covars_for(name)
    prepared <- get_prepared(cv)
    model_label <- phewas_model_label(name)
    message("[PheWAS] START ",name," | ",model_label,
            " | matched complete covariate N=",nrow(prepared$dat))
    context<-bm_context(prepared,name,panel,"scalar")
    rows<-bm_run_group(context,file.path(cache,paste0("scalar_",name)),
                       requested,REPLACE,engine)
    r<-publish(rows,"scalar",model_label)
    bm_atomic_rds(r,file.path(rawdir,paste0("phewas_",name,".rds")))
    scalar<-bind_rows(scalar,r)
    bm_atomic_rds(scalar,file.path(rawdir,"phewas.rds"));write_csv2(scalar,"phewas.csv")
    plot_phewas_figures(scalar,cohort=cohort,phenotype_ids=panel$ids,xmap=xmap,only=name)
    message("[PheWAS] DONE ",name)
    rm(context,rows,prepared);invisible(gc())
  }

  # Joint birth-month / interaction tests use the primary birth-context design.
  prepared_primary <- get_prepared(covars)
  groups<-list(
    c("month11df","cycle2df"),
    "north_interaction2df",
    "east_interaction2df",
    "move_interaction2df"
  )
  for(names in groups){
    context<-bm_context(prepared_primary,names,panel,"joint")
    rows<-bm_run_group(context,file.path(cache,paste(names,collapse="_")),
                       requested,REPLACE,engine)
    for(name in names)
      joint<-bind_rows(joint,publish(rows[rows$exposure==name,,drop=FALSE],
                                    "joint","primary_D4_birthContext"))
    bm_atomic_rds(joint,file.path(rawdir,"phewas_joint.rds"));write_csv2(joint,"phewas_joint.csv")
    plot_phewas_figures(scalar,joint,cohort=cohort,phenotype_ids=panel$ids,xmap=xmap,only=names)
    message("[PheWAS] DONE ",paste(names,collapse=", "))
    rm(context,rows);invisible(gc())
  }
}
run_bmonth_phewas()
