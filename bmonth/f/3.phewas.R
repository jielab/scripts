source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
suppressPackageStartupMessages(library(PheWAS))
engine <- file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R")
source(engine)

run_bmonth_phewas <- function() {
  input <- phewas_input_path(validate=FALSE)
  provenance <- phewas_provenance(input)
  if(!requireNamespace("brglm2",quietly=TRUE))
    stop("Install brglm2 before PheWAS: all unstable eligible fits use the prespecified AS_mean fallback")
  cohort_file <- file.path(rawdir,"bmonth.cohort.rds")
  if(!file.exists(cohort_file))stop("Run bmonth.sh prep first")
  if(!file.exists(input))stop("Missing ",input)

  signature <- list(
    format=6L,
    fit_settings=BM_FIT_SETTINGS,engine_version=BM_ENGINE_VERSION,
    phenotype_provenance=provenance,
    inputs=c(unname(tools::md5sum(cohort_file)),provenance$panel_md5),
    engine=unname(tools::md5sum(c(engine,file.path(Sys.getenv("BMONTH_FDIR"),"3.phewas.R"),
                                  file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"),
                                  file.path(Sys.getenv("BMONTH_FDIR"),"0.analysis_contract.R")))),
    R=R.version.string,PheWAS=as.character(packageVersion("PheWAS")),
    brglm2=as.character(packageVersion("brglm2")),
    primary_exposures=EXPOSURES,
    context_exposures=CONTEXT_EXPOSURES,
    exploratory_traits=EXPLORATORY_TRAITS,
    family_size_sensitivity_traits=FAMILY_SIZE_SENS_TRAITS,
    primary_covariates=covars,
    context_covariates=context_covars,
    exploratory_covariates=exploratory_covars,
    exploratory_family_covariates=exploratory_family_covars,
    contrasts=getOption("contrasts")
  )
  runtime_code <- signature$engine
  signature <- bm_checkpoint_signature(signature,
    file.path(Sys.getenv("BMONTH_FDIR"),"..","config","phewas_runtime_compat.tsv"))
  cache<-file.path(rawdir,"phewas_cache_v6",bm_hash(signature))
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  bm_atomic_rds(signature,file.path(cache,"signature.rds"))
  bm_atomic_rds(list(code_md5=runtime_code,scheduler="dynamic_batches_v1"),
                file.path(cache,"runtime_signature.rds"))

  panel<-phewas_packed_panel(provenance)
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
  diagnostics<-list()
  publish_diagnostics<-function(r,kind){
    q<-r %>% mutate(test_kind=kind) %>%
      group_by(test_kind,exposure,covariate_model,status,method,eligible) %>%
      summarise(n=n(),.groups="drop")
    diagnostics[[paste(kind,unique(r$exposure),unique(r$covariate_model),sep="|")]]<<-q
    write_csv2(bind_rows(diagnostics),"phewas_fit_status.csv")
  }
  dictionary<-annotate_phecodes(data.frame(phenotype=panel$phenotypes))
  xmap<-make_phewas_xmap(dictionary)

  publish<-function(rows,kind,model_label){
    r<-annotate_phecodes(rows)
    r$covariate_model <- model_label
    r$analysis_track <- analysis_track(r$exposure)
    if(kind=="scalar") r$exposure_unit <- exposure_unit(r$exposure)
    r<-bm_adjust_family(r)
    r$provenance_status<-if(!is.null(provenance$validation_status))
      as.character(provenance$validation_status) else "unknown"
    r$engine_version<-BM_ENGINE_VERSION
    r$phenotype_panel_md5<-provenance$panel_md5
    r$cohort_md5<-signature$inputs[1]
    publish_diagnostics(r,kind)
    if(kind=="joint") r$beta<-NA_real_ # joint tests do not have one effect direction
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
    plot_phewas_figures(scalar,cohort=cohort,phenotype_ids=panel$ids,xmap=xmap,only=name,provenance=provenance,cohort_md5=signature$inputs[1])
    message("[PheWAS] DONE ",name)
    rm(context,rows,prepared);invisible(gc())
  }

  # Family-size sensitivity for birth-order PheWAS.
  # Kept separate from primary PheWAS results and significance counts.
  fam_all <- scalar[0,]
  if(length(FAMILY_SIZE_SENS_TRAITS)){
    prepared_fam <- get_prepared(exploratory_family_covars)
    for(name in FAMILY_SIZE_SENS_TRAITS){
      message("[PheWAS-family-size] START ",name,
              " | D2 + categorical sibship size | N=",nrow(prepared_fam$dat))
      context <- bm_context(prepared_fam,name,panel,"scalar")
      rows <- bm_run_group(context,file.path(cache,paste0("family_size_",name)),
                           requested,REPLACE,engine)
      r <- publish(rows,"scalar","exploratory_D2_plus_family_size")
      fam_all <- bind_rows(fam_all,r)
      bm_atomic_rds(r,file.path(rawdir,paste0("phewas_family_size_",name,".rds")))
      rm(context,rows); invisible(gc())
    }
  }
  # Always publish the current selection, including an empty table when disabled.
  bm_atomic_rds(fam_all,file.path(rawdir,"phewas_family_size_sensitivity.rds"))
  write_csv2(fam_all,"phewas_family_size_sensitivity.csv")

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
    plot_phewas_figures(scalar,joint,cohort=cohort,phenotype_ids=panel$ids,xmap=xmap,only=names,provenance=provenance,cohort_md5=signature$inputs[1])
    message("[PheWAS] DONE ",paste(names,collapse=", "))
    rm(context,rows);invisible(gc())
  }
}
run_bmonth_phewas()
