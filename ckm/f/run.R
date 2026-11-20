#!/usr/bin/env Rscript
# CLI orchestration. Scientific modules have no top-level data side effects.
argv<-commandArgs(trailingOnly=TRUE)
if(length(argv)%%2L)stop("Each R CLI option requires a value")
arg<-list()
for(i in seq.int(1L,length(argv),by=2L)) {
  key<-sub("^--","",argv[i]);if(!startsWith(argv[i],"--"))stop("Invalid option ",argv[i])
  arg[[gsub("-","_",key)]]<-argv[i+1L]
}
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
fdir<-dirname(normalizePath(script));root<-dirname(fdir)
for(f in c("shared.R","common.R","risk.R","phenotype.R","clock.R","associations.R","markers.R","demo.R"))source(file.path(fdir,f))
ckm_require()
ckm_load_shared(file.path(root,"..","0f","phenotype.R"))
if(!is.null(arg$phenotype_dir))arg$ukb_dir<-arg$phenotype_dir
cfg<-ckm_config(arg$config %||% file.path(root,"config","ckm.json"),arg)
module<-arg$module %||% "all"
demo<-module=="demo"
if(demo) {cfg<-ckm_make_demo(cfg,as.integer(arg$demo_n %||% 6000L));module<-"all"}
data.table::setDTthreads(cfg$threads)
ph_file<-file.path(cfg$out_dir,"ckm.rds")
if(module=="phe"||demo) {
  message("Building CKM phenotype from ",cfg$paths$phe)
  d<-ckm_build_phenotype(cfg)
  attr(d,"ckm_build_id")<-paste(format(Sys.time(),"%Y%m%dT%H%M%OS6"),Sys.getpid(),sep="-")
  attr(d,"ckm_demo")<-isTRUE(cfg$demo)
  ckm_atomic_rds(d,ph_file)
  summary<-ckm_phenotype_summary(d)
  ckm_write(summary$stage_counts,file.path(cfg$out_dir,"stage.counts.tsv"))
  ckm_write(summary$missingness,file.path(cfg$out_dir,"missingness.tsv"))
} else {
  if(!file.exists(ph_file))stop("Run ./1.phe.sh first; phenotype not found: ",ph_file)
  d<-readRDS(ph_file)
  if(is.null(attr(d,"ckm_build_id")))stop("ckm.rds belongs to a different schema; rebuild with ./1.phe.sh")
}
cfg$demo<-isTRUE(attr(d,"ckm_demo"))
build_id<-attr(d,"ckm_build_id")
ph_meta<-attr(d,"ckm_metadata")
if(!identical(as.character(ph_meta$definition),as.character(cfg$phenotype$definition)) ||
   !identical(as.character(ph_meta$stage3_method),as.character(cfg$phenotype$stage3_method)))
  stop("Configured staging differs from ckm.rds. Rebuild the phenotype first, or use the matching config.")
if(!identical(as.character(ph_meta$date_follow_end),as.character(cfg$shared$date_follow_end)))
  stop("ckm.rds uses a different or unrecorded date_follow_end. Run ./1.phe.sh to rebuild with the current shared cutoff: ",cfg$shared$date_follow_end)
code_hash<-unname(tools::md5sum(c(file.path(fdir,c("shared.R","common.R","risk.R","phenotype.R","clock.R","associations.R","markers.R")),cfg$shared$file)))
cache_settings<-function(name) {
  settings<-switch(name,clock=cfg$clock,assoc=cfg$association,marker=cfg$markers,enrichment=cfg$enrichment,list())
  settings$workers<-NULL
  list(scientific=settings,exclude_met_lipid_block=isTRUE(cfg$exclude_met_lipid_block))
}
save_cache<-function(x,name) {
  obj<-list(phenotype_id=build_id,generated=format(Sys.time(),"%Y-%m-%dT%H:%M:%SZ",tz="UTC"),
      settings=cache_settings(name),domains=cfg$domains,code_hash=code_hash,data=x)
  ckm_atomic_rds(obj,file.path(cfg$out_dir,paste0(name,".rds")));invisible(obj)
}
load_cache<-function(name) {
  path<-file.path(cfg$out_dir,paste0(name,".rds"));if(!file.exists(path))return(NULL)
  x<-readRDS(path)
  if(!identical(x$phenotype_id,build_id)||!identical(x$settings,cache_settings(name))||
     !identical(x$domains,cfg$domains)||!identical(x$code_hash,code_hash))return(NULL)
  x
}
modules<-if(module=="all")c("clock","assoc","marker","enrich") else module
domains<-list()
if(any(modules %in% c("clock","assoc","marker"))) {
  for(nm in cfg$domains) {
    if(nm=="lifestyle" && !"assoc" %in% modules)next
    message("Reading ",nm," from its configured baseline source")
    domains[[nm]]<-ckm_domain(d,cfg,nm)
  }
}
if("clock" %in% modules) {
  result<-ckm_run_clock(d,cfg,domains[intersect(names(domains),c("prot","met"))])
  result$grades<-result$predictions[is.finite(risk)&is.finite(ipcw)&!is.na(grade),
    .(n=.N,events=sum(status==1),predicted=mean(risk),observed=sum(ipcw*target)/.N),
    by=.(cohort,model,stage,grade,horizon)]
  models<-result$models;result$models<-NULL
  # Stored prediction models need coefficients, not copies of participant records.
  for(nm in names(models))if(is.null(models[[nm]]$error)) {
    models[[nm]]$training_n<-length(models[[nm]]$training_eids);models[[nm]]$training_eids<-NULL
    models[[nm]]$prep$x<-NULL;models[[nm]]$death_prep$x<-NULL
    for(cause in c("event","death"))if(!is.null(models[[nm]][[cause]])) {
      ff<-models[[nm]][[cause]]
      # Formula/terms environments can otherwise serialize the complete training frame.
      models[[nm]][[cause]]<-structure(list(coefficients=coef(ff),var=vcov(ff),n=ff$n,
        nevent=ff$nevent,method=ff$method),class="coxph")
    }
  }
  ckm_atomic_rds(list(models=models,phenotype_id=build_id,definition=ph_meta,clock=cfg$clock,
    units="absolute risk in 0..1; position=stage+training-reference percentile",version=cfg$version),file.path(cfg$out_dir,"model.rds"))
  save_cache(result,"clock")
  for(nm in c("metrics","calibration","comparison","grades","selection","fit_status"))
    ckm_write(result[[nm]],file.path(cfg$out_dir,paste0("clock.",nm,".tsv")))
  # clock.rds retains person-level OOF predictions; the HTML receives aggregates only.
}
if("assoc" %in% modules) {
  pieces<-list();run_cfg<-cfg;run_cfg$out_dir<-NULL
  for(nm in names(domains)) {
    message("Associations: ",nm)
    dd<-d[d$eid %in% domains[[nm]]$data$eid,]
    pieces[[nm]]<-ckm_run_associations(dd,run_cfg,domains[[nm]])
  }
  keys<-setdiff(unique(unlist(lapply(pieces,names))),"settings")
  result<-stats::setNames(lapply(keys,function(k)data.table::rbindlist(lapply(pieces,`[[`,k),fill=TRUE)),keys)
  result$settings<-cfg$association
  save_cache(result,"assoc")
  for(nm in keys)ckm_write(result[[nm]],file.path(cfg$out_dir,paste0(nm,".tsv")))
  # Explicitly labelled sensitivity: only participants with known core staging components.
  sens<-list()
  if("stage_core_complete" %in% names(d))for(nm in names(domains)) {
    dd<-d[d$eid %in% domains[[nm]]$data$eid & d$stage_core_complete %in% TRUE,]
    if(nrow(dd)<100L)next
    sc<-run_cfg;sc$association$lag_years<-0;sc$association$yang<-FALSE;sc$association$yang_death<-FALSE
    rr<-ckm_run_associations(dd,sc,domains[[nm]])$assoc
    rr[,sensitivity:="complete_core_staging"];sens[[nm]]<-rr
  }
  ckm_write(data.table::rbindlist(sens,fill=TRUE),file.path(cfg$out_dir,"assoc.complete.tsv"))
}
if("marker" %in% modules) {
  pieces<-list()
  for(nm in intersect(names(domains),c("prot","met"))) {
    message("Cross-sectional stage markers: ",nm)
    pieces[[nm]]<-ckm_run_markers(d[d$eid %in% domains[[nm]]$data$eid,],cfg,domains[[nm]])
  }
  result<-list(marker=data.table::rbindlist(lapply(pieces,`[[`,"marker"),fill=TRUE),
    stage=data.table::rbindlist(lapply(pieces,`[[`,"stage"),fill=TRUE))
  save_cache(result,"marker")
  ckm_write(result$marker,file.path(cfg$out_dir,"marker.tsv"));ckm_write(result$stage,file.path(cfg$out_dir,"marker.stage.tsv"))
}
if("enrich" %in% modules) {
  a<-load_cache("assoc");result<-ckm_run_enrichment(if(is.null(a))list() else a$data,cfg)
  save_cache(result,"enrichment");ckm_write(result,file.path(cfg$out_dir,"enrichment.tsv"))
}
# Refresh report after every successful module, showing only matching cached phases.
clk<-load_cache("clock");asc<-load_cache("assoc");mark<-load_cache("marker");enr<-load_cache("enrichment")
clock_summary<-if(is.null(clk))list(notes="Clock has not been run for this phenotype/configuration.") else clk$data[setdiff(names(clk$data),c("predictions","folds"))]
if(!is.null(mark))clock_summary$markers<-mark$data$marker
assoc_summary<-if(is.null(asc))list() else asc$data[setdiff(names(asc$data),"settings")]
if(!is.null(asc)) {
  assoc_summary$global<-asc$data$assoc
  assoc_summary$pooled<-asc$data$pooled
}
payload<-list(title=if(cfg$demo)"CKM · SYNTHETIC DEMO / 合成数据演示" else "CKM · 连续分期与阶段特异因素",
  version=cfg$version,generated=format(Sys.time(),"%Y-%m-%d %H:%M UTC",tz="UTC"),
  phase=if(is.null(clk)&&is.null(asc))"phenotype_only" else "analysis",
  input_notes=c(if(cfg$demo)"SYNTHETIC TEST DATA. All numbers and gene sets are simulated, not UK Biobank findings." else "UK Biobank observational research outputs.",
    paste("Phenotype:",ph_meta$definition,"/",ph_meta$stage3_method),
    paste("Follow-up end:",ph_meta$date_follow_end,"from 0f/phenotype.R::date_follow_end; verified earlier coverage/loss/death may shorten individual follow-up."),
    "Stage comparison reference model includes age spline, sex and categorical stage.",
    "Clinical CVD includes CHD, heart failure, stroke, peripheral arterial disease and atrial fibrillation.",
    "Clock accuracy must be judged from held-out metrics and calibration; no accuracy gain is guaranteed.",
    ph_meta$score2_scope),
  phenotype=ckm_phenotype_summary(d),clock=clock_summary,associations=assoc_summary,
  enrichment=if(is.null(enr))data.frame(status="not_run") else enr$data,
  validation=list(phenotype_schema=ph_meta$schema,staging_definition=ph_meta$definition,
    date_follow_end=ph_meta$date_follow_end,followup_end_source="0f/phenotype.R::date_follow_end",shared_helpers=cfg$shared,
    clinical_assay_source=ph_meta$clinical_assay_source,statin_source=ph_meta$statin_source,
    R=R.version.string,packages=as.list(stats::setNames(vapply(c("data.table","survival","jsonlite"),function(p)as.character(utils::packageVersion(p)),character(1)),c("data.table","survival","jsonlite"))),
    source_reference_commit="971a61b44b84500ab46e2f5592478505a0eae50e",
    missing_incident_domains=ph_meta$missing_incident_domains,
    temporal_scope=ph_meta$temporal_scope,stage_uncertainty=ph_meta$stage_uncertainty,
    clock_status=if(is.null(clk))"not_run_or_changed_inputs" else "out_of_fold_results_available"))
report_json<-file.path(cfg$out_dir,"tmp","report.json")
jsonlite::write_json(payload,report_json,auto_unbox=TRUE,dataframe="rows",null="null",na="null",digits=10)
python<-Sys.getenv("CKM_PYTHON","python3")
exit<-system2(python,c(shQuote(file.path(fdir,"report.py")),"--input",shQuote(report_json),"--output",shQuote(file.path(cfg$out_dir,"ckm.html"))))
if(exit!=0L)stop("HTML report generation failed")
jsonlite::write_json(list(version=cfg$version,definition=ph_meta,settings=cfg),file.path(cfg$out_dir,"tmp","settings.json"),auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
message("Completed: ",module,"; ",file.path(cfg$out_dir,"ckm.html"))
