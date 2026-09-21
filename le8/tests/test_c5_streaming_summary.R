# Rscript tests/test_c5_streaming_summary.R
suppressPackageStartupMessages({library(dplyr);library(purrr);library(tibble);library(survival)})
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))[1]
fdir<-normalizePath(file.path(dirname(script),"..","f"))
for(f in c("c5_revision_validation.R","c5_joint_helpers.R","c5_joint_extensions.R",
           "c5_factorial.R","c5_joint_storage.R","c5_joint_outputs.R"))source(file.path(fdir,f))
# Exercise all numerical outputs; rendering/annotation are independent of storage.
c5_joint_plots<-c5_factorial_plots<-c5_final_panels<-c5_cell_output<-c5_evidence_output<-function(...)NULL
source<-function(file,...){
  if(!basename(file)%in%c("c5_factorial_plots.R","c5_panels_final.R"))base::source(file,...)
}
write_raw_csv<-function(x,file,rawdir){
  if(!ncol(x))x<-data.frame(note="No rows")
  write.csv(x,file.path(rawdir,file),row.names=FALSE)
}
expect_error<-function(expr,pattern){
  e<-tryCatch({force(expr);NULL},error=identity)
  stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e)))
}
run_tests<-function(){
  root<-tempfile("c5-stream-test-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  dir.create(file.path(root,"_c5_cache"));dir.create(file.path(root,"_c5_private"))
  set.seed(59);n<-300L;K<-3L
  yin<-data.frame(eid=as.character(seq_len(n)),.time=c(rep(2,30),rep(7,90),rep(12,180)),
    .event=c(rep(1,120),rep(0,180)),.fold=rep(seq_len(K),length.out=n),.entry=0,.risk=runif(n,.05,.4))
  design<-c5_factorial_design(FALSE);design<-design[design$available,]
  models<-c(design$model,"Clinical","Clinical_ProtMet_NS","Clinical_ProtMet_YS_YinYang")
  map<-data.frame(feature=c("prot__X","met__X"),assay=c("PX","MX"),layer=c("prot","met"))
  paths<-character();results<-list();arms<-c("all_assays","omit_GDF15_natriuretic")
  signature<-"synthetic-fixed-fits"
  availability<-bind_rows(lapply(c(0,5),function(L)data.frame(landmark=L,N=nrow(c5_landmark(yin,L,Inf)),available=TRUE)))
  for(L in c(0,5))for(fd in seq_len(K)){
    d<-c5_landmark(yin[yin$.fold==fd,],L,Inf)
    pred<-bind_rows(lapply(arms,function(a)bind_rows(lapply(models,function(mm){
      budget<-if(mm%in%c("F_C","Clinical"))0 else 5
      tibble(eid=d$eid,group=d$eid,fold=fd,time=d$.time,event=d$.event,model=mm,budget,
        arm=a,landmark=L,horizon=10-L,risk=d$.risk,lp=log(risk),
        high_training_q75=as.numeric(d$eid)>225,inflammatory_burden="unavailable")
    }))))
    members<-expand.grid(model=models,budget=5,arm=arms,fold=fd,landmark=L,feature=map$feature)
    r<-list(predictions=pred,members=as_tibble(members),diagnostics=tibble(model=models,status="ok"),
      marker_groups=tibble(eid=d$eid,fold=fd,landmark=L,GDF15_IL6_median="unavailable"))
    paths<-c(paths,file.path(root,"_c5_cache",paste0("L",L,".fold",fd,".rds")))
    c5_write_joint_checkpoint(list(signature=signature,result=r),tail(paths,1))
    results[[length(results)+1L]]<-r
  }
  roles<-yin[,c("eid",".fold")]
  write.csv(data.frame(Y="test",outer_folds=K,signature),file.path(root,"c5.cohort.csv"),row.names=FALSE)
  write.csv(availability,file.path(root,"c5.landmark_availability.csv"),row.names=FALSE)
  write.csv(map,file.path(root,"c5.assay_inventory.csv"),row.names=FALSE)
  write.csv(design,file.path(root,"c5.factorial_design.csv"),row.names=FALSE)
  write.csv(data.frame(status="synthetic"),file.path(root,"c5.prs_provenance.csv"),row.names=FALSE)
  writeLines("unavailable",file.path(root,"c5.omic_PGS_status.txt"))
  data.table::fwrite(roles,file.path(root,"_c5_private","roles.csv.gz"),compress="gzip")
  before<-tools::md5sum(paths)
  stage<-file.path(root,"stage")
  store<-c5_stage_joint_results(paths,signature,stage)
  all_pred<-bind_rows(lapply(results,`[[`,"predictions"))
  for(L in c(0,5))for(a in arms){
    expected<-all_pred[all_pred$landmark==L&all_pred$arm==a,]
    stopifnot(identical(c5_store_predictions(store,L,a),expected))
  }
  ids<-c5_validate_cached_cohorts(store$cohorts,availability,roles,K)
  for(L in c(0,5))stopifnot(identical(c5_output_ids(ids,L),c5_output_ids(yin,L)))
  expect_error(c5_validate_cached_cohorts(store$cohorts[-1,],availability,roles,K),"Incomplete")
  expect_error(c5_stage_joint_results(paths,"wrong-signature",file.path(root,"bad")),"signature mismatch")
  expect_error(c5_stage_joint_results(c(paths[1],paths[1]),signature,file.path(root,"duplicate")),"Duplicate")
  corrupt<-file.path(root,"corrupt.rds");writeLines("interrupted write",corrupt)
  expect_error(c5_stage_joint_results(corrupt,signature,file.path(root,"corrupt")),"Unreadable")
  expect_error(c5_joint_outputs(paths,yin,map,root,K,20,27,5,5,"wrong-signature"),"signature mismatch")
  stopifnot(!length(list.files(file.path(root,"_c5_private"),pattern="^[.]summary-",all.files=TRUE)))
  saved_plot<-c5_joint_plots
  c5_joint_plots<<-function(...)stop("simulated rendering failure")
  expect_error(c5_joint_outputs(paths,yin,map,root,K,0,27,5,5,signature),"simulated rendering failure")
  c5_joint_plots<<-saved_plot
  stopifnot(identical(tools::md5sum(paths),before))
  fresh<-c5_joint_outputs(paths,yin,map,root,K,20,27,5,5,signature)
  stopifnot(!any(file.exists(paths)),nrow(read.csv(file.path(root,"c5.fold_cache_cleanup.csv")))==length(paths))
  # Restore synthetic inputs to exercise recovery independently of the fresh run.
  for(i in seq_along(paths))c5_write_joint_checkpoint(list(signature=signature,result=results[[i]]),paths[i])
  stopifnot(identical(tools::md5sum(paths),before))
  unlink(file.path(root,"c5.res.rds"))
  frozen<-c5_resume_joint_outputs(root,"test",20,27,5,5)
  stopifnot(isTRUE(all.equal(fresh,frozen,tolerance=1e-12)),!any(file.exists(paths)))
  repeated<-c5_resume_joint_outputs(root,"test",20,27,5,5)
  stopifnot(isTRUE(all.equal(repeated,frozen,tolerance=1e-12)))
  ordinary<-c5_completed_joint_outputs(root,"test",20,27,5,5)
  stopifnot(isTRUE(all.equal(ordinary,frozen,tolerance=1e-12)),
    is.null(c5_completed_joint_outputs(root,"test",20,27,5,5,replace=TRUE)),
    is.null(c5_completed_joint_outputs(root,"other-trait",20,27,5,5)),
    is.null(c5_completed_joint_outputs(root,"test",21,27,5,5)))
  final_path<-file.path(root,"c5.res.rds");saved<-readRDS(final_path)
  failed<-saved;failed$summary_complete<-FALSE;saveRDS(failed,final_path)
  stopifnot(is.null(c5_completed_joint_outputs(root,"test",20,27,5,5)))
  failed<-saved;failed$signature<-"unrelated-fit";saveRDS(failed,final_path)
  stopifnot(is.null(c5_completed_joint_outputs(root,"test",20,27,5,5)))
  saveRDS(saved,final_path)
  prediction_path<-file.path(root,"_c5_private","out_of_fold_predictions.csv.gz")
  stopifnot(file.rename(prediction_path,paste0(prediction_path,".test")))
  stopifnot(is.null(c5_completed_joint_outputs(root,"test",20,27,5,5)))
  stopifnot(file.rename(paste0(prediction_path,".test"),prediction_path))
  stopifnot(all(frozen$paired_contrasts$delta_AUC==0),all(frozen$paired_contrasts$delta_Brier==0),
    all(frozen$factorial_contrasts$estimate==0))
  # Independent direct calculation verifies that blocking preserves absolute metrics.
  d<-all_pred[all_pred$model=="Clinical"&all_pred$landmark==0&all_pred$arm=="all_assays",]
  direct<-le8_evaluate_risk(d$time,d$event,d$risk,10,"Clinical",0,"nested outer validation",B=0)$metrics
  actual<-fresh$metrics|>filter(model=="Clinical",landmark==0,arm=="all_assays")
  stopifnot(isTRUE(all.equal(actual$AUC,direct$AUC)),isTRUE(all.equal(actual$Brier,direct$Brier)))
  export<-as.data.frame(data.table::fread(file.path(root,"_c5_private","out_of_fold_predictions.csv.gz")))
  stopifnot(nrow(export)==nrow(all_pred),sum(export$event)==sum(all_pred$event))
  stopifnot(!length(list.files(file.path(root,"_c5_private"),pattern="^[.]summary-",all.files=TRUE)))
  # Optionally compare every table with a saved pre-refactor implementation.
  previous<-Sys.getenv("LE8_TEST_PREVIOUS_C5_OUTPUTS","")
  if(nzchar(previous)){
    env<-new.env(parent=.GlobalEnv);sys.source(previous,envir=env)
    for(nm in c("c5_joint_plots","c5_cell_output","c5_evidence_output"))env[[nm]]<-function(...)NULL
    old<-env$c5_joint_outputs(results,yin,map,root,K,20,27,5,5,signature)
    canonical<-function(x){
      x<-as.data.frame(x);if(nrow(x)&&ncol(x))x<-x[do.call(order,lapply(x,as.character)),,drop=FALSE]
      rownames(x)<-NULL;x
    }
    for(nm in names(old))stopifnot(isTRUE(all.equal(canonical(old[[nm]]),canonical(fresh[[nm]]),tolerance=1e-12)))
  }
  unlink(paths[1]);expect_error(c5_resume_joint_outputs(root,"test",0,27,5,5),"every completed fold")
  cat("C5 streaming, frozen-cache recovery, integrity and numerical checks passed.\n")
}
run_tests()
