# Rscript tests/test_c5_storage_cleanup.R
suppressPackageStartupMessages({library(dplyr);library(tibble);library(survival)})
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))[1]
fdir<-normalizePath(file.path(dirname(script),"..","f"))
source(file.path(fdir,"c5_revision_validation.R"))
source(file.path(fdir,"c5_joint_storage.R"))
run_tests<-function(){
  root<-tempfile("c5-cleanup-test-");dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  private<-file.path(root,"_c5_private");cache<-file.path(root,"_c5_cache")
  dir.create(private);dir.create(cache)
  checkpoint<-file.path(cache,"L0.fold1.rds")
  disposable<-file.path(private,c("L0.fold1.all_assays.Clinical.k0.rds",
    "L0.fold1.all_assays.proxy.Yin.prot.bmi.rds"))
  newer<-file.path(private,"L0.fold1.all_assays.NewFit.k5.rds")
  outside<-file.path(root,"outside.rds");writeLines("original",outside)
  link<-file.path(private,"L0.fold1.all_assays.Linked.k5.rds");file.symlink(outside,link)
  protected<-file.path(private,c("roles.csv.gz","out_of_fold_predictions.csv.gz","notes.rds",
    "L0.fold2.all_assays.Clinical.k0.rds","L0.fold1.all_assays.Unknown.k0.rds"))
  for(p in c(disposable,newer,protected))writeLines("fixture",p)
  r<-list(predictions=data.frame(eid=letters[1:4],fold=1L,landmark=0),
    marker_groups=data.frame(eid=letters[1:4],fold=1L,landmark=0),
    coefficients=data.frame(variable="age",beta=.1),
    diagnostics=data.frame(model=c("Clinical","NewFit","Linked"),arm="all_assays",budget=c(0,5,5),status="ok"),
    proxy_validation=data.frame(arm="all_assays",cohort="Yin",modality="prot",target="bmi"))
  z<-list(signature="finished-fold",result=r)
  c5_write_joint_checkpoint(z,checkpoint);digest<-tools::md5sum(checkpoint)
  for(p in disposable)Sys.setFileTime(p,Sys.time()-10)
  Sys.setFileTime(newer,Sys.time()+10)
  plan<-c5_prune_fold_models(checkpoint,"finished-fold",dry_run=TRUE)
  stopifnot(setequal(plan$file,disposable),all(file.exists(disposable)))
  stopifnot(!nrow(c5_prune_fold_models(checkpoint,"wrong-signature")))
  incomplete<-z;incomplete$result$predictions$fold[1]<-NA_integer_
  stopifnot(!nrow(c5_prune_fold_models(checkpoint,"finished-fold",object=incomplete)))
  stopifnot(!nrow(c5_prune_fold_models(file.path(cache,"L0.fold2.rds"))))
  bad<-file.path(cache,"L9.fold1.rds");writeLines("interrupted",bad)
  stopifnot(!nrow(c5_prune_fold_models(bad)))
  removed<-c5_prune_fold_models(checkpoint,"finished-fold")
  stopifnot(all(removed$action=="removed"),!any(file.exists(disposable)),
    all(file.exists(c(protected,newer,outside,link))),identical(tools::md5sum(checkpoint),digest),
    nrow(read.csv(file.path(root,"c5.temporary_cleanup.csv")))==2L)
  # Fold checkpoints are retired only after matching final outputs exist;
  # changed checkpoints from another run are never removed.
  inventory<-file.info(checkpoint);inventory$file<-checkpoint
  used<-list(checkpoints=inventory[,c("file","size","mtime")])
  result_path<-file.path(root,"c5.res.rds")
  saveRDS(list(signature="finished-fold",summary_complete=FALSE),result_path)
  stopifnot(!nrow(c5_retire_fold_checkpoints(used,root,"finished-fold")),file.exists(checkpoint))
  saveRDS(list(signature="other-run",summary_complete=TRUE),result_path)
  stopifnot(!nrow(c5_retire_fold_checkpoints(used,root,"finished-fold")),file.exists(checkpoint))
  saveRDS(list(signature="finished-fold",summary_complete=TRUE),result_path)
  Sys.setFileTime(checkpoint,Sys.time()+10)
  kept<-c5_retire_fold_checkpoints(used,root,"finished-fold")
  stopifnot(file.exists(checkpoint),all(kept$action!="removed"))
  inventory<-file.info(checkpoint);inventory$file<-checkpoint;used$checkpoints<-inventory[,c("file","size","mtime")]
  retired<-c5_retire_fold_checkpoints(used,root,"finished-fold")
  stopifnot(all(retired$action=="removed"),!file.exists(checkpoint),all(file.exists(protected)))
  # Live owners, unrelated directories and symlinks must survive scavenging.
  live<-c5_summary_directory(private)
  abandoned<-c5_summary_directory(private)
  o<-readRDS(file.path(abandoned,".owner.rds"));o$pid<-2147483647
  saveRDS(o,file.path(abandoned,".owner.rds"))
  reused<-c5_summary_directory(private);o$pid<-Sys.getpid();o$token<-"previous-process"
  saveRDS(o,file.path(reused,".owner.rds"))
  unknown<-file.path(private,".summary-unknown");dir.create(unknown)
  other_host<-c5_summary_directory(private);o$host<-"unrelated-host"
  saveRDS(o,file.path(other_host,".owner.rds"))
  linked_dir<-file.path(private,".summary-link");file.symlink(live,linked_dir)
  c5_cleanup_stale_summaries(private)
  stopifnot(!dir.exists(abandoned),!dir.exists(reused),
    all(dir.exists(c(live,unknown,other_host,linked_dir))))
  # Release only consumed prediction shards, keeping metadata and other blocks.
  shards<-file.path(root,paste0("shard",1:3,".rds"));for(p in shards)saveRDS(1,p)
  store<-list(blocks=data.frame(landmark=c(0,0,5),arm=c("a","b","a"),file=shards))
  c5_release_predictions(store,0,"a")
  stopifnot(!file.exists(shards[1]),all(file.exists(shards[-1])))
  # Unused omic columns must never be serialized through a Cox formula environment.
  set.seed(108);n<-900L
  d<-data.frame(age=rnorm(n),x=rnorm(n),time=rexp(n,1/8),event=rbinom(n,1,.4))
  tr<-d[1:700,];te<-d[701:n,]
  wide<-tr;wide$unused_omics<-I(matrix(rnorm(nrow(tr)*500L),nrow(tr)))
  a<-le8_fit_budget_model(wide,te,"age","x","time","event","cox")
  b<-le8_fit_budget_model(tr,te,"age","x","time","event","cox")
  stopifnot(a$status=="ok",is.null(a[["fit"]]),isTRUE(all.equal(a,b)))
  serialized<-serialize(a,NULL)
  stopifnot(length(serialized)<150000L)
  restored<-unserialize(serialized)
  for(h in c(0,2,5,10,30))stopifnot(identical(le8_risk_at(a,h),le8_risk_at(restored,h)))
  # The optional baseline-hazard audit keeps jumps only without changing risk.
  lean<-a;lean$baseline_hazard<-a$baseline_hazard[!duplicated(a$baseline_hazard$hazard),]
  for(h in unique(c(0,a$baseline_hazard$time,50)))stopifnot(identical(le8_risk_at(a,h),le8_risk_at(lean,h)))
  previous<-Sys.getenv("LE8_TEST_PREVIOUS_MODEL","")
  if(nzchar(previous)){
    env<-new.env(parent=.GlobalEnv);sys.source(previous,envir=env)
    old<-env$le8_fit_budget_model(wide,te,"age","x","time","event","cox")
    before<-length(serialize(old,NULL));old$fit<-NULL
    stopifnot(isTRUE(all.equal(old,a)),before>10*length(serialized))
    cat("Serialized Cox result:",before,"->",length(serialized),"bytes; numerical fields unchanged\n")
  }
  cat("C5 temporary cleanup, live-run protection and compact-model regression tests passed.\n")
}
run_tests()
