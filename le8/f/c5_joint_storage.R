# Disk-backed aggregation: retain only one fold while staging and one
# landmark/ablation block while evaluating. No participant rows are dropped.
# Only known model-copy filenames for a durable, readable fold are eligible.
# Never remove fold checkpoints, roles, predictions, results, or newer files
# from a concurrently recomputing fold. Unknown files and symlinks stay intact.
c5_prune_fold_models<-function(checkpoint,signature=NULL,object=NULL,dry_run=FALSE){
  empty<-data.frame(file=character(),bytes=numeric(),action=character())
  if(!file.exists(checkpoint))return(empty)
  if(basename(dirname(checkpoint))!="_c5_cache"||nzchar(Sys.readlink(checkpoint)))return(empty)
  info<-file.info(checkpoint)
  if(is.na(info$mtime)||is.na(info$size)||info$isdir)return(empty)
  z<-if(is.null(object))tryCatch(readRDS(checkpoint),error=function(e)NULL)else object
  if(!is.list(z)||!is.character(z$signature)||length(z$signature)!=1L||is.na(z$signature)||
     !nzchar(z$signature)||(!is.null(signature)&&!identical(z$signature,signature)))return(empty)
  r<-z$result
  if(!is.list(r)||!is.data.frame(r$predictions)||!nrow(r$predictions)||
     !is.data.frame(r$coefficients)||!nrow(r$coefficients)||
     !is.data.frame(r$marker_groups)||!nrow(r$marker_groups))return(empty)
  mg<-r$marker_groups
  if(!all(c("eid","landmark","fold")%in%names(mg)))return(empty)
  L<-unique(mg$landmark);fd<-unique(mg$fold)
  if(length(L)!=1L||length(fd)!=1L||anyNA(c(L,fd))||anyNA(mg$eid)||anyDuplicated(mg$eid)||
     !all(c("eid","landmark","fold")%in%names(r$predictions))||
     anyNA(r$predictions[c("eid","landmark","fold")])||
     !all(r$predictions$eid%in%mg$eid)||!all(r$predictions$fold==fd)||!all(r$predictions$landmark==L)||
     basename(checkpoint)!=paste0("L",L,".fold",fd,".rds"))return(empty)
  prefix<-paste0("L",L,".fold",fd,".");names<-character()
  d<-r$diagnostics
  if(is.data.frame(d)&&all(c("model","arm","budget","status")%in%names(d))){
    d<-d[d$status%in%"ok"&is.finite(d$budget)&!is.na(d$model)&!is.na(d$arm),,drop=FALSE]
    names<-paste0(prefix,d$arm,".",d$model,".k",d$budget,".rds")
  }
  d<-r$proxy_validation
  if(is.data.frame(d)&&nrow(d)&&all(c("arm","cohort","modality","target")%in%names(d)))
    names<-c(names,paste0(prefix,d$arm,".proxy.",d$cohort,".",d$modality,".",d$target,".rds"))
  names<-unique(names[!grepl("[/\\\\]",names)])
  root<-dirname(dirname(checkpoint));private<-file.path(root,"_c5_private")
  if(nzchar(Sys.readlink(private)))return(empty)
  files<-file.path(private,names);fi<-file.info(files)
  keep<-!is.na(fi$size)&!is.na(fi$mtime)&!fi$isdir&fi$mtime<=info$mtime&!nzchar(Sys.readlink(files))
  files<-files[keep];fi<-fi[keep,,drop=FALSE]
  if(!length(files))return(empty)
  ans<-data.frame(file=files,bytes=fi$size,action=if(dry_run)"eligible"else"removed")
  if(dry_run)return(ans)
  # Recheck at deletion time so recently replaced/open-for-writing copies stay.
  deleted<-vapply(seq_along(files),function(i){
    now<-file.info(files[i])
    durable<-file.info(checkpoint)
    unchanged<-isTRUE(now$size==fi$size[i]&&now$mtime==fi$mtime[i]&&
      durable$size==info$size&&durable$mtime==info$mtime&&!nzchar(Sys.readlink(files[i])))
    unchanged&&unlink(files[i])==0L&&!file.exists(files[i])
  },logical(1))
  ans$action[!deleted]<-"kept: changed or removal failed"
  audit<-transform(ans,checkpoint=basename(checkpoint),time=format(Sys.time(),"%F %T %z"))
  log<-file.path(root,"c5.temporary_cleanup.csv")
  write.table(audit,log,sep=",",row.names=FALSE,col.names=!file.exists(log),append=file.exists(log))
  message("C5 cleanup: removed ",sum(deleted)," redundant model copies; ",
    sprintf("%.2f",sum(ans$bytes[deleted])/1024^3)," GiB; retained ",basename(checkpoint))
  ans
}

c5_process_token<-function(pid=Sys.getpid()){
  path<-file.path("/proc",as.character(pid),"stat")
  if(!file.exists(path))return(NA_character_)
  tryCatch({
    fields<-strsplit(sub("^.*[)] ","",suppressWarnings(readLines(path,warn=FALSE))[1])," +")[[1]]
    if(length(fields)<20L)NA_character_ else fields[20]
  },error=function(e)NA_character_)
}
c5_boot_id<-function(){
  tryCatch(suppressWarnings(readLines("/proc/sys/kernel/random/boot_id",warn=FALSE)[1]),error=function(e)NA_character_)
}
c5_summary_directory<-function(private){
  dir.create(private,recursive=TRUE,showWarnings=FALSE)
  d<-tempfile(paste0(".summary-",Sys.getpid(),"-"),tmpdir=private)
  dir.create(d,mode="0700")
  saveRDS(list(pid=Sys.getpid(),host=Sys.info()[["nodename"]],token=c5_process_token(),boot=c5_boot_id()),file.path(d,".owner.rds"))
  d
}
c5_cleanup_stale_summaries<-function(private){
  if(!dir.exists("/proc")||!dir.exists(private))return(invisible(NULL))
  dirs<-list.files(private,pattern="^[.]summary-",full.names=TRUE,all.files=TRUE)
  for(d in dirs){
    if(!dir.exists(d)||nzchar(Sys.readlink(d)))next
    o<-tryCatch(suppressWarnings(readRDS(file.path(d,".owner.rds"))),error=function(e)NULL)
    if(!is.list(o)||!identical(o$host,Sys.info()[["nodename"]])||length(o$pid)!=1L||
       !is.numeric(o$pid)||!is.finite(o$pid)||o$pid<1||length(o$token)!=1L||!is.character(o$token)||is.na(o$token))next
    current<-c5_process_token(o$pid)
    # An unreadable live /proc entry is ambiguous; keep it. A changed start
    # token identifies PID reuse and must not preserve a dead run's scratch.
    dead<-!dir.exists(file.path("/proc",as.character(o$pid)))
    boot<-c5_boot_id();reboot<-is.character(o$boot)&&length(o$boot)==1L&&!is.na(o$boot)&&!is.na(boot)&&!identical(boot,o$boot)
    if(dead||reboot||(!is.na(current)&&!identical(current,o$token))){
      unlink(d,recursive=TRUE);message("C5 cleanup: removed abandoned summary scratch ",basename(d))
    }
  }
  invisible(NULL)
}

c5_write_joint_checkpoint<-function(object,path){
  tmp<-tempfile(".fold-",tmpdir=dirname(path))
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(object,tmp)
  if(!file.rename(tmp,path))stop("Cannot publish C5 fold checkpoint: ",path)
}

# A completed summary, including plots/annotations and the published prediction
# export, replaces fold checkpoints. Failed summaries never reach this cleanup.
c5_retire_fold_checkpoints<-function(store,root,signature){
  empty<-data.frame(file=character(),bytes=numeric(),action=character())
  final<-tryCatch(readRDS(file.path(root,"c5.res.rds")),error=function(e)NULL)
  if(!isTRUE(final$summary_complete)||!identical(final$signature,signature)||
     !file.exists(file.path(root,"_c5_private","out_of_fold_predictions.csv.gz")))return(empty)
  inventory<-store$checkpoints
  if(!is.data.frame(inventory)||!nrow(inventory))return(empty)
  cache<-normalizePath(file.path(root,"_c5_cache"),mustWork=FALSE)
  rows<-lapply(seq_len(nrow(inventory)),function(i){
    p<-inventory$file[i];now<-file.info(p)
    safe<-identical(normalizePath(dirname(p),mustWork=FALSE),cache)&&
      grepl("^L[0-9.]+[.]fold[0-9]+[.]rds$",basename(p))&&!nzchar(Sys.readlink(p))&&
      !nzchar(Sys.readlink(dirname(p)))&&
      isTRUE(now$size==inventory$size[i]&&now$mtime==inventory$mtime[i])
    removed<-safe&&unlink(p)==0L&&!file.exists(p)
    data.frame(file=p,bytes=inventory$size[i],action=if(removed)"removed"else"kept: changed or unavailable")
  })
  ans<-do.call(rbind,rows)
  log<-file.path(root,"c5.fold_cache_cleanup.csv")
  write.table(transform(ans,time=format(Sys.time(),"%F %T %z")),log,sep=",",row.names=FALSE,
    col.names=!file.exists(log),append=file.exists(log))
  message("C5 cleanup: summary outputs published; removed ",sum(ans$action=="removed")," consumed fold checkpoints")
  ans
}

c5_stage_joint_results<-function(paths,signature,directory){
  dir.create(directory,recursive=TRUE,showWarnings=FALSE,mode="0700")
  checkpoint_info<-file.info(paths)
  checkpoint_info$file<-paths
  tables<-blocks<-cohorts<-list()
  seen<-character()
  for(i in seq_along(paths)){
    message("C5 summary: read fold ",i,"/",length(paths)," (",basename(paths[i]),")")
    z<-tryCatch(readRDS(paths[i]),error=function(e)stop("Unreadable C5 checkpoint: ",paths[i],": ",conditionMessage(e)))
    if(!identical(z$signature,signature))stop("C5 checkpoint signature mismatch: ",paths[i])
    r<-z$result;p<-r$predictions;mg<-r$marker_groups
    if(!is.data.frame(p)||!nrow(p)||!is.data.frame(mg)||!nrow(mg)||
       !all(c("eid","fold","landmark")%in%names(mg)))stop("Incomplete C5 checkpoint: ",paths[i])
    L<-unique(mg$landmark);fd<-unique(mg$fold)
    if(length(L)!=1L||length(fd)!=1L||anyNA(mg$eid)||anyDuplicated(mg$eid)||
       !all(p$fold==fd)||!all(p$landmark==L)||!all(p$eid%in%mg$eid))
      stop("Inconsistent C5 checkpoint participant/fold inventory: ",paths[i])
    key<-paste(L,fd,sep="/")
    if(key%in%seen)stop("Duplicate C5 landmark/fold: ",key)
    seen<-c(seen,key);cohorts[[i]]<-mg[,c("eid","fold","landmark"),drop=FALSE]
    for(nm in setdiff(names(r),"predictions")){
      file<-file.path(directory,paste0("fold-",i,"-",nm,".rds"))
      saveRDS(r[[nm]],file);tables[[nm]]<-c(tables[[nm]],file)
    }
    for(a in unique(p$arm)){
      file<-file.path(directory,paste0("prediction-",length(blocks)+1L,".rds"))
      saveRDS(p[p$arm==a,,drop=FALSE],file)
      blocks[[length(blocks)+1L]]<-data.frame(landmark=L,arm=a,fold=fd,file=file)
    }
    c5_prune_fold_models(paths[i],signature,object=z)
    rm(z,r,p,mg);invisible(gc())
  }
  if(!length(blocks))stop("No C5 predictions to summarize")
  list(tables=tables,blocks=bind_rows(blocks),cohorts=bind_rows(cohorts),
    checkpoints=checkpoint_info[,c("file","size","mtime"),drop=FALSE])
}

c5_store_table<-function(store,name){
  bind_rows(lapply(store$tables[[name]],readRDS))
}
c5_store_predictions<-function(store,landmark,arm){
  b<-store$blocks
  bind_rows(lapply(b$file[b$landmark==landmark&b$arm==arm],readRDS))
}
c5_release_predictions<-function(store,landmark,arm){
  b<-store$blocks
  unlink(b$file[b$landmark==landmark&b$arm==arm])
  invisible(NULL)
}
c5_output_ids<-function(yin,L){
  if(inherits(yin,"c5_cached_cohorts"))return(as.character(yin$eid[yin$landmark==L]))
  as.character(c5_landmark(yin,L,Inf)$eid)
}

c5_validate_cached_cohorts<-function(cohorts,availability,roles,K){
  if(!all(c("eid",".fold")%in%names(roles))||anyNA(roles$eid)||anyDuplicated(roles$eid))
    stop("Invalid saved C5 outer-fold roles")
  if(!setequal(unique(cohorts$landmark),availability$landmark))stop("Missing C5 landmark checkpoints")
  for(i in seq_len(nrow(availability))){
    a<-availability[i,];d<-cohorts[cohorts$landmark==a$landmark,,drop=FALSE]
    j<-match(d$eid,roles$eid)
    if(nrow(d)!=a$N||anyDuplicated(d$eid)||anyNA(j)||
       !setequal(unique(d$fold),seq_len(K))||any(d$fold!=roles$.fold[j]))
      stop("Incomplete or inconsistent saved C5 cohort at landmark ",a$landmark)
  }
  cohorts<-cohorts[order(match(cohorts$landmark,availability$landmark),match(cohorts$eid,roles$eid)),,drop=FALSE]
  class(cohorts)<-c("c5_cached_cohorts",class(cohorts));cohorts
}

# Completed numerical results outlive the disposable fold checkpoints. Normal
# runs and explicit recovery use the same check before loading raw omics.
c5_completed_joint_outputs<-function(root,trait,B,seed,primary_budget,primary_landmark,replace=FALSE){
  if(isTRUE(replace))return(NULL)
  required<-file.path(root,c("c5.cohort.csv","c5.res.rds","c5.metrics.csv",
    "c5.summary_provenance.csv","_c5_private/out_of_fold_predictions.csv.gz"))
  fi<-file.info(required)
  if(anyNA(fi$size)||any(fi$isdir)||any(fi$size<=0))return(NULL)
  cohort<-tryCatch(read.csv(required[1],stringsAsFactors=FALSE),error=function(e)NULL)
  if(!is.data.frame(cohort)||nrow(cohort)!=1L||!all(c("Y","signature")%in%names(cohort))||
     !identical(cohort$Y,trait)||is.na(cohort$signature)||!nzchar(cohort$signature))return(NULL)
  final<-tryCatch(readRDS(required[2]),error=function(e)NULL)
  settings<-list(bootstrap=B,seed=seed,primary_budget=primary_budget,primary_landmark=primary_landmark)
  if(!is.list(final)||!isTRUE(final$summary_complete)||!identical(final$signature,cohort$signature)||
     !isTRUE(all.equal(final$summary_settings,settings))||!is.list(final$tables)||
     !is.data.frame(final$tables$metrics)||!nrow(final$tables$metrics))return(NULL)
  message("C5: reuse completed joint outputs; no raw-data loading or model retraining. ",
    "Use --replace TRUE or a new output directory to recompute after changing inputs/model settings.")
  final$tables
}

# Explicit recovery uses the frozen cache's cohort, folds and signature. It
# never relabels old fits using newly loaded all.rds/omics or new training code.
c5_resume_joint_outputs<-function(root,trait,B,seed,primary_budget,primary_landmark){
  completed<-c5_completed_joint_outputs(root,trait,B,seed,primary_budget,primary_landmark)
  if(!is.null(completed))return(invisible(completed))
  required<-c("c5.cohort.csv","c5.landmark_availability.csv","c5.assay_inventory.csv",
    "c5.factorial_design.csv","c5.prs_provenance.csv","c5.omic_PGS_status.txt","_c5_private/roles.csv.gz")
  missing<-required[!file.exists(file.path(root,required))]
  if(length(missing))stop("C5 summary-only recovery needs saved metadata: ",paste(missing,collapse=", "))
  cohort<-read.csv(file.path(root,"c5.cohort.csv"),stringsAsFactors=FALSE)
  if(nrow(cohort)!=1L||cohort$Y!=trait||is.na(cohort$signature)||!nzchar(cohort$signature))stop("Invalid saved C5 cohort provenance")
  K<-as.integer(cohort$outer_folds)
  if(!is.finite(K)||K<3L)stop("Invalid saved C5 fold count")
  a<-read.csv(file.path(root,"c5.landmark_availability.csv"));a<-a[a$available%in%TRUE,,drop=FALSE]
  if(!nrow(a)||!primary_landmark%in%a$landmark)stop("Primary landmark unavailable in saved C5 run")
  paths<-unlist(lapply(a$landmark,function(L)file.path(root,"_c5_cache",paste0("L",L,".fold",seq_len(K),".rds"))),use.names=FALSE)
  if(any(!file.exists(paths)))stop("C5 summary-only requires every completed fold; missing: ",paste(basename(paths[!file.exists(paths)]),collapse=", "))
  roles<-as.data.frame(data.table::fread(file.path(root,"_c5_private","roles.csv.gz"),colClasses=c(eid="character")))
  map<-read.csv(file.path(root,"c5.assay_inventory.csv"),stringsAsFactors=FALSE)
  message("C5: summary-only recovery of ",length(paths)," frozen fold checkpoints; no model retraining")
  validation<-list(availability=a,roles=roles,K=K)
  c5_joint_outputs(paths,NULL,map,root,K,B,seed,primary_budget,primary_landmark,cohort$signature,
    cached_validation=validation)
}
