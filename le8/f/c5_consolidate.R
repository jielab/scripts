# FINAL C5 orchestrator: retain per-layer analyses, add joint validation,
# and always consolidate available C1–C5 evidence at analysis/le8/[Y]/.
fdir<-Sys.getenv("LE8_FDIR",file.path(Sys.getenv("DIRSCRIPT"),"f"))
source(file.path(fdir,"comm.f.R"))
LE8_JOB<-"c5_consolidate"
stage_log<-list()
run_stage<-function(name,file,env=character()){
  started<-le8_stage_start(paste0("C5/",name))
  logfile<-file.path(out.base,paste0("c5.",name,".log"))
  status<-tryCatch(system2(file.path(R.home("bin"),"Rscript"),shQuote(file.path(fdir,file)),
    stdout=logfile,stderr=logfile,env=env),error=function(e){writeLines(conditionMessage(e),logfile);1L})
  if(status==0)le8_stage_done(paste0("C5/",name),started)else message("[LE8] FAIL C5/",name," exit=",status," log=",logfile)
  stage_log[[length(stage_log)+1L]]<<-data.frame(stage=name,status=if(status==0)"complete"else"failed",
    detail=paste("Exit",status,";",basename(logfile)))
}
skip<-function(stage,detail)stage_log[[length(stage_log)+1L]]<<-data.frame(stage,status="unavailable",detail)
if(truthy(Sys.getenv("C5_RUN_REFERENCE","TRUE"))){
  for(layer in c(if(prot_DO)"prot",if(met_DO)"met")){
    d<-if(layer=="prot")out.prot else out.met
    if(file.exists(file.path(le8_job_dir(d,"c1_correlate"),"c1.res.rds")))
      run_stage(paste0("reference_",layer),"c5_legacy.R",env=c("LE8_C5_REFERENCE=TRUE",paste0("BIOM=",layer),
        paste0("PROT_DO=",if(layer=="prot")"TRUE"else"FALSE"),paste0("MET_DO=",if(layer=="met")"TRUE"else"FALSE")))
    else skip(paste0("reference_",layer),"C1 RDS unavailable; other layers/stages continue")
  }
}
if(truthy(Sys.getenv("C5_RUN_JOINT","TRUE"))){
  if(prot_DO&&met_DO)run_stage("joint","c5_joint.R")
  else skip("joint","Select --biom prot,met for a common-cohort joint comparison; per-layer modules remain available")
}
code<-file.path(fdir,"c5_systematic.py")
atlas_started<-le8_stage_start("C5/systematic_atlas")
status<-tryCatch(system2(Sys.getenv("PYTHON_BIN","python3"),vapply(c(code,"--root",out.base,"--cell-annotation"),shQuote,character(1))),error=function(e)1L)
if(status==0)le8_stage_done("C5/systematic_atlas",atlas_started)else message("[LE8] FAIL C5/systematic_atlas exit=",status)
stage_log[[length(stage_log)+1L]]<-data.frame(stage="systematic_atlas",status=if(status==0)"complete"else"failed",detail=paste("Exit",status))
log<-bind_rows(stage_log)
write.csv(log,file.path(out.base,"c5.pipeline_status.csv"),row.names=FALSE)
if(any(log$status=="failed"))stop("One or more C5 stages failed. Completed stages and upstream analyses are preserved; see c5.pipeline_status.csv")
