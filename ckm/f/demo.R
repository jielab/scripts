# Entirely synthetic integration fixture, never a scientific CKM result.
ckm_make_demo <- function(cfg,n=6000L) {
  n<-as.integer(n);if(n<2000L)stop("Demo requires at least 2000 participants")
  set.seed(cfg$seed);input<-file.path(cfg$out_dir,"test","input")
  dir.create(input,recursive=TRUE,showWarnings=FALSE)
  s<-sample(0:4,n,TRUE,c(.15,.2,.3,.25,.1));base<-as.Date("2010-01-01")+sample(0:365,n,TRUE)
  sex<-sample(0:1,n,TRUE);age<-runif(n,45,68)
  diet<-rnorm(n);shared<-rnorm(n);switch<-rnorm(n)
  alcohol<-sample(c("never","previous","current"),n,TRUE,c(.25,.15,.60))
  smoking<-sample(c("never","previous","current"),n,TRUE,c(.5,.3,.2))
  d<-data.table::data.table(eid=as.character(seq_len(n)+100000),date_attend=base,
    birth_date=base-round(age*365.25),age=age,sex=sex,center=paste0("demo_center_",sample(1:15,n,TRUE)),
    tdi=rnorm(n),PC1=rnorm(n),PC2=rnorm(n),ethnic.c="White",
    bmi=ifelse(s==0,runif(n,20,24),runif(n,26,33)),waist=ifelse(s==0,runif(n,65,77),runif(n,90,105)),
    sbp=ifelse(s<2,runif(n,105,118),runif(n,135,160)),dbp=ifelse(s<2,runif(n,65,76),runif(n,80,95)),
    hba1c=ifelse(s<2,runif(n,31,38),runif(n,38,52)),glucose=runif(n,4,5.4),
    tg=runif(n,.7,1.3),hdl=runif(n,1.5,2),total_chol=runif(n,4,5.6),
    egfr=ifelse(s==3,runif(n,32,42),runif(n,65,110)),uacr=ifelse(s==3,runif(n,35,120),runif(n,2,20)),
    creatinine=runif(n,60,90),fasting_hours=10,
    bp_tx=FALSE,statin=FALSE,dm_tx=FALSE,smoke_status=match(smoking,c("never","previous","current"))-1L,
    prevent_10y=ifelse(s==3,.25,.04),risk_index_date=base,risk_provenance="SYNTHETIC TEST VALUE; not PREVENT formula output",
    subclinical_ascvd=FALSE,pre_hf=FALSE,subclinical_ascvd_date=base,pre_hf_date=base,
    diet.pts=pmin(100,pmax(0,70+15*diet)),pa.pts=sample(seq(0,100,20),n,TRUE),
    smoke.c=factor(smoking,levels=c("never","previous","current")),sleep_duration=runif(n,5,9),
    alcohol.c=factor(alcohol,levels=c("never","previous","current")))
  rate<-.006*exp(.52*pmin(s,3)+.035*(age-55)+.35*shared+ifelse(s==1,.5,ifelse(s==3,-.4,.1))*switch-.2*diet)
  tcvd<-rexp(n,rate);td<-rexp(n,.006*exp(.05*(age-55)+.3*pmin(s,4)))
  lost<-ifelse(runif(n)<.15,runif(n,4,13),Inf);admin<-ckm_shared()$date_follow_end
  maxfu<-as.numeric(admin-base)/365.25
  d[,date_lost:=as.Date(ifelse(is.finite(lost)&lost<maxfu,as.numeric(base)+round(lost*365.25),NA_real_),origin="1970-01-01")]
  d[,date_death:=as.Date(ifelse(td<=pmin(lost,maxfu),as.numeric(base)+round(td*365.25),NA_real_),origin="1970-01-01")]
  cvd<-as.Date(ifelse(tcvd<=pmin(td,lost,maxfu),as.numeric(base)+pmax(1,round(tcvd*365.25)),NA_real_),origin="1970-01-01")
  cvd[s==4]<-base[s==4]-sample(1:2500,sum(s==4),TRUE)
  names_dx<-c("cvd_cad","cvd_hfail","cvd_stroke_i","cvd_stroke_ih","cvd_stroke_sh","cvd_pad","cvd_afib")
  ty<-sample(seq_along(names_dx),n,TRUE)
  for(j in seq_along(names_dx))d[[paste0("fod_icd10_",names_dx[j])]]<-as.Date(ifelse(ty==j,as.numeric(cvd),NA_real_),origin="1970-01-01")
  d[,`:=`(fod_icd10_t2dm=as.Date(NA),fod_icd10_cvd_htn=as.Date(NA))]
  d[sample(which(s<3),floor(n*.02)),uacr:=NA_real_]
  saveRDS(d,file.path(input,"all.rds"))
  prot<-data.table::data.table(eid=d$eid,SHARED_RISK=shared,STAGE_SWITCH=switch)
  for(j in 3:32)prot[[sprintf("PROT%03d",j)]]<-rnorm(n)+.12*s
  prot<-prot[sample(n,round(.7*n))]
  for(nm in setdiff(names(prot),"eid")) {j<-sample(nrow(prot),round(.03*nrow(prot)));data.table::set(prot,i=j,j=nm,value=NA_real_)}
  met<-data.table::data.table(eid=d$eid,MET_SHARED=exp(.4*shared+rnorm(n,.5,.4)),MET_SWITCH=exp(.4*switch+rnorm(n,.5,.4)))
  for(j in 3:18)met[[sprintf("MET%03d",j)]]<-exp(rnorm(n)+.08*s)
  saveRDS(prot,file.path(input,"prot.raw.rds"));saveRDS(met,file.path(input,"met.raw.rds"))
  ann<-data.table::rbindlist(list(data.table::data.table(source="prot",feature=setdiff(names(prot),"eid")),
    data.table::data.table(source="met",feature=setdiff(names(met),"eid"))))
  ann[,gene:=paste0("SYNTHETIC_GENE_",seq_len(.N))];ann[,lipid_block:=source=="met" & seq_len(.N)%%3L==0L]
  data.table::fwrite(ann,file.path(input,"annotation.tsv"),sep="\t")
  writeLines(c(paste(c("SYNTHETIC_SET_A","Demo only",ann$gene[1:15]),collapse="\t"),
    paste(c("SYNTHETIC_SET_B","Demo only",ann$gene[10:30]),collapse="\t"),
    paste(c("SYNTHETIC_SET_C","Demo only",tail(ann$gene,15)),collapse="\t")),file.path(input,"sets.gmt"))
  cfg$paths$ukb_dir<-input;cfg$paths$phe<-file.path(input,"all.rds")
  cfg$paths$biom<-cfg$paths$baseline<-cfg$paths$raw_vip<-cfg$paths$components<-cfg$paths$censor<-""
  cfg$paths$prot<-file.path(input,"prot.raw.rds");cfg$paths$met<-file.path(input,"met.raw.rds")
  cfg$paths$annotation<-file.path(input,"annotation.tsv");cfg$paths$gmt<-file.path(input,"sets.gmt")
  cfg$phenotype$allow_standardized_input<-TRUE
  cfg$clock$folds<-3L;cfg$clock$bootstrap<-20L;cfg$clock$max_omics<-4L
  cfg$clock$min_events<-30L;cfg$clock$min_at_horizon<-20L;cfg$clock$events_per_parameter<-5L
  cfg$association$min_n<-100L;cfg$association$min_events<-10L
  cfg$demo<-TRUE;cfg$version<-paste0(cfg$version," SYNTHETIC DEMO")
  jsonlite::write_json(cfg,file.path(input,"demo.config.json"),auto_unbox=TRUE,pretty=TRUE,null="null")
  cfg
}
