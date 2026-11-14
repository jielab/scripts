#!/usr/bin/env Rscript
# Reproducible scientific regression tests; uses only synthetic data.
args0<-commandArgs(FALSE)
this<-dirname(normalizePath(sub('^--file=','',grep('^--file=',args0,value=TRUE)[1])))
suppressPackageStartupMessages(library(data.table))
for(f in c('shared.R','common.R','risk.R','phenotype.R','clock.R','associations.R','markers.R'))source(file.path(this,f))
ckm_load_shared(file.path(this,'..','..','0f','phenotype.R'))
run_case<-function(label,f){cat('\nTEST: ',label,'\n',sep='');f();invisible(TRUE)}

run_case('Shared helper loading, configuration and random stream', function(){
  set.seed(31);before<-get('.Random.seed',envir=.GlobalEnv)
  helper<-file.path(this,'..','..','0f','phenotype.R')
  ckm_load_shared(helper)
  stopifnot(identical(before,get('.Random.seed',envir=.GlobalEnv)))
  stopifnot(!exists('t2e',envir=.GlobalEnv,inherits=FALSE),
    !exists('date_follow_end',envir=.GlobalEnv,inherits=FALSE))
  scratch<-tempfile('ckm-shared-');on.exit(unlink(scratch,recursive=TRUE))
  cfg<-ckm_config(file.path(this,'..','config','ckm.json'),list(out_dir=scratch))
  stopifnot(identical(cfg$association$covariates,ckm_shared()$vars.basic),
    identical(cfg$markers$covariates,ckm_shared()$vars.basic),
    identical(cfg$shared$date_follow_end,as.character(ckm_shared()$date_follow_end)),
    is.null(cfg$phenotype$admin_end))
  vals<-factor(c('2020-01-01','1900-01-01',NA))
  stopifnot(identical(ckm_ph_date(vals),ckm_shared()$.phenotype_date(vals,special_to_na=TRUE)))
})

run_case('Shared follow-up cutoff and shorter observation coverage', function(){
  shared<-ckm_shared();cutoff<-shared$date_follow_end
  on.exit(shared$date_follow_end<-cutoff)
  d<-data.table(eid=c('a','b','c'),baseline_date=cutoff-3650)
  cfg<-list(paths=list(),phenotype=list())
  default<-ckm_ph_censor(d,cfg)
  stopifnot(all(default$end==cutoff))
  shared$date_follow_end<-cutoff+40
  stopifnot(all(ckm_ph_censor(d,cfg)$end==cutoff+40))
  shared$date_follow_end<-cutoff
  file<-tempfile(fileext='.tsv');on.exit(unlink(file),add=TRUE)
  coverage<-data.table(eid=d$eid,censor_date=c(cutoff+100,cutoff-100,as.Date(NA)))
  fwrite(coverage,file,sep='\t');cfg$paths$censor<-file
  x<-ckm_ph_censor(d,cfg)
  stopifnot(identical(x$end,c(cutoff,cutoff-100,as.Date(NA))))
  coverage[,coverage_start:=c(cutoff-4000,cutoff-1000,cutoff-4000)]
  fwrite(coverage,file,sep='\t');x<-ckm_ph_censor(d,cfg)
  stopifnot(x$end[1]==cutoff,is.na(x$end[2]),is.na(x$end[3]))
  cfg$paths$censor<-'';cfg$phenotype$trust_input_censor<-TRUE
  d[,censor_date:=c(cutoff+100,cutoff-100,as.Date(NA))]
  stopifnot(identical(ckm_ph_censor(d,cfg)$end,c(cutoff,cutoff-100,as.Date(NA))))
})

run_case('Phenotype timing and stage rules', function(){
suppressPackageStartupMessages(library(data.table))
root <- tempfile('ckm-phenotype-'); dir.create(root)
on.exit(unlink(root, recursive=TRUE))
n <- 16L; base <- as.Date('2020-01-01')
d <- data.table(eid=as.character(1:n), date_attend=rep(base,n), birth_date=as.Date('1960-01-01'),
    age=60, sex=0, bmi=22, waist=70, sbp=110, dbp=70, hba1c=35,
    glucose=5, tg=1, hdl=1.6, total_chol=4, creatinine=70, egfr=90,
    uacr=10, bp_tx=FALSE,statin=FALSE,dm_tx=FALSE,smoke_status=0,
    fasting_hours=10, prevent_10y=.01,risk_index_date=base,risk_provenance='synthetic fixture',
    subclinical_ascvd=FALSE,pre_hf=FALSE,subclinical_ascvd_date=base,pre_hf_date=base,
    fod_icd10_cvd_cad=as.Date(NA),fod_icd10_cvd_hfail=as.Date(NA),
    fod_icd10_stroke=as.Date(NA),fod_icd10_cvd_pad=as.Date(NA),fod_icd10_cvd_afib=as.Date(NA),
    fod_icd10_t2dm=as.Date(NA),fod_icd10_cvd_htn=as.Date(NA), date_death=as.Date(NA))
d[2,bmi:=27]; d[3,sbp:=140]; d[4,prevent_10y:=.25]
d[5,fod_icd10_cvd_cad:=base-1];d[6,fod_icd10_cvd_cad:=base]
d[7:9,sbp:=140];d[7,fod_icd10_cvd_cad:=as.Date('2022-01-01')]
d[8,`:=`(fod_icd10_cvd_cad=as.Date('2023-01-01'),date_death=as.Date('2022-01-01'))]
d[9,`:=`(fod_icd10_cvd_cad=as.Date('2022-01-01'),date_death=as.Date('2022-01-01'))]
d[10,`:=`(egfr=40,uacr=NA_real_)];d[11,`:=`(egfr=25,uacr=NA_real_)]
d[12,`:=`(sbp=140,subclinical_ascvd=TRUE,subclinical_ascvd_date=base+100)]
d[13,fod_icd10_t2dm:=base+100]
d[,cnt_undated_icd10_cvd_cad:=0];d[14,cnt_undated_icd10_cvd_cad:=1]
d[15,fod_icd10_cvd_cad:=as.Date('1900-01-01')]
d[16,`:=`(date_death=as.Date('2022-01-01'),fod_icd10_cvd_cad=as.Date('2024-01-01'))]
saveRDS(d,file.path(root,'all.rds'))
cfg <- list(paths=list(phe=file.path(root,'all.rds'),ukb_dir=root,biom='',baseline='',raw_vip=''),
    phenotype=list(allow_standardized_input=TRUE,definition='aha2026',stage3_method='PREVENT'))
z <- ckm_build_phenotype(cfg)
stopifnot(identical(z$stage[1:7],c(0L,1L,2L,3L,4L,4L,2L)),
    z$status[7]==1L, is.na(z$status[8]), z$status[9]==1L,
    z$stage[10]==2L,z$stage_upper[10]==3L,!z$stage_core_complete[10],
    z$stage[11]==3L,z$stage_core_complete[11],z$stage[12]==2L,
    z$stage[13]==0L,!z$diabetes[13],is.na(z$status[14]),is.na(z$status[15]),
    all(is.na(z$status[5:6])),all(z$prevalent_cvd[5:6]),
    is.na(z$status[16]))
# Pure competing death remains in the risk set when there is no impossible CVD date.
d[8,fod_icd10_cvd_cad:=as.Date(NA)];saveRDS(d,file.path(root,'all.rds'))
z <- ckm_build_phenotype(cfg);stopifnot(z$status[8]==2L,z$event[8]==0L,z$event_death[8]==1L)
# Exact UACR units and missing UACR cannot become normal A1.
d[,uacr:=NA_real_];d[,`:=`(ualb=3.4,ucre=1000)];d[10,ualb:=NA_real_]
saveRDS(d,file.path(root,'all.rds'));z<-ckm_build_phenotype(cfg)
stopifnot(abs(z$uacr[1]-30.056)<1e-9,is.na(z$uacr[10]))
# A later component or risk-input date must fail before published staging.
d[1,risk_index_date:=base+1];saveRDS(d,file.path(root,'all.rds'))
stopifnot(inherits(try(ckm_build_phenotype(cfg),silent=TRUE),'try-error'))
# Raw clinical thresholds overwrite any imputed counterpart.
d <- d[1:2];d[,`:=`(risk_index_date=base,bmi=22,uacr=10,ualb=NA_real_,ucre=NA_real_,tg=9)]
saveRDS(d,file.path(root,'all.rds'))
raw <- data.table(eid=d$eid,bb_HBA1C=35,bb_GLU=5,bb_TG=1,bb_HDL=1.6,bb_TC=4,bb_CRE=70,bb_UALB=1,bb_UCRE=1000)
saveRDS(raw,file.path(root,'bbc.raw.rds'));cfg$paths$biom<-file.path(root,'bbc.raw.rds')
z <- ckm_build_phenotype(cfg);stopifnot(all(z$tg==1),all(z$stage==0L))
# Lost-to-follow-up cannot repair missing or post-baseline observation coverage.
d[,date_lost:=as.Date('2024-01-01')];saveRDS(d,file.path(root,'all.rds'))
coverage<-data.table(eid=d$eid,censor_date=as.Date('2025-01-01'),coverage_start=c(base,base+1))
fwrite(coverage,file.path(root,'coverage.tsv'),sep='\t');cfg$paths$censor<-file.path(root,'coverage.tsv')
z<-ckm_build_phenotype(cfg);stopifnot(z$incident_eligible[1],!z$incident_eligible[2],is.na(z$status[2]))
coverage[1,censor_date:=as.Date(NA)];fwrite(coverage,file.path(root,'coverage.tsv'),sep='\t')
z<-ckm_build_phenotype(cfg);stopifnot(all(is.na(z$status)))
# Original SCORE2 reference examples use no diabetic or CKD extrapolation.
s <- data.frame(age=c(50,50),sex=c('M','F'),sbp=140,total_chol=5.5,hdl=1.3,smoking=TRUE,diabetes=FALSE,egfr=90,uacr=10,prevalent_cvd=FALSE)
score <- ckm_score2_risk(s)$score2_10y
stopifnot(abs(score[1]-.059)<.002,abs(score[2]-.042)<.002)
s$diabetes[1]<-TRUE;stopifnot(is.na(ckm_score2_risk(s)$score2_10y[1]))
cat('Phenotype contracts: PASS (staging, same-day, death, future evidence, unknown timing, UACR, raw assay, SCORE2).\n')
})

run_case('Stage effects, shared equivalence and FDR', function(){
set.seed(20261113)
n <- 15000L
d <- data.frame(eid = as.character(seq_len(n)), stage = rep(1:3, each = n/3), age = runif(n, 45, 70),
                sex = factor(sample(c('F','M'),n,TRUE)), tdi = rnorm(n), PC1 = rnorm(n), PC2 = rnorm(n))
x <- rnorm(n); stage_x <- rnorm(n)
cat_x <- sample(c('Never','Previous','Current'),n,TRUE,prob=c(.5,.2,.3))
rate <- .06 * exp(.5*x + c(.8,0,-.8)[d$stage]*stage_x + .015*(d$age-55))
te <- rexp(n,rate); td <- rexp(n,.012); tc <- runif(n,7,13)
d$time <- pmin(te,td,tc); d$status <- ifelse(te<=td & te<=tc,1L,ifelse(td<te & td<=tc,2L,0L))
d$time_death <- pmin(td,tc);d$event_death <- as.integer(td<=tc)
extra <- d[1:900,];extra$eid <- paste0('Y',extra$eid);extra$stage <- 4L
d <- rbind(d,extra)
xx <- data.frame(eid=d$eid,shared=c(x,x[1:900]),specific=c(stage_x,stage_x[1:900]),
                 smoking=c(cat_x,cat_x[1:900]),all_missing=NA_real_)
m <- data.frame(feature=c('shared','specific','smoking','all_missing','not_present'),
                source='lifestyle',column=c('shared','specific','smoking','all_missing','not_present'),
                type=c('continuous','continuous','categorical','continuous','continuous'),
                ref=c('','','Never','',''))
cfg <- list(out_dir=NULL,association=list(covariates=c('age','sex'),stages=1:3,lag_years=c(0,2),
                   min_n=50,min_events=20,workers=1,ph_test=TRUE,equivalence_hr=1.2))
r <- ckm_run_associations(d,cfg,list(manifest=m,data=xx))
stopifnot(nrow(r$assoc)==5,nrow(r$sensitivity)==5,nrow(r$pooled)==6)
stopifnot(r$assoc[feature=='not_present',status]=='missing_exposure_column')
stopifnot(r$assoc[feature=='all_missing',status]=='constant_or_unavailable_exposure')
stopifnot(all(r$assoc$q_shared_planned_tests==5),all(r$assoc$q_interaction_planned_tests==5))
stopifnot(all.equal(r$assoc$q_shared, p.adjust(r$assoc$p_shared,'BH',n=5)))
stopifnot(r$assoc[feature=='specific',q_interaction]<1e-12)
stopifnot(r$stage[feature=='specific' & stage==1,beta]>0,r$stage[feature=='specific' & stage==3,beta]<0)
stopifnot(r$assoc[feature=='specific',pattern]=='stage_dependent_opposite_supported')
stopifnot(r$assoc[feature=='smoking',reference]=='Never',nrow(r$stage[feature=='smoking'])==6)
stopifnot(r$assoc[feature=='shared',events]==sum(d$status[d$stage %in%1:3]==1))
stopifnot(r$sensitivity[feature=='shared',n_before_complete]==sum(d$stage %in%1:3 & d$time>2))
stopifnot(all(r$yang[feature=='shared',status]=='ok'))
stopifnot(r$assoc[feature=='shared',ph_status]=='ok')
stopifnot(all(r$stage$q_stage_planned_tests==5))

# A rare third stage preserves estimable stage 1/2 effects and never invents
# an interaction P value from 'significant here, nonsignificant there'.
ds <- d[d$stage!=3 | seq_len(nrow(d)) %in% tail(which(d$stage==3),15),]
ms <- m[1,,drop=FALSE]
cs <- cfg;cs$association$lag_years <- 0;cs$association$yang <- FALSE;cs$association$yang_death <- FALSE
s <- ckm_run_associations(ds,cs,list(manifest=ms,data=xx))
stopifnot(is.na(s$assoc$p_interaction),s$assoc$interaction_status=='insufficient_stage_sample_or_events')
stopifnot(all(s$stage[stage%in%c(1,2),status]=='ok'),s$stage[stage==3,status]=='insufficient_stage_sample_or_events')
stopifnot(all(s$stage$estimation=='separate_stage_model'))

# Covariate errors are actionable configuration failures, not automatic
# changes to the scientific adjustment set.
cb <- cs;cb$association$covariates <- c('age','missing_covariate')
err <- try(ckm_run_associations(ds,cb,list(manifest=ms,data=xx)),silent=TRUE)
stopifnot(inherits(err,'try-error'),grepl('Configured adjustment columns are missing',as.character(err)))

# All endpoints and the FDR family are identical with forked execution.
cp <- cs;cp$association$workers <- 2
mp <- m[1:2,,drop=FALSE]
serial <- ckm_run_associations(d,cs,list(manifest=mp,data=xx))
parallel <- ckm_run_associations(d,cp,list(manifest=mp,data=xx))
stopifnot(identical(serial$assoc,parallel$assoc))

# A common adjusted effect can be equivalent within the prespecified margin
# even when a large sample detects a small stage difference.
d$known <- c(c(.8,0,-.8)[rep(1:3,each=n/3)]*stage_x,.8*stage_x[1:900])
ct <- cs;ct$association$covariates <- c('age','sex','known')
eq <- ckm_run_associations(d,ct,list(manifest=m[1,,drop=FALSE],data=xx))
stopifnot(eq$assoc$pattern=='shared_equivalent_within_margin',eq$assoc$q_equivalence<.05)
print(r$assoc[,c('feature','n','events','p_shared','q_interaction','p_equivalence','pattern','ph_status')])
cat('ASSOCIATION_SCIENTIFIC_TESTS_OK\n')
})

run_case('Markers and enrichment universe', function(){
set.seed(1113)
n <- 5000
d <- data.frame(eid=as.character(seq_len(n)),stage=rep(0:4,each=1000),age=rnorm(n,55,7),sex=factor(sample(c('F','M'),n,TRUE)),tdi=rnorm(n),PC1=rnorm(n),PC2=rnorm(n))
y <- .3*d$stage + .03*(d$age-55) + rnorm(n)
x <- data.frame(eid=d$eid,P1=y)
m <- data.frame(feature=c('P1','P_missing'),source='prot',column=c('P1','P_missing'))
cfg <- list(association=list(covariates=c('age','sex','tdi','PC1','PC2')),markers=list(min_stage_n=20))
r <- ckm_run_markers(d,cfg,list(manifest=m,data=x))
stopifnot(nrow(r$marker)==2,nrow(r$stage)==8,r$marker[feature=='P1',status]=='ok')
stopifnot(r$marker[feature=='P1',p]<1e-30,all(r$marker$q_planned_tests==2))
stopifnot(r$marker[feature=='P_missing',status]=='missing_exposure_column')
stopifnot(abs(r$stage[feature=='P1' & stage==4,beta]-1.2/sd(y))<.15)
nr <- ckm_run_markers(d[d$stage!=0,],cfg,list(manifest=m,data=x))
stopifnot(nr$marker[feature=='P1',status]=='reference_stage_unavailable_or_sparse')
sp <- d[d$stage!=3 | seq_len(n)%in%head(which(d$stage==3),5),]
sr <- ckm_run_markers(sp,cfg,list(manifest=m,data=x))
stopifnot(sr$marker[feature=='P1',status]=='ok',sr$stage[feature=='P1' & stage==3,status]=='sparse_stage')
stopifnot(all(sr$stage[feature=='P1' & stage%in%c(1,2,4),status]=='ok'))
stopifnot(ckm_run_enrichment(list(assoc=NULL,marker=NULL),cfg)$status=='not_configured')
ann <- data.frame(source='prot',feature=paste0('P',1:21),gene=paste0('G',1:21))
ann <- rbind(ann,ann[1,],data.frame(source='prot',feature='P2',gene='G1'))
annotation <- tempfile(fileext='.tsv');gmt <- tempfile(fileext='.gmt')
data.table::fwrite(ann,annotation,sep='\t')
writeLines(c(paste(c('A','Path A',paste0('G',1:6)),collapse='\t'),paste(c('B','Path B',paste0('G',10:15)),collapse='\t')),gmt)
cfg$paths <- list(annotation=annotation,gmt=gmt)
aa <- data.frame(source='prot',feature=paste0('P',1:21),status=c(rep('ok',20),'failed'),p_shared=c(rep(.01,20),NA),q_shared=c(rep(.02,4),rep(.9,16),NA),
                 p_interaction=c(rep(.8,20),NA),p_equivalence=c(rep(.001,20),NA),
                 pattern=c(rep('shared_equivalent_within_margin',4),rep('no_FDR_supported_association',17)),
                 direction=c(rep('protective_associated',4),rep('risk_associated',17)))
e <- ckm_run_enrichment(list(assoc=aa,pooled=NULL,marker=NULL),cfg)
t <- e[source=='prot' & set=='shared_protective' & status=='ok']
stopifnot(nrow(t)==2,all(t$background_genes==20),all(t$selected_genes==4),all(t$q_planned_tests==2))
stopifnot(t[pathway=='A',overlap]==4,t[pathway=='B',overlap]==0,t[pathway=='B',p]==1)
stopifnot(all.equal(t[pathway=='A',p],phyper(3,6,14,4,lower.tail=FALSE)))
unlink(c(annotation,gmt))
print(r$marker);print(t)
cat('MARKERS_AND_ENRICHMENT_TESTS_OK\n')
})


run_case('PREVENT exact threshold and reference example',function(){
 d<-data.frame(age=50,sex='F',sbp=160,bp_tx=TRUE,total_chol=200/38.67,hdl=45/38.67,statin=FALSE,diabetes=TRUE,smoking=FALSE,egfr=90,bmi=35,prevalent_cvd=FALSE,baseline_date=as.Date('2020-01-01'))
 stopifnot(abs(ckm_prevent_risk(d)$prevent_10y-.147)<.0005)
 d<-d[c(1,1),];d$age<-60
 for(j in 1:2)d$sbp[j]<-uniroot(function(v){q<-d[j,,drop=FALSE];q$sbp<-v;ckm_prevent_risk(q)$prevent_10y-c(.1998,.2002)[j]},c(90,180))$root
 stopifnot(identical(ckm_prevent_risk(d)$prevent_10y>=.2,c(FALSE,TRUE)))
})
run_case('Competing risk CIF and exact-horizon censoring',function(){
 grid<-seq(0,10,length.out=241)
 b1<-list(grid=grid,increment=rep(.03*10/240,240));b2<-list(grid=grid,increment=rep(.02*10/240,240))
 r<-ckm_cif(0,0,b1,b2)
 stopifnot(abs(r$risk-.03/.05*(1-exp(-.05*10)))<1e-10)
 stopifnot(abs(r$event_free_years-(1-exp(-.05*10))/.05)<1e-10)
 d<-data.frame(stage=rep(0:3,each=30),time=10,status=0L)
 w<-ckm_censor_weights(d,d,10)$weight
 stopifnot(all(w==1))
})
cat('\nALL SCIENTIFIC REGRESSION TESTS PASSED\n')
