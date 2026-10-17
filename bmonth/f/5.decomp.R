source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
suppressPackageStartupMessages(library(PheWAS))
dat<-load_cohort()
pf<-file.path(indir,"Rdata","PheWAS.rds")
phenos<-validate_ids(readRDS(pf),"id","PheWAS.rds")
phenos[setdiff(names(phenos),"id")]<-lapply(
  phenos[setdiff(names(phenos),"id")],function(v)if(is.logical(v))v else as_num(v)==1)

pick_candidates<-function(){
  out<-list()
  sf<-file.path(rawdir,"phewas.rds")
  if(file.exists(sf)){
    s<-as_tibble(readRDS(sf))
    for(x in intersect(unique(s$exposure),PHEWAS_EXPOSURES)){
      d<-s%>%filter(exposure==x,is.finite(p))%>%arrange(p)
      thr<-unique(d$bonferroni_threshold);thr<-thr[is.finite(thr)][1]
      hit<-if(length(thr))d%>%filter(p<thr) else d[0,]
      # Keep all Bonferroni hits plus top 10 for transparent attenuation checks.
      out[[x]]<-unique(c(as.character(hit$phenotype),as.character(head(d$phenotype,10))))
    }
  }
  # Explicitly preserve alcohol phenotypes because OLD geography PheWAS showed
  # exceptionally strong associations and they are the key regression test.
  for(x in intersect(c("birth_north_100km","birth_east_100km"),PHEWAS_EXPOSURES))
    out[[x]] <- unique(c(out[[x]],"317","317.1","317.11"))

  jf<-file.path(rawdir,"phewas_joint.rds")
  if(file.exists(jf)){
    j<-as_tibble(readRDS(jf))
    for(x in intersect(c("month11df","cycle2df"),unique(j$exposure))){
      d<-j%>%filter(exposure==x,is.finite(p))%>%arrange(p)
      thr<-unique(d$bonferroni_threshold);thr<-thr[is.finite(thr)][1]
      hit<-if(length(thr))d%>%filter(p<thr) else d[0,]
      out[[x]]<-unique(c(as.character(hit$phenotype),as.character(head(d$phenotype,10))))
    }
  }
  out
}
cand<-pick_candidates()
if(!length(cand)) stop("No PheWAS results found. Run bmonth.sh phewas first.")

fit_scalar_stage<-function(exposure,phenotype){
  if(!phenotype%in%names(phenos)) return(tibble())
  ytab<-phenos[,c("id",phenotype),drop=FALSE];names(ytab)<-c("eid",".y")
  z<-inner_join(dat,ytab,by="eid")

  # IMPORTANT: fixed common D5 sample -> changes across stages are adjustment
  # effects, not changing-sample artifacts.
  common_cov<-setdiff(covars_stage$D5_migration,exposure)
  z<-complete_model_data(z,unique(c(".y",exposure,common_cov)))
  if(nrow(z)<500||sum(z$.y)<20||sum(!z$.y)<20||length(unique(z[[exposure]]))<2)
    return(tibble())

  rows<-lapply(names(covars_stage),function(stage){
    cv<-setdiff(covars_stage[[stage]],exposure)
    fit<-model_try(glm(reformulate(c(exposure,cv),response=".y"),data=z,family=binomial()))
    if(is.null(fit$fit))return(tibble(stage=stage,N=nrow(z),cases=sum(z$.y),
      beta=NA_real_,SE=NA_real_,OR=NA_real_,p=NA_real_))
    sm<-coef(summary(fit$fit))
    if(!exposure%in%rownames(sm))return(tibble(stage=stage,N=nrow(z),cases=sum(z$.y),
      beta=NA_real_,SE=NA_real_,OR=NA_real_,p=NA_real_))
    b<-sm[exposure,"Estimate"];se<-sm[exposure,"Std. Error"]
    tibble(stage=stage,N=nrow(z),cases=sum(z$.y),beta=b,SE=se,OR=exp(b),
           p=sm[exposure,"Pr(>|z|)"])
  })
  r<-bind_rows(rows);r$exposure<-exposure;r$phenotype<-phenotype
  b1<-r$beta[match("D1_core",r$stage)]
  b2<-r$beta[match("D2_adultSES",r$stage)]
  b3<-r$beta[match("D3_geography",r$stage)]
  r$attenuation_vs_D1<-if(is.finite(b1)&&abs(b1)>1e-8)1-r$beta/b1 else NA_real_
  r$attenuation_vs_D2<-if(is.finite(b2)&&abs(b2)>1e-8)1-r$beta/b2 else NA_real_
  r$attenuation_vs_D3<-if(is.finite(b3)&&abs(b3)>1e-8)1-r$beta/b3 else NA_real_
  r
}

scalar_exposures<-intersect(names(cand),ALL_EXPOSURES)
scalar<-bind_rows(lapply(scalar_exposures,function(x)
  bind_rows(lapply(cand[[x]],function(y)fit_scalar_stage(x,y)))))
if (!nrow(scalar)) scalar <- tibble(stage=character(),N=integer(),cases=integer(),
  beta=double(),SE=double(),OR=double(),p=double(),exposure=character(),phenotype=character(),
  attenuation_vs_D1=double(),attenuation_vs_D2=double(),attenuation_vs_D3=double())
scalar$analysis_track <- analysis_track(scalar$exposure)
scalar$exposure_unit <- exposure_unit(scalar$exposure)
write_csv2(scalar,"decomposition_scalar.csv")
saveRDS(scalar,file.path(rawdir,"decomposition_scalar.rds"))

# Focus table for the old geography-alcohol finding.
geo_alcohol <- scalar %>%
  filter(exposure %in% c("birth_north_100km","birth_east_100km"),
         phenotype %in% c("317","317.1","317.11")) %>%
  arrange(phenotype,exposure,match(stage,names(covars_stage)))
write_csv2(geo_alcohol,"audit_geography_alcohol.csv")

fit_cycle_stage<-function(phenotype){
  if(!phenotype%in%names(phenos)) return(tibble())
  ytab<-phenos[,c("id",phenotype),drop=FALSE];names(ytab)<-c("eid",".y")
  z<-inner_join(dat,ytab,by="eid")
  z<-complete_model_data(z,unique(c(".y","bmonth_cos","bmonth_sin",covars_stage$D5_migration)))
  if(nrow(z)<500||sum(z$.y)<20||sum(!z$.y)<20)return(tibble())
  bind_rows(lapply(names(covars_stage),function(stage){
    cv<-covars_stage[[stage]]
    a<-model_try(glm(reformulate(cv,response=".y"),z,family=binomial()))
    b<-model_try(glm(reformulate(c(cv,"bmonth_cos","bmonth_sin"),response=".y"),z,family=binomial()))
    if(is.null(a$fit)||is.null(b$fit))
      return(tibble(stage=stage,N=nrow(z),cases=sum(z$.y),chisq=NA_real_,df=NA_real_,p=NA_real_))
    lr<-anova(a$fit,b$fit,test="LRT")
    tibble(stage=stage,N=nrow(z),cases=sum(z$.y),
           chisq=lr$Deviance[2],df=lr$Df[2],p=lr$`Pr(>Chi)`[2])
  }))%>%mutate(phenotype=phenotype)
}
cycle_candidates<-unique(c(cand$cycle2df,cand$month11df))
cycle<-bind_rows(lapply(cycle_candidates,fit_cycle_stage))
if(nrow(cycle)){
  base<-cycle%>%filter(stage=="D1_core")%>%select(phenotype,chisq_D1=chisq)
  cycle<-left_join(cycle,base,by="phenotype")%>%
    mutate(chisq_attenuation_vs_D1=if_else(is.finite(chisq_D1)&chisq_D1>0,1-chisq/chisq_D1,NA_real_))
}
write_csv2(cycle,"decomposition_cycle.csv")
saveRDS(cycle,file.path(rawdir,"decomposition_cycle.rds"))
message("[decomp] saved five-stage common-sample decomposition and audit_geography_alcohol.csv")
