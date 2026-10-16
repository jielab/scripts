source(file.path(Sys.getenv("BMONTH_FDIR"),"0.common.R"))
dat<-load_cohort()
tf<-file.path(rawdir,"bmonth_temperature.rds")
if(!file.exists(tf))stop("Run `bmonth.sh temp --temp-dir ...` first.")
temp<-validate_ids(readRDS(tf),"eid","bmonth_temperature.rds")
if(!identical(attr(temp,"temperature_metadata")$cohort_md5,
              unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds")))))
  stop("Temperature file does not match this cohort; rerun bmonth.sh temp first")
tvars<-intersect(c("temp_conception_mean","temp_birth_mean",
                   "temp_conception_anom","temp_birth_anom"),names(temp))
tvars<-tvars[vapply(temp[tvars],function(x)sum(is.finite(as_num(x)))>=1000,logical(1))]
if(!length(tvars))stop("No usable temperature variables.")
dat<-left_join(dat,temp[,c("eid",tvars),drop=FALSE],by="eid")

pf<-file.path(indir,"Rdata","PheWAS.rds")
phenos<-validate_ids(readRDS(pf),"id","PheWAS.rds")
phenos[setdiff(names(phenos),"id")]<-lapply(
  phenos[setdiff(names(phenos),"id")],function(v)if(is.logical(v))v else as_num(v)==1)

pick_temp_candidates<-function(){
  if(length(TEMPMED_PHENOS))return(TEMPMED_PHENOS)
  out<-character()
  sf<-file.path(rawdir,"phewas.rds")
  if(file.exists(sf)){
    s<-as_tibble(readRDS(sf))%>%
      filter(exposure%in%c("bmonth_cos","bmonth_sin","bmonth_winter_summer",
                           "bmonth_spring_autumn"),is.finite(p))%>%arrange(p)
    if(nrow(s)){
      hit<-s%>%filter(p<bonferroni_threshold)
      out<-c(out,as.character(hit$phenotype),as.character(head(s$phenotype,10)))
    }
  }
  jf<-file.path(rawdir,"phewas_joint.rds")
  if(file.exists(jf)){
    j<-as_tibble(readRDS(jf))%>%
      filter(exposure%in%c("cycle2df","month11df"),is.finite(p))%>%arrange(p)
    if(nrow(j)){
      hit<-j%>%filter(p<bonferroni_threshold)
      out<-c(out,as.character(hit$phenotype),as.character(head(j$phenotype,10)))
    }
  }
  unique(out)
}
candidates<-pick_temp_candidates()
candidates<-intersect(candidates,setdiff(names(phenos),"id"))
if(!length(candidates))stop("No birth-month PheWAS phenotypes available for temperature decomposition.")

# 2-df climate attenuation: how much annual-cycle evidence remains after adding temperature.
cycle_one<-function(yname,m){
  ytab<-phenos[,c("id",yname),drop=FALSE];names(ytab)<-c("eid",".y")
  z<-inner_join(dat,ytab,by="eid")
  z<-complete_model_data(z,c(".y","bmonth_cos","bmonth_sin",m,covars))
  if(nrow(z)<500||sum(z$.y)<20||sum(!z$.y)<20)return(tibble())
  f0<-model_try(glm(reformulate(covars,response=".y"),z,family=binomial()))
  f1<-model_try(glm(reformulate(c(covars,"bmonth_cos","bmonth_sin"),response=".y"),z,family=binomial()))
  fm0<-model_try(glm(reformulate(c(covars,m),response=".y"),z,family=binomial()))
  fm1<-model_try(glm(reformulate(c(covars,m,"bmonth_cos","bmonth_sin"),response=".y"),z,family=binomial()))
  if(any(vapply(list(f0$fit,f1$fit,fm0$fit,fm1$fit),is.null,logical(1))))return(tibble())
  a<-anova(f0$fit,f1$fit,test="LRT");b<-anova(fm0$fit,fm1$fit,test="LRT")
  chi0<-a$Deviance[2];chi1<-b$Deviance[2]
  tibble(phenotype=yname,mediator=m,N=nrow(z),cases=sum(z$.y),
         cycle_chisq_before=chi0,cycle_p_before=a$`Pr(>Chi)`[2],
         cycle_chisq_after=chi1,cycle_p_after=b$`Pr(>Chi)`[2],
         chisq_attenuation=if(is.finite(chi0)&&chi0>0)1-chi1/chi0 else NA_real_)
}
cycle<-bind_rows(lapply(candidates,function(y)
  bind_rows(lapply(tvars,function(m)cycle_one(y,m)))))
write_csv2(cycle,"temperature_cycle_attenuation.csv")
saveRDS(cycle,file.path(rawdir,"temperature_cycle_attenuation.rds"))

# Scalar mediation-like coefficient product for sine/cosine components.
B<-as.integer(Sys.getenv("BMONTH_MED_BOOT","100"));set.seed(12345)
scalar_one<-function(yname,x,m){
  ytab<-phenos[,c("id",yname),drop=FALSE];names(ytab)<-c("eid",".y")
  z<-inner_join(dat,ytab,by="eid")
  z<-complete_model_data(z,c(".y",x,m,covars))
  if(nrow(z)<1000||sum(z$.y)<50||sum(!z$.y)<50)return(tibble())
  if(length(unique(z[[x]]))<2L||length(unique(z[[m]]))<2L)return(tibble())
  z$.x<-as.numeric(scale(z[[x]]));z$.m<-as.numeric(scale(z[[m]]))
  af<-reformulate(c(".x",covars),response=".m")
  tf<-reformulate(c(".x",covars),response=".y")
  df<-reformulate(c(".x",".m",covars),response=".y")
  fa<-model_try(lm(af,z));ft<-model_try(glm(tf,z,family=binomial()))
  fd<-model_try(glm(df,z,family=binomial()))
  if(is.null(fa$fit)||is.null(ft$fit)||is.null(fd$fit))return(tibble())
  sa<-coef(summary(fa$fit));st<-coef(summary(ft$fit));sd<-coef(summary(fd$fit))
  if(!".x"%in%rownames(sa)||!".x"%in%rownames(st)||
     !all(c(".x",".m")%in%rownames(sd)))return(tibble())
  a<-sa[".x","Estimate"];ase<-sa[".x","Std. Error"]
  b<-sd[".m","Estimate"];bse<-sd[".m","Std. Error"]
  total<-st[".x","Estimate"];direct<-sd[".x","Estimate"];indirect<-a*b
  sobelse<-sqrt(b^2*ase^2+a^2*bse^2)
  pb<-rep(NA_real_,B)
  for(i in seq_len(B)){
    q<-z[sample.int(nrow(z),replace=TRUE),,drop=FALSE]
    aa<-coef(model_try(lm(af,q))$fit)[".x"]
    bb<-coef(model_try(glm(df,q,family=binomial()))$fit)[".m"]
    if(length(aa)!=1L||length(bb)!=1L)next
    if(is.finite(aa)&&is.finite(bb))pb[i]<-aa*bb
  }
  bb<-pb[is.finite(pb)]
  ok<-B>=50&&length(bb)>=max(50,ceiling(.8*B))
  tibble(phenotype=yname,exposure=x,mediator=m,N=nrow(z),cases=sum(z$.y),
         a_beta=a,b_beta=b,total_beta=total,direct_beta=direct,
         indirect_beta=indirect,p_sobel=if(sobelse>0)2*pnorm(abs(indirect/sobelse),lower.tail=FALSE)else NA_real_,
         bootstrap_valid=length(bb),
         p_boot=if(ok)min(1,2*min((sum(bb<=0)+1)/(length(bb)+1),(sum(bb>=0)+1)/(length(bb)+1)))else NA_real_,
         indirect_lo=if(ok)quantile(bb,.025,names=FALSE)else NA_real_,
         indirect_hi=if(ok)quantile(bb,.975,names=FALSE)else NA_real_,
         attenuation_beta=if(is.finite(total)&&abs(total)>1e-8)1-direct/total else NA_real_)
}
scalar<-bind_rows(lapply(candidates,function(y)
  bind_rows(lapply(c("bmonth_cos","bmonth_sin"),function(x)
    bind_rows(lapply(tvars,function(m)scalar_one(y,x,m)))))))
if(nrow(scalar))scalar$FDR_indirect<-p.adjust(scalar$p_boot,"BH")
write_csv2(scalar,"temperature_mediation.csv")
saveRDS(scalar,file.path(rawdir,"temperature_mediation.rds"))
writeLines(c(
  "Temperature analysis is mechanistic attenuation / mediation-like decomposition, not proof of causal mediation.",
  "Temperature is reconstructed from month/year of birth and birthplace, so it is strongly structured by the exposure.",
  "The 2-df cycle analysis compares annual-cycle likelihood-ratio evidence before versus after temperature adjustment.",
  "Scalar sine/cosine analyses use a linear mediator model and logistic outcome model; a*b is on the outcome log-odds scale.",
  "Migration distance is not included in the primary covariates because it occurs after birth; it is examined separately in D4."
),file.path(rawdir,"temperature_methods.txt"))
message("[tempmed] saved temperature cycle attenuation and scalar mediation-like results.")
