options(stringsAsFactors=FALSE,width=180)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(tibble);library(survival)})
dir0<-ifelse(Sys.info()[['sysname']]=='Windows','D:','/mnt/d')
indir<-Sys.getenv('BMONTH_INDIR',file.path(dir0,'data','ukb','phe'))
outdir<-Sys.getenv('BMONTH_OUTDIR',file.path(dir0,'analysis','bmonth'))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
rawdir<-file.path(outdir,'raw');dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
figdir<-file.path(outdir,'fig');dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
csv_env<-function(name,default=''){z<-trimws(strsplit(Sys.getenv(name,default),',',fixed=TRUE)[[1]]);unique(z[nzchar(z)])}
EXPOSURES<-csv_env('BMONTH_EXPOSURES','bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn,birth_north_100km,birth_east_100km')
OUTCOMES<-csv_env('BMONTH_OUTCOMES','cvd_cad,ra')
REPLACE<-toupper(Sys.getenv('BMONTH_REPLACE','FALSE'))%in%c('TRUE','1','YES')
ALL_EXPOSURES<-c('bmonth_cos','bmonth_sin','bmonth_winter_summer','bmonth_spring_autumn','birth_north_100km','birth_east_100km','bmonth_rank')
if(length(setdiff(EXPOSURES,ALL_EXPOSURES)))stop('Invalid BMONTH_EXPOSURES: ',paste(setdiff(EXPOSURES,ALL_EXPOSURES),collapse=', '))
validate_ids<-function(d,id,label){if(!is.data.frame(d)||!id%in%names(d))stop(label,': missing ',id);d<-as.data.frame(d);d[[id]]<-as.character(d[[id]]);if(anyNA(d[[id]])||any(!nzchar(d[[id]]))||anyDuplicated(d[[id]]))stop(label,': invalid IDs');d}
as_num<-function(x){if(is.numeric(x)||is.logical(x))return(as.numeric(x));suppressWarnings(as.numeric(as.character(x)))}
as_date<-function(x){if(inherits(x,'Date'))return(x);if(inherits(x,'POSIXt'))return(as.Date(x));if(is.factor(x))x<-as.character(x);if(is.character(x))return(suppressWarnings(as.Date(x)));z<-suppressWarnings(as.numeric(unclass(x)));z[!is.finite(z)]<-NA_real_;as.Date(z,origin='1970-01-01')}
pmin_date<-function(...){z<-lapply(list(...),as_date);x<-do.call(pmin,c(z,list(na.rm=TRUE)));x<-as_date(x);x[!is.finite(as.numeric(x))]<-as.Date(NA);x}
read_all<-function(){f<-file.path(indir,'Rdata','all.rds');if(!file.exists(f))stop('Missing ',f);validate_ids(readRDS(f),'eid','all.rds')}
clean_month<-function(x){z<-as_num(x);z[!z%in%1:12]<-NA_real_;z}
clean_coord<-function(x){z<-as_num(x);z[!is.finite(z)|z<0]<-NA_real_;z}
make_birth_month<-function(dat){
  req<-c('birth_month','birth_north','birth_east');miss<-setdiff(req,names(dat));if(length(miss))stop('Missing columns: ',paste(miss,collapse=', '))
  m<-clean_month(dat$birth_month);n<-clean_coord(dat$birth_north);e<-clean_coord(dat$birth_east)
  theta<-2*pi*m/12
  dat %>% mutate(
    birth_month=m,birth_north=n,birth_east=e,
    # Cyclic coding: bmonth_cos=-1 for Dec, +1 for Jun (same convention used in prior UKB season-of-birth work)
    bmonth_cos=-cos(theta),
    bmonth_sin=sin(theta),
    # Interpretable extreme-season contrasts; other months omitted from these two binary analyses
    bmonth_winter_summer=case_when(m%in%c(12,1,2)~0,m%in%c(6,7,8)~1,TRUE~NA_real_),
    bmonth_spring_autumn=case_when(m%in%c(3,4,5)~0,m%in%c(9,10,11)~1,TRUE~NA_real_),
    # Geographic coordinates scaled to 100 km. UKB field 129/130 are rounded north/east coordinates.
    birth_north_100km=n/100000,
    birth_east_100km=e/100000,
    # Secondary sensitivity only; months are circular, so Jan=1..Dec=12 is not a biologically preferred primary scale.
    bmonth_rank=m,
    # Joint seasonal interaction terms used in custom LRT PheWAS; do not interpret without their main effects.
    bmonth_cos_x_north=(-cos(theta))*(n/100000),
    bmonth_sin_x_north=sin(theta)*(n/100000),
    bmonth_cos_x_east=(-cos(theta))*(e/100000),
    bmonth_sin_x_east=sin(theta)*(e/100000)
  )
}
prep_covars<-function(dat){req<-c('age','sex','edu.sco','tdi','center');miss<-setdiff(req,names(dat));if(length(miss))stop('Missing covariates: ',paste(miss,collapse=', '));dat%>%mutate(age=as_num(age),sex=factor(sex),edu.sco=as_num(edu.sco),tdi=as_num(tdi),center=factor(center))}
covars<-c('age','sex','edu.sco','tdi','center')
outcome_date_columns<-function(dat,y){p<-paste0('fod_icd10_',y);if(p%in%names(dat))return(p);n<-names(dat);n[grepl('^(fod_|date_)',n)&endsWith(n,paste0('_',y))]}
default_follow_end_path<-file.path(dir0,'scripts','0f','phenotype.R')
follow_end<-function(){v<-Sys.getenv('BMONTH_FOLLOW_END','');if(nzchar(v))return(as.Date(v));if(!file.exists(default_follow_end_path))stop('Missing follow-up source; use --follow-end');ln<-readLines(default_follow_end_path,n=1,warn=FALSE);m<-regmatches(ln,regexpr('[0-9]{4}-[0-9]{2}-[0-9]{2}',ln));if(!nzchar(m))stop('Cannot parse date_follow_end');as.Date(m)}
ensure_outcome<-function(dat,y){tv<-paste0(y,'.t2e');ev<-paste0(y,'.Yt2e');if(all(c(tv,ev)%in%names(dat))){dat[[tv]]<-as_num(dat[[tv]]);dat[[ev]]<-as_num(dat[[ev]]);bad<-!is.finite(dat[[tv]])|dat[[tv]]<=0|!dat[[ev]]%in%c(0,1);dat[[tv]][bad]<-NA;dat[[ev]][bad]<-NA;return(dat)};if(!'date_attend'%in%names(dat))stop('Cannot reconstruct ',y);dc<-outcome_date_columns(dat,y);if(!length(dc))stop('No dated columns for ',y);event_date<-do.call(pmin_date,dat[dc]);base<-as_date(dat$date_attend);lost<-if('date_lost'%in%names(dat))as_date(dat$date_lost)else rep(as.Date(NA),nrow(dat));death<-if('date_death'%in%names(dat))as_date(dat$date_death)else rep(as.Date(NA),nrow(dat));censor<-pmin_date(rep(follow_end(),nrow(dat)),lost,death);prev<-!is.na(event_date)&!is.na(base)&event_date<=base;event<-!prev&!is.na(event_date)&event_date>base&event_date<=censor;exit<-as.Date(ifelse(event,as.numeric(event_date),as.numeric(censor)),origin='1970-01-01');dat[[tv]]<-as.numeric(exit-base)/365.25;dat[[ev]]<-as.integer(event);bad<-prev|!is.finite(dat[[tv]])|dat[[tv]]<=0;dat[[tv]][bad]<-NA;dat[[ev]][bad]<-NA;dat}
complete_model_data<-function(dat,vars){miss<-setdiff(vars,names(dat));if(length(miss))stop('Missing model columns: ',paste(miss,collapse=', '));d<-as.data.frame(dat[,vars,drop=FALSE]);ok<-complete.cases(d);for(v in names(d)[vapply(d,is.numeric,logical(1))])ok<-ok&is.finite(d[[v]]);droplevels(d[ok,,drop=FALSE])}
model_try<-function(expr){notes<-character();value<-tryCatch(withCallingHandlers(expr,warning=function(w){notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){notes<<-c(notes,conditionMessage(e));NULL});if(any(grepl('converg|infinite|singular',notes,ignore.case=TRUE)))value<-NULL;list(fit=value,note=paste(unique(notes),collapse='; '))}
load_cohort<-function(){f<-file.path(rawdir,'bmonth.cohort.rds');if(!file.exists(f))stop('Run bmonth.sh prep first');validate_ids(readRDS(f),'eid','bmonth cohort')}
fit_cox<-function(dat,exposure,y){dat<-ensure_outcome(dat,y);tv<-paste0(y,'.t2e');ev<-paste0(y,'.Yt2e');d<-complete_model_data(dat,unique(c(tv,ev,exposure,covars)));row<-tibble(outcome=y,exposure=exposure,N=nrow(d),events=sum(d[[ev]]==1),beta=NA_real_,se=NA_real_,HR=NA_real_,conf.low=NA_real_,conf.high=NA_real_,p=NA_real_,status='skipped',note='');if(nrow(d)<500||row$events<20||length(unique(d[[exposure]]))<2)return(row);f<-reformulate(c(exposure,covars),response=paste0('Surv(`',tv,'`,`',ev,'`)'));fit<-model_try(coxph(f,data=d,ties='efron'));if(is.null(fit$fit)){row$status<-'failed';row$note<-fit$note;return(row)};sm<-coef(summary(fit$fit));if(!exposure%in%rownames(sm))return(row);b<-sm[exposure,'coef'];s<-sm[exposure,'se(coef)'];row$beta<-b;row$se<-s;row$HR<-exp(b);row$conf.low<-exp(b-qnorm(.975)*s);row$conf.high<-exp(b+qnorm(.975)*s);row$p<-sm[exposure,'Pr(>|z|)'];row$status<-'ok';row}
write_csv2<-function(x,name)write.csv(x,file.path(rawdir,name),row.names=FALSE,na='')
# PheWAS::addPhecodeInfo uses an inner join. Annotation must not drop tests
# when the input dictionary contains a code absent from the installed package.
annotate_phecodes<-function(r){
  r<-as.data.frame(r);r$phenotype<-as.character(r$phenotype)
  info<-as.data.frame(PheWAS::addPhecodeInfo(data.frame(phenotype=unique(r$phenotype))))
  names(info)[names(info)=='group']<-'category'
  r<-left_join(r,info,by='phenotype')
  r$category[is.na(r$category)|!nzchar(r$category)]<-'other'
  missing_description<-is.na(r$description)|!nzchar(r$description)
  r$description[missing_description]<-r$phenotype[missing_description]
  r
}
