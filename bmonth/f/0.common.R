options(stringsAsFactors=FALSE,width=180)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(survival)
})

dir0 <- ifelse(Sys.info()[["sysname"]]=="Windows","D:","/mnt/d")
indir <- Sys.getenv("BMONTH_INDIR",file.path(dir0,"data","ukb","phe"))
outdir <- Sys.getenv("BMONTH_OUTDIR",file.path(dir0,"analysis","bmonth"))
rawdir <- file.path(outdir,"raw"); figdir <- file.path(outdir,"fig")
dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)

csv_env <- function(name,default=""){
  z <- trimws(strsplit(Sys.getenv(name,default),",",fixed=TRUE)[[1]])
  unique(z[nzchar(z)])
}
EXPOSURES <- csv_env("BMONTH_EXPOSURES",
  "bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn")
OUTCOMES <- csv_env("BMONTH_OUTCOMES","cvd_cad,ra")
TEMPMED_PHENOS <- csv_env("BMONTH_TEMPMED_PHENOS","")
REPLACE <- toupper(Sys.getenv("BMONTH_REPLACE","FALSE")) %in% c("TRUE","1","YES")

ALL_EXPOSURES <- c(
  "bmonth_cos","bmonth_sin","bmonth_winter_summer","bmonth_spring_autumn",
  "bmonth_rank","birth_north_100km","birth_east_100km","move_log1p")
if(length(setdiff(EXPOSURES,ALL_EXPOSURES)))
  stop("Invalid BMONTH_EXPOSURES: ",paste(setdiff(EXPOSURES,ALL_EXPOSURES),collapse=", "))

as_num <- function(x){
  if(is.numeric(x)||is.logical(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(as.character(x)))
}
as_date <- function(x){
  if(inherits(x,"Date")) return(x)
  if(inherits(x,"POSIXt")) return(as.Date(x))
  if(is.factor(x)) x<-as.character(x)
  if(is.character(x)) return(suppressWarnings(as.Date(x)))
  z<-suppressWarnings(as.numeric(unclass(x))); z[!is.finite(z)]<-NA_real_
  as.Date(z,origin="1970-01-01")
}
pmin_date <- function(...){
  z<-lapply(list(...),as_date)
  x<-do.call(pmin,c(z,list(na.rm=TRUE))); x<-as_date(x)
  x[!is.finite(as.numeric(x))]<-as.Date(NA); x
}
validate_ids <- function(d,id,label){
  if(!is.data.frame(d)||!id%in%names(d)) stop(label,": missing ",id)
  d<-as.data.frame(d); d[[id]]<-as.character(d[[id]])
  if(anyNA(d[[id]])||any(!nzchar(d[[id]]))||anyDuplicated(d[[id]]))
    stop(label,": invalid/duplicate IDs")
  if(anyDuplicated(names(d))) stop(label,": duplicate column names")
  d
}
read_all <- function(){
  p<-file.path(indir,"Rdata","all.rds")
  if(!file.exists(p)) stop("Missing ",p)
  validate_ids(readRDS(p),"eid","all.rds")
}
clean_month <- function(x){z<-as_num(x);z[!z%in%1:12]<-NA_real_;z}
clean_coord <- function(x){z<-as_num(x);z[!is.finite(z)|z<0]<-NA_real_;z}

birth_year_alias <- function(dat){
  cand <- c("birth_year","year_birth","year_of_birth","yob","p34")
  hit <- cand[cand%in%names(dat)]
  if(length(hit)) hit[1] else NA_character_
}
derive_birth_year <- function(dat){
  nm <- birth_year_alias(dat)
  if(!is.na(nm)){
    y<-as_num(dat[[nm]])
    y[!is.finite(y)|y<1900|y>2100|y!=floor(y)]<-NA_real_
    return(list(value=y,source=nm,approx=FALSE))
  }
  # Fallback only for climate linkage when no explicit year-of-birth column exists.
  # UKB age is integer age at assessment, so this can be off by one year.
  if(all(c("age","date_attend")%in%names(dat))){
    yr<-as_num(format(as_date(dat$date_attend),"%Y"))-as_num(dat$age)
    yr[!is.finite(yr)|yr<1900|yr>2100]<-NA_real_
    return(list(value=yr,source="assessment_year-age (approximate)",approx=TRUE))
  }
  list(value=rep(NA_real_,nrow(dat)),source="unavailable",approx=TRUE)
}

merge_move <- function(dat,required=TRUE){
  p<-file.path(indir,"Rdata","move.rds")
  if(!file.exists(p)){
    if(required) stop("Missing ",p,
      ". Run the ukb/f/phe.R move/tdi-at-birth block first.")
    return(dat)
  }
  mv<-validate_ids(readRDS(p),"eid","move.rds")
  keep<-intersect(c("eid","move","tdi.birth","tdi.birth.nn",
                    "lon.birth","lat.birth","lon.home","lat.home"),names(mv))
  if(!all(c("eid","move","tdi.birth")%in%keep))
    stop("move.rds must contain eid, move, and tdi.birth")
  mv<-mv[,keep,drop=FALSE]
  # all.rds may already contain these variables. Refresh from move.rds without
  # creating .x/.y names; retain existing values when the move table has no value.
  overlap<-intersect(setdiff(keep,"eid"),names(dat))
  out<-left_join(dat,mv,by="eid",suffix=c("",".move"))
  for(nm in overlap){
    out[[nm]]<-coalesce(as_num(out[[paste0(nm,".move")]]),as_num(out[[nm]]))
    out[[paste0(nm,".move")]]<-NULL
  }
  out
}

make_birth_month <- function(dat){
  req<-c("birth_month","birth_north","birth_east","move","tdi.birth")
  miss<-setdiff(req,names(dat)); if(length(miss))
    stop("Missing birth-month/location columns: ",paste(miss,collapse=", "))
  m<-clean_month(dat$birth_month); n<-clean_coord(dat$birth_north); e<-clean_coord(dat$birth_east)
  theta<-2*pi*m/12
  by<-derive_birth_year(dat)
  dat$birth_year_bm <- by$value
  if(!"tdi.birth.nn"%in%names(dat)) dat$tdi.birth.nn <- NA_real_
  out <- dat %>% mutate(
    birth_month=m,birth_north=n,birth_east=e,
    # Circular coding. -1 around Dec, +1 around Jun.
    bmonth_cos=-cos(theta),
    bmonth_sin=sin(theta),
    bmonth_winter_summer=case_when(m%in%c(12,1,2)~0,m%in%c(6,7,8)~1,TRUE~NA_real_),
    bmonth_spring_autumn=case_when(m%in%c(3,4,5)~0,m%in%c(9,10,11)~1,TRUE~NA_real_),
    bmonth_rank=m, # sensitivity only
    birth_north_100km=n/100000,
    birth_east_100km=e/100000,
    move_km=clean_coord(move),
    move_log1p=log1p(move_km),
    tdi.birth=as_num(tdi.birth),
    tdi.birth.nn=as_num(tdi.birth.nn),
    bmonth_cos_x_north=bmonth_cos*birth_north_100km,
    bmonth_sin_x_north=bmonth_sin*birth_north_100km,
    bmonth_cos_x_east=bmonth_cos*birth_east_100km,
    bmonth_sin_x_east=bmonth_sin*birth_east_100km,
    bmonth_cos_x_move=bmonth_cos*move_log1p,
    bmonth_sin_x_move=bmonth_sin*move_log1p
  )
  attr(out,"birth_year_source") <- by$source
  attr(out,"birth_year_approx") <- by$approx
  out
}

prep_covars <- function(dat){
  req<-c("age","sex","edu.sco","tdi","center","tdi.birth")
  miss<-setdiff(req,names(dat)); if(length(miss))
    stop("Missing covariates: ",paste(miss,collapse=", "))
  dat %>% mutate(
    age=as_num(age), sex=factor(sex), edu.sco=as_num(edu.sco),
    tdi=as_num(tdi), tdi.birth=as_num(tdi.birth), center=factor(center)
  )
}

# Four prespecified covariate stages.
# D4 is a post-birth migration sensitivity model, not a pure confounder model.
covars_stage <- list(
  D1_core=c("age","sex","center"),
  D2_adultSES=c("age","sex","center","edu.sco","tdi"),
  D3_birthContext=c("age","sex","center","edu.sco","tdi",
                    "birth_north_100km","birth_east_100km","tdi.birth"),
  D4_migration=c("age","sex","center","edu.sco","tdi",
                 "birth_north_100km","birth_east_100km","tdi.birth","move_log1p")
)
# Main birth-month analysis: isolate seasonality from birth geography/deprivation.
covars <- covars_stage$D3_birthContext

outcome_date_columns <- function(dat,y){
  p<-paste0("fod_icd10_",y); if(p%in%names(dat)) return(p)
  n<-names(dat); n[grepl("^(fod_|date_)",n)&endsWith(n,paste0("_",y))]
}
default_follow_end_path<-file.path(dir0,"scripts","0f","phenotype.R")
follow_end<-function(){
  v<-Sys.getenv("BMONTH_FOLLOW_END","")
  if(nzchar(v)) return(as.Date(v))
  if(!file.exists(default_follow_end_path)) stop("Missing follow-up source; use --follow-end")
  ln<-readLines(default_follow_end_path,n=1,warn=FALSE)
  m<-regmatches(ln,regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}",ln))
  if(!nzchar(m)) stop("Cannot parse date_follow_end"); as.Date(m)
}
ensure_outcome<-function(dat,y){
  tv<-paste0(y,".t2e"); ev<-paste0(y,".Yt2e")
  if(all(c(tv,ev)%in%names(dat))){
    dat[[tv]]<-as_num(dat[[tv]]); dat[[ev]]<-as_num(dat[[ev]])
    bad<-!is.finite(dat[[tv]])|dat[[tv]]<=0|!dat[[ev]]%in%c(0,1)
    dat[[tv]][bad]<-NA; dat[[ev]][bad]<-NA; return(dat)
  }
  if(!"date_attend"%in%names(dat)) stop("Cannot reconstruct ",y)
  dc<-outcome_date_columns(dat,y); if(!length(dc)) stop("No dated columns for ",y)
  event_date<-do.call(pmin_date,dat[dc]); base<-as_date(dat$date_attend)
  lost<-if("date_lost"%in%names(dat))as_date(dat$date_lost)else rep(as.Date(NA),nrow(dat))
  death<-if("date_death"%in%names(dat))as_date(dat$date_death)else rep(as.Date(NA),nrow(dat))
  censor<-pmin_date(rep(follow_end(),nrow(dat)),lost,death)
  prev<-!is.na(event_date)&!is.na(base)&event_date<=base
  event<-!prev&!is.na(event_date)&event_date>base&event_date<=censor
  exit<-as.Date(ifelse(event,as.numeric(event_date),as.numeric(censor)),origin="1970-01-01")
  dat[[tv]]<-as.numeric(exit-base)/365.25; dat[[ev]]<-as.integer(event)
  bad<-prev|!is.finite(dat[[tv]])|dat[[tv]]<=0
  dat[[tv]][bad]<-NA;dat[[ev]][bad]<-NA;dat
}
complete_model_data<-function(dat,vars){
  vars<-unique(vars)
  miss<-setdiff(vars,names(dat));if(length(miss))stop("Missing model columns: ",paste(miss,collapse=", "))
  d<-as.data.frame(dat[,vars,drop=FALSE]); ok<-complete.cases(d)
  for(v in names(d)[vapply(d,is.numeric,logical(1))]) ok<-ok&is.finite(d[[v]])
  droplevels(d[ok,,drop=FALSE])
}
model_try<-function(expr){
  notes<-character()
  value<-tryCatch(withCallingHandlers(expr,warning=function(w){
    notes<<-c(notes,conditionMessage(w));invokeRestart("muffleWarning")
  }),error=function(e){notes<<-c(notes,conditionMessage(e));NULL})
  if(any(grepl("converg|infinite|singular",notes,ignore.case=TRUE))) value<-NULL
  list(fit=value,note=paste(unique(notes),collapse="; "))
}
load_cohort<-function(){
  p<-file.path(rawdir,"bmonth.cohort.rds")
  if(!file.exists(p)) stop("Run bmonth.sh prep first")
  d<-validate_ids(readRDS(p),"eid","bmonth cohort")
  needed<-unique(c(ALL_EXPOSURES,unlist(covars_stage),"birth_year_bm"))
  if(length(setdiff(needed,names(d))) ||
     !identical(attr(d,"bmonth_metadata")$covars_primary,covars))
    stop("Cohort uses an older/incompatible model; rerun bmonth.sh prep first")
  d
}
fit_cox<-function(dat,exposure,y,covar_names=covars){
  covar_names<-setdiff(covar_names,exposure)
  dat<-ensure_outcome(dat,y);tv<-paste0(y,".t2e");ev<-paste0(y,".Yt2e")
  d<-complete_model_data(dat,unique(c(tv,ev,exposure,covar_names)))
  row<-tibble(outcome=y,exposure=exposure,N=nrow(d),events=sum(d[[ev]]==1),
    beta=NA_real_,se=NA_real_,HR=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
    p=NA_real_,status="skipped",note="")
  if(nrow(d)<500||row$events<20||length(unique(d[[exposure]]))<2) return(row)
  f<-reformulate(c(exposure,covar_names),response=paste0("Surv(`",tv,"`,`",ev,"`)"))
  fit<-model_try(coxph(f,data=d,ties="efron"))
  if(is.null(fit$fit)){row$status<-"failed";row$note<-fit$note;return(row)}
  sm<-coef(summary(fit$fit));if(!exposure%in%rownames(sm))return(row)
  b<-sm[exposure,"coef"];s<-sm[exposure,"se(coef)"]
  row$beta<-b;row$se<-s;row$HR<-exp(b);row$conf.low<-exp(b-qnorm(.975)*s)
  row$conf.high<-exp(b+qnorm(.975)*s);row$p<-sm[exposure,"Pr(>|z|)"];row$status<-"ok";row
}
annotate_phecodes <- function(x){
  # addPhecodeInfo performs an inner join: annotate a dictionary, then left join
  # so unmapped tests remain in the multiple-testing denominator and outputs.
  out<-as.data.frame(x);out$phenotype<-as.character(out$phenotype)
  info<-as.data.frame(PheWAS::addPhecodeInfo(data.frame(phenotype=unique(out$phenotype))))
  names(info)[names(info)=="group"]<-"category"
  out<-left_join(out,info,by="phenotype")
  if(!"description"%in%names(out))out$description<-out$phenotype
  missing_description<-is.na(out$description)|!nzchar(out$description)
  out$description[missing_description]<-out$phenotype[missing_description]
  if(!"category"%in%names(out)) out$category <- "other"
  out$category[is.na(out$category)|!nzchar(out$category)] <- "other"
  out
}
write_csv2<-function(x,name)write.csv(x,file.path(rawdir,name),row.names=FALSE,na="")
