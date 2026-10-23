options(stringsAsFactors=FALSE,width=180)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(survival)
})

dir0 <- ifelse(Sys.info()[["sysname"]]=="Windows","D:","/mnt/d")
indir <- Sys.getenv("BMONTH_INDIR",file.path(dir0,"data","ukb"))
outdir <- Sys.getenv("BMONTH_OUTDIR",file.path(dir0,"analysis","bmonth"))
rawdir <- file.path(outdir,"raw"); figdir <- file.path(outdir,"fig")
dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)

csv_env <- function(name,default=""){
  z <- trimws(strsplit(Sys.getenv(name,default),",",fixed=TRUE)[[1]])
  unique(z[nzchar(z)])
}

# Main scientific exposures: birth month / annual seasonality.
EXPOSURES <- csv_env("BMONTH_EXPOSURES",
  "bmonth_cos,bmonth_sin,bmonth_winter_summer,bmonth_spring_autumn")

# Context exposures: retained to reproduce and investigate the strong OLD
# birthplace geography -> alcohol PheWAS results. These are not birth-month effects.
CONTEXT_EXPOSURES <- csv_env("BMONTH_CONTEXT_EXPOSURES",
  "birth_north_100km,birth_east_100km")
EXPLORATORY_ALL <- c("border12", "border1n", "border_rank", "border_relative")
EXPLORATORY_TRAITS <- csv_env("BMONTH_EXPLORATORY_TRAITS", paste(EXPLORATORY_ALL, collapse = ","))
FAMILY_SIZE_SENS_TRAITS <- intersect(c("border1n","border_rank","border_relative"), EXPLORATORY_TRAITS)
if (identical(EXPLORATORY_TRAITS, "none")) EXPLORATORY_TRAITS <- character()
if (length(setdiff(EXPLORATORY_TRAITS, EXPLORATORY_ALL)))
	stop("Invalid BMONTH_EXPLORATORY_TRAITS; use birth-order trait names or 'none'")
ANALYSIS_EXPOSURES <- unique(c(EXPOSURES, EXPLORATORY_TRAITS))
PHEWAS_EXPOSURES <- unique(c(EXPOSURES, CONTEXT_EXPOSURES, EXPLORATORY_TRAITS))

OUTCOMES <- csv_env("BMONTH_OUTCOMES","cvd_cad,ra")
if(identical(OUTCOMES,"none"))OUTCOMES<-character()
TEMPMED_PHENOS <- csv_env("BMONTH_TEMPMED_PHENOS","")
REPLACE <- toupper(Sys.getenv("BMONTH_REPLACE","FALSE")) %in% c("TRUE","1","YES")

BASE_EXPOSURES <- c(
  "bmonth_cos","bmonth_sin","bmonth_winter_summer","bmonth_spring_autumn",
  "bmonth_rank","birth_north_100km","birth_east_100km","move_log1p")
ALL_EXPOSURES <- c(BASE_EXPOSURES, EXPLORATORY_ALL)
bad <- setdiff(unique(c(EXPOSURES,CONTEXT_EXPOSURES)),BASE_EXPOSURES)
if(length(bad)) stop("Invalid bmonth exposure(s): ",paste(bad,collapse=", "))

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


# 🚩 Exploratory birth-order traits, derived locally from the shared UKB inputs
clean_count <- function(x) {
	z <- as_num(x)
	z[!is.finite(z) | z < 0 | z != floor(z)] <- NA_real_
	z
}

make_birth_order <- function(dat) {
	required <- c("cnt_older", "cnt_brother", "cnt_sister", "adopted")
	missing <- setdiff(required, names(dat))
	if (length(missing)) {
		if (length(EXPLORATORY_TRAITS)) stop("Missing exploratory birth-order columns: ",
			paste(missing, collapse = ", "), "; use --exploratory-traits none to disable")
		return(dat)
	}
	dat %>% mutate(
		cnt_older = clean_count(cnt_older),
		cnt_brother = clean_count(cnt_brother),
		cnt_sister = clean_count(cnt_sister),
		adopted_num = if_else(as_num(adopted) %in% c(0, 1), as_num(adopted), NA_real_),
		sibling_n = cnt_brother + cnt_sister,
		sibship_size = sibling_n + 1,
		sibship_size_cap = pmin(sibship_size, 6),
		sibship_size_cat = factor(sibship_size_cap, levels = 2:6,
			labels = c("2","3","4","5","6+")),
		birth_rank = cnt_older + 1,
		impossible_rank = !is.na(cnt_older) & !is.na(sibling_n) & cnt_older > sibling_n,
		eligible_border = !is.na(cnt_older) & !is.na(cnt_brother) & !is.na(cnt_sister) &
			!impossible_rank & (is.na(adopted_num) | adopted_num != 1),
		border12 = case_when(
			!eligible_border | sibship_size != 2 ~ NA_real_,
			birth_rank == 1 ~ 0, birth_rank == 2 ~ 1, TRUE ~ NA_real_),
		border1n = case_when(
			!eligible_border | sibship_size < 2 ~ NA_real_,
			birth_rank == 1 ~ 0, birth_rank == sibship_size ~ 1, TRUE ~ NA_real_),
		border_rank = if_else(eligible_border & sibship_size >= 2, birth_rank, NA_real_),
		border_relative = if_else(eligible_border & sibship_size >= 2,
			(birth_rank - 1) / (sibship_size - 1), NA_real_)
	)
}

birth_order_qc <- function(dat) {
	if (!"eligible_border" %in% names(dat))
		return(tibble(metric = "birth-order derivation disabled: input columns unavailable", value = nrow(dat)))
	bind_rows(
		tibble(metric = "complete three sibling counts", value = sum(complete.cases(dat[, c("cnt_older", "cnt_brother", "cnt_sister")]))),
		tibble(metric = "adopted==1 excluded from birth-order traits", value = sum(dat$adopted_num == 1, na.rm = TRUE)),
		tibble(metric = "adoption status unknown (retained if otherwise eligible)", value = sum(is.na(dat$adopted_num))),
		tibble(metric = "impossible cnt_older > siblings", value = sum(dat$impossible_rank)),
		tibble(metric = "eligible birth-order N (including singletons)", value = sum(dat$eligible_border)),
		bind_rows(lapply(EXPLORATORY_ALL, function(x)
			tibble(metric = paste0(x, " nonmissing"), value = sum(!is.na(dat[[x]]))))),
		bind_rows(lapply(2:6, function(k)
			tibble(metric = paste0("sibship size ", ifelse(k == 6, "6+", k)),
				value = sum(dat$sibship_size_cap == k, na.rm = TRUE))))
	)
}

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
    bmonth_cos=-cos(theta),                 # exactly -1 in December, +1 in June
    bmonth_sin=sin(theta),
    bmonth_winter_summer=case_when(m%in%c(12,1,2)~0,m%in%c(6,7,8)~1,TRUE~NA_real_),
    bmonth_spring_autumn=case_when(m%in%c(3,4,5)~0,m%in%c(9,10,11)~1,TRUE~NA_real_),
    bmonth_rank=m,                          # sensitivity only; month is circular
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

# Five prespecified covariate stages.
# Separating D3 geography from D4 tdi.birth is deliberate: it identifies whether
# attenuation is caused by geography, birthplace deprivation proxy, or sample loss.
covars_stage <- list(
  D1_core=c("age","sex","center"),
  D2_adultSES=c("age","sex","center","edu.sco","tdi"),
  D3_geography=c("age","sex","center","edu.sco","tdi",
                 "birth_north_100km","birth_east_100km"),
  D4_birthContext=c("age","sex","center","edu.sco","tdi",
                    "birth_north_100km","birth_east_100km","tdi.birth"),
  D5_migration=c("age","sex","center","edu.sco","tdi",
                 "birth_north_100km","birth_east_100km","tdi.birth","move_log1p")
)

# Primary birth-month analysis.
covars <- covars_stage$D4_birthContext

# OLD-compatible model for birthplace-context PheWAS.
context_covars <- covars_stage$D2_adultSES
exploratory_covars <- covars_stage$D2_adultSES
exploratory_family_covars <- c(exploratory_covars, "sibship_size_cat")

analysis_covars_for <- function(exposure) {
	if (exposure %in% EXPLORATORY_ALL) exploratory_covars else covars
}
analysis_model_label <- function(exposure) {
	if (exposure %in% EXPLORATORY_ALL) "exploratory_D2_birth_order" else "primary_D4_birthContext"
}
analysis_track <- function(exposure) {
	ifelse(exposure %in% EXPLORATORY_ALL, "exploratory_birth_order",
		ifelse(exposure %in% c("birth_north_100km", "birth_east_100km"), "birthplace_context",
			ifelse(exposure == "move_log1p", "postbirth_sensitivity", "birth_month")))
}
exposure_unit <- function(exposure) {
	units <- c(border12 = "second-born vs first-born (two-child families)",
		border1n = "youngest vs oldest (middle children excluded)",
		border_rank = "per one-position later birth rank",
		border_relative = "youngest (1) vs oldest (0) on relative scale")
	ifelse(exposure %in% names(units), unname(units[exposure]), "per one-unit increase")
}

phewas_covars_for <- function(exposure){
  if(exposure %in% EXPLORATORY_ALL) return(exploratory_covars)
  if(exposure %in% c("birth_north_100km","birth_east_100km")) return(context_covars)
  if(exposure=="move_log1p") return(covars_stage$D4_birthContext)
  covars
}
phewas_model_label <- function(exposure){
  if(exposure %in% EXPLORATORY_ALL) return("exploratory_D2_birth_order")
  if(exposure %in% c("birth_north_100km","birth_east_100km"))
    return("context_D2_old_compatible")
  if(exposure=="move_log1p") return("postbirth_D4_sensitivity")
  "primary_D4_birthContext"
}

outcome_date_columns <- function(dat,y){
  p<-paste0("fod_icd10_",y); if(p%in%names(dat)) return(p)
  n<-names(dat); n[grepl("^(fod_|date_)",n)&endsWith(n,paste0("_",y))]
}
default_follow_end_path<-if(nzchar(Sys.getenv("BMONTH_FDIR","")))
  file.path(dirname(dirname(Sys.getenv("BMONTH_FDIR"))),"0f","phenotype.R") else
  file.path(dir0,"scripts","0f","phenotype.R")
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
  if(any(grepl("converg|infinite|singular|numerically 0 or 1|iteration limit",notes,ignore.case=TRUE))) value<-NULL
  list(fit=value,note=paste(unique(notes),collapse="; "))
}
load_cohort<-function(){
  p<-file.path(rawdir,"bmonth.cohort.rds")
  if(!file.exists(p)) stop("Run bmonth.sh prep first")
  d<-validate_ids(readRDS(p),"eid","bmonth cohort")
  needed<-unique(c(BASE_EXPOSURES,EXPLORATORY_TRAITS,unlist(covars_stage),
                   "birth_year_bm",if(length(FAMILY_SIZE_SENS_TRAITS))
                     c("sibship_size","sibship_size_cat")))
  if(length(setdiff(needed,names(d))) ||
     !identical(attr(d,"bmonth_metadata")$code_version,BMONTH_CODE_VERSION) ||
     !identical(attr(d,"bmonth_metadata")$covars_primary,covars))
    stop("Cohort uses an older/incompatible model; rerun bmonth.sh prep first")
  meta <- attr(d,"bmonth_metadata")
  requested <- unique(c(csv_env("BMONTH_EXTRA_COVARS"),csv_env("BMONTH_FOLLOWUP_EXTRA_COVARS"),
                        csv_env("BMONTH_OMICS_COVARIATES")))
  retained <- if(is.null(meta$extra_covariates)) character() else meta$extra_covariates
  if(!setequal(requested,retained))
    stop("Cohort covariate configuration changed; rerun bmonth.sh prep in a separate output directory")
  supplemental <- Sys.getenv("BMONTH_SUPPLEMENTAL_COVARIATES","")
  previous <- if(is.null(meta$supplemental_input)) "" else meta$supplemental_input
  if(nzchar(supplemental)!=nzchar(previous) ||
     (nzchar(supplemental) && (!file.exists(supplemental) ||
       !identical(unname(tools::md5sum(supplemental)),meta$supplemental_md5))))
    stop("Supplemental covariate input changed; rerun bmonth.sh prep in a separate output directory")
  if(nzchar(Sys.getenv("BMONTH_FOLLOW_END","")) && !is.null(meta$follow_end) &&
     !identical(as.character(follow_end()),meta$follow_end))
    stop("Follow-up end changed; rerun bmonth.sh prep in a separate output directory")
  d
}
empty_associations <- function() tibble(
  outcome=character(),exposure=character(),N=integer(),events=integer(),
  analysis_track=character(),covariate_model=character(),exposure_unit=character(),
  beta=double(),se=double(),HR=double(),conf.low=double(),conf.high=double(),p=double(),
  eligible=logical(),status=character(),note=character(),FDR=double(),
  FDR_all_prespecified=double(),bonferroni_all_prespecified=double(),
  multiplicity_note=character(),FDR_all_sensitivity=double())

read_survival_result <- function(path, required, empty=tibble()) {
  if(!length(OUTCOMES) || !file.exists(path)) return(empty)
  x <- if(grepl("[.]csv$",path)) read.csv(path,check.names=FALSE) else readRDS(path)
  if(!is.data.frame(x)) stop("Invalid survival result (expected table): ",path)
  if(!nrow(x)) return(empty)
  missing <- setdiff(c("outcome",required),names(x))
  if(length(missing)) stop("Invalid nonempty survival result: ",path,
    "; missing required columns: ",paste(missing,collapse=", "))
  as_tibble(x[x$outcome %in% OUTCOMES,,drop=FALSE])
}

fit_cox<-function(dat,exposure,y,covar_names=analysis_covars_for(exposure),
                   model_label=analysis_model_label(exposure)){
  covar_names<-setdiff(covar_names,exposure)
  dat<-ensure_outcome(dat,y);tv<-paste0(y,".t2e");ev<-paste0(y,".Yt2e")
  d<-complete_model_data(dat,unique(c(tv,ev,exposure,covar_names)))
  row<-tibble(outcome=y,exposure=exposure,N=nrow(d),events=sum(d[[ev]]==1),
    analysis_track=analysis_track(exposure),covariate_model=model_label,
    exposure_unit=exposure_unit(exposure),
    beta=NA_real_,se=NA_real_,HR=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
    p=NA_real_,eligible=FALSE,status="skipped",note="")
  if(nrow(d)<500||row$events<20||length(unique(d[[exposure]]))<2) return(row)
  row$eligible<-TRUE
  if(!exists("bm_prune_design",mode="function"))
    source(file.path(Sys.getenv("BMONTH_FDIR"),"phewas_engine.R"))
  # Include the constant in the nuisance space: Cox does not identify an
  # intercept, so an affine relationship (age = 50 + exposure) also aliases X.
  x <- bm_formula_matrix(d,c(exposure,covar_names),tested=exposure)
  constant <- attr(x,"dropped_terms")
  design <- bm_prune_design(x,tested=exposure)
  row$note <- paste(c(if(length(constant)) paste("Constant nuisance terms removed:",paste(constant,collapse=",")),
    if(length(design$dropped)) paste("Redundant nuisance columns removed:",paste(design$dropped,collapse=","))),collapse="; ")
  if(!design$ok) {
    row$status <- design$status
    row$note <- paste(row$note,"Tested exposure is not independently estimable given the nuisance space including a constant",sep="; ")
    return(row)
  }
  nuisance <- design$x[,setdiff(colnames(design$x),c("(Intercept)",exposure)),drop=FALSE]
  colnames(nuisance) <- sprintf(".cox_Z%03d",seq_len(ncol(nuisance)))
  fit_data <- cbind(d[,c(tv,ev,exposure),drop=FALSE],nuisance)
  f<-reformulate(c(exposure,colnames(nuisance)),response=paste0("Surv(`",tv,"`,`",ev,"`)"))
  fit<-model_try(coxph(f,data=fit_data,ties="efron",control=coxph.control(iter.max=50)))
  row$note <- paste(Filter(nzchar,c(row$note,fit$note)),collapse="; ")
  if(is.null(fit$fit)){row$status<-"failed";return(row)}
  sm<-coef(summary(fit$fit))
  if(!exposure%in%rownames(sm)){
    row$status<-"failed_missing_tested_term";return(row)
  }
  b<-sm[exposure,"coef"];s<-sm[exposure,"se(coef)"]
  est<-c(b,s,exp(b),exp(b-qnorm(.975)*s),exp(b+qnorm(.975)*s),sm[exposure,"Pr(>|z|)"])
  if(any(!is.finite(est))||s<=0||any(est[3:5]<=0)||est[6]<0||est[6]>1){
    row$status<-"invalid_numeric";row$note<-paste(row$note,"Non-finite/zero coefficient uncertainty or exponentiated estimate",sep="; ");return(row)
  }
  row$beta<-b;row$se<-s;row$HR<-exp(b);row$conf.low<-exp(b-qnorm(.975)*s)
  row$conf.high<-exp(b+qnorm(.975)*s);row$p<-sm[exposure,"Pr(>|z|)"];row$status<-"ok"
  row
}
annotate_phecodes <- function(x){
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
write_csv2<-function(x,name)bm_write_csv_atomic(x,file.path(rawdir,name))


# 🚩 Keep figure data alongside each rendered figure
remove_result_figure <- function(name) {
	unlink(file.path(figdir, c(paste0(name, ".png"), paste0(name, ".xlsx"))))
}
save_result_figure <- function(plot, name, results, width, height, dpi = 300) {
	remove_result_figure(name)
	if (is.null(plot) || !nrow(results)) return(invisible(NULL))
	ggplot2::ggsave(file.path(figdir, paste0(name, ".png")), plot,
		width = width, height = height, dpi = dpi, bg = "white")
	writexl::write_xlsx(list(results = as.data.frame(results)), file.path(figdir, paste0(name, ".xlsx")))
}

source(file.path(Sys.getenv("BMONTH_FDIR"),"0.analysis_contract.R"))
