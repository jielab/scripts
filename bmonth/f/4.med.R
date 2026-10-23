source(file.path(Sys.getenv("BMONTH_FDIR"), "0.common.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"), "mediation.R"))
source(file.path(Sys.getenv("BMONTH_FDIR"), "med_checkpoint.R"))
if(!length(OUTCOMES)) {
  med_atomic(empty_mediation(),file.path(rawdir,"mediation.rds"))
  write_csv2(empty_mediation(),"mediation.csv")
  write_csv2(data.frame(module="mediation",status="skipped_no_outcomes",
    detail="Survival outcomes disabled; no candidate, survival or omics inputs read"),
    "mediation_completion.csv")
  message("[mediation] skipped_no_outcomes")
  quit(save="no",status=0L)
}
base <- load_cohort()
B <- as.integer(Sys.getenv("BMONTH_MED_BOOT","100"))
if(length(B)!=1L||is.na(B)||B<0L)stop("BMONTH_MED_BOOT must be a nonnegative integer")
mode <- Sys.getenv("BMONTH_MED_MODE","screen")
if(!mode %in% c("screen","bootstrap"))stop("BMONTH_MED_MODE must be screen or bootstrap")
layers <- csv_env("BMONTH_MED_LAYERS","prot,met")
if(!length(layers)||length(setdiff(layers,c("prot","met"))))stop("BMONTH_MED_LAYERS must contain prot and/or met")
candidate_file <- Sys.getenv("BMONTH_MED_CANDIDATES","")
locked_candidates <- NULL
if(mode=="bootstrap") {
  if(!nzchar(candidate_file)||!file.exists(candidate_file))stop("Bootstrap requires an explicit BMONTH_MED_CANDIDATES CSV (layer,exposure,outcome,mediator). Run screening and review diagnostics first.")
  locked_candidates<-read.csv(candidate_file,stringsAsFactors=FALSE,colClasses="character")
  if(!all(c("layer","exposure","outcome","mediator")%in%names(locked_candidates)))stop("Invalid mediator candidate CSV schema")
  key_columns<-c("layer","exposure","outcome","mediator")
  if(anyNA(locked_candidates[,key_columns])||any(!nzchar(as.matrix(locked_candidates[,key_columns]))))stop("Missing mediator selection entries")
  if(anyDuplicated(locked_candidates[,key_columns]))stop("Duplicate mediator selections")
  if(length(setdiff(locked_candidates$layer,c("prot","met"))))stop("Unknown layer in mediator candidate CSV")
  if(length(setdiff(locked_candidates$exposure,ALL_EXPOSURES)))stop("Unknown exposure in mediator candidate CSV")
  known_outcomes<-sub("[.]Yt2e$","",grep("[.]Yt2e$",names(base),value=TRUE))
  if(length(setdiff(locked_candidates$outcome,unique(c(OUTCOMES,known_outcomes)))))stop("Unknown/unprepared outcome in mediator candidate CSV; prepare the requested outcomes first")
}
disposition<-if(is.null(locked_candidates))data.frame(layer=character(),exposure=character(),outcome=character(),mediator=character())else locked_candidates
disposition$status<-rep("not_run_selection",nrow(disposition))
disposition$status[disposition$layer%in%layers & disposition$exposure%in%ANALYSIS_EXPOSURES & disposition$outcome%in%OUTCOMES]<-"pending"
write_csv2(disposition,"mediation_candidate_disposition.csv")
if(B<50) message("[mediation] Fewer than 50 replicates: bootstrap intervals/P values will be NA")
set.seed(12345)
# Check contents rather than timestamps; regenerated identical cohorts remain reusable.
code_files <- file.path(Sys.getenv("BMONTH_FDIR"),c("0.common.R","0.analysis_contract.R","phewas_engine.R","mediation.R","4.med.R","med_checkpoint.R"))
med_signature <- list(version=2L,cohort=unname(tools::md5sum(file.path(rawdir,"bmonth.cohort.rds"))),
  code=unname(tools::md5sum(code_files)),B=B,exposures=ANALYSIS_EXPOSURES,outcomes=OUTCOMES,
  mode=mode,layers=layers,candidates=locked_candidates,
  center=Sys.getenv("BMONTH_MED_CENTER","strata"),extra_covariates=csv_env("BMONTH_OMICS_COVARIATES"),
  follow_end=if(any(!vapply(OUTCOMES,function(y) all(paste0(y,c(".t2e",".Yt2e")) %in% names(base)),logical(1)))) as.character(follow_end()) else NULL,
  R=R.version.string,packages=vapply(c("survival","dplyr","tibble"),function(p) as.character(packageVersion(p)),character(1)))
# Include BOTH layers so changing prot cannot silently reuse met's old RNG stream.
omics_files <- file.path(indir,"Rdata",paste0(layers,".rds"))
if(any(!file.exists(omics_files))) stop("Missing omics input")
med_signature$omics <- unname(tools::md5sum(omics_files))
cache <- file.path(rawdir,"mediation_cache_v2",med_hash(med_signature))
med_atomic(med_signature,file.path(cache,"signature.rds"))
message("[mediation] checkpoints: ",cache)
screen_audit <- list()
completed_parts <- list()
completion <- expand.grid(layer=layers,exposure=ANALYSIS_EXPOSURES,outcome=OUTCOMES,stringsAsFactors=FALSE)
completion$mode <- mode;completion$status <- "pending";completion$selected_candidates <- 0L;completion$completed_rows <- 0L
write_csv2(completion,"mediation_completion.csv")
med_atomic(empty_mediation(),file.path(rawdir,"mediation.rds"))
write_csv2(empty_mediation(),"mediation.csv")
cached <- function(key,compute) med_cached(cache,key,compute,REPLACE)
scan_layer <- function(layer) {
	mf <- file.path(indir,"Rdata",paste0(layer,".rds"))
	if(!file.exists(mf)) stop("Missing ",mf)
	om <- validate_ids(readRDS(mf),"eid",paste0(layer,".rds"))
	features <- setdiff(names(om),"eid")
	features <- features[vapply(om[features],is.numeric,logical(1))]
  if(mode=="bootstrap"){
    unknown<-setdiff(disposition$mediator[disposition$layer==layer & disposition$status=="pending"],features)
    if(length(unknown))stop("Unknown/non-numeric mediators in ",layer,": ",paste(unknown,collapse=", "))
  }
	if(length(intersect(features,names(base)))) stop("Omics feature names overlap cohort columns")
	d0 <- inner_join(base,om[,c("eid",features),drop=FALSE],by="eid")
	if(!nrow(d0) || !length(features)) stop("No shared participants/numeric features for ",layer)
	message("[mediation] Loaded ",layer,": ",nrow(d0)," shared participants, ",length(features),
					" features; bootstrap replicates per candidate=",B)
	allres <- list()
	for(x in ANALYSIS_EXPOSURES) {
		message("[mediation] ",layer," / ",x,": X-to-M screen (",length(features)," features)")
		# Stage A does not depend on outcome; fit only once per exposure/layer.
		xa <- cached(list("Xa",layer,x),function() bind_rows(lapply(features,function(m) screen_a(d0,x,m))))$value
		xa$FDR_a <- bm_planned_adjust(xa$p_a,xa$eligible)
		xa$n_tests_planned <- sum(xa$eligible)
		xa$n_tests_unresolved <- sum(xa$eligible & !is.finite(xa$p_a))
		xa$analysis_track <- rep(analysis_track(x),nrow(xa))
		xa$covariate_model <- rep(omics_model_label(x),nrow(xa))
		cand_a <- xa$mediator[which(xa$FDR_a<.05)]
		write_csv2(xa,paste0("med_",layer,"_",x,"_Xa.csv"))
		message("[mediation] ",layer," / ",x,": X-to-M BH FDR < 0.05: ",length(cand_a),"/",length(features)," selected")
		for(y in OUTCOMES) {
      completion_i<-which(completion$layer==layer & completion$exposure==x & completion$outcome==y)
      completion$status[completion_i]<<-"running"
      write_csv2(completion,"mediation_completion.csv")
			message("[mediation] ",layer," / ",x," / ",y,": M-to-Y screen (",length(cand_a)," FDR-selected features)")
			d<-ensure_outcome(d0,y)
			xb<-cached(list("Mb",layer,x,y),function() {
        if(length(cand_a)) bind_rows(lapply(cand_a,function(m) screen_b(d,x,m,y))) else
          tibble(mediator=character(),N_b=integer(),b=double(),se_b=double(),p_b=double(),eligible=logical(),status_b=character(),note_b=character())
      })$value
      invalid <- sum(!is.finite(xb$p_b))
      screen_audit[[paste(layer,x,y,sep="/")]] <<- tibble(layer=layer,exposure=x,outcome=y,
        selected_Xa=length(cand_a),valid_Mb=sum(is.finite(xb$p_b)),invalid_Mb=invalid,
        candidates=sum(xb$p_b<.05,na.rm=TRUE),
        status=if(!nrow(xb)) "no_Xa_candidates" else if(invalid==nrow(xb)) "all_Mb_invalid" else if(invalid) "some_Mb_invalid" else "ok")
      write_csv2(bind_rows(screen_audit),"mediation_screen_audit.csv")
      if(invalid) message("[mediation] WARNING: ",layer," / ",x," / ",y,": ",invalid,"/",nrow(xb),
        " M-to-Y fits have no valid P value; inspect note_b. These are NOT negative association results.")
			xb$FDR_b<-bm_planned_adjust(xb$p_b,xb$eligible)
			xb$analysis_track<-rep(analysis_track(x),nrow(xb))
			xb$covariate_model<-rep(omics_model_label(x),nrow(xb))
			write_csv2(xb,paste0("med_",layer,"_",x,"_",y,"_Mb.csv"))
			candidates<-xb$mediator[which(xb$p_b<.05)]
      if(mode=="screen") candidates<-character()
      else candidates<-intersect(candidates,locked_candidates$mediator[
        locked_candidates$layer==layer & locked_candidates$exposure==x & locked_candidates$outcome==y])
      requested_i<-which(disposition$layer==layer & disposition$exposure==x & disposition$outcome==y)
      if(length(requested_i))for(ii in requested_i){
        m<-disposition$mediator[ii];ai<-match(m,xa$mediator);bi<-match(m,xb$mediator)
        status<-if(is.na(ai))"unknown_mediator" else if(!xa$eligible[ai])"insufficient_Xa_data" else
          if(!is.finite(xa$p_a[ai]))"failed_Xa_model" else if(is.na(xa$FDR_a[ai])||xa$FDR_a[ai]>=.05)"not_selected_Xa_FDR" else
          if(is.na(bi))"Mb_not_evaluated" else if(!xb$eligible[bi])"insufficient_Mb_data" else
          if(!is.finite(xb$p_b[bi]))"failed_Mb_model" else if(xb$p_b[bi]>=.05)"not_selected_Mb_P" else "selected_for_bootstrap"
        disposition$status[ii]<<-status
      }
      write_csv2(disposition,"mediation_candidate_disposition.csv")
			# Only X-to-M BH-significant features can reach the mediation fits.
			stopifnot(all(candidates %in% cand_a))
			n_candidates<-length(candidates)
			message("[mediation] ",layer," / ",x," / ",y,": M-to-Y P < 0.05: ",sum(xb$p_b<.05,na.rm=TRUE),
							"/",length(cand_a),"; mode=",mode,"; selected coefficient products=",n_candidates," (bootstrap=",B,")")
			started<-proc.time()[["elapsed"]]
			# Report completed candidates, never individual bootstrap iterations.
			progress_every<-max(1L,min(5L,ceiling(n_candidates/10)))
			det<-if(n_candidates) bind_rows(lapply(seq_along(candidates),function(i) {
				checkpoint<-cached(list("candidate",layer,x,y,candidates[[i]]),
          function() mediate_feature(d,x,candidates[[i]],y,B))
        r<-checkpoint$value
        message("[mediation] ",if(checkpoint$reused) "SKIP cached: " else "SAVED: ",
          layer," / ",x," / ",y," / ",candidates[[i]])
				if(i %% progress_every == 0L || i == n_candidates)
					message("[mediation] ",layer," / ",x," / ",y,": completed ",i,"/",n_candidates,
									" candidates (",round(100*i/n_candidates),"%); elapsed ",
									round((proc.time()[["elapsed"]]-started)/60,1)," min")
				r
			})) else empty_mediation()
			det<-det %>% left_join(select(xa,mediator,a,p_a,FDR_a),by="mediator") %>%
				left_join(select(xb,mediator,b,p_b,FDR_b),by="mediator")
			det$layer<-rep(layer,nrow(det));det$FDR_indirect<-bm_planned_adjust(det$p_boot,rep(TRUE,nrow(det)))
      det$inference_label<-rep("exploratory selected coefficient products; no post-selection correction",nrow(det))
			allres[[paste(x,y)]]<-det
      completed_parts[[paste(layer,x,y,sep="/")]]<<-det
      if(length(requested_i))for(ii in requested_i){
        jj<-match(disposition$mediator[ii],det$mediator)
        if(!is.na(jj))disposition$status[ii]<<-if(det$status[jj]!="ok")"failed_final_model" else
          if(det$bootstrap_sufficient[jj])"bootstrap_complete" else "bootstrap_insufficient_valid_replicates"
      }
      write_csv2(disposition,"mediation_candidate_disposition.csv")
      published<-bind_rows(completed_parts)
      med_atomic(published,file.path(rawdir,"mediation.rds"));write_csv2(published,"mediation.csv")
      med_atomic(det,file.path(rawdir,paste0("mediation_",layer,"_",x,"_",y,".rds")))
      completion$status[completion_i]<<-if(mode=="screen")"screen_complete" else "bootstrap_complete"
      completion$selected_candidates[completion_i]<<-length(candidates)
      completion$completed_rows[completion_i]<<-nrow(det)
      write_csv2(completion,"mediation_completion.csv")
		}
	}
	bind_rows(allres)
}
res<-bind_rows(lapply(layers,scan_layer))
write_csv2(res,"mediation.csv");med_atomic(res,file.path(rawdir,"mediation.rds"))
message("[mediation] Saved ",nrow(res)," candidate results: ",file.path(rawdir,"mediation.rds"))
writeLines(c("Exploratory linear-by-Cox coefficient products; not causal natural indirect effects.",
	"X-to-M BH FDR < 0.05 within exposure/layer, then conditional M-to-Y P < 0.05.",
	"FDR_indirect is BH among selected candidates within exposure/outcome/layer; selection is not accounted for.",
	"p_boot is an empirical bootstrap sign-tail summary, not a selection-adjusted confirmatory P value.",
	"Bootstrap CI/P require at least max(50,80% of requested) valid replicates.",
	"No bootstrap repeats feature selection. A descriptive proportion is suppressed when signs disagree or the total Wald CI includes zero; reporting it does not establish mediation.",
	"Continuous exposures are per SD in each complete-case mediation sample; Cox associations use natural units.",
	"Binary season contrasts, border12 and border1n remain 0/1; other exposures are standardized.",
	"Variable-family-size birth-order traits use D2 plus categorical sibship size; border12 keeps fixed family size.",
	"Cox models are prespecified as centre-stratified by default; coefficient-centre sensitivity must use a separate output directory.",
	"X-to-M models use centre as a nuisance covariate; extra omics covariates are explicitly configured and never silently omitted.",
	"mediation_completion.csv distinguishes partial, screen-only, and completed bootstrap outputs.",
	"Default mode is screen. Bootstrap requires an explicit reviewed candidate CSV; selection is still not accounted for in inference.",
	"A full 2-df seasonal cycle has no single mediation coefficient; cyclic components are screened separately.",
	"Default 100 replicates are exploratory: minimum attainable two-sided p_boot is 2/101."),file.path(rawdir,"mediation_methods.txt"))
