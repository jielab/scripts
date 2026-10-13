source(file.path(Sys.getenv("BORDER_FDIR"), "0.common.R"))
source(file.path(Sys.getenv("BORDER_FDIR"), "mediation.R"))
base <- load_cohort()
B <- as.integer(Sys.getenv("BORDER_MED_BOOT","100"))
if(B<50) message("[mediation] Fewer than 50 replicates: bootstrap intervals/P values will be NA")
set.seed(12345)
scan_layer <- function(layer) {
  mf <- file.path(indir,"Rdata",paste0(layer,".rds"))
  if(!file.exists(mf)) stop("Missing ",mf)
  om <- validate_ids(readRDS(mf),"eid",paste0(layer,".rds"))
  features <- setdiff(names(om),"eid")
  features <- features[vapply(om[features],is.numeric,logical(1))]
  if(length(intersect(features,names(base)))) stop("Omics feature names overlap cohort columns")
  d0 <- inner_join(base,om[,c("eid",features),drop=FALSE],by="eid")
  if(!nrow(d0) || !length(features)) stop("No shared participants/numeric features for ",layer)
  message("[mediation] Loaded ",layer,": ",nrow(d0)," shared participants, ",length(features),
          " features; bootstrap replicates per candidate=",B)
  allres <- list()
  for(x in EXPOSURES) {
    message("[mediation] ",layer," / ",x,": X-to-M screen (",length(features)," features)")
    # Stage A does not depend on outcome; fit only once per exposure/layer.
    xa <- bind_rows(lapply(features,function(m) screen_a(d0,x,m)))
    xa$FDR_a <- p.adjust(xa$p_a,"BH")
    cand_a <- xa$mediator[which(xa$FDR_a<.05)]
    write_csv2(xa,paste0("med_",layer,"_",x,"_Xa.csv"))
    message("[mediation] ",layer," / ",x,": X-to-M BH FDR < 0.05: ",length(cand_a),"/",length(features)," selected")
    for(y in OUTCOMES) {
      message("[mediation] ",layer," / ",x," / ",y,": M-to-Y screen (",length(cand_a)," FDR-selected features)")
      d<-ensure_outcome(d0,y)
      xb<-if(length(cand_a)) bind_rows(lapply(cand_a,function(m) screen_b(d,x,m,y))) else
        tibble(mediator=character(),N_b=integer(),b=double(),p_b=double(),note_b=character())
      xb$FDR_b<-p.adjust(xb$p_b,"BH")
      write_csv2(xb,paste0("med_",layer,"_",x,"_",y,"_Mb.csv"))
      candidates<-xb$mediator[which(xb$p_b<.05)]
      # Only X-to-M BH-significant features can reach the mediation fits.
      stopifnot(all(candidates %in% cand_a))
      n_candidates<-length(candidates)
      message("[mediation] ",layer," / ",x," / ",y,": M-to-Y P < 0.05: ",n_candidates,
              "/",length(cand_a),"; starting mediation (bootstrap=",B,")")
      started<-proc.time()[["elapsed"]]
      # Report completed candidates, never individual bootstrap iterations.
      progress_every<-max(1L,min(5L,ceiling(n_candidates/10)))
      det<-if(n_candidates) bind_rows(lapply(seq_along(candidates),function(i) {
        r<-mediate_feature(d,x,candidates[[i]],y,B)
        if(i %% progress_every == 0L || i == n_candidates)
          message("[mediation] ",layer," / ",x," / ",y,": completed ",i,"/",n_candidates,
                  " candidates (",round(100*i/n_candidates),"%); elapsed ",
                  round((proc.time()[["elapsed"]]-started)/60,1)," min")
        r
      })) else empty_mediation()
      det<-det %>% left_join(select(xa,mediator,a,p_a,FDR_a),by="mediator") %>%
        left_join(select(xb,mediator,b,p_b,FDR_b),by="mediator")
      det$layer<-rep(layer,nrow(det));det$FDR_indirect<-p.adjust(det$p_boot,"BH")
      allres[[paste(x,y)]]<-det
    }
  }
  bind_rows(allres)
}
res<-bind_rows(scan_layer("prot"),scan_layer("met"))
write_csv2(res,"mediation.csv");saveRDS(res,file.path(rawdir,"mediation.rds"))
message("[mediation] Saved ",nrow(res)," candidate results: ",file.path(rawdir,"mediation.rds"))
writeLines(c("Exploratory linear-by-Cox coefficient products; not causal natural indirect effects.",
  "X-to-M BH FDR < 0.05 within exposure/layer, then conditional M-to-Y P < 0.05.",
  "FDR_indirect is BH among selected candidates within exposure/outcome/layer; selection is not accounted for.",
  "p_boot is an empirical bootstrap sign-tail summary, not a selection-adjusted confirmatory P value.",
  "Bootstrap CI/P require at least max(50,80% of requested) valid replicates.",
  "No bootstrap repeats feature selection. Proportion mediated can be unstable/outside [0,1].",
  "Continuous exposures are per SD in each complete-case mediation sample; Cox associations use natural units.",
  "Default 100 replicates are exploratory: minimum attainable two-sided p_boot is 2/101."),file.path(rawdir,"mediation_methods.txt"))
