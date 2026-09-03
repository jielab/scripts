# C1: prospective/prevalent associations, trajectories, clusters, and enrichment.

suppressPackageStartupMessages({
  fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
  source(file.path(fdir, "comm.f.R"))
  source(file.path(fdir, "c1_pgs.R"))
})
LE8_JOB <- "c1_correlate"
C1_CODE_VERSION     <- "2026-09-02.pgs_parallel1"
C1_SCAN_VERSION     <- "2026-08-31.riskset1" # observed-level models are unchanged

TOP_N              <- as.integer(Sys.getenv("C1_TOP_N", unset = "30"))
YY_TOP              <- as.integer(Sys.getenv("C1_YY_TOP", unset = "6"))
GRADIENT_TOP        <- as.integer(Sys.getenv("C1_GRADIENT_TOP", unset = "10"))
CLUSTER_TOP         <- as.integer(Sys.getenv("C1_CLUSTER_TOP", unset = "100"))
ATTENUATION_TOP     <- as.integer(Sys.getenv("C1_ATTENUATION_TOP", unset = "500"))
YY_BINS             <- as.integer(Sys.getenv("C1_YY_BINS", unset = "28"))
YY_MAX_YEAR         <- as.numeric(Sys.getenv("C1_YY_MAX_YEAR", unset = "16"))
GRADIENT_STEP       <- as.numeric(Sys.getenv("C1_GRADIENT_STEP", unset = "1"))
CLUSTER_STEP        <- as.numeric(Sys.getenv("C1_CLUSTER_STEP", unset = "0.5"))
MIN_BIN_N           <- as.integer(Sys.getenv("C1_MIN_BIN_N", unset = "20"))
CLUSTER_K_MAX       <- min(4L, as.integer(Sys.getenv("C1_CLUSTER_K_MAX", unset = "4")))
CLUSTER_GAP_B       <- as.integer(Sys.getenv("C1_CLUSTER_GAP_B", unset = "30"))
CLUSTER_STABILITY_B <- as.integer(Sys.getenv("C1_CLUSTER_STABILITY_B", unset = "40"))
C1_MIN_EVENT        <- 20L
C1_YY_SPAR          <- as.numeric(Sys.getenv("C1_YY_SPAR", unset = "0.62"))
C1_DIRECTION_ANCHORS<- unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
  unset="PCSK9,LPA,GDF15,NTPROBNP,MMP12"),",",fixed=TRUE)[[1]]))
C1_LANDMARK_YEARS   <- sort(unique(as.numeric(strsplit(Sys.getenv("C1_LANDMARK_YEARS",unset="0.5,1,2,5,10"),",",fixed=TRUE)[[1]])))
C1_LANDMARK_YEARS   <- C1_LANDMARK_YEARS[is.finite(C1_LANDMARK_YEARS)&C1_LANDMARK_YEARS>0]
C1_LANDMARK_TOP     <- as.integer(Sys.getenv("C1_LANDMARK_TOP",unset="500"))
C1_RISK_TOP         <- as.integer(Sys.getenv("C1_RISK_TOP",unset="16"))
C1_RISK_CUTS        <- sort(unique(as.numeric(strsplit(Sys.getenv("C1_RISK_CUTS",
  unset="0,0.5,1,2,5,10,16"),",",fixed=TRUE)[[1]])))
C1_RISK_CUTS        <- C1_RISK_CUTS[is.finite(C1_RISK_CUTS)&C1_RISK_CUTS>=0]
if(length(C1_RISK_CUTS)<2L)C1_RISK_CUTS<-c(0,.5,1,2,5,10,16)

# Publication typography. This overrides the shared theme only inside C1.
theme_5c <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size * 1.14, hjust = 0),
      plot.subtitle = element_text(face = "bold", size = base_size * .94, color = "grey30"),
      axis.title = element_text(face = "bold", size = base_size * 1.08),
      axis.text = element_text(face = "bold", size = base_size, color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size * 1.02),
      panel.grid.major.y = element_line(color = "grey91", linewidth = .25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(9, 13, 9, 13)
    )
}
forest_theme <- function(base_size = 10) theme_5c(base_size) + theme(panel.grid.major.y = element_blank())

assoc_has_results <- function(x) is.data.frame(x) && nrow(x) > 0 && any(is.finite(x$p.value))

assoc_empty_message <- function(x, case_label, min_event = C1_MIN_EVENT) {
  nt <- suppressWarnings(max(x$N_total, na.rm = TRUE)); ne <- suppressWarnings(max(x$N_event, na.rm = TRUE))
  if (!is.finite(nt)) nt <- NA_integer_; if (!is.finite(ne)) ne <- NA_integer_
  paste0("Complete-case N = ", format(nt, big.mark = ","), "\n",
         case_label, " cases = ", format(ne, big.mark = ","),
         if (is.finite(ne) && ne < min_event) paste0(" (<", min_event, " required)") else "")
}

assoc_blank_plot <- function(title, message) {
  ggplot() + annotate("text", x = 0, y = .16, label = title, fontface = "bold", size = 7) +
    annotate("text", x = 0, y = -.12, label = message, fontface = "bold", color = "grey25", size = 5) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void(base_size = 16)
}


# 🚩 Baseline prevalent association: logistic regression
make_prevalent_status <- function(dat, outcome = Y) {
  ydate <- paste0("fod_icd10_", outcome)
  if (!all(c(ydate, "date_attend") %in% names(dat))) stop("Cannot construct baseline prevalent status.", call. = FALSE)
  yd <- as.Date(dat[[ydate]]); ba <- as.Date(dat$date_attend)
  # Same-day diagnoses are deliberately excluded from the prevalent comparison.
  case_when(!is.na(yd) & !is.na(ba) & yd < ba ~ 1,
            is.na(yd) ~ 0,
            !is.na(yd) & !is.na(ba) & yd > ba ~ 0,
            TRUE ~ NA_real_)
}

logistic_scan <- function(dat, xs, covars, y, scale_x = TRUE, min_n = 500, min_case = 20) {
  xs <- intersect(xs, names(dat)); covars <- intersect(covars, names(dat))
  bind_rows(parallel_map(xs, function(x) {
    d <- dat[, unique(c(y, x, covars)), drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    nc <- sum(d[[y]] == 1, na.rm = TRUE)
    empty <- tibble(term = x, estimate = NA_real_, beta = NA_real_, std.error = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_, statistic = NA_real_, p.value = NA_real_,
                    N_total = nrow(d), N_event = nc)
    if (nrow(d) < min_n || nc < min_case || length(unique(d[[y]])) < 2) return(empty)
    d[[x]] <- suppressWarnings(as.numeric(d[[x]])); sx <- sd(d[[x]], na.rm = TRUE)
    if (!is.finite(sx) || sx <= 0) return(empty)
    if (scale_x) d[[x]] <- as.numeric(scale(d[[x]]))
    f <- reformulate(c(x, covars), response = y)
    fit <- tryCatch(glm(f, data = d, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(empty)
    sm <- coef(summary(fit)); if (!x %in% rownames(sm)) return(empty)
    b <- sm[x, "Estimate"]; se <- sm[x, "Std. Error"]
    tibble(term = x, estimate = exp(b), beta = b, std.error = se,
           conf.low = exp(b - 1.96 * se), conf.high = exp(b + 1.96 * se),
           statistic = b / se, p.value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
           N_total = nrow(d), N_event = nc)
  })) |> mutate(FDR = p.adjust(p.value, "BH")) |> arrange(p.value)
}

# Among baseline-prevalent cases, estimate whether current protein abundance
# varies with elapsed time since diagnosis.  This explicitly tests the
# disease -> protein interpretation that a time-agnostic case/control model
# cannot address.
prevalent_duration_scan <- function(dat,xs,bvar,covars,min_n=80L) {
  xs<-intersect(xs,names(dat));covars<-intersect(covars,names(dat))
  bind_rows(parallel_map(xs,function(x){
    d<-dat|>filter(is.finite(.data[[bvar]]),.data[[bvar]]<0)|>
      transmute(.protein=suppressWarnings(as.numeric(.data[[x]])),.duration=log1p(-.data[[bvar]]),across(all_of(covars)))
    d<-d[complete.cases(d),,drop=FALSE]
    empty<-tibble(term=x,beta=NA_real_,std.error=NA_real_,statistic=NA_real_,p.value=NA_real_,N_total=nrow(d))
    if(nrow(d)<min_n||sd(d$.protein)<=0||sd(d$.duration)<=0)return(empty)
    d$.protein<-as.numeric(scale(d$.protein));d$.duration<-as.numeric(scale(d$.duration))
    fit<-tryCatch(lm(reformulate(c(".duration",covars),response=".protein"),d),error=function(e)NULL)
    if(is.null(fit))return(empty);sm<-coef(summary(fit));if(!".duration"%in%rownames(sm))return(empty)
    tibble(term=x,beta=sm[".duration","Estimate"],std.error=sm[".duration","Std. Error"],
           statistic=sm[".duration","t value"],p.value=sm[".duration","Pr(>|t|)"],N_total=nrow(d))
  }))|>mutate(FDR=p.adjust(p.value,"BH"))|>arrange(p.value)
}

landmark_incident_scan <- function(dat,xs,tvar,evar,covars,landmarks=C1_LANDMARK_YEARS){
  xs<-intersect(xs,names(dat));covars<-intersect(covars,names(dat))
  map_dfr(landmarks,function(h){
    # A true landmark risk set: participants must be event-free and observed at
    # h years. Early cases and people censored before h are not controls.
    d<-dat|>filter(is.finite(.data[[tvar]]),!is.na(.data[[evar]]),.data[[tvar]]>h)|>
      mutate(.landmark_time=.data[[tvar]]-h,.landmark_event=.data[[evar]])
    z<-cox_scan(d,xs,covars,Y,time_var=".landmark_time",event_var=".landmark_event")
    z|>mutate(landmark_years=h,N_risk=max(N_total,na.rm=TRUE),events_after_landmark=max(N_event,na.rm=TRUE))
  })
}

# Formal diagnosis-window associations. Incident cases in [lo, hi) are
# compared with people known to remain disease-free and observed through hi;
# later cases therefore serve as valid controls for earlier windows. Baseline-
# prevalent duration bands use all non-prevalent participants as the reference.
risk_window_scan <- function(dat,xs,tvar,evar,bvar,covars,cuts=C1_RISK_CUTS){
  xs<-intersect(xs,names(dat));covars<-intersect(covars,names(dat))
  if(!length(xs)||length(cuts)<2L)return(tibble())
  rows<-map_dfr(seq_len(length(cuts)-1L),function(i){
    lo<-cuts[[i]];hi<-cuts[[i+1L]]
    inc_case<-dat[[evar]]==1&is.finite(dat[[tvar]])&dat[[tvar]]>=lo&dat[[tvar]]<hi
    inc_ctrl<-is.finite(dat[[tvar]])&dat[[tvar]]>=hi
    prev_case<-dat$prevalent_status==1&is.finite(dat[[bvar]])&(-dat[[bvar]])>=lo&(-dat[[bvar]])<hi
    prev_ctrl<-dat$prevalent_status==0
    one_side<-function(side,is_case,is_control,signed_mid){
      keep<-replace_na(is_case,FALSE)|replace_na(is_control,FALSE)
      base<-dat[keep,,drop=FALSE];base$.window_case<-as.integer(replace_na(is_case[keep],FALSE))
      bind_rows(parallel_map(xs,function(x){
        need<-unique(c(".window_case",x,covars));d<-base[,need,drop=FALSE];d<-d[complete.cases(d),,drop=FALSE]
        nc<-sum(d$.window_case==1);nn<-sum(d$.window_case==0)
        empty<-tibble(term=x,side,window_lo=lo,window_hi=hi,time=signed_mid,
          beta=NA_real_,std.error=NA_real_,conf.low=NA_real_,conf.high=NA_real_,
          p.value=NA_real_,N_case=nc,N_control=nn)
        if(nc<C1_MIN_EVENT||nn<100)return(empty)
        d[[x]]<-suppressWarnings(as.numeric(d[[x]]));sx<-sd(d[[x]],na.rm=TRUE)
        if(!is.finite(sx)||sx<=0)return(empty);d[[x]]<-as.numeric(scale(d[[x]]))
        fit<-tryCatch(glm(reformulate(c(x,covars),".window_case"),d,family=binomial()),error=function(e)NULL)
        if(is.null(fit))return(empty);sm<-coef(summary(fit));if(!x%in%rownames(sm))return(empty)
        b<-sm[x,"Estimate"];se<-sm[x,"Std. Error"]
        tibble(term=x,side,window_lo=lo,window_hi=hi,time=signed_mid,
          beta=b,std.error=se,conf.low=b-1.96*se,conf.high=b+1.96*se,
          p.value=sm[x,"Pr(>|z|)"],N_case=nc,N_control=nn)
      }))
    }
    bind_rows(one_side("Pre-baseline prevalent",prev_case,prev_ctrl,-(lo+hi)/2),
      one_side("Post-baseline incident",inc_case,inc_ctrl,(lo+hi)/2))
  })
  rows|>group_by(side,window_lo,window_hi)|>mutate(FDR=p.adjust(p.value,"BH"))|>ungroup()
}

plot_risk_window_scan <- function(z,anchors=C1_DIRECTION_ANCHORS){
  if(!nrow(z)||!any(is.finite(z$beta)))return(blank_plot("Risk-set diagnosis-window associations",
    "No diagnosis window had enough cases and eligible controls"))
  chosen<-unique(c(anchors,z|>filter(is.finite(p.value))|>group_by(term)|>
    summarise(best_p=min(p.value),.groups="drop")|>arrange(best_p)|>slice_head(n=C1_RISK_TOP)|>pull(term)))
  d<-z|>filter(term%in%chosen,is.finite(beta),is.finite(std.error))|>
    mutate(term=factor(term,levels=rev(chosen)),window=sprintf("%g–%g y",window_lo,window_hi))
  ggplot(d,aes(time,beta,color=side,fill=side,group=side))+
    geom_hline(yintercept=0,color="grey72")+geom_vline(xintercept=0,linetype=2,color="grey45")+
    geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=.10,color=NA)+
    geom_line(linewidth=.70)+geom_point(aes(size=N_case),alpha=.88)+
    facet_wrap(~term,scales="free_y",ncol=4)+
    scale_x_continuous(limits=c(-16,16),breaks=c(-16,-10,-5,-2,-1,0,1,2,5,10,16))+
    scale_color_manual(values=c(`Pre-baseline prevalent`="#8C6BB1",`Post-baseline incident`="#2C7FB8"))+
    scale_fill_manual(values=c(`Pre-baseline prevalent`="#8C6BB1",`Post-baseline incident`="#2C7FB8"))+
    labs(title="Risk-set associations across recorded diagnosis time",
      subtitle="Incident controls remain observed and disease-free through each window; coefficients are adjusted log odds ratios per 1-SD biomarker",
      x="Years relative to baseline blood draw",y="Adjusted beta (95% CI)",color=NULL,fill=NULL,size="Cases")+
    theme_5c(8)+theme(legend.position="bottom")
}

build_directionality_table <- function(incident,prevalent,duration,landmark=tibble(),birthline=tibble(),reverse=tibble()) {
  slim<-function(z,prefix){
    if(!nrow(z)||!all(c("term","beta","p.value","FDR")%in%names(z)))return(tibble(term=character()))
    z|>select(term,beta,p.value,FDR)|>
      rename(!!paste0("beta_",prefix):=beta,!!paste0("p_",prefix):=p.value,!!paste0("FDR_",prefix):=FDR)
  }
  landmark_target<-if(nrow(landmark)){z<-sort(unique(landmark$landmark_years[is.finite(landmark$landmark_years)]));
    if(any(z<=5))max(z[z<=5])else if(length(z))min(z)else NA_real_}else NA_real_
  lm5<-if(nrow(landmark)&&is.finite(landmark_target))landmark|>filter(landmark_years==landmark_target)|>
    select(term,beta_landmark5=beta,p_landmark5=p.value,FDR_landmark5=FDR) else
    tibble(term=character(),beta_landmark5=numeric(),p_landmark5=numeric(),FDR_landmark5=numeric())
  full_join(slim(incident,"incident"),slim(prevalent,"prevalent"),by="term")|>
    full_join(slim(duration,"duration"),by="term")|>full_join(lm5,by="term")|>
    full_join(slim(birthline,"birthline"),by="term")|>full_join(slim(reverse,"reverse"),by="term")|>
    mutate(incident_support=is.finite(FDR_incident)&FDR_incident<.05,
           prevalent_support=is.finite(FDR_prevalent)&FDR_prevalent<.05,
           duration_support=is.finite(FDR_duration)&FDR_duration<.05,
           distal_support=is.finite(FDR_landmark5)&FDR_landmark5<.05,
           birthline_support=is.finite(FDR_birthline)&FDR_birthline<.05,
           reactive_compatible=prevalent_support&(duration_support|(incident_support&!distal_support)),
           direction_class=case_when(
             distal_support&prevalent_support~"distal + disease-state (mixed)",
             distal_support&!prevalent_support~"distal antecedent supported",
             reactive_compatible~"near-diagnosis / reactive-compatible",
             !incident_support&prevalent_support~"established-disease associated",
             incident_support~"incident association only",
             TRUE~"unresolved"),
           incident_score=sign(beta_incident)*pmin(12,-log10(pmax(p_incident,1e-300))),
           prevalent_score=sign(beta_prevalent)*pmin(12,-log10(pmax(p_prevalent,1e-300))),
           duration_score=sign(beta_duration)*pmin(12,-log10(pmax(p_duration,1e-300))),
           landmark5_score=sign(beta_landmark5)*pmin(12,-log10(pmax(p_landmark5,1e-300))),
           birthline_score=sign(beta_birthline)*pmin(12,-log10(pmax(p_birthline,1e-300))),
           reverse_score=sign(beta_reverse)*pmin(12,-log10(pmax(p_reverse,1e-300))))|>
    arrange(match(direction_class,c("near-diagnosis / reactive-compatible","distal + disease-state (mixed)",
      "distal antecedent supported","incident association only","established-disease associated","unresolved")),p_incident)
}

plot_directionality_triage <- function(z,anchors=C1_DIRECTION_ANCHORS) {
  if(!nrow(z))return(blank_plot("Temporal directionality triage","No estimable protein results"))
  pal<-c("near-diagnosis / reactive-compatible"="#B45A4D","distal + disease-state (mixed)"="#8C6BB1",
         "distal antecedent supported"="#2B8CBE","incident association only"="#74A9CF",
         "established-disease associated"="#B38A3E","unresolved"="grey75")
  lab<-z|>filter(term%in%anchors)|>bind_rows(z|>filter(reactive_compatible)|>slice_min(p_incident,n=6))|>distinct(term,.keep_all=TRUE)
  pa<-ggplot(z,aes(beta_incident,beta_prevalent,color=direction_class))+geom_hline(yintercept=0,color="grey85")+geom_vline(xintercept=0,color="grey85")+
    geom_point(alpha=.65,size=1.8)+geom_text_repel(data=lab,aes(label=term),size=3,fontface="bold",max.overlaps=Inf)+
    scale_color_manual(values=pal)+labs(title="a. Incident versus baseline-prevalent association",
      subtitle="Prevalent association is disease-state compatible, not proof that disease caused the biomarker",
      x="Incident Cox beta",y="Prevalent logistic beta",color=NULL)+theme_5c(9)
  pbdat<-z|>filter(is.finite(beta_landmark5),is.finite(beta_duration))
  pb<-if(!nrow(pbdat))blank_plot("b. Distal versus established-disease evidence","Landmark or case-only duration result was unavailable")else
    ggplot(pbdat,aes(beta_landmark5,beta_duration,color=direction_class))+geom_hline(yintercept=0,color="grey85")+geom_vline(xintercept=0,color="grey85")+
      geom_point(alpha=.7,size=1.8)+geom_text_repel(data=pbdat|>filter(term%in%anchors),aes(label=term),size=3.2,fontface="bold",max.overlaps=Inf)+
      scale_color_manual(values=pal)+labs(title="b. Five-year landmark versus case-only duration",
        subtitle="Far-horizon persistence argues against a purely near-diagnosis signal; duration remains survivor/treatment sensitive",
        x="Incident Cox beta after 5-year landmark",y="Years-since-diagnosis slope",color=NULL)+theme_5c(9)
  (pa|pb)+plot_layout(guides="collect") &
    theme(legend.position="bottom",legend.box="vertical")
}

plot_directionality_supplement <- function(z,anchors=C1_DIRECTION_ANCHORS) {
  if(!nrow(z))return(blank_plot("Direction-of-time evidence detail","No estimable biomarker results"))
  pal<-c("near-diagnosis / reactive-compatible"="#B45A4D","distal + disease-state (mixed)"="#8C6BB1",
         "distal antecedent supported"="#2B8CBE","incident association only"="#74A9CF",
         "established-disease associated"="#B38A3E","unresolved"="grey75")
  chosen<-unique(c(anchors,z|>filter(reactive_compatible|distal_support)|>slice_min(p_incident,n=14)|>pull(term)))
  heat<-z|>filter(term%in%chosen)|>select(term,incident_score,birthline_score,landmark5_score,prevalent_score,duration_score)|>
    pivot_longer(-term,names_to="analysis",values_to="signed_evidence")|>
    mutate(analysis=factor(analysis,levels=c("incident_score","birthline_score","landmark5_score","prevalent_score","duration_score"),
                           labels=c("Incident Cox","Attained-age Cox","5-y landmark","Prevalent logistic","Duration slope")))
  pc<-ggplot(heat,aes(analysis,factor(term,levels=rev(chosen)),fill=signed_evidence))+geom_tile(color="white")+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-12,12),name="signed\n-log10(P)")+
    labs(title="a. Signed evidence by temporal model",x=NULL,y=NULL)+theme_5c(8)+theme(axis.text.x=element_text(angle=35,hjust=1))
  pd<-z|>count(direction_class,name="biomarkers")|>mutate(direction_class=factor(direction_class,levels=names(pal)))|>
    ggplot(aes(biomarkers,direction_class,fill=direction_class))+geom_col(width=.68)+geom_text(aes(label=biomarkers),hjust=-.12,fontface="bold")+
    scale_x_continuous(expand=expansion(mult=c(0,.18)))+scale_fill_manual(values=pal,guide="none")+
    labs(title="b. Full-screen evidence classes",x="Biomarkers",y=NULL)+theme_5c(9)
  pc|pd
}

plot_landmark_and_attained_age <- function(landmark,incident,birthline,anchors=C1_DIRECTION_ANCHORS) {
  top_landmark<-if(nrow(landmark))landmark|>filter(is.finite(p.value))|>group_by(term)|>
    summarise(best_p=min(p.value),.groups="drop")|>slice_min(best_p,n=8,with_ties=FALSE)|>pull(term) else character()
  chosen<-unique(c(anchors,top_landmark));lm<-landmark|>filter(term%in%chosen,is.finite(beta),is.finite(std.error))|>
    mutate(term=factor(term,levels=rev(chosen)),lo=beta-1.96*std.error,hi=beta+1.96*std.error)
  pa<-if(!nrow(lm))blank_plot("a. Landmark persistence","No landmark estimate was available")else
    ggplot(lm,aes(landmark_years,beta,color=term,group=term))+geom_hline(yintercept=0,linetype=3,color="grey65")+
      geom_ribbon(aes(ymin=lo,ymax=hi,fill=term),alpha=.08,color=NA)+geom_line(linewidth=.85)+geom_point(size=1.8)+
      scale_x_continuous(breaks=C1_LANDMARK_YEARS)+labs(title="a. Landmark persistence",
        subtitle="Risk sets begin 0.5, 1, 2, 5 or 10 years after baseline; estimates remain prospective",
        x="Event-free landmark after baseline (years)",y="Log HR per 1-SD baseline biomarker",color=NULL,fill=NULL)+
      theme_5c(9)+theme(legend.position="bottom")
  z<-incident|>select(term,beta_baseline=beta,p_baseline=p.value)|>
    inner_join(birthline|>select(term,beta_attained=beta,p_attained=p.value),by="term")|>
    filter(is.finite(beta_baseline),is.finite(beta_attained))
  rr<-if(nrow(z)>2)suppressWarnings(cor(z$beta_baseline,z$beta_attained,use="complete.obs"))else NA_real_
  lab<-z|>filter(term%in%anchors)|>bind_rows(z|>mutate(delta=abs(beta_attained-beta_baseline))|>slice_max(delta,n=6))|>distinct(term,.keep_all=TRUE)
  pb<-if(!nrow(z))blank_plot("b. Time-scale sensitivity","Attained-age result was unavailable")else
    ggplot(z,aes(beta_baseline,beta_attained))+geom_abline(slope=1,intercept=0,linetype=2,color="grey55")+
      geom_hline(yintercept=0,color="grey88")+geom_vline(xintercept=0,color="grey88")+geom_point(alpha=.55,color="#3F78A8")+
      ggrepel::geom_text_repel(data=lab,aes(label=term),size=2.8,fontface="bold",max.overlaps=Inf,seed=91)+
      annotate("text",x=Inf,y=-Inf,label=sprintf("r = %.3f",rr),hjust=1.08,vjust=-.5,fontface="bold",size=3.2)+
      labs(title="b. Baseline-time versus attained-age Cox",
        subtitle="Attained-age model uses delayed entry at blood draw; this is not a birth-cohort omics measurement",
        x="Baseline-time Cox beta",y="Attained-age delayed-entry Cox beta")+theme_5c(9)
  pa|pb
}

plot_reverse_time_exploratory <- function(reverse,prevalent,duration,anchors=C1_DIRECTION_ANCHORS) {
  z<-reverse|>select(term,beta_reverse=beta,p_reverse=p.value)|>
    left_join(prevalent|>select(term,beta_prevalent=beta),by="term")|>
    left_join(duration|>select(term,beta_duration=beta),by="term")
  lab<-z|>filter(term%in%anchors)|>bind_rows(z|>filter(is.finite(p_reverse))|>slice_min(p_reverse,n=8))|>distinct(term,.keep_all=TRUE)
  one<-function(y,title,ylab){
    d<-z|>filter(is.finite(beta_reverse),is.finite(.data[[y]]));ll<-lab|>semi_join(d,by="term")
    if(!nrow(d))return(blank_plot(title,"No exploratory reverse-time estimate was available"))
    ggplot(d,aes(beta_reverse,.data[[y]]))+geom_hline(yintercept=0,color="grey88")+geom_vline(xintercept=0,color="grey88")+
      geom_point(alpha=.58,color="#8C6BB1")+ggrepel::geom_text_repel(data=ll,aes(label=term),size=2.8,fontface="bold",max.overlaps=Inf,seed=92)+
      labs(title=title,x="Legacy reverse-time Cox beta",y=ylab)+theme_5c(9)
  }
  (one("beta_prevalent","a. Reverse-time versus prevalence","Prevalent logistic beta")|
     one("beta_duration","b. Reverse-time versus case duration","Years-since-diagnosis slope"))+
    plot_annotation(title="Exploratory only: incompatible reverse-time risk origins",
      subtitle="Cases use time since diagnosis, controls use age since birth; these panels are excluded from directionality grading")
}

# Same-N attenuation: both Cox models use the complete adj2 sample for each biomarker.
same_sample_attenuation <- function(dat, features, tvar, evar, covs_basic, covs_adj2) {
  min_basic_beta <- as.numeric(Sys.getenv("C1_ATTENUATION_MIN_ABS_BETA", unset = "0.02"))
  if(!is.finite(min_basic_beta)||min_basic_beta<=0)min_basic_beta<-0.02
  features <- intersect(features, names(dat))
  bind_rows(parallel_map(features, function(x) {
    need <- unique(c(tvar, evar, x, covs_adj2)); d <- dat[, need, drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    ne <- sum(d[[evar]] == 1, na.rm = TRUE)
    if (nrow(d) < 500 || ne < C1_MIN_EVENT) return(tibble(term=x,N=nrow(d),events=ne,
      beta_basic_sameN=NA_real_,p_basic_sameN=NA_real_,beta_adj2_sameN=NA_real_,p_adj2_sameN=NA_real_,
      effect_ratio=NA_real_,log2_effect_ratio=NA_real_,sign_flip=NA,attenuation_pct=NA_real_))
    d[[x]] <- as.numeric(scale(suppressWarnings(as.numeric(d[[x]]))))
    fit_one <- function(covars) {
      f <- as.formula(paste0("Surv(", bt(tvar), ",", bt(evar), ") ~ ", paste(bt(c(x, covars)), collapse=" + ")))
      fit <- tryCatch(coxph(f, d, ties="efron"), error=function(e) NULL)
      if (is.null(fit)) return(c(beta=NA_real_,p=NA_real_))
      sm <- coef(summary(fit)); if (!x %in% rownames(sm)) return(c(beta=NA_real_,p=NA_real_))
      c(beta=sm[x,"coef"], p=sm[x,"Pr(>|z|)"])
    }
    a <- fit_one(covs_basic); b <- fit_one(covs_adj2)
    stable<-is.finite(a[["beta"]])&&abs(a[["beta"]])>=min_basic_beta
    ratio<-if(stable)abs(b[["beta"]])/abs(a[["beta"]]) else NA_real_
    tibble(term=x,N=nrow(d),events=ne,beta_basic_sameN=a[["beta"]],p_basic_sameN=a[["p"]],
           beta_adj2_sameN=b[["beta"]],p_adj2_sameN=b[["p"]],
           effect_ratio=ratio,log2_effect_ratio=ifelse(is.finite(ratio)&ratio>0,log2(ratio),NA_real_),
           sign_flip=ifelse(is.finite(a[["beta"]])&is.finite(b[["beta"]]),sign(a[["beta"]])!=sign(b[["beta"]]),NA),
           attenuation_pct=ifelse(is.finite(ratio),100*(1-ratio),NA_real_))
  })) |> arrange(p_basic_sameN)
}


# 🚩 Yin-Yang: baseline = 0, prevalent < 0, incident > 0
residualize_z <- function(d, x, covars) {
  covars <- intersect(covars, names(d)); z <- suppressWarnings(as.numeric(d[[x]]))
  if (!length(covars)) return(as.numeric(scale(z)))
  dd <- d[, c(x, covars), drop=FALSE]; dd[[x]] <- z
  f <- reformulate(covars, response=x)
  fit <- tryCatch(lm(f, dd, na.action=na.exclude), error=function(e) NULL)
  if (is.null(fit)) return(as.numeric(scale(z)))
  as.numeric(scale(residuals(fit)))
}

global_trajectory_p <- function(time, z, df = 3) {
  ok <- is.finite(time) & is.finite(z)
  if (sum(ok) < 80 || length(unique(time[ok])) < 8) return(NA_real_)
  d <- data.frame(time=time[ok], z=z[ok])
  f0 <- tryCatch(lm(z ~ 1, d), error=function(e) NULL)
  f1 <- tryCatch(lm(z ~ splines::ns(time, df=df), d), error=function(e) NULL)
  if (is.null(f0) || is.null(f1)) return(NA_real_)
  a <- tryCatch(anova(f0,f1), error=function(e) NULL)
  if (is.null(a) || nrow(a) < 2) NA_real_ else as.numeric(a$`Pr(>F)`[2])
}

make_yy_panel <- function(dat, features, bvar, covars, title, ylab="Mean biomarker z-score",
                          bins=YY_BINS, max_year=YY_MAX_YEAR, min_bin_n=MIN_BIN_N,
                          smooth_lines=TRUE) {
  features <- intersect(features, names(dat)); if (!length(features)) return(blank_plot(title,"No feature available"))
  d0 <- dat |> filter(is.finite(.data[[bvar]]))
  lim <- min(max_year, max(abs(d0[[bvar]]), na.rm=TRUE)); if (!is.finite(lim) || lim <= 0) lim <- max_year
  br <- seq(-lim, lim, length.out=bins+1); mids <- (head(br,-1)+tail(br,-1))/2
  hist0 <- hist(d0[[bvar]], breaks=br, plot=FALSE)
  h <- tibble(mid=mids, xmin=head(br,-1), xmax=tail(br,-1), n=hist0$counts,
              side=ifelse(mids<0,"Prevalent (Yang)","Incident (Yin)"))
  rows <- map_dfr(features, function(x) {
    need <- unique(c(bvar,x,covars)); dd <- d0[, need, drop=FALSE]; dd <- dd[complete.cases(dd),,drop=FALSE]
    if (nrow(dd)<80) return(tibble())
    dd$z <- residualize_z(dd,x,covars)
    dd$bin <- cut(dd[[bvar]], br, include.lowest=TRUE, labels=FALSE)
    pg <- global_trajectory_p(dd[[bvar]],dd$z)
    sm <- dd |> filter(is.finite(z),!is.na(bin)) |> group_by(bin) |>
      summarise(mean=mean(z),sd=sd(z),se=sd/sqrt(n()),N=n(),.groups="drop") |>
      filter(N>=min_bin_n) |> mutate(time=mids[bin],feature=x,p_global=pg,N_total=nrow(dd))
    # Smooth only the observed bin summaries.  A weighted smoothing spline is
    # deterministic, does not extrapolate, and prevents small extreme bins from
    # dominating the visual trajectory.  Raw bin means remain in the output.
    if (nrow(sm) >= 5) {
      fit <- tryCatch(stats::smooth.spline(sm$time, sm$mean, w=sqrt(sm$N), spar=C1_YY_SPAR),
                      error=function(e) NULL)
      sm$mean_smooth <- if (is.null(fit)) sm$mean else
        as.numeric(predict(fit, x=sm$time)$y)
    } else sm$mean_smooth <- sm$mean
    sm
  })
  if (!nrow(rows)) return(blank_plot(title,"No bin reached the minimum sample size"))
  labs0 <- rows |> group_by(feature) |> summarise(p_global=first(p_global),N_total=first(N_total),.groups="drop") |>
    mutate(label=sprintf("%s (N=%s, Pglobal=%s)",feature,format(N_total,big.mark=","),
                         ifelse(is.finite(p_global),format.pval(p_global,digits=2,eps=1e-99),"NA")))
  rows <- rows |> left_join(labs0 |> select(feature,label),by="feature")
  yr <- range(rows$mean,na.rm=TRUE); if(!all(is.finite(yr))||diff(yr)<.2) yr <- c(-1,1)
  # Put the actual case-count histogram on the primary scale and biomarker
  # trajectories on the secondary scale, matching yy.R Fig1 panel a.
  pad <- .14*diff(yr); yr <- yr+c(-pad,pad); hmax <- max(h$n,1)
  plot_max <- hmax*1.12
  scale_factor <- plot_max/diff(yr)
  z_to_count <- function(z)(z-yr[1])*scale_factor
  rows <- rows |> mutate(line_value=if(smooth_lines)mean_smooth else mean,
                         y_plot=z_to_count(line_value), y_raw=z_to_count(mean))
  ggplot() +
    geom_col(data=h,aes(mid,n,fill=side),width=diff(br)[1]*.98,color="white",alpha=.55,inherit.aes=FALSE) +
    geom_vline(xintercept=0,linetype=2,color="grey38",linewidth=.75) +
    geom_hline(yintercept=z_to_count(0),linetype=3,color="grey62") +
    geom_line(data=rows,aes(time,y_plot,color=label,group=label),linewidth=1.05) +
    geom_point(data=rows,aes(time,y_raw,color=label),size=1.65,
               alpha=if(smooth_lines).42 else .90) +
    scale_fill_manual(values=c(`Prevalent (Yang)`="grey58",`Incident (Yin)`="#78B7E5"),guide="none") +
    scale_y_continuous(name="Patient count",limits=c(0,plot_max),
      sec.axis=sec_axis(~./scale_factor+yr[1],name=ylab)) +
    labs(title=title,
         subtitle=if(smooth_lines)
           "Weighted smoothing-spline lines; circles are observed bin means"
         else "Original display: solid circles are observed bin means joined without smoothing",
         x="Years relative to baseline",color=NULL) +
    coord_cartesian(xlim=c(-lim,lim),clip="off") + theme_5c(11) +
    guides(color=guide_legend(ncol=2,byrow=TRUE,override.aes=list(linewidth=1.1,size=2.8))) +
    theme(legend.position=c(.985,.985),legend.justification=c(1,1),legend.direction="horizontal",legend.box="horizontal",
          legend.background=element_rect(fill=scales::alpha("white",.82),color="grey75",linewidth=.25),
          legend.text=element_text(size=8.2,face="bold"),legend.margin=margin(3,4,3,4))
}


# 🚩 Fine temporal gradient (no extrapolated smoothed curve)
make_gradient_panel <- function(dat, features, bvar, side=c("incident","prevalent"), step=.25,
                                min_bin_n=20, title=NULL, covars=character(),
                                max_year=YY_MAX_YEAR) {
  side <- match.arg(side); features <- intersect(features,names(dat))
  # Selection differs by panel (incident top ten versus prevalent top ten), but
  # both panels deliberately show the complete -16..+16 diagnosis-time domain.
  # This is a cross-sectional pseudo-trajectory from one baseline sample per
  # person, not longitudinal within-person change.
  dd <- dat |> filter(is.finite(.data[[bvar]]),abs(.data[[bvar]])<=max_year) |>
    mutate(.plot_time=.data[[bvar]])
  domain_note <- if(side=="incident")"Features ranked by incident Cox P" else "Features ranked by baseline-prevalent logistic P"
  rows <- map_dfr(features,function(x){
    x0 <- suppressWarnings(as.numeric(dd[[x]])); if(sum(is.finite(x0))<50)return(tibble())
    z <- residualize_z(dd,x,covars)
    tb <- pmax(-max_year+step/2,pmin(max_year-step/2,
      floor((dd$.plot_time+max_year)/step)*step-max_year+step/2))
    tibble(time=tb,z=z) |> filter(is.finite(time),is.finite(z)) |> group_by(time) |>
      summarise(mean=mean(z),N=n(),.groups="drop") |>
      mutate(feature=x,reliable=N>=min_bin_n)
  })
  if(!nrow(rows)) return(blank_plot(title %||% "Temporal gradient","No sufficiently populated time bins"))
  ord <- rev(features[features %in% unique(rows$feature)]); rows$feature <- factor(rows$feature,levels=ord)
  ggplot(rows,aes(time,feature)) +
    geom_vline(xintercept=0,linetype=2,color="grey45") +
    geom_point(aes(size=pmin(N,200),fill=mean,alpha=reliable),shape=21,color="grey55",stroke=.15) +
    scale_fill_gradient2(low="#3B78A8",mid="#F7F7F7",high="#C86B4A",midpoint=0,name="Mean z") +
    scale_size_continuous(range=c(.6,4),name="Bin N") +
    scale_alpha_manual(values=c(`TRUE`=.92,`FALSE`=.35),guide="none")+
    scale_x_continuous(limits=c(-max_year,max_year),breaks=seq(-max_year,max_year,by=4),expand=expansion(mult=c(.01,.01))) +
    labs(title=title,subtitle=paste0(domain_note,"; adj2-residualized; faint circles have N < ",min_bin_n),
         x="Years of recorded diagnosis relative to baseline blood draw",y=NULL) + theme_5c(10) +
    theme(legend.position="right",axis.text.y=element_text(size=8.5,face="bold"))
}


# 🚩 Publication-style YY circle matrix
# Used as the first row of Fig4. No line joins adjacent circles: each dot is an
# observed time-bin summary.
make_gradient_circle_panel <- function(dat,features,bvar,step=GRADIENT_STEP,min_bin_n=MIN_BIN_N,
                                       max_year=YY_MAX_YEAR,title="a. Temporal biomarker gradients") {
  features<-intersect(features,names(dat)); d<-dat|>filter(is.finite(.data[[bvar]]),abs(.data[[bvar]])<=max_year)
  rows<-map_dfr(features,function(x){
    x0<-suppressWarnings(as.numeric(d[[x]])); if(sum(is.finite(x0))<50)return(tibble())
    z<-as.numeric(scale(x0)); tb<-floor(d[[bvar]]/step)*step+step/2
    tibble(time=tb,z=z)|>filter(is.finite(time),is.finite(z))|>group_by(time)|>
      summarise(mean=mean(z),N=n(),.groups="drop")|>filter(N>=min_bin_n)|>mutate(feature=x)
  })
  if(!nrow(rows))return(blank_plot(title,"No sufficiently populated time bins"))
  # Keep strongest proteins at the top, matching the conventional circle-matrix style.
  ord<-rev(features[features%in%unique(rows$feature)]);rows$feature<-factor(rows$feature,levels=ord)
  ggplot(rows,aes(time,feature))+
    geom_vline(xintercept=0,linetype=2,color="grey38",linewidth=.7)+
    geom_point(aes(size=N,fill=mean),shape=21,color="grey72",stroke=.22,alpha=.95)+
    scale_fill_gradient2(low="#2C7FB8",mid="white",high="#D7301F",midpoint=0,name="Mean z")+
    scale_size_continuous(range=c(.8,4.6),name="Bin N")+
    scale_x_continuous(breaks=pretty(c(-max_year,max_year),n=7))+
    labs(title=title,x="Years relative to baseline",y=NULL)+theme_5c(10)+
    theme(legend.position="right",axis.text.y=element_text(size=8.6,face="bold"))
}


# 🚩 Data-driven clustering: silhouette + gap + temporal-bin bootstrap stability
adjusted_rand <- function(a,b) {
  tab <- table(a,b); n <- sum(tab); if(n<2)return(NA_real_)
  c2 <- function(x) x*(x-1)/2
  nij <- sum(c2(tab)); ai <- sum(c2(rowSums(tab))); bj <- sum(c2(colSums(tab))); tot <- c2(n)
  expected <- ai*bj/tot; denom <- .5*(ai+bj)-expected
  if(!is.finite(denom)||denom==0) return(NA_real_)
  (nij-expected)/denom
}

mean_silhouette <- function(distmat, cl) {
  n <- length(cl); if(n<3 || length(unique(cl))<2)return(NA_real_)
  s <- vapply(seq_len(n),function(i){
    same <- which(cl==cl[i] & seq_len(n)!=i); a <- if(length(same))mean(distmat[i,same]) else 0
    oth <- setdiff(unique(cl),cl[i]); b <- min(vapply(oth,function(g)mean(distmat[i,cl==g]),numeric(1)))
    if(max(a,b)==0)0 else (b-a)/max(a,b)
  },numeric(1))
  mean(s,na.rm=TRUE)
}

within_ss <- function(x,cl) {
  sum(vapply(unique(cl),function(g){z<-x[cl==g,,drop=FALSE]; if(nrow(z)<2)return(0);sum(rowSums((z-matrix(colMeans(z),nrow(z),ncol(z),byrow=TRUE))^2))},numeric(1)))
}

trajectory_matrix <- function(dat,features,bvar,step=.5,min_bin_n=15,max_year=YY_MAX_YEAR) {
  features<-intersect(features,names(dat)); d<-dat|>filter(is.finite(.data[[bvar]]),abs(.data[[bvar]])<=max_year)
  grid<-seq(-max_year+step/2,max_year-step/2,by=step)
  raw<-map_dfr(features,function(x){
    z<-as.numeric(scale(suppressWarnings(as.numeric(d[[x]])))); tb<-floor((d[[bvar]]+max_year)/step)*step-max_year+step/2
    tibble(time=tb,z=z)|>filter(is.finite(z),is.finite(time))|>group_by(time)|>summarise(mean=mean(z),N=n(),.groups="drop")|>
      filter(N>=min_bin_n)|>mutate(feature=x)
  })
  if(!nrow(raw))return(NULL)
  wide<-raw|>select(feature,time,mean)|>pivot_wider(names_from=time,values_from=mean)
  rn<-wide$feature; X<-as.matrix(wide[,-1,drop=FALSE]);rownames(X)<-rn
  # Internal interpolation only for clustering. Plotted raw trajectories are never extrapolated.
  for(i in seq_len(nrow(X))){
    ok<-which(is.finite(X[i,])); if(length(ok)>=2){
      miss<-which(!is.finite(X[i,]) & seq_len(ncol(X))>=min(ok) & seq_len(ncol(X))<=max(ok))
      if(length(miss))X[i,miss]<-approx(ok,X[i,ok],xout=miss,rule=1)$y
    }
    m<-mean(X[i,],na.rm=TRUE); if(!is.finite(m))m<-0; X[i,!is.finite(X[i,])]<-m
  }
  X<-t(scale(t(X)));X[!is.finite(X)]<-0
  list(X=X,raw=raw,grid=grid)
}

choose_cluster_k <- function(X,kmax=4,gap_B=30,stability_B=40,seed=SEED) {
  n<-nrow(X); ks<-2:min(kmax,n-1); if(!length(ks))return(list(k=1,metrics=tibble()))
  D<-as.matrix(dist(X)); hc<-hclust(as.dist(D),method="ward.D2")
  sil<-vapply(ks,function(k)mean_silhouette(D,cutree(hc,k)),numeric(1))
  # Gap statistic on a uniform reference box with the same feature dimensions.
  set.seed(seed); obsW<-vapply(ks,function(k)within_ss(X,cutree(hc,k)),numeric(1)); obsW<-pmax(obsW,1e-12)
  ref_logW<-matrix(NA_real_,nrow=gap_B,ncol=length(ks))
  lo<-apply(X,2,min);hi<-apply(X,2,max)
  for(b in seq_len(gap_B)){
    XR<-sweep(matrix(runif(n*ncol(X)),nrow=n),2,hi-lo,"*");XR<-sweep(XR,2,lo,"+")
    hR<-hclust(dist(XR),method="ward.D2")
    ref_logW[b,]<-vapply(ks,function(k)log(pmax(within_ss(XR,cutree(hR,k)),1e-12)),numeric(1))
  }
  gap<-colMeans(ref_logW,na.rm=TRUE)-log(obsW); gap_se<-apply(ref_logW,2,sd,na.rm=TRUE)*sqrt(1+1/gap_B)
  # Stability: resample time dimensions, cluster again, compare with original labels by ARI.
  stab<-vapply(ks,function(k){
    orig<-cutree(hc,k); z<-rep(NA_real_,stability_B)
    for(b in seq_len(stability_B)){
      # Sparse outcomes can leave only one or two populated temporal bins.
      # Never request more columns than exist when bootstrapping stability.
      n_take<-min(ncol(X),max(1L,ceiling(.8*ncol(X))))
      cols<-sort(sample(seq_len(ncol(X)),n_take,replace=FALSE))
      cb<-cutree(hclust(dist(X[,cols,drop=FALSE]),method="ward.D2"),k)
      z[b]<-adjusted_rand(orig,cb)
    }
    median(z,na.rm=TRUE)
  },numeric(1))
  min_cluster_n<-vapply(ks,function(k)min(table(cutree(hc,k))),numeric(1))
  min_allowed<-max(3L,ceiling(.05*n))
  met<-tibble(k=ks,silhouette=sil,gap=gap,gap_se=gap_se,stability=stab,
              min_cluster_n=min_cluster_n,passes_min_cluster=min_cluster_n>=min_allowed)
  rank_metric<-function(x)rank(-ifelse(is.finite(x),x,-Inf),ties.method="min")
  met<-met|>mutate(rank_sum=rank_metric(silhouette)+rank_metric(gap)+rank_metric(stability))
  eligible<-met|>filter(passes_min_cluster)
  if(!nrow(eligible))eligible<-met|>filter(k==min(k))
  best<-eligible|>arrange(rank_sum,desc(silhouette),desc(stability),desc(gap))|>slice(1)|>pull(k)
  list(k=best,metrics=met,hc=hc)
}

make_cluster_figure <- function(dat,features,bvar,step=CLUSTER_STEP,max_year=YY_MAX_YEAR,
                                side=c("incident","prevalent"),panel_label="b",fixed_k=NULL) {
  side<-match.arg(side)
  train_dat<-dat|>filter(if(side=="incident") .data[[bvar]]>0 else .data[[bvar]]<0)
  tm<-trajectory_matrix(train_dat,features,bvar,step=step,min_bin_n=max(10L,floor(MIN_BIN_N*.75)),max_year=max_year)
  if(is.null(tm)||nrow(tm$X)<3){
    p0<-blank_plot(paste0(panel_label,". ",ifelse(side=="incident","Incident","Baseline-prevalent")," trajectory clusters"),
                   "Insufficient trajectories")
    return(list(plot=p0,cluster_plot=p0,diagnostic_plot=p0,cluster=tibble(),metrics=tibble(),raw=tibble(),k=1L))
  }
  ck<-choose_cluster_k(tm$X,CLUSTER_K_MAX,CLUSTER_GAP_B,CLUSTER_STABILITY_B)
  if(!is.null(fixed_k))ck$k<-min(max(1L,as.integer(fixed_k)),nrow(tm$X))
  cl<-tibble(feature=rownames(tm$X),cluster=cutree(ck$hc,ck$k))
  full_tm<-trajectory_matrix(dat,cl$feature,bvar,step=step,min_bin_n=max(10L,floor(MIN_BIN_N*.75)),max_year=max_year)
  lines<-(full_tm$raw%||%tm$raw)|>inner_join(cl,by="feature")
  cnt<-cl|>count(cluster,name="proteins")
  lines<-lines|>left_join(cnt,by="cluster")|>mutate(panel=factor(cluster,levels=sort(unique(cluster)),labels=paste0("Cluster ",sort(unique(cluster))," (N=",cnt$proteins[match(sort(unique(cluster)),cnt$cluster)],")")))
  centroid<-lines|>group_by(cluster,panel,time)|>
    summarise(estimate=weighted.mean(mean,pmax(N,1),na.rm=TRUE),lo=quantile(mean,.10,na.rm=TRUE),hi=quantile(mean,.90,na.rm=TRUE),N=sum(N),.groups="drop")|>
    group_by(cluster,panel)|>group_modify(~{
      q<-.x|>filter(is.finite(time),is.finite(estimate))|>arrange(time);if(nrow(q)<4)return(q)
      fit<-tryCatch(smooth.spline(q$time,q$estimate,w=sqrt(pmax(q$N,1)),spar=C1_YY_SPAR),error=function(e)NULL)
      if(!is.null(fit))q$estimate<-predict(fit,q$time)$y;q
    })|>ungroup()
  pA<-ggplot(lines,aes(time,mean,group=feature))+
    geom_vline(xintercept=0,linetype=2,color="grey45")+geom_line(color="grey55",alpha=.16,linewidth=.36)+
    geom_ribbon(data=centroid,aes(x=time,ymin=lo,ymax=hi,group=cluster,fill=factor(cluster)),inherit.aes=FALSE,alpha=.16,color=NA)+
    geom_line(data=centroid,aes(time,estimate,group=cluster,color=factor(cluster)),inherit.aes=FALSE,linewidth=1.45)+
    facet_wrap(~panel,nrow=1,scales="free_y")+scale_color_brewer(palette="Dark2",guide="none")+
    scale_fill_brewer(palette="Dark2",guide="none")+
    labs(title=paste0(panel_label,". ",ifelse(side=="incident","Incident","Baseline-prevalent"),
                      " trajectory clusters (top ",nrow(tm$X)," biomarkers; k = ",ck$k,")"),
         subtitle="Clusters are learned on the named side; both Yin and Yang data are displayed. Thick lines are weighted smoothing splines; ribbons are 10th-90th percentiles.",
         x="Years relative to baseline",y="Mean biomarker z-score")+theme_5c(10)
  # Do NOT min-max rescale the three criteria together. A constant stability
  # series can legitimately equal 1, but rescaling formerly collapsed the
  # display into horizontal 0/1 lines and hid the silhouette/gap structure.
  md<-ck$metrics|>select(k,silhouette,gap,stability)|>pivot_longer(-k,names_to="criterion",values_to="value")|>
    mutate(criterion=factor(criterion,levels=c("silhouette","gap","stability"),
      labels=c("Silhouette","Gap statistic","Bootstrap stability (ARI)")))
  bestpts<-md|>filter(k==ck$k)
  pB<-ggplot(md,aes(k,value,group=criterion))+geom_vline(xintercept=ck$k,linetype=2,color="grey40")+
    geom_line(linewidth=.85,color="#3F6F9F")+geom_point(size=2,color="#3F6F9F")+
    geom_point(data=bestpts,size=3.1,shape=21,fill="white",stroke=.9,color="#B2182B")+
    facet_wrap(~criterion,nrow=1,scales="free_y")+scale_x_continuous(breaks=unique(md$k))+
    labs(title="c. Cluster-number diagnostics",x="Number of clusters (k)",y="Native criterion value")+
    theme_5c(9)+theme(legend.position="none")
  list(plot=pA/pB+plot_layout(heights=c(3.2,1)),cluster_plot=pA,diagnostic_plot=pB,
       cluster=cl,metrics=ck$metrics,raw=tm$raw,k=ck$k)
}


# 🚩 Adjusted quintile analysis
plot_quantile_top <- function(dat, features, tvar, evar, covars, title_prefix) {
  features <- intersect(features,names(dat)); covars<-intersect(covars,names(dat))
  map(features,function(x){
    need<-unique(c(tvar,evar,x,covars)); d<-dat[,need,drop=FALSE];d<-d[complete.cases(d),,drop=FALSE]
    if(nrow(d)<500||sum(d[[evar]]==1)<20)return(blank_plot(x,"Insufficient complete-case follow-up"))
    d$value<-suppressWarnings(as.numeric(d[[x]]));d$quintile<-factor(ntile(d$value,5),levels=1:5,labels=paste0("Q",1:5))
    sf<-tryCatch(survfit(Surv(d[[tvar]],d[[evar]])~d$quintile),error=function(e)NULL)
    if(is.null(sf))return(blank_plot(x,"Survival model failed"))
    ss<-summary(sf); sdf<-tibble(time=ss$time,surv=ss$surv,lo=ss$lower,hi=ss$upper,strata=as.character(ss$strata))|>
      mutate(quintile=str_extract(strata,"Q[1-5]"))
    f<-as.formula(paste0("Surv(",bt(tvar),",",bt(evar),") ~ relevel(quintile,'Q1') + ",paste(bt(covars),collapse=" + ")))
    fit<-tryCatch(coxph(f,d,ties="efron"),error=function(e)NULL);sm<-if(is.null(fit))NULL else coef(summary(fit))
    rn<-if(is.null(sm))character() else rownames(sm);i5<-grep("Q5",rn,fixed=TRUE)[1]
    ann<-if(is.na(i5)||!length(i5))"" else sprintf("Adjusted Q5 vs Q1: HR %.2f (%.2f–%.2f), P=%s",
      exp(sm[i5,"coef"]),exp(sm[i5,"coef"]-1.96*sm[i5,"se(coef)"]),exp(sm[i5,"coef"]+1.96*sm[i5,"se(coef)"]),
      format.pval(sm[i5,"Pr(>|z|)"],digits=2,eps=1e-99))
    ggplot(sdf,aes(time,1-surv,color=quintile,fill=quintile))+geom_step(linewidth=.7)+
      geom_ribbon(aes(ymin=1-hi,ymax=1-lo),alpha=.07,color=NA)+annotate("text",x=Inf,y=Inf,label=ann,hjust=1.03,vjust=1.5,size=3,fontface="bold")+
      scale_color_brewer(palette="Spectral",direction=-1)+scale_fill_brewer(palette="Spectral",direction=-1)+
      labs(title=x,x="Years after baseline",y="Cumulative incidence",color=NULL,fill=NULL)+
      theme_5c(10)+theme(legend.position="top")
  })
}


# 🚩 Functional enrichment
# Protein enrichment prefers reproducible offline GO analysis with the assayed
# proteome as universe, then falls back to g:Profiler.
run_offline_go_enrichment <- function(sig,background){
  # A direct hypergeometric ORA avoids making clusterProfiler a hard dependency.
  # AnnotationDbi + org.Hs.eg.db are sufficient; GO.db is used only for labels.
  need<-c("AnnotationDbi","org.Hs.eg.db")
  if(!all(vapply(need,requireNamespace,logical(1),quietly=TRUE)))return(tibble())
  mp<-tryCatch(suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,keys=background,keytype="SYMBOL",columns=c("GO","ONTOLOGY"))),
    error=function(e)tibble())|>as_tibble()|>
    filter(!is.na(SYMBOL),!is.na(GO),ONTOLOGY%in%c("BP","MF","CC"))|>distinct(SYMBOL,GO,ONTOLOGY)
  if(!nrow(mp))return(tibble())
  nm<-if(requireNamespace("GO.db",quietly=TRUE))tryCatch(suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,keys=unique(mp$GO),keytype="GOID",columns="TERM")),error=function(e)tibble())|>
    as_tibble()|>transmute(GO=GOID,term_name=TERM)|>distinct(GO,.keep_all=TRUE) else
    tibble(GO=unique(mp$GO),term_name=unique(mp$GO))
  map_dfr(c("BP","MF","CC"),function(ont){
    mo<-mp|>filter(ONTOLOGY==ont);if(!nrow(mo))return(tibble())
    universe<-intersect(background,unique(mo$SYMBOL));query<-intersect(sig,universe)
    if(length(universe)<20||length(query)<2)return(tibble())
    mo|>filter(SYMBOL%in%universe)|>group_by(GO)|>
      summarise(term_size=n_distinct(SYMBOL),intersection_size=n_distinct(SYMBOL[SYMBOL%in%query]),.groups="drop")|>
      filter(term_size>=5,intersection_size>=1)|>
      mutate(p_raw=phyper(intersection_size-1,term_size,length(universe)-term_size,length(query),lower.tail=FALSE),
             adjusted_p=p.adjust(p_raw,"BH"),source=paste0("GO:",ont))|>
      left_join(nm,by="GO")|>transmute(source,term_id=GO,term_name=coalesce(term_name,GO),adjusted_p,
        intersection_size=as.integer(intersection_size),term_size=as.integer(term_size))
  })
}

run_functional_enrichment <- function(assoc, layer, background=assoc$term) {
  sig<-assoc|>filter(is.finite(p.value),p.value*nrow(assoc)<.05)|>pull(term)|>unique()
  background<-intersect(unique(as.character(background)),unique(as.character(assoc$term)))
  empty<-tibble(source=character(),term_id=character(),term_name=character(),adjusted_p=numeric(),
    intersection_size=integer(),term_size=integer(),query_size=integer(),background_n=integer(),enrichment_ratio=numeric())
  if(!length(sig))return(empty)
  # Metabolites do not have a GO-style gene ontology.  Use the curated UKB
  # met.lst hierarchy as the assay-specific universe and perform the same
  # over-representation test at super-group, group and subgroup levels.
  if(layer=="metabolite"){
    ann<-read_met_annotation(background)|>filter(trait%in%background)
    if(!nrow(ann))return(empty)
    one_level<-function(col,source_name){
      z<-ann|>transmute(trait,term=coalesce(na_if(as.character(.data[[col]]),""),"Other"))|>distinct()
      map_dfr(sort(unique(z$term)),function(tt){
        members<-z|>filter(term==tt)|>pull(trait)|>unique();m<-length(intersect(members,background));k<-length(intersect(members,sig))
        if(m<2||k<1)return(tibble())
        tibble(source=source_name,term_id=paste(source_name,tt,sep=":"),term_name=tt,
          p_raw=phyper(k-1,m,length(background)-m,length(sig),lower.tail=FALSE),intersection_size=k,term_size=m)
      })
    }
    ans<-bind_rows(one_level("super_group","Metabolite super-group"),one_level("group","Metabolite group"),
                   one_level("subgroup","Metabolite subgroup"))
    if(!nrow(ans))return(empty)
    return(ans|>group_by(source)|>mutate(adjusted_p=p.adjust(p_raw,"BH"))|>ungroup()|>
      mutate(query_size=length(sig),background_n=length(background),
        enrichment_ratio=(intersection_size/query_size)/(term_size/background_n))|>
      select(source,term_id,term_name,adjusted_p,intersection_size,term_size,query_size,background_n,enrichment_ratio)|>arrange(adjusted_p))
  }
  if(layer!="protein")return(empty)
  ans<-run_offline_go_enrichment(sig,background);used_custom_bg<-nrow(ans)>0
  if(!nrow(ans)&&requireNamespace("gprofiler2",quietly=TRUE)){
    gp<-tryCatch(gprofiler2::gost(sig,organism="hsapiens",correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),custom_bg=background,domain_scope="custom",user_threshold=1),error=function(e)NULL)
    if(is.null(gp)||is.null(gp$result)||!nrow(gp$result)){
      used_custom_bg<-FALSE
      gp<-tryCatch(gprofiler2::gost(sig,organism="hsapiens",correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),domain_scope="annotated",user_threshold=1),error=function(e)NULL)
    }
    if(!is.null(gp)&&!is.null(gp$result)&&nrow(gp$result))ans<-as_tibble(gp$result)|>
      transmute(source,term_id=native,term_name=name,adjusted_p=p_value,
        intersection_size=as.integer(intersection_size),term_size=as.integer(term_size),query_size=as.integer(query_size))
  }
  if(!nrow(ans)){
    py<-Sys.which("python3");script<-file.path(Sys.getenv("LE8_FDIR"),"c1_enrich.py")
    if(nzchar(py)&&file.exists(script)){
      sf<-tempfile(fileext=".txt");bf<-tempfile(fileext=".txt");of<-tempfile(fileext=".csv");on.exit(unlink(c(sf,bf,of)),add=TRUE)
      writeLines(sig,sf);writeLines(background,bf);status<-suppressWarnings(system2(py,c(script,sf,bf,of),stdout=FALSE,stderr=FALSE))
      if(status==0&&file.exists(of)){ans<-as_tibble(data.table::fread(of,showProgress=FALSE));used_custom_bg<-TRUE}
    }
  }
  # A network-free, non-GO fallback prevents a blank figure when Bioconductor
  # annotation is unavailable.  It is explicitly labelled as protein assay-group
  # enrichment and therefore is not confused with pathway enrichment.
  if(!nrow(ans)){
    ann<-tryCatch(read_prot_bed(background)|>filter(protein%in%background),error=function(e)tibble())
    if(nrow(ann)){
      ans<-map_dfr(sort(unique(ann$group)),function(g){
        members<-ann|>filter(group==g)|>pull(protein)|>unique();m<-length(members);k<-length(intersect(members,sig))
        if(m<2||k<1)return(tibble())
        tibble(source="protein assay group",term_id=paste0("PROT:",g),term_name=g,
          adjusted_p=phyper(k-1,m,length(background)-m,length(sig),lower.tail=FALSE),
          intersection_size=k,term_size=m,query_size=length(sig),background_n=length(background))
      })
      if(nrow(ans))ans<-ans|>mutate(adjusted_p=p.adjust(adjusted_p,"BH"));used_custom_bg<-TRUE
    }
  }
  if(!nrow(ans))return(empty)
  for(nm in c("term_size","query_size"))if(!nm%in%names(ans))ans[[nm]]<-NA_integer_
  ans|>mutate(query_size=coalesce(as.integer(query_size),length(sig)),
    background_n=if(used_custom_bg)length(background) else NA_integer_,
    enrichment_ratio=ifelse(is.finite(term_size)&term_size>0&is.finite(background_n)&background_n>0,
      (intersection_size/query_size)/(term_size/background_n),NA_real_))|>
    select(source,term_id,term_name,adjusted_p,intersection_size,term_size,query_size,background_n,enrichment_ratio)|>arrange(adjusted_p)
}

plot_functional_enrichment <- function(enrich,n_sig,assoc=NULL,layer="protein",panel_prefix=NULL) {
  add_prefix <- function(x) if(is.null(panel_prefix)||!nzchar(panel_prefix)) x else paste0(panel_prefix,"\n",x)
  if(!nrow(enrich)){
    return(blank_plot(add_prefix(paste0("Functional enrichment (",n_sig," significant biomarkers)")),
      if(layer=="protein")"No GO/pathway result was available; biomarker P values are not shown as pathway enrichment"
      else "No metabolite hierarchy contained an eligible enriched class"))
  }
  levels0<-c("GO:BP","GO:MF","GO:CC","KEGG","REAC","protein assay group","Metabolite super-group","Metabolite group","Metabolite subgroup")
  labs0<-c(`GO:BP`="BP",`GO:MF`="MF",`GO:CC`="CC",KEGG="KEGG",REAC="Reactome",`protein assay group`="protein group",
    `Metabolite super-group`="Super-group",`Metabolite group`="Group",`Metabolite subgroup`="Subgroup")
  is_met_enrichment<-any(str_detect(enrich$source,"^Metabolite"))
  levels0<-levels0[levels0%in%unique(enrich$source)]
  d<-enrich|>filter(is.finite(adjusted_p))|>group_by(source)|>slice_min(adjusted_p,n=6,with_ties=FALSE)|>ungroup()|>
    mutate(source=factor(source,levels=levels0,labels=unname(labs0[levels0])),score=-log10(pmax(adjusted_p,.Machine$double.xmin)),
      term_label=str_wrap(term_name,34),term_key=paste(source,term_id,sep="___"))|>
    arrange(source,adjusted_p)|>mutate(term_key=factor(term_key,levels=rev(unique(term_key))))
  if(!nrow(d))return(blank_plot(add_prefix("Functional enrichment"),"No finite enrichment result was available"))
  use_ratio<-all(is.finite(d$enrichment_ratio))
  d$x_value<-if(use_ratio)d$enrichment_ratio else d$score
  ggplot(d,aes(x_value,term_key,color=score,size=intersection_size))+
    {if(use_ratio)geom_vline(xintercept=1,linetype=3,color="grey70")else geom_blank()}+
    geom_point(alpha=.9)+facet_grid(source~.,scales="free_y",space="free_y")+
    scale_y_discrete(labels=setNames(d$term_label,as.character(d$term_key)))+
    scale_color_viridis_c(option="C",direction=-1,name=expression(-log[10](FDR)))+
    scale_size_continuous(range=c(2.2,7),name="Overlapping biomarkers")+
    labs(title=add_prefix(paste0("Functional enrichment of ",n_sig," Bonferroni-significant biomarkers")),
      subtitle=paste0(if(is_met_enrichment)"Assay-specific metabolite hierarchy" else if(any(d$source=="protein assay group"))"protein assay-group fallback (not GO/pathway evidence)" else if(any(is.na(d$background_n)))"Annotated-genome background" else "All assayed biomarkers as background",
        "; ",sum(d$adjusted_p<.05)," displayed terms have FDR < 0.05"),
      x=if(use_ratio)"Fold enrichment"else expression(-log[10](FDR)),y=NULL)+
    theme_5c(9)+theme(legend.position="right",panel.grid.major.y=element_line(color="grey93",linewidth=.25),
      strip.text.y=element_text(angle=0),plot.title.position="plot",plot.margin=margin(7,13,7,13))
}

plot_sameN_attenuation <- function(att, assoc_adj2, top_n=30L) {
  all_d <- att |> filter(is.finite(beta_basic_sameN),is.finite(beta_adj2_sameN)) |>
    left_join(assoc_adj2 |> select(term,p_adj2=p.value),by="term") |>
    mutate(direction=ifelse(beta_adj2_sameN>=0,"Positive","Inverse"),
      log2_ratio_plot=pmax(-2,pmin(2,log2_effect_ratio)))
  if(!nrow(all_d)) return(blank_plot("Same-sample covariate attenuation","No same-sample comparison was estimable"))
  pA<-ggplot(all_d,aes(beta_basic_sameN,beta_adj2_sameN,color=log2_ratio_plot,shape=sign_flip))+
    geom_abline(slope=1,intercept=0,linetype=2,color="grey50")+geom_hline(yintercept=0,color="grey88")+geom_vline(xintercept=0,color="grey88")+
    geom_point(alpha=.68,size=1.8)+
    scale_color_gradient2(low="#2C7FB8",mid="grey78",high="#D7301F",midpoint=0,
      breaks=c(-2,-1,0,1,2),labels=c("0.25x","0.5x","1x","2x","4x"),name="Adj2/basic\n|beta| ratio")+
    scale_shape_manual(values=c(`FALSE`=16,`TRUE`=4),na.value=1,name="Direction flip")+
    labs(title="a. Same-participant basic versus adj2 effects",
      subtitle="Colour is omitted when |basic beta| < 0.02, where percentage attenuation is unstable",
      x="Basic Cox beta",y="Adj2 Cox beta")+theme_5c(9)
  d<-all_d|>arrange(p_adj2)|>slice_head(n=min(top_n,24L))|>mutate(term=factor(term,levels=rev(term)),
    attenuation_label=case_when(is.na(effect_ratio)~"basic beta near 0",
      sign_flip~sprintf("%.2fx; flip",effect_ratio),TRUE~sprintf("%.2fx",effect_ratio)))
  xr<-range(c(d$beta_basic_sameN,d$beta_adj2_sameN),na.rm=TRUE);span<-diff(xr);if(!is.finite(span)||span==0)span<-1
  ann_x<-xr[2]+.10*span
  pB<-ggplot(d,aes(y=term))+
    geom_vline(xintercept=0,color="grey82")+
    geom_segment(aes(x=beta_basic_sameN,xend=beta_adj2_sameN,yend=term),color="grey72",linewidth=.8)+
    geom_point(aes(x=beta_basic_sameN,color="Basic"),size=2.5)+
    geom_point(aes(x=beta_adj2_sameN,color="Adj2"),size=2.8)+
    geom_text(aes(x=ann_x,label=attenuation_label),hjust=0,size=3,fontface="bold",color="grey35")+
    scale_color_manual(values=c(Basic="grey55",Adj2="#3F78A8"),name=NULL)+
    scale_x_continuous(limits=c(xr[1]-.04*span,xr[2]+.25*span))+
    labs(title="b. Strongest adj2 associations",
      subtitle="Basic and adj2 Cox models use the same complete-case sample; labels show the stable |adj2 beta| / |basic beta| ratio",
      x="Log hazard ratio per 1-SD biomarker",y=NULL)+theme_5c(9)+theme(legend.position="top")
  pA|pB
}

c1_panel_n <- function(x) {
  n <- suppressWarnings(as.integer(x$N_total))
  n <- sort(unique(n[is.finite(n)]))
  if (!length(n)) return("NA")
  if (length(n) == 1L) return(format(n, big.mark = ",", scientific = FALSE, trim = TRUE))
  paste(format(range(n), big.mark = ",", scientific = FALSE, trim = TRUE), collapse = "–")
}

plot_c1_fig12 <- function(layer, features_all, pgs_incident, pgs_prevalent,
                          pgs_attained_age, assoc_adj2, prevalent_adj2,
                          birthline_adj2, outdir) {
  omic_name <- if (layer == "protein") "protein" else "metabolite"
  pgs_title <- function(panel, analysis, x, include_outcome = FALSE) {
    paste0(panel, ". ", analysis, if (include_outcome) paste0(" ", Y) else "",
      " — inherited ", omic_name, " PGS, N = ", c1_panel_n(x))
  }
  circle_pgs_title <- function(panel, analysis, x) {
    str_replace(pgs_title(panel,analysis,x,TRUE),", N = ",",\nN = ")
  }

  # The inherited-score panels use everyone with a matched PGS and phenotype
  # record. They must not be restricted to participants with the corresponding
  # adult omics assay; that restriction discards most of the genetic cohort.
  if (layer == "protein") {
    p1a<-plot_pwas_manhattan(pgs_incident,features_all,pgs_title("a","Incident",pgs_incident,TRUE))
    p1b<-plot_pwas_manhattan(pgs_prevalent,features_all,pgs_title("b","Baseline-prevalent",pgs_prevalent,TRUE))
    p1c<-plot_pwas_manhattan(assoc_adj2,features_all,paste0("c. Incident ",Y," — adj2 Cox"))
    p1d<-plot_pwas_manhattan(prevalent_adj2,features_all,paste0("d. Baseline-prevalent ",Y," — adj2 logistic"))
    p1e<-plot_pwas_manhattan(pgs_attained_age,features_all,pgs_title("e","Attained-age",pgs_attained_age,TRUE))
    p1f<-plot_pwas_manhattan(birthline_adj2,features_all,paste0("f. Attained-age ",Y," — measured adult level, adj2"))
    save_plot((p1a|p1c)/plot_spacer()/(p1b|p1d)/plot_spacer()/(p1e|p1f)+
      plot_layout(heights=c(1,.06,1,.06,1)),"c1.Fig1.mh.png",24,19,outdir=outdir)
    p2a<-plot_volcano(pgs_incident,"beta","term",pgs_title("a","Incident",pgs_incident),.05/max(1,nrow(pgs_incident)),TOP_N)
    p2b<-plot_volcano(pgs_prevalent,"beta","term",pgs_title("b","Prevalent",pgs_prevalent),.05/max(1,nrow(pgs_prevalent)),TOP_N)
    p2c<-plot_volcano(assoc_adj2,"beta","term","c. Incident ProtWAS — adj2",.05/max(1,nrow(assoc_adj2)),TOP_N)
    p2d<-plot_volcano(prevalent_adj2,"beta","term","d. Prevalent ProtWAS — adj2 logistic",.05/max(1,nrow(prevalent_adj2)),TOP_N)
    p2e<-plot_volcano(pgs_attained_age,"beta","term",pgs_title("e","Attained-age",pgs_attained_age),.05/max(1,nrow(pgs_attained_age)),TOP_N)
    p2f<-plot_volcano(birthline_adj2,"beta","term","f. Attained-age ProtWAS — measured adult level, adj2",.05/max(1,nrow(birthline_adj2)),TOP_N)
    save_plot((p2a|p2c)/plot_spacer()/(p2b|p2d)/plot_spacer()/(p2e|p2f)+
      plot_layout(heights=c(1,.06,1,.06,1)),"c1.Fig2.vc.png",19,19,outdir=outdir)
  } else {
    ann<-read_met_annotation(features_all);annot<-function(z)z|>mutate(trait=term)|>left_join(ann,by="trait")|>mutate(label=coalesce(label,term),group=coalesce(group,"Other"))
    circ_data<-list(annot(pgs_incident),annot(pgs_prevalent),annot(assoc_adj2),annot(prevalent_adj2),annot(pgs_attained_age),annot(birthline_adj2))
    shared_beta_lim<-as.numeric(quantile(abs(unlist(map(circ_data,~.x$beta))),.98,na.rm=TRUE,names=FALSE))
    if(!is.finite(shared_beta_lim)||shared_beta_lim<=0)shared_beta_lim<-max(abs(unlist(map(circ_data,~.x$beta))),na.rm=TRUE)
    p1a<-plot_met_circle(circ_data[[1]],circle_pgs_title("a","Incident",pgs_incident),label_n=14,beta_limit=shared_beta_lim)
    p1b<-plot_met_circle(circ_data[[2]],circle_pgs_title("b","Prevalent",pgs_prevalent),label_n=14,beta_limit=shared_beta_lim)
    p1c<-plot_met_circle(circ_data[[3]],paste0("c. Incident ",Y," MWAS — adj2"),label_n=14,beta_limit=shared_beta_lim)
    p1d<-plot_met_circle(circ_data[[4]],paste0("d. Prevalent ",Y," MWAS — adj2 logistic"),label_n=14,beta_limit=shared_beta_lim)
    p1e<-plot_met_circle(circ_data[[5]],circle_pgs_title("e","Attained-age",pgs_attained_age),label_n=14,beta_limit=shared_beta_lim)
    p1f<-plot_met_circle(circ_data[[6]],paste0("f. Attained-age ",Y," — measured adult level, adj2"),label_n=14,beta_limit=shared_beta_lim)
    # Paper_t2dm Fig2 style: large radial panels and a dedicated shared-key row.
    # A guide_area prevents four independent legends from shrinking the circles.
    p_circle<-(p1a|p1c)/plot_spacer()/(p1b|p1d)/plot_spacer()/(p1e|p1f)/guide_area()+
      plot_layout(heights=c(1,.06,1,.06,1,.10),guides="collect") &
      theme(legend.position="bottom",legend.direction="horizontal")
    save_plot(p_circle,"c1.Fig1.circular.png",24,20,outdir=outdir)
    p2a<-plot_volcano(pgs_incident,"beta","term",pgs_title("a","Incident",pgs_incident),.05/max(1,nrow(pgs_incident)),12,x_quantile=.985)
    p2b<-plot_volcano(pgs_prevalent,"beta","term",pgs_title("b","Prevalent",pgs_prevalent),.05/max(1,nrow(pgs_prevalent)),12,x_quantile=.985)
    p2c<-plot_volcano(assoc_adj2,"beta","term","c. Incident MWAS — adj2",.05/max(1,nrow(assoc_adj2)),12,x_quantile=.985)
    p2d<-plot_volcano(prevalent_adj2,"beta","term","d. Prevalent MWAS — adj2 logistic",.05/max(1,nrow(prevalent_adj2)),12,x_quantile=.985)
    p2e<-plot_volcano(pgs_attained_age,"beta","term",pgs_title("e","Attained-age",pgs_attained_age),.05/max(1,nrow(pgs_attained_age)),12,x_quantile=.985)
    p2f<-plot_volcano(birthline_adj2,"beta","term","f. Attained-age MWAS — measured adult level, adj2",.05/max(1,nrow(birthline_adj2)),12,x_quantile=.985)
    save_plot((p2a|p2c)/plot_spacer()/(p2b|p2d)/plot_spacer()/(p2e|p2f)+
      plot_layout(heights=c(1,.06,1,.06,1)),"c1.Fig2.vc.png",19,19,outdir=outdir)
  }
}


# 🚩 Main C1
run_c1_layer <- function(layer=c("protein","metabolite")) {
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir)
  rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  scan_cache<-file.path(rawdir,paste0("c1.scan.",C1_SCAN_VERSION,".rds"));selected_cache<-file.path(rawdir,"c1.res.rds")
  if (truthy(Sys.getenv("C1_FIG12_ONLY", unset = "FALSE"))) {
    if (!file.exists(selected_cache) || file.size(selected_cache) <= 0)
      stop("C1_FIG12_ONLY requires an existing C1 result: ", selected_cache, call. = FALSE)
    old<-readRDS(selected_cache)
    needed<-c("pgs_incident","pgs_prevalent","pgs_attained_age","association_adj2","prevalent_adj2","birthline_adj2")
    missing<-setdiff(needed,names(old));if(length(missing))
      stop("Existing C1 result lacks fields required for Fig1/2: ",paste(missing,collapse=", "),call.=FALSE)
    features_all<-old$input_feature_annotation_audit$feature %||%
      unique(c(old$pgs_incident$term,old$association_adj2$term))
    plot_c1_fig12(layer,features_all,old$pgs_incident,old$pgs_prevalent,
      old$pgs_attained_age,old$association_adj2,old$prevalent_adj2,
      old$birthline_adj2,outdir)
    message("Completed C1 Fig1/2-only redraw: ",outdir)
    return(invisible(old))
  }
  pgs_signature<-c1_pgs_signature(layer)
  if(cache_valid(selected_cache)){
    old<-tryCatch(readRDS(selected_cache),error=function(e)NULL)
    if(!is.null(old)&&identical(old$meta$code_version%||%NA_character_,C1_CODE_VERSION)&&
       identical(old$meta$pgs_signature%||%NA_character_,pgs_signature)&&
       all(c("pgs_incident","pgs_prevalent","pgs_attained_age")%in%names(old))){
      cache_message(paste0("C1/",layer),selected_cache);return(old)
    }
    message("C1/",layer,": cache or PGS input predates ",C1_CODE_VERSION,
      "; recomputing inherited-score and observed-level outputs")
  }

  biom0<-if(layer=="protein")read_prot() else read_met();features_all<-setdiff(names(biom0),"eid")
  scan<-if(cache_valid(scan_cache))tryCatch(readRDS(scan_cache),error=function(e)NULL) else NULL
  reuse<-!is.null(scan)

  # If the expensive scans exist, load only the features needed for figures + same-N table.
  if(reuse){
    a0<-scan$association_adj2 %||% scan$association_basic; p0<-scan$prevalent_adj2 %||% scan$prevalent_basic
    fig_features<-unique(c(
      C1_DIRECTION_ANCHORS,
      a0|>filter(is.finite(p.value))|>slice_min(p.value,n=max(CLUSTER_TOP,ATTENUATION_TOP,TOP_N),with_ties=FALSE)|>pull(term),
      p0|>filter(is.finite(p.value))|>slice_min(p.value,n=max(GRADIENT_TOP,TOP_N),with_ties=FALSE)|>pull(term)))
    biom<-biom0[,intersect(c("eid",fig_features),names(biom0)),drop=FALSE]
  } else biom<-biom0

  need<-unique(c("eid","ethnic.c",vars.basic,vars.le8,"birth_date","date_attend","date_lost","date_death",paste0("fod_icd10_",Y)))
  # Restrict C1 to participants with the omics assay.
  dat<-read_all(need)|>inner_join(biom,by="eid")|>filter_analysis_cohort()|>make_outcome(Y)|>add_attained_age_time(Y)
  rm(biom,biom0);invisible(gc())
  features<-intersect(features_all,names(dat));covs_basic<-intersect(vars.basic,names(dat));covs_adj2<-intersect(vars.adj2,names(dat))
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");bvar<-paste0(Y,".b2e");bivar<-paste0(Y,".bi2e")
  rvar<-paste0(Y,".r2e");revar<-paste0(Y,".Yr2e")
  # Attained age is already the analysis time, so an additional linear age
  # covariate would double-adjust age.  Sex and the remaining covariates stay.
  age_covs<-unique(c("age",grep("^age($|[._])",unique(c(covs_basic,covs_adj2)),value=TRUE,ignore.case=TRUE)))
  birth_covs_basic<-setdiff(covs_basic,age_covs);birth_covs_adj2<-setdiff(covs_adj2,age_covs)
  dat$prevalent_status<-make_prevalent_status(dat,Y)

  if(!reuse){
    message("C1/",layer,": incident Cox scans: basic + adj2")
    assoc_basic<-cox_scan(dat,features,covs_basic,Y,time_var=tvar,event_var=evar)
    assoc_adj2<-cox_scan(dat,features,covs_adj2,Y,time_var=tvar,event_var=evar)
    message("C1/",layer,": attained-age Cox scans with delayed entry at baseline")
    birthline_basic<-cox_scan_delayed_entry(dat,features,birth_covs_basic,Y)
    birthline_adj2<-cox_scan_delayed_entry(dat,features,birth_covs_adj2,Y)
    message("C1/",layer,": baseline-prevalent logistic scans: basic + adj2")
    prevalent_basic<-logistic_scan(dat,features,covs_basic,"prevalent_status")
    prevalent_adj2<-logistic_scan(dat,features,covs_adj2,"prevalent_status")
    message("C1/",layer,": case-only diagnosis-duration and landmark incident scans")
    # Retained at the user's request as an exploratory legacy analysis.  Its
    # cases and controls have incompatible time origins, so it is never used in
    # directionality classification or C5 evidence grading.
    reverse_adj2<-if(all(c(rvar,revar)%in%names(dat)))
      cox_scan(dat,features,covs_adj2,Y,time_var=rvar,event_var=revar)|>
        mutate(analysis_status="exploratory only: incompatible time origins") else
      tibble(term=features,beta=NA_real_,p.value=NA_real_,FDR=NA_real_,analysis_status="not estimable")
    duration_adj2<-prevalent_duration_scan(dat,features,bvar,covs_adj2)
    landmark_features<-unique(c(C1_DIRECTION_ANCHORS,
      assoc_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=C1_LANDMARK_TOP,with_ties=FALSE)|>pull(term),
      prevalent_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=C1_LANDMARK_TOP,with_ties=FALSE)|>pull(term)))
    landmark_adj2<-landmark_incident_scan(dat,landmark_features,tvar,evar,covs_adj2,C1_LANDMARK_YEARS)
    risk_features<-unique(c(C1_DIRECTION_ANCHORS,
      assoc_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=C1_RISK_TOP,with_ties=FALSE)|>pull(term),
      prevalent_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=C1_RISK_TOP,with_ties=FALSE)|>pull(term)))
    risk_window_adj2<-risk_window_scan(dat,risk_features,tvar,evar,bvar,covs_adj2,C1_RISK_CUTS)
    # Same-N attenuation only needs the strongest features; set C1_ATTENUATION_TOP=0 for all.
    att_features<-assoc_basic|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)
    if(ATTENUATION_TOP>0)att_features<-head(att_features,ATTENUATION_TOP)
    attenuation_sameN<-same_sample_attenuation(dat,att_features,tvar,evar,covs_basic,covs_adj2)
    scan<-list(association_basic=assoc_basic,association_adj2=assoc_adj2,
               birthline_basic=birthline_basic,birthline_adj2=birthline_adj2,
               prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,reverse_adj2=reverse_adj2,
               duration_adj2=duration_adj2,landmark_adj2=landmark_adj2,risk_window_adj2=risk_window_adj2,
               attenuation_sameN=attenuation_sameN)
    saveRDS(scan,scan_cache,compress="xz")
  } else {
    assoc_basic<-scan$association_basic;assoc_adj2<-scan$association_adj2
    birthline_basic<-scan$birthline_basic%||%cox_scan_delayed_entry(dat,features,birth_covs_basic,Y)
    birthline_adj2<-scan$birthline_adj2%||%cox_scan_delayed_entry(dat,features,birth_covs_adj2,Y)
    prevalent_basic<-scan$prevalent_basic;prevalent_adj2<-scan$prevalent_adj2
    reverse_adj2<-scan$reverse_adj2%||%tibble(term=features,beta=NA_real_,p.value=NA_real_,FDR=NA_real_)
    duration_adj2<-scan$duration_adj2%||%tibble(term=features,beta=NA_real_,p.value=NA_real_,FDR=NA_real_)
    landmark_adj2<-scan$landmark_adj2%||%tibble()
    risk_features<-unique(c(C1_DIRECTION_ANCHORS,
      assoc_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=C1_RISK_TOP,with_ties=FALSE)|>pull(term),
      prevalent_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=C1_RISK_TOP,with_ties=FALSE)|>pull(term)))
    risk_window_adj2<-scan$risk_window_adj2%||%risk_window_scan(dat,risk_features,tvar,evar,bvar,covs_adj2,C1_RISK_CUTS)
    attenuation_sameN<-scan$attenuation_sameN %||% tibble()
    message("C1/",layer,": reuse association scans and regenerate figures")
  }

  assoc<-if(covs_use_name=="adj2")assoc_adj2 else assoc_basic
  assoc_prevalent<-if(covs_use_name=="adj2")prevalent_adj2 else prevalent_basic
  pgs<-run_c1_pgs_scan(layer,features_all,covs_basic,Y,rawdir,overlap_eids=dat$eid)
  pgs_incident<-pgs$incident;pgs_prevalent<-pgs$prevalent;pgs_attained_age<-pgs$attained_age
  pgs_incident_same<-pgs$incident_same_omic;pgs_prevalent_same<-pgs$prevalent_same_omic
  pgs_attained_age_same<-pgs$attained_age_same_omic
  pgs_concordance<-build_pgs_actual_concordance(
    list(incident=pgs_incident_same,prevalent=pgs_prevalent_same,attained_age=pgs_attained_age_same),
    list(incident=assoc_adj2,prevalent=prevalent_adj2,attained_age=birthline_adj2))
  prefix<-ifelse(layer=="protein","pwas","mwas")
  input_feature_audit<-layer_annotation_audit(layer,features_all)
  write_raw_csv(input_feature_audit,"c1.input_feature_annotation_audit.csv",rawdir)
  write_raw_csv(assoc_basic,paste0(prefix,"_incident_basic.csv"),rawdir);write_raw_csv(assoc_adj2,paste0(prefix,"_incident_adj2.csv"),rawdir)
  write_raw_csv(birthline_basic,paste0(prefix,"_birthline_attained_age_basic.csv"),rawdir)
  write_raw_csv(birthline_adj2,paste0(prefix,"_birthline_attained_age_adj2.csv"),rawdir)
  write_raw_csv(prevalent_basic,paste0(prefix,"_prevalent_basic.csv"),rawdir);write_raw_csv(prevalent_adj2,paste0(prefix,"_prevalent_adj2.csv"),rawdir)
  write_raw_csv(reverse_adj2,paste0(prefix,"_prevalent_reverse_cox_adj2.csv"),rawdir)
  write_raw_csv(duration_adj2,paste0(prefix,"_prevalent_duration_adj2.csv"),rawdir)
  write_raw_csv(landmark_adj2,paste0(prefix,"_incident_landmark_adj2.csv"),rawdir)
  write_raw_csv(risk_window_adj2,paste0(prefix,"_diagnosis_window_riskset_adj2.csv"),rawdir)
  write_raw_csv(attenuation_sameN,"adj2_attenuation_sameN.csv",rawdir)
  write_raw_csv(pgs$status,"c1.pgs_status.csv",rawdir)
  write_raw_csv(pgs_incident,paste0(prefix,"_pgs_incident_full_genetic.csv"),rawdir)
  write_raw_csv(pgs_prevalent,paste0(prefix,"_pgs_prevalent_full_genetic.csv"),rawdir)
  write_raw_csv(pgs_attained_age,paste0(prefix,"_pgs_attained_age_full_genetic.csv"),rawdir)
  write_raw_csv(pgs_incident_same,paste0(prefix,"_pgs_incident_same_omic.csv"),rawdir)
  write_raw_csv(pgs_prevalent_same,paste0(prefix,"_pgs_prevalent_same_omic.csv"),rawdir)
  write_raw_csv(pgs_attained_age_same,paste0(prefix,"_pgs_attained_age_same_omic.csv"),rawdir)
  write_raw_csv(pgs_concordance,"c1.pgs_actual_concordance.csv",rawdir)

  cohort<-tibble(layer,N_omics=nrow(dat),incident_events=sum(dat[[evar]]==1,na.rm=TRUE),
                 prevalent_cases=sum(dat$prevalent_status==1,na.rm=TRUE),features=length(features_all),
                 annotation_matched=sum(input_feature_audit$annotation_matched,na.rm=TRUE),
                 annotation_unmatched=sum(!input_feature_audit$annotation_matched,na.rm=TRUE),
                 attained_age_N=sum(is.finite(dat$.attained_entry)&is.finite(dat$.attained_exit)),
                 bi2e_available=sum(is.finite(dat[[bivar]])),
                 PGS_matched=length(pgs$score_map),PGS_file_signature=pgs_signature,
                 LE8_components=length(intersect(vars.le8,names(dat))),covs_use=covs_use_name)
  write_raw_csv(cohort,"c1.cohort.csv",rawdir)

  top<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=TOP_N,with_ties=FALSE)|>pull(term)
  top6<-head(top,YY_TOP);top_inc<-assoc_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=GRADIENT_TOP,with_ties=FALSE)|>pull(term)
  top_prev<-prevalent_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=GRADIENT_TOP,with_ties=FALSE)|>pull(term)
  cluster_features<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=CLUSTER_TOP,with_ties=FALSE)|>pull(term)

  # Fig1 + Fig2: the third row is an attained-age sensitivity analysis with
  # delayed entry at baseline. The left column is a fixed-at-conception omic
  # PGS in the full genetic cohort; the right column is the measured adult omic
  # level with adj2 adjustment. Basic observed-level analyses remain in the
  # workbook and raw files but are no longer displayed here.
  plot_c1_fig12(layer,features_all,pgs_incident,pgs_prevalent,pgs_attained_age,
    assoc_adj2,prevalent_adj2,birthline_adj2,outdir)

  # Fig3: same participants in basic and adj2 residualized YY panels.
  yy_need<-unique(c(bvar,top6,covs_adj2));yy_dat<-dat[complete.cases(dat[,intersect(yy_need,names(dat)),drop=FALSE]) & is.finite(dat[[bvar]]),,drop=FALSE]
  omic_name<-ifelse(layer=="protein","protein","metabolite")
  p3a<-make_yy_panel(yy_dat,top6,bvar,covs_basic,paste0("a. Yin–Yang ",omic_name," patterns — basic adjustment"),
                     ifelse(layer=="protein","Mean protein z-score","Mean metabolite z-score"),smooth_lines=FALSE)
  p3b<-make_yy_panel(yy_dat,top6,bvar,covs_adj2,"b. Yin–Yang trajectories — basic + LE8 adjustment",
                     ifelse(layer=="protein","Mean protein z-score","Mean metabolite z-score"),smooth_lines=TRUE)
  save_plot(p3a/p3b+plot_layout(heights=c(1,1)),"c1.Fig3.yy_top.png",16,13.5,outdir=outdir)

  # Fig4 shows adjusted diagnosis-timed cross-sectional omics gradients.
  p4i<-make_gradient_panel(dat,top_inc,bvar,"incident",GRADIENT_STEP,MIN_BIN_N,
    paste0("a. Incident (Yin): top ",length(top_inc)," adj2 biomarkers"),covars=covs_adj2)
  p4p<-make_gradient_panel(dat,top_prev,bvar,"prevalent",GRADIENT_STEP,MIN_BIN_N,
    paste0("b. Baseline-prevalent (Yang): top ",length(top_prev)," adj2 biomarkers"),covars=covs_adj2)
  save_plot((p4i/p4p+plot_layout(guides="collect"))&theme(legend.position="right"),
            "c1.Fig4.gradient.png",17,9.6,outdir=outdir)
  gradient_rank<-bind_rows(
    assoc_adj2|>filter(term%in%top_inc)|>transmute(analysis="incident_adj2",feature=term,p_value=p.value)|>arrange(p_value)|>mutate(rank=row_number()),
    prevalent_adj2|>filter(term%in%top_prev)|>transmute(analysis="prevalent_adj2",feature=term,p_value=p.value)|>arrange(p_value)|>mutate(rank=row_number()))
  write_raw_csv(gradient_rank,"c1.gradient_top10_provenance.csv",rawdir)

  # Fig5: separately selected and clustered incident and prevalent adj2 sets.
  cluster_prev<-prevalent_adj2|>filter(is.finite(p.value))|>slice_min(p.value,n=CLUSTER_TOP,with_ties=FALSE)|>pull(term)
  cli<-make_cluster_figure(dat,cluster_features,bvar,CLUSTER_STEP,YY_MAX_YEAR,"incident","a",fixed_k=NULL)
  clp<-make_cluster_figure(dat,cluster_prev,bvar,CLUSTER_STEP,YY_MAX_YEAR,"prevalent","b",fixed_k=NULL)
  # Keep the five top-level patchwork rows explicit.  Nesting the three-row
  # cluster_main object here is flattened by patchwork, so a three-entry outer
  # heights vector leaves room for only three of the resulting five panels.
  cluster_figure<-cli$cluster_plot/plot_spacer()/clp$cluster_plot/plot_spacer()/
    (cli$diagnostic_plot|clp$diagnostic_plot)+
    plot_layout(heights=c(1,.10,1,.05,.38))
  save_plot(cluster_figure,
            "c1.Fig5.gradient_cluster.png",18,10.5,outdir=outdir)
  cl_members<-bind_rows(cli$cluster|>mutate(analysis="incident_adj2"),clp$cluster|>mutate(analysis="prevalent_adj2"))
  cl_metrics<-bind_rows(cli$metrics|>mutate(analysis="incident_adj2"),clp$metrics|>mutate(analysis="prevalent_adj2"))
  write_raw_csv(cl_members,"c1.cluster_membership.csv",rawdir);write_raw_csv(cl_metrics,"c1.cluster_selection.csv",rawdir)

  # Fig6: adjusted Q5-vs-Q1 HRs on the same adj2 complete-case sample; no subtitle.
  qplots<-plot_quantile_top(dat,top6,tvar,evar,covs_adj2,paste0("Adjusted for ",paste(covs_adj2,collapse=", ")))
  save_plot(wrap_plots(qplots,ncol=3),"c1.Fig6.quantile_top.png",16,10,outdir=outdir)

  # Fig7: quantify rather than infer confounding from different P values.
  save_plot(plot_sameN_attenuation(attenuation_sameN,assoc_adj2),
            "c1.Fig7.attenuation_sameN.png",16,8.5,outdir=outdir)

  # Fig8: functional coherence of incident and baseline-prevalent adj2 scans.
  enrich<-run_functional_enrichment(assoc_adj2,layer,features_all)
  enrich_prev<-run_functional_enrichment(prevalent_adj2,layer,features_all)
  n_sig<-sum(is.finite(assoc_adj2$p.value)&assoc_adj2$p.value*nrow(assoc_adj2)<.05)
  n_sig_prev<-sum(is.finite(prevalent_adj2$p.value)&prevalent_adj2$p.value*nrow(prevalent_adj2)<.05)
  write_raw_csv(enrich,"c1.enrichment_incident_sig.csv",rawdir);write_raw_csv(enrich_prev,"c1.enrichment_prevalent_sig.csv",rawdir)
  pe_i<-plot_functional_enrichment(enrich,n_sig,assoc_adj2,layer,"a. Incident (Yin), adj2")
  pe_p<-plot_functional_enrichment(enrich_prev,n_sig_prev,prevalent_adj2,layer,"b. Baseline-prevalent (Yang), adj2")
  save_plot((pe_i/pe_p+plot_layout(guides="collect"))&theme(legend.position="right"),
            "c1.Fig8.enrich_sig.png",18,10.5,outdir=outdir)

  # Fig9 separates distal antecedent prediction from disease-state evidence.
  # GDF15 and PCSK9 are anchors by default, but the table is generated for all
  # assayed proteins/metabolites.
  directionality<-build_directionality_table(assoc_adj2,prevalent_adj2,duration_adj2,landmark_adj2,birthline_adj2,reverse_adj2)
  write_raw_csv(directionality,"c1.directionality_triage.csv",rawdir)
  save_plot(plot_directionality_triage(directionality,C1_DIRECTION_ANCHORS),
            "c1.Fig9.directionality_triage.png",16,9,outdir=outdir)
  save_plot(plot_landmark_and_attained_age(landmark_adj2,assoc_adj2,birthline_adj2,C1_DIRECTION_ANCHORS),
            "c1.Fig10.landmark_birthline_sensitivity.png",16,9,outdir=outdir)
  save_plot(plot_risk_window_scan(risk_window_adj2,C1_DIRECTION_ANCHORS),
            "c1.Fig11.diagnosis_window_riskset.png",18,11,outdir=outdir)
  save_plot(plot_directionality_supplement(directionality,C1_DIRECTION_ANCHORS),
            "c1.Fig12.directionality_detail.png",14,7.5,outdir=outdir)
  save_plot(plot_reverse_time_exploratory(reverse_adj2,prevalent_adj2,duration_adj2,C1_DIRECTION_ANCHORS),
            "c1.Fig13.reverse_time_exploratory.png",15.5,8,outdir=outdir)
  save_plot(plot_pgs_actual_concordance(pgs_concordance,C1_DIRECTION_ANCHORS),
            "c1.Fig14.pgs_actual_concordance.png",18,10.5,outdir=outdir)

  out<-list(meta=module_meta(layer,extra=list(N=nrow(dat),events=sum(dat[[evar]]==1,na.rm=TRUE),
              covs_use=covs_use_name,code_version=C1_CODE_VERSION,pgs_signature=pgs_signature,
              pgs_scan_signature=pgs$signature,
              pgs_matched=length(pgs$score_map),pgs_interpretation="fixed-at-conception score; not a biomarker measured at birth")),
            association=assoc,association_basic=assoc_basic,association_adj2=assoc_adj2,association_LE8=assoc_adj2,
            birthline_basic=birthline_basic,birthline_adj2=birthline_adj2,
            prevalent=assoc_prevalent,prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,
            reverse_prevalent=reverse_adj2,prevalent_duration=duration_adj2,landmark_incident=landmark_adj2,
            diagnosis_window_riskset=risk_window_adj2,directionality=directionality,
            attenuation=attenuation_sameN,enrichment=enrich,enrichment_prevalent=enrich_prev,
            pgs_status=pgs$status,pgs_incident=pgs_incident,pgs_prevalent=pgs_prevalent,
            pgs_attained_age=pgs_attained_age,pgs_incident_same_omic=pgs_incident_same,
            pgs_prevalent_same_omic=pgs_prevalent_same,pgs_attained_age_same_omic=pgs_attained_age_same,
            pgs_actual_concordance=pgs_concordance,
            input_feature_annotation_audit=input_feature_audit,top_features=top,gradient_top10_provenance=gradient_rank,clusters=cl_members,cluster_selection=cl_metrics)
  if(layer=="protein"){out$pwas_incident<-assoc;out$pwas_prevalent<-assoc_prevalent;out$top_proteins<-top}
  else {out$MWAS<-assoc;out$MWAS_prevalent<-assoc_prevalent}
  write_xlsx2(list(cohort=cohort,input_feature_audit=input_feature_audit,association=assoc,prevalent=assoc_prevalent,incident_basic=assoc_basic,incident_adj2=assoc_adj2,
                   birthline_basic=birthline_basic,birthline_adj2=birthline_adj2,
                   prevalent_basic=prevalent_basic,prevalent_adj2=prevalent_adj2,attenuation_sameN=attenuation_sameN,
                   reverse_prevalent=reverse_adj2,prevalent_duration=duration_adj2,incident_landmark=landmark_adj2,
                   diagnosis_window_riskset=risk_window_adj2,directionality=directionality,
                   pgs_status=pgs$status,pgs_incident=pgs_incident,pgs_prevalent=pgs_prevalent,
                   pgs_attained_age=pgs_attained_age,pgs_incident_same_omic=pgs_incident_same,
                   pgs_prevalent_same_omic=pgs_prevalent_same,pgs_attained_same_omic=pgs_attained_age_same,
                   pgs_actual_concordance=pgs_concordance,
                   gradient_top10_provenance=gradient_rank,
                   cluster_membership=cl_members,cluster_selection=cl_metrics,
                   enrichment_incident=enrich,enrichment_prevalent=enrich_prev),"c1.out.xlsx")
  finalize_outputs(LE8_JOB,outdir);saveRDS(out,selected_cache,compress="xz");out
}

if(prot_DO){invisible(run_c1_layer("protein"));gc(full=TRUE)}
if(met_DO){invisible(run_c1_layer("metabolite"));gc(full=TRUE)}
