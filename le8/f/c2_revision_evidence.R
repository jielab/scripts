# C2 result interpretation and score-weight provenance.
grade_c2_evidence <- function(mr,layer,mrlink2=list(results=tibble())) {
  primary<-if(layer=="protein")"cis"else"local"
  z<-mr|>filter(analysis==primary)
  if(!nrow(z))return(z)
  z<-z|>mutate(
    heterogeneous=ifelse(n_IV>1&is.finite(Q_p),Q_p<.05,NA),
    heterogeneity_tested=n_IV>1&is.finite(Q_p),
    egger_tested=n_IV>=3&is.finite(egger_intercept_p),
    egger_warning=ifelse(egger_tested,egger_intercept_p<.05,NA),
    weighted_median_concordant=ifelse(n_IV>=3&is.finite(tsmr_weighted_median_b),
      sign(tsmr_weighted_median_b)==sign(b),NA),
    steiger_ok=NA, many_IV_warning=n_IV>10,
    instrument_architecture=case_when(n_IV==0~"unavailable",n_IV==1~"single IV",
      n_IV<=5~"oligogenic (2-5 IVs)",n_IV<=10~"multi-IV (6-10)",TRUE~"many-IV (>10)"),
    evidence_grade=case_when(
      !is.finite(pval)~"U: unavailable",
      FDR_all<.05&n_IV==1~"S: single-IV support",
      FDR_all<.05&heterogeneity_tested&egger_tested&!heterogeneous&!egger_warning~"A: diagnostics-compatible",
      FDR_all<.05&coalesce(heterogeneous,FALSE)&coalesce(weighted_median_concordant,FALSE)~"B: heterogeneous",
      FDR_all<.05~"C: sensitivity unresolved",
      pval<.05~"D: nominal only",TRUE~"E: not detected"),
    grade_reason=case_when(
      n_IV==1~"Wald estimate; Egger, weighted-median and heterogeneity tests unavailable",
      evidence_grade=="A: diagnostics-compatible"~"FDR support; no detected Q/Egger warning, not proof of instrument validity",
      evidence_grade=="B: heterogeneous"~"FDR support with heterogeneity; weighted-median direction agrees",
      evidence_grade=="U: unavailable"~"Not tested / missing valid input",
      TRUE~"Provisional statistical evidence; requires same-region/signal corroboration"))
  rr<-as_tibble(mrlink2$results%||%tibble())
  z$MRlink2_p<-NA_real_;z$MRlink2_FDR<-NA_real_
  if(nrow(rr)) {
    tc<-pick_col_local(names(rr),c("^TRAIT$","^EXPOSURE$"));pc<-pick_col_local(names(rr),c("^P\\(ALPHA\\)$","^P_ALPHA$"))
    if(!is.na(tc)&&!is.na(pc)) {
      v<-tibble(exposure=as.character(rr[[tc]]),p=as.numeric(rr[[pc]]))|>
        filter(is.finite(p))|>group_by(exposure)|>
        summarise(MRlink2_p=pmin(1,min(p)*n()),.groups="drop")|>
        mutate(MRlink2_FDR=p.adjust(MRlink2_p,"BH"))
      i<-match(z$exposure,v$exposure);z$MRlink2_p<-v$MRlink2_p[i];z$MRlink2_FDR<-v$MRlink2_FDR[i]
    }
  }
  arrange(z,pval)
}
plot_c2_evidence_grades <- function(z,layer) {
  if(!nrow(z))return(blank_plot("C2 evidence grades","No estimable result"))
  chosen<-unique(c(C2_FIXED_TOP,head(z$exposure[order(z$pval)],24)))
  d<-z|>filter(exposure%in%chosen)|>mutate(exposure=factor(exposure,levels=rev(chosen)),
    lo=b-1.96*se,hi=b+1.96*se)
  pa<-ggplot(d,aes(b,exposure,color=evidence_grade))+geom_vline(xintercept=0)+
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=.12)+geom_point(aes(shape=instrument_architecture),size=2)+
    labs(title="a. Marginal-effect MR and uncertainty",x="MR effect (GWAS outcome scale)",y=NULL,color=NULL,shape=NULL)+theme_5c(9)
  pb<-z|>count(evidence_grade)|>ggplot(aes(n,reorder(evidence_grade,n),fill=evidence_grade))+
    geom_col(show.legend=FALSE)+labs(title="b. Missing diagnostics are not passed tests",x="Biomarkers",y=NULL)+theme_5c(9)
  pa|pb
}
build_genetic_score_manifest <- function(iv_list,layer) {
  imap_dfr(iv_list,function(obj,feature) {
    z<-as_tibble(obj$score_instruments%||%tibble());if(!nrow(z))return(tibble())
    z|>transmute(feature,SNP,CHR,POS,effect_allele=EA,other_allele=NEA,weight=BETA,SE,P,EAF,N,
      component=paste0("COJO PGS / ",analysis),effect_type="joint PGS weight; not MR input",
      interpretation="Existing prot.pgs.rds/met.pgs.rds unchanged. Not a protein concentration measured at birth.")
  })
}

if(!exists(".le8_original_run_mrlink2_step",envir=.GlobalEnv,inherits=FALSE))
  .le8_original_run_mrlink2_step<-run_mrlink2_step
run_mrlink2_step <- function(...) {
  tryCatch(.le8_original_run_mrlink2_step(...),error=function(e) {
    warning("Optional MR-link-2 failed: ",conditionMessage(e),call.=FALSE);invisible(1L)
  })
}
if(!exists(".le8_original_integrate_c2",envir=.GlobalEnv,inherits=FALSE))
  .le8_original_integrate_c2<-integrate_c2_directionality
integrate_c2_directionality <- function(mr,c1,dan=list(),reverse_mr=tibble()) {
  # Do not pick whichever of cis or trans has the smaller P as causal evidence.
  z<-.le8_original_integrate_c2(mr[mr$analysis%in%c("cis","local"),,drop=FALSE],c1,dan,reverse_mr)
  if(nrow(z)) {
    if("causal_interpretation"%in%names(z))
      z$causal_interpretation<-gsub("no causal support","not established / not detected",z$causal_interpretation,fixed=TRUE)
    z$MR_scope<-"cis/local primary; trans/distal remains in its own evidence table"
  }
  z
}

if(!exists(".le8_original_instrument_availability",envir=.GlobalEnv,inherits=FALSE))
  .le8_original_instrument_availability<-instrument_availability_audit
instrument_availability_audit <- function(...) {
  z<-.le8_original_instrument_availability(...)
  if(nrow(z)&&"n_independent"%in%names(z)) {
    z$n_COJO_selected<-z$n_independent;z$n_independent<-NA_integer_
    z$LD_interpretation<-"Input availability is not LD verification; final pruned counts and fallback status are in c2.MR_all.csv"
  }
  z
}
