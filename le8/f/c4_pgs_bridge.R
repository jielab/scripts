# Matched biomarker PGS -> measured omic. This is calibration/association,
# not a disease PRS -> all proteins scan and not identification of mediation.
le8_pgs_bridge <- function(dat, features, fold, covars, disease, membership, layer, rawdir) {
  empty<-list(scan=tibble(),bridges=tibble(),status=tibble(status="PGS input unavailable"))
  emit<-function(out) {
    for(nm in names(out))write_raw_csv(out[[nm]],paste0("c4.matched_PGS_",nm,".csv"),rawdir)
    out
  }
  f<-Sys.getenv("C4_PGS_FILE",unset=find_c1_pgs_file(layer))
  if(is.na(f)||!nzchar(f)||!file.exists(f))return(emit(empty))
  message("C4: load matched omic PGS: ",f)
  scores<-read_c1_pgs(f)
  if(anyDuplicated(scores$eid)||anyDuplicated(dat$eid))stop("Duplicate eid in omic/PGS input")
  mapping<-map_c1_pgs_columns(features,setdiff(names(scores),"eid"))
  if(!length(mapping))return(emit(empty))
  idx<-match(as.character(dat$eid),scores$eid)
  scores<-scores[idx,unname(mapping),drop=FALSE]
  covars<-intersect(covars,names(dat));rows<-list()
  for(h in unique(fold))for(feature in names(mapping)) {
    y<-suppressWarnings(as.numeric(dat[[feature]]));g<-as.numeric(scores[[mapping[[feature]]]])
    keep<-fold==h&is.finite(y)&is.finite(g)&complete.cases(dat[,covars,drop=FALSE])
    n<-sum(keep);r<-p<-NA_real_
    if(n>=100) {
      # Subset covariates only: copying thousands of omic columns for every
      # QR fit made a genome-wide bridge needlessly slow and memory intensive.
      M<-model.matrix(reformulate(covars),dat[keep,covars,drop=FALSE]);q<-qr(M)
      yr<-qr.resid(q,y[keep]);gr<-qr.resid(q,g[keep]);den<-sqrt(sum(yr^2)*sum(gr^2));df<-n-q$rank-1
      if(is.finite(den)&&den>0&&df>2){r<-max(-.999999,min(.999999,sum(yr*gr)/den));p<-2*pt(abs(r)*sqrt(df/(1-r*r)),df,lower.tail=FALSE)}
    }
    rows[[length(rows)+1L]]<-tibble(feature,score_column=mapping[[feature]],split=h,n,r,p.value=p)
    if(length(rows)%%500L==0L)message("C4 matched PGS: ",length(rows)," / ",length(unique(fold))*length(mapping)," fits")
  }
  rm(scores);invisible(gc())
  scan<-bind_rows(rows)|>group_by(split)|>mutate(FDR=p.adjust(p.value,"BH"))|>ungroup()
  a<-scan|>filter(split=="discovery")|>select(feature,score_column,n_disc=n,r_disc=r,p_disc=p.value,FDR_disc=FDR)
  b<-scan|>filter(split=="replication")|>select(feature,n_rep=n,r_rep=r,p_rep=p.value,FDR_rep=FDR)
  z<-left_join(a,b,by="feature")|>mutate(replicated=coalesce(FDR_disc<.05&FDR_rep<.05&r_disc*r_rep>0,FALSE),
    partial_R2_rep=r_rep^2,edge="biomarker PGS -> measured biomarker",
    causal_mediation_identified=FALSE,source_file=f)
  if(nrow(membership))z<-left_join(z,membership|>select(any_of(c("feature","primary_component","set","strict_YS"))),by="feature")
  if(nrow(disease))z<-left_join(z,disease|>transmute(feature=term,observed_disease_p=p.value),by="feature")
  c1file<-file.path(dirname(rawdir),"c1_correlate","c1.res.rds")
  if(file.exists(c1file)) {
    c1<-readRDS(c1file);pg<-c1$pgs_incident%||%tibble()
    if(nrow(pg))z<-left_join(z,pg|>transmute(feature=term,PGS_disease_beta=beta,PGS_disease_p=p.value,PGS_disease_FDR=FDR),by="feature")
  }
  status<-tibble(status="ok",matched_scores=length(mapping),replicated_edges=sum(z$replicated),
    note="Split-replicated partial association; discovery GWAS overlap and score transportability require separate assessment")
  for(nm in c("scan","bridges","status")) {
    value<-switch(nm,scan=scan,bridges=z,status=status)
    write_raw_csv(value,paste0("c4.matched_PGS_",nm,".csv"),rawdir)
  }
  list(scan=scan,bridges=z,status=status)
}
