# Fixed-panel comparisons of Clinical (C), proteome (P), metabolome (M),
# and the endpoint's disease PRS (G). Biomarker PGS are a separate analysis.
c5_factorial_design<-function(has_prs=TRUE){
  z<-expand.grid(C=c(FALSE,TRUE),P=c(FALSE,TRUE),M=c(FALSE,TRUE),G=c(FALSE,TRUE))
  z<-z[rowSums(z)>0,,drop=FALSE]
  z$model<-apply(z,1,function(r)paste0("F_",paste(c("C","P","M","G")[as.logical(r)],collapse="")))
  z$label<-apply(z[c("C","P","M","G")],1,function(r)
    paste(c("Clinical","Protein","Metabolite","Disease PRS")[as.logical(r)],collapse=" + "))
  z$available<-!z$G|has_prs
  z$status<-ifelse(z$available,"scheduled","endpoint PRS unavailable; model not fitted")
  z$budget_definition<-"k assays PER included measured layer; P+M uses 2k, G is one precomputed endpoint score"
  z[order(rowSums(z[c("C","P","M","G")]),z$model),,drop=FALSE]
}
c5_factorial_edges<-function(design){
  z<-design[design$available,,drop=FALSE];ans<-list()
  for(i in seq_len(nrow(z)))for(j in seq_len(nrow(z))){
    a<-as.logical(unlist(z[i,c("C","P","M","G")],use.names=FALSE));b<-as.logical(unlist(z[j,c("C","P","M","G")],use.names=FALSE))
    if(sum(a&!b)==1&&!any(b&!a))ans[[length(ans)+1L]]<-data.frame(model=z$model[i],reference=z$model[j],
      added_layer=c("C","P","M","G")[a&!b],clinical_background=b[1],focus=z$model[i]=="F_CPMG"&&z$model[j]%in%c("F_CPM","F_CPG","F_CMG"))
  }
  do.call(rbind,ans)
}
c5_factorial_cutpoints<-function(){
  s<-Sys.getenv("C5_RISK_CUTS","");if(!nzchar(trimws(s)))return(numeric())
  v<-suppressWarnings(as.numeric(strsplit(s,",",fixed=TRUE)[[1]]))
  if(!length(v)||any(!is.finite(v)|v<=0|v>=1)||anyDuplicated(v))stop("C5_RISK_CUTS requires unique probabilities strictly between 0 and 1")
  sort(v)
}
# Compare ONLY within-fold pairs, because separate Cox models can have
# different LP scales. Administratively restrict all C-index calculations
# to the same prediction horizon as AUC/Brier.
c5_within_fold_c<-function(time,event,risk,fold,horizon){
  count<-c(0,0,0)
  for(f in unique(fold)){
    ii<-which(fold==f);if(length(ii)<2)next
    d<-data.frame(tt=pmin(time[ii],horizon),ee=as.integer(event[ii]==1&time[ii]<=horizon),rr=risk[ii])
    obj<-tryCatch(survival::concordance(survival::Surv(tt,ee)~rr,data=d,reverse=TRUE),error=function(e)NULL)
    if(!is.null(obj))count<-count+as.numeric(obj$count[c("concordant","discordant","tied.x")])
  }
  if(any(!is.finite(count))||sum(count)<=0)return(NA_real_)
  (count[1]+.5*count[3])/sum(count)
}
c5_factorial_reclassification<-function(new,old,y,w,cuts=numeric()){
  meanpart<-function(x,case){ii<-y==case&w>0;if(!any(ii)||sum(w[ii])<=0)return(NA_real_);sum(w[ii]*x[ii])/sum(w[ii])}
  direction<-sign(new-old);ne<-meanpart(direction,1);nn<--meanpart(direction,0)
  out<-c(IDI=meanpart(new-old,1)-meanpart(new-old,0),continuous_NRI=ne+nn,NRI_event=ne,NRI_nonevent=nn)
  if(length(cuts)){
    change<-sign(findInterval(new,cuts)-findInterval(old,cuts))
    ce<-meanpart(change,1);cn<--meanpart(change,0)
    out<-c(out,categorical_NRI=ce+cn,categorical_NRI_event=ce,categorical_NRI_nonevent=cn)
  }
  out
}
c5_factorial_shapley<-function(values){
  # All 3! orderings of P/M/G added to C; exact layer-subset decomposition.
  players<-c("P","M","G");permutations<-list(c("P","M","G"),c("P","G","M"),c("M","P","G"),c("M","G","P"),c("G","P","M"),c("G","M","P"))
  needed<-c("F_C","F_CP","F_CM","F_CG","F_CPM","F_CPG","F_CMG","F_CPMG")
  if(!all(needed%in%names(values))||any(!is.finite(values[needed])))return(setNames(rep(NA_real_,3),players))
  out<-setNames(numeric(3),players)
  for(order in permutations){
    included<-character();previous<-"F_C"
    for(p in order){
      included<-c(included,p);current<-paste0("F_C",paste(players[players%in%included],collapse=""))
      out[p]<-out[p]+(values[[current]]-values[[previous]])/length(permutations);previous<-current
    }
  }
  out
}
c5_factorial_stat<-function(d,risk,edges,horizon,cuts,lp=risk){
  iw<-le8_ipcw(d$time,d$event,horizon)
  model<-matrix(NA_real_,ncol(risk),3,dimnames=list(colnames(risk),c("AUC","C_index","Brier")))
  for(j in seq_len(ncol(risk))){
    r<-risk[,j]
    if(iw$status=="ok"){
      model[j,"AUC"]<-le8_weighted_auc(r,iw$y,iw$w)
      model[j,"Brier"]<-mean(iw$w*(iw$y-r)^2)
    }
    model[j,"C_index"]<-c5_within_fold_c(d$time,d$event,lp[,j],d$fold,horizon)
  }
  out<-numeric()
  for(nm in rownames(model))for(mt in colnames(model))out[paste("model",nm,mt,sep="|")]<-model[nm,mt]
  reclass_names<-names(c5_factorial_reclassification(c(.1,.2),c(.1,.2),c(0,1),c(1,1),cuts))
  for(i in seq_len(nrow(edges))){
    a<-edges$model[i];b<-edges$reference[i]
    for(mt in colnames(model))out[paste("delta",a,b,mt,sep="|")]<-model[a,mt]-model[b,mt]
    re<-if(iw$status=="ok")c5_factorial_reclassification(risk[,a],risk[,b],iw$y,iw$w,cuts)else setNames(rep(NA_real_,length(reclass_names)),reclass_names)
    for(mt in names(re))out[paste("delta",a,b,mt,sep="|")]<-re[[mt]]
  }
  for(mt in colnames(model)){
    value<-model[,mt];if(mt=="Brier")value<--value
    sh<-c5_factorial_shapley(value)
    for(p in names(sh))out[paste("shapley",p,mt,sep="|")]<-sh[[p]]
  }
  out
}
c5_factorial_data<-function(pred,design,k,L,arm,expected_ids,K){
  models<-design$model[design$available];lst<-list();audit<-list()
  for(nm in models){
    ds<-design[design$model==nm,];budget<-if(ds$P||ds$M)k else 0
    x<-pred[pred$model==nm&pred$budget==budget&pred$landmark==L&pred$arm==arm,,drop=FALSE]
    ok<-nrow(x)==length(expected_ids)&&!anyDuplicated(x$eid)&&setequal(x$eid,expected_ids)&&
      length(unique(x$fold))==K&&all(is.finite(x$risk))&&all(is.finite(x$lp))
    audit[[nm]]<-data.frame(model=nm,status=if(ok)"complete"else"missing/failed outer-fold predictions; excluded from comparisons",N=nrow(x))
    if(ok)lst[[nm]]<-x[match(expected_ids,x$eid),,drop=FALSE]
  }
  if(!length(lst))return(list(audit=do.call(rbind,audit),data=NULL,risk=NULL))
  d<-lst[[1]]
  for(x in lst)if(!all(vapply(c("eid","fold","time","event","group"),function(v)identical(x[[v]],d[[v]]),logical(1))))
    stop("Factorial models do not share the same outcomes/folds/participants")
  risk<-do.call(cbind,lapply(lst,`[[`,"risk"));colnames(risk)<-names(lst)
  lp<-do.call(cbind,lapply(lst,`[[`,"lp"));colnames(lp)<-names(lst)
  list(data=d,risk=risk,lp=lp,models=lst,audit=do.call(rbind,audit))
}
