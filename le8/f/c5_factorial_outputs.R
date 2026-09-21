# Paired layer increments and uncertainty from frozen out-of-fold predictions.
c5_factorial_outputs<-function(pred,yin,design,K,B,seed,primary_budget,primary_landmark){
  cuts<-c5_factorial_cutpoints();edges_all<-c5_factorial_edges(design)
  metrics<-contrasts<-shapley<-status<-strata<-correlations<-roc<-list()
  all_boot<-tolower(Sys.getenv("C5_FACTORIAL_BOOT_ALL","FALSE"))%in%c("true","1","yes")
  budgets<-sort(unique(pred$budget[pred$budget>0&startsWith(pred$model,"F_")]))
  for(L in sort(unique(pred$landmark)))for(arm in unique(pred$arm))for(k in budgets){
    ids<-if(inherits(yin,"c5_cached_cohorts"))as.character(yin$eid[yin$landmark==L])else c5_landmark(yin,L,Inf)$eid
    z<-c5_factorial_data(pred,design,k,L,arm,ids,K)
    id<-data.frame(budget_per_layer=k,landmark=L,arm=arm)
    status[[length(status)+1L]]<-cbind(id,z$audit)
    if(is.null(z$data))next
    d<-z$data;r<-z$risk;h<-d$horizon[1]
    edges<-edges_all[edges_all$model%in%colnames(r)&edges_all$reference%in%colnames(r),,drop=FALSE]
    est<-c5_factorial_stat(d,r,edges,h,cuts,z$lp)
    nboot<-if(all_boot||(k==primary_budget&&arm=="all_assays"))B else 0L
    message("C5 layer comparisons: L=",L,"; ",arm,"; k/layer=",k,"; N=",nrow(d),"; models=",ncol(r),"; paired bootstrap=",nboot)
    bs<-matrix(NA_real_,length(est),nboot,dimnames=list(names(est),NULL))
    if(nboot>0){
      g<-split(seq_len(nrow(d)),d$group);ng<-length(g)
      set.seed(seed+as.integer(100*L)+k)
      for(b in seq_len(nboot)){
        ix<-unlist(g[sample.int(ng,ng,replace=TRUE)],use.names=FALSE)
        bs[,b]<-c5_factorial_stat(d[ix,,drop=FALSE],r[ix,,drop=FALSE],edges,h,cuts,z$lp[ix,,drop=FALSE])[names(est)]
      }
    }
    intervals<-function(nm){
      v<-bs[nm,];v<-v[is.finite(v)]
      ci<-if(length(v)>=max(20,ceiling(.8*nboot)))quantile(v,c(.025,.975),names=FALSE)else c(NA_real_,NA_real_)
      data.frame(estimate=unname(est[nm]),lo=ci[1],hi=ci[2],bootstrap_success=length(v),bootstrap_requested=nboot,
        status=if(is.finite(est[nm]))"estimated"else"not estimable")
    }
    for(nm in names(est)){
      pieces<-strsplit(nm,"|",fixed=TRUE)[[1]];row<-cbind(id,horizon=h,intervals(nm))
      row$N<-nrow(d);row$events_by_horizon<-sum(d$event==1&d$time<=h)
      row$uncertainty<-"paired group bootstrap; fixed trained models; screening/training uncertainty excluded"
      if(pieces[1]=="model"){
        ds<-design[design$model==pieces[2],]
        row$model<-pieces[2];row$metric<-pieces[3];row$label<-ds$label
        row$protein_assays<-as.integer(ds$P)*k;row$metabolite_assays<-as.integer(ds$M)*k
        row$total_assays<-row$protein_assays+row$metabolite_assays;row$disease_PRS_scores<-as.integer(ds$G)
        row$clinical<-ds$C;metrics[[length(metrics)+1L]]<-row
      }else if(pieces[1]=="delta"){
        ee<-edges[edges$model==pieces[2]&edges$reference==pieces[3],]
        row$model<-pieces[2];row$reference<-pieces[3];row$metric<-pieces[4]
        row$added_layer<-ee$added_layer;row$clinical_background<-ee$clinical_background
        row$focus<-ee$focus&k==primary_budget&L==primary_landmark&arm=="all_assays"
        row$cutpoints<-paste(cuts,collapse=",")
        contrasts[[length(contrasts)+1L]]<-row
      }else{
        row$layer<-pieces[2];row$metric<-if(pieces[3]=="Brier")"Brier_reduction"else pieces[3]
        row$interpretation<-"mean predictive increment across all layer-addition orders conditional on Clinical; not causal contribution"
        shapley[[length(shapley)+1L]]<-row
      }
    }
    if(k!=primary_budget)next
    iw<-le8_ipcw(d$time,d$event,h)
    if(all(c("F_P","F_M","F_G")%in%names(z$models))&&iw$status=="ok"){
      flags<-data.frame(P=z$models$F_P$high_training_q75,M=z$models$F_M$high_training_q75,G=z$models$F_G$high_training_q75)
      labels<-apply(flags,1,function(x)paste0("P",as.integer(x[1])," M",as.integer(x[2])," G",as.integer(x[3])))
      for(s in sort(unique(labels))){
        ii<-which(labels==s);obs<-sum(iw$w[ii]*iw$y[ii])/sum(iw$w[ii])
        q<-tryCatch(summary(survival::survfit(survival::Surv(time,event)~1,data=d[ii,,drop=FALSE]),times=h,extend=FALSE),error=function(e)NULL)
        km<-lo<-hi<-NA_real_
        if(!is.null(q)&&length(q$surv)==1){km<-1-q$surv;lo<-1-q$upper;hi<-1-q$lower}
        strata[[length(strata)+1L]]<-cbind(id,horizon=h,stratum=s,N=length(ii),events=sum(iw$y[ii]),
          observed_ipcw=if(is.finite(obs))obs else NA_real_,KM_net_risk=km,KM_lo=lo,KM_hi=hi,
          interpretation="1 = above fold-training 75th percentile of the layer-only predicted LP; net risk, no causal subtype")
      }
    }
    # Correlations within folds avoid mixing differently scaled risk scores.
    single<-intersect(c("F_P","F_M","F_G"),colnames(r))
    if(length(single)>1)for(a in seq_len(length(single)-1L))for(b in seq.int(a+1,length(single)))for(fd in sort(unique(d$fold))){
      ii<-d$fold==fd;rho<-suppressWarnings(cor(z$lp[ii,single[a]],z$lp[ii,single[b]],method="spearman"))
      correlations[[length(correlations)+1L]]<-cbind(id,fold=fd,score_a=single[a],score_b=single[b],N=sum(ii),Spearman_r=rho)
    }
    if(iw$status=="ok"&&arm=="all_assays")for(nm in colnames(r)){
      p<-r[,nm];th<-c(Inf,sort(unique(quantile(p,seq(0,1,length.out=101),names=FALSE)),decreasing=TRUE),-Inf)
      aa<-vapply(th,function(t)sum(iw$w*iw$y*(p>=t))/sum(iw$w*iw$y),numeric(1))
      bb<-vapply(th,function(t)sum(iw$w*(1-iw$y)*(p>=t))/sum(iw$w*(1-iw$y)),numeric(1))
      roc[[length(roc)+1L]]<-cbind(id,horizon=h,model=nm,data.frame(threshold=th,sensitivity=aa,false_positive_rate=bb))
    }
  }
  list(factorial_metrics=bind_rows(metrics),factorial_contrasts=bind_rows(contrasts),factorial_shapley=bind_rows(shapley),
    factorial_status=bind_rows(status),factorial_strata=bind_rows(strata),factorial_score_correlations=bind_rows(correlations),factorial_ROC=bind_rows(roc))
}
