# LE8-first selection and transparent comparison of distinct assay panels.
# Neither function below accesses disease outcomes or validation measurements.
focus_heterogeneity<-function(contrasts,boot){
  if(!nrow(boot)||!nrow(contrasts))return(tibble())
  lo<-"Low baseline inflammation";hi<-"High baseline inflammation"
  specs<-contrasts|>filter(stratum%in%c(lo,hi))|>distinct(model,reference,landmark,horizon)
  out<-map_dfr(seq_len(nrow(specs)),function(i){
    s<-specs[i,];q<-contrasts|>filter(model==s$model,reference==s$reference,landmark==s$landmark,stratum%in%c(lo,hi))
    if(!all(c(lo,hi)%in%q$stratum))return(tibble())
    a<-boot|>filter(model==s$model,landmark==s$landmark,stratum%in%c(lo,hi))
    b<-boot|>filter(model==s$reference,landmark==s$landmark,stratum%in%c(lo,hi))
    z<-inner_join(a,b,by=c("replicate","stratum","landmark","horizon"),suffix=c("_a","_b"))|>
      mutate(delta=AUC_a-AUC_b)|>select(replicate,stratum,delta)|>pivot_wider(names_from=stratum,values_from=delta)
    if(!all(c(lo,hi)%in%names(z)))return(tibble())
    v<-z[[hi]]-z[[lo]];v<-v[is.finite(v)]
    if(length(v)<20)return(tibble())
    est<-q$delta_AUC[match(hi,q$stratum)]-q$delta_AUC[match(lo,q$stratum)]
    se<-sd(v);p<-if(is.finite(se)&&se>0)2*pnorm(abs(est/se),lower.tail=FALSE)else NA_real_
    bind_cols(s,tibble(delta_high_minus_low=est,lo=unname(quantile(v,.025)),hi=unname(quantile(v,.975)),p_heterogeneity=p,
      inference="Frozen-fit exploratory difference in delta AUC; independent within-stratum bootstrap; not inflammatory causation"))
  })
  if(nrow(out))out$FDR<-p.adjust(out$p_heterogeneity,"BH")
  out
}

focus_balanced_panel <- function(membership,k,components) {
  z<-membership|>filter(selected%in%TRUE,is.finite(r1),is.finite(r2))|>
    mutate(strength=pmin(abs(r1),abs(r2)))|>
    arrange(desc(strength),feature)|>distinct(feature,.keep_all=TRUE)
  queues<-lapply(components,function(cmp)z$feature[z$component==cmp])
  names(queues)<-components;panel<-character()
  while(length(panel)<k&&any(lengths(queues)>0)) {
    for(cmp in components)if(length(queues[[cmp]])&&length(panel)<k) {
      panel<-c(panel,queues[[cmp]][1]);queues[[cmp]]<-queues[[cmp]][-1]
    }
  }
  unique(panel)
}

focus_panel_audit <- function(designs,membership,components) {
  coverage<-map_dfr(designs,function(ds) {
    mm<-membership|>filter(cohort==ds$cohort,selected%in%TRUE,feature%in%ds$features)
    counts<-table(factor(mm$component,levels=components))
    tibble(model=ds$name,component=components,assays=as.integer(counts),
      assay_budget=length(ds$features),selection_cohort=ds$cohort,
      scope="Replicated partial association labels, not proof of cell origin or intervention responsiveness")
  })
  overlap<-map_dfr(seq_along(designs),function(i) {
    if(i==1)return(tibble())
    map_dfr(seq_len(i-1L),function(j) {
      a<-designs[[i]];b<-designs[[j]]
      if(a$budget!=b$budget||a$budget==0)return(tibble())
      tibble(model=a$name,reference=b$name,budget=a$budget,
        intersection=length(intersect(a$features,b$features)),
        union=length(union(a$features,b$features)),
        Jaccard=length(intersect(a$features,b$features))/length(union(a$features,b$features)),
        identical_panel=setequal(a$features,b$features))
    })
  })
  list(panel_coverage=coverage,panel_overlap=overlap)
}
