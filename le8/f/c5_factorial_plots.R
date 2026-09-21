# Two compact six-panel figures; all numerical outputs remain in root CSVs.
c5_factorial_plots<-function(t,design,root,k,L){
  blank<-function(title)ggplot()+annotate("text",x=0,y=0,label="Not estimable / input unavailable")+theme_void()+labs(title=title)
  style<-theme_bw(base_size=10)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
  short<-function(x)gsub("([CPMG])(?=[CPMG])","\\1 + ",sub("^F_","",x),perl=TRUE)
  forest<-function(d,title,xlab){
    if(!nrow(d))return(blank(title))
    ggplot(d,aes(estimate,reorder(label,estimate)))+
      geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.15,na.rm=TRUE)+geom_point()+
      labs(title=title,x=xlab,y=NULL)+style
  }
  save<-function(p,nm){for(ext in c("png","pdf"))ggsave(file.path(root,paste0(nm,".",ext)),p,width=18,height=12,dpi=240,limitsize=FALSE)}
  m<-t$factorial_metrics;d<-t$factorial_contrasts;s<-t$factorial_shapley
  if(!nrow(m)){save(blank("Factorial comparison unavailable"),"c5.Fig23.PRS_prot_met_prediction");return(invisible(NULL))}
  m$label<-short(m$model)
  q<-m|>filter(budget_per_layer==k,landmark==L,arm=="all_assays",status=="estimated")
  a<-forest(q|>filter(metric=="C_index"),"A. All available model combinations","Within-fold, horizon-restricted Harrell C")
  b<-forest(q|>filter(metric=="AUC"),"B. Fixed-horizon discrimination","Out-of-fold IPCW AUC")
  x<-if(nrow(d))d|>filter(budget_per_layer==k,landmark==L,arm=="all_assays",model=="F_CPMG",reference%in%c("F_CPM","F_CPG","F_CMG"),metric=="C_index")|>
    mutate(label=paste0("Add ",added_layer," to ",short(reference)))else data.frame()
  c<-forest(x,"C. Each layer added after the other two","Paired delta C (conditional 95% CI)")
  if(nrow(x))c<-c+geom_vline(xintercept=0,linetype=2)
  clinical_models<-design$model[design$C&design$available]
  x<-t$calibration
  if(nrow(x))x<-x|>filter(model%in%clinical_models,arm=="all_assays",landmark==L,budget%in%c(0,k))
  e<-if(!nrow(x))blank("D. Calibration")else ggplot(x,aes(predicted,observed_ipcw,color=short(model)))+
    geom_abline(slope=1,intercept=0,linetype=2)+geom_line()+geom_point(size=1)+
    labs(title="D. Calibration on held-out predictions",x="Predicted net risk",y="Observed IPCW risk",color=NULL)+style
  x<-t$decision_curves
  if(nrow(x))x<-x|>filter(model%in%clinical_models,arm=="all_assays",landmark==L,budget%in%c(0,k))
  f<-if(!nrow(x))blank("E. Decision curves")else ggplot(x,aes(threshold,net_benefit,color=short(model)))+
    geom_line()+geom_line(aes(y=treat_all),color="grey50",linetype=2)+geom_hline(yintercept=0,linetype=3)+
    labs(title="E. Exploratory net benefit",x="Risk threshold",y="Net benefit",color=NULL)+style
  x<-if(nrow(s))s|>filter(budget_per_layer==k,landmark==L,arm=="all_assays",metric=="AUC",status=="estimated")|>mutate(label=layer)else data.frame()
  g<-forest(x,"F. Average contribution across addition orders","AUC gain allocated across P / M / G")
  if(nrow(x))g<-g+geom_vline(xintercept=0,linetype=2)
  caption<-paste0("C=clinical; P=proteins; M=metabolites; G=disease PRS, distinct from biomarker PGS. k=",k,
    " per measured layer; P+M uses ",2*k," assays. Shared participants and folds. Landmark ",L,
    ". Conditional bootstrap CI excludes training uncertainty. G discovery overlap: see c5.prs_provenance.csv. Death censored; net risks, not competing-risk CIFs.")
  caption<-paste(strwrap(caption,width=160),collapse="\n")
  annotation_theme<-theme(plot.caption=element_text(hjust=0,size=9),plot.title=element_text(size=15,face="bold"))
  save((a|b|c)/(e|f|g)+plot_annotation(title=paste(Y,"— disease PRS, proteome and metabolome"),caption=caption,theme=annotation_theme),"c5.Fig23.PRS_prot_met_prediction")
  x<-t$factorial_ROC
  a<-blank("A. Time-dependent ROC")
  if(nrow(x)){
    x<-x|>filter(landmark==L,model%in%clinical_models)
    if(nrow(x))a<-ggplot(x,aes(false_positive_rate,sensitivity,color=short(model)))+geom_line()+
      geom_abline(slope=1,intercept=0,linetype=2)+coord_equal()+labs(title="A. IPCW ROC on the common cohort",x="False positive rate",y="Sensitivity",color=NULL)+style
  }
  x<-t$factorial_score_correlations
  b<-blank("B. Layer-score correlations")
  if(nrow(x)){
    x<-x|>filter(landmark==L,arm=="all_assays")|>group_by(score_a,score_b)|>summarise(r=mean(Spearman_r,na.rm=TRUE),.groups="drop")
    if(nrow(x))b<-ggplot(x,aes(short(score_a),short(score_b),fill=r))+geom_tile()+geom_text(aes(label=sprintf("%.2f",r)))+
      scale_fill_gradient2(low="#326da8",mid="white",high="#ba493e",limits=c(-1,1))+
      labs(title="B. Mean within-fold Spearman correlation",x=NULL,y=NULL,fill="r")+style
  }
  x<-t$factorial_strata
  c<-blank("C. Three-layer risk strata")
  if(nrow(x)){
    x<-x|>filter(landmark==L,arm=="all_assays")
    if(nrow(x))c<-ggplot(x,aes(stratum,KM_net_risk))+geom_col(fill="#437a9c")+
      geom_errorbar(aes(ymin=KM_lo,ymax=KM_hi),width=.15,na.rm=TRUE)+
      geom_text(aes(label=paste0(events,"/",N)),vjust=-.4,size=2.5)+
      labs(title="C. High/lower layer-score combinations",x="1 = above training 75th percentile; labels = events / N",y="KM net risk (95% CI)")+style+
      theme(axis.text.x=element_text(angle=45,hjust=1))
  }
  x<-m|>filter(budget_per_layer==k,arm=="all_assays",metric=="C_index",model%in%clinical_models,status=="estimated")
  e<-if(!nrow(x))blank("D. Lead time")else ggplot(x,aes(landmark,estimate,color=short(model)))+
    geom_line()+geom_point()+geom_errorbar(aes(ymin=lo,ymax=hi),width=.1,na.rm=TRUE)+
    labs(title="D. Refit at each event-free landmark",x="Years since sampling",y="Horizon-restricted C-index",color=NULL)+style
  x<-if(nrow(d))d|>filter(budget_per_layer==k,landmark==L,arm=="all_assays",model=="F_CPMG",reference%in%c("F_CPM","F_CPG","F_CMG"),
    metric%in%c("IDI","continuous_NRI","categorical_NRI"))|>mutate(label=paste("Add",added_layer))else data.frame()
  f<-if(!nrow(x))blank("E. Reclassification")else ggplot(x,aes(estimate,label))+
    geom_vline(xintercept=0,linetype=2)+geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.1,na.rm=TRUE)+geom_point()+
    facet_wrap(~metric,scales="free_x")+labs(title="E. Paired reclassification (exploratory)",x="Estimate and conditional 95% CI",y=NULL)+style
  x<-m|>filter(landmark==L,arm=="all_assays",metric=="C_index",model%in%c("F_CP","F_CM","F_CPM","F_CPMG"),status=="estimated")
  g<-if(!nrow(x))blank("F. Assay requirements")else ggplot(x,aes(total_assays,estimate,color=short(model)))+geom_line()+geom_point()+
    labs(title="F. Performance versus actual assay count",x="Total measured assays (PGS reported separately)",y="Horizon-restricted C-index",color=NULL)+style
  save((a|b|c)/(e|f|g)+plot_annotation(title=paste(Y,"— complementarity, lead time and assay requirements"),caption=caption,theme=annotation_theme),"c5.Fig24.PRS_prot_met_complementarity")
}
