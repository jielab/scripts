# Compact multi-panel figures supplement (and retain) the existing figures.
c5_final_panels<-function(t,root,k){
  blank<-function(title)ggplot()+annotate("text",x=0,y=0,label="Not estimable / input unavailable")+theme_void()+labs(title=title)
  style<-theme_bw(base_size=9)+theme(legend.position="bottom",plot.title=element_text(face="bold"))
  save<-function(p,name){for(ext in c("png","pdf"))ggsave(file.path(root,paste0(name,".",ext)),p,width=18,height=12,dpi=240,limitsize=FALSE)}
  short<-function(x)gsub("Clinical_","C + ",gsub("ProtMet","P+M",gsub("Matched","",x)))
  m<-t$metrics|>filter(status=="ok",arm=="all_assays",budget%in%c(0,k))
  keep<-c("Clinical","Clinical_PRS","Clinical_Protein_NS","Clinical_Metabolite_NS","Clinical_ProtMet_NS","Clinical_ProtMet_PRS_NS")
  a<-if(!nrow(m))blank("A. Joint modality comparisons")else ggplot(m|>filter(model%in%keep),aes(AUC,short(model),color=factor(landmark)))+
    geom_point()+labs(title="A. Same participants; declared assay budgets",x="Out-of-fold IPCW AUC",y=NULL,color="Landmark")+style
  d<-t$paired_contrasts|>filter(budget==k,arm=="all_assays",grepl("PGS|sharedBudget",paste(model,reference)))
  b<-if(!nrow(d))blank("B. Incremental PGS / modality information")else ggplot(d,aes(delta_AUC,paste(short(model),"vs",short(reference)),color=factor(landmark)))+
    geom_vline(xintercept=0,linetype=2)+geom_errorbar(aes(xmin=AUC_lo,xmax=AUC_hi),orientation="y",width=.15)+geom_point()+
    labs(title="B. Paired increments; measured reference matches PGS coverage",x="Delta AUC with conditional 95% CI",y=NULL,color="Landmark")+style
  d<-t$pgs_reconstruction
  c<-if(!nrow(d))blank("C. Matched PGS capture of measured omics")else ggplot(d|>filter(arm=="all_assays"),aes(sub("__.*$","",feature),R2,color=factor(landmark)))+
    geom_hline(yintercept=0,linetype=2)+geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=.1,height=0,seed=2026),alpha=.45,size=.7)+
    labs(title="C. Adult measurement predicted from its PGS",x=NULL,y="Held-out reconstruction R²",color="Landmark")+style
  d<-m|>filter(grepl("PGS_captured|PGS_remainder|MatchedMeasured|MatchedPGS",model))
  e<-if(!nrow(d))blank("D. Genetic-score captured / remaining information")else ggplot(d,aes(AUC,short(model),color=factor(landmark)))+
    geom_point()+labs(title="D. Statistical decomposition; not causal partition",x="Out-of-fold IPCW AUC",y=NULL,color="Landmark")+style
  d<-t$cross_omic_links
  f<-if(!nrow(d))blank("E. Protein–metabolite connections")else {
    z<-d|>filter(arm=="all_assays")|>group_by(protein,metabolite)|>summarise(r=mean(r_validation,na.rm=TRUE),.groups="drop")
    ggplot(z,aes(sub("met__","",metabolite),sub("prot__","",protein),fill=r))+geom_tile()+
      scale_fill_gradient2(low="#4374a9",mid="white",high="#ba514e",limits=c(-1,1))+
      labs(title="E. Held-out partial correlations, selected pairs",x=NULL,y=NULL,fill="Mean r")+style+
      theme(axis.text.x=element_text(angle=60,hjust=1,size=6),axis.text.y=element_text(size=6))
  }
  d<-t$marker_heterogeneity
  g<-if(!nrow(d))blank("F. Inflammation-stratum heterogeneity")else ggplot(d|>filter(arm=="all_assays"),aes(difference,paste(definition,short(model)),color=factor(landmark)))+
    geom_vline(xintercept=0,linetype=2)+geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.15)+geom_point()+
    labs(title="F. Difference in improvement: higher minus lower burden",x="Difference in delta AUC",y=NULL,color="Landmark")+style
  save((a|b|c)/(e|f|g)+plot_annotation(title=paste(Y,"— measured proteome, metabolome and matched biomarker PGS"),
    caption="P=proteins; M=metabolites; C=clinical. Disease PRS is separate from biomarker PGS. Automatic UKB-derived PGS are exploratory until discovery overlap is resolved. All contrasts and failures are retained in CSV outputs."),"c5.Fig20.integrated_omics_PGS")
  d<-t$domain_support
  a<-if(!nrow(d))blank("A. Supported domains")else {
    z<-d|>filter(arm=="all_assays")|>group_by(domain,cohort)|>summarise(n=median(eligible),.groups="drop")
    ggplot(z,aes(cohort,domain,fill=log1p(n)))+geom_tile()+geom_text(aes(label=n))+labs(title="A. Training-replicated domain candidates",x=NULL,y=NULL,fill="log(1+n)")+style
  }
  d<-t$proxy_validation
  b<-if(!nrow(d))blank("B. Proxy reconstruction")else ggplot(d|>filter(arm=="all_assays"),aes(modality,R2_vs_training_mean,color=cohort))+
    geom_point(position=position_jitter(width=.12,height=0,seed=2026))+facet_wrap(~target)+labs(title="B. Reconstruct measured domains in held-out participants",x=NULL,y="R² vs training mean")+style
  d<-m|>filter(grepl("Proxy|^Clinical$",model))
  c<-if(!nrow(d))blank("C. Proxy addition versus replacement")else ggplot(d,aes(AUC,short(model),color=factor(landmark)))+geom_point()+
    labs(title="C. Disease benefit is tested separately",x="Out-of-fold IPCW AUC",y=NULL,color="Landmark")+style
  d<-t$panel_counts|>filter(budget==k,arm=="all_assays",model%in%c("Clinical_ProtMet_NS","Clinical_ProtMet_YS_YinYang","Clinical_ProtMet_YSplus_YinYang"))
  e<-if(!nrow(d))blank("D. Actual measured assay requirements")else ggplot(d,aes(short(model),total_assays,color=factor(landmark)))+
    geom_point(position=position_jitter(width=.1,height=0,seed=2026))+labs(title="D. Score count is not assay count",x=NULL,y="Measured assays",color="Landmark")+style+theme(axis.text.x=element_text(angle=20,hjust=1))
  d<-t$coefficients|>filter(budget==k,arm=="all_assays",fold==1,model=="Clinical_ProtMet_YS_YinYang")
  f<-if(!nrow(d))blank("E. Transparent model coefficients")else ggplot(d,aes(beta,variable,color=factor(landmark)))+geom_vline(xintercept=0,linetype=2)+geom_point()+
    labs(title="E. Standardized coefficients, outer fold 1",x="Log-hazard coefficient (training scale)",y=NULL,color="Landmark")+style
  cf<-file.path(root,"c5.cell.coverage.csv");d<-if(file.exists(cf))read.csv(cf)else data.frame()
  g<-blank("F. External cell-type interpretation")
  # Coverage is displayed only after inspecting the annotation output schema.
  if(nrow(d)&&all(c("model","cell_labelled_genes","unique_genes")%in%names(d))){
    d<-d|>filter(grepl(paste0("[|]",k,"[|]all_assays[|]"),model))|>
      mutate(method=sub("[|].*$","",model))|>filter(method%in%c("Clinical_ProtMet_NS","Clinical_ProtMet_YS_YinYang","Clinical_ProtMet_YSplus_YinYang"))
    if(nrow(d))g<-ggplot(d,aes(short(method),cell_labelled_genes/pmax(1,unique_genes)))+
      geom_boxplot()+geom_point(position=position_jitter(width=.1,height=0,seed=2026))+
      labs(title="F. External cell-label coverage across folds",x=NULL,y="Labelled / unique panel genes")+style+theme(axis.text.x=element_text(angle=20,hjust=1))
  }
  save((a|b|c)/(e|f|g)+plot_annotation(title=paste(Y,"— supported LE8 domains and interpretable panels"),
    caption="One to eight candidate domains; unsupported domains contribute no forced proxies. BMI/non-HDL reconstruction does not demonstrate a better clinical measure. No inflammatory-driven subtype is inferred from marker thresholds."),"c5.Fig21.domains_interpretation")
  keep<-c("Clinical","Clinical_ProtMet_NS","Clinical_ProtMet_YS_YinYang")
  d<-t$metrics|>filter(status=="ok",model%in%keep,budget%in%c(0,k))
  a<-if(!nrow(d))blank("A. Lead time and marker omission")else ggplot(d,aes(landmark,AUC,color=short(model),linetype=arm))+
    geom_line()+geom_point()+labs(title="A. Retraining at each landmark",x="Event-free landmark after sampling (years)",y="IPCW AUC",color=NULL,linetype=NULL)+style
  d<-t$calibration
  b<-if(!nrow(d))blank("B. Calibration")else ggplot(d|>filter(model%in%keep,budget%in%c(0,k),arm=="all_assays"),aes(predicted,observed_ipcw,color=short(model)))+
    geom_abline(slope=1,intercept=0,linetype=2)+geom_line()+geom_point()+facet_wrap(~landmark)+
    labs(title="B. Out-of-fold calibration",x="Predicted net risk",y="Observed IPCW net risk",color=NULL)+style
  d<-t$decision_curves
  c<-if(!nrow(d))blank("C. Decision curves")else ggplot(d|>filter(model%in%keep,budget%in%c(0,k),arm=="all_assays"),aes(threshold,net_benefit,color=short(model)))+
    geom_hline(yintercept=0,linetype=2)+geom_line()+facet_wrap(~landmark)+
    labs(title="C. Exploratory net benefit",x="Risk threshold",y="Net benefit",color=NULL)+style
  d<-t$biomarker_PGS_omics_strata
  e<-if(!nrow(d))blank("D. Measured score × biomarker PGS")else ggplot(d|>filter(arm=="all_assays"),aes(stratum,observed_risk,fill=stratum))+
    geom_col()+geom_text(aes(label=paste0("N=",N)),vjust=-.3,size=2.7)+facet_wrap(~landmark)+
    labs(title="D. Complementary measured and matched genetic scores",x=NULL,y="Observed IPCW net risk")+style+
    theme(axis.text.x=element_text(angle=35,hjust=1,size=7),legend.position="none")
  save((a|b)/(c|e)+plot_annotation(title=paste(Y,"— temporal validation and genetic-score strata"),
    caption="Each landmark uses a different eligible risk set and a shorter horizon to the fixed baseline-year endpoint. Cutoffs are learned in training. Death is censored: these displays estimate net risk, not competing-risk incidence."),"c5.Fig22.temporal_calibration_genetic_strata")
}
