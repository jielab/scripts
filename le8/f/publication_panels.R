# All plot data are exported with figure/panel provenance. Missing panels are
# omitted, and each page contains no more than six actual axes.
pub_theme<-function() theme_classic(base_size=11)+theme(
  plot.title=element_text(face='bold',size=12),plot.subtitle=element_text(size=9,color='#526172'),
  axis.text=element_text(color='#26364A'),legend.position='bottom',legend.title=element_blank(),
  panel.grid.major.y=element_line(color='#E9EDF1',linewidth=.25),plot.margin=margin(9,12,9,9))
pub_cols<-c('#286B8B','#BC5366','#6D629A','#C58D37','#388B78','#687784','#9C674E','#6493B0')
pub_ok<-function(d,cols) is.data.frame(d)&&nrow(d)>0&&all(cols%in%names(d))
pub_num<-function(x) suppressWarnings(as.numeric(x))
pub_forest<-function(d,label,beta,lo,hi,title,xlab='Effect estimate',null=0,n=12) {
  if(!pub_ok(d,c(label,beta,lo,hi)))return(NULL)
  d<-d[is.finite(d[[beta]])&is.finite(d[[lo]])&is.finite(d[[hi]]),,drop=FALSE]
  d<-head(d,n);if(!nrow(d))return(NULL)
  d$.label<-factor(d[[label]],levels=rev(unique(d[[label]])))
  ggplot(d,aes(x=.data[[beta]],y=.label))+geom_vline(xintercept=null,color='grey65',linetype=2)+
    geom_segment(aes(x=.data[[lo]],xend=.data[[hi]],yend=.label),color=pub_cols[1])+
    geom_point(color=pub_cols[1],size=2.2)+labs(title=title,x=xlab,y=NULL)+pub_theme()
}
pub_count<-function(d,col,title) {
  if(!pub_ok(d,col))return(NULL)
  z<-d|>count(.data[[col]],name='n');z$.label<-stringr::str_wrap(as.character(z[[col]]),32)
  ggplot(z,aes(n,reorder(.label,n)))+geom_col(fill=pub_cols[1],width=.68)+
    geom_text(aes(label=n),hjust=-.15,size=3)+scale_x_continuous(expand=expansion(mult=c(0,.15)))+
    labs(title=title,x='Biomarkers',y=NULL)+pub_theme()
}
pub_enrich<-function(d,title,asset_dir) {
  if(!pub_ok(d,c('term_name','adjusted_p','source')))return(NULL)
  d<-d|>filter(is.finite(adjusted_p),adjusted_p<.05)|>arrange(adjusted_p)|>
    group_by(source)|>slice_head(n=3)|>ungroup()|>arrange(adjusted_p)|>slice_head(n=10)
  if(!nrow(d))return(NULL)
  f<-file.path(asset_dir,'go_terms.tsv')
  if(file.exists(f)){labs<-fread(f);i<-match(d$term_name,labs$term);hit<-!is.na(i);d$term_name[hit]<-labs$TERM[i[hit]]}
  d$.label<-factor(stringr::str_wrap(d$term_name,40),levels=rev(unique(stringr::str_wrap(d$term_name,40))))
  ggplot(d,aes(-log10(pmax(adjusted_p,1e-300)),.label,color=source))+geom_point(size=3)+
    scale_color_manual(values=pub_cols)+labs(title=title,x=expression(-log[10](FDR)),y=NULL)+pub_theme()
}
publication_run<-function(trait,layer,analysis_root) {
  stopifnot(layer%in%c('prot','met'))
  root<-file.path(analysis_root,trait,layer)
  if(!dir.exists(root)){message('FINAL: no analysis directory: ',root);return(invisible(NULL))}
  modules<-c(c1='c1_correlate',c2='c2_cause',c3='c3_coloc',c4='c4_connect',focus='c4_focus',c5='c5_consolidate',extra='final_supp')
  audit<-list();tables<-list();provenance<-list();figures<-list();captions<-character()
  read<-function(key,module,file) {
    if(layer=='met')file<-sub('^pwas_','mwas_',file)
    path<-file.path(root,modules[[module]],file)
    d<-if(file.exists(path))tryCatch(as_tibble(fread(path,showProgress=FALSE)),error=function(e)tibble())else tibble()
    audit[[key]]<<-tibble(table=key,source=path,rows=nrow(d),bytes=if(file.exists(path))file.info(path)$size else NA_real_,
      modified=if(file.exists(path))as.character(file.info(path)$mtime)else NA_character_,status=if(nrow(d))'available'else'absent or no rows')
    tables[[key]]<<-d;d
  }
  cohort<-read('cohort','c1','c1.cohort.csv')
  a<-read('incident','c1','pwas_incident_adj2.csv');prev<-read('prevalent','c1','pwas_prevalent_adj2.csv')
  temporal<-read('temporal','c1','c1.directionality_triage.csv')
  pgs<-read('pgs','c1','c1.pgs_actual_concordance.csv')
  paired<-read('paired','c1','c1.paired_pgs_measured.csv')
  enrich<-read('enrich','c1','c1.enrichment_incident_sig.csv')
  enrich_prev<-read('enrich_prev','c1','c1.enrichment_prevalent_sig.csv')
  mr<-read('mr','c2','c2.MR_all.csv');grades<-read('grades','c2','c2.evidence_grades.csv')
  dan<-read('dandelion_audit','c2','c2.dandelion_input_audit.csv')
  co<-read('coloc','c3','c3.coloc_summary.csv')
  # Reuse the project's strict overlap helper without sourcing any model worker.
  helper<-parse(file.path(Sys.getenv('LE8_FDIR',unset='.'),'c0_revision_core.R'))
  for(ex in helper)if(is.call(ex)&&identical(ex[[1]],as.name('<-'))&&identical(ex[[2]],as.name('le8_same_locus_evidence')))eval(ex)
  same<-le8_same_locus_evidence(mr,co,if(layer=='prot')'protein'else'metabolite')
  tables$same_region<-same
  membership<-read('membership','c4','c4.proxy_membership_YS_YSP_NS.csv')
  med<-read('mediation','c4','c4.mediation_all.csv')
  bridge<-read('bridges','c4','c4.genetic_omic_disease_bridges.csv')
  matched_pgs<-read('matched_pgs','c4','c4.matched_PGS_bridges.csv')
  cigma<-read('cigma_status','c3','c3.CIGMA_status.csv')
  read('cigma_annotation','c3','c3.CIGMA_annotation.csv')
  fm<-read('focus_metrics','focus','c4.focus.metrics.csv')
  fc<-read('focus_contrasts','focus','c4.focus.contrasts.csv')
  fp<-read('focus_proxy_accuracy','focus','c4.focus.proxy_accuracy.csv')
  fpc<-read('focus_pillars','focus','c4.focus.pillar_counts.csv')
  read('focus_design','focus','c4.focus.design.csv')
  read('focus_members','focus','c4.focus.panel_members.csv')
  imaging<-read('imaging_status','c4','c4.imaging_status.csv')
  perf<-read('performance','c5','c5.prediction_summary.csv')
  budget<-read('budget','c5','c5.review_budget_metrics.csv')
  cal<-read('calibration','c5','c5.review_calibration.csv')
  delta<-read('paired_delta','c5','c5.review_paired_delta_CI.csv')
  design<-read('prediction_design','c5','c5.review_design.csv')
  lead<-read('leadtime','c5','c5.leadtime_discrimination.csv')
  ev<-read('legacy_evidence','c5','c5.evidence_consolidation.csv')
  panel_members<-read('panel_members','c5','c5.review_panel_members.csv')
  assets<-file.path(Sys.getenv('LE8_FDIR',unset='.'),'assets')
  add<-function(fig,p,keys,label) {
    if(is.null(p))return(invisible(NULL))
    figures[[fig]][[length(figures[[fig]])+1L]]<<-p
    provenance[[length(provenance)+1L]]<<-tibble(figure=fig,panel=LETTERS[length(figures[[fig]])],
      title=label,tables=paste(keys,collapse=';'),selection='Prespecified panel recipe; display ranking only, not prediction selection')
  }
  # Figure 1: four panels; counts, volcano, anchors and prevalent concordance.
  if(pub_ok(a,c('term','beta','p.value','estimate','conf.low','conf.high'))) {
    z<-a|>filter(is.finite(beta),is.finite(p.value))|>mutate(support=ifelse(p.value<.05/nrow(a),'Bonferroni','Other'))
    p<-ggplot(z,aes(beta,-log10(pmax(p.value,1e-300)),color=support))+geom_point(alpha=.65,size=1.1)+
      geom_hline(yintercept=-log10(.05/nrow(a)),linetype=2,color='grey55')+
      scale_color_manual(values=c(Bonferroni=pub_cols[2],Other='#B6C1CC'))+
      labs(title='Incident association scan',x='Log hazard ratio per SD',y=expression(-log[10](P)))+pub_theme()
    add('Fig1',p,'incident','Incident association scan')
    counts<-tibble(analysis=c('Incident: tested','Incident: Bonferroni','Prevalent: tested','Prevalent: Bonferroni'),
      n=c(nrow(a),sum(a$p.value<.05/nrow(a),na.rm=TRUE),nrow(prev),if(nrow(prev))sum(prev$p.value<.05/nrow(prev),na.rm=TRUE)else 0))
    tables$association_counts<-counts
    add('Fig1',ggplot(counts,aes(n,reorder(analysis,n)))+geom_col(fill=pub_cols[1],width=.65)+geom_text(aes(label=n),hjust=-.12,size=3.5)+
      scale_x_continuous(expand=expansion(mult=c(0,.18)))+labs(title='Association yield',x='Biomarkers',y=NULL)+pub_theme(),'association_counts','Association yield')
    aa<-a|>arrange(p.value)
    add('Fig1',pub_forest(aa,'term','estimate','conf.low','conf.high','Leading incident associations','Hazard ratio per SD',1), 'incident','Leading incident associations')
    if(pub_ok(prev,c('term','beta'))) {
      z<-inner_join(a|>select(term,incident=beta),prev|>select(term,prevalent=beta),by='term')
      add('Fig1',ggplot(z,aes(incident,prevalent))+geom_hline(yintercept=0,color='grey85')+geom_vline(xintercept=0,color='grey85')+
        geom_point(color=pub_cols[1],alpha=.5,size=1)+labs(title='Incident versus prevalent associations',x='Incident log HR',y='Prevalent log OR')+pub_theme(),c('incident','prevalent'),'Incident versus prevalent')
    }
  }
  captions['Fig1']<-'Baseline prevalent cases contribute to descriptive and prevalent analyses; incident models exclude baseline disease. Bonferroni correction is within each omics scan. Prevalent odds ratios and incident hazard ratios estimate different quantities.'
  # Figure 2: inherited propensity, cis/local MR, colocalization, and diagnostic coverage.
  if(pub_ok(paired,c('feature','model','component','beta','std.error'))) {
    z<-paired|>filter(grepl('^separate;',model),is.finite(beta),is.finite(std.error))
    shown<-head(unique(z$feature),10)
    z<-z|>filter(feature%in%shown)|>mutate(feature=factor(feature,levels=rev(shown)),lo=beta-1.96*std.error,hi=beta+1.96*std.error)
    if(nrow(z))add('Fig2',ggplot(z,aes(beta,feature,color=component))+geom_vline(xintercept=0,color='grey70',linetype=2)+
      geom_errorbar(aes(xmin=lo,xmax=hi),orientation='y',width=.2,position=position_dodge(width=.5))+
      geom_point(position=position_dodge(width=.5))+scale_color_manual(values=pub_cols)+
      labs(title='Paired measured and inherited propensity',subtitle='Identical participants and covariates within each biomarker',x='Log HR per own SD (95% CI)',y=NULL,color=NULL)+pub_theme(),'paired','Paired measured and inherited propensity')
  }
  primary<-if(layer=='prot')'cis'else'local'
  if(pub_ok(mr,c('analysis','b','se','pval','exposure'))) {
    mm<-mr|>filter(analysis==primary)|>arrange(pval)|>mutate(lo=b-1.96*se,hi=b+1.96*se)
    add('Fig2',pub_forest(mm,'exposure','b','lo','hi',paste(primary,'MR estimates'),'Effect on GWAS outcome scale'),'mr','Primary MR estimates')
  }
  if(pub_ok(co,c('PP.H4','PP.H4_robust_min','locus_class','status'))) {
    z<-co|>filter(status=='ok',is.finite(PP.H4),is.finite(PP.H4_robust_min))
    if(nrow(z))add('Fig2',ggplot(z,aes(PP.H4,PP.H4_robust_min,color=locus_class))+geom_point(alpha=.65)+
      geom_hline(yintercept=.7,linetype=2)+scale_color_manual(values=pub_cols)+
      labs(title='Colocalization and prior sensitivity',x='PP.H4, default prior',y='Minimum PP.H4 across priors')+pub_theme(),'coloc','Colocalization prior sensitivity')
  }
  add('Fig2',pub_count(grades,'evidence_grade','MR diagnostic coverage'),'grades','MR diagnostic coverage')
  captions['Fig2']<-'PGS is an inherited genetic score, not a protein concentration measured at birth. Paired measured and PGS estimates use identical people and covariates, standardized to their own SD; their units are not interchangeable. Cis/local MR is primary; trans/distal associations remain exploratory. Colocalization supports a shared region, not proof of causal direction.'
  # Figure 3: four panels capturing LE8 supervision, reproducibility and mediation.
  if(pub_ok(membership,c('in_YS','in_YSP_plus','in_NS','primary_component','r_disc','r_rep'))) {
    z<-membership|>filter(in_YS%in%TRUE)
    pillars<-tibble(primary_component=c('diet','pa','smoke','sleep','bmi','nonhdl','hba1c','bp'))|>left_join(z|>count(primary_component),by='primary_component')|>mutate(n=coalesce(n,0L))
    tables$pillar_counts<-pillars
    add('Fig3',ggplot(pillars,aes(factor(primary_component,levels=primary_component),n))+geom_col(fill=pub_cols[1],width=.65)+geom_text(aes(label=n),vjust=-.4)+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(title='LE8 pillars represented among YS proxies',x=NULL,y='Biomarkers')+pub_theme(),'pillar_counts','LE8 pillar coverage')
    counts<-tibble(set=c('YS','YSP plus','NS'),n=c(sum(membership$in_YS%in%TRUE),sum(membership$in_YSP_plus%in%TRUE),sum(membership$in_NS%in%TRUE)))
    tables$set_counts<-counts
    add('Fig3',ggplot(counts,aes(set,n,fill=set))+geom_col(width=.6,show.legend=FALSE)+geom_text(aes(label=n),vjust=-.4)+
      scale_fill_manual(values=pub_cols)+scale_y_continuous(expand=expansion(mult=c(0,.15)))+labs(title='Proxy and comparator sets',x=NULL,y='Biomarkers (sets may overlap)')+pub_theme(),'set_counts','Proxy and comparator sets')
    if(nrow(z))add('Fig3',ggplot(z,aes(r_disc,r_rep,color=primary_component))+geom_abline(slope=1,intercept=0,color='grey65',linetype=2)+
      geom_point(alpha=.6,size=1.4)+scale_color_manual(values=pub_cols)+labs(title='Discovery and replication profiles',x='Discovery partial correlation',y='Replication partial correlation')+pub_theme(),'membership','Replication profiles')
  }
  if(pub_ok(med,c('component','feature','indirect_beta','indirect_lo','indirect_hi','FDR_indirect'))) {
    z<-med|>arrange(FDR_indirect)|>mutate(label=paste(component,feature,sep=' / '))
    add('Fig3',pub_forest(z,'label','indirect_beta','indirect_lo','indirect_hi','Exploratory mediation paths','Product-of-coefficients indirect effect',0,10),'mediation','Mediation paths')
  }
  captions['Fig3']<-'LE8 → omics → disease ← omics ← genetics is the organizing hypothesis. YS and YSP-plus are supervised and additional proxies; NS is an independently defined unsupervised comparator and may overlap them. Baseline mediation is exploratory and does not establish intervention effects. An empty genetic bridge table means that connection was not demonstrated.'
  # Figure 4: held-out discrimination, assay budgets, calibration, paired uncertainty.
  if(pub_ok(perf,c('biom_set','model','AUC','C_index'))) {
    z<-perf|>filter(is.finite(AUC))
    add('Fig4',ggplot(z,aes(AUC,biom_set,color=model))+geom_point(position=position_dodge(width=.5),size=2)+
      scale_color_manual(values=pub_cols)+scale_y_discrete(labels=function(x)stringr::str_wrap(x,25))+
      labs(title='Held-out discrimination',x='IPCW AUC (recorded horizon)',y=NULL)+pub_theme(),'performance','Held-out discrimination')
  }
  if(pub_ok(budget,c('ablation','horizon','actual_n_assays','AUC','paradigm','status'))&&all(LE8_ASSAY_BUDGETS%in%budget$actual_n_assays)) {
    z<-budget|>filter(ablation=='none',horizon==10,status=='ok')
    if(nrow(z))add('Fig4',ggplot(z,aes(actual_n_assays,AUC,color=paradigm,group=paradigm))+geom_line(na.rm=TRUE)+geom_point(size=2)+
      scale_color_manual(values=pub_cols)+labs(title='Performance at explicit assay budgets',x='Actual fitted assays',y='10-year IPCW AUC')+pub_theme(),'budget','Assay budget comparison')
  }
  if(pub_ok(cal,c('horizon','ablation','budget','predicted','observed_ipcw','model'))) {
    z<-cal|>filter(horizon==10,ablation=='none',budget%in%c(0,10))
    if(nrow(z))add('Fig4',ggplot(z,aes(predicted,observed_ipcw,color=model))+geom_abline(slope=1,intercept=0,linetype=2,color='grey65')+
      geom_line()+geom_point(size=1.4)+scale_color_manual(values=pub_cols)+labs(title='Calibration at a fixed 10-assay budget',x='Predicted 10-year risk',y='Observed IPCW risk')+pub_theme(),'calibration','Calibration at 10 assays')
  }
  if(pub_ok(delta,c('horizon','model','delta_AUC_lo','delta_AUC_hi'))) {
    z<-delta|>filter(horizon==10,is.finite(delta_AUC_lo),is.finite(delta_AUC_hi),
      !grepl('_[0-9]+$',model)|pub_num(sub('.*_','',model))%in%LE8_ASSAY_BUDGETS)|>slice_head(n=12)
    if(nrow(z))add('Fig4',ggplot(z,aes(y=reorder(model,delta_AUC_lo)))+geom_vline(xintercept=0,linetype=2,color='grey55')+
      geom_segment(aes(x=delta_AUC_lo,xend=delta_AUC_hi,yend=reorder(model,delta_AUC_lo)),linewidth=1.2,color=pub_cols[1])+
      labs(title='Paired improvement versus clinical model',x='95% bootstrap interval for ΔAUC (10 years)',y=NULL)+pub_theme(),'paired_delta','Paired AUC intervals')
  }
  captions['Fig4']<-'Frozen models evaluated on the held-out sample. The clinical comparator is defined in the design sheet; it is not automatically SCORE2. Assay budgets and 10-year horizon are prespecified in this presentation. Bootstrap intervals describe frozen fits, not the uncertainty of the entire model-selection process. Deaths are censored; these risks are not competing-risk cumulative incidences.'
  # Prefer the direct test of the article's hypothesis once it exists, including
  # null/adverse results. Do not choose panels by significance or best AUC.
  if(pub_ok(fpc,c('cohort','component','n'))&&pub_ok(fp,c('component','model','delta_R2'))&&length(figures[['Fig3']])>=4) {
    figures[['Fig3']][[1]]<-ggplot(fpc,aes(component,n,fill=cohort))+geom_col(position='dodge')+
      labs(title='LE8 proxy discovery: Yin and Yin + Yang',x=NULL,y='Replicated proxies')+pub_theme()+theme(axis.text.x=element_text(angle=35,hjust=1))
    display_budget<-max(fm$budget,na.rm=TRUE)
    z<-fp|>filter(grepl(paste0('_',display_budget,'$'),model))
    figures[['Fig3']][[2]]<-ggplot(z,aes(component,model,fill=delta_R2))+geom_tile()+
      scale_fill_gradient2(low=pub_cols[2],mid='white',high=pub_cols[1],midpoint=0)+
      labs(title=paste('Held-out LE8 reconstruction at',display_budget,'assays'),x=NULL,y=NULL,fill='Incremental R²')+pub_theme()+theme(axis.text.x=element_text(angle=35,hjust=1))
    for(i in seq_along(provenance))if(provenance[[i]]$figure=='Fig3'&&provenance[[i]]$panel%in%c('A','B')) {
      provenance[[i]]$tables<-if(provenance[[i]]$panel=='A')'focus_pillars'else'focus_proxy_accuracy'
      provenance[[i]]$title<-if(provenance[[i]]$panel=='A')'Yin/Yang proxy discovery'else'Held-out LE8 reconstruction'
    }
    captions['Fig3']<-paste(captions['Fig3'],'Panels A/B use only outer-training donors for discovery and fit; Yang donors contribute to proxy learning. Incremental R² is measured beyond basic covariates on held-out incident participants, and is not intervention responsiveness.')
  }
  if(pub_ok(fm,c('stratum','landmark','budget','AUC','paradigm'))&&pub_ok(fc,c('stratum','landmark','model','delta_AUC','delta_lo','delta_hi'))) {
    figures[['Fig4']]<-list();provenance<-Filter(function(x)x$figure!='Fig4',provenance)
    z<-fm|>filter(stratum=='All',landmark==0,budget>0)
    clinical_auc<-fm$AUC[fm$model=='Clinical'&fm$stratum=='All'&fm$landmark==0][1]
    add('Fig4',ggplot(z,aes(budget,AUC,color=paradigm))+geom_line()+geom_point()+
      geom_hline(yintercept=clinical_auc,linetype=2,color='grey50')+
      scale_x_continuous(breaks=sort(unique(z$budget)))+
      labs(title='Matched assay budgets',subtitle=sprintf('Dashed line: clinical covariates alone (AUC %.3f)',clinical_auc),x='Measured assays',y='10-year IPCW AUC')+pub_theme(),'focus_metrics','Matched assay budgets')
    z<-fc|>filter(stratum=='All',landmark==0)
    add('Fig4',ggplot(z,aes(delta_AUC,reorder(model,delta_AUC),color=contrast))+geom_vline(xintercept=0,linetype=2)+
      geom_errorbar(aes(xmin=delta_lo,xmax=delta_hi),orientation='y',width=.2)+geom_point()+
      labs(title='Supervision and added-Yang contrasts',x='Paired ΔAUC with 95% bootstrap interval',y=NULL)+pub_theme(),'focus_contrasts','Supervision and Yang contrasts')
    z<-fm|>filter(stratum=='All',budget%in%c(0,max(fm$budget,na.rm=TRUE)))
    add('Fig4',ggplot(z,aes(landmark,AUC,color=paradigm))+geom_line()+geom_point()+
      labs(title='Conditional prediction to baseline year 10',x='Disease-free landmark (years)',y='IPCW AUC after landmark')+pub_theme(),'focus_metrics','Landmark validation')
    z<-fm|>filter(stratum!='All',landmark==0,budget%in%c(0,max(fm$budget,na.rm=TRUE)))
    if(nrow(z))add('Fig4',ggplot(z,aes(AUC,paradigm,color=stratum))+geom_point(position=position_dodge(width=.4))+
      labs(title='Baseline inflammation strata',x='10-year IPCW AUC',y=NULL)+pub_theme(),'focus_metrics','Baseline inflammation strata')
    else add('Fig4',pub_count(tables$focus_members,'model','Assays in each evaluated panel'),'focus_members','Panel composition')
    captions['Fig4']<-'Same incident test cohort, same assay budget, same Cox model and clinical covariates. Yang donors enter only proxy discovery. YSplus reserves 80% of slots (rounded up) for YS by default; see focus_design. Low baseline inflammation is an outcome-independent marker stratum, not a validated CAD subtype. These exploratory comparisons were designed after inspecting prior results; bootstrap intervals condition on frozen fits and require external confirmation. Death is censored.'
  }
  if(pub_ok(matched_pgs,c('feature','r_disc','r_rep','PGS_disease_FDR'))&&length(figures[['Fig3']])>=3) {
    disease_label<-paste0('PGS-',trait,' FDR < 0.05')
    z<-matched_pgs|>filter(is.finite(r_disc),is.finite(r_rep))|>
      mutate(disease_support=case_when(!is.finite(PGS_disease_FDR)~'Unavailable',PGS_disease_FDR<.05~disease_label,TRUE~'Not detected'))
    figures[['Fig3']][[3]]<-ggplot(z,aes(r_disc,r_rep,color=disease_support))+
      geom_abline(slope=1,intercept=0,linetype=2,color='grey70')+geom_point(alpha=.55,size=1.2)+
      scale_color_manual(values=pub_cols)+labs(title='Matched inherited propensity → measured omic',
        subtitle='Each biomarker is paired with its own PGS',x='Discovery partial correlation',y='Replication partial correlation',color=NULL)+pub_theme()
    for(i in seq_along(provenance))if(provenance[[i]]$figure=='Fig3'&&provenance[[i]]$panel=='C') {
      provenance[[i]]$tables<-'matched_pgs';provenance[[i]]$title<-'Matched PGS calibration'
    }
    captions['Fig3']<-sub('An empty genetic bridge table means that connection was not demonstrated.','',captions['Fig3'],fixed=TRUE)
    captions['Fig3']<-paste(captions['Fig3'],'Panel C uses two training halves. Source-GWAS overlap may inflate PGS calibration; shared associations do not identify a causal mediation chain.')
  }
  # Conservative evidence matrix built afresh, never endorses legacy C5 grades.
  if(pub_ok(a,c('term','p.value','FDR'))) {
    strict<-a|>transmute(feature=term,observed=ifelse(is.finite(FDR),FDR<.05,NA),p=p.value)
    pg<-if(pub_ok(pgs,c('analysis','feature','pgs_FDR'))) pgs|>filter(analysis=='Incident')|>distinct(feature,.keep_all=TRUE)|>transmute(feature,PGS=ifelse(is.finite(pgs_FDR),pgs_FDR<.05,NA)) else tibble(feature=character(),PGS=logical())
    me<-if(pub_ok(mr,c('analysis','exposure','FDR_all')))mr|>filter(analysis==primary)|>group_by(exposure)|>summarise(MR=if(all(!is.finite(FDR_all)))NA else any(FDR_all<.05,na.rm=TRUE),.groups='drop')|>rename(feature=exposure)else tibble(feature=character(),MR=logical())
    ys<-if(pub_ok(membership,c('feature','in_YS')))membership|>distinct(feature,.keep_all=TRUE)|>transmute(feature,YS=in_YS%in%TRUE)else tibble(feature=character(),YS=logical())
    sl<-if(nrow(same))same|>group_by(feature)|>summarise(same_region=any(eligible%in%TRUE),.groups='drop')else tibble(feature=character(),same_region=logical())
    strict<-strict|>left_join(pg,by='feature')|>left_join(me,by='feature')|>left_join(sl,by='feature')|>left_join(ys,by='feature')|>
      mutate(role=case_when(same_region%in%TRUE~'Cis/local region-supported candidate',MR%in%TRUE~'MR support; locus unresolved',YS%in%TRUE~'LE8 proxy',TRUE~'Association / unresolved'))|>arrange(desc(same_region%in%TRUE),p)
    tables$publication_evidence<-strict
    top<-head(strict,18);hm<-top|>pivot_longer(c(observed,PGS,MR,same_region,YS),names_to='domain',values_to='support')|>
      mutate(domain=factor(domain,levels=c('observed','PGS','MR','same_region','YS'),labels=c('Observed','PGS','MR','Same region','YS')),feature=factor(feature,levels=rev(top$feature)),state=case_when(is.na(support)~'Unavailable',support~'Supported',TRUE~'Not detected'))
    add('Fig5',ggplot(hm,aes(domain,feature,fill=state))+geom_tile(color='white',linewidth=.5)+
      scale_fill_manual(values=c(Supported=pub_cols[1],`Not detected`='#E5EAF0',Unavailable='#F4D9B3'))+
      labs(title='Conservative cross-domain evidence',x=NULL,y=NULL)+pub_theme(),c('publication_evidence','same_region'),'Cross-domain evidence')
    add('Fig5',pub_count(strict,'role','Evidence roles, without causal proof'),'publication_evidence','Evidence roles')
  }
  add('Fig5',pub_enrich(enrich,'Functional context of incident associations',assets),'enrich','Incident enrichment')
  if(pub_ok(lead,c('horizon','AUC','kind','method'))) {
    z<-lead|>filter(kind=='Omic score',method%in%c('C4 YS','C4 YSplus','C4 NS','Pradeep-style / glmnet','Yu-style / LightGBM'),is.finite(AUC))
    if(nrow(z))add('Fig5',ggplot(z,aes(horizon,AUC,color=method))+geom_line()+geom_point(size=1.5)+
      scale_color_manual(values=pub_cols)+labs(title='Discrimination with minimum lead time',x='Minimum years before diagnosis',y='Case/control AUC')+pub_theme(),'leadtime','Minimum lead-time discrimination')
  }
  captions['Fig5']<-'Cis/local MR and conservative colocalization must overlap retained instruments in the same region. This is region-level corroboration, not signal-resolved causality; missing tests are distinct from negative tests. The minimum-lead-time case/control AUC is not the IPCW AUC in Fig4. Annotation enrichment is contextual evidence.'
  # The generic association landscape is retained as a supplement. Main
  # Figure 1 now joins C1 matched scores, calibrated components and C3 loci.
  old_fig1<-figures[['Fig1']]
  if(length(old_fig1))pgs_save_panels(old_fig1,file.path(root,'final_supp'),
    'final.association_landscape',tables[c('incident','prevalent','cohort')],
    paste(trait,toupper(layer),'association context'),captions[['Fig1']])
  figures[['Fig1']]<-list()
  provenance<-Filter(function(z)!identical(z$figure,'Fig1'),provenance)
  pf<-file.path(root,'c1_correlate','c1.pgs_focus.rds')
  focus_pgs<-if(file.exists(pf))tryCatch(readRDS(pf),error=function(e)list())else list()
  if(pgs_ok(focus_pgs$paired,c('feature','evidence_pattern'))) {
    triangulation<-make_c3_pgs_integration(coloc_summary=co,focus=focus_pgs)
    loci<-attr(triangulation,'loci')
    panels<-pgs_main_panels(focus_pgs,triangulation,loci)
    for(nm in names(focus_pgs))if(is.data.frame(focus_pgs[[nm]]))tables[[paste0('PGS_',nm)]]<-focus_pgs[[nm]]
    tables$PGS_loci<-loci;tables$PGS_triangulation<-as.data.frame(triangulation)
    keys<-grep('^PGS_',names(tables),value=TRUE)
    for(nm in names(panels))add('Fig1',panels[[nm]],keys,paste('PGS discordance:',nm))
    audit$PGS_focus<-tibble(table='PGS_focus',source=pf,rows=nrow(focus_pgs$paired),
      bytes=file.info(pf)$size,modified=as.character(file.info(pf)$mtime),status='available')
  }else audit$PGS_focus<-tibble(table='PGS_focus',source=pf,rows=0L,bytes=NA_real_,modified=NA_character_,
    status='Run pgs_focus; legacy unequal-sample directions are not used as main evidence')
  captions['Fig1']<-'Identical-sample measured/PGS associations and cross-fitted captured (G) versus remaining (R) components. Joint contrasts include covariance; conditional CIs and refitted-bootstrap intervals are distinguished. R is not a pure lifestyle fraction. Colocalization is locus-specific, with prior sensitivity retained. Candidate selection is exploratory; opposite directions do not prove antagonistic pleiotropy.'
  manifest<-if(length(provenance))bind_rows(provenance)else tibble(figure=character(),panel=character(),title=character(),tables=character(),selection=character())
  manifest$analysis<-manifest$figure
  source_audit<-bind_rows(audit)
  # Preserve old main/supplement versions before assembling the new edition.
  old<-list.files(root,pattern='^(Fig[0-9]+|FigS[0-9]+)\\.(png|xlsx)$',full.names=TRUE)
  le8_archive_figure_files(old,root)
  status<-list();figure_number<-0L;rendered_manifest<-list()
  dpi<-as.integer(Sys.getenv('LE8_FINAL_DPI',unset='320'))
  for(fig in paste0('Fig',1:5)) {
    expanded<-lapply(figures[[fig]],le8_expand_facets)
    pp<-unlist(expanded,recursive=FALSE);n<-length(pp)
    if(!n){status[[fig]]<-tibble(figure=NA_character_,analysis=fig,panels=0,status='withheld: no usable panels');next}
    keys<-unique(unlist(strsplit(manifest$tables[manifest$analysis==fig],';',fixed=TRUE)))
    source_rows<-manifest[manifest$analysis==fig,]
    panel_rows<-map_dfr(seq_along(expanded),function(i) {
      if(!length(expanded[[i]]))return(tibble())
      map_dfr(expanded[[i]],function(p)source_rows[i,]|>mutate(source_panel=panel,
        title=paste(as.character(p$labels$title),collapse=' ')))
    })
    if(n>6L) {
      pgs_save_panels(pp[7:n],file.path(root,'final_supp'),paste0('final.overflow.',fig),
        tables[keys],paste(trait,toupper(layer),fig,'additional panels'),captions[[fig]])
      pp<-pp[1:6];panel_rows<-panel_rows[1:6,];n<-6L
    }
    actual<-character()
    for(start in seq(1L,n,by=6L)) {
      page<-pp[start:min(n,start+5L)];figure_number<-figure_number+1L;dest<-fig;actual<-c(actual,dest)
      page_rows<-panel_rows[start:min(n,start+5L),]|>mutate(figure=dest,panel=LETTERS[seq_along(page)])
      rendered_manifest[[length(rendered_manifest)+1L]]<-page_rows
      sheets<-c(list(provenance=page_rows,caption=tibble(caption=captions[[fig]])),tables[keys])
      p<-wrap_plots(page,ncol=if(length(page)==1)1 else 2,widths=c(1,1))+plot_annotation(title=paste(trait,toupper(layer),'|',switch(fig,Fig1='Measured versus genetically predicted biomarkers',Fig2='Genetic triangulation',Fig3='LE8 connections',Fig4='Prediction and assay budget',Fig5='Evidence consolidation')),
        caption=stringr::str_wrap(gsub('Fig4','the prediction figure',captions[[fig]],fixed=TRUE),155),tag_levels='A',theme=theme(plot.title=element_text(size=17,face='bold'),plot.caption=element_text(size=9,hjust=0),plot.tag=element_text(face='bold')))
      ggsave(file.path(root,paste0(dest,'.png')),p,width=19,height=ceiling(length(page)/2)*6+.9,dpi=dpi,bg='white',limitsize=FALSE)
      pub_workbook(sheets,file.path(root,paste0(dest,'.xlsx')))
    }
    status[[fig]]<-tibble(figure=paste(actual,collapse=';'),analysis=fig,panels=n,status='rendered')
  }
  manifest<-if(length(rendered_manifest))bind_rows(rendered_manifest)else
    tibble(figure=character(),panel=character(),title=character(),tables=character(),selection=character(),analysis=character(),source_panel=character())
  # Retain every available module figure in final supplements. Manuscript
  # selection must not silently delete a view from the reviewable output set.
  # No raster tiling or stretching: preserve original resolution and aspect ratio.
  supp<-list();sn<-0L
  for(code in names(modules)) {
    module<-modules[[code]];mf<-file.path(root,module,'figure_manifest.csv')
    if(!file.exists(mf))next
    fm<-as_tibble(fread(mf));files<-file.path(root,module,fm$file[fm$panels<=6])
    for(src in files[file.exists(files)]) {
    sn<-sn+1L;dest<-paste0('FigS',sn,'.png')
    if(!file.copy(src,file.path(root,dest),overwrite=TRUE))stop('Cannot copy supplement: ',src)
    sx<-sub('[.]png$','.xlsx',src)
    if(file.exists(sx)&&!file.copy(sx,file.path(root,sub('[.]png$','.xlsx',dest)),overwrite=TRUE))stop('Cannot copy supplementary workbook: ',sx)
    supp[[sn]]<-tibble(figure=dest,source=src,module=module,review='Retained module figure; original source layout and previous editions are preserved')
    }
  }
  pub_workbook(c(list(sources=source_audit,figure_status=bind_rows(status),panels=manifest,supplements=bind_rows(supp)),
    tables[c('cohort','prediction_design','dandelion_audit','imaging_status','bridges','matched_pgs','cigma_status','cigma_annotation','focus_metrics','focus_contrasts','focus_proxy_accuracy','focus_pillars','focus_design','focus_members','paired','same_region','publication_evidence','panel_members','enrich_prev')]),file.path(root,'TablesS.xlsx'))
  fwrite(manifest,file.path(root,'final.figure_manifest.csv'))
  supp_table<-if(length(supp))bind_rows(supp)else tibble(figure=character(),source=character(),module=character(),review=character())
  fwrite(supp_table,file.path(root,'final.supplement_manifest.csv'))
  fwrite(bind_rows(status),file.path(root,'final.status.csv'))
  pub_report(root,trait,layer,tables,source_audit,bind_rows(status),captions)
  message('FINAL: ',root,'; ',sum(vapply(status,function(z)z$status=='rendered',logical(1))),' main figures; ',sn,' supplements')
  invisible(status)
}
pub_workbook<-function(sheets,file) {
  wb<-createWorkbook();sheets<-sheets[!vapply(sheets,is.null,logical(1))]
  names(sheets)<-make.unique(substr(names(sheets),1,28))
  for(nm in names(sheets)) {
    z<-as.data.frame(sheets[[nm]]);if(!ncol(z))z<-data.frame(status='No available result')
    addWorksheet(wb,nm);if(nrow(z))writeDataTable(wb,nm,z,tableStyle='TableStyleMedium2')else writeData(wb,nm,z)
    freezePane(wb,nm,firstRow=TRUE);setColWidths(wb,nm,cols=seq_len(ncol(z)),widths=18)
  }
  saveWorkbook(wb,file,overwrite=TRUE)
}
pub_report<-function(root,trait,layer,t,source_audit,status,captions) {
  lines<-c(paste('#',trait,'/',layer,'— 结果与写作清单'),'',paste('生成时间：',Sys.time()),'',
    '本摘要读取模块内的聚合结果；缺失结果不会写成阴性发现。论文数值以本次输出的来源表和实际运行设计为准。','')
  if(pub_ok(t$cohort,c('N_omics','incident_events','prevalent_cases','features'))) {
    d<-t$cohort[1,];lines<-c(lines,sprintf('Omics cohort：%s 人；incident events %s；baseline prevalent cases %s；%s 个 biomarker。各模型 complete-case N 需单独报告。',d$N_omics,d$incident_events,d$prevalent_cases,d$features),'')
  }
  for(key in c('incident','prevalent')) {
    d<-t[[key]];if(!pub_ok(d,c('p.value','FDR')))next
    lines<-c(lines,sprintf('%s：测试 %d 项，Bonferroni 显著 %d 项，BH-FDR <0.05 为 %d 项。',key,nrow(d),sum(d$p.value<.05/nrow(d),na.rm=TRUE),sum(d$FDR<.05,na.rm=TRUE)))
  }
  if(pub_ok(t$PGS_paired,c('evidence_pattern','feature'))) {
    d<-t$PGS_paired;opposite<-d$feature[d$evidence_pattern=='Both supported: opposite']
    lines<-c(lines,'',sprintf('同样本 measured/PGS：%d 项可匹配指标；两边均 BH-FDR <0.05 且方向相反 %d 项。',nrow(d),length(opposite)),
      paste('方向相反候选：',if(length(opposite))paste(opposite,collapse=', ')else'无达到该证据标准的候选'),
      'G/R 的效应比较、校准能力、随访窗口及位点证据见 Fig1.xlsx；R 不等于纯生活方式来源。')
  }
  if(pub_ok(t$performance,c('model','biom_set','AUC','C_index'))) {
    d<-t$performance|>filter(model=='Combined',is.finite(AUC))
    lines<-c(lines,'','验证集 Combined 模型（逐模型报告，不按验证集结果挑选最终模型）：','')
    for(i in seq_len(nrow(d)))lines<-c(lines,sprintf('- %s：C-index %.3f；%s 年 IPCW AUC %.3f；%s 个 biomarker。',d$biom_set[i],d$C_index[i],d$AUC_horizon_years[i],d$AUC[i],d$n_selected[i]))
  }
  if(pub_ok(t$same_region,'eligible'))lines<-c(lines,'',sprintf('同一 cis/local 区域 MR + 保守共定位支持：%d 个 biomarker（区域层面的候选，不等于因果证明）。',n_distinct(t$same_region$feature[t$same_region$eligible%in%TRUE])))
  if(pub_ok(t$mediation,'FDR_indirect'))lines<-c(lines,sprintf('探索性 mediation：%d 条路径，FDR <0.05 为 %d 条。',nrow(t$mediation),sum(t$mediation$FDR_indirect<.05,na.rm=TRUE)))
  lines<-c(lines,'','## 尚未建立的结论与写作边界','',
    '- PGS 反映遗传倾向，不能当作出生时实测 omics；需交代来源 GWAS、样本重叠、权重和预测能力。',
    '- 诊断前时间模式保留用于文章比较，方法中注明为不同人的一次基线采样。',
    '- 预测有效性与病因作用分开评价；无 MR 支持不等于没有临床预测价值。',
    '- 当前 clinical model 和论文中的 SCORE2/协变量不一定一致，不能仅凭 AUC 接近宣称严格复现。',
    '- C5 legacy evidence 的最小 P cis/trans 混选、不同区域合并和启发式权重不作为最终因果分级；Fig5/同区域表采用独立的保守规则。',
    '- 同期 LE8—omics mediation 不足以证明 lifestyle intervention 的因果效果。',
    '- 外部验证、锁定小 panel 后的完整选择流程验证、竞争风险及绝对风险可迁移性仍需另行建立。')
  if(pub_ok(t$matched_pgs,'replicated'))lines<-c(lines,sprintf('- 匹配 biomarker PGS 的遗传—实测 omic bridge：%d 项检验，%d 项两半样本复制；这不等同于因果中介链。',nrow(t$matched_pgs),sum(t$matched_pgs$replicated%in%TRUE)))
  else if(is.null(t$bridges)||!nrow(t$bridges))lines<-c(lines,'- 遗传—omics—疾病 bridge 表没有结果：四维 connection 尚未由此模块贯通。')
  if(pub_ok(t$focus_metrics,c('stratum','landmark','model','AUC'))) {
    z<-t$focus_metrics|>filter(stratum=='All',landmark==0)
    lines<-c(lines,'','## 同等 assay 预算的监督比较','',paste0('- ',z$model,': AUC=',sprintf('%.4f',z$AUC)),
      '','这些比较属于已查看旧结果之后的探索性分析；不能将最高 AUC 当作预先锁定的成功模型。Yin/Yang 的差异和 NS 对照见 focus_contrasts。')
  }
  if(pub_ok(t$coloc,'susie_status')&&all(t$coloc$susie_status!='ok',na.rm=TRUE))lines<-c(lines,'- 当前没有成功的 SuSiE 信号级 fine-mapping；ABF 的 H4 条件 SNP 后验不能当作 trait fine-mapping credible set。')
  if(pub_ok(t$dandelion_audit,c('metric','value'))&&any(t$dandelion_audit$metric=='primary_eligible'&t$dandelion_audit$value=='FALSE'))lines<-c(lines,'- DANDELION 当前标记 primary_eligible=FALSE，仅作为已记录的敏感性分析，不宣称完整复现原始方法。')
  missing<-source_audit|>filter(status!='available')
  if(nrow(missing))lines<-c(lines,'','缺失或空表：',paste0('- ',missing$table,'：',missing$source))
  lines<-c(lines,'','## 图表与文章结构','',paste0('- ',status$figure,'：',status$status,'（',status$panels,' panels）'),'',
    '正文可围绕“LE8-supervised、遗传证据分层且限制 assay 数量的疾病风险评估”组织；遗传、可干预解释和预测收益必须分别有对应结果，不能由统计关联推导 intervention promise。','',
    'Fig1 实测与遗传预测的分歧 → Fig2 遗传位点证据 → Fig3 LE8 connections → Fig4 预测与 assay budget → Fig5 证据整合。',
    'Fig*.xlsx 保存来源数据和 panel 映射；TablesS.xlsx 保存设计、缺口、同区域证据和入选附图路径。','')
  writeLines(lines,file.path(root,'final.results.md'),useBytes=TRUE)
  writeLines(c('# Figure legends','',unlist(lapply(names(captions),function(n)c(paste('##',n),'',captions[[n]],'')))),file.path(root,'final.legends.md'))
}
