c5_joint_outputs<-function(results,yin,map,root,K,B,seed,primary_budget,primary_landmark,signature,cached_validation=NULL){
  private<-file.path(root,"_c5_private")
  c5_cleanup_stale_summaries(private)
  directory<-c5_summary_directory(private)
  on.exit(unlink(directory,recursive=TRUE),add=TRUE)
  store<-c5_stage_joint_results(results,signature,directory)
  if(is.null(yin)){
    if(is.null(cached_validation))stop("Cached C5 summaries require saved cohort validation")
    yin<-do.call(c5_validate_cached_cohorts,c(list(cohorts=store$cohorts),cached_validation))
  }
  combine<-function(nm)c5_store_table(store,nm)
  members<-combine("members")|>filter(feature%in%map$feature)|>distinct()
  diag<-combine("diagnostics")
  keys<-c("model","budget","arm","landmark","horizon")
  metrics<-cal<-dec<-foldmetrics<-pairs<-strata<-burden<-heterogeneity<-list()
  marker_groups<-combine("marker_groups")
  source(file.path(fdir,"c5_factorial_outputs.R"))
  design<-read.csv(file.path(root,"c5.factorial_design.csv"),stringsAsFactors=FALSE)
  factorial<-genetic_strata<-list()
  blocks<-unique(store$blocks[c("landmark","arm")])
  prediction_export<-file.path(directory,"out_of_fold_predictions.csv.gz")
  for(block in seq_len(nrow(blocks))){
  L0<-blocks$landmark[block];a0<-blocks$arm[block]
  pred<-c5_store_predictions(store,L0,a0)
  message("C5 summary: block ",block,"/",nrow(blocks),"; landmark=",L0,"; arm=",a0,"; prediction rows=",nrow(pred))
  groups<-split(seq_len(nrow(pred)),do.call(interaction,c(pred[keys],list(drop=TRUE,lex.order=TRUE))))
  for(ii in groups){
    d<-pred[ii,,drop=FALSE]
    id<-d[1,keys,drop=FALSE];expected<-length(c5_output_ids(yin,id$landmark))
    if(nrow(d)!=expected||anyDuplicated(d$eid)||length(unique(d$fold))!=K){
      metrics[[length(metrics)+1L]]<-cbind(id,status="incomplete outer predictions; do not compare",N=nrow(d));next
    }
    z<-le8_evaluate_risk(d$time,d$event,d$risk,id$horizon,id$model,id$budget,"nested outer validation",B=0)
    z$metrics$arm<-id$arm;z$metrics$landmark<-id$landmark
    metrics[[length(metrics)+1L]]<-z$metrics
    cal[[length(cal)+1L]]<-z$calibration|>mutate(arm=id$arm,landmark=id$landmark)
    dec[[length(dec)+1L]]<-z$decision|>mutate(arm=id$arm,landmark=id$landmark,budget=id$budget)
    for(fd in seq_len(K)){
      q<-d[d$fold==fd,,drop=FALSE]
      q$.c_time<-pmin(q$time,id$horizon);q$.c_event<-as.integer(q$event==1&q$time<=id$horizon)
      cc<-tryCatch(survival::concordance(survival::Surv(.c_time,.c_event)~lp,data=q,reverse=TRUE),error=function(e)NULL)
      foldmetrics[[length(foldmetrics)+1L]]<-cbind(id,fold=fd,N=nrow(q),events=sum(q$event),
        Harrell_C=if(is.null(cc))NA_real_ else cc$concordance,variance=if(is.null(cc))NA_real_ else cc$var,
        scope="within-fold Harrell C restricted to the stated prediction horizon")
    }
    if(id$budget==primary_budget&&id$model%in%c("Clinical_ProtMet_NS","Clinical_ProtMet_YS_YinYang","Clinical_ProtMet_PRS_NS")){
      for(bg in unique(d$inflammatory_burden)){
        q<-d[d$inflammatory_burden==bg,,drop=FALSE]
        zz<-le8_evaluate_risk(q$time,q$event,q$risk,id$horizon,id$model,id$budget,"exploratory marker-burden subgroup",B=0)
        burden[[length(burden)+1L]]<-zz$metrics|>mutate(burden=bg,arm=id$arm,landmark=id$landmark)
      }
    }
  }
  # Prespecified contrasts. No selection of the best validation budget/model.
  contrasts<-data.frame(model=c("Clinical_ProtMet_PRS_NS","Clinical_ProtMet_PRS_NS","Clinical_ProtMet_YS_YinYang",
      "Clinical_ProtMet_YS_YinYang","Clinical_ProtMet_YSplus_YinYang","Clinical_ProtMet_PRS_YS_YinYang",
      "Clinical_AddProxy_joint_YinYang","Clinical_ReplaceProxy_joint_YinYang"),
    reference=c("Clinical_ProtMet_NS","Clinical_PRS","Clinical_ProtMet_NS","Clinical_ProtMet_YS_Yin",
      "Clinical_ProtMet_NS","Clinical_ProtMet_PRS_NS","Clinical","Clinical"))
  contrasts<-bind_rows(contrasts,data.frame(
    model=c("Clinical_MatchedMeasured_PGS_NS","Clinical_MatchedMeasured_PGS_NS","Clinical_MatchedMeasured_PGS_PRS_NS",
      "Clinical_ProtMet_NS","Clinical_ProtMet_NS","Clinical_ProtMet_YSplus_YinYang"),
    reference=c("Clinical_MatchedMeasured_NS","Clinical_MatchedPGS_NS","Clinical_MatchedMeasured_PGS_NS",
      "Clinical_Protein_sharedBudget_NS","Clinical_Metabolite_sharedBudget_NS","Clinical_ProtMet_YSplus_Yin")))
  contrasts<-bind_rows(contrasts,data.frame(
    model=c("Clinical_ProtMet_PGSselected_NS","Clinical_ProtMet_PGSselected_NS","Clinical_ProtMet_PGSselected_PRS_NS","Clinical_ProtMet_YS_PGSselected_YinYang"),
    reference=c("Clinical_ProtMet_NS","Clinical_PGSselected_NS","Clinical_ProtMet_PGSselected_NS","Clinical_ProtMet_PGSselected_NS")))
  for(L in sort(unique(pred$landmark)))for(a in unique(pred$arm))for(k in sort(unique(pred$budget[pred$budget>0])))for(i in seq_len(nrow(contrasts))){
    co<-contrasts[i,];x<-pred[pred$model==co$model&pred$landmark==L&pred$arm==a&pred$budget==k,,drop=FALSE]
    kb<-if(co$reference%in%c("Clinical","Clinical_PRS"))0 else k
    y<-pred[pred$model==co$reference&pred$landmark==L&pred$arm==a&pred$budget==kb,,drop=FALSE]
    expected<-length(c5_output_ids(yin,L))
    if(nrow(x)!=expected||nrow(y)!=expected||anyDuplicated(x$eid)||anyDuplicated(y$eid))next
    h<-x$horizon[1]
    delta<-c5_paired_delta(x,y,h,B,seed+as.integer(L*100)+k)
    if(k==primary_budget&&co$reference=="Clinical_ProtMet_NS"&&co$model%in%c("Clinical_ProtMet_YS_YinYang","Clinical_ProtMet_YSplus_YinYang")){
      mg<-marker_groups[marker_groups$landmark==L,,drop=FALSE];mg<-mg[,setdiff(names(mg),c("fold","landmark")),drop=FALSE]
      he<-c5_marker_heterogeneity(x,y,mg,h,B,seed+as.integer(L*100)+k)
      if(nrow(he))heterogeneity[[length(heterogeneity)+1L]]<-he|>mutate(model=co$model,reference=co$reference,landmark=L,arm=a,budget=k)
    }
    pairs[[length(pairs)+1L]]<-delta|>mutate(model=co$model,reference=co$reference,budget=k,landmark=L,arm=a,
      primary=(k==primary_budget&L==primary_landmark&a=="all_assays"&
        co$model=="Clinical_ProtMet_YS_YinYang"&co$reference=="Clinical_ProtMet_NS"))
  }
  # Frozen training 75th-percentile thresholds in each fold. Show observed risk
  # in the four PRS x omic strata; no cutpoints optimized on held-out outcomes.
  for(L in sort(unique(pred$landmark)))for(a in unique(pred$arm)){
    p<-pred[pred$model=="ProtMet_only_NS"&pred$budget==primary_budget&pred$landmark==L&pred$arm==a,,drop=FALSE]
    g<-pred[pred$model=="PRS_only"&pred$landmark==L&pred$arm==a,,drop=FALSE]
    if(!nrow(p)||!nrow(g))next
    d<-merge(p,g[,c("eid","high_training_q75")],by="eid",suffixes=c(".omics",".PRS"))
    d$stratum<-paste0("omics ",ifelse(d$high_training_q75.omics,"high","lower")," / PRS ",ifelse(d$high_training_q75.PRS,"high","lower"))
    iw<-le8_ipcw(d$time,d$event,d$horizon[1])
    if(iw$status!="ok")next
    d$w<-iw$w;d$y<-iw$y
    strata[[length(strata)+1L]]<-d|>group_by(stratum)|>summarise(N=n(),events=sum(y),
      observed_risk_ipcw=sum(w*y)/sum(w),.groups="drop")|>mutate(landmark=L,arm=a,
        interpretation="Cause-specific net-risk display; high defined from training 75th percentile, not validated clinical cutoffs")
  }
  genetic_strata[[block]]<-c5_genetic_measured_strata(pred,primary_budget)
  factorial[[block]]<-c5_factorial_outputs(pred,yin,design,K,B,seed,primary_budget,primary_landmark)
  data.table::fwrite(pred,prediction_export,compress="gzip",append=block>1L)
  c5_release_predictions(store,L0,a0)
  rm(pred,groups);invisible(gc())
  }
  overlap<-members|>select(model,budget,arm,fold,landmark,feature)
  overlap_rows<-list()
  for(L in unique(overlap$landmark))for(a in unique(overlap$arm))for(k in unique(overlap$budget))for(fd in seq_len(K)){
    z<-overlap[overlap$landmark==L&overlap$arm==a&overlap$budget==k&overlap$fold==fd,,drop=FALSE]
    for(nm in c("Clinical_ProtMet_YS_YinYang","Clinical_ProtMet_YSplus_YinYang")){
      u<-z$feature[z$model==nm];v<-z$feature[z$model=="Clinical_ProtMet_YS_Yin"]
      if(length(u)&&length(v))overlap_rows[[length(overlap_rows)+1L]]<-tibble(model=nm,reference="Clinical_ProtMet_YS_Yin",
        fold=fd,landmark=L,arm=a,budget=k,identical_panel=setequal(u,v),Jaccard=length(intersect(u,v))/length(union(u,v)))
    }
  }
  het<-bind_rows(heterogeneity);if(nrow(het))het$FDR<-p.adjust(het$p,"BH")
  tables<-list(biomarker_PGS_omics_strata=bind_rows(genetic_strata),marker_heterogeneity=het,marker_definitions=combine("marker_audit"),
    domain_support=combine("domain_status"),pgs_reconstruction=combine("pgs_reconstruction"),pgs_decomposition_coefficients=combine("pgs_coefficients"),
    cross_omic_links=combine("cross_omic"),predictor_inventory=combine("members"),metrics=bind_rows(metrics),paired_contrasts=bind_rows(pairs),calibration=bind_rows(cal),decision_curves=bind_rows(dec),
    fold_C_index=bind_rows(foldmetrics),panel_members=left_join(members,map,by="feature"),panel_overlap=bind_rows(overlap_rows),
    panel_counts=members|>left_join(map,by="feature")|>group_by(model,budget,arm,fold,landmark)|>
      summarise(total_assays=n_distinct(feature),protein_assays=n_distinct(feature[layer=="prot"]),metabolite_assays=n_distinct(feature[layer=="met"]),.groups="drop"),
    panel_stability=members|>group_by(model,budget,arm,landmark,feature)|>summarise(selected_folds=n_distinct(fold),selection_frequency=n_distinct(fold)/K,.groups="drop"),
    fit_diagnostics=diag,coefficients=combine("coefficients"),model_preprocessing=combine("risk_preprocessing"),
    baseline_hazards=combine("baseline_hazards"),proxy_validation=combine("proxy_validation"),
    proxy_coefficients=combine("proxy_weights"),proxy_preprocessing=combine("proxy_preprocessing"),domain_associations=combine("domain_associations"),
    training_screens=combine("screen"),genetic_training_screens=combine("genetic_screen"),PRS_omics_strata=bind_rows(strata),inflammatory_burden=bind_rows(burden))
  factorial_tables<-setNames(lapply(names(factorial[[1]]),function(nm)bind_rows(lapply(factorial,`[[`,nm))),names(factorial[[1]]))
  tables<-c(tables,factorial_tables)
  provenance<-read.csv(file.path(root,"c5.prs_provenance.csv"),stringsAsFactors=FALSE)
  prs_scope<-paste(provenance$status,collapse="; ")
  for(nm in names(tables)[startsWith(names(tables),"factorial_")])
    if(nrow(tables[[nm]]))tables[[nm]]$disease_PRS_scope<-prs_scope
  for(nm in names(tables))write_raw_csv(tables[[nm]],paste0("c5.",nm,".csv"),root)
  primary_rows<-tables$paired_contrasts
  primary_ok<-nrow(primary_rows)>0&&any(primary_rows$primary%in%TRUE&primary_rows$status=="ok")
  write.csv(data.frame(budget=primary_budget,landmark=primary_landmark,
    comparison="Clinical_ProtMet_YS_YinYang vs Clinical_ProtMet_NS",
    status=if(primary_ok)"estimated; inspect effect and uncertainty"else"not estimable; no alternative primary selected"),
    file.path(root,"c5.primary_status.csv"),row.names=FALSE)
  if(!file.rename(prediction_export,file.path(root,"_c5_private","out_of_fold_predictions.csv.gz")))
    stop("Cannot publish C5 out-of-fold predictions")
  write.csv(data.frame(signature=signature,mode=if(is.null(cached_validation))"current run"else"frozen-checkpoint recovery",
    folds=K,bootstrap=B,seed=seed,primary_budget=primary_budget,primary_landmark=primary_landmark,
    prediction_blocks=nrow(blocks)),file.path(root,"c5.summary_provenance.csv"),row.names=FALSE)
  writeLines(c("C5 FINAL — internal, common-cohort, grouped outer validation.",
    paste("Primary: YS YinYang vs NS; total budget",primary_budget,"; landmark",primary_landmark),
    "Fixed end at baseline year 10 by default: landmark 5 evaluates years 5–10 among those still event-free and observed at year 5.",
    "Each landmark retrains and reselects using eligible training participants. This does not identify causal biomarkers.",
    "Clinical includes continuous BMI and non-HDL. Proxy replacement and incremental addition are separate tests.",
    "Proxy preprocessing/selection/fit are cross-fitted; missing OOF values never become in-sample predictions.",
    "Natriuretic/GDF15 omission excludes candidates before equal-budget re-selection; it is not an estimate of causal status.",
    "No automatic inflammatory/non-inflammatory disease subtype is inferred from these markers.",
    "IPCW assumes independent censoring; death is censored. Risks are net risks, not real-world competing-risk CIFs.",
    "Paired bootstrap conditions on trained models; it does not include full training/selection uncertainty.",
    "Secondary contrasts, strata, domain reconstructions and cell enrichments are exploratory.",
    "Clinical deployment claims require external validation, recalibration and competing-risk assessment."),file.path(root,"c5.interpretation.txt"))
  tables$metrics$validation_scope<-"Internal measured-omics validation; no external confirmation"
  pgsmodels<-grepl("PGS",tables$metrics$model)
  tables$metrics$validation_scope[pgsmodels]<-paste("Biomarker PGS:",paste(readLines(file.path(root,"c5.omic_PGS_status.txt")),collapse=" "))
  prsmodels<-grepl("PRS",tables$metrics$model)|tables$metrics$model%in%design$model[design$G]
  tables$metrics$validation_scope[prsmodels]<-paste(tables$metrics$validation_scope[prsmodels],"Disease PRS:",prs_scope)
  write_raw_csv(tables$metrics,"c5.metrics.csv",root)
  source(file.path(fdir,"c5_factorial_plots.R"));c5_factorial_plots(tables,design,root,primary_budget,primary_landmark)
  c5_joint_plots(tables,root,primary_budget)
  c5_cell_output(tables$panel_members,map,root)
  source(file.path(fdir,"c5_panels_final.R"));c5_final_panels(tables,root,primary_budget)
  c5_evidence_output(tables$panel_members,root)
  c5_write_joint_checkpoint(list(signature=signature,tables=tables,complete=all(diag$status=="ok"),summary_complete=TRUE,
    summary_settings=list(bootstrap=B,seed=seed,primary_budget=primary_budget,primary_landmark=primary_landmark)),file.path(root,"c5.res.rds"))
  c5_retire_fold_checkpoints(store,root,signature)
  invisible(tables)
}
c5_evidence_output<-function(panels,root){
  # Post-validation annotation only. Never reuse these full-cohort statistics
  # for screening, tuning or claiming independent replication of prediction.
  ledger<-evidence<-list()
  for(layer in c("prot","met")){
    candidates<-unique(panels$assay[panels$layer==layer])
    sources<-data.frame(module=c("c1_correlate","c2_cause","c2_cause","c3_coloc","c4_connect"),
      file=c(if(layer=="prot")"pwas_incident_adj2.csv"else"mwas_incident_adj2.csv","c2.MR_all.csv",
        "c2.dandelion_targets_all.csv","c3.coloc_summary.csv","c4.mediation_all.csv"),
      key=c("term","exposure","feature","feature","feature"))
    for(i in seq_len(nrow(sources))){
      s<-sources[i,];path<-file.path(root,layer,s$module,s$file)
      optionfile<-file.path(dirname(path),"analysis_options.rds")
      current<-file.exists(optionfile)&&identical(readRDS(optionfile),le8_analysis_options())
      status<-if(!file.exists(path))"unavailable"else if(current)"current baseline contract"else"legacy or unverified baseline contract; rerun before inference"
      ledger[[length(ledger)+1L]]<-tibble(layer,module=s$module,file=s$file,status,used_for_prediction=FALSE)
      if(!file.exists(path))next
      z<-as.data.frame(data.table::fread(path));if(!s$key%in%names(z))next
      z<-z[as.character(z[[s$key]])%in%candidates,,drop=FALSE]
      if(!nrow(z))next
      z$assay<-as.character(z[[s$key]]);z$record<-seq_len(nrow(z))
      names(z)[names(z)=="layer"]<-"source_layer"
      fields<-setdiff(names(z),c("assay","record"));z[fields]<-lapply(z[fields],as.character)
      evidence[[length(evidence)+1L]]<-tidyr::pivot_longer(z,cols=all_of(fields),names_to="statistic",values_to="value")|>
        mutate(layer,module=s$module,source_file=s$file,baseline_contract=status,
          interpretation="Descriptive triangulation; correlated evidence, no causal voting score")
    }
  }
  write_raw_csv(bind_rows(ledger),"c5.evidence_availability.csv",root)
  write_raw_csv(bind_rows(evidence),"c5.evidence_long.csv",root)
}
c5_joint_plots<-function(tables,root,k){
  savefig<-function(p,n,w=11,h=7)ggplot2::ggsave(file.path(root,n),p,width=w,height=h,dpi=200)
  m<-tables$metrics|>filter(status=="ok",budget%in%c(0,k),arm=="all_assays")
  keep<-c("Clinical","Clinical_PRS","Clinical_Protein_NS","Clinical_Metabolite_NS","Clinical_ProtMet_NS","Clinical_ProtMet_PRS_NS","Clinical_ProtMet_YS_YinYang")
  if(nrow(m))savefig(ggplot(m|>filter(model%in%keep),aes(AUC,reorder(model,AUC),color=factor(landmark)))+
    geom_point(size=2)+labs(x="Out-of-fold IPCW AUC",y=NULL,color="Landmark (years)",title="Joint modality comparisons on the same participants")+theme_bw(),"c5.Fig1.joint_modalities.png")
  p<-tables$paired_contrasts
  if(nrow(p))savefig(ggplot(p|>filter(budget==k,arm=="all_assays"),aes(delta_AUC,paste(model,"vs",reference),color=factor(landmark)))+
    geom_vline(xintercept=0,linetype=2)+geom_errorbar(aes(xmin=AUC_lo,xmax=AUC_hi),orientation="y",width=.15)+geom_point()+
    labs(x="Paired delta AUC (conditional 95% bootstrap CI)",y=NULL,color="Landmark")+theme_bw(),"c5.Fig2.paired_increment.png",13,8)
  q<-tables$proxy_validation
  if(nrow(q))savefig(ggplot(q|>filter(arm=="all_assays"),aes(modality,R2_vs_training_mean,color=cohort))+
    geom_point(position=position_jitter(width=.08,height=0))+facet_grid(target~landmark)+
    labs(y="Held-out reconstruction R²",x=NULL,title="BMI and lipid proxies: reconstruction is distinct from disease benefit")+theme_bw(),"c5.Fig3.domain_proxies.png")
  q<-tables$PRS_omics_strata
  if(nrow(q))savefig(ggplot(q|>filter(arm=="all_assays"),aes(stratum,observed_risk_ipcw,fill=stratum))+
    geom_col()+facet_wrap(~landmark)+labs(y="Observed IPCW net risk",x=NULL,title="Complementary PRS and omic information")+
    theme_bw()+theme(axis.text.x=element_text(angle=25,hjust=1),legend.position="none"),"c5.Fig4.PRS_omics_strata.png")
  q<-tables$calibration|>filter(model%in%keep,budget%in%c(0,k),arm=="all_assays")
  if(nrow(q))savefig(ggplot(q,aes(predicted,observed_ipcw,color=model))+geom_abline(slope=1,intercept=0,linetype=2)+
    geom_line()+geom_point()+facet_wrap(~landmark)+theme_bw()+labs(x="Predicted net risk",y="Observed IPCW risk"),"c5.Fig5.calibration.png")
  q<-tables$metrics|>filter(status=="ok",budget==k,model%in%c("Clinical_ProtMet_NS","Clinical_ProtMet_YS_YinYang"))
  if(nrow(q))savefig(ggplot(q,aes(landmark,AUC,color=arm,linetype=model))+geom_line()+geom_point()+theme_bw()+
    labs(x="Retraining landmark (years)",y="IPCW AUC",title="Longer lead-time and prespecified marker omission"),"c5.Fig6.leadtime_omission.png")
}
c5_cell_output<-function(panels,map,root){
  prot<-map[map$layer=="prot",,drop=FALSE]
  # Exact symbols only; NTPROBNP is a fragment of NPPB, not a separate gene.
  prot$gene<-toupper(prot$assay);prot$gene[prot$gene=="NTPROBNP"]<-"NPPB"
  mf<-Sys.getenv("C5_ASSAY_GENE_MAP","")
  if(nzchar(mf)){
    g<-read.csv(mf,stringsAsFactors=FALSE)
    if(!all(c("assay","gene")%in%names(g))||anyDuplicated(g$assay))stop("Assay-gene map requires unique assay rows")
    prot$gene<-g$gene[match(prot$assay,g$assay)]
  }
  universe<-data.frame(assay=prot$feature,gene=prot$gene)
  p<-panels|>filter(layer=="prot")|>mutate(model=paste(model,budget,arm,landmark,fold,sep="|"))|>select(model,feature)|>distinct()
  write.csv(universe,file.path(root,"c5.cell_assay_universe.csv"),row.names=FALSE)
  write.csv(p,file.path(root,"c5.cell_panels.csv"),row.names=FALSE)
  atlas<-Sys.getenv("C5_CELL_ATLAS",file.path(dirname(fdir),"data/cellage/atlas.csv"))
  code<-file.path(fdir,"c3_cell_annotation.py")
  status<-system2(Sys.getenv("PYTHON_BIN","python3"),vapply(c(code,"--universe",file.path(root,"c5.cell_assay_universe.csv"),
    "--atlas",atlas,"--panels",file.path(root,"c5.cell_panels.csv"),"--outdir",root,"--prefix","c5.cell"),shQuote,character(1)))
  if(status!=0)stop("Cell annotation failed; inspect Python dependencies and assay map")
}
