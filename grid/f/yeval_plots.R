# Report layer. Scores and performance estimates are computed by Yeval.R.
colours <- c('GRID-tuned'='#278C78',GRID_shared='#76A88C',GRID_posterior='#5D9B91',GRID_matched='#86A9AA',
 COJO='#7C8798','PRS-CSx-auto-meta'='#E1A63B','PRS-CSx'='#3275B4','DiscoDivas-untuned'='#CD8B95',
 'DiscoDivas-tuned'='#D45260','PRS-CSx-fixed-meta'='#7755A2','csx.AFR'='#C6813B','csx.EAS'='#60A67A','csx.EUR'='#6593CC','csx.SAS'='#9E83BD')
ancestry_colours<-c(EUR='#437CB3',AFR='#D59421',EAS='#269B78',SAS='#AE6BA5',OTH='#8B9098',UNASSIGNED='#BABEC5')
extra<-setdiff(groups,names(ancestry_colours))
if(length(extra))ancestry_colours<-c(ancestry_colours,setNames(scales::hue_pal()(length(extra)),extra))
metric_label<-switch(type,ct='Prediction R²',dt='AUC (covariates + PRS)',t2e='Harrell C-index (covariates + PRS)')
fmt_n<-function(x)format(x,big.mark=',',scientific=FALSE,trim=TRUE)
fmt<-function(x,digits=4)ifelse(is.finite(x),formatC(x,digits=digits,format='f'),'NA')
ci_text<-paste0(nfold,'-fold out-of-fold predictions; ',if(nboot>0)paste0(nboot,' paired subject bootstrap resamples') else 'intervals disabled')
counts<-d[eligible==TRUE,.(N=.N,events=if(type=='ct')NA_integer_ else sum(outcome)),by=target]
target_label<-setNames(vapply(pops,function(g){
 z<-performance[target==g][1]
 if(!nrow(z)||is.na(z$N))return(paste0(g,'\nUnavailable'))
 paste0(g,'\nN = ',fmt_n(z$N),if(type!='ct')paste0('\n',if(type=='t2e')'Events' else 'Cases',' = ',fmt_n(z$events)))
},character(1)),pops)
axis_methods<-function(z){
 labels<-c(COJO='COJO','PRS-CSx-auto-meta'='Auto\nmeta','PRS-CSx-fixed-meta'='Fixed\nmeta','PRS-CSx'='Four-score\nfit',
           'DiscoDivas-tuned'='Disco\ntuned','DiscoDivas-untuned'='Disco\nuntuned','GRID-tuned'='GRID\ntuned')
 ans<-unname(labels[z]);ans[is.na(ans)]<-z[is.na(ans)];ans
}
plot_theme<-theme_classic(base_size=12)+theme(legend.position='bottom',legend.title=element_text(face='bold'),
 strip.background=element_rect(fill='#F2F4F7',colour=NA),strip.text=element_text(face='bold',margin=margin(8,5,8,5)),
 plot.title=element_text(face='bold',size=17),plot.subtitle=element_text(colour='#526071',size=11,lineheight=1.15),
 plot.caption=element_text(hjust=0,colour='#626A76',size=10),panel.spacing=grid::unit(1.2,'lines'),plot.margin=margin(12,14,10,12))
primary_caption<-switch(type,
 ct='Prediction R² = squared correlation of covariate-residualized phenotype and PRS prediction, evaluated on held-out people.',
 dt='AUC evaluates held-out case-control discrimination; covariates are included. ΔAUC and Brier are saved separately.',
 t2e='C-index measures ranking of event times, including covariates. Dashed lines show covariates-only C; 0.5 is chance ranking.')
make_comparison<-function(methods,title,subtitle,caption) {
 z<-merge(CJ(target=pops,method=methods,unique=TRUE),performance[target%in%pops],by=c('target','method'),all.x=TRUE)
 z[,method:=factor(method,levels=methods)];z[,target:=factor(target,levels=pops)]
 p<-ggplot(z,aes(method,estimate,fill=method))+
  geom_hline(yintercept=if(type!='ct').5 else 0,colour='#AEB7C2',linewidth=.35)+
  geom_col(width=.66,na.rm=TRUE,alpha=.95)+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.16,linewidth=.5,na.rm=TRUE)+
  geom_text(data=z[is.na(estimate)],aes(x=method,y=if(type=='t2e').51 else 0,label='Unavailable'),
            inherit.aes=FALSE,angle=90,hjust=0,size=2.8,colour='#777777')+
  facet_grid(.~target,labeller=labeller(target=target_label),drop=FALSE)+
  scale_fill_manual(name='Predictor',values=colours,limits=methods,breaks=methods,drop=FALSE)+
  scale_x_discrete(labels=axis_methods,drop=FALSE)+scale_y_continuous(expand=expansion(mult=c(.04,.08)))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+plot_theme+
  theme(axis.text.x=element_text(size=10,lineheight=.95),legend.text=element_text(size=11))+
  labs(title=title,subtitle=subtitle,x=NULL,y=metric_label,caption=caption)
 if(type=='t2e') {
  b<-unique(z[is.finite(baseline_C),.(target,baseline_C)])
  p<-p+geom_hline(data=b,aes(yintercept=baseline_C),linetype=2,colour='#303842',linewidth=.6)+
   coord_cartesian(ylim=c(min(.5,z$lower95,z$estimate,z$baseline_C,na.rm=TRUE),
                          min(1,max(.82,z$upper95,z$estimate,z$baseline_C,na.rm=TRUE)+.025)))
 }
 p
}
p1<-make_comparison(main_methods,paste(Y,'| Prediction by target ancestry'),
 'Overall benchmark: COJO, auto-meta, fitted PRS-CSx and DiscoDivas within each ancestry.\nThe same participants are used for every available predictor in a panel.',primary_caption)
p2<-make_comparison(c('PRS-CSx-auto-meta','PRS-CSx-fixed-meta','PRS-CSx',if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned'),
 paste(Y,'| Combined-score comparison'),
 'How CSx scores are combined: posterior meta-analysis (auto/fixed phi), four-score regression (PRS-CSx),\nor genetic-distance interpolation (DiscoDivas). COJO is omitted; shared bars repeat the first figure.',
 if(type=='t2e')primary_caption else 'auto/fixed meta combine SNP effects. Four-score regression learns target-ancestry weights in training folds.')

# Four distinct questions: ancestry, empirical prediction, geometry, posterior.
set.seed(seed+1L)
landscape<-d[eligible==TRUE,.SD[sample.int(.N,min(.N,2000L))],by=target]
shown_groups<-groups[groups%in%unique(landscape$target)]
category<-merge(data.table(target=shown_groups),performance[method=='PRS-CSx'],by='target',all.x=TRUE)
category[,target:=factor(target,levels=shown_groups)]
small_theme<-plot_theme+theme(plot.title=element_text(size=13,face='bold'),plot.subtitle=element_text(size=10),
 legend.text=element_text(size=10),axis.title=element_text(size=11),plot.tag=element_text(face='bold',size=18))
ancestry_scale<-function()scale_colour_manual(name='Ancestry',values=ancestry_colours,limits=shown_groups,drop=FALSE)
pa<-ggplot(landscape,aes(proj_PC1,proj_PC2,colour=target))+geom_point(size=.65,alpha=.55)+
 ancestry_scale()+small_theme+coord_equal()+
 guides(colour=guide_legend(override.aes=list(size=3,alpha=1),nrow=1))+
 labs(title='Ancestry groups',subtitle=paste0('Original labels: ',gc),x='Projected PC1',y='Projected PC2')

category[,`:=`(value=estimate,lo=lower95,hi=upper95)]
if(type=='t2e')category[,`:=`(value=delta_C,lo=delta_C_lower95,hi=delta_C_upper95)]
if(type=='dt')category[,`:=`(value=delta_AUC,lo=delta_AUC_lower95,hi=delta_AUC_upper95)]
category_labels<-setNames(vapply(shown_groups,function(g){
 z<-counts[target==g]
 paste0(g,'\nN=',fmt_n(z$N),if(type!='ct')paste0('\nEvents=',fmt_n(z$events)))
},character(1)),shown_groups)
gain_label<-switch(type,ct='Prediction R²',dt='ΔAUC from PRS',t2e='ΔC from PRS')
pb<-ggplot(category,aes(target,value,colour=target))+
 geom_hline(yintercept=0,colour='#CDD3DA',linewidth=.4)+
 geom_errorbar(aes(ymin=lo,ymax=hi),width=.12,linewidth=.6,na.rm=TRUE)+geom_point(size=3.2,na.rm=TRUE)+
 geom_text(aes(label=ifelse(is.finite(value),fmt(value,3),'Unavailable')),vjust=-1.3,size=3,na.rm=TRUE)+
 ancestry_scale()+scale_x_discrete(labels=category_labels,drop=FALSE)+
 scale_y_continuous(expand=expansion(mult=c(.12,.24)))+small_theme+theme(legend.position='none')+
 labs(title='PRS-CSx prediction by group',subtitle=if(type=='ct')'Covariate-adjusted prediction in held-out participants' else 'Improvement beyond covariates in held-out participants',x=NULL,y=gain_label)

# c uses a small number of real people to explain geometry, not another cloud
# recoloured by distance. All four source centres are shown, with no EUR default.
examples<-landscape[, .SD[which.min(abs(distance.analysis-median(distance.analysis)))],by=target]
examples[,point_id:=paste0('i',seq_len(.N))]
gc_plot<-copy(geometry$centers)
gc_plot[,label:=paste0(POP,if(geometry$status=='discovery')' training' else ' reference')]
segments<-rbindlist(lapply(seq_len(nrow(examples)),function(i){
 data.table(x=gc_plot$PC1,y=gc_plot$PC2,xend=examples$proj_PC1[i],yend=examples$proj_PC2[i],weight=gc_plot$mixture_weight)
}))
pcp<-ggplot()+
 geom_segment(data=segments,aes(x=x,y=y,xend=xend,yend=yend,linewidth=weight),colour='#AAB2BD',alpha=.65,
              arrow=grid::arrow(length=grid::unit(.07,'inches'),type='closed'))+
 scale_linewidth_continuous(range=c(.25,1.1),guide='none')+
 geom_point(data=gc_plot,aes(PC1,PC2),shape=17,size=3.3,colour='#273646')+
 geom_text(data=gc_plot,aes(PC1,PC2,label=label),nudge_y=.04*diff(range(landscape$proj_PC2)),size=3,check_overlap=TRUE)+
 geom_point(data=examples,aes(proj_PC1,proj_PC2,colour=target),size=3.2)+
 geom_text(data=examples,aes(proj_PC1,proj_PC2,label=point_id,colour=target),nudge_y=-.04*diff(range(landscape$proj_PC2)),size=3.4,show.legend=FALSE)+
 ancestry_scale()+scale_x_continuous(expand=expansion(mult=.18))+scale_y_continuous(expand=expansion(mult=.18))+small_theme+coord_equal()+theme(legend.position='none')+
 labs(title=if(geometry$status=='discovery')'Distance to GWAS training groups' else 'Reference geometry (exploratory)',
      subtitle=paste0('Triangles: centres; dots: example people; distances use ',npc,' PCs'),x='Projected PC1',y='Projected PC2')

individual_metric<-arg('individual-metric',if(type=='ct'&&arg('posterior-mode','required')!='off')'reliability' else 'sd')
if(!individual_metric%in%c('auto','reliability','sd'))stop('--individual-metric must be auto, reliability or sd')
if(individual_metric=='auto')individual_metric<-if(nrow(individual)&&any(is.finite(individual$individual_R2)))'reliability' else 'sd'
if(individual_metric=='reliability'&&(!nrow(individual)||!any(is.finite(individual$individual_R2))))stop('Individual reliability requires a valid --genetic-variance-file on the correct outcome scale')
if(nrow(individual)){
 individual[,plot_value:=if(individual_metric=='reliability')individual_R2 else posterior_sd]
 plotted<-individual[is.finite(plot_value)]
 plot_limit<-intarg('individual-max-points',5000,100)
 set.seed(seed+2L)
 dots<-plotted[,.SD[sample.int(.N,min(.N,plot_limit))],by=target]
 ylab<-if(individual_metric=='reliability')'Individual model-based R²' else switch(type,
     ct='Posterior SD of PRS prediction',dt='Posterior SD of PRS log-odds',t2e='Posterior SD of PRS log-hazard')
 pd<-ggplot(dots,aes(distance.analysis,plot_value,colour=target))+geom_point(size=.55,alpha=.3)+
  ancestry_scale()+scale_x_continuous(trans='sqrt',labels=scales::label_number())+
  small_theme+theme(legend.position='none')+
  labs(title=if(individual_metric=='reliability')'Individual PRS-CSx reliability' else 'Individual PRS-CSx uncertainty',
       subtitle='Each dot is one person; no bins or imposed decay curve',x=geometry$label,y=ylab)
 if(individual_metric=='reliability'&&min(dots$plot_value)<0)pd<-pd+geom_hline(yintercept=0,colour='#9DA8B4',linewidth=.4)
 missing_prior<-nrow(individual)-nrow(plotted)
 posterior_caption<-paste0('d: ',fmt_n(nrow(plotted)),' individual estimates; up to ',fmt_n(plot_limit),' dots per group displayed. ',
  if(individual_metric=='reliability')paste0('Model-based reliability is distinct from empirical prediction R² in b; ',missing_prior,' people lack a variance scale.') else 'Posterior SD is shown; lower values mean less uncertainty. This is not a prediction R².')
}else{
 pd<-ggplot()+theme_void()+labs(title='Individual posterior analysis disabled',subtitle='Run 1csx.sh --posterior TRUE, then Yeval with --posterior-mode required')+small_theme
 posterior_caption<-'Individual posterior estimates disabled explicitly.'
}
distance_caption<-paste0('a: original ancestry labels. c: lines illustrate PC1/PC2 only; numerical distances use all ',npc,' PCs.\n',
 if(geometry$status=='discovery')'Multi-training distance: sqrt(sum(N_group / N_total × distance_to_group²)).' else 'Reference proxy: equal weights over AFR/EAS/EUR/SAS 1KG centres; these are not the GWAS training centres.',
 '\n',posterior_caption)
p_distance<-(pa+pb)/(pcp+pd)+plot_annotation(title=paste(Y,'| Prediction along genetic distance'),
 subtitle='PRS-CSx: observed prediction by ancestry and posterior information for each person.',tag_levels='a',caption=distance_caption,
 theme=theme(plot.title=element_text(face='bold',size=18),plot.subtitle=element_text(size=12,colour='#526071'),plot.caption=element_text(hjust=0,size=9.5),plot.margin=margin(12,14,12,12)))

# The empirical distance bins remain a separate validation plot, not fake
# individual observations. Survival panels emphasize incremental discrimination.
binplot<-copy(distperf)
if(nrow(binplot)){
 binplot[,`:=`(value=estimate,lo=lower95,hi=upper95)]
 if(type=='t2e')binplot[,`:=`(value=delta_C,lo=delta_C_lower95,hi=delta_C_upper95)]
 if(type=='dt')binplot[,`:=`(value=delta_AUC,lo=delta_AUC_lower95,hi=delta_AUC_upper95)]
 p_bins<-ggplot(binplot,aes(distance,value,colour=target))+geom_hline(yintercept=0,colour='#ADB5C0')+
  geom_errorbar(aes(ymin=lo,ymax=hi),width=0,alpha=.6)+geom_point(size=2)+
  ancestry_scale()+facet_wrap(~target,scales='free_x',nrow=1)+scale_x_continuous(trans='sqrt')+plot_theme+
  labs(title=paste(Y,'| Empirical prediction across distance bins'),subtitle='Each point summarizes held-out participants in a distance bin.',x=geometry$label,y=gain_label,
       caption='These group estimates validate prediction on observed outcomes; they are separate from individual posterior reliability or uncertainty.')
}else p_bins<-ggplot()+theme_void()+labs(title='Insufficient samples/events for empirical distance bins')
if(nrow(comparison)){
 comparison[,target:=factor(target,levels=groups)]
 p_paired<-ggplot(comparison,aes(target,difference))+geom_hline(yintercept=0,linetype=2,colour='#777777')+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.12,colour=colours[disco_method],linewidth=.7)+
  geom_point(size=3.5,colour=colours[disco_method])+plot_theme+
  labs(title=paste(Y,'|',disco_method,'versus PRS-CSx'),subtitle='Positive differences favour DiscoDivas; negative differences favour PRS-CSx.',
       x='Target ancestry',y=paste0('Difference in ',switch(type,ct='prediction R²',dt='AUC',t2e='C-index')))
}else p_paired<-ggplot()+theme_void()+labs(title='DiscoDivas versus PRS-CSx',subtitle='Comparison unavailable or intervals disabled')
figures<-list(comparison=p1,combined_scores=p2,distance_performance=p_distance,paired_improvement=p_paired,distance_bins=p_bins)
if(capabilities('cairo'))grDevices::cairo_pdf(file.path(out,'plots.pdf'),width=14,height=11,onefile=TRUE) else pdf(file.path(out,'plots.pdf'),width=14,height=11,useDingbats=FALSE)
for(p in figures)print(p)
invisible(dev.off())
for(nm in names(figures))ggsave(file.path(out,paste0(nm,'.png')),figures[[nm]],width=14,height=if(nm=='distance_performance')11 else 6.8,dpi=170,bg='white')

methods_text<-c(paste0('# ',Y,' PRS evaluation'),'',paste('Outcome:',outcome_definition),paste('Covariates:',paste(covars,collapse=', ')),
 paste('Evaluation:',ci_text),'',
 '## Prediction metrics',
 '- Continuous: Prediction R² = cor(y - baseline_OOF, full_OOF - baseline_OOF)^2. This is covariate-adjusted squared prediction correlation. All residualization, score standardization and combination fitting use training folds.',
 '- OLS full_OOF - baseline_OOF equals the score-weighted, training-covariate-residualized score. No regression is fitted within a held-out distance bin.',
 '- Squared correlation does not assess calibration or direction. prediction_r, RMSE, full_R2, baseline_R2, delta_R2 and the former SSE_partial_R2 remain in performance.tsv.',
 '- Binary: main metric is logistic AUC; baseline_AUC, delta_AUC, paired delta intervals and Brier are also saved. No unvalidated liability conversion is applied.',
 '- Survival: main benchmark is covariates + PRS Harrell C. Comparisons use within-fold comparable pairs. The ancestry and empirical distance panels show delta_C to isolate PRS increment.',
 '- These metrics follow common PRS reporting conventions but are not a numerical reproduction of any publication, because training data, covariates, splits and outcomes differ.','',
 '## Four panels',
 '- a: original ancestry labels in PC space. b: empirical group prediction (R² or delta discrimination). c: separate geometric illustration with four source centres and representative real people. d: one individual per posterior point.',
 '- The first bar chart is the overall method benchmark including COJO. The second isolates combined-score strategies, omits COJO and adds fixed-meta; repeated bars are identical.',
 '- The DiscoDivas difference chart follows the four-panel figure. Empirical bins are a separate validation figure after that.','',
 '## Posterior model and interpretation',
 '- All four population beta draws share a retained iteration within a chromosome. Individual score covariance includes SNP LD and all cross-population covariance terms.',
 '- Scores are centred at harmonized discovery EAF. Missing genotypes are mean-imputed at the same EAF. The posterior means replace the four old uncentred scores for fitting, so means, scales and covariance refer to one predictor.',
 '- Chromosome means and covariance matrices are summed under PRS-CSx chromosome independence. Independent chromosomes are not artificially coupled by matching iteration numbers.',
 '- With training-fold weights w = regression_coefficient / training_score_SD, individual prediction variance is w^T Sigma_i w. Both diagonal-only and cross-population contributions are saved.',
 '- If an external genetic variance Vg on the correct scale is supplied, model-based individual R² = 1 - w^T Sigma_i w / Vg. For a continuous phenotype, supplied residual SNP h² is multiplied by training-fold covariate-residual phenotype variance.',
 '- This stacked-CSx reliability is an exploratory extension conditional on fitted combination/covariate weights and the supplied Vg, not an established equivalence to LDpred2 reliability. It requires model calibration and compatible priors/scales. It is not empirical phenotype prediction R².',
 '- Negative model-based reliability is retained as a diagnostic of variance scaling, uncertainty or model mismatch. No clipping, artificial dots or forced decay is applied.',
 '- If genetic variance is absent, d reports posterior SD explicitly; it is uncertainty, not accuracy. For survival this is SD of the PRS log-hazard contribution; for binary outcomes SD of log-odds. Liability h² cannot be used in their place.',
 '- Bootstrap intervals condition on fixed OOF fits. Posterior covariance conditions on GWAS summary statistics, LD reference and fitted combination weights; neither accounts for all sources of model misspecification.','',
 '## Genetic distance',geometry$description,
 '- Discovery centres must be means in the same PCA coordinate system, with source and pca_space identifiers. Multiple source groups use N_GWAS-weighted RMS distance, retaining distances to every group in the individual table.',
 '- A mixture centroid alone can conceal large distance from every actual training group. RMS distance avoids that cancellation, but is a descriptive extension, not a proven predictor of multi-ancestry accuracy.',
 '- Without discovery centres, a clearly labelled equal-weight four-reference proxy is used. It is not the former EUR-only distance and is never called a GWAS training distance.',
 paste0('- At most ',nbins,' empirical quantile bins per ancestry; minimum ',minn,' people, and for dt/t2e ',min_bin_events,' events/cases plus non-events/controls. Sparse bins are not used to infer precise trends.'),'',
 '## References',
 '- PRS-CS: https://doi.org/10.1038/s41467-019-09718-5',
 '- PRS-CSx: https://doi.org/10.1038/s41588-022-01054-7',
 '- Individual accuracy: https://doi.org/10.1038/s41586-023-06079-4',
 '- Source method: https://github.com/yidingdd/individual-pgs-accuracy/tree/main/pgs-accuracy')
writeLines(methods_text,file.path(out,'methods.md'))
escape<-function(x){x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
html_table<-function(z)paste0('<div class="table-wrap"><table><tr>',paste0('<th>',escape(names(z)),'</th>',collapse=''),'</tr>',
 paste(apply(as.data.frame(z),1,function(r)paste0('<tr>',paste0('<td>',escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table></div>')
csx_summary<-performance[method=='PRS-CSx']
if(type=='ct'){
 overview<-csx_summary[,.(Ancestry=target,N=fmt_n(N),`Prediction R²`=fmt(estimate),`Prediction r`=fmt(prediction_r),`SSE-based partial R²`=fmt(SSE_partial_R2),RMSE=fmt(RMSE,3))]
 metric_explanation<-'<p><strong>Prediction R²</strong> is the squared correlation between covariate-residualized observed and predicted phenotypes in held-out participants. It differs from the previous SSE-based partial R², which remains available for comparison.</p>'
}else if(type=='t2e'){
 overview<-csx_summary[,.(Ancestry=target,N=fmt_n(N),Events=fmt_n(events),`Full C`=fmt(estimate),`Covariates C`=fmt(baseline_C),`ΔC from PRS`=fmt(delta_C),`ΔC 95% CI`=paste0(fmt(delta_C_lower95),'–',fmt(delta_C_upper95)))]
 metric_explanation<-'<p><strong>C-index</strong> measures discrimination of event times. The main comparison bars show covariates + PRS; panel b and empirical distance bins show <strong>ΔC</strong> from adding PRS. Posterior SD in d measures uncertainty of the individual PRS log-hazard contribution; there is no individual C-index.</p>'
}else{
 overview<-csx_summary[,.(Ancestry=target,N=fmt_n(N),Cases=fmt_n(events),AUC=fmt(estimate),`Covariates AUC`=fmt(baseline_AUC),`ΔAUC`=fmt(delta_AUC),Brier=fmt(Brier))]
 metric_explanation<-'<p><strong>AUC</strong> measures held-out case-control discrimination. Panel b and empirical distance bins show ΔAUC from adding PRS. Posterior SD measures uncertainty of individual PRS log-odds.</p>'
}
shown<-copy(performance);for(nm in names(shown))if(is.numeric(shown[[nm]]))set(shown,j=nm,value=signif(shown[[nm]],5))
writeLines(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>',escape(Y),' PRS prediction</title>',
 '<style>body{font:16px system-ui;max-width:1440px;margin:32px auto;padding:0 24px;color:#263446;line-height:1.55}img{width:100%;height:auto;margin:18px 0}.table-wrap{overflow-x:auto}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:right;white-space:nowrap}th{background:#f2f4f7}td:first-child,th:first-child{text-align:left}.metric{background:#f2f6fa;padding:16px 22px;border-radius:8px}pre{white-space:pre-wrap;font:14px system-ui}a{color:#245f9b}details{margin:18px 0}</style>',
 '<h1>',escape(Y),' | PRS prediction comparison</h1><p><strong>Outcome:</strong> ',escape(outcome_definition),' · <strong>Covariates:</strong> ',escape(paste(covars,collapse=', ')),' · ',nfold,' held-out folds.</p>',
 '<div class="metric">',metric_explanation,'<p>',escape(geometry$description),'</p></div>',
 '<p><a href="plots.pdf">Figures PDF</a> · <a href="performance.tsv">Performance</a> · <a href="distance_performance.tsv">Empirical bins</a> · <a href="methods.md">Methods</a>',
 if(nrow(individual))' · <a href="individual_posterior.tsv.gz">Individual posterior estimates</a>' else '', '</p>',
 '<h2>PRS-CSx overview</h2>',html_table(overview),
 paste0('<img src="',names(figures),'.png" alt="',names(figures),'">',collapse=''),
 '<h2>All predictors</h2>',html_table(shown),
 if(length(skips))paste0('<details><summary>Skipped evaluations</summary>',html_table(rbindlist(skips)),'</details>') else '',
 '<details><summary>Evaluation inputs</summary>',html_table(manifest),'</details>',
 '<details><summary>Score definitions</summary>',html_table(model_map),'</details>',
 '<details><summary>Methods and references</summary><pre>',escape(paste(methods_text,collapse='\n')),'</pre></details></html>'),file.path(out,'report.html'))
