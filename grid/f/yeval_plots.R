# Report layer. Scores and performance estimates are computed by Yeval.R.
colours <- c('GRID-tuned'='#278C78',GRID_shared='#76A88C',GRID_posterior='#5D9B91',GRID_matched='#86A9AA',
 COJO='#7C8798','PRS-CSx-auto-meta'='#E1A63B','PRS-CSx'='#3275B4','DiscoDivas-untuned'='#CD8B95',
 'DiscoDivas-tuned'='#D45260','PRS-CSx-fixed-meta'='#7755A2','csx.AFR'='#C6813B','csx.EAS'='#60A67A','csx.EUR'='#6593CC','csx.SAS'='#9E83BD')
ancestry_colours<-c(EUR='#437CB3',AFR='#D59421',EAS='#269B78',SAS='#AE6BA5',OTH='#8B9098',UNASSIGNED='#BABEC5')
extra<-setdiff(groups,names(ancestry_colours))
if(length(extra))ancestry_colours<-c(ancestry_colours,setNames(scales::hue_pal()(length(extra)),extra))
metric_label<-switch(type,ct='Partial R²',dt='Observed-scale partial R²',t2e='Harrell C-index (covariates + PRS)')
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
 ct='Partial R² = 1 - SSE(covariates + PRS) / SSE(covariates), evaluated on held-out people.',
 dt='Partial R² uses linear probability predictions on the observed scale; logistic AUC/Brier are reported separately.',
 t2e='C-index measures ranking of event times, including covariates. Dashed lines show covariates-only C; 0.5 is chance ranking.')
make_comparison<-function(methods,title,subtitle,caption) {
 z<-merge(CJ(target=pops,method=methods,unique=TRUE),performance[target%in%pops],by=c('target','method'),all.x=TRUE)
 z[,method:=factor(method,levels=methods)];z[,target:=factor(target,levels=pops)]
 p<-ggplot(z,aes(method,estimate,fill=method))+
  geom_hline(yintercept=if(type=='t2e').5 else 0,colour='#AEB7C2',linewidth=.35)+
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

# Same people and coordinates in a/c; retain the original input ancestry labels.
# No artificial cluster boundaries, individual R² values or fitted decay curve.
set.seed(seed+1L)
landscape<-d[eligible==TRUE,.SD[sample.int(.N,min(.N,2000L))],by=target]
shown_groups<-groups[groups%in%unique(landscape$target)]
landscape[,target:=factor(target,levels=shown_groups)]
category<-merge(data.table(target=shown_groups),performance[method=='PRS-CSx'],by='target',all.x=TRUE)
category[,target:=factor(target,levels=shown_groups)]
binplot<-copy(distperf[method=='PRS-CSx'])
binplot[,target:=factor(target,levels=shown_groups)]
small_theme<-plot_theme+theme(plot.title=element_text(size=13,face='bold'),plot.subtitle=element_text(size=10),
                            legend.text=element_text(size=10),axis.title=element_text(size=11),plot.tag=element_text(face='bold',size=18))
pa<-ggplot(landscape,aes(proj_PC1,proj_PC2,colour=target))+geom_point(size=.65,alpha=.55)+
 scale_colour_manual(values=ancestry_colours,limits=shown_groups,drop=FALSE)+small_theme+
 guides(colour=guide_legend(override.aes=list(size=3,alpha=1),nrow=1))+
 labs(title='Original ancestry groups',subtitle=paste0('Input labels: ',gc),x='Reference-projected PC1',y='Reference-projected PC2',colour='Ancestry')
category_labels<-setNames(vapply(shown_groups,function(g){
 c0<-counts[target==g]
 paste0(g,'\nN=',fmt_n(c0$N),if(type!='ct')paste0('\n',if(type=='t2e')'E=' else 'Cases=',fmt_n(c0$events)))
},character(1)),shown_groups)
y_values<-c(category$estimate,category$lower95,category$upper95,binplot$estimate,binplot$lower95,binplot$upper95,
             if(type=='t2e')category$baseline_C else 0)
y_range<-range(y_values[is.finite(y_values)],if(type=='t2e').5 else 0)
y_padding<-max(.015,diff(y_range)*.12)
y_limits<-c(y_range[1]-y_padding,y_range[2]+y_padding)
if(type=='t2e')y_limits<-pmax(0,pmin(1,y_limits))
pb<-ggplot(category,aes(target,estimate,colour=target))+
 geom_hline(yintercept=if(type=='t2e').5 else 0,colour='#CDD3DA',linewidth=.35)+
 geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.14,linewidth=.65,na.rm=TRUE)+geom_point(size=3.3,na.rm=TRUE)+
 geom_text(aes(label=ifelse(is.finite(estimate),fmt(estimate,3),'')),vjust=-1.05,size=3.1,show.legend=FALSE,na.rm=TRUE)+
 geom_text(data=category[is.na(estimate)],aes(y=y_limits[1]+y_padding,label='Unavailable'),size=3,show.legend=FALSE)+
 scale_colour_manual(values=ancestry_colours,limits=shown_groups,drop=FALSE)+scale_x_discrete(labels=category_labels,drop=FALSE)+
 coord_cartesian(ylim=y_limits)+small_theme+theme(legend.position='none')+
 labs(title='PRS-CSx by ancestry group',subtitle=if(type=='t2e')'Full-model C; black ticks = covariates-only C' else 'One held-out estimate per original group',x=NULL,y=metric_label)
if(type=='t2e')pb<-pb+geom_point(aes(y=baseline_C),shape=95,size=7,colour='#303842',na.rm=TRUE)
eu_center<-data.table(proj_PC1=centers[POP=='EUR',PC1],proj_PC2=centers[POP=='EUR',PC2])
pcp<-ggplot(landscape,aes(proj_PC1,proj_PC2,colour=distance.EUR))+geom_point(size=.65,alpha=.6)+
 geom_point(data=eu_center,aes(proj_PC1,proj_PC2),inherit.aes=FALSE,shape=4,size=3.8,stroke=1.2,colour='black')+
 scale_colour_viridis_c(option='C',trans='sqrt')+small_theme+
 guides(colour=guide_colourbar(barwidth=grid::unit(6,'cm'),barheight=grid::unit(.35,'cm')))+
 labs(title='Continuous genetic distance',subtitle=paste0('Same people; distance uses ',npc,' PCs; cross = EUR centre'),
      x='Reference-projected PC1',y='Reference-projected PC2',colour='Distance to 1KG EUR')
if(nrow(binplot)) {
 pd<-ggplot(binplot,aes(distance,estimate,colour=target,group=target))+
  geom_hline(yintercept=if(type=='t2e').5 else 0,colour='#CDD3DA',linewidth=.35)+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=0,linewidth=.5,alpha=.6,na.rm=TRUE)+
  geom_line(linewidth=.65,na.rm=TRUE)+geom_point(size=2.4,na.rm=TRUE)+
  scale_colour_manual(values=ancestry_colours,limits=shown_groups,drop=FALSE)+
  scale_x_continuous(trans='sqrt',labels=scales::label_number())+
  coord_cartesian(ylim=y_limits)+small_theme+theme(legend.position='none')+
  labs(title='PRS-CSx along genetic distance',subtitle='Each point = a distance bin; colours = original ancestry groups',
       x='Distance to 1KG EUR centre (square-root axis)',y=metric_label)
} else pd<-ggplot()+theme_void()+small_theme+labs(title='PRS-CSx along genetic distance',subtitle='Insufficient samples/events for at least two bins per group')
distance_caption<-paste0('PCA displays at most 2,000 eligible participants per group; a/c use exactly the same people. Panels b/d use all evaluable people.\n',
 'Panel d: up to ',nbins,' quantile bins/group; at least ',minn,' people/bin',
 if(type!='ct')paste0(' and ',min_bin_events,' events/cases plus ',min_bin_events,' non-events/controls/bin') else '',
 '. Sparse groups have fewer or no bins.\n',
 'Distances are to a 1KG EUR reference centre, not the multi-ancestry discovery GWAS. Points estimate group/bin performance, not individual accuracy.')
p_distance<-(pa+pb)/(pcp+pd)+plot_annotation(title=paste(Y,'| Prediction along genetic distance'),
 subtitle='PRS-CSx only: original ancestry groups above, the continuous distance view below.',tag_levels='a',caption=distance_caption,
 theme=theme(plot.title=element_text(face='bold',size=18),plot.subtitle=element_text(size=12,colour='#526071'),plot.caption=element_text(hjust=0,size=10),plot.margin=margin(12,14,12,12)))

if(nrow(comparison)) {
 comparison[,target:=factor(target,levels=groups)]
 p_paired<-ggplot(comparison,aes(target,difference))+geom_hline(yintercept=0,linetype=2,colour='#777777')+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.12,colour=colours[disco_method],linewidth=.7)+
  geom_point(size=3.5,colour=colours[disco_method])+plot_theme+
  labs(title=paste(Y,'|',disco_method,'versus PRS-CSx'),
       subtitle='Additional comparison: positive values favour DiscoDivas; negative values favour PRS-CSx.',
       x='Target ancestry',y=paste0('Difference in ',if(type=='t2e')'Harrell C-index' else metric_label),
       caption='Differences use the same participants and held-out folds for both predictors.')
} else p_paired<-ggplot()+theme_void()+labs(title='DiscoDivas versus PRS-CSx',subtitle='Comparison unavailable (or bootstrap disabled).')
figures<-list(comparison=p1,combined_scores=p2,distance_performance=p_distance,paired_improvement=p_paired)
# Cairo preserves R²/Greek characters in the PDF where available.
if(capabilities('cairo'))grDevices::cairo_pdf(file.path(out,'plots.pdf'),width=14,height=10,onefile=TRUE) else
 pdf(file.path(out,'plots.pdf'),width=14,height=10,useDingbats=FALSE)
for(p in figures)print(p)
dev.off()
for(nm in names(figures))ggsave(file.path(out,paste0(nm,'.png')),figures[[nm]],width=14,
                              height=if(nm=='distance_performance')11 else 6.8,dpi=160,bg='white')

metric_methods<-switch(type,
 ct=c('- Partial R² = 1 - SSE(full)/SSE(covariates) = (full_R2 - baseline_R2)/(1 - baseline_R2). Both models predict held-out participants.',
      '- delta_R2 uses total phenotype variance as its denominator; partial R² uses residual prediction error after covariates. Neither is adjusted R² for model complexity.',
      '- baseline_SSE, full_SSE, baseline_R2, full_R2, delta_R2 and RMSE are saved so the calculation can be checked. Negative held-out estimates are retained.',
      '- The PRS-CS paper describes covariate-adjusted observed-versus-predicted R²; its text alone does not establish numerical identity with this OOF SSE metric. Squared correlation and predictive R² can differ when calibration differs.'),
 dt=c('- Main metric: observed-scale partial R² from held-out linear probability models. Logistic models supply AUC, baseline_AUC, delta_AUC and Brier.',
      '- Partial R² is not liability-scale R². Prevalence K is descriptive and is not used to convert partial R².',
      '- Default t2dm dt derives baseline ICD10 status from Yr2e/Yt2e. Use --phenotype-col for a different explicit definition.'),
 t2e=c('- Main metric: covariates + PRS Harrell C-index, weighted by comparable pairs within each held-out fold. Censored observations are retained.',
       '- baseline_C uses the covariates-only model on the same held-out people and folds. delta_C = full C - baseline C; its intervals use paired bootstrap differences.',
       '- C measures ranking, not explained variance or the probability of future disease. It must not be compared numerically with partial R², AUC or liability R².',
       '- Event fraction is the observed fraction with an event during follow-up, not a fixed-horizon risk or population prevalence. Follow-up time is used exactly as supplied.',
       '- The incident-risk cohort and baseline covariates must be defined upstream; this evaluation does not derive absolute risk or choose a prediction horizon.'))
methods_text<-c(paste0('# ',Y,' UKB prediction evaluation'),'',paste('Outcome:',outcome_definition),
 paste('Covariates:',paste(covars,collapse=', ')),paste('Original ancestry label column:',gc),paste('Evaluation:',ci_text),'',
 '## How to read the figures',
 '- Prediction by target ancestry: overall comparison including COJO and the principal combined predictors.',
 '- Combined-score comparison: focuses on how CSx evidence is combined; adds fixed-meta and omits COJO. Repeated predictors have exactly the same estimates as in the first plot.',
 '- Prediction along genetic distance: four panels, PRS-CSx only. a: original ancestry labels in PCA space; b: ancestry-level performance; c: the same people coloured by distance; d: performance in distance bins.',
 '- DiscoDivas versus PRS-CSx is a supplementary comparison placed after the distance figure. Error bars in all figures show 95% bootstrap intervals when enabled.','',
 '## Score definitions',
 '- COJO: ancestry-matched SNP/refA/bJ score (or explicit --pt-effect b); not a P+T grid search.',
 '- PRS-CSx-auto-meta: learned phi and posterior meta-analysis of SNP effects.',
 '- PRS-CSx-fixed-meta: fixed phi and posterior meta-analysis of SNP effects.',
 '- PRS-CSx: four CSx scores jointly regressed on the outcome with covariates within each target ancestry training fold. No phi-grid selection. A separate combination is fit for each ancestry/fold.',
 '- DiscoDivas-tuned: four ancestry anchor models fitted in training folds, followed by distance interpolation and training-only calibration. DiscoDivas-untuned uses the saved score from 2disco.sh.',
 '- Every full model includes the specified covariates. Auto/fixed-meta still receive training-fold outcome calibration.','',
 '## Metrics',metric_methods,'',
 '## Distance and uncertainty',
 paste0('- Euclidean distance to the 1KG EUR reference median, using projected PC1..PC',npc,'. It is not a distance to the discovery-GWAS centre. Colour and x-axis square-root transforms change display spacing only; ticks retain raw distance units.'),
 '- Original groups are taken from --group-col; they are not reassigned using the plotted distance and need not be self-reported ethnicity. OTH/UNASSIGNED labels are retained and do not imply proven admixture.',
 paste0('- Up to ',nbins,' equal-count distance bins within each group, requiring ',minn,' people/bin. For dt/t2e, reduce the bin count until every bin has at least ',min_bin_events,' events/cases and ',min_bin_events,' non-events/controls. Groups without two valid bins are omitted from panel d.'),
 '- Bin boundaries use distances; event counts only determine whether bins are sufficiently populated. Predictions stay fixed, without within-bin refitting or performance-based bin selection.',
 '- distance_performance.tsv saves bin bounds, medians, N, event counts, performance and intervals. Its metric is the same as in the ancestry panels.',
 '- Bootstrap resamples people jointly across methods on fixed OOF predictions. It does not include discovery/model-fitting uncertainty or cluster relatives. Distance-bin bootstrap uses at most 100 resamples.',
 '- Predictions are saved only with --write-predictions TRUE. A successful report has a SUCCESS marker.','',
 '## References',
 '- PRS-CS: https://doi.org/10.1038/s41467-019-09718-5',
 '- PRS-CSx: https://doi.org/10.1038/s41588-022-01054-7',
 '- Four-panel presentation inspired by Ding et al., Fig. 1 (schematic); this report displays empirical group/bin performance: https://doi.org/10.1038/s41586-023-06079-4',
 '- DiscoDivas: https://doi.org/10.1016/j.ajhg.2026.05.006')
writeLines(methods_text,file.path(out,'methods.md'))
escape<-function(x){x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
html_table<-function(z)paste0('<div class="table-wrap"><table><tr>',paste0('<th>',escape(names(z)),'</th>',collapse=''),'</tr>',
 paste(apply(as.data.frame(z),1,function(r)paste0('<tr>',paste0('<td>',escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table></div>')
csx_summary<-copy(performance[method=='PRS-CSx'])
csx_summary<-csx_summary[order(match(target,groups))]
if(type=='ct'){
 overview<-csx_summary[,.(Ancestry=target,N=fmt_n(N),`Partial R²`=fmt(estimate),`Covariates R²`=fmt(baseline_R2),`Full R²`=fmt(full_R2),`Incremental R²`=fmt(delta_R2),RMSE=fmt(RMSE,3))]
 ex<-csx_summary[target=='EUR'];if(!nrow(ex))ex<-head(csx_summary,1)
 metric_explanation<-paste0('<p><strong>How partial R² is calculated.</strong> First predict the phenotype from covariates; then add PRS. ',
  'Both predictions are made outside the fitting fold. Partial R² is the fraction of the baseline model’s squared prediction error removed by adding PRS.</p>',
  '<p class="formula">Partial R² = 1 − SSE(full) / SSE(covariates) = (R²(full) − R²(covariates)) / (1 − R²(covariates))</p>',
  '<p><strong>',escape(ex$target),' PRS-CSx:</strong> (',fmt(ex$full_R2,6),' − ',fmt(ex$baseline_R2,6),') / (1 − ',fmt(ex$baseline_R2,6),') = <strong>',fmt(ex$estimate,4),'</strong>. ',
  'This removes ',fmt(100*ex$estimate,2),'% of the error remaining after covariates; the gain relative to total phenotype variance is ',fmt(ex$delta_R2,4),'.</p>')
 if(!nrow(ex))metric_explanation<-'<p>Partial R² = 1 − SSE(full)/SSE(covariates), using held-out predictions. A worked PRS-CSx example is unavailable because the four-score model could not be evaluated.</p>'
} else if(type=='t2e'){
 overview<-csx_summary[,.(Ancestry=target,N=fmt_n(N),Events=fmt_n(events),`Event fraction`=paste0(fmt(100*P,2),'%'),
   `Full C`=fmt(estimate),`Full C 95% CI`=paste0(fmt(lower95),'–',fmt(upper95)),`Covariates C`=fmt(baseline_C),
   `ΔC from PRS`=fmt(delta_C),`ΔC 95% CI`=paste0(fmt(delta_C_lower95),'–',fmt(delta_C_upper95))) ]
 metric_explanation<-paste0('<p><strong>Time-to-event analysis.</strong> The main metric is Harrell’s C-index for <strong>covariates + PRS</strong>. ',
  'C = 0.5 means chance ranking and C = 1 means perfect ranking of comparable event times. It is not R² or absolute disease risk.</p>',
  '<p>The dashed baseline in the bar plots and black ticks in panel b show covariates-only C. <strong>ΔC = full C − covariates-only C</strong> isolates the improvement in ranking after adding PRS. ',
  'N, event counts and observed event fractions are shown below; event fractions are not fixed-horizon disease risks.</p>')
} else {
 overview<-csx_summary[,.(Ancestry=target,N=fmt_n(N),Cases=fmt_n(events),`Case fraction`=paste0(fmt(100*P,2),'%'),
   `Observed partial R²`=fmt(estimate),AUC=fmt(AUC),`Covariates AUC`=fmt(baseline_AUC),`ΔAUC`=fmt(delta_AUC),Brier=fmt(Brier))]
 metric_explanation<-'<p><strong>Binary outcome.</strong> Main plots use observed-scale partial R². The table also shows logistic AUC and Brier. These R² values are not on the liability scale; the sample case fraction is not population prevalence.</p>'
}
shown<-copy(performance);for(nm in names(shown))if(is.numeric(shown[[nm]]))set(shown,j=nm,value=signif(shown[[nm]],5))
writeLines(paste0('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>',escape(Y),' UKB prediction</title>',
 '<style>body{font:16px system-ui;max-width:1440px;margin:32px auto;padding:0 24px;color:#263446;line-height:1.55}img{width:100%;height:auto;margin:16px 0}.table-wrap{overflow-x:auto}table{border-collapse:collapse;font-size:13px;width:100%}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:right;white-space:nowrap}th{background:#f2f4f7}td:first-child,th:first-child{text-align:left}.metric{background:#f2f6fa;padding:16px 22px;border-radius:8px}.formula{font-family:ui-monospace,monospace}pre{white-space:pre-wrap;font:14px system-ui}a{color:#245f9b}details{margin:18px 0}</style>',
 '<h1>',escape(Y),' | UKB prediction comparison</h1><p><strong>Outcome:</strong> ',escape(outcome_definition),
 ' · <strong>Covariates:</strong> ',escape(paste(covars,collapse=', ')),' · <strong>Evaluation:</strong> ',nfold,' held-out folds.</p>',
 if(length(missing))paste0('<p><strong>Partial report. Unavailable inputs: ',escape(paste(missing,collapse=', ')),'.</strong></p>') else '',
 '<p><a href="plots.pdf">Figures PDF</a> · <a href="performance.tsv">Performance</a> · <a href="distance_performance.tsv">Distance bins</a> · <a href="methods.md">Methods</a> · <a href="cohort.tsv">Cohort</a></p>',
 '<div class="metric">',metric_explanation,'</div><h2>PRS-CSx overview</h2>',html_table(overview),
 paste0('<img src="',names(figures),'.png" alt="',c('Overall prediction benchmark','Score combination strategies','Four-panel ancestry and genetic-distance performance','DiscoDivas difference from PRS-CSx'),'">',collapse=''),
 '<h2>All predictors</h2>',html_table(shown),
 if(length(skips))paste0('<details><summary>Skipped evaluations</summary>',html_table(rbindlist(skips)),'</details>') else '',
 '<details><summary>Evaluation inputs</summary>',html_table(manifest),'</details>',
 '<details><summary>Score definitions</summary>',html_table(model_map),'</details>',
 '<details><summary>Methods and references</summary><pre>',escape(paste(methods_text,collapse='\n')),'</pre></details></html>'),file.path(out,'report.html'))
