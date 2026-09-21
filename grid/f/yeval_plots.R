# Plot/report layer, sourced by Yeval.R after metrics have been saved.
colours <- c(PT='#7C8798',`PRS-CS-multi`='#E1A63B',`PRS-CSX`='#3275B4',
 DiscoDivas='#D45260',`CSx-meta`='#7755A2',`csx.AFR`='#C6813B',`csx.EAS`='#60A67A',`csx.EUR`='#6593CC',`csx.SAS`='#9E83BD')
ancestry_colours<-c(EUR='#437CB3',AFR='#DBA13B',EAS='#55A589',SAS='#AE6BA5',OTH='#9A9DA4')
metric_label<-switch(type,ct='Prediction R²',dt='Liability R²',t2e='Harrell C-index')
ci_text<-paste0(nfold,'-fold out-of-fold predictions; paired bootstrap 95% intervals')
target_label<-setNames(paste0('Target = ',pops,'\nN = ',vapply(pops,function(g){v<-performance[target==g]$N;if(length(v))format(v[1],big.mark=',',scientific=FALSE) else 'unavailable'},character(1))),pops)
label_methods<-function(z)sub('PRS-CS-multi','PRS-CS-\nmulti',z,fixed=TRUE)
plot_theme<-theme_classic(base_size=12)+theme(legend.position='bottom',strip.background=element_rect(fill='#F2F4F7',colour=NA),
 strip.text=element_text(face='bold',margin=margin(8,5,8,5)),plot.title=element_text(face='bold',size=17),
 plot.subtitle=element_text(colour='#526071',size=10),axis.text.x=element_text(angle=35,hjust=1),
 plot.caption=element_text(hjust=0,colour='#626A76',size=9),panel.spacing=grid::unit(1.2,'lines'))
make_comparison<-function(methods,title,caption) {
 grid<-CJ(target=pops,method=methods,unique=TRUE)
 z<-merge(grid,performance[target%in%pops],by=c('target','method'),all.x=TRUE)
 z[,method:=factor(method,levels=methods)];z[,target:=factor(target,levels=pops)]
 ggplot(z,aes(method,estimate,fill=method))+geom_hline(yintercept=if(type=='t2e').5 else 0,colour='#AEB7C2',linewidth=.3)+
  geom_col(width=.62,na.rm=TRUE,alpha=.85)+geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.16,linewidth=.55,na.rm=TRUE)+
  geom_point(shape=21,size=2.1,na.rm=TRUE)+
  geom_text(data=z[is.na(estimate)],aes(x=method,y=0,label='Unavailable'),inherit.aes=FALSE,angle=90,hjust=0,size=2.6,colour='#777777')+
  facet_grid(.~target,labeller=labeller(target=target_label),drop=FALSE)+
  scale_fill_manual(values=colours,drop=FALSE)+scale_x_discrete(labels=label_methods,drop=FALSE)+
  scale_y_continuous(expand=expansion(mult=c(.06,.13)))+plot_theme+theme(legend.position='none')+
  labs(title=title,subtitle=paste('UKB |',ci_text),x=NULL,y=metric_label,caption=caption)
}
primary_caption<-if(type=='ct')'R² = 1 - SSE(covariates + score) / SSE(covariates), using held-out predictions.' else if(type=='dt')
 'Liability R²: Lee et al. transformation of held-out linear-model partial R²; K and sample P are recorded separately.' else
 'Concordance uses comparable pairs within each held-out fold; censoring is retained.'
p1<-make_comparison(main_methods,paste(Y,'| Prediction by target ancestry'),primary_caption)
p2<-make_comparison(c('PRS-CS-multi','CSx-meta','PRS-CSX','DiscoDivas'),paste(Y,'| Combined-score comparison'),
 'PRS-CS-multi = csx.auto; CSx-meta = csx.meta. Labels follow the requested comparison; see methods.md for exact provenance.')
if(nrow(comparison)) {
 comparison[,target:=factor(target,levels=groups)]
 p3<-ggplot(comparison,aes(target,difference))+geom_hline(yintercept=0,linetype=2,colour='#777777')+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.12,colour=colours['DiscoDivas'],linewidth=.7)+
  geom_point(size=3.5,colour=colours['DiscoDivas'])+plot_theme+
  labs(title=paste(Y,'| DiscoDivas versus PRS-CSX'),subtitle=ci_text,x='Target ancestry',
       y=paste0('Difference in ',metric_label),caption='Positive values favour DiscoDivas. Same participants, folds and bootstrap resamples for both methods.')
} else p3<-ggplot()+theme_void()+labs(title='Paired comparison unavailable')
# Show all genetically labelled samples, including OTH, without making up an
# individual R2: performance is measured only in held-out groups/distance bins.
set.seed(seed+1L)
landscape<-d[!is.na(target),.SD[sample.int(.N,min(.N,2000L))],by=target]
landscape[,target:=factor(target,levels=groups)]
pa<-ggplot(landscape,aes(proj_PC1,proj_PC2,colour=target))+geom_point(size=.55,alpha=.5)+
 scale_colour_manual(values=ancestry_colours,na.value='grey70')+plot_theme+theme(axis.text.x=element_text())+
 labs(title='a  Genetic ancestry',x='Reference-projected PC1',y='Reference-projected PC2',colour='Target')
pb<-ggplot(landscape,aes(proj_PC1,proj_PC2,colour=distance.EUR))+geom_point(size=.55,alpha=.55)+
 scale_colour_viridis_c(option='C')+plot_theme+theme(axis.text.x=element_text())+
 labs(title='b  Continuous genetic distance',x='Reference-projected PC1',y='Reference-projected PC2',colour='Distance to EUR')
pcp<-ggplot(landscape,aes(distance.EUR,colour=target,fill=target))+geom_density(alpha=.08,linewidth=.7)+
 scale_colour_manual(values=ancestry_colours,na.value='grey70')+scale_fill_manual(values=ancestry_colours,na.value='grey70')+
 plot_theme+theme(axis.text.x=element_text())+labs(title='c  Distance distributions within ancestry',x=paste0('Euclidean distance to EUR reference center (',npc,' PCs)'),y='Density',colour='Target',fill='Target')
p4<-(pa+pb)/pcp+plot_annotation(title=paste(Y,'| The genetic ancestry continuum'),
 caption='PCA/density display: at most 2,000 people per ancestry. Distances use the shared DiscoDivas reference PC space.\nEUR is a reference center, not the discovery-GWAS center; these plots do not estimate individual-level theoretical accuracy.',theme=theme(plot.title=element_text(face='bold',size=17),plot.caption=element_text(hjust=0,size=9)))
if(nrow(distperf)) {
 distperf[,target:=factor(target,levels=groups)]
 distperf[,panel:=paste(target,ifelse(distance_axis=='distance.EUR','EUR reference','Nearest reference'),sep=' | ')]
 distperf[,panel:=factor(panel,levels=unlist(lapply(c('EUR reference','Nearest reference'),function(a)paste(groups,a,sep=' | '))))]
 p5<-ggplot(distperf,aes(distance,estimate,colour=method,group=method))+geom_hline(yintercept=if(type=='t2e').5 else 0,colour='#CCCCCC')+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=0,alpha=.45,linewidth=.45)+geom_line(linewidth=.65)+geom_point(size=1.9)+
  facet_wrap(~panel,ncol=length(unique(distperf$target)),scales='free_x')+scale_colour_manual(values=colours)+plot_theme+theme(axis.text.x=element_text(angle=25,hjust=1))+
  labs(title=paste(Y,'| Prediction along genetic distance'),subtitle=paste(nbins,'distance quantiles within each ancestry; points at bin medians'),x='Euclidean genetic distance',y=metric_label,colour=NULL,
       caption='Each point evaluates fixed out-of-fold predictions on a distance bin; no outcome-based binning or refitting.\nSparse bins are omitted. These are empirical local performance estimates, not individual R².')
} else p5<-ggplot()+theme_void()+labs(title='Distance performance unavailable: insufficient samples/cases in bins')
figures<-list(comparison=p1,combined_scores=p2,paired_improvement=p3,genetic_landscape=p4,distance_performance=p5)
pdf(file.path(out,'plots.pdf'),width=13,height=6.8,useDingbats=FALSE)
for(p in figures)print(p)
dev.off()
for(nm in names(figures)) {
 ggsave(file.path(out,paste0(nm,'.png')),figures[[nm]],width=13,height=if(nm=='genetic_landscape')8 else 6.8,dpi=160,bg='white')
}
methods_text<-c(
 paste0('# ',Y,' UKB prediction evaluation'),'',
 paste('Outcome:',outcome_definition),paste('Covariates:',paste(covars,collapse=', ')),
 paste('Folds:',nfold,'; seed:',seed,'; bootstrap:',nboot),
 '', '## Score definitions',
 '- PT: ancestry-matched .jma.cojo SNP/refA/bJ scores (or the explicitly selected --pt-effect). This is the requested COJO-based comparator, not the manuscript’s clumping/P-threshold grid search. Match rsID, Chr/bp and effect allele; duplicate rsIDs are scored only when exactly one genotype record matches. Ambiguous records are excluded and documented in pt.pgs.gz.matches.tsv beside the PT score cache. Chromosomal SCORE_SUM values are added; averages are never added. t2dm.AFA maps to AFR.',
 '- PRS-CS-multi: display label requested for csx.auto. The current upstream pipeline generates csx.auto with PRS-CSx automatic phi and posterior meta-analysis. It is not the independently fitted PRS-CS-mult algorithm in the 2022 paper.',
 '- PRS-CSX: simultaneous regression on training-fold-standardized csx.AFR, csx.EAS, csx.EUR and csx.SAS, separately within each target ancestry. Standardization and score coefficients are fitted using the training fold only. Only the available fixed-phi scores are used; no additional phi grid selection is claimed.',
 '- DiscoDivas: the existing saved disco.pgs.gz from 2disco.sh. Its current pipeline uses the official distance-matrix correction and PC-residualized ancestry scores, 1000 Genomes reference medians, and quality factors of one. It does not fine-tune the input PRS separately on ancestry-specific phenotype training sets as in the full published pipeline. Its phenotype calibration is fitted in each training fold here. The comparison evaluates this saved implementation; it does not claim a full reproduction of the manuscript.',
 '- CSx-meta and individual csx.* scores are also evaluated with the same folds and complete-case cohort. csx.meta is shown in the second figure and all individual scores are in performance.tsv.',
 '', '## Prediction metrics and uncertainty',
 '- Continuous: primary metric is covariate-adjusted out-of-fold partial R² = 1 - SSE(full)/SSE(covariates). Full and baseline total R² and RMSE are also saved. Negative held-out values are retained.',
 '- Binary: logistic models provide OOF probabilities, AUC and Brier scores. Separate OOF linear probability models provide observed-scale partial R² for the Lee et al. (2012) liability transformation. Logistic pseudo-R² and AUC are not relabelled as liability R².',
 '- Liability transformation: t=qnorm(1-K), z=dnorm(t), C=K²(1-K)²/[z²P(1-P)], a=(z/K)(P-K)/(1-K); R²_liability=C R²_observed/[1+C a(a-t) R²_observed]. P is the evaluated ancestry’s sample case fraction. Nonpositive denominators are undefined. Negative predictive R² is retained as a diagnostic extension, not interpreted as a negative biological variance.',
 '- --prevalence cohort uses the full available phenotype cohort within each genetic ancestry before score/covariate filtering, as explicitly requested. This is a UKB cohort-based working assumption, not a representative general-population or lifetime prevalence. K, source, numerator and denominator are shown in the report’s Prevalence assumptions table for binary outcomes. Supply externally justified K for population inference.',
 '- The default derived t2dm binary outcome is dated baseline ICD10 T2D. Incident cases are non-cases at baseline. Unknown dates/invalid outcomes remain missing. This source may under-ascertain diabetes; it is not the broader algorithmic/HbA1c definition used by some UKB publications.',
 '- Survival: fold-wise comparable-pair-weighted Harrell C with censoring. No liability conversion is applied to incident-event indicators.',
 '- Confidence intervals: paired resampling of individuals on fixed OOF predictions, stratified by case/event for binary/survival outcomes. They are conditional on fitted models and do not include discovery-GWAS or model-fitting uncertainty. Relatives are not clustered; independent discovery and unrelated samples must be established upstream.',
 '- Model-specific missing values never become zero. The available methods share complete cases within each target. A partial run labels missing methods explicitly. All phenotype-informed score weights are learned outside the evaluated fold.',
 '', '## Genetic distance analysis',
 paste0('- Euclidean distances use PC1..PC',npc,' in the same reference-projected space and the four reference centers used by the saved DiscoDivas pipeline.'),
 '- Both distance to EUR reference and distance to the nearest reference are examined. Neither is presented as distance to the actual GWAS discovery sample, whose center is unavailable.',
 '- Quantile bins are defined by genetic distance within each target. Each bin evaluates the existing OOF predictions. Local R²/CI uses the bin’s residual errors; binary liability conversions keep the target-level K and use the bin case fraction. Bins with too few cases/controls are omitted.',
 '- The Nature 2023 figure is a conceptual illustration. Here PCA and distance distributions are empirical; binned performance is not the paper’s theoretical individual r_i² estimator, which needs additional discovery genotype/heritability information.',
 '', '## References',
 '- PRS-CSx: https://doi.org/10.1038/s41588-022-01054-7',
 '- DiscoDivas: https://doi.org/10.1016/j.ajhg.2026.05.006',
 '- Genetic ancestry continuum: https://doi.org/10.1038/s41586-023-06079-4',
 '- Liability R²: Lee et al. (2012), https://doi.org/10.1002/gepi.21614 ; implementation reference https://cnsgenomics.com/data/teaching/GNGWS23/module5/Practical2_accuracy.html',
 '- UKB prevalence context: https://pmc.ncbi.nlm.nih.gov/articles/PMC6936483/ reports 4.4% in White European and 16.4% in South Asian participants under a different baseline T2D definition. These rates are not silently substituted for genetic-ancestry-specific rates in this cohort.'
)
writeLines(methods_text,file.path(out,'methods.md'))
escape<-function(x) {x<-gsub('&','&amp;',as.character(x),fixed=TRUE);x<-gsub('<','&lt;',x,fixed=TRUE);gsub('>','&gt;',x,fixed=TRUE)}
html_table<-function(z)paste0('<table><tr>',paste0('<th>',escape(names(z)),'</th>',collapse=''),'</tr>',paste(apply(as.data.frame(z),1,function(r)paste0('<tr>',paste0('<td>',escape(r),'</td>',collapse=''),'</tr>')),collapse=''),'</table>')
shown<-copy(performance);for(nm in names(shown))if(is.numeric(shown[[nm]]))set(shown,j=nm,value=signif(shown[[nm]],4))
writeLines(paste0('<!doctype html><meta charset="utf-8"><title>',escape(Y),' UKB prediction</title>',
 '<style>body{font:16px system-ui;max-width:1400px;margin:36px auto;padding:0 24px;color:#263446}img{width:100%}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border-bottom:1px solid #ddd}pre{white-space:pre-wrap}</style>',
 '<h1>',escape(Y),' | UKB prediction comparison</h1><p>',escape(ci_text),'. ',escape(primary_caption),'</p>',
 if(length(missing))paste0('<p><strong>PARTIAL REPORT. Unavailable: ',escape(paste(missing,collapse=', ')),'.</strong></p>') else '',
 '<p><a href="plots.pdf">All figures PDF</a> · <a href="performance.tsv">Performance table</a> · <a href="methods.md">Methods and references</a> · <a href="cohort.tsv">Cohort audit</a></p>',
 if(type=='dt')paste0('<h2>Prevalence assumptions</h2>',html_table(prevalence)) else '',
 paste0('<img src="',names(figures),'.png" alt="',names(figures),'">',collapse=''),
 '<h2>Performance</h2>',html_table(shown),
 if(length(skips))paste0('<h2>Skipped evaluations</h2>',html_table(rbindlist(skips))) else '',
 '<details><summary>Evaluation inputs</summary>',html_table(manifest),'</details>',
 '<details><summary>Score definitions</summary>',html_table(model_map),'</details>',
 '<h2>Methods</h2><pre>',escape(paste(methods_text,collapse='\n')),'</pre>'),file.path(out,'report.html'))
