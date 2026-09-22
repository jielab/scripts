# Plot/report layer, sourced by Yeval.R after metrics have been saved.
colours <- c(`GRID-tuned`='#278C78',GRID_shared='#76A88C',GRID_posterior='#5D9B91',GRID_matched='#86A9AA',COJO='#7C8798',`PRS-CSx-auto-meta`='#E1A63B',`PRS-CSx`='#3275B4',
 `DiscoDivas-untuned`='#CD8B95',`DiscoDivas-tuned`='#D45260',`PRS-CSx-fixed-meta`='#7755A2',`csx.AFR`='#C6813B',`csx.EAS`='#60A67A',`csx.EUR`='#6593CC',`csx.SAS`='#9E83BD')
ancestry_colours<-c(EUR='#437CB3',AFR='#DBA13B',EAS='#55A589',SAS='#AE6BA5',OTH='#9A9DA4')
metric_label<-switch(type,ct='Partial R²',dt='Observed-scale partial R²',t2e='Harrell C-index')
ci_text<-paste0(nfold,'-fold out-of-fold predictions; ',if(nboot>0)'paired bootstrap 95% intervals' else 'no confidence intervals (bootstrap=0)')
target_label<-setNames(paste0('Target = ',pops,'\nN = ',vapply(pops,function(g){v<-performance[target==g]$N;if(length(v))format(v[1],big.mark=',',scientific=FALSE) else 'unavailable'},character(1))),pops)
label_methods<-function(z){z<-sub('PRS-CSx-auto-meta','PRS-CSx\nauto-meta',z,fixed=TRUE);z<-sub('PRS-CSx-fixed-meta','PRS-CSx\nfixed-meta',z,fixed=TRUE);sub('DiscoDivas-','DiscoDivas\n',z,fixed=TRUE)}
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
 'Binary R² uses held-out linear probability models on the observed scale; AUC/Brier use logistic models.' else
 'Concordance uses comparable pairs within each held-out fold; censoring is retained.'
p1<-make_comparison(main_methods,paste(Y,'| Prediction by target ancestry'),primary_caption)
p2<-make_comparison(c('PRS-CSx-auto-meta','PRS-CSx-fixed-meta','PRS-CSx',if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned'),paste(Y,'| Combined-score comparison'),
 'auto-meta: learned phi; fixed-meta: fixed phi. PRS-CSx combines four scores within the training fold.')
if(nrow(comparison)) {
 comparison[,target:=factor(target,levels=groups)]
 p3<-ggplot(comparison,aes(target,difference))+geom_hline(yintercept=0,linetype=2,colour='#777777')+
  geom_errorbar(aes(ymin=lower95,ymax=upper95),width=.12,colour=colours[disco_method],linewidth=.7)+
  geom_point(size=3.5,colour=colours[disco_method])+plot_theme+
  labs(title=paste(Y,'|',disco_method,'versus PRS-CSx'),subtitle=ci_text,x='Target ancestry',
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
 paste('Folds:',nfold,'; seed:',seed,'; bootstrap:',nboot),'',
 '## Score definitions',
 '- COJO: ancestry-matched SNP/refA/bJ score (or explicit --pt-effect b). This is not a clumping/P-threshold grid search.',
 '- PRS-CSx-auto-meta: PRS-CSx learned phi plus official posterior meta-analysis; not independently fitted PRS-CS-mult.',
 '- PRS-CSx-fixed-meta: posterior meta-analysis using the configured fixed phi.',
 '- PRS-CSx: joint regression of four training-standardized csx.AFR/EAS/EUR/SAS scores within the target ancestry. No additional phi grid tuning.',
 '- DiscoDivas-untuned: saved 2disco.sh output using official geometry and raw ancestry-specific scores. Its PCA residualization uses the input scoring cohort (transductive, no outcomes). It is a diagnostic comparator.',
 '- DiscoDivas-tuned: each outer training fold fits four ancestry-specific C+4PRS anchor models. Only their genetic contributions enter interpolation. Centers are PC medians of the actual anchor training people; PC residualization/scaling uses a balanced subset of training anchors. Calibration is also fit outside the test fold. Quality factors are fixed, not phenotype-selected. This adapts DiscoDivas to CSx inputs; it is not a reproduction of the paper’s LDpred2 experiment.',
 '', '## Metrics and uncertainty',
 '- Continuous: held-out partial R² = 1-SSE(C+PRS)/SSE(C). Baseline and full total R² and RMSE are additional metrics; negative estimates are retained.',
 '- Binary: the main R² is observed-scale partial R² from held-out linear probability models. Logistic models supply AUC/Brier. The old direct Lee conversion of covariate-adjusted partial R² has been removed: the residual denominator is not the raw binary variance P(1-P). These outputs must not be labelled liability R².',
 '- Prevalence K is retained only as descriptive/context information. Cohort K is not general-population or lifetime prevalence.',
 '- Default t2dm dt explicitly derives baseline ICD10 status from Yr2e/Yt2e. Use --phenotype-col for another definition. Survival uses the provided positive follow-up time and 0/1 event; prevalent-case exclusions must be correct upstream.',
 '- Survival: fold-wise comparable-pair-weighted Harrell C; folds do not share an arbitrary Cox baseline. No liability conversion.',
 '- Bootstrap resamples the same individuals for every method on fixed OOF predictions. It does not include GWAS/model-fitting uncertainty and does not cluster relatives. Use an unrelated evaluation cohort.',
 '- Missing scores never become zero. Models within an ancestry share complete cases and folds. UNASSIGNED means missing ancestry labels, not proven admixture.',
 '- Predictions are saved only with --write-predictions TRUE. Successful reports have a SUCCESS marker. No real-UKB performance improvement is implied by software tests.',
 '', '## Distance analysis',
 paste0('- Reference distances use PC1..PC',npc,'. They are descriptive distances to 1KG centers; tuned Disco uses its own training centers.'),
 '- Distance bins evaluate already-fixed OOF predictions; they do not refit/select models or estimate individual heritability.',
 '', '## References',
 '- PRS-CSx: https://doi.org/10.1038/s41588-022-01054-7',
 '- DiscoDivas: https://doi.org/10.1016/j.ajhg.2026.05.006',
 '- PLINK scoring: https://www.cog-genomics.org/plink/2.0/score'
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
