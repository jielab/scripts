# Publication-style companions. Always plot this run's results, never paper data.
le8_mock_volcano<-function(a,layer='protein') {
  z<-a|>filter(is.finite(beta),is.finite(p.value))|>mutate(HR=exp(beta),significance=case_when(p.value<.05/nrow(a)&beta>0~'Risk',p.value<.05/nrow(a)&beta<0~'Protective',TRUE~'Not significant'))
  lab<-z|>arrange(p.value)|>slice_head(n=15)
  ggplot(z,aes(HR,-log10(pmax(p.value,1e-300)),color=significance))+geom_point(size=1.5,alpha=.7)+
    geom_hline(yintercept=-log10(.05/nrow(a)),linetype=2)+geom_vline(xintercept=1,color='grey75')+
    ggrepel::geom_text_repel(data=lab,aes(label=term),seed=SEED,max.overlaps=Inf,size=3)+
    scale_color_manual(values=c(Risk='#CE4965',Protective='#7463AB',`Not significant`='grey80'))+
    labs(title=paste('A. Plasma',if(layer=='protein')'proteins'else'metabolites','and incident',Y),
      x=paste('Hazard ratio per 1-SD',if(layer=='protein')'protein'else'metabolite'),y=expression(-log[10](P)),color=NULL)+theme_5c(10)
}
le8_mock_trajectories<-function(dat,a,covars,tvar,evar,layer='protein') {
  features<-a|>filter(is.finite(p.value),p.value<.05/nrow(a))|>arrange(p.value)|>pull(term)|>intersect(names(dat))
  empty<-list(trajectories=tibble(),clusters=tibble(),status=tibble(status='No significant biomarkers or insufficient matched observations',layer=layer))
  if(!length(features))return(empty)
  if(!requireNamespace('MatchIt',quietly=TRUE)){empty$status$status<-'Install MatchIt for matched-control trajectories';return(empty)}
  z<-as.data.frame(dat[,unique(c('eid',covars,tvar,evar,features)),drop=FALSE])
  z<-z[complete.cases(z[,c(covars,tvar,evar),drop=FALSE])&z[[tvar]]>0,,drop=FALSE]
  z$.case<-z[[evar]];if(sum(z$.case==1)<20||sum(z$.case==0)<20)return(empty)
  set.seed(SEED)
  match_covars<-covars[vapply(z[covars],function(v)length(unique(v))>1L,logical(1))]
  form<-if(length(match_covars))reformulate(match_covars,response='.case')else as.formula('.case ~ 1')
  m<-tryCatch(MatchIt::matchit(form,z,method='nearest',ratio=max(1L,min(10L,floor(sum(z$.case==0)/sum(z$.case==1)))),replace=FALSE),error=identity)
  if(inherits(m,'error')){empty$status$status<-paste('Matching unavailable:',conditionMessage(m));return(empty)}
  md<-MatchIt::match.data(m);ca<-md[md$.case==1,,drop=FALSE];co<-md[md$.case==0,,drop=FALSE]
  grid<-seq(-15,0,by=.1)
  tr<-bind_rows(lapply(features,function(x){
    mu<-mean(co[[x]],na.rm=TRUE);s<-sd(co[[x]],na.rm=TRUE)
    if(!is.finite(s)||s<=0)return(tibble())
    d<-data.frame(years=-ca[[tvar]],z=(ca[[x]]-mu)/s);d<-d[complete.cases(d)&d$years>=-15,,drop=FALSE]
    if(nrow(d)<20||length(unique(d$years))<5)return(tibble())
    # Only fitted means are used; omit unused standard-error calculations.
    fit<-tryCatch(loess(z~years,d,span=.75,control=loess.control(statistics='none')),error=function(e)NULL);if(is.null(fit))return(tibble())
    tibble(feature=x,years=grid,z=as.numeric(predict(fit,data.frame(years=grid))))
  }))
  if(!nrow(tr))return(empty)
  eligible<-tr|>group_by(feature)|>summarise(cross=any(abs(z)>.25,na.rm=TRUE),.groups='drop')|>filter(cross)|>pull(feature)
  mat<-tr|>filter(feature%in%eligible)|>pivot_wider(names_from=years,values_from=z)
  cl<-tibble(feature=character(),cluster=integer())
  if(nrow(mat)>=3){v<-as.matrix(mat[,-1]);ok<-colSums(is.finite(v))==nrow(v)
    if(sum(ok)>5){hc<-hclust(dist(v[,ok,drop=FALSE]),method='ward.D2');cl<-tibble(feature=mat$feature,cluster=unname(cutree(hc,k=3)))}
  }
  list(trajectories=tr,clusters=cl,status=tibble(status='ok',layer=layer,matching_covariates=paste(match_covars,collapse=','),cases=nrow(ca),controls=nrow(co),biomarkers=length(features),
    proteins=if(layer=='protein')length(features)else NA_integer_,metabolites=if(layer=='metabolite')length(features)else NA_integer_,
    interpretation='One baseline sample per person; diagnosis-timed LOESS, not longitudinal within-person change',reference='Mean/SD of matched controls; absolute Z > 0.25 for clustering'))
}
le8_plot_mock2<-function(z,outdir,layer='protein') {
  tr<-z$trajectories;cl<-z$clusters
  omic<-if(layer=='protein')'protein'else'metabolite'
  write_xlsx2(z,le8_artifact_path('c1.Fig7.profiles.xlsx',outdir))
  if(!nrow(tr)){save_plot(blank_plot(paste('Diagnosis-timed',omic,'trajectories'),z$status$status[1]),'c1.Fig7.diagnosis_timed_profiles.png',14,9,outdir=outdir);return()}
  ord<-if(nrow(cl))c(cl$feature[order(cl$cluster)],setdiff(unique(tr$feature),cl$feature))else unique(tr$feature)
  h<-tr|>mutate(feature=factor(feature,levels=rev(ord)),display=ifelse(abs(z)>.25,z,0))
  p<-ggplot(h,aes(years,feature,fill=display))+geom_tile()+scale_fill_gradient2(low='#3878B9',mid='white',high='#CF4050',midpoint=0,na.value='grey90')+
    guides(fill=guide_colourbar(barwidth=grid::unit(6,'cm'),barheight=grid::unit(.3,'cm'),title.position='top'))+
    labs(title=paste0('A. ',tools::toTitleCase(omic),' levels before diagnosis'),x='Years before diagnosis',y=NULL,fill='Z score')+theme_5c(8)+theme(panel.grid=element_blank())
  if(length(ord)>150)p<-p+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(y=paste(length(ord),paste0(omic,'s; row names in CSV')))
  lines<-lapply(1:3,function(k){d<-inner_join(tr,cl,by='feature')|>filter(cluster==k)
    if(!nrow(d))return(blank_plot(paste('Cluster',k),'Fewer than three estimable clusters'))
    avg<-d|>group_by(years)|>summarise(z=mean(z,na.rm=TRUE),.groups='drop')
    ggplot(d,aes(years,z,group=feature))+geom_line(color=c('#4C8CB5','#D87463','#7C67A4')[k],alpha=.3,na.rm=TRUE)+
      geom_line(data=avg,aes(group=1),linewidth=1,color='black',na.rm=TRUE)+geom_hline(yintercept=0,linetype=2,color='grey65')+
      labs(title=paste0(LETTERS[k+1],'. Cluster ',k,' (n=',n_distinct(d$feature),')'),x='Years before diagnosis',y='Z score')+theme_5c(9)
  })
  caption<-paste0(z$status$interpretation[1],'. Heatmap: ',length(ord),' incident-associated ',omic,'s; ',nrow(cl),
    ' enter Ward clustering (k = 3; |Z| > 0.25 at any supported time). Thin lines: individual ',omic,'s; black: cluster mean. White heatmap cells: |Z| ≤ 0.25; gray: unsupported LOESS range.')
  save_plot(wrap_plots(c(list(p),lines),design='AB\nAC\nAD',widths=c(1.3,1),heights=c(1,1,1))+plot_annotation(caption=caption),
    'c1.Fig7.diagnosis_timed_profiles.png',22,12,outdir=outdir)
}
le8_mock_bar<-function(z,title) {
  if(!nrow(z))return(blank_plot(title,'No FDR-significant terms / annotation unavailable'))
  z<-z|>filter(is.finite(adjusted_p),adjusted_p<.05)|>arrange(adjusted_p)|>slice_head(n=12)
  if(!nrow(z))return(blank_plot(title,'No FDR-significant terms'))
  ggplot(z,aes(-log10(pmax(adjusted_p,1e-300)),reorder(term_name,-adjusted_p),fill=source))+geom_col(width=.7)+
    scale_fill_manual(values=c('GO:BP'='#BC5366','GO:CC'='#286B8B','GO:MF'='#6D629A',KEGG='#388B78',MGI='#C58D37',TF='#687784'))+
    scale_y_discrete(labels=function(x)stringr::str_wrap(x,45))+labs(title=title,x=expression(-log[10]('FDR')),y=NULL,fill=NULL)+theme_5c(8)+theme(axis.text.y=element_text(size=8))
}
le8_mock_network<-function(edges,title,values=NULL) {
  if(!nrow(edges))return(blank_plot(title,'No supported edges / annotation unavailable'))
  edges<-edges|>distinct(from,to,.keep_all=TRUE);g<-igraph::graph_from_data_frame(as.data.frame(edges[,c('from','to')]),directed=FALSE)
  set.seed(SEED);xy<-igraph::layout_with_fr(g);n<-tibble(node=igraph::V(g)$name,x=xy[,1],y=xy[,2],value=as.numeric(igraph::degree(g)))
  if(!is.null(values))n$value<-values[n$node]
  label_nodes<-names(sort(igraph::degree(g),decreasing=TRUE))[seq_len(min(18,nrow(n)))]
  n$label<-ifelse(n$node%in%label_nodes,n$node,'')
  if(!'score'%in%names(edges))edges$score<-1
  e<-edges|>left_join(n|>select(from=node,x,y),by='from')|>left_join(n|>select(to=node,xend=x,yend=y),by='to')
  ggplot()+geom_segment(data=e,aes(x,y,xend=xend,yend=yend,linewidth=score),color='grey75')+
    geom_point(data=n,aes(x,y,color=value),size=3)+ggrepel::geom_text_repel(data=n,aes(x,y,label=label),max.overlaps=Inf,seed=SEED,size=2.3)+
    scale_linewidth_continuous(range=c(.2,1.1),guide='none')+scale_color_gradient(low='#AEC9DB',high='#9E2847',na.value='grey60')+labs(title=title,subtitle='Labels: up to 18 highest-degree nodes; all edges retained',color=if(is.null(values))'Degree' else '-log10(P)')+theme_void()+theme(plot.title=element_text(face='bold'))
}
le8_mock_enrich<-function(a,enrich,outdir,enrich_prev=tibble()) {
  rawdir<-le8_job_dir(outdir,'c1_correlate');genes<-a$term[a$p.value<.05/nrow(a)&is.finite(a$p.value)]
  sf<-file.path(rawdir,'c1.mock_gene_universe.csv');write_raw_csv(a|>transmute(gene=term,selected=term%in%genes),'c1.mock_gene_universe.csv',rawdir)
  cached<-all(file.exists(file.path(rawdir,c('c1.mock_function_terms.csv','c1.mock_tf_edges.csv','c1.mock_ppi_edges.csv'))))
  # Keep public annotation downloads alongside this outcome's C1 outputs.
  py<-Sys.which('python3');status<-if(LE8_REUSE_RESULTS&&cached)0L else if(nzchar(py))system2(py,shQuote(c(file.path(Sys.getenv('LE8_FDIR'),'mock_annotations.py'),sf,rawdir,file.path(rawdir,'le8_annotations'))),timeout=600)else 1L
  rd<-function(n){f<-file.path(rawdir,n);if(status==0L&&file.exists(f))as_tibble(data.table::fread(f))else tibble()}
  terms<-rd('c1.mock_function_terms.csv');tf<-rd('c1.mock_tf_edges.csv');ppi<-rd('c1.mock_ppi_edges.csv')
  if(!'source'%in%names(terms))terms<-tibble(source=character(),term_name=character(),adjusted_p=numeric())
  go<-enrich|>filter(source%in%c('GO:BP','GO:CC','GO:MF'))
  label_file<-file.path(Sys.getenv('LE8_FDIR'),'assets/go_terms.tsv')
  if(nrow(go)&&file.exists(label_file)){
    labels<-data.table::fread(label_file);ix<-match(go$term_name,labels$term)
    hit<-!is.na(ix);go$term_name[hit]<-labels$TERM[ix[hit]]
  }
  pathway<-bind_rows(go,terms|>filter(source=='KEGG'))|>group_by(source)|>slice_min(adjusted_p,n=3,with_ties=FALSE)|>ungroup()
  values<-setNames(-log10(pmax(a$p.value,1e-300)),a$term)
  # Six readable panels integrate prevalent and incident enrichment with annotations.
  gp<-enrich_prev
  if(nrow(gp)&&file.exists(label_file)) {
    labels<-data.table::fread(label_file);ix<-match(gp$term_name,labels$term)
    hit<-!is.na(ix);gp$term_name[hit]<-labels$TERM[ix[hit]]
  }
  top_sources<-function(d) if(nrow(d)) d|>group_by(source)|>slice_min(adjusted_p,n=3,with_ties=FALSE)|>ungroup() else d
  panels<-list(le8_mock_bar(top_sources(go),'A. Incident GO enrichment'),
    le8_mock_bar(top_sources(gp),'B. Prevalent enrichment'),
    le8_mock_bar(terms|>filter(source=='KEGG'),'C. Incident pathways'),
    le8_mock_bar(terms|>filter(source=='MGI'),'D. Mammalian phenotype'),
    le8_mock_network(tf,'E. Transcription-factor targets',values),
    le8_mock_network(ppi,'F. Physical protein interactions'))
  save_plot(wrap_plots(panels,ncol=2)+plot_annotation(caption='Bonferroni-selected associations; library-wise FDR and assayed background. Networks are annotation support, not causal evidence.'),
    'c1.Fig8.enrich_sig.png',20,18,outdir=outdir)
  write_xlsx2(list(incident=enrich,prevalent=enrich_prev,annotations=terms,TF_edges=tf,PPI_edges=ppi,annotation_status=data.frame(exit_code=status)),le8_artifact_path('c1.Fig8.enrichment.xlsx',outdir))
  list(terms=terms,tf_edges=tf,ppi_edges=ppi,annotation_exit=status)
}
le8_mock_c1<-function(dat,a,enrich,layer,covars,tvar,evar,outdir,enrich_prev=tibble()) {
  save_plot(le8_mock_volcano(a,layer),'c1.Fig17.measured_volcano.png',10,7,outdir=outdir)
  z<-le8_mock_trajectories(dat,a,covars,tvar,evar,layer)
  for(n in names(z))write_raw_csv(z[[n]],paste0('c1.mock_',n,'.csv'),le8_job_dir(outdir,'c1_correlate'))
  le8_plot_mock2(z,outdir,layer)
  e<-if(layer=='protein')le8_mock_enrich(a,enrich,outdir,enrich_prev)else list()
  list(temporal=z,functional=e)
}

le8_mock_c5<-function(obj,outdir) {
  imp<-as_tibble(obj$scores$Yu$importance%||%tibble())
  if(!nrow(imp)&&!is.null(obj$scores$Yu$fit))imp<-tryCatch(as_tibble(lightgbm::lgb.importance(obj$scores$Yu$fit)),error=function(e)tibble())
  pa<-if(nrow(imp))imp|>arrange(desc(Gain))|>slice_head(n=15)|>ggplot(aes(reorder(Feature,Gain),Gain))+geom_col(fill='#BA5272')+coord_flip()+labs(title='A. Training LightGBM feature importance',x=NULL,y='Split gain')+theme_5c(8) else blank_plot('A. Feature importance','LightGBM importance unavailable')
  seq<-obj$sequential$log%||%tibble()
  curve<-if(nrow(seq))ggplot(seq,aes(step,auc))+geom_line(color='#356C96')+geom_point()+labs(title='B. Training-only forward selection',x='Number of predictors',y='Inner holdout AUC')+theme_5c(9) else blank_plot('B. Forward selection','No estimable training selection curve')
  pred<-obj$prediction;label<-unique(pred$biom_set[grepl('LightGBM',pred$biom_set)])
  if(!length(label))label<-unique(pred$biom_set)[1]
  d<-pred|>filter(biom_set==label[1]);if('split'%in%names(d))d<-d|>filter(split=='validation')
  full_h<-as.numeric(Sys.getenv('C5_MOCK_HORIZON','15'));windows<-list(c(0,full_h),c(0,10),c(10,full_h))
  curves<-list();metrics<-list();plots<-lapply(seq_along(windows),function(i){
    lo<-windows[[i]][1];hi<-windows[[i]][2];z<-d|>filter(time>lo)|>mutate(time=time-lo)
    rr<-z|>group_by(model)|>group_modify(~{
      q<-ipcw_roc_curve(.x$time,.x$event,.x$score,hi-lo)
      if(nrow(q))q$AUC<-weighted_time_auc(.x$time,.x$event,.x$score,hi-lo)
      q
    })|>ungroup()
    title<-paste0(LETTERS[i+2],'. ',if(lo==0)paste0('Within ',hi,' years')else paste0('Years ',lo,'–',hi,'; event-free at ',lo))
    if(!nrow(rr))return(blank_plot(title,'Insufficient cases/controls with valid predictions'))
    rr$window<-paste(lo,hi,sep='-');curves[[i]]<<-rr
    mm<-rr|>group_by(model)|>summarise(AUC=first(AUC),.groups='drop')|>mutate(window=paste(lo,hi,sep='-'));metrics[[i]]<<-mm
    text<-paste(sprintf('%s: %.3f',mm$model,mm$AUC),collapse='\n')
    ggplot(rr,aes(fpr,tpr,color=model))+geom_abline(slope=1,intercept=0,linetype=2,color='grey70')+geom_step(linewidth=.8)+
      annotate('text',x=.98,y=.05,label=text,hjust=1,vjust=0,size=2.7)+coord_equal()+
      scale_color_manual(values=c(`Model 0`='grey45',Biomarkers='#B64C69',Combined='#367F9B'))+
      labs(title=title,x='False positive rate',y='True positive rate',color=NULL)+theme_5c(9)
  })
  save_plot((pa|curve)/wrap_plots(plots,nrow=1)+
    plot_annotation(caption='A: LightGBM importance; B: separate forward-selection path. Held-out predictions; IPCW ROC accounting for censoring. Fixed baseline models, no validation-based feature selection. The >10-year panel is a landmark evaluation.'),
    'c5.Fig14.heldout_ROC.png',20,11.5,outdir=outdir)
  rawdir<-le8_job_dir(outdir,'c5_consolidate')
  write_raw_csv(imp,'c5.mock_importance.csv',rawdir);write_raw_csv(bind_rows(metrics),'c5.mock_auc.csv',rawdir)
  write_raw_csv(bind_rows(curves),'c5.mock_roc.csv',rawdir)
  write_xlsx2(list(importance=imp,selection=seq,auc=bind_rows(metrics),ROC=bind_rows(curves)),le8_artifact_path('c5.Fig14.heldout_ROC.xlsx',outdir))
  list(importance=imp,auc=bind_rows(metrics))
}

le8_mock_brain<-function(counts,atlas_name,title) {
  if(!requireNamespace('sf',quietly=TRUE))return(blank_plot(title,'Install sf to render atlas polygons'))
  env<-new.env();load(file.path(Sys.getenv('LE8_FDIR'),'assets/ggseg',paste0(atlas_name,'.rda')),envir=env)
  a<-get(atlas_name,env);d<-sf::st_as_sf(as.data.frame(a$data))
  d$n_sig<-counts$n_sig[match(d$label,counts$atlas_label)]
  if(atlas_name=='dk'){
    # Reposition the four real atlas views in a compact 2x2 grid (translation only).
    group<-interaction(d$hemi,d$side,drop=TRUE)
    for(g in levels(group)){ii<-which(group==g);bb<-sf::st_bbox(d[ii,]);dx<-if(d$hemi[ii[1]]=='right')400 else 0;dy<-if(d$side[ii[1]]=='lateral')240 else 0
      sf::st_geometry(d)[ii]<-sf::st_geometry(d)[ii]+c(dx-bb['xmin'],dy-bb['ymin'])}
  }else d<-d[d$side=='coronal',]
  # These bundled MULTIPOLYGON atlases are planar drawings. Render their native
  # vertices without coord_sf's geographic transformations or PROJ database.
  # Keep separate polygons and rings (including holes) within each atlas region.
  vertices<-as.data.frame(sf::st_coordinates(d))
  vertices$n_sig<-d$n_sig[vertices$L3]
  vertices$polygon<-interaction(vertices$L3,vertices$L2,drop=TRUE)
  vertices$ring<-interaction(vertices$L3,vertices$L2,vertices$L1,drop=TRUE)
  ggplot(vertices,aes(X,Y,group=polygon,subgroup=ring,fill=n_sig))+
    geom_polygon(rule='evenodd',color='grey45',linewidth=.12)+coord_equal()+scale_fill_gradient(low='#FFF3DC',high='#BD3651',na.value='grey88',limits=c(0,max(1,counts$n_sig,na.rm=TRUE)))+
    labs(title=title,fill='Significant\nbiomarkers')+theme_void()+theme(plot.title=element_text(face='bold',size=10),legend.position='right')
}
le8_mock_c4<-function(a,outdir) {
  if(!nrow(a))return(invisible(NULL))
  global<-a|>filter(measure%in%c('img_total_cortical_gm','img_p26517','img_p24486'))|>
    mutate(structure=recode(measure,img_total_cortical_gm='CGV',img_p26517='SGV',img_p24486='WMH'),mark=case_when(bonferroni<.05~'**',FDR_all<.05~'*',TRUE~''))
  top_global<-global|>group_by(feature)|>summarise(pmin=min(p,na.rm=TRUE),.groups='drop')|>arrange(pmin)|>slice_head(n=35)|>pull(feature)
  global<-global|>filter(feature%in%top_global)|>mutate(feature=factor(feature,levels=rev(top_global)))
  pa<-ggplot(global,aes(structure,feature,fill=beta))+geom_tile()+geom_text(aes(label=mark),size=2)+
    scale_fill_gradient2(low='#3678B0',mid='white',high='#CC4656',midpoint=0,na.value='grey88')+
    labs(title='A. Global structure: top 35 features',x=NULL,y=NULL,fill='Standardized\nbeta')+theme_5c(7)+theme(panel.grid=element_blank())
  if(n_distinct(global$feature)>150)pa<-pa+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(y=paste(n_distinct(global$feature),'features; row names in CSV'))
  counts<-a|>group_by(family,label,measure)|>summarise(n_sig=sum(FDR_all<.05,na.rm=TRUE),tested=sum(is.finite(p)),.groups='drop')
  counts<-counts|>mutate(atlas_label=case_when(
    family%in%c('Cortical area','Cortical volume')~sub('^aparc-Desikan_([lr]h)_(area|volume)_','\\1_',label),
    family=='Subcortical volume'~sub('^aseg_rh_volume_','Right-',sub('^aseg_lh_volume_','Left-',label)),TRUE~NA_character_))
  pb<-le8_mock_brain(counts|>filter(family=='Cortical area'),'dk','B. Cortical surface area')
  pc<-le8_mock_brain(counts|>filter(family=='Cortical volume'),'dk','C. Cortical gray matter volume')
  pd<-le8_mock_brain(counts|>filter(family=='Subcortical volume'),'aseg','D. Subcortical gray matter volume')
  wm<-counts|>filter(family=='White matter FA/MD')|>mutate(metric=ifelse(grepl('_FA_',label),'FA','MD'),tract=sub('^.*_(FA|MD)_','',label))
  pe<-ggplot(wm,aes(metric,tract,fill=n_sig))+geom_tile(color='white')+geom_text(aes(label=n_sig),size=2.4)+scale_fill_gradient(low='#FFF3DC',high='#BD3651')+
    labs(title='E. White matter tracts',x=NULL,y=NULL,fill='Significant\nbiomarkers')+theme_5c(8)+theme(panel.grid=element_blank())
  save_plot((pa|wrap_plots(pb,pc,pd,pe,ncol=2))+plot_layout(widths=c(.7,2))+
    plot_annotation(caption='A: * global BH < 0.05; ** Bonferroni < 0.05. B–E: number of biomarkers passing global BH. Gray atlas regions: not mapped/tested. DK/aseg atlas; FA/MD shown by measured tract.'),
    'c4.Fig14.imaging_atlas.png',21,12,outdir=outdir)
  write_raw_csv(counts,'c4.mock_brain_region_counts.csv',le8_job_dir(outdir,'c4_connect'))
  write_xlsx2(list(global_top35=global,region_counts=counts),le8_artifact_path('c4.Fig14.imaging_atlas.xlsx',outdir))
  invisible(counts)
}
