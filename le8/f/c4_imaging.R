# C4 brain structure connections: baseline omics -> imaging measurements.
# Complete-case standardized linear models; no image imputation or causal claim.
c4_img_manifest <- function(indir) {
  p<-file.path(indir,'Rdata/img.dictionary.csv')
  if(!file.exists(p))return(tibble())
  m<-as_tibble(data.table::fread(p))
  m<-m|>mutate(family=case_when(
    grepl('^aparc-Desikan_.*_area_',label)&!grepl('TotalSurface',label)~'Cortical area',
    grepl('^aparc-Desikan_.*_volume_',label)&!grepl('Total',label)~'Cortical volume',
    grepl('^aseg_[lr]h_volume_(Thalamus-Proper|Caudate|Putamen|Pallidum|Hippocampus|Amygdala|Accumbens-area|VentralDC)$',label)~'Subcortical volume',
    grepl('dMRI_ProbtrackX_(FA|MD)_',label)~'White matter FA/MD',
    field=='p24486'~'White matter hyperintensity',
    field=='p26517'~'Global grey matter',TRUE~NA_character_))|>filter(!is.na(family))
  m|>mutate(adjust_TIV=!family%in%c('Cortical area','White matter FA/MD'))
}
c4_img_fit <- function(dat,feature,measure,covars,min_n=100L) {
  columns<-unique(c(feature,measure,covars));z<-as.data.frame(dat[,columns,drop=FALSE])
  ok<-complete.cases(z)&is.finite(z[[feature]])&is.finite(z[[measure]])
  z<-z[ok,,drop=FALSE];n<-nrow(z)
  empty<-tibble(feature=feature,measure=measure,N=n,beta=NA_real_,SE=NA_real_,p=NA_real_,status='insufficient observations or variance')
  if(n<min_n||!is.finite(sd(z[[feature]]))||sd(z[[feature]])==0||sd(z[[measure]])==0)return(empty)
  cv<-covars[vapply(z[covars],function(x)length(unique(x))>1,logical(1))]
  z[[feature]]<-as.numeric(scale(z[[feature]]));z[[measure]]<-as.numeric(scale(z[[measure]]))
  tryCatch({
    fit<-lm(reformulate(c(feature,cv),response=measure),z)
    s<-coef(summary(fit));if(!feature%in%rownames(s))return(empty)
    tibble(feature=feature,measure=measure,N=n,beta=s[feature,1],SE=s[feature,2],p=s[feature,4],status='ok')
  },error=function(e){empty$status<-conditionMessage(e);empty})
}
run_c4_imaging <- function(layer,disease=NULL,outdir=if(layer=='protein')out.prot else out.met) {
  rawdir<-le8_job_dir(outdir,'c4_connect');dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  enabled<-truthy(Sys.getenv('C4_IMG_ENABLED',unset='TRUE'))
  missing_result<-function(reason){
    status<-tibble(status='not run',detail=reason)
    write_raw_csv(status,'c4.imaging_status.csv',rawdir)
    if(enabled)save_plot(blank_plot('Omics and brain structure',reason),'c4.Fig14.imaging_atlas.png',14,9,outdir=outdir)
    list(status=status,associations=tibble(),fields=tibble())
  }
  if(!enabled)return(missing_result('C4_IMG_ENABLED=FALSE'))
  imgfile<-Sys.getenv('C4_IMG_FILE',unset=file.path(indir,'Rdata/img.raw.rds'))
  if(!file.exists(imgfile))return(missing_result('Prepare unimputed img.raw.rds with ukb/f/img.f.R (phe.sh biom step)'))
  if(is.null(disease)){
    f<-file.path(le8_job_dir(outdir,'c1_correlate'),'c1.res.rds')
    if(!file.exists(f))return(missing_result('C1 outcome association table unavailable'))
    c1<-readRDS(f);le8_check_options(c1);disease<-as_tibble(c1$association%||%c1$pwas_incident%||%c1$MWAS)
  }
  max_features<-as.integer(Sys.getenv('C4_IMG_MAX_FEATURES',unset='100'))
  explicit<-trimws(strsplit(Sys.getenv('C4_IMG_FEATURES',unset=''),',',fixed=TRUE)[[1]])
  explicit<-explicit[nzchar(explicit)]
  if(length(explicit))features<-unique(explicit) else {
    stopifnot(all(c('term','p.value')%in%names(disease)))
    features<-disease|>filter(is.finite(p.value),p.value<.05/nrow(disease))|>arrange(p.value)|>pull(term)|>unique()|>head(max_features)
  }
  if(!length(features))return(missing_result('No C1 Bonferroni-significant feature; no data-driven threshold relaxation'))
  covars<-trimws(strsplit(Sys.getenv('C4_IMG_COVARS',unset=paste(if(le8_custom_adjustment())le8_custom_covars else vars.basic,collapse=',')),',',fixed=TRUE)[[1]])
  datecol<-Sys.getenv('C4_IMG_OUTCOME_DATE',unset=le8_y_date())
  files<-c(imgfile,file.path(indir,'Rdata/img.dictionary.csv'),file.path(indir,'Rdata/img.visits.rds'),
    file.path(indir,'Rdata/all.rds'),if(layer=='protein')file.path(indir,'Rdata/prot.rds')else file.path(indir,'Rdata/met.rds'))
  signature<-list(version=3L,features=features,covars=covars,datecol=datecol,white=Sys.getenv('LE8_WHITE_ONLY',unset='TRUE'),
    files=file.info(files)[,c('size','mtime')],paths=files)
  cache<-file.path(rawdir,'c4.imaging.rds')
  if(file.exists(cache)&&!truthy(Sys.getenv('C4_IMG_REPLACE',unset='FALSE'))){
    old<-readRDS(cache);if(identical(old$signature,signature))return(old)
  }
  message('C4 imaging: ',length(features),' ',layer,' features; complete-case models')
  all<-read_all(unique(c('eid',covars,'center','ethnic.c','date_attend',datecol)))|>filter_analysis_cohort()
  required<-unique(c('eid',covars,'center','date_attend',datecol))
  if(length(setdiff(required,names(all))))return(missing_result(paste('Missing covariates/date:',paste(setdiff(required,names(all)),collapse=','))))
  # Keep the outcome-free baseline cohort without requiring complete LE8 scores.
  all<-all[!is.na(all$date_attend)&(is.na(all[[datecol]])|all[[datecol]]>all$date_attend),,drop=FALSE]
  all$eid<-as.character(all$eid);all$center<-factor(all$center)
  if(layer=='protein'){
    bio<-read_prot();omic_input<-'Existing prepared protein matrix (same input as C1)'
  }else{bio<-read_met();omic_input<-'Existing prepared metabolite matrix'}
  bio$eid<-as.character(bio$eid);features<-intersect(features,names(bio))
  if(!length(features))return(missing_result('Selected features absent from omics input'))
  img<-readRDS(imgfile);img$eid<-as.character(img$eid)
  stopifnot(!anyDuplicated(all$eid),!anyDuplicated(bio$eid),!anyDuplicated(img$eid))
  fields<-c4_img_manifest(indir)
  if(all(c('img_p26552','img_p26583')%in%names(img))){
    img$img_total_cortical_gm<-img$img_p26552+img$img_p26583
    fields<-bind_rows(fields,tibble(field='derived',label='Total cortical grey matter volume',variable='img_total_cortical_gm',family='Global grey matter',adjust_TIV=TRUE))
  }
  fields<-fields|>filter(variable%in%names(img))
  if(!nrow(fields))return(missing_result('No supported brain measurements found in dictionary'))
  vfile<-file.path(indir,'Rdata/img.visits.rds')
  if(file.exists(vfile)){
    visits<-readRDS(vfile);visits$eid<-as.character(visits$eid)
    all<-merge(all,visits,by='eid',all.x=TRUE,sort=FALSE)
    if('p54_i2'%in%names(all))all$center<-factor(all$p54_i2)
    if('p53_i2'%in%names(all))all<-all[!is.na(all$p53_i2)&all$p53_i2>=all$date_attend,,drop=FALSE]
  }
  covars<-unique(c(covars,'center'))
  measures<-unique(c(fields$variable,'img_p26521'));measures<-intersect(measures,names(img))
  dat<-merge(merge(all,bio[,c('eid',features),drop=FALSE],by='eid'),img[,c('eid',measures),drop=FALSE],by='eid')
  rm(all,bio,img);invisible(gc())
  if(!'img_p26521'%in%names(dat)&&any(fields$adjust_TIV))return(missing_result('TIV (26521) missing; volume models not fit without required adjustment'))
  rows<-parallel_map(seq_len(nrow(fields)),function(i){
    f<-fields[i,];cv<-unique(c(covars,if(f$adjust_TIV)'img_p26521'))
    bind_rows(lapply(features,function(x)c4_img_fit(dat,x,f$variable,cv)))
  })
  associations<-bind_rows(rows)|>left_join(fields,by=c('measure'='variable'))|>
    mutate(FDR_all=p.adjust(p,'BH'),bonferroni=pmin(1,p*n()))|>group_by(family)|>mutate(FDR_family=p.adjust(p,'BH'))|>ungroup()
  status<-tibble(status='ok',participants=nrow(dat),features=length(features),measurements=nrow(fields),tested=sum(is.finite(associations$p)),
    input=omic_input,covariates=paste(covars,collapse=','),outcome_date=datecol,
    interpretation='Baseline omics associated with later brain structure; no causal or within-person change estimate')
  result<-list(signature=signature,status=status,associations=associations,fields=fields)
  saveRDS(result,cache);write_raw_csv(status,'c4.imaging_status.csv',rawdir)
  write_raw_csv(associations,'c4.imaging_associations.csv',rawdir);write_raw_csv(fields,'c4.imaging_fields.csv',rawdir)
  plot_c4_imaging(associations, features, outdir)
  result
}

plot_c4_imaging <- function(associations, features, outdir) {
  le8_mock_c4(associations,outdir)
  top<-associations|>filter(is.finite(p))|>arrange(p)|>slice_head(n=20)|>
    mutate(pair=stringr::str_wrap(paste(feature,label,sep=' → '),48),pair=factor(pair,levels=rev(unique(pair))))
  if(nrow(top)){
    g<-ggplot(top,aes(beta,pair,color=family))+geom_vline(xintercept=0,color='grey70')+
      geom_errorbarh(aes(xmin=beta-1.96*SE,xmax=beta+1.96*SE),height=.15)+geom_point(aes(shape=FDR_all<.05))+
      labs(title='Omics and brain structure',subtitle='Top 20 associations; complete-case coefficients; exploratory brain context',x='Standardized beta (95% CI)',y=NULL,shape='Global FDR < 0.05',color=NULL)+theme_5c(9)
    save_plot(g,'c4.Fig12.imaging_associations.png',15,13,outdir=outdir)
    summary<-associations|>group_by(feature,family)|>summarise(tested=sum(is.finite(p)),significant=sum(FDR_all<.05,na.rm=TRUE),.groups='drop')
    show<-summary|>group_by(feature)|>summarise(n=sum(significant),.groups='drop')|>arrange(desc(n),feature)|>slice_head(n=30)|>pull(feature)
    summary<-summary|>filter(feature%in%show)
    g<-ggplot(summary,aes(family,feature,fill=significant))+geom_tile(color='white')+scale_fill_viridis_c()+
      labs(title='Brain associations across measurement families',subtitle='Top 30 features by significant-pair count; FDR across all tested pairs',x=NULL,y=NULL,fill='Significant')+theme_5c(8)+theme(axis.text.x=element_text(angle=25,hjust=1))
    save_plot(g,'c4.Fig13.imaging_overview.png',13,max(7,length(features)*.15),outdir=outdir)
  }
}

attach_c4_imaging <- function(obj,layer,outdir,disease=NULL) {
  obj$imaging<-run_c4_imaging(layer,disease,outdir)
  wbfile<-le8_artifact_path('c4.out.xlsx',outdir)
  if(file.exists(wbfile)){
    wb<-openxlsx::loadWorkbook(wbfile)
    for(nm in c('status','associations','fields')){
      sheet<-paste0('imaging_',nm);if(sheet%in%names(wb))openxlsx::removeWorksheet(wb,sheet)
      openxlsx::addWorksheet(wb,sheet);x<-obj$imaging[[nm]]
      if(nrow(x))openxlsx::writeData(wb,sheet,x)
    }
    openxlsx::saveWorkbook(wb,wbfile,overwrite=TRUE)
  }
  obj
}
