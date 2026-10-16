# FINAL C4 multi-route extensions.
# Adds, without replacing native C4:
#   1) cross-fitted LE8-pillar reconstruction by NS / lasso / elastic-net / ridge;
#   2) global omics biological-age clocks;
#   3) brain MRI clock when the existing UKB imaging inputs are available;
#   4) auto-detected phenotype-system clocks (cardiovascular/metabolic/liver/kidney/lung/immune)
#      and optional user manifest clocks;
#   5) Sleep-Chart-style nonlinear pillar -> BAG -> incident disease connections.
# All disease analyses are explicitly associational. No route is selected by test-Y performance.

c4_final_group_folds <- function(groups,k=5L,seed=2026L) {
  groups <- as.character(groups); ug <- unique(groups)
  if(length(ug)<k) k <- max(2L,length(ug))
  set.seed(seed); ug <- sample(ug)
  mp <- setNames(rep(seq_len(k),length.out=length(ug)),ug)
  unname(mp[groups])
}

c4_final_prepare_matrix <- function(train,test,features,max_features=200L,target=NULL) {
  features <- intersect(features,intersect(names(train),names(test)))
  if(!length(features)) return(NULL)
  # Supervised screening is strictly training-only.
  if(!is.null(target) && length(features)>max_features) {
    rr <- vapply(features,function(v) suppressWarnings(cor(as.numeric(train[[v]]),target,use='pairwise.complete.obs')),numeric(1))
    rr[!is.finite(rr)] <- 0
    features <- names(sort(abs(rr),decreasing=TRUE))[seq_len(min(max_features,length(rr)))]
  }
  A <- as.matrix(data.frame(lapply(train[,features,drop=FALSE],function(z)suppressWarnings(as.numeric(z))),check.names=FALSE))
  B <- as.matrix(data.frame(lapply(test[,features,drop=FALSE],function(z)suppressWarnings(as.numeric(z))),check.names=FALSE))
  keep <- vapply(seq_len(ncol(A)),function(j) sum(is.finite(A[,j]))>=max(50,.2*nrow(A)) && is.finite(sd(A[,j],na.rm=TRUE)) && sd(A[,j],na.rm=TRUE)>0,logical(1))
  A<-A[,keep,drop=FALSE];B<-B[,keep,drop=FALSE];features<-features[keep]
  if(!ncol(A)) return(NULL)
  med <- apply(A,2,function(z){v<-median(z[is.finite(z)],na.rm=TRUE);ifelse(is.finite(v),v,0)})
  for(j in seq_along(med)){A[!is.finite(A[,j]),j]<-med[j];B[!is.finite(B[,j]),j]<-med[j]}
  mu<-colMeans(A);ss<-apply(A,2,sd);ss[!is.finite(ss)|ss==0]<-1
  A<-sweep(sweep(A,2,mu,'-'),2,ss,'/');B<-sweep(sweep(B,2,mu,'-'),2,ss,'/')
  list(train=A,test=B,features=features)
}

c4_final_proxy_multiroute <- function(dat,features,components,rawdir,outdir=NULL,seed=2026L) {
  if(!requireNamespace('glmnet',quietly=TRUE)) {
    z<-tibble::tibble(status='unavailable',reason='glmnet missing');write_raw_csv(z,'c4.final_proxy_status.csv',rawdir);return(list(status=z,metrics=tibble::tibble(),predictions=tibble::tibble()))
  }
  components<-intersect(components,names(dat));features<-intersect(features,names(dat))
  if(!length(components)||length(features)<2) return(list(status=tibble::tibble(status='unavailable',reason='No LE8 target/features'),metrics=tibble::tibble(),predictions=tibble::tibble()))
  maxf<-as.integer(Sys.getenv('C4_FINAL_PROXY_MAX_FEATURES','200'));K<-as.integer(Sys.getenv('C4_FINAL_PROXY_FOLDS','5'))
  groups<-if('.group'%in%names(dat))dat$.group else dat$eid; fold<-c4_final_group_folds(groups,K,seed)
  preds<-list();audit<-list();kk<-0L;ka<-0L
  routes<-c(NS=NA_real_,Lasso=1,ElasticNet=.5,Ridge=0)
  for(cmp in components) {
    yall<-suppressWarnings(as.numeric(dat[[cmp]]));if(sum(is.finite(yall))<300||sd(yall,na.rm=TRUE)==0)next
    for(fd in sort(unique(fold))) {
      tr<-which(fold!=fd & is.finite(yall));te<-which(fold==fd & is.finite(yall));if(length(tr)<200||length(te)<30)next
      pp<-c4_final_prepare_matrix(dat[tr,,drop=FALSE],dat[te,,drop=FALSE],features,maxf,yall[tr]);if(is.null(pp))next
      rr<-vapply(seq_len(ncol(pp$train)),function(j)suppressWarnings(cor(pp$train[,j],yall[tr],use='complete.obs')),numeric(1));rr[!is.finite(rr)]<-0
      bestj<-which.max(abs(rr));xtr<-pp$train[,bestj];xte<-pp$test[,bestj]
      lmfit<-lm.fit(cbind(1,xtr),yall[tr]);ns<-drop(cbind(1,xte)%*%lmfit$coefficients)
      kk<-kk+1L;preds[[kk]]<-tibble::tibble(eid=as.character(dat$eid[te]),component=cmp,fold=fd,route='NS',observed=yall[te],predicted=ns,selected_features=1L,lead_feature=pp$features[bestj])
      # Group-safe inner fold IDs for cv.glmnet.
      inner_groups<-as.character(groups[tr]);ifold<-c4_final_group_folds(inner_groups,min(5L,length(unique(inner_groups))),seed+fd+match(cmp,components)*100)
      for(rn in c('Lasso','ElasticNet','Ridge')) {
        alpha<-routes[[rn]]
        fit<-tryCatch(glmnet::cv.glmnet(pp$train,yall[tr],family='gaussian',alpha=alpha,foldid=ifold,standardize=FALSE,keep=FALSE),error=function(e)NULL)
        if(is.null(fit))next
        pr<-as.numeric(predict(fit,pp$test,s='lambda.1se'))
        bb<-as.matrix(coef(fit,s='lambda.1se'));nsel<-sum(abs(bb[-1,1])>1e-10)
        kk<-kk+1L;preds[[kk]]<-tibble::tibble(eid=as.character(dat$eid[te]),component=cmp,fold=fd,route=rn,observed=yall[te],predicted=pr,selected_features=nsel,lead_feature=NA_character_)
      }
      ka<-ka+1L;audit[[ka]]<-tibble::tibble(component=cmp,fold=fd,training_N=length(tr),test_N=length(te),screened_features=ncol(pp$train),NS_feature=pp$features[bestj],NS_training_r=rr[bestj])
    }
  }
  pr<-dplyr::bind_rows(preds);au<-dplyr::bind_rows(audit)
  if(!nrow(pr)){z<-tibble::tibble(status='unavailable',reason='No cross-fitted proxy predictions');write_raw_csv(z,'c4.final_proxy_status.csv',rawdir);return(list(status=z,metrics=tibble::tibble(),predictions=pr))}
  met<-pr|>dplyr::group_by(.data$component,.data$route)|>dplyr::summarise(N=dplyr::n(),R=cor(.data$observed,.data$predicted,use='complete.obs'),R2=1-sum((.data$observed-.data$predicted)^2)/sum((.data$observed-mean(.data$observed))^2),RMSE=sqrt(mean((.data$observed-.data$predicted)^2)),median_selected=median(.data$selected_features,na.rm=TRUE),.groups='drop')
  ns<-met|>dplyr::filter(.data$route=='NS')|>dplyr::select(.data$component,NS_R2=.data$R2,NS_R=.data$R)
  met<-met|>dplyr::left_join(ns,by='component')|>dplyr::mutate(delta_R2_vs_NS=.data$R2-.data$NS_R2,delta_absR_vs_NS=abs(.data$R)-abs(.data$NS_R))
  winners<-met|>dplyr::filter(.data$route!='NS')|>dplyr::group_by(.data$component)|>dplyr::slice_max(.data$R2,n=1,with_ties=FALSE)|>dplyr::ungroup()|>dplyr::mutate(YS_beats_NS=.data$delta_R2_vs_NS>0)
  write_raw_csv(pr,'c4.final_proxy_oof_predictions.csv',rawdir);write_raw_csv(met,'c4.final_proxy_route_metrics.csv',rawdir);write_raw_csv(winners,'c4.final_proxy_route_winners.csv',rawdir);write_raw_csv(au,'c4.final_proxy_fold_audit.csv',rawdir)
  st<-tibble::tibble(status='completed',components=dplyr::n_distinct(pr$component),routes=paste(unique(pr$route),collapse=';'),interpretation='Cross-fitted route exploration; route winner is not selected using disease Y')
  write_raw_csv(st,'c4.final_proxy_status.csv',rawdir)
  if(!is.null(outdir)&&requireNamespace('ggplot2',quietly=TRUE)){
    p<-ggplot2::ggplot(met,ggplot2::aes(.data$route,.data$R2,fill=.data$route))+ggplot2::geom_col()+ggplot2::facet_wrap(~component,scales='free_y')+ggplot2::geom_hline(data=ns,ggplot2::aes(yintercept=.data$NS_R2),linetype=2,color='grey35')+ggplot2::labs(title='C4 multi-route LE8 proxy reconstruction',subtitle='All values are outer group-cross-fitted; dashed line is the single-biomarker NS route',x=NULL,y='OOF R²',fill=NULL)+theme_5c(9)+ggplot2::theme(legend.position='bottom')
    save_plot(p,'c4.Fig18.final_proxy_routes.png',16,10,outdir=outdir)
  }
  list(status=st,metrics=met,winners=winners,predictions=pr,audit=au)
}

c4_final_clock_oof <- function(d,features,clock,source='omics',groups=NULL,seed=2026L,max_features=500L) {
  if(!requireNamespace('glmnet',quietly=TRUE)||!'age'%in%names(d))return(tibble::tibble())
  features<-intersect(features,names(d));features<-setdiff(features,c('age','eid','.group'));if(length(features)<3)return(tibble::tibble())
  age<-suppressWarnings(as.numeric(d$age));ok<-is.finite(age);if(sum(ok)<500)return(tibble::tibble())
  if(is.null(groups))groups<-if('.group'%in%names(d))d$.group else d$eid
  K<-min(5L,length(unique(groups[ok])));if(K<3)return(tibble::tibble())
  fd<-c4_final_group_folds(groups,K,seed);rows<-list();k<-0L
  for(f in sort(unique(fd))) {
    tr<-which(fd!=f&ok);te<-which(fd==f&ok);if(length(tr)<300||length(te)<50)next
    pp<-c4_final_prepare_matrix(d[tr,,drop=FALSE],d[te,,drop=FALSE],features,max_features,age[tr]);if(is.null(pp)||ncol(pp$train)<3)next
    ifold<-c4_final_group_folds(groups[tr],min(5L,length(unique(groups[tr]))),seed+f)
    for(route in c('ridge','elasticnet','lasso')){
      alpha<-if(route=='ridge')0 else if(route=='lasso')1 else .5
      fit<-tryCatch(glmnet::cv.glmnet(pp$train,age[tr],family='gaussian',alpha=alpha,foldid=ifold,standardize=FALSE),error=function(e)NULL);if(is.null(fit))next
      pred<-as.numeric(predict(fit,pp$test,s='lambda.1se'));k<-k+1L
      rows[[k]]<-tibble::tibble(eid=as.character(d$eid[te]),age=age[te],clock=clock,source=source,route=route,predicted_age=pred,fold=f,feature_N=ncol(pp$train))
    }
  }
  z<-dplyr::bind_rows(rows);if(!nrow(z))return(z)
  # Age-bias correction uses OOF predictions only; this correction never sees disease Y.
  z<-z|>dplyr::group_by(.data$clock,.data$route)|>dplyr::group_modify(function(g,key){g$raw_gap<-g$predicted_age-g$age;fit<-lm(raw_gap~age,data=g);g$BAG<-residuals(fit);g})|>dplyr::ungroup()
  z
}

c4_final_auto_system_features <- function(nms) {
  nms0<-tolower(nms);pick<-function(rx)nms[grepl(rx,nms0,perl=TRUE)]
  list(
    cardiovascular=pick('(^|[._])(sbp|dbp|pulse|resting.?heart|arterial.?stiff)'),
    metabolic=pick('(^|[._])(bmi|waist|glucose|hba1c|cholesterol|ldl|hdl|trig|apob|apoa)'),
    liver=pick('(^|[._])(alt|ast|ggt|bilirubin|albumin|alp)([._]|$)'),
    kidney=pick('creatinine|cystatin|egfr|(^|[._])urea([._]|$)'),
    lung=pick('fev1|fvc|pef|spirom'),
    immune=pick('(^|[._])(crp|wbc|white.?blood|neutroph|lymph|monocyte|platelet)')
  )
}

c4_final_clock_connections <- function(clock_rows,phen,time_col,event_col,components,covars,rawdir,outdir=NULL) {
  if(!nrow(clock_rows)||!requireNamespace('mgcv',quietly=TRUE))return(list())
  # Retain a single age column and adjust the disease models for age as requested.
  phen<-phen[,setdiff(names(phen),setdiff(names(clock_rows),'eid')),drop=FALSE]
  d<-dplyr::left_join(clock_rows,phen,by='eid');components<-intersect(components,names(d));covars<-setdiff(intersect(covars,names(d)),'sex')
  s<-mgcv::s
  curves<-list();tests<-list();disease<-list();kc<-kt<-kd<-0L
  for(key in split(seq_len(nrow(d)),interaction(d$clock,d$route,drop=TRUE))) {
    q<-d[key,,drop=FALSE];clock<-q$clock[1];route<-q$route[1]
    if(sum(is.finite(q$BAG))<300)next
    for(cmp in components){
      x<-suppressWarnings(as.numeric(q[[cmp]]));use<-is.finite(x)&is.finite(q$BAG);if(sum(use)<250||length(unique(x[use]))<5)next
      a<-q[use,,drop=FALSE];a$.x<-x[use]
      rhs<-c("s(.x,k=5,bs='cr')",if('sex'%in%names(a))'factor(sex)',setdiff(covars,cmp))
      gm<-tryCatch(mgcv::gam(as.formula(paste('BAG ~',paste(rhs,collapse=' + '))),data=a,method='REML'),error=function(e)NULL);if(is.null(gm))next
      st<-summary(gm)$s.table;pnon<-if(nrow(st))st[1,'p-value']else NA_real_;edf<-if(nrow(st))st[1,'edf']else NA_real_
      grid<-seq(quantile(a$.x,.02,na.rm=TRUE),quantile(a$.x,.98,na.rm=TRUE),length.out=101);nd<-data.frame(.x=grid)
      for(cv in c(if('sex'%in%names(a))'sex',covars)){
        if(is.factor(a[[cv]]))nd[[cv]]<-factor(levels(a[[cv]])[1],levels=levels(a[[cv]])) else if(cv=='sex')nd[[cv]]<-a[[cv]][which(!is.na(a[[cv]]))[1]] else nd[[cv]]<-median(as.numeric(a[[cv]]),na.rm=TRUE)
      }
      pr<-tryCatch(predict(gm,newdata=nd,se.fit=TRUE),error=function(e)NULL);if(is.null(pr))next;j<-which.min(pr$fit)
      kc<-kc+1L;curves[[kc]]<-tibble::tibble(clock=clock,route=route,component=cmp,x=grid,BAG=as.numeric(pr$fit),lo=as.numeric(pr$fit-1.96*pr$se.fit),hi=as.numeric(pr$fit+1.96*pr$se.fit))
      # Sex interaction is an explicit sensitivity route.
      pint<-NA_real_
      if('sex'%in%names(a)&&length(unique(a$sex))>=2){
        a$.sex<-factor(a$sex);gmi<-tryCatch(mgcv::gam(as.formula(paste('BAG ~ .sex + s(.x,by=.sex,k=5,bs="cr")',if(length(covars))paste0(' + ',paste(covars,collapse=' + '))else'')),data=a,method='REML'),error=function(e)NULL)
        if(!is.null(gmi))pint<-tryCatch(anova(gm,gmi,test='Chisq')$`Pr(>Chi)`[2],error=function(e)NA_real_)
      }
      kt<-kt+1L;tests[[kt]]<-tibble::tibble(clock=clock,route=route,component=cmp,N=nrow(a),nonlinear_p=pnon,edf=edf,optimum=grid[j],BAG_min=pr$fit[j],sex_interaction_p=pint)
      if(all(c(time_col,event_col)%in%names(a))){
        cc<-complete.cases(a[,unique(c(time_col,event_col,'BAG','.x',covars,if('sex'%in%names(a))'sex')),drop=FALSE]);b<-a[cc,,drop=FALSE]
        if(nrow(b)>=300&&sum(b[[event_col]]==1)>=30){
          cv<-c(if('sex'%in%names(b))'factor(sex)',setdiff(covars,cmp));tail<-if(length(cv))paste0(' + ',paste(cv,collapse=' + '))else''
          fbag<-as.formula(paste0('survival::Surv(',time_col,',',event_col,') ~ scale(BAG)',tail))
          f0<-as.formula(paste0('survival::Surv(',time_col,',',event_col,') ~ scale(.x)',tail));f1<-as.formula(paste0('survival::Surv(',time_col,',',event_col,') ~ scale(.x) + scale(BAG)',tail))
          mb<-tryCatch(survival::coxph(fbag,b),error=function(e)NULL);m0<-tryCatch(survival::coxph(f0,b),error=function(e)NULL);m1<-tryCatch(survival::coxph(f1,b),error=function(e)NULL)
          if(!is.null(mb)&&!is.null(m0)&&!is.null(m1)){
            s0<-coef(summary(m0));s1<-coef(summary(m1));sb<-coef(summary(mb));b0<-s0['scale(.x)','coef'];b1<-s1['scale(.x)','coef']
            kd<-kd+1L;disease[[kd]]<-tibble::tibble(clock=clock,route=route,component=cmp,N=nrow(b),events=sum(b[[event_col]]==1),BAG_HR=exp(sb['scale(BAG)','coef']),BAG_p=sb['scale(BAG)','Pr(>|z|)'],pillar_beta_without_BAG=b0,pillar_beta_with_BAG=b1,attenuation_fraction=ifelse(abs(b0)>1e-12,(b0-b1)/b0,NA_real_),interpretation='Associational attenuation; not causal mediation')
          }
        }
      }
    }
  }
  cu<-dplyr::bind_rows(curves);te<-dplyr::bind_rows(tests);di<-dplyr::bind_rows(disease)
  if(nrow(te)){te$FDR_nonlinear<-p.adjust(te$nonlinear_p,'BH');te$FDR_sex_interaction<-p.adjust(te$sex_interaction_p,'BH')}
  if(nrow(di))di$FDR_BAG_disease<-p.adjust(di$BAG_p,'BH')
  write_raw_csv(te,'c4.clock_pillar_tests.csv',rawdir);write_raw_csv(cu,'c4.clock_pillar_curves.csv',rawdir);write_raw_csv(di,'c4.clock_disease_connections.csv',rawdir)
  if(!is.null(outdir)&&nrow(te)&&requireNamespace('ggplot2',quietly=TRUE)){
    p<-ggplot2::ggplot(te,ggplot2::aes(.data$component,.data$clock,fill=pmin(-log10(pmax(.data$FDR_nonlinear,1e-300)),10)))+ggplot2::geom_tile(color='white')+ggplot2::facet_wrap(~route)+ggplot2::labs(title='C4 biological-age clock connections',subtitle='Nonlinear LE8 pillar → biological-age-gap associations; FDR across tested clock–pillar pairs',x=NULL,y=NULL,fill='-log10(FDR)')+theme_5c(8.5)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))
    save_plot(p,'c4.Fig20.clock_connection_heatmap.png',16,10,outdir=outdir)
  }
  list(tests=te,curves=cu,disease=di)
}

run_c4_final_extensions <- function(layer=c('protein','metabolite')) {
  layer<-match.arg(layer);outdir<-if(layer=='protein')out.prot else out.met;rawdir<-le8_job_dir(outdir,'c4_connect');dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  biom<-if(layer=='protein')read_prot() else read_met();biom$eid<-as.character(biom$eid);features<-setdiff(names(biom),'eid')
  all0<-read_all();need<-unique(c('eid','age','sex','ethnic.c',Sys.getenv('LE8_GROUP_COLUMN',Sys.getenv('PGS_GROUP_COLUMN','')),le8_custom_covars,vars.basic,vars.le8,'birth_date','date_attend','date_lost','date_death',paste0('fod_icd10_',Y)))
  ph<-all0[,intersect(need,names(all0)),drop=FALSE]|>filter_analysis_cohort()|>make_outcome(Y);ph$eid<-as.character(ph$eid);ph$.group<-le8_participant_groups(ph);rm(all0);invisible(gc())
  dat<-dplyr::inner_join(ph,biom,by='eid');rm(biom);invisible(gc())
  tvar<-paste0(Y,'.t2e');evar<-paste0(Y,'.Yt2e');comps<-intersect(vars.le8,names(dat));covs<-intersect(if(le8_custom_adjustment())le8_custom_covars else vars.basic,names(dat))
  proxy<-c4_final_proxy_multiroute(dat,features,comps,rawdir,outdir,SEED+4100L)

  clocks<-list();kc<-0L
  # Global omics clocks require no extra external data.
  om<-c4_final_clock_oof(dat,features,paste0(if(layer=='protein')'proteomic' else 'metabolomic','_global'),source=layer,groups=dat$.group,seed=SEED+4200L,max_features=as.integer(Sys.getenv('C4_CLOCK_MAX_OMICS','500')))
  if(nrow(om)){kc<-kc+1L;clocks[[kc]]<-om}

  # Auto phenotype-system clocks from existing all.rds variables. These are only
  # created when >=3 prespecified-pattern variables are actually present.
  base<-ph[,intersect(c('eid','age','.group'),names(ph)),drop=FALSE]
  if(toupper(Sys.getenv('C4_FINAL_SYSTEM_CLOCKS','TRUE'))%in%c('TRUE','1','YES')){
    aa<-read_all();aa$eid<-as.character(aa$eid);sys<-c4_final_auto_system_features(names(aa))
    for(nm in names(sys)){
      vv<-unique(setdiff(sys[[nm]],c('eid','age',vars.le8)));if(length(vv)<3)next
      ds<-dplyr::inner_join(base,aa[,c('eid',vv),drop=FALSE],by='eid')
      cr<-c4_final_clock_oof(ds,vv,paste0('phenotype_',nm),source='all.rds',groups=ds$.group,seed=SEED+4300L+match(nm,names(sys)),max_features=80L)
      if(nrow(cr)){kc<-kc+1L;clocks[[kc]]<-cr}
    }
    rm(aa);invisible(gc())
  }

  # Brain MRI age clock uses the existing C4 imaging files if already prepared.
  imgfile<-Sys.getenv('C4_IMG_FILE',unset=file.path(indir,'Rdata/img.raw.rds'))
  if(file.exists(imgfile)){
    img<-tryCatch(readRDS(imgfile),error=function(e)NULL)
    if(!is.null(img)&&'eid'%in%names(img)){
      img$eid<-as.character(img$eid);ff<-tryCatch(c4_img_manifest(indir),error=function(e)tibble::tibble());vv<-intersect(unique(ff$variable),names(img))
      if(length(vv)>=3){ds<-dplyr::inner_join(ph[,intersect(c('eid','age','.group'),names(ph)),drop=FALSE],img[,c('eid',vv),drop=FALSE],by='eid');cr<-c4_final_clock_oof(ds,vv,'MRI_brain',source='img.raw.rds',groups=ds$.group,seed=SEED+4400L,max_features=300L);if(nrow(cr)){kc<-kc+1L;clocks[[kc]]<-cr}}
      rm(img);invisible(gc())
    }
  }

  # Optional explicit organ/system manifest: clock,source,variable. This lets a
  # future UKB phenotype/imaging set be added without changing C4 code.
  mf<-Sys.getenv('C4_CLOCK_MANIFEST','')
  if(nzchar(mf)&&file.exists(mf)){
    mm<-tryCatch(as.data.frame(data.table::fread(mf)),error=function(e)NULL)
    if(!is.null(mm)&&all(c('clock','source','variable')%in%names(mm))){
      for(nm in unique(as.character(mm$clock))){
        m<-mm[as.character(mm$clock)==nm,,drop=FALSE];src<-unique(as.character(m$source));vv<-unique(as.character(m$variable))
        if(length(src)!=1||!src%in%c('all','img','omics')||length(vv)<3)next
        if(src=='omics'){
          ds<-dat;vv<-intersect(vv,names(ds))
        } else if(src=='all'){
          aa2<-read_all();aa2$eid<-as.character(aa2$eid);vv<-intersect(vv,names(aa2));ds<-dplyr::inner_join(base,aa2[,c('eid',vv),drop=FALSE],by='eid');rm(aa2);invisible(gc())
        } else {
          if(!file.exists(imgfile))next;im2<-tryCatch(readRDS(imgfile),error=function(e)NULL);if(is.null(im2)||!'eid'%in%names(im2))next;im2$eid<-as.character(im2$eid);vv<-intersect(vv,names(im2));ds<-dplyr::inner_join(base,im2[,c('eid',vv),drop=FALSE],by='eid');rm(im2);invisible(gc())
        }
        if(length(vv)>=3){cr<-c4_final_clock_oof(ds,vv,paste0('manifest_',nm),source=paste0('manifest:',src),groups=ds$.group,seed=SEED+4500L+match(nm,unique(as.character(mm$clock))),max_features=500L);if(nrow(cr)){kc<-kc+1L;clocks[[kc]]<-cr}}
      }
    }
  }

  clockrows<-dplyr::bind_rows(clocks)
  if(nrow(clockrows)){
    cm<-clockrows|>dplyr::group_by(.data$clock,.data$route,.data$source)|>dplyr::summarise(N=dplyr::n(),R=cor(.data$predicted_age,.data$age,use='complete.obs'),R2=1-sum((.data$predicted_age-.data$age)^2)/sum((.data$age-mean(.data$age))^2),MAE=mean(abs(.data$predicted_age-.data$age)),BAG_SD=sd(.data$BAG),median_features=median(.data$feature_N),.groups='drop')
    write_raw_csv(clockrows,'c4.clock_oof_predictions.csv',rawdir);write_raw_csv(cm,'c4.clock_performance.csv',rawdir)
    conn<-c4_final_clock_connections(clockrows,ph,tvar,evar,comps,covs,rawdir,outdir)
    if(requireNamespace('ggplot2',quietly=TRUE)){
      p<-ggplot2::ggplot(cm,ggplot2::aes(.data$R2,reorder(paste(.data$clock,.data$route,sep=' | '),.data$R2),fill=.data$source))+ggplot2::geom_col()+ggplot2::labs(title='C4 cross-fitted biological-age clocks',x='OOF age-prediction R²',y=NULL,fill=NULL,caption='Global omics uses existing prot/met; phenotype/brain clocks run only when corresponding UKB variables already exist.')+theme_5c(9)+ggplot2::theme(legend.position='bottom');save_plot(p,'c4.Fig19.biological_age_clocks.png',14,9,outdir=outdir)
    }
  } else {cm<-tibble::tibble();conn<-list();write_raw_csv(tibble::tibble(status='unavailable',reason='No age clock had sufficient existing variables'),'c4.clock_status.csv',rawdir)}
  status<-tibble::tibble(status='completed',proxy_routes=nrow(proxy$metrics%||%tibble::tibble()),clocks=nrow(cm),note='FINAL extensions are additive; native C4 outputs remain primary')
  write_raw_csv(status,'c4.final_extensions_status.csv',rawdir)
  list(proxy=proxy,clock_performance=cm,clock_connections=conn,status=status)
}
