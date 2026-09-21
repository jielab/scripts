# Common-cohort prediction for any Y. No C1/C2/C3 selection enters validation.
c5_csv<-function(name,default)Filter(nzchar,trimws(strsplit(Sys.getenv(name,default),",",fixed=TRUE)[[1]]))
c5_assert_ids<-function(x){
  if(!"eid"%in%names(x)||anyNA(x$eid)||anyDuplicated(x$eid))stop("Missing/duplicate participant IDs")
  x$eid<-as.character(x$eid);x
}
c5_group_folds<-function(group,k,seed){
  if(anyNA(group)||any(!nzchar(as.character(group))))stop("Missing split group")
  g<-sort(unique(as.character(group)));if(length(g)<k)stop("Too few independent groups")
  set.seed(seed);f<-sample(rep(seq_len(k),length.out=length(g)));f[match(group,g)]
}
c5_landmark<-function(d,L,end){
  covered<-if(".entry"%in%names(d))is.finite(d$.entry)&d$.entry<=L else rep(TRUE,nrow(d))
  d<-d[covered&is.finite(d$.time)&d$.time>L&d$.event%in%0:1,,drop=FALSE]
  # Preserve held-out follow-up beyond the evaluation horizon: dynamic AUC
  # controls require time > horizon, not administratively truncated equality.
  d$.time<-d$.time-L;d
}
c5_prs_manifest_input<-function(Y){
  f<-Sys.getenv("C5_PRS_MANIFEST","")
  if(!nzchar(f))stop("C5 requires C5_PRS_MANIFEST: explicit disease PRS for each Y; omic PGS is not a disease PRS")
  m<-read.csv(f,stringsAsFactors=FALSE,check.names=FALSE)
  required<-c("Y","file","column","source","build","ancestry","sample_overlap")
  if(!all(required%in%names(m)))stop("PRS manifest requires: ",paste(required,collapse=","))
  m<-m[m$Y==Y,,drop=FALSE]
  if(nrow(m)!=1||anyNA(m[,required])||any(!nzchar(as.matrix(m[,required]))))stop("Exactly one complete PRS manifest row is required for ",Y)
  if(m$sample_overlap!="none")stop("Use an external disease PRS whose discovery excludes this cohort. Internally derived PRS requires end-to-end outer-fold GWAS/scoring and is not supported by a single precomputed score column.")
  if(!grepl("^(/|[A-Za-z]:)",m$file))m$file<-file.path(dirname(f),m$file)
  d<-if(grepl("[.]rds$",m$file,ignore.case=TRUE))readRDS(m$file)else data.table::fread(m$file,data.table=FALSE)
  d<-c5_assert_ids(as.data.frame(d))
  if(!m$column%in%names(d)||!is.numeric(d[[m$column]]))stop("PRS column missing or nonnumeric")
  d<-data.frame(eid=d$eid,disease_PRS=d[[m$column]])
  list(data=d[is.finite(d$disease_PRS),,drop=FALSE],manifest=m)
}
# Exact endpoint naming, never a fuzzy match to biomarker PGS or another trait.
# An explicit manifest takes precedence. Unknown discovery overlap is labelled,
# not silently treated as independent; no participant's missing PRS is imputed.
c5_prs_input<-function(Y,ph){
  column<-paste0(Y,".pgs")
  unavailable<-function(status,m=NULL){
    if(is.null(m))m<-data.frame(Y=Y,file="all.rds",column=column,source="not supplied",
      build="unknown",ancestry="unknown",sample_overlap="unknown")
    m$status<-status;m$enabled<-FALSE
    list(data=data.frame(eid=as.character(ph$eid)),manifest=m,enabled=FALSE)
  }
  if(nzchar(Sys.getenv("C5_PRS_MANIFEST",""))){
    obj<-c5_prs_manifest_input(Y)
    obj$manifest$status<-"explicit manifest; discovery independence declared by user"
  }else{
    if(!column%in%names(ph))return(unavailable(paste("absent:",column,"; non-PRS comparisons continue")))
    if(!is.numeric(ph[[column]]))return(unavailable(paste("unusable nonnumeric column:",column)))
    overlap<-Sys.getenv("C5_DISEASE_PRS_OVERLAP","unknown")
    if(!overlap%in%c("none","unknown","yes"))stop("C5_DISEASE_PRS_OVERLAP must be none, unknown or yes")
    m<-data.frame(Y=Y,file="all.rds",column=column,
      source=Sys.getenv("C5_DISEASE_PRS_SOURCE","all.rds exact endpoint column; provenance not supplied"),
      build=Sys.getenv("C5_DISEASE_PRS_BUILD","unknown"),ancestry=Sys.getenv("C5_DISEASE_PRS_ANCESTRY","unknown"),
      sample_overlap=overlap)
    if(overlap=="yes")return(unavailable("known discovery overlap; precomputed PRS excluded from validation",m))
    g<-ph[[column]];obj<-list(data=data.frame(eid=as.character(ph$eid[is.finite(g)]),disease_PRS=g[is.finite(g)]),manifest=m)
    obj$manifest$status<-if(overlap=="none")"exact all.rds endpoint PRS; independence declared by user"else
      "exploratory endpoint PRS; discovery overlap unknown; outer CV does not resolve PRS discovery leakage"
  }
  z<-obj$data$disease_PRS
  if(length(z)<2L||!is.finite(sd(z))||sd(z)<=0)return(unavailable("PRS has insufficient finite values or zero variance",obj$manifest))
  obj$enabled<-TRUE;obj$manifest$enabled<-TRUE
  obj$manifest$N_available<-length(z);obj
}
c5_join_omics<-function(ph,prot,met,prs){
  ph<-c5_assert_ids(ph);map<-list()
  for(layer in c("prot","met")){
    b<-c5_assert_ids(if(layer=="prot")prot else met)
    ff<-setdiff(names(b),"eid")
    if(any(grepl("([.]pgs$|[.]prs$|^fod_|^date_|[.]t2e$|[.]Yt2e$)",ff,ignore.case=TRUE)))stop("Non-assay columns in ",layer," matrix")
    if(!all(vapply(b[ff],is.numeric,logical(1))))stop("Omics matrices require numeric assay columns only")
    map[[layer]]<-data.frame(feature=paste0(layer,"__",ff),assay=ff,layer=layer)
    names(b)[match(ff,names(b))]<-map[[layer]]$feature
    ph<-merge(ph,b,by="eid",sort=FALSE)
  }
  ph<-merge(ph,c5_assert_ids(prs),by="eid",sort=FALSE)
  list(data=ph[order(ph$eid),,drop=FALSE],map=bind_rows(map))
}
c5_attach_omic_pgs<-function(dat,map){
  mf<-Sys.getenv("C5_OMIC_PGS_MANIFEST","")
  if(!nzchar(mf))return(list(data=dat,map=map,inputs=character(),status="not supplied; disease PRS and biomarker PGS are distinct"))
  m<-read.csv(mf,stringsAsFactors=FALSE)
  req<-c("layer","file","map_file","source","sample_overlap")
  if(!all(req%in%names(m))||anyNA(m[,req])||anyDuplicated(m$layer)||!setequal(m$layer,c("prot","met"))||any(m$sample_overlap!="none"))
    stop("Omic-PGS manifest requires prot/met rows and independent discovery (sample_overlap=none)")
  inputs<-mf;map$pgs_feature<-NA_character_
  for(i in seq_len(nrow(m))){
    path<-function(p)if(grepl("^(/|[A-Za-z]:)",p))p else file.path(dirname(mf),p)
    f<-path(m$file[i]);mp<-path(m$map_file[i]);inputs<-c(inputs,f,mp)
    lookup<-read.csv(mp,stringsAsFactors=FALSE)
    if(!all(c("assay","column")%in%names(lookup))||anyDuplicated(lookup$assay)||anyDuplicated(lookup$column))stop("Omic-PGS map needs one unique column per assay")
    g<-c5_assert_ids(as.data.frame(if(grepl("[.]rds$",f,ignore.case=TRUE))readRDS(f)else data.table::fread(f,data.table=FALSE)))
    hit<-which(map$layer==m$layer[i]&map$assay%in%lookup$assay)
    cols<-lookup$column[match(map$assay[hit],lookup$assay)]
    if(!all(cols%in%names(g))||!all(vapply(g[cols],is.numeric,logical(1))))stop("PGS mapped columns missing or nonnumeric")
    names_new<-paste0("pgs__",map$feature[hit]);map$pgs_feature[hit]<-names_new
    g<-g[,c("eid",cols),drop=FALSE];names(g)<-c("eid",names_new)
    dat<-merge(dat,g,by="eid",sort=FALSE)
  }
  list(data=dat[order(dat$eid),,drop=FALSE],map=map,inputs=inputs,status="independent matched biomarker PGS available")
}
# Outcome-free, covariate-adjusted domain ranking. This is association strength,
# not a mediation estimate. A feature may be relevant to both domains.
c5_domain_rank<-function(train,features,target,basic){
  d<-as.data.frame(train);ok<-is.finite(d[[target]])
  if(sum(ok)<100)stop("Insufficient observed training target: ",target)
  d<-d[ok,,drop=FALSE];z<-le8_prepare_prediction_matrix(d,d,basic)
  q<-qr(cbind(1,z$train));y<-qr.resid(q,d[[target]]);sy<-sqrt(sum(y*y))
  if(!is.finite(sy)||sy<=0)stop("Invariant domain target")
  blocks<-split(features,ceiling(seq_along(features)/64))
  out<-map_dfr(blocks,function(bb){
    X<-as.matrix(d[,bb,drop=FALSE]);storage.mode(X)<-"double"
    observed<-colSums(is.finite(X))
    for(j in seq_len(ncol(X))){v<-X[,j];med<-median(v[is.finite(v)],na.rm=TRUE)
      if(!is.finite(med))med<-0;v[!is.finite(v)]<-med;X[,j]<-v}
    X<-qr.resid(q,X);den<-sqrt(colSums(X*X))*sy
    r<-as.numeric(crossprod(X,y))/den;r[observed<100|!is.finite(r)]<-NA_real_
    tibble(feature=bb,domain=target,r=r,N_observed=observed)
  })
  out|>filter(is.finite(r))|>arrange(desc(abs(r)),feature)
}
c5_balance<-function(queues,k){
  ans<-character()
  while(length(ans)<k&&any(lengths(queues)>0))for(j in seq_along(queues)){
    queues[[j]]<-setdiff(queues[[j]],ans)
    if(length(queues[[j]])&&length(ans)<k){ans<-c(ans,queues[[j]][1]);queues[[j]]<-queues[[j]][-1]}
  }
  ans
}
c5_ys_panel<-function(queues,k,ns=NULL){
  # Match the NS allocation as well as total assay budget.
  panel<-character()
  for(layer in c("prot","met")){
    kl<-if(layer=="prot")ceiling(k/2)else floor(k/2)
    q<-lapply(queues,function(v)v[startsWith(v,paste0(layer,"__"))])
    required<-if(is.null(ns))kl else max(1,floor(.8*kl))
    keep<-c5_balance(q,required)
    if(length(keep)<required)return(character())
    if(!is.null(ns))keep<-head(c(keep,setdiff(ns[startsWith(ns,paste0(layer,"__"))],keep)),kl)
    panel<-c(panel,keep)
  }
  panel
}
c5_gaussian_proxy<-function(train,test,target,features,seed){
  ok<-is.finite(train[[target]]);tr<-train[ok,,drop=FALSE]
  if(sum(ok)<100)stop("Too few observed proxy targets")
  X<-le8_prepare_prediction_matrix(tr,test,features)
  if(!ncol(X$train))stop("No usable proxy assays")
  if(ncol(X$train)==1){
    fit<-lm.fit(cbind(1,X$train),tr[[target]])
    pred<-as.numeric(cbind(1,X$test)%*%fit$coefficients);b<-fit$coefficients
  }else{
    fd<-c5_group_folds(tr$.group,5,seed)
    fit<-glmnet::cv.glmnet(X$train,tr[[target]],family="gaussian",alpha=0,
      foldid=fd,standardize=TRUE,type.measure="mse")
    pred<-as.numeric(predict(fit,X$test,s="lambda.1se"));b<-as.numeric(coef(fit,s="lambda.1se"))
  }
  if(any(!is.finite(pred)))stop("Non-finite proxy prediction; no in-sample fallback")
  list(pred=pred,features=features,preprocess=X$audit,
    coefficient=data.frame(variable=c("(Intercept)",colnames(X$train)),beta=b))
}
# Every training proxy prediction excludes that participant's entire group from
# ranking, proxy fitting and lambda tuning. Yang never includes test groups.
c5_crossfit_proxy<-function(train,test,yang,features,target,k,basic,seed){
  fd<-c5_group_folds(train$.group,5,seed);oof<-rep(NA_real_,nrow(train))
  fitone<-function(a,b,extra,s){
    learn<-bind_rows(a,extra);rank<-c5_domain_rank(learn,features,target,basic)
    panel<-head(rank$feature,k)
    if(length(panel)!=k)stop("Insufficient eligible proxy panel")
    c5_gaussian_proxy(learn,b,target,panel,s)
  }
  for(j in seq_len(5)){
    hold<-train$.group[fd==j];extra<-yang[!yang$.group%in%hold,,drop=FALSE]
    oof[fd==j]<-fitone(train[fd!=j,,drop=FALSE],train[fd==j,,drop=FALSE],extra,seed+j)$pred
  }
  obj<-fitone(train,test,yang,seed+20)
  if(any(!is.finite(oof)))stop("Incomplete cross-fitted proxy; cannot fit disease model")
  list(train=oof,test=obj$pred,final=obj)
}
c5_inflammation_score<-function(train,test,map){
  a<-toupper(map$assay);hit<-map$feature[map$layer=="prot"&a%in%c("GDF15","IL6","IL-6")]
  if(length(hit)<2)return(NULL)
  X<-le8_prepare_prediction_matrix(train,test,hit)
  if(ncol(X$train)!=length(hit))return(NULL)
  ztr<-rowMeans(X$train);zte<-rowMeans(X$test)
  threshold<-unname(quantile(ztr,.75))
  list(group=ifelse(zte>=threshold,"higher inflammatory-marker burden","lower inflammatory-marker burden"),threshold=threshold)
}
c5_paired_delta<-function(a,b,horizon,B,seed){
  z<-merge(a,b,by=c("eid","fold","time","event"),suffixes=c(".a",".b"))
  if(nrow(z)!=nrow(a)||nrow(z)!=nrow(b))return(tibble(status="incomplete common predictions"))
  stat<-function(ix){w<-le8_ipcw(z$time[ix],z$event[ix],horizon)
    if(w$status!="ok")return(c(AUC=NA_real_,Brier=NA_real_))
    c(AUC=le8_weighted_auc(z$risk.a[ix],w$y,w$w)-le8_weighted_auc(z$risk.b[ix],w$y,w$w),
      Brier=mean(w$w*((w$y-z$risk.a[ix])^2-(w$y-z$risk.b[ix])^2)))}
  est<-stat(seq_len(nrow(z)));set.seed(seed)
  # Resample independent groups, retaining repeated relatives together.
  groups<-split(seq_len(nrow(z)),z$group.a);ng<-length(groups)
  bs<-if(B>0)replicate(B,stat(unlist(groups[sample.int(ng,ng,replace=TRUE)],use.names=FALSE)))else matrix(NA_real_,2,0)
  ci<-function(v,p){v<-v[is.finite(v)];if(length(v)>=max(20,.8*B))unname(quantile(v,p))else NA_real_}
  tibble(status=if(all(is.finite(est)))"ok"else"insufficient supported outcome follow-up",delta_AUC=est[1],AUC_lo=ci(bs[1,],.025),AUC_hi=ci(bs[1,],.975),
    delta_Brier=est[2],Brier_lo=ci(bs[2,],.025),Brier_hi=ci(bs[2,],.975),
    uncertainty="paired group bootstrap of frozen out-of-fold predictions; training uncertainty excluded")
}
