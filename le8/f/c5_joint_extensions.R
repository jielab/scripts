# FINAL additions to joint validation. No full-cohort C1 statistics enter here.

c5_endpoint_coverage<-function(dat,Y){
  dat$.entry<-0;mf<-Sys.getenv("LE8_ENDPOINT_MANIFEST","")
  status<-"No endpoint coverage manifest supplied; baseline coverage must be verified"
  if(nzchar(mf)){
    m<-read.csv(mf,stringsAsFactors=FALSE)
    if(!all(c("Y","coverage_start","coverage_end","prebaseline_capture")%in%names(m)))stop("Incomplete endpoint manifest")
    m<-m[m$Y==Y,,drop=FALSE]
    if(nrow(m)>1)stop("Duplicate endpoint manifest rows")
    if(nrow(m)==1){
      a<-as.Date(m$coverage_start);b<-as.Date(m$coverage_end)
      if(is.na(a)||is.na(b)||b<=a)stop("Verified registry start/end dates are required")
      baseline<-as.Date(dat$date_attend)
      dat$.entry<-pmax(0,as.numeric(a-baseline)/365.25)
      end<-as.numeric(b-baseline)/365.25
      dat$.event<-as.integer(dat$.event==1&dat$.time<=end)
      dat$.time<-pmin(dat$.time,end)
      status<-paste("Only participants already covered at each landmark enter that risk set; prebaseline capture",m$prebaseline_capture)
    }
  }
  write.csv(data.frame(Y,status),file.path(out.base,"c5.endpoint_coverage.csv"),row.names=FALSE)
  dat
}

c5_load_biomarker_pgs<-function(dat,map){
  mf<-Sys.getenv("C5_OMIC_PGS_MANIFEST","")
  if(nzchar(mf)){
    # Explicit provenance is preferred; keep independent inputs on their
    # original strict path. Automatic existing UKB PGS is exploratory below.
    result<-c5_attach_omic_pgs(dat,map)
    write.csv(read.csv(mf,stringsAsFactors=FALSE),file.path(out.base,"c5.biomarker_PGS_provenance.csv"),row.names=FALSE)
    return(result)
  }
  source(file.path(fdir,"c1_pgs.R"))
  inputs<-character();audit<-list();map$pgs_feature<-NA_character_
  for(layer in c("prot","met")){
    f<-file.path(indir,"Rdata",paste0(layer,".pgs.rds"))
    if(!file.exists(f)){audit[[layer]]<-data.frame(layer,status="not supplied",mapped=0);next}
    g<-c5_assert_ids(as.data.frame(read_c1_pgs(f)))
    g<-g[g$eid%in%dat$eid,,drop=FALSE]
    mi<-which(map$layer==layer);lookup<-map_c1_pgs_columns(map$assay[mi],names(g))
    hit<-mi[map$assay[mi]%in%names(lookup)];cols<-unname(lookup[map$assay[hit]])
    if(!length(cols)){audit[[layer]]<-data.frame(layer,status="no matched scores",mapped=0);next}
    if(anyDuplicated(cols)||!all(vapply(g[cols],is.numeric,logical(1))))stop("Ambiguous or nonnumeric matched PGS")
    nf<-paste0("pgs__",map$feature[hit]);map$pgs_feature[hit]<-nf
    g<-g[,c("eid",cols),drop=FALSE];names(g)<-c("eid",nf)
    # Retain the measured common cohort. PGS missingness never becomes zero.
    dat<-merge(dat,g,by="eid",all.x=TRUE,sort=FALSE)
    inputs<-c(inputs,f);audit[[layer]]<-data.frame(layer,status="existing biomarker PGS; discovery overlap unverified; exploratory",mapped=length(hit))
  }
  if(all(is.na(map$pgs_feature)))map$pgs_feature<-NULL
  write.csv(bind_rows(audit),file.path(out.base,"c5.biomarker_PGS_provenance.csv"),row.names=FALSE)
  list(data=dat[order(dat$eid),,drop=FALSE],map=map,inputs=inputs,
    status=if(length(inputs))"exploratory: automatic biomarker PGS; discovery independence unverified"else"biomarker PGS unavailable")
}

c5_replicated_domains<-function(train,features,targets,basic,seed){
  # Stable donor hashing: adding Yang must not move existing Yin donors
  # between the discovery and replication halves.
  stable_half<-function(g){v<-as.double(seed);for(x in utf8ToInt(as.character(g)))v<-(v*131+x)%%2147483629;as.integer(v%%2)+1L}
  h<-vapply(train$.group,stable_half,integer(1))
  ranks<-list();audit<-list();minr<-as.numeric(Sys.getenv("C5_DOMAIN_MIN_R","0.10"))
  if(!is.finite(minr)||minr<0||minr>1)stop("C5_DOMAIN_MIN_R must be in [0,1]")
  for(tg in targets){
    if(!tg%in%names(train)){audit[[tg]]<-tibble(domain=tg,status="missing target",eligible=0L);next}
    one<-function(hh){
      z<-tryCatch(c5_domain_rank(train[h==hh,,drop=FALSE],features,tg,basic),error=function(e)tibble())
      if(!nrow(z))return(z)
      df<-pmax(1,z$N_observed-length(basic)-2)
      z$p<-2*pt(abs(z$r)*sqrt(df/pmax(1e-12,1-z$r^2)),df=df,lower.tail=FALSE)
      z$q<-p.adjust(z$p,"BH",n=length(features));z
    }
    a<-one(1);b<-one(2)
    if(!nrow(a)||!nrow(b)){audit[[tg]]<-tibble(domain=tg,status="insufficient training target",eligible=0L);next}
    z<-inner_join(a,b,by=c("feature","domain"),suffix=c("1","2"))|>
      mutate(eligible=q1<.05&q2<.05&sign(r1)==sign(r2)&pmin(abs(r1),abs(r2))>=minr,
        strength=pmin(abs(r1),abs(r2)))|>arrange(desc(strength),feature)
    ranks[[tg]]<-z;audit[[tg]]<-tibble(domain=tg,status=if(any(z$eligible))"supported in both training halves"else"no eligible proxies",eligible=sum(z$eligible))
  }
  list(ranks=lapply(ranks,function(z)z[z$eligible,,drop=FALSE]),all=bind_rows(ranks),audit=bind_rows(audit))
}

# Orthogonalization is a statistical decomposition, not causal partitioning.
# The genetic-predicted part uses only that biomarker's own PGS. The remainder
# includes acquired biology, untagged genetic effects, confounding and error.
c5_pgs_decompose<-function(train,test,features,map,seed){
  tr<-train;te<-test;gp<-resid<-character();aud<-coef<-list()
  folds<-c5_group_folds(train$.group,5,seed)
  for(f in features){
    g<-map$pgs_feature[match(f,map$feature)]
    if(length(g)!=1||is.na(g)||!g%in%names(train))next
    observed<-is.finite(train[[f]])&is.finite(train[[g]])
    if(sum(observed)<100||sd(train[[g]][observed])<=0)next
    fit_predict<-function(ii,dd){
      ok<-observed&ii
      if(sum(ok)<50)return(NULL)
      fit<-lm.fit(cbind(1,train[[g]][ok]),train[[f]][ok]);b<-fit$coefficients
      if(any(!is.finite(b)))return(NULL)
      list(pred=b[1]+b[2]*dd[[g]],coef=b)
    }
    oof<-rep(NA_real_,nrow(train))
    for(h in seq_len(5)){
      z<-fit_predict(folds!=h,train[folds==h,,drop=FALSE])
      if(!is.null(z))oof[folds==h]<-z$pred
    }
    obj<-fit_predict(rep(TRUE,nrow(train)),test)
    if(is.null(obj)||any(!is.finite(oof[observed])))next
    gn<-paste0("genpart__",f);rn<-paste0("remainder__",f)
    tr[[gn]]<-oof;te[[gn]]<-obj$pred
    tr[[rn]]<-train[[f]]-oof;te[[rn]]<-test[[f]]-obj$pred
    gp<-c(gp,gn);resid<-c(resid,rn)
    ok<-is.finite(test[[f]])&is.finite(obj$pred)
    obs<-test[[f]][ok];pr<-obj$pred[ok];den<-sum((obs-mean(train[[f]][observed]))^2)
    aud[[f]]<-tibble(feature=f,pgs=g,N_training=sum(observed),N_validation=sum(ok),
      R2=if(sum(ok)>2&&den>0)1-sum((obs-pr)^2)/den else NA_real_,
      measured_PGS_r=if(sum(ok)>2)cor(obs,pr)else NA_real_,
      interpretation="Prediction of adult measurement from matched PGS; remainder is not an environmental or causal effect")
    coef[[f]]<-tibble(feature=f,pgs=g,intercept=obj$coef[1],slope=obj$coef[2])
  }
  list(train=tr,test=te,genetic=gp,remainder=resid,audit=bind_rows(aud),coefficient=bind_rows(coef))
}

c5_cross_omic_links<-function(train,test,features,basic){
  pp<-features[startsWith(features,"prot__")];mm<-features[startsWith(features,"met__")]
  if(!length(pp)||!length(mm))return(tibble())
  getr<-function(d){
    cv<-le8_prepare_prediction_matrix(d,d,basic)$train;q<-qr(cbind(1,cv))
    # Pairwise observed assays, with basic-covariate residuals; training and
    # held-out cohorts are kept separate and all requested pairs are retained.
    residual<-function(v){ok<-is.finite(v);ans<-rep(NA_real_,length(v))
      if(sum(ok)>length(basic)+5)ans[ok]<-lm.fit(cbind(1,cv[ok,,drop=FALSE]),v[ok])$residuals;ans}
    a<-sapply(d[,pp,drop=FALSE],residual);b<-sapply(d[,mm,drop=FALSE],residual)
    if(is.null(dim(a)))a<-matrix(a,ncol=1,dimnames=list(NULL,pp))
    if(is.null(dim(b)))b<-matrix(b,ncol=1,dimnames=list(NULL,mm))
    map_dfr(pp,function(p)map_dfr(mm,function(m){
      ok<-is.finite(a[,p])&is.finite(b[,m]);n<-sum(ok)
      r<-if(n>20)cor(a[ok,p],b[ok,m])else NA_real_
      df<-max(1,n-length(basic)-2);pv<-if(is.finite(r))2*pt(abs(r)*sqrt(df/max(1e-12,1-r*r)),df,lower.tail=FALSE)else NA_real_
      tibble(protein=p,metabolite=m,r=r,N=n,p=pv)
    }))
  }
  z<-inner_join(getr(train),getr(test),by=c("protein","metabolite"),suffix=c("_training","_validation"))
  z$q_validation<-p.adjust(z$p_validation,"BH")
  z$replicated_direction<-sign(z$r_training)==sign(z$r_validation)
  z$interpretation<-"Partial correlation; no directionality or protein-to-metabolite mediation established"
  z
}

c5_marker_strata<-function(train,test,map){
  definitions<-list(GDF15_IL6_median=list(assay=c("GDF15","IL6"),q=.5),
    GDF15_IL6_q75=list(assay=c("GDF15","IL6"),q=.75),
    GlycA_median=list(assay="GLYCA",q=.5))
  ans<-data.frame(eid=test$eid);audit<-list()
  for(nm in names(definitions)){
    spec<-definitions[[nm]];hit<-map$feature[match(spec$assay,toupper(map$assay))]
    if(anyNA(hit)){ans[[nm]]<-"unavailable";audit[[nm]]<-tibble(definition=nm,status="required marker absent");next}
    mu<-vapply(train[hit],mean,numeric(1),na.rm=TRUE);ss<-vapply(train[hit],sd,numeric(1),na.rm=TRUE)
    if(any(!is.finite(ss)|ss<=0)){ans[[nm]]<-"unavailable";next}
    sc<-function(d)rowMeans(sweep(sweep(as.matrix(d[hit]),2,mu),2,ss,"/"))
    tr<-sc(train);te<-sc(test);cut<-unname(quantile(tr,spec$q,na.rm=TRUE))
    ans[[nm]]<-ifelse(!is.finite(te),"unavailable",ifelse(te<=cut,"lower","higher"))
    audit[[nm]]<-tibble(definition=nm,assay=hit,center=mu,sd=ss,cutoff=cut,status="training threshold; complete-marker score")
  }
  list(groups=ans,audit=bind_rows(audit))
}

c5_marker_heterogeneity<-function(a,b,groups,horizon,B,seed){
  z<-merge(a,b,by=c("eid","fold","time","event"),suffixes=c(".a",".b"))
  if(nrow(z)!=nrow(a)||nrow(z)!=nrow(b))return(tibble())
  z<-merge(z,groups,by="eid");out<-list()
  for(nm in setdiff(names(groups),"eid")){
    d<-z[z[[nm]]%in%c("lower","higher"),,drop=FALSE]
    if(!all(c("lower","higher")%in%d[[nm]]))next
    stat<-function(ix){sapply(c("lower","higher"),function(g){
      ii<-ix[d[[nm]][ix]==g];w<-le8_ipcw(d$time[ii],d$event[ii],horizon)
      if(w$status!="ok")return(NA_real_)
      le8_weighted_auc(d$risk.a[ii],w$y,w$w)-le8_weighted_auc(d$risk.b[ii],w$y,w$w)})}
    est<-stat(seq_len(nrow(d)));set.seed(seed)
    gg<-split(seq_len(nrow(d)),d$group.a)
    bs<-if(B>0)replicate(B,stat(unlist(gg[sample.int(length(gg),length(gg),replace=TRUE)],use.names=FALSE)))else matrix(NA_real_,2,0)
    dd<-bs[2,]-bs[1,];dd<-dd[is.finite(dd)];diff<-unname(est[2]-est[1]);se<-sd(dd)
    out[[nm]]<-tibble(definition=nm,N=nrow(d),delta_lower=est[1],delta_higher=est[2],difference=diff,
      lo=if(length(dd)>=20)unname(quantile(dd,.025))else NA_real_,hi=if(length(dd)>=20)unname(quantile(dd,.975))else NA_real_,
      p=if(length(dd)>=20&&is.finite(se)&&se>0)2*pnorm(abs(diff/se),lower.tail=FALSE)else NA_real_,
      interpretation="Exploratory difference in delta AUC, group bootstrap of frozen fits; marker burden is not a causal subtype")
  }
  bind_rows(out)
}

c5_genetic_measured_strata<-function(pred,budget){
  rows<-list()
  for(L in unique(pred$landmark))for(arm in unique(pred$arm)){
    a<-pred[pred$model=="MatchedMeasured_only_NS"&pred$budget==budget&pred$landmark==L&pred$arm==arm,,drop=FALSE]
    b<-pred[pred$model=="MatchedPGS_only_NS"&pred$budget==budget&pred$landmark==L&pred$arm==arm,,drop=FALSE]
    if(!nrow(a)||!nrow(b)||!setequal(a$eid,b$eid))next
    z<-merge(a,b[,c("eid","high_training_q75")],by="eid",suffixes=c(".measured",".PGS"))
    z$stratum<-paste0("measured ",ifelse(z$high_training_q75.measured,"high","lower")," / PGS ",ifelse(z$high_training_q75.PGS,"high","lower"))
    for(g in unique(z$stratum)){
      d<-z[z$stratum==g,,drop=FALSE];w<-le8_ipcw(d$time,d$event,d$horizon[1])
      rows[[length(rows)+1L]]<-tibble(stratum=g,landmark=L,arm=arm,N=nrow(d),events=w$N_case,
        observed_risk=if(w$status=="ok")sum(w$w*w$y)/sum(w$w)else NA_real_,status=w$status,
        interpretation="Training upper-quartile thresholds; same matched biomarkers; acquired or genetic causation is not identified")
    }
  }
  bind_rows(rows)
}
