suppressPackageStartupMessages({library(data.table);library(ggplot2)})
a=commandArgs(trailingOnly=TRUE)
arg=function(k,d=NULL){i=match(k,a);if(is.na(i))d else a[i+1]}
trait=tolower(arg('--trait')); phefile=arg('--phe'); ancestryfile=arg('--ancestry'); pcafile=arg('--pca',''); scoredir=arg('--score-dir'); outdir=arg('--out-dir'); phenocol=arg('--phenotype-col','auto'); covspec=arg('--covariates','auto'); B=as.integer(arg('--bootstrap','100')); K=as.integer(arg('--folds','5')); ming=as.integer(arg('--min-group','100')); tuned=toupper(arg('--tuned','TRUE'))=='TRUE'; writepred=toupper(arg('--write-predictions','FALSE'))=='TRUE'; seed=as.integer(arg('--seed','20260904'));set.seed(seed)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
read_any=function(f){if(!nzchar(f)||!file.exists(f))return(NULL);if(grepl('\\.rds$',f,ignore.case=TRUE))return(readRDS(f));fread(f)}
find_df=function(x){
 if(inherits(x,c('data.frame','data.table')))return(as.data.table(x))
 if(is.list(x)){z=lapply(x,find_df);z=z[!vapply(z,is.null,logical(1))];if(length(z))return(z[[which.max(vapply(z,nrow,integer(1)))]])}
 NULL
}
norm_id=function(d){id=intersect(c('eid','IID','#IID','ID_2','f.eid','id'),names(d));if(!length(id))return(NULL);setnames(d,id[1],'eid');d[,eid:=as.character(eid)];d}
phe=norm_id(find_df(readRDS(phefile)));if(is.null(phe))stop('Cannot identify participant table/ID in ',phefile)
# Outcome selection deliberately writes candidates before failing.
select_outcome=function(d,trait,requested){
 if(requested!='auto'){if(!requested%in%names(d))stop('Requested phenotype column missing: ',requested);return(requested)}
 aliases=switch(trait,
  height=c('height','standing_height','standing.height','f.50.0.0','f_50_0_0','50-0.0'),
  ldl=c('ldl','ldl_c','ldl.direct','ldl_direct','direct_ldl','f.30780.0.0','f_30780_0_0','30780-0.0'),
  t2dm=c('t2dm','t2d','type2_diabetes','type_2_diabetes','dm2','diabetes_type2','incident_t2dm','prevalent_t2dm'),character())
 low=tolower(names(d));for(x in aliases){j=which(low==tolower(x));if(length(j))return(names(d)[j[1]])}
 pat=switch(trait,height='(^|[._])height($|[._])|standing.*height',ldl='(^|[._])ldl($|[._])|low.*density.*lipoprotein',t2dm='t2d|type.?2.*diabet|diabet.*type.?2',trait)
 cand=names(d)[grepl(pat,names(d),ignore.case=TRUE)];fwrite(data.table(candidate=cand),file.path(outdir,'phenotype_candidates.tsv'),sep='\t')
 if(length(cand)==1)return(cand);stop('Could not uniquely identify ',trait,' outcome. See phenotype_candidates.tsv and rerun with --phenotype-col.')
}
ycol=select_outcome(phe,trait,phenocol)
# Normalize binary T2DM; quantitative traits remain numeric.
if(trait=='t2dm'){
 x=phe[[ycol]]
 if(is.logical(x))y=as.integer(x) else if(is.factor(x)||is.character(x)){z=tolower(as.character(x));y=ifelse(z%in%c('1','yes','y','case','true','t2dm','type 2 diabetes'),1,ifelse(z%in%c('0','no','n','control','false','none'),0,NA))} else {y=as.numeric(x);u=sort(unique(y[is.finite(y)]));if(!all(u%in%c(0,1)))y=as.integer(y>0)}
 phe[,outcome:=y];family_type='binary'
}else{phe[,outcome:=as.numeric(get(ycol))];family_type='quantitative'}
# Merge ancestry and PCA/covariates.
anc=norm_id(as.data.table(read_any(ancestryfile)));if(is.null(anc))stop('No ancestry table')
pc=find_df(read_any(pcafile));if(!is.null(pc))pc=norm_id(pc)
d=merge(phe,anc,by='eid',all=FALSE,suffixes=c('','.anc'));if(!is.null(pc))d=merge(d,pc,by='eid',all.x=TRUE,suffixes=c('','.pc'))
# Merge score files.
scorefiles=c('csx.tsv.gz','csx_mixed.tsv.gz','disco_zero.tsv.gz','grid.tsv.gz')
for(f in scorefiles){p=file.path(scoredir,f);if(file.exists(p)){s=norm_id(fread(p));newcols=setdiff(names(s),names(d));d=merge(d,s[,c('eid',newcols),with=FALSE],by='eid',all.x=TRUE)}}
# Create matched CSX if mixed file was not present.
pops=c('AFR','EAS','EUR','SAS');qcols=paste0('posterior_',pops);csxcols=paste0('CSX_',pops)
if(all(csxcols%in%names(d))&&!('CSX_posterior'%in%names(d))){
 if(all(qcols%in%names(d)))Q=as.matrix(d[,..qcols]) else {ac=intersect(c('ancestry','genetic_ancestry','predicted_ancestry'),names(d));if(!length(ac))stop('No ancestry label');Q=sapply(pops,function(p)as.numeric(toupper(d[[ac[1]]])==p))}
 Q[!is.finite(Q)]=0;den=rowSums(Q);Q=Q/pmax(den,1);Q[den==0,]=.25;S=as.matrix(d[,..csxcols]);d[,CSX_posterior:=rowSums(S*Q,na.rm=TRUE)];d[,CSX_matched:=S[cbind(seq_len(.N),max.col(Q,ties.method='first'))]]
}
ac=intersect(c('ancestry','genetic_ancestry','predicted_ancestry'),names(d));if(length(ac))d[,eval_ancestry:=toupper(as.character(get(ac[1])))] else d[,eval_ancestry:='ALL']
# Automatic covariates: age, sex, genotyping array/center, and available PCs.
if(covspec!='auto'){covars=intersect(strsplit(covspec,'[,; ]+')[[1]],names(d))}else{
 exact=c('age','age_enroll','baseline_age','f.21003.0.0','sex','f.31.0.0','assessment_centre','center','genotyping_array','array')
 covars=unique(c(names(d)[tolower(names(d))%in%tolower(exact)],grep('^PC([1-9]|1[0-9]|20)$',names(d),value=TRUE,ignore.case=TRUE)))
 covars=setdiff(covars,c(ycol,'outcome'))
}
# Remove unusable/duplicated covariates.
covars=covars[vapply(d[,..covars],function(x)length(unique(x[!is.na(x)]))>1,logical(1))];covars=unique(covars)
fwrite(data.table(trait=trait,outcome_column=ycol,family=family_type,covariates=paste(covars,collapse=','),N=nrow(d)),file.path(outdir,'evaluation_manifest.tsv'),sep='\t')
# OOF phenotype-tuned stacking comparator. It is labelled separately and never used to construct GRID.
if(tuned&&all(csxcols%in%names(d))){
 ok=complete.cases(d[,c('outcome',covars,csxcols),with=FALSE]);idx=which(ok);fold=rep(NA_integer_,nrow(d));
 if(family_type=='binary'){for(v in 0:1){ii=idx[d$outcome[idx]==v];fold[ii]=sample(rep(seq_len(K),length.out=length(ii)))}}else fold[idx]=sample(rep(seq_len(K),length.out=length(idx)))
 gs=rep(NA_real_,nrow(d))
 for(k in seq_len(K)){tr=idx[fold[idx]!=k];te=idx[fold[idx]==k];if(!length(te))next
  fb=as.formula(paste('outcome~',if(length(covars))paste(covars,collapse='+') else '1'));ff=as.formula(paste('outcome~',paste(c(covars,csxcols),collapse='+')))
  if(family_type=='binary'){mb=glm(fb,d[tr],family=binomial());mf=glm(ff,d[tr],family=binomial());gs[te]=predict(mf,d[te],type='link')-predict(mb,d[te],type='link')}else{mb=lm(fb,d[tr]);mf=lm(ff,d[tr]);gs[te]=predict(mf,d[te])-predict(mb,d[te])}
 }
 d[,CSX_tuned_stack:=gs]
}
methods=intersect(c('CSX_matched','CSX_posterior','DiscoDivas_zero','GRID_shared','GRID_matched','GRID_posterior','CSX_tuned_stack'),names(d))
if(!length(methods))stop('No evaluable PRS methods found in ',scoredir)
auc_rank=function(y,p){ok=is.finite(y)&is.finite(p);y=y[ok];p=p[ok];n1=sum(y==1);n0=sum(y==0);if(n1<2||n0<2)return(NA_real_);(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)}
metric_one=function(x,m){
 z=x[complete.cases(x[,c('outcome',covars,m),with=FALSE])];if(nrow(z)<ming)return(NULL);z[,score_z:=(as.numeric(get(m))-mean(as.numeric(get(m))))/sd(as.numeric(get(m)))]
 if(!is.finite(sd(z$score_z)))return(NULL);fb=as.formula(paste('outcome~',if(length(covars))paste(covars,collapse='+') else '1'));ff=as.formula(paste('outcome~',paste(c(covars,'score_z'),collapse='+')))
 if(family_type=='quantitative'){
  b=lm(fb,z);f=lm(ff,z);sm=summary(f);co=coef(sm)['score_z',];pred=predict(f,z);data.table(N=nrow(z),cases=NA_integer_,metric='incremental_R2',estimate=summary(f)$r.squared-summary(b)$r.squared,beta=co[1],se=co[2],p=co[4],secondary_metric='RMSE',secondary_value=sqrt(mean((z$outcome-pred)^2)))
 }else{
  b=glm(fb,z,family=binomial());f=glm(ff,z,family=binomial());co=coef(summary(f))['score_z',];pb=predict(b,z,type='response');pf=predict(f,z,type='response');data.table(N=nrow(z),cases=sum(z$outcome==1),metric='delta_AUC',estimate=auc_rank(z$outcome,pf)-auc_rank(z$outcome,pb),beta=co[1],se=co[2],p=co[4],secondary_metric='AUC_full',secondary_value=auc_rank(z$outcome,pf),OR_per_SD=exp(co[1]),Brier=mean((z$outcome-pf)^2))
 }
}
groups=c('ALL',sort(intersect(pops,unique(d$eval_ancestry))));res=list();boots=list()
for(g in groups){x=if(g=='ALL')d else d[eval_ancestry==g];if(nrow(x)<ming)next
 for(m in methods){r=metric_one(x,m);if(is.null(r))next;r[,`:=`(group=g,method=m)];res[[length(res)+1]]=r
  if(B>0){bv=rep(NA_real_,B);for(b in seq_len(B)){xx=x[sample.int(nrow(x),replace=TRUE)];rr=try(metric_one(xx,m),silent=TRUE);if(!inherits(rr,'try-error')&&!is.null(rr))bv[b]=rr$estimate[1]};boots[[length(boots)+1]]=data.table(group=g,method=m,boot=seq_len(B),estimate=bv)}
 }
}
perf=rbindlist(res,fill=TRUE);if(!nrow(perf))stop('No group/method passed minimum sample size')
if(length(boots)){bt=rbindlist(boots);ci=bt[,.(ci_low=quantile(estimate,.025,na.rm=TRUE),ci_high=quantile(estimate,.975,na.rm=TRUE)),by=.(group,method)];perf=merge(perf,ci,by=c('group','method'),all.x=TRUE);fwrite(bt,file.path(outdir,'performance_bootstrap.tsv.gz'),sep='\t')}
fwrite(perf,file.path(outdir,'performance.tsv'),sep='\t')
# Figure 1: main zero-shot/tuned comparison.
p=ggplot(perf,aes(x=estimate,y=reorder(method,estimate)))+geom_vline(xintercept=0,lty=2)+geom_point()+geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.15,na.rm=TRUE)+facet_wrap(~group,scales='free_y')+labs(x=unique(perf$metric),y=NULL,title=toupper(trait),subtitle='Zero-shot methods and phenotype-tuned stacking are labelled separately')+theme_bw(base_size=11)
ggsave(file.path(outdir,'eval.Fig1.performance.png'),p,width=11,height=7,dpi=300)
# Figure 2: non-EUR performance relative to EUR, when possible.
e=perf[group%in%pops,.(group,method,estimate)];if(nrow(e)){eur=e[group=='EUR',.(method,EUR=estimate)];e=merge(e,eur,by='method');e[,relative_to_EUR:=estimate/EUR];p2=ggplot(e[group!='EUR'],aes(x=group,y=relative_to_EUR,group=method,shape=method))+geom_hline(yintercept=1,lty=2)+geom_point(position=position_dodge(width=.4))+labs(x=NULL,y='Performance / EUR performance',title=paste(toupper(trait),'cross-ancestry relative performance'))+theme_bw(base_size=11);ggsave(file.path(outdir,'eval.Fig2.relative.png'),p2,width=9,height=6,dpi=300)}
if(requireNamespace('openxlsx',quietly=TRUE)){wb=openxlsx::createWorkbook();openxlsx::addWorksheet(wb,'performance');openxlsx::writeData(wb,'performance',perf);openxlsx::addWorksheet(wb,'manifest');openxlsx::writeData(wb,'manifest',fread(file.path(outdir,'evaluation_manifest.tsv')));openxlsx::saveWorkbook(wb,file.path(outdir,'GRID_evaluation.xlsx'),overwrite=TRUE)}
if(writepred)fwrite(d[,c('eid','outcome','eval_ancestry',methods),with=FALSE],file.path(outdir,'predictions.tsv.gz'),sep='\t')
cat('outcome=',ycol,' N=',nrow(d),' methods=',paste(methods,collapse=','),'\n',sep='')
