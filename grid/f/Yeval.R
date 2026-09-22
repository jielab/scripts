# Paired prediction evaluation. SNP weights are fixed; all outcome fitting is OOF.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(survival);library(pROC);library(patchwork)})
setDTthreads(4)
args <- commandArgs(TRUE)
allowed <- c('trait','type','method','score-dir','pgs-file','disco-file','pt-file','pheno-file',
 'ancestry-file','group-col','covar-name','phenotype-col','event-col','time-col','prevalence',
 'pca-file','med-file','distance-pcs','distance-bins','folds','seed','bootstrap','min-n',
 'remove','out-root','dir-gwas','dir-gen','pt-effect','threads','check','allow-missing-scores','disco-tune','disco-a','min-anchor','write-predictions','run-dir','grid-file')
opt <- list(); i <- 1L
while(i <= length(args)) {
 key <- sub('^--','',args[i]); if(!key %in% allowed)stop('Unknown option: ',args[i])
 if(key %in% c('check','allow-missing-scores')) {opt[[key]] <- TRUE;i <- i+1L} else {
  if(i==length(args))stop('Missing value: ',args[i]);opt[[key]] <- args[i+1L];i <- i+2L
 }
}
arg <- function(k,default=NULL)if(is.null(opt[[k]]))default else opt[[k]]
Y <- arg('trait'); type <- arg('type')
if(is.null(Y)||!type %in% c('ct','dt','t2e'))stop('--trait and --type ct|dt|t2e required')
if(!arg('method','all') %in% c('all','csx','disco'))stop('Invalid --method')
out <- arg('run-dir',file.path(arg('out-root','/mnt/d/analysis/grid/Yeval'),Y))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,name)fwrite(x,file.path(out,name),sep='\t',na='NA')
intarg <- function(k,v,minimum) {z<-suppressWarnings(as.integer(arg(k,v)));if(is.na(z)||z<minimum)stop('Invalid --',k);z}
nfold <- intarg('folds',5,2); seed <- intarg('seed',20260904,0)
nboot <- intarg('bootstrap',200,0); minn <- intarg('min-n',100,20)
npc <- intarg('distance-pcs',10,2); nbins <- intarg('distance-bins',5,2)
tune <- toupper(arg('disco-tune','TRUE'))
if(!tune %in% c('TRUE','FALSE'))stop('--disco-tune must be TRUE/FALSE')
tune<-tune=='TRUE'; min_anchor<-intarg('min-anchor',100,10)
quality<-as.numeric(strsplit(arg('disco-a','1,1,1,1'),',',fixed=TRUE)[[1]])
if(length(quality)!=4||any(!is.finite(quality)|quality<0)||!any(quality>0))stop('Invalid --disco-a (AFR,EAS,EUR,SAS order)')
partial <- isTRUE(arg('allow-missing-scores',FALSE))
pops <- c('EUR','AFR','EAS','SAS'); base_scores <- paste0('csx.',c('AFR','EAS','EUR','SAS'))
quality <- setNames(quality,c('AFR','EAS','EUR','SAS'))[pops]
main_methods <- c('COJO','PRS-CSx-auto-meta','PRS-CSx',if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned')
score_dir <- file.path(arg('score-dir','/mnt/d/data/ukb/pgs'),Y)
files <- c(phenotype=arg('pheno-file','/mnt/d/data/ukb/phe/Rdata/all.rds'),
 csx=arg('pgs-file',file.path(score_dir,'csx.pgs.gz')),
 disco=arg('disco-file',file.path(score_dir,'disco.pgs.gz')),
 pt=arg('pt-file',file.path(score_dir,'pt.pgs.gz')),
 ancestry=arg('ancestry-file','/mnt/d/data/ukb/pca_proj/ukb.ancestry.auto.tsv.gz'),
 pca=arg('pca-file','/mnt/d/data/ukb/pca_proj/ukb.discodivas.pca.tsv.gz'),
 centers=arg('med-file','/mnt/d/files/DiscoDivas/med.g1000.4pop.tsv'))
read_table <- function(f) {
 if(!file.exists(f))stop('Missing input: ',f)
 if(grepl('\\.rds$',f,ignore.case=TRUE))return(as.data.table(readRDS(f)))
 header<-names(fread(f,nrows=0))
 fread(f,colClasses=list(character=intersect(header,c('eid','IID','#IID','FID','#FID'))),na.strings=c('NA','NaN',''))
}
ids <- function(d) {
 nm<-intersect(c('eid','IID','#IID','ID_2'),names(d));if(!length(nm))stop('Missing sample ID')
 if(nm[1]!='eid')setnames(d,nm[1],'eid')
 d[,eid:=as.character(eid)];if(anyNA(d$eid)||anyDuplicated(d$eid))stop('Missing/duplicate sample IDs')
 d
}
num <- function(x,name) {y<-suppressWarnings(as.numeric(as.character(x)));if(any(!is.na(x)&is.na(y)))stop('Non-numeric: ',name);y}
cat('Reading full phenotype cohort: ',files['phenotype'],'\n',sep='')
phe <- ids(read_table(files['phenotype']))
gc <- arg('group-col','genetic_ancestry')
if(!gc %in% names(phe)) {
 anc<-ids(read_table(files['ancestry']));if(!gc %in% names(anc))stop('Missing ancestry column: ',gc)
 phe<-merge(phe,anc[,c('eid',gc),with=FALSE],by='eid',all.x=TRUE)
}
phe[,target:=as.character(get(gc))]
phe[is.na(target)|!nzchar(target),target:='UNASSIGNED']
rem <- arg('remove','/mnt/d/files/ukb.exclude.id'); excluded <- character()
if(nzchar(rem)) {
 if(!file.exists(rem))stop('Missing withdrawal file: ',rem)
 if(file.info(rem)$size>0) {r<-fread(rem,header=FALSE,colClasses='character');excluded<-r[[min(2L,ncol(r))]]}
}
phe <- phe[!startsWith(eid,'-')&!eid %in% excluded]
covars <- trimws(strsplit(arg('covar-name','age,sex,PC1,PC2'),',',fixed=TRUE)[[1]])
covars[covars=='PC']<-'PC1';covars<-unique(covars[nzchar(covars)&covars!='none'])
yc <- arg('phenotype-col',Y); ec <- arg('event-col',paste0(Y,'.Yt2e')); tc <- arg('time-col',paste0(Y,'.t2e'))
outcome_definition <- if(type=='t2e')paste(ec,tc) else yc
# Yt2e is incident-only. Baseline disease includes dated prevalent cases, and
# treats subsequently diagnosed people as baseline non-cases, never as missing.
if(type=='dt'&&Y=='t2dm'&&is.null(arg('phenotype-col'))) {
 needed<-c('t2dm.Yr2e','t2dm.Yt2e');if(!all(needed%in%names(phe)))stop('Need explicit 0/1 --phenotype-col or t2dm.Yr2e/Yt2e')
 phe[,t2dm:=fifelse(get('t2dm.Yr2e')==1,1,
                   fifelse(get('t2dm.Yt2e')%in%c(0,1)|get('t2dm.Yr2e')==0,0,NA_real_),na=NA_real_)]
 # fifelse NA condition above must not erase incident cases (Yr2e=NA).
 phe[is.na(get('t2dm.Yr2e')) & get('t2dm.Yt2e')%in%c(0,1),t2dm:=0]
 outcome_definition <- 'Baseline ICD10 T2D: Yr2e=1 case; valid nonprevalent Yt2e=0/1 control; undated/invalid NA'
}
required<-unique(c(if(type=='t2e')c(ec,tc) else yc,covars))
if(!all(required%in%names(phe)))stop('Missing phenotype/covariates: ',paste(setdiff(required,names(phe)),collapse=','))
if(length(intersect(covars,c(yc,ec,tc,base_scores,'csx.auto','csx.meta','disco'))))stop('Outcome/score cannot be a covariate')
phe[,outcome:=num(get(if(type=='t2e')ec else yc),'outcome')]
if(type!='ct'&&any(!is.na(phe$outcome)&!phe$outcome%in%c(0,1)))stop('Binary/event outcome must be 0/1')
if(type=='t2e')phe[,time:=num(get(tc),tc)] else phe[,time:=NA_real_]
# K is estimated before score matching and covariate complete-case selection.
prevalence <- phe[!is.na(target)&is.finite(outcome),.(population_N=.N,population_cases=sum(outcome==1),K=mean(outcome)),by=target]
prevalence[,source:='Full phenotype cohort, before PGS/covariate filtering (UKB cohort assumption)']
if(type=='dt') {
 spec<-arg('prevalence','cohort')
 if(spec!='cohort') {
  if(grepl('=',spec,fixed=TRUE)) {
   kv<-strsplit(strsplit(spec,',',fixed=TRUE)[[1]],'=',fixed=TRUE)
   if(any(lengths(kv)!=2))stop('Invalid prevalence mapping')
   ks<-setNames(vapply(kv,function(z)as.numeric(z[2]),numeric(1)),vapply(kv,`[`,character(1),1))
   if(anyDuplicated(names(ks)))stop('Duplicate prevalence target')
   prevalence[,K:=ks[target]]
  } else prevalence[,K:=as.numeric(spec)]
  prevalence[,source:='User-specified population prevalence']
 }
 if(any(!is.finite(prevalence[target%in%pops]$K)|prevalence[target%in%pops]$K<=0|prevalence[target%in%pops]$K>=1))stop('Each target needs prevalence 0<K<1')
}
keep <- unique(c('eid','target','outcome','time',covars))
d <- phe[,..keep]; rm(phe); invisible(gc(verbose=FALSE))
for(v in covars)if(is.character(d[[v]])||is.factor(d[[v]])||v=='sex')set(d,j=v,value=factor(d[[v]]))
audit <- list(data.table(stage='phenotype_after_withdrawal',target=d$target)[,.N,by=.(stage,target)])
pgs <- ids(read_table(files['csx']))
missing <- setdiff(c(base_scores,'csx.auto','csx.meta'),names(pgs))
if(length(missing)&&!partial)stop('Missing CSx columns: ',paste(missing,collapse=', '),'. Complete upstream scores or explicitly use --allow-missing-scores.')
available <- intersect(c(base_scores,'csx.auto','csx.meta'),names(pgs))
d <- merge(d,pgs[,c('eid',available),with=FALSE],by='eid');rm(pgs)
for(kind in c('disco','pt')) {
 sc<-if(kind=='disco')'disco' else paste0('pt.',pops)
 if(file.exists(files[kind])) {
  z<-ids(read_table(files[kind]));bad<-setdiff(sc,names(z))
  if(length(bad)&&!partial)stop('Missing ',kind,' columns: ',paste(bad,collapse=','))
  use<-intersect(sc,names(z));available<-c(available,use);missing<-c(missing,bad)
  d<-merge(d,z[,c('eid',use),with=FALSE],by='eid',all.x=TRUE)
 } else {
  missing<-c(missing,sc)
  if(!partial&&!isTRUE(arg('check',FALSE))&&!(kind=='disco'&&tune))stop('Missing ',files[kind])
 }
}
grid_scores<-character()
if(!is.null(arg('grid-file'))){
 z<-ids(read_table(arg('grid-file')));grid_scores<-intersect(c(paste0('GRID_',c('AFR','EAS','EUR','SAS')),'GRID_shared','GRID_posterior','GRID_matched'),names(z))
 if(!length(grid_scores))stop('No recognized GRID columns')
 d<-merge(d,z[,c('eid',grid_scores),with=FALSE],by='eid',all.x=TRUE);available<-c(available,grid_scores)
}
for(s in available)set(d,j=s,value=num(d[[s]],s))
# Use reference-projected PCs only for genetic distance; phenotype PCs remain covariates.
pc <- paste0('PC',seq_len(npc)); pca <- ids(read_table(files['pca'])); centers <- read_table(files['centers'])
if(!all(pc%in%names(pca))||!all(c('POP',pc)%in%names(centers)))stop('Missing reference-projected PCs/centers')
centers<-centers[match(pops,POP)]
if(anyNA(centers$POP)||anyDuplicated(read_table(files['centers'])$POP))stop('Need one center for each EUR/AFR/EAS/SAS')
pm<-as.matrix(pca[,..pc]); cm<-as.matrix(centers[,..pc])
if(any(!is.finite(pm))||any(!is.finite(cm)))stop('Nonfinite projected PCs or centers')
pc_names<-paste0('ancPC',seq_len(npc))
gd <- data.table(eid=pca$eid,proj_PC1=pm[,1],proj_PC2=pm[,2])
for(j in seq_len(npc))gd[,(pc_names[j]):=pm[,j]]
for(j in seq_along(pops))gd[,(paste0('distance.',pops[j])):=sqrt(rowSums(sweep(pm,2,cm[j,],'-')^2))]
gd[,nearest_distance:=do.call(pmin,.SD),.SDcols=paste0('distance.',pops)]
d<-merge(d,gd,by='eid');rm(gd,pca,pm);invisible(gc(verbose=FALSE))
models <- list(COJO='pt.TARGET',`PRS-CSx-auto-meta`='csx.auto',`PRS-CSx`=base_scores,`DiscoDivas-untuned`='disco',`PRS-CSx-fixed-meta`='csx.meta')
if(tune){models[['DiscoDivas-tuned']]<-'disco.cv';available<-c(available,'disco.cv');d[,disco.cv:=0]}
for(s in base_scores)models[[s]]<-s
if(length(grid_scores)){
 g4<-paste0('GRID_',c('AFR','EAS','EUR','SAS'))
 if(all(g4%in%grid_scores)){models[['GRID-tuned']]<-g4;main_methods<-c(main_methods,'GRID-tuned')}
 for(g in intersect(c('GRID_shared','GRID_posterior','GRID_matched'),grid_scores))models[[g]]<-g
}
model_map<-rbindlist(lapply(names(models),function(m)data.table(method=m,scores=paste(models[[m]],collapse=','))))
manifest <- data.table(field=c('trait','type',names(files),'covariates','outcome','folds','seed','bootstrap','distance_PCs','missing_scores','disco_tuned','grid_file','uncertainty'),
 value=c(Y,type,files,paste(covars,collapse=','),outcome_definition,nfold,seed,nboot,npc,paste(missing,collapse=','),tune,arg('grid-file','not supplied'),'Paired subject bootstrap conditional on fixed OOF fits; no discovery/fit uncertainty'))
if(length(missing))cat('Unavailable scores: ',paste(missing,collapse=', '),'\n',sep='')
cat('Ancestry-matched candidates:\n');print(d[,.N,by=target])
if(isTRUE(arg('check',FALSE))) {
 print(manifest);print(model_map);print(rbindlist(audit))
 cat('CHECK complete. Missing PT can be generated by a normal run.\n');quit(status=0)
}
quote_name<-function(x)paste0('`',x,'`')
formula_for <- function(cv,ns=0L,surv=FALSE) {
 terms<-c(quote_name(cv),if(ns)paste0('z',seq_len(ns)))
 as.formula(paste(if(surv)'Surv(time,outcome)' else 'outcome','~',if(length(terms))paste(terms,collapse='+') else '1'))
}
fit_model <- function(f,x,kind) {
 m<-switch(kind,ct=lm(f,data=x),dt=glm(f,data=x,family=binomial()),t2e=coxph(f,data=x,ties='efron'))
 if(kind=='dt'&&!m$converged)stop('Logistic model did not converge')
 if(any(!is.finite(coef(m))))stop('Singular/nonfinite fit; check covariates and score collinearity')
 m
}
predict_model<-function(m,x,kind)as.numeric(if(kind=='ct')predict(m,x) else predict(m,x,type=if(kind=='dt')'response' else 'lp'))
liability <- function(r2,K,P) {
 t<-qnorm(1-K);z<-dnorm(t);i<-z/K;C<-K^2*(1-K)^2/(z^2*P*(1-P));a<-i*(P-K)/(1-K)
 den<-1+C*a*(a-t)*r2
 ifelse(is.finite(den)&den>0,C*r2/den,NA_real_)
}
cindex <- function(y,time,p,fold) {
 num<-den<-0
 for(k in unique(fold)) {
  take<-which(fold==k);z<-data.frame(y=y[take],time=time[take],p=p[take])
  if(nrow(z)<2||!any(z$y==1 & z$time<max(z$time)))next
  cc<-survival::concordancefit(Surv(z$time,z$y),z$p,reverse=TRUE,std.err=FALSE);np<-sum(cc$count[1:3])
  if(np>0&&is.finite(cc$concordance)){num<-num+np*cc$concordance;den<-den+np}
 }
 if(den>0)num/den else NA_real_
}
# All methods share bootstrap indices. Negative held-out R2 is retained.
summarize_predictions <- function(x,pred,base,linear,linear_base,Kpop=NA_real_,B=nboot) {
 nm<-colnames(pred);n<-nrow(x);P<-mean(x$outcome)
 if(type!='t2e') {
  err<-sweep(linear,1,x$outcome,'-')^2;err0<-(linear_base-x$outcome)^2
  stat<-function(ix) {counts<-tabulate(ix,nbins=n);v<-1-as.numeric(crossprod(counts,err))/sum(counts*err0);names(v)<-nm;v}
 } else stat<-function(ix)vapply(nm,function(m)cindex(x$outcome[ix],x$time[ix],pred[ix,m],x$fold[ix]),numeric(1))
 point<-stat(seq_len(n));boots<-matrix(NA_real_,B,length(nm),dimnames=list(NULL,nm))
 for(b in seq_len(B)) {
  ix<-if(type=='ct')sample.int(n,n,replace=TRUE) else unlist(lapply(split(seq_len(n),x$outcome),function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)
  boots[b,]<-stat(ix)
  if(B>=100L&&b%%max(1L,B%/%4L)==0L){cat('  bootstrap ',b,'/',B,'\n',sep='');flush.console()}
 }
 ci<-function(v)if(sum(is.finite(v))>=max(10,.8*B))quantile(v,c(.025,.975),na.rm=TRUE,names=FALSE) else c(NA_real_,NA_real_)
 intervals<-if(B>0)t(apply(boots,2,ci)) else matrix(NA_real_,length(nm),2);res<-data.table(method=nm,estimate=as.numeric(point),lower95=intervals[,1],upper95=intervals[,2],N=n,events=if(type=='ct')NA_integer_ else sum(x$outcome),K=Kpop,P=if(type=='ct')NA_real_ else P)
 res[,metric:=switch(type,ct='OOF_partial_R2',dt='OOF_observed_partial_R2',t2e='OOF_Harrell_C')]
 list(performance=res,bootstrap=boots)
}
set.seed(seed)
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'yeval_disco.R'))
groups<-c(pops,setdiff(sort(unique(na.omit(d$target))),pops))
d[,`:=`(fold=0L,eligible=FALSE,row_index=.I)]
for(g in groups) {
 ix<-which(d$target==g);mm<-lapply(models,function(z)sub('TARGET',g,z,fixed=TRUE))
 mm<-mm[vapply(mm,function(z)all(z%in%available),logical(1))]
 req<-unique(c('outcome',if(type=='t2e')'time',covars,unlist(mm)))
 valid<-complete.cases(d[ix,..req]);for(v in req)if(is.numeric(d[[v]]))valid<-valid&is.finite(d[[v]][ix])
 if(type=='t2e')valid<-valid&d$time[ix]>0
 good<-ix[valid];d$eligible[good]<-TRUE
 strata<-if(type=='ct')list(good) else split(good,d$outcome[good])
 for(j in strata)if(length(j))d$fold[j]<-sample(rep(seq_len(nfold),length.out=length(j)))
}
disco_folds<-if(tune)build_disco_folds(d,d$eligible,pc_names,quality,min_anchor) else NULL
all_predictions<-list()
perf<-differences<-distance_results<-skips<-list()
groups<-c(pops,setdiff(sort(unique(na.omit(d$target))),pops))
for(g in groups) {
 x<-copy(d[target==g]);n_start<-nrow(x)
 mm<-lapply(models,function(s)sub('TARGET',g,s,fixed=TRUE))
 ok_models<-vapply(mm,function(s)all(s%in%available),logical(1));mm<-mm[ok_models]
 for(m in names(models)[!ok_models])skips[[length(skips)+1L]]<-data.table(target=g,method=m,reason='Missing score input or no ancestry-matched PT')
 used<-unique(unlist(mm));required<-c('outcome',if(type=='t2e')'time',covars,used)
 valid<-complete.cases(x[,..required]);for(v in required)if(is.numeric(x[[v]]))valid<-valid&is.finite(x[[v]])
 if(type=='t2e')valid<-valid&x$time>0
 x<-droplevels(x[valid]);n<-nrow(x)
 audit[[length(audit)+1L]]<-data.table(stage=c('matched_before_complete_cases','paired_complete_cases'),target=g,N=c(n_start,n))
 if(n<minn||!length(mm)||(type!='ct'&&min(table(factor(x$outcome,levels=0:1)))<nfold*2)) {
  skips[[length(skips)+1L]]<-data.table(target=g,method='ALL',reason=paste('Insufficient samples or cases/controls:',n));next
 }
 cv<-covars[vapply(x[,..covars],function(v)uniqueN(v)>1,logical(1))]
 if(any(x$fold==0L))stop('Internal fold assignment mismatch')
 cat('START ',g,': N=',n,'; folds=',nfold,'\n',sep='');flush.console()
 pred<-lin<-matrix(NA_real_,n,length(mm),dimnames=list(NULL,names(mm)))
 base<-lb<-rep(NA_real_,n)
 for(k in seq_len(nfold)) {
  tr<-copy(x[fold!=k]);te<-copy(x[fold==k]);it<-which(x$fold==k)
  if(tune){tr[,disco.cv:=disco_folds[[k]][row_index]];te[,disco.cv:=disco_folds[[k]][row_index]]}
  f0<-formula_for(cv,surv=type=='t2e');bm<-fit_model(f0,tr,type);base[it]<-predict_model(bm,te,type)
  lb[it]<-if(type=='dt')predict_model(fit_model(formula_for(cv),tr,'ct'),te,'ct') else base[it]
  for(m in names(mm)) {
   sc<-mm[[m]];mu<-vapply(sc,function(s)mean(tr[[s]]),numeric(1));ss<-vapply(sc,function(s)sd(tr[[s]]),numeric(1))
   if(any(!is.finite(ss)|ss<=0))stop('Constant training score: ',g,' ',m)
   for(j in seq_along(sc)){tr[,(paste0('z',j)):=(get(sc[j])-mu[j])/ss[j]];te[,(paste0('z',j)):=(get(sc[j])-mu[j])/ss[j]]}
   f<-formula_for(cv,length(sc),type=='t2e');fm<-fit_model(f,tr,type);pred[it,m]<-predict_model(fm,te,type)
   lmfit<-if(type=='dt')fit_model(formula_for(cv,length(sc)),tr,'ct') else fm
   lin[it,m]<-if(type=='dt')predict_model(lmfit,te,'ct') else pred[it,m]
  }
  cat('  fold ',k,'/',nfold,' DONE\n',sep='');flush.console()
 }
 if(toupper(arg('write-predictions','FALSE'))=='TRUE')all_predictions[[g]]<-cbind(x[,.(eid,target,fold,outcome,time)],data.table(baseline=base),as.data.table(pred))
 if(any(!is.finite(pred))||any(!is.finite(lin)))stop('Nonfinite held-out prediction: ',g)
 kp<-if(type=='dt')prevalence[target==g]$K else NA_real_
 if(type=='dt'&&(length(kp)!=1||!is.finite(kp)||kp<=0||kp>=1)){skips[[length(skips)+1L]]<-data.table(target=g,method='ALL',reason='Missing valid K');next}
 sm<-summarize_predictions(x,pred,base,lin,lb,kp);pp<-sm$performance;pp[,target:=g]
 # Additional metrics retain the same cohort and fold assignments.
 for(m in names(mm)) {
  if(type=='ct') {
   den<-sum((x$outcome-mean(x$outcome))^2)
   pp[method==m,`:=`(baseline_R2=1-sum((x$outcome-base)^2)/den,full_R2=1-sum((x$outcome-pred[,m])^2)/den,RMSE=sqrt(mean((x$outcome-pred[,m])^2)))]
  } else if(type=='dt') {
   auc<-function(p)as.numeric(pROC::auc(pROC::roc(x$outcome,p,levels=0:1,direction='<',quiet=TRUE)))
   pp[method==m,`:=`(AUC=auc(pred[,m]),baseline_AUC=auc(base),Brier=mean((x$outcome-pred[,m])^2),baseline_Brier=mean((x$outcome-base)^2))]
  }
 }
 perf[[length(perf)+1L]]<-pp
 disco_method<-if(tune)'DiscoDivas-tuned' else 'DiscoDivas-untuned'
 if(nboot>0&&all(c('PRS-CSx',disco_method)%in%names(mm))) {
  bs<-sm$bootstrap;dd<-bs[,disco_method]-bs[,'PRS-CSx'];ci<-quantile(dd,c(.025,.975),na.rm=TRUE)
  reference<-pp[method=='PRS-CSx',estimate];val<-pp[method==disco_method,estimate]-reference
  rel<-if(reference>0&&mean(bs[,'PRS-CSx']>0,na.rm=TRUE)>=.975)100*dd/bs[,'PRS-CSx'] else rep(NA_real_,nboot)
  ri<-if(any(is.finite(rel)))quantile(rel,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
  differences[[length(differences)+1L]]<-data.table(target=g,metric=pp$metric[1],difference=val,lower95=ci[1],upper95=ci[2],relative_percent=if(reference>0)100*val/reference else NA_real_,relative_lower95=ri[1],relative_upper95=ri[2],N=n)
 }
 # Fixed predictions only: bins never choose scores, fit models or use outcomes.
 for(axis in c('distance.EUR','nearest_distance')) {
  breaks<-unique(quantile(x[[axis]],seq(0,1,length.out=nbins+1L),na.rm=TRUE))
  if(length(breaks)<3)next
  bins<-cut(x[[axis]],breaks,include.lowest=TRUE,labels=FALSE)
  selected<-intersect(c('PRS-CSx-auto-meta','PRS-CSx','DiscoDivas-tuned','DiscoDivas-untuned','PRS-CSx-fixed-meta'),names(mm))
  if(!length(selected))next
  for(bin in sort(unique(bins))) {
   ix<-which(bins==bin)
   if(length(ix)<minn||(type!='ct'&&min(table(factor(x$outcome[ix],levels=0:1)))<5))next
   zz<-summarize_predictions(x[ix],pred[ix,selected,drop=FALSE],base[ix],lin[ix,selected,drop=FALSE],lb[ix],kp,B=min(nboot,100L))$performance
   zz[,`:=`(target=g,distance_axis=axis,bin=bin,distance=median(x[[axis]][ix]),distance_min=min(x[[axis]][ix]),distance_max=max(x[[axis]][ix]))]
   distance_results[[length(distance_results)+1L]]<-zz
  }
 }
 cat('Evaluated ',g,': N=',n,'; methods=',paste(names(mm),collapse=', '),'\n',sep='');flush.console()
}
if(!length(perf))stop('No evaluable ancestry groups')
if(length(all_predictions))write_tsv(rbindlist(all_predictions,fill=TRUE),'predictions.tsv.gz')
performance<-rbindlist(perf,fill=TRUE);write_tsv(performance,'performance.tsv')
write_tsv(rbindlist(audit,fill=TRUE),'cohort.tsv')
comparison<-rbindlist(differences)
distperf<-rbindlist(distance_results)
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'yeval_plots.R'))
cat('DONE: ',file.path(out,'report.html'),'\n',sep='')
