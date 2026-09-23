# Joint posterior propagation and distance definitions; loaded by Yeval.R.
posterior_pops <- c('AFR','EAS','EUR','SAS')
cov_pairs <- t(combn(1:4,2))
cov_columns <- unlist(lapply(1:4,function(j)paste0('cov.',posterior_pops[j],'.',posterior_pops[j:4])))

load_posterior <- function() {
 mode<-arg('posterior-mode','required')
 if(!mode%in%c('required','off'))stop('--posterior-mode must be required or off')
 if(mode=='off')return(NULL)
 f<-arg('posterior-file',file.path(score_dir,'csx.posterior.tsv.gz'))
 if(!file.exists(f))stop('Missing individual posterior: ',f,'. Run 1csx.sh --posterior TRUE; --posterior-mode off is only for a benchmark without individual posterior panels.')
 meta<-read_table(paste0(f,'.metadata.tsv'))
 md<-setNames(meta$value,meta$field)
 if(md[['schema']]!='grid_csx_moments_v1'||md[['centering']]!='discovery_EAF'||
    !identical(unname(md[['md5']]),unname(tools::md5sum(f))))stop('Invalid/stale posterior metadata or centering')
 if(arg('allow-chromosome-subset','FALSE')!='TRUE' && !setequal(strsplit(md[['chromosomes']],',',fixed=TRUE)[[1]],as.character(1:22)))stop('Posterior lacks all 22 autosomes; use --allow-chromosome-subset TRUE only for deliberate subset analyses')
 z<-ids(read_table(f));req<-c(paste0('csx.',posterior_pops),cov_columns)
 if(!all(req%in%names(z))||any(!is.finite(as.matrix(z[,..req]))))stop('Invalid posterior means/covariance')
 # Every individual covariance must be positive semidefinite. Cholesky-like
 # eigenvalue checks on 4x4 matrices are inexpensive relative to scoring.
 mat<-as.matrix(z[,..cov_columns]);diag_ix<-c(1L,5L,8L,10L)
 if(any(mat[,diag_ix]<0))stop('Negative posterior variance')
 for(j in 1:4)for(k in j:4){
  v<-z[[paste0('cov.',posterior_pops[j],'.',posterior_pops[k])]]
  bound<-sqrt(z[[paste0('cov.',posterior_pops[j],'.',posterior_pops[j])]]*z[[paste0('cov.',posterior_pops[k],'.',posterior_pops[k])]])
  if(any(abs(v)>bound+1e-9*pmax(1,bound)))stop('Posterior covariance violates Cauchy-Schwarz')
 }
 # The producer constructs PSD covariance from synchronized scores. Full PSD
 # checks on a deterministic spread of rows detect malformed imported files.
 for(i in unique(round(seq(1,nrow(z),length.out=min(1000,nrow(z)))))){
  s<-matrix(0,4,4);at<-0L
  for(j in 1:4)for(k in j:4){at<-at+1L;s[j,k]<-s[k,j]<-mat[i,at]}
  if(min(eigen(s,symmetric=TRUE,only.values=TRUE)$values)< -1e-8*max(1,max(diag(s))))stop('Non-PSD individual covariance')
 }
 z
}

training_geometry <- function(gd,pm,pc,reference_centers) {
 path<-arg('training-centers',file.path(score_dir,'training_centers.tsv'))
 distance_source<-arg('distance-source','discovery')
 if(!distance_source%in%c('discovery','reference'))stop('--distance-source must be discovery or reference')
 if(distance_source=='discovery'&&!file.exists(path))stop('Discovery centres missing: ',path,'. Build them with f/training_geometry.py. A reference-only exploratory plot requires explicit --distance-source reference.')
 if(distance_source=='discovery'){
  ct<-read_table(path);need<-c('POP',pc,'N_GWAS','pca_space','source','kind')
  if(!all(need%in%names(ct))||anyDuplicated(ct$POP)||!all(posterior_pops%in%ct$POP))stop('Training centres need unique POP, PC1.., N_GWAS, pca_space, source, kind')
  ct<-ct[match(posterior_pops,POP)]
  if(any(ct$kind!='discovery')||any(!nzchar(ct$source))||anyNA(ct$source))stop('Training centres must explicitly document discovery provenance')
  space<-arg('pca-space','')
  if(!nzchar(space)||any(ct$pca_space!=space))stop('--pca-space must match discovery centres and the actual target projection basis')
  if(any(!is.finite(ct$N_GWAS)|ct$N_GWAS<=0)||any(!is.finite(as.matrix(ct[,..pc]))))stop('Invalid discovery centres/N_GWAS')
  weights<-ct$N_GWAS/sum(ct$N_GWAS)
  status<-'discovery';label<-'RMS distance to GWAS training groups'
  description<-'Sample-size-weighted RMS distance to discovery-group means in the supplied common PC space; a multi-training extension, not the original single-training distance.'
 }else{
  ct<-copy(reference_centers)[match(posterior_pops,POP)]
  weights<-rep(.25,4);status<-'reference_proxy';label<-'RMS distance to four 1KG reference groups'
  description<-'Equal-weight four-reference distance. Discovery centres were not supplied; this axis is a reference proxy and must not be called distance to GWAS training data.'
 }
 distmat<-vapply(seq_len(4),function(j)sqrt(rowSums(sweep(pm,2,as.numeric(ct[j,..pc]),'-')^2)),numeric(nrow(pm)))
 for(j in 1:4)gd[,(paste0('distance.training.',posterior_pops[j])):=distmat[,j]]
 gd[,distance.analysis:=sqrt(as.numeric(distmat^2%*%weights))]
 ct[,mixture_weight:=weights]
 list(gd=gd,centers=ct,status=status,label=label,description=description)
}

read_variance <- function() {
 f<-arg('genetic-variance-file',file.path(score_dir,'genetic_variance.tsv'))
 if(!file.exists(f))return(NULL)
 v<-read_table(f)
 if(!all(c('trait','target','scale','source')%in%names(v)))stop('Genetic variance table needs trait,target,scale,source, and h2 or genetic_variance')
 v<-v[trait==Y]
 if(anyDuplicated(v$target)||anyNA(v$source)||any(!nzchar(v$source)))stop('Duplicate variance target or missing provenance')
 if(!'h2'%in%names(v))v[,h2:=NA_real_]
 if(!'genetic_variance'%in%names(v))v[,genetic_variance:=NA_real_]
 if(any(is.finite(v$h2)&(!is.finite(v$h2)|v$h2<=0|v$h2>=1)))stop('Require 0 < SNP h2 < 1')
 if(any(is.finite(v$genetic_variance)&v$genetic_variance<=0))stop('Genetic variance must be positive')
 if(any(!is.finite(v$h2)&!is.finite(v$genetic_variance)))stop('Every supplied variance row needs a finite h2 or genetic_variance; template NA values must be filled')
 if(any(is.finite(v$h2)&is.finite(v$genetic_variance)))stop('Supply h2 OR genetic_variance per row, not both')
 expected<-switch(type,ct='residual_phenotype',dt='log_odds',t2e='log_hazard')
 if(any(v$scale!=expected))stop('Variance scale must be ',expected,'; liability h2 cannot be substituted for a log-hazard/log-odds variance')
 if(type!='ct'&&any(is.finite(v$h2)))stop('For dt/t2e supply absolute genetic_variance on log-odds/log-hazard scale, not h2')
 v
}

posterior_variance <- function(z,w) {
 ans<-numeric(nrow(z));diagonal<-numeric(nrow(z))
 for(j in 1:4)for(k in j:4){
  term<-w[j]*w[k]*z[[paste0('cov.',posterior_pops[j],'.',posterior_pops[k])]]
  ans<-ans+if(j==k)term else 2*term
  if(j==k)diagonal<-diagonal+term
 }
 if(any(ans< -1e-8*pmax(1,diagonal)))stop('Negative combined variance: non-PSD input covariance')
 list(total=pmax(ans,0),diagonal=diagonal)
}

individual_posterior <- function(te,tr,fit,baseline_fit,ss,g,k) {
 w<-as.numeric(coef(fit)[paste0('z',1:4)])/ss
 vv<-posterior_variance(te,w)
 prior<-NA_real_;source<-'not supplied'
 if(!is.null(variance_spec)){
  row<-variance_spec[target==g]
  if(nrow(row)){
   prior<-row$genetic_variance;source<-row$source
   if(type=='ct'&&is.finite(row$h2))prior<-row$h2*var(residuals(baseline_fit))
  }
 }
 r2<-if(is.finite(prior)&&prior>0)1-vv$total/prior else rep(NA_real_,nrow(te))
 # Negative reliability is diagnostic of scale/prior/calibration problems;
 # preserve it. Never clip or force a decreasing curve.
 out<-te[,c('eid','target','fold','proj_PC1','proj_PC2','distance.analysis',paste0('distance.training.',posterior_pops)),with=FALSE]
 out[,`:=`(posterior_variance=vv$total,posterior_sd=sqrt(vv$total),
           diagonal_only_variance=vv$diagonal,cross_population_covariance=vv$total-vv$diagonal,
           prior_genetic_variance=prior,individual_R2=r2,variance_source=source)]
 for(j in 1:4)out[,(paste0('weight.',posterior_pops[j])):=w[j]]
 out
}
