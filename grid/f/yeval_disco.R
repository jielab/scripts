# Phenotype-tuned DiscoDivas within each OUTER training fold.
# The four input anchor models each combine all four raw PRS-CSx scores.
disco_weights <- function(pc,centers,quality=rep(1,4)) {
 G<-as.matrix(dist(centers))
 if(!is.finite(kappa(G))||kappa(G)>1e12)stop('Singular Disco anchor geometry')
 correction<-as.numeric(solve(G,rep(1,nrow(G))))*quality
 dist<-vapply(seq_len(nrow(centers)),function(j)sqrt(rowSums(sweep(pc,2,centers[j,],'-')^2)),numeric(nrow(pc)))
 eps<-max(G,1)*1e-12
 w<-sweep(1/pmax(dist,eps),2,correction,'*');den<-rowSums(w)
 if(any(abs(den)<=1e-12*pmax(rowSums(abs(w)),1e-300)))stop('Undefined Disco interpolation denominator')
 w<-w/den
 exact<-which(rowSums(dist<=eps)>0)
 for(i in exact){hit<-which(dist[i,]<=eps & correction!=0);if(length(hit)==1){w[i,]<-0;w[i,hit]<-1}}
 w
}
build_disco_folds <- function(d,eligible,pc_names,quality,min_anchor=100L) {
 answer<-vector('list',nfold)
 for(k in seq_len(nfold)) {
  training<-which(eligible & d$fold!=k)
  models<-matrix(NA_real_,nrow(d),4,dimnames=list(NULL,pops))
  centers<-matrix(NA_real_,4,length(pc_names),dimnames=list(pops,pc_names))
  anchors<-list()
  for(j in seq_along(pops)) {
   ix<-training[d$target[training]==pops[j]]
   if(length(ix)<min_anchor)stop('Too few Disco anchor training samples: ',pops[j],'; use --disco-tune FALSE for the saved-score diagnostic')
   tr<-copy(d[ix]);sc<-base_scores
   mu<-vapply(tr[,..sc],mean,numeric(1));ss<-vapply(tr[,..sc],sd,numeric(1))
   if(any(!is.finite(ss)|ss<=0))stop('Constant anchor PRS: ',pops[j])
   for(q in seq_along(sc))tr[,(paste0('z',q)):=(get(sc[q])-mu[q])/ss[q]]
   cv<-covars[vapply(tr[,..covars],function(z)uniqueN(z)>1,logical(1))]
   fit<-fit_model(formula_for(cv,4L,type=='t2e'),tr,type)
   beta<-coef(fit)[paste0('z',1:4)]
   models[,j]<-as.numeric(sweep(sweep(as.matrix(d[,..sc]),2,mu,'-'),2,ss,'/')%*%beta)
   centers[j,]<-vapply(tr[,..pc_names],median,numeric(1))
   # Balanced TRAINING reference for ancestry residualization and scaling.
   anchors[[j]]<-ix
  }
  balanced_n<-min(10000L,lengths(anchors))
  harmonize<-unlist(lapply(anchors,function(ix)sample(ix,balanced_n)),use.names=FALSE)
  pc<-as.matrix(d[,..pc_names]);hx<-cbind(1,pc[harmonize,,drop=FALSE]);allx<-cbind(1,pc)
  hf<-lm.fit(hx,models[harmonize,,drop=FALSE])
  if(any(!is.finite(hf$coefficients)))stop('Singular training PC harmonization')
  residual<-models-allx%*%hf$coefficients
  mu<-colMeans(residual[harmonize,,drop=FALSE]);ss<-apply(residual[harmonize,,drop=FALSE],2,sd)
  if(any(!is.finite(ss)|ss<=0))stop('Zero residual PRS variance')
  scaled<-sweep(sweep(residual,2,mu,'-'),2,ss,'/')
  weights<-disco_weights(pc,centers,quality)
  answer[[k]]<-rowSums(scaled*weights)
  cat('Disco tuned: fold ',k,'/',nfold,' DONE; training=',length(training),'\n',sep='');flush.console()
 }
 answer
}
