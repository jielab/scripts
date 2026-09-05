suppressPackageStartupMessages({library(data.table)})
a=commandArgs(trailingOnly=TRUE)
getarg=function(x,d=NULL){i=match(x,a);if(is.na(i))d else a[i+1]}
csx=getarg('--csx'); out=getarg('--out'); pca_file=getarg('--pca',''); centers_file=getarg('--centers',''); ancestry_file=getarg('--ancestry',''); npc=as.integer(getarg('--n-pcs','10'))
if(is.null(csx)||is.null(out))stop('--csx and --out required')
read_any=function(f){
 if(!nzchar(f)||!file.exists(f))return(NULL)
 if(grepl('\\.rds$',f,ignore.case=TRUE)){x=readRDS(f);return(as.data.table(x))}
 fread(f)
}
norm_id=function(d){z=intersect(c('eid','IID','#IID','ID_2','id'),names(d));if(!length(z))stop('No ID column in ',deparse(substitute(d)));setnames(d,z[1],'eid');d[,eid:=as.character(eid)];d}
y=norm_id(fread(csx)); pops=c('AFR','EAS','EUR','SAS'); sc=paste0('CSX_',pops)
if(!all(sc%in%names(y)))stop('Missing CSX columns: ',paste(setdiff(sc,names(y)),collapse=','))
P=NULL; method=''
pca=read_any(pca_file); cen=read_any(centers_file)
if(!is.null(pca)&&!is.null(cen)){
 pca=norm_id(pca); pc=intersect(paste0('PC',seq_len(npc)),names(pca)); if(length(pc)<2)pc=grep('^PC[0-9]+$',names(pca),value=TRUE)[seq_len(min(npc,length(grep('^PC[0-9]+$',names(pca),value=TRUE))))]
 popc=intersect(c('pop','POP','super_pop','ancestry','population'),names(cen)); if(length(pc)>=2&&length(popc)){
   setnames(cen,popc[1],'pop'); cen[,pop:=toupper(as.character(pop))]; cen=cen[pop%in%pops]
   yy=merge(y,pca[,c('eid',pc),with=FALSE],by='eid',all.x=TRUE); q=matrix(NA_real_,nrow(yy),length(pops),dimnames=list(NULL,paste0('q_',pops)))
   for(k in seq_along(pops)){cc=cen[pop==pops[k]][1]; if(nrow(cc))q[,k]=sqrt(rowSums((as.matrix(yy[,..pc])-matrix(as.numeric(cc[,..pc]),nrow(yy),length(pc),byrow=TRUE))^2,na.rm=TRUE))}
   q=1/(q^2+1e-8); q[!is.finite(q)]=0; den=rowSums(q); q=sweep(q,1,pmax(den,1e-300),'/'); q[den<=0,]=1/length(pops); P=as.data.table(q); y=yy; method='PCA_inverse_squared_distance'
 }
}
if(is.null(P)){
 anc=read_any(ancestry_file); if(is.null(anc))stop('No usable PCA/centers and no ancestry posterior file')
 anc=norm_id(anc); y=merge(y,anc,by='eid',all.x=TRUE)
 qn=paste0('posterior_',pops); if(all(qn%in%names(y)))P=as.data.table(as.matrix(y[,..qn])) else {
   ac=intersect(c('ancestry','genetic_ancestry','predicted_ancestry'),names(y)); if(!length(ac))stop('No posterior or ancestry column')
   P=as.data.table(sapply(pops,function(p)as.numeric(toupper(y[[ac[1]]])==p)))
 }
 setnames(P,paste0('q_',pops)); method='ancestry_posterior_fallback'
}
for(j in seq_along(sc)){
 v=as.numeric(y[[sc[j]]]); pcs=grep('^PC[0-9]+$',names(y),value=TRUE); pcs=head(pcs,npc)
 if(length(pcs)>=2){fit=lm(v~.,data=as.data.frame(y[,..pcs]),na.action=na.exclude);v=residuals(fit)}
 s=sd(v,na.rm=TRUE); y[[paste0(sc[j],'_z')]]=if(is.finite(s)&&s>0)(v-mean(v,na.rm=TRUE))/s else v
}
Q=as.matrix(P); Q[!is.finite(Q)]=0; den=rowSums(Q); Q=sweep(Q,1,pmax(den,1e-300),'/'); Q[den<=0,]=1/length(pops); Z=as.matrix(y[,paste0(sc,'_z'),with=FALSE])
outd=data.table(eid=y$eid,DiscoDivas_zero=rowSums(Z*Q)); outd=cbind(outd,as.data.table(Q)); setnames(outd,3:6,paste0('q_',pops)); outd[,disco_weight_source:=method]
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE); fwrite(outd,out,sep='\t')
cat('rows=',nrow(outd),' method=',method,'\n',sep='')
