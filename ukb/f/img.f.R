# Preserve numeric imaging measurements without imputation or missingness filters.
# Metadata is kept separately from measured brain phenotypes.
prepare_ukb_imaging <- function(indir) {
  file<-file.path(indir,'rap/img.tab.gz')
  x<-data.table::fread(file,showProgress=FALSE)
  x[,eid:=as.character(eid)];stopifnot(!anyDuplicated(x$eid))
  extra<-file.path(indir,'rap/img.extra.tab.gz')
  if(file.exists(extra)){
    e<-data.table::fread(extra,showProgress=FALSE);e[,eid:=as.character(eid)]
    missing<-setdiff(names(e),names(x))
    if(length(missing))x<-merge(x,e[,c('eid',missing),with=FALSE],by='eid',all=TRUE,sort=FALSE)
  }
  map<-data.table::fread(file.path(indir,'common/img.lst'),header=FALSE,fill=TRUE)[,1:2]
  data.table::setnames(map,c('field','label'))
  cols<-grep('_i2$',names(x),value=TRUE)
  nums<-cols[vapply(x[,..cols],is.numeric,logical(1))]
  nums<-setdiff(nums,c('p54_i2'))
  img<-x[,c('eid',nums),with=FALSE]
  data.table::setnames(img,nums,paste0('img_',sub('_i2$','',nums)))
  img<-img[rowSums(!is.na(img[,-1]))>0]
  dir.create(file.path(indir,'Rdata'),recursive=TRUE,showWarnings=FALSE)
  saveRDS(as.data.frame(img),file.path(indir,'Rdata/img.raw.rds'),compress=FALSE)
  visit_cols<-intersect(c('p53_i2','p53_i3','p54_i2','p54_i3'),names(x))
  visits<-x[,c('eid',visit_cols),with=FALSE]
  for(cc in intersect(c('p53_i2','p53_i3'),names(visits)))visits[,(cc):=as.Date(get(cc))]
  saveRDS(as.data.frame(visits),file.path(indir,'Rdata/img.visits.rds'))
  map<-map[field%in%sub('^img_','',names(img))];map[,variable:=paste0('img_',field)]
  data.table::fwrite(map,file.path(indir,'Rdata/img.dictionary.csv'))
  message('Unimputed imaging: ',nrow(img),' participants, ',ncol(img)-1,' variables')
  as.data.frame(img)
}
