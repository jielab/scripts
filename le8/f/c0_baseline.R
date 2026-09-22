# Reconstruct time-safe LE8 variables in memory, before ANY phenotype selection.
# Unsuffixed source assays follow ukb/f/phe.R's baseline convention. Custom
# datasets must supply LE8_BASELINE_MAP (canonical,column,unit; see runbook).
LE8_BASELINE_VERSION <- "2026-09-21.baseline-med-source-v3"
.le8_med_cache <- new.env(parent=emptyenv())
le8_attach_baseline_medication <- function(x) {
  # Raw i0 responses retain -7 (none of the above), which phe.R may have
  # converted to NA. Never infer a negative response from that processed NA.
  root<-get0("indir",ifnotfound=Sys.getenv("UKB_PHE","/mnt/d/data/ukb/phe"))
  path<-Sys.getenv("LE8_BASELINE_MED_FILE",file.path(root,"rap/vip.tab.gz"))
  if(!file.exists(path))return(x)
  explicit<-nzchar(Sys.getenv("LE8_BASELINE_MED_FILE",""))
  key<-paste(path,file.info(path)$size,as.numeric(file.info(path)$mtime),sep="|")
  if(!identical(.le8_med_cache$key,key)) {
    hdr<-if(grepl("[.]rds$",path,ignore.case=TRUE))names(readRDS(path))else names(data.table::fread(path,nrows=0,showProgress=FALSE))
    cols<-grep("^(p)?(6153|6177)(_i0(_a[0-9]+)?|[.]0[.][0-9]+|-0[.][0-9]+)$",hdr,value=TRUE)
    ids<-intersect(c("eid","f.eid","IID"),hdr)
    if(!length(cols)||!length(ids)) {
      if(explicit)stop("LE8_BASELINE_MED_FILE needs eid and explicit baseline p6153/p6177 fields")
      return(x)
    }
    med<-if(grepl("[.]rds$",path,ignore.case=TRUE))as.data.frame(readRDS(path))[,c(ids[1],cols),drop=FALSE]else
      as.data.frame(data.table::fread(path,select=c(ids[1],cols),showProgress=FALSE))
    names(med)[1]<-"eid";med$eid<-as.character(med$eid)
    if(anyNA(med$eid)||anyDuplicated(med$eid))stop("Raw baseline medication input has missing/duplicate eid")
    names(med)[-1]<-paste0(".le8_raw_med_i0_",seq_along(cols))
    .le8_med_cache$data<-med;.le8_med_cache$key<-key;.le8_med_cache$source<-paste(path,paste(cols,collapse=","),sep=" : ")
  }
  med<-.le8_med_cache$data;i<-match(as.character(x$eid),med$eid)
  for(v in setdiff(names(med),"eid"))x[[v]]<-med[[v]][i]
  attr(x,"baseline_med_source")<-.le8_med_cache$source;x
}
le8_rebuild_baseline <- function(dat, audit_dir=NULL) {
  x<-le8_attach_baseline_medication(as.data.frame(dat));n<-nrow(x)
  getnum<-function(v)if(v%in%names(x))suppressWarnings(as.numeric(as.character(x[[v]])))else rep(NA_real_,n)
  dt<-function(z)if(inherits(z,"Date"))as.Date(z)else if(is.numeric(z))as.Date(z,origin="1970-01-01")else as.Date(as.character(z))
  if(!"date_attend"%in%names(x))stop("Baseline date_attend is required")
  baseline<-dt(x$date_attend)
  old_names<-intersect(c("drug.lipid","drug.dm","drug.htn","dm.yes","htn.yes","nonhdl.pts","hba1c.pts","bp.pts"),names(x))
  old<-x[,old_names,drop=FALSE]
  mapfile<-Sys.getenv("LE8_BASELINE_MAP","")
  if(nzchar(mapfile)){
    m<-read.csv(mapfile,stringsAsFactors=FALSE)
    if(!all(c("canonical","column","unit")%in%names(m))||anyDuplicated(m$canonical))stop("Invalid baseline map")
    if(any(!m$column%in%names(x)))stop("Baseline mapping refers to missing columns")
    expected<-c(bmi="kg/m2",bb_TC="mmol/L",bb_HDL="mmol/L",bb_HBA1C="mmol/mol",sbp="mmHg",dbp="mmHg")
    if(any(!m$canonical%in%names(expected))||any(m$unit!=expected[m$canonical]))stop("Convert baseline units before mapping")
    for(j in seq_len(nrow(m)))x[[m$canonical[j]]]<-x[[m$column[j]]]
  }
  # Only explicitly tagged first-visit responses, or explicitly supplied columns.
  medcols<-trimws(strsplit(Sys.getenv("LE8_BASELINE_MED_COLUMNS",""),",",fixed=TRUE)[[1]])
  medcols<-medcols[nzchar(medcols)]
  if(!length(medcols))medcols<-grep("^([.]le8_raw_med_i0_|drug[.]big3.*[_.]i0([_.]|$)|(p)?(6153|6177)_i0(_a[0-9]+)?$)",names(x),value=TRUE)
  if(any(!medcols%in%names(x)))stop("Unknown baseline medication column")
  meds<-matrix(NA_integer_,n,3,dimnames=list(NULL,c("drug.lipid","drug.htn","drug.dm")))
  if(length(medcols)) {
    patterns<-do.call(paste,c(lapply(x[medcols],as.character),sep=" "))
    keys<-unique(patterns)
    lookup<-t(vapply(keys,function(s){
      tokens<-suppressWarnings(as.integer(regmatches(s,gregexpr("-?[0-9]+",s))[[1]]))
      known<-any(tokens%in%c(-7L,1:5))&&!any(tokens%in%c(-1L,-3L))&&
        !(any(tokens==-7L,na.rm=TRUE)&&any(tokens%in%1:5))
      if(known)as.integer(1:3%in%tokens)else rep(NA_integer_,3)
    },integer(3)))
    meds[,]<-lookup[match(patterns,keys),,drop=FALSE]
  }
  for(v in colnames(meds))x[[v]]<-meds[,v]
  baseline_diagnosis<-function(suffix){
    cc<-grep(paste0("^fod_(srd|icd10|ref)_",suffix,"$"),names(x),value=TRUE)
    hit<-rep(FALSE,n)
    for(v in cc){d<-dt(x[[v]]);hit<-hit|(!is.na(d)&!is.na(baseline)&d<=baseline)}
    # No recorded diagnosis is not proof of absence; this indicator describes records.
    hit
  }
  x$nonhdl<- (getnum("bb_TC")-getnum("bb_HDL"))*38.67
  x$nonhdl[!is.finite(x$nonhdl)|x$nonhdl<0]<-NA_real_
  x$hba1c_ngsp<-getnum("bb_HBA1C")*.0915+2.15
  x$hba1c_ngsp[getnum("bb_HBA1C")<=0]<-NA_real_
  for(v in c("sbp","dbp")){
    cc<-grep(paste0("^",v,"_.*_i0([_.]|$)"),names(x),value=TRUE)
    if(length(cc)){
      z<-as.matrix(x[,cc,drop=FALSE]);storage.mode(z)<-"double";z[z<=0]<-NA_real_
      x[[v]]<-rowMeans(z,na.rm=TRUE)
    }
    if(!v%in%names(x))x[[v]]<-rep(NA_real_,n)
    x[[v]][!is.finite(x[[v]])|x[[v]]<=0]<-NA_real_
  }
  x$dm.yes<-ifelse(baseline_diagnosis("t2dm")|x$drug.dm==1|x$hba1c_ngsp>=6.5,1,0)
  # Unknown medication + no positive evidence remains unknown.
  x$htn.yes<-ifelse(baseline_diagnosis("cvd_htn")|x$drug.htn==1,1,0)
  cutscore<-function(v,br,sc)as.numeric(as.character(cut(v,br,labels=sc,right=FALSE)))
  x$bmi.pts<-cutscore(getnum("bmi"),c(0,25,30,35,40,Inf),c(100,70,30,15,0))
  z<-cutscore(x$nonhdl,c(0,130,160,190,220,Inf),c(100,60,40,20,0))
  x$nonhdl.pts<-ifelse(z==0,0,pmax(z-20*x$drug.lipid,0))
  h<-x$hba1c_ngsp
  x$hba1c.pts<-ifelse(x$dm.yes==0&h<5.7,100,ifelse(x$dm.yes==0&h<6.5,60,
    ifelse(h>=6.5|x$dm.yes==1,cutscore(h,c(0,7,8,9,10,Inf),c(40,30,20,10,0)),NA_real_)))
  s<-x$sbp;d<-x$dbp
  z<-ifelse(s>=160|d>=100,0,ifelse(s>=140|d>=90,25,ifelse(s>=130|d>=80,50,ifelse(s>=120,75,100))))
  x$bp.pts<-ifelse(z==0,0,pmax(z-20*x$drug.htn,0))
  x$.le8_baseline_version<-LE8_BASELINE_VERSION
  audit<-do.call(rbind,lapply(old_names,function(v){a<-as.character(old[[v]]);b<-as.character(x[[v]])
    data.frame(variable=v,N=n,changed=sum((is.na(a)!=is.na(b))|(!is.na(a)&!is.na(b)&a!=b)),
      old_missing=sum(is.na(a)),new_missing=sum(is.na(b)),version=LE8_BASELINE_VERSION)}))
  if(is.null(audit))audit<-data.frame(variable=character(),N=integer())
  attr(x,"baseline_audit")<-audit
  if(!is.null(audit_dir)){
    dir.create(audit_dir,recursive=TRUE,showWarnings=FALSE)
    write.csv(audit,file.path(audit_dir,"c0.baseline_rebuild_audit.csv"),row.names=FALSE)
    writeLines(c(LE8_BASELINE_VERSION,paste("baseline medication columns:",paste(medcols,collapse=",")),
      paste("raw medication source:",if(is.null(attr(x,"baseline_med_source")))"not available"else attr(x,"baseline_med_source")),
      paste("mapping:",if(nzchar(mapfile))mapfile else "ukb/f/phe.R baseline raw-assay convention"),
      "Diagnosis requires date <= date_attend; unresolved medication status stays NA.",
      "Raw assays and other LE8 components still require source visit provenance review."),file.path(audit_dir,"c0.baseline_provenance.txt"))
  }
  x
}
