# Minimal, read-only adapter for the existing LE8 RDS layout. No source() side effects.
a <- commandArgs(trailingOnly=TRUE)
stopifnot(length(a)==8)
phe_file<-a[1]; omic_file<-a[2]; out<-a[3]; cols<-strsplit(a[4],",",fixed=TRUE)[[1]]
max_n<-as.integer(a[5]); seed<-as.integer(a[6]); id<-a[7]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
p<-readRDS(phe_file)
if (!is.data.frame(p) || !id %in% names(p)) stop("Phenotype must be a data.frame with ID")
missing<-setdiff(cols,names(p)); if(length(missing))stop("Missing phenotype columns: ",paste(missing,collapse=", "))
p<-as.data.frame(p)[,unique(c(id,cols)),drop=FALSE]
p[[id]]<-as.character(p[[id]])
if(anyNA(p[[id]])||anyDuplicated(p[[id]]))stop("Invalid/duplicate phenotype IDs")
x<-if(grepl("\\.rds$",omic_file))as.data.frame(readRDS(omic_file)) else {
  if(!requireNamespace("data.table",quietly=TRUE))stop("R package data.table needed for text omics")
  as.data.frame(data.table::fread(omic_file,check.names=FALSE))
}
if(a[8]=="prot")names(x)[names(x)!=id]<-toupper(names(x)[names(x)!=id])
if(!id %in% names(x))stop("Omics table needs explicit ID column")
x[[id]]<-as.character(x[[id]])
if(anyNA(x[[id]])||anyDuplicated(x[[id]]))stop("Invalid/duplicate omics IDs: select one baseline visit first")
x<-x[x[[id]] %in% p[[id]],,drop=FALSE]
if(max_n>0 && nrow(x)>max_n){set.seed(seed);x<-x[sort(sample.int(nrow(x),max_n)),,drop=FALSE]}
p<-p[match(x[[id]],p[[id]]),,drop=FALSE]
features<-setdiff(names(x),id)
if(anyDuplicated(features))stop("Duplicate feature names")
if(!all(vapply(x[,features,drop=FALSE],is.numeric,logical(1))))stop("All omics columns except ID must be numeric; remove metadata first")
# Dates must be ISO text, not ambiguous numeric days.
for(nm in names(p))if(inherits(p[[nm]],"Date")||inherits(p[[nm]],"POSIXt"))p[[nm]]<-format(p[[nm]],"%Y-%m-%d")
write.csv(p,file.path(out,"phenotype.csv"),row.names=FALSE,na="")
writeLines(features,file.path(out,"features.txt"))
writeLines(as.character(c(nrow(x),length(features))),file.path(out,"shape.txt"))
con<-file(file.path(out,"omics.f32"),"wb")
# Bounded temporary allocation, row-major little endian float32 for numpy memmap.
for(start in seq.int(1,nrow(x),by=1000L)){
  ix<-start:min(nrow(x),start+999L)
  writeBin(as.numeric(t(as.matrix(x[ix,features,drop=FALSE]))),con,size=4,endian="little")
}
close(con)
cat("Exported",nrow(x),"people x",length(features),"features\n")
