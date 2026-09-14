# Serve only validated local references; Shiny's static handler ignores Range.
bplot_local_igv_references <- function(root) {
  result<-list(references=list(),files=list(),indexes=list())
  for(genome in c("hg19","hg38")) {
    name<-if(genome=="hg19")"GRCH37.fasta"else "GRCH38.fasta"
    path<-file.path(root,name)
    if(!file.exists(path) || !file.exists(paste0(path,".fai")))next
    idx<-tryCatch(read.delim(paste0(path,".fai"),header=FALSE,colClasses="character"),error=function(e)data.frame())
    expected<-if(genome=="hg19")249250621 else 248956422
    if(ncol(idx)!=5L || !any(sub("^chr","",idx[[1]])=="1" & idx[[2]]==expected))next
    # Display aliases only: offsets and line widths still address the original
    # FASTA bytes. UCSC gene annotations and cytobands use chr-prefixed names.
    chr<-sub("^chr","",idx[[1]])
    canonical<-chr %in% c(as.character(1:22),"X","Y","M","MT")
    chr[chr=="MT"]<-"M"
    idx[[1]][canonical]<-paste0("chr",chr[canonical])
    # The shared GWAS viewport covers 1-22/X, so IGV uses the same chromosomes.
    idx<-idx[chr %in% c(as.character(1:22),"X"),,drop=FALSE]
    fasta_url<-paste0("bplot_reference/",name)
    index_url<-paste0("bplot_reference/",genome,".fai")
    result$references[[genome]]<-list(fastaURL=fasta_url,indexURL=index_url)
    result$files[[paste0("/",fasta_url)]]<-path
    result$indexes[[paste0("/",index_url)]]<-charToRaw(paste0(paste(apply(idx,1,paste,collapse="\t"),collapse="\n"),"\n"))
  }
  result
}

bplot_igv_reference_response <- function(req,resources) {
  route<-req$PATH_INFO
  if(!route %in% c(names(resources$files),names(resources$indexes)))return(NULL)
  if(!req$REQUEST_METHOD %in% c("GET","HEAD"))return(shiny::httpResponse(405,headers=list(Allow="GET, HEAD")))
  index<-resources$indexes[[route]]
  if(!is.null(index))return(shiny::httpResponse(200,"text/plain",if(req$REQUEST_METHOD=="HEAD")raw()else index,
    headers=list("Content-Length"=as.character(length(index)))))
  path<-resources$files[[route]];size<-file.info(path)$size
  headers<-list("Accept-Ranges"="bytes")
  if(req$REQUEST_METHOD=="HEAD") {
    headers$`Content-Length`<-format(size,scientific=FALSE,trim=TRUE)
    return(shiny::httpResponse(200,"application/octet-stream",raw(),headers))
  }
  range<-req$HTTP_RANGE
  match<-if(is.null(range))character()else regmatches(range,regexec("^bytes=([0-9]+)-([0-9]*)$",range))[[1]]
  if(length(match)==3L) {
    start<-as.numeric(match[2]);end<-if(nzchar(match[3]))min(as.numeric(match[3]),size-1)else min(size-1,start+16*1024^2-1)
    if(is.finite(start) && is.finite(end) && start>=0 && end>=start && start<size && end-start<16*1024^2) {
      con<-file(path,"rb");on.exit(close(con));seek(con,start,origin="start")
      data<-readBin(con,"raw",n=as.integer(end-start+1))
      headers$`Content-Range`<-sprintf("bytes %.0f-%.0f/%.0f",start,end,size)
      return(shiny::httpResponse(206,"application/octet-stream",data,headers))
    }
  }
  headers$`Content-Range`<-sprintf("bytes */%.0f",size)
  shiny::httpResponse(416,"text/plain","A valid byte range of at most 16 MiB is required.",headers)
}
