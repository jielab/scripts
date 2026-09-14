# CIGMA is an orthogonal expression-variance annotation, never a PP(H4).
le8_c3_cigma <- function(layer,outdir) {
  rd<-le8_job_dir(outdir,"c3_coloc");dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  manifest<-Sys.getenv("C3_CIGMA_MANIFEST",unset="")
  result<-Sys.getenv("C3_CIGMA_RESULTS",unset="")
  status<-tibble(status="unavailable",detail="Provide C3_CIGMA_MANIFEST with donor pseudobulk/noise, proportions and kinship, or C3_CIGMA_RESULTS from that runner; omic PGS and GWAS alone are insufficient")
  if(nzchar(manifest)&&nzchar(result))stop("Configure CIGMA manifest OR results, not both")
  if(nzchar(manifest)) {
    if(!file.exists(manifest))stop("Missing CIGMA manifest: ",manifest)
    od<-file.path(rd,"cigma_native");dir.create(od,recursive=TRUE,showWarnings=FALSE)
    installed<-file.path(path.expand("~"),".local/share/le8/cigma-1.1.0/bin/python")
    py<-Sys.getenv("C3_CIGMA_PYTHON",unset=if(file.exists(installed))installed else "python3")
    args<-c(file.path(Sys.getenv("LE8_FDIR"),"c3_cigma.py"),"--manifest",normalizePath(manifest),"--outdir",od)
    rc<-system2(py,vapply(args,shQuote,character(1)),stdout=file.path(od,"runner.log"),stderr=file.path(od,"runner.stderr.log"))
    status<-tibble(status=if(rc==0)"native_run_ok"else"native_run_failed",detail=paste("exit",rc,"; see",od))
    result<-file.path(od,"cigma.results.csv")
    if(rc!=0)result<-"" # keep partial outputs local; do not promote a failed run
  }
  annotation<-tibble()
  if(nzchar(result)) {
    if(!file.exists(result))stop("Missing CIGMA results: ",result)
    z<-as_tibble(data.table::fread(result));need<-c("gene","tissue","specific_p","specific_FDR_manifest","specificity")
    if(!all(need%in%names(z)))stop("CIGMA results missing fields: ",paste(setdiff(need,names(z)),collapse=","))
    if(anyDuplicated(paste(z$gene,z$tissue)))stop("Duplicate CIGMA gene/tissue rows")
    if(nrow(z)) {
      features<-if(layer=="protein"){
        cf<-file.path(rd,"c3.coloc_summary.csv")
        if(file.exists(cf))unique(data.table::fread(cf)$feature)else unique(z$gene)
      }else character()
      map<-tibble(feature=features,gene=unname(le8_assay_genes(features)))
      annotation<-inner_join(map,z,by="gene")|>mutate(source_file=normalizePath(result),
        interpretation="Cell-specific expression regulation annotation; no disease colocalization or protein mediation is established")
      status<-bind_rows(status|>filter(status!="unavailable"),tibble(status="annotated",detail=paste(nrow(annotation),"assay/tissue rows;",nrow(z),"external gene/tissue rows")))
    }
  }
  write_raw_csv(annotation,"c3.CIGMA_annotation.csv",rd)
  write_raw_csv(status,"c3.CIGMA_status.csv",rd)
  openxlsx::write.xlsx(list(status=status,annotation=annotation),file.path(rd,"c3.CIGMA.xlsx"),overwrite=TRUE)
  list(annotation=annotation,status=status)
}
