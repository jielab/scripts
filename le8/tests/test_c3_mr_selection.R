# Run with Rscript tests/test_c3_mr_selection.R
suppressPackageStartupMessages(library(dplyr))
script<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))[1]
fdir<-file.path(dirname(normalizePath(script)),"..","f")
load_function<-function(file,name) {
  for(expr in parse(file))if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&
    identical(expr[[2]],as.name(name)))eval(expr,.GlobalEnv)
}
for(n in c("read_stage_cache","write_stage_cache","cache_valid"))load_function(file.path(fdir,"comm.f.R"),n)
load_function(file.path(fdir,"c0_revision_core.R"),"le8_hash_object")
load_function(file.path(fdir,"c3_revision_locus.R"),"select_nonoverlapping_leads")
source(file.path(fdir,"c3_mr_selection.R"))
LE8_REPLACE<-FALSE;WINDOW_BP<-500000;MIN_SNPS<-50
C3_P12<-c(conservative=1e-6,default=1e-5,liberal=1e-4)
truthy<-function(x)tolower(x)%in%c("true","1","yes")
mr<-tibble(exposure=c("A","A","B","C"),analysis=c("local","distal","distal","unknown"),
  pval=c(.001,.001,.002,.001),FDR_all=c(.01,.1,.02,.01),
  instrument_snps=c("rs1","rs2","rs3;rs4","rs5"),instrument_positions=c("1:100","2:200","3:300;3:2000300","4:400"))
sig<-c3_select_mr(mr,"metabolite")
stopifnot(identical(sig$exposure,c("A","B")),nrow(c3_select_mr(mr[0,],"metabolite"))==0,
  !"rs2"%in%c3_mr_instruments(sig)$SNP,
  inherits(try(c3_select_mr(select(mr,-FDR_all),"metabolite"),silent=TRUE),"try-error"))
root<-tempfile("c3-test-");dir.create(root)
qfile<-file.path(root,"qtl.tsv");writeLines("fixture",qfile)
reads<-0L
read_sumstat_snps<-function(file,snps) {
  reads<<-reads+1L
  tibble(SNP=c("rs1","rs3","rs4"),CHR=c("1","3","3"),POS=c(100,300,2000300),P=c(1e-9,1e-8,1e-7),
    BETA=c(.1,.2,.3),SE=.01) |> filter(SNP%in%snps)
}
iv<-c3_read_mr_qtl(sig,qfile,file.path(root,"qtl_cache"))
stopifnot(reads==1L,nrow(iv)==3L,identical(iv$analysis,c("local","distal","distal")))
stopifnot(identical(c3_read_mr_qtl(sig,qfile,file.path(root,"qtl_cache")),iv),reads==1L)
# A changed source must invalidate an otherwise identical cached extraction.
writeLines("changed fixture",qfile)
invisible(c3_read_mr_qtl(sig,qfile,file.path(root,"qtl_cache")))
stopifnot(reads==2L)
# Reordered candidate list cannot change locus identity; only MR-used SNPs enter.
leads<-select_nonoverlapping_leads(iv |> filter(analysis=="distal"),3L,500000)
stopifnot(nrow(leads)==2L,all(leads$analysis=="distal"))
missing_mr<-sig;missing_mr$instrument_snps[1]<-"absent"
stopifnot(inherits(try(c3_read_mr_qtl(missing_mr,qfile,root),silent=TRUE),"try-error"))
# Old numeric cache indices are irrelevant, but class and geometry must match.
yfile<-file.path(root,"y.tsv");writeLines("outcome",yfile)
oldfile<-file.path(root,"0042_A_chr1_100.rds")
z<-list(summary=tibble(layer="metabolite",feature="A",chr="1",lead_pos=100,start=1,end=500100,
  locus_class="local",case_fraction_source="not required for beta/varbeta cc ABF"),variants=tibble(),regional=tibble())
write_stage_cache(z,oldfile)
Sys.setFileTime(oldfile,Sys.time()+2)
newfile<-file.path(root,"new.rds")
stopifnot(is.null(c3_cached_locus(newfile,oldfile,"A","metabolite","1",100,"distal",qfile,yfile,"cc",NA_real_)))
stopifnot(identical(c3_cached_locus(newfile,oldfile,"A","metabolite","1",100,"local",qfile,yfile,"cc",NA_real_),z))
# New caches still load when the legacy file no longer exists (resumption).
unlink(oldfile)
stopifnot(identical(c3_cached_locus(newfile,character(),"A","metabolite","1",100,"local",qfile,yfile,"cc",NA_real_),z))
unlink(root,recursive=TRUE)
cat("C3 selection, QTL caching, missing-input guard and locus migration tests passed\n")

# Exercise the actual layer driver without running plots or GPU jobs.
# An old aggregate result must not bypass MR selection; a new result can resume.
load_function(file.path(fdir,"c3_coloc.R"),"run_c3_layer")
`%||%`<-function(x,y)if(is.null(x))y else x
sandbox<-tempfile("c3-driver-");dir.create(sandbox)
out.prot<-file.path(sandbox,"trait","prot");out.met<-file.path(sandbox,"trait","met")
dir.create(file.path(out.met,"c2_cause"),recursive=TRUE)
qfile<-file.path(sandbox,"qtl.tsv");yfile<-file.path(sandbox,"y.tsv")
writeLines("QTL",qfile);writeLines("Y",yfile)
saveRDS(list(MR=mr),file.path(out.met,"c2_cause","c2.res.rds"))
dir.create(file.path(out.met,"c3_coloc"))
saveRDS(list(summary=tibble(feature="STALE"),regional=tibble(),variants=tibble()),
  file.path(out.met,"c3_coloc","c3.res.rds"))
LE8_REUSE_RESULTS<-FALSE;LE8_JOB<-"c3_coloc";C3_MR_FDR<-.05;MAX_FEATURES<-200L
MAX_LOCI_PER_FEATURE<-3L;H4_STRONG<-.7;C3_CODE_VERSION<-C3_SELECTION_VERSION
Y<-"trait";dir.met.gwas<-sandbox;dir.X<-sandbox
le8_c3_cigma<-function(...)NULL
le8_load_revision<-le8_finish_revision<-setwd2<-function(...)NULL
le8_job_dir<-function(outdir,job)file.path(outdir,job)
get_y_gwas_file<-function(...)yfile
find_qtl_files<-function(...)list(full=qfile)
infer_outcome_type<-function(...)"cc"
get_case_fraction<-function(...)NA_real_
write_raw_csv<-write_raw_tsv<-write_xlsx2<-finalize_outputs<-function(...)NULL
le8_stage_start<-function(...)Sys.time()
le8_stage_done<-function(...)NULL
le8_stage<-function(label,expr,...)force(expr)
plot_coloc_results<-plot_gpu_coloc_validation<-plot_c3_pgs_integration<-function(...)NULL
read_c3_pgs_integration<-function(...)tibble()
credible_set_audit<-function(...)list(overall=tibble(),by_locus=tibble())
le8_c3_sets<-function(...)list(Causal_Tier1=character(),Causal_Tier2plus=character(),Causal_any=character())
module_meta<-function(layer,extra)extra
calls<-0L;gpu_manifests<-list()
coloc_one_locus<-function(feature,layer,qtl_file,lead_chr,lead_pos,ygwas_file,outcome_type,case_frac,locus_class=NULL) {
 calls<<-calls+1L
 stopifnot(locus_class%in%c("local","distal"))
 list(summary=tibble(feature=feature,locus=paste0("chr",lead_chr,":",lead_pos),status="ok",PP.H4=.8),
   variants=tibble(),regional=tibble())
}
run_gpu_coloc_step<-function(layer,rawdir,manifest,ygfile,outtype,gpudir) {
 gpu_manifests[[length(gpu_manifests)+1L]]<<-manifest
 stopifnot(grepl("/gpu_coloc/runs/",gpudir))
 0L
}
read_gpu_coloc_results<-function(...)list(results=tibble(),status=tibble())
first<-run_c3_layer("metabolite")
stopifnot(calls==3L,nrow(first$summary)==3L,!"STALE"%in%first$summary$feature,
 identical(unique(first$summary$feature),c("A","B")),nrow(gpu_manifests[[1]])==3L,
 !is.null(first$meta$selection_signature))
second<-run_c3_layer("metabolite")
stopifnot(calls==3L,length(gpu_manifests)==1L)
# Tightening MR evidence changes the aggregate/GPU key but reuses A's locus.
mr$FDR_all[mr$exposure=="B"]<-.5
saveRDS(list(MR=mr),file.path(out.met,"c2_cause","c2.res.rds"))
third<-run_c3_layer("metabolite")
stopifnot(calls==3L,nrow(third$summary)==1L,nrow(gpu_manifests[[2]])==1L,
 !identical(first$meta$selection_signature,third$meta$selection_signature))
mr$FDR_all<-1
saveRDS(list(MR=mr),file.path(out.met,"c2_cause","c2.res.rds"))
empty<-run_c3_layer("metabolite")
stopifnot(nrow(empty$summary)==0L,calls==3L,length(gpu_manifests)==2L)
unlink(sandbox,recursive=TRUE)
cat("C3 layer driver: stale aggregate, exact GPU manifest, resume, changed MR and empty selection passed\n")

load_function(file.path(fdir,"comm.f.R"),"read_sumstat_header")
header_file<-tempfile(fileext=".gz")
con<-gzfile(header_file,"wt");writeLines(c("SNP\tCHR\tPOS","rs1\t1\t100"),con);close(con)
before<-nrow(showConnections(all=TRUE))
for(i in seq_len(150))stopifnot(read_sumstat_header(header_file)=="SNP\tCHR\tPOS")
stopifnot(nrow(showConnections(all=TRUE))==before)
unlink(header_file)
cat("Repeated compressed header reads close their connections\n")
