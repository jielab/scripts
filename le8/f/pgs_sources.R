# Source-specific score sensitivity. Existing full PGSs cannot be algebraically
# relabelled as cis/non-MHC scores: retained alleles must actually be rescored.
pgs_shared_mhc <- function(build) {
  if(!build%in%c("37","38"))stop("Verified build 37/38 required for MHC annotation")
  get<-function(n){v<-Sys.getenv(n,"");if(!nzchar(v)){x<-get0(n,ifnotfound=NULL);if(!is.null(x))v<-as.character(x)};suppressWarnings(as.numeric(v))}
  a<-get(paste0("MHC_START_b",build));b<-get(paste0("MHC_END_b",build))
  if(length(a)!=1||length(b)!=1||!is.finite(a)||!is.finite(b)||a>=b)
    stop("Shared MHC_START_b",build,"/MHC_END_b",build," unavailable: run through le8.sh (sources PHE_F). No substitute interval is used.")
  c(start=a,end=b)
}
pgs_source_masks <- function(w,build,regions=data.frame()) {
  mhc<-pgs_shared_mhc(build);ch<-sub("^chr","",as.character(w$CHR));pos<-as.numeric(w$POS)
  if(any(!is.finite(pos))||anyNA(ch))stop("Missing source-score variant coordinates")
  in_mhc<-ch=="6"&pos>=mhc[1]&pos<=mhc[2]
  kind<-sub("^.* / ","",w$component)
  masks<-list(full_recomputed=rep(TRUE,nrow(w)),no_MHC=!in_mhc)
  # Protein gene-cis and metabolite lead-local remain explicitly distinct.
  for(k in intersect(c("cis","trans","local","distal"),unique(kind))) {
    masks[[k]]<-kind%in%k
    if(k%in%c("trans","distal"))masks[[paste0(k,"_no_MHC")]]<-kind%in%k&!in_mhc
  }
  if(nrow(regions))for(i in seq_len(nrow(regions))) {
    r<-regions[i,];if(!all(c("name","chr","start","end","build")%in%names(r)))stop("Invalid PGS_EXCLUDE_REGIONS schema")
    if(as.character(r$build)!=build)stop("Exclusion-region build differs from score weights")
    if(!is.finite(r$start)||!is.finite(r$end)||r$start>r$end)stop("Invalid exclusion interval")
    masks[[paste0("without_",make.names(r$name))]]<-!(ch==sub("^chr","",as.character(r$chr))&pos>=r$start&pos<=r$end)
  }
  masks
}
pgs_plink <- function(args,log) {
  exe<-Sys.getenv("PLINK2",Sys.which("plink2"));if(!nzchar(exe))stop("plink2 unavailable")
  rc<-system2(exe,vapply(args,shQuote,character(1)),stdout=log,stderr=log)
  if(rc!=0)stop("PLINK2 failed; see ",log)
}
pgs_pfile_files <- function(prefix)c(paste0(prefix,c(".pgen",".psam")),if(file.exists(paste0(prefix,".pvar")))paste0(prefix,".pvar")else paste0(prefix,".pvar.zst"))
pgs_pfile_args <- function(prefix)c("--pfile",prefix,if(!file.exists(paste0(prefix,".pvar")))"vzs")
pgs_prepare_genotype <- function(w,ids,pfiledir,cache) {
  allfiles<-unlist(lapply(unique(w$CHR),function(ch)pgs_pfile_files(file.path(pfiledir,paste0("chr",ch)))))
  if(any(!file.exists(allfiles)))stop("Source-score genotype files unavailable: ",paste(head(allfiles[!file.exists(allfiles)],3),collapse=";"))
  key<-pgs_hash(list(inputs=pgs_stamp(allfiles),variants=sort(unique(w$SNP)),ids=sort(ids)))
  rd<-file.path(cache,key);dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  keep<-file.path(rd,"keep.txt");data.table::fwrite(data.frame(`#IID`=ids,check.names=FALSE),keep,sep="\t")
  for(ch in unique(w$CHR)) {
    prefix<-file.path(rd,paste0("chr",ch));done<-paste0(prefix,".complete")
    if(file.exists(done)&&all(file.exists(pgs_pfile_files(prefix))))next
    extract<-paste0(prefix,".variants.txt");writeLines(unique(w$SNP[w$CHR==ch]),extract)
    pgs_plink(c(pgs_pfile_args(file.path(pfiledir,paste0("chr",ch))),"--extract",extract,
      "--keep",keep,"--make-pgen","--threads",Sys.getenv("N_CORES","2"),"--out",prefix),paste0(prefix,".console.log"))
    file.create(done)
  }
  rd
}
pgs_score_sources <- function(w,masks,pfiledir,cache) {
  # Exact ID + effect-allele matching; no strand guessing, no duplicate-ID skip.
  if(anyDuplicated(w$SNP)||any(!grepl("^[ACGT]+$",toupper(w$effect_allele)))||any(!is.finite(w$weight)))stop("Invalid/duplicate score alleles or weights")
  scopes<-names(masks);tot<-NULL;used<-list();status<-list()
  for(ch in unique(w$CHR)) {
    ix<-which(as.character(w$CHR)==as.character(ch));a<-w[ix,];prefix<-file.path(cache,paste0("chr",ch));dir.create(cache,recursive=TRUE,showWarnings=FALSE)
    weight<-data.frame(ID=a$SNP,EA=toupper(a$effect_allele))
    for(n in scopes)weight[[n]]<-a$weight*as.numeric(masks[[n]][ix])
    sf<-paste0(prefix,".weights.tsv");data.table::fwrite(weight,sf,sep="\t")
    pgs_plink(c(pgs_pfile_args(file.path(pfiledir,paste0("chr",ch))),"--score",sf,"1","2","header-read","no-mean-imputation",
      "list-variants","cols=+scoresums","--score-col-nums",paste0("3-",ncol(weight)),"--threads",Sys.getenv("N_CORES","2"),"--out",prefix),paste0(prefix,".console.log"))
    vf<-paste0(prefix,".sscore.vars");if(!file.exists(vf))stop("PLINK did not retain scored-variant list")
    matched<-readLines(vf,warn=FALSE);hit<-a$SNP%in%matched
    aa<-a;aa$matched_by_plink<-hit;used[[length(used)+1L]]<-aa
    z<-as.data.frame(data.table::fread(paste0(prefix,".sscore"),showProgress=FALSE));names(z)<-sub("^#","",names(z))
    if(!all(c("IID",paste0(scopes,"_SUM"))%in%names(z)))stop("Unexpected PLINK multi-score output schema")
    z$IID<-as.character(z$IID);if(anyDuplicated(z$IID))stop("Duplicate IID in genotype score output")
    sc<-data.frame(eid=z$IID);for(n in scopes)sc[[n]]<-z[[paste0(n,"_SUM")]]
    # With no-mean-imputation, flag participants with any missing dosage in
    # the union. The primary source comparison uses this exact common set.
    sc$.complete_genotype<-z$ALLELE_CT>=2*sum(hit)
    if(is.null(tot))tot<-sc else {
      i<-match(tot$eid,sc$eid);if(anyNA(i)||nrow(sc)!=nrow(tot))stop("Genotype sample sets differ across chromosomes")
      for(n in scopes)tot[[n]]<-tot[[n]]+sc[[n]][i]
      tot$.complete_genotype<-tot$.complete_genotype&sc$.complete_genotype[i]
    }
  }
  u<-pgs_bind(used)
  for(n in scopes) {
    mask<-masks[[n]][match(u$SNP,w$SNP)]
    status[[n]]<-data.frame(scope=n,requested=sum(mask),matched=sum(mask&u$matched_by_plink),
      match_fraction=if(sum(mask))mean(u$matched_by_plink[mask])else NA_real_)
  }
  list(scores=tot,variants=u,status=pgs_bind(status))
}
pgs_source_analysis <- function(layer,outdir,candidates,base,biom,im,covars,paired,getdata) {
  rd<-file.path(outdir,"c1_correlate");wf<-Sys.getenv(if(layer=="protein")"PGS_PROT_WEIGHTS"else"PGS_MET_WEIGHTS",
    file.path(outdir,"c2_cause","c2.genetic_score_weights.tsv"))
  status<-list();counts<-list();variants<-list();effects<-list();summaries<-list();components<-list();calibration<-list();region_audit<-list();contrasts<-list()
  fail<-function(detail)list(source_status=data.frame(status="unavailable",detail=detail))
  if(toupper(Sys.getenv("PGS_BUILD_SOURCES","AUTO"))%in%c("FALSE","NO","0"))return(fail("Explicitly disabled PGS_BUILD_SOURCES"))
  if(!file.exists(wf))return(fail(paste("COJO weight manifest missing:",wf)))
  w<-as.data.frame(data.table::fread(wf,showProgress=FALSE))
  if(!all(c("feature","SNP","CHR","POS","effect_allele","weight","component")%in%names(w)))return(fail("Invalid COJO weight manifest schema"))
  w$CHR<-sub("^chr","",as.character(w$CHR))
  source_default<-unique(c(intersect(candidates,c("PCSK9","LPA","CCL19","CCL21","CXCL13","L_VLDL_TG.pct","GlycA","Lactate")),
    head(paired$feature[paired$evidence_pattern=="Both supported: opposite"],3)))
  chosen<-intersect(pgs_csv("PGS_SOURCE_FEATURES",paste(source_default,collapse=",")),candidates)
  if(!length(chosen))return(fail("No source-score candidate among measured/PGS matches"))
  private_cache<-file.path(Sys.getenv("PGS_SOURCE_CACHE_DIR",file.path(indir,".le8_pgs_source_cache")),layer)
  dir.create(private_cache,recursive=TRUE,showWarnings=FALSE)
  project<-if(layer=="protein")dir.X else dir.met.gwas
  # The C2 manifest retains QTL coordinates; verify against its GWAS QC build.
  build_for<-function(f) {
    declared<-Sys.getenv("PGS_WEIGHTS_BUILD","")
    qcf<-file.path(project,"common",f,"qc",paste0(f,".grch"))
    qc<-if(file.exists(qcf))trimws(readLines(qcf,n=1,warn=FALSE))else""
    qc<-sub("^(GRCh|b)","",qc,ignore.case=TRUE)
    if(nzchar(declared)&&nzchar(qc)&&declared!=qc)stop("PGS_WEIGHTS_BUILD conflicts with QTL QC metadata")
    b<-if(nzchar(declared))declared else qc
    if(!b%in%c("37","38"))stop("Weight build unverified; supply PGS_WEIGHTS_BUILD or QTL qc/<feature>.grch")
    b
  }
  regions_file<-Sys.getenv("PGS_EXCLUDE_REGIONS","")
  regions<-if(nzchar(regions_file))as.data.frame(data.table::fread(regions_file))else data.frame()
  for(f in chosen) {
    tryCatch({
      a<-w[w$feature==f,,drop=FALSE];if(!nrow(a))stop("No COJO weights for candidate")
      if(any(!a$CHR%in%as.character(1:22)))stop("Autosomal source analysis requires autosomal weights")
      build<-build_for(f)
      annotation_status<-"Metabolite lead-local/distal; not gene cis"
      if(layer=="protein") {
        bed<-get0("prot_bed_file",ifnotfound="")
        bbuild<-Sys.getenv("LE8_PROT_BED_BUILD","")
        if(!nzchar(bbuild)&&grepl("(b|[.])38[.]",basename(bed)))bbuild<-"38"
        if(!nzchar(bbuild)&&grepl("(b|[.])37[.]",basename(bed)))bbuild<-"37"
        a$component<-"COJO PGS / unknown"
        annotation_status<-"cis/trans withheld: annotation build unverified or differs from QTL build"
        if(file.exists(bed)&&identical(build,bbuild)) {
          ann<-as.data.frame(data.table::fread(bed,header=FALSE,showProgress=FALSE))
          hit<-which(as.character(ann[[4]])==f)
          if(length(hit)==1&&all(is.finite(as.numeric(ann[hit,2:3])))) {
            pad<-pgs_num("C2_CIS_WINDOW_BP",1e6)
            cis<-a$CHR==sub("^chr","",as.character(ann[hit,1]))&a$POS>=as.numeric(ann[hit,2])+1-pad&a$POS<=as.numeric(ann[hit,3])+pad
            a$component<-paste0("COJO PGS / ",ifelse(cis,"cis","trans"))
            annotation_status<-paste("Gene annotation verified in build",build)
          }
        }
      }
      local_regions<-regions
      # A leave-hub-out sensitivity uses observed weight-table coordinates,
      # never a hard-coded SH2B3/MHC interval or another build's coordinates.
      for(tag in pgs_csv("PGS_HUB_VARIANTS","rs3184504")) {
        hit<-w[w$SNP==tag,,drop=FALSE]
        if(nrow(hit)) {
          same<-vapply(hit$feature,function(ff)identical(tryCatch(build_for(ff),error=function(e)""),build),logical(1))
          loc<-unique(hit[same,c("CHR","POS"),drop=FALSE])
          if(nrow(loc)==1) {
            pad<-pgs_num("PGS_HUB_WINDOW_BP",1e6)
            r<-data.frame(name=paste0(tag,"_region"),chr=loc$CHR,start=max(1,loc$POS-pad),end=loc$POS+pad,build=build)
            overlap<-a$CHR==as.character(r$chr)&a$POS>=r$start&a$POS<=r$end
            if(any(overlap))local_regions<-pgs_bind(list(local_regions,r))
          }
        }
      }
      if(nrow(local_regions)){rr<-local_regions;rr$feature<-f;region_audit[[f]]<-rr}
      masks<-pgs_source_masks(a,build,local_regions)
      pd<-Sys.getenv("PGS_PFILE_DIR",file.path("/mnt/e/ukbGen",build,"imp"))
      genotypebuild<-Sys.getenv("PGS_PFILE_BUILD",if(identical(pd,file.path("/mnt/e/ukbGen",build,"imp")))build else "")
      if(!identical(genotypebuild,build))stop("Custom PGS_PFILE_DIR requires matching PGS_PFILE_BUILD; liftover is not guessed")
      inputfiles<-unlist(lapply(unique(a$CHR),function(ch)pgs_pfile_files(file.path(pd,paste0("chr",ch)))))
      sig<-pgs_hash(list(version=PGS_FOCUS_VERSION,w=a,masks=masks,build=build,pfiles=pgs_stamp(inputfiles),ids=sort(base$eid)))
      cache<-file.path(private_cache,"scores",sig);dir.create(cache,recursive=TRUE,showWarnings=FALSE)
      file<-file.path(cache,"scores.rds")
      res<-if(file.exists(file)&&!LE8_REPLACE)readRDS(file)else {
        # Cached subset makes scoring multiple source columns cheap. Each
        # feature is isolated so a missing chromosome cannot erase other candidates.
        pp<-pgs_prepare_genotype(a,base$eid,pd,file.path(private_cache,"genotypes"))
        z<-pgs_score_sources(a,masks,pp,cache);saveRDS(z,file);z
      }
      res$status$feature<-f;res$status$build<-build;res$status$MHC_interval<-paste(pgs_shared_mhc(build),collapse="-")
      counts[[f]]<-res$status;res$variants$feature<-f;res$variants$build<-build;variants[[f]]<-res$variants
      dd<-base;dd$.m<-as.numeric(biom[[f]][im]);j<-match(dd$eid,res$scores$eid)
      good<-!is.na(j)&res$scores$.complete_genotype[j]%in%TRUE
      dd<-dd[good,,drop=FALSE];j<-j[good]
      validscopes<-res$status$scope[res$status$matched>0&res$status$match_fraction>=pgs_num("PGS_MIN_VARIANT_MATCH",.9)]
      existing<-getdata(f);dd$.g<-existing$.g[match(dd$eid,existing$eid)]
      common<-pgs_complete(dd,c(".m",".g",covars))
      for(n in validscopes)common<-common&is.finite(res$scores[[n]][j])
      dd<-dd[common,,drop=FALSE];j<-j[common]
      r0<-pgs_pair(dd,covars,f,scope="existing_full_common_genotype")
      summaries[[paste(f,"existing")]]<-r0$summary;effects[[paste(f,"existing")]]<-r0$effects
      fullcor<-if(sum(is.finite(dd$.g))>2)cor(dd$.g,res$scores$full_recomputed[j],use="complete.obs")else NA_real_
      for(n in validscopes) {
        d<-dd;d$.g<-res$scores[[n]][j]
        r<-pgs_pair(d,covars,f,scope=n);effects[[paste(f,n)]]<-r$effects;summaries[[paste(f,n)]]<-r$summary
        ca<-pgs_calibrate(d,covars);ca$audit$feature<-f;ca$audit$scope<-n;calibration[[paste(f,n)]]<-ca$audit
        co<-pgs_component_models(ca$data,covars,f,scope=n);components[[paste(f,n)]]<-co$effects;contrasts[[paste(f,n)]]<-co$contrast
      }
      status[[f]]<-data.frame(feature=f,status=if(length(validscopes))"ok"else"no adequately matched scopes",build=build,
        annotation_status=annotation_status,full_score_correlation=fullcor,N_source_common=nrow(dd),N_genotype_complete=sum(good),N_genotype_incomplete=sum(!good),detail="Existing PGS untouched; common sample across recomputed scopes. Cis/trans labels inherit C2 annotation; check its build.")
    },error=function(e){status[[f]]<<-data.frame(feature=f,status="unavailable",detail=conditionMessage(e))})
  }
  list(source_status=pgs_bind(status),source_regions=pgs_bind(region_audit),source_counts=pgs_bind(counts),source_variants=pgs_bind(variants),
    source_models=pgs_adjust(pgs_bind(effects),groups=c("scope","model","term")),source_paired=pgs_bind(summaries),
    source_components=pgs_adjust(pgs_bind(components),groups=c("scope","model","term")),source_calibration=pgs_bind(calibration),
    source_contrasts=pgs_adjust(pgs_bind(contrasts),groups="scope"))
}
