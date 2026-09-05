# LE8 5C revision 2026-09-05. No UKB data are bundled.
# Shared integrity, outcome and reproducibility functions.
LE8_REVISION <- "2026-09-05.5c-audit-v1"
.le8_revision_state <- new.env(parent=emptyenv())
.le8_revision_state$hashes <- new.env(parent=emptyenv())
.le8_revision_state$deps <- character()

le8_num_env <- function(name, default) {
  z <- suppressWarnings(as.numeric(Sys.getenv(name, unset=as.character(default))))
  if(length(z)!=1L || !is.finite(z)) stop("Invalid numeric setting: ",name,call.=FALSE)
  z
}
le8_csv_env <- function(name, default="") {
  z <- trimws(strsplit(Sys.getenv(name,unset=default),",",fixed=TRUE)[[1]])
  unique(z[nzchar(z)])
}
le8_hash_object <- function(x) {
  f <- tempfile(); on.exit(unlink(f),add=TRUE)
  base::saveRDS(x,f,version=2,compress=FALSE)
  unname(tools::md5sum(f))
}
le8_file_hash <- function(path) {
  if(length(path)!=1L || is.na(path) || !nzchar(path) || !file.exists(path)) return("missing")
  path <- normalizePath(path,winslash="/",mustWork=TRUE)
  fi <- file.info(path)
  if(isTRUE(fi$isdir)) return("directory")
  key <- paste(path,fi$size,as.numeric(fi$mtime),as.numeric(fi$ctime),sep="|")
  h <- .le8_revision_state$hashes
  if(exists(key,envir=h,inherits=FALSE)) return(get(key,envir=h))
  # Persist full-content digests; file ctime is part of the invalidation key.
  root <- Sys.getenv("LE8_ANALYSIS_ROOT",unset=file.path(tempdir(),"le8"))
  hd <- file.path(root,".revision_hashes");dir.create(hd,recursive=TRUE,showWarnings=FALSE)
  cf <- file.path(hd,paste0(le8_hash_object(path),".rds"))
  old <- if(file.exists(cf)) tryCatch(readRDS(cf),error=function(e)NULL) else NULL
  if(!is.null(old) && identical(old$key,key)) value <- old$value else {
    value <- unname(tools::md5sum(path))
    if(is.na(value)) stop("Cannot fingerprint input: ",path,call.=FALSE)
    base::saveRDS(list(key=key,value=value),cf)
  }
  assign(key,value,envir=h);value
}
le8_dependency <- function(path) {
  z <- path[!is.na(path)&nzchar(path)&file.exists(path)]
  if(!is.null(.le8_revision_state$rawdir)) {
    own<-paste0(normalizePath(.le8_revision_state$rawdir,winslash="/",mustWork=FALSE),"/")
    z<-z[!startsWith(normalizePath(z,winslash="/",mustWork=FALSE),own)]
  }
  .le8_revision_state$deps <- unique(c(.le8_revision_state$deps,
    normalizePath(z,winslash="/",mustWork=FALSE)))
  invisible(path)
}
le8_table <- function(x) if(is.null(x)) tibble::tibble() else tibble::as_tibble(x)
le8_append_workbook <- function(file,tables,prefix="review_") {
  if(!length(tables)) return(invisible(NULL))
  wb <- if(file.exists(file)) openxlsx::loadWorkbook(file) else openxlsx::createWorkbook()
  for(nm in names(tables)) {
    d <- as.data.frame(tables[[nm]])
    if(!ncol(d)) d <- data.frame(status="No estimable result / unavailable")
    sh <- substr(paste0(prefix,nm),1,31)
    if(sh %in% names(wb)) openxlsx::removeWorksheet(wb,sh)
    openxlsx::addWorksheet(wb,sh)
    if(nrow(d)) openxlsx::writeDataTable(wb,sh,d) else openxlsx::writeData(wb,sh,d)
    openxlsx::freezePane(wb,sh,firstRow=TRUE)
    openxlsx::setColWidths(wb,sh,cols=seq_len(ncol(d)),widths=18)
  }
  openxlsx::saveWorkbook(wb,file,overwrite=TRUE);invisible(file)
}

# This function is also inserted verbatim into 0f/0phe.f.R by the installer.
t2e <- function(dat,domain,Y_date,birth_date,date_attend,date_lost,date_death,
                date_end=date_follow_end,prefix=NA,time_unit="year") {
  n <- nrow(dat)
  getdate <- function(x) {
    if(inherits(x,"Date")) {
      if(length(x)==1L) return(rep(x,n))
      if(length(x)!=n) stop("Date vector length mismatch")
      return(x)
    }
    if(is.null(x)||length(x)!=1L||is.na(x)||!x%in%names(dat)) return(rep(as.Date(NA),n))
    z <- dat[[x]]
    if(is.numeric(z)&&!inherits(z,"Date")) as.Date(z,origin="1970-01-01") else as.Date(z)
  }
  ba<-getdate(date_attend);bd<-getdate(birth_date);yd<-getdate(Y_date)
  lost<-getdate(date_lost);dead<-getdate(date_death)
  end<-as.Date(date_end);if(length(end)==1L)end<-rep(end,n)
  if(length(end)!=n||anyNA(end))stop("Administrative end date must be valid",call.=FALSE)
  censor<-pmin(lost,dead,end,na.rm=TRUE)
  valid<-!is.na(ba)&!is.na(censor)&censor>=ba
  # Same-day diagnosis is baseline disease (never an incident non-case).
  prevalent<-!is.na(yd)&!is.na(ba)&yd<=ba
  incident<-valid&!is.na(yd)&yd>ba&yd<=censor
  beyond<-valid&!is.na(yd)&yd>censor
  dirty<-rep(FALSE,n)
  if(!is.null(domain)&&length(domain)==1L&&!is.na(domain)) {
    dv<-if(domain%in%names(dat))domain else paste0("icd10Ct_",domain)
    if(dv%in%names(dat)) {
      zz<-dat[[dv]]
      evidence<-if(inherits(zz,"Date"))!is.na(zz) else suppressWarnings(as.numeric(zz))>0
      dirty<-is.na(yd)&!is.na(evidence)&evidence
    }
  }
  yt<-ifelse(valid&!prevalent&!dirty,as.numeric(incident),NA_real_)
  yr<-ifelse(valid&!dirty&!incident,as.numeric(prevalent),NA_real_)
  stopdate<-censor;stopdate[incident]<-yd[incident]
  forward<-ifelse(!is.na(yt),as.numeric(stopdate-ba),NA_real_)
  # Backward clock remains descriptive; non-cases' r2e equals baseline age.
  reverse<-ifelse(!is.na(yr),ifelse(prevalent,as.numeric(ba-yd),as.numeric(ba-bd)),NA_real_)
  signed<-ifelse(prevalent,as.numeric(yd-ba),ifelse(incident,forward,NA_real_))
  exitdate<-stopdate;exitdate[prevalent]<-yd[prevalent]
  ageexit<-ifelse((prevalent|!is.na(yt))&!is.na(bd),as.numeric(exitdate-bd),NA_real_)
  divisor<-if(time_unit=="year")365.25 else if(time_unit%in%c("day","days"))1 else
    stop("time_unit must be year or day",call.=FALSE)
  vals<-list(Yt2e=yt,Yr2e=yr,t2e=forward/divisor,r2e=reverse/divisor,
             b2e=signed/divisor,bi2e=ageexit/divisor)
  nm<-names(vals);if(length(prefix)==1L&&!is.na(prefix))nm<-paste0(prefix,".",nm)
  dat[nm]<-vals
  attr(dat,"le8_outcome_audit")<-data.frame(metric=c("N","prevalent_including_same_day",
    "same_day","incident_within_followup","diagnosis_after_censor","invalid_followup","unknown_date_with_disease_evidence"),
    value=c(n,sum(prevalent),sum(!is.na(yd)&!is.na(ba)&yd==ba),sum(incident),
      sum(beyond),sum(!valid),sum(dirty)))
  attr(dat,"le8_t2e_version")<-"2026-09-05.censor-inclusive-baseline-v1"
  dat
}

make_outcome <- function(dat,outcome=Y) {
  req<-c("birth_date","date_attend",paste0("fod_icd10_",outcome))
  if(length(setdiff(req,names(dat))))stop("Missing outcome fields: ",paste(setdiff(req,names(dat)),collapse=", "))
  # Always rebuild: old phenotype RDS may contain pre-fix t2e variables.
  for(v in c("date_lost","date_death"))if(!v%in%names(dat))dat[[v]]<-as.Date(NA)
  ans<-t2e(dat,NA,paste0("fod_icd10_",outcome),"birth_date","date_attend",
    "date_lost","date_death",date_follow_end,outcome,"year")
  if(!is.null(.le8_revision_state$rawdir)) {
    au<-attr(ans,"le8_outcome_audit")
    au$outcome<-outcome
    data.table::fwrite(au,file.path(.le8_revision_state$rawdir,
      paste0("outcome_audit.",nrow(ans),".csv")))
  }
  ans
}
make_prevalent_status <- function(dat,outcome=Y) {
  y<-as.Date(dat[[paste0("fod_icd10_",outcome)]]);b<-as.Date(dat$date_attend)
  ifelse(is.na(b),NA_real_,as.numeric(!is.na(y)&y<=b))
}

le8_register_reads <- function() {
  # Wrap only project readers, not base readRDS. Preserve public interfaces.
  for(nm in c("read_all","read_prot","read_met","read_table_auto","read_sumstat","read_sumstat_snps","read_sumstat_region")) {
    if(!exists(nm,envir=.GlobalEnv,inherits=FALSE))next
    old<-get(nm,envir=.GlobalEnv)
    if(isTRUE(attr(old,"revision_wrapped")))next
    wrapper<-local({fn<-old;name<-nm;function(...) {
      args<-list(...)
      path<-switch(name,read_all=file.path(indir,"Rdata/all.rds"),
        read_prot=file.path(indir,"Rdata/prot.rds"),read_met=file.path(indir,"Rdata/met.rds"),
        if(length(args))args[[1]]else character())
      if(is.character(path))le8_dependency(path)
      fn(...)
    }})
    attr(wrapper,"revision_wrapped")<-TRUE;assign(nm,wrapper,envir=.GlobalEnv)
  }
}
le8_revision_begin <- function(layer,module) {
  outdir<-if(layer=="protein")out.prot else out.met
  rawdir<-file.path(outdir,module);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE)
  .le8_revision_state$rawdir<-rawdir
  .le8_revision_state$deps<-character()
  fdir<-Sys.getenv("LE8_FDIR",unset=file.path(dir0,"scripts/le8/f"))
  cf<-c(list.files(fdir,pattern="\\.(R|py|sh)$",full.names=TRUE),
    file.path(dir0,"scripts/0f",c("0phe.f.R","pred.f.R","0phe.f.sh")))
  inputs<-file.path(indir,"Rdata",c("all.rds",if(layer=="protein")c("prot.rds","prot.pgs.rds") else c("met.rds","met.pgs.rds")))
  inputs<-unique(c(inputs,prot_bed_file,met_list_file,
    Sys.getenv(c("C2_DANDELION_GENE_P_FILE","C2_MAGMA_GENE_FILE","C2_DANDELION_GENE_ANNOTATION",
      "C2_DANDELION_SNP_FILE","C2_DANDELION_SNP_GENE_MAP","C2_PROT_HERITABILITY_FILE",
      "C2_MET_HERITABILITY_FILE","LE8_GWAS_MANIFEST","LE8_LD_MANIFEST"),unset="")))
  prior<-file.path(rawdir,"revision_manifest.rds")
  old<-if(file.exists(prior))tryCatch(readRDS(prior),error=function(e)NULL)else NULL
  inputs<-unique(c(inputs,old$dependencies))
  # Upstream content, not merely 'available', enters each downstream fingerprint.
  up<-switch(module,c2_cause="c1_correlate",c3_coloc=c("c1_correlate","c2_cause"),
    c4_connect="c1_correlate",c5_consolidate=c("c1_correlate","c2_cause","c3_coloc","c4_connect"),character())
  for(u in up)inputs<-c(inputs,list.files(file.path(outdir,u),pattern="\\.res\\.rds$",full.names=TRUE))
  e<-Sys.getenv();e<-e[grepl("^(C[1-5]_|LE8_|RUN_|MRLINK|COLOC_|Y$|SEED$|N_CORES$|DATE_FOLLOW_END$)",names(e))]
  e<-e[!grepl("REPLACE|FDIR|JOB|REVISION",names(e))]
  inputs<-sort(unique(inputs[!is.na(inputs)&nzchar(inputs)]))
  base_sig<-list(version=LE8_REVISION,layer=layer,module=module,Y=Y,
    code=setNames(vapply(cf,le8_file_hash,character(1)),cf),settings=e[order(names(e))],
    covariates=list(basic=vars.basic,adj2=vars.adj2,le8=vars.le8),end=date_follow_end)
  digest<-le8_hash_object(list(base_sig,files=setNames(vapply(inputs,le8_file_hash,character(1)),inputs)))
  stale<-is.null(old)||!identical(old$signature,digest)
  existing<-list.files(rawdir,all.files=FALSE,full.names=TRUE)
  if(stale&&length(existing)) {
    stamp<-format(Sys.time(),"%Y%m%d-%H%M%S")
    hist<-file.path(outdir,"_history",paste0(stamp,"-",module))
    dir.create(hist,recursive=TRUE,showWarnings=FALSE)
    # Same-filesystem rename retains old RDS/CSV/figures without duplication.
    for(f in existing)if(!file.rename(f,file.path(hist,basename(f))))
      stop("Cannot archive stale output: ",f,call.=FALSE)
    prefix<-sub("_.*$","",module)
    pubs<-list.files(outdir,pattern=paste0("^",prefix,"[.].*\\.(png|xlsx|pdf|svg)$"),full.names=TRUE)
    for(f in pubs)if(!file.rename(f,file.path(hist,basename(f))))stop("Cannot archive figure: ",f)
    message("LE8 ",module,": incompatible cache archived at ",hist)
  }
  .le8_revision_state$base_sig<-base_sig
  .le8_revision_state$inputs<-inputs
  .le8_revision_state$manifest<-prior
  base::saveRDS(list(signature=digest,dependencies=inputs,version=LE8_REVISION),prior)
  invisible(NULL)
}
le8_revision_end <- function() {
  if(is.null(.le8_revision_state$manifest))return(invisible(NULL))
  inputs<-sort(unique(c(.le8_revision_state$inputs,.le8_revision_state$deps)))
  sig<-le8_hash_object(list(.le8_revision_state$base_sig,
    files=setNames(vapply(inputs,le8_file_hash,character(1)),inputs)))
  base::saveRDS(list(signature=sig,dependencies=inputs,version=LE8_REVISION),.le8_revision_state$manifest)
  invisible(NULL)
}
le8_load_revision <- function(layer,module) {
  fdir<-Sys.getenv("LE8_FDIR",unset=file.path(dir0,"scripts/le8/f"))
  files<-switch(module,c1_correlate="c1_revision_temporal.R",
    c2_cause=c("c1_revision_temporal.R","c2_revision_decompose.R","c2_revision_dandelion.R","c2_revision_evidence.R"),
    c3_coloc="c3_revision_locus.R",c4_connect="c4_revision_connection.R",
    c5_consolidate="c5_revision_validation.R",character())
  for(f in files)sys.source(file.path(fdir,f),envir=.GlobalEnv)
  # Shared formula keeps genetic association separate from downstream LE8.
  le8_register_reads()
  # Legacy drug.lipid may include post-baseline visits: do not silently use it.
  if(!truthy(Sys.getenv("LE8_TREATMENT_BASELINE_CONFIRMED",unset="FALSE"))) {
    for(nm in c("C1_TREATMENT_VARS"))if(exists(nm,.GlobalEnv,inherits=FALSE)) {
      z<-get(nm,.GlobalEnv)
      if("drug.lipid"%in%z) {
        warning("drug.lipid excluded: baseline-only provenance unconfirmed; use a baseline-only field or LE8_TREATMENT_BASELINE_CONFIRMED=TRUE after verification",call.=FALSE)
        assign(nm,setdiff(z,"drug.lipid"),.GlobalEnv)
      }
    }
    Sys.setenv(C2_TREATMENT_VARS=paste(setdiff(le8_csv_env("C2_TREATMENT_VARS",Sys.getenv("C1_TREATMENT_VARS")),"drug.lipid"),collapse=","))
  }
  le8_revision_begin(layer,module)
  invisible(NULL)
}

# An assay can differ from its encoding-gene symbol. User mappings override
# this small nomenclature map. Gene coordinates still come from annotation.
le8_assay_genes <- function(features) {
  out<-setNames(as.character(features),as.character(features))
  alias<-c(NTPROBNP="NPPB",NTproBNP="NPPB",`NT-proBNP`="NPPB")
  hit<-intersect(names(out),names(alias));out[hit]<-alias[hit]
  f<-Sys.getenv("LE8_ASSAY_GENE_MAP",unset="")
  if(nzchar(f)) {
    if(!file.exists(f))stop("LE8_ASSAY_GENE_MAP is missing: ",f)
    le8_dependency(f);m<-data.table::fread(f)
    if(!all(c("assay","gene")%in%names(m)))stop("Assay map needs assay and gene columns")
    if(anyDuplicated(m$assay))stop("Assay map has duplicate assay identifiers")
    hit<-intersect(names(out),m$assay);out[hit]<-m$gene[match(hit,m$assay)]
  }
  out
}
le8_same_locus_evidence <- function(mr,coloc,layer) {
  empty<-tibble::tibble(feature=character(),locus=character(),locus_class=character(),
    MR_FDR=numeric(),coloc_robust_min=numeric(),instrument_overlap=logical(),
    eligible=logical(),claim=character())
  if(!is.data.frame(mr)||!nrow(mr)||!is.data.frame(coloc)||!nrow(coloc))return(empty)
  need_m<-c("exposure","analysis","FDR_all","instrument_chr","instrument_pos_min","instrument_pos_max")
  need_c<-c("feature","locus","chr","start","end","status","PP.H4_robust_min","locus_class")
  if(!all(need_m%in%names(mr))||!all(need_c%in%names(coloc)))return(empty)
  primary<-if(layer=="protein")"cis"else"local"
  m<-mr[mr$analysis==primary&is.finite(mr$FDR_all),,drop=FALSE]
  c<-coloc[coloc$locus_class==primary&coloc$status=="ok",,drop=FALSE]
  if(!nrow(m)||!nrow(c))return(empty)
  z<-dplyr::inner_join(c,m,by=c("feature"="exposure"))
  # The interval must overlap actual retained MR instruments, not merely the
  # same feature. Single-signal ABF remains region-level, not signal-resolved.
  z$instrument_overlap<-vapply(seq_len(nrow(z)),function(i) {
    same<-as.character(z$chr[i])%in%strsplit(as.character(z$instrument_chr[i]),";",fixed=TRUE)[[1]]
    if(!same)return(FALSE)
    if("instrument_positions"%in%names(z)) {
      keys<-strsplit(as.character(z$instrument_positions[i]),";",fixed=TRUE)[[1]]
      pp<-strsplit(keys,":",fixed=TRUE)
      return(any(vapply(pp,function(q)length(q)==2L&&q[1]==as.character(z$chr[i])&&
        is.finite(suppressWarnings(as.numeric(q[2])))&&as.numeric(q[2])>=z$start[i]&&as.numeric(q[2])<=z$end[i],logical(1))))
    }
    # Old min/max ranges do not prove a particular retained SNP overlaps.
    FALSE
  },logical(1))
  z|>dplyr::transmute(feature,locus,locus_class,MR_FDR=FDR_all,
    coloc_robust_min=PP.H4_robust_min,instrument_overlap,
    eligible=instrument_overlap&is.finite(FDR_all)&FDR_all<.05&
      is.finite(PP.H4_robust_min)&PP.H4_robust_min>=.7,
    claim="Same cis/local region MR + ABF support; not a resolved causal signal or causal proof")
}
le8_same_locus_candidates <- function(mr,coloc,layer) {
  z<-le8_same_locus_evidence(mr,coloc,layer)
  unique(z$feature[z$eligible%in%TRUE])
}
# Optional additions fail locally and leave the already-produced core outputs.
le8_optional <- function(label,expr) {
  tryCatch(force(expr),error=function(e) {
    warning(label,": ",conditionMessage(e),call.=FALSE)
    rd<-.le8_revision_state$rawdir
    if(!is.null(rd))data.table::fwrite(data.frame(module=label,status="failed",
      message=conditionMessage(e)),file.path(rd,paste0(gsub("[^A-Za-z0-9_]","_",label),".status.csv")))
    list(status=tibble::tibble(status="failed",message=conditionMessage(e)))
  })
}
# Honor an explicit administrative end; an unconfigured run retains the
# project's source-defined date instead of guessing a more recent cutoff.
if(nzchar(Sys.getenv("DATE_FOLLOW_END",unset=""))) {
  .le8_end<-as.Date(Sys.getenv("DATE_FOLLOW_END"))
  if(is.na(.le8_end))stop("DATE_FOLLOW_END must be YYYY-MM-DD")
  date_follow_end<-.le8_end
}

le8_gwas_metadata <- function(file) {
  out<-list(N=NA_real_,N_case=NA_real_,N_control=NA_real_,sdY=NA_real_,
    type=NA_character_,beta_scale=NA_character_,build=NA_character_,
    source="unprovided",discovery_overlap="unknown")
  if(length(file)!=1L||is.na(file)||!nzchar(file))return(out)
  f<-Sys.getenv("LE8_GWAS_MANIFEST",unset="")
  if(!nzchar(f)||!file.exists(f))return(out)
  le8_dependency(f);m<-data.table::fread(f,showProgress=FALSE)
  if(!"file"%in%names(m))stop("LE8_GWAS_MANIFEST needs a file column")
  normal<-function(x)normalizePath(x,winslash="/",mustWork=FALSE)
  idx<-which(vapply(as.character(m$file),normal,character(1))==normal(file))
  if(length(idx)>1)stop("Duplicate exact file in GWAS manifest: ",file)
  if(!length(idx))return(out)
  for(n in intersect(names(out),names(m)))out[[n]]<-m[[n]][idx]
  for(n in c("N","N_case","N_control","sdY"))out[[n]]<-suppressWarnings(as.numeric(out[[n]]))
  if(!is.finite(out$N)&&is.finite(out$N_case)&&is.finite(out$N_control))out$N<-out$N_case+out$N_control
  out$source<-f;out
}
get_case_fraction <- function(outcome=Y) {
  explicit<-suppressWarnings(as.numeric(Sys.getenv("COLOC_CASE_FRAC",unset="NA")))
  if(is.finite(explicit)&&explicit>0&&explicit<1)return(explicit)
  f<-get_y_gwas_file(outcome,FALSE)
  if(is.na(f))return(NA_real_)
  m<-le8_gwas_metadata(f)
  if(is.finite(m$N_case)&&is.finite(m$N_control)&&m$N_case>0&&m$N_control>0)
    m$N_case/(m$N_case+m$N_control) else NA_real_
}

if(!exists(".le8_original_layer_annotation",envir=.GlobalEnv,inherits=FALSE)&&exists("layer_annotation",envir=.GlobalEnv))
  .le8_original_layer_annotation<-layer_annotation
layer_annotation <- function(layer,features) {
  z<-.le8_original_layer_annotation(layer,features)
  if(layer!="protein")return(z)
  genes<-le8_assay_genes(features);z$encoding_gene<-unname(genes[z$feature])
  aliases<-names(genes)[genes!=names(genes)]
  if(length(aliases)) {
    a<-.le8_original_layer_annotation(layer,unique(unname(genes[aliases])))
    for(i in which(z$feature%in%aliases&(!is.finite(z$start)|!is.finite(z$end)|is.na(z$chr)))) {
      j<-match(z$encoding_gene[i],a$feature)
      if(!is.na(j))for(nm in intersect(c("chr","start","end","pos"),names(z)))z[[nm]][i]<-a[[nm]][j]
    }
  }
  z
}
le8_finish_revision <- function(layer,module,env) {
  on.exit(le8_revision_end(),add=TRUE)
  # An early cache return has already-complete additions. A failed core run
  # must not be disguised by an expensive on-exit analysis.
  if(!exists("out",envir=env,inherits=FALSE))return(invisible(NULL))
  obj<-get("out",env);outdir<-if(layer=="protein")out.prot else out.met
  get0local<-function(n,default=NULL)if(exists(n,envir=env,inherits=FALSE))get(n,env)else default
  tables<-le8_optional(paste0(module,"_review_additions"),switch(module,
    c1_correlate=le8_c1_additions(get0local("dat"),layer,get0local("covs_adj2"),get0local("tvar"),get0local("evar")),
    c2_cause={
      dan<-obj$DANDELION%||%list();de<-obj$individual_decomposition%||%list()
      le8_dandelion_plot_bundle(dan,outdir)
      list(dan_targets_all=dan$targets_all%||%tibble(),dan_input_audit=dan$input_audit%||%tibble(),
        decomp_folds=de$folds%||%tibble(),decomp_summary=de$summary%||%tibble())
    },
    c3_coloc=le8_c3_additions(get0local("res",obj$summary),get0local("mr",tibble()),layer,outdir),
    c4_connect=le8_c4_additions(obj,outdir),
    c5_consolidate=le8_c5_additions(get0local("dat"),get0local("biom_vars"),get0local("ranked"),
      get0local("training_screen"),get0local("distalset",character()),get0local("geneticset",character()),
      get0local("clinical"),get0local("tvar"),get0local("evar"),get0local("split"),layer,outdir),list()))
  workbook<-file.path(outdir,if(module=="c5_consolidate")"c5.prediction_panels.xlsx"else paste0(sub("_.*$","",module),".out.xlsx"))
  le8_optional(paste0(module,"_review_workbook"),le8_append_workbook(workbook,tables))
  # Attach aggregate additions to the same result RDS so they are not orphaned
  # from C1-C5. Do not overwrite unknown files or serialized individual data.
  cache<-get0local("selected_cache",get0local("final_cache",get0local("cache")))
  if(is.character(cache)&&length(cache)==1L&&file.exists(cache)&&grepl("(res[.]rds)$",cache)) {
    obj$review<-tables;obj$meta$revision<-LE8_REVISION
    base::saveRDS(obj,cache,compress="xz")
  }
  le8_optional(paste0(module,"_output_index"),finalize_outputs(module,outdir))
  invisible(NULL)
}
