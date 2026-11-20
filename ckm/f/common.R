# Shared CKM IO. No data are read or written when this file is sourced.
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

ckm_require <- function() {
  for (p in c("data.table", "survival", "jsonlite"))
    if (!requireNamespace(p, quietly=TRUE)) stop("Missing R package: ", p, ". Run install.sh first.")
  data.table::setDTthreads(as.integer(Sys.getenv("CKM_THREADS", "4")))
}

ckm_expand <- function(s) {
  if (!is.character(s)) return(s)
  vapply(s, function(z) {
    m <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", z)[[1]]
    if (m[1] > 0L) {
      tokens <- regmatches(z, list(m))[[1]]
      for (tok in unique(tokens)) {
        key <- substr(tok, 3, nchar(tok)-1); value <- Sys.getenv(key, unset=NA_character_)
        if (is.na(value)) stop("Unset path variable ", tok, "; run through ckm.sh or set this path explicitly.")
        z <- gsub(tok, value, z, fixed=TRUE)
      }
    }
    path.expand(z)
  }, character(1), USE.NAMES=FALSE)
}

ckm_config <- function(file, overrides=list()) {
  if (!file.exists(file)) stop("Config not found: ", file)
  cfg <- jsonlite::read_json(file, simplifyVector=TRUE)
  cfg$config_file <- normalizePath(file)
  cfg$config_dir <- dirname(cfg$config_file)
  cfg$paths <- lapply(cfg$paths, ckm_expand)
  for (nm in c("manifest", "annotation", "gmt")) {
    p <- cfg$paths[[nm]]
    if (!is.null(p) && nzchar(p) && !grepl("^(/|[A-Za-z]:)", p)) cfg$paths[[nm]] <- file.path(cfg$config_dir,p)
  }
  cfg$out_dir <- overrides$out_dir %||% ckm_expand(cfg$out_dir)
  if (!grepl("^(/|[A-Za-z]:)",cfg$out_dir)) cfg$out_dir <- file.path(getwd(),cfg$out_dir)
  cfg$out_dir <- normalizePath(cfg$out_dir,mustWork=FALSE)
  cfg$phenotype <- cfg$phenotype %||% list()
  for(nm in c("definition","stage3_method")) if (!is.null(overrides[[nm]])) cfg$phenotype[[nm]] <- overrides[[nm]]
  if (!is.null(overrides$admin_end) || (!is.null(cfg$phenotype$admin_end) && any(nzchar(cfg$phenotype$admin_end))))
    stop("CKM v2 uses date_follow_end from 0f/phenotype.R. Remove the old admin_end override; change the shared variable when updating the study cutoff.")
  cfg$phenotype$admin_end <- NULL
  cfg$shared <- ckm_shared_info()
  for (section in c("association", "markers")) {
    covariates <- cfg[[section]]$covariates
    if (is.null(covariates) || identical(unname(unlist(covariates)), "vars.basic"))
      cfg[[section]]$covariates <- ckm_shared()$vars.basic
  }
  if(!is.null(overrides$ukb_dir)) {
    old <- cfg$paths$ukb_dir
    cfg$paths$ukb_dir <- overrides$ukb_dir
    for(nm in names(cfg$paths)) if(nm!="ukb_dir" && is.character(cfg$paths[[nm]]) && length(cfg$paths[[nm]])==1 && startsWith(cfg$paths[[nm]],paste0(old,"/")))
      cfg$paths[[nm]] <- paste0(overrides$ukb_dir,substring(cfg$paths[[nm]],nchar(old)+1))
  }
  if(!is.null(overrides$domains)) cfg$domains <- strsplit(overrides$domains,",",fixed=TRUE)[[1]]
  cfg$threads <- as.integer(overrides$threads %||% cfg$threads %||% 4L)
  cfg$association$workers <- cfg$threads
  cfg$seed <- as.integer(cfg$seed %||% 20261113L)
  if(any(!cfg$domains %in% c("lifestyle","prot","met"))) stop("domains must be lifestyle,prot,met")
  dir.create(cfg$out_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(cfg$out_dir,"tmp"), recursive=TRUE, showWarnings=FALSE)
  cfg
}

ckm_atomic_rds <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temp <- tempfile(".ckm-", tmpdir=dirname(path)); on.exit(unlink(temp), add=TRUE)
  saveRDS(x,temp,compress=FALSE)
  if(!file.rename(temp,path)) stop("Could not install ",path)
  invisible(path)
}

ckm_write <- function(x,path) {
  if(is.null(x)) x <- data.frame(status="not_run")
  x <- data.table::as.data.table(x)
  if(!ncol(x)) x <- data.table::data.table(status=character())
  temp <- tempfile(".ckm-",tmpdir=dirname(path));on.exit(unlink(temp),add=TRUE)
  data.table::fwrite(x,temp,sep="\t",na="NA")
  if(!file.rename(temp,path))stop("Could not install ",path)
  invisible(path)
}

ckm_table <- function(path,label="input") {
  if(is.null(path)||!nzchar(path)||!file.exists(path))stop("Missing ",label,": ",path)
  x <- if(grepl("\\.rds$",path,ignore.case=TRUE))readRDS(path) else data.table::fread(path)
  if(is.matrix(x))x<-as.data.frame(x)
  if(!is.data.frame(x))stop(label," must be a table: ",path)
  x<-data.table::as.data.table(x)
  if(!"eid" %in% names(x))stop(label," needs an eid column")
  x[,eid:=as.character(eid)]
  if(anyNA(x$eid)||any(!nzchar(x$eid))||anyDuplicated(x$eid))stop(label," has missing or duplicate participant IDs")
  x
}

ckm_domain <- function(d,cfg,domain) {
  manifest <- data.table::fread(cfg$paths$manifest,na.strings=c("NA",""))
  required<-c("feature","source","column","type","transform")
  if(length(setdiff(required,names(manifest))))stop("Exposure manifest requires ",paste(required,collapse=", "))
  manifest<-manifest[source==domain]
  if(domain=="lifestyle") {
    x<-if(nzchar(cfg$paths$lifestyle %||% ""))ckm_table(cfg$paths$lifestyle,"lifestyle") else data.table::copy(d)
    needed<-intersect(manifest$column,names(x))
    # Missing planned columns remain in the manifest, so failed tests stay in the FDR family.
    x<-x[,c("eid",needed),with=FALSE]
  } else {
    x<-ckm_table(cfg$paths[[domain]],domain)
    wildcard<-manifest[column=="*"]
    if(nrow(wildcard)) {
      cols<-setdiff(names(x),"eid")
      numeric_cols<-cols[vapply(x[,..cols],is.numeric,logical(1))]
      cols<-numeric_cols[!grepl("^(date_|age$|sex$|batch$|plate$|instance$|sample_id$)",numeric_cols)]
      if(!length(cols))stop("No numeric assays in ",domain)
      auto<-wildcard[rep(1L,length(cols))]
      auto[,`:=`(feature=cols,column=cols)]
      manifest<-data.table::rbindlist(list(manifest[column!="*"],auto),fill=TRUE)
    }
    if(anyDuplicated(manifest$feature))stop("Repeated feature in ",domain," manifest")
    x<-x[eid %in% d$eid,c("eid",intersect(manifest$column,names(x))),with=FALSE]
    keep<-setdiff(names(x),"eid")
    if(length(keep))x<-x[rowSums(!is.na(as.matrix(x[,..keep])))>0L]
  }
  if(!nrow(manifest))stop("No manifest entries for requested domain ",domain)
  # Optional annotation supplies the lipid-block flag; no name heuristic is used.
  if(nzchar(cfg$paths$annotation %||% "") && file.exists(cfg$paths$annotation)) {
    ann<-data.table::fread(cfg$paths$annotation)
    if(all(c("source","feature","lipid_block") %in% names(ann))) {
      ann<-unique(ann[,.(source,feature,lipid_block)])
      if(anyDuplicated(ann[,.(source,feature)]))stop("Conflicting lipid_block annotation")
      flag<-ann$lipid_block[match(paste(manifest$source,manifest$feature),paste(ann$source,ann$feature))]
      if(!"lipid_block" %in% names(manifest))manifest[,lipid_block:=FALSE]
      ii<-which(!is.na(flag));manifest[ii,lipid_block:=as.logical(flag[ii])]
    }
  }
  if(domain=="met" && isTRUE(cfg$exclude_met_lipid_block)) {
    if(!nzchar(cfg$paths$annotation %||% "") || !file.exists(cfg$paths$annotation))
      stop("exclude_met_lipid_block requires an explicit metabolite annotation file")
    ann<-data.table::fread(cfg$paths$annotation)
    if(!all(c("source","feature","lipid_block") %in% names(ann)))stop("Annotation needs source,feature,lipid_block")
    flag<-ann$lipid_block[match(paste(manifest$source,manifest$feature),paste(ann$source,ann$feature))]
    flag<-as.logical(flag)
    if(anyNA(flag))stop("Every planned metabolite needs an explicit TRUE/FALSE lipid_block annotation for this sensitivity")
    manifest<-manifest[!flag]
    x<-x[,c("eid",intersect(manifest$column,names(x))),with=FALSE]
    if(!nrow(manifest))stop("No metabolites remain after the configured lipid-block exclusion")
  }
  list(manifest=manifest,data=x)
}

ckm_phenotype_summary <- function(d) {
  a<-data.table::as.data.table(d)
  stages<-a[,.(n=.N,events=sum(status==1,na.rm=TRUE),deaths=sum(event_death==1,na.rm=TRUE),competing_deaths=sum(status==2,na.rm=TRUE)),by=stage][order(stage)]
  flags<-intersect(c("stage_core_complete","incident_eligible","riskset_eligible","stage_exact","clinical_cvd_without_current_risk"),names(a))
  complete<-data.table::rbindlist(lapply(flags,function(k)data.table::data.table(name=k,n=sum(a[[k]] %in% TRUE),percent=100*mean(a[[k]] %in% TRUE))),fill=TRUE)
  vars<-intersect(c("age","sex","bmi","sbp","dbp","hba1c","glucose","tg","hdl","total_chol","egfr","uacr","prevent_10y","score2_10y"),names(a))
  miss<-data.table::rbindlist(lapply(vars,function(k)data.table::data.table(variable=k,missing=sum(is.na(a[[k]])),n=nrow(a))))
  list(n_total=nrow(a),stage_counts=stages,completeness=complete,missingness=miss,
       stage_comparison=if(all(c("stage_prevent","stage_score2") %in% names(a)))a[,.(n=.N),by=.(stage_prevent,stage_score2)] else data.frame(),
       notes=c("Baseline stages use available pre-baseline evidence. Unknown components remain unknown.","Stage 4 is operational clinical CVD; baseline stage 4 is excluded from first-incident-CVD models.","Deaths counts all observed follow-up deaths, including deaths after CVD. Competing deaths counts deaths before first incident CVD.","A single eGFR/UACR assessment establishes a CKD-risk phenotype, not confirmed chronicity."))
}
