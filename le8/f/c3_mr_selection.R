# C3 candidates are significant C2 analyses and their actual retained instruments.
# This avoids re-running COJO coordinate matching (including PGS-only preparation).
C3_SELECTION_VERSION <- "2026-09-19.mr-fdr-v1"
c3_file_stamp <- function(paths) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  info <- file.info(paths)
  data.frame(path=normalizePath(paths,mustWork=FALSE),size=info$size,
             mtime=as.numeric(info$mtime),ctime=as.numeric(info$ctime))
}
c3_select_mr <- function(mr,layer,fdr=.05,max_features=200L) {
  required <- c("exposure","analysis","pval","FDR_all","instrument_snps","instrument_positions")
  if(!all(required %in% names(mr)))
    stop("C3 requires C2 MR results with FDR_all and retained instrument IDs/positions; rerun C2.",call.=FALSE)
  classes <- if(layer=="protein")c("cis","trans")else c("local","distal")
  z <- mr |> filter(is.finite(FDR_all),FDR_all<fdr,is.finite(pval),
                    !is.na(exposure),nzchar(exposure),analysis %in% classes) |>
    arrange(FDR_all,pval,exposure,analysis)
  if(any(is.na(z$instrument_snps)|!nzchar(z$instrument_snps)))
    stop("Significant C2 MR rows lack retained instrument IDs; rerun C2.",call.=FALSE)
  features <- unique(z$exposure)
  if(length(features)>max_features)
    message("C3: ",length(features)," MR-significant features; C3_MAX_FEATURES retains ",max_features)
  z |> filter(exposure %in% head(features,max_features))
}
c3_mr_instruments <- function(mr) {
  bind_rows(lapply(seq_len(nrow(mr)),function(i) {
    ids <- trimws(strsplit(mr$instrument_snps[i],";",fixed=TRUE)[[1]])
    pos <- strsplit(ifelse(is.na(mr$instrument_positions[i]),"",mr$instrument_positions[i]),";",fixed=TRUE)[[1]]
    if(length(pos)!=length(ids))pos<-rep(NA_character_,length(ids))
    tibble(SNP=ids,analysis=mr$analysis[i],mr_position=pos)
  })) |> distinct()
}
c3_read_mr_qtl <- function(mr,qfile,cache_root) {
  wanted <- c3_mr_instruments(mr)
  key <- le8_hash_object(list(version=C3_SELECTION_VERSION,file=c3_file_stamp(qfile),
    wanted=wanted,metadata=c3_file_stamp(Sys.getenv("LE8_GWAS_MANIFEST","")),N=Sys.getenv("C2_SUMSTAT_N","100000")))
  cache <- file.path(cache_root,paste0(key,".rds"))
  iv <- read_stage_cache(cache)
  if(is.data.frame(iv)) {
    message("  QTL instrument cache hit: ",nrow(iv)," records")
    return(iv)
  }
  # Positions in C2 are from marginal QTL records, already on the QTL build.
  # Query indexed intervals first; missing/legacy coordinates get one ID scan.
  positions <- unique(wanted$mr_position)
  positions <- positions[!is.na(positions)&grepl("^[^:]+:[0-9]+$",positions)]
  tabix_ok <- nzchar(Sys.which("tabix")) && any(vapply(paste0(qfile,c(".tbi",".csi")),
    function(f)file.exists(f)&&file.info(f)$mtime>=file.info(qfile)$mtime,logical(1)))
  q <- tibble(SNP=character())
  if(tabix_ok && length(positions)) {
    message("  QTL indexed lookup: ",length(positions)," MR instrument positions")
    q <- bind_rows(lapply(positions,function(p) {
      s<-strsplit(p,":",fixed=TRUE)[[1]]
      read_sumstat_region(qfile,s[1],as.numeric(s[2]),as.numeric(s[2]))
    }))
  }
  missing <- setdiff(wanted$SNP,q$SNP)
  if(length(missing)) {
    message("  QTL ID scan: ",length(missing)," instruments absent from indexed lookup")
    q <- bind_rows(q,read_sumstat_snps(qfile,missing))
  }
  if(!nrow(q))q<-tibble(SNP=character(),CHR=character(),POS=numeric(),P=numeric())
  q <- q |> distinct(SNP,.keep_all=TRUE)
  iv <- wanted |> inner_join(q,by="SNP")
  lost <- setdiff(wanted$SNP,iv$SNP)
  if(length(lost))stop("C3 cannot recover retained MR instruments from ",qfile,": ",
    paste(lost,collapse=", "),"; check C2/QTL provenance.",call.=FALSE)
  if(any(!is.finite(iv$POS)|is.na(iv$CHR)|!is.finite(iv$P)))
    stop("C3 retained MR instruments have invalid QTL coordinates/P values: ",qfile,call.=FALSE)
  write_stage_cache(iv,cache)
  iv
}
c3_locus_settings <- function() list(window=WINDOW_BP,min_snps=MIN_SNPS,p12=C3_P12,
  susie=Sys.getenv("C3_RUN_SUSIE","FALSE"),susie_max=Sys.getenv("C3_SUSIE_MAX_SNPS","5000"),
  ld=c3_file_stamp(Sys.getenv("LE8_LD_MANIFEST","")),
  metadata=c3_file_stamp(Sys.getenv("LE8_GWAS_MANIFEST","")),N=Sys.getenv("C2_SUMSTAT_N","100000"))
c3_cached_locus <- function(cache,legacy_files,feature,layer,chr,pos,cls,qfile,yfile,outtype,sfrac) {
  z<-read_stage_cache(cache)
  if(!is.null(z))return(z)
  # Legacy caches have no input signature. Migrate only default ABF settings,
  # unchanged input files and exactly matching geometry/class; never a SuSiE run.
  legacy_ok<-!LE8_REPLACE && WINDOW_BP==500000 && MIN_SNPS==50 &&
    identical(unname(C3_P12),c(1e-6,1e-5,1e-4)) && !truthy(Sys.getenv("C3_RUN_SUSIE","FALSE")) &&
    Sys.getenv("C2_SUMSTAT_N","100000")=="100000" && !nzchar(Sys.getenv("LE8_GWAS_MANIFEST","")) && !is.finite(sfrac)
  if(!legacy_ok)return(NULL)
  suffix<-paste0("_",gsub("[^A-Za-z0-9._-]","_",feature),"_chr",chr,"_",format(pos,scientific=FALSE,trim=TRUE),".rds")
  paths<-legacy_files[endsWith(basename(legacy_files),suffix)]
  for(f in paths) {
    if(any(file.info(c(qfile,yfile))$mtime>file.info(f)$mtime))next
    old<-read_stage_cache(f);s<-old$summary
    if(is.data.frame(s)&&nrow(s)==1L&&all(c("locus_class","case_fraction_source")%in%names(s))&&
       identical(as.character(s$feature),feature)&&identical(as.character(s$layer),layer)&&
       identical(as.character(s$locus_class),cls)&&as.character(s$chr)==as.character(chr)&&
       s$lead_pos==pos&&s$start==max(1,pos-WINDOW_BP)&&s$end==pos+WINDOW_BP&&
       s$case_fraction_source=="not required for beta/varbeta cc ABF"&&
       all(c("variants","regional")%in%names(old))) {
      write_stage_cache(old,cache);return(old)
    }
  }
  NULL
}
