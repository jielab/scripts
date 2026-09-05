# MR uses marginal QTL and outcome effects. COJO bJ remains a PGS weight only.
# Save the original estimator once; it supplies the existing output schema.
if(!exists(".le8_original_run_mr",envir=.GlobalEnv,inherits=FALSE))
  .le8_original_run_mr <- run_mr

le8_ld_read <- function(path) {
  if(length(path)!=1L||is.na(path)||!file.exists(path))return(NULL)
  le8_dependency(path)
  z<-tryCatch(if(grepl("\\.rds$",path))readRDS(path)else
    read.table(path,header=TRUE,check.names=FALSE,stringsAsFactors=FALSE),error=function(e)NULL)
  if(is.list(z)&&!is.data.frame(z)&&!is.matrix(z))z<-z$R
  if(is.data.frame(z)) {
    if(ncol(z)==nrow(z)+1L){rn<-as.character(z[[1]]);z<-as.matrix(z[,-1,drop=FALSE]);rownames(z)<-rn}
    else z<-as.matrix(z)
  }
  if(!is.matrix(z)||nrow(z)!=ncol(z)||is.null(rownames(z))||is.null(colnames(z)))return(NULL)
  storage.mode(z)<-"double"
  if(anyDuplicated(rownames(z))||anyDuplicated(colnames(z))||
     !setequal(rownames(z),colnames(z)))return(NULL)
  z<-z[,rownames(z),drop=FALSE]
  if(any(!is.finite(z))||max(abs(z-t(z)))>1e-5||
     max(abs(diag(z)-1))>.01||max(abs(z))>1.001)return(NULL)
  z
}
le8_ld_greedy <- function(iv,R,r2=0.001) {
  stopifnot(is.matrix(R),r2>=0,r2<1)
  iv<-iv[order(iv$P,iv$SNP),,drop=FALSE]
  candidates<-intersect(iv$SNP,rownames(R));selected<-character()
  for(s in candidates)if(!length(selected)||all(R[s,selected]^2<=r2))selected<-c(selected,s)
  iv[match(selected,iv$SNP),,drop=FALSE]
}
le8_ld_for_iv <- function(iv,feature) {
  explicit<-Sys.getenv("C2_MR_LD_DIR",unset="")
  sf<-unique(iv$source_file);sf<-sf[!is.na(sf)&nzchar(sf)]
  dirs<-unique(c(explicit,dirname(sf)))
  paths<-unique(unlist(lapply(dirs,function(d)file.path(d,c(paste0(feature,".ld.rds"),
    paste0(feature,".ldr.cojo"),paste0(feature,".jma.ldr"),paste0(feature,".jma.cojo.ldr"))))))
  for(p in paths){R<-le8_ld_read(p);if(!is.null(R)&&sum(iv$SNP%in%rownames(R))>=2L)
    return(list(R=R,source=p))}
  list(R=NULL,source="unavailable")
}
le8_prepare_mr_iv <- function(iv,feature) {
  iv<-iv[is.finite(iv$BETA)&is.finite(iv$SE)&iv$SE>0&is.finite(iv$P),,drop=FALSE]
  if(!nrow(iv))return(list(iv=iv,status="no eligible variants",source="none",n_input=0L))
  # Even COJO-conditional discoveries must satisfy marginal relevance here.
  iv<-iv[iv$P<=le8_num_env("C2_MR_P",5e-8)&(iv$BETA/iv$SE)^2>=le8_num_env("C2_MR_MIN_F",10),,drop=FALSE]
  if(!nrow(iv))return(list(iv=iv,status="no marginally strong variants",source="none",n_input=0L))
  n0<-nrow(iv);ld<-le8_ld_for_iv(iv,feature)
  if(!is.null(ld$R)) {
    known<-iv[iv$SNP%in%rownames(ld$R),,drop=FALSE]
    kept<-le8_ld_greedy(known,ld$R,le8_num_env("C2_MR_LD_R2",.001))
    list(iv=kept,status="LD-pruned marginal IVs",source=ld$source,n_input=n0,
      n_without_LD=n0-nrow(known),max_retained_r2=if(nrow(kept)>1)
        max(ld$R[kept$SNP,kept$SNP][upper.tri(ld$R[kept$SNP,kept$SNP])]^2)else 0)
  }else {
    # Never assume COJO conditional signals are mutually uncorrelated.
    kept<-iv[order(iv$P,iv$SNP)[1L],,drop=FALSE]
    list(iv=kept,status=if(n0==1)"single marginal IV"else"single-IV fallback; LD unavailable",
      source="none",n_input=n0,n_without_LD=n0,max_retained_r2=NA_real_)
  }
}
read_qtl_instruments <- function(feature,base_dir,layer=c("protein","metabolite"),annotation=NULL) {
  layer<-match.arg(layer);fs<-find_qtl_files(feature,base_dir,layer)
  le8_dependency(unlist(fs[c("joint","full","cis")]))
  empty<-list(instruments=tibble(),score_instruments=tibble(),files=fs)
  if(is.na(fs$joint)||is.na(fs$full))return(empty)
  j<-read_sumstat(fs$joint,joint=TRUE)
  if(!nrow(j))return(empty)
  q<-read_sumstat_snps(fs$full,j$SNP)
  if(nrow(q)<nrow(j)) {
    jm<-tryCatch(read_sumstat_matched(fs$joint,fs$full,joint=TRUE),error=function(e)tibble())
    if(nrow(jm)){j<-jm;q<-read_sumstat_snps(fs$full,j$SNP)}
  }
  # Alleles, BETA, SE, P and N all come from the same marginal QTL record.
  iv<-q[q$SNP%in%j$SNP&!is.na(q$EA)&nzchar(q$EA)&!is.na(q$NEA)&nzchar(q$NEA),,drop=FALSE]
  if(!nrow(iv))return(empty)
  classify<-function(z) {
    if(layer=="protein") {
      a<-if(!is.null(annotation))annotation[annotation$feature==feature,,drop=FALSE]else NULL
      if(!is.null(a)&&nrow(a)&&all(is.finite(c(a$start[1],a$end[1])))&&!is.na(a$chr[1])) {
        pad<-le8_num_env("C2_CIS_WINDOW_BP",1e6)
        z$analysis<-ifelse(as.character(z$CHR)==as.character(a$chr[1])&
          z$POS>=a$start[1]-pad&z$POS<=a$end[1]+pad,"cis","trans")
        z$cis_annotation_status<-"gene coordinate defined"
      }else {z$analysis<-"unknown";z$cis_annotation_status<-"missing annotation; never inferred from lead SNP"}
    }else {
      lead<-iv[which.min(iv$P),,drop=FALSE]
      z$analysis<-ifelse(z$CHR==lead$CHR&abs(z$POS-lead$POS)<=le8_num_env("C2_LOCAL_WINDOW_BP",1e6),"local","distal")
      z$cis_annotation_status<-"metabolite local locus; not gene cis"
    }
    z
  }
  iv<-classify(iv);iv$effect_type<-"marginal"
  score<-tryCatch(recover_qtl_alleles(j,fs$full),error=function(e)j)
  score<-classify(score);score$effect_type<-"COJO_joint_PGS_only"
  list(instruments=iv,score_instruments=score,files=fs)
}
run_mr <- function(iv,ygwas,exposure,analysis) {
  # Also repairs reverse-MR callers that passed a raw COJO disease set.
  if(nrow(iv)&&"joint"%in%names(iv)&&any(iv$joint%in%TRUE)) {
    sources<-unique(iv$source_file);sources<-sources[!is.na(sources)]
    full<-if(length(sources))sub("\\.jma\\.cojo$",".gz",sources[1])else NA_character_
    marginal<-if(!is.na(full)&&file.exists(full))read_sumstat_snps(full,iv$SNP)else tibble()
    if(!nrow(marginal)) {
      warning("MR ",exposure,": joint effects cannot be repaired; no estimate",call.=FALSE)
      iv<-iv[0,,drop=FALSE]
    }else iv<-marginal
  }
  # Harmonize before LD pruning/fallback, so a missing top outcome SNP does
  # not hide another valid overlapping instrument.
  if(nrow(iv)) {
    h<-harmonize_sumstats(iv,ygwas)
    iv<-iv[iv$SNP%in%h$SNP,,drop=FALSE]
  }
  pp<-le8_prepare_mr_iv(iv,exposure)
  if(analysis=="unknown")pp$iv<-pp$iv[0,,drop=FALSE]
  ans<-.le8_original_run_mr(pp$iv,ygwas,exposure,analysis)
  ans$effect_type<-"marginal/marginal"
  ans$LD_status<-pp$status;ans$LD_source<-pp$source
  ans$n_IV_before_LD<-pp$n_input;ans$n_without_LD<-pp$n_without_LD%||%NA_integer_
  ans$max_retained_r2<-pp$max_retained_r2%||%NA_real_
  ans$steiger_method<-"per-IV t-statistic heuristic; not formal liability-scale Steiger"
  ans$steiger_formal_tested<-FALSE
  ans$steiger_r2_sum_exposure_audit<-ans$steiger_r2_exposure
  ans$steiger_r2_sum_outcome_audit<-ans$steiger_r2_outcome
  ans$steiger_r2_exposure<-NA_real_;ans$steiger_r2_outcome<-NA_real_
  ans$steiger_support<-NA # Do not let a heuristic masquerade as a formal test.
  ans$tsmr_verified<-is.finite(ans$tsmr_ivw_b)&is.finite(ans$b)&
    abs(ans$tsmr_ivw_b-ans$b)<1e-8&is.finite(ans$tsmr_ivw_p)&
    abs(ans$tsmr_ivw_p-ans$pval)<1e-5
  if(nrow(pp$iv)) {
    ans$instrument_snps<-paste(pp$iv$SNP,collapse=";")
    ans$instrument_chr<-paste(unique(pp$iv$CHR),collapse=";")
    ans$instrument_positions<-paste(paste(pp$iv$CHR,pp$iv$POS,sep=":"),collapse=";")
    ans$instrument_pos_min<-min(pp$iv$POS);ans$instrument_pos_max<-max(pp$iv$POS)
  }else {ans$instrument_snps<-"";ans$instrument_chr<-"";ans$instrument_positions<-"";ans$instrument_pos_min<-NA_real_;ans$instrument_pos_max<-NA_real_}
  ans
}
