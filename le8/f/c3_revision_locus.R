# Region-coherent colocalization; ABF-H4 conditional sets are not fine-mapping.
# Optional SuSiE requires dense, explicitly allele-labelled LD and actual N.
select_nonoverlapping_leads <- function(iv,max_loci=MAX_LOCI_PER_FEATURE,window_bp=WINDOW_BP) {
  iv<-iv|>filter(is.finite(POS),!is.na(CHR),is.finite(P))
  if(!nrow(iv))return(iv)
  # Reserve the first slot for a genuine cis/local signal, never a trans proxy.
  iv<-iv|>mutate(.priority=ifelse(analysis%in%c("cis","local"),0L,1L))|>arrange(.priority,P,SNP)
  keep<-integer()
  for(i in seq_len(nrow(iv))) {
    overlap<-length(keep)&&any(iv$CHR[keep]==iv$CHR[i]&abs(iv$POS[keep]-iv$POS[i])<=2*window_bp)
    if(!overlap)keep<-c(keep,i)
    if(length(keep)>=max_loci)break
  }
  iv[keep,,drop=FALSE]|>select(-.priority)
}
le8_locus_class <- function(feature,layer,chr,pos) {
  if(layer=="metabolite") {
    fs<-find_qtl_files(feature,dir.met.gwas,layer)
    if(is.na(fs$joint))return("unknown")
    iv<-read_qtl_instruments(feature,dir.met.gwas,layer)$instruments
    if(!nrow(iv))return("unknown")
    lead<-iv[which.min(iv$P),,drop=FALSE]
    return(if(chr==lead$CHR&&abs(pos-lead$POS)<=le8_num_env("C2_LOCAL_WINDOW_BP",1e6))"local"else"distal")
  }
  a<-layer_annotation(layer,feature)
  if(!nrow(a)||is.na(a$chr[1])||!is.finite(a$start[1])||!is.finite(a$end[1]))return("unknown")
  pad<-le8_num_env("C2_CIS_WINDOW_BP",1e6)
  if(as.character(chr)==as.character(a$chr[1])&&pos>=a$start[1]-pad&&pos<=a$end[1]+pad)"cis"else"trans"
}
le8_coloc_dataset <- function(d,which,type,file,case_frac=NA_real_) {
  suffix<-if(which==1L)"x"else"y";m<-le8_gwas_metadata(file)
  beta<-d[[paste0("BETA_",suffix)]];se<-d[[paste0("SE_",suffix)]]
  z<-list(beta=beta,varbeta=se^2,snp=d$SNP,position=d[[paste0("POS_",suffix)]],type=type)
  n<-d[[paste0("N_",suffix)]];n<-n[is.finite(n)&n>0]
  N<-if(length(n))median(n)else m$N
  if(is.finite(N)&&N>0)z$N<-N
  maf<-pmin(d[[paste0("EAF_",suffix)]],1-d[[paste0("EAF_",suffix)]])
  if(all(is.finite(maf)&maf>0&maf<=.5))z$MAF<-maf
  if(type=="cc") {
    # beta/varbeta case-control ABF does not require s. Never substitute UKB
    # incident prevalence for the discovery GWAS case fraction.
    if(is.finite(case_frac)&&case_frac>0&&case_frac<1)z$s<-case_frac
  }else {
    if(is.finite(m$sdY)&&m$sdY>0)z$sdY<-m$sdY
    if(is.null(z$sdY)&&(is.null(z$N)||is.null(z$MAF)))
      stop("Quantitative coloc needs documented sdY or real N plus MAF: ",file)
  }
  z
}
le8_dense_ld <- function(trait,chr,start,end,snps,EA,NEA) {
  mf<-Sys.getenv("LE8_LD_MANIFEST",unset="")
  if(!nzchar(mf)||!file.exists(mf))return(list(status="dense LD manifest unavailable"))
  le8_dependency(mf);m<-data.table::fread(mf,showProgress=FALSE)
  req<-c("trait","chr","start","end","file")
  if(!all(req%in%names(m)))return(list(status="LD manifest requires trait,chr,start,end,file"))
  ii<-which(m$trait==trait&as.character(m$chr)==as.character(chr)&m$start<=start&m$end>=end)
  if(!length(ii))return(list(status=paste("no covering dense LD for",trait)))
  j<-ii[which.min(m$end[ii]-m$start[ii])];f<-as.character(m$file[j]);le8_dependency(f)
  z<-tryCatch(readRDS(f),error=function(e)NULL)
  if(!is.list(z)||is.null(z$R)||is.null(z$alleles))return(list(status="LD RDS must contain R and alleles"))
  R<-z$R;a<-as.data.frame(z$alleles)
  if(!all(c("SNP","EA","NEA")%in%names(a))||anyDuplicated(a$SNP)||
    !is.matrix(R)||nrow(R)!=ncol(R)||is.null(rownames(R))||!identical(rownames(R),colnames(R)))
    return(list(status="invalid LD names/alleles"))
  if(!all(snps%in%rownames(R))||!all(snps%in%a$SNP))return(list(status="incomplete dense LD coverage; no silent SNP thinning"))
  a<-a[match(snps,a$SNP),,drop=FALSE];R<-R[snps,snps,drop=FALSE]
  same<-toupper(a$EA)==toupper(EA)&toupper(a$NEA)==toupper(NEA)
  flip<-toupper(a$EA)==toupper(NEA)&toupper(a$NEA)==toupper(EA)
  if(anyNA(same|flip)||any(!(same|flip)))return(list(status="LD/summary allele mismatch"))
  sign<-ifelse(same,1,-1);R<-R*outer(sign,sign)
  if(any(!is.finite(R))||max(abs(R-t(R)))>1e-6||max(abs(diag(R)-1))>1e-4)
    return(list(status="invalid dense LD correlations"))
  # chol permits positive semidefinite external LD only after a tiny numeric
  # ridge. No shrinkage is silently used to rescue a materially indefinite R.
  ev<-tryCatch(min(eigen(R,symmetric=TRUE,only.values=TRUE)$values),error=function(e)NA_real_)
  if(!is.finite(ev)||ev< -1e-5)return(list(status="dense LD is not positive semidefinite"))
  list(status="ok",R=R,file=f,n_ref=z$n_ref%||%NA_integer_,build=z$build%||%NA_character_,
    ancestry=z$ancestry%||%NA_character_)
}
le8_run_susie <- function(dx,dy,d,feature,chr,start,end,locus) {
  no<-list(status="disabled",summary=tibble(),pip=tibble())
  if(!truthy(Sys.getenv("C3_RUN_SUSIE",unset="FALSE")))return(no)
  if(!requireNamespace("susieR",quietly=TRUE))return(modifyList(no,list(status="susieR missing")))
  if(is.null(dx$N)||is.null(dy$N))return(modifyList(no,list(status="actual discovery N required for each SuSiE trait")))
  if(length(dx$snp)>le8_num_env("C3_SUSIE_MAX_SNPS",5000))
    return(modifyList(no,list(status="locus exceeds configured dense-LD memory cap; not thinned")))
  lx<-le8_dense_ld(feature,chr,start,end,d$SNP,d$EA_x,d$NEA_x)
  ly<-le8_dense_ld(Y,chr,start,end,d$SNP,d$EA_x,d$NEA_x) # harmonization orients both betas to EA_x
  if(lx$status!="ok"||ly$status!="ok")return(modifyList(no,list(status=paste(lx$status,ly$status,sep="; "))))
  dx$LD<-lx$R;dy$LD<-ly$R
  ans<-tryCatch({
    fx<-coloc::runsusie(dx,maxit=1000,estimate_residual_variance=FALSE)
    fy<-coloc::runsusie(dy,maxit=1000,estimate_residual_variance=FALSE)
    if(!isTRUE(fx$converged)||!isTRUE(fy$converged))stop("SuSiE did not converge")
    fits<-lapply(C3_P12,function(p12)coloc::coloc.susie(fx,fy,p1=1e-4,p2=1e-4,p12=p12))
    sm<-imap_dfr(fits,function(ff,nm)as_tibble(ff$summary)|>mutate(prior=nm,feature=feature,locus=locus))
    pp<-tibble(feature,locus,SNP=d$SNP,QTL_PIP=fx$pip,outcome_PIP=fy$pip)
    rd<-.le8_revision_state$rawdir;dir.create(file.path(rd,"susie"),showWarnings=FALSE)
    stem<-paste0(gsub("[^A-Za-z0-9_.-]","_",paste(feature,locus)),".rds")
    saveRDS(list(QTL=fx,outcome=fy,coloc=fits,summary=sm,pip=pp,LD_QTL=lx$file,LD_outcome=ly$file),file.path(rd,"susie",stem))
    list(status="ok; signal-pair posterior available",summary=sm,pip=pp)
  },error=function(e)modifyList(no,list(status=conditionMessage(e))))
  ans
}
coloc_one_locus <- function(feature,layer,qtl_file,lead_chr,lead_pos,ygwas_file,outcome_type,case_frac) {
  start<-max(1,lead_pos-WINDOW_BP);end<-lead_pos+WINDOW_BP
  locus<-paste0("chr",lead_chr,":",floor(start),"-",ceiling(end))
  cls<-le8_locus_class(feature,layer,lead_chr,lead_pos)
  s<-tibble(layer,feature,locus,chr=as.character(lead_chr),start,end,lead_pos,locus_class=cls,
    status="failed",n_snps=0L,PP.H0=NA_real_,PP.H1=NA_real_,PP.H2=NA_real_,PP.H3=NA_real_,PP.H4=NA_real_,
    PP.H4_p12_conservative=NA_real_,PP.H4_p12_default=NA_real_,PP.H4_p12_liberal=NA_real_,PP.H4_robust_min=NA_real_,
    lead_shared=NA_character_,lead_shared_pp=NA_real_,credible_set_n=NA_integer_,message="",
    credible_set_interpretation="SNP posterior CONDITIONAL ON H4; not a trait fine-mapping credible set",
    susie_status="not attempted",case_fraction_source=if(is.finite(case_frac))"configured discovery GWAS"else"not required for beta/varbeta cc ABF")
  fail<-function(msg){s$message<-msg;list(summary=s,variants=tibble(),regional=tibble())}
  le8_dependency(c(qtl_file,ygwas_file))
  q<-read_sumstat_region(qtl_file,lead_chr,start,end);y<-read_sumstat_region(ygwas_file,lead_chr,start,end)
  d<-harmonize_sumstats(q,y)
  if(!nrow(d))return(fail("no allele-aligned overlapping SNPs"))
  d<-d|>filter(is.finite(BETA_x),is.finite(SE_x),SE_x>0,is.finite(BETA_y),is.finite(SE_y),SE_y>0)
  s$n_snps<-nrow(d)
  if(nrow(d)<MIN_SNPS)return(fail("insufficient dense overlap"))
  mq<-le8_gwas_metadata(qtl_file);my<-le8_gwas_metadata(ygwas_file)
  if(!is.na(mq$build)&&!is.na(my$build)&&as.character(mq$build)!=as.character(my$build))return(fail("declared discovery builds differ; harmonize coordinates before C3"))
  datasets<-tryCatch(list(x=le8_coloc_dataset(d,1L,"quant",qtl_file),
    y=le8_coloc_dataset(d,2L,outcome_type,ygwas_file,case_frac)),error=function(e)e)
  if(inherits(datasets,"condition"))return(fail(conditionMessage(datasets)))
  fits<-lapply(C3_P12,function(p12)tryCatch({
    f<-NULL;invisible(capture.output(f<-coloc::coloc.abf(datasets$x,datasets$y,p1=1e-4,p2=1e-4,p12=p12)));f
  },error=function(e)e))
  f<-fits[["default"]];if(inherits(f,"condition"))return(fail(conditionMessage(f)))
  for(i in 0:4)s[[paste0("PP.H",i)]]<-as.numeric(f$summary[[paste0("PP.H",i,".abf")]])
  pp<-vapply(fits,function(x)if(inherits(x,"condition"))NA_real_ else as.numeric(x$summary[["PP.H4.abf"]]),numeric(1))
  s$PP.H4_p12_conservative<-pp["conservative"];s$PP.H4_p12_default<-pp["default"];s$PP.H4_p12_liberal<-pp["liberal"]
  s$PP.H4_robust_min<-if(all(is.finite(pp)))min(pp)else NA_real_ # missing prior is not a pass
  v<-as_tibble(f$results)|>arrange(desc(SNP.PP.H4))|>
    mutate(cum_pp=cumsum(SNP.PP.H4),credible95=lag(cum_pp,default=0)<.95,
      feature=feature,layer=layer,locus=locus,
      interpretation="Conditional shared-variant posterior under H4")
  s$lead_shared<-as.character(v$snp[1]);s$lead_shared_pp<-v$SNP.PP.H4[1];s$credible_set_n<-sum(v$credible95)
  s$status<-"ok"
  su<-le8_run_susie(datasets$x,datasets$y,d,feature,lead_chr,start,end,locus);s$susie_status<-su$status
  r<-bind_rows(d|>transmute(layer,feature,locus,dataset="QTL",SNP,position=POS_x,p=P_x,beta=BETA_x,se=SE_x),
    d|>transmute(layer,feature,locus,dataset=Y,SNP,position=POS_y,p=P_y,beta=BETA_y,se=SE_y))
  list(summary=s,variants=v,regional=r,susie=su)
}
le8_c3_sets <- function(res,mr,layer) {
  e<-le8_same_locus_evidence(mr,res,layer);g<-unique(e$feature[e$eligible%in%TRUE])
  # Compatibility keys retained; names are NOT a claim of identified causation.
  list(Causal_Tier1=character(),Causal_Tier2plus=g,Causal_any=g)
}
le8_c3_additions <- function(res,mr,layer,outdir) {
  e<-le8_same_locus_evidence(mr,res,layer);rd<-le8_job_dir(outdir,"c3_coloc")
  write_raw_csv(e,"c3.same_locus_evidence.csv",rd)
  ss<-list.files(file.path(rd,"susie"),pattern="\\.rds$",full.names=TRUE)
  su<-bind_rows(lapply(ss,function(f){z<-readRDS(f);as_tibble(z$summary)}))
  pi<-bind_rows(lapply(ss,function(f){z<-readRDS(f);as_tibble(z$pip)}))
  write_raw_csv(su,"c3.susie_signal_pairs.csv",rd);write_raw_csv(pi,"c3.susie_trait_PIP.csv",rd)
  status<-if("susie_status"%in%names(res))res|>count(susie_status)else tibble(status="no loci")
  list(same_locus=e,susie_pairs=su,susie_PIP=pi,susie_status=status)
}
