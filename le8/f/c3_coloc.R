# C3: Colocalization — locus-specific QTL/outcome colocalization and fine-mapping evidence.
suppressPackageStartupMessages({
  fdir <- Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f"))
  source(file.path(fdir, "comm.f.R"))
  source(file.path(fdir, "c3_pgs.R"))
})
LE8_JOB <- "c3_coloc"
C3_CODE_VERSION <- "2026-09-02.pgs_triangulation1"
suppressPackageStartupMessages(pacman::p_load(coloc))
MAX_FEATURES <- as.integer(Sys.getenv("C3_MAX_FEATURES", unset = "200"))
MAX_LOCI_PER_FEATURE <- as.integer(Sys.getenv("C3_MAX_LOCI_PER_FEATURE", unset = "3"))
WINDOW_BP <- as.numeric(Sys.getenv("COLOC_WINDOW_BP", unset = "500000"))
MIN_SNPS <- as.integer(Sys.getenv("C3_MIN_NSNP", unset = "50"))
H4_STRONG <- as.numeric(Sys.getenv("C3_H4", unset = "0.70"))
C3_P12 <- c(conservative=1e-6,default=1e-5,liberal=1e-4)

run_gpu_coloc_step <- function(layer, rawdir, manifest, ygfile, outtype) {
  gpudir <- file.path(rawdir, "gpu_coloc")
  dir.create(gpudir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(gpudir, "gpu_coloc.results.tsv")
  if (!truthy(Sys.getenv("RUN_GPU_COLOC", unset = "TRUE"))) return(invisible(0L))
  if (cache_valid(result_file)) {
    cache_message(paste0("GPU-coloc/", layer), result_file)
    return(invisible(0L))
  }
  if (!nrow(manifest)) return(invisible(0L))
  sh <- file.path(Sys.getenv("LE8_FDIR"), "c3_coloc_GPU.sh")
  manifest_file <- file.path(rawdir, "qtl_cad_manifest.tsv")
  write_raw_tsv(manifest, basename(manifest_file), rawdir)
  status <- system2("bash", c(sh, "--qtl-manifest", manifest_file, "--cad-gwas", ygfile,
    "--outdir", gpudir, "--outcome-type", outtype, "--H4", as.character(H4_STRONG)))
  if (status != 0) warning("GPU-coloc runner returned status ", status)
  invisible(status)
}

read_gpu_coloc_results <- function(rawdir){
  gpudir<-file.path(rawdir,"gpu_coloc");rf<-file.path(gpudir,"gpu_coloc.results.tsv");sf<-file.path(gpudir,"signal_preparation_status.tsv")
  ans<-list(results=tibble(),status=tibble())
  if(!truthy(Sys.getenv("RUN_GPU_COLOC",unset="TRUE")))return(ans)
  if(file.exists(rf)&&file.size(rf)>0)ans$results<-tryCatch(as_tibble(data.table::fread(rf,showProgress=FALSE,check.names=FALSE)),error=function(e)tibble())
  if(file.exists(sf)&&file.size(sf)>0)ans$status<-tryCatch(as_tibble(data.table::fread(sf,showProgress=FALSE,check.names=FALSE)),error=function(e)tibble())
  d<-ans$results
  if(nrow(d)){
    pp<-names(d)[toupper(names(d))=="PP.H4"][1]
    if(is.na(pp))pp<-names(d)[str_detect(toupper(names(d)),"PP.*H4")][1]
    sigcols<-names(d)[vapply(d,function(x)any(str_detect(as.character(x),"QTL__"),na.rm=TRUE),logical(1))]
    if("trait"%in%names(d))d$feature<-as.character(d$trait) else if(length(sigcols)){
      sig<-as.character(d[[sigcols[1]]]);trait<-str_match(sig,"^QTL__[^_]+__(.*?)__chr")
      d$feature<-if(ncol(trait)>=2)trait[,2]else NA_character_
    } else d$feature<-NA_character_
    d$GPU_PP.H4<-if(length(pp)&&!is.na(pp))suppressWarnings(as.numeric(d[[pp]]))else NA_real_
    ans$results<-d
  }
  ans
}

plot_gpu_coloc_validation <- function(gpu,abf,outdir){
  d<-gpu$results;st<-gpu$status
  if(nrow(d)&&"GPU_PP.H4"%in%names(d)){
    if("region"%in%names(d)&&"locus"%in%names(abf)){
      z<-d|>filter(is.finite(GPU_PP.H4))|>group_by(feature,region)|>slice_max(GPU_PP.H4,n=1,with_ties=FALSE)|>ungroup()|>
        left_join(abf|>filter(status=="ok")|>transmute(feature,region=locus,ABF_PP.H4=PP.H4),by=c("feature","region"))
    }else{
      z<-d|>filter(is.finite(GPU_PP.H4))|>group_by(feature)|>slice_max(GPU_PP.H4,n=1,with_ties=FALSE)|>ungroup()|>
        left_join(abf|>filter(status=="ok")|>group_by(feature)|>summarise(ABF_PP.H4=max(PP.H4,na.rm=TRUE),.groups="drop"),by="feature")|>
        mutate(region=NA_character_)
    }
    z<-z|>arrange(desc(GPU_PP.H4))|>slice_head(n=40)|>
      mutate(label=ifelse(is.na(region)|!nzchar(region),feature,paste(feature,region,sep=" | ")),label=factor(label,levels=rev(label)))
    pa<-if(!nrow(z))blank_plot("a. GPU-coloc regional results")else ggplot(z,aes(GPU_PP.H4,label))+
      geom_vline(xintercept=H4_STRONG,linetype=2,color="grey55")+geom_segment(aes(x=0,xend=GPU_PP.H4,yend=label),color="grey75")+
      geom_point(aes(color=GPU_PP.H4>=H4_STRONG),size=2.2)+scale_color_manual(values=c(`TRUE`="#1B9E77",`FALSE`="grey65"),guide="none")+
      scale_x_continuous(limits=c(0,1),labels=label_percent())+
      labs(title="a. GPU-coloc regional marginal-signal evidence",subtitle="Exact manifest QTL-outcome pairs only",x="GPU-coloc PP(H4)",y=NULL)+theme_5c(8)
    pb<-if(!nrow(z)||!any(is.finite(z$ABF_PP.H4)))blank_plot("b. ABF versus GPU-coloc","No feature could be aligned across methods")else ggplot(z,aes(ABF_PP.H4,GPU_PP.H4))+
      geom_abline(slope=1,intercept=0,linetype=2,color="grey55")+geom_hline(yintercept=H4_STRONG,linetype=3,color="grey70")+geom_vline(xintercept=H4_STRONG,linetype=3,color="grey70")+
      geom_point(color="#6A3D9A",size=2)+ggrepel::geom_text_repel(aes(label=label),size=2.5,max.overlaps=15,seed=63)+
      coord_equal(xlim=c(0,1),ylim=c(0,1))+
      labs(title="b. coloc.abf versus GPU-coloc",subtitle="Same marginal locus and p12; disagreement flags input or prior-scale differences",x="coloc.abf PP(H4)",y="GPU-coloc PP(H4)")+theme_5c(9)
  } else {pa<-blank_plot("a. GPU-coloc regional results","GPU-coloc was disabled, unavailable, or produced no PP(H4) column");pb<-blank_plot("b. ABF versus GPU-coloc")}
  pc<-if(!nrow(st)||!"status"%in%names(st))blank_plot("c. GPU-coloc preparation audit","No status file")else st|>count(status)|>
    ggplot(aes(n,fct_reorder(status,n),fill=status))+geom_col()+geom_text(aes(label=n),hjust=-.1,fontface="bold")+scale_x_continuous(expand=expansion(mult=c(0,.2)))+
      labs(title="c. Signal-preparation audit",x="Loci",y=NULL,fill=NULL)+theme_5c(9)+theme(legend.position="none")
  save_plot((pa|pb)/pc+plot_layout(heights=c(1,.55)),"c3.Fig6.gpu_coloc_validation.png",15.5,11,outdir=outdir)
}

infer_outcome_type <- function(y) {
  explicit <- tolower(Sys.getenv("C3_OUTCOME_TYPE", unset = "auto")); if (explicit %in% c("cc","quant")) return(explicit)
  if (str_detect(tolower(y), "bmi|height|weight|bp$|ldl|hdl|hba1c|glucose|age")) "quant" else "cc"
}

select_nonoverlapping_leads <- function(iv, max_loci = MAX_LOCI_PER_FEATURE, window_bp = WINDOW_BP) {
  iv <- iv |> filter(is.finite(POS), !is.na(CHR), is.finite(P)) |> arrange(P)
  if (!nrow(iv)) return(iv)
  keep <- logical(nrow(iv)); chosen_chr <- character(); chosen_pos <- numeric()
  for (i in seq_len(nrow(iv))) {
    overlap <- any(chosen_chr == iv$CHR[[i]] & abs(chosen_pos - iv$POS[[i]]) <= 2 * window_bp)
    if (!overlap) {
      keep[[i]] <- TRUE; chosen_chr <- c(chosen_chr, iv$CHR[[i]]); chosen_pos <- c(chosen_pos, iv$POS[[i]])
      if (sum(keep) >= max_loci) break
    }
  }
  iv[keep, , drop = FALSE]
}

coloc_one_locus <- function(feature, layer, qtl_file, lead_chr, lead_pos, ygwas_file, outcome_type, case_frac) {
  start <- max(1, lead_pos - WINDOW_BP); end <- lead_pos + WINDOW_BP
  q <- read_sumstat_region(qtl_file, lead_chr, start, end)
  y <- read_sumstat_region(ygwas_file, lead_chr, start, end)
  d <- harmonize_sumstats(q, y)
  locus <- paste0("chr", lead_chr, ":", floor(start), "-", ceiling(end))
  empty <- list(summary = tibble(layer, feature, locus, chr=lead_chr,start,end,lead_pos,status="failed",n_snps=nrow(d),
                                 PP.H0=NA_real_,PP.H1=NA_real_,PP.H2=NA_real_,PP.H3=NA_real_,PP.H4=NA_real_,
                                 PP.H4_p12_conservative=NA_real_,PP.H4_p12_default=NA_real_,PP.H4_p12_liberal=NA_real_,
                                 PP.H4_robust_min=NA_real_,lead_shared=NA_character_,lead_shared_pp=NA_real_,
                                 credible_set_n=NA_integer_,message="Insufficient aligned SNPs"),
                variants=tibble(), regional=tibble())
  if (nrow(d) < MIN_SNPS) return(empty)
  # coloc requires allele frequencies. Use the matched dataset as a fallback, but never invent MAF=0.5.
  d <- d |> mutate(EAF_q = coalesce(EAF_x, EAF_y), EAF_y2 = coalesce(EAF_y, EAF_x)) |>
    filter(is.finite(EAF_q), EAF_q > 0, EAF_q < 1, is.finite(EAF_y2), EAF_y2 > 0, EAF_y2 < 1,
           is.finite(BETA_x), is.finite(SE_x), SE_x > 0, is.finite(BETA_y), is.finite(SE_y), SE_y > 0)
  empty$summary$n_snps <- nrow(d)
  if (nrow(d) < MIN_SNPS) { empty$summary$message <- "Insufficient aligned SNPs with usable allele frequencies"; return(empty) }
  mafq <- pmin(d$EAF_q, 1-d$EAF_q); mafy <- pmin(d$EAF_y2, 1-d$EAF_y2)
  nx <- median(d$N_x,na.rm=TRUE); ny <- median(d$N_y,na.rm=TRUE)
  if(!is.finite(nx)) nx <- as.numeric(Sys.getenv("C2_SUMSTAT_N",unset="100000")); if(!is.finite(ny)) ny <- as.numeric(Sys.getenv("C2_SUMSTAT_N",unset="100000"))
  dx <- list(beta=d$BETA_x,varbeta=d$SE_x^2,snp=d$SNP,position=d$POS_x,MAF=mafq,N=nx,type="quant")
  dy <- list(beta=d$BETA_y,varbeta=d$SE_y^2,snp=d$SNP,position=d$POS_y,MAF=mafy,N=ny,type=outcome_type)
  if (outcome_type=="cc") dy$s <- case_frac
  # coloc.abf() prints the five hypothesis posteriors for every locus. Capture
  # that routine stdout to keep the console readable; warnings and errors are
  # intentionally left visible.
  fits <- lapply(C3_P12,function(p12)tryCatch({
    fit_value <- NULL
    invisible(utils::capture.output(
      fit_value <- coloc::coloc.abf(dx,dy,p1=1e-4,p2=1e-4,p12=p12),type="output"))
    fit_value
  },error=function(e)e))
  fit <- fits[["default"]]
  if (inherits(fit,"error")) { empty$summary$message <- conditionMessage(fit); return(empty) }
  s <- fit$summary; vr <- as_tibble(fit$results)
  ppcol <- grep("SNP.PP.H4",names(vr),value=TRUE)[1]
  if (is.na(ppcol)) vr$SNP.PP.H4 <- NA_real_ else names(vr)[names(vr)==ppcol] <- "SNP.PP.H4"
  vr <- vr |> arrange(desc(SNP.PP.H4)) |> mutate(cum_pp=cumsum(coalesce(SNP.PP.H4,0)),credible95=cum_pp<=.95|lag(cum_pp,default=0)<.95)
  lead_shared <- if(nrow(vr) && "snp" %in% names(vr)) as.character(vr$snp[[1]]) else NA_character_
  pp4_sens<-vapply(fits,function(ff)if(inherits(ff,"error"))NA_real_ else as.numeric(ff$summary[["PP.H4.abf"]]),numeric(1))
  summ <- tibble(layer,feature,locus,chr=lead_chr,start,end,lead_pos,status="ok",n_snps=nrow(d),
                 PP.H0=as.numeric(s[["PP.H0.abf"]]),PP.H1=as.numeric(s[["PP.H1.abf"]]),PP.H2=as.numeric(s[["PP.H2.abf"]]),
                 PP.H3=as.numeric(s[["PP.H3.abf"]]),PP.H4=as.numeric(s[["PP.H4.abf"]]),
                 PP.H4_p12_conservative=pp4_sens[["conservative"]],PP.H4_p12_default=pp4_sens[["default"]],
                 PP.H4_p12_liberal=pp4_sens[["liberal"]],
                 PP.H4_robust_min=if(any(is.finite(pp4_sens)))min(pp4_sens[is.finite(pp4_sens)])else NA_real_,
                 lead_shared=lead_shared,lead_shared_pp=if(nrow(vr))vr$SNP.PP.H4[[1]]else NA_real_,
                 credible_set_n=sum(vr$credible95,na.rm=TRUE),message="")
  reg <- bind_rows(d |> transmute(layer,feature,locus,dataset="QTL",SNP,position=POS_x,p=P_x,beta=BETA_x,se=SE_x),
                   d |> transmute(layer,feature,locus,dataset=Y,SNP,position=POS_y,p=P_y,beta=BETA_y,se=SE_y))
  list(summary=summ,variants=vr |> mutate(layer,feature,locus),regional=reg)
}

# Retained for reproducibility of the original four-panel output.  The
# publication function below is the active implementation.
plot_coloc_results_legacy <- function(res, regional, variants, layer, outdir) {
  ok <- res |> filter(status=="ok") |> arrange(desc(PP.H4))
  if (!nrow(ok)) {
    blank_files <- c("c3.Fig1.posterior_lollipop.png", "c3.Fig2.posterior_hypotheses.png",
                     "c3.Fig3.regional_top_loci.png", "c3.Fig4.coloc_network.png")
    for (i in seq_along(blank_files)) save_plot(blank_plot(paste0("C3 Figure ", i),
      "No locus passed the minimum aligned-SNP requirement"), blank_files[[i]], 9, 5, outdir = outdir)
    return(invisible(NULL))
  }
  top <- ok |> slice_head(n=50) |> mutate(label=paste(feature,locus,sep=" | "),label=factor(label,levels=rev(label)),tier=case_when(PP.H4>=.8~"Strong",PP.H4>=H4_STRONG~"Probable",PP.H4>=.5~"Suggestive",TRUE~"Weak"))
  p1 <- ggplot(top,aes(PP.H4,label))+geom_segment(aes(x=0,xend=PP.H4,y=label,yend=label,color=tier),linewidth=.75)+geom_point(aes(color=tier,size=n_snps),alpha=.9)+
    geom_vline(xintercept=H4_STRONG,linetype=2,color="grey45")+scale_color_manual(values=c(Strong="#1B9E77",Probable="#66A61E",Suggestive="#E6AB02",Weak="grey70"))+
    scale_x_continuous(limits=c(0,1),labels=label_percent())+labs(title="Colocalization posterior probability",subtitle="Each row is a QTL locus tested against the disease GWAS",x="PP(H4): shared causal signal",y=NULL,color=NULL,size="Aligned SNPs")+theme_5c(9)+theme(legend.position="bottom")
  save_plot(p1,"c3.Fig1.posterior_lollipop.png",10,10,outdir=outdir)

  post <- top |> select(label,starts_with("PP.H")) |> pivot_longer(starts_with("PP.H"),names_to="hypothesis",values_to="posterior")|>
    mutate(hypothesis=factor(hypothesis,levels=paste0("PP.H",0:4)))
  p2 <- ggplot(post,aes(label,posterior,fill=hypothesis))+geom_col(width=.78)+coord_flip()+scale_fill_brewer(palette="Set2")+
    scale_y_continuous(labels=label_percent())+labs(title="Competing colocalization hypotheses",x=NULL,y="Posterior probability",fill="Hypothesis")+theme_5c(9)+theme(legend.position="bottom")
  save_plot(p2,"c3.Fig2.posterior_hypotheses.png",10,10,outdir=outdir)

  loci <- ok |> slice_head(n=6) |> pull(locus); rr <- regional |> filter(locus%in%loci,is.finite(p),p>0)|>mutate(y=-log10(pmax(p,1e-300)),position_mb=position/1e6)
  p3 <- ggplot(rr,aes(position_mb,y,color=dataset))+geom_point(size=.8,alpha=.72)+geom_vline(data=ok|>filter(locus%in%loci),aes(xintercept=lead_pos/1e6),inherit.aes=FALSE,linetype=2,color="grey55")+
    facet_grid(dataset~feature+locus,scales="free_x",space="free_x")+scale_color_manual(values=c(QTL="#4C78A8",setNames("#E45756",Y)))+
    labs(title="Regional QTL–disease comparison at top loci",x="Position (Mb)",y=expression(-log[10](P)),color=NULL)+theme_5c(8)+theme(legend.position="bottom",axis.text.x=element_text(angle=45,hjust=1))
  save_plot(p3,"c3.Fig3.regional_top_loci.png",17,7.8,outdir=outdir)

  # Evidence network: feature nodes on the outside, genomic loci inside; edge width is PP4.
  net <- ok |> filter(PP.H4>=.5) |> distinct(feature,locus,PP.H4,chr,lead_pos)
  if (!nrow(net)) p4 <- blank_plot("Colocalization evidence network","No locus had PP(H4) >= 0.5") else {
    fnode <- tibble(name=unique(net$feature),type="feature",angle=seq(0,2*pi,length.out=n_distinct(net$feature)+1)[-(n_distinct(net$feature)+1)],r=1)
    lnode <- tibble(name=unique(net$locus),type="locus",angle=seq(0,2*pi,length.out=n_distinct(net$locus)+1)[-(n_distinct(net$locus)+1)]+pi/10,r=.48)
    nodes <- bind_rows(fnode,lnode)|>mutate(x=r*cos(angle),y=r*sin(angle))
    edges <- net |> left_join(nodes|>select(feature=name,x1=x,y1=y),by="feature") |> left_join(nodes|>select(locus=name,x2=x,y2=y),by="locus")
    p4 <- ggplot()+geom_curve(data=edges,aes(x=x1,y=y1,xend=x2,yend=y2,linewidth=PP.H4,color=PP.H4),curvature=.15,alpha=.65)+
      geom_point(data=nodes,aes(x,y,shape=type),size=3,fill="white")+ggrepel::geom_text_repel(data=nodes,aes(x,y,label=name),size=2.5,max.overlaps=40,seed=4)+
      scale_linewidth(range=c(.4,2.5),guide="none")+scale_color_gradient(low="#FEE8C8",high="#B30000",limits=c(.5,1),name="PP(H4)")+coord_equal()+theme_void()+labs(title="Shared-locus evidence network")+theme(plot.title=element_text(face="bold"),legend.position="bottom")
  }
  save_plot(p4,"c3.Fig4.coloc_network.png",11,9,outdir=outdir)
}

# Publication figures cover evidence, posterior, regional, and prior diagnostics.
plot_coloc_results <- function(res,regional,variants,layer,outdir){
  ok<-res|>filter(status=="ok")|>arrange(desc(PP.H4))
  files<-c("c3.Fig1.evidence_triage.png","c3.Fig2.posterior_diagnostics.png","c3.Fig3.regional_top_loci.png",
           "c3.Fig4.credible_sets.png","c3.Fig5.prior_sensitivity.png")
  if(!nrow(ok)){
    for(i in seq_along(files))save_plot(blank_plot(paste0("C3 Figure ",i),"No locus passed the aligned-SNP requirement"),files[[i]],9,5,outdir=outdir)
    return(invisible(NULL))
  }
  for(nm in c("PP.H4_p12_conservative","PP.H4_p12_default","PP.H4_p12_liberal","PP.H4_robust_min","lead_shared_pp"))
    if(!nm%in%names(ok))ok[[nm]]<-if(nm=="PP.H4_p12_default")ok$PP.H4 else NA_real_
  ok<-ok|>arrange(desc(coalesce(PP.H4_robust_min,PP.H4)),desc(PP.H4))
  top<-ok|>slice_head(n=30)|>mutate(label=paste(feature,str_replace(locus,"^chr","chr"),sep=" | "),
    label=factor(label,levels=rev(label)),robust=coalesce(PP.H4_robust_min,PP.H4)>=H4_STRONG,
    lo=pmin(PP.H4_p12_conservative,PP.H4_p12_default,PP.H4_p12_liberal,na.rm=TRUE),
    hi=pmax(PP.H4_p12_conservative,PP.H4_p12_default,PP.H4_p12_liberal,na.rm=TRUE))
  p1<-ggplot(top,aes(PP.H4,label))+geom_vline(xintercept=H4_STRONG,linetype=2,color="grey50")+
    geom_errorbarh(aes(xmin=lo,xmax=hi,color=robust),height=.10,linewidth=.8)+geom_point(aes(size=lead_shared_pp,color=robust),alpha=.9)+
    scale_color_manual(values=c(`TRUE`="#1B9E77",`FALSE`="#D95F02"),labels=c(`TRUE`="Robust to p12",`FALSE`="Prior-sensitive"))+
    scale_x_continuous(limits=c(0,1),labels=label_percent())+scale_size_continuous(range=c(1.8,5),name="Lead SNP PP")+
    labs(title="Colocalization evidence triage",subtitle="Point: default PP(H4); interval: p12 = 1e-6 to 1e-4",
      x="Posterior probability of a shared signal",y=NULL,color=NULL)+theme_5c(8)+theme(legend.position="bottom")
  save_plot(p1,files[[1]],11.5,10,outdir=outdir)

  t25<-top|>slice_head(n=25);post<-t25|>select(label,PP.H0,PP.H1,PP.H2,PP.H3,PP.H4)|>
    pivot_longer(starts_with("PP.H"),names_to="hypothesis",values_to="posterior")|>
    mutate(hypothesis=factor(hypothesis,levels=paste0("PP.H",0:4)))
  pa<-ggplot(post,aes(hypothesis,label,fill=posterior))+geom_tile(color="white",linewidth=.25)+
    scale_fill_viridis_c(limits=c(0,1),labels=label_percent())+
    labs(title="a. Competing hypotheses",x=NULL,y=NULL,fill="Posterior")+theme_5c(8)
  pb<-ggplot(t25,aes(credible_set_n,n_snps,color=PP.H4,size=lead_shared_pp))+geom_point(alpha=.82)+
    ggrepel::geom_text_repel(aes(label=as.character(label)),size=2.3,seed=52,max.overlaps=12)+
    scale_x_log10()+scale_y_log10()+scale_color_viridis_c(limits=c(0,1),name="PP(H4)")+
    labs(title="b. Shared-signal posterior resolution",subtitle="Small sets and high lead-SNP posterior are preferred, conditional on the coloc single-signal model",
      x="95% shared-signal posterior-set size",y="Aligned SNPs",size="Lead SNP PP")+theme_5c(8)
  save_plot(pa|pb,files[[2]],16,9.5,outdir=outdir)

  loc<-ok|>slice_head(n=4)|>select(feature,locus,lead_pos);rr<-regional|>semi_join(loc,by=c("feature","locus"))|>
    filter(is.finite(p),p>0)|>mutate(y=-log10(pmax(p,1e-300)),position_mb=position/1e6,panel=paste(feature,locus,sep=" | "))
  p3<-if(!nrow(rr))blank_plot("Regional QTL–disease comparisons")else ggplot(rr,aes(position_mb,y,color=dataset))+
    geom_line(linewidth=.35,alpha=.5)+geom_point(size=.85,alpha=.72)+
    geom_vline(data=loc|>mutate(panel=paste(feature,locus,sep=" | ")),aes(xintercept=lead_pos/1e6),inherit.aes=FALSE,linetype=2,color="grey45")+
    facet_wrap(~panel,ncol=2,scales="free")+scale_color_manual(values=c(QTL="#4C78A8",setNames("#E45756",Y)))+
    labs(title="Regional QTL–disease comparisons",subtitle="Four strongest loci; dashed line is the selected QTL lead",
      x="Position (Mb)",y=expression(-log[10](P)),color=NULL)+theme_5c(8)+theme(legend.position="bottom")
  save_plot(p3,files[[3]],14.5,10,outdir=outdir)

  vv<-variants|>semi_join(loc,by=c("feature","locus"))
  if(nrow(vv)&&all(c("SNP.PP.H4","snp")%in%names(vv))){
    vv<-vv|>group_by(feature,locus)|>arrange(desc(SNP.PP.H4),.by_group=TRUE)|>slice_head(n=40)|>ungroup()|>
      mutate(rank=ave(-SNP.PP.H4,interaction(feature,locus),FUN=rank),panel=paste(feature,locus,sep=" | "),
        credible95=coalesce(credible95,FALSE),posterior_plot=pmax(coalesce(SNP.PP.H4,0),1e-10))
    p4a<-ggplot(vv,aes(rank,posterior_plot,color=credible95))+geom_segment(aes(xend=rank,y=1e-10,yend=posterior_plot),linewidth=.40)+geom_point(size=1.55)+
      facet_wrap(~panel,ncol=2,scales="free_x")+scale_y_log10(limits=c(1e-10,1),breaks=10^c(-10,-8,-6,-4,-2,0))+
      scale_color_manual(values=c(`TRUE`="#1B9E77",`FALSE`="grey70"))+
      labs(title="a. Variant-level shared-signal posterior (log scale)",
        subtitle="A posterior near 1 with all alternatives near 0 is visible rather than appearing as a blank panel",
        x="Variant rank within locus",y="SNP posterior under H4",color="95% credible set")+theme_5c(8)+theme(legend.position="bottom")
  }else p4a<-blank_plot("a. Variant-level shared-signal posterior","No variant posterior was available")
  cs<-ok|>filter(is.finite(credible_set_n),credible_set_n>0,is.finite(PP.H4))|>
    mutate(robust=coalesce(PP.H4_robust_min,PP.H4)>=H4_STRONG,
      degenerate=credible_set_n==1&coalesce(lead_shared_pp,0)>.999,
      label=ifelse(feature%in%loc$feature|degenerate&PP.H4>=H4_STRONG,paste(feature,locus,sep=" | "),NA_character_))
  p4b<-if(!nrow(cs))blank_plot("b. Credible-set resolution","No credible-set summary was available")else
    ggplot(cs,aes(credible_set_n,coalesce(PP.H4_robust_min,PP.H4),color=degenerate,size=n_snps))+
      geom_hline(yintercept=H4_STRONG,linetype=2,color="grey55")+geom_point(alpha=.62)+
      ggrepel::geom_text_repel(aes(label=label),size=2.2,max.overlaps=15,seed=72)+scale_x_log10()+
      scale_color_manual(values=c(`TRUE`="#D7301F",`FALSE`="#3F78A8"),labels=c(`TRUE`="1-SNP posterior concentration",`FALSE`="Other"))+
      scale_y_continuous(limits=c(0,1),labels=label_percent())+
      labs(title="b. Resolution across every tested locus",x="95% credible-set size (log scale)",y="Robust minimum PP(H4)",color=NULL,size="Aligned SNPs")+theme_5c(8)+theme(legend.position="bottom")
  p4c<-if(!nrow(cs))blank_plot("c. Credible-set size distribution")else
    ggplot(cs,aes(credible_set_n,color=robust))+stat_ecdf(linewidth=.9)+scale_x_log10()+
      scale_color_manual(values=c(`TRUE`="#1B9E77",`FALSE`="grey55"),labels=c(`TRUE`="Robust shared signal",`FALSE`="Not robust"))+
      labs(title="c. Shared-signal resolution audit",
        subtitle=paste0(sum(cs$degenerate)," / ",nrow(cs)," loci have a one-SNP set with lead PP > 0.999; verify LD and harmonization"),
        x="95% credible-set size (log scale)",y="Cumulative fraction of loci",color=NULL)+theme_5c(8)+theme(legend.position="bottom")
  p4<-p4a/plot_spacer()/(p4b|p4c)+plot_layout(heights=c(1.25,.07,1))
  save_plot(p4,files[[4]],17,10,outdir=outdir)

  sens<-top|>mutate(highlight=row_number()<=10)|>select(label,highlight,conservative=PP.H4_p12_conservative,default=PP.H4_p12_default,liberal=PP.H4_p12_liberal)|>
    pivot_longer(c(conservative,default,liberal),names_to="prior",values_to="PP4")|>mutate(prior=factor(prior,levels=c("conservative","default","liberal"),labels=c("p12=1e-6","p12=1e-5","p12=1e-4")))
  p5<-ggplot(sens,aes(prior,PP4,group=label,color=highlight))+geom_hline(yintercept=H4_STRONG,linetype=2,color="grey55")+
    geom_line(alpha=.48)+geom_point(size=1.7)+scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="grey72"),guide="none")+
    scale_y_continuous(limits=c(0,1),labels=label_percent())+
    labs(title="Shared-signal prior sensitivity",subtitle="Highlighted lines are the ten strongest default-prior loci",
      x="Prior probability that one variant affects both traits",y="PP(H4)")+theme_5c(10)
  save_plot(p5,files[[5]],11.5,7.5,outdir=outdir)
}

credible_set_audit <- function(res,variants=tibble()){
  ok<-res|>filter(status=="ok")
  locus<-if(!nrow(variants)||!all(c("feature","locus","SNP.PP.H4")%in%names(variants)))tibble()else
    variants|>group_by(feature,locus)|>arrange(desc(SNP.PP.H4),.by_group=TRUE)|>
      summarise(variant_rows=n(),posterior_sum=sum(SNP.PP.H4,na.rm=TRUE),lead_pp=max(SNP.PP.H4,na.rm=TRUE),
        second_pp=ifelse(n()>1,nth(SNP.PP.H4,2),NA_real_),zero_pp_fraction=mean(coalesce(SNP.PP.H4,0)<=1e-12),.groups="drop")
  by_locus<-ok|>select(feature,locus,n_snps,PP.H4,PP.H4_robust_min,lead_shared_pp,credible_set_n)|>
    left_join(locus,by=c("feature","locus"))|>
    mutate(robust=coalesce(PP.H4_robust_min,PP.H4)>=H4_STRONG,
      one_snp_concentration=credible_set_n==1&coalesce(lead_shared_pp,lead_pp,0)>.999,
      audit_flag=case_when(one_snp_concentration~"verify LD/harmonization: one-SNP posterior concentration",
        !is.finite(credible_set_n)|credible_set_n<1~"missing credible set",TRUE~"ok"))
  overall<-tibble(metric=c("tested loci","robust shared-signal loci","median credible-set size","one-SNP posterior concentration","missing credible set"),
    value=c(nrow(by_locus),sum(by_locus$robust,na.rm=TRUE),median(by_locus$credible_set_n,na.rm=TRUE),
      sum(by_locus$one_snp_concentration,na.rm=TRUE),sum(!is.finite(by_locus$credible_set_n)|by_locus$credible_set_n<1)))
  list(overall=overall,by_locus=by_locus)
}

run_c3_layer <- function(layer=c("protein","metabolite")) {
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE);cache<-file.path(rawdir,"c3.res.rds")
  if(cache_valid(cache)){
    old<-tryCatch(readRDS(cache),error=function(e)NULL)
    if(!is.null(old)&&all(c("summary","regional","variants")%in%names(old))){
      message("C3/",layer,": reuse locus results and regenerate ",C3_CODE_VERSION," figures")
      plot_coloc_results(old$summary,old$regional,old$variants,layer,outdir)
      gpu<-old$GPU_coloc%||%read_gpu_coloc_results(rawdir);plot_gpu_coloc_validation(gpu,old$summary,outdir)
      aud<-credible_set_audit(old$summary,old$variants)
      tri<-read_c3_pgs_integration(layer,outdir,old$summary);plot_c3_pgs_integration(tri,outdir)
      write_raw_csv(aud$overall,"c3.credible_set_audit.csv",rawdir);write_raw_csv(aud$by_locus,"c3.credible_set_by_locus.csv",rawdir)
      write_raw_csv(tri,"c3.pgs_observed_coloc_triangulation.csv",rawdir)
      old$credible_set_audit<-aud;old$pgs_triangulation<-tri;old$meta$code_version<-C3_CODE_VERSION
      lists<-old$causal_lists%||%list();gpu<-old$GPU_coloc%||%list(results=tibble(),status=tibble())
      write_xlsx2(list(coloc_summary=old$summary,credible_set_audit=aud$overall,credible_set_by_locus=aud$by_locus,
        variant_posteriors=old$variants,regional=old$regional,GPU_results=gpu$results%||%tibble(),GPU_status=gpu$status%||%tibble(),
        GPU_manifest=old$manifest%||%tibble(),causal_sets=if(length(lists))stack(lists)else tibble(),
        pgs_observed_coloc=tri),"c3.out.xlsx")
      saveRDS(old,cache,compress="xz");finalize_outputs(LE8_JOB,outdir);return(old)
    }
  }
  base0<-if(layer=="protein")dir.X else dir.met.gwas;ygfile0<-get_y_gwas_file(Y,TRUE)
  c1f<-file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds");c2f<-file.path(le8_job_dir(outdir,"c2_cause"),"c2.res.rds")
  if(!file.exists(c1f))stop("Run C1 first.",call.=FALSE);c1<-readRDS(c1f);assoc<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble()
  c2<-if(file.exists(c2f))tryCatch(readRDS(c2f),error=function(e)e) else list()
  if(inherits(c2,"condition")){
    warning("C3/",layer,": optional C2 result is unreadable; candidates will be ranked from C1 only: ",
      conditionMessage(c2),call.=FALSE);c2<-list()
  }
  mr<-c2$MR%||%tibble()
  mr_ranked<-if(nrow(mr)&&all(c("exposure","pval")%in%names(mr)))mr|>
    filter(is.finite(pval))|>arrange(pval)|>pull(exposure)else character()
  assoc_ranked<-if(nrow(assoc)&&all(c("term","p.value")%in%names(assoc)))assoc|>
    filter(is.finite(p.value))|>arrange(p.value)|>pull(term)else character()
  candidates<-head(unique(c(mr_ranked,assoc_ranked)),MAX_FEATURES)
  base<-base0;ann<-layer_annotation(layer,candidates);ygfile<-ygfile0
  outtype<-infer_outcome_type(Y);sfrac<-if(outtype=="cc")get_case_fraction(Y) else NA_real_
  if(outtype=="cc"&&(!is.finite(sfrac)||sfrac<=0||sfrac>=1))stop("C3 could not determine the case fraction. Set COLOC_CASE_FRAC.",call.=FALSE)
  rows<-list();vrows<-list();rrows<-list();manifest<-list();k<-0L
  locus_cache_dir<-file.path(rawdir,"locus_cache");dir.create(locus_cache_dir,recursive=TRUE,showWarnings=FALSE)
  for(feature in candidates){
    q<-read_qtl_instruments(feature,base,layer,ann);iv<-q$instruments;qf<-if(!is.na(q$files$full)) q$files$full else q$files$cis
    if(!nrow(iv)||is.na(qf)||!file.exists(qf))next
    leads<-select_nonoverlapping_leads(iv,MAX_LOCI_PER_FEATURE,WINDOW_BP)
    for(i in seq_len(nrow(leads))){
      k<-k+1
      cache_name<-sprintf("%04d_%s_chr%s_%s.rds",k,gsub("[^A-Za-z0-9._-]","_",feature),leads$CHR[i],format(leads$POS[i],scientific=FALSE,trim=TRUE))
      locus_cache<-file.path(locus_cache_dir,cache_name);z<-read_stage_cache(locus_cache)
      if(is.null(z)){
        z<-coloc_one_locus(feature,layer,qf,leads$CHR[i],leads$POS[i],ygfile,outtype,sfrac)
        write_stage_cache(z,locus_cache)
      }
      rows[[k]]<-z$summary;vrows[[k]]<-z$variants;rrows[[k]]<-z$regional
      manifest[[k]]<-tibble(omics=layer,trait=feature,file=qf,type="quant",region=z$summary$locus)}
  }
  res<-bind_rows(rows);variants<-bind_rows(vrows);regional<-bind_rows(rrows);mani<-bind_rows(manifest)
  if(!nrow(res)){
    message("C3/",layer,": no QTL locus could be constructed; write auditable empty outputs and blank panels")
    res<-tibble(layer=character(),feature=character(),locus=character(),chr=character(),
      start=numeric(),end=numeric(),lead_pos=numeric(),status=character(),n_snps=integer(),
      PP.H0=numeric(),PP.H1=numeric(),PP.H2=numeric(),PP.H3=numeric(),PP.H4=numeric(),
      PP.H4_p12_conservative=numeric(),PP.H4_p12_default=numeric(),PP.H4_p12_liberal=numeric(),
      PP.H4_robust_min=numeric(),lead_shared=character(),lead_shared_pp=numeric(),
      credible_set_n=integer(),message=character(),PP4_rank_fraction=numeric(),tier=character())
    aud<-credible_set_audit(res,variants)
    write_raw_csv(res,"c3.coloc_summary.csv",rawdir);write_raw_csv(variants,"c3.variant_posteriors.csv",rawdir)
    write_raw_csv(regional,"c3.regional_rows.csv",rawdir);write_raw_tsv(mani,"qtl_cad_manifest.tsv",rawdir)
    write_raw_csv(aud$overall,"c3.credible_set_audit.csv",rawdir);write_raw_csv(aud$by_locus,"c3.credible_set_by_locus.csv",rawdir)
    plot_coloc_results(res,regional,variants,layer,outdir)
    gpu<-list(results=tibble(),status=tibble());plot_gpu_coloc_validation(gpu,res,outdir)
    tri<-read_c3_pgs_integration(layer,outdir,res);plot_c3_pgs_integration(tri,outdir)
    write_raw_csv(tri,"c3.pgs_observed_coloc_triangulation.csv",rawdir)
    lists<-list(Causal_Tier1=character(),Causal_Tier2plus=character(),Causal_any=character())
    out<-list(meta=module_meta(layer,extra=list(outcome_type=outtype,case_fraction=sfrac,
      code_version=C3_CODE_VERSION,status="no QTL locus constructed")),summary=res,variants=variants,
      regional=regional,manifest=mani,GPU_coloc=gpu,causal_lists=lists,credible_set_audit=aud,
      pgs_triangulation=tri)
    saveRDS(out,cache,compress="xz")
    write_xlsx2(list(coloc_summary=res,credible_set_audit=aud$overall,
      credible_set_by_locus=aud$by_locus,variant_posteriors=variants,regional=regional,
      GPU_results=gpu$results,GPU_status=gpu$status,GPU_manifest=mani,causal_sets=stack(lists),
      pgs_observed_coloc=tri),"c3.out.xlsx")
    finalize_outputs(LE8_JOB,outdir);return(out)
  }
  res<-res|>mutate(PP4_rank_fraction=ifelse(status=="ok"&sum(status=="ok")>0,rank(-PP.H4,ties.method="min",na.last="keep")/sum(status=="ok"),NA_real_),tier=case_when(status!="ok"~"Not tested",PP.H4>=.8~"Tier 1",PP.H4>=H4_STRONG~"Tier 2",PP.H4>=.5~"Tier 3",TRUE~"Tier 4"))
  aud<-credible_set_audit(res,variants)
  write_raw_csv(res,"c3.coloc_summary.csv",rawdir);write_raw_csv(variants,"c3.variant_posteriors.csv",rawdir);write_raw_csv(regional,"c3.regional_rows.csv",rawdir);write_raw_tsv(mani,"qtl_cad_manifest.tsv",rawdir)
  write_raw_csv(aud$overall,"c3.credible_set_audit.csv",rawdir);write_raw_csv(aud$by_locus,"c3.credible_set_by_locus.csv",rawdir)
  run_gpu_coloc_step(layer,rawdir,mani,ygfile,outtype)
  plot_coloc_results(res,regional,variants,layer,outdir)
  gpu<-read_gpu_coloc_results(rawdir);plot_gpu_coloc_validation(gpu,res,outdir)
  tri<-read_c3_pgs_integration(layer,outdir,res);plot_c3_pgs_integration(tri,outdir)
  write_raw_csv(tri,"c3.pgs_observed_coloc_triangulation.csv",rawdir)
  lists<-list(Causal_Tier1=res|>filter(tier=="Tier 1")|>pull(feature)|>unique(),Causal_Tier2plus=res|>filter(tier%in%c("Tier 1","Tier 2"))|>pull(feature)|>unique(),Causal_any=res|>filter(PP.H4>=.5)|>pull(feature)|>unique())
  saveRDS(lists,file.path(rawdir,paste0("c3.causal_",layer,"_lists.rds")),compress="xz")
  if(layer=="protein")saveRDS(lists,file.path(rawdir,"c3.causal_protein_lists.rds"),compress="xz")
  out<-list(meta=module_meta(layer,extra=list(outcome_type=outtype,case_fraction=sfrac,code_version=C3_CODE_VERSION)),
    summary=res,variants=variants,regional=regional,manifest=mani,GPU_coloc=gpu,causal_lists=lists,
    credible_set_audit=aud,pgs_triangulation=tri)
  saveRDS(out,cache,compress="xz");write_xlsx2(list(coloc_summary=res,credible_set_audit=aud$overall,
    credible_set_by_locus=aud$by_locus,variant_posteriors=variants,regional=regional,GPU_results=gpu$results,
    GPU_status=gpu$status,GPU_manifest=mani,causal_sets=stack(lists),pgs_observed_coloc=tri),
    "c3.out.xlsx");finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO)run_c3_layer("protein")
if(met_DO)run_c3_layer("metabolite")
