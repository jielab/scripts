# C3: Colocalization — locus-specific QTL/outcome colocalization and fine-mapping evidence.
suppressPackageStartupMessages({
  source(file.path(Sys.getenv("LE8_FDIR", unset = file.path(Sys.getenv("DIRSCRIPT"), "f")), "comm.f.R"))
})
LE8_JOB <- "c3_coloc"
suppressPackageStartupMessages(pacman::p_load(coloc))
MAX_FEATURES <- as.integer(Sys.getenv("C3_MAX_FEATURES", unset = "200"))
MAX_LOCI_PER_FEATURE <- as.integer(Sys.getenv("C3_MAX_LOCI_PER_FEATURE", unset = "3"))
WINDOW_BP <- as.numeric(Sys.getenv("COLOC_WINDOW_BP", unset = "500000"))
MIN_SNPS <- as.integer(Sys.getenv("C3_MIN_NSNP", unset = "50"))
H4_STRONG <- as.numeric(Sys.getenv("C3_H4", unset = "0.70"))

migrate_c3_figure_names <- function(outdir) {
  for (i in 1:4) {
    old <- file.path(outdir, sprintf("c3.Fig%02d.%s.png", i,
      c("posterior_lollipop","posterior_hypotheses","regional_top_loci","coloc_network")[[i]]))
    new <- file.path(outdir, sprintf("c3.Fig%d.%s.png", i,
      c("posterior_lollipop","posterior_hypotheses","regional_top_loci","coloc_network")[[i]]))
    if (file.exists(old)) {
      if (file.exists(new)) unlink(old, force=TRUE)
      else if (!file.rename(old,new)) {
        if (file.copy(old,new,overwrite=TRUE)) unlink(old,force=TRUE)
      }
    }
  }
}

run_gpu_coloc_step <- function(layer, rawdir, manifest, ygfile, outtype) {
  gpudir <- file.path(rawdir, "gpu_coloc")
  dir.create(gpudir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(gpudir, "gpu_coloc.results.tsv")
  if (cache_valid(result_file)) {
    cache_message(paste0("GPU-coloc/", layer), result_file)
    return(invisible(0L))
  }
  if (!truthy(Sys.getenv("RUN_GPU_COLOC", unset = "TRUE")) || !nrow(manifest)) return(invisible(0L))
  manifest_file <- file.path(rawdir, "qtl_cad_manifest.tsv")
  if (!file.exists(manifest_file)) write_raw_tsv(manifest, basename(manifest_file), rawdir)
  sh <- file.path(Sys.getenv("LE8_FDIR"), "c3_coloc_GPU.sh")
  status <- system2("bash", c(sh, "--qtl-manifest", manifest_file, "--cad-gwas", ygfile,
    "--outdir", gpudir, "--outcome-type", outtype, "--H4", as.character(H4_STRONG)))
  if (status != 0) warning("GPU-coloc runner returned status ", status)
  invisible(status)
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
                                 lead_shared=NA_character_,credible_set_n=NA_integer_,message="Insufficient aligned SNPs"),
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
  fit <- tryCatch({
    fit_value <- NULL
    invisible(utils::capture.output(
      fit_value <- coloc::coloc.abf(dx,dy),
      type = "output"
    ))
    fit_value
  },error=function(e)e)
  if (inherits(fit,"error")) { empty$summary$message <- conditionMessage(fit); return(empty) }
  s <- fit$summary; vr <- as_tibble(fit$results)
  ppcol <- grep("SNP.PP.H4",names(vr),value=TRUE)[1]
  if (is.na(ppcol)) vr$SNP.PP.H4 <- NA_real_ else names(vr)[names(vr)==ppcol] <- "SNP.PP.H4"
  vr <- vr |> arrange(desc(SNP.PP.H4)) |> mutate(cum_pp=cumsum(coalesce(SNP.PP.H4,0)),credible95=cum_pp<=.95|lag(cum_pp,default=0)<.95)
  lead_shared <- if(nrow(vr) && "snp" %in% names(vr)) as.character(vr$snp[[1]]) else NA_character_
  summ <- tibble(layer,feature,locus,chr=lead_chr,start,end,lead_pos,status="ok",n_snps=nrow(d),
                 PP.H0=as.numeric(s[["PP.H0.abf"]]),PP.H1=as.numeric(s[["PP.H1.abf"]]),PP.H2=as.numeric(s[["PP.H2.abf"]]),
                 PP.H3=as.numeric(s[["PP.H3.abf"]]),PP.H4=as.numeric(s[["PP.H4.abf"]]),
                 lead_shared=lead_shared,credible_set_n=sum(vr$credible95,na.rm=TRUE),message="")
  reg <- bind_rows(d |> transmute(layer,feature,locus,dataset="QTL",SNP,position=POS_x,p=P_x,beta=BETA_x,se=SE_x),
                   d |> transmute(layer,feature,locus,dataset=Y,SNP,position=POS_y,p=P_y,beta=BETA_y,se=SE_y))
  list(summary=summ,variants=vr |> mutate(layer,feature,locus),regional=reg)
}

plot_coloc_results <- function(res, regional, variants, layer, outdir) {
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

run_c3_layer <- function(layer=c("protein","metabolite")) {
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir);rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE);cache<-file.path(rawdir,"c3.res.rds")
  migrate_c3_figure_names(outdir)
  if(cache_valid(cache)){
    cache_message(paste0("C3/",layer),cache)
    cached <- readRDS(cache); ygfile <- get_y_gwas_file(Y,TRUE); outtype <- infer_outcome_type(Y)
    run_gpu_coloc_step(layer,rawdir,cached$manifest %||% tibble(),ygfile,outtype)
    return(cached)
  }
  c1f<-file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds");c2f<-file.path(le8_job_dir(outdir,"c2_cause"),"c2.res.rds")
  if(!file.exists(c1f))stop("Run C1 first.",call.=FALSE);c1<-readRDS(c1f);assoc<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble()
  c2<-if(file.exists(c2f))readRDS(c2f) else list();mr<-c2$MR%||%tibble()
  candidates<-unique(c(mr|>filter(is.finite(pval))|>arrange(pval)|>pull(exposure),assoc|>arrange(p.value)|>pull(term)))|>head(MAX_FEATURES)
  base<-if(layer=="protein")dir.X else dir.met.gwas;ann<-layer_annotation(layer,candidates);ygfile<-get_y_gwas_file(Y,TRUE)
  outtype<-infer_outcome_type(Y);sfrac<-if(outtype=="cc")get_case_fraction(Y) else NA_real_
  if(outtype=="cc"&&(!is.finite(sfrac)||sfrac<=0||sfrac>=1))stop("C3 could not determine the case fraction. Set COLOC_CASE_FRAC.",call.=FALSE)
  rows<-list();vrows<-list();rrows<-list();manifest<-list();k<-0L
  for(feature in candidates){
    q<-read_qtl_instruments(feature,base,layer,ann);iv<-q$instruments;qf<-if(!is.na(q$files$full)) q$files$full else q$files$cis
    if(!nrow(iv)||is.na(qf)||!file.exists(qf))next
    leads<-select_nonoverlapping_leads(iv,MAX_LOCI_PER_FEATURE,WINDOW_BP)
    for(i in seq_len(nrow(leads))){k<-k+1;z<-coloc_one_locus(feature,layer,qf,leads$CHR[i],leads$POS[i],ygfile,outtype,sfrac);rows[[k]]<-z$summary;vrows[[k]]<-z$variants;rrows[[k]]<-z$regional
      manifest[[k]]<-tibble(omics=layer,trait=feature,file=qf,type="quant",region=z$summary$locus)}
  }
  res<-bind_rows(rows);variants<-bind_rows(vrows);regional<-bind_rows(rrows);mani<-bind_rows(manifest)
  if(!nrow(res))stop("C3/",layer,": no QTL locus could be constructed.",call.=FALSE)
  res<-res|>mutate(PP4_rank_fraction=ifelse(status=="ok"&sum(status=="ok")>0,rank(-PP.H4,ties.method="min",na.last="keep")/sum(status=="ok"),NA_real_),tier=case_when(status!="ok"~"Not tested",PP.H4>=.8~"Tier 1",PP.H4>=H4_STRONG~"Tier 2",PP.H4>=.5~"Tier 3",TRUE~"Tier 4"))
  write_raw_csv(res,"c3.coloc_summary.csv",rawdir);write_raw_csv(variants,"c3.variant_posteriors.csv",rawdir);write_raw_csv(regional,"c3.regional_rows.csv",rawdir);write_raw_tsv(mani,"qtl_cad_manifest.tsv",rawdir)
  plot_coloc_results(res,regional,variants,layer,outdir)
  run_gpu_coloc_step(layer,rawdir,mani,ygfile,outtype)
  lists<-list(Causal_Tier1=res|>filter(tier=="Tier 1")|>pull(feature)|>unique(),Causal_Tier2plus=res|>filter(tier%in%c("Tier 1","Tier 2"))|>pull(feature)|>unique(),Causal_any=res|>filter(PP.H4>=.5)|>pull(feature)|>unique())
  saveRDS(lists,file.path(rawdir,paste0("c3.causal_",layer,"_lists.rds")),compress="xz")
  if(layer=="protein")saveRDS(lists,file.path(rawdir,"c3.causal_protein_lists.rds"),compress="xz")
  out<-list(meta=module_meta(layer,extra=list(outcome_type=outtype,case_fraction=sfrac)),summary=res,variants=variants,regional=regional,manifest=mani,causal_lists=lists)
  saveRDS(out,cache,compress="xz");write_xlsx2(list(coloc_summary=res,variant_posteriors=variants,regional=regional,GPU_manifest=mani,causal_sets=stack(lists)),"c3.out.xlsx");finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO)run_c3_layer("protein")
if(met_DO)run_c3_layer("metabolite")
