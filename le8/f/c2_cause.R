# C2: Causation — cis/local and trans/distal MR + DANDELION trans-regulatory
# disease-driver prioritization + instrument architecture.
#
# DANDELION implementation in this file follows the package's SNP-based design:
#   disease-associated SNP (gene1) -> trans-regulated protein/gene (gene2)
#                            -> gene-level disease association (MAGMA genes.out)
#
# IMPORTANT: MAGMA *gene-level* results (*.genes.out) are required by DANDELION.
# MAGMA gene-set output (*.gsa.out) is useful downstream for pathway enrichment,
# but it is not the p.wes / gene-level input required by med_gene().

suppressPackageStartupMessages({
  source(file.path(Sys.getenv("LE8_FDIR", unset=file.path(Sys.getenv("DIRSCRIPT"),"f")),"comm.f.R"))
})
LE8_JOB <- "c2_cause"
MAX_FEATURES <- as.integer(Sys.getenv("C2_MAX_FEATURES", unset="500"))
TOP_FOREST <- as.integer(Sys.getenv("C2_TOP_FOREST", unset="32"))
N_DEFAULT <- as.numeric(Sys.getenv("C2_SUMSTAT_N", unset="100000"))
RUN_DANDELION <- truthy(Sys.getenv("RUN_DANDELION", unset="TRUE"))
DANDELION_FDR <- as.numeric(Sys.getenv("C2_DANDELION_FDR", unset="0.10"))
DANDELION_CIS_BP <- as.numeric(Sys.getenv("C2_DANDELION_CIS_BP", unset="5000000"))
DANDELION_GWS <- as.numeric(Sys.getenv("C2_DANDELION_GWS", unset="5e-8"))
DANDELION_LEAD_BP <- as.numeric(Sys.getenv("C2_DANDELION_LEAD_BP", unset="5000000"))
DANDELION_MAX_SNPS <- as.integer(Sys.getenv("C2_DANDELION_MAX_SNPS", unset="100"))
DANDELION_MAX_GENE2 <- as.integer(Sys.getenv("C2_DANDELION_MAX_GENE2", unset="0")) # 0 = all available

# Publication typography for C2 only.
theme_5c <- function(base_size=12){
  theme_classic(base_size=base_size)+theme(
    plot.title=element_text(face="bold",size=base_size*1.14,hjust=0),
    plot.subtitle=element_text(face="bold",size=base_size*.94,color="grey30"),
    axis.title=element_text(face="bold",size=base_size*1.08),
    axis.text=element_text(face="bold",size=base_size,color="black"),
    legend.title=element_text(face="bold"),legend.text=element_text(face="bold"),
    strip.background=element_blank(),strip.text=element_text(face="bold"),
    panel.grid.major.y=element_line(color="grey91",linewidth=.25),panel.grid.minor=element_blank(),
    plot.margin=margin(9,13,9,13))
}
forest_theme <- function(base_size=10)theme_5c(base_size)+theme(panel.grid.major.y=element_blank())

run_mrlink2_step <- function(layer,rawdir,jobs,ygfile){
  linkdir<-file.path(rawdir,"mrlink2");dir.create(linkdir,recursive=TRUE,showWarnings=FALSE)
  jobs_file<-file.path(linkdir,"c2.mrlink2.jobs.tsv");if(!file.exists(jobs_file)&&nrow(jobs))write_raw_tsv(jobs,basename(jobs_file),linkdir)
  complete_file<-file.path(linkdir,"mrlink2.complete")
  if(cache_valid(complete_file))return(invisible(0L))
  if(!truthy(Sys.getenv("RUN_MRLINK2",unset="TRUE"))||!nrow(jobs))return(invisible(0L))
  sh<-file.path(Sys.getenv("LE8_FDIR"),"c2_mr_link2.sh");status<-system2("bash",c(sh,"--jobs",jobs_file,"--cad-gwas",ygfile,"--outdir",linkdir))
  if(status!=0)warning("MR-link-2 runner returned status ",status);invisible(status)
}

plot_c2_fig1 <- function(mr,assoc,layer){
  local_name<-if(layer=="protein")"cis"else"local"
  best<-mr|>filter(is.finite(pval))|>group_by(exposure)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()|>
    left_join(assoc|>select(exposure=term,obs_beta=beta,obs_se=std.error,obs_p=p.value),by="exposure")|>
    mutate(concordance=case_when(FDR_all<.05&sign(b)==sign(obs_beta)~"FDR-significant, concordant",
                                 FDR_all<.05&sign(b)!=sign(obs_beta)~"FDR-significant, discordant",TRUE~"No FDR MR support"))
  if(!nrow(best))return(list(plot=blank_plot("C2 MR evidence","No harmonized MR result"),bars=tibble(),forest=tibble(),best=best))
  bars<-bind_rows(best|>transmute(panel="All genetically tested traits",concordance),
                  best|>filter(obs_p<.05/max(1,nrow(assoc)))|>transmute(panel="Observationally associated traits",concordance),
                  best|>filter(analysis==local_name)|>transmute(panel=paste0("Primary ",local_name," evidence"),concordance))|>
    count(panel,concordance)|>group_by(panel)|>mutate(percent=100*n/sum(n))|>ungroup()
  pA<-ggplot(bars,aes(panel,percent,fill=concordance))+geom_col(width=.7,color="white")+geom_text(aes(label=paste0("n=",n)),position=position_stack(vjust=.5),size=3,fontface="bold")+coord_flip()+
    labs(title="a. Concordance of observational and genetically informed evidence",x=NULL,y="Percent",fill=NULL)+theme_5c(10)+theme(legend.position="top")
  top<-best|>arrange(pval)|>slice_head(n=TOP_FOREST)|>mutate(exposure=factor(exposure,levels=rev(exposure)))
  fl<-bind_rows(top|>transmute(exposure,evidence="Observational",b=obs_beta,se=obs_se,p=obs_p),top|>transmute(exposure,evidence=paste0("MR (",analysis,")"),b,se,p=pval))|>mutate(lo=b-1.96*se,hi=b+1.96*se)
  pB<-ggplot(fl,aes(b,exposure))+geom_vline(xintercept=0,color="grey55")+geom_errorbarh(aes(xmin=lo,xmax=hi),height=0,color="grey50",na.rm=TRUE)+
    geom_point(aes(color=p<.05),size=1.9,na.rm=TRUE)+facet_grid(.~evidence,scales="free_x")+scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="grey55"),guide="none")+
    labs(title="b. Top traits",x="Effect per 1-SD higher omic trait",y=NULL)+forest_theme(9)
  list(plot=pA/pB+plot_layout(heights=c(.55,1.45)),bars=bars,forest=fl,best=best)
}

plot_c2_fig2 <- function(mr,assoc,layer){
  local_name<-if(layer=="protein")"cis"else"local";distal_name<-if(layer=="protein")"trans"else"distal"
  wide<-mr|>select(exposure,analysis,b,se,pval)|>pivot_wider(names_from=analysis,values_from=c(b,se,pval))|>
    left_join(assoc|>select(exposure=term,obs_beta=beta,obs_se=std.error,obs_p=p.value),by="exposure")
  bl<-paste0("b_",local_name);sl<-paste0("se_",local_name);pl<-paste0("pval_",local_name);bd<-paste0("b_",distal_name);sd<-paste0("se_",distal_name);pd<-paste0("pval_",distal_name)
  a<-wide|>filter(is.finite(.data[[bl]]),is.finite(obs_beta))|>mutate(label=ifelse(min_rank(pmin(obs_p,.data[[pl]],na.rm=TRUE))<=12,exposure,NA_character_))
  pA<-if(!nrow(a))blank_plot(paste0("a. Observational versus ",local_name," MR"),"No paired estimate") else ggplot(a,aes(obs_beta,.data[[bl]]))+geom_hline(yintercept=0,color="grey70")+geom_vline(xintercept=0,color="grey70")+
    geom_smooth(method="lm",se=TRUE,color="grey35",fill="grey85",linewidth=.7)+geom_point(aes(color=.data[[pl]]<.05),size=2)+ggrepel::geom_text_repel(aes(label=label),size=2.8,fontface="bold",seed=1,max.overlaps=20,na.rm=TRUE)+
    scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="grey70"),guide="none")+labs(title=paste0("a. Observational versus ",local_name,"-MR effects"),x="Observational effect",y=paste0(local_name," MR effect"))+theme_5c(11)
  b<-wide|>filter(is.finite(.data[[bl]]),is.finite(.data[[bd]]))|>mutate(label=ifelse(min_rank(pmin(.data[[pl]],.data[[pd]],na.rm=TRUE))<=12,exposure,NA_character_))
  pB<-if(!nrow(b))blank_plot(paste0("b. ",local_name," versus ",distal_name),"No trait had both instrument classes") else ggplot(b,aes(.data[[bl]],.data[[bd]]))+geom_abline(slope=1,intercept=0,linetype=2,color="grey45")+
    geom_hline(yintercept=0,color="grey70")+geom_vline(xintercept=0,color="grey70")+geom_smooth(method="lm",se=TRUE,color="grey35",fill="grey85",linewidth=.7)+
    geom_point(aes(color=sign(.data[[bl]])==sign(.data[[bd]])),size=2)+ggrepel::geom_text_repel(aes(label=label),size=2.8,fontface="bold",seed=2,max.overlaps=20,na.rm=TRUE)+
    scale_color_manual(values=c(`TRUE`="#1B9E77",`FALSE`="#D95F02"),guide="none")+labs(title=paste0("b. ",local_name," versus ",distal_name," MR"),x=paste0(local_name," effect"),y=paste0(distal_name," effect"))+theme_5c(11)
  list(plot=pA|pB,wide=wide)
}

# ----------------------------------------------------------------------------
# DANDELION helpers
# ----------------------------------------------------------------------------
pick_col_local <- function(nms,patterns){
  nl<-toupper(nms)
  for(p in patterns){i<-grep(p,nl,perl=TRUE)[1];if(!is.na(i))return(nms[[i]])}
  NA_character_
}

find_magma_gene_file <- function(ygfile){
  explicit<-Sys.getenv("C2_MAGMA_GENE_FILE",unset="");if(nzchar(explicit)&&file.exists(explicit))return(normalizePath(explicit,winslash="/",mustWork=FALSE))
  d<-dirname(ygfile);trait<-sub("\\.gz$","",basename(ygfile),ignore.case=TRUE)
  # Standard gwas_post.sh layout is <GWAS root>/magma/<trait>/<trait>.genes.out,
  # while older runs sometimes placed the result beside the cleaned GWAS.
  clean_root<-dirname(d);gwas_root<-dirname(clean_root)
  dirs<-unique(c(d,file.path(d,paste0(trait,".magma")),file.path(gwas_root,"magma",trait)))
  explicit_candidates<-c(file.path(gwas_root,"magma",trait,paste0(trait,".genes.out")),file.path(d,paste0(trait,".genes.out")))
  discovered<-unlist(lapply(dirs,function(path){if(dir.exists(path))list.files(path,pattern="(genes\\.out|gene\\.out)$",full.names=TRUE,ignore.case=TRUE)else character()}),use.names=FALSE)
  fs<-unique(c(explicit_candidates,discovered));fs<-fs[file.exists(fs)&!grepl("(gsa|sets)\\.out$",fs,ignore.case=TRUE)]
  if(length(fs))return(normalizePath(fs[[which.max(file.info(fs)$mtime)]],winslash="/",mustWork=FALSE))
  NA_character_
}

map_magma_gene_ids <- function(ids,target_symbols){
  ids<-as.character(ids);target_symbols<-unique(as.character(target_symbols));over<-mean(ids%in%target_symbols,na.rm=TRUE)
  if(is.finite(over)&&over>.05)return(ids)
  mapfile<-Sys.getenv("C2_MAGMA_GENE_MAP",unset="")
  if(nzchar(mapfile)&&file.exists(mapfile)){
    m<-data.table::fread(mapfile,showProgress=FALSE);if(ncol(m)>=2){mp<-setNames(as.character(m[[2]]),as.character(m[[1]]));z<-unname(mp[ids]);z[is.na(z)]<-ids[is.na(z)];return(z)}
  }
  if(requireNamespace("org.Hs.eg.db",quietly=TRUE)&&requireNamespace("AnnotationDbi",quietly=TRUE)){
    z<-suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,keys=ids,column="SYMBOL",keytype="ENTREZID",multiVals="first"));z<-as.character(z);z[is.na(z)]<-ids[is.na(z)];return(z)
  }
  ids
}

read_magma_gene_p <- function(file,target_symbols){
  if(is.na(file)||!file.exists(file))return(numeric())
  d<-data.table::fread(file,showProgress=FALSE,check.names=FALSE);nms<-names(d)
  gc<-pick_col_local(nms,c("^GENE$","^GENE_ID$","^ID$","^SYMBOL$"));pc<-pick_col_local(nms,c("^P$","^PVAL$","^P_VALUE$"))
  if(is.na(gc)||is.na(pc))stop("MAGMA gene file needs a gene identifier column and P column: ",file,call.=FALSE)
  gn<-map_magma_gene_ids(d[[gc]],target_symbols);p<-suppressWarnings(as.numeric(d[[pc]]));ok<-!is.na(gn)&nzchar(gn)&is.finite(p)&p>0&p<=1
  p<-p[ok];names(p)<-gn[ok];p<-tapply(p,names(p),min,na.rm=TRUE);as.numeric(p)|>setNames(names(p))
}

read_lead_file <- function(file,ygfile){
  if(is.na(file)||!file.exists(file))return(tibble())
  d<-data.table::fread(file,showProgress=FALSE,check.names=FALSE,fill=TRUE);nms<-names(d)
  sc<-pick_col_local(nms,c("^SNP$","^RSID$","^ID$"));cc<-pick_col_local(nms,c("^CHR$","^CHROM$"));bc<-pick_col_local(nms,c("^BP$","^POS$","^POSITION$"));pc<-pick_col_local(nms,c("^P$","^PVAL$","^P_VALUE$","^P_BOLT_LMM$"))
  if(is.na(sc))return(tibble())
  z<-tibble(SNP=as.character(d[[sc]]),CHR=if(!is.na(cc))as.character(d[[cc]]) else NA_character_,POS=if(!is.na(bc))as.numeric(d[[bc]]) else NA_real_,P=if(!is.na(pc))as.numeric(d[[pc]]) else NA_real_)
  miss<-z$SNP[!is.finite(z$P)|!is.finite(z$POS)|is.na(z$CHR)];if(length(miss)){
    full<-read_sumstat_snps(ygfile,miss,N_DEFAULT);if(nrow(full))z<-z|>left_join(full|>select(SNP,CHR2=CHR,POS2=POS,P2=P),by="SNP")|>mutate(CHR=coalesce(CHR,CHR2),POS=coalesce(POS,POS2),P=coalesce(P,P2))|>select(-ends_with("2"))
  }
  z|>filter(!is.na(SNP),SNP!="",!is.na(CHR),CHR!="",is.finite(POS))|>distinct(SNP,.keep_all=TRUE)
}

stream_gws_hits <- function(ygfile,threshold=DANDELION_GWS){
  hdr<-tryCatch(readLines(if(grepl("\\.gz$",ygfile))gzfile(ygfile,"rt") else file(ygfile,"rt"),n=1),error=function(e)"")
  if(!nzchar(hdr))return(tibble());sep<-if(grepl("\\t",hdr))"\t"else" ";nms<-strsplit(trimws(hdr),"[[:space:]]+")[[1]]
  sc<-pick_col_local(nms,c("^SNP$","^RSID$","^ID$","^VARIANT_ID$"));cc<-pick_col_local(nms,c("^CHR$","^CHROM$"));bc<-pick_col_local(nms,c("^BP$","^POS$","^POSITION$"));pc<-pick_col_local(nms,c("^P$","^PVAL$","^P_VALUE$","^P_BOLT_LMM$"))
  idx<-match(c(sc,cc,bc,pc),nms);if(any(!is.finite(idx)))return(tibble())
  if(.Platform$OS.type!="windows"&&nzchar(Sys.which("awk"))){
    dec<-if(grepl("\\.gz$",ygfile))paste("gzip -cd",shQuote(ygfile))else paste("cat",shQuote(ygfile))
    awk<-sprintf("awk 'BEGIN{OFS=\"\\t\"} NR==1{next} $%d<=%g {print $%d,$%d,$%d,$%d}'",idx[4],threshold,idx[1],idx[2],idx[3],idx[4])
    z<-tryCatch(data.table::fread(cmd=paste(dec,"|",awk),col.names=c("SNP","CHR","POS","P"),showProgress=FALSE),error=function(e)NULL)
    if(!is.null(z))return(as_tibble(z)|>mutate(SNP=as.character(SNP),CHR=as.character(CHR),POS=as.numeric(POS),P=as.numeric(P))|>
      filter(!is.na(SNP),SNP!="",!is.na(CHR),CHR!="",is.finite(POS),is.finite(P),P<=threshold))
  }
  d<-data.table::fread(ygfile,select=unique(idx),showProgress=FALSE,check.names=FALSE)
  tibble(SNP=as.character(d[[sc]]),CHR=as.character(d[[cc]]),POS=as.numeric(d[[bc]]),P=as.numeric(d[[pc]]))|>
    filter(!is.na(SNP),SNP!="",!is.na(CHR),CHR!="",is.finite(POS),is.finite(P),P<=threshold)
}

distance_prune_leads <- function(z,window=DANDELION_LEAD_BP,max_n=DANDELION_MAX_SNPS){
  if(!nrow(z))return(z)
  window<-suppressWarnings(as.numeric(window));if(length(window)!=1L||!is.finite(window)||window<0)stop("DANDELION lead-pruning window must be a finite non-negative number.",call.=FALSE)
  max_n<-suppressWarnings(as.integer(max_n));if(length(max_n)!=1L||is.na(max_n)||max_n<0)stop("DANDELION_MAX_SNPS must be a non-negative integer.",call.=FALSE)
  z<-z|>mutate(CHR=str_remove(as.character(CHR),"^chr"),POS=suppressWarnings(as.numeric(POS)),P=suppressWarnings(as.numeric(P)))|>
    filter(!is.na(SNP),SNP!="",!is.na(CHR),CHR!="",is.finite(POS),is.finite(P))|>arrange(P)
  if(!nrow(z))return(z)
  keep<-logical(nrow(z));picked<-list()
  for(i in seq_len(nrow(z))){
    chr<-z$CHR[[i]];pos<-z$POS[[i]];old<-picked[[chr]]%||%numeric()
    if(!length(old)||all(abs(pos-old)>window)){keep[[i]]<-TRUE;picked[[chr]]<-c(old,pos)}
  }
  z<-z[keep,,drop=FALSE];if(max_n>0)z<-head(z,max_n);z
}

find_disease_leads <- function(ygfile){
  explicit<-Sys.getenv("C2_DANDELION_SNP_FILE",unset="");cand<-character()
  if(nzchar(explicit))cand<-c(cand,explicit)
  d<-dirname(ygfile);cand<-c(cand,list.files(d,pattern="(jma\\.cojo$|\\.clumped$|lead.*\\.(txt|tsv|csv)$|indep.*\\.(txt|tsv|csv)$)",full.names=TRUE,ignore.case=TRUE))
  cand<-unique(cand[file.exists(cand)])
  if(length(cand)){
    for(f in cand){z<-read_lead_file(f,ygfile);if(nrow(z)>=2){attr(z,"source")<-f;return(head(z|>arrange(P),DANDELION_MAX_SNPS))}}
  }
  warning("DANDELION: no independent/lead-SNP file found beside outcome GWAS; using distance-pruned genome-wide significant SNPs. For the primary analysis, provide an LD-independent file via C2_DANDELION_SNP_FILE.",call.=FALSE)
  z<-distance_prune_leads(stream_gws_hits(ygfile,DANDELION_GWS));attr(z,"source")<-"distance-pruned GWS from outcome GWAS";z
}

build_trans_p_matrix <- function(gene2,lead,base){
  snps<-lead$SNP
  rows<-parallel_map(gene2,function(g){
    fs<-find_qtl_files(g,base,"protein");f<-fs$full
    if(is.na(f)||!file.exists(f))return(rep(NA_real_,length(snps)))
    z<-tryCatch(read_sumstat_snps(f,snps,N_DEFAULT),error=function(e)tibble());p<-rep(NA_real_,length(snps));names(p)<-snps
    if(nrow(z))p[z$SNP]<-z$P;p
  })
  M<-do.call(rbind,rows);rownames(M)<-gene2;colnames(M)<-snps;storage.mode(M)<-"double";M
}

# A comprehensive gene annotation is strongly preferred for SNP -> cis-gene
# interpretation.  If none is supplied, the assayed-protein annotation is used
# as a conservative fallback; DANDELION itself can then infer the nearest gene
# among those available in ref.table.keep.
read_dandelion_gene_ref <- function(assayed_ref){
  f<-Sys.getenv("C2_DANDELION_GENE_ANNOTATION",unset="")
  if(!nzchar(f)||!file.exists(f)){
    warning("DANDELION: C2_DANDELION_GENE_ANNOTATION was not supplied; SNP-to-cis-gene labels will be inferred only from assayed protein genes. Supply a genome-wide gene annotation table for the primary disease-driver analysis.",call.=FALSE)
    return(assayed_ref)
  }
  d<-data.table::fread(f,showProgress=FALSE,check.names=FALSE);nms<-names(d)
  gc<-pick_col_local(nms,c("^GENE_NAME$","^GENE$","^SYMBOL$","^GENESYMBOL$","^GENE_SYMBOL$"))
  tc<-pick_col_local(nms,c("^TYPE$","^GENE_TYPE$","^BIOTYPE$"))
  cc<-pick_col_local(nms,c("^CHROMOSOME$","^CHR$","^CHROM$"))
  sc<-pick_col_local(nms,c("^START$","^GENE_START$","^START_POS$"))
  ec<-pick_col_local(nms,c("^END$","^GENE_END$","^END_POS$"))
  if(any(is.na(c(gc,cc,sc,ec))))stop("C2_DANDELION_GENE_ANNOTATION needs gene/symbol, chromosome, start and end columns: ",f,call.=FALSE)
  z<-tibble(gene_name=as.character(d[[gc]]),type=if(!is.na(tc))as.character(d[[tc]]) else "protein_coding",
    Chromosome=paste0("chr",str_remove(as.character(d[[cc]]),regex("^chr",ignore_case=TRUE))),
    start=suppressWarnings(as.numeric(d[[sc]])),end=suppressWarnings(as.numeric(d[[ec]])))|>
    filter(!is.na(gene_name),gene_name!="",is.finite(start),is.finite(end))|>
    mutate(type=coalesce(na_if(type,""),"protein_coding"),type=case_when(str_detect(str_to_lower(type),"protein")~"protein_coding",str_detect(str_to_lower(type),"linc|lnc")~"lincRNA",TRUE~type))|>distinct(gene_name,.keep_all=TRUE)
  # Always retain exact coordinates for the assayed proteins used as gene2.
  bind_rows(assayed_ref,z|>filter(!gene_name%in%assayed_ref$gene_name))|>distinct(gene_name,.keep_all=TRUE)
}

read_dandelion_snp_gene_map <- function(lead){
  f<-Sys.getenv("C2_DANDELION_SNP_GENE_MAP",unset="")
  empty<-data.frame(SNP=character(),GeneSymbol=character(),stringsAsFactors=FALSE)
  if(!nzchar(f)||!file.exists(f))return(empty)
  d<-data.table::fread(f,showProgress=FALSE,check.names=FALSE);nms<-names(d)
  sc<-pick_col_local(nms,c("^SNP$","^RSID$","^ID$"));gc<-pick_col_local(nms,c("^GENESYMBOL$","^GENE_SYMBOL$","^SYMBOL$","^GENE$"))
  if(is.na(sc)||is.na(gc))stop("C2_DANDELION_SNP_GENE_MAP needs SNP/rsID and GeneSymbol/gene columns: ",f,call.=FALSE)
  d<-tibble(SNP=as.character(d[[sc]]),GeneSymbol=as.character(d[[gc]]))|>
    filter(SNP%in%lead$SNP,!is.na(GeneSymbol),GeneSymbol!="")|>distinct(SNP,.keep_all=TRUE)
  as.data.frame(d)
}

safe_min_finite <- function(x){x<-x[is.finite(x)];if(length(x))min(x)else NA_real_}

plot_dandelion_results <- function(tg,pd,gene_pair,outdir){
  if(!nrow(tg)){
    save_plot(blank_plot("DANDELION disease-driver prioritization","No SNP–gene2 pair passed the target FDR"),"c2.Fig6.dandelion.png",10,6,outdir=outdir)
    return(invisible(NULL))
  }
  top<-tg|>slice_head(n=30)|>mutate(gene2=factor(gene2,levels=rev(gene2)))
  pA<-ggplot(top,aes(-log10(pmax(DANDELION_p,1e-300)),gene2))+
    geom_segment(aes(x=0,xend=-log10(pmax(DANDELION_p,1e-300)),yend=gene2),color="grey76")+
    geom_point(aes(size=n_distal_loci,color=-log10(pmax(MAGMA_p,1e-300))))+
    labs(title="a. DANDELION-prioritized disease-proximal genes",x=expression(-log[10](DANDELION~P)),y=NULL,
         size="Disease loci",color=expression(-log[10](gene-level~P)))+theme_5c(10)
  pB<-ggplot(pd,aes(-log10(pmax(trans_p,1e-300)),-log10(pmax(magma_p,1e-300))))+
    geom_point(aes(size=-log10(pmax(DANDELION_p,1e-300)),color=-log10(pmax(DANDELION_p,1e-300))),alpha=.72)+
    ggrepel::geom_text_repel(data=pd|>arrange(DANDELION_p)|>slice_head(n=20),aes(label=gene2),size=2.6,fontface="bold",seed=12,max.overlaps=Inf)+
    labs(title="b. Trans-regulatory + disease-gene evidence",x=expression(-log[10](trans-regulatory~P)),
         y=expression(-log[10](gene-level~disease~P)),size=expression(-log[10](DANDELION~P)),color=expression(-log[10](DANDELION~P)))+
    scale_color_viridis_c(option="C",end=.9)+theme_5c(10)
  gp<-as_tibble(gene_pair)
  if(nrow(gp)&&all(c("region","gene2")%in%names(gp))){
    if(!"DANDELION_p"%in%names(gp))gp<-gp|>left_join(pd|>select(gene2,DANDELION_p)|>group_by(gene2)|>summarise(DANDELION_p=safe_min_finite(DANDELION_p),.groups="drop"),by="gene2")
    keepg<-tg|>slice_head(n=20)|>pull(gene2);gp<-gp|>filter(gene2%in%keepg)|>distinct(region,gene2,.keep_all=TRUE)|>slice_head(n=100)
    pC<-if(!nrow(gp))blank_plot("c. Disease locus → disease-proximal gene map") else ggplot(gp,aes(fct_reorder(region,-log10(pmax(DANDELION_p,1e-300)),max),fct_reorder(gene2,-log10(pmax(DANDELION_p,1e-300)),max)))+
      geom_point(aes(size=-log10(pmax(DANDELION_p,1e-300)),color=-log10(pmax(DANDELION_p,1e-300))),alpha=.78)+coord_flip()+
      scale_color_viridis_c(option="D",end=.88)+labs(title="c. Disease locus → disease-proximal gene map",x="Distal disease locus / cis gene",y="Prioritized gene2",size=expression(-log[10](P)),color=expression(-log[10](P)))+theme_5c(9)
  } else pC<-blank_plot("c. Disease locus → disease-proximal gene map","No interpretable SNP-to-cis-gene mapping was available")
  save_plot((pA|pB)/pC+plot_layout(heights=c(1.05,.95)),"c2.Fig6.dandelion.png",16,13,outdir=outdir)
}

run_dandelion_step <- function(layer,assoc,ann,base,ygfile,rawdir,outdir){
  empty<-list(status="not_run",pairs=tibble(),gene_pairs=tibble(),targets=tibble(),lead_snps=tibble(),snp_gene_map=tibble(),magma_file=NA_character_)
  if(!RUN_DANDELION)return(modifyList(empty,list(status="disabled")))
  if(layer!="protein")return(modifyList(empty,list(status="protein_only")))
  if(!requireNamespace("DANDELION",quietly=TRUE)){
    warning("DANDELION R package is not installed; run Rscript setup_dandelion.R or remotes::install_github('mxxptian/DANDELION').",call.=FALSE)
    return(modifyList(empty,list(status="package_missing")))
  }
  magma<-find_magma_gene_file(ygfile)
  if(is.na(magma)){
    gsa<-list.files(dirname(ygfile),pattern="gsa\\.out$",full.names=TRUE,ignore.case=TRUE)
    msg<-if(length(gsa))"MAGMA gene-set (*.gsa.out) found, but DANDELION requires gene-level (*.genes.out) P values." else "No MAGMA gene-level (*.genes.out) result found."
    warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg)))
  }
  ref_assayed<-ann|>filter(!is.na(feature),is.finite(start),is.finite(end),!is.na(chr))|>
    transmute(gene_name=feature,type="protein_coding",Chromosome=paste0("chr",str_remove(chr,"^chr")),start=start,end=end)|>distinct(gene_name,.keep_all=TRUE)
  ref_all<-read_dandelion_gene_ref(ref_assayed)
  p_gene<-read_magma_gene_p(magma,ref_assayed$gene_name);common<-intersect(ref_assayed$gene_name,names(p_gene))
  if(length(common)<20){
    msg<-paste0("Only ",length(common)," MAGMA genes matched assayed protein symbols. If MAGMA uses Entrez IDs, install org.Hs.eg.db or set C2_MAGMA_GENE_MAP.")
    warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,magma_file=magma)))
  }
  lead<-find_disease_leads(ygfile)|>filter(SNP!="",is.finite(POS))|>distinct(SNP,.keep_all=TRUE)
  if(nrow(lead)<2){msg<-"Fewer than two independent disease-associated lead SNPs available.";warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,magma_file=magma,lead_snps=lead)))}
  gene2<-sort(common);if(DANDELION_MAX_GENE2>0&&length(gene2)>DANDELION_MAX_GENE2){set.seed(SEED);gene2<-sort(sample(gene2,DANDELION_MAX_GENE2))}
  message("C2/DANDELION: ",length(lead$SNP)," disease SNPs x ",length(gene2)," assayed gene2 proteins; gene-level disease evidence=",basename(magma))
  ptrans<-build_trans_p_matrix(gene2,lead,base);okrow<-rowSums(is.finite(ptrans))>=2;ptrans<-ptrans[okrow,,drop=FALSE];pwes<-p_gene[rownames(ptrans)]
  if(nrow(ptrans)<20){msg<-"Too few assayed proteins had trans-pQTL data at the selected disease SNPs.";warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,magma_file=magma,lead_snps=lead)))}
  snpref<-lead|>transmute(SNP,SNPPos=POS,SNPChr=str_remove(as.character(CHR),"^chr"))
  # ref.table may contain genome-wide genes; med_gene() uses it to remove local
  # cis targets.  Only gene2 present in p.trans/p.wes can become DPGs.
  ref_for_med<-ref_all|>filter(type%in%c("lincRNA","protein_coding"),!Chromosome%in%c("chrM","chrX","chrY"))|>distinct(gene_name,.keep_all=TRUE)
  uniq_snp<-read_dandelion_snp_gene_map(lead)
  med<-tryCatch(DANDELION::med_gene(p.trans=ptrans,p.wes=pwes,ref.table=ref_for_med,gene1.list=colnames(ptrans),target.fdr=DANDELION_FDR,
    dist=DANDELION_CIS_BP,gene1.type="SNP",SNP.ref=snpref,n.cores=N_CORES,verbose=TRUE),error=function(e)e)
  if(inherits(med,"condition")){msg<-conditionMessage(med);warning("DANDELION failed: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,magma_file=magma,lead_snps=lead,snp_gene_map=uniq_snp)))}
  pairs<-tryCatch(DANDELION::calc_pair.snp(mat.sig=med$mat.sig,mat.p=med$mat.p,p.wes=pwes,gene1=med$gene1,uniq_snp=uniq_snp,
    ref.table.keep=ref_for_med,eta.wgs=.05/max(1,length(pwes)),SNP.ref=snpref,verbose=TRUE),error=function(e)e)
  if(inherits(pairs,"condition")){warning("DANDELION post-processing failed: ",conditionMessage(pairs),call.=FALSE);pairs<-NULL}
  pd<-if(is.null(pairs))tibble() else as_tibble(pairs$pairs_dact)
  gp<-if(is.null(pairs))tibble() else as_tibble(pairs$gene.pair)
  if(nrow(pd)){
    pd<-pd|>mutate(trans_p=map2_dbl(gene2,rsid,~if(.x%in%rownames(ptrans)&&.y%in%colnames(ptrans))ptrans[.x,.y]else NA_real_),magma_p=pwes[gene2])
    tg<-pd|>group_by(gene2)|>summarise(n_distal_loci=n_distinct(rsid),DANDELION_p=safe_min_finite(DANDELION_p),MAGMA_p=safe_min_finite(magma_p),best_trans_p=safe_min_finite(trans_p),.groups="drop")|>arrange(DANDELION_p)
  } else tg<-tibble()
  # If no explicit SNP-to-gene file was supplied, expose the mapping inferred by
  # calc_pair.snp() so the distal-locus interpretation remains auditable.
  inferred_map<-if(nrow(pd)&&all(c("rsid","gene1")%in%names(pd)))pd|>transmute(SNP=rsid,GeneSymbol=gene1)|>filter(!is.na(GeneSymbol),GeneSymbol!="")|>distinct() else tibble()
  map_out<-bind_rows(as_tibble(uniq_snp),inferred_map)|>distinct(SNP,.keep_all=TRUE)
  write_raw_csv(lead,"c2.dandelion_lead_snps.csv",rawdir);write_raw_csv(map_out,"c2.dandelion_snp_to_cis_gene.csv",rawdir)
  write_raw_csv(pd,"c2.dandelion_pairs.csv",rawdir);write_raw_csv(gp,"c2.dandelion_gene_pairs.csv",rawdir);write_raw_csv(tg,"c2.dandelion_targets.csv",rawdir)
  plot_dandelion_results(tg,pd,gp,outdir)
  list(status="ok",result=med,pairs=pd,gene_pairs=gp,targets=tg,lead_snps=lead,snp_gene_map=map_out,magma_file=magma,
       ptrans_dimensions=dim(ptrans),target_fdr=DANDELION_FDR,cis_bp=DANDELION_CIS_BP,
       gene_annotation=Sys.getenv("C2_DANDELION_GENE_ANNOTATION",unset=""),snp_gene_map_file=Sys.getenv("C2_DANDELION_SNP_GENE_MAP",unset=""))
}


# ----------------------------------------------------------------------------
# Main C2
# ----------------------------------------------------------------------------
restore_c2_figures <- function(cached,layer,outdir){
  mr<-as_tibble(cached$MR%||%tibble());assoc<-as_tibble(cached$observational%||%tibble())
  if(!nrow(mr)||!nrow(assoc))return(invisible(FALSE))
  fig1<-plot_c2_fig1(mr,assoc,layer);save_plot(fig1$plot,"c2.Fig1.prots.top.png",13.5,9.8,outdir=outdir)
  fig2<-plot_c2_fig2(mr,assoc,layer)
  top_r2<-mr|>filter(is.finite(r2_bounded),is.finite(pval),n_IV>0)|>arrange(pval)|>group_by(exposure)|>slice(1)|>ungroup()|>slice_head(n=50)|>
    transmute(exposure,r2_total=100*(1-pmax(1-r2_bounded,1e-12)^(1/pmax(n_IV,1))),mean_F,n_IV,analysis,pval)|>arrange(r2_total)|>mutate(exposure=factor(exposure,levels=exposure))
  p3<-if(!nrow(top_r2))blank_plot("Variance explained by QTL instruments","No instrument set") else ggplot(top_r2,aes(r2_total,exposure,fill=analysis))+geom_col(width=.72)+geom_text(aes(label=sprintf("%.1f%%",r2_total)),hjust=-.1,size=2.7,fontface="bold")+
    scale_x_continuous(labels=label_percent(scale=1),limits=c(0,max(1,max(top_r2$r2_total,na.rm=TRUE)*1.18)))+labs(title="Mean variance explained per independent QTL instrument",x="Mean per-instrument explained variance",y=NULL,fill="Instrument class")+theme_5c(10)+theme(legend.position="bottom")
  save_plot(p3,"c2.Fig3.pQTL_R2.png",8.5,9,outdir=outdir)
  best<-cached$MR_best%||%fig1$best
  if(!"y"%in%names(best))best<-best|>mutate(y=-log10(pmax(pval,1e-300)),direction=case_when(FDR_all<.05&b>0~"Positive",FDR_all<.05&b<0~"Inverse",TRUE~"NS"),label=ifelse(min_rank(pval)<=25,exposure,NA_character_))
  p4<-ggplot(best,aes(b,y))+geom_hline(yintercept=-log10(.05/max(1,nrow(best))),linetype=2,color="grey50")+geom_vline(xintercept=0,color="grey65")+geom_point(aes(color=direction,size=pmin(n_IV,8)),alpha=.85)+
    ggrepel::geom_text_repel(aes(label=label),size=2.8,fontface="bold",seed=5,max.overlaps=Inf,na.rm=TRUE)+scale_color_manual(values=c(Positive="#D7301F",Inverse="#2C7FB8",NS="grey78"))+labs(title=paste0("Genetically predicted ",layer," levels and ",Y),x="MR effect",y=expression(-log[10](P)),color=NULL,size="Number of IVs")+theme_5c(12)+theme(legend.position="bottom")
  arch<-cached$architecture%||%tibble()
  p5<-if(!nrow(arch))blank_plot("QTL instrument-strength and pleiotropy atlas","No instrument metric")else ggplot(arch,aes(mean_F,100*r2_bounded))+geom_vline(xintercept=10,linetype=2,color="grey55")+geom_point(aes(size=n_IV,color=evidence),alpha=.8)+ggrepel::geom_text_repel(aes(label=label),size=2.7,fontface="bold",seed=8,max.overlaps=20,na.rm=TRUE)+scale_x_log10()+labs(title="QTL instrument-strength and pleiotropy atlas",x="Mean F statistic (log scale)",y="Bounded explained variance (%)",size="IV count",color=NULL)+theme_5c(11)+theme(legend.position="bottom")
  save_plot(fig2$plot/(p4|p5)+plot_layout(heights=c(1,1.08)),"c2.Fig2.concordance.png",16,13,outdir=outdir)
  dan<-cached$DANDELION%||%list();if(layer=="protein")plot_dandelion_results(dan$targets%||%tibble(),dan$pairs%||%tibble(),dan$gene_pairs%||%tibble(),outdir)
  invisible(TRUE)
}

read_c1_cache_for_c2 <- function(outdir){
  c1file <- file.path(le8_job_dir(outdir,"c1_correlate"),"c1.res.rds")
  wait_seconds <- suppressWarnings(as.numeric(Sys.getenv("C2_C1_WAIT_SECONDS",unset="60")))
  if(!is.finite(wait_seconds)||wait_seconds<0)wait_seconds<-0
  deadline <- Sys.time()+wait_seconds
  announced <- FALSE
  repeat{
    if(file.exists(c1file)&&isTRUE(file.size(c1file)>0)){
      c1 <- tryCatch(readRDS(c1file),error=function(e)e)
      if(!inherits(c1,"condition"))return(c1)
      last_error <- conditionMessage(c1)
    } else last_error <- "file does not exist or is empty"
    if(Sys.time()>=deadline){
      stop("C2 requires an existing, readable C1 result: ",c1file,
           " (",last_error,"). C2 does not run C1 automatically.",call.=FALSE)
    }
    if(!announced){
      message("C2: waiting up to ",wait_seconds," seconds for existing C1 result: ",c1file)
      announced <- TRUE
    }
    Sys.sleep(1)
  }
}

run_c2_layer <- function(layer=c("protein","metabolite")){
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir)
  rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE);cache<-file.path(rawdir,"c2.res.rds")
  unlink(file.path(outdir,c("c2.Fig4.genetic_protein_cvd_cad.png","c2.Fig4.genetic_metabolite_cvd_cad.png",
    paste0("c2.Fig4.genetic_protein_",Y,".png"),paste0("c2.Fig4.genetic_metabolite_",Y,".png"),
    "c2.Fig5.instrument_architecture.png","c2.Fig5.protein_LE8_variance.png","c2.Fig5.metabolite_LE8_variance.png",
    "c2.Fig6.instrument_architecture.png","c2.Fig7.dandelion.png")),force=TRUE)
  if(cache_valid(cache)){
    cached<-readRDS(cache)
    if(!RUN_DANDELION||!is.null(cached$DANDELION)){
      required<-c("c2.Fig1.prots.top.png","c2.Fig2.concordance.png","c2.Fig3.pQTL_R2.png")
      if(layer=="protein")required<-c(required,"c2.Fig6.dandelion.png")
      # Figure code evolves independently of the expensive MR/DANDELION cache;
      # always redraw from cached numerical results on a resume run.
      message("C2: valid result cache found; regenerating figures from cached results.")
      restore_c2_figures(cached,layer,outdir);finalize_outputs(LE8_JOB,outdir)
      run_mrlink2_step(layer,rawdir,cached$MRLink2_jobs%||%tibble(),get_y_gwas_file(Y,TRUE));return(cached)
    } else message("C2: cached result lacks required DANDELION output; rebuilding C2 outputs.")
  }
  c1<-read_c1_cache_for_c2(outdir);assoc<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble();if(!"beta"%in%names(assoc))assoc<-assoc|>mutate(beta=safe_log(estimate))
  base<-if(layer=="protein")dir.X else dir.met.gwas;ann<-layer_annotation(layer,unique(assoc$term));ranked<-assoc|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)
  mapped<-ranked[vapply(ranked,function(x){fs<-find_qtl_files(x,base,layer);any(!is.na(unlist(fs[c("joint","full","cis","trans")])))},logical(1))]
  candidates<-head(mapped,MAX_FEATURES);if(!length(candidates))stop("C2/",layer,": no QTL files mapped under ",base,call.=FALSE)
  write_raw_csv(tibble(feature=candidates),"c2.input_features.csv",rawdir);message("C2/",layer,": reading instruments for ",length(candidates)," traits")
  iv_list<-setNames(map(candidates,~read_qtl_instruments(.x,base,layer,ann)),candidates);all_snps<-unique(unlist(map(iv_list,~.$instruments$SNP)))
  ygfile<-get_y_gwas_file(Y,TRUE);ygwas<-read_sumstat_snps(ygfile,all_snps,N_DEFAULT);if(!nrow(ygwas))stop("No outcome-GWAS variants overlapped QTL instruments.",call.=FALSE)
  mr<-imap_dfr(iv_list,function(obj,feature){iv<-obj$instruments;if(!nrow(iv))return(tibble());map_dfr(unique(iv$analysis),~run_mr(iv|>filter(analysis==.x),ygwas,feature,.x))})
  if(!nrow(mr))stop("C2/",layer,": no harmonized MR estimate.",call.=FALSE)
  mr<-mr|>group_by(analysis)|>mutate(FDR_analysis=p.adjust(pval,"BH"))|>ungroup()|>mutate(FDR_all=p.adjust(pval,"BH"));write_raw_csv(mr,"c2.MR_all.csv",rawdir)
  fig1<-plot_c2_fig1(mr,assoc,layer);save_plot(fig1$plot,"c2.Fig1.prots.top.png",13.5,9.8,outdir=outdir)
  fig2<-plot_c2_fig2(mr,assoc,layer)

  r2<-mr|>filter(is.finite(r2_bounded),n_IV>0)|>mutate(analysis=factor(analysis));top_r2<-r2|>filter(is.finite(pval))|>arrange(pval)|>group_by(exposure)|>slice(1)|>ungroup()|>slice_head(n=50)|>
    transmute(exposure,r2_total=100*(1-pmax(1-r2_bounded,1e-12)^(1/pmax(n_IV,1))),mean_F,n_IV,analysis=as.character(analysis),pval)|>arrange(r2_total)|>mutate(exposure=factor(exposure,levels=exposure))
  p3<-if(!nrow(top_r2))blank_plot("Variance explained by QTL instruments","No instrument set") else ggplot(top_r2,aes(r2_total,exposure,fill=analysis))+geom_col(width=.72)+geom_text(aes(label=sprintf("%.1f%%",r2_total)),hjust=-.1,size=2.7,fontface="bold")+
    scale_x_continuous(labels=label_percent(scale=1),limits=c(0,max(1,max(top_r2$r2_total,na.rm=TRUE)*1.18)))+labs(title="Mean variance explained per independent QTL instrument",x="Mean per-instrument explained variance",y=NULL,fill="Instrument class")+theme_5c(10)+theme(legend.position="bottom")
  save_plot(p3,"c2.Fig3.pQTL_R2.png",8.5,9,outdir=outdir)

  best<-fig1$best|>mutate(y=-log10(pmax(pval,1e-300)),direction=case_when(FDR_all<.05&b>0~"Positive",FDR_all<.05&b<0~"Inverse",TRUE~"NS"),label=ifelse(min_rank(pval)<=25,exposure,NA_character_))
  p4<-ggplot(best,aes(b,y))+geom_hline(yintercept=-log10(.05/max(1,nrow(best))),linetype=2,color="grey50")+geom_vline(xintercept=0,color="grey65")+geom_point(aes(color=direction,size=pmin(n_IV,8)),alpha=.85)+
    ggrepel::geom_text_repel(aes(label=label),size=2.8,fontface="bold",seed=5,max.overlaps=Inf,na.rm=TRUE)+scale_color_manual(values=c(Positive="#D7301F",Inverse="#2C7FB8",NS="grey78"))+
    labs(title=paste0("Genetically predicted ",layer," levels and ",Y),x="MR effect",y=expression(-log[10](P)),color=NULL,size="Number of IVs")+theme_5c(12)+theme(legend.position="bottom")

  # LE8-to-omics connection analysis belongs to C4, not C2. Remove
  # the former C2 LE8-variance panel to keep the five-C separation explicit.

  arch<-mr|>filter(is.finite(mean_F),is.finite(r2_bounded))|>mutate(evidence=case_when(FDR_all<.05~"MR FDR < 0.05",Q_p<.05~"Heterogeneous",TRUE~"Other"),label=ifelse(min_rank(pval)<=15,exposure,NA_character_))
  p6<-if(!nrow(arch))blank_plot("QTL instrument-strength and pleiotropy atlas","No instrument metric")else ggplot(arch,aes(mean_F,100*r2_bounded))+geom_vline(xintercept=10,linetype=2,color="grey55")+geom_point(aes(size=n_IV,color=evidence),alpha=.8)+
    ggrepel::geom_text_repel(aes(label=label),size=2.7,fontface="bold",seed=8,max.overlaps=20,na.rm=TRUE)+scale_x_log10()+labs(title="QTL instrument-strength and pleiotropy atlas",x="Mean F statistic (log scale)",y="Bounded explained variance (%)",size="IV count",color=NULL)+theme_5c(11)+theme(legend.position="bottom")
  save_plot(fig2$plot/(p4|p6)+plot_layout(heights=c(1,1.08)),"c2.Fig2.concordance.png",16,13,outdir=outdir)

  # DANDELION is complementary to MR: it asks whether a distal disease SNP has
  # trans-regulatory evidence to a gene/protein that itself has disease-genetic evidence.
  dandelion<-run_dandelion_step(layer,assoc,ann,base,ygfile,rawdir,outdir)

  jobs<-imap_dfr(iv_list,function(obj,feature){iv<-obj$instruments;if(!nrow(iv))return(tibble());primary<-if(layer=="protein")"cis"else"local";lead_pool<-iv|>filter(analysis==primary);if(!nrow(lead_pool))lead_pool<-iv;lead<-lead_pool|>arrange(P)|>slice(1);qf<-if(!is.na(obj$files$full))obj$files$full else obj$files$joint;tibble(omics=layer,trait=feature,exposure=qf,outcome=ygfile,region=paste0(lead$CHR,":",max(1,lead$POS-1e7),"-",lead$POS+1e7))})|>filter(!is.na(exposure),file.exists(exposure))
  linkdir<-file.path(rawdir,"mrlink2");dir.create(linkdir,recursive=TRUE,showWarnings=FALSE);write_raw_tsv(jobs,"c2.mrlink2.jobs.tsv",linkdir);run_mrlink2_step(layer,rawdir,jobs,ygfile)

  out<-list(meta=module_meta(layer),MR=mr,MR_best=best,observational=assoc,Fig1_summary=fig1$bars,Fig1_forest=fig1$forest,Fig2=fig2$wide,R2_QTL=top_r2,architecture=arch,DANDELION=dandelion,MRLink2_jobs=jobs)
  saveRDS(out,cache,compress="xz")
  write_xlsx2(list(MR_all=mr,MR_best=best,concordance_summary=fig1$bars,concordance_forest=fig1$forest,cis_local_vs_trans_distal=fig2$wide,QTL_R2=top_r2,instrument_architecture=arch,
                   DANDELION_lead_snps=dandelion$lead_snps%||%tibble(),DANDELION_snp_gene_map=dandelion$snp_gene_map%||%tibble(),DANDELION_pairs=dandelion$pairs%||%tibble(),DANDELION_gene_pairs=dandelion$gene_pairs%||%tibble(),DANDELION_targets=dandelion$targets%||%tibble(),MRLink2_jobs=jobs),"c2.out.xlsx")
  finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO)invisible(run_c2_layer("protein"))
if(met_DO)invisible(run_c2_layer("metabolite"))
