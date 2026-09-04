# C2: cis/trans MR and DANDELION driver prioritization using gene-level disease evidence.

suppressPackageStartupMessages({
  .c2_fdir<-Sys.getenv("LE8_FDIR", unset=file.path(Sys.getenv("DIRSCRIPT"),"f"))
  source(file.path(.c2_fdir,"comm.f.R"))
  source(file.path(.c2_fdir,"c2_genetic_decompose.R"))
})
LE8_JOB <- "c2_cause"
C2_CODE_VERSION <- "2026-09-03.genetic_components_le4_3"
C2_STAGE_VERSION <- "2026-08-31.final1" # MR/DANDELION estimators are unchanged
MAX_FEATURES <- as.integer(Sys.getenv("C2_MAX_FEATURES", unset="500"))
TOP_FOREST <- as.integer(Sys.getenv("C2_TOP_FOREST", unset="32"))
N_DEFAULT <- as.numeric(Sys.getenv("C2_SUMSTAT_N", unset="100000"))
normalize_c2_mode <- function(x,name){
  x<-str_to_lower(str_trim(as.character(x)[1]));x<-case_when(
    x=="top"~"Top",x%in%c("all","true","1","yes","y")~"All",
    x%in%c("none","false","0","no","n")~"None",TRUE~NA_character_)
  if(is.na(x))stop(name," must be Top, All, or None.",call.=FALSE)
  x
}
# These names are intentionally case-sensitive.  RUN_MRLINK2 and
# RUN_DANDELION are retired so logs and batch scripts expose one contract.
RUN_MRlink2 <- normalize_c2_mode(Sys.getenv("RUN_MRlink2",unset="Top"),"RUN_MRlink2")
RUN_Dandelion <- normalize_c2_mode(Sys.getenv("RUN_Dandelion",unset="Top"),"RUN_Dandelion")
C2_FIXED_TOP <- unique(trimws(strsplit(Sys.getenv("C2_TOP_CANDIDATES",
  unset="PCSK9,LPA,GDF15,NTPROBNP,MMP12"),",",fixed=TRUE)[[1]]))
C2_TOP_MAX <- suppressWarnings(as.integer(Sys.getenv("C2_TOP_MAX",unset="50")))
if(!is.finite(C2_TOP_MAX)||C2_TOP_MAX<5L)C2_TOP_MAX<-50L
DANDELION_FDR <- as.numeric(Sys.getenv("C2_DANDELION_FDR", unset="0.10"))
DANDELION_CIS_BP <- as.numeric(Sys.getenv("C2_DANDELION_CIS_BP", unset="5000000"))
DANDELION_GWS <- as.numeric(Sys.getenv("C2_DANDELION_GWS", unset="5e-8"))
DANDELION_LEAD_BP <- as.numeric(Sys.getenv("C2_DANDELION_LEAD_BP", unset="5000000"))
DANDELION_MAX_SNPS <- as.integer(Sys.getenv("C2_DANDELION_MAX_SNPS", unset="100"))
DANDELION_MAX_GENE2 <- as.integer(Sys.getenv("C2_DANDELION_MAX_GENE2", unset="0")) # 0 = all available
DANDELION_ALLOW_MAGMA <- truthy(Sys.getenv("C2_DANDELION_ALLOW_MAGMA",unset="TRUE"))
DANDELION_MAX_TARGET_FRACTION <- as.numeric(Sys.getenv("C2_DANDELION_MAX_TARGET_FRACTION",unset="0.25"))

c2_method_scope_audit <- function(layer){
  tibble(method=c("cis/local MR","trans/distal MR","MR-link-2","DANDELION","CIGMA"),
    role=c("causal estimation","pleiotropy-sensitive causal estimation","regional LD-aware causal estimation",
      "distal regulatory driver prioritization","cell-type-shared/specific eQTL variance decomposition"),
    eligible_with_current_inputs=c(TRUE,TRUE,layer=="protein",layer=="protein",FALSE),
    decision=c("core C2","core C2 with stricter interpretation","optional C2 regional sensitivity",
      ifelse(layer=="protein","primary only with valid disease gene/rare-variant evidence; MAGMA fallback is sensitivity","not defined for metabolites"),
      "do not add to C2 causation core; optional cell-type annotation after external CIGMA analysis"),
    required_input=c("omic QTL + disease GWAS","omic QTL + disease GWAS","dense cis/local QTL + disease GWAS + LD reference",
      "distal QTL plus valid disease-gene evidence",
      "population-scale single-cell RNA-seq h5ad, donor/cell metadata, pseudobulk and kinship"),
    reference=c(NA_character_,NA_character_,NA_character_,NA_character_,
      "Nature 2026, doi:10.1038/s41586-026-10577-6"))
}
if(!is.finite(DANDELION_MAX_TARGET_FRACTION)||DANDELION_MAX_TARGET_FRACTION<=0||DANDELION_MAX_TARGET_FRACTION>1)
  DANDELION_MAX_TARGET_FRACTION<-.25

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

select_c2_top_candidates <- function(layer,assoc,mr,universe){
  primary<-if(layer=="protein")"cis"else"local"
  fixed<-C2_FIXED_TOP[C2_FIXED_TOP%in%universe]
  obs<-as_tibble(assoc)|>filter(term%in%universe,is.finite(p.value))|>
    arrange(p.value)|>slice_head(n=20)|>pull(term)
  gen<-as_tibble(mr)|>filter(exposure%in%universe,analysis==primary,is.finite(pval))|>
    arrange(pval)|>slice_head(n=25)|>pull(exposure)
  x<-head(unique(c(fixed,gen,obs)),C2_TOP_MAX)
  tibble(feature=x,top_rank=seq_along(x),top_reason=case_when(
    feature%in%fixed~"prespecified anchor",
    feature%in%gen~paste0("top ",primary," MR"),TRUE~"top C1 association"))
}

run_mrlink2_step <- function(layer,rawdir,jobs,ygfile,mode=RUN_MRlink2){
  linkdir<-file.path(rawdir,"mrlink2",str_to_lower(mode));dir.create(linkdir,recursive=TRUE,showWarnings=FALSE)
  complete_file<-file.path(linkdir,"mrlink2.complete")
  if(mode=="None"){
    message("C2/",layer,": MR-link-2 not run (RUN_MRlink2=None)")
    return(invisible(0L))
  }
  if(cache_valid(complete_file)){
    cache_message(paste0("MR-link-2/",layer),complete_file)
    return(invisible(0L))
  }
  if(!nrow(jobs))return(invisible(0L))
  jobs_file<-file.path(linkdir,"c2.mrlink2.jobs.tsv");write_raw_tsv(jobs,basename(jobs_file),linkdir)
  ref_dir<-Sys.getenv("MRLINK2_REF_PFILE_DIR",unset="");ref_bed<-Sys.getenv("MRLINK2_REF_BED",unset="")
  ref_pop<-toupper(Sys.getenv("MRLINK2_REF_POP",unset="EUR"))
  ref_id_dir<-Sys.getenv("MRLINK2_REF_ID_DIR",unset="")
  ref_samples<-Sys.getenv("MRLINK2_REF_SAMPLES",unset="")
  if(!nzchar(ref_id_dir)&&nzchar(ref_dir))ref_id_dir<-file.path(dirname(ref_dir),"id")
  if(!nzchar(ref_samples)&&nzchar(ref_dir))ref_samples<-file.path(dirname(ref_dir),"samples.txt")
  ref_keep<-if(!nzchar(ref_bed)&&nzchar(ref_id_dir)&&ref_pop!="ALL")file.path(ref_id_dir,paste0(ref_pop,".id.2col"))else ""
  ref_bfile_dir<-if(!nzchar(ref_bed)&&nzchar(ref_dir))file.path(dirname(ref_dir),"bfile",ref_pop)else ""
  if(nzchar(ref_bed)){ref_keep<-"";ref_samples<-""}
  sh<-file.path(Sys.getenv("LE8_FDIR"),"c2_mr_link2.sh")
  status<-system2("bash",c(sh,"--jobs",jobs_file,"--cad-gwas",ygfile,"--outdir",linkdir))
  if(status!=0)warning("MR-link-2 runner returned status ",status);invisible(status)
}

# MR-link-2 needs dense marginal summary statistics from one biologically
# defined cis/local region.  The cleaning pipeline's *.cis.gz file is therefore
# authoritative when it exists; its observed bounds also preserve the exact
# flank used by gwas_post.sh instead of silently replacing it with a 20-Mb
# lead-SNP window.
mrlink2_file_bounds <- function(file, preferred_chr=NA_character_){
  empty<-tibble(chr=character(),start=numeric(),end=numeric(),n_region_variants=integer())
  if(length(file)!=1L||is.na(file)||!file.exists(file)||file.size(file)<=0)return(empty)
  con<-if(grepl("\\.gz$",file,ignore.case=TRUE))gzfile(file,"rt")else base::file(file,"rt")
  hdr<-tryCatch(readLines(con,n=1,warn=FALSE),error=function(e)"")
  try(close(con),silent=TRUE)
  if(!length(hdr)||!nzchar(hdr))return(empty)
  nms<-strsplit(hdr,"\t",fixed=TRUE)[[1]]
  cchr<-.pick_col(nms,c("CHR","CHROM","CHROMOSOME"));cpos<-.pick_col(nms,c("POS","BP","POSITION","BASE_PAIR_LOCATION"))
  if(is.na(cchr)||is.na(cpos))return(empty)
  d<-tryCatch(data.table::fread(file,select=unique(c(cchr,cpos)),showProgress=FALSE,check.names=FALSE),error=function(e)NULL)
  if(is.null(d)||!nrow(d))return(empty)
  z<-tibble(chr=str_remove(as.character(d[[cchr]]),regex("^chr",ignore_case=TRUE)),
            pos=suppressWarnings(as.numeric(d[[cpos]])))|>
    filter(!is.na(chr),nzchar(chr),is.finite(pos),pos>=1)
  if(!nrow(z))return(empty)
  pref<-str_remove(as.character(preferred_chr)[1],regex("^chr",ignore_case=TRUE))
  if(!is.na(pref)&&nzchar(pref)&&any(z$chr==pref))z<-z|>filter(chr==pref)
  z|>group_by(chr)|>summarise(start=floor(min(pos)),end=ceiling(max(pos)),n_region_variants=n(),.groups="drop")|>
    arrange(desc(n_region_variants),start)|>slice(1)
}

build_mrlink2_jobs <- function(iv_list,annotation,layer,ygfile){
  primary<-if(layer=="protein")"cis"else"local"
  fallback_bp<-as.numeric(Sys.getenv(if(layer=="protein")"C2_CIS_WINDOW_BP"else"C2_LOCAL_WINDOW_BP",unset="1000000"))
  jobs<-imap_dfr(iv_list,function(obj,feature){
    iv<-obj$instruments
    if(!nrow(iv))return(tibble())
    lead_pool<-iv|>filter(analysis==primary,is.finite(POS),!is.na(CHR),is.finite(P))
    # MR-link-2 is a cis/local method.  Do not relabel a trans/distal-only
    # exposure as cis merely to force a job through the runner.
    if(!nrow(lead_pool))return(tibble())
    lead<-lead_pool|>arrange(P)|>slice(1)
    cis_file<-obj$files$cis
    use_cis<-length(cis_file)==1L&&!is.na(cis_file)&&file.exists(cis_file)&&file.size(cis_file)>0
    bounds<-if(use_cis)mrlink2_file_bounds(cis_file,lead$CHR)else tibble()
    if(nrow(bounds)){
      qf<-cis_file;chr<-bounds$chr[[1]];start<-bounds$start[[1]];end<-bounds$end[[1]]
      source<-"cis.gz observed bounds";n_region<-bounds$n_region_variants[[1]]
    }else{
      qf<-obj$files$full
      if(length(qf)!=1L||is.na(qf)||!file.exists(qf))return(tibble())
      chr<-as.character(lead$CHR[[1]]);start<-max(1,lead$POS[[1]]-fallback_bp);end<-lead$POS[[1]]+fallback_bp
      source<-paste0("full summary; ",primary," lead +/-",format(fallback_bp,scientific=FALSE)," bp fallback")
      n_region<-NA_integer_
      if(layer=="protein"&&!is.null(annotation)&&nrow(annotation)){
        a<-annotation|>filter(.data$feature==.env$feature)|>slice(1)
        if(nrow(a)&&!is.na(a$chr[[1]])&&is.finite(a$start[[1]])&&is.finite(a$end[[1]])){
          chr<-as.character(a$chr[[1]]);start<-max(1,a$start[[1]]-fallback_bp);end<-a$end[[1]]+fallback_bp
          source<-paste0("full summary; gene +/-",format(fallback_bp,scientific=FALSE)," bp fallback")
        }
      }
    }
    tibble(omics=layer,trait=feature,exposure=qf,outcome=ygfile,
      region=paste0(chr,":",floor(start),"-",ceiling(end)),region_source=source,
      region_variants=n_region,lead_snp=lead$SNP[[1]],lead_p=lead$P[[1]])
  })
  if(!nrow(jobs)||!"exposure"%in%names(jobs))return(tibble())
  jobs|>filter(!is.na(exposure),file.exists(exposure))
}

venn_membership_counts <- function(a,b,c){
  u<-sort(unique(c(a,b,c)));z<-tibble(id=u,A=u%in%a,B=u%in%b,C=u%in%c)
  tibble(region=c("A","B","C","AB","AC","BC","ABC"),
         n=c(sum(z$A&!z$B&!z$C),sum(!z$A&z$B&!z$C),sum(!z$A&!z$B&z$C),
             sum(z$A&z$B&!z$C),sum(z$A&!z$B&z$C),sum(!z$A&z$B&z$C),sum(z$A&z$B&z$C)))
}

plot_venn_evidence <- function(a,b,c,names3,title,highlight="A"){
  # Circle identity is fixed across panels: local/cis = left, distal/trans =
  # top, observational = right.  Only the panel's reference set changes fill.
  theta<-seq(0,2*pi,length.out=240)
  cent<-tibble(set=c("A","B","C"),cx=c(-.42,0,.42),cy=c(-.15,.33,-.15))
  circ<-cent|>group_by(set)|>group_modify(~tibble(x=.x$cx+.66*cos(theta),y=.x$cy+.66*sin(theta)))
  pos<-tibble(region=c("A","B","C","AB","AC","BC","ABC"),
    x=c(-.72,0,.72,-.25,0,.25,0),y=c(-.30,.70,-.30,.19,-.40,.19,.02))|>
    left_join(venn_membership_counts(a,b,c),by="region")
  lab<-tibble(x=c(-.88,0,.88),y=c(-.94,1.08,-.94),
    label=paste0(names3,"\n(n=",c(length(unique(a)),length(unique(b)),length(unique(c))),")"))
  ggplot(circ,aes(x,y,group=set,color=set))+
    geom_polygon(data=circ|>filter(set==highlight),aes(x,y,group=set),inherit.aes=FALSE,
                 fill="#FFD54F",color=NA,alpha=.88)+
    geom_path(linewidth=1.25)+geom_text(data=pos,aes(x,y,label=n),inherit.aes=FALSE,fontface="bold",size=4)+
    geom_text(data=lab,aes(x,y,label=label),inherit.aes=FALSE,fontface="bold",size=3)+
    scale_color_manual(values=c(A="#D95F02",B="#1B9E77",C="#4C78A8"),guide="none")+
    coord_equal(xlim=c(-1.30,1.30),ylim=c(-1.10,1.20),clip="off")+theme_void()+
    labs(title=title,subtitle="Yellow circle = reference set")+
    theme(plot.title=element_text(face="bold",size=11,hjust=.5),
          plot.subtitle=element_text(face="bold",size=8.5,hjust=.5,color="#8A5A00"),
          plot.margin=margin(5,16,5,16))
}

plot_c2_fig1 <- function(mr,assoc,layer){
  local_name<-if(layer=="protein")"cis"else"local";distal_name<-if(layer=="protein")"trans"else"distal"
  mr0<-mr|>filter(is.finite(pval));obs_sig<-assoc|>filter(is.finite(p.value),p.value<.05/max(1,nrow(assoc)))|>pull(term)
  local<-mr0|>filter(analysis==local_name);distal<-mr0|>filter(analysis==distal_name)
  local_sig<-local|>filter(FDR_analysis<.05)|>pull(exposure);distal_sig<-distal|>filter(FDR_analysis<.05)|>pull(exposure)
  topn<-min(TOP_FOREST,max(10L,floor(sqrt(max(1,nrow(mr0)))*2)))
  top_local<-local|>slice_min(pval,n=topn,with_ties=FALSE)|>pull(exposure)
  top_distal<-distal|>slice_min(pval,n=topn,with_ties=FALSE)|>pull(exposure)
  top_obs<-assoc|>filter(is.finite(p.value))|>slice_min(p.value,n=topn,with_ties=FALSE)|>pull(term)
  va<-plot_venn_evidence(top_local,distal_sig,obs_sig,
    c(paste0("Top ",local_name),paste0(distal_name," FDR"),"Observational Bonf."),paste0("a. Anchored on ",local_name," MR"),"A")
  vb<-plot_venn_evidence(local_sig,top_distal,obs_sig,
    c(paste0(local_name," FDR"),paste0("Top ",distal_name),"Observational Bonf."),paste0("b. Anchored on ",distal_name," MR"),"B")
  vc<-plot_venn_evidence(local_sig,distal_sig,top_obs,
    c(paste0(local_name," FDR"),paste0(distal_name," FDR"),"Top observational"),"c. Anchored on observational evidence","C")
  top_ids<-unique(c(head(top_local,12),head(top_distal,12),head(top_obs,12)))
  ev<-bind_rows(
    assoc|>filter(term%in%top_ids)|>transmute(exposure=term,evidence="Observational",b=beta,se=std.error,p=p.value),
    mr0|>filter(exposure%in%top_ids,analysis%in%c(local_name,distal_name))|>transmute(exposure,evidence=paste0("MR (",analysis,")"),b,se,p=pval))|>
    complete(exposure=top_ids,evidence=c("Observational",paste0("MR (",local_name,")"),paste0("MR (",distal_name,")")))|>
    mutate(available=is.finite(b)&is.finite(se),lo=b-1.96*se,hi=b+1.96*se,b_plot=ifelse(available,b,0),
           exposure=factor(exposure,levels=rev(top_ids)))
  pforest<-ggplot(ev,aes(b_plot,exposure))+geom_vline(xintercept=0,color="grey60")+
    geom_errorbarh(data=ev|>filter(available),aes(xmin=lo,xmax=hi),height=0,color="grey55")+
    geom_point(data=ev|>filter(available),aes(color=p<.05),size=2)+
    geom_point(data=ev|>filter(!available),shape=4,size=2.2,color="grey72")+
    facet_grid(.~evidence,scales="free_x")+scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="#4C78A8"),guide="none")+
    labs(title="d. Effect estimates; × denotes no valid instrument result",x="Effect per 1-SD higher omic trait",y=NULL)+forest_theme(8)
  best<-mr0|>group_by(exposure)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()|>
    left_join(assoc|>select(exposure=term,obs_beta=beta,obs_se=std.error,obs_p=p.value),by="exposure")|>
    mutate(concordance=case_when(FDR_all<.05&sign(b)==sign(obs_beta)~"FDR-significant, concordant",
      FDR_all<.05&sign(b)!=sign(obs_beta)~"FDR-significant, discordant",TRUE~"No FDR MR support"))
  counts<-bind_rows(tibble(anchor=local_name,venn_membership_counts(top_local,distal_sig,obs_sig)),
                    tibble(anchor=distal_name,venn_membership_counts(local_sig,top_distal,obs_sig)),
                    tibble(anchor="observational",venn_membership_counts(local_sig,distal_sig,top_obs)))
  composed<-(va|vb|vc)/pforest+plot_layout(heights=c(.82,1.18))+
    plot_annotation(caption=paste0(toupper(substr(local_name,1,1)),substr(local_name,2,nchar(local_name)),
      " is always left; ",distal_name," is always top; observational is always right. Yellow = panel reference set."),
                    theme=theme(plot.caption=element_text(size=8,color="grey40",hjust=1)))
  list(plot=composed,bars=counts,forest=ev,best=best)
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
  ev<-bind_rows(
    assoc|>transmute(exposure=term,evidence="Observational",effect=beta,p=p.value),
    mr|>filter(analysis%in%c(local_name,distal_name))|>
      transmute(exposure,evidence=ifelse(analysis==local_name,paste0("MR: ",local_name),paste0("MR: ",distal_name)),effect=b,p=pval))|>
    filter(is.finite(p),is.finite(effect))|>group_by(evidence)|>mutate(q=p.adjust(p,"BH"))|>ungroup()
  top_ev<-ev|>group_by(exposure)|>summarise(best_p=min(p),.groups="drop")|>slice_min(best_p,n=28,with_ties=FALSE)|>pull(exposure)
  ev<-ev|>filter(exposure%in%top_ev)|>mutate(score=stable_neglog10_p(p),
      signed_score=sign(effect)*pmin(score,12),significant=q<.05,
      evidence=factor(evidence,levels=c("Observational",paste0("MR: ",local_name),paste0("MR: ",distal_name))),
      exposure=factor(exposure,levels=rev(top_ev)))
  pC<-if(!nrow(ev))blank_plot("c. Evidence matrix")else ggplot(ev,aes(evidence,exposure))+
    geom_tile(aes(fill=signed_score),color="white",linewidth=.35)+
    geom_point(data=ev|>filter(significant),shape=8,size=2.1,color="black")+
    scale_fill_gradient2(low="#3F78A8",mid="white",high="#C86B4A",midpoint=0,limits=c(-12,12),
      name="signed\n-log10(P)")+
    labs(title="c. Signed cross-evidence matrix",
      subtitle="Blue = protective; orange = risk-increasing; asterisk = within-analysis FDR < 0.05; scale capped at 12",
      x=NULL,y=NULL)+theme_5c(8)+theme(legend.position="right")
  list(plot=(pA|pB)/pC+plot_layout(heights=c(1,.9)),wide=wide)
}

plot_qtl_variance <- function(mr,layer){
  cls<-if(layer=="protein")c("cis","trans")else c("local","distal")
  ps<-map(cls,function(cc){
    d<-mr|>filter(analysis==cc,is.finite(r2_median),n_IV>0)|>
      mutate(r2_med=100*r2_median,r2_lo=100*r2_q25,r2_hi=100*r2_q75,r2_p90_plot=100*r2_p90)|>
      slice_max(r2_p90_plot,n=20,with_ties=FALSE)|>arrange(r2_med)|>mutate(exposure=factor(exposure,levels=exposure))
    if(!nrow(d))return(blank_plot(paste0(toupper(substr(cc,1,1)),substr(cc,2,nchar(cc))," instruments"),"No valid instrument set"))
    ggplot(d,aes(r2_med,exposure))+
      geom_segment(aes(x=r2_med,xend=r2_p90_plot,yend=exposure),color="grey72",linewidth=.7)+
      geom_errorbarh(aes(xmin=r2_lo,xmax=r2_hi),height=.16,color=ifelse(cc%in%c("cis","local"),"#D95F02","#12AEB5"),linewidth=.8)+
      geom_point(aes(size=pmin(n_IV,500)),color=ifelse(cc%in%c("cis","local"),"#D95F02","#12AEB5"))+
      geom_text(aes(x=r2_p90_plot,label=sprintf("median %.3f%%; K=%d",r2_med,n_IV)),hjust=-.05,size=2.35,fontface="bold")+
      scale_x_continuous(expand=expansion(mult=c(.02,.34)))+scale_size_continuous(range=c(1.6,4.8),name="IV count")+
      labs(title=paste0(cc," QTL instruments"),subtitle="Per-IV partial R²: point = median; interval = IQR; grey tail = 90th percentile",
           x="Per-instrument partial phenotypic R² (%)",y=NULL)+theme_5c(8)
  })
  (ps[[1]]|ps[[2]])+
    plot_annotation(title="QTL instrument-level variance architecture",
      subtitle="This figure does not sum R² and does not use 1 - product(1 - R²); values cannot saturate merely because K is large",
      caption="Instrument-level R² distribution",
      theme=theme(plot.title=element_text(face="bold",size=15),
                  plot.subtitle=element_text(size=10,color="grey30"),
                  plot.caption=element_text(size=8,color="grey40",hjust=1)))
}

# C2 has one output contract.  Standard figure names are overwritten whenever
# figures are regenerated from a valid numerical cache.
save_c2_plot <- function(p,file,width,height,outdir){
  save_plot(p,file,width,height,outdir=outdir)
}

plot_c2_fig4 <- function(mr,layer){
  for(nm in c("egger_intercept_p","steiger_support","steiger_support_fraction","steiger_n","median_F","r2_median"))if(!nm%in%names(mr))mr[[nm]]<-NA
  d<-mr|>filter(is.finite(pval),is.finite(median_F),is.finite(r2_median))|>
    mutate(class=analysis,significant=FDR_analysis<.05,label=ifelse(min_rank(pval)<=8,exposure,NA_character_),
      q_score=stable_neglog10_p(Q_p),egger_score=stable_neglog10_p(egger_intercept_p),
      q_plot=compress_extreme_tail(q_score,10,4),egger_plot=compress_extreme_tail(egger_score,10,4))
  if(!nrow(d))return(blank_plot("MR sensitivity and instrument architecture","No complete instrument metrics"))
  pa<-ggplot(d,aes(median_F,100*r2_median,color=class,size=pmin(n_IV,500)))+geom_vline(xintercept=10,linetype=2,color="grey55")+
    geom_point(alpha=.78)+ggrepel::geom_text_repel(aes(label=label),size=2.35,seed=31,max.overlaps=10,na.rm=TRUE)+scale_x_log10()+
    labs(title="a. Instrument strength and typical variance",subtitle="Median per-IV metrics avoid the many-IV product-to-one artefact",
         x="Median F statistic (log scale)",y="Median per-IV partial R² (%)",color=NULL,size="IV count")+theme_5c(9)
  ds<-d|>filter(is.finite(q_plot),is.finite(egger_plot))
  pb<-if(!nrow(ds))blank_plot("b. Heterogeneity and directional pleiotropy","Egger diagnostics require at least three valid IVs")else {
    sx<-compressed_tail_scale(ds$q_score,10,4);sy<-compressed_tail_scale(ds$egger_score,10,4)
    ggplot(ds,aes(q_plot,egger_plot,color=class))+geom_vline(xintercept=-log10(.05),linetype=2,color="grey60")+
      geom_hline(yintercept=-log10(.05),linetype=2,color="grey60")+geom_point(aes(size=pmin(n_IV,500)),alpha=.78)+
      scale_x_continuous(breaks=sx$breaks,labels=sx$labels)+scale_y_continuous(breaks=sy$breaks,labels=sy$labels)+
      labs(title="b. Heterogeneity and directional pleiotropy",subtitle="Axes show -log10(P); values above 10 are monotonically compressed",
           x=expression(-log[10](P[Cochran~Q])),y=expression(-log[10](P[Egger~intercept])),color=NULL,size="IV count")+theme_5c(9)
  }
  dz<-d|>filter(is.finite(steiger_support_fraction))
  pc<-if(!nrow(dz))blank_plot("c. Steiger directionality support","No per-IV directionality comparison was estimable")else
    ggplot(dz,aes(class,steiger_support_fraction,color=class))+geom_hline(yintercept=.5,linetype=2,color="grey55")+
      geom_violin(aes(fill=class),alpha=.12,color=NA,trim=FALSE)+geom_boxplot(width=.18,outlier.shape=NA,alpha=.55)+
      geom_jitter(aes(size=pmin(steiger_n,500)),width=.08,alpha=.48)+scale_y_continuous(limits=c(0,1),labels=label_percent())+
      labs(title="c. Steiger directionality support",subtitle="Fraction of harmonized IVs with R²(exposure) > R²(outcome); >50% supports exposure → outcome",
           x=NULL,y="Supporting IV fraction",color=NULL,fill=NULL,size="IV count")+theme_5c(9)+theme(legend.position="none")
  (pa|pb)/pc+plot_layout(heights=c(1,.72))
}

read_mrlink2_results <- function(rawdir,mode=RUN_MRlink2){
  linkdir<-file.path(rawdir,"mrlink2",str_to_lower(mode))
  f<-file.path(linkdir,"mrlink2.all.tsv");s<-file.path(linkdir,"mrlink2.status.tsv")
  ans<-list(results=tibble(),status=tibble(),mode=mode)
  if(file.exists(f)&&file.size(f)>0)ans$results<-tryCatch(as_tibble(data.table::fread(f,showProgress=FALSE,check.names=FALSE)),error=function(e)tibble())
  if(file.exists(s)&&file.size(s)>0)ans$status<-tryCatch(as_tibble(data.table::fread(s,showProgress=FALSE,check.names=FALSE)),error=function(e)tibble())
  ans
}

plot_mrlink2_results <- function(x,outdir){
  d<-x$results;st<-x$status
  if(nrow(d)){
    ac<-pick_col_local(names(d),c("^ALPHA$"));sc<-pick_col_local(names(d),c("^SE\\(ALPHA\\)$","^SE_ALPHA$"));pc<-pick_col_local(names(d),c("^P\\(ALPHA\\)$","^P_ALPHA$"))
    yc<-pick_col_local(names(d),c("^SIGMA_Y$","^SIGMAY$"));ypc<-pick_col_local(names(d),c("^P\\(SIGMA_Y\\)$","^P_SIGMA_Y$","^PSIGMAY$"))
    nc<-pick_col_local(names(d),c("^M_SNPS_OVERLAP$","^MOVERLAP$","^N_SNPS$"));rc<-pick_col_local(names(d),c("^REGION$"))
    if(all(!is.na(c(ac,sc,pc)))){
      num_or_na<-function(nm)if(is.na(nm))rep(NA_real_,nrow(d))else suppressWarnings(as.numeric(d[[nm]]))
      char_or_na<-function(nm)if(is.na(nm))rep(NA_character_,nrow(d))else as.character(d[[nm]])
      trait_col<-pick_col_local(names(d),c("^TRAIT$","^EXPOSURE$"))
      z_all<-tibble(trait=char_or_na(trait_col),region=char_or_na(rc),alpha=num_or_na(ac),se=num_or_na(sc),p=num_or_na(pc),
                    sigma_y=num_or_na(yc),p_sigma_y=num_or_na(ypc),n_overlap=num_or_na(nc))|>
        filter(!is.na(trait),trait!="",is.finite(alpha),is.finite(se),se>0,is.finite(p))|>
        group_by(trait)|>slice_min(p,n=1,with_ties=FALSE)|>ungroup()|>arrange(p)
      if(!nrow(z_all)){
        pa<-blank_plot("a. MR-link-2 regional causal effects","Mapped columns contained no finite estimates")
        pb<-blank_plot("b. MR-link-2 causal evidence","Mapped columns contained no finite estimates")
        pc<-blank_plot("c. Regional pleiotropy audit","Mapped columns contained no finite estimates")
      }else{
        z<-z_all|>slice_head(n=32)|>mutate(lo=alpha-1.96*se,hi=alpha+1.96*se)
        robust_lim<-as.numeric(quantile(abs(c(z$lo,z$hi)),.98,na.rm=TRUE,names=FALSE))
        if(!is.finite(robust_lim)||robust_lim<=0)robust_lim<-max(abs(c(z$lo,z$hi)),na.rm=TRUE)
        z<-z|>mutate(lo_plot=pmax(lo,-robust_lim),hi_plot=pmin(hi,robust_lim),clipped=lo< -robust_lim|hi>robust_lim,
                     trait=factor(trait,levels=rev(trait)))
        pa<-ggplot(z,aes(alpha,trait))+geom_vline(xintercept=0,color="grey60")+
          geom_errorbarh(aes(xmin=lo_plot,xmax=hi_plot),height=.08)+geom_point(aes(color=p<.05/max(1,nrow(z))),size=2)+
          geom_point(data=z|>filter(clipped),aes(x=ifelse(alpha>=0,robust_lim,-robust_lim),y=trait),shape=17,size=2.1,color="grey35",inherit.aes=FALSE)+
          scale_color_manual(values=c(`TRUE`="#D95F02",`FALSE`="#4C78A8"),guide="none")+
          coord_cartesian(xlim=c(-robust_lim,robust_lim))+
          labs(title="a. MR-link-2 regional causal effects",subtitle="LD-aware estimate per trait; triangles mark CIs clipped at the robust display limit",
               x=expression(alpha),y=NULL)+forest_theme(8)
        zb<-z_all|>mutate(score=stable_neglog10_p(p),score_plot=compress_extreme_tail(score,10,4),label=ifelse(min_rank(p)<=14,trait,NA_character_))
        sy<-compressed_tail_scale(zb$score,10,4)
        pb<-ggplot(zb,aes(alpha,score_plot))+geom_vline(xintercept=0,color="grey60")+
          geom_hline(yintercept=-log10(.05/max(1,nrow(zb))),linetype=2,color="grey55")+
          geom_point(color="#6A3D9A",size=2,alpha=.78)+
          ggrepel::geom_text_repel(aes(label=label),size=2.45,seed=41,max.overlaps=15,na.rm=TRUE)+
          scale_y_continuous(breaks=sy$breaks,labels=sy$labels)+
          labs(title="b. MR-link-2 causal evidence",subtitle="-log10(P) above 10 is monotonically compressed",
               x=expression(alpha),y=expression(-log[10](P[alpha])))+theme_5c(9)
        audit_txt<-if(nrow(st)&&"status"%in%names(st))paste(paste0(names(table(st$status)),"=",as.integer(table(st$status))),collapse="; ")else"status file unavailable"
        zc<-z_all|>filter(is.finite(sigma_y))|>
          mutate(label=ifelse(min_rank(p)<=12,trait,NA_character_),
                 pleiotropy=ifelse(is.finite(p_sigma_y)&p_sigma_y<.05,"P(sigma_y) < 0.05","No sigma_y support"),
                 n_overlap=ifelse(is.finite(n_overlap),pmax(1,n_overlap),1))
        pc<-if(!nrow(zc))blank_plot("c. Regional pleiotropy audit","sigma_y fields were not returned by this MR-link-2 build")else
          ggplot(zc,aes(alpha,sigma_y,color=pleiotropy,size=n_overlap))+geom_hline(yintercept=0,color="grey70")+geom_vline(xintercept=0,color="grey70")+
            geom_point(alpha=.75)+ggrepel::geom_text_repel(aes(label=label),size=2.4,seed=42,max.overlaps=12,na.rm=TRUE)+
            scale_color_manual(values=c(`P(sigma_y) < 0.05`="#D95F02",`No sigma_y support`="#4C78A8"))+
            scale_size_continuous(range=c(1.5,5),name="Overlapping SNPs")+
            labs(title="c. Regional causal effect and residual pleiotropy",subtitle=paste0("Execution audit: ",audit_txt),
                 x=expression(alpha),y=expression(sigma[y]),color=NULL)+theme_5c(9)
      }
    }else{
      pa<-blank_plot("a. MR-link-2 regional causal effects","Output columns could not be mapped")
      pb<-blank_plot("b. MR-link-2 causal evidence","Output columns could not be mapped")
      pc<-blank_plot("c. Regional pleiotropy audit","Output columns could not be mapped")
    }
  } else {
    why<-if(identical(x$mode,"None"))"RUN_MRlink2=None; no regional jobs were requested" else
      paste0("RUN_MRlink2=",x$mode%||%"unknown","; no completed regional estimate")
    pa<-blank_plot("MR-link-2 execution status",why)
    pb<-blank_plot("Candidate selection",why)
    pc<-blank_plot("Regional pleiotropy audit",why)
  }
  save_plot(pa/(pb|pc)+plot_layout(heights=c(1.2,.8)),"c2.Fig5.mrlink2.png",14.5,11,outdir=outdir)
}


# 🚩 DANDELION helpers
pick_col_local <- function(nms,patterns){
  nl<-toupper(nms)
  for(p in patterns){i<-grep(p,nl,perl=TRUE)[1];if(!is.na(i))return(nms[[i]])}
  NA_character_
}

find_magma_gene_file <- function(ygfile){
  explicit<-Sys.getenv("C2_MAGMA_GENE_FILE",unset="");if(nzchar(explicit)&&file.exists(explicit))return(normalizePath(explicit,winslash="/",mustWork=FALSE))
  trait<-sub("\\.gz$","",basename(ygfile),ignore.case=TRUE)
  # gwas_post contract: <project>/common/<trait>/{gwas,magma}.
  magma_dir<-gwas_magma_dir_from_clean_file(ygfile)
  explicit_candidates<-file.path(magma_dir,paste0(trait,".genes.out"))
  discovered<-if(!is.na(magma_dir)&&dir.exists(magma_dir))list.files(magma_dir,pattern="(genes\\.out|gene\\.out)$",full.names=TRUE,ignore.case=TRUE)else character()
  fs<-unique(c(explicit_candidates,discovered));fs<-fs[file.exists(fs)&!grepl("(gsa|sets)\\.out$",fs,ignore.case=TRUE)]
  if(length(fs))return(normalizePath(fs[[which.max(file.info(fs)$mtime)]],winslash="/",mustWork=FALSE))
  NA_character_
}

find_dandelion_gene_file <- function(ygfile){
  explicit<-Sys.getenv("C2_DANDELION_GENE_P_FILE",unset="")
  if(nzchar(explicit)){
    if(!file.exists(explicit)||file.size(explicit)<=0)stop("DANDELION gene-level disease file is missing: ",explicit,call.=FALSE)
    return(list(path=normalizePath(explicit,winslash="/",mustWork=FALSE),
                evidence_type="user-supplied gene-level disease association (WES burden preferred)",
                primary_eligible=TRUE,analysis_class="primary candidate"))
  }
  if(!DANDELION_ALLOW_MAGMA)return(list(path=NA_character_,
    evidence_type="not run: explicit gene-level rare-variant/WES input required",
    primary_eligible=FALSE,analysis_class="not run"))
  magma<-find_magma_gene_file(ygfile)
  list(path=magma,evidence_type="MAGMA gene-level GWAS (adapted DANDELION; not WES burden)",
       primary_eligible=FALSE,analysis_class="MAGMA sensitivity")
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

read_gene_level_p <- function(file,target_symbols){
  if(is.na(file)||!file.exists(file))return(numeric())
  d<-data.table::fread(file,showProgress=FALSE,check.names=FALSE);nms<-names(d)
  gc<-pick_col_local(nms,c("^GENE$","^GENE_ID$","^GENEID$","^ID$","^SYMBOL$","^GENE_SYMBOL$"))
  pc<-pick_col_local(nms,c("^P$","^PVAL$","^P_VALUE$","^PVALUE$","^PVAL_BURDEN$","^P_VALUE_BURDEN$","^BURDEN_P$","^P_BURDEN$"))
  if(is.na(gc)||is.na(pc))stop("Gene-level disease file needs a gene identifier column and P column: ",file,call.=FALSE)
  gn<-map_magma_gene_ids(d[[gc]],target_symbols);p<-suppressWarnings(as.numeric(d[[pc]]));ok<-!is.na(gn)&nzchar(gn)&is.finite(p)&p>0&p<=1
  p<-p[ok];names(p)<-gn[ok];p<-tapply(p,names(p),min,na.rm=TRUE);as.numeric(p)|>setNames(names(p))
}

read_lead_file <- function(file,ygfile){
  if(is.na(file)||!file.exists(file))return(tibble())
  d<-if(grepl("jma\\.cojo$",file,ignore.case=TRUE))match_GRCH_table(file,ygfile) else data.table::fread(file,showProgress=FALSE,check.names=FALSE,fill=TRUE);nms<-names(d)
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

find_disease_lead_candidates <- function(ygfile){
  explicit<-Sys.getenv("C2_DANDELION_SNP_FILE",unset="");cand<-character()
  if(nzchar(explicit))cand<-c(cand,explicit)
  d<-dirname(ygfile);cand<-c(cand,list.files(d,pattern="(jma\\.cojo$|\\.clumped$|lead.*\\.(txt|tsv|csv)$|indep.*\\.(txt|tsv|csv)$)",full.names=TRUE,ignore.case=TRUE))
  unique(cand[file.exists(cand)])
}

find_disease_leads <- function(ygfile){
  cand<-find_disease_lead_candidates(ygfile)
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
    if(is.na(f)||!file.exists(f))return(list(p=setNames(rep(NA_real_,length(snps)),snps),file=NA_character_,n=0L,error="full QTL missing"))
    err<-NA_character_;z<-tryCatch(read_sumstat_snps(f,snps,N_DEFAULT),error=function(e){err<<-conditionMessage(e);tibble()})
    p<-rep(NA_real_,length(snps));names(p)<-snps
    if(nrow(z))p[z$SNP]<-z$P
    list(p=p,file=f,n=sum(is.finite(p)),error=err)
  })
  M<-do.call(rbind,lapply(rows,`[[`,"p"));rownames(M)<-gene2;colnames(M)<-snps;storage.mode(M)<-"double"
  attr(M,"qtl_audit")<-tibble(gene2=gene2,full_qtl=vapply(rows,function(x)x$file%||%NA_character_,character(1)),
    n_disease_snps_found=vapply(rows,function(x)as.integer(x$n%||%0L),integer(1)),
    read_error=vapply(rows,function(x)x$error%||%NA_character_,character(1)))
  M
}

# Fall back to assayed-protein annotations when no gene annotation is supplied.
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
  if(!"dpg_class"%in%names(tg))tg<-tg|>mutate(dpg_class="DANDELION prioritized target")
  sensitivity_note<-if("gene_evidence_type"%in%names(tg)&&any(str_detect(str_to_lower(tg$gene_evidence_type),"magma|adapted"),na.rm=TRUE))
    "MAGMA-adapted sensitivity analysis; not counted as primary rare-variant DANDELION evidence" else
    "Gene-level input and selection diagnostics are reported in the input-audit table"
  # Panel A combines target significance and hub evidence. A target supported
  # by several independent distal disease loci is the biologically interesting
  # convergence pattern; it is not a generic gene1 -> gene2 -> gene3 chain.
  top<-tg|>arrange(DANDELION_p,desc(n_distal_loci))|>slice_head(n=24)|>
    mutate(path_score=stable_neglog10_p(DANDELION_p),path_plot=compress_extreme_tail(path_score,10,4),
           disease_score=stable_neglog10_p(gene_level_p),gene2=factor(gene2,levels=rev(gene2)))
  sxA<-compressed_tail_scale(top$path_score,10,4)
  pA<-ggplot(top,aes(path_plot,gene2))+
    geom_segment(aes(x=0,xend=path_plot,yend=gene2),color="grey76")+
    geom_point(aes(size=n_distal_loci,color=disease_score,shape=dpg_class),alpha=.88)+
    geom_text(aes(label=ifelse(n_distal_loci>1,paste0(n_distal_loci," loci"),"")),nudge_x=.22,hjust=0,size=2.45,fontface="bold")+
    scale_x_continuous(breaks=sxA$breaks,labels=sxA$labels,expand=expansion(mult=c(0,.18)))+
    scale_color_viridis_c(option="C",end=.9)+
    labs(title="a. Prioritized targets and convergent hubs",subtitle=sensitivity_note,
         x=expression(-log[10](DANDELION~P)),y=NULL,size="Upstream loci",shape="Target class",color=expression(-log[10](gene-level~P)))+theme_5c(9)

  # The 0--10 range is kept linear because this is where evidential grades are
  # distinguishable. More extreme values remain ordered but consume less ink.
  evid<-pd|>filter(is.finite(trans_p),is.finite(gene_level_p),is.finite(DANDELION_p))|>
    mutate(trans_score=stable_neglog10_p(trans_p),disease_score=stable_neglog10_p(gene_level_p),path_score=stable_neglog10_p(DANDELION_p),
           trans_plot=compress_extreme_tail(trans_score,10,4),disease_plot=compress_extreme_tail(disease_score,10,4))
  if(nrow(evid)){
    label_rows<-evid|>group_by(gene2)|>slice_min(DANDELION_p,n=1,with_ties=FALSE)|>ungroup()|>arrange(DANDELION_p)|>slice_head(n=18)
    sxB<-compressed_tail_scale(evid$trans_score,10,4);syB<-compressed_tail_scale(evid$disease_score,10,4)
    pB<-ggplot(evid,aes(trans_plot,disease_plot))+
      geom_point(aes(size=path_score,color=path_score),alpha=.68)+
      ggrepel::geom_text_repel(data=label_rows,aes(label=gene2),size=2.45,fontface="bold",seed=12,max.overlaps=Inf)+
      scale_x_continuous(breaks=sxB$breaks,labels=sxB$labels)+scale_y_continuous(breaks=syB$breaks,labels=syB$labels)+
      scale_color_viridis_c(option="C",end=.9)+
      labs(title="b. Trans-regulatory and disease-gene evidence",subtitle="Both axes are linear from 0 to 10 and monotonically compressed thereafter",
           x=expression(-log[10](trans-regulatory~P)),y=expression(-log[10](gene-level~disease~P)),
           size=expression(-log[10](DANDELION~P)),color=expression(-log[10](DANDELION~P)))+theme_5c(9)
  }else pB<-blank_plot("b. Trans-regulatory and disease-gene evidence","No finite pair-level evidence")

  # DANDELION paths link disease loci to trans-regulated disease-proximal genes.
  path<-pd
  if(nrow(path)){
    src<-if("gene1"%in%names(path))as.character(path$gene1)else if("region"%in%names(path))as.character(path$region)else if("rsid"%in%names(path))as.character(path$rsid)else rep(NA_character_,nrow(path))
    fallback_src<-if("rsid"%in%names(path))as.character(path$rsid)else if("region"%in%names(path))as.character(path$region)else rep(NA_character_,nrow(path))
    preferred_targets<-tg|>arrange(DANDELION_p,desc(n_distal_loci))|>slice_head(n=14)|>pull(gene2)
    path<-path|>mutate(source=ifelse(is.na(src)|src=="",fallback_src,src),target=as.character(gene2),
      path_score=stable_neglog10_p(DANDELION_p),trans_score=stable_neglog10_p(trans_p),disease_score=stable_neglog10_p(gene_level_p))|>
      filter(target%in%preferred_targets,!is.na(source),source!="")|>arrange(DANDELION_p)|>
      distinct(source,target,.keep_all=TRUE)|>group_by(target)|>slice_head(n=3)|>ungroup()|>slice_head(n=30)
  }
  if(!nrow(path)){
    pC<-blank_plot("c. Distal-regulatory disease-driver network","No interpretable upstream-to-target path")
  }else{
    target_order<-unique(path$target)
    target_nodes<-tibble(target=target_order,target_y=rev(seq_along(target_order)))|>
      left_join(tg|>transmute(target=as.character(gene2),n_distal_loci,DANDELION_p),by="target")|>
      mutate(path_score=stable_neglog10_p(DANDELION_p),n_distal_loci=replace_na(as.numeric(n_distal_loci),1))
    source_order<-unique(path$source)
    source_nodes<-tibble(source=source_order,source_y=seq(max(target_nodes$target_y),min(target_nodes$target_y),length.out=length(source_order)))
    disease_y<-mean(range(target_nodes$target_y))
    edge_up<-path|>left_join(source_nodes,by="source")|>left_join(target_nodes|>select(target,target_y),by="target")|>
      mutate(edge_strength=compress_extreme_tail(trans_score,10,4))
    edge_down<-target_nodes|>mutate(disease_y=disease_y,edge_strength=compress_extreme_tail(path_score,10,4))
    pC<-ggplot()+
      geom_curve(data=edge_up,aes(x=0,y=source_y,xend=1,yend=target_y,alpha=edge_strength),curvature=.08,color="#4C78A8",
                 arrow=grid::arrow(length=grid::unit(.055,"inches"),type="closed"),linewidth=.55)+
      geom_curve(data=edge_down,aes(x=1,y=target_y,xend=2,yend=disease_y,alpha=edge_strength),curvature=.06,color="#E45756",
                 arrow=grid::arrow(length=grid::unit(.055,"inches"),type="closed"),linewidth=.75)+
      geom_point(data=source_nodes,aes(0,source_y),shape=24,size=3.1,fill="#D9D9D9",color="grey25")+
      geom_point(data=target_nodes,aes(1,target_y,size=n_distal_loci,fill=path_score),shape=21,color="grey20")+
      annotate("point",x=2,y=disease_y,shape=22,size=5,fill="#E45756",color="grey20")+
      geom_text(data=source_nodes,aes(0,source_y,label=source),hjust=1.12,size=2.5)+
      geom_text(data=target_nodes,aes(1,target_y,label=target),hjust=-.16,size=2.65,fontface="bold")+
      annotate("text",x=2.08,y=disease_y,label=Y,hjust=0,size=3.1,fontface="bold")+
      annotate("text",x=0,y=max(target_nodes$target_y)+.8,label="Upstream locus / cis gene",fontface="bold",size=3.1)+
      annotate("text",x=1,y=max(target_nodes$target_y)+.8,label="Trans target / driver hub",fontface="bold",size=3.1)+
      annotate("text",x=2,y=max(target_nodes$target_y)+.8,label="Disease",fontface="bold",size=3.1)+
      scale_alpha_continuous(range=c(.18,.82),guide="none")+scale_size_continuous(range=c(2.8,7),name="Upstream loci")+
      scale_fill_viridis_c(option="C",end=.9,name=expression(-log[10](DANDELION~P)))+
      coord_cartesian(xlim=c(-.55,2.35),ylim=c(min(target_nodes$target_y)-.8,max(target_nodes$target_y)+1.25),clip="off")+
      labs(title="c. Distal-regulatory disease-driver network",subtitle="Blue: distal trans regulation; red: target-level integration with disease evidence. Only data-supported hubs are shown.")+theme_void(base_size=9)+
      theme(plot.title=element_text(face="bold",size=12,hjust=0),plot.subtitle=element_text(color="grey35",size=8.5),
            legend.position="right",plot.margin=margin(8,80,8,90))
  }
  save_plot(pA|pB,"c2.Fig6.dandelion.png",17,8.5,outdir=outdir)
  save_plot(pC,"c2.Fig9.dandelion_network.png",15,9,outdir=outdir)
}

plot_dandelion_landscape <- function(evidence,tg,outdir,evidence_type="gene-level disease association"){
  if(!nrow(evidence)){
    save_plot(blank_plot("DANDELION tested-pair landscape","No finite trans-regulatory/gene-level pair was available"),
              "c2.Fig7.dandelion_evidence.png",11,7,outdir=outdir);return(invisible(NULL))
  }
  d<-evidence|>filter(is.finite(trans_p),is.finite(gene_level_p),is.finite(DANDELION_p))|>
    mutate(trans_score=stable_neglog10_p(trans_p),gene_score=stable_neglog10_p(gene_level_p),
           dandelion_score=stable_neglog10_p(DANDELION_p),selected=as.logical(significant))
  labels<-d|>filter(selected)|>arrange(DANDELION_p)|>distinct(gene2,.keep_all=TRUE)|>slice_head(n=14)
  sx<-compressed_tail_scale(d$trans_score,10,4);sy<-compressed_tail_scale(d$gene_score,10,4)
  pA<-ggplot(d,aes(compress_extreme_tail(trans_score,10,4),compress_extreme_tail(gene_score,10,4)))+
    geom_point(data=d|>filter(!selected),color="grey80",alpha=.22,size=.55)+
    geom_point(data=d|>filter(selected),aes(color=dandelion_score),alpha=.88,size=1.65)+
    ggrepel::geom_text_repel(data=labels,aes(label=gene2),size=2.45,fontface="bold",seed=72,max.overlaps=Inf)+
    scale_x_continuous(breaks=sx$breaks,labels=sx$labels)+scale_y_continuous(breaks=sy$breaks,labels=sy$labels)+
    scale_color_viridis_c(option="C",end=.9)+
    labs(title="a. Complete two-component evidence landscape",
      subtitle=paste0("Grey = tested pairs; colour = DANDELION-selected paths. Gene evidence: ",evidence_type),
      x=expression(-log[10](trans-regulatory~P)),y=expression(-log[10](gene-level~disease~P)),
      color=expression(-log[10](DANDELION~P)))+theme_5c(9)

  selected<-d|>filter(selected)
  top_targets<-tg|>arrange(DANDELION_p,desc(n_distal_loci))|>slice_head(n=18)|>pull(gene2)
  top_loci<-selected|>filter(gene2%in%top_targets)|>group_by(rsid)|>summarise(best=min(DANDELION_p),.groups="drop")|>
    arrange(best)|>slice_head(n=25)|>pull(rsid)
  heat<-selected|>filter(gene2%in%top_targets,rsid%in%top_loci)|>
    mutate(gene2=factor(gene2,levels=rev(top_targets)),rsid=factor(rsid,levels=top_loci),score=stable_neglog10_p(DANDELION_p))
  pB<-if(!nrow(heat))blank_plot("b. Disease-locus to DPG map","No selected pair to display") else
    ggplot(heat,aes(rsid,gene2,fill=compress_extreme_tail(score,10,4)))+geom_tile(color="white",linewidth=.25)+
      scale_fill_viridis_c(option="C",end=.9)+
      labs(title="b. Disease-locus to disease-proximal-gene map",
           subtitle="Only selected DANDELION paths are shown; rows/columns are ordered by strongest integrated evidence",
           x="Independent disease SNP",y="Candidate disease-proximal gene",fill=expression(-log[10](P)))+
      theme_5c(8)+theme(axis.text.x=element_text(angle=55,hjust=1,vjust=1))
  save_plot(pA/pB+plot_layout(heights=c(1,1.05)),"c2.Fig7.dandelion_evidence.png",16.5,14,outdir=outdir)
}

plot_dandelion_mr_integration <- function(dan,mr,assoc,outdir){
  tg<-as_tibble(dan$targets%||%tibble())
  if(!nrow(tg)){
    save_plot(blank_plot("DANDELION / MR integration","No DANDELION target was available"),
              "c2.Fig8.dandelion_mr_integration.png",10,6,outdir=outdir);return(tibble())
  }
  mrbest<-mr|>filter(is.finite(pval))|>group_by(exposure,analysis)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()
  obs<-assoc|>filter(is.finite(p.value))|>group_by(term)|>slice_min(p.value,n=1,with_ties=FALSE)|>ungroup()
  if(!"consolidation_eligible"%in%names(tg)){
    primary_input<-if("gene_evidence_type"%in%names(tg))
      !str_detect(str_to_lower(coalesce(tg$gene_evidence_type,"")),"magma|adapted") else FALSE
    tested_n<-suppressWarnings(as.numeric((dan$ptrans_dimensions%||%c(NA_real_))[[1]]))
    target_fraction<-if(is.finite(tested_n)&&tested_n>0)nrow(tg)/tested_n else NA_real_
    broad<-is.finite(target_fraction)&&target_fraction>DANDELION_MAX_TARGET_FRACTION
    tg$consolidation_eligible<-primary_input&!broad
  }
  wide<-tg|>transmute(feature=as.character(gene2),DANDELION_p,n_distal_loci,gene_level_p,
      consolidation_eligible=as.logical(consolidation_eligible))|>
    left_join(mrbest|>select(feature=exposure,analysis,mr_b=b,mr_p=pval,mr_fdr=FDR_all)|>
      pivot_wider(names_from=analysis,values_from=c(mr_b,mr_p,mr_fdr)),by="feature")|>
    left_join(obs|>transmute(feature=term,obs_b=beta,obs_p=p.value),by="feature")
  local_name<-if("cis"%in%mrbest$analysis)"cis"else"local";distal_name<-if("trans"%in%mrbest$analysis)"trans"else"distal"
  for(nm in c(paste0(c("mr_b_","mr_p_","mr_fdr_"),local_name),paste0(c("mr_b_","mr_p_","mr_fdr_"),distal_name),"obs_b","obs_p"))
    if(!nm%in%names(wide))wide[[nm]]<-NA_real_
  patt<-wide|>transmute(feature,
    local_MR=.data[[paste0("mr_fdr_",local_name)]]<.05,
    distal_MR=.data[[paste0("mr_fdr_",distal_name)]]<.05,
    observational=obs_p<.05/max(1,nrow(assoc)))|>
    mutate(across(-feature,~replace_na(.x,FALSE)),independent_support=as.integer(local_MR)+as.integer(distal_MR)+as.integer(observational))
  wide<-wide|>left_join(patt|>select(feature,local_MR,distal_MR,observational,independent_support),by="feature")
  top_features<-wide|>arrange(desc(independent_support),DANDELION_p,desc(n_distal_loci))|>slice_head(n=32)|>pull(feature)
  long<-bind_rows(
    wide|>transmute(feature,evidence="DANDELION",p=DANDELION_p,effect=NA_real_,supported=consolidation_eligible),
    wide|>transmute(feature,evidence=paste0("MR: ",local_name),p=.data[[paste0("mr_p_",local_name)]],effect=.data[[paste0("mr_b_",local_name)]],supported=.data[[paste0("mr_fdr_",local_name)]]<.05),
    wide|>transmute(feature,evidence=paste0("MR: ",distal_name),p=.data[[paste0("mr_p_",distal_name)]],effect=.data[[paste0("mr_b_",distal_name)]],supported=.data[[paste0("mr_fdr_",distal_name)]]<.05),
    wide|>transmute(feature,evidence="Observational",p=obs_p,effect=obs_b,supported=obs_p<.05/max(1,nrow(assoc))))|>
    filter(feature%in%top_features)|>
    mutate(score=stable_neglog10_p(p),direction=case_when(evidence=="DANDELION"~"DPG path",effect>0~"Positive",effect<0~"Inverse",TRUE~"Unavailable"),
           feature=factor(feature,levels=rev(top_features)),
           evidence=factor(evidence,levels=c("DANDELION",paste0("MR: ",local_name),paste0("MR: ",distal_name),"Observational")))
  pA<-ggplot(long|>filter(is.finite(score)),aes(evidence,feature))+
    geom_point(aes(size=compress_extreme_tail(score,10,4),fill=direction,alpha=supported),shape=21,color="grey35")+
    scale_fill_manual(values=c(`DPG path`="#6A3D9A",Positive="#D95F02",Inverse="#2C7FB8",Unavailable="grey80"))+
    scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=.32),guide="none")+
    scale_size_continuous(range=c(1.2,7),name=expression(-log[10](P)~"(tail compressed)"))+
    labs(title="a. Top 32 DANDELION targets across independent C2 evidence",
         subtitle="Ranked first by the number of independent supports, then DANDELION P; faint DANDELION points are sensitivity-only",
         x=NULL,y=NULL,fill=NULL)+theme_5c(8)
  tiers<-patt|>mutate(tier=factor(case_when(independent_support==0~"DANDELION only",independent_support==1~"+ 1 independent support",TRUE~"+ 2–3 independent supports"),
      levels=c("DANDELION only","+ 1 independent support","+ 2–3 independent supports")))|>
    count(tier,name="n_targets",.drop=FALSE)|>mutate(percent=n_targets/sum(n_targets),label=sprintf("%s (%.1f%%)",format(n_targets,big.mark=","),100*percent))
  pB<-ggplot(tiers,aes(n_targets,tier,fill=tier))+geom_col(width=.68)+geom_text(aes(label=label),hjust=-.08,fontface="bold")+
    scale_fill_manual(values=c("DANDELION only"="grey72","+ 1 independent support"="#6BAED6","+ 2–3 independent supports"="#2171B5"),guide="none")+
    scale_x_continuous(expand=expansion(mult=c(0,.30)))+
    labs(title="b. Independent support among all DANDELION targets",
      subtitle="Independent support = cis/local MR FDR, trans/distal MR FDR, or observational Bonferroni; DANDELION itself is not double-counted",
      x="Number of targets",y=NULL)+theme_5c(9)
  save_plot(pA|pB,"c2.Fig8.dandelion_mr_integration.png",17,10.5,outdir=outdir)
  wide
}

write_dandelion_package_network <- function(gene_pair,p_gene,rawdir,outdir){
  if(!nrow(gene_pair)||!all(c("region","gene2")%in%names(gene_pair)))return(NA_character_)
  netdir<-file.path(rawdir,"dandelion_network");dir.create(netdir,recursive=TRUE,showWarnings=FALSE)
  src<-tryCatch(quiet_package_call(
    DANDELION::gen_fig(gene.pair=as.data.frame(gene_pair),p.wes=p_gene,
      eta.wgs=.05/max(1,length(p_gene)),pic_dir=netdir)
  ),error=function(e){warning("DANDELION package network failed: ",conditionMessage(e),call.=FALSE);NULL})
  if(is.null(src)||!file.exists(src))return(NA_character_)
  # Preserve the package-native graph as a raw audit artifact; the sequential
  # publication figure is c2.Fig9.dandelion_network.png.
  dst<-file.path(netdir,"dandelion_package_network.pdf")
  if(!file.copy(src,dst,overwrite=TRUE))return(NA_character_)
  normalizePath(dst,winslash="/",mustWork=FALSE)
}

run_dandelion_step <- function(layer,assoc,ann,base,ygfile,rawdir,outdir,
  top_candidates=tibble(),mode=RUN_Dandelion){
  empty<-list(status="not_run",pairs=tibble(),gene_pairs=tibble(),targets=tibble(),lead_snps=tibble(),snp_gene_map=tibble(),
    evidence_plot=tibble(),exposure_qc=tibble(),qtl_audit=tibble(),input_audit=tibble(),integration=tibble(),
    sig_gene2=character(),non_sig_gene2=character(),gene_level_file=NA_character_,gene_evidence_type=NA_character_,network_file=NA_character_,
    primary_eligible=FALSE,analysis_class="not run",target_fraction=NA_real_,broad_selection_warning=FALSE,consolidation_eligible=FALSE)
  if(mode=="None")return(modifyList(empty,list(status="mode_none",analysis_class="not run by RUN_Dandelion=None")))
  if(layer!="protein")return(modifyList(empty,list(status="protein_only")))
  if(!requireNamespace("DANDELION",quietly=TRUE)){
    warning("DANDELION R package is not installed.",call.=FALSE)
    return(modifyList(empty,list(status="package_missing")))
  }
  gene_input<-find_dandelion_gene_file(ygfile);gene_file<-gene_input$path;gene_evidence_type<-gene_input$evidence_type
  primary_eligible<-isTRUE(gene_input$primary_eligible);analysis_class<-gene_input$analysis_class%||%"sensitivity"
  if(is.na(gene_file)){
    msg<-if(!DANDELION_ALLOW_MAGMA)
      "No explicit gene-level rare-variant/WES file was supplied; MAGMA sensitivity is disabled."
      else "No eligible gene-level disease-association file was found."
    warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,
      gene_evidence_type=gene_evidence_type,primary_eligible=primary_eligible,analysis_class=analysis_class)))
  }
  ref_assayed<-ann|>filter(!is.na(feature),is.finite(start),is.finite(end),!is.na(chr))|>
    transmute(gene_name=feature,type="protein_coding",Chromosome=paste0("chr",str_remove(chr,"^chr")),start=start,end=end)|>distinct(gene_name,.keep_all=TRUE)
  ref_all<-read_dandelion_gene_ref(ref_assayed)
  p_gene<-read_gene_level_p(gene_file,ref_assayed$gene_name);common<-intersect(ref_assayed$gene_name,names(p_gene))
  if(length(common)<20){
    msg<-paste0("Only ",length(common)," gene-level disease results matched assayed protein symbols. Check gene identifiers or set C2_MAGMA_GENE_MAP.")
    warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,gene_level_file=gene_file,gene_evidence_type=gene_evidence_type)))
  }
  lead0<-find_disease_leads(ygfile);lead_source<-attr(lead0,"source")%||%NA_character_
  lead<-lead0|>filter(SNP!="",is.finite(POS))|>distinct(SNP,.keep_all=TRUE)
  if(nrow(lead)<2){msg<-"Fewer than two independent disease-associated lead SNPs available.";warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,gene_level_file=gene_file,gene_evidence_type=gene_evidence_type,lead_snps=lead)))}
  gene2<-sort(common)
  if(mode=="Top")gene2<-intersect(as.character(top_candidates$feature%||%character()),gene2)
  if(DANDELION_MAX_GENE2>0&&length(gene2)>DANDELION_MAX_GENE2)
    gene2<-head(gene2[order(p_gene[gene2],na.last=TRUE)],DANDELION_MAX_GENE2)
  min_gene2<-if(mode=="Top")5L else 20L
  if(length(gene2)<min_gene2){
    msg<-paste0("RUN_Dandelion=",mode," retained only ",length(gene2),
      " gene2 proteins; need at least ",min_gene2,".")
    warning("DANDELION: ",msg,call.=FALSE)
    return(modifyList(empty,list(status=msg,gene_level_file=gene_file,gene_evidence_type=gene_evidence_type)))
  }
  message("C2/DANDELION [",mode,"]: ",length(lead$SNP)," independent disease SNPs x ",
    length(gene2)," assayed gene2 proteins; gene evidence=",basename(gene_file))
  ptrans<-build_trans_p_matrix(gene2,lead,base);qtl_audit<-as_tibble(attr(ptrans,"qtl_audit")%||%tibble())
  okrow<-rowSums(is.finite(ptrans))>=2;ptrans<-ptrans[okrow,,drop=FALSE];pwes<-p_gene[rownames(ptrans)]
  if(nrow(ptrans)<min_gene2){msg<-"Too few selected proteins had trans-pQTL data at the disease SNPs.";warning("DANDELION: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,gene_level_file=gene_file,gene_evidence_type=gene_evidence_type,lead_snps=lead,qtl_audit=qtl_audit)))}
  snpref<-lead|>transmute(SNP,SNPPos=POS,SNPChr=str_remove(as.character(CHR),"^chr"))
  # ref.table may contain genome-wide genes; med_gene() uses it to remove local
  # cis targets.  Only gene2 present in p.trans/p.wes can become DPGs.
  ref_for_med<-ref_all|>filter(type%in%c("lincRNA","protein_coding"),!Chromosome%in%c("chrM","chrX","chrY"))|>distinct(gene_name,.keep_all=TRUE)
  uniq_snp<-read_dandelion_snp_gene_map(lead)
  med<-tryCatch(quiet_package_call(
    DANDELION::med_gene(p.trans=ptrans,p.wes=pwes,ref.table=ref_for_med,gene1.list=colnames(ptrans),target.fdr=DANDELION_FDR,
      dist=DANDELION_CIS_BP,gene1.type="SNP",SNP.ref=snpref,n.cores=N_CORES,verbose=FALSE)
  ),error=function(e)e)
  if(inherits(med,"condition")){msg<-conditionMessage(med);warning("DANDELION failed: ",msg,call.=FALSE);return(modifyList(empty,list(status=msg,gene_level_file=gene_file,gene_evidence_type=gene_evidence_type,lead_snps=lead,snp_gene_map=uniq_snp,qtl_audit=qtl_audit)))}
  pairs<-tryCatch(quiet_package_call(
    DANDELION::calc_pair.snp(mat.sig=med$mat.sig,mat.p=med$mat.p,p.wes=pwes,gene1=med$gene1,uniq_snp=uniq_snp,
      ref.table.keep=ref_for_med,eta.wgs=.05/max(1,length(pwes)),SNP.ref=snpref,verbose=FALSE)
  ),error=function(e)e)
  if(inherits(pairs,"condition")){warning("DANDELION post-processing failed: ",conditionMessage(pairs),call.=FALSE);pairs<-NULL}
  pd<-if(is.null(pairs))tibble() else as_tibble(pairs$pairs_dact)
  gp<-if(is.null(pairs))tibble() else as_tibble(pairs$gene.pair)
  sig_gene2<-if(is.null(pairs))character()else as.character(pairs$sig_gene2%||%character())
  non_sig_gene2<-if(is.null(pairs))character()else as.character(pairs$non_sig.gene2%||%character())
  if(nrow(pd)){
    pd<-pd|>mutate(trans_p=map2_dbl(gene2,rsid,~if(.x%in%rownames(ptrans)&&.y%in%colnames(ptrans))ptrans[.x,.y]else NA_real_),
      gene_level_p=pwes[gene2],global_pair_BH=p.adjust(DANDELION_p,"BH"),gene_evidence_type=gene_evidence_type)
    tg<-pd|>group_by(gene2)|>summarise(n_selected_pairs=n(),n_distal_loci=n_distinct(rsid),DANDELION_p=safe_min_finite(DANDELION_p),
      gene_level_p=safe_min_finite(gene_level_p),best_trans_p=safe_min_finite(trans_p),.groups="drop")|>
      mutate(target_BH=p.adjust(DANDELION_p,"BH"),gene_evidence_type=gene_evidence_type,
        gene_level_bonferroni=gene_level_p<.05/max(1,length(pwes)),
        dpg_class=case_when(gene2%in%sig_gene2~"DANDELION + conventional gene evidence",
                            gene2%in%non_sig_gene2~"DANDELION-only prioritized target",
                            TRUE~"DANDELION prioritized target"))|>arrange(DANDELION_p)
  } else tg<-tibble()
  target_fraction<-nrow(tg)/max(1,nrow(ptrans))
  broad_selection_warning<-is.finite(target_fraction)&&target_fraction>DANDELION_MAX_TARGET_FRACTION
  consolidation_eligible<-primary_eligible&&!broad_selection_warning
  selection_warning<-if(broad_selection_warning)
    paste0("Selected ",scales::percent(target_fraction,accuracy=.1)," of tested gene2 targets; exceeds ",
      scales::percent(DANDELION_MAX_TARGET_FRACTION,accuracy=1)," audit threshold") else NA_character_
  if(nrow(pd))pd<-pd|>mutate(primary_eligible=primary_eligible,analysis_class=analysis_class,
    target_fraction=target_fraction,broad_selection_warning=broad_selection_warning,
    consolidation_eligible=consolidation_eligible,selection_warning=selection_warning)
  if(nrow(tg))tg<-tg|>mutate(primary_eligible=primary_eligible,analysis_class=analysis_class,
    target_fraction=target_fraction,broad_selection_warning=broad_selection_warning,
    consolidation_eligible=consolidation_eligible,selection_warning=selection_warning)
  # If no explicit SNP-to-gene file was supplied, expose the mapping inferred by
  # calc_pair.snp() so the distal-locus interpretation remains auditable.
  inferred_map<-if(nrow(pd)&&all(c("rsid","gene1")%in%names(pd)))pd|>transmute(SNP=rsid,GeneSymbol=gene1)|>filter(!is.na(GeneSymbol),GeneSymbol!="")|>distinct() else tibble()
  map_out<-bind_rows(as_tibble(uniq_snp),inferred_map)|>distinct(SNP,.keep_all=TRUE)
  matp<-med$mat.p[rownames(ptrans),colnames(ptrans),drop=FALSE];mats<-med$mat.sig[rownames(ptrans),colnames(ptrans),drop=FALSE]
  evidence<-tibble(gene2=rep(rownames(ptrans),times=ncol(ptrans)),rsid=rep(colnames(ptrans),each=nrow(ptrans)),
    trans_p=as.vector(ptrans),gene_level_p=rep(as.numeric(pwes),times=ncol(ptrans)),
    DANDELION_p=as.vector(matp),significant=replace_na(as.vector(mats)!=0,FALSE))|>
    filter(is.finite(trans_p),is.finite(gene_level_p),is.finite(DANDELION_p))
  sig_i<-which(evidence$significant);other_i<-which(!evidence$significant)
  if(length(other_i)>50000L)other_i<-other_i[unique(as.integer(round(seq(1,length(other_i),length.out=50000L))))]
  evidence_plot<-evidence[sort(unique(c(sig_i,other_i))),,drop=FALSE]
  exposure_qc<-evidence|>group_by(rsid)|>summarise(n_targets_tested=n(),n_selected=sum(significant),
    selected_fraction=n_selected/n_targets_tested,min_trans_p=safe_min_finite(trans_p),
    min_DANDELION_p=safe_min_finite(DANDELION_p),.groups="drop")|>
    mutate(broad_regulator=n_targets_tested>=min_gene2&selected_fraction>DANDELION_MAX_TARGET_FRACTION)|>
    left_join(lead|>select(rsid=SNP,CHR,POS,outcome_P=P),by="rsid")
  input_audit<-tibble(metric=c("RUN_Dandelion","outcome_GWAS","independent_disease_SNP_source","gene_level_disease_file","gene_evidence_type","analysis_class","primary_eligible","method_boundary",
    "protein_QTL_project","n_independent_disease_SNPs","n_gene2_with_gene_evidence","n_gene2_with_at_least_two_trans_tests","n_finite_pairs","n_selected_pairs","n_selected_DPGs","selected_target_fraction","broad_selection_warning","consolidation_eligible","target_FDR","cis_exclusion_bp"),
    value=as.character(c(mode,ygfile,lead_source,gene_file,gene_evidence_type,
      analysis_class,primary_eligible,
      "trans-pQTL adaptation; the Cell reference used trans-eQTL plus gene-level rare-variant/WES evidence",base,nrow(lead),length(common),nrow(ptrans),nrow(evidence),sum(evidence$significant),nrow(tg),target_fraction,broad_selection_warning,consolidation_eligible,DANDELION_FDR,DANDELION_CIS_BP)))
  write_raw_csv(lead,"c2.dandelion_lead_snps.csv",rawdir);write_raw_csv(map_out,"c2.dandelion_snp_to_cis_gene.csv",rawdir)
  write_raw_csv(pd,"c2.dandelion_pairs.csv",rawdir);write_raw_csv(gp,"c2.dandelion_gene_pairs.csv",rawdir);write_raw_csv(tg,"c2.dandelion_targets.csv",rawdir)
  write_raw_csv(qtl_audit,"c2.dandelion_qtl_coverage.csv",rawdir);write_raw_csv(exposure_qc,"c2.dandelion_exposure_qc.csv",rawdir)
  write_raw_csv(input_audit,"c2.dandelion_input_audit.csv",rawdir)
  data.table::fwrite(evidence,file.path(rawdir,"c2.dandelion_tested_pairs.csv.gz"),sep=",",na="NA",compress="gzip")
  plot_dandelion_results(tg,pd,gp,outdir);plot_dandelion_landscape(evidence_plot,tg,outdir,gene_evidence_type)
  network_file<-write_dandelion_package_network(gp,pwes,rawdir,outdir)
  list(status="ok",mode=mode,result=med,pairs=pd,gene_pairs=gp,targets=tg,lead_snps=lead,snp_gene_map=map_out,
       sig_gene2=sig_gene2,non_sig_gene2=non_sig_gene2,
       gene_level_file=gene_file,gene_evidence_type=gene_evidence_type,evidence_plot=evidence_plot,
       exposure_qc=exposure_qc,qtl_audit=qtl_audit,input_audit=input_audit,network_file=network_file,
       ptrans_dimensions=dim(ptrans),target_fdr=DANDELION_FDR,cis_bp=DANDELION_CIS_BP,
       primary_eligible=primary_eligible,analysis_class=analysis_class,target_fraction=target_fraction,
       broad_selection_warning=broad_selection_warning,consolidation_eligible=consolidation_eligible,
       gene_annotation=Sys.getenv("C2_DANDELION_GENE_ANNOTATION",unset=""),snp_gene_map_file=Sys.getenv("C2_DANDELION_SNP_GENE_MAP",unset=""))
}


# 🚩 Main C2
instrument_availability_audit <- function(iv_list,ygwas,layer){
  classes<-if(layer=="protein")c("cis","trans")else c("local","distal")
  imap_dfr(iv_list,function(obj,feature){
    iv<-obj$instruments%||%tibble();fs<-obj$files%||%list(joint=NA_character_,full=NA_character_)
    map_dfr(classes,function(cc){
      z<-if(nrow(iv))iv|>filter(analysis==cc)else tibble();ov<-if(nrow(z))sum(z$SNP%in%ygwas$SNP)else 0L
      hz<-if(nrow(z)&&ov>0)tryCatch(nrow(harmonize_sumstats(z,ygwas)),error=function(e)0L)else 0L
      reason<-case_when(nrow(z)>0&&hz>0~"Available",
        is.na(fs$joint)||!file.exists(fs$joint)~"No independent clump/COJO instrument file",
        !nrow(iv)~"No allele-complete independent instrument",
        !nrow(z)~paste0("No ",cc," instrument after genomic classification"),
        ov==0~"No SNP overlap with outcome GWAS",hz==0~"Alleles could not be harmonized",TRUE~"Unavailable")
      tibble(exposure=feature,analysis=cc,joint_file=fs$joint%||%NA_character_,full_file=fs$full%||%NA_character_,
        n_independent=nrow(iv),n_class=nrow(z),n_outcome_overlap=ov,n_harmonized=hz,status=reason)
    })
  })
}

restore_c2_figures <- function(cached,layer,outdir,rawdir=NULL){
  mr<-as_tibble(cached$MR%||%tibble());assoc<-as_tibble(cached$observational%||%tibble())
  if(!nrow(mr)||!nrow(assoc))return(invisible(FALSE))
  fig1<-plot_c2_fig1(mr,assoc,layer);save_c2_plot(fig1$plot,"c2.Fig1.prots.top.png",17,13.5,outdir=outdir)
  save_c2_plot(plot_qtl_variance(mr,layer),"c2.Fig2.pQTL_R2.png",15.5,10.5,outdir=outdir)
  fig3<-plot_c2_fig2(mr,assoc,layer);save_plot(fig3$plot,"c2.Fig3.effect_concordance.png",16,12.5,outdir=outdir)
  save_plot(plot_c2_fig4(mr,layer),"c2.Fig4.sensitivity_architecture.png",15.5,11.5,outdir=outdir)
  if(!is.null(rawdir))plot_mrlink2_results(read_mrlink2_results(rawdir),outdir)
  dan<-cached$DANDELION%||%list()
  if(layer=="protein"){
    plot_dandelion_results(dan$targets%||%tibble(),dan$pairs%||%tibble(),dan$gene_pairs%||%tibble(),outdir)
    plot_dandelion_landscape(dan$evidence_plot%||%tibble(),dan$targets%||%tibble(),outdir,dan$gene_evidence_type%||%"gene-level disease association")
    integration<-plot_dandelion_mr_integration(dan,mr,assoc,outdir)
    if(!is.null(rawdir)&&nrow(dan$gene_pairs%||%tibble())){
      p_gene<-tryCatch(read_gene_level_p(dan$gene_level_file%||%NA_character_,unique(as.character(dan$gene_pairs$gene2))),error=function(e)numeric())
      if(length(p_gene))write_dandelion_package_network(dan$gene_pairs, p_gene, rawdir, outdir)
    }
    if(nrow(integration)&&!is.null(rawdir))write_raw_csv(integration,"c2.dandelion_mr_integration.csv",rawdir)
  }
  dirint<-cached$directionality_integration%||%tibble()
  save_plot(plot_c2_directionality(dirint),"c2.Fig10.directionality_causal.png",16,12,outdir=outdir)
  if(nrow(dirint)&&!is.null(rawdir))write_raw_csv(dirint,"c2.directionality_causal.csv",rawdir)
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

run_reverse_mr_stage <- function(iv_list,ygfile,layer){
  # Disease-liability -> omic MR uses independent disease GWAS instruments as
  # exposure and each full pQTL/mQTL summary as outcome. It complements (and
  # does not replace) clinical lead-time sensitivity because liability is not
  # the same estimand as established disease or treatment.
  leads<-find_disease_leads(ygfile);if(!nrow(leads))return(list(MR=tibble(),audit=tibble()))
  yiv<-read_sumstat_snps(ygfile,leads$SNP,N_DEFAULT)|>filter(SNP%in%leads$SNP)
  if(!nrow(yiv))return(list(MR=tibble(),audit=tibble(feature=names(iv_list),status="No disease-IV beta/SE in outcome GWAS")))
  rows<-parallel_map(names(iv_list),function(feature){
    fs<-iv_list[[feature]]$files
    f<-fs$full%||%NA_character_
    if(length(f)!=1L||is.na(f)||!file.exists(f))
      return(list(result=tibble(),audit=tibble(feature=feature,status="Full QTL summary missing",n_disease_IV=nrow(yiv),n_QTL_overlap=0L)))
    qtl<-tryCatch(read_sumstat_snps(f,yiv$SNP,N_DEFAULT),error=function(e)tibble())
    if(!nrow(qtl))return(list(result=tibble(),audit=tibble(feature=feature,status="No disease-IV overlap in QTL",n_disease_IV=nrow(yiv),n_QTL_overlap=0L)))
    rr<-run_mr(yiv,qtl,Y,"disease_liability_to_omic")|>mutate(feature=feature,.before=1)
    list(result=rr,audit=tibble(feature=feature,status=ifelse(rr$n_IV>0,"Available","Harmonization failed"),n_disease_IV=nrow(yiv),n_QTL_overlap=nrow(qtl)))
  })
  mr<-bind_rows(lapply(rows,`[[`,"result"));if(nrow(mr))mr<-mr|>mutate(FDR_reverse=p.adjust(pval,"BH"))|>arrange(pval)
  list(MR=mr,audit=bind_rows(lapply(rows,`[[`,"audit")),disease_instruments=yiv)
}

plot_bidirectional_mr <- function(forward,reverse,layer){
  fw<-forward|>filter(is.finite(pval))|>group_by(exposure)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()|>
    transmute(feature=exposure,forward_beta=b,forward_p=pval,forward_FDR=FDR_all)
  rv<-if(nrow(reverse)&&all(c("feature","b","pval","FDR_reverse","n_IV")%in%names(reverse)))reverse|>
    filter(is.finite(pval))|>transmute(feature,reverse_beta=b,reverse_p=pval,reverse_FDR=FDR_reverse,n_reverse_IV=n_IV)else
    tibble(feature=character(),reverse_beta=numeric(),reverse_p=numeric(),reverse_FDR=numeric(),n_reverse_IV=integer())
  z<-full_join(fw,rv,by="feature")|>mutate(class=case_when(
    is.finite(forward_FDR)&forward_FDR<.05&is.finite(reverse_FDR)&reverse_FDR<.05~"Bidirectional genetic support",
    is.finite(forward_FDR)&forward_FDR<.05~"Omic -> disease only",
    is.finite(reverse_FDR)&reverse_FDR<.05~"Disease liability -> omic only",TRUE~"No FDR support"),
    x=sign(forward_beta)*-log10(pmax(forward_p,1e-300)),y=sign(reverse_beta)*-log10(pmax(reverse_p,1e-300)),
    label_p=pmin(forward_p,reverse_p,na.rm=TRUE),label_p=ifelse(is.finite(label_p),label_p,NA_real_))
  pa<-if(!nrow(z))blank_plot("Bidirectional MR","No reverse-direction estimate")else ggplot(z,aes(x,y,color=class))+
    geom_hline(yintercept=0,color="grey80")+geom_vline(xintercept=0,color="grey80")+geom_point(alpha=.72,size=2)+
    geom_text_repel(data=z|>filter(class!="No FDR support",is.finite(label_p))|>slice_min(label_p,n=20),aes(label=feature),size=2.8,max.overlaps=Inf)+
    labs(title="a. Forward and reverse-direction MR",subtitle="Axes are signed -log10(P); reverse MR estimates disease liability, not post-diagnosis treatment effects",
      x="Omic -> disease MR",y="Disease liability -> omic MR",color=NULL)+theme_5c(9)+theme(legend.position="top")
  cnt<-z|>count(class,name="biomarkers")
  pb<-if(!nrow(cnt))blank_plot("b. Direction classes")else ggplot(cnt,aes(biomarkers,fct_reorder(class,biomarkers),fill=class))+
    geom_col(width=.68,show.legend=FALSE)+geom_text(aes(label=biomarkers),hjust=-.12,fontface="bold")+
    scale_x_continuous(expand=expansion(mult=c(0,.20)))+labs(title="b. Genetic direction classes",x="Biomarkers",y=NULL)+theme_5c(9)
  pa|pb
}

build_genetic_score_manifest <- function(iv_list,layer){
  primary<-if(layer=="protein")"cis"else"local"
  imap_dfr(iv_list,function(obj,feature){
    z<-as_tibble(obj$instruments%||%tibble());if(!nrow(z))return(tibble())
    z|>filter(analysis==primary,is.finite(BETA),is.finite(SE),!is.na(EA),EA!="")|>
      transmute(feature,SNP,effect_allele=EA,other_allele=NEA,weight=BETA,SE,P,EAF,N,
        component=paste0("genetically predicted ",layer," (",primary," instruments)"),
        interpretation="Instrument-weight manifest; individual scores require UKB genotypes and cross-fitted/external pQTL weights")
  })
}

integrate_c2_directionality <- function(mr,c1,dan=list(),reverse_mr=tibble()){
  d<-as_tibble(c1$directionality%||%tibble());if(!nrow(d))return(tibble())
  mb<-mr|>filter(is.finite(pval))|>group_by(exposure)|>slice_min(pval,n=1,with_ties=FALSE)|>ungroup()|>
    transmute(term=exposure,MR_analysis=analysis,MR_beta=b,MR_p=pval,MR_FDR=FDR_all)
  tg<-as_tibble(dan$targets%||%tibble())
  if(nrow(tg)&&!"dpg_class"%in%names(tg))tg$dpg_class<-"DANDELION prioritized target"
  if(nrow(tg)&&!"consolidation_eligible"%in%names(tg)){
    primary_input<-if("gene_evidence_type"%in%names(tg))
      !str_detect(str_to_lower(coalesce(tg$gene_evidence_type,"")),"magma|adapted") else FALSE
    tested_n<-suppressWarnings(as.numeric((dan$ptrans_dimensions%||%c(NA_real_))[[1]]))
    frac<-if(is.finite(tested_n)&&tested_n>0)nrow(tg)/tested_n else NA_real_
    tg$consolidation_eligible<-primary_input&!(is.finite(frac)&&frac>DANDELION_MAX_TARGET_FRACTION)
  }
  dg<-if(nrow(tg))tg|>transmute(term=as.character(gene2),DANDELION_p,
      DANDELION_class=dpg_class%||%"DANDELION prioritized target",
      DANDELION_primary=as.logical(consolidation_eligible))else
    tibble(term=character(),DANDELION_p=numeric(),DANDELION_class=character(),DANDELION_primary=logical())
  rv<-if(nrow(reverse_mr))reverse_mr|>transmute(term=feature,reverse_MR_beta=b,reverse_MR_p=pval,reverse_MR_FDR=FDR_reverse)else
    tibble(term=character(),reverse_MR_beta=numeric(),reverse_MR_p=numeric(),reverse_MR_FDR=numeric())
  d|>left_join(mb,by="term")|>left_join(rv,by="term")|>left_join(dg,by="term")|>
    mutate(MR_support=is.finite(MR_FDR)&MR_FDR<.05,
           reverse_MR_support=is.finite(reverse_MR_FDR)&reverse_MR_FDR<.05,
           DANDELION_sensitivity=is.finite(DANDELION_p),
           DANDELION_support=DANDELION_sensitivity&replace_na(DANDELION_primary,FALSE),
           reactive_warning=coalesce(reactive_compatible,FALSE),
           distal_antecedent=coalesce(distal_support,FALSE),
           causal_interpretation=case_when(
             MR_support&reactive_warning~"Forward MR + reactive-compatible (mixed)",
             MR_support&distal_antecedent~"Forward MR + distal antecedent",
             reverse_MR_support&!MR_support~"Disease liability -> omic MR supported",
             MR_support~"Forward MR supported",
             !MR_support&DANDELION_support&reactive_warning~"DANDELION + reactive-compatible",
             !MR_support&DANDELION_support~"DANDELION-prioritized",
             !MR_support&DANDELION_sensitivity~"DANDELION sensitivity only",
             reactive_warning~"Observational reactive-compatible",
             distal_antecedent~"Observational distal antecedent",
             TRUE~"no causal support"))|>arrange(desc(MR_support),desc(DANDELION_support),p_incident)
}

plot_c2_directionality <- function(z,anchors=unique(trimws(strsplit(Sys.getenv("C1_DIRECTION_ANCHORS",
  unset="PCSK9,LPA,GDF15,NTPROBNP,MMP12"),",",fixed=TRUE)[[1]]))){
  if(!nrow(z))return(blank_plot("Causal evidence with temporal directionality","C1 directionality table was unavailable"))
  if(!"DANDELION_sensitivity"%in%names(z))z$DANDELION_sensitivity<-z$DANDELION_support
  chosen<-unique(c(anchors,z|>filter(MR_support|DANDELION_sensitivity)|>
    arrange(pmin(MR_p,DANDELION_p,na.rm=TRUE))|>slice_head(n=24)|>pull(term)))
  long<-z|>filter(term%in%chosen)|>transmute(term,`Incident Cox`=incident_score,`5-y landmark`=landmark5_score,
      `Prevalent logistic`=prevalent_score,`Duration slope`=duration_score,
      `Forward MR`=sign(MR_beta)*pmin(12,-log10(pmax(MR_p,1e-300))),
      `Reverse MR`=sign(reverse_MR_beta)*pmin(12,-log10(pmax(reverse_MR_p,1e-300))),
      `DANDELION`=pmin(12,-log10(pmax(DANDELION_p,1e-300))))|>
    pivot_longer(-term,names_to="evidence",values_to="signed_score")
  pa<-ggplot(long,aes(evidence,factor(term,levels=rev(chosen)),fill=signed_score))+geom_tile(color="white")+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-12,12),name="signed\n-log10(P)")+
    labs(title="a. Genetic evidence beside Yin-Yang temporal evidence",x=NULL,y=NULL)+theme_5c(8)+theme(axis.text.x=element_text(angle=35,hjust=1))
  sc<-z|>filter(is.finite(MR_beta),is.finite(beta_incident))
  pb<-if(!nrow(sc))blank_plot("b. MR and incident association","No overlapping estimate")else
    ggplot(sc,aes(MR_beta,beta_incident,color=causal_interpretation))+geom_hline(yintercept=0,color="grey85")+geom_vline(xintercept=0,color="grey85")+
      geom_point(size=2,alpha=.75)+geom_text_repel(data=sc|>filter(term%in%anchors|reactive_warning),aes(label=term),size=3,max.overlaps=20,fontface="bold")+
      labs(title="b. Genetic causality and reactive signal may coexist",subtitle="Prevalent/duration evidence is a disease-state warning, not a veto on valid cis-MR",x="MR beta",y="Incident Cox beta",color=NULL)+theme_5c(9)
  pc<-z|>count(causal_interpretation,name="proteins")|>ggplot(aes(proteins,reorder(causal_interpretation,proteins),fill=causal_interpretation))+
    geom_col(width=.72)+guides(fill="none")+labs(title="c. Integrated interpretation",x="Proteins",y=NULL)+theme_5c(9)
  pa/(pb|pc)+plot_layout(heights=c(1.25,1))
}

grade_c2_evidence <- function(mr,layer,mrlink2=list(results=tibble())){
  primary<-if(layer=="protein")"cis"else"local"
  z<-mr|>filter(analysis==primary,is.finite(pval))|>mutate(
    heterogeneous=n_IV>1&is.finite(Q_p)&Q_p<.05,
    instrument_architecture=case_when(n_IV<=1~"single IV",n_IV<=5~"oligogenic (2-5 IVs)",
      n_IV<=10~"multi-IV (6-10)",TRUE~"many-IV (>10)"),
    many_IV_warning=n_IV>10,
    egger_warning=n_IV>=3&is.finite(egger_intercept_p)&egger_intercept_p<.05,
    weighted_median_concordant=n_IV==1|(is.finite(tsmr_weighted_median_b)&
      is.finite(tsmr_weighted_median_p)&tsmr_weighted_median_p<.05&sign(tsmr_weighted_median_b)==sign(b)),
    steiger_ok=is.finite(steiger_support_fraction)&steiger_support_fraction>.5,
    evidence_grade=case_when(
      FDR_all<.05&steiger_ok&!egger_warning&!heterogeneous~"A: robust",
      FDR_all<.05&steiger_ok&!egger_warning&heterogeneous&weighted_median_concordant~"B: heterogeneous",
      FDR_all<.05~"C: sensitivity warning",
      pval<.05~"D: nominal only",TRUE~"E: unsupported"),
    grade_reason=case_when(
      evidence_grade=="A: robust"~"FDR + Steiger; no Egger warning or Cochran-Q heterogeneity",
      evidence_grade=="B: heterogeneous"~"FDR + Steiger; heterogeneous, but weighted-median direction is concordant",
      evidence_grade=="C: sensitivity warning"~"FDR supported, but Steiger or Egger sensitivity warning",
      evidence_grade=="D: nominal only"~"Nominal primary MR only",TRUE~"No nominal primary MR support"))
  rr<-as_tibble(mrlink2$results%||%tibble())
  if(nrow(rr)){
    tc<-pick_col_local(names(rr),c("^TRAIT$","^EXPOSURE$"))
    pc<-pick_col_local(names(rr),c("^P\\(ALPHA\\)$","^P_ALPHA$"))
    if(!is.na(tc)&&!is.na(pc)){
      ml<-tibble(exposure=as.character(rr[[tc]]),MRlink2_p=suppressWarnings(as.numeric(rr[[pc]])))|>
        filter(!is.na(exposure),is.finite(MRlink2_p))|>group_by(exposure)|>
        slice_min(MRlink2_p,n=1,with_ties=FALSE)|>ungroup()|>
        mutate(MRlink2_FDR=p.adjust(MRlink2_p,"BH"))
      z<-z|>left_join(ml,by="exposure")
    }
  }
  if(!"MRlink2_p"%in%names(z))z<-z|>mutate(MRlink2_p=NA_real_,MRlink2_FDR=NA_real_)
  z|>arrange(factor(evidence_grade,levels=c("A: robust","B: heterogeneous",
    "C: sensitivity warning","D: nominal only","E: unsupported")),pval)
}

plot_c2_evidence_grades <- function(z,layer){
  if(!nrow(z))return(blank_plot("C2 evidence grades","No cis/local MR estimate"))
  chosen<-unique(c(C2_FIXED_TOP,z|>arrange(pval)|>slice_head(n=24)|>pull(exposure)))
  d<-z|>filter(exposure%in%chosen)|>mutate(
    exposure=factor(exposure,levels=rev(unique(exposure))),lo=b-1.96*se,hi=b+1.96*se)
  pa<-ggplot(d,aes(b,exposure,color=evidence_grade))+geom_vline(xintercept=0,color="grey75")+
    geom_errorbarh(aes(xmin=lo,xmax=hi),height=.10)+geom_point(aes(shape=heterogeneous),size=2.5)+
    labs(title=paste0("a. ",if(layer=="protein")"cis"else"local"," MR evidence grade"),
      subtitle="Significance is separated from heterogeneity, Egger and Steiger diagnostics",
      x="MR log-odds effect",y=NULL,color=NULL,shape="Cochran Q P<0.05")+
    theme_5c(9)+theme(legend.position="bottom")
  cnt<-z|>count(evidence_grade,name="biomarkers")|>mutate(evidence_grade=factor(evidence_grade,
    levels=c("A: robust","B: heterogeneous","C: sensitivity warning","D: nominal only","E: unsupported")))
  pb<-ggplot(cnt,aes(biomarkers,evidence_grade,fill=evidence_grade))+geom_col(width=.70,show.legend=FALSE)+
    geom_text(aes(label=biomarkers),hjust=-.12,fontface="bold")+
    scale_x_continuous(expand=expansion(mult=c(0,.22)))+
    labs(title="b. Evidence-grade distribution",
      subtitle="C2 grades remain provisional until C3 colocalization/fine-mapping",
      x="Biomarkers",y=NULL)+theme_5c(9)
  pa|pb
}

run_c2_layer <- function(layer=c("protein","metabolite")){
  layer<-match.arg(layer);outdir<-if(layer=="protein")out.prot else out.met;setwd2(outdir)
  rawdir<-le8_job_dir(outdir,LE8_JOB);dir.create(rawdir,recursive=TRUE,showWarnings=FALSE);cache<-file.path(rawdir,"c2.res.rds")
  method_scope<-c2_method_scope_audit(layer);write_raw_csv(method_scope,"c2.method_scope_audit.csv",rawdir)
  h2_file0<-find_c2_heritability_file(layer)
  h2_signature<-if(is.na(h2_file0))"h2=missing"else paste0("h2=",h2_file0,
    ";size=",file.size(h2_file0),";mtime=",format(file.info(h2_file0)$mtime,tz="UTC",usetz=TRUE))
  component_covariate_signature<-paste0("LE4=",Sys.getenv("C2_LE4_COVARS",
    unset="diet.pts,pa.pts,smoke.pts,sleep.pts"),";treatment=",
    Sys.getenv("C2_TREATMENT_VARS",unset=Sys.getenv("C1_TREATMENT_VARS",unset="")))
  mode_signature<-paste0("MRlink2=",RUN_MRlink2,";Dandelion=",RUN_Dandelion,
    ";",h2_signature,";",component_covariate_signature)
  pgs_file0<-find_c2_score_file(layer)
  pgs_check<-length(pgs_file0)==1L&&!is.na(pgs_file0)&&nzchar(pgs_file0)&&
    file.exists(pgs_file0)&&file.size(pgs_file0)>0
  pgs_expected<-file.path(indir,"Rdata",if(layer=="protein")"prot.pgs.rds"else"met.pgs.rds")
  pgs_signature<-if(is.na(pgs_file0))"PGS=missing"else paste0("PGS=",pgs_file0,";size=",file.size(pgs_file0),
    ";mtime=",format(file.info(pgs_file0)$mtime,tz="UTC",usetz=TRUE))
  if(!pgs_check)message("C2/",layer,": PGS check: missing ",pgs_expected,
    "; skip individual genetic decomposition and lead-time analysis")
  if(cache_valid(cache)){
    old<-tryCatch(readRDS(cache),error=function(e)NULL)
    if(!is.null(old)&&identical(old$meta$code_version%||%NA_character_,C2_CODE_VERSION)&&
       identical(old$meta$mode_signature%||%NA_character_,mode_signature)&&
       identical(old$meta$pgs_signature%||%NA_character_,pgs_signature)){
      cache_message(paste0("C2/",layer),cache);return(old)
    }
    message("C2/",layer,": cache code/mode/PGS contract differs from ",C2_CODE_VERSION," / ",mode_signature,
      "; reusing compatible stage caches and rebuilding final outputs")
  }
  base0<-if(layer=="protein")dir.X else dir.met.gwas
  c1<-read_c1_cache_for_c2(outdir);assoc<-(c1$association%||%c1$pwas_incident%||%c1$MWAS)|>as_tibble();if(!"beta"%in%names(assoc))assoc<-assoc|>mutate(beta=safe_log(estimate))
  base<-base0;ann<-layer_annotation(layer,unique(assoc$term));ranked<-assoc|>filter(is.finite(p.value))|>arrange(p.value)|>pull(term)
  # MR requires an independent COJO set plus the corresponding full summary
  # statistics for allele recovery. A lone cis extract or LD (.ldr.cojo) file
  # must not consume a candidate slot.
  mapped<-ranked[vapply(ranked,function(x){fs<-find_qtl_files(x,base,layer);!is.na(fs$joint)&&!is.na(fs$full)},logical(1))]
  fixed_assayed<-intersect(C2_FIXED_TOP,unique(assoc$term))
  fixed_mapped<-fixed_assayed[vapply(fixed_assayed,function(x){fs<-find_qtl_files(x,base,layer);!is.na(fs$joint)&&!is.na(fs$full)},logical(1))]
  mr_candidates<-unique(c(fixed_mapped,head(mapped,MAX_FEATURES)))
  audit_features<-unique(c(fixed_assayed,head(ranked,MAX_FEATURES),mr_candidates))
  if(!length(mr_candidates))stop("C2/",layer,": no QTL files mapped under ",base,call.=FALSE)
  input_table<-tibble(feature=audit_features,top_observational=feature%in%head(ranked,MAX_FEATURES),MR_candidate=feature%in%mr_candidates)
  write_raw_csv(input_table,"c2.input_features.csv",rawdir);message("C2/",layer,": auditing ",length(audit_features)," traits; ",length(mr_candidates)," have mapped QTL files")
  iv_cache<-file.path(rawdir,paste0("c2.instruments.",C2_STAGE_VERSION,".rds"));iv_list<-read_stage_cache(iv_cache)
  if(is.null(iv_list)){
    message("C2/",layer,": build independent instrument stage")
    iv_list<-setNames(map(audit_features,~read_qtl_instruments(.x,base,layer,ann)),audit_features)
    write_stage_cache(iv_list,iv_cache)
  }else message("C2/",layer,": reuse independent instrument stage")
  ygfile<-get_y_gwas_file(Y,TRUE)
  mr_cache<-file.path(rawdir,paste0("c2.mr_stage.",C2_STAGE_VERSION,".rds"));mr_stage<-read_stage_cache(mr_cache)
  if(is.null(mr_stage)){
    all_snps<-unique(unlist(map(iv_list,~.$instruments$SNP)))
    ygwas<-read_sumstat_snps(ygfile,all_snps,N_DEFAULT);if(!nrow(ygwas))stop("No outcome-GWAS variants overlapped QTL instruments.",call.=FALSE)
    availability<-instrument_availability_audit(iv_list,ygwas,layer)
    mr<-imap_dfr(iv_list,function(obj,feature){iv<-obj$instruments;if(!nrow(iv))return(tibble());map_dfr(unique(iv$analysis),~run_mr(iv|>filter(analysis==.x),ygwas,feature,.x))})
    if(!nrow(mr))stop("C2/",layer,": no harmonized MR estimate.",call.=FALSE)
    mr<-mr|>group_by(analysis)|>mutate(FDR_analysis=p.adjust(pval,"BH"))|>ungroup()|>mutate(FDR_all=p.adjust(pval,"BH"))
    mr_stage<-list(MR=mr,availability=availability);write_stage_cache(mr_stage,mr_cache)
  }else{mr<-mr_stage$MR;availability<-mr_stage$availability;message("C2/",layer,": reuse harmonization/MR stage")}
  write_raw_csv(availability,"c2.instrument_availability.csv",rawdir);write_raw_csv(mr,"c2.MR_all.csv",rawdir)
  reverse_cache<-file.path(rawdir,paste0("c2.reverse_mr_stage.",C2_STAGE_VERSION,".rds"));reverse_stage<-read_stage_cache(reverse_cache)
  if(is.null(reverse_stage)){
    message("C2/",layer,": disease-liability -> omic reverse MR")
    reverse_stage<-run_reverse_mr_stage(iv_list,ygfile,layer);write_stage_cache(reverse_stage,reverse_cache)
  }else message("C2/",layer,": reuse reverse-direction MR stage")
  reverse_mr<-reverse_stage$MR%||%tibble();reverse_audit<-reverse_stage$audit%||%tibble()
  write_raw_csv(reverse_mr,"c2.reverse_MR_all.csv",rawdir);write_raw_csv(reverse_audit,"c2.reverse_MR_audit.csv",rawdir)
  top_candidates<-select_c2_top_candidates(layer,assoc,mr,names(iv_list))
  write_raw_csv(top_candidates,"c2.top_candidates.csv",rawdir)
  genetic_score_manifest<-build_genetic_score_manifest(iv_list,layer)
  write_raw_tsv(genetic_score_manifest,"c2.genetic_score_weights.tsv",rawdir)
  individual_decomposition<-if(pgs_check)
    run_individual_genetic_decomposition(layer,outdir,top_candidates) else
    list(status=tibble(status="not run",detail=paste0("PGS check failed; score file not found: ",pgs_expected)),
      summary=tibble(),trajectory=tibble())
  write_raw_csv(individual_decomposition$status%||%tibble(),"c2.individual_genetic_decomposition_status.csv",rawdir)
  write_raw_csv(individual_decomposition$heritability_status%||%tibble(),"c2.heritability_status.csv",rawdir)
  save_plot(plot_individual_genetic_decomposition(individual_decomposition$summary%||%tibble()),
    "c2.Fig12.genetic_decomposition.png",17,13,outdir=outdir)
  save_plot(plot_component_leadtime(individual_decomposition$trajectory%||%tibble()),
    "c2.Fig13.genetic_component_leadtime.png",17,11,outdir=outdir)
  save_plot(plot_bidirectional_mr(mr,reverse_mr,layer),"c2.Fig11.bidirectional_mr.png",16,8.5,outdir=outdir)
  fig1<-plot_c2_fig1(mr,assoc,layer);save_c2_plot(fig1$plot,"c2.Fig1.prots.top.png",17,13.5,outdir=outdir)
  save_c2_plot(plot_qtl_variance(mr,layer),"c2.Fig2.pQTL_R2.png",15.5,10.5,outdir=outdir)
  fig3<-plot_c2_fig2(mr,assoc,layer);save_plot(fig3$plot,"c2.Fig3.effect_concordance.png",16,12.5,outdir=outdir)
  save_plot(plot_c2_fig4(mr,layer),"c2.Fig4.sensitivity_architecture.png",15.5,11.5,outdir=outdir)
  best<-fig1$best|>mutate(y=-log10(pmax(pval,1e-300)),direction=case_when(FDR_all<.05&b>0~"Positive",FDR_all<.05&b<0~"Inverse",TRUE~"NS"),label=ifelse(min_rank(pval)<=25,exposure,NA_character_))
  top_r2<-mr|>filter(is.finite(r2_median),n_IV>0)|>
    mutate(r2_median_pct=100*r2_median,r2_q25_pct=100*r2_q25,r2_q75_pct=100*r2_q75,
           r2_p90_pct=100*r2_p90,r2_max_pct=100*r2_max)|>arrange(desc(r2_p90_pct))
  arch<-mr|>filter(is.finite(median_F),is.finite(r2_median))|>
    mutate(evidence=case_when(FDR_all<.05~"MR FDR < 0.05",Q_p<.05~"Heterogeneous",TRUE~"Other"),
           label=ifelse(min_rank(pval)<=15,exposure,NA_character_))

  # DANDELION is complementary to MR: it asks whether a distal disease SNP has
  # trans-regulatory evidence to a gene/protein that itself has disease-genetic evidence.
  dan_cache<-file.path(rawdir,paste0("c2.dandelion_stage.",str_to_lower(RUN_Dandelion),".",C2_STAGE_VERSION,".rds"));dandelion<-read_stage_cache(dan_cache)
  if(is.null(dandelion)){
    dandelion<-run_dandelion_step(layer,assoc,ann,base,ygfile,rawdir,outdir,
      top_candidates,RUN_Dandelion)
    write_stage_cache(dandelion,dan_cache)
  }else message("C2/",layer,": reuse DANDELION stage")

  if(layer=="protein"){
    plot_dandelion_results(dandelion$targets%||%tibble(),dandelion$pairs%||%tibble(),dandelion$gene_pairs%||%tibble(),outdir)
    plot_dandelion_landscape(dandelion$evidence_plot%||%tibble(),dandelion$targets%||%tibble(),outdir,dandelion$gene_evidence_type%||%"gene-level disease association")
    dandelion$integration<-plot_dandelion_mr_integration(dandelion,mr,assoc,outdir)
    write_raw_csv(dandelion$integration%||%tibble(),"c2.dandelion_mr_integration.csv",rawdir)
    network_file<-dandelion$network_file%||%NA_character_
    network_missing<-length(network_file)!=1L||is.na(network_file)||!nzchar(network_file)||!file.exists(network_file)
    if(nrow(dandelion$gene_pairs%||%tibble())&&network_missing){
      p_gene<-tryCatch(read_gene_level_p(dandelion$gene_level_file%||%NA_character_,unique(as.character(dandelion$gene_pairs$gene2))),error=function(e)numeric())
      if(length(p_gene))dandelion$network_file<-write_dandelion_package_network(dandelion$gene_pairs,p_gene,rawdir,outdir)
    }
    if(!identical(dandelion$status%||%"","ok")){
      why<-paste0("DANDELION unavailable: ",dandelion$status%||%"no usable result")
      for(i in 6:9)save_plot(blank_plot(paste0("C2 Figure ",i),why),
        c("c2.Fig6.dandelion.png","c2.Fig7.dandelion_evidence.png",
          "c2.Fig8.dandelion_mr_integration.png","c2.Fig9.dandelion_network.png")[[i-5]],10,6,outdir=outdir)
    }else if(!nrow(dandelion$gene_pairs%||%tibble())){
      save_plot(blank_plot("DANDELION network","No significant gene-pair network was available"),
        "c2.Fig9.dandelion_network.png",10,6,outdir=outdir)
    }
  }else{
    why<-"DANDELION is a gene/protein regulatory-network analysis and is not defined for metabolites"
    for(i in 6:9)save_plot(blank_plot(paste0("C2 Figure ",i),why),
      c("c2.Fig6.dandelion.png","c2.Fig7.dandelion_evidence.png",
        "c2.Fig8.dandelion_mr_integration.png","c2.Fig9.dandelion_network.png")[[i-5]],10,6,outdir=outdir)
  }

  directionality_integration<-integrate_c2_directionality(mr,c1,dandelion,reverse_mr)
  write_raw_csv(directionality_integration,"c2.directionality_causal.csv",rawdir)
  save_plot(plot_c2_directionality(directionality_integration),"c2.Fig10.directionality_causal.png",16,12,outdir=outdir)

  jobs_all<-build_mrlink2_jobs(iv_list,ann,layer,ygfile)
  jobs_audit<-if(nrow(jobs_all))jobs_all|>
    left_join(top_candidates|>select(trait=feature,top_rank,top_reason),by="trait")|>
    mutate(mode=RUN_MRlink2,selected=case_when(RUN_MRlink2=="All"~TRUE,
      RUN_MRlink2=="Top"~is.finite(top_rank),TRUE~FALSE),
      selection_reason=case_when(selected&RUN_MRlink2=="All"~"RUN_MRlink2=All",
        selected~top_reason,RUN_MRlink2=="None"~"RUN_MRlink2=None",TRUE~"outside Top set")) else
    tibble(trait=character(),mode=character(),selected=logical(),selection_reason=character())
  write_raw_csv(jobs_audit,"c2.mrlink2_job_audit.csv",rawdir)
  jobs<-if(nrow(jobs_all))jobs_audit|>filter(selected)|>select(all_of(names(jobs_all)))else jobs_all
  run_mrlink2_step(layer,rawdir,jobs,ygfile,RUN_MRlink2)
  mrlink2<-read_mrlink2_results(rawdir,RUN_MRlink2);plot_mrlink2_results(mrlink2,outdir)
  # Older completed runs sometimes retained the combined estimates but not the
  # generated job manifest. Recover a minimal provenance audit instead of
  # publishing a non-empty result figure beside an empty audit.
  if(!nrow(jobs_audit)&&nrow(mrlink2$results%||%tibble())){
    trait_col<-pick_col_local(names(mrlink2$results),c("^TRAIT$","^EXPOSURE$"))
    region_col<-pick_col_local(names(mrlink2$results),c("^REGION$"))
    trait<-if(is.na(trait_col))rep(NA_character_,nrow(mrlink2$results))else as.character(mrlink2$results[[trait_col]])
    region<-if(is.na(region_col))rep(NA_character_,nrow(mrlink2$results))else as.character(mrlink2$results[[region_col]])
    jobs_audit<-tibble(trait=trait,region=region,mode=RUN_MRlink2,selected=TRUE,
      selection_reason="Recovered from completed MR-link-2 output; original job manifest absent")|>
      filter(!is.na(trait),nzchar(trait))|>distinct()
    write_raw_csv(jobs_audit,"c2.mrlink2_job_audit.csv",rawdir)
  }
  evidence_grades<-grade_c2_evidence(mr,layer,mrlink2)
  write_raw_csv(evidence_grades,"c2.evidence_grades.csv",rawdir)
  save_plot(plot_c2_evidence_grades(evidence_grades,layer),
    "c2.Fig14.evidence_grades.png",17,9,outdir=outdir)

  out<-list(meta=module_meta(layer,extra=list(code_version=C2_CODE_VERSION,mode_signature=mode_signature,
            pgs_signature=pgs_signature,h2_signature=h2_signature,
            RUN_MRlink2=RUN_MRlink2,RUN_Dandelion=RUN_Dandelion)),MR=mr,MR_reverse=reverse_mr,MR_reverse_audit=reverse_audit,genetic_score_manifest=genetic_score_manifest,individual_decomposition=individual_decomposition,MR_best=best,observational=assoc,instrument_availability=availability,
            Fig1_summary=fig1$bars,Fig1_forest=fig1$forest,Fig3=fig3$wide,R2_QTL=top_r2,
            architecture=arch,top_candidates=top_candidates,evidence_grades=evidence_grades,
            DANDELION=dandelion,directionality_integration=directionality_integration,
            MRLink2_job_audit=jobs_audit,MRLink2_jobs=jobs,MRLink2=mrlink2,method_scope=method_scope)
  saveRDS(out,cache,compress="xz")
  write_xlsx2(list(MR_all=mr,MR_reverse=reverse_mr,MR_reverse_audit=reverse_audit,
                   method_scope=method_scope,evidence_grades=evidence_grades,top_candidates=top_candidates,
                   genetic_score_manifest=genetic_score_manifest,
                   genetic_decomp_status=individual_decomposition$status%||%tibble(),
                   heritability_status=individual_decomposition$heritability_status%||%tibble(),
                   genetic_decomp_summary=individual_decomposition$summary%||%tibble(),
                   genetic_leadtime=individual_decomposition$trajectory%||%tibble(),
                   MR_best=best,instrument_availability=availability,evidence_overlap=fig1$bars,effect_forest=fig1$forest,cis_local_vs_trans_distal=fig3$wide,QTL_R2=top_r2,instrument_architecture=arch,
                   DANDELION_input_audit=dandelion$input_audit%||%tibble(),DANDELION_QTL_coverage=dandelion$qtl_audit%||%tibble(),DANDELION_exposure_QC=dandelion$exposure_qc%||%tibble(),
                   DANDELION_lead_snps=dandelion$lead_snps%||%tibble(),DANDELION_snp_gene_map=dandelion$snp_gene_map%||%tibble(),DANDELION_pairs=dandelion$pairs%||%tibble(),DANDELION_gene_pairs=dandelion$gene_pairs%||%tibble(),DANDELION_targets=dandelion$targets%||%tibble(),DANDELION_MR_integration=dandelion$integration%||%tibble(),
                   directionality_causal=directionality_integration,
                   MRLink2_job_audit=jobs_audit,MRLink2_jobs=jobs,
                   MRLink2_results=mrlink2$results%||%tibble(),
                   MRLink2_status=mrlink2$status%||%tibble()),"c2.out.xlsx")
  finalize_outputs(LE8_JOB,outdir);out
}

if(prot_DO)invisible(run_c2_layer("protein"))
if(met_DO)invisible(run_c2_layer("metabolite"))
