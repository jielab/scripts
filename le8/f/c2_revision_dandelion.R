# DANDELION: disease locus -> trans-regulated target gene -> gene-level disease evidence.
# Official DACT output, conservative MaxP sensitivity, and other 5C evidence
# remain separate. No MR-derived gene P is accepted as independent burden evidence.
le8_safe_call <- function(fun,args) {
  f<-names(formals(fun));if(!"..."%in%f)args<-args[names(args)%in%f]
  do.call(fun,args)
}
find_dandelion_gene_file <- function(ygfile) {
  f<-Sys.getenv("C2_DANDELION_GENE_P_FILE",unset="")
  type<-Sys.getenv("C2_DANDELION_GENE_EVIDENCE",unset="unverified")
  allowed<-c("WES_disease","MAGMA_disease","other_disease","unverified")
  if(!type%in%allowed)stop("C2_DANDELION_GENE_EVIDENCE must be one of: ",paste(allowed,collapse=", "))
  if(nzchar(f)) {
    if(!file.exists(f))return(list(path=NA_character_,evidence_type="configured gene file missing",primary_eligible=FALSE,analysis_class="unavailable"))
    le8_dependency(f)
    matches<-identical(Sys.getenv("C2_DANDELION_GENE_OUTCOME",unset=""),Y)
    unfiltered<-identical(Sys.getenv("C2_DANDELION_GENE_SCOPE",unset="unknown"),"unfiltered")
    primary<-type=="WES_disease"&&matches&&unfiltered
    return(list(path=f,evidence_type=type,primary_eligible=primary,
      analysis_class=if(primary)"WES-supported trans-pQTL adaptation"else"declared/unknown-source sensitivity"))
  }
  f<-if(DANDELION_ALLOW_MAGMA)find_magma_gene_file(ygfile)else NA_character_
  le8_dependency(f)
  list(path=f,evidence_type="MAGMA_disease",primary_eligible=FALSE,
    analysis_class="MAGMA trans-pQTL sensitivity; common-variant evidence is not independent WES evidence")
}
read_gene_level_p <- function(file,target_symbols) {
  if(is.na(file)||!file.exists(file))return(numeric())
  le8_dependency(file);d<-data.table::fread(file,showProgress=FALSE,check.names=FALSE)
  gc<-pick_col_local(names(d),c("^GENE$","^GENE_ID$","^SYMBOL$","^GENE_SYMBOL$","^GENEID$"))
  pc<-Sys.getenv("C2_DANDELION_GENE_P_COLUMN",unset="")
  if(!nzchar(pc))pc<-pick_col_local(names(d),c("^P$","^PVAL$","^P_VALUE$","^BURDEN_P$","^P_BURDEN$","^PVAL_BURDEN$"))
  if(is.na(gc)||is.na(pc)||!pc%in%names(d))stop("Gene-level disease input needs gene and explicit P columns")
  # Optional prespecified test/mask choice is applied before multiplicity handling.
  test<-Sys.getenv("C2_DANDELION_GENE_TEST",unset="")
  if(nzchar(test)) {
    tc<-pick_col_local(names(d),c("^TEST$","^MASK$","^ANNOTATION$"))
    if(is.na(tc))stop("GENE_TEST configured but no TEST/MASK column")
    d<-d[as.character(d[[tc]])==test,,drop=FALSE]
  }
  z<-tibble(gene=map_magma_gene_ids(d[[gc]],target_symbols),p=as.numeric(d[[pc]]))|>
    filter(!is.na(gene),nzchar(gene),is.finite(p),p>=0,p<=1)|>
    mutate(p=pmax(p,.Machine$double.xmin))|>group_by(gene)|>
    summarise(n_tests=n(),p_raw_min=min(p),p=pmin(1,min(p)*n()),.groups="drop")
  # Bonferroni across masks/tests, not an unadjusted minimum.
  if(!is.null(.le8_revision_state$rawdir))write_raw_csv(z,"c2.dandelion_gene_p_aggregation.csv",.le8_revision_state$rawdir)
  ans<-setNames(z$p,z$gene)
  attr(ans,"n_gene_tests")<-length(unique(as.character(d[[gc]])[!is.na(d[[gc]])]))
  ans
}
le8_dandelion_family <- function(evidence,universe,alpha=.1) {
  d<-evidence
  valid<-is.finite(d$DANDELION_p)&d$DANDELION_p>=0&d$DANDELION_p<=1
  d$global_pair_BH<-d$global_pair_BY<-NA_real_
  d$global_pair_BH[valid]<-p.adjust(d$DANDELION_p[valid],"BH")
  d$global_pair_BY[valid]<-p.adjust(d$DANDELION_p[valid],"BY")
  d$maxP<-pmax(d$trans_p,d$gene_level_p)
  d$maxP_pair_BH<-NA_real_;ii<-is.finite(d$maxP);d$maxP_pair_BH[ii]<-p.adjust(d$maxP[ii],"BH")
  rule<-toupper(Sys.getenv("C2_DANDELION_MULTIPLICITY",unset="BY"))
  if(!rule%in%c("BH","BY"))stop("C2_DANDELION_MULTIPLICITY must be BH or BY")
  qp<-if(rule=="BY")d$global_pair_BY else d$global_pair_BH
  d$significant<-is.finite(qp)&qp<=alpha
  # DACT unavailable is NA, not a non-significant result. MaxP still uses all
  # eligible edges and remains explicitly a different sensitivity method.
  tg<-map_dfr(universe,function(g){
    all<-d[d$gene2==g,,drop=FALSE];z<-all[is.finite(all$DANDELION_p),,drop=FALSE];m<-nrow(z)
    ps<-sort(z$DANDELION_p)
    tibble(gene2=g,n_tested_pairs=m,n_eligible_trans_pairs=nrow(all),
      DANDELION_p=if(m)min(1,min(ps)*m)else NA_real_,
      loo_p=if(m>1)min(1,ps[2]*(m-1))else if(m==1)1 else NA_real_,
      n_selected_pairs=sum(z$significant),n_distal_loci=n_distinct(z$locus_id[z$significant]),
      gene_level_p=if(nrow(all))all$gene_level_p[1]else NA_real_,
      best_trans_p=if(nrow(all))min(all$trans_p)else NA_real_,
      maxP_target_p=if(nrow(all))min(1,min(all$maxP)*nrow(all))else NA_real_)
  })|>mutate(target_BH=p.adjust(DANDELION_p,"BH",n=length(universe)),
    target_BY=p.adjust(DANDELION_p,"BY",n=length(universe)),
    maxP_target_BH=p.adjust(maxP_target_p,"BH",n=length(universe)),
    leave_best_edge_BH=p.adjust(loo_p,"BH",n=length(universe)),
    selected=if(rule=="BY")is.finite(target_BY)&target_BY<=alpha else is.finite(target_BH)&target_BH<=alpha,
    selection_rule=rule,
    note="Target P is Bonferroni-min across ALL valid edges; target family includes every target; leave-best-edge is descriptive sensitivity")
  list(pairs=d,targets=tg)
}

le8_dandelion_plot_bundle <- function(dan,outdir) {
  tg<-dan$targets_all%||%tibble();e<-dan$evidence_plot%||%tibble()
  caption<-paste(dan$analysis_class%||%"unavailable","; no causal effect size or mediation fraction is estimated")
  if(!nrow(tg)) {
    for(pair in list(c("c2.Fig6.dandelion.png","DANDELION targets"),
      c("c2.Fig7.dandelion_evidence.png","DANDELION evidence"),c("c2.Fig9.dandelion_network.png","DANDELION network")))
      save_plot(blank_plot(pair[2],dan$status%||%"unavailable"),pair[1],12,7,outdir=outdir)
    return(invisible(NULL))
  }
  z<-tg|>arrange(target_BH)|>slice_head(n=24)|>mutate(gene2=factor(gene2,levels=rev(gene2)))
  a<-ggplot(z,aes(-log10(pmax(target_BH,1e-300)),gene2))+geom_point(aes(size=n_distal_loci,shape=selected))+
    labs(title="a. Target-level evidence (all targets corrected)",x="-log10(target BH)",y=NULL,size="Distinct upstream loci",shape="Selected")+theme_5c(9)
  b<-ggplot(z,aes(-log10(pmax(gene_level_p,1e-300)),gene2))+geom_point()+
    labs(title="b. Gene-level disease evidence",x="-log10(gene disease P)",y=NULL)+theme_5c(9)
  c<-ggplot(tg,aes(-log10(pmax(target_BH,1e-300)),-log10(pmax(maxP_target_BH,1e-300))))+
    geom_point(aes(shape=selected))+geom_abline(slope=1,intercept=0,linetype=2)+
    labs(title="c. DACT and conservative MaxP sensitivity",x="-log10(DACT target BH)",y="-log10(MaxP target BH)")+theme_5c(9)
  d<-ggplot(z,aes(-log10(pmax(leave_best_edge_BH,1e-300)),gene2))+geom_point()+
    labs(title="d. Remove the strongest edge (aggregation sensitivity)",x="-log10(leave-best-edge BH)",y=NULL)+theme_5c(9)
  save_plot((a|b)/(c|d)+plot_annotation(caption=caption),"c2.Fig6.dandelion.png",17,12,outdir=outdir)
  if(nrow(e)) {
    a<-ggplot(e,aes(-log10(pmax(trans_p,1e-300)),-log10(pmax(gene_level_p,1e-300))))+
      geom_point(aes(color=significant),alpha=.4)+labs(title="a. Regulatory and gene-disease evidence are separate",x="-log10(trans-pQTL P)",y="-log10(gene disease P)")+theme_5c(9)
    b<-dan$exposure_qc|>ggplot(aes(n_targets_tested,n_selected))+geom_point()+
      labs(title="b. Broad-regulator audit (not an automatic exclusion)",x="Tested targets per locus",y="Selected edges")+theme_5c(9)
    save_plot((a|b)+plot_annotation(caption=caption),"c2.Fig7.dandelion_evidence.png",16,8,outdir=outdir)
  }
  edges<-dan$pairs%||%tibble();edges<-edges|>filter(gene2%in%as.character(z$gene2))
  if(nrow(edges)) {
    genes<-unique(edges$gene2);loci<-unique(edges$locus_id)
    edges$yg<-match(edges$gene2,genes)/max(1,length(genes));edges$yl<-match(edges$locus_id,loci)/max(1,length(loci))
    a<-ggplot(edges)+geom_segment(aes(x=0,xend=1,y=yl,yend=yg,alpha=-log10(pmax(global_pair_BH,1e-300))))+
      geom_text(data=unique(edges[,c("locus_id","yl")]),aes(x=0,y=yl,label=locus_id),hjust=1,size=2.6)+
      geom_text(data=unique(edges[,c("gene2","yg")]),aes(x=1,y=yg,label=gene2),hjust=0,size=2.8)+
      xlim(-.8,1.8)+theme_void()+labs(title="Disease loci → target genes",subtitle="SNP-to-upstream-gene labels are annotations, not proven effector genes",alpha="Edge evidence")
  }else a<-blank_plot("DANDELION network","No edge passed the configured full-family threshold")
  save_plot(a+plot_annotation(caption=caption),"c2.Fig9.dandelion_network.png",15,10,outdir=outdir)
}
run_dandelion_step <- function(layer,assoc,ann,base,ygfile,rawdir,outdir,
                              top_candidates=tibble(),mode=RUN_Dandelion) {
  empty<-list(status="not_run",pairs=tibble(),targets=tibble(),targets_all=tibble(),gene_pairs=tibble(),
    lead_snps=tibble(),snp_gene_map=tibble(),evidence_plot=tibble(),exposure_qc=tibble(),qtl_audit=tibble(),
    input_audit=tibble(),integration=tibble(),sig_gene2=character(),non_sig_gene2=character(),
    primary_eligible=FALSE,consolidation_eligible=FALSE,broad_selection_warning=FALSE,target_fraction=NA_real_,
    gene_level_file=NA_character_,gene_evidence_type="unavailable",network_file=NA_character_,analysis_class="not run")
  if(mode=="None"||layer!="protein")return(modifyList(empty,list(status=if(layer!="protein")"protein-only adaptation"else"RUN_Dandelion=None")))
  input<-find_dandelion_gene_file(ygfile)
  if(is.na(input$path))return(modifyList(empty,list(status="gene-level disease input unavailable",analysis_class=input$analysis_class)))
  gene_map<-le8_assay_genes(ann$feature)
  ref<-ann|>mutate(gene_name=unname(gene_map[feature]))|>
    filter(!is.na(chr),is.finite(start),is.finite(end))|>
    arrange(feature)|>distinct(gene_name,.keep_all=TRUE) # prespecified lexicographic representative, never minimum P
  p_gene<-read_gene_level_p(input$path,ref$gene_name)
  ref<-ref|>filter(gene_name%in%names(p_gene))
  if(mode=="Top")ref<-ref|>filter(feature%in%top_candidates$feature)
  if(DANDELION_MAX_GENE2>0)ref<-head(ref,DANDELION_MAX_GENE2) # not ranked by gene-disease P
  if(nrow(ref)<5)return(modifyList(empty,list(status="fewer than five matched annotated targets",gene_level_file=input$path)))
  lead<-find_disease_leads(ygfile);src<-attr(lead,"source")%||%"unknown"
  le8_dependency(src)
  lead<-lead|>filter(is.finite(P),P<=DANDELION_GWS,is.finite(POS),!is.na(CHR))|>distinct(SNP,.keep_all=TRUE)
  # Attempt explicit LD validation of the selected disease variants.
  fake<-lead;fake$source_file<-ygfile;fake$BETA<-1;fake$SE<-1
  ld<-le8_ld_for_iv(fake,sub("\\.gz$","",basename(ygfile)))
  ld_verified<-!is.null(ld$R)&&all(lead$SNP%in%rownames(ld$R))
  if(ld_verified)lead<-le8_ld_greedy(lead,ld$R,le8_num_env("C2_MR_LD_R2",.001))
  if(nrow(lead)<2)return(modifyList(empty,list(status="fewer than two eligible disease loci",gene_level_file=input$path)))
  lead$locus_id<-paste0("chr",sub("^chr","",lead$CHR),":",lead$POS)
  if(!ld_verified) {
    # Distance groups provide descriptive counts only; not called LD-independent.
    for(ch in unique(lead$CHR)) {
      ix<-which(lead$CHR==ch);ix<-ix[order(lead$POS[ix])];g<-cumsum(c(TRUE,diff(lead$POS[ix])>DANDELION_LEAD_BP))
      for(k in unique(g)){jj<-ix[g==k];lead$locus_id[jj]<-paste0("chr",sub("^chr","",ch),":",min(lead$POS[jj]),"-",max(lead$POS[jj]))}
    }
  }
  pt<-build_trans_p_matrix(ref$feature,lead,base);qa<-attr(pt,"qtl_audit")
  rownames(pt)<-ref$gene_name[match(rownames(pt),ref$feature)]
  refm<-ref|>transmute(gene_name,type="protein_coding",Chromosome=paste0("chr",sub("^chr","",chr)),start,end)
  pw<-p_gene[rownames(pt)]
  # Explicit local exclusion before all tests; missing tests remain NA, not P=1.
  for(i in seq_len(nrow(pt))){r<-refm[match(rownames(pt)[i],refm$gene_name),]
    local<-paste0("chr",sub("^chr","",lead$CHR))==r$Chromosome&lead$POS>=r$start-DANDELION_CIS_BP&lead$POS<=r$end+DANDELION_CIS_BP
    pt[i,local]<-NA_real_}
  med<-NULL;package_status<-"package unavailable; MaxP sensitivity only"
  if(requireNamespace("DANDELION",quietly=TRUE)) {
    med<-tryCatch(le8_safe_call(DANDELION::med_gene,list(p.trans=pt,p.wes=pw,ref.table=as.data.frame(refm),
      gene1.list=colnames(pt),gene1.type="SNP",SNP.ref=as.data.frame(lead|>transmute(SNP,SNPPos=POS,SNPChr=CHR)),
      target.fdr=DANDELION_FDR,dist=DANDELION_CIS_BP,n.cores=N_CORES,verbose=FALSE)),error=function(e)e)
    package_status<-if(inherits(med,"condition"))conditionMessage(med)else"ok"
    if(inherits(med,"condition"))med<-NULL
  }
  dact<-matrix(NA_real_,nrow(pt),ncol(pt),dimnames=dimnames(pt));native<-matrix(FALSE,nrow(pt),ncol(pt),dimnames=dimnames(pt))
  if(!is.null(med)) {
    rr<-intersect(rownames(pt),rownames(med$mat.p));cc<-intersect(colnames(pt),colnames(med$mat.p))
    dact[rr,cc]<-med$mat.p[rr,cc];native[rr,cc]<-med$mat.sig[rr,cc]!=0
  }
  evidence<-tibble(gene2=rep(rownames(pt),times=ncol(pt)),rsid=rep(colnames(pt),each=nrow(pt)),
    trans_p=as.vector(pt),gene_level_p=rep(as.numeric(pw),times=ncol(pt)),DANDELION_p=as.vector(dact),
    package_selected=as.vector(native))|>filter(is.finite(trans_p),is.finite(gene_level_p))|>
    left_join(lead|>select(rsid=SNP,locus_id,CHR,POS,outcome_P=P),by="rsid")
  fam<-le8_dandelion_family(evidence,rownames(pt),DANDELION_FDR);ev<-fam$pairs;alltg<-fam$targets
  primary<-input$primary_eligible&&mode=="All"&&ld_verified&&package_status=="ok"
  class<-paste(input$analysis_class,if(mode=="Top")"; selected Top family is exploratory"else"; all eligible target genes",
    if(!ld_verified)"; disease-locus LD unverified"else"; disease-locus LD verified")
  alltg<-alltg|>left_join(ref|>select(gene2=gene_name,feature),by="gene2")|>
    mutate(gene_evidence_type=input$evidence_type,primary_eligible=primary,analysis_class=class,
      consolidation_eligible=primary&selected,gene_level_bonferroni=gene_level_p<.05/max(1,attr(p_gene,"n_gene_tests")%||%length(p_gene)),
      dpg_class=case_when(!is.finite(DANDELION_p)~"DACT not tested / unavailable",
        selected&gene_level_bonferroni~"DACT + gene-disease evidence",
        selected~"DACT-prioritized hypothesis",gene_level_bonferroni~"Gene-disease evidence only",
        TRUE~"Not selected"))
  tg<-alltg|>filter(selected);fraction<-nrow(tg)/max(1,nrow(alltg));broad<-fraction>DANDELION_MAX_TARGET_FRACTION
  # Broad selection triggers an audit flag, not a data-dependent veto threshold.
  for(nm in c("alltg","tg")){z<-get(nm);z$target_fraction<-fraction;z$broad_selection_warning<-broad;assign(nm,z)}
  pd<-ev|>filter(significant)|>left_join(alltg|>select(gene2,feature,primary_eligible,consolidation_eligible,analysis_class,gene_evidence_type),by="gene2")
  mp<-read_dandelion_snp_gene_map(lead)
  if(nrow(mp))pd<-left_join(pd,as_tibble(mp)|>group_by(SNP)|>summarise(gene1=paste(unique(GeneSymbol),collapse=";"),.groups="drop"),by=c("rsid"="SNP"))
  if(!"gene1"%in%names(pd))pd$gene1<-pd$rsid
  pd$gene1[is.na(pd$gene1)]<-pd$rsid[is.na(pd$gene1)]
  gp<-pd|>mutate(region=locus_id) # locus is the upstream unit; labels do not establish a causal gene.
  qc<-ev|>group_by(rsid,locus_id)|>summarise(n_targets_tested=n(),n_selected=sum(significant),selected_fraction=mean(significant),.groups="drop")
  audit<-tibble(metric=c("mode","gene_file","gene_evidence","analysis_class","primary_eligible","disease_LD_verified",
    "DACT_package_status","target_universe","valid_trans_tests","selected_targets","target_fraction","multiplicity","scope_boundary"),
    value=as.character(c(mode,input$path,input$evidence_type,class,primary,ld_verified,package_status,nrow(alltg),nrow(ev),nrow(tg),fraction,
      paste0("All finite DACT pairs; target Bonferroni-min then BH/BY over all targets; selection=",Sys.getenv("C2_DANDELION_MULTIPLICITY",unset="BY")),
      "trans-pQTL adaptation; neither causal direction nor protein mediation fraction is identified")))
  ans<-modifyList(empty,list(status=if(package_status=="ok")"ok"else"DACT unavailable; see MaxP sensitivity",mode=mode,
    result=med,pairs=pd,targets=tg,targets_all=alltg,gene_pairs=gp,lead_snps=lead,snp_gene_map=mp,evidence_plot=ev,
    exposure_qc=qc,qtl_audit=qa,input_audit=audit,sig_gene2=tg$gene2[tg$gene_level_bonferroni%in%TRUE],
    non_sig_gene2=tg$gene2[!tg$gene_level_bonferroni%in%TRUE],gene_level_file=input$path,gene_evidence_type=input$evidence_type,
    analysis_class=class,primary_eligible=primary,consolidation_eligible=primary,target_fraction=fraction,
    broad_selection_warning=broad,ptrans_dimensions=dim(pt)))
  tabs<-list(pairs=pd,targets=tg,targets_all=alltg,gene_pairs=gp,lead_snps=lead,snp_to_cis_gene=mp,
    exposure_qc=qc,qtl_coverage=qa,input_audit=audit)
  for(nm in names(tabs))write_raw_csv(tabs[[nm]],paste0("c2.dandelion_",nm,".csv"),rawdir)
  data.table::fwrite(ev,file.path(rawdir,"c2.dandelion_tested_pairs.csv.gz"),compress="gzip")
  le8_dandelion_plot_bundle(ans,outdir)
  invisible(tabs) # appended after the legacy workbook writer by the revision hook
  ans
}
# Existing C2 invokes the old plot functions again. Keep the canonical figures
# produced by le8_dandelion_plot_bundle rather than overwrite them with old p-values.
plot_dandelion_results <- function(tg,pd,gene_pair,outdir) invisible(NULL)
plot_dandelion_landscape <- function(evidence_plot,tg,outdir,gene_evidence_type) invisible(NULL)
write_dandelion_package_network <- function(gene_pair,p_gene,rawdir,outdir) NA_character_
plot_dandelion_mr_integration <- function(dandelion,mr,assoc,outdir) {
  tg<-dandelion$targets_all%||%tibble();if(!nrow(tg)) {
    save_plot(blank_plot("DANDELION and MR","DANDELION target tests unavailable"),"c2.Fig8.dandelion_mr_integration.png",14,8,outdir=outdir)
    return(tibble())
  }
  m<-mr|>filter(analysis=="cis")|>select(feature=exposure,MR_beta=b,MR_p=pval,MR_FDR=FDR_all)
  z<-tg|>left_join(m,by="feature")|>left_join(assoc|>select(feature=term,observed_beta=beta,observed_p=p.value),by="feature")
  a<-ggplot(z,aes(-log10(pmax(target_BH,1e-300)),-log10(pmax(MR_FDR,1e-300))))+
    geom_point(aes(shape=primary_eligible))+labs(title="a. Orthogonal evidence shown side by side",x="-log10(DACT target BH)",y="-log10(cis-MR BH)",shape="WES/All/LD-eligible")+theme_5c(9)
  b<-z|>arrange(target_BH)|>slice_head(n=20)|>ggplot(aes(observed_beta,reorder(gene2,observed_beta)))+
    geom_point()+geom_vline(xintercept=0)+labs(title="b. Observed association does not determine regulatory direction",x="Observed incident log HR",y=NULL)+theme_5c(9)
  save_plot(a|b,"c2.Fig8.dandelion_mr_integration.png",16,9,outdir=outdir);z
}

read_dandelion_snp_gene_map <- function(lead) {
  f<-Sys.getenv("C2_DANDELION_SNP_GENE_MAP",unset="")
  if(!nzchar(f)||!file.exists(f))return(tibble(SNP=character(),GeneSymbol=character()))
  le8_dependency(f);d<-data.table::fread(f,showProgress=FALSE)
  sc<-pick_col_local(names(d),c("^SNP$","^RSID$"));gc<-pick_col_local(names(d),c("^GENESYMBOL$","^GENE_SYMBOL$","^GENE$","^SYMBOL$"))
  if(is.na(sc)||is.na(gc))stop("SNP gene map needs SNP and GeneSymbol columns")
  tibble(SNP=as.character(d[[sc]]),GeneSymbol=as.character(d[[gc]]))|>
    filter(SNP%in%lead$SNP,!is.na(GeneSymbol),nzchar(GeneSymbol))|>distinct()
}

# Keep an invalid optional input local to DANDELION; preserve completed MR.
.le8_dandelion_impl<-run_dandelion_step
run_dandelion_step <- function(layer,assoc,ann,base,ygfile,rawdir,outdir,
                              top_candidates=tibble(),mode=RUN_Dandelion) {
  tryCatch(.le8_dandelion_impl(layer,assoc,ann,base,ygfile,rawdir,outdir,top_candidates,mode),
    error=function(e) {
      msg<-conditionMessage(e);warning("DANDELION input/method failure: ",msg,call.=FALSE)
      au<-tibble(status="failed",message=msg)
      write_raw_csv(au,"c2.dandelion_failure.csv",rawdir)
      list(status=paste("failed:",msg),pairs=tibble(),targets=tibble(),targets_all=tibble(),
        gene_pairs=tibble(),evidence_plot=tibble(),exposure_qc=tibble(),qtl_audit=tibble(),
        input_audit=au,lead_snps=tibble(),snp_gene_map=tibble(),integration=tibble(),
        gene_level_file=NA_character_,gene_evidence_type="unavailable",network_file=NA_character_,
        primary_eligible=FALSE,analysis_class="Failed optional analysis; not a negative finding")
    })
}
