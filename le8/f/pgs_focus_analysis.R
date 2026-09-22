# Called by c1_pgs.R; independent pgs_focus entry avoids rerunning MR/ML.
pgs_focus_anchors <- function(layer) {
  default<-if(layer=="protein")"PCSK9,LPA,GDF15,NTPROBNP,MMP12,CCL19,CCL21,CXCL13,AGER,IL15,TNFRSF4" else
    "L_VLDL_TG.pct,L_VLDL_TG,Total_TG,ApoB,VLDL_size,GlycA,Albumin,Phe,Lactate,DHA,LA.pct"
  pgs_csv(if(layer=="protein")"PGS_PROT_ANCHORS"else"PGS_MET_ANCHORS",default)
}
pgs_focus_selection <- function(summary,layer) {
  d<-summary;ranked<-if(nrow(d))d$feature[order(!(d$evidence_pattern=="Both supported: opposite"),
    pmax(d$measured_FDR,d$pgs_FDR),d$measured_p,na.last=TRUE)]else character()
  unique(c(intersect(pgs_focus_anchors(layer),d$feature),head(ranked,pgs_num("PGS_DEEP_MAX",24))))
}
pgs_focus_export <- function(tables,root,prefix="c1.pgs_focus") {
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  for(nm in names(tables)) {
    d<-tables[[nm]];if(!is.data.frame(d))next
    if(!ncol(d))d<-data.frame(status="unavailable / not estimated")
    data.table::fwrite(d,file.path(root,paste0(prefix,".",nm,".csv")))
  }
  wb<-openxlsx::createWorkbook()
  for(nm in names(tables)) {
    d<-tables[[nm]];if(!is.data.frame(d))next
    if(!ncol(d))d<-data.frame(status="unavailable / not estimated")
    # Excel has a hard row limit; CSV retains all rows if a large family exceeds it.
    if(nrow(d)>1048500)d<-data.frame(status="Full table exceeds Excel row limit; see companion CSV",rows=nrow(d))
    sh<-substr(nm,1,31);openxlsx::addWorksheet(wb,sh)
    if(nrow(d))openxlsx::writeDataTable(wb,sh,d)else openxlsx::writeData(wb,sh,d)
    openxlsx::freezePane(wb,sh,firstRow=TRUE);openxlsx::setColWidths(wb,sh,seq_len(ncol(d)),18)
  }
  openxlsx::saveWorkbook(wb,file.path(root,paste0(prefix,".xlsx")),overwrite=TRUE)
}
pgs_focus_deep <- function(d,covars,feature,layer,boot=0L) {
  k<-as.integer(pgs_num("PGS_FOLDS",5));seed<-as.integer(pgs_num("SEED",2026))
  cal<-pgs_calibrate(d,covars,k,seed);cal$audit$feature<-feature;cal$folds$feature<-rep(feature,nrow(cal$folds))
  main<-pgs_component_models(cal$data,covars,feature)
  cuts<-sort(unique(as.numeric(pgs_csv("PGS_TIME_CUTS","0,2,5,10,Inf"))))
  if(length(cuts)<2||anyNA(cuts)||cuts[1]!=0||any(diff(cuts)<=0))stop("Invalid PGS_TIME_CUTS")
  windows<-lapply(seq_len(length(cuts)-1L),function(i)pgs_component_models(cal$data,covars,feature,landmark=cuts[i],end=cuts[i+1]))
  landmarks<-lapply(c(2,5),function(L)pgs_component_models(cal$data,covars,feature,landmark=L))
  sens<-list();sensitivity_sets<-list(LE4=unique(c(covars,pgs_csv("PGS_LE4_COVARS","diet.pts,pa.pts,smoke.pts,sleep.pts"))),
    LE8=unique(c(covars,get0("vars.le8",ifnotfound=character()))),
    treatment=unique(c(covars,pgs_csv("PGS_TREATMENT_COVARS","drug.lipid,drug.htn,drug.dm"))))
  for(nm in names(sensitivity_sets)) {
    cv<-sensitivity_sets[[nm]];z<-d[pgs_complete(d,c(".m",".g",cv)),,drop=FALSE]
    if(identical(cv,covars))next
    # Both models use the SAME restricted sample; no silent covariate dropping.
    a<-pgs_pair(z,covars,feature,paste0(nm,"_restricted_basic"))
    b<-pgs_pair(z,cv,feature,paste0(nm,"_adjusted"))
    sens[[nm]]<-pgs_bind(list(a$summary,b$summary))
  }
  # Outlier influence sensitivity is separate from the primary preprocessing.
  w<-d;ii<-w$.prev%in%0&is.finite(w$.m)
  if(sum(ii)>200){lim<-quantile(w$.m[ii],c(.01,.99));w$.m<-pmin(pmax(w$.m,lim[1]),lim[2]);sens$winsor<-pgs_pair(w,covars,feature,"measured_winsor_1_99")$summary}
  life<-pgs_lifestyle(cal$data,covars,intersect(get0("vars.le8",ifnotfound=character()),names(d)),feature)
  list(calibration=cal$audit,folds=cal$folds,components=main$effects,contrasts=main$contrast,
    prevalent=pgs_prevalent(cal$data,covars,feature),
    windows=pgs_bind(lapply(windows,`[[`,"effects")),window_contrasts=pgs_bind(lapply(windows,`[[`,"contrast")),
    landmarks=pgs_bind(lapply(landmarks,`[[`,"effects")),landmark_contrasts=pgs_bind(lapply(landmarks,`[[`,"contrast")),
    sensitivity=pgs_bind(sens),lifestyle=life,
    bootstrap=pgs_bootstrap(d,covars,feature,B=boot,k=k,seed=seed))
}
run_c1_pgs_focus <- function(layer,outdir) {
  rd<-file.path(outdir,"c1_correlate");dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  sf<-find_c1_pgs_file(layer)
  if(is.na(sf)) {
    out<-list(status=data.frame(status="unavailable",detail="Biomarker PGS input missing"))
    pgs_focus_export(out,rd);saveRDS(out,file.path(rd,"c1.pgs_focus.rds"));stop("Biomarker PGS input missing for ",layer,"; focused analysis stopped")
  }
  primary<-pgs_csv("PGS_BASIC_COVARS",paste(if(length(get0("le8_custom_covars",ifnotfound=character())))le8_custom_covars else vars.basic,collapse=","))
  # A supplied genetic score is not a disease PRS, even when its name is similar.
  scores<-pgs_ids(read_c1_pgs(sf),"PGS");biom<-pgs_ids(if(layer=="protein")read_prot()else read_met(),"Measured omics")
  feats<-setdiff(names(biom),"eid");mp<-map_c1_pgs_columns(feats,names(scores))
  only<-pgs_csv("PGS_FEATURES");if(length(only))mp<-mp[intersect(names(mp),only)]
  if(!length(mp))stop("No exact measured/PGS feature matches for ",layer)
  ph<-pgs_ids(read_all(),"Phenotypes")|>filter_analysis_cohort()|>make_outcome(Y)
  if(length(setdiff(primary,names(ph))))stop("Missing primary PGS covariates: ",paste(setdiff(primary,names(ph)),collapse=","))
  ph$.prev<-make_prevalent_status(ph,Y);ph$.time<-ph[[paste0(Y,".t2e")]];ph$.event<-ph[[paste0(Y,".Yt2e")]]
  # A diagnosis count without an onset date is an unresolved baseline state,
  # not a healthy control for calibration or prevalent analysis.
  dc<-if(Y%in%names(ph))Y else paste0("icd10Ct_",Y)
  yd<-paste0("fod_icd10_",Y)
  if(all(c(dc,yd)%in%names(ph))) {
    positive<-if(inherits(ph[[dc]],"Date"))!is.na(ph[[dc]])else suppressWarnings(as.numeric(ph[[dc]]))>0
    unknown<-is.na(ph[[yd]])&positive%in%TRUE;ph$.prev[unknown]<-NA_real_
  }
  groupcol<-Sys.getenv("PGS_GROUP_COLUMN","")
  if(nzchar(groupcol)&&!groupcol%in%names(ph))stop("PGS_GROUP_COLUMN missing: ",groupcol)
  ph$.group<-if(nzchar(groupcol))as.character(ph[[groupcol]])else ph$eid
  if(anyNA(ph$.group)||any(!nzchar(ph$.group)))stop("Missing PGS group IDs; do not silently split related participants")
  need<-unique(c("eid",".prev",".time",".event",".group",primary,vars.le8,
    pgs_csv("PGS_TREATMENT_COVARS","drug.lipid,drug.htn,drug.dm")))
  base<-ph[,intersect(need,names(ph)),drop=FALSE];rm(ph);invisible(gc())
  ids<-intersect(intersect(base$eid,biom$eid),scores$eid)
  base<-base[match(ids,base$eid),,drop=FALSE]
  im<-match(ids,biom$eid);ig<-match(ids,scores$eid)
  base$.fold<-pgs_folds(base$.group,as.integer(pgs_num("PGS_FOLDS",5)),as.integer(pgs_num("SEED",2026)))
  inputs<-c(sf,file.path(indir,"Rdata",c("all.rds",if(layer=="protein")"prot.rds"else"met.rds")),
    Sys.getenv("LE8_BASELINE_MED_FILE",file.path(indir,"rap/vip.tab.gz")))
  sig<-pgs_hash(list(version=PGS_FOCUS_VERSION,files=pgs_stamp(inputs),Y=Y,primary=primary,
    baseline=LE8_BASELINE_VERSION,options=le8_analysis_options(),features=names(mp),group=groupcol,
    end=as.character(date_follow_end),folds=base$.fold,settings=Sys.getenv(c("PGS_FOLDS","PGS_DEEP_MAX","PGS_MIN_N","PGS_MIN_EVENTS","PGS_BOOT","PGS_BOOT_FEATURES","PGS_TIME_CUTS","PGS_LE4_COVARS","PGS_TREATMENT_COVARS","LE8_BASELINE_MED_COLUMNS","LE8_BASELINE_MAP","PGS_PROT_ANCHORS","PGS_MET_ANCHORS","PGS_MIN_PARTIAL_R2","PGS_GWAS_SOURCE","PGS_GWAS_OVERLAP"))))
  cache<-file.path(rd,".pgs_focus_cache",sig);dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  getdata<-function(f){d<-base;d$.m<-suppressWarnings(as.numeric(biom[[f]][im]));d$.g<-suppressWarnings(as.numeric(scores[[mp[[f]]]][ig]));d}
  message("[LE8] START PGS/",layer," matched scan | ",length(mp)," features; ",nrow(base)," matched participants")
  scanfile<-file.path(cache,"scan.rds");pairs<-if(file.exists(scanfile)&&!LE8_REPLACE)readRDS(scanfile)else {
    ans<-lapply(names(mp),function(f)pgs_pair(getdata(f),primary,f));saveRDS(ans,scanfile,compress=FALSE);ans}
  paired<-pgs_bind(lapply(pairs,`[[`,"summary"));effects<-pgs_bind(lapply(pairs,`[[`,"effects"));rm(pairs)
  paired$measured_FDR<-p.adjust(paired$measured_p,"BH");paired$pgs_FDR<-p.adjust(paired$pgs_p,"BH")
  paired$pgs_joint_FDR<-p.adjust(paired$pgs_joint_p,"BH")
  # Intersection-union test for opposite directions, with two possible
  # orientations corrected before BH over the full matched scan.
  zm<-paired$measured_beta/paired$measured_se;zg<-paired$pgs_beta/paired$pgs_se
  paired$opposite_conjunction_p<-pmin(1,2*pmin(pmax(pnorm(-zm),pnorm(zg)),pmax(pnorm(zm),pnorm(-zg))))
  paired$opposite_conjunction_FDR<-p.adjust(paired$opposite_conjunction_p,"BH")
  paired$opposite_conjunction_supported<-is.finite(paired$opposite_conjunction_FDR)&paired$opposite_conjunction_FDR<.05
  paired$evidence_pattern<-pgs_classify(paired$measured_beta,paired$pgs_beta,paired$measured_FDR,paired$pgs_FDR)
  effects<-pgs_adjust(effects,groups=c("model","term"))
  candidates<-pgs_focus_selection(paired,layer);paired$deep_selected<-paired$feature%in%candidates
  selection<-data.frame(feature=candidates,selection=ifelse(candidates%in%pgs_focus_anchors(layer),"declared anchor","exploratory same-sample ranking"),
    caveat="Same data used for candidate selection; candidate-only FDR is not independent confirmation")
  defaultboot<-if(layer=="protein")"PCSK9,LPA,CCL19,CCL21"else"L_VLDL_TG.pct,GlycA,Lactate"
  bootfeatures<-pgs_csv("PGS_BOOT_FEATURES",defaultboot);B<-as.integer(pgs_num("PGS_BOOT",100))
  message("[LE8] DONE PGS/",layer," matched scan | opposite=",sum(paired$evidence_pattern=="Both supported: opposite"))
  message("[LE8] START PGS/",layer," component analysis | ",length(candidates)," candidates; bootstrap=",B," for ",length(intersect(candidates,bootfeatures))," anchors")
  deep<-lapply(candidates,function(f){
    cf<-file.path(cache,paste0("feature_",pgs_hash(f),".rds"))
    if(file.exists(cf)&&!LE8_REPLACE)return(readRDS(cf))
    z<-pgs_focus_deep(getdata(f),primary,f,layer,if(f%in%bootfeatures)B else 0L)
    tmp<-paste0(cf,".tmp");saveRDS(z,tmp,compress=FALSE);if(!file.rename(tmp,cf))stop("Cannot save PGS feature checkpoint");z
  })
  nm<-unique(unlist(lapply(deep,names)));out<-setNames(lapply(nm,function(n)pgs_bind(lapply(deep,`[[`,n))),nm)
  out$paired<-paired;out$matched_models<-effects;out$selection<-selection
  reference_file<-file.path(rd,paste0(if(layer=="protein")"pwas"else"mwas","_pgs_incident_full_genetic.csv"))
  if(file.exists(reference_file)) {
    ref<-as.data.frame(data.table::fread(reference_file,showProgress=FALSE))
    if(nrow(ref)){ref$comparison_role<-"Historical full-genetic-cohort reference; not matched to current measured analysis; see original adjustment metadata"
      out$full_cohort_reference<-ref}
  }
  if(layer=="metabolite")out$composition<-pgs_composition_analysis(base,biom,im,mp,getdata,primary)
  for(n in intersect(c("components","windows","landmarks","prevalent"),names(out)))out[[n]]<-pgs_adjust(out[[n]],groups=intersect(c("scope","model","term","landmark","end"),names(out[[n]])))
  for(n in intersect(c("contrasts","window_contrasts","landmark_contrasts"),names(out)))out[[n]]<-pgs_adjust(out[[n]],groups=intersect(c("scope","landmark","end"),names(out[[n]])))
  out$lifestyle<-pgs_adjust(out$lifestyle,groups=intersect(c("part"),names(out$lifestyle)))
  out$status<-data.frame(status="ok",version=PGS_FOCUS_VERSION,signature=sig,layer,outcome=Y,N_matched=nrow(base),features=length(mp),estimable_pairs=sum(is.finite(paired$measured_p)&is.finite(paired$pgs_p)),deep_features=length(candidates),
    adjustment=paste(primary,collapse=";"),PGS_source=sf,GWAS_source=Sys.getenv("PGS_GWAS_SOURCE","unknown"),
    GWAS_overlap=Sys.getenv("PGS_GWAS_OVERLAP","unknown"),group_column=if(nzchar(groupcol))groupcol else"eid; family structure not supplied",
    interpretation="Association decomposition, not causal partition or biomarker concentration at birth")
  # Source-specific PGSs are built only with verified build/alleles/genotypes.
  if(exists("pgs_source_analysis",mode="function")) {
    src<-pgs_source_analysis(layer,outdir,candidates,base,biom,im,primary,paired,getdata)
    out<-c(out,src)
  }
  pgs_focus_export(out,rd);saveRDS(out,file.path(rd,"c1.pgs_focus.rds"),compress="gzip")
  message("[LE8] DONE PGS/",layer," component analysis")
  out
}

pgs_composition_analysis <- function(base,biom,im,mp,getdata,covars) {
  targets<-intersect(pgs_csv("PGS_COMPOSITION_TARGETS","L_VLDL_TG.pct"),names(mp))
  controls<-pgs_csv("PGS_COMPOSITION_CONTROLS","L_VLDL_TG,Total_TG,ApoB,VLDL_size")
  rows<-list()
  for(f in targets)for(v in controls)for(kind in c("measured","pgs")) {
    tag<-paste(f,v,kind)
    if(!v%in%names(biom)||(kind=="pgs"&&!v%in%names(mp))) {
      rows[[tag]]<-data.frame(feature=f,control=v,control_type=kind,status="control unavailable");next
    }
    d<-getdata(f);d$.burden<-if(kind=="measured")as.numeric(biom[[v]][im])else getdata(v)$.g
    d<-d[pgs_complete(d,c(".m",".g",".burden",covars)),,drop=FALSE]
    a<-pgs_pair(d,covars,f,"burden_restricted_basic")$summary
    b<-pgs_pair(d,c(covars,".burden"),f,"burden_conditioned")$summary
    z<-pgs_bind(list(a,b));z$control<-v;z$control_type<-kind
    z$status<-"Conditional sensitivity; potential mediator/collider adjustment, not a causal direct effect"
    rows[[tag]]<-z
  }
  pgs_bind(rows)
}
