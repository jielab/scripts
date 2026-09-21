# Focused, prospective test of LE8 supervision and use of prevalent (Yang) donors.
# All risk models use incident training participants only. Yang enters proxy learning.
suppressPackageStartupMessages({
  fdir<-Sys.getenv("LE8_FDIR",unset="/mnt/d/scripts/le8/f")
  source(file.path(fdir,"comm.f.R"))
  source(file.path(fdir,"c5_revision_validation.R"))
  source(file.path(fdir,"c1_pgs.R"))
  source(file.path(fdir,"c4_pgs_bridge.R"))
  source(file.path(fdir,"c4_focus_extensions.R"))
})
LE8_JOB<-"c4_focus"
# Edit f/le8_budget_config.R; C4_FOCUS_BUDGETS=5,10,50 overrides one run.
C4_FOCUS_BUDGETS<-LE8_ASSAY_BUDGETS

# Screening and proxy discovery do not depend on assay budgets/bootstrap count.
# Keep their cache keys separate so changing c(5,10,50) only refits the panels.
focus_stage_signature <- function(name,layer,inputs) {
  fi<-file.info(inputs)
  funcs<-c("cox_scan","focus_split","make_outcome","t2e","filter_analysis_cohort",
    "le8_training_connection_set","le8_select_phenotypes")
  le8_hash_object(list(stage=name,layer=layer,outcome=Y,inputs=inputs,size=fi$size,mtime=fi$mtime,
    options=le8_analysis_options(),clinical=covs_use,basic=vars.basic,components=vars.le8,
    follow_end=date_follow_end,seed=SEED,proxy_assignment="fixed Yin halves v1",
    requested_components=Sys.getenv("C4_FOCUS_COMPONENTS",paste(vars.le8,collapse=",")),
    clinical_default=unique(c(vars.basic,"smoke.pts","sbp","hba1c_ngsp","bmi","nonhdl")),
    algorithms=lapply(funcs,function(n)deparse(get(n,mode="function"))),
    cap=Sys.getenv("C5_CONNECTION_MAX_N",unset="60000")))
}

focus_split <- function(dat,evar,seed=SEED) {
  set.seed(seed)
  strata<-split(seq_len(nrow(dat)),as.character(dat[[evar]]),drop=TRUE)
  val<-unlist(lapply(strata,function(ii)sample(ii,max(1L,round(length(ii)*.2)))),use.names=FALSE)
  ifelse(seq_len(nrow(dat))%in%val,"validation","training")
}

focus_proxy_accuracy <- function(train,test,features,components,covars,model) {
  map_dfr(components,function(cmp) {
    tr<-train[is.finite(train[[cmp]]),,drop=FALSE];te<-test[is.finite(test[[cmp]]),,drop=FALSE]
    if(nrow(tr)<100||nrow(te)<50)return(tibble())
    predict_one<-function(vars) {
      x<-le8_prepare_prediction_matrix(tr,te,vars)
      xx<-cbind(Intercept=1,x$train);xt<-cbind(Intercept=1,x$test)
      b<-lm.fit(xx,tr[[cmp]])$coefficients
      if(anyNA(b))return(rep(NA_real_,nrow(te)))
      drop(xt%*%b)
    }
    y<-te[[cmp]];mu<-mean(tr[[cmp]]);den<-sum((y-mu)^2)
    base<-predict_one(covars);full<-predict_one(unique(c(covars,features)))
    only<-predict_one(features)
    tibble(model,component=cmp,N=nrow(te),R2_omics=1-sum((y-only)^2)/den,
      R2_basic=1-sum((y-base)^2)/den,R2_basic_omics=1-sum((y-full)^2)/den,
      delta_R2=(sum((y-base)^2)-sum((y-full)^2))/den,
      note="Common Yin-trained held-out LE8 reconstruction; not intervention response")
  })
}

focus_inflammation <- function(train,test,features) {
  markers<-intersect(le8_csv_env("C4_FOCUS_INFLAMMATION","GDF15,IL6,CRP"),features)
  score<-function(d,ref) {
    x<-sapply(markers,function(v)(as.numeric(d[[v]])-ref[v,"center"])/ref[v,"sd"])
    if(is.null(dim(x)))x<-matrix(x,ncol=length(markers))
    ans<-rowMeans(x);ans[rowSums(is.finite(x))!=length(markers)]<-NA_real_;ans
  }
  if(length(markers)<2L)return(list(group=rep("unavailable",nrow(test)),audit=tibble(status="fewer than two prespecified inflammatory markers")))
  ref<-t(vapply(markers,function(v)c(center=mean(train[[v]],na.rm=TRUE),sd=sd(train[[v]],na.rm=TRUE)),numeric(2)))
  if(any(!is.finite(ref))||any(ref[,"sd"]<=0))return(list(group=rep("unavailable",nrow(test)),audit=tibble(status="invalid training marker distribution")))
  tr<-score(train,ref);te<-score(test,ref);cut<-median(tr[is.finite(tr)])
  group<-ifelse(!is.finite(te),"unavailable",ifelse(te<=cut,"Low baseline inflammation","High baseline inflammation"))
  list(group=group,audit=tibble(marker=markers,center=ref[,"center"],sd=ref[,"sd"],cutoff=cut,
    definition="Complete-marker mean z score; training median threshold; NOT a non-inflammatory CAD subtype"))
}

focus_panel <- function(ranked,ys,k,kind,ys_fraction=.8) {
  r<-ranked[ranked%in%ys]
  if(kind=="NS")return(head(ranked,k))
  if(kind=="YS")return(head(r,k))
  # Fix the supervised allocation before seeing validation performance.
  nys<-ceiling(k*ys_fraction)
  if(length(r)<nys)return(character())
  unique(c(head(r,nys),head(setdiff(ranked,ys),k-nys)))
}

focus_contrasts <- function(metrics,boot,budgets) {
  empty<-tibble(stratum=character(),landmark=numeric(),horizon=numeric(),
    AUC_a=numeric(),AUC_b=numeric(),model=character(),reference=character(),
    contrast=character(),delta_AUC=numeric(),inference=character())
  if(!nrow(metrics))return(empty)
  specs<-list()
  for(k in budgets) {
    for(cohort in c("Yin","YinYang"))for(kind in c("YS","YSplus","YSbalanced"))
      specs[[length(specs)+1L]]<-c(paste(kind,cohort,k,sep="_"),paste("NS",k,sep="_"),"Supervision vs NS")
    for(kind in c("YS","YSplus","YSbalanced"))specs[[length(specs)+1L]]<-c(paste(kind,"YinYang",k,sep="_"),paste(kind,"Yin",k,sep="_"),"Added Yang for proxy learning")
  }
  q<-function(x,p)if(sum(is.finite(x))>=20)as.numeric(quantile(x,p,na.rm=TRUE))else NA_real_
  result<-map_dfr(specs,function(s) {
    a<-metrics|>filter(model==s[1]);b<-metrics|>filter(model==s[2])
    if(!nrow(a)||!nrow(b))return(tibble())
    z<-inner_join(a|>select(stratum,landmark,horizon,AUC_a=AUC),b|>select(stratum,landmark,horizon,AUC_b=AUC),by=c("stratum","landmark","horizon"))
    if(nrow(boot)) {
      ba<-boot|>filter(model==s[1]);bb<-boot|>filter(model==s[2])
      zz<-inner_join(ba,bb,by=c("replicate","stratum","landmark","horizon"),suffix=c("_a","_b"))|>
        mutate(delta=AUC_a-AUC_b)|>group_by(stratum,landmark,horizon)|>
        summarise(delta_lo=q(delta,.025),delta_hi=q(delta,.975),valid=sum(is.finite(delta)),.groups="drop")
      z<-left_join(z,zz,by=c("stratum","landmark","horizon"))
    }
    z|>mutate(model=s[1],reference=s[2],contrast=s[3],delta_AUC=AUC_a-AUC_b,
      inference="Exploratory paired validation bootstrap; multiple budgets/strata, no external confirmation")
  })
  bind_rows(empty,result)
}

focus_plot <- function(metrics,contrasts,proxy,pillars,rd) {
  if(!nrow(metrics))return(invisible(NULL))
  a<-metrics|>filter(stratum=="All",landmark==0,budget>0)|>ggplot(aes(budget,AUC,color=paradigm))+
    geom_line()+geom_point()+scale_x_continuous(breaks=sort(unique(metrics$budget[metrics$budget>0])))+
    labs(title="a. Same assay budget and incident validation cohort",x="Measured proteins/metabolites",y="10-year IPCW AUC",color=NULL)+theme_5c(9)
  clinical_auc<-metrics$AUC[metrics$model=="Clinical"&metrics$stratum=="All"&metrics$landmark==0]
  if(length(clinical_auc))a<-a+geom_hline(yintercept=clinical_auc[1],linetype=2,color="grey45")+
    labs(subtitle=sprintf("Dashed line: clinical covariates alone (AUC %.3f)",clinical_auc[1]))
  # Older saved results can contain a zero-column tibble when no pair exists.
  z<-if(nrow(contrasts))contrasts|>filter(stratum=="All",landmark==0) else contrasts
  if(nrow(z)) {
  b<-ggplot(z,aes(delta_AUC,reorder(model,delta_AUC),color=contrast))+geom_vline(xintercept=0,linetype=2)+
    geom_point()+labs(title="b. Paired comparison with NS and with Yin",x="Change in AUC",y=NULL,color=NULL)+theme_5c(9)
  if(all(c("delta_lo","delta_hi")%in%names(z)))b<-b+geom_errorbar(aes(xmin=delta_lo,xmax=delta_hi),orientation="y",width=.2)
  } else {
    b<-ggplot()+annotate("text",x=0,y=0,label="No eligible paired comparisons\nSee design_status and fit_diagnostics")+
      labs(title="b. Paired comparison with NS and with Yin")+theme_void()
  }
  c<-proxy|>ggplot(aes(component,model,fill=delta_R2))+geom_tile()+
    scale_fill_gradient2(low="#A64E59",mid="white",high="#24748A",midpoint=0)+
    labs(title="c. LE8 reconstruction in held-out participants",x=NULL,y=NULL,fill="Incremental R²")+
    theme_5c(8)+theme(axis.text.x=element_text(angle=40,hjust=1))
  d<-pillars|>ggplot(aes(component,n,fill=cohort))+geom_col(position="dodge")+
    labs(title="d. Replicated proxies by pillar, learned in training",x=NULL,y="Eligible features",fill=NULL)+theme_5c(9)+theme(axis.text.x=element_text(angle=40,hjust=1))
  p<-(a|b)/(c|d)+plot_annotation(caption="Yang is used for LE8-proxy discovery only. All risk coefficients are fitted in incident training participants. No mechanistic subtype is inferred.")
  le8_queue_figure(p,file.path(rd,"c4_focus.Fig1.equal_budget.png"),19,13,220)
  le8_flush_figures(rd)
}

focus_write_outputs<-function(tables,rd) {
  for(nm in setdiff(names(tables),"PGS"))write_raw_csv(tables[[nm]],paste0("c4.focus.",nm,".csv"),rd)
  focus_plot(tables$metrics,tables$contrasts,tables$proxy_accuracy,tables$pillar_counts,rd)
  sheets<-tables[setdiff(names(tables),"PGS")]
  if(is.list(tables$PGS))for(nm in names(tables$PGS))sheets[[paste0("PGS_",nm)]]<-tables$PGS[[nm]]
  openxlsx::write.xlsx(sheets,file.path(rd,"c4.focus.xlsx"),overwrite=TRUE)
}

run_c4_focus <- function(layer) {
  outdir<-if(layer=="protein")out.prot else out.met
  rd<-le8_job_dir(outdir,"c4_focus");dir.create(rd,recursive=TRUE,showWarnings=FALSE)
  budget_values<-suppressWarnings(as.numeric(le8_csv_env("C4_FOCUS_BUDGETS",paste(C4_FOCUS_BUDGETS,collapse=","))))
  if(any(!is.finite(budget_values))||any(budget_values!=floor(budget_values)))stop("Assay budgets must be whole numbers")
  budgets<-sort(unique(as.integer(budget_values)))
  solver<-Sys.getenv("C4_FOCUS_SOLVER",unset="ridge")
  if(!solver%in%c("cox","ridge"))stop("C4_FOCUS_SOLVER must be cox or ridge")
  B<-as.integer(le8_num_env("C4_FOCUS_BOOT",200));fraction<-le8_num_env("C4_FOCUS_YS_FRACTION",.8)
  if(anyNA(budgets)||!length(budgets)||any(budgets<1|budgets>500)||B<0||fraction<=0||fraction>1)stop("Invalid C4_FOCUS settings (budgets: 1..500)")
  biomfile<-file.path(indir,"Rdata",if(layer=="protein")"prot.rds"else"met.rds")
  inputs<-c(file.path(indir,"Rdata/all.rds"),biomfile)
  fi<-file.info(inputs)
  signature<-le8_hash_object(list(inputs=inputs,size=fi$size,mtime=fi$mtime,options=le8_analysis_options(),
    code=tools::md5sum(file.path(fdir,c("c4_focus.R","c4_focus_extensions.R","c5_revision_validation.R","c4_pgs_bridge.R","c1_pgs.R","comm.f.R","c0_revision_core.R"))),
    settings=Sys.getenv()[grepl("^(C4_FOCUS|C5_CONNECTION|C4_PGS|C1_PGS|DATE_FOLLOW_END)",names(Sys.getenv()))],budgets=budgets,B=B,seed=SEED,solver=solver))
  cache<-file.path(rd,"c4.res.rds")
  if(!LE8_REPLACE&&file.exists(cache)) {
    old<-readRDS(cache)
    if(identical(old$signature,signature)){
      message("C4 focus: matching completed result reused; regenerate presentation files")
      focus_write_outputs(old$tables,rd);return(invisible(old))
    }
  }
  stage<-function(name,expr) {
    f<-file.path(rd,paste0(name,".rds"));old<-if(file.exists(f))readRDS(f)else NULL
    key<-focus_stage_signature(name,layer,inputs)
    if(!LE8_REPLACE&&identical(old$signature,key))return(old$value)
    value<-force(expr);saveRDS(list(signature=key,value=value),f);value
  }
  message("C4 focus: loading current phenotypes and ",layer)
  biom<-if(layer=="protein")read_prot()else read_met();features<-setdiff(names(biom),"eid")
  clinical<-if(le8_custom_adjustment())le8_custom_covars else unique(c(vars.basic,"smoke.pts","sbp","hba1c_ngsp","bmi","nonhdl"))
  requested_components<-le8_csv_env("C4_FOCUS_COMPONENTS",paste(vars.le8,collapse=","))
  if(length(requested_components)<1||length(requested_components)>8)stop("Request 1–8 candidate LE8 domains; unsupported domains need not contribute assays")
  need<-unique(c("eid","ethnic.c",clinical,vars.basic,vars.le8,requested_components,"birth_date","date_attend","date_lost","date_death",le8_y_date()))
  ph<-read_all(need)|>filter_analysis_cohort()|>make_outcome(Y)
  ph$eid<-as.character(ph$eid);biom$eid<-as.character(biom$eid)
  if(anyDuplicated(ph$eid)||anyDuplicated(biom$eid))stop("Duplicate phenotype/omic eid")
  dat<-inner_join(ph,biom,by="eid");rm(ph,biom);invisible(gc())
  tvar<-paste0(Y,".t2e");evar<-paste0(Y,".Yt2e");bvar<-paste0(Y,".b2e")
  yang<-dat[is.finite(dat[[bvar]])&dat[[bvar]]<=0,,drop=FALSE]
  yin<-dat[is.finite(dat[[tvar]])&dat[[tvar]]>0&dat[[evar]]%in%c(0,1),,drop=FALSE]
  if(length(intersect(yin$eid,yang$eid)))stop("Yin/Yang sets overlap")
  fold<-focus_split(yin,evar);train<-yin[fold=="training",,drop=FALSE];test<-yin[fold=="validation",,drop=FALSE]
  # The same Yin donor remains in the same discovery/replication half when
  # Yang donors are added; resplitting Yin would confound the comparison.
  set.seed(SEED+411)
  train$.le8_proxy_half<-sample(rep(1:2,length.out=nrow(train)))
  yang$.le8_proxy_half<-sample(rep(1:2,length.out=nrow(yang)))
  rm(dat,yin);invisible(gc())
  if(nrow(train)<1000||nrow(test)<100)stop("Insufficient omic participants")
  if(length(setdiff(clinical,names(train))))stop("Clinical covariates missing")
  missing_domains<-setdiff(requested_components,names(train))
  write_raw_csv(tibble(component=requested_components,status=ifelse(requested_components%in%missing_domains,"unavailable","candidate; support determined within training")),"c4.focus.domain_availability.csv",rd)
  components<-intersect(requested_components,names(train));basic<-intersect(vars.basic,names(train))
  message("C4 focus: training=",nrow(train),", validation=",nrow(test),", Yang=",nrow(yang))
  data.table::fwrite(bind_rows(tibble(eid=train$eid,role="incident_training"),tibble(eid=test$eid,role="incident_validation"),tibble(eid=yang$eid,role="Yang_proxy_training")),file.path(rd,"c4.focus.roles.csv.gz"),compress="gzip")
  screen<-stage("training_screen",cox_scan(train,features,clinical,Y,time_var=tvar,event_var=evar))
  ranked<-screen|>filter(is.finite(p.value))|>arrange(p.value,term)|>pull(term)
  memberships<-list();sets<-list()
  for(cohort in c("Yin","YinYang")) {
    message("C4 focus: learn ",cohort," LE8 proxies")
    subdir<-file.path(rd,cohort);dir.create(subdir,showWarnings=FALSE)
    learn<-if(cohort=="Yin")train else bind_rows(train,yang)
    sets[[cohort]]<-stage(paste0("proxy_",cohort),le8_training_connection_set(learn,features,components,basic,subdir))
    f<-file.path(subdir,"c5.connection_membership_training_only.csv")
    if(file.exists(f))memberships[[cohort]]<-as_tibble(fread(f))|>mutate(cohort=cohort)
  }
  membership<-bind_rows(memberships)
  if(!ncol(membership))membership<-tibble(feature=character(),component=character(),cohort=character(),selected=logical(),r1=numeric(),r2=numeric(),FDR1=numeric(),FDR2=numeric(),specificity=numeric())
  pillars<-membership|>filter(selected)|>count(cohort,component)|>
    complete(cohort=c("Yin","YinYang"),component=components,fill=list(n=0L))
  designs<-list(list(name="Clinical",features=character(),budget=0,paradigm="Clinical",cohort="Yin"))
  design_status<-list()
  for(k in budgets)for(kind in c("NS","YS","YSplus"))for(cohort in if(kind=="NS")"Yin"else c("Yin","YinYang")) {
    fs<-focus_panel(ranked,sets[[cohort]],k,kind,fraction)
    nm<-if(kind=="NS")paste(kind,k,sep="_")else paste(kind,cohort,k,sep="_")
    design_status[[length(design_status)+1L]]<-tibble(model=nm,budget=k,available=length(fs),status=if(length(fs)==k)"ready"else"insufficient eligible assays")
    if(length(fs)==k)designs[[length(designs)+1L]]<-list(name=nm,features=fs,budget=k,paradigm=if(kind=="NS")kind else paste(kind,cohort),cohort=cohort)
  }
  # A separate LE8-first arm: balance available pillars, rank by replicated
  # proxy strength only. Do not force unvalidated sleep/activity proxies.
  for(k in budgets)for(cohort in c("Yin","YinYang")) {
    fs<-focus_balanced_panel(membership|>filter(.data$cohort==.env$cohort),k,components)
    nm<-paste("YSbalanced",cohort,k,sep="_")
    design_status[[length(design_status)+1L]]<-tibble(model=nm,budget=k,available=length(fs),
      status=if(length(fs)==k)"ready"else"insufficient replicated proxies")
    if(length(fs)==k)designs[[length(designs)+1L]]<-list(name=nm,features=fs,budget=k,
      paradigm=paste("YSbalanced",cohort),cohort=cohort)
  }
  panel_audit<-focus_panel_audit(designs,membership,components)
  inflammatory<-focus_inflammation(train,test,features)
  stratum<-list(All=seq_len(nrow(test)))
  if(layer=="protein")for(g in c("Low baseline inflammation","High baseline inflammation"))stratum[[g]]<-which(inflammatory$group==g)
  checkpoint_file<-file.path(rd,"c4.focus.prediction_checkpoint.rds")
  checkpoint<-if(file.exists(checkpoint_file))readRDS(checkpoint_file)else NULL
  if(!LE8_REPLACE&&identical(checkpoint$signature,signature)) {
    message("C4 focus: reuse completed prediction checkpoint")
    tables<-checkpoint$tables;met<-tables$metrics;contrasts<-tables$contrasts;pr<-tables$proxy_accuracy
  } else {
  metrics<-boot<-proxy<-members<-calibration<-decision<-coefficients<-preprocessing<-diagnostics<-baseline_hazards<-list();mi<-0L
  for(ds in designs) {
    message("C4 focus: fit ",ds$name)
    obj<-le8_fit_budget_model(train,test,clinical,ds$features,tvar,evar,solver=solver,seed=SEED+27)
    diagnostics[[length(diagnostics)+1L]]<-tibble(model=ds$name,status=obj$status,
      fit_method=obj$fit_method%||%solver,condition_number=obj$condition_number%||%NA_real_,
      lambda=obj$lambda%||%NA_real_,warnings=obj$warnings%||%"")
    members[[length(members)+1L]]<-tibble(model=ds$name,feature=if(length(ds$features))ds$features else NA_character_,budget=ds$budget,status=obj$status)
    if(obj$status!="ok")next
    coefficients[[length(coefficients)+1L]]<-obj$coefficient|>mutate(model=ds$name,lp_center=obj$lp_center)
    preprocessing[[length(preprocessing)+1L]]<-obj$preprocess|>mutate(model=ds$name)
    baseline_hazards[[length(baseline_hazards)+1L]]<-as_tibble(obj$baseline_hazard)|>mutate(model=ds$name,lp_center=obj$lp_center)
    # Common reconstruction training cohort isolates panel selection from the
    # effect of fitting reconstruction coefficients in a different population.
    proxy[[length(proxy)+1L]]<-focus_proxy_accuracy(train,test,ds$features,components,basic,ds$name)
    for(g in names(stratum))for(L in c(0,2,5)) {
      ii<-stratum[[g]];ii<-ii[test[[tvar]][ii]>L]
      if(length(ii)<100)next
      # Conditional risk from landmark L to baseline year 10, with a frozen fit.
      bh<-obj$baseline_hazard;H<-function(t){i<-which(bh$time<=t);if(length(i))bh$hazard[max(i)]else 0}
      risk<--expm1(-max(0,H(10)-H(L))*exp(pmin(30,obj$lp[ii])))
      ev<-le8_evaluate_risk(test[[tvar]][ii]-L,test[[evar]][ii],risk,10-L,ds$name,ds$budget,ds$paradigm,
        B=B,seed=SEED+1000*match(g,names(stratum))+as.integer(L*10))
      mi<-mi+1L;metrics[[mi]]<-ev$metrics|>mutate(stratum=g,landmark=L,actual_assays=obj$N_selected,fit_method=obj$fit_method)
      if(nrow(ev$bootstrap))boot[[mi]]<-ev$bootstrap|>mutate(stratum=g,landmark=L)
      if(nrow(ev$calibration))calibration[[mi]]<-ev$calibration|>mutate(stratum=g,landmark=L)
      if(nrow(ev$decision))decision[[mi]]<-ev$decision|>mutate(stratum=g,landmark=L)
    }
  }
  met<-bind_rows(metrics);bo<-bind_rows(boot);pr<-bind_rows(proxy)
  contrasts<-focus_contrasts(met,bo,budgets)
  heterogeneity<-focus_heterogeneity(contrasts,bo)
  tables<-list(metrics=met,contrasts=contrasts,heterogeneity=heterogeneity,proxy_accuracy=pr,panel_members=bind_rows(members),
    pillar_counts=pillars,membership=membership,inflammation_definition=inflammatory$audit,
    design_status=bind_rows(design_status),training_screen=screen,
    calibration=bind_rows(calibration),decision=bind_rows(decision),
    model_coefficients=bind_rows(coefficients),preprocessing=bind_rows(preprocessing),baseline_hazards=bind_rows(baseline_hazards),
    fit_diagnostics=bind_rows(diagnostics),panel_overlap=panel_audit$panel_overlap,
    panel_coverage=panel_audit$panel_coverage,
    design=tibble(item=c("training","YinYang","test","YSP allocation","inference"),value=c(
      paste(nrow(train),"incident participants"),paste(nrow(yang),"additional prevalent donors for proxy learning only"),
      paste(nrow(test),"incident participants; shared 80/20 split algorithm with C5"),
      paste0("ceil(budget * ",fraction,") YS; remainder disease-ranked outside YS"),
      "Post hoc hypothesis development after inspecting earlier results; frozen-fit internal bootstrap, requires new external validation")))
  for(nm in names(tables))write_raw_csv(tables[[nm]],paste0("c4.focus.",nm,".csv"),rd)
  data.table::fwrite(bo,file.path(rd,"c4.focus.bootstrap.csv.gz"),compress="gzip")
  saveRDS(list(signature=signature,tables=tables),checkpoint_file)
  }
  if(truthy(Sys.getenv("C4_FOCUS_PGS",unset="TRUE"))) {
    # Save bridge results under the established C4 directory as well as this run.
    pd<-bind_rows(train,yang);set.seed(SEED+811);half<-sample(rep(c("discovery","replication"),length.out=nrow(pd)))
    mm<-membership|>filter(cohort=="YinYang")|>transmute(feature,primary_component=component,strict_YS=selected)
    prd<-le8_job_dir(outdir,"c4_connect");dir.create(prd,recursive=TRUE,showWarnings=FALSE)
    tables$PGS<-le8_pgs_bridge(pd,features,half,basic,screen,mm,layer,prd)
  }
  out<-list(signature=signature,meta=module_meta(layer,extra=list(status="ok")),tables=tables)
  saveRDS(out,cache);focus_write_outputs(tables,rd);invisible(out)
}

if(!truthy(Sys.getenv("LE8_FOCUS_FUNCTIONS_ONLY",unset="FALSE"))) {
  if(prot_DO)le8_stage("C4/focus/protein", run_c4_focus("protein"))
  if(met_DO)le8_stage("C4/focus/metabolite", run_c4_focus("metabolite"))
}
