# C5 v2: one joint incident cohort, one validation scheme, any prespecified Y.
fdir<-Sys.getenv("LE8_FDIR",file.path(Sys.getenv("DIRSCRIPT"),"f"))
source(file.path(fdir,"comm.f.R"))
source(file.path(fdir,"c5_revision_validation.R"))
source(file.path(fdir,"c5_joint_helpers.R"))
source(file.path(fdir,"c5_joint_outputs.R"))
source(file.path(fdir,"c5_joint_extensions.R"))
source(file.path(fdir,"c5_factorial.R"))
source(file.path(fdir,"c5_joint_storage.R"))
LE8_JOB<-"c5_consolidate"
if(!requireNamespace("glmnet",quietly=TRUE))stop("C5 requires glmnet for cross-fitted proxy learning")
risk_solver<-Sys.getenv("C5_RISK_SOLVER","ridge")
if(!risk_solver%in%c("cox","ridge"))stop("C5_RISK_SOLVER must be cox or ridge")
if(!prot_DO||!met_DO)stop("Joint C5 requires --biom prot,met; all comparisons use their common cohort")
if(nzchar(Sys.getenv("C5_INC_PROT",""))||nzchar(Sys.getenv("C5_INC_MET","")))stop("Outcome-derived include lists cannot be reused for joint validation. Supply full assay matrices.")
K<-as.integer(Sys.getenv("C5_OUTER_FOLDS","5"));B<-as.integer(Sys.getenv("C5_BOOT","200"))
budgets<-as.numeric(c5_csv("C5_ASSAY_BUDGETS","5,10,50"))
landmarks<-as.numeric(c5_csv("C5_LANDMARKS","0,2,5"));end<-as.numeric(Sys.getenv("C5_END_YEARS","10"))
if(!is.finite(K)||K<3||!is.finite(B)||B<0||any(!is.finite(budgets)|budgets<2|budgets>100|budgets!=floor(budgets))||
  any(!is.finite(landmarks)|landmarks<0|landmarks>=end))stop("Invalid fold, budget or landmark settings")
budgets<-sort(unique(budgets));landmarks<-sort(unique(landmarks))
primary_budget<-as.integer(Sys.getenv("C5_PRIMARY_BUDGET","10"))
primary_landmark<-as.numeric(Sys.getenv("C5_PRIMARY_LANDMARK","5"))
if(!primary_budget%in%budgets||!primary_landmark%in%landmarks)stop("Prespecified primary budget/landmark must be included")
if(truthy(Sys.getenv("C5_JOINT_SUMMARY_ONLY","FALSE"))){
  if(LE8_REPLACE)stop("C5_JOINT_SUMMARY_ONLY cannot be combined with --replace TRUE")
  c5_resume_joint_outputs(out.base,Y,B,SEED,primary_budget,primary_landmark)
  quit(save="no",status=0)
}
completed<-c5_completed_joint_outputs(out.base,Y,B,SEED,primary_budget,primary_landmark,replace=LE8_REPLACE)
if(!is.null(completed))quit(save="no",status=0)
rm(completed)
clinical<-c5_csv("C5_CLINICAL_VARS",paste(unique(c(vars.basic,"smoke.pts","sbp","hba1c_ngsp","bmi","nonhdl")),collapse=","))
basic<-vars.basic
targets<-c5_csv("C5_DOMAINS","bmi,nonhdl,hba1c_ngsp,sbp,smoke.pts,diet.pts,pa.pts,sleep.pts")
proxy_targets<-c5_csv("C5_PROXY_TARGETS","bmi,nonhdl")
if(length(targets)<1||length(targets)>8||length(proxy_targets)<1||length(proxy_targets)>8)stop("Choose 1–8 candidate domains / proxy targets")
if(primary_budget<length(proxy_targets))stop("Primary assay budget must cover requested proxy targets")
if(!all(c("bmi","nonhdl")%in%clinical))stop("Continuous BMI and non-HDL are required clinical comparators")
if(any(grepl("(^fod_|^date_|[.]Yt2e$|[.]t2e$|[.]b2e$|^drug[.]|^dm[.]|^htn[.])",clinical)))stop("Unapproved diagnosis/medication/outcome-derived predictor; use explicit baseline fields")
if(any(grepl("[.](pgs|prs)$|^disease_PRS$",c(clinical,basic,targets,proxy_targets),ignore.case=TRUE)))
  stop("Keep disease PRS / biomarker PGS out of the clinical and LE8-domain comparator fields")
groupcol<-Sys.getenv("C5_GROUP_COLUMN","")
ph<-read_all(unique(c("eid","ethnic.c",clinical,basic,targets,proxy_targets,"birth_date","date_attend","date_lost","date_death",le8_y_date(),groupcol,paste0(Y,".pgs"))))|>
  filter_analysis_cohort()|>make_outcome(Y)
if(length(setdiff(c(clinical,basic),names(ph))))stop("Missing clinical baseline fields")
prs<-c5_prs_input(Y,ph);has_prs<-isTRUE(prs$enabled)
prs$manifest$file[prs$manifest$file=="all.rds"]<-file.path(indir,"Rdata/all.rds")
# First count the common measured cohort. Never impute a missing disease PRS.
joined<-c5_join_omics(ph,read_prot(),read_met(),data.frame(eid=as.character(ph$eid)))
prs_flow<-data.frame(stage="clinical_prot_met_intersection",N=nrow(joined$data))
if(has_prs){
  joined$data<-merge(joined$data,prs$data,by="eid",all.x=TRUE,sort=FALSE)
  prs_flow<-rbind(prs_flow,data.frame(stage="finite_endpoint_PRS_in_intersection",N=sum(is.finite(joined$data$disease_PRS))))
}
pgs<-c5_load_biomarker_pgs(joined$data,joined$map)
dat<-pgs$data;map<-pgs$map;pgs_inputs<-pgs$inputs;pgs_status<-pgs$status
rm(ph,joined,pgs);invisible(gc())
dat$.time<-dat[[paste0(Y,".t2e")]];dat$.event<-dat[[paste0(Y,".Yt2e")]]
dat<-c5_endpoint_coverage(dat,Y)
dat$.baseline_disease<-is.finite(dat[[paste0(Y,".b2e")]])&dat[[paste0(Y,".b2e")]]<=0
if(nzchar(groupcol)&&!groupcol%in%names(dat))stop("Group column missing")
dat$.group<-if(nzchar(groupcol))as.character(dat[[groupcol]])else dat$eid
yang<-dat[dat$.baseline_disease,,drop=FALSE]
yin<-dat[!dat$.baseline_disease&is.finite(dat$.time)&dat$.time>dat$.entry&dat$.event%in%0:1,,drop=FALSE]
prs_flow<-rbind(prs_flow,data.frame(stage="incident_measured_cohort_before_PRS",N=nrow(yin)))
if(has_prs){
  eligible<-is.finite(yin$disease_PRS)
  if(sum(eligible)>=1000&&sum(yin$.event[eligible])>=100&&sd(yin$disease_PRS[eligible])>0){
    yin<-yin[eligible,,drop=FALSE]
    yang<-yang[is.finite(yang$disease_PRS),,drop=FALSE]
    prs_flow<-rbind(prs_flow,data.frame(stage="common_incident_cohort_for_ALL_models_with_and_without_PRS",N=nrow(yin)))
  }else{
    has_prs<-FALSE;prs$enabled<-FALSE;prs$manifest$enabled<-FALSE
    prs$manifest$status<-paste(prs$manifest$status,"; insufficient common incident PRS cohort, non-PRS suite continues")
  }
}
write.csv(prs_flow,file.path(out.base,"c5.prs_cohort_flow.csv"),row.names=FALSE)
if(length(intersect(yin$eid,yang$eid)))stop("Yin/Yang overlap")
if(nrow(yin)<1000||sum(yin$.event)<100)stop("Insufficient common incident cohort; do not substitute a different cohort per modality")
yin$.fold<-c5_group_folds(yin$.group,K,SEED)
private<-file.path(out.base,"_c5_private");dir.create(private,recursive=TRUE,showWarnings=FALSE)
cache<-file.path(out.base,"_c5_cache");dir.create(cache,recursive=TRUE,showWarnings=FALSE)
writeLines(c("*","!.gitignore"),file.path(private,".gitignore"))
writeLines(c("*","!.gitignore"),file.path(cache,".gitignore"))
input_files<-c(file.path(indir,"Rdata",c("all.rds","prot.rds","met.rds")),prs$manifest$file,Sys.getenv("C5_PRS_MANIFEST"),Sys.getenv("LE8_ENDPOINT_MANIFEST"),pgs_inputs)
fi<-file.info(input_files)
signature<-le8_hash_object(list(version="2026-09-16.factorial-prs",Y=Y,inputs=data.frame(path=input_files,size=fi$size,mtime=as.character(fi$mtime)),
  code=tools::md5sum(file.path(fdir,c("c5_joint.R","c5_factorial.R","c5_factorial_outputs.R","c5_factorial_plots.R","c5_joint_extensions.R","c5_panels_final.R","c5_joint_helpers.R","c1_pgs.R","c5_joint_outputs.R","c5_joint_storage.R","c5_revision_validation.R","c0_baseline.R","c0_revision_core.R","comm.f.R"))),
  options=le8_analysis_options(),settings=Sys.getenv()[grepl("^(C5_|DATE_FOLLOW_END)",names(Sys.getenv()))],
  ids=yin$eid,fold=yin$.fold,clinical=clinical,budgets=budgets,landmarks=landmarks,end=end,seed=SEED))
write.csv(prs$manifest,file.path(out.base,"c5.prs_provenance.csv"),row.names=FALSE)
factorial_design<-c5_factorial_design(has_prs)
write.csv(factorial_design,file.path(out.base,"c5.factorial_design.csv"),row.names=FALSE)
write.csv(data.frame(Y=Y,N_joint=nrow(dat),N_Yin=nrow(yin),N_Yang=nrow(yang),incident_events=sum(yin$.event),
  outer_folds=K,grouping=if(nzchar(groupcol))groupcol else "eid; relatedness not controlled",signature=signature),file.path(out.base,"c5.cohort.csv"),row.names=FALSE)
data.table::fwrite(yin[,c("eid",".fold",".group")],file.path(private,"roles.csv.gz"),compress="gzip")
write.csv(map,file.path(out.base,"c5.assay_inventory.csv"),row.names=FALSE)
writeLines(pgs_status,file.path(out.base,"c5.omic_PGS_status.txt"))
# yin/yang already own the needed rows; the full joined matrix is no longer used.
rm(dat);invisible(gc())
drop_assays<-toupper(c5_csv("C5_EXCLUDE_ASSAYS","GDF15,NTPROBNP,NPPB"))
features<-map$feature
arms<-list(all_assays=features,omit_GDF15_natriuretic=map$feature[!toupper(map$assay)%in%drop_assays])
if(truthy(Sys.getenv("C5_INCLUDE_LIPID_ABLATION","TRUE"))){
  direct_lipid<-toupper(c5_csv("C5_DIRECT_LIPID_ASSAYS","Non_HDL_C,NonHDL_C,Total_C,HDL_C"))
  arms$omit_direct_lipid_reconstruction<-map$feature[!(map$layer=="met"&toupper(map$assay)%in%direct_lipid)]
}
write.csv(data.frame(arm=names(arms),available_assays=lengths(arms)),file.path(out.base,"c5.ablation_inventory.csv"),row.names=FALSE)
all_results<-character()
landmark_availability<-map_dfr(landmarks,function(L){
  d<-c5_landmark(yin,L,end);ne<-vapply(seq_len(K),function(fd)sum(d$.fold==fd&d$.event==1&d$.time<=end-L),numeric(1))
  ok<-nrow(d)>=1000&&all(ne>=10)&&all(sum(ne)-ne>=50)
  tibble(landmark=L,N=nrow(d),events=sum(ne),min_test_events=min(ne),available=ok,
    status=if(ok)"estimable"else"insufficient covered event-free participants/events; no substitute landmark chosen")
})
write.csv(landmark_availability,file.path(out.base,"c5.landmark_availability.csv"),row.names=FALSE)
landmarks_to_run<-landmark_availability$landmark[landmark_availability$available]
if(!length(landmarks_to_run))stop("No prespecified landmark is estimable; see c5.landmark_availability.csv")
for(L in landmarks_to_run)for(fd in seq_len(K)){
  checkpoint<-file.path(cache,paste0("L",L,".fold",fd,".rds"))
  old<-if(file.exists(checkpoint)&&!LE8_REPLACE)tryCatch(readRDS(checkpoint),error=function(e)NULL)else NULL
  reusable<-identical(old$signature,signature)&&is.list(old$result)
  if(reusable)c5_prune_fold_models(checkpoint,signature,object=old)
  rm(old);invisible(gc())
  if(reusable){
    message("C5: reuse landmark ",L,", outer fold ",fd,"/",K)
    all_results<-c(all_results,checkpoint);next
  }
  tr<-c5_landmark(yin[yin$.fold!=fd,,drop=FALSE],L,end)
  te<-c5_landmark(yin[yin$.fold==fd,,drop=FALSE],L,end)
  tr$.event<-as.integer(tr$.event==1&tr$.time<=end-L)
  tr$.time<-pmin(tr$.time,end-L)
  yy<-yang[!yang$.group%in%yin$.group[yin$.fold==fd],,drop=FALSE]
  if(length(intersect(tr$.group,te$.group))||length(intersect(yy$.group,te$.group)))stop("Group leakage")
  if(sum(tr$.event)<50||sum(te$.event)<10)stop("Too few events at landmark/fold; reduce prespecified model complexity or obtain more data")
  message(Y,": joint C5 landmark ",L,", outer fold ",fd,"/",K,"; events ",sum(tr$.event)," / ",sum(te$.event))
  # Identical screening reference for all panels; screening fits use outer train only.
  screen<-cox_scan(tr,features,clinical,Y,time_var=".time",event_var=".event")
  ranked<-screen|>filter(is.finite(p.value))|>arrange(p.value,term)|>pull(term)
  genetic_screen<-tibble();genetic_ranked<-character()
  if("pgs_feature"%in%names(map)){
    gs<-map$pgs_feature[!is.na(map$pgs_feature)]
    if(length(gs)){
      genetic_screen<-cox_scan(tr,gs,clinical,Y,time_var=".time",event_var=".event")
      genetic_ranked<-genetic_screen|>filter(is.finite(p.value))|>arrange(p.value,term)|>pull(term)
    }
  }
  preds<-members<-diagnostics<-coefficients<-proxy_validation<-proxy_weights<-proxy_preprocessing<-domain_tables<-list()
  risk_preprocessing<-baseline_hazards<-list()
  pgs_reconstruction<-pgs_coefficients<-cross_omic<-domain_status<-list()
  marker<-c5_marker_strata(tr,te,map)
  inflam<-c5_inflammation_score(tr,te,map)
  fit_cache<-new.env(parent=emptyenv())
  addfit<-function(model,k,arm,cv,ff,train=tr,test=te){
    if(!has_prs&&"disease_PRS"%in%c(cv,ff))return(invisible(NULL))
    cacheable<-missing(train)&&missing(test)
    ck<-paste(paste(cv,collapse="|"),paste(ff,collapse="|"),sep="::")
    obj<-if(cacheable&&exists(ck,fit_cache,inherits=FALSE))get(ck,fit_cache)else
      le8_fit_budget_model(train,test,cv,ff,".time",".event",solver=risk_solver,seed=SEED+fd)
    if(cacheable)assign(ck,obj,fit_cache)
    transformed<-ff[grepl("^(genpart__|remainder__)",ff)]
    inherited_dependencies<-if(length(transformed))paste0("pgs__",sub("^(genpart__|remainder__)","",transformed))else character()
    measured_dependencies<-unique(c(ff[ff%in%map$feature],sub("^remainder__","",ff[startsWith(ff,"remainder__")])))
    id<-length(diagnostics)+1L
    diagnostics[[id]]<<-tibble(model,budget=k,arm,fold=fd,landmark=L,status=obj$status,
      requested_predictors=length(ff),fitted_predictors=obj$N_selected%||%NA_integer_,
      direct_measured_assays=length(measured_dependencies),
      disease_PRS_predictors=as.integer("disease_PRS"%in%c(cv,ff)),
      genetic_score_predictors=length(unique(c(ff[startsWith(ff,"pgs__")],inherited_dependencies))))
    if(length(ff))members[[length(members)+1L]]<<-tibble(model,budget=k,arm,fold=fd,landmark=L,feature=ff)
    if("disease_PRS"%in%cv)members[[length(members)+1L]]<<-tibble(model,budget=k,arm,fold=fd,landmark=L,feature="disease_PRS")
    dependencies<-sub("^remainder__","",ff[startsWith(ff,"remainder__")])
    if(length(dependencies))members[[length(members)+1L]]<<-tibble(model,budget=k,arm,fold=fd,landmark=L,feature=dependencies)
    if(obj$status!="ok")return(invisible(NULL))
    risk<-le8_risk_at(obj,end-L)
    tx<-le8_prepare_prediction_matrix(train,train,unique(c(cv,ff)))
    lp_training<-drop(tx$train%*%obj$coefficient$beta)-obj$lp_center
    q75<-unname(quantile(lp_training,.75))
    preds[[length(preds)+1L]]<<-tibble(eid=test$eid,group=test$.group,fold=fd,time=test$.time,event=test$.event,
      model,budget=k,arm,landmark=L,horizon=end-L,risk=risk,lp=obj$lp,high_training_q75=obj$lp>=q75,
      inflammatory_burden=if(is.null(inflam))"unavailable"else inflam$group)
    coefficients[[length(coefficients)+1L]]<<-obj$coefficient|>mutate(model,budget=k,arm,fold=fd,landmark=L,
      lambda=obj$lambda,condition_number=obj$condition_number)
    risk_preprocessing[[length(risk_preprocessing)+1L]]<<-obj$preprocess|>mutate(model,budget=k,arm,fold=fd,landmark=L)
    bh<-obj$baseline_hazard
    baseline_hazards[[length(baseline_hazards)+1L]]<<-bh[!duplicated(bh$hazard),,drop=FALSE]|>
      mutate(model,budget=k,arm,fold=fd,landmark=L,lp_center=obj$lp_center)
    # Fold tables already retain these values; never persist per-model RDS copies.
    invisible(obj)
  }
  for(arm in names(arms)){
    pool<-arms[[arm]];rr<-ranked[ranked%in%pool]
    addfit("Clinical",0,arm,clinical,character())
    addfit("PRS_only",0,arm,"disease_PRS",character())
    addfit("Clinical_PRS",0,arm,c(clinical,"disease_PRS"),character())
    for(i in which(!factorial_design$P&!factorial_design$M&factorial_design$available)){
      ds<-factorial_design[i,];cv<-c(if(ds$C)clinical,if(ds$G)"disease_PRS")
      addfit(ds$model,0,arm,cv,character())
    }
    domains<-list()
    for(cohort in c("Yin","YinYang")){
      learn<-if(cohort=="Yin")tr else bind_rows(tr,yy)
      dom<-c5_replicated_domains(learn,pool,targets,basic,SEED+fd+411)
      domains[[cohort]]<-dom$ranks
      domain_tables[[length(domain_tables)+1L]]<-dom$all|>mutate(cohort,arm,fold=fd,landmark=L)
      domain_status[[length(domain_status)+1L]]<-dom$audit|>mutate(cohort,arm,fold=fd,landmark=L)
    }
    for(k in budgets){
      pp<-head(rr[startsWith(rr,"prot__")],k);mm<-head(rr[startsWith(rr,"met__")],k)
      # Joint NS budget is TOTAL assays: ceil(k/2) protein, floor(k/2) metabolite.
      joint<-c(head(pp,ceiling(k/2)),head(mm,floor(k/2)))
      # Fixed layer panels: adding M retains ALL k proteins (k+k assays).
      # The legacy equal-total-budget comparison below remains unchanged.
      for(i in which((factorial_design$P|factorial_design$M)&factorial_design$available)){
        ds<-factorial_design[i,];ff<-c(if(ds$P)pp,if(ds$M)mm)
        cv<-c(if(ds$C)clinical,if(ds$G)"disease_PRS")
        if((ds$P&&length(pp)!=k)||(ds$M&&length(mm)!=k))next
        addfit(ds$model,k,arm,cv,ff)
      }
      for(mode in c("Protein","Metabolite","ProtMet")){
        ff<-switch(mode,Protein=pp,Metabolite=mm,ProtMet=joint)
        if(length(ff)!=k)stop("Insufficient measured assays for requested budget")
        addfit(paste0(mode,"_only_NS"),k,arm,character(),ff)
        addfit(paste0("Clinical_",mode,"_NS"),k,arm,clinical,ff)
      }
      addfit("Clinical_Protein_sharedBudget_NS",k,arm,clinical,joint[startsWith(joint,"prot__")])
      addfit("Clinical_Metabolite_sharedBudget_NS",k,arm,clinical,joint[startsWith(joint,"met__")])
      addfit("Clinical_ProtMet_PRS_NS",k,arm,c(clinical,"disease_PRS"),joint)
      gp<-character()
      if(length(genetic_ranked)){
        eligible<-map$pgs_feature[map$feature%in%pool&!is.na(map$pgs_feature)]
        gr<-genetic_ranked[genetic_ranked%in%eligible]
        gp<-c(head(gr[startsWith(gr,"pgs__prot__")],ceiling(k/2)),head(gr[startsWith(gr,"pgs__met__")],floor(k/2)))
        if(length(gp)==k){
          addfit("Clinical_PGSselected_NS",k,arm,clinical,gp)
          addfit("Clinical_ProtMet_PGSselected_NS",k,arm,clinical,c(joint,gp))
          addfit("Clinical_ProtMet_PGSselected_PRS_NS",k,arm,c(clinical,"disease_PRS"),c(joint,gp))
          genetically_ranked_assays<-map$feature[match(gp,map$pgs_feature)]
          addfit("Clinical_GeneticSelectedMeasured_NS",k,arm,clinical,genetically_ranked_assays)
        }
      }
      if("pgs_feature"%in%names(map)){
        gg<-map$pgs_feature[match(joint,map$feature)]
        valid<-vapply(gg,function(g)!is.na(g)&&g%in%names(tr)&&sum(is.finite(tr[[g]]))>=100,logical(1))
        matched<-joint[valid];gg<-gg[valid]
        if(length(gg)>=2){
          # The reduced matched subset gets its OWN measured reference.
          addfit("MatchedPGS_only_NS",k,arm,character(),gg)
          addfit("MatchedMeasured_only_NS",k,arm,character(),matched)
          addfit("Clinical_MatchedMeasured_NS",k,arm,clinical,matched)
          addfit("Clinical_MatchedPGS_NS",k,arm,clinical,gg)
          addfit("Clinical_MatchedMeasured_PGS_NS",k,arm,clinical,c(matched,gg))
          addfit("Clinical_MatchedMeasured_PGS_PRS_NS",k,arm,c(clinical,"disease_PRS"),c(matched,gg))
          if(k==primary_budget){
            decomp<-c5_pgs_decompose(tr,te,matched,map,SEED+fd+312)
            pgs_reconstruction[[length(pgs_reconstruction)+1L]]<-decomp$audit|>mutate(fold=fd,landmark=L,arm)
            pgs_coefficients[[length(pgs_coefficients)+1L]]<-decomp$coefficient|>mutate(fold=fd,landmark=L,arm)
            if(length(decomp$genetic)>=2){
              addfit("Clinical_PGS_captured",k,arm,clinical,decomp$genetic,decomp$train,decomp$test)
              addfit("Clinical_PGS_remainder",k,arm,clinical,decomp$remainder,decomp$train,decomp$test)
              addfit("Clinical_PGS_captured_remainder",k,arm,clinical,c(decomp$genetic,decomp$remainder),decomp$train,decomp$test)
            }
          }
        }
      }
      for(cohort in c("Yin","YinYang")){
        queues<-lapply(domains[[cohort]],function(z)z$feature)
        ys<-c5_ys_panel(queues,k);plus<-c5_ys_panel(queues,k,rr)
        if(length(ys)!=k||length(plus)!=k){
          diagnostics[[length(diagnostics)+1L]]<-tibble(model=paste0("YS_",cohort),budget=k,arm,fold=fd,landmark=L,
            status="insufficient replicated domain proxies; no forced pillar or NS substitution")
          next
        }
        addfit(paste0("Clinical_ProtMet_YS_",cohort),k,arm,clinical,ys)
        addfit(paste0("Clinical_ProtMet_YSplus_",cohort),k,arm,clinical,plus)
        addfit(paste0("Clinical_ProtMet_PRS_YS_",cohort),k,arm,c(clinical,"disease_PRS"),ys)
        if(length(gp)==k)addfit(paste0("Clinical_ProtMet_YS_PGSselected_",cohort),k,arm,clinical,c(ys,gp))
        if(k==primary_budget&&cohort=="YinYang")cross_omic[[length(cross_omic)+1L]]<-
          c5_cross_omic_links(tr,te,union(joint,ys),basic)|>mutate(fold=fd,landmark=L,arm)
      }
    }
    # Domain proxy models at one prespecified budget. Report union assay count;
    # two predicted targets do not mean a two-assay clinical panel.
    pk<-primary_budget
    for(cohort in c("Yin","YinYang")){
      extra<-if(cohort=="Yin")yy[FALSE,,drop=FALSE]else yy
      for(modality in c("prot","met","joint")){
        candidates<-if(modality=="joint")pool else pool[startsWith(pool,paste0(modality,"__"))]
        trp<-tr;tep<-te;union_panel<-character();proxy_ok<-TRUE
        for(tg in proxy_targets){
          kp<-pk%/%length(proxy_targets)+as.integer(match(tg,proxy_targets)<=pk%%length(proxy_targets))
          obj<-tryCatch(c5_crossfit_proxy(tr,te,extra,candidates,tg,kp,basic,SEED+fd+111),error=function(e){
            diagnostics[[length(diagnostics)+1L]]<<-tibble(model=paste0("Proxy_",tg,"_",modality,"_",cohort),budget=pk,arm,fold=fd,landmark=L,status=conditionMessage(e))
            NULL})
          if(is.null(obj)){proxy_ok<-FALSE;next}
          nv<-paste0("proxy_",tg);trp[[nv]]<-obj$train;tep[[nv]]<-obj$test
          union_panel<-union(union_panel,obj$final$features)
          ok<-is.finite(te[[tg]]);pred<-obj$test[ok];obs<-te[[tg]][ok]
          denom<-sum((obs-mean(tr[[tg]],na.rm=TRUE))^2)
          proxy_validation[[length(proxy_validation)+1L]]<-tibble(fold=fd,landmark=L,arm,cohort,modality,target=tg,
            N=sum(ok),R2_vs_training_mean=1-sum((obs-pred)^2)/denom,RMSE=sqrt(mean((obs-pred)^2)),
            correlation=cor(obs,pred),scope="Held-out reconstruction of continuous measurement; not evidence of intervention responsiveness")
          proxy_weights[[length(proxy_weights)+1L]]<-obj$final$coefficient|>mutate(fold=fd,landmark=L,arm,cohort,modality,target=tg)
          proxy_preprocessing[[length(proxy_preprocessing)+1L]]<-obj$final$preprocess|>mutate(fold=fd,landmark=L,arm,cohort,modality,target=tg)
        }
        if(!proxy_ok)next
        for(action in c("Replace","Add")){
          if(action=="Replace"&&!all(proxy_targets%in%clinical))next
          cv<-if(action=="Replace")setdiff(clinical,proxy_targets)else clinical
          label<-paste0("Clinical_",action,"Proxy_",modality,"_",cohort)
          addfit(label,pk,arm,cv,paste0("proxy_",proxy_targets),trp,tep)
          # Correct the feature inventory: measured inputs, not two score names.
          members[[length(members)+1L]]<-tibble(model=label,budget=pk,arm,fold=fd,landmark=L,feature=union_panel)
        }
      }
    }
  }
  result<-list(predictions=bind_rows(preds),members=bind_rows(members),diagnostics=bind_rows(diagnostics),
    risk_preprocessing=bind_rows(risk_preprocessing),baseline_hazards=bind_rows(baseline_hazards),
    coefficients=bind_rows(coefficients),proxy_validation=bind_rows(proxy_validation),proxy_weights=bind_rows(proxy_weights),proxy_preprocessing=bind_rows(proxy_preprocessing),
    domain_associations=bind_rows(domain_tables),domain_status=bind_rows(domain_status),
    pgs_reconstruction=bind_rows(pgs_reconstruction),pgs_coefficients=bind_rows(pgs_coefficients),cross_omic=bind_rows(cross_omic),
    marker_groups=marker$groups|>mutate(fold=fd,landmark=L),marker_audit=marker$audit|>mutate(fold=fd,landmark=L),
    genetic_screen=genetic_screen|>mutate(fold=fd,landmark=L),
    screen=screen|>mutate(fold=fd,landmark=L))
  c5_write_joint_checkpoint(list(signature=signature,result=result),checkpoint)
  c5_prune_fold_models(checkpoint,signature,object=list(signature=signature,result=result))
  all_results<-c(all_results,checkpoint)
  rm(result,preds,members,diagnostics,coefficients,proxy_validation,proxy_weights,proxy_preprocessing,
    domain_tables,pgs_reconstruction,pgs_coefficients,cross_omic,domain_status,fit_cache,risk_preprocessing,baseline_hazards)
  rm(list=intersect(c("tr","te","yy","learn","trp","tep","decomp","obj","extra","marker","dom","domains","screen","genetic_screen"),ls()))
  invisible(gc())
}
# Summaries need only eligibility/IDs; release the assay matrices and final
# fold's working copies before reading any saved predictions.
yin<-yin[,intersect(c("eid",".time",".event",".entry"),names(yin)),drop=FALSE]
rm(list=intersect(c("dat","yang","tr","te","yy","learn","trp","tep","decomp","obj","extra","marker","dom","domains","screen","genetic_screen"),ls()))
invisible(gc())
c5_joint_outputs(all_results,yin,map,out.base,K,B,SEED,primary_budget,primary_landmark,signature)
